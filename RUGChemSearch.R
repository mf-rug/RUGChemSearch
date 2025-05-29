library(httr)
library(tidyverse)
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)
library(DT)
library(googlesheets4)

smiles_to_pubchem_url <- function(smiles, cid_list, max_display = 10000, wait_seconds = 2, substructure='identity', department=TRUE) {
  # Step 1: Start the substructure search job
  # if (substructure) {
  #   url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/JSON"
  # } else {
  #   url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/JSON"
  # }
  url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/', substructure, '/smiles/JSON')
  
  init_res <- POST(
    url,
    body = list(smiles = smiles),
    encode = "form"
  )
  
  if (!init_res$status_code %in% c(202,200)) {
    stop("Substructure search initiation failed.")
  }
  
  listkey <- content(init_res, as = "parsed")$Waiting$ListKey
  if (!is.null(listkey)) {
    cat('listkey is -->', listkey, '<--\n')
    
    # Step 2: Poll for completion
    Sys.sleep(wait_seconds)  # wait before first try
    repeat {
      poll_res <- GET(paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/",
        listkey,
        "/cids/JSON"
      ))
      
      if (poll_res$status_code == 200) {
        cat('Finished', '\n')
        break
      } else if (poll_res$status_code == 202) {
        cat('Running', '\n')
        Sys.sleep(wait_seconds)  # still processing
      } else {
        stop(paste("Unexpected error while polling for search result:", poll_res))
      }
    }
    cids_found <- content(poll_res, as = "parsed")$IdentifierList$CID %>% unlist()
  } else {
    cat('no listkey, result directly returned\n')
    result <- content(init_res, as = "parsed")
    cids_found <- result$PC_Compounds[[1]]$id$id$cid
  }
  if (is.null(cids_found)) {
    cat('No results', '\n')
    return(NULL)
  }
  cat(length(cids_found), ' results found', '\n')
  
  # Step 3: Intersect with your list
  cid_list <- as.integer(cid_list)
  if (department) {
    match_cids <- intersect(cids_found, cid_list)
  }
  else {
    match_cids <- cids_found
  }
  
  cat(length(match_cids), ' matches found', '\n')
  if (length(match_cids) == 0) return('No results')
  
  # Step 4: Create PubChem web URL (limit output)
  if (length(match_cids) > 1) {
    url <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/#query=",
      paste0(head(match_cids, max_display), collapse = ","),
      "&collection=compound"
    )   
  } else {
    url <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/compound/",
      match_cids
    )
  }
  return(url)
}


ui <- fluidPage(
  useShinyjs(),
  titlePanel("Chemicals Search"),
  sidebarLayout(
    sidebarPanel(
      tags$head(
        tags$script(src = "https://unpkg.com/kekule/dist/kekule.min.js"),
        tags$link(rel = "stylesheet", href = "https://unpkg.com/kekule/dist/themes/default/kekule.css")
      ),
      tags$h3("Draw a molecule"),
      tags$div(id = "editor", style = "width: 500px; height: 400px; border: 1px solid #ccc;"),
      verbatimTextOutput("smiles"),
      
      tags$script(HTML("
var editor;

function initKekule() {
  var container = document.getElementById('editor');
  if (!container) {
    console.error('Kekule editor container not found.');
    return;
  }

  editor = new Kekule.Editor.Composer(container);

  editor.addEventListener('change', function() {
    var mol = editor.getChemObj();
    if (mol) {
      var smiles = Kekule.IO.saveFormatData(mol, 'smi');
      if (smiles && smiles.trim() !== '') {
        Shiny.setInputValue('from_kekule', smiles, {priority: 'event'});
      }
    }
  });
}

// ensure DOM and Shiny are ready
if (window.Shiny && document.readyState === 'complete') {
  initKekule();
} else {
  document.addEventListener('DOMContentLoaded', function() {
    if (window.Shiny) {
      initKekule();
    } else {
      document.addEventListener('shiny:connected', initKekule);
    }
  });
}

")),
      br(), 
      textInput("smiles", "Draw above or directly enter SMILES", ""),
      radioGroupButtons('search', NULL, c('identity', 'similarity', 'substructure', 'superstructure', '3d similarity')),
      materialSwitch(
        'department',
        HTML('Restrict to <a href="#" onclick="Shiny.setInputValue(\'main_tabs\', \'browse\')">department chemicals</a>'),
        value = TRUE,
        status = 'primary',
        right = TRUE
      ),
      
      actionButton("submit", "Search"),
    ),
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel('search',
                 withSpinner(uiOutput("result_url"))
        ),
        tabPanel('browse',
                 br(),
                 fluidRow(materialSwitch('ordered', 'update list with chemicals from the Chemical Order Form (googlesheet)', status = 'primary', width = '100%', right = TRUE)),
                 fluidRow(hidden(materialSwitch('subset', 'Only show the last search hits', status = 'primary', value = TRUE, right = TRUE))),
                 withSpinner(DTOutput("dep_chemicals"))
        )
      )
    )
  )
)


safe_read <- possibly(function(name) {
  read_sheet(order_gsheet, sheet = name)
}, otherwise = NULL)

cas_regex <- "^\\d{2,7}-\\d{2}-\\d$"

# Function to check if a column contains mostly CAS-like entries
is_mostly_cas_column <- function(column) {
  if (!is.character(column)) column <- as.character(column)
  matches <- str_detect(column, cas_regex)
  valid <- sum(!is.na(column))
  if (valid == 0) return(FALSE)
  mean(matches, na.rm = TRUE) >= 0.5
}

# Function to check if a data frame contains at least one mostly-CAS column
has_mostly_cas_column <- function(df) {
  any(map_lgl(df, is_mostly_cas_column))
}


server <- function(input, output, session) {
  search_trigger <- reactiveVal(NULL)
  cas_hits <- reactiveVal(NULL)
  
  cid_list <- readLines("bc_pubchem.txt")
  
  dep_chems <- read_csv('biotechnology_chemicals.csv', show_col_types = FALSE, col_names = TRUE)
  pubchem <- read_csv('biochemistry_all_cas_and_pubchem_ids.csv', show_col_types = FALSE, col_names = TRUE)
  dep_chems <- dep_chems %>%
    left_join(pubchem, by = "Cas-nr") %>%
    rename(`PubChem CID` = `PubChem CID`)  # Optional, just keeps column label
  
  
  result_url <- eventReactive(input$submit, {
    req(input$smiles)
    runjs("$('#main_tabs li a[data-value=\"search\"]').tab('show');")
    smiles_to_pubchem_url(input$smiles, cid_list, substructure = input$search, department=input$department)
  })
  
  # When the tab becomes active and table is ready
  observeEvent({
    input$main_tabs
    input$dep_table_ready
  }, {
    if (input$main_tabs == "browse" && isTRUE(input$dep_table_ready) && !is.null(search_trigger())) {
      term <- search_trigger()
      runjs(sprintf("$('#dep_chemicals input[type=search]').val('%s').trigger('input');", term))
      search_trigger(NULL)
    }
  })
  
  # CC(=O)C is acetone
  observeEvent(result_url(), {
    req(result_url())
    updateTabsetPanel(session, "main_tabs", selected = 'search')
    
    # Extract CID and CAS
    url <- result_url()
    cid <- str_extract(url, '[0-9]+$')
    if (!is.na(cid)) {
      cas <- dep_chems[dep_chems$`PubChem CID` == cid, "Cas-nr"] %>% unique()
      cas <- cas[!is.na(cas)][1]
      cas_hits(cas)
      search_trigger(cas)
    } else {
      cid_pattern <- str_extract(url, '(?:[0-9]+[,&])+') %>% str_replace_all(., ',', '|') %>% str_remove(., '&')
      cas <- dep_chems[which(str_detect(as.character(dep_chems$`PubChem CID`), cid_pattern)), "Cas-nr"] %>% unlist() %>% as.vector()
      cas_hits(cas)
    }
  })
  
  output$result_url <- renderUI({
    url <- tryCatch(result_url(), error = function(e) NULL)
    if (is.null(url)) {
      print('url is null')
      tags$div(tags$br(),
               tags$p(style = 'color:grey;', tags$i("Perform a search using the menu on the left.")))
    } else {
      if (url == 'No results') {
        tags$p(style = 'color:grey;', tags$i("No matching compounds found!"))
      } else{
        tags$div(
          tags$br(),
          tags$a(href = url, target='_blank', url),
          tags$iframe(src = url, width = "100%", height = "800px", style = "border-width: 1px; border-style: solid;")
        )
      }
    }
  })
  
  ordered_chems <- eventReactive(input$ordered, {
    #include google sheets with ordered chemicals
    order_gsheet <- 'https://docs.google.com/spreadsheets/d/1bhcbzRX7Y_mJAx11_-TvE7ksrw5xSUORBgOQsUHOXAA/edit?gid=709388515#gid=709388515'
    gs4_deauth()  # skip auth for public sheets
    tabs <-  sheet_names(order_gsheet)
    tabs <- tabs[!str_detect(tabs, '(Status)|(Currency)|(penalty)|(How to)')]
    sheets <- map(tabs, safe_read)
    names(sheets) <- tabs
    
    # Drop failed (NULL) reads
    sheets <- compact(sheets)
    sheets <- sheets[map_lgl(sheets, has_mostly_cas_column)]
    sheets
  })
  
  
  output$dep_chemicals <- renderDT({
    if (!is.null(cas_hits())) {
      show('subset')
    } else {
      hide('subset')
    }
    df <- dep_chems
    
    if (input$subset && !is.null(cas_hits())) {
      # browser()
      df <- df[df$`Cas-nr` %in% cas_hits(),]
    }
    if (input$ordered) {
      df$Status <- 'in stock at last inventory'
      dfo <- ordered_chems()
      dfo <- bind_rows(dfo, .id = "sheet_name")
      dfo$Status <- paste0(dfo$Status, ' (', str_extract(dfo$sheet_name, '[2][0-9]{3}'), ')')
      dfo <- dfo[,c("Product description", "CAS number (required for chemicals)", "Amount", "Location", "Status")]
      colnames(dfo) <- c('Name', 'Cas-nr', 'Quantity', 'Location', 'Status')
      dfo$`Cas-nr` <- vapply(dfo$`Cas-nr`, function(x) {
        if (length(x) == 0 || is.null(x)) NA_character_ else as.character(x)
      }, character(1))
      
      dfo <- dfo %>%
        left_join(pubchem, by = "Cas-nr") %>%
        rename(`PubChem CID` = `PubChem CID`)  # Optional, just keeps column label
      df <- bind_rows(df, dfo)
    }
    df$`PubChem CID` <- ifelse(
      is.na(df$`PubChem CID`),
      NA,
      paste0(
        '<a href="https://pubchem.ncbi.nlm.nih.gov/compound/',
        df$`PubChem CID`,
        '" target="_blank">',
        df$`PubChem CID`,
        '</a>'
      )
    )
    if (input$ordered) {
      df <- df[,c("Name", "Cas-nr", "Location", "Quantity", "PubChem CID", "Status")]
    } else {
      df <- df[,c("Name", "Cas-nr", "Location", "Quantity", "PubChem CID")]
    }
    df$`Cas-nr` <- ifelse(
      is.na(df$`Cas-nr`),
      NA,
      paste0(
        '<a href="https://pubchem.ncbi.nlm.nih.gov/#query=',
        df$`Cas-nr`,
        '" target="_blank">',
        df$`Cas-nr`,
        '</a>'
      )
    )
    datatable(df, escape = FALSE, rownames = FALSE, 
              options = list(pageLength = 20,
                             initComplete = JS("function(settings, json) { Shiny.setInputValue('dep_table_ready', true, {priority: 'event'}); }")))
  })
  
  observeEvent(input$from_kekule, {
    updateTextInput(session, "smiles", value = input$from_kekule)
  })
  
  # observeEvent(input$main_tabs, {
  #   updateTabsetPanel(session, "main_tabs", selected = input$main_tabs)
  # })
  
}

shinyApp(ui, server)

