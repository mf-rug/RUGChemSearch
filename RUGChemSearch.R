library(httr)
library(tidyverse)
library(shiny)
library(shinyWidgets)
library(shinycssloaders)

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
  if (length(match_cids) == 0) return(NULL)
  
  # Step 4: Create PubChem web URL (limit output)
  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/#query=",
    paste0(head(match_cids, max_display), collapse = ","),
    "&collection=compound"
  )
  
  return(url)
}


ui <- fluidPage(
  titlePanel("SMILES to PubChem Filtered Search"),
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
    editor = new Kekule.Editor.Composer(document.getElementById('editor'));
  }

  // Ensure Shiny is ready
  if (window.Shiny) {
    initKekule();
  } else {
    document.addEventListener('shiny:connected', function() {
      initKekule();
    });
  }

  document.addEventListener('click', function(e) {
    if (e.target && e.target.id === 'get_smiles') {
      if (editor) {
        var mol = editor.getChemObj();
        var smiles = Kekule.IO.saveFormatData(mol, 'smi');
        Shiny.setInputValue('from_kekule', smiles, {priority: 'event'});
      }
    }
  });
")),
      
      br(),
      actionButton("get_smiles", "Get SMILES"),
      br(), br(),
      textInput("smiles", "Draw above or directly enter SMILES", ""),
      radioGroupButtons('search', NULL, c('identity', 'similarity', 'substructure', 'superstructure', '3d similarity')),
      materialSwitch('department', 'Restrict to department chemicals',value = TRUE),
      actionButton("submit", "Search")
    ),
    mainPanel(
      withSpinner(uiOutput("result_url"))
    )
  )
)

server <- function(input, output, session) {
  cid_list <- readLines("/Users/max/Documents/Code/bc_pubchem.txt")
  
  result_url <- eventReactive(input$submit, {
    req(input$smiles)
    smiles_to_pubchem_url(input$smiles, cid_list, substructure = input$search, department=input$department)
  })
  
  output$result_url <- renderUI({
    url <- result_url()
    if (is.null(url)) {
      tags$p("No matching compounds found.")
    } else {
      tags$div(
        style = "word-wrap: break-word; white-space: normal; width: 100%;",
        tags$a(href = url, url, target = "_blank")
      )
    }
  })
  observeEvent(input$from_kekule, {
    print('sdfa')
    updateTextInput(session, "smiles", value = input$from_kekule)
  })
  
}

shinyApp(ui, server)

