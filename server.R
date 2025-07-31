max_result_display <- 10000
max_wait_rounds <- 3

time_since_last_call <- local({
  last_time <- Sys.time()
  function() {
    now <- Sys.time()
    elapsed <- as.numeric(difftime(now, last_time, units = "secs"))
    last_time <<- now
    cat(sprintf("Elapsed: %.3f seconds\n", elapsed))
  }
})

extract_ghs_codes <- function(urls) {
  if (is.null(urls) || all(is.na(urls))) return(character(0))
  sub(".*(GHS[0-9]{2})\\.svg.*", "\\1", urls)
}

get_ghs_pictograms <- function(cid, cache_file = "ghs_cache.json") {
  # Load existing cache (or initialize empty list)
  if (file.exists(cache_file)) {
    cache <- fromJSON(cache_file, simplifyVector = FALSE)
  } else {
    cache <- list()
  }
  
  # Return cached result if available
  if (!is.null(cache[[as.character(cid)]])) {
    message(sprintf("Getting GHS cache result for CID %s", cid))
    return(cache[[as.character(cid)]])
  }
  
  # Not cached: fetch from API
  message(sprintf("Fetching CID %s", cid))
  url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON/?heading=GHS+Classification", cid)
  resp <- GET(url)
  
  if (status_code(resp) != 200) {
    cache[[as.character(cid)]] <- NA
    write_json(cache, cache_file, auto_unbox = TRUE, pretty = TRUE)
    return(NA)
  }
  
  json <- content(resp, as = "parsed", type = "application/json")
  
  # Extract GHS pictogram URLs
  pictogram_urls <- str_extract_all(
    paste(unlist(json$Record$Section), collapse = '\n'),
    'https://pubchem.ncbi.nlm.nih.gov/images/ghs/GHS[0-9][0-9].svg'
  )[[1]] %>% unique() %>% sort()
  
  # pictogram_urls <- tryCatch({
  #   ghs_sections <- json$Record$Section %>%
  #     keep(~ .x$TOCHeading == "Safety and Hazards") %>%
  #     map("Section", .default = NULL) %>% compact() %>% purrr::flatten() %>%
  #     keep(~ .x$TOCHeading == "Hazards Identification") %>%
  #     map("Section", .default = NULL) %>% compact() %>% purrr::flatten() %>%
  #     keep(~ .x$TOCHeading == "GHS Classification")
  #   
  #   ghs_sections %>%
  #     map("Information", .default = NULL) %>% compact() %>% purrr::flatten() %>%
  #     keep(~ .x$Name == "Pictogram(s)") %>%
  #     map("Value") %>%
  #     map("StringWithMarkup", .default = NULL) %>% compact() %>% purrr::flatten() %>%
  #     map("Markup", .default = NULL) %>% compact() %>% purrr::flatten() %>%
  #     keep(~ .x$Type == "Icon") %>%
  #     map_chr("URL") %>%
  #     unique()
  # }, error = function(e) NA)
  
  # Update cache and save
  cache[[as.character(cid)]] <- pictogram_urls
  write_json(cache, cache_file, auto_unbox = TRUE, pretty = TRUE)
  
  pictogram_urls
}

pubchem_cid_to_url <- function(cid, max_display=max_result_display) {
  # Create PubChem web URL (limit output)
  if (length(cid) > 1) {
    paste0(
      "https://pubchem.ncbi.nlm.nih.gov/#query=",
      paste0(head(cid, max_display), collapse = ","),
      "&collection=compound"
    )   
  } else {
    paste0(
      "https://pubchem.ncbi.nlm.nih.gov/compound/",
      cid
    )
  }
}

smiles_to_pubchem_url <- function(smiles, cid_list, max_display = max_result_display, wait_seconds = 2, search_type='identity', department=TRUE) {
  # Step 1: Start the search job
  url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/', search_type, '/smiles/JSON')
  setProgress(value = 0.2, detail = paste('initiating a ', search_type, 'search via', url))
  
  init_res <- POST(
    url,
    body = list(smiles = smiles),
    encode = "form"
  )
  
  if (!init_res$status_code %in% c(202,200)) {
    stop("Search initiation failed.")
  }
  
  listkey <- content(init_res, as = "parsed")$Waiting$ListKey
  setProgress(value = 0.3, detail = paste('got listkey', listkey))
  if (!is.null(listkey)) {
    cat('listkey is -->', listkey, '<--\n')
    
    # Step 2: Poll for completion
    Sys.sleep(wait_seconds)  # wait before first try
    i <- 1
    repeat {
      poll_res <- GET(paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/",
        listkey,
        "/cids/JSON"
      ))    
      if (poll_res$status_code == 200) {
        cat('Finished', '\n')
        break
      } else if (poll_res$status_code == 202 && i <= max_wait_rounds) {
        setProgress(value = 0.3 + (i / 100), detail = paste('waiting for result of listkey', listkey))
        i <- i + 1
        cat('Running', '\n')
        Sys.sleep(wait_seconds)  # still processing
      } else if (i > max_wait_rounds) {
        message(paste("PubChem search timed out:", poll_res))
        return('PubChem search timed out')
      } else {
        message(paste("Unexpected error while polling for search result:", poll_res))
        return('PubChem search fail')
      }
    }
    cids_found <- content(poll_res, as = "parsed")$IdentifierList$CID %>% unlist()
  } else {
    message('no listkey, result directly returned\n')
    result <- content(init_res, as = "parsed")
    cids_found <- result$PC_Compounds[[1]]$id$id$cid
  }
  if (is.null(cids_found)) {
    message('No results', '\n')
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
  if (length(match_cids) == 0) {
    print('0 matches, returning unfiltered result url')
    return(list('No matches found', pubchem_cid_to_url(cids_found)))
  }
  return(list('Matches found', pubchem_cid_to_url(match_cids)))
}

get_smiles <- function(cid, cache_file = "all_smiles.json") {
  if (is.na(cid) || cid == 'NA') {
    return(NA)
  }
  # Load existing cache (or initialize empty list)
  if (file.exists(cache_file)) {
    cache <- fromJSON(cache_file, simplifyVector = FALSE)
  } else {
    cache <- list()
  }
  # Return cached result if available
  if (!is.null(cache[[as.character(cid)]]) && is.character(cache[[as.character(cid)]])) {
    message(sprintf("Getting smiles cache result for CID %s", cid))
    # if (!is.character(cache[[as.character(cid)]])) {
    #   return(NA)
    # } else {
      return(cache[[as.character(cid)]])
    # }
  }
  
  # Not cached: fetch from API
  message(sprintf("Fetching CID %s", cid))
  url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/smiles/JSON/", cid)
  resp <- GET(url)
  
  if (status_code(resp) != 200) {
    cache[[as.character(cid)]] <- NA
    write_json(cache, cache_file, auto_unbox = TRUE, pretty = TRUE)
    return(NA)
  }
  
  json <- content(resp, as = "parsed", type = "application/json")
  
  # Extract GHS pictogram URLs
  smiles <- tryCatch({
    json$PropertyTable$Properties[[1]]$SMILES
  }, error = function(e) NA)
  
  # Update cache and save
  cache[[as.character(cid)]] <- smiles
  write_json(cache, cache_file, auto_unbox = TRUE, pretty = TRUE)
  
  smiles
}

safe_read <- possibly(function(order_gsheet, name) {
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



# server 
server <- function(input, output, session) {
  
  search_trigger <- reactiveVal(NULL)
  cas_hits <- reactiveVal(NULL)
  
  cid_list <- readLines("bc_pubchem.txt")
  
  # Load data
  # print('ho')
  # client_secret_json = "client_secret_746595120067-ci7drg9hmtegjgp7hucd3kceup39p2mo.apps.googleusercontent.com.json"
  # print('hu')
  # # gs4_auth_configure(path = client_secret_json)  # from GCP
  # # gs4_auth(cache = FALSE, email = NULL)
  # gs4_auth(path = client_secret_json)
  # print('ha')
  # # online csv is https://docs.google.com/spreadsheets/d/12Ow84zk6B3lhcui2NXxSd6k7BttRUrjqqZgL3vm7iPE/edit?usp=sharing
  # dep_chems <- read_sheet('https://docs.google.com/spreadsheets/d/12Ow84zk6B3lhcui2NXxSd6k7BttRUrjqqZgL3vm7iPE/edit?usp=sharing')
  # 
  dep_chems <- read_csv('biotechnology_chemicals_fixed.csv', show_col_types = FALSE, col_names = TRUE)
  pubchem <- read_csv('biochemistry_all_cas_and_pubchem_ids.csv', show_col_types = FALSE, col_names = TRUE)
  dep_chems <- dep_chems %>%
    left_join(pubchem, by = "Cas-nr") 
  
  # Load SMILES map
  smiles_map <- jsonlite::fromJSON("all_smiles.json", simplifyVector = FALSE)
  
  # Add SMILES from JSON map (safe lookup)
  dep_chems$Smiles <- sapply(
    as.character(dep_chems$`PubChem CID`),
    function(cid) {
      val <- smiles_map[[cid]]
      if (!is.null(val) && length(val) == 1) val else NA_character_
    }
  )
  
  # Fill missing SMILES via fallback function
  missing <- is.na(dep_chems$Smiles) & !is.na(dep_chems$`PubChem CID`)
  if (any(missing)) {
    dep_chems$Smiles[missing] <- sapply(
      as.numeric(dep_chems$`PubChem CID`[missing]),
      get_smiles
    )
  }
  
  # deprecated by pubchem and doesnt seem correct data
  # dep_chems$canSmiles <- map_chr(dep_chems$`PubChem CID`, get_can_smiles)
  
  # dep_chems$canSmiles <- map_chr(dep_chems$Smiles, function(x) {
  #   if (!is.character(x) || is.na(x)) { return(NA) }
  #   mol <- parse.smiles(x)[[1]]
  #   get.smiles(mol, smiles.flavors(c("Canonical")))
  # })
  
  # Load can_smiles map
  can_smiles_map <- jsonlite::fromJSON("all_can_smiles.json", simplifyVector = FALSE)
  
  # Add can_smiles from JSON map (safe lookup)
  dep_chems$canSmiles <- sapply(
    as.character(dep_chems$`PubChem CID`),
    function(cid) {
      val <- can_smiles_map[[cid]]
      if (!is.null(val) && length(val) == 1) val else NA_character_
    }
  )
  
  # Fill missing can_smiles via fallback function
  missing <- (is.na(dep_chems$canSmiles) | dep_chems$canSmiles == '') & !is.na(dep_chems$`PubChem CID`)

  dep_chems$canSmiles[missing] <- map_chr(dep_chems$Smiles[missing], function(x) {
    # if (!is.character(x) || is.na(x)) { return(NA) }
    mol <- parse.smiles(x)[[1]]
    get.smiles(mol, smiles.flavors(c("Canonical")))
  })
  
  result_url <- eventReactive(input$submit, {
    req(input$smiles)
    runjs("$('#main_tabs li a[data-value=\"search\"]').tab('show');")
    if (input$search == 'identity' && input$department) {
      'smiles'
    } else {
      withProgress(smiles_to_pubchem_url(input$smiles, cid_list, search_type = input$search, department=input$department),
                   value = 0.1, message = 'Running search structure in PubChem')
    }
  })
  
  # When the tab becomes active and table is ready
  observeEvent({
    input$main_tabs
    input$dep_table_ready
  }, {
    if (input$main_tabs == "browse department chemicals" && isTRUE(input$dep_table_ready) && !is.null(search_trigger())) {
      term <- search_trigger()
      runjs(sprintf("$('#dep_chemicals input[type=search]').val('%s').trigger('input');", term))
      search_trigger(NULL)
    }
  })
  
  # CC(=O)C is acetone
  observeEvent(result_url(), {
    req(result_url())
    if (length(result_url()) > 1) {
      outcome <- result_url()[[1]]
      
      if (outcome == 'Matches found') {
        url <- result_url()[[2]]
        updateTabsetPanel(session, "main_tabs", selected = 'search')
        
        # Extract CID and CAS
        cid <- str_extract(url, '[0-9]+$')
        # browser()
        if (!is.na(cid)) {
          cas <- dep_chems[dep_chems$`PubChem CID` == cid, "Cas-nr"] %>% unique()
          cas <- cas[!is.na(cas)][1]
          cas_hits(cas)
          search_trigger(cas)
        } else {
          # cid_pattern <- str_extract(url, '(?:[0-9]+[,&])+') %>% str_replace_all(., ',', '|') %>% str_remove(., '&')
          # cas <- dep_chems[which(str_detect(as.character(dep_chems$`PubChem CID`), cid_pattern)), "Cas-nr"] %>% unlist() %>% as.vector()
          
          cids <- str_extract(url, '(?:[0-9]+[,&])+') %>% str_remove(., '&') %>% str_split(.,',')
          cids <- cids[[1]]
          cas <- dep_chems[which(as.character(dep_chems$`PubChem CID`) %in% cids), "Cas-nr"] %>% unlist() %>% as.vector()
          cas_hits(cas)
        }
      } else {
        cas_hits(NULL)
      }
    } else {
      if (input$search == 'identity' && input$department) {
        mol <- parse.smiles(input$smiles)[[1]]
        tryCatch({
          canSmiles <- get.smiles(mol, smiles.flavors(c("Canonical")))
          hits.df <- dep_chems[dep_chems$canSmiles == canSmiles & !is.na(dep_chems$canSmiles),]
          cas_hits(hits.df$`Cas-nr`)
        }, error = function(e) {
          cas_hits(NULL)
          NA
        })
      } else {
        cas_hits(NULL)
      }
    }
  })
  
  output$update_with_secret <- renderUI({
    last_update <- file.mtime('biotechnology_chemicals_fixed.csv') 
    if (difftime(Sys.time(), last_update, units = 'days') >1) {
      actionButton('update_with_secret_btn', 
                   paste0(
                     'Initiate update of chemicals list (last update: ',
                     last_update %>% format(., format = '%d.%m.%Y %H:%M'),
                     ')'
                   ))
    }
  })
  
  observeEvent(input$update_with_secret_btn,{
    showModal(modalDialog(
      title = "This action requires extra privileges",
      passwordInput("user_pw", "Enter password:"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("submit_text", "Submit")
      )
    ))
  })
  
  
  observeEvent(input$submit_text,{
    removeModal()
    ug_pw <- retrieve_secret_via_pw(input$user_pw, 'secret.json')
    print('ug_pw')
  })
  
  
  output$result_url <- renderUI({
    req(result_url())
    if (length(result_url()) > 1) {
      url_result <- tryCatch(result_url(), error = function(e) NULL)
      outcome <- url_result[[1]]
      url <- url_result[[2]]
    } else {
      cid <- as.numeric(unique(dep_chems[dep_chems$`Cas-nr` %in% cas_hits(),"PubChem CID"])[1])
      if (!is.na(cid) && is.numeric(cid)) {
        outcome <- 'Matches found'
        url <- pubchem_cid_to_url(cid)
      } else {
        if (str_detect(result_url(), "^PubChem search")) {
          outcome <- result_url()
        } else if (length(result_url()) == 1 && result_url() == 'smiles') {
          outcome <- 'No matches found'
          url <- paste0('https://pubchem.ncbi.nlm.nih.gov/#query=', input$smiles)
        } else {
          outcome <- 'No matches found'
        }
      }
    }
    if (is.null(url)) {
      tags$div(tags$br(),
               tags$p(style = 'color:grey;', tags$i("Perform a search using the menu on the left.")))
    } else {
      # browser()
      if (length(result_url()) > 1) {
        result_url <- result_url()[[1]]
      } else {
        result_url <- result_url()
      }
      if (outcome == 'Matches found' || !input$department) {
        tags$div(
          tags$br(),
          tags$div(
            style = "word-wrap: break-word; white-space: normal; width: 100%;",
            tags$a(href = url, ifelse(str_count(url) > 200, paste0(str_extract(url, '^.{200}'), '...'), url), target = "_blank")
          ),
          tags$iframe(src = url, width = "100%", height = "800px", style = "border-width: 1px; border-style: solid;")
        )
      } else if (str_detect(result_url, "^PubChem search")) {
        tags$div(
          tags$br(),
          tags$div(
            style = 'color:grey;',
            HTML(paste0(
              "<i>", result_url, "<br>",
              "A common reason is that you the search is not specific enough and would return an enormous number of molecules, e.g. if you do a substructure search with a very common functional group.",
              '</i>'
            ))
          )
        )
      } else {
        tags$div(
          tags$br(),
          tags$div(
            style = 'color:grey;',
            HTML(paste0(
              "<i>No matching compounds found in ",
              if (isolate(input$department)) {
                "the department's chemical inventory!<br/>Uncheck 'Restrict to department chemicals' in the left menu to see PubChem result."
              } else {
                "all of PubChem."
              },
              '</i>'
            ))
          )
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
    sheets <- map(tabs, function(x) safe_read(order_gsheet, x))
    names(sheets) <- tabs
    
    # Drop failed (NULL) reads
    sheets <- compact(sheets)
    sheets <- sheets[map_lgl(sheets, has_mostly_cas_column)]
    sheets
  })
  
  
  observe({
    if (input$department) {
      if (length(cas_hits()) > 0) {
        shinyjs::show('subset')
        updateMaterialSwitch(session, 'subset', TRUE)
      } else {
        shinyjs::hide('subset')
        updateMaterialSwitch(session, 'subset', FALSE)
      }
    } else {
      updateMaterialSwitch(session, 'subset', FALSE)
      shinyjs::hide('subset')
    }
  })
  
  output$dep_chemicals <- renderDT({
    
    df <- dep_chems
    
    if (input$ordered) {
      df$Status <- 'in stock at last inventory'
      dfo <- ordered_chems()
      dfo <- bind_rows(dfo, .id = "sheet_name")
      if (is.data.frame(dfo) && nrow(dfo) > 0) {
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
    }
    
    if (input$subset && !is.null(cas_hits())) {
      df <- df[df$`Cas-nr` %in% cas_hits(),]
    }
    
    # Load GHS map
    if (file.exists("ghs_cache.json")) {
      ghs_map <- jsonlite::fromJSON("ghs_cache.json", simplifyVector = FALSE)
    } else {
      ghs_map <- NULL
    }
    # Add GHS from JSON map (safe lookup)
    hazard_img_urls <- lapply(
      as.character(df$`PubChem CID`),
      function(cid) {
        val <- ghs_map[[cid]]
        if (!is.null(val)) unlist(val) else NA_character_
      }
    )
    
    # Fill missing GHS links via fallback function
    missing <- is.na(hazard_img_urls) & !is.na(df$`PubChem CID`)
    
    if (any(missing)) {
      hazard_img_urls[missing] <- lapply(
        as.numeric(df$`PubChem CID`[missing]),
        get_ghs_pictograms
      )
    }
    
    # data for title of pictograms
    ghs_data <- jsonlite::fromJSON("ghs_data.json")
    ghs_labels <- setNames(ghs_data$name, ghs_data$code)
    hazard_codes <- lapply(hazard_img_urls, extract_ghs_codes)
    hazard_labels <- lapply(hazard_codes, function(codes) {
      ghs_labels[codes]
    })
    
    img_html <- mapply(function(urls, labels) {
      if (length(urls) == 0 || all(is.na(urls))) return("")
      paste0(
        mapply(function(url, label) {
          sprintf('<img src="%s" width="30" height="30" style="margin-right:4px; margin-top:-4px; margin-bottom:-4px;" title="%s">', url, label)
        }, urls, labels, USE.NAMES = FALSE),
        collapse = ""
      )
    }, hazard_img_urls, hazard_labels, USE.NAMES = FALSE)
    
    df$Hazards <- img_html
    
    strc_img_urls <- paste0('https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=', df$`PubChem CID`, '&t=l')
    df$Structure <- sprintf(
      '<div class="hover-container" style="position: relative; display: inline-block;">
         <div class="thumb-wrapper" style="display: inline-block;">
           <img src="%s" class="thumb" style="width:50px;filter: brightness(1.06);">
         </div>
         <img class="hover-preview" src="%s">
       </div>',
      strc_img_urls, strc_img_urls
    )
    
    df$Structure[is.na(df$`PubChem CID`)] <- NA
    
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
    if (input$ordered && is.data.frame(dfo) && nrow(dfo) > 0) {
      df <- df[,c("Name", "Cas-nr", "Location", "Hazards", "Quantity", "PubChem CID", 'Structure', "Status")]
    } else {
      df <- df[,c("Name", "Cas-nr", "Location", "Hazards", "Quantity", 'Structure', "PubChem CID")]
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
                             search = list(regex = input$regex),
                             initComplete = JS("function(settings, json) { Shiny.setInputValue('dep_table_ready', true, {priority: 'event'}); }"))
    )  %>%
      formatStyle(0, target = 'row', lineHeight='75%') 
  })
  
  observeEvent(input$from_kekule, {
    updateTextInput(session, "smiles", value = input$from_kekule)
  })
}
