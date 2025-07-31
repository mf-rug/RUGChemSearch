library(tidyverse)
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)
library(DT)
library(googlesheets4)
library(gargle)
library(rcdk)
library(httr)
library(jsonlite)


# source("secret_retrieval_function.R")

if (Sys.getenv("SHINY_PORT") != "") {
  # We're on shinyapps.io
  env_name <- "venv_secret_retrieval_shiny"
  envs <- reticulate::virtualenv_list()
  if (!(env_name %in% envs)) {
    reticulate::virtualenv_create(env_name)
    # reticulate::virtualenv_install(env_name, packages = c("cryptography"))
  }
  reticulate::use_virtualenv(env_name, required = TRUE)
  
} else {
  # We're running locally
    reticulate::use_python("/opt/homebrew/bin/python3", required = TRUE)
}


ui <- fluidPage(
  useShinyjs(),
  titlePanel("Chemicals Search"),
  tags$head(
    tags$style("
  .hover-container {
    position: relative;
    display: inline-block;
    margin-top:-8px;
    margin-bottom:-8px;
  }
  .hover-container:hover .thumb {
    outline: 1px solid cyan;
    filter: brightness(1.2);
    transition: all 0.2s ease;
  }
  .hover-container .hover-preview {
    display: none;
    background-color: white !important;
    position: absolute;
    top: 50%;
    transform: translateY(-50%);
    right: 110%;
    width: 300px;
    z-index: 10;
    filter: brightness(1.06);
    border: 1px solid #ccc;
  }
  .hover-container:hover .hover-preview {
    display: block;
  }
    /* Default button styling */
    .btn-default {
        --bs-btn-bg: #cceef8 !important;
        padding-right: 8px;
        padding-left: 8px;
    }
    /* Default button styling */
    .btn-run {
        background: #200640;
        color: aliceblue;
        font-size: 13pt;
        box-shadow: #9700ff 0px 0px 6px 1px inset !important;
        padding-right: 20px;
        padding-left: 20px;
    }
               "),
    tags$script(src = "https://cdn.jsdelivr.net/npm/kekule/dist/kekule.min.js"),
    tags$link(rel = "stylesheet", href = "https://cdn.jsdelivr.net/npm/kekule/dist/themes/default/kekule.css")
  ),
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

      tabsetPanel(
        id = "main_tabs",
        tabPanel('browse department chemicals',
                 br(),
                 fluidRow(column(
                   6,
                   materialSwitch(
                     'ordered',
                     'include chemicals from the Chemical Order Form (googlesheet)',
                     status = 'primary',
                     width = '100%',
                     right = TRUE
                   )
                 ),
                 column(
                   6,
                   materialSwitch(
                     'regex',
                     'enable regular expressions search',
                     status = 'primary',
                     width = '100%',
                     right = TRUE
                   )
                 )),
                 # fluidRow(column(12, uiOutput('update_with_secret'))), br(),
                 fluidRow(column(12,
                   hidden(
                     materialSwitch(
                       'subset',
                       'Only show the last search hits',
                       status = 'primary',
                       value = TRUE,
                       right = TRUE
                     )
                 ))),
                 withSpinner(DTOutput("dep_chemicals")),
                 br(), br(),br()
        ),
        
        tabPanel('search by structure',
                 sidebarLayout(
                   sidebarPanel(width = 4,
                                tags$p(style = "font-size:13pt;", tags$strong("Draw a molecule")),
                                tags$div(id = "editor", style = "width: 100%; height: 50vh; border: 1px solid #ccc;"),
                                verbatimTextOutput("smiles"),
                                
                                
                                br(), 
                                textInput("smiles", "Draw above or directly enter SMILES", ""),
                                radioGroupButtons('search', NULL, c('identity', 'similarity', 'substructure', 'superstructure', '3d similarity')),
                                materialSwitch(
                                  'department',
                                  HTML('Restrict to <a href="#" onclick="Shiny.setInputValue(\'main_tabs\', \'browse department chemicals\')">department chemicals</a>'),
                                  value = TRUE,
                                  status = 'primary',
                                  right = TRUE
                                ),
                                div(style = "text-align: center;",
                                  actionButton(inputId = "submit",label =  "Search", class = 'btn-run')
                                ),
                   ),
                   mainPanel(withSpinner(uiOutput("result_url")))
      )
    )
  )
)

