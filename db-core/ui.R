library(shiny)
library(bslib)
library(htmltools)
library(ggplot2)
library(dplyr)
library(DT)
library(shinyjqui)
library(shinyWidgets)
library(networkD3)
options(shiny.maxRequestSize=1024*1024^2)

# readxl::read_excel('scdb_core.xlsx',sheet = 'Sheet1') %>% 
#     tidyr::separate_rows(Tissue,sep = ';') %>% pull(Tissue) |> unique() |> length()
#Total tissue number
coretable <- readxl::read_excel('scdb_core.xlsx',sheet = 'Sheet3')

ui <- 
    navbarPage(
    includeCSS("www/style.css"), 
    title = div(img(src = "logo-removebg-preview.png",width = '80px'),class = 'logo'),
    bg = "#F9F7F3",
    id = "inTabset",
    theme = bslib::bs_theme(),
    ###Home-----
    tabPanel(title = "Home",
              icon = icon('home',lib="glyphicon"),
              div(class = 'home',
              fluidPage(
                  card(
                      card_header("Introduction",class = "Introduction", style = "font-size: 24px;"),
                      status = "primary",
                      width = 12,
                      height = NULL,
                      card_body(
                          uiOutput("HOME_output_text")
                      )

                  ),
                  layout_columns(
                      fill = FALSE,
                      value_box(
                          "Total cell",
                          '>4 Milions',
                          showcase = bsicons::bs_icon("Egg fried", size = NULL),
                          style = 'background-color: #5092D0!important;'
                      ),
                      value_box(
                          "Total tissue",
                          '28',
                          showcase = bsicons::bs_icon("person-hearts", size = NULL),
                          style = 'background-color: #459943!important;'
                      ),
                      value_box(
                          "Total disease",
                          '48',
                          showcase = bsicons::bs_icon("universal-access", size = NULL),
                          style = 'background-color: #ea9c9d!important;'
                      )
                    
                  ),
                  tags$br(),
                  card(
                      card_image(file = "www/DRCTDB workflow_A4.jpg" ),
                      style = 'background:#FFFFFF'
                  ),
                  style = "font-size:150%;width:80%;"
              )
              )
        ),
     ###Search-----
    tabPanel(title = "Search",
               icon = icon('search',lib = 'glyphicon'),
               fluidPage(
                   div(
                   card(
                       card_header("Dataset overview"),
                       div(DT::dataTableOutput('coretable',width = "100%"))
                   ),
                   card(
                       actionButton("go_to_panel", "Expore dataset....",class = 'jump')
                       #textOutput('test')
                   )

                    ),style = "font-size:120%;width:80%;")
               ),


    ###Results-----
    tabPanel(title = "Results",
             value =  "Results",
              icon = icon('chart-line',lib="font-awesome"),
              fluidPage(
                  layout_column_wrap(
                      width = 1/2,
                      height = NULL,
                      navset_card_tab(
                          title = textOutput('study_name'),
                          id = "umap", height = "800px",
                          full_screen = TRUE,
                          #static umap
                          # nav_panel("raw_uamp",
                          #          div(imageOutput("show_umap1"),style = "margin-left: auto; margin-right: auto;background:#FFFFFF;")
                          # ),
                          #dynamic umap
                          nav_panel("UMAP",
                                    div(imageOutput("show_umap"),style = "margin-left: auto; margin-right: auto;")
                          ),
                          nav_panel("DRCT",
                                    div(DT::dataTableOutput('show_DRCT_table',width = "100%"))
                          ),
                          nav_panel("DAR",
                                    div(DT::dataTableOutput('show_DER_table',width = "100%"))
                          )
                      ),
                      navset_card_tab(
                          title = div(class = 'select disease',
                                        uiOutput('dropdown')
                                      ),
                          id = "tabset1", height = "800px",
                          full_screen = TRUE,
                          nav_panel("CCC",
                                        switchInput(
                                            inputId = "switchccc",
                                            label = '-To-',
                                            value = TRUE,
                                            onLabel = "Figure",
                                            offLabel = "Table",
                                            size = 'normal'
                                        ),
                                    div(uiOutput('show_ccc'),style = "margin-left: auto; margin-right: auto;")
                                   
                          ),
                          nav_panel("ATAC",
                                    div(uiOutput('show_atac_enrich'),style = "margin-left: auto; margin-right: auto;")
                          ),
                          nav_panel('RNA',
                                    div(uiOutput('show_rna_enrich'),style = "margin-left: auto; margin-right: auto;")
                          ),
                          nav_panel("GRN",
                                   div(forceNetworkOutput("grn_output"))
                          )

                      )
                  ),
                  style = "font-size:150%;width:100%;"
              )

              ),

    ###Online tools-----

    tabPanel(title = "Tutorials", 
              icon = icon('bookmark',lib = 'glyphicon'),
              p("Second tab content.")),
    tabPanel(title = "Download", 
              icon = icon('download',lib = 'glyphicon'),
             fluidPage(
                 div(
                     card(
                         card_header("Dataset overview"),
                         div(DT::dataTableOutput('download_table',width = "100%"))
                     )
                     
                 ),style = "font-size:120%;width:80%;")
             )
             ,
    tabPanel(title = "Contact",
              icon =  icon('envelope',lib = 'glyphicon'),
             fluidPage(
                 card(
                     card_header("Contact us",class = "Introduction", style = "font-size:200%;line-height: 150%; padding: 20px;"),
                     status = "primary",
                     width = 10,
                     height = NULL,
                     card_body(
                         uiOutput("Contact_text")
                     )
                     
                 ),style = "font-size:120%;width:80%;")
             ),
    nav_spacer(),
    nav_item(tags$a(shiny::icon("github"), "DRCTDB", href = "https://github.com/jiang-junyao/DRCTdb", target = "_blank"))
)
