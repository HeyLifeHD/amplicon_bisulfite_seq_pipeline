# library(magrittr)
# library(shinythemes)
# library(shinycssloaders)
# library(shinydashboard)
# library(shinybusy)
# library(shiny)
# library(shinyBS)
# library(shinyjs)
# Libraries are now specified in the .Rprofile of the shiny user for a faster start up

ui <- navbarPage(strong("AmpliconDesign"),id="page",
                 windowTitle = "AmpliconDesign | Primer design for DNA methylation assays",


                 tabPanel("Home",
                          mainPanel(
                            h1(strong("AmpliconDesign"),align="center"),
                            column(8,includeMarkdown("home.Rmd"),offset=2,align="center")
                          ,width = 11)),

                 tabPanel("MassArray",
                          add_busy_spinner(spin = "fading-circle",position="full-page",timeout = 1000),
                          fluidPage(theme=shinytheme("cosmo"),
                                    sidebarLayout(
                                      sidebarPanel(
                                        ##general setting for the whole app###
                                        shiny::tags$head(
                                          shiny::tags$style(
                                            HTML(".shiny-notification {
                                                 position:fixed;
                                                 top: calc(50%);;
                                                 left: calc(50%);;
                                                 width: 500px;;
                                                 height: 300px;;
                                                 }
                                                 "
                                            )
                                            )
                                            ),
                                        shiny::tags$style(type="text/css",
                                                          ".shiny-output-error { visibility: hidden; }",
                                                          ".shiny-output-error:before { visibility: hidden; }"
                                        ),
                                        shiny::tags$style(type = "text/css",
                                                          "mark{
                                                          background-color: #ccccff
                                                          }"),

# UI Settings -------------------------------------------------------------

                    selectInput("selectGenome","Select genome",
                                c("Human GRCh38/hg38" = "hg38",
                                  "Human GRCh37/hg19" = "hg19",
                                  "Mouse GRCh38/mm10" = "mm10"),width ="100%"),
                    textInput("textRegion", "Enter genomic coordinate"),
                    actionButton("SubmitBTNRegion","Submit"),
                    br(),


# Conditional Panel with Primer Design Parameters -------------------------

                  conditionalPanel(condition="input.MassArrayTabset == 'Primer Design' && output.clicked_ma=='yes'",
                                   br(),br(),
                                   h4(strong("Primer Design Parameter")),
                                   numericInput("MaPrimS", "Optimal Primer Size",value="20"),
                                   numericInput("MaMinPrimS", "Minimal Primer Size",value="15"),
                                   numericInput("MaMaxPrimS", "Maximal Primer Size",value="24"),
                                   numericInput("MaPrimT", "Optimal Primer melting Temparature",value="55"),
                                   numericInput("MaMinPrimT", "Minimal Primer melting Temparature",value="50"),
                                   numericInput("MaMaxPrimT", "Maximal Primer melting Temparature",value="62"),
                                   numericInput("MaAmpMin", "Minimal Amplicon Size",value="150"),
                                   numericInput("MaAmpMax", "Maximal Amplicon Size",value="300"),
                                   #numericInput("MaPrimer","Numer of Primer to design per region","10"),
                                   checkboxInput("MaExcludeCG","Exclude CGs from Primer",value=TRUE,width="100%"),
                                   uiOutput("cgincluded"),
                                   br(),
                                   actionButton("MA_primerDesign","Submit")

                  ),
                    width=2),


# Main panel --------------------------------------------------------------

                    mainPanel(
                      tabsetPanel(id="MassArrayTabset",

# Home Tab ----------------------------------------------------------------

                                  tabPanel("Home",
                                           textOutput("warning"),
                                           shiny::tags$head(shiny::tags$style(HTML("#warning{color:red;}"))),

                                           h1(strong("Home"),align="center"),
                                           includeMarkdown("MA_main.Rmd"),
                                           actionButton("MA_sampleData","Load Sample Data")

                                           ),

# Input Region Panel ------------------------------------------------------

                                  tabPanel("Input Region",value="inputRegion_tabset",
                                           h1(strong("Input Region"),align="center"),
                                           conditionalPanel(condition="output.clicked_ma=='no'",
                                                            br(),
                                                            column(12,align="center",p("Information about the input region will be displayed here."))),
                                           conditionalPanel(condition="output.clicked_ma=='yes'",
                                                  shiny::tags$head(shiny::tags$style("#inputRegion,#textfw,#textRC{
                                                                         font-size:100%
                                                                         }")),
                                                 br(),
                                                 column(12,align="center",tableOutput("inputMA")),
                                                 br(),
                                                 br(),
                                                 column(12,align="center",textOutput("inputRegion")),
                                                 shiny::tags$head(shiny::tags$style("#inputRegion{color:red}")),
                                                  br(),
                                                  column(12,align="center",textOutput("sliderInfo"),
                                                         br(),
                                                         uiOutput("coordinateSlider")),
                                                  shiny::tags$head(shiny::tags$style(HTML("#coordinateSlider{color:gray;width:100%}"))),
      
                                                  br(),
                                                  br(),
                                                  textOutput("textfw"), #%>% withSpinner(color="#0dc5c1"),
                                                  br(),
                                                  uiOutput("forwardSeq",inline =FALSE),
                                                  # splitLayout(cellWidths = c("90%","50%") ,
                                                  #             #htmlOutput("letterNumberForw",inline=FALSE),
                                                  #             uiOutput("forwardSeq",inline =FALSE),
                                                  #             plotOutput("description",width = "40%")),
                
                                                  shiny::tags$head(shiny::tags$style("#letterNumberForw,letterNumberBack,{
                                                                                     width:100%;
                                                                                     background-color:white;
                                                                                     }")),
                                                  br(),
                                                  textOutput("textRC"),
                                                  br(),
                                                  uiOutput("revC",inline =FALSE),
                                                  br(),
                                                  htmlOutput("textColoring"),
                                                  br(),
                                                  # splitLayout(cellWidths = c("7%","90%","45%"),
                                                  #             htmlOutput("letterNumberBack"),
                                                  #             htmlOutput("revC"),
                                                  #             htmlOutput("suggestionAmplicon")),
                                                  shiny::tags$head(shiny::tags$style("#forwardSeq,#revC {
                                                                                     background-color: white;
                                                                                     width :100%;
                                                                                     } ")))),




# Genomic Features Plot ---------------------------------------------------

                            tabPanel("Genome Annotations",
                                     h1(strong("Genome Annotations"),align="center"),
                                     conditionalPanel(condition="output.clicked_ma=='no'",
                                                  br(),
                                                  column(12,align="center",p("Annotations for the chosen region will be displayed here."))),
                                     conditionalPanel(condition="output.clicked_ma=='yes'",
                                                 plotOutput("PlotOutputs",height ="500px"),#%>% withSpinner(color="#666666"),
            
                                                 br(),br(),
                                                 splitLayout(cellWidths = c("35%","45%"),
                                                             textOutput("CGFWtext")#,
                                                            # textOutput("CGinRCtext")
                                                 ),
                                                 splitLayout(cellWidths = c("35%","35%"),
                                                             tableOutput("CGinFWseq")#,
                                                           #  tableOutput("CGinRC")
                                                 ),
                                                 shiny::tags$head(shiny::tags$style(HTML("#CGinFWseq,#CGinRC,#CGFWasGenomic{;height:500px}"))),
                                                 br(),br()
                                  )),


# Amplicon Prediction Plot ------------------------------------------------
                            tabPanel("Amplicon Prediction",
                                     h1(strong("Amplicon Prediction"),align="center"),
                                     conditionalPanel(condition="output.clicked_ma=='no'",
                                                      br(),
                                                      column(12,align="center",p("Amplicon prediction plots will be displayed here."))),
                                     conditionalPanel(condition="output.clicked_ma=='yes'",
                                                 column(12,align="center",
                                                        textOutput("numberCG"),
                                                        br(),
                                                        p("The following table shows if the CpGs can be assayed by MassArray:"),
                                                        tableOutput("assayCG"),
                                                        textOutput("desciptionPlot"),
                                              plotOutput("amlikonPrediction",width = "90%", height = "600px"), #%>% withSpinner(color="#666666"),
                                              plotOutput("legend",width ="90%")))),


# Primer Design -----------------------------------------------------------
                            tabPanel("Primer Design",
                                  h1(strong("Manual Primer Design"),align="center"),
                                  conditionalPanel(condition="output.clicked_ma=='no'",
                                                   br(),
                                                   column(12,align="center",p("Designed primer pairs will be displayed here."))),
                                  br(),
                                  br(),
                                  conditionalPanel(condition="input.SubmitBTNRegion!=0 && output.clicked_ma=='yes'",
                                  h5(strong("Suggestions for MassArray Primer Design")),
                                  p("The primer melting temperature should be between 52°C and 60°C.",br(),
                                    "The primer should not overlap with CpGs or SNPs. If unavoidable than place those in the 5' position.",br(),
                                    "Design primers in a way that the amplicons contain at least 1 CpG and have a size of 100-500 bp (100-200 for FFPE DNA)", br(),
                                    "After primer design, the forward primer must be preceded at 5’-end by AGGAAGAGAG and the reverse primer at 5’-end by CAGTAATACGACTCACTATAGGGAGAAGCT before ordering"),
                                  br(),
                                  br(),
                                  column(12,align="center",selectInput("standSelector","Select Strand on which to design the primer",c("Watson","Crick")),
                                         br()),
                                  column(12,align="center",uiOutput("slider1"),br()),
                                  column(12,align="center",uiOutput("seqTextForward"),br()),
                                  column(12,align="center",textOutput("temperatureForw"),br(),br()),
                                  column(12,align="center",uiOutput("slider2")),
                                  #column(12,align="center",uiOutput("showRCorC"),br()),
                                  column(12,align="center",uiOutput("seqTextBackward"),br()),
                                  column(12,align="center",textOutput("temperatureBack"),br(),br(),br()),
                                  column(12,align="center",tableOutput("primerInfo")), br(),br(),
                                  column(12,align="center",textOutput("ampsize"),br(),br(),br()),
                                  shiny::tags$head(shiny::tags$style(HTML("#slider2,#slider1{;color:gray;}"))),
                                  shiny::tags$head(shiny::tags$style("#suggestionPrimer{
                                                                     width :80%;height:80%;
                                                                     } ")),
                                  shiny::tags$head(shiny::tags$style("#suggestionAmplicon {
                                                                     width :80%;height:80%;
                                                                     } ")),


# Automated Primer Design by calling Primer3 ------------------------------

                                  h1(strong("Automatic Primer Design"),align="center"),
                                  br(),
                                  textOutput("textPrimerDesign"),
                                  br(),
                                  #tableOutput("primerDesign"),
                                  withSpinner(DT::dataTableOutput("primerDesign"),color="#666666") %>% withSpinner(color="#0dc5c1"),
                                  uiOutput("MassArrayPopup"),
                                  br(),br(),
                                  uiOutput("downloadB"))
                                                  ))))),width=11),



           #####BISULFITAMPLICONSEQUECING#######
           tabPanel("AmpliconBisulfiteSequencing", #changed from BisulfiteAmpliconSequencing
                    add_busy_spinner(spin = "fading-circle",position="full-page"),
                    sidebarLayout(
                      sidebarPanel(
                        selectInput("selectGenomeBisulfit","Select genome",
                                    c("Human GRCh37/hg19" = "hg19",
                                      "Human GRCh38/hg38" = "hg38",
                                      "Mouse GRCm38/mm10" = "mm10"
                                    ),width ="200px"),



                          uiOutput("uiInputSelection"),



                        conditionalPanel(
                          condition = "input.uiInputChoiceBisulfit == 'fasta'",
                          checkboxInput("uiConvert","Bisulfite Convert",value=TRUE,width="100%")#,
                          #numericInput("uiTargetC", paste("Target CG", "(Position of C)", sep="\n"),value="250")
                        ),
                        conditionalPanel(
                          condition = "input.uiInputChoiceBisulfit == 'cpgID'",
                          selectInput("uiEPIC","Select Array type",
                                      c("EPIC-Array" = "EPIC",
                                        "450k-Array" = "450k"
                                      ),width ="200px")
                          #numericInput("uiTargetC", paste("Target CG", "(Position of C)", sep="\n"),value="250")
                        ),


                        #Decisions for Primer Design
                        numericInput("uiPrimS", "Optimal Primer Size",value="20"),
                        numericInput("uiMinPrimS", "Minimal Primer Size",value="15"),
                        numericInput("uiMaxPrimS", "Maximal Primer Size",value="24"),
                        numericInput("uiPrimT", "Optimal Primer melting Temparature",value="55"),
                        numericInput("uiMinPrimT", "Minimal Primer melting Temparature",value="50"),
                        numericInput("uiMaxPrimT", "Maximal Primer melting Temparature",value="62"),
                        numericInput("uiAmpMin", "Minimal Amplicon Size",value="150"),
                        numericInput("uiAmpMax", "Maximal Amplicon Size",value="300"),
                        numericInput("nPrimer","Numer of Primer to design per region","10",min=1,max=999),
                        checkboxInput("uiExcludeCG","Exclude CGs from Primer",value=TRUE,width="100%"),

                        uiOutput("uiBisulfit"),

                        actionButton("advSet","Use advanced settings",width="100%"),
                        br(),
                        br(),


                        actionButton("SubmitBTNRegionBisulfit","Submit")


                      ,width=2),
                      mainPanel(
                        tabsetPanel(id="BSampTabset",

# Home Panel --------------------------------------------------------------
                          tabPanel("Home",
                          h1(strong("Home"),align="center"),
                          includeMarkdown("AMP_main.Rmd"),
                          actionButton("AMP_sampleData","Load Sample Data")
                          ),
                          tabPanel("Designed Primer",value="designedPrimer",
                                   h1(strong("Designed Primer"),align="center"),
                                   conditionalPanel(condition="input.SubmitBTNRegionBisulfit==0",
                                                    br(),
                                                    column(12,align="center",p("Designes primer pairs will be displayed here."))),
                                   #withSpinner(DT::dataTableOutput("test"),color="#666666"),
                                   DT::dataTableOutput("test"),
                                   #add_busy_spinner(spin = "fading-circle",position="full-page"),
                                   uiOutput("popup"),
                                   br(),
                                   uiOutput("downloadBut"),
                                   br(),
                                   br()),
                          tabPanel("Used Settings",
                                   h1(strong("Used Settings"),align="center"),
                                   br(),
                                   column(12,align="center",tableOutput("settings")))
                          )
                        )
                    )
                  ),


#BisBlast Section ------------------------------------------------------------
tabPanel("BisBlast",
         add_busy_spinner(spin = "fading-circle",position="full-page"),
         sidebarLayout(
           sidebarPanel(
           selectInput("selectGenome_blast","Select genome",
                       c("Human GRCh37/hg19" = "hg19",
                         "Human GRCh38/hg38" = "hg38",
                         "Mouse GRCm38/mm10" = "mm10"
                       ),width ="200px"),
           textInput("blast_in","Primer Input"),
           actionButton("submit_blast","Submit"),
           width=2

           ),
           mainPanel(
             tabsetPanel(
               id="bisblast",
               tabPanel("Home",
                 h1(strong("BisBlast"),align="center"),
                 includeMarkdown("blast_main.Rmd")
                 ),
               tabPanel("Query Results",value="blast_query",
                 h1(strong("Results"),align="center"),
                 conditionalPanel(condition="input.submit_blast==0",
                                  br(),
                                  column(12,align="center",p("BisBlast query results will be displayed here."))),
                 DT::dataTableOutput("blast_table")
               )

             )
           )
           )),


# Analysis Pipeline ----------------------------------------------------------

tabPanel("Analysis Pipeline",
         add_busy_spinner(spin = "fading-circle",position="full-page",timeout = 1000),

         ### Specify input options
         sidebarLayout(
           sidebarPanel(#style = "position:fixed;width:inherit;",
             fileInput("input_analysis_sample", "Upload Sample Sheet",
                       multiple = FALSE,accept = ".txt"),
             fileInput("input_analysis_regions", "Upload Regions",
                       multiple = FALSE,accept = ".bed"),
             fileInput("input_analysis_files", "Upload .cov Files",
                       multiple = TRUE,accept = ".cov"),
             numericInput("cutoff","Coverage Cut-Off",value=10,min=1),
             actionButton("submit_analysis","Submit"),
             width = 2),
         mainPanel(
           tabsetPanel(
             id="analysispanel",

             ### Home Panel
             tabPanel("Home",
           column(12,h1(strong("Snakemake Pipeline"),align="center"),
                  includeMarkdown("snakemake.Rmd"),
                  br(),
                  downloadButton("Pipeline_sampleData","Download Pipeline Sample Data"))
             ),

            ### Analysis Panel
           tabPanel("Analysis",value="analysis_selected",
                    column(12,h1(strong("Analyzed AmpBS-Seq Results"),align="center"),
                           conditionalPanel(condition="output.clicked_ampBS=='no'",
                                            br(),
                                            column(12,align="center",p("The AmpBS-Seq analysis results will be displayed here."))
                                            ),

                           conditionalPanel(condition="input.submit_analysis!=0 && output.clicked_ampBS=='yes'",
                                            br(),
                                            column(12,align="center",
                                                   selectInput("select_analysis",
                                                               "",
                                                               choices=c("Overview","Quality Control"="QC",
                                                                         "QC Filtered"="CF",
                                                                         "PCA","Heatmap","Region")))),


                           conditionalPanel(condition="input.submit_analysis!=0 && input.select_analysis=='Overview'",
                                            br(),
                                            column(12,align="center",
                                                   tableOutput("dt_overview"),
                                                   tableOutput("cov_sites"))),

                           conditionalPanel(condition="input.submit_analysis!=0 && input.select_analysis=='QC'",
                                            br(),
                                            column(12,align="center",
                                                   tags$b("Barplot of the total reads per sample (log10 scale):"),
                                                   plotOutput("qc_totalReads"),
                                                   br(),
                                                   tags$b("Barplot showing the number of covered CpG sites per sample:"),
                                                   plotOutput("qc_covSites"),
                                                   br(),
                                                   tags$b("Correlation between coverage and detected CpG sites:"),
                                                   br(),
                                                   HTML("The number of total reads per sample (log10 scale) is plotted on the x-axis against the number of detected CpG sites (y-axis)."),
                                                   plotOutput("qc_cor"),
                                                   br(),
                                                   tags$b("Boxplot showing the number of reads per CpG site (log10 scale):"),
                                                   plotOutput("qc_covCG"))),

                           conditionalPanel(condition="input.submit_analysis!=0 && input.select_analysis=='CF'",
                                            br(),
                                            column(12,align="center",
                                                   tags$b("Barplot showing the number of covered CpG sites per sample after coverage filtering:"),
                                                   plotOutput("cf_covSites"),
                                                   br(),
                                                   tags$b("Correlation between coverage and detected CpG sites after coverage filtering:"),
                                                   br(),
                                                   HTML("The number of total reads per sample (log10 scale) is plotted on the x-axis against the number of detected CpG sites (y-axis)."),
                                                   plotOutput("cf_cor"),
                                                   br(),
                                                   tags$b("Boxplot showing the number of reads per CpG site (log10 scale) after coverage filtering:"),
                                                   plotOutput("cf_covCG"))),

                           conditionalPanel(condition="input.submit_analysis!=0 && input.select_analysis=='PCA'",
                                            br(),
                                            column(12,align="center",
                                                   HTML("Please select which sample annotation column should be used for coloring the samples."),
                                                   uiOutput("select_pca"),
                                                   tags$b("The principal component analysis (PCA) will be computed based on all CpG sites which are covered in all samples:"),
                                                   br(),
                                                   plotOutput("pca_ampBS"))),

                           conditionalPanel(condition="input.submit_analysis!=0 && input.select_analysis=='Heatmap'",
                                            br(),
                                            column(12,align="center",
                                                   uiOutput("select_heatmap"),
                                                   div(style="display: inline-block;vertical-align:top;",uiOutput("heatmap_coloring")),
                                                   br(),
                                                   div(style="display: inline-block;vertical-align:top;",uiOutput("heatmap_rows")),
                                                   div(style="display: inline-block;vertical-align:top;",uiOutput("heatmap_cols")),
                                                   div(style="display: inline-block;vertical-align:top;",uiOutput("heatmap_colNames")),
                                                   div(style="display: inline-block;vertical-align:top;",uiOutput("heatmap_rowNames")),
                                                   br(),
                                                   tags$b("The heatmap will be computed based on all CpG sites which are covered in all samples:"),
                                                   br(),
                                                   plotOutput("heatmap_ampBS"))),
                           
                           conditionalPanel(condition="input.submit_analysis!=0 && input.select_analysis=='Region'",
                                            br(),
                                            column(12,align="center",
                                                   uiOutput("select_region"),
                                                   div(style="display: inline-block;vertical-align:top;",uiOutput("co_region_ui")),
                                                   br(),
                                                   tags$b("Heatmap of all CpG sites in the analyzed regions. By applying a coverage cut-off, CpGs with a low coverage might be discarded:"),
                                                   br(),
                                                   plotOutput("heatmap_region")))


                           )),
           tabPanel("Download",value="download_ampBS_analysis",
                    column(12,h1(strong("Download AmpBS-Seq Results"),align="center"),
                           conditionalPanel(condition="output.clicked_ampBS=='no'",
                                            br(),
                                            column(12,align="center",p("The AmpBS-Seq download report will be available here."))),
                           conditionalPanel(condition="output.clicked_ampBS=='yes'",
                                            br(),
                                            column(12,align="center",
                                                   uiOutput("downloadBut_ampBS")
                                                   ))


                    ))

         )))),


#Help Section ------------------------------------------------------------
      tabPanel("Help",
               sidebarLayout(
                 sidebarPanel(#style = "position:fixed;width:inherit;",
                   selectInput("help_input","Method",choices=c("MassArray"="a","AmpBS-Seq"="b")),width=2),
               mainPanel(

# Help for MassArray ------------------------------------------------------
              conditionalPanel(condition = ("input.help_input=='a'"),
                                column(12,
                                       h1(strong("Tutorial for MassArray Primer Design"),align="center"),
                                       includeMarkdown("MA_intro.Rmd"),
                                       actionButton("interactive_MA","Interactive Sample Data"))
                               ),
              conditionalPanel(condition = ("input.help_input=='b'"),
                               column(12,h1(strong("Tutorial for AmpBS-Seq Primer Design"),align="center"),
                                      includeMarkdown("AMP_intro.Rmd"),
                                      useShinyjs(),
                                      actionButton("interactive_AMP","Interactive Sample Data"))
                               )
                             )

)),





# Write the Impressum -----------------------------------------------------
      tabPanel("Impressum/About",
                 mainPanel(
                      column(12,h1(strong("About AmpliconDesign"),align="center"),
                             includeMarkdown("about.Rmd"))

               ))


)
