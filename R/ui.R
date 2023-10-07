# Paradox of converging evidence application ----
# Mediation analyses

# # Illinois colours
IlliniOrange <- "#FF552E"
IlliniDkBlue <- "#13294B"
IlliniBlue   <- "#7FC3E1"
IlliniAltOrg <- "#DD3403"
IlliniArches <- "#009FD4"



## Define UI for app ----
ui <- fluidPage(

    tags$head(
        includeCSS("./www/CSS.css")
    ),

    p(),

    fluidRow(column(4,
                    img(src = "University-Wordmark-Full-Color-RGB.png", height = 70),
    ),
    column(8,

           ## App title ----
           titlePanel(HTML("What does a popuation-level mediation reveal about individual people? <br> Computing Pattern Probabilities")),
    )
    ),
    p(),

    # Sidebar layout with instructions for the App
    sidebarLayout(

        ## Sidebar panel for inputs ----
        sidebarPanel(width = 4,


                     helpText("Select the tab correspoding to the structure of your",
                              "mediation analysis: Two- or three-step mediation.",
                              p(),
                              "Set the correlations between the variables by",
                              "using the slide bars or introducing them into the",
                              "respective text boxes.",
                              style = paste('color:', IlliniDkBlue)),

                     br(),

                     ### Tabset to separate two-step and three-step mediation parts of app ----
                     tabsetPanel(id = "analysisTabs", type = "tabs", selected = 3,

                                 #### Two-Step Panel ----
                                 tabPanel(title = HTML("<center>Two-Step Mediation <br> X &rarr; M &rarr; Y</center>"),
                                          value = 1,
                                          h3("Correlations"),

                                          fluidRow(
                                              # # Correlations: Input slide bars
                                              column(7,
                                                     tags$div(title = "Correlation between X and M",
                                                     textInput(inputId = "XM2ScorText",
                                                               label = "X-M",
                                                               value = 0,
                                                               width = "40%"),
                                                     ),

                                                     sliderInput(inputId = "XM2Scor",
                                                                 label = NULL,
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     tags$div(title = "Correlation between M and Y",
                                                     textInput(inputId = "MY2ScorText",
                                                               label = "M-Y",
                                                               value = 0,
                                                               width = "40%"),
                                                     ),

                                                     sliderInput(inputId = "MY2Scor",
                                                                 label = NULL,
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     tags$div(title = "Correlation between X and Y",
                                                     textInput(inputId = "XY2ScorText",
                                                               label = "X-Y",
                                                               value = 0,
                                                               width = "40%"),
                                                     ),

                                                     sliderInput(inputId = "XY2Scor",
                                                                 label = NULL,
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     tags$div(title = "Resets all correlations to zero",
                                                     actionButton("resetInput2", "Reset correlations"),
                                                     ),

                                              ),

                                              # Input: Parameters for how to cut the distribution and how to plot the sample from the joint distribution.
                                              column(5,
                                                     tags$div(title = "If population correlations are known, set at zero (0).\nIf computed from a random sample, write the sample size.\nMarginal confidence intervals will be computed for the pattern proportions\nwhen sample size is greater than zero.",
                                                     textInput(inputId = "n2SText",
                                                               label = "Sample size:",
                                                               value = 0,
                                                               width = "70%"),
                                                     ),

                                                     p(),

                                                     tags$div(title = "Number of points to simulate for the figure on the right.\nUsed for illustration purposes only.",
                                                     textInput(inputId = "d2SText",
                                                               label = "Number of points:",
                                                               value = 100,
                                                               width = "70%"),
                                                     ),

                                                     p(),

                                                     tags$div(title = "Percentile of the cutoff point for 'high' groups.",
                                                     textInput(inputId = "cutOff2SText",
                                                               label = "Splitting percentile",
                                                               value = 50,
                                                               width = "50%"),
                                                     ),

                                                     sliderInput("cutOff2S",
                                                                 label = NULL,
                                                                 value = 50,
                                                                 min = 50,
                                                                 max = 99),

                                                     p(),

                                                     sliderInput("digits2S",
                                                                 "Number of decimal places for proportions:",
                                                                 value = 3,
                                                                 min = 2,
                                                                 max = 5),

                                                     p(),

                                                     tags$div(title = "Check to plot the ellipse that envelopes 95% of the multivariate normal distribution.",
                                                     checkboxInput("ellipsoids2S",
                                                                   "Plot ellipsoids"
                                                     ),
                                                     ),

                                                     p(),

                                                     tags$div(title = "Check to plot the cubes enclosing the No violation patterns.",
                                                     checkboxInput("cubes2S",
                                                                   "Plot cubes",
                                                                   value = TRUE
                                                     )
                                                     ),

                                              ),
                                          )#,


                                 ),

                                 #### Three-Step Panel ----
                                 tabPanel(title = HTML("<center>Three-Step Mediation <br> X &rarr; M<sub>1</sub> &rarr; M<sub>2</sub> &rarr; Y</center>"),
                                          value = 2,


                                          h3("Correlations"),

                                          fluidRow(
                                              # # Correlations: Input slide bars
                                              column(6,
                                                     tags$div(title = HTML("Correlation between X and M1"),
                                                     textInput(inputId = "XM13ScorText",
                                                               label = "X-M1",
                                                               value = 0,
                                                               width = "40%"),
                                                     ),

                                                     sliderInput(inputId = "XM13Scor",
                                                                 label = NULL,
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     tags$div(title = HTML("Correlation between M1 and Y"),
                                                     textInput(inputId = "M1Y3ScorText",
                                                               label = "M1-Y",
                                                               value = 0,
                                                               width = "40%"),
                                                     ),

                                                     sliderInput(inputId = "M1Y3Scor",
                                                                 label = NULL,
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     tags$div(title = HTML("Correlation between M1 and M2"),
                                                     textInput(inputId = "M1M23ScorText",
                                                               label = "M1-M2",
                                                               value = 0,
                                                               width = "40%"),
                                                     ),

                                                     sliderInput(inputId = "M1M23Scor",
                                                                 label = NULL,
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),
                                              ),

                                              column(6,
                                                     tags$div(title = HTML("Correlation between X and M2"),
                                                     textInput(inputId = "XM23ScorText",
                                                               label = "X-M2",
                                                               value = 0,
                                                               width = "40%"),
                                                     ),

                                                     sliderInput(inputId = "XM23Scor",
                                                                 label = NULL,
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     tags$div(title = HTML("Correlation between M2 and Y"),
                                                     textInput(inputId = "M2Y3ScorText",
                                                               label = "M2-Y",
                                                               value = 0,
                                                               width = "40%"),
                                                     ),

                                                     sliderInput(inputId = "M2Y3Scor",
                                                                 label = NULL,
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     tags$div(title = HTML("Correlation between X and Y"),
                                                     textInput(inputId = "XY3ScorText",
                                                               label = "X-Y",
                                                               value = 0,
                                                               width = "40%"),
                                                     ),

                                                     sliderInput(inputId = "XY3Scor",
                                                                 label = NULL,
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),
                                              ),

                                              # Input: Parameters for how to cut the distribution and how to plot the sample from the joint distribution.
                                              column(6,
                                                     tags$div(title = "If population correlations are known, set at zero (0).\nIf computed from a random sample, write the sample size.\nMarginal confidence intervals will be computed for the pattern proportions\nwhen sample size is greater than zero.",
                                                     textInput(inputId = "n3SText",
                                                               label = "Sample size:",
                                                               value = 0,
                                                               width = "70%"),
                                                     ),

                                                     p(),

                                                     tags$div(title = "Number of points to simulate for the figure on the right.\nUsed for illustration purposes only.",
                                                     textInput(inputId = "d3SText",
                                                               label = "Number of points:",
                                                               value = 100,
                                                               width = "70%"),
                                                     ),

                                                     p(),

                                                     sliderInput("digits3S",
                                                                 "Number of decimal places for proportions:",
                                                                 value = 3,
                                                                 min = 2,
                                                                 max = 5),
                                              ),

                                              column(5,

                                                     tags$div(title = "Percentile of the cutoff point for 'high' groups.",
                                                     textInput(inputId = "cutOff3SText",
                                                               label = "Splitting percentile",
                                                               value = 50,
                                                               width = "50%"),
                                                     ),

                                                     sliderInput("cutOff3S",
                                                                 label = NULL,
                                                                 value = 50,
                                                                 min = 50,
                                                                 max = 99),

                                                     p(),


                                                     tags$div(title = "Check to plot the ellipse that envelopes 95% of the multivariate normal distribution.",
                                                     checkboxInput("ellipsoids3S",
                                                                   "Plot ellipsoids"
                                                     ),
                                                     ),

                                                     p(),

                                                     tags$div(title = "Check to plot the cubes enclosing the No violation patterns.",
                                                     checkboxInput("cubes3S",
                                                                   "Plot cubes",
                                                                   value = TRUE
                                                     ),
                                                     ),

                                                     tags$div(title = "Resets all correlations to zero",
                                                     actionButton("resetInput3", "Reset correlations"),
                                                     ),

                                              ),
                                          )#,

                                 ),

                                 # # About Panel
                                 tabPanel(title = HTML("<center>About</center><br>"),
                                          value = 3,

                                          helpText(p(),
                                                   "This app is a companion to",
                                                   "'Bridging the model-theory gap in mediation analysis:",
                                                   "What does a population-level mediation reveal about individual people?'",
                                                   "(Bogdan, Cervantes & Regenwetter, 2023).",
                                                   p(),
                                                   "You may use it to compute the probabilities for patterns",
                                                   "of high/low values across three variables involved in a",
                                                   "mediation analysis (on tab Two-Step Mediation)",
                                                   "of across four variables in a mediation chaing",
                                                   "(on tab Three-Step Mediation)."
                                          )

                                 )

                     )

        ),

        # Main panel
        mainPanel(width = 8,

                  conditionalPanel(
                      condition = "input.analysisTabs == 1",

                      fluidRow(
                          column(6, align = "left",
                                 p(),

                                 uiOutput("twoStepPartialMessage"),

                                 p(),
                                 br(),
                                 p(),

                                 uiOutput("twoStepXMProbMessage"),

                                 p(),

                                 uiOutput("twoStepMYProbMessage"),

                                 p(),

                                 uiOutput("twoStepXYProbMessage"),

                                 p(),
                                 br(),
                                 p(),

                                 conditionalPanel("output.twoStepMatrixCaption != ''",
                                 selectInput('tableId_2S', 'Choose proportion table to show',
                                             c("Actual proportions",
                                               "Implied by Full Mediation",
                                               "Implied by No Mediation")),
                                 ),

                                 conditionalPanel("input.tableId_2S == 'Actual proportions'",
                                   tableOutput("twoStepTableAct")
                                 ),

                                 conditionalPanel("input.tableId_2S == 'Implied by Full Mediation'",
                                   tableOutput("twoStepTableFull")
                                 ),

                                 conditionalPanel("input.tableId_2S == 'Implied by No Mediation'",
                                                  tableOutput("twoStepTableNo")
                                 ),

                          ),

                          column(6, align = "right",
                                 rgl::rglwidgetOutput("twoStepPlot", width = 550),
                          )
                      ),

                      fluidRow(
                          p(),
                          column(12, align = "center", htmlOutput("twoStepMatrixCaption")),

                          column(4, align = "center",
                                 htmlOutput("twoStepCorMatrixCaption"),
                                 plotOutput("twoStepCorMatrix", width = 220)
                          ),
                          column(4, align = "center",
                                 htmlOutput("twoStepFullMatrixCaption"),
                                 plotOutput("twoStepFullMedMatrix", width = 220)
                          ),
                          column(4, align = "center",
                                 htmlOutput("twoStepNoMatrixCaption"),
                                 plotOutput("twoStepNoMedMatrix", width = 220)
                          ),
                          p(),
                          column(12, align = "center",
                                 textOutput("ellipsoid2SMessage"),
                          )
                      )
                  ),


        conditionalPanel(
            condition = "input.analysisTabs == 2",
            fluidRow(
                column(6, align = "left",
                       p(),

                       uiOutput("threeStepPartialMessage"),

                       p(),
                       br(),
                       p(),

                       uiOutput("threeStepXM1ProbMessage"),

                       p(),

                       uiOutput("threeStepM1M2ProbMessage"),

                       p(),

                       uiOutput("threeStepM2YProbMessage"),

                       p(),

                       uiOutput("threeStepXYProbMessage"),

                       p(),
                       br(),
                       p(),

                       conditionalPanel("output.twoStepMatrixCaption != ''",
                                        selectInput('tableId_3S', 'Choose proportion table to show',
                                                    c("Actual proportions",
                                                      "Implied by Full Mediation",
                                                      "Implied by No Mediation")),
                       ),

                        conditionalPanel("input.tableId_3S == 'Actual proportions'",
                                         tableOutput("threeStepTableAct")
                        ),

                        conditionalPanel("input.tableId_3S == 'Implied by Full Mediation'",
                                         tableOutput("threeStepTableFull")
                        ),

                       conditionalPanel("input.tableId_3S == 'Implied by No Mediation'",
                                        tableOutput("threeStepTableNo")
                       ),

                ),

                column(6, align = "right",
                       rgl::rglwidgetOutput("threeStepM1Plot", height = 350, width = 550),
                       rgl::rglwidgetOutput("threeStepM2Plot", height = 350, width = 550),

                )
            ),

            fluidRow(
                          p(),
                          column(12, align = "center", htmlOutput("threeStepMatrixCaption")),
            ),

            fluidRow(
                          column(4, align = "center",
                                 htmlOutput("threeStepCorMatrixCaption"),
                                 plotOutput("threeStepCorMatrix", width = 220)
                          ),
                          column(4, align = "center",
                                 htmlOutput("threeStepFullMatrixCaption"),
                                 plotOutput("threeStepFullMedMatrix", width = 220)
                          ),
                          column(4, align = "center",
                                 htmlOutput("threeStepNoMatrixCaption"),
                                 plotOutput("threeStepNoMedMatrix", width = 220)
                          ),
                          p(),
                 column(12, align = "center",
                       textOutput("ellipsoid3SMessage")
                )
            ),

        ),



    )
)
)

