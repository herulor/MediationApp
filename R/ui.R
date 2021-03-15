
# # Paradox of converging evidence application
# $ Mediation analyses

# # Illinois colours
IlliniOrange <- "#FF552E"
IlliniDkBlue <- "#13294B"
IlliniBlue   <- "#7FC3E1"
IlliniAltOrg <- "#DD3403"
IlliniArches <- "#009FD4"


# Define UI for random distribution app ----
ui <- fluidPage(

    tags$head(
        includeCSS("www/CSS.css")
    ),

    shinyWidgets::setSliderColor(rep(IlliniOrange, 13), seq(0, 13)),

    p(),
    fluidRow(column(4,
                    img(src = "University-Wordmark-Full-Color-RGB.png", height = 70),
    ),
    column(8,

           # App title ----
           titlePanel(HTML("What does mediation reveal about individual people? <br> Computing Pattern Probabilities")),
    )
    ),
    p(),

    # Sidebar layout with instructions for the App
    sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(width = 4,


                     helpText("Select the tab correspoding to the structure of your",
                              "mediation analysis: Two- or three-step mediation.",
                              p(),
                              "Set the correlations between the variables by",
                              "using the slide bars.",
                              style = paste('color:', IlliniDkBlue)),

                     br(),

                     # Tabset to separate two-step and three-step mediation parts of app
                     tabsetPanel(id = "analysisTabs", type = "tabs", selected = 3,

                                 # # Two-Step Panel
                                 tabPanel(title = HTML("<center>Two-Step Mediation <br> X &rarr; M &rarr; Y</center>"),
                                          value = 1,
                                          h3("Correlations"),

                                          fluidRow(
                                              # # Correlations: Input slide bars
                                              column(7,
                                                     sliderInput(inputId = "XM2Scor",
                                                                 label = "X-M",
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     sliderInput(inputId = "MY2Scor",
                                                                 label = "M-Y",
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     sliderInput(inputId = "XY2Scor",
                                                                 label = "X-Y",
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),
                                              ),

                                              # Input: Parameters for how to cut the distribution and how to plot the sample from the joint distribution.
                                              column(5,
                                                     sliderInput("cutOff2S",
                                                                 "Splitting percentile",
                                                                 value = 50,
                                                                 min = 50,
                                                                 max = 99),

                                                     p(),

                                                     sliderInput("n2S",
                                                                 "Number of points:",
                                                                 value = 50,
                                                                 min = 10,
                                                                 max = 300),

                                                     p(),

                                                     checkboxInput("ellipsoids2S",
                                                                   "Plot ellipsoids"
                                                     ),

                                                     p(),

                                                     checkboxInput("cubes2S",
                                                                   "Plot cubes",
                                                                   value = TRUE
                                                     )

                                              ),
                                          )#,


                                 ),

                                 # # Three-Step Panel
                                 tabPanel(title = HTML("<center>Three-Step Mediation <br> X &rarr; M<sub>1</sub> &rarr; M<sub>2</sub> &rarr; Y</center>"),
                                          value = 2,


                                          h3("Correlations"),

                                          fluidRow(
                                              # # Correlations: Input slide bars
                                              column(6,
                                                     sliderInput(inputId = "XM13Scor",
                                                                 label = "X-M1",
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     sliderInput(inputId = "M1Y3Scor",
                                                                 label = "M1-Y",
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     sliderInput(inputId = "M1M23Scor",
                                                                 label = "M1-M2",
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),
                                              ),

                                              column(6,
                                                     sliderInput(inputId = "XM23Scor",
                                                                 label = "X-M2",
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     sliderInput(inputId = "M2Y3Scor",
                                                                 label = "M2-Y",
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),

                                                     sliderInput(inputId = "XY3Scor",
                                                                 label = "X-Y",
                                                                 min = -1,
                                                                 max =  1,
                                                                 value = 0,
                                                                 step = .01, round = FALSE),
                                              ),

                                              # Input: Parameters for how to cut the distribution and how to plot the sample from the joint distribution.
                                              column(6,
                                                     sliderInput("cutOff3S",
                                                                 "Splitting percentile",
                                                                 value = 50,
                                                                 min = 50,
                                                                 max = 99),

                                                     p(),

                                                     sliderInput("n3S",
                                                                 "Number of points:",
                                                                 value = 50,
                                                                 min = 10,
                                                                 max = 300),
                                              ),

                                              column(5,
                                                     checkboxInput("ellipsoids3S",
                                                                   "Plot ellipsoids"
                                                     ),

                                                     p(),

                                                     checkboxInput("cubes3S",
                                                                   "Plot cubes",
                                                                   value = TRUE
                                                     )

                                              ),
                                          )#,

                                 ),

                                 # # Three-Step Panel
                                 tabPanel(title = HTML("<center>About</center><br>"),
                                          value = 3,

                                          helpText(p(),
                                        "This app is a companion to",
                                          "'Bridging the method-theory gap in mediation analysis:",
                                          "What does mediation reveal about individual people?'",
                                          "(Bogdan, Cervantes & Regenwetter, 2021).",
                                          p(),
                                          ""
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

                                 uiOutput("twoStepXMProbMessage"),

                                 p(),

                                 uiOutput("twoStepMYProbMessage"),

                                 p(),

                                 uiOutput("twoStepXYProbMessage"),

                                 p(),
                                 br(),
                                 p(),

                                 tableOutput("twoStepTable")

                          ),

                          column(6, align = "center",
                                 rgl::rglwidgetOutput("twoStepPlot", width = 550)
                          )
                      ),

                      fluidRow(
                          column(12, align = "center",
                                 textOutput("ellipsoid2SMessage")
                          )
                      ),



                  ),

                  conditionalPanel(
                      condition = "input.analysisTabs == 2",
                      fluidRow(
                          column(6, align = "left",
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

                                 tableOutput("threeStepTable")

                          ),

                          column(6, align = "center",
                                 rgl::rglwidgetOutput("threeStepM1Plot", height = 350, width = 550),
                                 rgl::rglwidgetOutput("threeStepM2Plot", height = 350, width = 550),

                          )
                      ),

                      fluidRow(
                          column(12, align = "center",
                                 textOutput("ellipsoid3SMessage")
                          )
                      ),

                  ),



        )
    )
)

