# Must be executed BEFORE rgl is loaded on headless devices.
options(rgl.useNULL = TRUE)

library(shiny)
library(shinyWidgets)
library(tidyverse)
library(rgl)
library(kableExtra)
library(corpcor)
library(mvtnorm)
library(corrplot)


# # Illinois colours
IlliniOrange <- "#FF552E"
IlliniDkBlue <- "#13294B"
IlliniBlue   <- "#7FC3E1"
IlliniAltOrg <- "#DD3403"
IlliniArches <- "#009FD4"


corMatricesBgPalette <- colorRampPalette(c(IlliniAltOrg, "white", IlliniBlue))(201)


# Define server logic for random distribution app ----
server <- function(input, output, session) {

  # Digits for output
  digits2S <- reactive(input$digits2S)
  format2S <- reactive(paste0("%0.", digits2S(), "f"))

  # Update labels with correlations according with the respective slideBars ----
  # Two-step
  observeEvent(input$XM2Scor, {
    label <- paste0("X-M = ", sprintf("%.2f", input$XM2Scor))
    updateSliderInput(session, "XM2Scor", label = label)
  })

  observeEvent(input$MY2Scor, {
    label <- paste0("M-Y = ", sprintf("%.2f", input$MY2Scor))
    updateSliderInput(session, "MY2Scor", label = label)
  })

  observeEvent(input$XY2Scor, {
    label <- paste0("X-Y = ", sprintf("%.2f", input$XY2Scor))
    updateSliderInput(session, "XY2Scor", label = label)
  })

  # Three-step
  observeEvent(input$XY3Scor, {
    label <- paste0("X-Y = ", sprintf("%.2f", input$XY3Scor))
    updateSliderInput(session, "XY3Scor", label = label)
  })

  observeEvent(input$XM13Scor, {
    label <- paste0("X-M1 = ", sprintf("%.2f", input$XM13Scor))
    updateSliderInput(session, "XM13Scor", label = label)
  })

  observeEvent(input$M1Y3Scor, {
    label <- paste0("M1-Y = ", sprintf("%.2f", input$M1Y3Scor))
    updateSliderInput(session, "M1Y3Scor", label = label)
  })

  observeEvent(input$XM23Scor, {
    label <- paste0("X-M2 = ", sprintf("%.2f", input$XM23Scor))
    updateSliderInput(session, "XM23Scor", label = label)
  })

  observeEvent(input$M2Y3Scor, {
    label <- paste0("M2-Y = ", sprintf("%.2f", input$M2Y3Scor))
    updateSliderInput(session, "M2Y3Scor", label = label)
  })

  observeEvent(input$M1M23Scor, {
    label <- paste0("M1-M2 = ", sprintf("%.2f", input$M1M23Scor))
    updateSliderInput(session, "M1M23Scor", label = label)
  })



  # # Two-Step mediation computations ----


  # Reactive reading of the correlations and the percentile for defining the top (bottom) groups
  #
  quantileGroup2S <- reactive(qnorm(input$cutOff2S / 100))
  corXM2S <- reactive(input$XM2Scor)
  corMY2S <- reactive(input$MY2Scor)
  corXY2S <- reactive(input$XY2Scor)


  # Reactive expression to check if the matrix is positive definite
  twoStepMatrix <- reactive({
    mat <- matrix(c(1,         corXM2S(), corXY2S(),
                    corXM2S(), 1,         corMY2S(),
                    corXY2S(), corMY2S(), 1),

                  nrow = 3)
  })

  twoStepFullMedMatrix <- reactive({
    matrix(c(1,                     corXM2S(), corXM2S() * corMY2S(),
             corXM2S(),             1,         corMY2S(),
             corXM2S() * corMY2S(), corMY2S(), 1),

           nrow = 3)
  })

  twoStepNoMedMatrix <- reactive({
    matrix(c(1,         0, corXY2S(),
             0,         1, 0,
             corXY2S(), 0, 1),

           nrow = 3)
  })

  twoStepPcor <- reactive(cor2pcor(twoStepMatrix()))

  pcorXY2S <- reactive(twoStepPcor()[1, 3])


  twoStepSimData <- reactive({
    mvtnorm::rmvnorm(n = input$n2S,
                     mean = c(0, 0, 0),
                     sigma = twoStepMatrix())
  })

  isTwoStepSemiPosDef <- reactive(all(eigen(twoStepMatrix())$values > -10 * .Machine$double.eps))
  isTwoStepPosDef     <- reactive(all(eigen(twoStepMatrix())$values >  10 * .Machine$double.eps))


  ### Proportions Two-Step ------

  propXM2S <- reactive({
    if (quantileGroup2S() == 0) {
      propXM2S <- .25 + (asin(abs(corXM2S()))) / (2 * pi)
    } else {
      propXM2S <- mvtnorm::pmvnorm(lower = c(quantileGroup2S(), quantileGroup2S()),
                                   upper = c(Inf, Inf),
                                   mean = c(0, 0),
                                   corr = abs(twoStepMatrix()[c(1, 2), c(1, 2)]))
    }
    return(propXM2S)
  })

  propMY2S <- reactive({
    if (quantileGroup2S() == 0) {
      propMY2S <- .25 + (asin(abs(corMY2S()))) / (2 * pi)
    } else {
      propMY2S <- mvtnorm::pmvnorm(lower = c(quantileGroup2S(), quantileGroup2S()),
                                   upper = c(Inf, Inf),
                                   mean = c(0, 0),
                                   corr = abs(twoStepMatrix()[c(2, 3), c(2, 3)]))
    }
    return(propMY2S)
  })

  corXY2SSign <- reactive({
    ifelse(sign(corXM2S()) * sign(corMY2S()) == 0,
           1,
           sign(corXM2S()) * sign(corMY2S())
    )
  })

  propXY2S <- reactive({
    if (quantileGroup2S() == 0) {
      propXY2S <- .25 + (asin(corXY2SSign() * (corXY2S()))) / (2 * pi)
    } else {
      propXY2S <- mvtnorm::pmvnorm(lower = c(quantileGroup2S(), quantileGroup2S()),
                                   upper = c(Inf, Inf),
                                   mean = c(0, 0),
                                   corr = matrix(c(1, corXY2SSign() * corXY2S(),
                                                   corXY2SSign() * corXY2S(), 1), nrow = 2))
    }
    return(propXY2S)
  })

  qropXY2S <- reactive({
    if (quantileGroup2S() == 0) {
      propXY2S <- .25 + (asin(-corXY2SSign() * (corXY2S()))) / (2 * pi)
    } else {
      propXY2S <- mvtnorm::pmvnorm(lower = c(quantileGroup2S(), quantileGroup2S()),
                                   upper = c(Inf, Inf),
                                   mean = c(0, 0),
                                   corr = matrix(c(1, -corXY2SSign() * corXY2S(),
                                                   -corXY2SSign() * corXY2S(), 1), nrow = 2))
    }
    return(propXY2S)
  })

  propXYStar2S <- reactive({
    if (quantileGroup2S() == 0) {
      propXY2S <- .25 + (asin(corXY2SSign() * (corXM2S() * corMY2S()))) / (2 * pi)
    } else {
      propXY2S <- mvtnorm::pmvnorm(lower = c(quantileGroup2S(), quantileGroup2S()),
                                   upper = c(Inf, Inf),
                                   mean = c(0, 0),
                                   corr = matrix(c(1, corXY2SSign() * corXM2S() * corMY2S(),
                                                   corXY2SSign() * corXM2S() * corMY2S(), 1), nrow = 2))
    }
    return(propXY2S)
  })


  # Message for partial correlation
  output$twoStepPartialMessage <- reactive({
    if (isTwoStepSemiPosDef()) {
      textPartial2S <- paste0("The partial correlation of X and Y is ",
                              strong(sprintf(format2S(), round(pcorXY2S(), digits2S()))))
    } else {
      validate("The correlation matrix is not semi-positive definite. Change some correlation values.")
    }

    return(textPartial2S)
  })


  # Messages for the proportions of 'satisfaction'
  #
  # X-M pair
  output$twoStepXMProbMessage <- reactive({

    if (corXM2S() >= 0) {
      textXM2S <- paste0("The proportion of observations that are on the top ",
                         100 - input$cutOff2S, "% on both X and M is ")
    } else {
      textXM2S <- paste0("The proportion of observations that are on the top ",
                         100 - input$cutOff2S, "% on X and the bottom ",
                         100 - input$cutOff2S, "% on M is ")
    }

    textXM2S <- paste0(textXM2S, strong(sprintf(format2S(), round(propXM2S(), digits2S()))))

    if (!isTwoStepSemiPosDef()) {
      textXM2S <- ""
    }

    return(textXM2S)
  })


  # M-Y pair
  output$twoStepMYProbMessage <- reactive({

    if (corMY2S() >= 0) {
      textMY2S <- paste0("The proportion of observations that are on the top ",
                         100 - input$cutOff2S, "% on both Y and M is ")
    } else {
      textMY2S <- paste0("The proportion of observations that are on the top ",
                         100 - input$cutOff2S, "% on M and the bottom ",
                         100 - input$cutOff2S, "% on Y is ")
    }

    textMY2S <- paste0(textMY2S, strong(sprintf(format2S(), round(propMY2S(), digits2S()))))

    if (!isTwoStepSemiPosDef()) {
      textMY2S <- ""
    }

    return(textMY2S)
  })


  # X-Y pair
  output$twoStepXYProbMessage <- reactive({

    if (corXY2SSign() >= 0) {
      textXY2S <- paste0("The proportion of observations that are on the top ",
                         100 - input$cutOff2S, "% on both X and Y is ")
    } else {
      textXY2S <- paste0("The proportion of observations that are on the top ",
                         100 - input$cutOff2S, "% on X and the bottom ",
                         100 - input$cutOff2S, "% on Y is ")
    }

    textXY2S <- paste0(textXY2S, strong(sprintf(format2S(), round(propXY2S(), digits2S()))))

    if (!isTwoStepSemiPosDef()) {
      textXY2S <- ""
    }

    return(textXY2S)
  })


  # Actual proportions in distribution -----
  twoStepDiagonal <- reactive({
    if (quantileGroup2S() == 0) {
      propMediationDiagonal <- propXM2S() + propMY2S() + propXY2S() - .5
    } else {
      signXM <- sign(corXM2S())
      signMY <- sign(corMY2S())

      lower1 <- c(-Inf, ifelse(signXM >= 0, -Inf, quantileGroup2S()), ifelse((signXM * signMY) >= 0, -Inf, quantileGroup2S()))
      upper1 <- c(-quantileGroup2S(), ifelse(signXM < 0, Inf, -quantileGroup2S()), ifelse((signXM * signMY) < 0, Inf, -quantileGroup2S()))

      propMediationDiagonal <- 2 * mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                                    mean = c(0, 0, 0),
                                                    corr = twoStepMatrix())[1]
    }

    if (!isTwoStepSemiPosDef()) {
      propMediationDiagonal <- NA
    }

    return(propMediationDiagonal)
  })

  twoStepBadXM <- reactive({
    if (quantileGroup2S() == 0) {
      propMediationDiagonal <- -propXM2S() + propMY2S() - propXY2S() + .5
    } else {
      signXM <- sign(corXM2S())
      signMY <- sign(corMY2S())

      lower1 <- c(-Inf, ifelse(signXM < 0, -Inf, quantileGroup2S()), ifelse((signXM * signMY) < 0, -Inf, quantileGroup2S()))
      upper1 <- c(-quantileGroup2S(), ifelse(signXM >= 0, Inf, -quantileGroup2S()), ifelse((signXM * signMY) >= 0, Inf, -quantileGroup2S()))

      propMediationDiagonal <- 2 * mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                                    mean = c(0, 0, 0),
                                                    corr = twoStepMatrix())[1]
    }

    if (!isTwoStepSemiPosDef()) {
      propMediationDiagonal <- NA
    }

    return(propMediationDiagonal)
  })

  twoStepBadMY <- reactive({
    if (quantileGroup2S() == 0) {
      propMediationDiagonal <- propXM2S() - propMY2S() - propXY2S() + .5
    } else {
      signXM <- sign(corXM2S())
      signMY <- sign(corMY2S())

      lower1 <- c(-Inf, ifelse(signXM >= 0, -Inf, quantileGroup2S()), ifelse((signXM * signMY) < 0, -Inf, quantileGroup2S()))
      upper1 <- c(-quantileGroup2S(), ifelse(signXM < 0, Inf, -quantileGroup2S()), ifelse((signXM * signMY) >= 0, Inf, -quantileGroup2S()))

      propMediationDiagonal <- 2 * mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                                    mean = c(0, 0, 0),
                                                    corr = twoStepMatrix())[1]
    }

    if (!isTwoStepSemiPosDef()) {
      propMediationDiagonal <- NA
    }

    return(propMediationDiagonal)
  })

  twoStepBadXMY <- reactive({
    if (quantileGroup2S() == 0) {
      propMediationDiagonal <- -propXM2S() - propMY2S() + propXY2S() + .5
    } else {
      signXM <- sign(corXM2S())
      signMY <- sign(corMY2S())

      lower1 <- c(-Inf, ifelse(signXM < 0, -Inf, quantileGroup2S()), ifelse((signXM * signMY) >= 0, -Inf, quantileGroup2S()))
      upper1 <- c(-quantileGroup2S(), ifelse(signXM >= 0, Inf, -quantileGroup2S()), ifelse((signXM * signMY) < 0, Inf, -quantileGroup2S()))

      propMediationDiagonal <- 2 * mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                                    mean = c(0, 0, 0),
                                                    corr = twoStepMatrix())[1]
    }

    if (!isTwoStepSemiPosDef()) {
      propMediationDiagonal <- NA
    }

    return(propMediationDiagonal)
  })


  # Proportion according to mediation
  twoStepMediation<- reactive({
    if (quantileGroup2S() == 0) {
      propMediation <- propXM2S() + propMY2S() + propXYStar2S() - .5
    } else {
      signXM <- sign(corXM2S())
      signMY <- sign(corMY2S())

      lower1 <- c(-Inf, ifelse(signXM >= 0, -Inf, quantileGroup2S()), ifelse((signXM * signMY) >= 0, -Inf, quantileGroup2S()))
      upper1 <- c(-quantileGroup2S(), ifelse(signXM < 0, Inf, -quantileGroup2S()), ifelse((signXM * signMY) < 0, Inf, -quantileGroup2S()))
      corMat <- twoStepMatrix()
      corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]

      propMediation <- 2 * mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                            mean = c(0, 0, 0),
                                            corr = corMat)[1]
    }

    if (!isTwoStepSemiPosDef()) {
      propMediation <- NA
    }

    return(propMediation)
  })

  twoStepMediationBadXM <- reactive({
    if (quantileGroup2S() == 0) {
      propMediation <- -propXM2S() + propMY2S() - propXYStar2S() + .5
    } else {
      signXM <- sign(corXM2S())
      signMY <- sign(corMY2S())

      lower1 <- c(-Inf, ifelse(signXM < 0, -Inf, quantileGroup2S()), ifelse((signXM * signMY) < 0, -Inf, quantileGroup2S()))
      upper1 <- c(-quantileGroup2S(), ifelse(signXM >= 0, Inf, -quantileGroup2S()), ifelse((signXM * signMY) >= 0, Inf, -quantileGroup2S()))
      corMat <- twoStepMatrix()
      corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]

      propMediation <- 2 * mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                            mean = c(0, 0, 0),
                                            corr = corMat)[1]
    }

    if (!isTwoStepSemiPosDef()) {
      propMediation <- NA
    }

    return(propMediation)
  })

  twoStepMediationBadMY <- reactive({
    if (quantileGroup2S() == 0) {
      propMediation <- propXM2S() - propMY2S() - propXYStar2S() + .5
    } else {
      signXM <- sign(corXM2S())
      signMY <- sign(corMY2S())

      lower1 <- c(-Inf, ifelse(signXM >= 0, -Inf, quantileGroup2S()), ifelse((signXM * signMY) < 0, -Inf, quantileGroup2S()))
      upper1 <- c(-quantileGroup2S(), ifelse(signXM < 0, Inf, -quantileGroup2S()), ifelse((signXM * signMY) >= 0, Inf, -quantileGroup2S()))
      corMat <- twoStepMatrix()
      corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]

      propMediation <- 2 * mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                            mean = c(0, 0, 0),
                                            corr = corMat)[1]
    }

    if (!isTwoStepSemiPosDef()) {
      propMediation <- NA
    }

    return(propMediation)
  })

  twoStepMediationBadXMY <- reactive({
    if (quantileGroup2S() == 0) {
      propMediation <- -propXM2S() - propMY2S() + propXYStar2S() + .5
    } else {
      signXM <- sign(corXM2S())
      signMY <- sign(corMY2S())

      lower1 <- c(-Inf, ifelse(signXM < 0, -Inf, quantileGroup2S()), ifelse((signXM * signMY) >= 0, -Inf, quantileGroup2S()))
      upper1 <- c(-quantileGroup2S(), ifelse(signXM >= 0, Inf, -quantileGroup2S()), ifelse((signXM * signMY) < 0, Inf, -quantileGroup2S()))
      corMat <- twoStepMatrix()
      corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]

      propMediation <- 2 * mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                            mean = c(0, 0, 0),
                                            corr = corMat)[1]
    }

    if (!isTwoStepSemiPosDef()) {
      propMediation <- NA
    }

    return(propMediation)
  })




  # Generate two-step output Table

  checkmark <- HTML('&#10003;')
  otimes    <- '&#8855;'

  output$twoStepTable <- function() ({

    table2S <- data.frame(Violations = c('No violation', 'No violation', 'XM', 'XM', 'MY', 'MY', 'XM; MY', 'XM; MY'),
                          X = rep(c(checkmark, otimes), 4),
                          M = rep(c(checkmark, otimes, otimes, checkmark), 2),
                          Y = c(checkmark, otimes, otimes, checkmark, otimes, checkmark, checkmark, otimes),
                          Actual = c(rep(sprintf(format2S(), round(.5 * twoStepDiagonal(), digits2S())), 2),
                                     rep(sprintf(format2S(), round(.5 * twoStepBadXM(), digits2S())), 2),
                                     rep(sprintf(format2S(), round(.5 * twoStepBadMY(), digits2S())), 2),
                                     rep(sprintf(format2S(), round(.5 * twoStepBadXMY(), digits2S())), 2)),
                          Complete = c(rep(sprintf(format2S(), round(.5 * twoStepMediation(), digits2S())), 2),
                                       rep(sprintf(format2S(), round(.5 * twoStepMediationBadXM(), digits2S())), 2),
                                       rep(sprintf(format2S(), round(.5 * twoStepMediationBadMY(), digits2S())), 2),
                                       rep(sprintf(format2S(), round(.5 * twoStepMediationBadXMY(), digits2S())), 2)),
                          NoM = c(rep(sprintf(format2S(), round((1 - (input$cutOff2S / 100)) * propXY2S(), digits2S())), 2),
                                  rep(sprintf(format2S(), round((1 - (input$cutOff2S / 100)) * qropXY2S(), digits2S())), 2),
                                  rep(sprintf(format2S(), round((1 - (input$cutOff2S / 100)) * qropXY2S(), digits2S())), 2),
                                  rep(sprintf(format2S(), round((1 - (input$cutOff2S / 100)) * propXY2S(), digits2S())), 2)))


    if (quantileGroup2S() != 0) {

      table2S[nrow(table2S) + 1, -c(1:4)] <- sprintf(format2S(),
                                                     1 - colSums(apply(table2S[, -c(1:4)], 2, as.numeric)))
      table2S[nrow(table2S), 1:4] <- c('Other', rep('-', 3))

    }

    table2S <- kable(table2S, align = 'c', format = 'html',
                     escape = FALSE,
                     col.names = c(
                       'Pattern Violations',
                       'X', 'M', 'Y',
                       'Actual Proportion',
                       'Full Mediation',
                       'No Mediation'
                     )
    ) %>%
      kable_styling(bootstrap_options = c("striped", "responsive"),
                    full_width = TRUE,
                    position = "left", font_size = 14)%>%
      column_spec(1, bold = TRUE) %>%
      collapse_rows(columns = 1, valign = "middle") ->
      table2S



    if (!isTwoStepSemiPosDef()) {
      table2S <- ""
    }

    return(table2S)
  }
  )



  # Generate plots: X -> M -> Y ----

  output$twoStepMatrixCaption <- reactive({
    if (!isTwoStepSemiPosDef()) {
      outText <- ""
    } else {
      outText <- as.character(h4(strong("Correlation matrices")))
    }
    return(outText)
  })

  output$twoStepCorMatrixCaption <- reactive({
    if (!isTwoStepSemiPosDef()) {
      outText <- ""
    } else {
      outText <- as.character(strong("Actual", br(), "correlations"))
    }
    return(outText)
  })

  output$twoStepFullMatrixCaption <- reactive({
    if (!isTwoStepSemiPosDef()) {
      outText <- ""
    } else {
      outText <- as.character(strong("Implied by", br(), "full mediation"))
    }
    return(outText)
  })

  output$twoStepNoMatrixCaption <- reactive({
    if (!isTwoStepSemiPosDef()) {
      outText <- ""
    } else {
      outText <- as.character(strong("Implied by", br(), "no mediation"))
    }
    return(outText)
  })

  output$twoStepCorMatrix <- renderPlot({
    if (!isTwoStepSemiPosDef()) {
    } else {
      mat <- twoStepMatrix()
      rownames(mat) <- colnames(mat) <- c("X", "M", "Y")
      corrplot::corrplot(corr = mat,
                         type = "lower", method = "color",
                         diag = FALSE,
                         number.cex = 1.5,
                         addCoef.col = IlliniDkBlue,
                         addgrid.col = IlliniDkBlue,
                         cl.pos = FALSE,
                         tl.cex = 1.4,
                         tl.offset = 1,
                         tl.srt = 0,
                         tl.col = IlliniDkBlue,
                         col = corMatricesBgPalette)
    }
  })

  output$twoStepFullMedMatrix <- renderPlot({
    if (!isTwoStepSemiPosDef()) {
    } else {
      mat <- twoStepFullMedMatrix()
      rownames(mat) <- colnames(mat) <- c("X", "M", "Y")
      corrplot::corrplot(corr = mat,
                         type = "lower", method = "color",
                         diag = FALSE,
                         number.cex = 1.5,
                         addCoef.col = IlliniDkBlue,
                         addgrid.col = IlliniDkBlue,
                         cl.pos = FALSE,
                         tl.cex = 1.4,
                         tl.offset = 1,
                         tl.srt = 0,
                         tl.col = IlliniDkBlue,
                         col = corMatricesBgPalette)
    }
   })

  output$twoStepNoMedMatrix <- renderPlot({
     if (!isTwoStepSemiPosDef()) {
    } else {
      mat <- twoStepNoMedMatrix()
      rownames(mat) <- colnames(mat) <- c("X", "M", "Y")
      corrplot::corrplot(corr = mat,
                         type = "lower", method = "color",
                         diag = FALSE,
                         number.cex = 1.5,
                         addCoef.col = IlliniDkBlue,
                         addgrid.col = IlliniDkBlue,
                         cl.pos = FALSE,
                         tl.cex = 1.4,
                         tl.offset = 1,
                         tl.srt = 0,
                         tl.col = IlliniDkBlue,
                         col = corMatricesBgPalette)
    }
  })


  output$twoStepPlot <- renderRglwidget({

    if (!isTwoStepSemiPosDef()) {
    } else {
      close3d()
      plotData <- twoStepSimData()

      if ((sign(corXM2S()) >= 0) & (sign(corMY2S()) >= 0)) {
        colored <- (plotData[, 1] > quantileGroup2S() &
                      plotData[, 2] > quantileGroup2S() &
                      plotData[, 3] > quantileGroup2S()) |
          (plotData[, 1] < -quantileGroup2S() &
             plotData[, 2] < -quantileGroup2S() &
             plotData[, 3] < -quantileGroup2S())
      } else if ((sign(corXM2S()) >= 0) & (sign(corMY2S()) < 0)) {
        colored <- (plotData[, 1] > quantileGroup2S() &
                      plotData[, 2] > quantileGroup2S() &
                      plotData[, 3] < -quantileGroup2S()) |
          (plotData[, 1] < -quantileGroup2S() &
             plotData[, 2] < -quantileGroup2S() &
             plotData[, 3] > quantileGroup2S())
      } else if ((sign(corXM2S()) < 0) & (sign(corMY2S()) >= 0)) {
        colored <- (plotData[, 1] > quantileGroup2S() &
                      plotData[, 2] < -quantileGroup2S() &
                      plotData[, 3] < -quantileGroup2S()) |
          (plotData[, 1] < -quantileGroup2S() &
             plotData[, 2] > quantileGroup2S() &
             plotData[, 3] > quantileGroup2S())
      } else if ((sign(corXM2S()) < 0) & (sign(corMY2S()) < 0)) {
        colored <- (plotData[, 1] > quantileGroup2S() &
                      plotData[, 2] < -quantileGroup2S() &
                      plotData[, 3] > quantileGroup2S()) |
          (plotData[, 1] < -quantileGroup2S() &
             plotData[, 2] > quantileGroup2S() &
             plotData[, 3] < -quantileGroup2S())
      }

      pointColors <- rep("black", times = input$n2S)
      pointColors[colored] <- IlliniAltOrg

      plot3d(plotData,
             col = pointColors,
             xlab = "X",
             ylab = "M",
             zlab = "Y",
             xlim = c(-3, 3),
             ylim = c(-3, 3),
             zlim = c(-3, 3)
      )

      if (input$ellipsoids2S) {
        if (!isTwoStepPosDef()) {
          output$ellipsoid2SMessage <- renderText("The correlation matrix is not positive definite. No ellipsoids will be drawn.")
        } else {
          output$ellipsoid2SMessage <- renderText("")

          plot3d(
            ellipse3d(x = twoStepMatrix(),
                      centre = c(0, 0, 0),
                      level = 0.95),
            xlab = "X",
            ylab = "M",
            zlab = "Y",
            col = IlliniBlue,
            alpha = 0.15
          )
        }
      }

      if (input$cubes2S) {
        cubeT <- cube3d(col = IlliniOrange, alpha = 0.2)
        cubeT$vb[cubeT$vb == -1] <- 0
        cubeT$vb[4,  ] <- 1 / (3 - quantileGroup2S())
        cubeT$vb[-4, ] <- cubeT$vb[-4, ] + (quantileGroup2S() / (3 - quantileGroup2S()))

        cubeB <- cube3d(col = IlliniArches, alpha = 0.25)
        cubeB$vb[cubeB$vb ==  1] <- 0
        cubeB$vb[4,  ] <- 1 / (3 - quantileGroup2S())
        cubeB$vb[-4, ] <- cubeB$vb[-4, ] - (quantileGroup2S() / (3 - quantileGroup2S()))

        if ((sign(corXM2S()) >= 0) & (sign(corMY2S()) < 0)) {
          cubeT$vb[3, ] <- -cubeT$vb[3, ]
          cubeB$vb[3, ] <- -cubeB$vb[3, ]
        } else if ((sign(corXM2S()) < 0) & (sign(corMY2S()) >= 0)) {
          cubeT$vb[1, ] <- -cubeT$vb[1, ]
          cubeB$vb[1, ] <- -cubeB$vb[1, ]
        } else if ((sign(corXM2S()) < 0) & (sign(corMY2S()) < 0)) {
          cubeT$vb[2, ] <- -cubeT$vb[2, ]
          cubeB$vb[2, ] <- -cubeB$vb[2, ]
        }

        shade3d(cubeT)
        shade3d(cubeB)

      }


      rglwidget()
    }

  })



  # # Three-Step mediation computations ----

  # Digits for output
  digits3S <- reactive(input$digits3S)
  format3S <- reactive(paste0("%0.", digits3S(), "f"))


  # Reactive reading of the correlations and the percentile for defining the top (bottom) groups
  #
  quantileGroup3S <- reactive(qnorm(input$cutOff3S / 100))
  corXM13S  <- reactive(input$XM13Scor)
  corM1Y3S  <- reactive(input$M1Y3Scor)
  corXM23S  <- reactive(input$XM23Scor)
  corM2Y3S  <- reactive(input$M2Y3Scor)
  corM1M23S <- reactive(input$M1M23Scor)
  corXY3S   <- reactive(input$XY3Scor)


  # Reactive expression to check if the matrix is positive definite
  threeStepMatrix <- reactive({
    matrix(c(1,          corXM13S(),  corXM23S(),  corXY3S(),
             corXM13S(), 1,           corM1M23S(), corM1Y3S(),
             corXM23S(), corM1M23S(), 1,           corM2Y3S(),
             corXY3S(),  corM1Y3S(),  corM2Y3S(),  1),

           nrow = 4)
  })

  threeStepFullMedMatrix <- reactive({
    matrix(c(1,                                     corXM13S(),               corXM13S() * corM1M23S(), corXM13S() * corM1M23S() * corM2Y3S(),
             corXM13S(),                            1,                        corM1M23S(),              corM1M23S() * corM2Y3S(),
             corXM13S() * corM1M23S(),              corM1M23S(),              1,                        corM2Y3S(),
             corXM13S() * corM1M23S() * corM2Y3S(), corM1M23S() * corM2Y3S(), corM2Y3S(),               1),

           nrow = 4)
  })

  threeStepNoMedMatrix <- reactive({
    matrix(c(1,         0,           0,           corXY3S(),
             0,         1,           corM1M23S(), 0,
             0,         corM1M23S(), 1,           0,
             corXY3S(), 0,           0,           1),

           nrow = 4)
  })

  threeStepPcor <- reactive(cor2pcor(threeStepMatrix()))

  pcorXY3S <- reactive(threeStepPcor()[1, 4])


  threeStepSimData <- reactive({
    mvtnorm::rmvnorm(n = input$n3S,
                     mean = c(0, 0, 0, 0),
                     sigma = threeStepMatrix())
  })

  isThreeStepSemiPosDef <- reactive(all(eigen(threeStepMatrix())$values > -10 * .Machine$double.eps))
  isThreeStepPosDef     <- reactive(all(eigen(threeStepMatrix())$values >  10 * .Machine$double.eps))


  ### Proportions Three-Step ------

  propXM13S <- reactive({
    if (quantileGroup3S() == 0) {
      propXM13S <- .25 + (asin(abs(corXM13S()))) / (2 * pi)
    } else {
      propXM13S <- mvtnorm::pmvnorm(lower = c(quantileGroup3S(), quantileGroup3S()),
                                    upper = c(Inf, Inf),
                                    mean = c(0, 0),
                                    corr = abs(threeStepMatrix()[c(1, 2), c(1, 2)]))
    }
    return(propXM13S)
  })


  propM1M23S <- reactive({
    if (quantileGroup3S() == 0) {
      propM1M23S <- .25 + (asin(abs(corM1M23S()))) / (2 * pi)
    } else {
      propM1M23S <- mvtnorm::pmvnorm(lower = c(quantileGroup3S(), quantileGroup3S()),
                                     upper = c(Inf, Inf),
                                     mean = c(0, 0),
                                     corr = abs(threeStepMatrix()[c(2, 3), c(2, 3)]))
    }
    return(propM1M23S)
  })

  qropM1M23S <- reactive({
    if (quantileGroup3S() == 0) {
      propM1M23S <- .25 + (asin(-abs(corM1M23S()))) / (2 * pi)
    } else {
      propM1M23S <- mvtnorm::pmvnorm(lower = c(quantileGroup3S(), quantileGroup3S()),
                                     upper = c(Inf, Inf),
                                     mean = c(0, 0),
                                     corr = matrix(c(1, -corM1M23S(),
                                                     -corM1M23S(), 1), nrow = 2))
    }
    return(propM1M23S)
  })

  propM2Y3S <- reactive({
    if (quantileGroup3S() == 0) {
      propM2Y3S <- .25 + (asin(abs(corM2Y3S()))) / (2 * pi)
    } else {
      propM2Y3S <- mvtnorm::pmvnorm(lower = c(quantileGroup3S(), quantileGroup3S()),
                                    upper = c(Inf, Inf),
                                    mean = c(0, 0),
                                    corr = abs(threeStepMatrix()[c(3, 4), c(3, 4)]))
    }
    return(propM2Y3S)
  })

  corXY3SSign <- reactive({
    ifelse(sign(corXM13S()) * sign(corM1M23S()) * sign(corM2Y3S()) == 0,
           1,
           sign(corXM13S()) * sign(corM1M23S()) * sign(corM2Y3S())
    )
  })

  propXY3S <- reactive({
    if (quantileGroup3S() == 0) {
      propXY3S <- .25 + (asin(corXY3SSign() * (corXY3S()))) / (2 * pi)
    } else {
      propXY3S <- mvtnorm::pmvnorm(lower = c(quantileGroup3S(), quantileGroup3S()),
                                   upper = c(Inf, Inf),
                                   mean = c(0, 0),
                                   corr = matrix(c(1, corXY3SSign() * corXY3S(),
                                                   corXY3SSign() * corXY3S(), 1), nrow = 2))
    }
    return(propXY3S)
  })

  qropXY3S <- reactive({
    if (quantileGroup3S() == 0) {
      propXY3S <- .25 + (asin(-corXY3SSign() * (corXY3S()))) / (2 * pi)
    } else {
      propXY3S <- mvtnorm::pmvnorm(lower = c(quantileGroup3S(), quantileGroup3S()),
                                   upper = c(Inf, Inf),
                                   mean = c(0, 0),
                                   corr = matrix(c(1, -corXY3SSign() * corXY3S(),
                                                   -corXY3SSign() * corXY3S(), 1), nrow = 2))
    }
    return(propXY3S)
  })



  # Message for partial correlation
  output$threeStepPartialMessage <- reactive({
    if (isThreeStepSemiPosDef()) {
      textPartial3S <- paste0("The partial correlation of X and Y is ",
                              strong(sprintf(format3S(), round(pcorXY3S(), digits3S()))))
    } else {
      validate("The correlation matrix is not semi-positive definite. Change some correlation values.")
    }

    return(textPartial3S)
  })



  # Messages for the proportions of 'satisfaction'
  #
  # X-M1 pair
  output$threeStepXM1ProbMessage <- reactive({

    if (corXM13S() >= 0) {
      textXM13S <- paste0("The proportion of observations that are on the top ",
                          100 - input$cutOff3S, "% on both X and M", tags$sub(1)," is ")
    } else {
      textXM13S <- paste0("The proportion of observations that are on the top ",
                          100 - input$cutOff3S, "% on X and the bottom ",
                          100 - input$cutOff3S, "% on M", tags$sub(1)," is ")
    }

    textXM13S <- paste0(textXM13S, strong(sprintf(format3S(), round(propXM13S(), digits3S()))))

    if (!isThreeStepSemiPosDef()) {
      textXM13S <- ""
    }

    return(textXM13S)
  })

  # X-M2 pair
  output$threeStepXM2ProbMessage <- reactive({

    if (corXM23S() >= 0) {
      textXM23S <- paste0("The proportion of observations that are on the top ",
                          100 - input$cutOff3S, "% on both X and M", tags$sub(2)," is ")
    } else {
      textXM23S <- paste0("The proportion of observations that are on the top ",
                          100 - input$cutOff3S, "% on X and the bottom ",
                          100 - input$cutOff3S, "% on M", tags$sub(2)," is ")
    }

    textXM23S <- paste0(textXM23S, strong(sprintf(format3S(), round(propXM23S(), digits3S()))))

    if (!isThreeStepSemiPosDef()) {
      textXM23S <- ""
    }

    return(textXM23S)
  })


  # M1-Y pair
  output$threeStepM1YProbMessage <- reactive({

    if (corM1Y3S() >= 0) {
      textM1Y3S <- paste0("The proportion of observations that are on the top ",
                          100 - input$cutOff3S, "% on both Y and M", tags$sub(1)," is ")
    } else {
      textM1Y3S <- paste0("The proportion of observations that are on the top ",
                          100 - input$cutOff3S, "% on M", tags$sub(1)," and the bottom ",
                          100 - input$cutOff3S, "% on Y is ")
    }

    textM1Y3S <- paste0(textM1Y3S, strong(sprintf(format3S(), round(propM1Y3S(), digits3S()))))

    if (!isThreeStepSemiPosDef()) {
      textM1Y3S <- ""
    }

    return(textM1Y3S)
  })


  # M2-Y pair
  output$threeStepM2YProbMessage <- reactive({

    if (corM2Y3S() >= 0) {
      textM2Y3S <- paste0("The proportion of observations that are on the top ",
                          100 - input$cutOff3S, "% on both Y and M", tags$sub(2)," is ")
    } else {
      textM2Y3S <- paste0("The proportion of observations that are on the top ",
                          100 - input$cutOff3S, "% on M", tags$sub(2)," and the bottom ",
                          100 - input$cutOff3S, "% on Y is ")
    }

    textM2Y3S <- paste0(textM2Y3S, strong(sprintf(format3S(), round(propM2Y3S(), digits3S()))))

    if (!isThreeStepSemiPosDef()) {
      textM2Y3S <- ""
    }

    return(textM2Y3S)
  })


  # M1-M2 pair
  output$threeStepM1M2ProbMessage <- reactive({

    if (corM1M23S() >= 0) {
      textM1M23S <- paste0("The proportion of observations that are on the top ",
                           100 - input$cutOff3S, "% on both M", tags$sub(1)," and M", tags$sub(2)," is ")
    } else {
      textM1M23S <- paste0("The proportion of observations that are on the top ",
                           100 - input$cutOff3S, "% on M", tags$sub(1)," and the bottom ",
                           100 - input$cutOff3S, "% on M", tags$sub(2)," is ")
    }

    textM1M23S <- paste0(textM1M23S, strong(sprintf(format3S(), round(propM1M23S(), digits3S()))))

    if (!isThreeStepSemiPosDef()) {
      textM1M23S <- ""
    }

    return(textM1M23S)
  })


  # X-Y pair
  output$threeStepXYProbMessage <- reactive({

    if (corXY3SSign() >= 0) {
      textXY3S <- paste0("The proportion of observations that are on the top ",
                         100 - input$cutOff3S, "% on both X and Y is ")
    } else {
      textXY3S <- paste0("The proportion of observations that are on the top ",
                         100 - input$cutOff3S, "% on X and the bottom ",
                         100 - input$cutOff3S, "% on Y is ")
    }

    textXY3S <- paste0(textXY3S, strong(sprintf(format3S(), round(propXY3S(), 2))))

    if (!isThreeStepSemiPosDef()) {
      textXY3S <- ""
    }

    return(textXY3S)
  })


  # Proportion in mediation diagonal -----
  #
  threeStepDiagonal <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, Inf, -quantileGroup3S()))

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = threeStepMatrix())[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepBadXM1 <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, Inf, -quantileGroup3S()))

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = threeStepMatrix())[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepBadXM1M1M2 <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, Inf, -quantileGroup3S()))

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = threeStepMatrix())[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepBadM1M2M2Y <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, Inf, -quantileGroup3S()))

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = threeStepMatrix())[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepBadM2Y <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, Inf, -quantileGroup3S()))

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = threeStepMatrix())[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepBadM1M2 <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, Inf, -quantileGroup3S()))

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = threeStepMatrix())[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepBadXM1M2Y <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, Inf, -quantileGroup3S()))

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = threeStepMatrix())[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepBadXM1M1M2M2Y <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, Inf, -quantileGroup3S()))

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = threeStepMatrix())[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })



  # Proportion according to mediation
  #
  threeStepMediation<- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, Inf, -quantileGroup3S()))
    corMat <- threeStepMatrix()
    corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]
    corMat[2, 4] <- corMat[4, 2] <- corMat[2, 3] * corMat[3, 4]
    corMat[1, 4] <- corMat[4, 1] <- corMat[1, 2] * corMat[2, 3] * corMat[3, 4]

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = corMat)[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepMediationBadXM1 <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, Inf, -quantileGroup3S()))
    corMat <- threeStepMatrix()
    corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]
    corMat[2, 4] <- corMat[4, 2] <- corMat[2, 3] * corMat[3, 4]
    corMat[1, 4] <- corMat[4, 1] <- corMat[1, 2] * corMat[2, 3] * corMat[3, 4]

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = corMat)[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepMediationBadXM1M1M2 <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, Inf, -quantileGroup3S()))
    corMat <- threeStepMatrix()
    corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]
    corMat[2, 4] <- corMat[4, 2] <- corMat[2, 3] * corMat[3, 4]
    corMat[1, 4] <- corMat[4, 1] <- corMat[1, 2] * corMat[2, 3] * corMat[3, 4]

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = corMat)[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepMediationBadM1M2M2Y <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, Inf, -quantileGroup3S()))
    corMat <- threeStepMatrix()
    corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]
    corMat[2, 4] <- corMat[4, 2] <- corMat[2, 3] * corMat[3, 4]
    corMat[1, 4] <- corMat[4, 1] <- corMat[1, 2] * corMat[2, 3] * corMat[3, 4]

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = corMat)[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepMediationBadM2Y <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, Inf, -quantileGroup3S()))
    corMat <- threeStepMatrix()
    corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]
    corMat[2, 4] <- corMat[4, 2] <- corMat[2, 3] * corMat[3, 4]
    corMat[1, 4] <- corMat[4, 1] <- corMat[1, 2] * corMat[2, 3] * corMat[3, 4]

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = corMat)[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepMediationBadM1M2 <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, Inf, -quantileGroup3S()))
    corMat <- threeStepMatrix()
    corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]
    corMat[2, 4] <- corMat[4, 2] <- corMat[2, 3] * corMat[3, 4]
    corMat[1, 4] <- corMat[4, 1] <- corMat[1, 2] * corMat[2, 3] * corMat[3, 4]

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = corMat)[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepMediationBadXM1M2Y <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, Inf, -quantileGroup3S()))
    corMat <- threeStepMatrix()
    corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]
    corMat[2, 4] <- corMat[4, 2] <- corMat[2, 3] * corMat[3, 4]
    corMat[1, 4] <- corMat[4, 1] <- corMat[1, 2] * corMat[2, 3] * corMat[3, 4]

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = corMat)[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })

  threeStepMediationBadXM1M1M2M2Y <- reactive({
    signXM1  <- sign(corXM13S())
    signM1M2 <- sign(corM1M23S())
    signM2Y  <- sign(corM2Y3S())

    lower1 <- c(-Inf,
                ifelse(signXM1 < 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2) >= 0, -Inf, quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) < 0, -Inf, quantileGroup3S()))
    upper1 <- c(-quantileGroup3S(),
                ifelse(signXM1 >= 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2) < 0, Inf, -quantileGroup3S()),
                ifelse((signXM1 * signM1M2 * signM2Y) >= 0, Inf, -quantileGroup3S()))
    corMat <- threeStepMatrix()
    corMat[1, 3] <- corMat[3, 1] <- corMat[1, 2] * corMat[2, 3]
    corMat[2, 4] <- corMat[4, 2] <- corMat[2, 3] * corMat[3, 4]
    corMat[1, 4] <- corMat[4, 1] <- corMat[1, 2] * corMat[2, 3] * corMat[3, 4]

    propDiag <- mvtnorm::pmvnorm(lower = lower1, upper = upper1,
                                 mean = c(0, 0, 0, 0),
                                 corr = corMat)[1]

    if (!isThreeStepSemiPosDef()) {
      propDiag <- NA
    }

    return(propDiag)
  })




  # Generate three-step output Table

  output$threeStepTable <- function() ({

    table3S <- data.frame(Violations = c('No violation', 'No violation',
                                         paste0('XM', tags$sub(1)), paste0('XM', tags$sub(1)),
                                         paste0('XM', tags$sub(1), '; M', tags$sub(1), 'M', tags$sub(2)), paste0('XM', tags$sub(1), '; M', tags$sub(1), 'M', tags$sub(2)),
                                         paste0('M', tags$sub(1), 'M', tags$sub(2), '; M', tags$sub(2), 'Y'), paste0('M', tags$sub(1), 'M', tags$sub(2), '; M', tags$sub(2), 'Y'),
                                         paste0('M', tags$sub(2), 'Y'), paste0('M', tags$sub(2), 'Y'),
                                         paste0('M', tags$sub(1), 'M', tags$sub(2)), paste0('M', tags$sub(1), 'M', tags$sub(2)),
                                         paste0('XM', tags$sub(1), '; M', tags$sub(2), 'Y'), paste0('XM', tags$sub(1), '; M', tags$sub(2), 'Y'),
                                         paste0('XM', tags$sub(1), '; M', tags$sub(1), 'M', tags$sub(2), '; M', tags$sub(2), 'Y'), paste0('XM', tags$sub(1), '; M', tags$sub(1), 'M', tags$sub(2), '; M', tags$sub(2), 'Y')
    ),
    X = rep(c(checkmark, otimes), 8),
    M1 = c(checkmark, otimes, otimes, checkmark, otimes, checkmark, checkmark, otimes,
           checkmark, otimes, checkmark, otimes, otimes, checkmark, otimes, checkmark),
    M2 = c(rep(c(checkmark, otimes, otimes, checkmark), 3),
           otimes, checkmark, checkmark, otimes),
    Y = c(checkmark, otimes, otimes, checkmark, checkmark, otimes, checkmark, otimes,
          otimes, checkmark, otimes, checkmark, checkmark, otimes, otimes, checkmark),
    Actual = c(rep(sprintf(format3S(), round(threeStepDiagonal(), digits3S())), 2),
               rep(sprintf(format3S(), round(threeStepBadXM1(), digits3S())), 2),
               rep(sprintf(format3S(), round(threeStepBadXM1M1M2(), digits3S())), 2),
               rep(sprintf(format3S(), round(threeStepBadM1M2M2Y(), digits3S())), 2),
               rep(sprintf(format3S(), round(threeStepBadM2Y(), digits3S())), 2),
               rep(sprintf(format3S(), round(threeStepBadM1M2(), digits3S())), 2),
               rep(sprintf(format3S(), round(threeStepBadXM1M2Y(), digits3S())), 2),
               rep(sprintf(format3S(), round(threeStepBadXM1M1M2M2Y(), digits3S())), 2)),
    Complete = c(rep(sprintf(format3S(), round(threeStepMediation(), digits3S())), 2),
                 rep(sprintf(format3S(), round(threeStepMediationBadXM1(), digits3S())), 2),
                 rep(sprintf(format3S(), round(threeStepMediationBadXM1M1M2(), digits3S())), 2),
                 rep(sprintf(format3S(), round(threeStepMediationBadM1M2M2Y(), digits3S())), 2),
                 rep(sprintf(format3S(), round(threeStepMediationBadM2Y(), digits3S())), 2),
                 rep(sprintf(format3S(), round(threeStepMediationBadM1M2(), digits3S())), 2),
                 rep(sprintf(format3S(), round(threeStepMediationBadXM1M2Y(), digits3S())), 2),
                 rep(sprintf(format3S(), round(threeStepMediationBadXM1M1M2M2Y(), digits3S())), 2)),
    NoM = c(rep(sprintf(format3S(), round(propM1M23S() * propXY3S(), digits3S())), 2),
            rep(sprintf(format3S(), round(propM1M23S() * qropXY3S(), digits3S())), 2),
            rep(sprintf(format3S(), round(qropM1M23S() * propXY3S(), digits3S())), 2),
            rep(sprintf(format3S(), round(qropM1M23S() * propXY3S(), digits3S())), 2),
            rep(sprintf(format3S(), round(propM1M23S() * qropXY3S(), digits3S())), 2),
            rep(sprintf(format3S(), round(qropM1M23S() * qropXY3S(), digits3S())), 2),
            rep(sprintf(format3S(), round(propM1M23S() * propXY3S(), digits3S())), 2),
            rep(sprintf(format3S(), round(qropM1M23S() * qropXY3S(), digits3S())), 2)))

    if (quantileGroup3S() != 0) {

      table3S[nrow(table3S) + 1, -c(1:5)] <- sprintf(format3S(),
                                                     1 - colSums(apply(table3S[, -c(1:5)], 2, as.numeric)))
      table3S[nrow(table3S), 1:5] <- c('Other', rep('-', 4))

    }

    kable(table3S, align = 'c', format = 'html',
          escape = FALSE,
          col.names = c(
            'Pattern Violations',
            'X', paste0('M', tags$sub(1)), paste0('M', tags$sub(2)), 'Y',
            'Actual Proportion',
            'Full Mediation',
            'No Mediation'
          )
    ) %>%
      kable_styling(bootstrap_options = c("striped", "responsive"),
                    full_width = TRUE,
                    position = "left", font_size = 14)%>%
      column_spec(1, bold = TRUE) %>%
      collapse_rows(columns = 1, valign = "middle") ->
      table3S

    if (!isThreeStepSemiPosDef()) {
      table3S <- ""
    }

    return(table3S)
  }
  )





  # Generate three-step plots: X -> M1 -> Y and X -> M2 -> Y----

  output$threeStepMatrixCaption <- reactive({
    if (!isThreeStepSemiPosDef()) {
      outText <- ""
    } else {
      outText <- as.character(h4(strong("Correlation matrices")))
    }
    return(outText)
  })

  output$threeStepCorMatrixCaption <- reactive({
    if (!isThreeStepSemiPosDef()) {
      outText <- ""
    } else {
      outText <- as.character(strong("Actual", br(), "correlations"))
    }
    return(outText)
  })

  output$threeStepFullMatrixCaption <- reactive({
    if (!isThreeStepSemiPosDef()) {
      outText <- ""
    } else {
      outText <- as.character(strong("Implied by", br(), "full mediation"))
    }
    return(outText)
  })

  output$threeStepNoMatrixCaption <- reactive({
    if (!isThreeStepSemiPosDef()) {
      outText <- ""
    } else {
      outText <- as.character(strong("Implied by", br(), "no mediation"))
    }
    return(outText)
  })

  output$threeStepCorMatrix <- renderPlot({
    if (!isThreeStepSemiPosDef()) {
    } else {
      mat <- threeStepMatrix()
      rownames(mat) <- colnames(mat) <- c("X", "M_1", "M_2", "Y")
      corrplot::corrplot(corr = mat,
                         type = "lower", method = "color",
                         diag = FALSE,
                         number.cex = 1.5,
                         addCoef.col = IlliniDkBlue,
                         addgrid.col = IlliniDkBlue,
                         cl.pos = FALSE,
                         tl.cex = 1.4,
                         tl.offset = 1,
                         tl.srt = 0,
                         tl.col = IlliniDkBlue,
                         col = corMatricesBgPalette)
    }
  })

  output$threeStepFullMedMatrix <- renderPlot({
    if (!isThreeStepSemiPosDef()) {
    } else {
      mat <- threeStepFullMedMatrix()
      rownames(mat) <- colnames(mat) <- c("X", "M_1", "M_2", "Y")
      corrplot::corrplot(corr = mat,
                         type = "lower", method = "color",
                         diag = FALSE,
                         number.cex = 1.5,
                         addCoef.col = IlliniDkBlue,
                         addgrid.col = IlliniDkBlue,
                         cl.pos = FALSE,
                         tl.cex = 1.4,
                         tl.offset = 1,
                         tl.srt = 0,
                         tl.col = IlliniDkBlue,
                         col = corMatricesBgPalette)
    }
  })

  output$threeStepNoMedMatrix <- renderPlot({
    if (!isThreeStepSemiPosDef()) {
    } else {
      mat <- threeStepNoMedMatrix()
      rownames(mat) <- colnames(mat) <- c("X", "M_1", "M_2", "Y")
      corrplot::corrplot(corr = mat,
                         type = "lower", method = "color",
                         diag = FALSE,
                         number.cex = 1.5,
                         addCoef.col = IlliniDkBlue,
                         addgrid.col = IlliniDkBlue,
                         cl.pos = FALSE,
                         tl.cex = 1.4,
                         tl.offset = 1,
                         tl.srt = 0,
                         tl.col = IlliniDkBlue,
                         col = corMatricesBgPalette)
    }
  })


  output$threeStepM1Plot <- renderRglwidget({

    if (!isThreeStepSemiPosDef()) {
    } else {
      close3d()
      plotData <- threeStepSimData()

      if ((sign(corXM13S()) >= 0) & (sign(corM1Y3S()) >= 0)) {
        colored <- (plotData[, 1] > quantileGroup3S() &
                      plotData[, 2] > quantileGroup3S() &
                      plotData[, 4] > quantileGroup3S()) |
          (plotData[, 1] < -quantileGroup3S() &
             plotData[, 2] < -quantileGroup3S() &
             plotData[, 4] < -quantileGroup3S())
      } else if ((sign(corXM13S()) >= 0) & (sign(corM1Y3S()) < 0)) {
        colored <- (plotData[, 1] > quantileGroup3S() &
                      plotData[, 2] > quantileGroup3S() &
                      plotData[, 4] < -quantileGroup3S()) |
          (plotData[, 1] < -quantileGroup3S() &
             plotData[, 2] < -quantileGroup3S() &
             plotData[, 4] > quantileGroup3S())
      } else if ((sign(corXM13S()) < 0) & (sign(corM1Y3S()) >= 0)) {
        colored <- (plotData[, 1] > quantileGroup3S() &
                      plotData[, 2] < -quantileGroup3S() &
                      plotData[, 4] < -quantileGroup3S()) |
          (plotData[, 1] < -quantileGroup3S() &
             plotData[, 2] > quantileGroup3S() &
             plotData[, 4] > quantileGroup3S())
      } else if ((sign(corXM13S()) < 0) & (sign(corM1Y3S()) < 0)) {
        colored <- (plotData[, 1] > quantileGroup3S() &
                      plotData[, 2] < -quantileGroup3S() &
                      plotData[, 4] > quantileGroup3S()) |
          (plotData[, 1] < -quantileGroup3S() &
             plotData[, 2] > quantileGroup3S() &
             plotData[, 4] < -quantileGroup3S())
      }

      pointColors <- rep("black", times = input$n3S)
      pointColors[colored] <- IlliniAltOrg

      plot3d(plotData[, c(1, 2, 4)],
             main = expression(X %->% M[1] %->% Y),
             col = pointColors,
             xlab = "X",
             ylab = "M1",
             zlab = "Y",
             xlim = c(-3, 3),
             ylim = c(-3, 3),
             zlim = c(-3, 3)
      )

      if (input$ellipsoids3S) {
        if (!isThreeStepPosDef()) {
          output$ellipsoid3SMessage <- renderText("The correlation matrix is not positive definite. No ellipsoids will be drawn.")
        } else {
          output$ellipsoid3SMessage <- renderText("")

          plot3d(
            ellipse3d(x = threeStepMatrix()[c(1, 2, 4), c(1, 2, 4)],
                      centre = c(0, 0, 0),
                      level = 0.95),
            xlab = "X",
            ylab = "M1",
            zlab = "Y",
            col = IlliniBlue,
            alpha = 0.15
          )
        }
      }

      if (input$cubes3S) {
        cubeT <- cube3d(col = IlliniOrange, alpha = 0.2)
        cubeT$vb[cubeT$vb == -1] <- 0
        cubeT$vb[4,  ] <- 1 / (3 - quantileGroup3S())
        cubeT$vb[-4, ] <- cubeT$vb[-4, ] + (quantileGroup3S() / (3 - quantileGroup3S()))

        cubeB <- cube3d(col = IlliniArches, alpha = 0.25)
        cubeB$vb[cubeB$vb ==  1] <- 0
        cubeB$vb[4,  ] <- 1 / (3 - quantileGroup3S())
        cubeB$vb[-4, ] <- cubeB$vb[-4, ] - (quantileGroup3S() / (3 - quantileGroup3S()))

        if ((sign(corXM13S()) >= 0) & (sign(corM1Y3S()) < 0)) {
          cubeT$vb[3, ] <- -cubeT$vb[3, ]
          cubeB$vb[3, ] <- -cubeB$vb[3, ]
        } else if ((sign(corXM13S()) < 0) & (sign(corM1Y3S()) >= 0)) {
          cubeT$vb[1, ] <- -cubeT$vb[1, ]
          cubeB$vb[1, ] <- -cubeB$vb[1, ]
        } else if ((sign(corXM13S()) < 0) & (sign(corM1Y3S()) < 0)) {
          cubeT$vb[2, ] <- -cubeT$vb[2, ]
          cubeB$vb[2, ] <- -cubeB$vb[2, ]
        }

        shade3d(cubeT)
        shade3d(cubeB)

      }

      rglwidget()
    }

  })

  output$threeStepM2Plot <- renderRglwidget({

    if (!isThreeStepSemiPosDef()) {
    } else {
      close3d()
      plotData <- threeStepSimData()

      if ((sign(corXM23S()) >= 0) & (sign(corM2Y3S()) >= 0)) {
        colored <- (plotData[, 1] > quantileGroup3S() &
                      plotData[, 3] > quantileGroup3S() &
                      plotData[, 4] > quantileGroup3S()) |
          (plotData[, 1] < -quantileGroup3S() &
             plotData[, 3] < -quantileGroup3S() &
             plotData[, 4] < -quantileGroup3S())
      } else if ((sign(corXM23S()) >= 0) & (sign(corM2Y3S()) < 0)) {
        colored <- (plotData[, 1] > quantileGroup3S() &
                      plotData[, 3] > quantileGroup3S() &
                      plotData[, 4] < -quantileGroup3S()) |
          (plotData[, 1] < -quantileGroup3S() &
             plotData[, 3] < -quantileGroup3S() &
             plotData[, 4] > quantileGroup3S())
      } else if ((sign(corXM23S()) < 0) & (sign(corM2Y3S()) >= 0)) {
        colored <- (plotData[, 1] > quantileGroup3S() &
                      plotData[, 3] < -quantileGroup3S() &
                      plotData[, 4] < -quantileGroup3S()) |
          (plotData[, 1] < -quantileGroup3S() &
             plotData[, 3] > quantileGroup3S() &
             plotData[, 4] > quantileGroup3S())
      } else if ((sign(corXM23S()) < 0) & (sign(corM2Y3S()) < 0)) {
        colored <- (plotData[, 1] > quantileGroup3S() &
                      plotData[, 3] < -quantileGroup3S() &
                      plotData[, 4] > quantileGroup3S()) |
          (plotData[, 1] < -quantileGroup3S() &
             plotData[, 3] > quantileGroup3S() &
             plotData[, 4] < -quantileGroup3S())
      }

      pointColors <- rep("black", times = input$n3S)
      pointColors[colored] <- IlliniAltOrg

      plot3d(plotData[, c(1, 3, 4)],
             main = expression(X %->% M[2] %->% Y),
             col = pointColors,
             xlab = "X",
             ylab = "M2",
             zlab = "Y",
             xlim = c(-3, 3),
             ylim = c(-3, 3),
             zlim = c(-3, 3)
      )

      if (input$ellipsoids3S) {

        plot3d(
          ellipse3d(x = threeStepMatrix()[c(1, 3, 4), c(1, 3, 4)],
                    centre = c(0, 0, 0),
                    level = 0.95),
          xlab = "X",
          ylab = "M2",
          zlab = "Y",
          col = IlliniBlue,
          alpha = 0.15
        )
      }

      if (input$cubes3S) {
        cubeT <- cube3d(col = IlliniOrange, alpha = 0.2)
        cubeT$vb[cubeT$vb == -1] <- 0
        cubeT$vb[4,  ] <- 1 / (3 - quantileGroup3S())
        cubeT$vb[-4, ] <- cubeT$vb[-4, ] + (quantileGroup3S() / (3 - quantileGroup3S()))

        cubeB <- cube3d(col = IlliniArches, alpha = 0.25)
        cubeB$vb[cubeB$vb ==  1] <- 0
        cubeB$vb[4,  ] <- 1 / (3 - quantileGroup3S())
        cubeB$vb[-4, ] <- cubeB$vb[-4, ] - (quantileGroup3S() / (3 - quantileGroup3S()))

        if ((sign(corXM23S()) >= 0) & (sign(corM2Y3S()) < 0)) {
          cubeT$vb[3, ] <- -cubeT$vb[3, ]
          cubeB$vb[3, ] <- -cubeB$vb[3, ]
        } else if ((sign(corXM23S()) < 0) & (sign(corM2Y3S()) >= 0)) {
          cubeT$vb[1, ] <- -cubeT$vb[1, ]
          cubeB$vb[1, ] <- -cubeB$vb[1, ]
        } else if ((sign(corXM23S()) < 0) & (sign(corM2Y3S()) < 0)) {
          cubeT$vb[2, ] <- -cubeT$vb[2, ]
          cubeB$vb[2, ] <- -cubeB$vb[2, ]
        }

        shade3d(cubeT)
        shade3d(cubeB)

      }

      rglwidget()
    }

  })





}