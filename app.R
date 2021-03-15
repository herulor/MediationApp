#
# This is a Shiny web application.
#
# Find out more about building applications with Shiny here:
#    http://shiny.rstudio.com/
#
#
# What does mediation reveal about individual people?
# Computing Pattern Probabilities
#
# Víctor H Cervantes
# 2021
# https://github.com/herulor/MediationApp
# https://herulor.shinyapps.io/MediationApp



source("R/ui.R")
source("R/server.R")

shinyApp(ui = ui, server = server)
