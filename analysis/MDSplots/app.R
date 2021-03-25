#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(ggplot2)
library(tidyverse)
library(shiny)

load("mdsplots.RData")


ui = fluidPage(
    titlePanel("Interactive Samples PCA"),
    fluidRow(
        selectInput("col1", "Color", list("CellLine","gRNA", "KO", "DIFF", "RAPA")),
        selectInput("shape1", "Shape", list("CellLine","gRNA", "KO", "DIFF", "RAPA"))),
    fluidRow(
        plotOutput("PCAplot")
    ),
    titlePanel("Interactive Samples MDS"),
    fluidRow(
        selectInput("col2", "Color", list("CellLine","gRNA", "KO", "DIFF", "RAPA")),
        selectInput("shape2", "Shape", list("CellLine","gRNA", "KO", "DIFF", "RAPA"))),
    fluidRow(
        plotOutput("MDSplot")
    )
)

server = function(input, output) {
    GPCA <- reactive(gpca.dat %>% mutate(col=!!sym(input$col1),
                                         shp=!!sym(input$shape1)))
    MDS <- reactive(mds%>% mutate(col=!!sym(input$col2),
                                  shp=!!sym(input$shape2)))

    output$PCAplot = renderPlot({
        ggplot(GPCA(), aes(x = dim1, y = dim2,
                             color = col, shape = shp)) +
            geom_point(size = 2)  + ggtitle("glmpca - Generalized PCA")})

    output$MDSplot = renderPlot({
        ggplot(MDS(), aes(x = `1`, y = `2`, color = col, shape = shp)) +
            geom_point(size = 2)  + ggtitle("MDS with log2 counts")})
}

# Run the application
shinyApp(ui = ui, server = server)


