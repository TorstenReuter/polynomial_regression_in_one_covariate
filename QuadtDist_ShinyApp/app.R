#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#####################################################################
#### Shiny App to visualize the D-optimal subsampling design for #### 
#### quadratic regression given a t-distribution                 ####
#####################################################################


library(shiny)
library(ggplot2)

dfmain <- readRDS("QuadtDistVis_main_data")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel(withMathJax("Optimal Subsampling Design for the \\(t\\)-Distribution for Quadratic Regression")),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("dof", "Degrees of freedom:",
                       c("5" = 5,
                         "6" = 6,
                         "7" = 7,
                         "8" = 8,
                         "9" = 9,
                         "30"= 30),
                       selected = 8),
            sliderInput("ga",
                        withMathJax("Percentage \\( \\alpha \\) to be sampled"),
                        min = 0.01,
                        max = 0.99,
                        step = 0.01,
                        value = 0.5,
                        animate = animationOptions(interval = 1200))
        ),
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        df <- dfmain[dfmain$dof == input$dof,]
        n <- 800
        eps <- 10^(-3)
        points <- c(seq(-4.5,-0.1,length.out = n/4), seq(-0.1, 0.1, length.out = n/2), seq(0.1,4.5,length.out = n/4))
        CDFphi <- dt(points, df = as.numeric(input$dof))
        idx <- abs(points) >= df[1000*input$ga,]$a | abs(points) <= df[1000*input$ga,]$b
        CDFg <- points
        CDFg[idx] <- CDFphi[idx] 
        CDFg[!idx] <- 0

        dfphi <- data.frame("x" = points, "CDF" = CDFphi)
        dfg <- data.frame("x" = points, "CDF" = CDFg)
        
        ggplot() +
            geom_area(data = dfg, mapping = aes(x=x, y = CDF, fill = "lightcoral")) +
            geom_line(data = dfphi, aes(x = x, y = CDF), size = 1.4, linetype = "dashed", color = "royalblue4") +
            geom_line(data = dfg, aes(x = x, y = CDF), size = 0.6, linetype = "solid", color = "firebrick3") +
            scale_x_continuous(limits = c(-4.5, 4.5)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  legend.position = "none") +
            xlab("x") + 
            ylab("Density")
    })
    

}

# Run the application 
shinyApp(ui = ui, server = server)
