library(shinythemes)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  #shinythemes::themeSelector(), # <--- Add this somewhere in the UI
  titlePanel("Simple Random Sampling Simulation"),
  navlistPanel(
    widths = c(2, 10),
    tabPanel(
      "Finite Population",
      sidebarPanel(
        fluidRow(column(7,
                        h4(
                          "Population Option"
                        )),
                 column(
                   5,
                   actionButton("action", "Generate", class = "btn-primary") #,style = "text-align:right; padding:15px"
                 )),
        numericInput("psize", "Population size N", value = 1000, min = 1),
        selectInput(
          "pdistr",
          "Population distribution",
          choices = list(
            "Uniform" = 1,
            "Bernoulli" = 2,
            "Normal" = 3
          )
        ),
        uiOutput("para")
      ),
      mainPanel(
        h4("Population, Y"),
        fluidRow(column(6,
                        tableOutput("summary1"))),
        h4("Historgram of finite population"),
        plotOutput("popplot")
      )
    ),
    tabPanel(
      "Sampling",
      sidebarPanel(
        fluidRow(column(7,
                        h4("Sampling Option")),
                 column(5,
                        actionButton(
                          "action2", "Generate"
                        ))),
        numericInput("ssize", "Sample size n", value = 100, min = 1),
        hr(),
        h4("Computation Option"),
        numericInput(
          "simsize",
          "Simulation size M",
          value = 10000,
          min = 1
        ),
        hr(),
        h4("Histogram Output Option"),
        selectInput(
          "samplot2",
          NULL,
          choices = list(
            "y bar" = 1,
            "sqrt(n)(y bar - Y bar)" = 2,
            "(y bar - Y bar) / V" = 3
          )
        ),
        uiOutput("makepdf")
      ),
      mainPanel(
        h4("Sample mean, y bar"),
        fluidRow(column(6,
                        tableOutput("summary2"))
                 ),
        h4("Histogram of sample"),
        plotOutput("sampplot")
      )
    ),
    tabPanel("Coverage", 
             sidebarPanel(
               numericInput("clevel", "Confidence Level(%)", value = 95, min = 1, max = 100)
               ),
             mainPanel(
               tableOutput("coverage")
             )
             )
  )
  
)

server <- function(input, output) {
  output$para <- renderUI({
    if (input$pdistr == 1) {
      tagList(strong("Unif(a,b)"),
              fluidRow(column(6,
                              numericInput(
                                "a", "a", value = 0
                              )),
                       column(6,
                              numericInput(
                                "b", "b", value = 1
                              ))))
      
    } else if (input$pdistr == 2) {
      tagList(strong("Bernoulli(p)"),
              numericInput("p", "p", value = 0.5))
    } else if (input$pdistr == 3) {
      tagList(strong("Normal(mu,sigma^2)"),
              fluidRow(column(
                6,
                numericInput("mu", "mu", value = 0)
              ),
              column(
                6,
                numericInput("sigma", "sigma", value = 1)
              )))
    }
  })
  reacpop <- eventReactive(input$action, {
    if (input$pdistr == 1) {
      pop <- runif(input$psize, min = input$a, max = input$b)
      var2 <- 1 / 12 * (input$b - input$a) ^ 2
    } else if (input$pdistr == 2) {
      pop <- rbinom(input$psize, 1, input$p)
      var2 <- input$p * (1 - input$p)
    } else if (input$pdistr == 3) {
      pop <- rnorm(input$psize, mean = input$mu, sd = input$sigma)
      var2 <- input$sigma ^ 2
    }
    return(list(
      pop = pop,
      meanpop = mean(pop),
      var = var2
    ))
  })
  reacsamp <- eventReactive(input$action2, {
    samp <- c()
    sampvar <- c()
    for (i in 1:input$simsize) {
      x <- sample(reacpop()$pop, input$ssize, replace = FALSE)
      samp <- c(samp, mean(x))
      sampvar <- c(sampvar, var(x))
    }
    return(list(samp = samp, sampvar = sampvar))
  })
  output$summary1 <- renderTable({
    table2 <- c(mean(reacpop()$pop), var(reacpop()$pop))
    names(table2) <- c("mean", "variance")
    return(table2)
  }, digits = 7, rownames = TRUE, colnames = FALSE)
  output$summary2 <- renderTable({
    table2 <- c(mean(reacsamp()$samp), var(reacsamp()$samp))
    names(table2) <- c("mean", "variance")
    return(table2)
  }, digits = 7, rownames = TRUE, colnames = FALSE)
  output$makepdf <- renderUI({
    if (input$samplot2 == 2) {
      checkboxInput("makepdf", "pdf of N(0, sigma^2)")
    } else if (input$samplot2 == 3) {
      checkboxInput("makepdf", "pdf of N(0, 1)")
    }
  })
  output$popplot <- renderPlot({
    hist(
      reacpop()$pop,
      main = expression(paste("Historgram of Y")),
      xlab = "population",
      freq = FALSE
    )
  })
  output$sampplot <- renderPlot({
    res2 <- NULL
    expr2 <- NULL
    if (input$samplot2 == 1) {
      res2 <- reacsamp()$samp
      expr2 <- expression(paste("Historgram of ", bar("y")))
    } else if (input$samplot2 == 2) {
      res2 <- (reacsamp()$samp - reacpop()$meanpop) * sqrt(input$ssize)
      expr2 <-
        expression(paste("Historgram of ", sqrt("n")(bar("y") - bar("Y"))))
    } else if (input$samplot2 == 3) {
      V <- (1 / input$ssize - 1 / input$psize) * reacsamp()$sampvar
      res2 <- (reacsamp()$samp - reacpop()$meanpop) / sqrt(V)
      expr2 <-
        expression(paste(
          "Historgram of ",
          (bar("y") - bar("Y")) / sqrt("V") * ", V" == (1 / "n"-1 / "N") * "s" ^
            2
        ))
    }
    hist(res2,
         main = expr2,
         xlab = "sample",
         freq = FALSE)
    if (!is.null(input$makepdf) & all(input$samplot2 == 2, input$makepdf)) {
      min1 <- min(res2)
      max1 <- max(res2)
      vec <- seq(min1, max1, length = 1000)
      lines(vec, dnorm(vec, mean = 0, sd = sqrt(reacpop()$var)), col = "red")
    } else if (!is.null(input$makepdf) & all(input$samplot2 == 3, input$makepdf)) {
      min1 <- min(res2)
      max1 <- max(res2)
      vec <- seq(min1, max1, length = 1000)
      lines(vec, dnorm(vec, mean = 0, sd = 1), col = "red")
    }
  })
  
  output$coverage <- renderTable({
    V <- (1 / input$ssize - 1 / input$psize) * reacsamp()$sampvar
    z <- qnorm(1 - (1 - input$clevel/100) / 2)
    sum2 <-
      sum ((reacpop()$meanpop > reacsamp()$samp - z * sqrt(V)) &
             (reacpop()$meanpop < reacsamp()$samp + z * sqrt(V))
      )
    res3 <- t(sum2 / input$simsize)
    row.names(res3) <- sprintf("Coverage of %d%% confidence interval", input$clevel)
    return(res3)
  }, rownames = TRUE, colnames = FALSE, digits = 5)
}

shinyApp(ui = ui, server = server)