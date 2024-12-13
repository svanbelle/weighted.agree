#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(sjPlot)
library(corrplot)
library(dplyr)
library(vtable)
library(vcd)
library(mvtnorm)
source("./script/functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

  navbarPage("User Interface:",
             tabPanel("Upload",titlePanel("Uploading Files"),
                      sidebarLayout(
                        sidebarPanel(
                          fileInput("file1", "Choose CSV File"),
                          
                          # Input: Checkbox if file has header ----
                          checkboxInput("header", "Header", TRUE),
                          
                          # Input: Select separator ----
                          radioButtons("sep", "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),
                                       selected = ","),
                          
                          strong("Note"),
                          helpText("Categories should be coded from 1 to K")),
                        
                        mainPanel(
                          strong("Dataset overview"),
                          tableOutput("contents"),
                          br(),
                          strong("Marginal probability distribution for the observers (counts)"),
                          verbatimTextOutput("summary"),
                          strong("Marginal probability distribution for the observers (proportions)"),
                          verbatimTextOutput("summary2")
                        ))), 
             
             tabPanel("Several observers",
                      titlePanel("Agreement between several observers"),
                      sidebarLayout(
                        sidebarPanel(
                          checkboxGroupInput("variables", label = ""), 
                          hr(),
                          selectInput("stat_many", strong("Select the agreement index"),
                                       c("Proportion of agreement" = 1, "Proportion of disagreement" = 2,
                                                      "Mean absolute error" = 3,"Mean squared error" = 4,"Kappa coefficient (linear weights, different marginals)" = 5,
                                                      "Kappa coefficient (linear weights, same marginals)" = 6,"Kappa coefficient (quadratic weights, different marginals)" = 7,
                                                      "Kappa coefficient (quadratic weights, same marginals)" = 8))),
                          
                          mainPanel(
                          conditionalPanel(
                            condition = "input.stat_many == 1",
                            strong("Proportion of agreement"),
                            htmlOutput("pom"),
                            br(),
                            strong("Proportion of agreement between pairs of observers"),
                            plotOutput("plotpair1",width = "75%")),
                          conditionalPanel(
                            condition = "input.stat_many == 2",
                            strong("Proportion of disagreement"),
                            htmlOutput("qom"),
                            br(),
                            strong("Proportion of disagreement between pairs of observers"),
                            plotOutput("plotpair2",width = "75%")),
                          conditionalPanel(
                            condition = "input.stat_many == 3",
                            strong("Mean absolute error (MAE)"),
                            htmlOutput("MAE"),
                            br(),
                            strong("Mean absolute error between pairs of observers"),
                            plotOutput("plotpair3",width = "75%")),
                          conditionalPanel(
                            condition = "input.stat_many == 4",
                            strong("Mean square error (MSE)"),
                            htmlOutput("MSE"),
                            br(),
                            strong("Mean square error between pairs of observers"),
                            plotOutput("plotpair4",width = "75%")),
                          conditionalPanel(
                            condition = "input.stat_many == 5",
                            strong("Linear weighted kappa (different marginals)"),
                            htmlOutput("LWK2"),
                            br(),
                            strong("Linear weighted kappa (different marginals) between pairs of observers"),
                            plotOutput("plotpair5",width = "75%")),
                          conditionalPanel(
                            condition = "input.stat_many == 6",
                            strong("Linear weighted kappa (same marginals)"),
                            htmlOutput("LWK3"),
                            br(),
                            strong("Linear weighted kappa (same marginals) between pairs of observers"),
                            plotOutput("plotpair6",width = "75%")),
                          conditionalPanel(
                            condition = "input.stat_many == 7",
                            strong("Quadratic weighted kappa (different marginals)"),
                            htmlOutput("QWK2"),
                            br(),
                            strong("Quadratic weighted kappa (different marginals) between pairs of observers"),
                            plotOutput("plotpair7",width = "75%")),
                          conditionalPanel(
                            condition = "input.stat_many == 8",
                            strong("Quadratic weighted kappa (same marginals)"),
                            htmlOutput("QWK3"),
                            br(),
                            strong("Quadratic weighted kappa (same marginals) between pairs of observers"),
                            plotOutput("plotpair8",width = "75%"))
                        ))),
             tabPanel("Two observers",
                      titlePanel("Agreement between two observers"),
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("row", "Select rater 1 (row)", character(0)),
                          selectInput("col", "Select rater 2 (column)", character(0)),
                          hr(),
                          selectInput("stat_two", strong("Select the agreement index"),
                                      c("Proportion of agreement" = 1, "Proportion of disagreement" = 2,
                                        "Mean absolute error" = 3,"Mean squared error" = 4,"Kappa coefficient (linear weights, different marginals)" = 5,
                                        "Kappa coefficient (linear weights, same marginals)" = 6,"Kappa coefficient (quadratic weights, different marginals)" = 7,
                                        "Kappa coefficient (quadratic weights, same marginals)" = 8))),
                        mainPanel(
                           strong("Contingency table (with cell %)"),
                           htmlOutput("contingency"),
                           br(),
                           strong("Agreement plot"),
                           plotOutput("plotagree",width = "75%"),
                           br(),
                           conditionalPanel(
                             condition = "input.stat_two == 1", 
                             strong("Proportion of agreement"),
                             htmlOutput("pom2")),
                           conditionalPanel(
                             condition = "input.stat_two == 2", 
                             strong("Proportion of disagreement"),
                             htmlOutput("qom2")),
                           conditionalPanel(
                             condition = "input.stat_two == 3", 
                             strong("Mean absolute error"),
                             htmlOutput("MAE2")),
                           conditionalPanel(
                             condition = "input.stat_two == 4", 
                             strong("Mean square error"),
                             htmlOutput("MSE2")),
                           conditionalPanel(
                             condition = "input.stat_two == 5", 
                             strong("Linear weighted kappa (different marginals)"),
                             htmlOutput("LWK22")),
                           conditionalPanel(
                             condition = "input.stat_two == 6", 
                             strong("Linear weighted kappa (same marginals)"),
                             htmlOutput("LWK32")),
                           conditionalPanel(
                             condition = "input.stat_two == 7", 
                             strong("Quadratic weighted kappa (different marginals)"),
                             htmlOutput("QWK22")),
                           conditionalPanel(
                             condition = "input.stat_two == 8", 
                             strong("Quadratic weighted kappa (same marginals)"),
                             htmlOutput("QWK32"))
                        ))),
             tabPanel("Sample size (confidence interval approach)",
                      titlePanel("What can I get with my sample sizes"),
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("ncat", label = "Number of categories:",
                                      choices = 2:10, selected = 2),
                          numericInput("category_1", "% of participants expected in category (1):",0.50,min=0, max=1,step=0.01),
                          uiOutput("new_categories"),
                          selectInput("nrater", label = "Number of raters:",
                                      choices = 2:30, selected = 2),
                          sliderInput("nsub", label = "Number of particpants:",
                                      min= 5, max= 1000, value = 30),
                          sliderInput("alpha", label = "level of uncertnainty:",
                                      min= 0, max= 1, value = 0.05),
                          sliderInput("nsim", label = "Number of simulated datasets:",
                                      min = 100, max = 10000, value = 500),
                          # sliderInput("ncover", label = "Expected width of the confidence interval:",
                          #             min = 0, max= 100, value = 100),
                          numericInput("seed", label = "Seed number:",
                                      min = 0,step=1, value = 1234),
                          selectInput("size_agree", strong("Select the agreement index"),
                                      c("Proportion of agreement" = 1, "Proportion of disagreement" = 2,
                                        "Mean absolute error" = 3,"Mean squared error" = 4,"Kappa coefficient (linear weights)" = 5,
                                        "Kappa coefficient (quadratic weights)" = 6)),
                          conditionalPanel("input.size_agree==1", uiOutput("CI_PA")),
                          conditionalPanel("input.size_agree==2", uiOutput("CI_PD")),
                          conditionalPanel("input.size_agree==3", uiOutput("CI_MAE")),
                          conditionalPanel("input.size_agree==4", uiOutput("CI_MSE")),
                          conditionalPanel("input.size_agree==5", uiOutput("CI_LINEAR")),
                          conditionalPanel("input.size_agree==6", uiOutput("CI_QUADRATIC")),
                          actionButton(inputId = "size_GO", label = "Calculate")
                          ),
                        mainPanel(
                          htmlOutput("size_CIagree"),
                          plotOutput("size_CIagree2")
                        ))),
             tabPanel("Sample size (testing approach)",
                      titlePanel("What can I get with my sample sizes"),
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("tncat", label = "Number of categories:",
                                      choices = 2:10, selected = 2),
                          numericInput("tcategory_1", "% of participants expected in category (1):",0.50,min=0, max=1,step=0.01),
                          uiOutput("tnew_categories"),
                          selectInput("tnrater", label = "Number of raters:",
                                      choices = 2:30, selected = 2),
                          sliderInput("tnsub", label = "Number of particpants:",
                                      min= 5, max= 1000, value = 30),
                          sliderInput("talpha", label = "level of uncertnainty:",
                                      min= 0, max= 1, value = 0.05),
                          sliderInput("tnsim", label = "Number of simulated datasets:",
                                      min = 100, max = 10000, value = 500),
                          sliderInput("tncover", label = "Percentile for power (%):",
                                      min = 0, max= 100, value = 100),
                          numericInput("tseed", label = "Seed number:",
                                       min = 0,step=1, value = 1234),
                          selectInput("tsize_agree", strong("Select the agreement index"),
                                      c("Proportion of agreement" = 1, "Proportion of disagreement" = 2,
                                        "Mean absolute error" = 3,"Mean squared error" = 4,"Kappa coefficient (linear weights)" = 5,
                                        "Kappa coefficient (quadratic weights)" = 6)),
                          conditionalPanel("input.tsize_agree==1", uiOutput("TEST_PA")),
                          conditionalPanel("input.tsize_agree==2", uiOutput("TEST_PD")),
                          conditionalPanel("input.tsize_agree==3", uiOutput("TEST_MAE")),
                          conditionalPanel("input.tsize_agree==4", uiOutput("TEST_MSE")),
                          conditionalPanel("input.tsize_agree==5", uiOutput("TEST_LINEAR")),
                          conditionalPanel("input.tsize_agree==6", uiOutput("TEST_QUADRATIC")),
                          actionButton(inputId = "tsize_GO", label = "Calculate")
                        ),
                        mainPanel(
                          htmlOutput("size_tagree"),
                          plotOutput("size_tagree2")
                        ))),
             tabPanel("Explanations",
                      titlePanel("Explanations"),
                      mainPanel(
                        h4("Agreement plot"),
                        p("The agreement plot is a visual summary of the cross-classification table.
                        The graph is divided in 4 rectangles according to the marginal probability of the two observers. 
                        For example, if observer 1 rates 60% of the participants as positive, an horizontal line will be drawn at 60.
                        For observer 2, a vertical line is drawn. The smallest side of the upper left rectangle represents 
                        the maximum proportion of participants on which the two obervers can agree for the category 0.
                        Similarly, the lower right rectangle is relative to category 1."),
                        p("The black rectangles represent the current proportion of agreement on category 0 (upper left) and category 1 (lower right).
                         The maximum possible agreement is obtained when one side of a black rectangle coincides with one side of the white rectangles."),
                        p("When the white rectangles cross on the red line, this means that the two observers rated the same proportion of 
                        participants in category 1 (i.e., they have the same marginal probability distribution). If the rectangles
                        cross below the red line, this means that a smaller proportion of participants was rated positive by the first observer
                        than by second observer. The inverse is true when the 4 rectangles cross above the red lines."),
                        br(),
                        h4("Correlogram"),
                        p("The correlogram gives the proportion of agreement (po), of positive agreement (pa) or Chamberlain positive agreement (pca)
                         between all possible pairs of observers. The size of a circle is proportional to the value of the agreement, given
                         in white in the circles. Different colors are also used depending on the agreement level."),
                        p("This plot permits to determine if one particular observer contributes very poorly to the overall agreement index
                         or if the agreement is homogeneous between all pairs of observers."),
                        br(),
                        h4("Agreement between several observers"),
                        p("The agreement indexes are weighted means of the pairwise indexes. Pairs with a higher proportion of participants 
                        rated positive have a higher weight in the multi-observer agreement index. They are not the mean of the pairwise agreement indexes, 
                        except when all pairs have the same propensity to classify participants as positive. In that case, all weights will be equal.")
                        
                      ))
  ))



# Define server logic required to draw a histogram
server <- function(input, output, session) {

  data <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    df1<-read.csv(inFile$datapath,header = input$header,sep = input$sep)
    updateCheckboxGroupInput(session, "variables", label = "Select the observers", choices = colnames(df1), selected = colnames(df1))
    return(df1)
  })
  
  
  observeEvent(data(), {
    updateSelectInput(session, "row", choices = names(data()))
    updateSelectInput(session, "col", choices = names(data()))
    
  })
  
  data2 <- reactive({
    req(data(),input$variables)
    df2 <- data() %>% select(input$variables)
    return(df2)
  })
  
  
  output$contents <- renderTable({
    req(data())
    head(data())
  })
  
  output$summary <- renderPrint({
    
    req(data())
    #compute unique levels in data frame
    lvls <- unique(unlist(data()))
    
    # apply the summation per value 
    sapply(data(),function(x) table(factor(x, levels = lvls,ordered = TRUE)))
    
  })
  
  output$summary2 <- renderPrint({
    
    req(data())
    #compute unique levels in data frame
    lvls <- unique(unlist(data()))
    
    # apply the summation per value 
    sapply(data(),function(x) round(prop.table(table(factor(x, levels = lvls,ordered = TRUE))),2))
    
  })
  
 #proportion of agreement
  
  output$plotpair1 <- renderPlot({
    
    req(data2())
    
    correlo1(data2())
    
  })
  
  output$pom <- renderUI({
    req(data2())
    
    po<-round(qow(data2(),"binary_agree"),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("po"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))))
    
  }) 
  
  #proportion of disagreement
  output$plotpair2 <- renderPlot({
    
    req(data2())
    
    correlo2(data2())
    
  })
  
  output$qom <- renderUI({
    req(data2())
    
    po<-round(qow(data2(),"binary_disagree"),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("qo"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))))
    
  }) 
  
  
  #MAE
  output$plotpair3 <- renderPlot({
    
    req(data2())
    
    correlo3(data2())
    
  })
  
  output$MAE <- renderUI({
    req(data2())
    
    po<-round(qow(data2(),"linear"),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("MAE"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))))
    
  }) 
  
  #MSE
  output$plotpair4 <- renderPlot({
    
    req(data2())
    
    correlo4(data2())
    
  })
  
  output$MSE <- renderUI({
    req(data2())
    
    po<-round(qow(data2(),"quadratic"),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("MSE"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))))
    
  }) 
  
  #Linear weighted kappa, different marginals
  output$plotpair5 <- renderPlot({
    
    req(data2())
    
    correlo5(data2())
    
  })
  
  output$LWK2 <- renderUI({
    req(data2())
    
    po<-round(kappaw(data2(),"linear",2),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Linear k (different marginals)"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))))
    
  }) 
  
  #Linear weighted kappa, same marginals
  output$plotpair6 <- renderPlot({
    
    req(data2())
    
    correlo6(data2())
    
  })
  
  output$LWK3 <- renderUI({
    req(data2())
    
    po<-round(kappaw(data2(),"linear",3),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Linear k (same marginals)"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))))
    
  }) 
  
  #quadratic weighted kappa, different marginals
  output$plotpair7 <- renderPlot({
    
    req(data2())
    
    correlo7(data2())
    
  })
  
  output$QWK2 <- renderUI({
    req(data2())
    
    po<-round(kappaw(data2(),"quadratic",2),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Quadratic k (different marginals)"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))))
    
  }) 
  
  
  #quadratic weighted kappa, same marginals
  output$plotpair8 <- renderPlot({
    
    req(data2())
    
    correlo8(data2())
    
  })
  
  output$QWK3 <- renderUI({
    req(data2())
    
    po<-round(kappaw(data2(),"quadratic",3),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Quadratic k (same marginals)"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))))
    
  }) 
  
  #outputs for two observers
  output$contingency <- renderUI({
    req(data())
    
    custom_table<-sjt.xtab(data()[[input$row]],data()[[input$col]],var.labels=c(paste(input$row),paste(input$col)),show.cell.prc=TRUE,show.summary=FALSE)
    HTML(custom_table$knitr)
  })
  
  output$plotagree <- renderPlot({
    
    req(data())
    
    agreementplot(table(data()[[input$row]],data()[[input$col]]),reverse_y=FALSE,)
    
  })
  
  
  output$pom2 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    po<-round(qow(data2,"binary_agree"),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("po"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))
      ))
  })
  
  output$qom2 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    po<-round(qow(data2,"binary_disagree"),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("qo"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))
      ))
  })
  
  output$MAE2 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    po<-round(qow(data2,"linear"),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("MAE"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))
      ))
  })
  
  output$MSE2 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    po<-round(qow(data2,"quadratic"),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("MSE"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))
      ))
  })
  
  output$LWK22 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    po<-round(kappaw(data2,"linear",2),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Linear weighted kappa (different marginals)"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))
      ))
  })
  
  output$LWK32 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    po<-round(kappaw(data2,"linear",3),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Linear weighted kappa (same marginals)"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))
      ))
  }) 
  
  output$QWK22 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    po<-round(kappaw(data2,"quadratic",2),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Quadratic weighted kappa (different marginals)"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))
      ))
  }) 
  
  output$QWK32 <- renderUI({
    req(data())
    
    data2<-as.data.frame(cbind(data()[[input$row]],data()[[input$col]]))
    
    po<-round(kappaw(data2,"quadratic",3),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("Value"),
        tags$th("SE"),
        tags$th("95% confidence interval")
      ),
      tags$tr(
        tags$td("Quadratic weighted kappa (same marginals)"),
        tags$td(paste(po[1])),
        tags$td(paste(po[2])),
        tags$td(paste(paste(paste(round(po[1]-1.96*po[2],3)),";")),paste(round(po[1]+1.96*po[2],3)))
      ))
  }) 
  
  output$CI_PA <- renderUI({
    tagList(
        sliderInput("nagree", label = "target agreement:",min = 0.1, value = 0.8, max = 1),
        sliderInput("nwidth_agree", label = "target width:",min = 0.05, value = 0.1, max = 0.3)
    )
  })
  output$CI_PD <- renderUI({
    tagList(
      sliderInput("ndisagree", label = "target disagreement:",min = 0.1, value = 0.8, max = 1),
      sliderInput("nwidth_disagree", label = "target width:",min = 0.05, value = 0.1, max = 0.3)
    )
  })
  
  output$CI_MAE <- renderUI({
    tagList(
      sliderInput("nMAE", label = "target MAE:",min = 0.0, value = 0.5, max = as.numeric(req(input$ncat))-1),
      sliderInput("nwidth_MAE", label = "target width:",min = 0.05, value = 0.1, max = (as.numeric(req(input$ncat))-1)/2)
    )
  })
  
  output$CI_MSE <- renderUI({
    tagList(
      sliderInput("nMSE", label = "target MSE:",min = 0.0, value = 0.5, max = as.numeric(req(input$ncat))-1),
      sliderInput("nwidth_MSE", label = "target width:",min = 0.05, value = 0.1, max = (as.numeric(req(input$ncat))-1)/2)
    )
  })
  
  output$CI_LINEAR <- renderUI({
    tagList(
      sliderInput("nlin", label = "target agreement:",min = 0.1, value = 0.8, max = 1),
      sliderInput("nwidth_lin", label = "target width:",min = 0.05, value = 0.1, max = 0.3)
    )
  })
  
  output$CI_QUADRATIC <- renderUI({
    tagList(
      sliderInput("nquad", label = "target agreement:",min = 0.1, value = 0.8, max = 1),
      sliderInput("nwidth_quad", label = "target width:",min = 0.05, value = 0.1, max = 0.3)
    )
  })
  
  
  
  new_categories<-reactive({
    n<-input$ncat 
    if(n>1){
      lapply(2:n, function(i) {
        br()
        numericInput(inputId = paste0("category_",i),label = paste0("Put a proportion for category (", i,"):"),0.5,min=0,max=1,step=0.1)
      }) # end of the function and lapply
    }     # end of the if(n>1)
    
  })     # end of the reactive
  
 
  
  output$new_categories <- renderUI({ new_categories() })
  
  
  operation <- function(which_coef){
    
      marginals<-vector()
      for (i in 1:req(input$ncat)){
      marginals[i]<-req(input[[paste0("category_",i)]])
      }
      
      
        if (which_coef==1){
        data_sim<-sizeCI_qow(wqo=as.numeric(req(input$nagree)),weight="binary_agree",nrat=as.numeric(req(input$nrater)),npar=as.numeric(req(input$nsub)),prob=as.numeric(marginals),alpha=as.numeric(req(input$alpha)),nsim=as.numeric(req(input$nsim)),seed=as.numeric(req(input$seed)))
        }
      if (which_coef==2){
        data_sim<-sizeCI_qow(wqo=as.numeric(req(input$ndisagree)),weight="binary_disagree",nrat=as.numeric(req(input$nrater)),npar=as.numeric(req(input$nsub)),prob=as.numeric(marginals),alpha=as.numeric(req(input$alpha)),nsim=as.numeric(req(input$nsim)),seed=as.numeric(req(input$seed)))
      }
      if (which_coef==3){
        data_sim<-sizeCI_qow(wqo=as.numeric(req(input$nMAE)),weight="linear",nrat=as.numeric(req(input$nrater)),npar=as.numeric(req(input$nsub)),prob=as.numeric(marginals),alpha=as.numeric(req(input$alpha)),nsim=as.numeric(req(input$nsim)),seed=as.numeric(req(input$seed)))
      }
      if (which_coef==4){
        data_sim<-sizeCI_qow(wqo=as.numeric(req(input$nMSE)),weight="quadratic",nrat=as.numeric(req(input$nrater)),npar=as.numeric(req(input$nsub)),prob=as.numeric(marginals),alpha=as.numeric(req(input$alpha)),nsim=as.numeric(req(input$nsim)),seed=as.numeric(req(input$seed)))
      }
      if (which_coef==5){
        data_sim<-sizeCI_kappaw(wkappa=as.numeric(req(input$nlin)),weight="linear",nrat=as.numeric(req(input$nrater)),npar=as.numeric(req(input$nsub)),prob=as.numeric(marginals),alpha=as.numeric(req(input$alpha)),nsim=as.numeric(req(input$nsim)),seed=as.numeric(req(input$seed)))
      }
      if (which_coef==6){
        data_sim<-sizeCI_kappaw(wkappa=as.numeric(req(input$nquad)),weight="quadratic",nrat=as.numeric(req(input$nrater)),npar=as.numeric(req(input$nsub)),prob=as.numeric(marginals),alpha=as.numeric(req(input$alpha)),nsim=as.numeric(req(input$nsim)),seed=as.numeric(req(input$seed)))
      }
       #return(c(marginals,input$nagree,input$nrater,input$nsub,input$alpha,input$nsim,input$seed))
      
      return(data_sim)
  }
  
  
  
  dat_sim<-eventReactive(input$size_GO,{
    
    
    withProgress(message = 'Simulations in progress', value = 0,{operation(req(input$size_agree))})
  })
  
  
  output$size_CIagree <- renderUI({
    
    req(dat_sim())
    
    res<-round(summary(dat_sim()),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("min"),
        tags$th("Q1"),
        tags$th("Median"),
        tags$th("Mean"),
        tags$th("Q3"),
        tags$th("Max")
      ),
      tags$tr(
        tags$td("Confidence interval width"),
        tags$td(paste(res[1])),
        tags$td(paste(res[2])),
        tags$td(paste(res[3])),
        tags$td(paste(res[4])),
        tags$td(paste(res[5])),
        tags$td(paste(res[6]))
      ))
  }) 
  
  
  output$size_CIagree3 <- renderUI({
    
    req(dat_sim())
    
    percent<-vect()
    
    if (req(input$size_agree)==1){percent<-mean(ifelse(dat_sim()>as.numeric(input$nwidthagree),1,0))}
    if (req(input$size_agree)==2){percent<-mean(ifelse(dat_sim()>as.numeric(input$nwidthagree),1,0))}
    if (req(input$size_agree)==3){percent<-mean(ifelse(dat_sim()>as.numeric(input$nwidthagree),1,0))}
    if (req(input$size_agree)==4){percent<-mean(ifelse(dat_sim()>as.numeric(input$nwidthagree),1,0))}
    if (req(input$size_agree)==5){percent<-mean(ifelse(dat_sim()>as.numeric(input$nwidthagree),1,0))}
    if (req(input$size_agree)==6){percent<-mean(ifelse(dat_sim()>as.numeric(input$nwidthagree),1,0))}
    
    
    
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("min"),
        tags$th("Q1"),
        tags$th("Median"),
        tags$th("Mean"),
        tags$th("Q3"),
        tags$th("Max")
      ),
      tags$tr(
        tags$td("Confidence interval width"),
        tags$td(paste(res[1])),
        tags$td(paste(res[2])),
        tags$td(paste(res[3])),
        tags$td(paste(res[4])),
        tags$td(paste(res[5])),
        tags$td(paste(res[6]))
      ))
  }) 
  
  output$size_CIagree2 <- renderPlot({
    
    req(dat_sim())
    
    plot(density(dat_sim()),main="Density plot of the width of the confidence interval", xlab="Width of the confidence interval",ylab="")
    
   
  }) 
  

############## TESTING APPROACH
  output$TEST_PA <- renderUI({
    tagList(
      sliderInput("agree0", label = "Proportion of agreement under H0:",min = 0, value = 0.7, max = 1),
      sliderInput("agree1", label = "Proportion of agreement under H1:",min = 0, value = 0.8, max = 1),
 
    )
  })
  output$TEST_PD <- renderUI({
    tagList(
      sliderInput("disagree0", label = "Proportion of disagreement under H0:",min = 0, value = 0.3, max = 1),
      sliderInput("disagree1", label = "Proportion of disagreement under H1:",min = 0, value = 0.2, max = 1),
    )
  })
  
  output$TEST_MAE <- renderUI({
    tagList(
      sliderInput("MAE0", label = "Mean absolute error under H0:",min = 0, value = 0.3, max = as.numeric(req(input$tncat))-1),
      sliderInput("MAE1", label = "Mean absolute error under H1:",min = 0, value = 0.2, max = as.numeric(req(input$tncat))-1),
    )
  })
  
  output$TEST_MSE <- renderUI({
    tagList(
      sliderInput("MSE0", label = "Mean squared error under H0:",min = 0, value = 0.3, max = as.numeric(req(input$tncat))-1),
      sliderInput("MSE1", label = "Mean squared error under H1:",min = 0, value = 0.2, max = as.numeric(req(input$tncat))-1),
    )
  })
  
  output$TEST_LINEAR <- renderUI({
    tagList(
      sliderInput("lin0", label = "Linear weighted kappa under H0:",min = 0, value = 0.7, max = 1),
      sliderInput("lin1", label = "Linear weighted kappa under H1:",min = 0, value = 0.8, max = 1),
    )
  })
  
  output$TEST_QUADRATIC <- renderUI({
    tagList(
      sliderInput("qua0", label = "Quadratic weighted kappa under H0:",min = 0, value = 0.7, max = 1),
      sliderInput("qua1", label = "Quadratic weighted kappa under H1:",min = 0, value = 0.8, max = 1),
    )
  })
  
  
  
  tnew_categories<-reactive({
    n<-input$tncat 
    if(n>1){
      lapply(2:n, function(i) {
        br()
        numericInput(inputId = paste0("tcategory_",i),label = paste0("Put a proportion for category (", i,"):"),0.5,min=0,max=1,step=0.1)
      }) # end of the function and lapply
    }     # end of the if(n>1)
    
  })     # end of the reactive
  
  
  
  output$tnew_categories <- renderUI({ tnew_categories() })
  
  
  toperation <- function(twhich_coef){
    
    tmarginals<-vector()
    for (i in 1:req(input$tncat)){
      tmarginals[i]<-req(input[[paste0("tcategory_",i)]])
    }
    
    
    if (twhich_coef==1){
      tdata_sim<-test_qow(wqo0=as.numeric(req(input$agree0)),wqo1=as.numeric(req(input$agree1)),weight="binary_agree",nrat=as.numeric(req(input$tnrater)),npar=as.numeric(req(input$tnsub)),prob=as.numeric(tmarginals),alpha=as.numeric(req(input$talpha)),nsim=as.numeric(req(input$tnsim)),seed=as.numeric(req(input$tseed)))
    }
    if (twhich_coef==2){
      tdata_sim<-test_qow(wqo0=as.numeric(req(input$disagree0)),wqo1=as.numeric(req(input$disagree1)),weight="binary_disagree",nrat=as.numeric(req(input$tnrater)),npar=as.numeric(req(input$tnsub)),prob=as.numeric(tmarginals),alpha=as.numeric(req(input$talpha)),nsim=as.numeric(req(input$tnsim)),seed=as.numeric(req(input$tseed)))
    }
    if (twhich_coef==3){
      tdata_sim<-test_qow(wqo0=as.numeric(req(input$MAE0)),wqo1=as.numeric(req(input$MAE1)),weight="linear",nrat=as.numeric(req(input$tnrater)),npar=as.numeric(req(input$tnsub)),prob=as.numeric(tmarginals),alpha=as.numeric(req(input$talpha)),nsim=as.numeric(req(input$tnsim)),seed=as.numeric(req(input$tseed)))
    }
    if (twhich_coef==4){
      tdata_sim<-test_qow(wqo0=as.numeric(req(input$MSE0)),wqo1=as.numeric(req(input$MSE1)),weight="quadratic",nrat=as.numeric(req(input$tnrater)),npar=as.numeric(req(input$tnsub)),prob=as.numeric(tmarginals),alpha=as.numeric(req(input$talpha)),nsim=as.numeric(req(input$tnsim)),seed=as.numeric(req(input$tseed)))
    }
    if (twhich_coef==5){
      tdata_sim<-test_kappaw(wkappa0=as.numeric(req(input$MSE0)),wkappa1=as.numeric(req(input$MSE1)),weight="linear",nrat=as.numeric(req(input$tnrater)),npar=as.numeric(req(input$tnsub)),prob=as.numeric(tmarginals),alpha=as.numeric(req(input$talpha)),nsim=as.numeric(req(input$tnsim)),seed=as.numeric(req(input$tseed)))
    }
    if (twhich_coef==6){
      tdata_sim<-test_kappaw(wkappa0=as.numeric(req(input$MSE0)),wkappa1=as.numeric(req(input$MSE1)),weight="quadratic",nrat=as.numeric(req(input$tnrater)),npar=as.numeric(req(input$tnsub)),prob=as.numeric(tmarginals),alpha=as.numeric(req(input$talpha)),nsim=as.numeric(req(input$tnsim)),seed=as.numeric(req(input$tseed)))
    }
    
    return(tdata_sim)
  }
  
  
  
  tdat_sim<-eventReactive(input$tsize_GO,{
    
    
    withProgress(message = 'Simulations in progress', value = 0,{toperation(req(input$tsize_agree))})
  })
  
  
  output$size_tagree <- renderUI({
    
    req(tdat_sim())
    
    tres<-round(summary(tdat_sim()),3)
    
    tags$table(
      style = "width:100%",
      tags$tr(
        tags$th(""),
        tags$th("min"),
        tags$th("Q1"),
        tags$th("Median"),
        tags$th("Mean"),
        tags$th("Q3"),
        tags$th("Max"),
        tags$th(paste0(req(input$tncover),"%"))
      ),
      tags$tr(
        tags$td("Confidence interval width"),
        tags$td(paste(tres[1])),
        tags$td(paste(tres[2])),
        tags$td(paste(tres[3])),
        tags$td(paste(tres[4])),
        tags$td(paste(tres[5])),
        tags$td(paste(tres[6])),
        tags$td(paste(round(quantile(tdat_sim(),req(input$tncover)/100),3)))
      ))
  }) 
  
  output$size_tagree2 <- renderPlot({
    
    req(tdat_sim())
    
    plot(density(tdat_sim()),main="Density plot of the empirical power", xlab="Empirical power",ylab="")
    
    
  }) 
  
    
                                          
}

# Run the application 
shinyApp(ui = ui, server = server)
