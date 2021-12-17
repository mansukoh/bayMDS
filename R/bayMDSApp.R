#' @title Shiny App for exploring the results of \code{bmds} function
#' 
#' @description Call Shiny to show the results of Bayesian analysis of
#' multidimensional scaling in a web-based application.
#' @usage bayMDSApp(out)
#' @param out an object of class \code{bmds}, the output of the \code{bmds} function
#' @return open Shiny app
#' @export
#' @examples
#' \donttest{
#' data(cityDIST)
#' out <- bmds(cityDIST, min_p=1, max_p=6 )
#' if(interactive()){bayMDSApp(out)}
#' }

bayMDSApp<-function(out){

  bayMDS.env<-new.env()
 
  Plist = out$min_p:out$max_p
  optP=Plist[which.min(MDSIC(out,plot=FALSE)$mdsic)]
  bayMDS_app=shiny::shinyApp(
    ui=shiny::fluidPage(
      theme=shinythemes::shinytheme("cerulean"),
      # Application title
      shiny::titlePanel(shiny::h1(shiny::strong("Explore bayMDS output"))),

      # Sidebar with a slider input for number of bins
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::radioButtons("type",
                       label = shiny::h3((shiny::strong("# of dimensions"))),
                       choices=as.character(Plist),selected=as.character(optP)),
          width=3),

        
        shiny::mainPanel(
          shiny::tabsetPanel(
            shiny::tabPanel(shiny::h4(shiny::strong("Object Configuration")),
                            shiny::plotOutput("show_obj")),            
            shiny::tabPanel(shiny::h4(shiny::strong("MDSIC")),
                            shiny::plotOutput("show_mdsic"),
                            shiny::tableOutput("show_mdsic_table")),
            shiny::tabPanel(shiny::h4(shiny::strong("Delta vs. DIST")),
                            shiny::plotOutput("plot_deldist")),            
            shiny::tabPanel(shiny::h4(shiny::strong("Trace-Delta")),
                            shiny::plotOutput("trace_delta")),
            shiny::tabPanel(shiny::h4(shiny::strong("Trace-Sigma")),
                            shiny::plotOutput("trace_sigma")),
            shiny::tabPanel(shiny::h4(shiny::strong("Trace-Lambda")),
                            shiny::plotOutput("trace_lambda")))
        )
      )), #end ui

    server<-function(input,output,session){

      output$show_mdsic <- shiny::renderPlot({
        MDSIC(out)
      })
      output$show_mdsic_table <- shiny::renderTable({
        temp=MDSIC(out,plot=FALSE)
        tempT = cbind(out$min_p:out$max_p,temp$mdsic,temp$llike,temp$penalty)
        colnames(tempT) = c("p","MDSIC","LRT","penalty")
        tempT = data.frame(tempT)
        tempT$p = as.character(tempT$p)
        tempT
      },digits=4)
      output$show_obj <- shiny::renderPlot({
        selp=as.numeric(input$type)
        plotObj(out$BMDSp[[selp]],pairs=TRUE)
      })
      output$plot_deldist <- shiny::renderPlot({
        selp=as.numeric(input$type)
        plotDelDist(out$BMDSp[[selp]])
      })      
      output$trace_delta <- shiny::renderPlot({
        selp=as.numeric(input$type)
        plotTrace(out$BMDSp[[selp]],para=c("del"),linecolor="blue")
      })
      output$trace_sigma <- shiny::renderPlot({
        selp=as.numeric(input$type)
        plotTrace(out$BMDSp[[selp]],para=c("sigma"),linecolor="blue")
      })
      output$trace_lambda <- shiny::renderPlot({
        selp=as.numeric(input$type)
        plotTrace(out$BMDSp[[selp]],para=c("lambda"),linecolor="blue")
      })      
    }#end server
  )#end App
  shiny::runApp(bayMDS_app,launch.browser=TRUE)
}



