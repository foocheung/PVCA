#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function() {
  fluidPage(
    titlePanel("PVCA Plot Generator"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file_exprs", "Upload Expression Data"),
        fileInput("file_metadata", "Upload Metadata"),
        sliderInput("width", "Plot Width:", min = 300, max = 1200, value = 600),
        sliderInput("height", "Plot Height:", min = 300, max = 800, value = 400),
        actionButton("submit", "Generate Plot")
      ),
      mainPanel(
        plotOutput("pca_plot")
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "PVCA"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
