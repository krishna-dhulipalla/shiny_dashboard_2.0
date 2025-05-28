# ui.R
source("global.R")

header <- dashboardHeader(
  title = "Gene Co-Expression Analysis",
  titleWidth = 250
)

sidebar <- dashboardSidebar(
  width = 250,
  sidebarMenu(
    id = "tabs",
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("Heatmap", tabName = "heatmap", icon = icon("th"))
  ),
  
  conditionalPanel(
    condition = "input.tabs == 'heatmap'",
    selectizeInput("file", "Select Dataset", 
                   choices = list(
                     "Tissue-specific (NE)" = c(
                       "Cauline" = "Cauline",
                       "flower" = "flower",
                       "leaf" = "leaf",
                       "root" = "root",
                       "shoot" = "shoot",
                       "silique" = "silique",
                       "stem" = "stem"
                     ),
                     "Kelsey" = c(
                       "Kelsey coexpression_bulk" = "Kelsey coexpression_bulk",
                       "Kelsey coexpression_all" = "Kelsey coexpression_all"
                     ),
                     "Tylor" = c(
                       "Tylor coexpression_bulk" = "Tylor coexpression_bulk",
                       "Tylor coexpression_all" = "Tylor coexpression_all"
                     ),
                     "ATH Leaf 2k" = c(
                       "ATH Leaf 2k Pearson" = "ATH Leaf 2k Pearson",
                       "ATH Leaf 2k Spearman" = "ATH Leaf 2k Spearman"
                     ),
                     "General Coexpression" = c(
                       "Coexpression_all" = "Coexpression_all",
                       "Coexpression_bulk" = "Coexpression_bulk"
                     )
                   ),
                   selected = NULL
    ),
    radioButtons(
      inputId = "label_type",
      label = "Display Labels As:",
      choices = c("GeneID", "GeneName"),
      selected = "GeneName",
      inline = TRUE
    ),
    numericInput("num_clusters", "Number of Clusters (cut tree):", value = 5, min = 5, max = 20),
    radioButtons("matrix_type", "Matrix Type:",
                 choices = c("Full Matrix", "Reduced Matrix", "Top N Genes"),
                 selected = "Full Matrix"),
    conditionalPanel(
      condition = "input.matrix_type == 'Reduced Matrix'",
      numericInput("matrix_size", "Submatrix Size", 
                   value = 50, min = 10, max = 10000, step = 50)
    ),
    conditionalPanel(
      condition = "input.matrix_type == 'Top N Genes'",
      numericInput("top_n_genes", "Top N Genes", value = 100, min = 10,step = 10)
    ),
    radioButtons("filter_mode", "Filter Mode:",
                 choices = c("None", "Inside Range", "Outside Range"),
                 selected = "None"),
    conditionalPanel(
      condition = "input.filter_mode == 'Inside Range'",
      sliderInput("inside_range", "Keep Only Correlations Between:",
                  min = -1, max = 1, value = c(-0.5, 0.5), step = 0.05)
    ),
    conditionalPanel(
      condition = "input.filter_mode == 'Outside Range'",
      sliderInput("outside_range", "Keep Only Correlations Outside:",
                  min = -1, max = 1, value = c(-0.5, 0.5), step = 0.05)
    ),
    selectInput(
      inputId = "highlight_category",
      label = "Highlight a Category (default: Top 5 shown):",
      choices = list( "Top5"),
      selected = "Top5"
    ),
    actionButton("submit", "Generate Heatmap"),
    tags$hr()
  )
)

body <- dashboardBody(
  tags$head(
    useShinyjs(),
    tags$style(HTML("

  .main-heatmap {
    overflow: auto;
    position: relative;
    padding: 15px;
    margin: -15px;
  }

  .control-ui-container {
    position: relative !important;
    z-index: 500;
    background: white;
    padding: 10px;
    border-bottom: 1px solid #eee;
  }

  .control_panel {
    max-width: 100% !important;
    margin-top: 10px !important;
    left: 15px !important;
    right: 15px !important;
  }

  .control_panel button,
  .control_panel input {
    margin: 3px;
    padding: 5px 10px;
    font-size: 0.9em;
  }

  .shiny-output-error { color: #ff4444; }

  .loading-message {
    font-size: 18px;
    padding: 20px;
    text-align: center;
  }

  .content-wrapper {
    overflow-x: auto;
  }

  .content {
    width: auto;
  }
  #heatmap_output_heatmap_resize {
    width: auto !important;
    height: auto !important;
  }
  #heatmap_output_sub_heatmap_resize {
    width: auto !important;
    height: auto !important;
  }
  #heatmap_output_heatmap_control, #heatmap_output_sub_heatmap_control {
    width:100% !important;
  }
  #heatmap_tooltip {
      position: fixed;
      background: rgba(0,0,0,0.85);
      color: white;
      padding: 6px 10px;
      border-radius: 5px;
      font-size: 13px;
      pointer-events: none;
      z-index: 9999;
      max-width: 300px;
      white-space: pre-line;
      display: none;
    }
    "))
  ),
  
  tabItems(
    tabItem(
      tabName = "home",
      fluidRow(
        box(
          width = 12,
          title = "Welcome to Gene Co-Expression Analysis",
          status = "primary",
          solidHeader = TRUE,
          tabBox(
            width = 12,
            tabPanel(
              "About",
              h3("Features:"),
              tags$ul(
                tags$li("Interactive heatmap exploration with zoom, brush, and tooltip support"),
                tags$li("Gene search with GeneID or GeneName toggle"),
                tags$li("Sub-heatmap for selected gene regions"),
                tags$li("DotPlot visualization for gene expression"),
                tags$li("Advanced filtering: top genes and correlation range selection"),
                tags$li("Export of selected gene table and heatmap images")
              ),
              tags$img(src = "welcome_page.png", height = "500px", width = "auto")
            ),
            tabPanel(
              "Datasets",
              h3("Available Datasets (NE genes only):"),
              div(class = "dataset-info",
                  tags$ul(
                    tags$li(strong("Cauline:"), "NE gene-specific co-expression data from cauline tissue"),
                    tags$li(strong("Flower:"), "NE gene-specific co-expression data from flower tissue"),
                    tags$li(strong("Leaf:"), "NE gene-specific co-expression data from leaf tissue"),
                    tags$li(strong("Root:"), "NE gene-specific co-expression data from root tissue"),
                    tags$li(strong("Shoot:"), "NE gene-specific co-expression data from shoot tissue"),
                    tags$li(strong("Silique:"), "NE gene-specific co-expression data from silique tissue"),
                    tags$li(strong("Stem:"), "NE gene-specific co-expression data from stem tissue")
                  )
              )
            ),
            tabPanel(
              "Quick Start",
              h3("Getting Started:"),
              tags$ol(
                tags$li("Navigate to the 'Heatmap' tab"),
                tags$li("Select a dataset focused on NE genes"),
                tags$li("Choose matrix size or top N genes if needed"),
                tags$li("Click 'Generate Heatmap' to visualize correlations"),
                tags$li("Use features: zoom, search, subheatmap, DotPlot")
              )
            )
          )
        )
      )
    ),
    tabItem(
      tabName = "heatmap",
      
      # Row 1: Full-width main heatmap
      fluidRow(
        box(
          title = "Main Heatmap",
          width = NULL,  
          solidHeader = TRUE,
          status = "primary",
          originalHeatmapOutput(
            "heatmap_output",
            #output_ui = FALSE,  # 🛑 prevents auto-initialization
            width = "100%",
            action = "hover",
            response = c("hover", "brush"),
            height = 700,
            containment = TRUE
          )
        ),
        div(id = "heatmap_tooltip")
      ),
      
      # Row 2: Sub-heatmap + Heatmap Info side by side
      fluidRow(
        column(width = 6,
               box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
                   plotlyOutput("sub_heatmaply", height = "500px")%>% withSpinner(type = 4)
               )
        ),
        column(width = 6,
               box(title = "Sub-heatmap Mode", width = NULL, solidHeader = TRUE, status = "primary",
                   radioButtons("subheatmap_mode", "Choose Sub-heatmap Mode:",
                                choices = c("Zoom", "Cluster"), selected = "Zoom", inline = TRUE),
                   conditionalPanel(
                     condition = "input.subheatmap_mode == 'Cluster'",
                     selectInput("selected_cluster", "Select Cluster", choices = NULL)
                   )
               ),
               box(title = "Heatmap Info", width = NULL, solidHeader = TRUE, status = "primary",
                   HeatmapInfoOutput("heatmap_output", title = NULL)
               ),
               box(title = "Result table of the selected genes", width = NULL, solidHeader = TRUE, status = "primary",
                   DT::DTOutput("selected_table")
               )
        )
      ),
      
      # Row 3: Dot plot
      fluidRow(
        box(
          title = "Gene Expression Dot Plot",
          width = 12,
          solidHeader = TRUE,
          status = "primary",
          plotOutput("dotplot", height = "500px")
        )
      )
    )
  ),
  tags$script(HTML("
  const tooltip = document.getElementById('heatmap_tooltip');
  const heatmapBox = document.getElementById('heatmap_output_heatmap_resize');
  let lastMouseY = 0;

  document.addEventListener('mousemove', function(e) {
    lastMouseY = e.clientY;
    if (tooltip) {
      tooltip.style.left = (e.clientX + 10) + 'px';
      tooltip.style.top = (e.clientY + 10) + 'px';
    }
  });

  // Hide tooltip when mouse leaves heatmap area
  document.addEventListener('DOMContentLoaded', function () {
    if (heatmapBox && tooltip) {
      heatmapBox.addEventListener('mouseleave', function () {
        tooltip.style.display = 'none';
      });
    }
  });
  
  // Additional safeguard: hide tooltip if page scrolls beyond heatmap
  document.addEventListener('scroll', function () {
    const scrollY = window.scrollY;
    const heatmapBottom = heatmapBox.getBoundingClientRect().bottom;
    // If the mouse is no longer over the heatmap area after scroll
    if (heatmapBottom < 0 || lastMouseY > heatmapBottom) {
      tooltip.style.display = 'none';
    }
  }, true);
"))
)
ui <- dashboardPage(header, sidebar, body)