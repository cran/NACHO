suppressPackageStartupMessages(invisible(
  sapply(
    X = c("shiny", "shinyWidgets", "utils", "ggplot2", "NACHO"),
    FUN = library, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE
  )
))

source("utils.R")

nacho_object <- shiny::getShinyOption("nacho_object", NULL)

ui <- shiny::tagList(
  shiny::tags$head(shiny::tags$style(shiny::HTML(
    ".navbar-nav { float: none !important; } .navbar-nav > li:nth-child(7) { float: right; }"
  ))),
  shiny::navbarPage(
    theme = "united-bootstrap.min.css",
    title = shiny::tags$span(shiny::tags$img(src = "nacho_hex.png", height = 18), "NACHO"),
    windowTitle = "NACHO",
    collapsible = TRUE,
    id = "main-menu",
    selected = "qc_metrics",
    shiny::tabPanel("Upload RCC Files", icon = shiny::icon("file-upload"), value = "upload-tab",
      shiny::fluidRow(
        shiny::column(width = 6,
          card(title = "Normalisation Settings", body = {
            shiny::radioButtons("norm_method",
              label = shiny::tags$span(
                "Normalisation Method",
                shiny::helpText("(See", shiny::tags$a(
                  href = "https://m.canouil.dev/NACHO/reference/load_rcc.html",
                  target = "_blank",
                  shiny::tags$code("load_rcc()")
                ), "for details and more options)")
              ),
              choices = c("GEO", "GLM"), selected = "GEO",
              inline = TRUE
            )
          })
        ),
        shiny::column(width = 6,
          card(title = "Upload RCC Files", body = {
            shiny::fileInput("rcc_files",
              label = shiny::tags$span(
                "Choose One or Several RCC Files and Optionally a CSV File",
                shiny::helpText(
                  "(The CSV file must contains \"IDFILE\",",
                  shiny::tags$i("i.e."),
                  "\"BASENAME.RCC\", and optionnally \"plexset_id\",",
                  shiny::tags$i("i.e."),
                  ", \"S1\", \"S2\", ...)"
                )
              ),
              multiple = TRUE,
              accept = c(".RCC", "application/zip", ".csv")
            )
          })
        )
      ),
      shiny::uiOutput("upload_ui")
    ),
    panelInputUI("qc_metrics", "QC Metrics",
      sidebar = {
        list(
          shiny::tags$h3("QC Thresholds"),
          shiny::radioButtons("qc_bd_metrics",
            shiny::tags$span("Binding Density",
              shiny::actionLink("about_bd", NULL, icon = shiny::icon("info-circle"), style = "text-decoration:none;")
            ),
            choiceNames = list(
              shiny::tags$span("MAX/FLEX", shiny::helpText("(Default: 0.1 - 2.25)")),
              shiny::tags$span("SPRINT", shiny::helpText("(Default: 0.1 - 1.8)"))
            ),
            choiceValues = list("MAX/FLEX", "SPRINT"),
            inline = TRUE
          ),
          shiny::sliderInput("qc_bd_thresh", NULL,
            min = 0, max = 2.5, value = c(0.1, 2.25), step = 0.05
          ),
          shiny::sliderInput("qc_fov_thresh",
            shiny::tags$span("Field of View",
              shiny::actionLink("about_fov", NULL, icon = shiny::icon("info-circle"), style = "text-decoration:none;"),
              shiny::helpText("(Default: 75)")
            ),
            min = 50, max = 100, value = 75
          ),
          shiny::sliderInput("qc_pcl_thresh",
            shiny::tags$span("Positive Control Linearity",
              shiny::actionLink("about_pcl", NULL, icon = shiny::icon("info-circle"), style = "text-decoration:none;"),
              shiny::helpText("(Default: 0.95)")
            ),
            min = 0.5, max = 1, value = 0.95
          ),
          shiny::sliderInput("qc_lod_thresh",
            shiny::tags$span("Limit of Detection",
              shiny::actionLink("about_lod", NULL, icon = shiny::icon("info-circle"), style = "text-decoration:none;"),
              shiny::helpText("(Default: 2)")
            ),
            min = 0, max = 30, value = 2
          )
        )
      },
      plotInputUI("Binding Density"),
      plotInputUI("Field of View", right = TRUE),
      plotInputUI("Positive Control Linearity"),
      plotInputUI("Limit of Detection", right = TRUE)
    ),
    panelInputUI("qc_control", "Control Probes",
      sidebar = NULL,
      plotInputUI("Positive Probes"),
      plotInputUI("Negative Probes", right = TRUE),
      plotInputUI("Housekeeping Genes Probes"),
      plotInputUI("Control Probe Expression", right = TRUE)
    ),
    panelInputUI("qc_count", "Counts",
      sidebar = NULL,
      plotInputUI("Average Count vs. Binding Density"),
      plotInputUI("Average Count vs. Median Count", right = TRUE),
      plotInputUI("Principal Component Analysis"),
      plotInputUI("Principal Component Analysis Inertia", right = TRUE)
    ),
    panelInputUI("norm", "Normalisation",
      sidebar = {
        list(
          shiny::sliderInput("qc_pf_thresh",
            shiny::tags$span("Positive Factor",
              shiny::actionLink("about_pf", NULL, icon = shiny::icon("info-circle"), style = "text-decoration:none;"),
              shiny::helpText("(Default: 0.25 - 4)")
            ),
            min = 0, max = 5, value = c(0.25, 4), step = 0.25
          ),
          shiny::sliderInput("qc_hgf_thresh",
            shiny::tags$span("Housekeeping Genes Factor",
              shiny::actionLink("about_hgf", NULL, icon = shiny::icon("info-circle"), style = "text-decoration:none;"),
              shiny::helpText("(Default: 0.09 - 11)")
            ),
            min = 0, max = 15, value = c(0.09, 11), step = 0.01
          )
        )
      },
      plotInputUI("Positive Factor vs. Negative Factor"),
      plotInputUI("Housekeeeping Genes Factor", right = TRUE),
      plotInputUI("Normalisation Result", right = TRUE)
    ),
    shiny::tabPanel(title = "Outliers", value = "outliers-tab",
      card(
        title = shiny::tags$h4("Outliers List"),
        body = list(shiny::uiOutput("outliers-thresholds"), shiny::tableOutput("outliers"))
      )
    ),
    shiny::tabPanel("About", icon = shiny::icon("info"), value = "about-tab",
      shiny::tags$p(shiny::includeMarkdown("www/about-nacho.md"))
    )
  )
)

server <- function(input, output, session) {
  # ---------------------------------------- Upload
  nacho_react <- shiny::reactive({
    if (inherits(nacho_object, "nacho")) return(nacho_object)

    targets <- shiny::req(input$rcc_files)
    if (nrow(targets) > 0) {
      targets <- do.call("rbind.data.frame", apply(
        X = targets[, c("name", "datapath", "type")],
        MARGIN = 1,
        FUN = function(.row) {
          name <- .row[1]
          datapath <- .row[2]
          type <- .row[3]
          if (type == "application/x-zip-compressed") {
            ex_dir <- file.path(dirname(datapath), sub(".zip$", "", name))
            utils::unzip(datapath, exdir = ex_dir)
            data.frame(
              name = file.path(sub(".zip$", "", name), list.files(ex_dir)),
              datapath = list.files(ex_dir, full.names = TRUE),
              type = type,
              IDFILE = file.path(sub(".zip$", "", name), list.files(ex_dir)),
              stringsAsFactors = FALSE
            )
          } else {
            file.rename(
              from = datapath,
              to = file.path(dirname(datapath), name)
            )
            data.frame(
              name = name,
              datapath = file.path(dirname(datapath), name),
              type = type,
              IDFILE = name,
              stringsAsFactors = FALSE
            )
          }
        }
      ))

      targets_ssheet <- targets[grep("\\.csv", targets[["name"]]), ]
      if (nrow(targets_ssheet) > 0) {
        targets <- targets[grep("\\.RCC", targets[["name"]]), ]
        ssheet_dt <- data.table::fread(targets_ssheet[["datapath"]])
      }

      check_multiplex <- all(sapply(
        X = targets$datapath,
        FUN = function(.x) any(grepl("Endogenous8s", readLines(.x)))
      ))
      save(list = ls(), file = "all.rdata")
      if (check_multiplex) {
        targets <- merge(
          x = targets,
          y = expand.grid(
            IDFILE = targets[["IDFILE"]],
            plexset_id = paste0("S", seq_len(8)),
            stringsAsFactors = FALSE
          ),
          by = "IDFILE"
        )
      }

      if (any(grepl("^IDFILE$", names(ssheet_dt)))) {
        if (check_multiplex) {
          if (any(grepl("^plexset_id$", names(ssheet_dt)))) {
            targets <- merge(
              x = targets,
              y = ssheet_dt,
              by = c("IDFILE", "plexset_id")
            )
          } else {
            warning(
              "[NACHO] Missing \"plexset_id\" column in sample sheet file!\n",
              "  Sample sheet file is discarded."
            )
          }
        } else {
          targets <- merge(
            x = targets,
            y = ssheet_dt,
            by = "IDFILE"
          )
        }
      } else {
        warning(
          "[NACHO] Missing \"IDFILE\" column in sample sheet file!\n",
          "  Sample sheet file is discarded."
        )
      }
      suppressMessages(
        x = NACHO::load_rcc(
          data_directory = unique(mapply(FUN = function(.x, .y) sub(.x, "", .y), targets$IDFILE, targets$datapath)),
          ssheet_csv = targets,
          id_colname = "IDFILE",
          normalisation_method = input[["norm_method"]]
        )
      )
    }
  })

  output$rcc_contents <- shiny::renderTable(shiny::req(input$rcc_files))
  output$rcc_contents_summary <- shiny::renderUI({
    rcc_size <- sum(input$rcc_files[, "size"])
    class(rcc_size) <- "object_size"
    shiny::tags$p(
      "A total of", shiny::tags$strong(length(input$rcc_files[, "name"])), "RCC files were succesfully uploaded,",
      "for a total amount of", shiny::tags$strong(format(rcc_size, units = "Mb")), "."
    )
  })
  output$upload_ui <- shiny::renderUI({
    if (is.null(input$rcc_files)) {
      shiny::fluidRow(style = "padding-top: 1em;",
        shiny::column(width = 12,
          card(title = "Summary", body = shiny::tags$p("No RCC files uploaded."))
        )
      )
    } else {
      list(
        shiny::fluidRow(style = "padding-top: 1em;",
          shiny::column(width = 12,
            card(title = "Summary", body = shiny::uiOutput("rcc_contents_summary"))
          )
        ),
        shiny::fluidRow(style = "padding-top: 1em;",
          shiny::column(width = 12,
            card(title = "RCC Files Uploaded", body = shiny::tableOutput("rcc_contents"))
          )
        )
      )
    }
  })

  # ---------------------------------------- UI / SERVER
  # Global UI input
  shiny::observe({
    nacho_tmp <- nacho_custom()
    lapply(
      X = c(
        "bd", "fov", "pcl", "lod",
        "pp", "np", "hgp", "cpe",
        "acvbd", "acvmc", "pca", "pcai",
        "pfvnf", "hgf", "nr"
      ),
      F = plotInput,
      nacho = nacho_tmp
    )
  })

  # QC metrics UI input
  ## Update UI with thresholds
  shiny::observe({
    switch(shiny::req(input$qc_bd_metrics),
      "MAX/FLEX" = {
        shiny::updateSliderInput(session, "qc_bd_thresh",
          value = min(shiny::isolate(input$qc_bd_thresh), 2.25)
        )
      },
      "SPRINT" = {
        shiny::updateSliderInput(session, "qc_bd_thresh",
          value = min(shiny::isolate(input$qc_bd_thresh), 1.8)
        )
      }
    )
    shiny::updateSliderInput(session, "qc_fov_thresh",
      value = shiny::isolate(input$qc_fov_thresh)
    )
    shiny::updateSliderInput(session, "qc_pcl_thresh",
      value = shiny::isolate(input$qc_pcl_thresh)
    )
    shiny::updateSliderInput(session, "qc_lod_thresh",
      value = shiny::isolate(input$qc_lod_thresh)
    )
  })

  ## Help for QC metrics
  lapply(
    X = c(
      "Binding Density",
      "Field of View",
      "Positive Control Linearity",
      "Limit of Detection",
      "Positive Factor",
      "Housekeeping Genes Factor"
    ),
    FUN = function(.x) {
      short_x <- tolower(sub("\\b(\\pL)\\pL|.", "\\U\\1", .x, perl = TRUE))
      shiny::observeEvent(input[[paste0("about_", short_x)]], {
        shiny::showModal(shiny::modalDialog(
          title = .x,
          shiny::tags$p(shiny::includeMarkdown(paste0("www/about-", short_x, ".md"))),
          easyClose = TRUE
        ))
      })
    }
  )

  # ---------------------------------------- Input
  # Get nacho object and update thresholds
  nacho_custom <- shiny::reactive({
    nacho <- shiny::req(nacho_react())

    nacho$outliers_thresholds[["BD"]] <- input$qc_bd_thresh %||%
      nacho$outliers_thresholds[["BD"]]
    nacho$outliers_thresholds[["FoV"]] <- input$qc_fov_thresh %||%
      nacho$outliers_thresholds[["FoV"]]
    nacho$outliers_thresholds[["LoD"]] <- input$qc_lod_thresh %||%
      nacho$outliers_thresholds[["LoD"]]
    nacho$outliers_thresholds[["PCL"]] <- input$qc_pcl_thresh %||%
      nacho$outliers_thresholds[["PCL"]]
    nacho$outliers_thresholds[["Positive_factor"]] <- input$qc_pf_thresh %||%
      nacho$outliers_thresholds[["Positive_factor"]]
    nacho$outliers_thresholds[["House_factor"]] <- input$qc_hgf_thresh %||%
      nacho$outliers_thresholds[["House_factor"]]

    NACHO::check_outliers(nacho)
  })
  observe({
    if (inherits(nacho_object, "nacho")) {
      nacho_object$outliers_thresholds <- nacho_custom()$outliers_thresholds
      message(
        '[NACHO] Updated "nacho_object" can be loaded with:\n',
        '  nacho_object <- readRDS("', tempdir(), '/nacho_object.rds")'
      )
      saveRDS(object = nacho_object, file = file.path(tempdir(), "nacho_object.rds"))
    }
  })

  # ---------------------------------------- Output
  outliers_list <- shiny::reactive({
    columns_qc <- intersect(c(
      "IDFILE", "CartridgeID", "BD", "FoV", "PCL", "LoD", "MC", "MedC",
      "Positive_factor", "House_factor"
    ), colnames(nacho_custom()$nacho))
    unique(
      data.table::as.data.table(nacho_custom()$nacho)[
        (is_outlier),
        .SD,
        .SDcols = columns_qc
      ]
    )
  })
  output[["outliers"]] <- shiny::renderTable(outliers_list())
  output[["outliers-thresholds"]] <- shiny::renderUI({
    ot <- lapply(nacho_custom()$outliers_thresholds, round, digits = 3)
    shiny::tags$div(
      shiny::tags$ul(
        shiny::tags$li(
          "Binding Density (", shiny::tags$code("BD"), ") <",
          shiny::tags$strong(min(ot[["BD"]])),
          "or Binding Density (", shiny::tags$code("BD"), ") >",
          shiny::tags$strong(max(ot[["BD"]]))
        ),
        shiny::tags$li(
          "Field of View (", shiny::tags$code("FoV"), ") <",
          shiny::tags$strong(ot[["FoV"]])
        ),
        shiny::tags$li(
          "Positive Control Linearity (", shiny::tags$code("PCL"), ") <",
          shiny::tags$strong(min(ot[["PCL"]]))
        ),
        shiny::tags$li(
          "Limit of Detection (", shiny::tags$code("LoD"), ") <",
          shiny::tags$strong(min(ot[["LoD"]]))
        ),
        shiny::tags$li(
          "Positive Normalisation Dactor (", shiny::tags$code("Positive_factor"),
          ") <", shiny::tags$strong(min(ot[["Positive_factor"]])),
          "or Positive Normalisation Dactor (", shiny::tags$code("Positive_factor"),
          ") >", shiny::tags$strong(max(ot[["Positive_factor"]]))),
        shiny::tags$li(
          "Housekeeping Normalisation Factor (", shiny::tags$code("house_factor"),
          ") <", shiny::tags$strong(min(ot[["House_factor"]])),
          "or Housekeeping Normalisation Dactor (", shiny::tags$code("house_factor"),
          ") >", shiny::tags$strong(max(ot[["House_factor"]]))
        )
      )
    )
  })

  # ---------------------------------------- Show / Hide tabs
  shiny::observe({
    if (!inherits(nacho_object, "nacho") & is.null(input$rcc_files)) {
      shiny::showTab("main-menu", target = "upload-tab", select = TRUE)
      lapply(
        X = paste0(c("qc_metrics", "qc_control", "qc_count", "norm", "outliers"), "-tab"),
        FUN = function(.x) shiny::hideTab("main-menu", target = .x)
      )
    }

    if (inherits(nacho_object, "nacho") & is.null(input$rcc_files)) {
      lapply(
        X = paste0(c("qc_metrics", "qc_control", "qc_count", "norm", "outliers"), "-tab"),
        FUN = function(.x) shiny::showTab("main-menu", target = .x, select = .x == "qc_metrics")
      )
      shiny::hideTab("main-menu", target = "upload-tab")
    }

    if (!is.null(input$rcc_files)) {
      lapply(
        X = paste0(c("upload", "qc_metrics", "qc_control", "qc_count", "norm", "outliers"), "-tab"),
        FUN = function(.x) shiny::showTab("main-menu", target = .x, select = .x == "qc_metrics")
      )
    }

    if (nrow(outliers_list()) == 0) {
      shiny::hideTab("main-menu", target = "outliers-tab")
    } else {
      shiny::showTab("main-menu", target = "outliers-tab")
    }
  })
}

shiny::shinyApp(ui = ui, server = server)
