#' CMMID Shiny for the efficacy of airport screening for infectious diseases
#'
#' @description Runs a local instance of a Shiny app developed by the Centre for
#' the Mathematical Modelling of Infectious Diseases at the London School of
#' Hygiene and tropical Medicine, to model the efficacy of implementing airport
#' screening for infectious diseases, as a part of the CMMID response to the
#' Covid-19 pandemic, in early 2020.
#' @return Runs the `{shiny}` app locally.
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # Choose options from the app GUI
#' \dontrun{
#' run_app()
#' }
run_app <- function() {
  waffle_description <- system.file(
    "info", "waffle_description.md",
    package = "airportscreening"
  )
  density_description <- system.file(
    "info", "density_description.md",
    package = "airportscreening"
  )
  assumptions <- system.file(
    "info", "assumptions.md",
    package = "airportscreening"
  )
  references <- system.file(
    "info", "references.md",
    package = "airportscreening"
  )

  pathogen_parameters <- airportscreening::pathogen_parameters

  ui <- list(
    shiny::tags$style(type = "text/css", "
  body { padding-top: 20px; padding-left: 20px; padding-right: 20px}
  .inline label.control-label, .inline .selectize-control.single {
    display: table-cell;
    text-align: left;
    vertical-align: middle;
  }
  .inline .form-group {
    display: table-row;
  }
  .inline .selectize-control.single div.item {
    padding-right: 15px;
  }
  .seroselect .selectize-control.single { width:10em; }
  .demoselect .selectize-control.single { width:15em; }
  .seroselect.inline label.control-label { width: 70px; }
  .demoselect.inline label.control-label { width: 70px; }
  .costs.inline label.control-label { width: 90px; }
  hr { margin:5px; border-top: 1px solid grey; }
"),
    shiny::fluidPage(
      shiny::titlePanel(
        "Effectiveness of airport screening for Bundibugyo ebolavirus (BDBV)"
      ),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::sliderInput("dur.flight",
            label = "Travel duration (hours) — DRC/Uganda to UK: ~12 h",
            min = 1, max = 20, step = 1, value = 12
          ),
          shiny::checkboxInput("do_exit",
            label = "Exit screening?",
            value = TRUE
          ),
          shiny::conditionalPanel(
            condition = "input.do_exit == true",
            shiny::sliderInput(
              inputId = "sens.exit",
              label = "Sensitivity of exit screening",
              value = 86, min = 0, max = 100, step = 1, post = " %"
            )
          ),
          shiny::checkboxInput("do_entry",
            label = "Entry screening?",
            value = TRUE
          ),
          shiny::conditionalPanel(
            condition = "input.do_entry == true",
            shiny::sliderInput(
              inputId = "sens.entry",
              label = "Sensitivity of entry screening",
              value = 86, min = 0, max = 100, step = 1, post = " %"
            )
          ),
          shiny::sliderInput(
            inputId = "prop.asy",
            label = "Proportion of cases that are asymptomatic",
            value = 17, min = 0, max = 100, step = 1, post = " %"
          ),
          shiny::checkboxInput("do_growth",
            label = "Adjust for epidemic growth phase?",
            value = FALSE
          ),
          shiny::conditionalPanel(
            condition = "input.do_growth == true",
            shiny::sliderInput(
              inputId = "doubling_time",
              label = shiny::HTML(
                "Epidemic doubling time (days)<br/>",
                "<small style='color:grey'>epiforecasts BDBV 2026, 90% CrI: 13.8\u201322.8 (initial), 4.5\u201343.6 (latest)</small>"
              ),
              min = 4.5, max = 44, step = 0.5, value = 18
            )
          ),
          shiny::div(
            class = "header",
            shiny::selectInput("pathogen",
              label = "Pathogen",
              choices = unique(pathogen_parameters$name),
              selected = pathogen_parameters$name[1]
            ),
            shiny::numericInput("mu_inc",
              "Days from infection to symptom onset (mean)",
              value = 9.0, min = 0.1, max = 30, step = 0.1
            ),
            shiny::numericInput("sigma_inc",
              "Days from infection to symptom onset (SD, days)",
              value = 5.0, min = 0.01, max = 30, step = 0.01
            ),
            shiny::numericInput("mu_inf",
              paste0(
                "Days from symptom onset to severe symptoms e.g.",
                " hospitalisation (mean)"
              ),
              value = 4.03, min = 0.1, max = 30, step = 0.1
            ),
            shiny::numericInput("sigma_inf",
              paste0(
                "Days from symptom onset to severe symptoms e.g.",
                " hospitalisation (SD, days)"
              ),
              value = 3.70, min = 0.01, max = 30, step = 0.01
            ),
            shiny::numericInput("mu_od",
              "Days from symptom onset to death/recovery (mean)",
              value = 11.71, min = 0.1, max = 60, step = 0.1
            ),
            shiny::numericInput("sigma_od",
              "Days from symptom onset to death/recovery (SD, days)",
              value = 6.54, min = 0.01, max = 30, step = 0.01
            )
          ),
          shiny::checkboxInput("uncert",
            label = "Show uncertainty (takes longer)",
            value = FALSE
          )
        ),
        shiny::mainPanel(
          shiny::tabsetPanel(
            type = "tabs",
            shiny::tabPanel(
              title = "Plot",
              (shiny::includeMarkdown(waffle_description)),
              shiny::fluidRow(
                shiny::uiOutput("waffle_plot"),
                shiny::tableOutput("detailed_estimates"),
                align = "center"
              ),
              shiny::includeMarkdown(density_description),
              shiny::fluidRow(shiny::uiOutput("density_plot"))
            ),
            shiny::tabPanel(
              title = "Model",
              shiny::fluidRow(
                shiny::includeMarkdown(assumptions)
              )
            ),
            shiny::tabPanel(
              title = "References",
              shiny::includeMarkdown(references)
            ),
            shiny::tabPanel(
              title = "About",
              shiny::tags$p(
                "This app was originally developed by Quilty et al. (2020) as part of the ",
                "CMMID response to the COVID-19 pandemic. It has been adapted here for ",
                "Bundibugyo ebolavirus (BDBV) using natural-history parameters from the ",
                "2012 Isiro outbreak."
              ),
              shiny::tags$p(shiny::tags$strong("If you use this app, please cite:")),
              shiny::tags$blockquote(
                "Quilty BJ, Clifford S, CMMID nCoV Working Group, Flasche S, Eggo RM. ",
                "Effectiveness of airport screening at detecting travellers infected with ",
                "novel coronavirus (2019-nCoV). ",
                shiny::tags$em("Euro Surveill."),
                " 2020;25(5):pii=2000080. ",
                shiny::tags$a(
                  href = "https://doi.org/10.2807/1560-7917.ES.2020.25.5.2000080",
                  "https://doi.org/10.2807/1560-7917.ES.2020.25.5.2000080"
                )
              )
            )
          )
        )
      )
    )
  )
  server <- function(input, output, session) {
    # Effective inputs: zero out sensitivity when screening is toggled off
    eff_input <- shiny::reactive({
      list(
        dur.flight = input$dur.flight,
        mu_inc = input$mu_inc,
        sigma_inc = input$sigma_inc^2,
        mu_inf = input$mu_inf,
        sigma_inf = input$sigma_inf^2,
        mu_od = input$mu_od,
        sigma_od = input$sigma_od^2,
        sens.exit = if (isTRUE(input$do_exit)) input$sens.exit else 0,
        sens.entry = if (isTRUE(input$do_entry)) input$sens.entry else 0,
        prop.asy = input$prop.asy,
        growth_rate = if (isTRUE(input$do_growth)) log(2) / input$doubling_time else 0
      )
    })

    shiny::observe({
      pathogen_input <- input$pathogen

      shiny::updateNumericInput(session, "prop.asy",
        value = pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$prop.asy
      )
      shiny::updateNumericInput(session, "mu_inc",
        value = pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$mu_inc
      )
      shiny::updateNumericInput(session, "sigma_inc",
        value = round(sqrt(pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$sigma_inc), 2)
      )
      shiny::updateNumericInput(session, "mu_inf",
        value = pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$mu_inf
      )
      shiny::updateNumericInput(session, "sigma_inf",
        value = round(sqrt(pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$sigma_inf), 2)
      )
      shiny::updateNumericInput(session, "mu_od",
        value = pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$mu_od
      )
      shiny::updateNumericInput(session, "sigma_od",
        value = round(sqrt(pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$sigma_od), 2)
      )
    })

    waffle_df <- shiny::reactive({
      travellers <- generate_travellers(eff_input(), i = rep(10000, 1))

      probs <- generate_probabilities(travellers)

      n_exit <- round(probs$prop_symp_at_exit[1] * 10)
      n_sev <- round(probs$prop_sev_at_entry[1] * 10)
      n_entry <- round(probs$prop_symp_at_entry[1] * 10)
      n_missed <- 1000L - (n_exit + n_sev + n_entry)

      waffle_labels <- data.frame(
        desc = factor(
          c(
            rep("detected at exit screening", n_exit),
            rep("detected as severe on flight", n_sev),
            rep("detected at entry screening", n_entry),
            rep("not detected", n_missed)
          ),
          ordered = TRUE,
          levels = c(
            "detected at exit screening",
            "detected as severe on flight",
            "detected at entry screening",
            "not detected"
          )
        )
      ) %>%
        dplyr::mutate(
          desc = factor(
            .data$desc,
            levels = c(
              "detected at exit screening",
              "detected as severe on flight",
              "detected at entry screening",
              "not detected"
            )
          )
        )


      waffle_counts <- dplyr::count(waffle_labels, .data$desc, .drop = FALSE) %>%
        dplyr::mutate(desc_comb = paste0(.data$desc, " (", .data$n, ")"))

      waffle_df <- expand.grid(y = seq_len(25), x = seq_len(40)) %>%
        tibble::as_tibble() %>%
        dplyr::bind_cols(waffle_labels) %>%
        dplyr::mutate(desc_comb = factor(
          .data$desc,
          levels = waffle_counts$desc,
          labels = waffle_counts$desc_comb
        ))
      return(waffle_df)
    })

    nat_hist_periods <- shiny::reactive({
      periods <- data.frame(
        inc_period = time_to_event(1e4, input$mu_inc, input$sigma_inc^2),
        inf_period = time_to_event(1e4, input$mu_inf, input$sigma_inf^2)
      )
      return(periods)
    })

    output$waffleplot <- shiny::renderPlot(expr = {
      waffle <- waffle_df()

      tile_colors <- c(
        "detected at exit screening"   = "#66C2A5",
        "detected as severe on flight" = "#E78AC3",
        "detected at entry screening"  = "#8DA0CB",
        "not detected"                 = "#FC8D62"
      )

      ggplot2::ggplot(
        waffle,
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$desc_comb)
      ) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
        ggplot2::scale_fill_manual(
          values = setNames(
            tile_colors[as.character(levels(waffle$desc))],
            levels(waffle$desc_comb)
          )
        ) +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::labs(
          title = "Out of 1000 infected travellers:",
          fill  = NULL
        ) +
        ggplot2::theme_void(base_size = 13) +
        ggplot2::theme(
          legend.position = "bottom",
          legend.text = ggplot2::element_text(size = 10),
          plot.title = ggplot2::element_text(
            size = 15, face = "bold",
            margin = ggplot2::margin(b = 8)
          ),
          plot.margin = ggplot2::margin(12, 16, 12, 16)
        )
    })

    output$densityplot <- shiny::renderPlot(expr = {
      period_plot_data <- nat_hist_periods()

      period_colors <- c(
        "Infection to onset"  = "#4E84C4",
        "Onset to severe"     = "#E84040",
        "Infection to severe" = "#52854C"
      )

      period_means <- period_plot_data %>%
        dplyr::mutate(severe_period = .data$inc_period + .data$inf_period) %>%
        tidyr::pivot_longer(
          names_to = "Period", cols = dplyr::everything(), values_to = "value"
        ) %>%
        dplyr::mutate(
          Period = factor(
            .data$Period,
            levels = c("inc_period", "inf_period", "severe_period"),
            labels = c("Infection to onset", "Onset to severe", "Infection to severe")
          )
        ) %>%
        dplyr::group_by(.data$Period) %>%
        dplyr::summarise(mean = mean(.data$value), .groups = "drop")

      period_plot_long <- period_plot_data %>%
        dplyr::mutate(severe_period = .data$inc_period + .data$inf_period) %>%
        tidyr::pivot_longer(
          names_to = "Period", cols = dplyr::everything(), values_to = "value"
        ) %>%
        dplyr::mutate(
          Period = factor(
            .data$Period,
            levels = c("inc_period", "inf_period", "severe_period"),
            labels = c("Infection to onset", "Onset to severe", "Infection to severe")
          )
        )

      ggplot2::ggplot(period_plot_long, ggplot2::aes(
        x = .data$value,
        fill = .data$Period,
        colour = .data$Period
      )) +
        ggplot2::geom_density(alpha = 0.25, linewidth = 0.8) +
        ggplot2::geom_vline(
          data = period_means,
          ggplot2::aes(xintercept = .data$mean, colour = .data$Period),
          linetype = "dashed", linewidth = 0.7
        ) +
        ggplot2::scale_fill_manual(values = period_colors) +
        ggplot2::scale_colour_manual(values = period_colors) +
        ggplot2::facet_wrap(~ .data$Period, scales = "free_x", nrow = 1) +
        ggplot2::labs(
          x     = "Time (days)",
          y     = "Density",
          title = "Natural history delay distributions (dashed lines = means)"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          legend.position = "none",
          strip.text = ggplot2::element_text(face = "bold"),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(
            hjust = 0.5, size = 11,
            colour = "grey30"
          )
        )
    }, execOnResize = FALSE)

    output$waffle_plot <- shiny::renderUI({
      shiny::plotOutput("waffleplot")
    })

    output$density_plot <- shiny::renderUI({
      shiny::plotOutput("densityplot", width = "100%", height = "2in")
    })

    output$detailed_estimates <- shiny::renderTable(
      if (input$uncert == TRUE) {
        travellers <- generate_travellers(eff_input(), i = rep(100, 200))
        probs <- generate_probabilities(travellers)

        unconditional_levels <- c(
          "prop_symp_at_exit",
          "prop_sev_at_entry",
          "prop_symp_at_entry",
          "prop_undetected"
        )
        unconditional_labels <- c(
          "Detected at exit screening",
          "Severely ill during flight",
          "Detected at entry screening",
          "Not detected"
        )
        conditional_levels <- c(
          "cond_sev_at_entry",
          "cond_symp_at_entry",
          "cond_undetected"
        )
        conditional_labels <- c(
          "— Severely ill during flight",
          "— Detected at entry screening",
          "— Not detected"
        )

        make_rows <- function(levels, labels) {
          data.frame(
            "Detection outcome" = labels,
            "Estimate (95% CI)" = apply(
              X = probs[, levels] * 10, MARGIN = 2, FUN = make_ci_label
            ),
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        }

        rbind(
          make_rows(unconditional_levels, unconditional_labels),
          data.frame("Detection outcome" = "Of every 1000 who flew:", "Estimate (95% CI)" = "", check.names = FALSE, stringsAsFactors = FALSE),
          make_rows(conditional_levels, conditional_labels)
        )
      } else {
        NULL
      },
      align = "lr"
    )
  }
  shiny::shinyApp(ui = ui, server = server)
}
