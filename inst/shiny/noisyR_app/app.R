# noisyR interactive app: count matrix approach
# Launch with noisyr::run_noisyr_app()
library(shiny)
library(bslib)
library(plotly)

# functions come from the installed package, or from an environment supplied
# through shinyOptions(noisyr_functions = ...) when developing from source
nr <- getShinyOption("noisyr_functions", asNamespace("noisyr"))
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

simulate_counts <- function(n.genes = 8000, n.samples = 6, seed = 42) {
  set.seed(seed)
  mu <- stats::rlnorm(n.genes, meanlog = 3, sdlog = 2)
  m <- vapply(seq_len(n.samples), function(j) {
    stats::rnbinom(n.genes, mu = mu * stats::runif(1, 0.7, 1.3), size = 2)
  }, numeric(n.genes))
  dimnames(m) <- list(paste0("gene", seq_len(n.genes)), paste0("sample", seq_len(n.samples)))
  m
}

read_counts <- function(path, name) {
  sep <- if (grepl("\\.(tsv|txt)(\\.gz)?$", name, ignore.case = TRUE)) "\t" else ","
  df <- utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.numeric(df[[1]])) {
    rownames(df) <- make.unique(as.character(df[[1]]))
    df <- df[, -1, drop = FALSE]
  }
  nr$cast_matrix_to_numeric(as.matrix(df))
}

methods.threshold <- nr$get_methods_calculate_noise_threshold()
methods.similarity <- nr$get_methods_correlation_distance()

step_card <- function(title, ...) card(card_header(title), ...)

ui <- page_navbar(
  title = "noisyR",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel(
    "1 · Data",
    layout_sidebar(
      sidebar = sidebar(
        fileInput("counts_file", "Count matrix (CSV/TSV, genes x samples)",
                  accept = c(".csv", ".tsv", ".txt", ".gz")),
        helpText("First column: gene identifiers. Remaining columns: raw or normalised counts."),
        actionButton("use_demo", "Use simulated demo data", class = "btn-outline-primary")
      ),
      layout_columns(
        value_box("Genes", textOutput("n_genes")),
        value_box("Samples", textOutput("n_samples")),
        value_box("Zero entries", textOutput("pct_zero"))
      ),
      layout_columns(
        step_card("Library sizes", plotlyOutput("libsize_plot", height = 300)),
        step_card("Expression distribution (log2 counts + 1)", plotlyOutput("density_plot", height = 300))
      )
    )
  ),
  nav_panel(
    "2 · Similarity",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("similarity_measure", "Similarity measure", methods.similarity,
                    selected = "correlation_pearson"),
        numericInput("window_length", "Genes per window", value = NA, min = 2),
        helpText("Leave empty for the default (10% of genes)."),
        sliderInput("step_fraction", "Step (fraction of window)", 0.01, 0.5, 0.05, 0.01),
        numericInput("n_cores", "Cores", value = 1, min = 1, max = parallel::detectCores()),
        actionButton("run_similarity", "Compute similarity", class = "btn-primary"),
        hr(),
        actionButton("run_optimise", "Suggest window length (slower)", class = "btn-outline-secondary")
      ),
      step_card("Similarity of each sliding window to the other samples",
                plotlyOutput("similarity_plot", height = 420)),
      step_card("Window length optimisation (Jensen-Shannon entropy)",
                textOutput("optimise_result"))
    )
  ),
  nav_panel(
    "3 · Threshold",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("method_chosen", "Threshold method", methods.threshold, selected = "Boxplot-IQR"),
        sliderInput("similarity_threshold", "Similarity threshold", 0, 1, 0.25, 0.01),
        numericInput("binsize", "Boxplot bin size (log2)", 0.1, min = 0.01, step = 0.05),
        helpText("The slider updates instantly: curves are computed once per method."),
        hr(),
        actionButton("run_stats", "Compare all methods", class = "btn-outline-primary"),
        helpText("Click a heatmap cell to select that method and threshold.")
      ),
      layout_columns(
        col_widths = c(8, 4),
        step_card("Threshold read-off per sample", plotlyOutput("threshold_plot", height = 420)),
        step_card("Noise thresholds", tableOutput("threshold_table"),
                  textOutput("threshold_cv"))
      ),
      step_card("Coefficient of variation of thresholds across samples (lower = more consistent)",
                plotlyOutput("stats_heatmap", height = 360))
    )
  ),
  nav_panel(
    "4 · Denoise",
    layout_sidebar(
      sidebar = sidebar(
        checkboxInput("add_threshold", "Add threshold to all entries (else floor at threshold)", TRUE),
        checkboxInput("average_threshold", "Use the average threshold across samples", TRUE),
        checkboxInput("remove_noisy", "Remove genes below the threshold in all samples", TRUE),
        downloadButton("download_matrix", "Denoised matrix (CSV)"),
        downloadButton("download_script", "Reproducible R script")
      ),
      layout_columns(
        value_box("Genes kept", textOutput("genes_kept")),
        value_box("Genes removed", textOutput("genes_removed")),
        value_box("Threshold used", textOutput("threshold_used"))
      ),
      step_card("Expression before and after denoising (log2 counts + 1)",
                plotlyOutput("denoise_plot", height = 380))
    )
  )
)

server <- function(input, output, session) {
  counts <- reactiveVal(NULL)
  observeEvent(input$counts_file, {
    m <- tryCatch(read_counts(input$counts_file$datapath, input$counts_file$name),
                  error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    counts(m)
  })
  observeEvent(input$use_demo, counts(simulate_counts()))

  output$n_genes <- renderText(if (is.null(counts())) "-" else format(nrow(counts()), big.mark = ","))
  output$n_samples <- renderText(if (is.null(counts())) "-" else ncol(counts()))
  output$pct_zero <- renderText(if (is.null(counts())) "-" else sprintf("%.1f%%", 100 * mean(counts() == 0)))

  output$libsize_plot <- renderPlotly({
    req(counts())
    plot_ly(x = colnames(counts()), y = colSums(counts()), type = "bar") |>
      layout(xaxis = list(title = ""), yaxis = list(title = "total counts"))
  })
  output$density_plot <- renderPlotly({
    req(counts())
    p <- plot_ly()
    for (j in seq_len(ncol(counts()))) {
      d <- stats::density(log2(counts()[, j] + 1))
      p <- add_lines(p, x = d$x, y = d$y, name = colnames(counts())[j])
    }
    layout(p, xaxis = list(title = "log2(count + 1)"), yaxis = list(title = "density"))
  })

  # Step 2: similarity -------------------------------------------------------
  summary_params <- reactiveVal(NULL)
  expression_summary <- eventReactive(input$run_similarity, {
    req(counts())
    window <- if (is.na(input$window_length)) NULL else input$window_length
    summary_params(list(similarity.measure = input$similarity_measure,
                        n.elements.per.window = window,
                        n.step.fraction = input$step_fraction))
    withProgress(message = "Computing window similarities...", {
      quiet(nr$calculate_expression_similarity_counts(
        counts(), similarity.measure = input$similarity_measure,
        n.elements.per.window = window, n.step.fraction = input$step_fraction,
        n.cores = input$n_cores))
    })
  })

  sample_names <- reactive({
    if (is.null(colnames(counts()))) paste0("sample", seq_len(ncol(counts()))) else colnames(counts())
  })

  output$similarity_plot <- renderPlotly({
    es <- expression_summary()
    p <- plot_ly()
    for (j in seq_len(ncol(es$expression.levels))) {
      p <- add_lines(p, x = log2(es$expression.levels[, j] + 1), y = es$expression.levels.similarity[, j],
                     name = sample_names()[j])
    }
    layout(p, xaxis = list(title = "window mean expression, log2(count + 1)"),
           yaxis = list(title = input$similarity_measure))
  })

  optimised <- eventReactive(input$run_optimise, {
    req(counts())
    withProgress(message = "Optimising window length...", {
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off())
      quiet(nr$optimise_window_length(counts(), similarity.measure = input$similarity_measure,
                                      iteration.number = 20, n.cores = input$n_cores))
    })
  })
  observeEvent(optimised(), updateNumericInput(session, "window_length", value = optimised()))
  output$optimise_result <- renderText({
    w <- optimised()
    if (is.null(w)) "No stable window length found; try fewer iterations or a finer search."
    else paste("Suggested window length:", w, "genes (filled in on the left, recompute to use it).")
  })

  # Step 3: thresholds -------------------------------------------------------
  # the expensive part depends on the method only, so the slider is instant
  curves <- reactive({
    quiet(nr$calculate_noise_threshold_curves(expression_summary(), method.chosen = input$method_chosen,
                                              binsize = input$binsize))
  })
  noise_thresholds <- reactive(quiet(nr$noise_thresholds_from_curves(curves(), input$similarity_threshold)))

  output$threshold_plot <- renderPlotly({
    es <- expression_summary()
    thr <- rep_len(noise_thresholds(), ncol(es$expression.levels))
    p <- plot_ly()
    for (j in seq_len(ncol(es$expression.levels))) {
      p <- add_lines(p, x = log2(es$expression.levels[, j] + 1), y = es$expression.levels.similarity[, j],
                     name = sample_names()[j], legendgroup = j, line = list(width = 1))
    }
    shapes <- c(
      list(list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                y0 = input$similarity_threshold, y1 = input$similarity_threshold,
                line = list(dash = "dash", color = "black"))),
      lapply(log2(thr + 1), function(x) list(type = "line", yref = "paper", y0 = 0, y1 = 1,
                                             x0 = x, x1 = x, line = list(color = "grey", width = 1)))
    )
    layout(p, shapes = shapes, xaxis = list(title = "log2(count + 1)"),
           yaxis = list(title = "similarity"))
  })
  output$threshold_table <- renderTable({
    thr <- noise_thresholds()
    data.frame(sample = if (length(thr) == 1) "all (fixed)" else sample_names(), threshold = thr)
  }, digits = 1)
  output$threshold_cv <- renderText({
    thr <- noise_thresholds()
    if (length(thr) > 1) sprintf("CV across samples: %.3f", stats::sd(thr) / mean(thr)) else ""
  })

  stats_table <- eventReactive(input$run_stats, {
    withProgress(message = "Evaluating all methods...", {
      quiet(nr$calculate_noise_threshold_method_statistics(
        expression_summary(), similarity.threshold.sequence = seq(0.05, 0.6, 0.05),
        method.chosen.sequence = methods.threshold[!grepl("^Density", methods.threshold)],
        binsize = input$binsize))
    })
  })
  output$stats_heatmap <- renderPlotly({
    st <- stats_table()
    st$method.full <- paste(st$approach, st$method, sep = "-")
    plot_ly(st, x = ~similarity.threshold, y = ~method.full, z = ~noise.threshold.coefficient.of.variation,
            type = "heatmap", colorscale = "Viridis", source = "stats",
            text = ~sprintf("mean threshold %.1f", noise.threshold.mean), hoverinfo = "x+y+z+text") |>
      layout(xaxis = list(title = "similarity threshold"), yaxis = list(title = "")) |>
      event_register("plotly_click")
  })
  observeEvent(event_data("plotly_click", source = "stats"), {
    click <- event_data("plotly_click", source = "stats")
    updateSelectInput(session, "method_chosen", selected = click$y)
    updateSliderInput(session, "similarity_threshold", value = click$x)
  })

  # Step 4: denoising --------------------------------------------------------
  denoised <- reactive({
    quiet(nr$remove_noise_from_matrix(counts(), noise_thresholds(), add.threshold = input$add_threshold,
                                      average.threshold = input$average_threshold,
                                      remove.noisy.features = input$remove_noisy))
  })
  output$genes_kept <- renderText(format(nrow(denoised()), big.mark = ","))
  output$genes_removed <- renderText(format(nrow(counts()) - nrow(denoised()), big.mark = ","))
  output$threshold_used <- renderText({
    thr <- noise_thresholds()
    if (input$average_threshold || length(thr) == 1) round(mean(thr)) else paste(round(range(thr)), collapse = "-")
  })
  output$denoise_plot <- renderPlotly({
    before <- stats::density(log2(as.vector(counts()) + 1))
    after <- stats::density(log2(as.vector(denoised()) + 1))
    plot_ly() |>
      add_lines(x = before$x, y = before$y, name = "before") |>
      add_lines(x = after$x, y = after$y, name = "after") |>
      layout(xaxis = list(title = "log2(count + 1)"), yaxis = list(title = "density"))
  })
  output$download_matrix <- downloadHandler(
    filename = function() "noisyR_denoised_counts.csv",
    content = function(file) utils::write.csv(denoised(), file)
  )
  output$download_script <- downloadHandler(
    filename = function() "noisyR_analysis.R",
    content = function(file) {
      p <- summary_params()
      writeLines(c(
        "library(noisyr)",
        "expression.matrix <- cast_matrix_to_numeric(as.matrix(read.csv(\"counts.csv\", row.names = 1)))",
        sprintf("expression.summary <- calculate_expression_similarity_counts(expression.matrix, similarity.measure = %s, n.elements.per.window = %s, n.step.fraction = %s)",
                deparse(p$similarity.measure), deparse(p$n.elements.per.window), p$n.step.fraction),
        sprintf("noise.thresholds <- calculate_noise_threshold(expression.summary, similarity.threshold = %s, method.chosen = %s, binsize = %s)",
                input$similarity_threshold, deparse(input$method_chosen), input$binsize),
        sprintf("expression.matrix.denoised <- remove_noise_from_matrix(expression.matrix, noise.thresholds, add.threshold = %s, average.threshold = %s, remove.noisy.features = %s)",
                input$add_threshold, input$average_threshold, input$remove_noisy)
      ), file)
    }
  )
}

shinyApp(ui, server)
