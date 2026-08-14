# Figure 2

# Functions

check_packages <- function(packages) {
  missing <- packages[!vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]

  if (length(missing) > 0L) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "), "\n",
      "Install the required packages listed in README.md."
    )
  }
}

ensure_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_required_csv <- function(path, required_columns = character()) {
  if (!file.exists(path)) {
    stop("Required input file not found: ", path)
  }

  data <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "N/A", ".")
  )

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "File ", path, " is missing column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  data
}

assert_unique_id <- function(data, id = "patient_id") {
  duplicated_ids <- unique(data[[id]][duplicated(data[[id]])])
  if (length(duplicated_ids) > 0L) {
    stop(
      "Expected one row per ", id, ". Duplicated values: ",
      paste(duplicated_ids, collapse = ", ")
    )
  }
  invisible(TRUE)
}

as_binary <- function(x, name = deparse(substitute(x))) {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.numeric(x)) {
    invalid <- !is.na(x) & !x %in% c(0, 1)
    if (any(invalid)) {
      stop(name, " must be coded as 0, 1, or NA.")
    }
    return(as.integer(x))
  }

  key <- tolower(trimws(as.character(x)))
  positive <- c("1", "yes", "positive", "pos", "present", "mut", "mutant", "cr", "event")
  negative <- c("0", "no", "negative", "neg", "absent", "wt", "wild-type", "wild type", "non-cr", "censored")

  result <- rep(NA_integer_, length(key))
  result[key %in% positive] <- 1L
  result[key %in% negative] <- 0L
  result[is.na(x)] <- NA_integer_

  invalid <- !is.na(x) & is.na(result)
  if (any(invalid)) {
    stop(
      "Unable to convert ", name, " to 0/1. Unrecognized value(s): ",
      paste(unique(x[invalid]), collapse = ", ")
    )
  }

  result
}

normalise_alteration <- function(x) {
  raw <- trimws(as.character(x))
  key <- toupper(gsub("[ -]+", "_", raw))

  aliases <- c(
    MISSENSE = "Missense",
    MISSENSE_MUTATION = "Missense",
    INFRAME_INDEL = "Inframe_indel",
    IN_FRAME_INDEL = "Inframe_indel",
    IN_FRAME_DEL = "Inframe_indel",
    IN_FRAME_INS = "Inframe_indel",
    CNV = "CNV",
    COPY_NUMBER_VARIATION = "CNV",
    AMPLIFICATION = "Amplification",
    AMP = "Amplification",
    DELETION = "Deletion",
    DEL = "Deletion",
    FRAMESHIFT = "Frameshift",
    FRAME_SHIFT = "Frameshift",
    FRAME_SHIFT_DEL = "Frameshift",
    FRAME_SHIFT_INS = "Frameshift",
    SPLICING = "Splicing",
    SPLICE_SITE = "Splicing",
    STOP_GAINED = "Stop_gained",
    NONSENSE = "Stop_gained",
    NONSENSE_MUTATION = "Stop_gained",
    STOP_LOST = "Stop_lost",
    NONSTOP_MUTATION = "Stop_lost",
    FUSION = "Fusion"
  )

  mapped <- unname(aliases[key])
  mapped[is.na(mapped)] <- raw[is.na(mapped)]
  mapped
}

wald_or_table <- function(model) {
  coefficient_table <- summary(model)$coefficients
  output <- data.frame(
    term = rownames(coefficient_table),
    estimate = coefficient_table[, "Estimate"],
    std.error = coefficient_table[, "Std. Error"],
    statistic = coefficient_table[, "z value"],
    p.value = coefficient_table[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )

  output$OR <- exp(output$estimate)
  output$conf.low <- exp(output$estimate - 1.96 * output$std.error)
  output$conf.high <- exp(output$estimate + 1.96 * output$std.error)
  output
}

loocv_glm_predictions <- function(data, formula) {
  needed <- all.vars(formula)
  analysis_data <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]

  if (nrow(analysis_data) < 5L) {
    stop("Too few complete cases for LOOCV model: ", deparse(formula))
  }

  outcome <- all.vars(formula)[1]
  if (length(unique(analysis_data[[outcome]])) != 2L) {
    stop("Outcome must contain two classes for: ", deparse(formula))
  }

  predictions <- rep(NA_real_, nrow(analysis_data))

  for (i in seq_len(nrow(analysis_data))) {
    training <- analysis_data[-i, , drop = FALSE]
    testing <- analysis_data[i, , drop = FALSE]

    fit <- suppressWarnings(stats::glm(
      formula,
      data = training,
      family = stats::binomial()
    ))

    predictions[i] <- suppressWarnings(stats::predict(
      fit,
      newdata = testing,
      type = "response"
    ))
  }

  if (any(!is.finite(predictions))) {
    stop(
      "LOOCV produced non-finite predictions for ", deparse(formula),
      ". Check separation, sparse categories, and variable coding."
    )
  }

  data.frame(
    patient_id = analysis_data$patient_id,
    observed = analysis_data[[outcome]],
    predicted = predictions,
    stringsAsFactors = FALSE
  )
}

safe_file_label <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

format_p_value <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f"))
  )
}

# Batch two-group tests for binary features.
batch_binary_tests <- function(
    data,
    group,
    variables,
    method = c("fisher", "chisq", "auto"),
    p_adjust_method = "BH") {
  method <- match.arg(method)
  group_values <- data[[group]]
  group_levels <- sort(unique(group_values[!is.na(group_values)]))

  if (length(group_levels) != 2L) {
    stop(group, " must contain exactly two non-missing groups.")
  }

  results <- lapply(variables, function(variable) {
    complete <- stats::complete.cases(group_values, data[[variable]])
    current_group <- factor(group_values[complete], levels = group_levels)
    current_feature <- as_binary(data[[variable]][complete], variable)
    contingency <- table(
      group = current_group,
      feature = factor(current_feature, levels = c(0, 1))
    )

    expected <- suppressWarnings(stats::chisq.test(contingency, correct = FALSE)$expected)
    selected_method <- if (method == "auto") {
      if (all(expected >= 5)) "chisq" else "fisher"
    } else {
      method
    }

    if (selected_method == "fisher") {
      test <- stats::fisher.test(contingency, alternative = "two.sided")
      odds_ratio <- if (length(test$estimate)) unname(test$estimate) else NA_real_
      conf_low <- if (length(test$conf.int)) unname(test$conf.int[1]) else NA_real_
      conf_high <- if (length(test$conf.int)) unname(test$conf.int[2]) else NA_real_
      statistic <- NA_real_
    } else {
      test <- suppressWarnings(stats::chisq.test(contingency, correct = FALSE))
      odds_ratio <- NA_real_
      conf_low <- NA_real_
      conf_high <- NA_real_
      statistic <- unname(test$statistic)
    }

    data.frame(
      variable = variable,
      group_1 = as.character(group_levels[1]),
      group_2 = as.character(group_levels[2]),
      group_1_negative = unname(contingency[1, "0"]),
      group_1_positive = unname(contingency[1, "1"]),
      group_2_negative = unname(contingency[2, "0"]),
      group_2_positive = unname(contingency[2, "1"]),
      method = selected_method,
      odds_ratio = odds_ratio,
      conf.low = conf_low,
      conf.high = conf_high,
      statistic = statistic,
      p.value = test$p.value,
      stringsAsFactors = FALSE
    )
  })

  output <- do.call(rbind, results)
  output$p.adjust <- stats::p.adjust(output$p.value, method = p_adjust_method)
  output[order(output$p.value), , drop = FALSE]
}

# Batch univariable Cox models.
batch_univariable_cox <- function(data, time, event, variables) {
  results <- lapply(variables, function(variable) {
    needed <- c(time, event, variable)
    complete <- data[stats::complete.cases(data[, needed, drop = FALSE]), needed, drop = FALSE]

    if (
      nrow(complete) < 5L ||
      sum(complete[[event]]) < 2L ||
      length(unique(complete[[variable]])) < 2L
    ) {
      return(NULL)
    }

    formula <- stats::as.formula(
      paste0("survival::Surv(", time, ", ", event, ") ~ ", variable)
    )
    fit <- tryCatch(
      survival::coxph(formula, data = complete, ties = "efron"),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }

    fit_summary <- summary(fit)
    coefficients <- fit_summary$coefficients
    confidence <- fit_summary$conf.int

    data.frame(
      feature = variable,
      term = rownames(coefficients),
      n = nrow(complete),
      events = sum(complete[[event]]),
      beta = coefficients[, "coef"],
      HR = confidence[, "exp(coef)"],
      conf.low = confidence[, "lower .95"],
      conf.high = confidence[, "upper .95"],
      p.value_Wald = coefficients[, "Pr(>|z|)"],
      p.value_likelihood = unname(fit_summary$logtest[["pvalue"]]),
      p.value_score = unname(fit_summary$sctest[["pvalue"]]),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })

  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0L) {
    return(data.frame())
  }
  do.call(rbind, results)
}
# Generate the baseline genomic Oncoprint (manuscript Figure 2A).
# This script begins with a processed mutation table and does not perform
# alignment, variant calling, or other upstream bioinformatics processing.

check_packages(c("ComplexHeatmap"))

input_file <- file.path("data", "mutations_long.csv")
gene_order_file <- file.path("data", "gene_order.csv")
output_dir <- file.path("results", "oncoprint")
output_pdf <- file.path(output_dir, "Figure_2A_oncoprint.pdf")
output_frequency <- file.path(output_dir, "gene_frequency.csv")

max_genes <- 30L
show_sample_names <- TRUE
figure_width <- 13
figure_height <- 10

alteration_colours <- c(
  Missense = "#377EB8",
  Inframe_indel = "#4DAF4A",
  CNV = "#984EA3",
  Amplification = "#E41A1C",
  Deletion = "#A65628",
  Frameshift = "#FF7F00",
  Splicing = "#F781BF",
  Stop_gained = "#000000",
  Stop_lost = "#999999",
  Fusion = "#FFD92F"
)

response_colours <- c(`0` = "#E76F51", `1` = "#2A9D8F")
arm_colours <- c(Consolidation = "#4C78A8", Induction = "#F58518")

mutations <- read_required_csv(
  input_file,
  required_columns = c("patient_id", "gene", "alteration", "arm", "CR")
)
mutations$patient_id <- trimws(as.character(mutations$patient_id))
mutations$gene <- toupper(trimws(as.character(mutations$gene)))
mutations$alteration <- normalise_alteration(mutations$alteration)
mutations$arm <- trimws(as.character(mutations$arm))
mutations$CR <- as_binary(mutations$CR, "CR")
mutations <- unique(mutations)

valid <- stats::complete.cases(mutations[, c("patient_id", "gene", "alteration", "arm", "CR")]) &
  mutations$patient_id != "" & mutations$gene != "" & mutations$alteration != ""
mutations <- mutations[valid, , drop = FALSE]

if (nrow(mutations) == 0L) {
  stop("No valid mutation records are available for the Oncoprint.")
}

unknown_alterations <- setdiff(unique(mutations$alteration), names(alteration_colours))
if (length(unknown_alterations) > 0L) {
  stop(
    "Define colours for the following alteration class(es): ",
    paste(unknown_alterations, collapse = ", ")
  )
}

patient_annotation <- unique(mutations[, c("patient_id", "arm", "CR")])
if (any(duplicated(patient_annotation$patient_id))) {
  stop("Each patient must have one arm and one CR value in mutations_long.csv.")
}

patient_annotation$input_order <- seq_len(nrow(patient_annotation))
patient_annotation <- patient_annotation[
  order(-patient_annotation$CR, patient_annotation$arm, patient_annotation$input_order),
  ,
  drop = FALSE
]
sample_order <- patient_annotation$patient_id

gene_frequency <- vapply(
  split(mutations$patient_id, mutations$gene),
  function(x) length(unique(x)),
  integer(1)
)
gene_frequency <- sort(gene_frequency, decreasing = TRUE)

frequency_table <- data.frame(
  gene = names(gene_frequency),
  altered_patients = unname(gene_frequency),
  frequency = unname(gene_frequency) / length(unique(mutations$patient_id)),
  stringsAsFactors = FALSE
)

if (file.exists(gene_order_file)) {
  requested_genes <- read_required_csv(gene_order_file, "gene")$gene
  gene_order <- unique(toupper(trimws(as.character(requested_genes))))
  gene_order <- gene_order[gene_order %in% mutations$gene]
  if (length(gene_order) == 0L) {
    stop("No gene in gene_order.csv was found in mutations_long.csv.")
  }
} else {
  gene_order <- head(names(gene_frequency), max_genes)
}

alterations_present <- names(alteration_colours)[
  names(alteration_colours) %in% unique(mutations$alteration)
]

template <- matrix(
  0L,
  nrow = length(gene_order),
  ncol = length(sample_order),
  dimnames = list(gene_order, sample_order)
)

matrix_list <- setNames(
  lapply(alterations_present, function(x) template),
  alterations_present
)

plot_records <- mutations[
  mutations$gene %in% gene_order & mutations$patient_id %in% sample_order,
  ,
  drop = FALSE
]

for (i in seq_len(nrow(plot_records))) {
  matrix_list[[plot_records$alteration[i]]][
    plot_records$gene[i],
    plot_records$patient_id[i]
  ] <- 1L
}

alter_fun <- function(x, y, w, h, v) {
  grid::grid.rect(
    x, y,
    width = w - grid::unit(0.5, "mm"),
    height = h - grid::unit(0.5, "mm"),
    gp = grid::gpar(fill = "#F2F2F2", col = NA)
  )

  active <- names(which(v))
  if (length(active) == 0L) {
    return(invisible(NULL))
  }

  cell_height <- h - grid::unit(0.5, "mm")
  for (i in seq_along(active)) {
    grid::grid.rect(
      x = x,
      y = y - cell_height / 2 + cell_height * (i - 0.5) / length(active),
      width = w - grid::unit(0.5, "mm"),
      height = cell_height / length(active),
      gp = grid::gpar(fill = alteration_colours[[active[i]]], col = NA)
    )
  }
}

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  CR = factor(patient_annotation$CR, levels = c(1, 0), labels = c("CR", "non-CR")),
  Arm = patient_annotation$arm,
  `Alterations per patient` = ComplexHeatmap::anno_oncoprint_barplot(),
  col = list(
    CR = c(CR = response_colours[["1"]], `non-CR` = response_colours[["0"]]),
    Arm = arm_colours
  )
)

oncoprint <- ComplexHeatmap::oncoPrint(
  matrix_list,
  alter_fun = alter_fun,
  alter_fun_is_vectorized = FALSE,
  col = alteration_colours[alterations_present],
  row_order = gene_order,
  column_order = sample_order,
  top_annotation = top_annotation,
  right_annotation = ComplexHeatmap::rowAnnotation(
    `Altered patients` = ComplexHeatmap::anno_oncoprint_barplot()
  ),
  show_column_names = show_sample_names,
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 9, fontface = "italic"),
  pct_gp = grid::gpar(fontsize = 8),
  remove_empty_rows = TRUE,
  remove_empty_columns = FALSE,
  heatmap_legend_param = list(
    title = "Alteration",
    at = alterations_present,
    labels = gsub("_", " ", alterations_present)
  )
)

ensure_directory(output_dir)
utils::write.csv(frequency_table, output_frequency, row.names = FALSE)

grDevices::pdf(output_pdf, width = figure_width, height = figure_height, family = "Helvetica")
ComplexHeatmap::draw(
  oncoprint,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
grDevices::dev.off()

message("Oncoprint saved to: ", output_pdf)

# Response-association analyses for Figures 2-3 and Tables S5-S7.

check_packages(c("DescTools", "dplyr", "ggplot2", "tidyr"))

clinical_file <- file.path("data", "analysis_data.csv")
mutation_file <- file.path("data", "mutations_long.csv")
output_dir <- file.path("results", "response_associations")
ensure_directory(output_dir)

data <- read_required_csv(
  clinical_file,
  required_columns = c("patient_id", "CR")
)
assert_unique_id(data)
data$CR <- as_binary(data$CR, "CR")

binary_columns <- intersect(
  c("cMRF", "cEMVI", "SMAD4_mut", "FLT3_amp", "ctDNA_T1", "ctDNA_T2", "ctDNA_T3", "ctDNA_T4", "ctDNA_T23", "ctDNA_T234"),
  names(data)
)
for (variable in binary_columns) {
  data[[variable]] <- as_binary(data[[variable]], variable)
}

safe_fisher <- function(x, y) {
  complete <- stats::complete.cases(x, y)
  contingency <- table(x[complete], y[complete])
  if (length(dim(contingency)) != 2L || any(dim(contingency) < 2L)) {
    return(NA_real_)
  }
  stats::fisher.test(contingency, alternative = "two.sided")$p.value
}

categorical_variables <- intersect(
  c("arm", "sex", "stage", binary_columns),
  names(data)
)

fisher_results <- do.call(rbind, lapply(categorical_variables, function(variable) {
  data.frame(
    variable = variable,
    n = sum(stats::complete.cases(data[, c(variable, "CR")])),
    p.value = safe_fisher(data[[variable]], data$CR),
    stringsAsFactors = FALSE
  )
}))
fisher_results$p.adjust_BH <- stats::p.adjust(fisher_results$p.value, method = "BH")

continuous_variables <- intersect(
  c("age", "baseline_CEA", "post_CEA", "TMB"),
  names(data)
)

welch_results <- do.call(rbind, lapply(continuous_variables, function(variable) {
  complete <- data[stats::complete.cases(data[, c(variable, "CR")]), c(variable, "CR")]
  test <- if (length(unique(complete$CR)) == 2L) {
    stats::t.test(complete[[variable]] ~ complete$CR, var.equal = FALSE)
  } else {
    NULL
  }
  data.frame(
    variable = variable,
    n = nrow(complete),
    mean_CR = mean(complete[complete$CR == 1, variable], na.rm = TRUE),
    mean_non_CR = mean(complete[complete$CR == 0, variable], na.rm = TRUE),
    p.value = if (is.null(test)) NA_real_ else test$p.value,
    stringsAsFactors = FALSE
  )
}))
welch_results$p.adjust_BH <- stats::p.adjust(welch_results$p.value, method = "BH")

mrtrg_variables <- intersect(c("mrTRG_T3", "mrTRG_T4"), names(data))
trend_results <- do.call(rbind, lapply(mrtrg_variables, function(variable) {
  complete <- data[stats::complete.cases(data[, c(variable, "CR")]), c(variable, "CR")]
  complete[[variable]] <- as.numeric(complete[[variable]])
  contingency <- table(complete[[variable]], factor(complete$CR, levels = c(0, 1)))
  test <- if (nrow(contingency) >= 2L && ncol(contingency) == 2L) {
    DescTools::CochranArmitageTest(contingency, alternative = "two.sided")
  } else {
    NULL
  }
  data.frame(
    variable = variable,
    n = nrow(complete),
    statistic = if (is.null(test)) NA_real_ else unname(test$statistic),
    p.value = if (is.null(test)) NA_real_ else test$p.value,
    stringsAsFactors = FALSE
  )
}))

# Gene-level mutation enrichment in CR versus non-CR.
mutations <- read_required_csv(
  mutation_file,
  required_columns = c("patient_id", "gene")
)
mutations$gene <- toupper(trimws(as.character(mutations$gene)))

gene_presence <- mutations |>
  dplyr::filter(!is.na(.data$gene), .data$gene != "") |>
  dplyr::distinct(.data$patient_id, .data$gene) |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(
    names_from = "gene",
    values_from = "present",
    values_fill = 0L
  )

gene_data <- data[, c("patient_id", "CR")] |>
  dplyr::left_join(gene_presence, by = "patient_id")
gene_columns <- setdiff(names(gene_data), c("patient_id", "CR"))
gene_data[gene_columns] <- lapply(gene_data[gene_columns], function(x) {
  x[is.na(x)] <- 0L
  x
})

gene_fisher_results <- batch_binary_tests(
  gene_data,
  group = "CR",
  variables = gene_columns,
  method = "fisher",
  p_adjust_method = "BH"
)
names(gene_fisher_results)[names(gene_fisher_results) == "variable"] <- "gene"
names(gene_fisher_results)[names(gene_fisher_results) == "p.adjust"] <- "p.adjust_BH"

# Univariable logistic regressions. Character variables are handled as factors by glm().
logistic_candidates <- intersect(
  c(
    "arm", "age", "sex", "stage", "cMRF", "cEMVI", "baseline_CEA",
    "post_CEA", "mrTRG_T3", "mrTRG_T4", "TMB", "SMAD4_mut", "FLT3_amp",
    "ctDNA_T1", "ctDNA_T2", "ctDNA_T3", "ctDNA_T4", "ctDNA_T23", "ctDNA_T234"
  ),
  names(data)
)

univariable_logistic <- lapply(logistic_candidates, function(variable) {
  complete <- data[stats::complete.cases(data[, c("CR", variable)]), c("CR", variable)]
  if (nrow(complete) < 5L || length(unique(complete$CR)) < 2L || length(unique(complete[[variable]])) < 2L) {
    return(NULL)
  }

  formula <- stats::as.formula(paste("CR ~", variable))
  fit <- suppressWarnings(stats::glm(formula, data = complete, family = stats::binomial()))
  table <- wald_or_table(fit)
  table <- table[table$term != "(Intercept)", , drop = FALSE]
  table$variable <- variable
  table$n <- nrow(complete)
  table[, c("variable", "term", "n", "OR", "conf.low", "conf.high", "p.value")]
})
univariable_logistic <- do.call(rbind, univariable_logistic[!vapply(univariable_logistic, is.null, logical(1))])

utils::write.csv(fisher_results, file.path(output_dir, "categorical_fisher_tests.csv"), row.names = FALSE)
utils::write.csv(welch_results, file.path(output_dir, "continuous_welch_tests.csv"), row.names = FALSE)
utils::write.csv(trend_results, file.path(output_dir, "mrTRG_cochran_armitage_tests.csv"), row.names = FALSE)
utils::write.csv(gene_fisher_results, file.path(output_dir, "gene_fisher_tests_BH.csv"), row.names = FALSE)
utils::write.csv(univariable_logistic, file.path(output_dir, "univariable_logistic_regression.csv"), row.names = FALSE)

# Manuscript-style response plots for selected features.
plot_binary_response <- function(variable, label) {
  plot_data <- data[stats::complete.cases(data[, c(variable, "CR")]), c(variable, "CR")]
  plot_data[[variable]] <- factor(plot_data[[variable]], levels = c(0, 1), labels = c("Negative/WT", "Positive/Mut"))

  summary_data <- aggregate(CR ~ group, data = transform(plot_data, group = plot_data[[variable]]), FUN = mean)

  ggplot2::ggplot(summary_data, ggplot2::aes(x = .data$group, y = .data$CR)) +
    ggplot2::geom_col(width = 0.65, fill = "#4C78A8") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = label, y = "Complete response rate") +
    ggplot2::theme_classic(base_size = 11)
}

check_packages(c("scales"))
selected_binary_plots <- intersect(c("FLT3_amp", "SMAD4_mut", "ctDNA_T234"), names(data))
for (variable in selected_binary_plots) {
  plot <- plot_binary_response(variable, variable)
  ggplot2::ggsave(
    file.path(output_dir, paste0("CR_rate_", variable, ".pdf")),
    plot = plot,
    width = 4,
    height = 4
  )
}

if ("arm" %in% names(data)) {
  arm_data <- data[stats::complete.cases(data[, c("arm", "CR")]), c("arm", "CR")]
  arm_summary <- aggregate(CR ~ arm, data = arm_data, FUN = mean)
  arm_plot <- ggplot2::ggplot(arm_summary, ggplot2::aes(x = .data$arm, y = .data$CR)) +
    ggplot2::geom_col(width = 0.65, fill = "#4C78A8") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = "Treatment arm", y = "Complete response rate") +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(file.path(output_dir, "CR_rate_arm.pdf"), arm_plot, width = 4.5, height = 4)
}

if ("TMB" %in% names(data)) {
  tmb_data <- data[stats::complete.cases(data[, c("TMB", "CR")]), ]
  tmb_data$Response <- factor(tmb_data$CR, levels = c(1, 0), labels = c("CR", "non-CR"))
  tmb_plot <- ggplot2::ggplot(tmb_data, ggplot2::aes(x = .data$Response, y = .data$TMB)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
    ggplot2::labs(x = NULL, y = "Tumor mutational burden (mut/Mb)") +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(file.path(output_dir, "TMB_by_response.pdf"), tmb_plot, width = 4, height = 4)
}

message("Response-association results saved to: ", output_dir)

