# Figure 3

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
# Multivariable logistic regression, LOOCV ROC curves, and bootstrap AUC tests.
# Corresponds to manuscript Figure 3D-I.

check_packages(c("ggplot2", "pROC"))

input_file <- file.path("data", "analysis_data.csv")
output_dir <- file.path("results", "logistic_models")
ensure_directory(output_dir)

set.seed(20260813)

data <- read_required_csv(
  input_file,
  required_columns = c("patient_id", "CR")
)
assert_unique_id(data)
data$CR <- as_binary(data$CR, "CR")

binary_columns <- intersect(
  c("SMAD4_mut", "ctDNA_T23", "ctDNA_T234"),
  names(data)
)
for (variable in binary_columns) {
  data[[variable]] <- as_binary(data[[variable]], variable)
}

# Each pair is evaluated in the same complete-case population so the ROC curves
# and paired bootstrap comparison are directly comparable.
model_comparisons <- list(
  T4_tissue = list(
    figure = "Figure_3E",
    base = c("mrTRG_T4"),
    integrated = c("mrTRG_T4", "SMAD4_mut", "TMB")
  ),
  T4_tissue_ctDNA = list(
    figure = "Figure_3G",
    base = c("mrTRG_T4"),
    integrated = c("mrTRG_T4", "SMAD4_mut", "TMB", "ctDNA_T234")
  ),
  T3_tissue = list(
    figure = "Figure_3H",
    base = c("mrTRG_T3"),
    integrated = c("mrTRG_T3", "SMAD4_mut", "TMB")
  ),
  T3_tissue_ctDNA = list(
    figure = "Figure_3I",
    base = c("mrTRG_T3"),
    integrated = c("mrTRG_T3", "SMAD4_mut", "TMB", "ctDNA_T23")
  )
)

missing_model_columns <- unique(unlist(lapply(model_comparisons, function(x) {
  setdiff(c(x$base, x$integrated), names(data))
})))
if (length(missing_model_columns) > 0L) {
  stop(
    "analysis_data.csv is missing model column(s): ",
    paste(missing_model_columns, collapse = ", ")
  )
}

auc_results <- list()
prediction_results <- list()
coefficient_results <- list()

for (comparison_name in names(model_comparisons)) {
  specification <- model_comparisons[[comparison_name]]
  complete_columns <- unique(c("patient_id", "CR", specification$integrated))
  comparison_data <- data[
    stats::complete.cases(data[, complete_columns, drop = FALSE]),
    complete_columns,
    drop = FALSE
  ]

  if (nrow(comparison_data) < 10L || length(unique(comparison_data$CR)) < 2L) {
    warning("Skipping ", comparison_name, ": insufficient complete cases or outcome classes.")
    next
  }

  base_formula <- stats::as.formula(
    paste("CR ~", paste(specification$base, collapse = " + "))
  )
  integrated_formula <- stats::as.formula(
    paste("CR ~", paste(specification$integrated, collapse = " + "))
  )

  base_predictions <- loocv_glm_predictions(comparison_data, base_formula)
  integrated_predictions <- loocv_glm_predictions(comparison_data, integrated_formula)

  if (!identical(base_predictions$patient_id, integrated_predictions$patient_id)) {
    stop("LOOCV patient ordering differs within comparison: ", comparison_name)
  }

  roc_base <- pROC::roc(
    response = base_predictions$observed,
    predictor = base_predictions$predicted,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  roc_integrated <- pROC::roc(
    response = integrated_predictions$observed,
    predictor = integrated_predictions$predicted,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  bootstrap_test <- pROC::roc.test(
    roc_base,
    roc_integrated,
    method = "bootstrap",
    boot.n = 1000,
    paired = TRUE,
    boot.stratified = TRUE
  )

  auc_base <- as.numeric(pROC::auc(roc_base))
  auc_integrated <- as.numeric(pROC::auc(roc_integrated))

  auc_results[[comparison_name]] <- data.frame(
    comparison = comparison_name,
    figure = specification$figure,
    n = nrow(comparison_data),
    events_CR = sum(comparison_data$CR == 1),
    AUC_base = auc_base,
    AUC_integrated = auc_integrated,
    AUC_difference = auc_integrated - auc_base,
    bootstrap_p_value = bootstrap_test$p.value,
    stringsAsFactors = FALSE
  )

  base_predictions$comparison <- comparison_name
  base_predictions$model <- "mrTRG only"
  integrated_predictions$comparison <- comparison_name
  integrated_predictions$model <- "Integrated"
  prediction_results[[paste0(comparison_name, "_base")]] <- base_predictions
  prediction_results[[paste0(comparison_name, "_integrated")]] <- integrated_predictions

  base_fit <- suppressWarnings(stats::glm(
    base_formula,
    data = comparison_data,
    family = stats::binomial()
  ))
  integrated_fit <- suppressWarnings(stats::glm(
    integrated_formula,
    data = comparison_data,
    family = stats::binomial()
  ))

  base_coefficients <- wald_or_table(base_fit)
  base_coefficients$comparison <- comparison_name
  base_coefficients$model <- "mrTRG only"
  integrated_coefficients <- wald_or_table(integrated_fit)
  integrated_coefficients$comparison <- comparison_name
  integrated_coefficients$model <- "Integrated"
  coefficient_results[[paste0(comparison_name, "_base")]] <- base_coefficients
  coefficient_results[[paste0(comparison_name, "_integrated")]] <- integrated_coefficients

  roc_plot <- pROC::ggroc(
    list(`mrTRG only` = roc_base, Integrated = roc_integrated),
    legacy.axes = TRUE,
    size = 1
  ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey60") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = specification$figure,
      x = "1 - Specificity",
      y = "Sensitivity",
      colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 11)

  ggplot2::ggsave(
    file.path(output_dir, paste0(specification$figure, "_LOOCV_ROC.pdf")),
    plot = roc_plot,
    width = 5,
    height = 5
  )
}

if (length(auc_results) > 0L) {
  utils::write.csv(
    do.call(rbind, auc_results),
    file.path(output_dir, "LOOCV_AUC_bootstrap_comparisons.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    do.call(rbind, prediction_results),
    file.path(output_dir, "LOOCV_patient_predictions.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    do.call(rbind, coefficient_results),
    file.path(output_dir, "multivariable_logistic_coefficients.csv"),
    row.names = FALSE
  )
}

message("Logistic-model results saved to: ", output_dir)

