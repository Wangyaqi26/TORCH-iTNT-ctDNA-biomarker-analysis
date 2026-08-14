# Supplementary analyses

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
# Summarize baseline and treatment characteristics by TORCH treatment arm.

input_file <- file.path("data", "analysis_data.csv")
output_dir <- file.path("results", "cohort_summary")
ensure_directory(output_dir)

data <- read_required_csv(
  input_file,
  required_columns = c("patient_id", "arm", "CR")
)
assert_unique_id(data)
data$CR <- as_binary(data$CR, "CR")

continuous_variables <- intersect(
  c("age", "baseline_CEA", "post_CEA", "TMB"),
  names(data)
)

categorical_variables <- intersect(
  c(
    "sex", "stage", "cMRF", "cEMVI", "mrTRG_T3", "mrTRG_T4",
    "treatment_choice", "CR"
  ),
  names(data)
)

arm_levels <- unique(stats::na.omit(data$arm))
if (length(arm_levels) != 2L) {
  stop("Expected exactly two treatment-arm levels; found: ", paste(arm_levels, collapse = ", "))
}

summarize_continuous <- function(variable) {
  values <- lapply(arm_levels, function(current_arm) {
    x <- data[data$arm == current_arm, variable]
    x <- x[is.finite(x)]
    data.frame(
      variable = variable,
      arm = current_arm,
      n = length(x),
      median = if (length(x) > 0L) stats::median(x) else NA_real_,
      minimum = if (length(x) > 0L) min(x) else NA_real_,
      maximum = if (length(x) > 0L) max(x) else NA_real_,
      mean = if (length(x) > 0L) mean(x) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  output <- do.call(rbind, values)
  complete <- data[stats::complete.cases(data[, c("arm", variable)]), c("arm", variable)]
  p_value <- if (length(unique(complete$arm)) == 2L) {
    stats::t.test(complete[[variable]] ~ complete$arm, var.equal = FALSE)$p.value
  } else {
    NA_real_
  }
  output$p.value_Welch <- p_value
  output
}

continuous_summary <- if (length(continuous_variables) > 0L) {
  do.call(rbind, lapply(continuous_variables, summarize_continuous))
} else {
  data.frame()
}

summarize_categorical <- function(variable) {
  complete <- data[stats::complete.cases(data[, c("arm", variable)]), c("arm", variable)]
  contingency <- table(complete[[variable]], complete$arm)

  p_value <- if (all(dim(contingency) >= 2L)) {
    stats::fisher.test(contingency)$p.value
  } else {
    NA_real_
  }

  output <- as.data.frame(contingency, stringsAsFactors = FALSE)
  names(output) <- c("level", "arm", "n")
  output$variable <- variable
  denominators <- table(complete$arm)
  output$percent <- 100 * output$n / as.numeric(denominators[as.character(output$arm)])
  output$p.value_Fisher <- p_value
  output[, c("variable", "level", "arm", "n", "percent", "p.value_Fisher")]
}

categorical_summary <- if (length(categorical_variables) > 0L) {
  do.call(rbind, lapply(categorical_variables, summarize_categorical))
} else {
  data.frame()
}

utils::write.csv(
  continuous_summary,
  file.path(output_dir, "continuous_variables_by_arm.csv"),
  row.names = FALSE
)
utils::write.csv(
  categorical_summary,
  file.path(output_dir, "categorical_variables_by_arm.csv"),
  row.names = FALSE
)

message("Cohort summaries saved to: ", output_dir)

# Optional representativeness analysis: biomarker sub-cohort versus the overall
# evaluable TORCH population (Additional Table S1).
population_file <- file.path("data", "torch_population.csv")
if (file.exists(population_file)) {
  population <- read_required_csv(
    population_file,
    required_columns = c("patient_id", "included_biomarker_cohort")
  )
  assert_unique_id(population)
  population$included_biomarker_cohort <- as_binary(
    population$included_biomarker_cohort,
    "included_biomarker_cohort"
  )

  population_continuous <- intersect(c("age"), names(population))
  population_categorical <- intersect(
    c("arm", "sex", "stage", "cMRF", "cEMVI"),
    names(population)
  )

  population_continuous_results <- lapply(population_continuous, function(variable) {
    complete <- population[stats::complete.cases(population[, c("included_biomarker_cohort", variable)]), ]
    test <- stats::t.test(
      complete[[variable]] ~ complete$included_biomarker_cohort,
      var.equal = FALSE
    )
    data.frame(
      variable = variable,
      n_included = sum(complete$included_biomarker_cohort == 1L),
      n_not_included = sum(complete$included_biomarker_cohort == 0L),
      median_included = stats::median(complete[complete$included_biomarker_cohort == 1L, variable]),
      median_not_included = stats::median(complete[complete$included_biomarker_cohort == 0L, variable]),
      p.value_Welch = test$p.value,
      stringsAsFactors = FALSE
    )
  })

  population_categorical_results <- lapply(population_categorical, function(variable) {
    complete <- population[stats::complete.cases(population[, c("included_biomarker_cohort", variable)]), ]
    contingency <- table(complete[[variable]], complete$included_biomarker_cohort)
    data.frame(
      variable = variable,
      n = nrow(complete),
      p.value_Fisher = if (all(dim(contingency) >= 2L)) stats::fisher.test(contingency)$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  if (length(population_continuous_results) > 0L) {
    utils::write.csv(
      do.call(rbind, population_continuous_results),
      file.path(output_dir, "Table_S1_continuous_representativeness.csv"),
      row.names = FALSE
    )
  }
  if (length(population_categorical_results) > 0L) {
    utils::write.csv(
      do.call(rbind, population_categorical_results),
      file.path(output_dir, "Table_S1_categorical_representativeness.csv"),
      row.names = FALSE
    )
  }
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
# Kaplan-Meier and Cox proportional-hazards analyses.
# Corresponds to Figures 4B-E, Figure 5, Figures S4-S6, and Tables S8-S11.

check_packages(c("ggplot2", "survival", "survminer"))

input_file <- file.path("data", "analysis_data.csv")
output_root <- file.path("results", "survival")
ensure_directory(output_root)

data <- read_required_csv(input_file, required_columns = c("patient_id"))
assert_unique_id(data)

binary_columns <- intersect(
  c(
    "CR", "cMRF", "cEMVI", "SMAD4_mut", "ctDNA_T1", "ctDNA_T2",
    "ctDNA_T3", "ctDNA_T4", "ctDNA_T23", "ctDNA_T234",
    "ctDNA_clearance_T2", "ctDNA_clearance_T3", "ctDNA_clearance_T4"
  ),
  names(data)
)
for (variable in binary_columns) {
  data[[variable]] <- as_binary(data[[variable]], variable)
}

endpoints <- list(
  DFS = c(time = "DFS_time", event = "DFS_event"),
  LRFS = c(time = "LRFS_time", event = "LRFS_event"),
  DMFS = c(time = "DMFS_time", event = "DMFS_event"),
  OS = c(time = "OS_time", event = "OS_event")
)

km_group_variables <- intersect(
  c(
    "CR", "arm", "SMAD4_mut", "ctDNA_T1", "ctDNA_T2", "ctDNA_T3",
    "ctDNA_T4", "ctDNA_T23", "ctDNA_T234", "ctDNA_clearance_T2",
    "ctDNA_clearance_T3", "ctDNA_clearance_T4"
  ),
  names(data)
)

cox_candidate_variables <- intersect(
  c(
    "age", "sex", "arm", "stage", "cMRF", "cEMVI", "baseline_CEA",
    "post_CEA", "mrTRG_T3", "mrTRG_T4", "TMB", "SMAD4_mut",
    "ctDNA_T1", "ctDNA_T2", "ctDNA_T3", "ctDNA_T4", "ctDNA_T23", "ctDNA_T234"
  ),
  names(data)
)

cox_table <- function(fit, feature = NA_character_) {
  summary_fit <- summary(fit)
  coefficients <- summary_fit$coefficients
  confidence <- summary_fit$conf.int

  data.frame(
    feature = feature,
    term = rownames(coefficients),
    HR = confidence[, "exp(coef)"],
    conf.low = confidence[, "lower .95"],
    conf.high = confidence[, "upper .95"],
    p.value = coefficients[, "Pr(>|z|)"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

for (endpoint_name in names(endpoints)) {
  time_column <- endpoints[[endpoint_name]][["time"]]
  event_column <- endpoints[[endpoint_name]][["event"]]

  if (!all(c(time_column, event_column) %in% names(data))) {
    warning("Skipping ", endpoint_name, ": missing time or event column.")
    next
  }

  data[[event_column]] <- as_binary(data[[event_column]], event_column)
  endpoint_dir <- file.path(output_root, endpoint_name)
  ensure_directory(endpoint_dir)

  # Kaplan-Meier analyses.
  logrank_results <- list()
  for (group_variable in km_group_variables) {
    complete <- data[stats::complete.cases(data[, c(time_column, event_column, group_variable)]), , drop = FALSE]
    if (nrow(complete) < 5L || length(unique(complete[[group_variable]])) < 2L) {
      next
    }

    formula <- stats::as.formula(
      paste0("survival::Surv(", time_column, ", ", event_column, ") ~ ", group_variable)
    )
    fit <- survival::survfit(formula, data = complete)
    difference <- survival::survdiff(formula, data = complete)
    p_value <- 1 - stats::pchisq(difference$chisq, df = length(difference$n) - 1L)

    logrank_results[[group_variable]] <- data.frame(
      endpoint = endpoint_name,
      group = group_variable,
      n = nrow(complete),
      events = sum(complete[[event_column]]),
      p.value_logrank = p_value,
      stringsAsFactors = FALSE
    )

    km_plot <- survminer::ggsurvplot(
      fit,
      data = complete,
      risk.table = TRUE,
      pval = TRUE,
      conf.int = FALSE,
      xlab = "Time (months)",
      ylab = paste0(endpoint_name, " probability"),
      legend.title = group_variable,
      ggtheme = ggplot2::theme_classic(base_size = 11)
    )

    output_pdf <- file.path(
      endpoint_dir,
      paste0("KM_", endpoint_name, "_by_", safe_file_label(group_variable), ".pdf")
    )
    grDevices::pdf(output_pdf, width = 6.5, height = 6.5)
    print(km_plot)
    grDevices::dev.off()
  }

  if (length(logrank_results) > 0L) {
    utils::write.csv(
      do.call(rbind, logrank_results),
      file.path(endpoint_dir, paste0(endpoint_name, "_logrank_tests.csv")),
      row.names = FALSE
    )
  }

  # Univariable Cox analyses.
  univariable_table <- batch_univariable_cox(
    data,
    time = time_column,
    event = event_column,
    variables = cox_candidate_variables
  )

  if (nrow(univariable_table) == 0L) {
    next
  }
  utils::write.csv(
    univariable_table,
    file.path(endpoint_dir, paste0(endpoint_name, "_univariable_Cox.csv")),
    row.names = FALSE
  )

  model_p_values <- tapply(
    univariable_table$p.value_likelihood,
    univariable_table$feature,
    function(x) x[1]
  )
  selected_variables <- names(model_p_values)[
    is.finite(model_p_values) & model_p_values < 0.05
  ]

  if (length(selected_variables) == 0L) {
    warning("No univariable Cox variable met P < 0.05 for ", endpoint_name, ".")
    next
  }

  multivariable_columns <- c(time_column, event_column, selected_variables)
  multivariable_data <- data[
    stats::complete.cases(data[, multivariable_columns, drop = FALSE]),
    multivariable_columns,
    drop = FALSE
  ]

  if (sum(multivariable_data[[event_column]]) <= length(selected_variables)) {
    warning(
      endpoint_name,
      " has no more events than selected variables; estimates may be unstable."
    )
  }

  multivariable_formula <- stats::as.formula(
    paste0(
      "survival::Surv(", time_column, ", ", event_column, ") ~ ",
      paste(selected_variables, collapse = " + ")
    )
  )
  multivariable_fit <- survival::coxph(
    multivariable_formula,
    data = multivariable_data,
    ties = "efron",
    x = TRUE
  )

  multivariable_table <- cox_table(multivariable_fit, feature = "multivariable")
  multivariable_table$n <- nrow(multivariable_data)
  multivariable_table$events <- sum(multivariable_data[[event_column]])
  utils::write.csv(
    multivariable_table,
    file.path(endpoint_dir, paste0(endpoint_name, "_multivariable_Cox.csv")),
    row.names = FALSE
  )

  proportional_hazards <- tryCatch(
    survival::cox.zph(multivariable_fit),
    error = function(e) NULL
  )
  if (!is.null(proportional_hazards)) {
    ph_table <- data.frame(
      term = rownames(proportional_hazards$table),
      proportional_hazards$table,
      row.names = NULL,
      check.names = FALSE
    )
    utils::write.csv(
      ph_table,
      file.path(endpoint_dir, paste0(endpoint_name, "_proportional_hazards_test.csv")),
      row.names = FALSE
    )
  }

  forest_plot <- survminer::ggforest(
    multivariable_fit,
    data = multivariable_data,
    main = paste0(endpoint_name, " multivariable Cox model")
  )
  ggplot2::ggsave(
    file.path(endpoint_dir, paste0(endpoint_name, "_multivariable_forest.pdf")),
    plot = forest_plot,
    width = 7,
    height = max(4.5, 1 + 0.55 * length(stats::coef(multivariable_fit)))
  )
}

message("Survival-analysis results saved to: ", output_root)
# mIHC and cytokine analyses for Supplementary Figures S1-S3.

check_packages(c("dplyr", "ggplot2", "tidyr"))

input_file <- file.path("data", "immune_cytokine_long.csv")
output_dir <- file.path("results", "immune_cytokine")
ensure_directory(output_dir)

data <- read_required_csv(
  input_file,
  required_columns = c("patient_id", "assay", "marker", "timepoint", "value", "CR")
)
data$patient_id <- trimws(as.character(data$patient_id))
data$assay <- trimws(as.character(data$assay))
data$marker <- trimws(as.character(data$marker))
data$value <- as.numeric(data$value)
data$CR <- as_binary(data$CR, "CR")

time_key <- tolower(trimws(as.character(data$timepoint)))
data$timepoint_standard <- ifelse(
  time_key %in% c("t1", "baseline", "pre", "pretreatment"),
  "Baseline",
  ifelse(time_key %in% c("t4", "post", "post-treatment", "posttreatment"), "Post", NA_character_)
)

if (any(is.na(data$timepoint_standard) & !is.na(data$timepoint))) {
  stop(
    "Unrecognized immune/cytokine timepoint(s): ",
    paste(unique(data$timepoint[is.na(data$timepoint_standard)]), collapse = ", ")
  )
}

data <- data[stats::complete.cases(data[, c("patient_id", "assay", "marker", "timepoint_standard", "value", "CR")]), ]

duplicate_measurements <- duplicated(data[, c("patient_id", "assay", "marker", "timepoint_standard")])
if (any(duplicate_measurements)) {
  stop("Duplicate patient/assay/marker/timepoint measurements were detected.")
}

paired_results <- list()
group_results <- list()
change_results <- list()
pair_index <- 1L
group_index <- 1L
change_index <- 1L

analysis_groups <- unique(data[, c("assay", "marker")])
for (i in seq_len(nrow(analysis_groups))) {
  current_assay <- analysis_groups$assay[i]
  current_marker <- analysis_groups$marker[i]
  subset <- data[data$assay == current_assay & data$marker == current_marker, ]

  wide <- tidyr::pivot_wider(
    subset[, c("patient_id", "CR", "timepoint_standard", "value")],
    id_cols = c("patient_id", "CR"),
    names_from = "timepoint_standard",
    values_from = "value"
  )

  if (all(c("Baseline", "Post") %in% names(wide))) {
    paired <- wide[stats::complete.cases(wide[, c("Baseline", "Post")]), ]

    if (nrow(paired) >= 3L) {
      paired_t <- tryCatch(
        stats::t.test(paired$Post, paired$Baseline, paired = TRUE),
        error = function(e) NULL
      )
      paired_wilcoxon <- tryCatch(
        suppressWarnings(stats::wilcox.test(
          paired$Post,
          paired$Baseline,
          paired = TRUE,
          exact = FALSE
        )),
        error = function(e) NULL
      )

      paired_results[[pair_index]] <- data.frame(
        assay = current_assay,
        marker = current_marker,
        n_pairs = nrow(paired),
        mean_baseline = mean(paired$Baseline),
        mean_post = mean(paired$Post),
        paired_t_p.value = if (is.null(paired_t)) NA_real_ else paired_t$p.value,
        paired_wilcoxon_p.value = if (is.null(paired_wilcoxon)) NA_real_ else paired_wilcoxon$p.value,
        stringsAsFactors = FALSE
      )
      pair_index <- pair_index + 1L

      paired$change <- paired$Post - paired$Baseline
      group_sizes <- table(paired$CR)
      if (length(group_sizes) == 2L && all(group_sizes >= 2L)) {
        change_test <- tryCatch(
          stats::t.test(change ~ CR, data = paired, var.equal = FALSE),
          error = function(e) NULL
        )
        change_results[[change_index]] <- data.frame(
          assay = current_assay,
          marker = current_marker,
          n = nrow(paired),
          mean_change_CR = mean(paired$change[paired$CR == 1]),
          mean_change_non_CR = mean(paired$change[paired$CR == 0]),
          p.value_Welch = if (is.null(change_test)) NA_real_ else change_test$p.value,
          stringsAsFactors = FALSE
        )
        change_index <- change_index + 1L
      }
    }
  }

  for (current_timepoint in unique(subset$timepoint_standard)) {
    group_data <- subset[subset$timepoint_standard == current_timepoint, ]
    response_counts <- table(group_data$CR)
    if (nrow(group_data) >= 4L && length(response_counts) == 2L && all(response_counts >= 2L)) {
      group_test <- tryCatch(
        stats::t.test(value ~ CR, data = group_data, var.equal = FALSE),
        error = function(e) NULL
      )
      group_results[[group_index]] <- data.frame(
        assay = current_assay,
        marker = current_marker,
        timepoint = current_timepoint,
        n = nrow(group_data),
        mean_CR = mean(group_data$value[group_data$CR == 1]),
        mean_non_CR = mean(group_data$value[group_data$CR == 0]),
        p.value_Welch = if (is.null(group_test)) NA_real_ else group_test$p.value,
        stringsAsFactors = FALSE
      )
      group_index <- group_index + 1L
    }
  }
}

if (length(paired_results) > 0L) {
  paired_table <- do.call(rbind, paired_results)
  paired_table$paired_t_p.adjust_BH <- stats::p.adjust(paired_table$paired_t_p.value, method = "BH")
  paired_table$paired_wilcoxon_p.adjust_BH <- stats::p.adjust(paired_table$paired_wilcoxon_p.value, method = "BH")
  utils::write.csv(
    paired_table,
    file.path(output_dir, "baseline_post_paired_tests.csv"),
    row.names = FALSE
  )
}

if (length(group_results) > 0L) {
  group_table <- do.call(rbind, group_results)
  group_table$p.adjust_BH <- stats::p.adjust(group_table$p.value_Welch, method = "BH")
  utils::write.csv(
    group_table,
    file.path(output_dir, "CR_nonCR_Welch_tests.csv"),
    row.names = FALSE
  )
}

if (length(change_results) > 0L) {
  change_table <- do.call(rbind, change_results)
  change_table$p.adjust_BH <- stats::p.adjust(change_table$p.value_Welch, method = "BH")
  utils::write.csv(
    change_table,
    file.path(output_dir, "change_from_baseline_by_response.csv"),
    row.names = FALSE
  )
}

data$Response <- factor(data$CR, levels = c(1, 0), labels = c("CR", "non-CR"))
for (current_assay in unique(data$assay)) {
  plot_data <- data[data$assay == current_assay, ]
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$timepoint_standard,
      y = .data$value,
      colour = .data$Response
    )
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = interaction(.data$timepoint_standard, .data$Response)),
      outlier.shape = NA,
      position = ggplot2::position_dodge(width = 0.7)
    ) +
    ggplot2::geom_jitter(
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
      size = 1.2,
      alpha = 0.7
    ) +
    ggplot2::facet_wrap(~marker, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c(CR = "#2A9D8F", `non-CR` = "#E76F51")) +
    ggplot2::labs(x = NULL, y = "Measurement", colour = "Response") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(
    file.path(output_dir, paste0("Supplementary_", safe_file_label(current_assay), "_plots.pdf")),
    plot = plot,
    width = 10,
    height = max(5, 2.5 * ceiling(length(unique(plot_data$marker)) / 4))
  )
}

message("Immune/cytokine results saved to: ", output_dir)

