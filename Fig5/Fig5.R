# Figure 5

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

