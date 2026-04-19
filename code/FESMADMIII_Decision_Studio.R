

suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(readxl)
  library(writexl)
  library(ggplot2)
  library(DT)
  library(magrittr)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

EPS <- 1e-12

restore_shape <- function(ref, y) {
  if (is.matrix(ref) || is.array(ref)) {
    y <- array(as.vector(y), dim = dim(ref), dimnames = dimnames(ref))
  } else if (!is.null(names(ref)) && length(y) == length(ref)) {
    names(y) <- names(ref)
  }
  y
}

as_num_vec <- function(x, allow_empty = FALSE) {
  if (is.null(x)) {
    if (allow_empty) return(numeric(0))
    stop("Encountered NULL where a numeric vector was required.")
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x) || is.array(x)) x <- as.vector(x)
  x <- suppressWarnings(as.numeric(unname(x)))
  if (!allow_empty && length(x) == 0) stop("Encountered an empty numeric vector during computation.")
  x
}

clamp01 <- function(x) restore_shape(x, pmin(1, pmax(0, x)))

safe_div <- function(num, den, eps = EPS) {
  num <- as_num_vec(num)
  den <- as_num_vec(den)
  if (length(den) == 1L && length(num) > 1L) den <- rep(den, length(num))
  if (length(num) == 1L && length(den) > 1L) num <- rep(num, length(den))
  if (length(num) != length(den)) stop("safe_div received incompatible lengths.")
  num / pmax(den, eps)
}

renorm_simplex <- function(x, eps = EPS) {
  x <- pmax(as_num_vec(x), 0)
  if (length(x) == 0) stop("Cannot normalize an empty simplex vector.")
  s <- sum(x)
  if (!is.finite(s) || s <= eps) return(rep(1 / length(x), length(x)))
  x / s
}

entropy_raw <- function(p, eps = EPS) {
  p <- pmax(as_num_vec(p), 0)
  s <- sum(p)
  if (!is.finite(s) || s <= eps) return(0)
  p <- p / s
  nz <- p > eps
  -sum(p[nz] * log2(p[nz]))
}

project_box_simplex <- function(target, lower, upper, total = 1, eps = EPS) {
  target <- as_num_vec(target)
  repaired <- repair_interval_bounds(lower, upper, total = total, eps = eps)
  lower <- repaired$lower
  upper <- repaired$upper

  if (length(target) != length(lower)) stop("Projection target and bounds must have the same length.")
  if (sum(lower) > total + 1e-10 || sum(upper) < total - 1e-10) {
    stop("Projection bounds do not define a nonempty box-constrained simplex.")
  }

  phi <- function(lambda) sum(pmin(pmax(target - lambda, lower), upper)) - total

  lo <- min(target - upper) - total - 1
  hi <- max(target - lower) + total + 1

  for (iter in seq_len(300)) {
    mid <- (lo + hi) / 2
    val <- phi(mid)
    if (abs(val) <= 1e-13) break
    if (val > 0) lo <- mid else hi <- mid
  }

  x <- pmin(pmax(target - (lo + hi) / 2, lower), upper)
  diff <- total - sum(x)

  if (abs(diff) > 1e-11) {
    if (diff > 0) {
      slack <- upper - x
      ord <- order(slack, decreasing = TRUE)
      for (i in ord) {
        add <- min(diff, slack[i])
        x[i] <- x[i] + add
        diff <- diff - add
        if (diff <= 1e-12) break
      }
    } else {
      diff <- -diff
      slack <- x - lower
      ord <- order(slack, decreasing = TRUE)
      for (i in ord) {
        sub <- min(diff, slack[i])
        x[i] <- x[i] - sub
        diff <- diff - sub
        if (diff <= 1e-12) break
      }
    }
  }

  x <- pmin(pmax(x, lower), upper)
  if (abs(sum(x) - total) > 1e-8) {
    stop("Projection failed to satisfy the simplex equality within tolerance.")
  }
  x
}

mid_interval <- function(lower, upper) {
  lower <- as_num_vec(lower)
  upper <- as_num_vec(upper)
  if (length(lower) != length(upper)) stop("mid_interval received incompatible lengths.")
  (lower + upper) / 2
}

fmt6 <- function(x) sprintf("%.6f", x)
fmt4 <- function(x) sprintf("%.4f", x)

interval_simplex_bounds <- function(lower, upper, eps = EPS) {
  lower <- as_num_vec(lower)
  upper <- as_num_vec(upper)
  if (length(lower) != length(upper)) stop("Lower and upper bounds must have the same length.")
  m <- length(lower)

  outL <- numeric(m)
  outU <- numeric(m)

  for (i in seq_len(m)) {
    denL <- lower[i] + sum(upper[-i])
    denU <- upper[i] + sum(lower[-i])
    outL[i] <- safe_div(lower[i], denL, eps = eps)
    outU[i] <- safe_div(upper[i], denU, eps = eps)
  }

  outL <- clamp01(outL)
  outU <- clamp01(pmax(outU, outL))
  list(lower = outL, upper = outU)
}

repair_interval_bounds <- function(lower, upper, total = 1, eps = EPS) {
  lower <- as_num_vec(lower)
  upper <- as_num_vec(upper)
  if (length(lower) != length(upper)) stop("Interval bounds must have the same length.")
  if (length(lower) == 0) stop("Interval bounds cannot be empty.")

  lower[!is.finite(lower)] <- 0
  upper[!is.finite(upper)] <- 0

  lower <- pmax(lower, 0)
  upper <- pmax(upper, lower)
  lower <- pmin(lower, total)
  upper <- pmin(upper, total)

  sL <- sum(lower)
  if (sL > total + 1e-10) {
    lower <- lower / sL * total
    upper <- pmax(upper, lower)
    upper <- pmin(upper, total)
  }

  sU <- sum(upper)
  if (sU < total - 1e-10) {
    room <- pmax(total - upper, 0)
    if (sum(room) <= eps) {
      upper <- rep(total / length(upper), length(upper))
    } else {
      upper <- upper + (total - sU) * room / sum(room)
      upper <- pmax(upper, lower)
      upper <- pmin(upper, total)
    }
  }

  list(lower = lower, upper = upper)
}

entropy_min_vector <- function(lower, upper, total = 1, eps = EPS) {
  repaired <- repair_interval_bounds(lower, upper, total = total, eps = eps)
  lower <- repaired$lower
  upper <- repaired$upper

  p <- lower
  residual <- total - sum(p)
  if (residual <= eps) return(renorm_simplex(p, eps = eps))

  ord <- order(upper, decreasing = TRUE)
  for (i in ord) {
    cap <- upper[i] - p[i]
    add <- min(residual, cap)
    p[i] <- p[i] + add
    residual <- residual - add
    if (residual <= eps) break
  }

  renorm_simplex(p, eps = eps)
}

entropy_max_vector <- function(lower, upper, total = 1, eps = EPS) {
  repaired <- repair_interval_bounds(lower, upper, total = total, eps = eps)
  lower <- repaired$lower
  upper <- repaired$upper

  f <- function(cval) {
    sum(pmin(pmax(cval, lower), upper)) - total
  }

  lo <- min(lower)
  hi <- max(upper)
  for (iter in seq_len(200)) {
    mid <- (lo + hi) / 2
    val <- f(mid)
    if (abs(val) <= 1e-12) break
    if (val > 0) hi <- mid else lo <- mid
  }

  p <- pmin(pmax((lo + hi) / 2, lower), upper)
  diff <- total - sum(p)

  if (abs(diff) > 1e-10) {
    if (diff > 0) {
      slack <- upper - p
      ord <- order(slack, decreasing = TRUE)
      for (i in ord) {
        add <- min(diff, slack[i])
        p[i] <- p[i] + add
        diff <- diff - add
        if (diff <= 1e-12) break
      }
    } else {
      diff <- -diff
      slack <- p - lower
      ord <- order(slack, decreasing = TRUE)
      for (i in ord) {
        sub <- min(diff, slack[i])
        p[i] <- p[i] - sub
        diff <- diff - sub
        if (diff <= 1e-12) break
      }
    }
  }

  renorm_simplex(p, eps = eps)
}

entropy_interval_from_bounds <- function(lower, upper, n_alts, eps = EPS) {
  p_min <- entropy_min_vector(lower, upper, total = 1, eps = eps)
  p_max <- entropy_max_vector(lower, upper, total = 1, eps = eps)

  H_lower_raw <- entropy_raw(p_min, eps = eps)
  H_upper_raw <- entropy_raw(p_max, eps = eps)

  if (n_alts <= 1) {
    h_lower <- 0
    h_upper <- 0
  } else {
    h_lower <- H_lower_raw / log2(n_alts)
    h_upper <- H_upper_raw / log2(n_alts)
  }

  list(
    p_min = as_num_vec(p_min),
    p_max = as_num_vec(p_max),
    H_lower_raw = H_lower_raw,
    H_upper_raw = H_upper_raw,
    h_lower = h_lower,
    h_upper = h_upper
  )
}

alpha_cut_bounds <- function(central, delta, alpha, eps = EPS) {
  central <- as.numeric(central)
  delta <- pmax(as.numeric(delta), 0)
  spread <- (1 - alpha) * delta
  lower <- central - spread
  upper <- central + spread
  upper <- pmax(upper, lower)
  list(lower = lower, upper = upper)
}

read_sheet_safe <- function(path, sheet_name) {
  out <- NULL
  tryCatch({
    out <- as.data.frame(read_excel(path, sheet = sheet_name), stringsAsFactors = FALSE)
    names(out) <- trimws(names(out))
  }, error = function(e) {
    out <<- NULL
  })
  out
}

extract_single_alpha <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  nums <- suppressWarnings(as.numeric(unlist(df)))
  nums <- nums[is.finite(nums)]
  if (length(nums) == 0) return(NULL)
  nums[1]
}

is_blank_cell <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

drop_empty_rows_cols <- function(df) {
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) return(df)
  keep_rows <- !apply(df, 1, function(r) all(is_blank_cell(r)))
  keep_cols <- !vapply(df, function(col) all(is_blank_cell(col)), logical(1))
  df[keep_rows, keep_cols, drop = FALSE]
}

clean_matrix_sheet <- function(df, sheet_name) {
  df <- drop_empty_rows_cols(df)
  if (is.null(df) || ncol(df) < 3) {
    stop(paste0(sheet_name, " must contain one criterion column and at least two alternative columns."))
  }

  crit <- trimws(as.character(df[[1]]))
  num_df <- as.data.frame(
    lapply(df[, -1, drop = FALSE], function(col) suppressWarnings(as.numeric(col))),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  keep <- !(is.na(crit) | crit == "") & rowSums(!is.na(num_df)) > 0
  df2 <- data.frame(
    Criterion = crit[keep],
    num_df[keep, , drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(df2) < 2) stop(paste0(sheet_name, " must contain at least two valid criteria rows."))
  if (anyDuplicated(df2$Criterion)) stop(paste0(sheet_name, " contains duplicate criterion names."))
  if (anyNA(as.matrix(df2[, -1, drop = FALSE]))) {
    stop(paste0(sheet_name, " contains non-numeric or empty values in the alternatives block."))
  }

  df2
}

clean_subjective_weights_sheet <- function(df) {
  df <- drop_empty_rows_cols(df)
  need <- c("Criterion", "Central", "Delta")
  if (!all(need %in% names(df))) stop("FuzzySubjectiveWeights must contain: Criterion, Central, Delta.")

  crit <- trimws(as.character(df$Criterion))
  central <- suppressWarnings(as.numeric(df$Central))
  delta <- suppressWarnings(as.numeric(df$Delta))
  delta[!is.finite(delta)] <- 0

  keep <- !(is.na(crit) | crit == "") & is.finite(central)
  out <- data.frame(
    Criterion = crit[keep],
    Central = central[keep],
    Delta = pmax(delta[keep], 0),
    stringsAsFactors = FALSE
  )

  if (nrow(out) < 2) stop("FuzzySubjectiveWeights must contain at least two valid rows.")
  if (anyDuplicated(out$Criterion)) stop("FuzzySubjectiveWeights contains duplicate criterion names.")
  out
}

clean_benefit_cost_sheet <- function(df) {
  df <- drop_empty_rows_cols(df)
  need <- c("Criterion", "Type")
  if (!all(need %in% names(df))) stop("BenefitCost must contain: Criterion and Type.")

  crit <- trimws(as.character(df$Criterion))
  type <- tolower(trimws(as.character(df$Type)))
  keep <- !(is.na(crit) | crit == "") & !(is.na(type) | type == "")
  out <- data.frame(Criterion = crit[keep], Type = type[keep], stringsAsFactors = FALSE)

  if (nrow(out) < 2) stop("BenefitCost must contain at least two valid rows.")
  if (anyDuplicated(out$Criterion)) stop("BenefitCost contains duplicate criterion names.")
  bad_types <- setdiff(unique(out$Type), c("benefit", "cost"))
  if (length(bad_types) > 0) stop("BenefitCost Type values must be only 'benefit' or 'cost'.")
  out
}

build_model_bundle <- function(dfMatrix, dfDelta, dfSbj, dfBC) {
  dfMatrix <- clean_matrix_sheet(dfMatrix, "DataMatrix")
  dfDelta <- clean_matrix_sheet(dfDelta, "FuzzyDeviations")
  dfSbj <- clean_subjective_weights_sheet(dfSbj)
  dfBC <- clean_benefit_cost_sheet(dfBC)

  if (!identical(names(dfMatrix), names(dfDelta))) stop("DataMatrix and FuzzyDeviations must have identical alternative columns.")

  critNames <- trimws(as.character(dfMatrix[[1]]))
  altNames <- trimws(names(dfMatrix)[-1])

  dataMat <- as.matrix(dfMatrix[, -1, drop = FALSE])
  deltaMat <- as.matrix(dfDelta[, -1, drop = FALSE])

  storage.mode(dataMat) <- "numeric"
  storage.mode(deltaMat) <- "numeric"

  rownames(dataMat) <- critNames
  colnames(dataMat) <- altNames
  rownames(deltaMat) <- trimws(as.character(dfDelta[[1]]))
  colnames(deltaMat) <- trimws(names(dfDelta)[-1])

  list(
    dataMat = dataMat,
    deltaMat = deltaMat,
    sbjDF = dfSbj,
    bcDF = dfBC
  )
}

validate_bundle <- function(dataMat, deltaMat, sbjDF, bcDF, alpha) {
  if (!identical(dim(dataMat), dim(deltaMat))) stop("DataMatrix and FuzzyDeviations must have identical dimensions.")
  if (!identical(rownames(dataMat), rownames(deltaMat))) stop("Criteria names must match between DataMatrix and FuzzyDeviations.")
  if (!identical(colnames(dataMat), colnames(deltaMat))) stop("Alternative names must match between DataMatrix and FuzzyDeviations.")
  if (nrow(dataMat) < 2) stop("At least two criteria are required.")
  if (ncol(dataMat) < 2) stop("At least two alternatives are required.")
  if (any(!is.finite(dataMat))) stop("DataMatrix contains non-numeric or missing values.")
  if (any(!is.finite(deltaMat))) stop("FuzzyDeviations contains non-numeric or missing values.")
  if (any(deltaMat < 0, na.rm = TRUE)) stop("FuzzyDeviations cannot contain negative values.")

  critNames <- rownames(dataMat)
  if (!all(c("Criterion", "Central", "Delta") %in% names(sbjDF))) {
    stop("FuzzySubjectiveWeights must contain: Criterion, Central, Delta.")
  }
  if (!all(c("Criterion", "Type") %in% names(bcDF))) {
    stop("BenefitCost must contain: Criterion and Type.")
  }

  miss_sbj <- setdiff(critNames, sbjDF$Criterion)
  miss_bc <- setdiff(critNames, bcDF$Criterion)
  if (length(miss_sbj) > 0) stop(paste("Missing subjective weights for:", paste(miss_sbj, collapse = ", ")))
  if (length(miss_bc) > 0) stop(paste("Missing Benefit/Cost type for:", paste(miss_bc, collapse = ", ")))

  sbj_row <- sbjDF[match(critNames, sbjDF$Criterion), , drop = FALSE]
  acut <- alpha_cut_bounds(sbj_row$Central, sbj_row$Delta, alpha)
  if (any(acut$lower < 0)) {
    stop("The selected alpha produces a negative lower bound in the subjective-weight alpha-cut. Adjust Central/Delta or alpha.")
  }

  TRUE
}

matrix_to_df <- function(mat, row_label = "Criterion") {
  mat <- as.matrix(mat)
  rn <- rownames(mat)
  if (is.null(rn)) rn <- paste0("Row_", seq_len(nrow(mat)))
  df <- data.frame(tmp = rn, stringsAsFactors = FALSE, check.names = FALSE)
  names(df)[1] <- row_label
  for (j in seq_len(ncol(mat))) {
    nm <- colnames(mat)[j] %||% paste0("V", j)
    df[[nm]] <- mat[, j]
  }
  df
}

interval_matrix_to_df <- function(lower, upper, row_label = "Criterion") {
  lower <- as.matrix(lower)
  upper <- as.matrix(upper)
  if (!identical(dim(lower), dim(upper))) stop("lower and upper interval matrices must have identical dimensions.")
  rn <- rownames(lower)
  if (is.null(rn)) rn <- paste0("Row_", seq_len(nrow(lower)))
  cn <- colnames(lower)
  if (is.null(cn)) cn <- paste0("V", seq_len(ncol(lower)))

  out <- data.frame(tmp = rn, stringsAsFactors = FALSE, check.names = FALSE)
  names(out)[1] <- row_label
  for (j in seq_len(ncol(lower))) {
    out[[paste0(cn[j], "_L")]] <- lower[, j]
    out[[paste0(cn[j], "_U")]] <- upper[, j]
  }
  out
}

compute_fesmadmii_v2 <- function(dataMat, deltaMat, sbjDF, bcDF, alpha = 0.5, eps = EPS) {
  alpha_num <- suppressWarnings(as.numeric(alpha)[1])
  if (!is.finite(alpha_num)) alpha_num <- 0.5
  alpha <- min(max(alpha_num, 0), 1)

  validate_bundle(dataMat, deltaMat, sbjDF, bcDF, alpha)

  critNames <- rownames(dataMat)
  altNames <- colnames(dataMat)
  M <- nrow(dataMat)
  N <- ncol(dataMat)
  logN <- log2(N)

  bc_map <- setNames(tolower(trimws(as.character(bcDF$Type))), as.character(bcDF$Criterion))
  critTypes <- bc_map[critNames]
  critTypes[is.na(critTypes)] <- "benefit"

  xi_lower <- matrix(0, nrow = M, ncol = N, dimnames = dimnames(dataMat))
  xi_upper <- matrix(0, nrow = M, ncol = N, dimnames = dimnames(dataMat))
  r_lower <- matrix(0, nrow = M, ncol = N, dimnames = dimnames(dataMat))
  r_upper <- matrix(0, nrow = M, ncol = N, dimnames = dimnames(dataMat))
  shifts <- numeric(M)
  names(shifts) <- critNames

  for (i in seq_len(M)) {
    acut <- alpha_cut_bounds(dataMat[i, ], deltaMat[i, ], alpha, eps = eps)
    rawL <- as.numeric(acut$lower)
    rawU <- as.numeric(acut$upper)

    minLower <- min(rawL, na.rm = TRUE)
    shiftVal <- if (is.finite(minLower) && minLower <= 0) abs(minLower) + eps else 0
    shifts[i] <- shiftVal

    rawL <- rawL + shiftVal
    rawU <- rawU + shiftVal

    rawL <- pmax(rawL, eps)
    rawU <- pmax(rawU, rawL + eps)

    xi_lower[i, ] <- rawL
    xi_upper[i, ] <- rawU

    if (critTypes[i] == "cost") {
      refVal <- min(rawL, na.rm = TRUE)
      r_lower[i, ] <- safe_div(refVal, rawU, eps = eps)
      r_upper[i, ] <- safe_div(refVal, rawL, eps = eps)
    } else {
      refVal <- max(rawU, na.rm = TRUE)
      r_lower[i, ] <- safe_div(rawL, refVal, eps = eps)
      r_upper[i, ] <- safe_div(rawU, refVal, eps = eps)
    }
  }

  r_lower <- clamp01(r_lower)
  r_upper <- clamp01(pmax(r_upper, r_lower))

  rowSumL <- rowSums(r_lower)
  rowSumU <- rowSums(r_upper)

  P_lower <- matrix(0, nrow = M, ncol = N, dimnames = dimnames(dataMat))
  P_upper <- matrix(0, nrow = M, ncol = N, dimnames = dimnames(dataMat))

  for (i in seq_len(M)) {
    P_lower[i, ] <- safe_div(r_lower[i, ], rowSumU[i], eps = eps)
    P_upper[i, ] <- safe_div(r_upper[i, ], rowSumL[i], eps = eps)
  }

  feasibility_table <- data.frame(
    Criterion = critNames,
    Sum_P_Lower = rowSums(P_lower),
    Sum_P_Upper = rowSums(P_upper),
    Feasible = (rowSums(P_lower) <= 1 + 1e-10) & (rowSums(P_upper) >= 1 - 1e-10),
    stringsAsFactors = FALSE
  )

  h_lower <- h_upper <- d_lower <- d_upper <- numeric(M)
  p_entropy_min <- matrix(0, nrow = M, ncol = N, dimnames = dimnames(dataMat))
  p_entropy_max <- matrix(0, nrow = M, ncol = N, dimnames = dimnames(dataMat))

  for (i in seq_len(M)) {
    eInt <- entropy_interval_from_bounds(P_lower[i, ], P_upper[i, ], n_alts = N, eps = eps)
    h_lower[i] <- max(0, min(1, eInt$h_lower))
    h_upper[i] <- max(h_lower[i], min(1, eInt$h_upper))
    d_lower[i] <- 1 - h_upper[i]
    d_upper[i] <- 1 - h_lower[i]
    p_entropy_min[i, ] <- eInt$p_min
    p_entropy_max[i, ] <- eInt$p_max
  }

  names(h_lower) <- names(h_upper) <- names(d_lower) <- names(d_upper) <- critNames

  sbj_row <- sbjDF[match(critNames, sbjDF$Criterion), , drop = FALSE]
  sbj_acut <- alpha_cut_bounds(sbj_row$Central, sbj_row$Delta, alpha, eps = eps)
  xSBJ_lower <- as_num_vec(sbj_acut$lower)
  xSBJ_upper <- as_num_vec(sbj_acut$upper)

  sbjInt <- interval_simplex_bounds(xSBJ_lower, xSBJ_upper, eps = eps)
  wSBJ_lower <- sbjInt$lower
  wSBJ_upper <- sbjInt$upper

  objInt <- interval_simplex_bounds(d_lower, d_upper, eps = eps)
  wOBJ_lower <- objInt$lower
  wOBJ_upper <- objInt$upper

  q_lower <- pmax(wSBJ_lower * wOBJ_lower, 0)
  q_upper <- pmax(wSBJ_upper * wOBJ_upper, q_lower)

  intInt <- interval_simplex_bounds(q_lower, q_upper, eps = eps)
  wINT_lower <- intInt$lower
  wINT_upper <- intInt$upper

  names(xSBJ_lower) <- names(xSBJ_upper) <- critNames
  names(wSBJ_lower) <- names(wSBJ_upper) <- critNames
  names(wOBJ_lower) <- names(wOBJ_upper) <- critNames
  names(q_lower) <- names(q_upper) <- critNames
  names(wINT_lower) <- names(wINT_upper) <- critNames

  P_mid <- (P_lower + P_upper) / 2
  P_star <- matrix(0, nrow = M, ncol = N, dimnames = dimnames(dataMat))
  for (i in seq_len(M)) {
    P_star[i, ] <- project_box_simplex(P_mid[i, ], P_lower[i, ], P_upper[i, ], total = 1, eps = eps)
  }

  wINT_mid <- (wINT_lower + wINT_upper) / 2
  wINT_star <- project_box_simplex(wINT_mid, wINT_lower, wINT_upper, total = 1, eps = eps)
  names(wINT_star) <- critNames

  scores <- as.numeric(wINT_star %*% P_star)
  names(scores) <- altNames

  S_mu_star <- numeric(M)
  for (i in seq_len(M)) S_mu_star[i] <- entropy_raw(P_star[i, ], eps = eps)
  names(S_mu_star) <- critNames

  ICI <- wINT_star * (1 - safe_div(S_mu_star, logN, eps = eps))
  names(ICI) <- critNames

  S_Y <- entropy_raw(scores, eps = eps)
  S_X <- entropy_raw(wINT_star, eps = eps)
  S_Y_given_X <- sum(wINT_star * S_mu_star)
  S_XY <- S_X + S_Y_given_X
  J_XY <- max(0, S_Y - S_Y_given_X)

  NMI <- if (S_Y <= eps) 1 else J_XY / S_Y
  NMI <- clamp01(NMI)
  CES <- clamp01(sum(ICI))
  ADI <- clamp01(1 - S_Y / logN)
  NMGI <- clamp01(mean(c(NMI, CES, ADI)))

  alternatives_table <- data.frame(
    Alternative = altNames,
    Score = scores,
    Rank = rank(-scores, ties.method = "min"),
    stringsAsFactors = FALSE
  )
  alternatives_table <- alternatives_table[order(alternatives_table$Rank, alternatives_table$Alternative), ]

  criteria_table <- data.frame(
    Criterion = critNames,
    Type = critTypes,
    xSBJ_L = xSBJ_lower,
    xSBJ_U = xSBJ_upper,
    wSBJ_L = wSBJ_lower,
    wSBJ_U = wSBJ_upper,
    h_L = h_lower,
    h_U = h_upper,
    d_L = d_lower,
    d_U = d_upper,
    wOBJ_L = wOBJ_lower,
    wOBJ_U = wOBJ_upper,
    q_L = q_lower,
    q_U = q_upper,
    wINT_L = wINT_lower,
    wINT_U = wINT_upper,
    wINT_star = wINT_star,
    S_mu_star = S_mu_star,
    ICI = ICI,
    stringsAsFactors = FALSE
  )

  diagnostics_table <- data.frame(
    Measure = c(
      "S_Y - Entropy of alternatives",
      "S_X - Entropy of criteria",
      "S_Y|X - Conditional entropy",
      "S_XY - Joint entropy",
      "J_X;Y - Mutual information",
      "NMI",
      "CES",
      "ADI",
      "NMGI"
    ),
    Value = c(S_Y, S_X, S_Y_given_X, S_XY, J_XY, NMI, CES, ADI, NMGI),
    stringsAsFactors = FALSE
  )

  entropy_table <- data.frame(
    Criterion = critNames,
    h_L = h_lower,
    h_U = h_upper,
    d_L = d_lower,
    d_U = d_upper,
    stringsAsFactors = FALSE
  )

  representative_prob_table <- matrix_to_df(P_star, row_label = "Criterion")
  representative_weight_table <- data.frame(
    Criterion = critNames,
    wINT_star = wINT_star,
    S_mu_star = S_mu_star,
    ICI = ICI,
    stringsAsFactors = FALSE
  )

  validation_flags <- data.frame(
    Check = c(
      "sum(wINT_star)=1",
      "each row of P_star sums to 1",
      "sum(scores)=1",
      "all diagnostics in [0,1]",
      "all P_star within interval bounds",
      "wINT_star within interval bounds"
    ),
    Result = c(
      abs(sum(wINT_star) - 1) <= 1e-8,
      all(abs(rowSums(P_star) - 1) <= 1e-8),
      abs(sum(scores) - 1) <= 1e-8,
      all(c(NMI, CES, ADI, NMGI) >= -1e-10 & c(NMI, CES, ADI, NMGI) <= 1 + 1e-10),
      all(P_star >= P_lower - 1e-10 & P_star <= P_upper + 1e-10),
      all(wINT_star >= wINT_lower - 1e-10 & wINT_star <= wINT_upper + 1e-10)
    ),
    stringsAsFactors = FALSE
  )

  summary_table <- data.frame(
    Item = c("Alpha", "Top Alternative", "Top Score", "NMI", "CES", "ADI", "NMGI"),
    Value = c(
      fmt6(alpha),
      alternatives_table$Alternative[1],
      fmt6(alternatives_table$Score[1]),
      fmt6(NMI),
      fmt6(CES),
      fmt6(ADI),
      fmt6(NMGI)
    ),
    stringsAsFactors = FALSE
  )

  list(
    alpha = alpha,
    M = M,
    N = N,
    critNames = critNames,
    altNames = altNames,
    critTypes = critTypes,
    shifts = shifts,
    xi_lower = xi_lower,
    xi_upper = xi_upper,
    r_lower = r_lower,
    r_upper = r_upper,
    rowSumL = rowSumL,
    rowSumU = rowSumU,
    P_lower = P_lower,
    P_upper = P_upper,
    p_entropy_min = p_entropy_min,
    p_entropy_max = p_entropy_max,
    h_lower = h_lower,
    h_upper = h_upper,
    d_lower = d_lower,
    d_upper = d_upper,
    xSBJ_lower = xSBJ_lower,
    xSBJ_upper = xSBJ_upper,
    wSBJ_lower = wSBJ_lower,
    wSBJ_upper = wSBJ_upper,
    wOBJ_lower = wOBJ_lower,
    wOBJ_upper = wOBJ_upper,
    q_lower = q_lower,
    q_upper = q_upper,
    wINT_lower = wINT_lower,
    wINT_upper = wINT_upper,
    P_star = P_star,
    wINT_star = wINT_star,
    scores = scores,
    S_mu_star = S_mu_star,
    ICI = ICI,
    S_Y = S_Y,
    S_X = S_X,
    S_Y_given_X = S_Y_given_X,
    S_XY = S_XY,
    J_XY = J_XY,
    NMI = NMI,
    CES = CES,
    ADI = ADI,
    NMGI = NMGI,
    feasibility_table = feasibility_table,
    validation_flags = validation_flags,
    entropy_table = entropy_table,
    criteria_table = criteria_table,
    alternatives_table = alternatives_table,
    diagnostics_table = diagnostics_table,
    representative_prob_table = representative_prob_table,
    representative_weight_table = representative_weight_table,
    summary_table = summary_table
  )
}

scenario_summary_row <- function(name, res, source = "Manual", origin = "") {
  data.frame(
    Scenario = name,
    Source = source,
    Origin = origin,
    Alpha = res$alpha,
    Winner = res$alternatives_table$Alternative[1],
    WinnerScore = res$alternatives_table$Score[1],
    NMI = res$NMI,
    CES = res$CES,
    ADI = res$ADI,
    NMGI = res$NMGI,
    stringsAsFactors = FALSE
  )
}

strip_extension <- function(x) sub("\\.[^.]+$", "", basename(x))

make_unique_scenario_name <- function(name, existing_names) {
  nm <- trimws(name)
  if (!nzchar(nm)) nm <- "Imported_Scenario"
  if (!(nm %in% existing_names)) return(nm)
  base <- nm
  k <- 1
  repeat {
    candidate <- paste0(base, "_", k)
    if (!(candidate %in% existing_names)) return(candidate)
    k <- k + 1
  }
}

validate_scenario_compatibility <- function(baseline_bundle, scenario_bundle) {
  if (!identical(rownames(baseline_bundle$dataMat), rownames(scenario_bundle$dataMat))) {
    stop("Imported scenario criteria do not match the baseline case study.")
  }
  if (!identical(colnames(baseline_bundle$dataMat), colnames(scenario_bundle$dataMat))) {
    stop("Imported scenario alternatives do not match the baseline case study.")
  }
  TRUE
}

compute_alpha_grid <- function(bundle, alpha_from, alpha_to, alpha_by) {
  if (alpha_by <= 0) stop("Alpha step must be positive.")
  if (alpha_to < alpha_from) stop("Alpha 'to' must be greater than or equal to alpha 'from'.")
  alphas <- seq(alpha_from, alpha_to, by = alpha_by)
  alphas <- round(pmin(pmax(alphas, 0), 1), 10)

  out <- lapply(alphas, function(a) {
    res <- compute_fesmadmii_v2(
      dataMat = bundle$dataMat,
      deltaMat = bundle$deltaMat,
      sbjDF = bundle$sbjDF,
      bcDF = bundle$bcDF,
      alpha = a
    )
    data.frame(
      Alpha = a,
      Winner = res$alternatives_table$Alternative[1],
      WinnerScore = res$alternatives_table$Score[1],
      NMI = res$NMI,
      CES = res$CES,
      ADI = res$ADI,
      NMGI = res$NMGI,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

safe_sheet_name <- function(x, used = character()) {
  x <- gsub("[\\\\/?*\\[\\]:]", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- ifelse(nchar(x) > 31, substr(x, 1, 31), x)
  base <- x
  k <- 1
  while (x %in% used) {
    suffix <- paste0("_", k)
    x <- paste0(substr(base, 1, max(1, 31 - nchar(suffix))), suffix)
    k <- k + 1
  }
  x
}

build_export_list <- function(res, prefix = "Baseline") {
  out <- list()
  out[[paste0(prefix, "_Summary")]] <- res$summary_table
  out[[paste0(prefix, "_AlphaCuts")]] <- interval_matrix_to_df(res$xi_lower, res$xi_upper, row_label = "Criterion")
  out[[paste0(prefix, "_NormalizedR")]] <- interval_matrix_to_df(res$r_lower, res$r_upper, row_label = "Criterion")
  out[[paste0(prefix, "_ProbIntervals")]] <- interval_matrix_to_df(res$P_lower, res$P_upper, row_label = "Criterion")
  out[[paste0(prefix, "_EntropyBounds")]] <- res$entropy_table
  out[[paste0(prefix, "_Criteria")]] <- res$criteria_table
  out[[paste0(prefix, "_RepresentativeP")]] <- res$representative_prob_table
  out[[paste0(prefix, "_RepresentativeW")]] <- res$representative_weight_table
  out[[paste0(prefix, "_Alternatives")]] <- res$alternatives_table
  out[[paste0(prefix, "_Diagnostics")]] <- res$diagnostics_table
  out[[paste0(prefix, "_Feasibility")]] <- res$feasibility_table
  out[[paste0(prefix, "_ValidationFlags")]] <- res$validation_flags
  out
}


export_summary_table <- function(res) {
  data.frame(
    Item = c("Alpha", "Top Alternative", "Top Score", "NMI", "CES", "ADI", "NMGI"),
    Value = c(
      sprintf("%.6f", res$alpha),
      as.character(res$alternatives_table$Alternative[1]),
      sprintf("%.12f", res$alternatives_table$Score[1]),
      sprintf("%.12f", res$NMI),
      sprintf("%.12f", res$CES),
      sprintf("%.12f", res$ADI),
      sprintf("%.12f", res$NMGI)
    ),
    stringsAsFactors = FALSE
  )
}

alternatives_export_df <- function(res) {
  df <- res$alternatives_table
  df[order(df$Alternative), , drop = FALSE]
}

alpha_grid_export_df <- function(alpha_grid) {
  if (is.null(alpha_grid) || nrow(alpha_grid) == 0) return(NULL)
  alpha_grid[order(alpha_grid$Alpha), , drop = FALSE]
}

build_validation_export_list <- function(base_res, alpha_grid = NULL) {
  out <- list(
    Baseline_Summary = export_summary_table(base_res),
    Baseline_AlphaCuts = interval_matrix_to_df(base_res$xi_lower, base_res$xi_upper, row_label = "Criterion"),
    Baseline_NormalizedR = interval_matrix_to_df(base_res$r_lower, base_res$r_upper, row_label = "Criterion"),
    Baseline_ProbIntervals = interval_matrix_to_df(base_res$P_lower, base_res$P_upper, row_label = "Criterion"),
    Baseline_EntropyBounds = base_res$entropy_table,
    Baseline_Criteria = base_res$criteria_table,
    Baseline_RepresentativeP = base_res$representative_prob_table,
    Baseline_RepresentativeW = base_res$representative_weight_table,
    Baseline_Alternatives = alternatives_export_df(base_res),
    Baseline_Diagnostics = base_res$diagnostics_table,
    Baseline_Feasibility = base_res$feasibility_table,
    Baseline_ValidationFlags = base_res$validation_flags
  )
  ag <- alpha_grid_export_df(alpha_grid)
  if (!is.null(ag)) out[["Alpha_Grid"]] <- ag
  out
}

build_full_export_list <- function(base_res, current_res, portfolio_summary, portfolio, alpha_grid,
                                   baseline_alpha, current_alpha) {
  export_list <- list(
    Run_Info = data.frame(
      Item = c("Export timestamp", "Baseline alpha", "Current scenario alpha", "Saved scenarios", "Alpha-grid rows"),
      Value = c(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        sprintf("%.6f", baseline_alpha),
        sprintf("%.6f", current_alpha),
        length(portfolio),
        if (is.null(alpha_grid)) 0 else nrow(alpha_grid)
      ),
      stringsAsFactors = FALSE
    )
  )

  tmp <- build_export_list(base_res, prefix = "Baseline")
  tmp[["Baseline_Summary"]] <- export_summary_table(base_res)
  tmp[["Baseline_Alternatives"]] <- alternatives_export_df(base_res)
  export_list <- c(export_list, tmp)

  tmp <- build_export_list(current_res, prefix = "CurrentScenario")
  tmp[["CurrentScenario_Summary"]] <- export_summary_table(current_res)
  tmp[["CurrentScenario_Alternatives"]] <- alternatives_export_df(current_res)
  export_list <- c(export_list, tmp)

  if (nrow(portfolio_summary) > 0) {
    export_list[["ScenarioPortfolio_Summary"]] <- portfolio_summary
    for (nm in names(portfolio)) {
      scen <- portfolio[[nm]]$result
      scen_prefix <- paste0("PF_", nm)
      export_list[[paste0(scen_prefix, "_Summary")]] <- data.frame(
        Item = c("Scenario", "Source", "Origin", "Alpha", "Winner", "WinnerScore", "NMI", "CES", "ADI", "NMGI"),
        Value = c(
          nm,
          portfolio[[nm]]$source %||% "Manual",
          portfolio[[nm]]$origin %||% "",
          sprintf("%.6f", scen$alpha),
          as.character(scen$alternatives_table$Alternative[1]),
          sprintf("%.12f", scen$alternatives_table$Score[1]),
          sprintf("%.12f", scen$NMI),
          sprintf("%.12f", scen$CES),
          sprintf("%.12f", scen$ADI),
          sprintf("%.12f", scen$NMGI)
        ),
        stringsAsFactors = FALSE
      )
      export_list[[paste0(scen_prefix, "_Alternatives")]] <- alternatives_export_df(scen)
      export_list[[paste0(scen_prefix, "_Criteria")]] <- scen$criteria_table
      export_list[[paste0(scen_prefix, "_Diagnostics")]] <- scen$diagnostics_table
    }
  }

  ag <- alpha_grid_export_df(alpha_grid)
  if (!is.null(ag)) {
    export_list[["Alpha_Grid"]] <- ag
    export_list[["Alpha_Grid_Best"]] <- ag[which.max(ag$NMGI), , drop = FALSE]
  }

  export_list
}


load_bundle_from_workbook <- function(path,
                                      sheet_matrix = "DataMatrix",
                                      sheet_delta = "FuzzyDeviations",
                                      sheet_weights = "FuzzySubjectiveWeights",
                                      sheet_bc = "BenefitCost",
                                      sheet_alpha = "GlobalAlpha",
                                      fallback_alpha = 0.60) {
  sheets <- readxl::excel_sheets(path)
  required <- c(sheet_matrix, sheet_delta, sheet_weights, sheet_bc)
  missing <- required[!required %in% sheets]
  if (length(missing) > 0) stop(paste("Missing required sheets:", paste(missing, collapse = ", ")))

  dfMatrix <- read_sheet_safe(path, sheet_matrix)
  dfDelta <- read_sheet_safe(path, sheet_delta)
  dfSbj <- read_sheet_safe(path, sheet_weights)
  dfBC <- read_sheet_safe(path, sheet_bc)
  alpha_df <- if (sheet_alpha %in% sheets) read_sheet_safe(path, sheet_alpha) else NULL

  alpha_global <- extract_single_alpha(alpha_df)
  alpha_final <- if (is.null(alpha_global)) fallback_alpha else alpha_global
  alpha_final <- min(max(alpha_final, 0), 1)

  bundle <- build_model_bundle(dfMatrix, dfDelta, dfSbj, dfBC)
  validate_bundle(bundle$dataMat, bundle$deltaMat, bundle$sbjDF, bundle$bcDF, alpha_final)
  list(bundle = bundle, alpha = alpha_final)
}

write_validation_workbook_from_input <- function(input_path,
                                                 output_path,
                                                 alpha_from = NULL,
                                                 alpha_to = NULL,
                                                 alpha_by = NULL,
                                                 fallback_alpha = 0.60) {
  loaded <- load_bundle_from_workbook(input_path, fallback_alpha = fallback_alpha)
  res <- compute_fesmadmii_v2(
    dataMat = loaded$bundle$dataMat,
    deltaMat = loaded$bundle$deltaMat,
    sbjDF = loaded$bundle$sbjDF,
    bcDF = loaded$bundle$bcDF,
    alpha = loaded$alpha
  )

  alpha_grid <- NULL
  if (!is.null(alpha_from) && !is.null(alpha_to) && !is.null(alpha_by)) {
    alpha_grid <- compute_alpha_grid(
      bundle = loaded$bundle,
      alpha_from = alpha_from,
      alpha_to = alpha_to,
      alpha_by = alpha_by
    )
  }

  export_list <- build_validation_export_list(res, alpha_grid = alpha_grid)
  writexl::write_xlsx(export_list, path = output_path)
  invisible(output_path)
}

choose_save_path <- function(default_name = "FESMADMIII_results.xlsx") {
  if (.Platform$OS.type == "windows" && requireNamespace("tcltk", quietly = TRUE)) {
    path <- tryCatch({
      tcltk::tclvalue(
        tcltk::tkgetSaveFile(
          defaultextension = ".xlsx",
          initialfile = default_name,
          filetypes = "{{Excel Workbook} {.xlsx}} {{All files} *}"
        )
      )
    }, error = function(e) "")
    path <- trimws(path)
    if (!nzchar(path)) return("")
    if (!grepl("\\.xlsx$", path, ignore.case = TRUE)) path <- paste0(path, ".xlsx")
    return(path)
  }

  path <- file.path(getwd(), default_name)
  if (!grepl("\\.xlsx$", path, ignore.case = TRUE)) path <- paste0(path, ".xlsx")
  path
}

find_chromium_app <- function() {
  candidates <- unique(c(
    Sys.which("chrome"),
    Sys.which("chrome.exe"),
    Sys.which("msedge"),
    Sys.which("msedge.exe"),
    "C:/Program Files/Google/Chrome/Application/chrome.exe",
    "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
    "C:/Program Files/Microsoft/Edge/Application/msedge.exe",
    "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
  ))
  candidates <- candidates[nzchar(candidates)]
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) return("")
  hits[1]
}

launch_standalone_window <- function(url) {
  browser_path <- find_chromium_app()

  if (nzchar(browser_path)) {
    args <- c(
      "--new-window",
      paste0("--app=", url),
      "--disable-session-crashed-bubble",
      "--disable-features=Translate"
    )
    system2(browser_path, args = args, wait = FALSE)
  } else {
    utils::browseURL(url)
  }
}

fes_dt <- function(df, pageLength = 10, scrollX = FALSE, dom = "tip", rownames = FALSE) {
  datatable(
    df,
    rownames = rownames,
    class = "compact stripe hover row-border order-column nowrap",
    options = list(
      pageLength = pageLength,
      autoWidth = TRUE,
      scrollX = scrollX,
      dom = dom
    )
  )
}

kpi_box <- function(title, value, note = NULL, fill = "#1F4E79") {
  div(
    class = "kpi-box",
    style = paste0("border-top: 4px solid ", fill, ";"),
    div(class = "kpi-title", title),
    div(class = "kpi-value", value),
    if (!is.null(note)) div(class = "kpi-note", note)
  )
}

ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(
    tags$style(HTML("
      body { background: #F6F9FC; }
      .main-title { font-weight: 800; color: #17324D; margin-bottom: 4px; }
      .sub-title { color: #52667A; margin-bottom: 18px; }
      .soft-card {
        background: #FFFFFF; border: 1px solid #E5EDF5; border-radius: 14px;
        padding: 16px; margin-bottom: 16px; box-shadow: 0 6px 18px rgba(20,40,80,0.06);
      }
      .kpi-box {
        background: #FFFFFF; border-radius: 14px; padding: 14px 16px; margin-bottom: 14px;
        box-shadow: 0 6px 18px rgba(20,40,80,0.06); border: 1px solid #E5EDF5;
        min-height: 104px;
      }
      .kpi-title {
        color: #607286; font-size: 12px; font-weight: 700; text-transform: uppercase;
        letter-spacing: 0.04em; margin-bottom: 6px;
      }
      .kpi-value {
        color: #17324D; font-size: 24px; font-weight: 800; line-height: 1.1;
      }
      .kpi-note {
        color: #66788A; font-size: 12px; margin-top: 6px;
      }
      .section-head {
        font-weight: 800; color: #17324D; margin-top: 0; margin-bottom: 12px;
      }
      .help-note { color: #61788B; font-size: 13px; line-height: 1.6; }
      .nav-tabs > li > a { font-weight: 700; }
      .well { background: #FFFFFF; border-radius: 14px; border: 1px solid #E5EDF5; box-shadow: none; }
      .control-label { font-weight: 700; color: #17324D; }
      .hero-card { padding-top: 18px; padding-bottom: 18px; }
      .creator-line { color: #476176; font-size: 13px; font-weight: 600; display: flex; align-items: center; gap: 8px; flex-wrap: wrap; margin-top: 10px; }
      .hero-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; }
      .hero-tile {
        background: linear-gradient(180deg, #F8FBFF 0%, #EEF5FC 100%);
        border: 1px solid #DEEAF6; border-radius: 14px; padding: 12px 12px 10px 12px;
        min-height: 118px; box-shadow: inset 0 1px 0 rgba(255,255,255,0.8);
      }
      .hero-icon-wrap {
        width: 42px; height: 42px; border-radius: 12px; display: flex; align-items: center; justify-content: center;
        background: #17324D; color: #FFFFFF; font-size: 18px; margin-bottom: 10px;
      }
      .hero-tile-title { color: #17324D; font-weight: 800; font-size: 13px; margin-bottom: 4px; }
      .hero-tile-note { color: #607286; font-size: 11.5px; line-height: 1.45; }
    "))
  ),

  div(class = "soft-card hero-card",
      fluidRow(
        column(
          8,
          h2(class = "main-title", "FES-MADM III Decision Studio"),
          p(class = "sub-title",
            "Validated computational studio for FES-MADM III: alpha-cuts, ratio normalization, interval conditional probabilities, entropy envelopes, subjective-objective-integrated weights, feasible Step-6 projections, scores, ICI, NMI, CES, ADI, NMGI, scenario portfolio, and alpha-grid analysis."),
          div(class = "creator-line",
              icon("user"),
              span("Creator: Lt Col (Ordnance) Dr Sideris Kiratsoudis"),
              tags$span(" | ", style = "margin: 0 6px;"),
              span("Postdoctoral Researcher, Department of Physics, Democritus University of Thrace")
          )
        ),
        column(
          4,
          div(class = "hero-grid",
              div(class = "hero-tile",
                  div(class = "hero-icon-wrap", icon("line-chart")),
                  div(class = "hero-tile-title", "Decision Analytics"),
                  div(class = "hero-tile-note", "Scoring, ranking, and scenario comparison")
              ),
              div(class = "hero-tile",
                  div(class = "hero-icon-wrap", icon("sitemap")),
                  div(class = "hero-tile-title", "Scenario Portfolio"),
                  div(class = "hero-tile-note", "Integrated multi-scenario storage and comparison")
              ),
              div(class = "hero-tile",
                  div(class = "hero-icon-wrap", icon("sliders")),
                  div(class = "hero-tile-title", "Alpha Control"),
                  div(class = "hero-tile-note", "Confidence-grid exploration and robustness")
              ),
              div(class = "hero-tile",
                  div(class = "hero-icon-wrap", icon("shield")),
                  div(class = "hero-tile-title", "Auditability"),
                  div(class = "hero-tile-note", "Traceable diagnostics and export-ready outputs")
              )
          )
        )
      )
  ),

  tabsetPanel(
    id = "mainTabs",

    tabPanel(
      "Overview",
      fluidRow(
        column(
          8,
          div(class = "soft-card",
              h3(class = "section-head", "Model scope"),
              p(class = "help-note", "The implementation follows the validated FES-MADM III compact formulation step by step: uncertain performances and subjective weights are transformed into alpha-cut intervals, normalized by criterion orientation, converted to interval conditional probabilities, processed through entropy-based weighting, and then aggregated through feasible Step-6 projections into representative scores and diagnostic indices.")
          ),
          uiOutput("overviewSummary")
        ),
        column(
          4,
          div(class = "soft-card",
              h3(class = "section-head", "Expected input workbook"),
              tags$ul(
                tags$li("DataMatrix"),
                tags$li("FuzzyDeviations"),
                tags$li("FuzzySubjectiveWeights"),
                tags$li("BenefitCost"),
                tags$li("GlobalAlpha (optional)")
              ),
              p(class = "help-note", "Multiple scenarios are handled inside the application through the scenario portfolio. Export uses a native save dialog and writes the workbook directly from the server side.")
          )
        )
      )
    ),

    tabPanel(
      "Data Import",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          fileInput("fileExcel", "Upload Excel workbook", accept = c(".xls", ".xlsx")),
          textInput("sheetMatrix", "DataMatrix sheet", "DataMatrix"),
          textInput("sheetDelta", "FuzzyDeviations sheet", "FuzzyDeviations"),
          textInput("sheetWeights", "FuzzySubjectiveWeights sheet", "FuzzySubjectiveWeights"),
          textInput("sheetBC", "BenefitCost sheet", "BenefitCost"),
          textInput("sheetAlphaGlobal", "Optional GlobalAlpha sheet", "GlobalAlpha"),
          numericInput("fallbackAlpha", "Fallback global alpha", value = 0.60, min = 0, max = 1, step = 0.05),
          actionButton("loadData", "Load workbook", class = "btn-primary"),
          hr(),
          p(class = "help-note", "After loading, the current workbook becomes the baseline model state. The scenario workspace starts as a copy of the baseline and can then be modified safely.")
        ),
        mainPanel(
          width = 8,
          div(class = "soft-card",
              h3(class = "section-head", "Workbook previews"),
              tabsetPanel(
                tabPanel("DataMatrix", br(), DTOutput("tableMatrixPreview")),
                tabPanel("FuzzyDeviations", br(), DTOutput("tableDeltaPreview")),
                tabPanel("FuzzySubjectiveWeights", br(), DTOutput("tableWeightsPreview")),
                tabPanel("BenefitCost", br(), DTOutput("tableBCPreview")),
                tabPanel("Alpha info", br(), DTOutput("tableAlphaInfo"))
              )
          )
        )
      )
    ),

    tabPanel(
      "Baseline",
      uiOutput("baselineKPIs"),
      tabsetPanel(
        tabPanel("Alternatives",
                 div(class = "soft-card", DTOutput("tableAltBaseline")),
                 div(class = "soft-card", plotOutput("plotAltBaseline", height = "420px"))),
        tabPanel("Criteria",
                 div(class = "soft-card", DTOutput("tableCritBaseline")),
                 div(class = "soft-card", plotOutput("plotWeightsBaseline", height = "420px")),
                 div(class = "soft-card", plotOutput("plotICIBaseline", height = "420px"))),
        tabPanel("Diagnostics",
                 div(class = "soft-card", DTOutput("tableDiagBaseline")),
                 div(class = "soft-card", plotOutput("plotDiagBaseline", height = "420px"))),
        tabPanel("Computation trace",
                 div(class = "soft-card", h4(class = "section-head", "Alpha-cut intervals"), DTOutput("tableAlphaCutsBaseline")),
                 div(class = "soft-card", h4(class = "section-head", "Normalized intervals"), DTOutput("tableNormBaseline")),
                 div(class = "soft-card", h4(class = "section-head", "Conditional probability intervals"), DTOutput("tableProbBaseline")),
                 div(class = "soft-card", h4(class = "section-head", "Entropy bounds"), DTOutput("tableEntropyBaseline")),
                 div(class = "soft-card", h4(class = "section-head", "Feasibility check"), DTOutput("tableFeasibilityBaseline")))
      )
    ),

    tabPanel(
      "Scenario Lab",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          numericInput("scenarioAlpha", "Scenario global alpha", value = 0.60, min = 0, max = 1, step = 0.05),
          hr(),
          h4("Modify one subjective weight"),
          selectInput("modWeightCrit", "Criterion", choices = NULL),
          numericInput("modWeightCentral", "New central weight", value = 0.25, step = 0.01),
          numericInput("modWeightDelta", "New delta", value = 0.02, step = 0.01),
          actionButton("applyWeightChange", "Apply weight change"),
          hr(),
          h4("Modify one data entry"),
          selectInput("modDataCrit", "Criterion", choices = NULL),
          selectInput("modDataAlt", "Alternative", choices = NULL),
          numericInput("modDataValue", "New central performance", value = 1, step = 0.1),
          numericInput("modDataDelta", "New fuzzy deviation", value = 0.1, step = 0.01),
          actionButton("applyDataChange", "Apply data/deviation change"),
          hr(),
          h4("Modify criterion type"),
          selectInput("modTypeCrit", "Criterion", choices = NULL),
          radioButtons("modTypeValue", "Type", choices = c("benefit", "cost"), inline = TRUE),
          actionButton("applyTypeChange", "Apply type change"),
          hr(),
          textInput("scenarioName", "Scenario name", "Scenario_1"),
          actionButton("saveScenario", "Save current scenario to portfolio", class = "btn-success"),
          actionButton("resetScenario", "Reset current scenario to baseline", class = "btn-warning")
        ),
        mainPanel(
          width = 8,
          uiOutput("scenarioKPIs"),
          tabsetPanel(
            tabPanel("Scenario outputs",
                     div(class = "soft-card", DTOutput("tableAltScenario")),
                     div(class = "soft-card", DTOutput("tableDiagScenario")),
                     div(class = "soft-card", plotOutput("plotAltScenarioCompare", height = "420px"))),
            tabPanel("Baseline vs current scenario",
                     div(class = "soft-card", DTOutput("tableCompScenario")),
                     div(class = "soft-card", plotOutput("plotCompScenario", height = "420px"))),
            tabPanel("Current scenario trace",
                     div(class = "soft-card", DTOutput("tableTraceScenario")))
          )
        )
      )
    ),

    tabPanel(
      "Scenario Import",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          fileInput("scenarioImportFiles", "Upload scenario Excel files", accept = c(".xls", ".xlsx"), multiple = TRUE),
          checkboxInput("importClearFirst", "Clear existing portfolio before this import", value = FALSE),
          radioButtons(
            "importDuplicateMode",
            "If a scenario name already exists",
            choices = c("Overwrite existing scenario" = "overwrite", "Keep both (auto-rename)" = "rename"),
            selected = "overwrite"
          ),
          actionButton("importScenarios", "Import scenarios to portfolio", class = "btn-primary"),
          hr(),
          p(class = "help-note",
            "Each imported workbook is treated as an independent scenario of the currently loaded baseline case study. The app validates that criteria and alternatives match the baseline, computes the full model, and stores the results directly in the shared Scenario Portfolio without overwriting earlier imported scenarios unless you explicitly choose overwrite.")
        ),
        mainPanel(
          width = 8,
          uiOutput("scenarioImportKPIs"),
          div(class = "soft-card", h3(class = "section-head", "Import log"), DTOutput("tableScenarioImportLog"))
        )
      )
    ),

    tabPanel(
      "Scenario Portfolio",
      fluidRow(
        column(
          4,
          div(class = "soft-card",
              h3(class = "section-head", "Scenario portfolio controls"),
              selectInput("portfolioRemoveName", "Saved scenario", choices = NULL),
              actionButton("removeScenario", "Remove selected scenario", class = "btn-danger"),
              actionButton("clearPortfolio", "Clear all scenarios"),
              p(class = "help-note", style = "margin-top:12px;",
                "Each saved scenario stores the full modified state (data, deviations, subjective weights, criterion types, and alpha). Scenarios may come either from the Scenario Lab or from imported Excel workbooks, and the comparison layer evaluates them together.")
          )
        ),
        column(
          8,
          uiOutput("portfolioKPIs")
        )
      ),
      div(class = "soft-card", DTOutput("tablePortfolio")),
      div(class = "soft-card", plotOutput("plotPortfolioNMGI", height = "420px")),
      div(class = "soft-card", plotOutput("plotPortfolioTopScore", height = "420px")),
      div(class = "soft-card", plotOutput("plotPortfolioHeatmap", height = "460px"))
    ),

    tabPanel(
      "Alpha Grid",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          numericInput("alphaFrom", "From alpha", value = 0.00, min = 0, max = 1, step = 0.05),
          numericInput("alphaTo", "To alpha", value = 1.00, min = 0, max = 1, step = 0.05),
          numericInput("alphaBy", "By", value = 0.10, min = 0.01, max = 1, step = 0.01),
          actionButton("runAlphaSweep", "Run alpha sweep", class = "btn-primary"),
          p(class = "help-note", style = "margin-top:12px;",
            "The alpha grid computes the model for every selected alpha and reports winner, winner score, NMI, CES, ADI, and NMGI.")
        ),
        mainPanel(
          width = 8,
          uiOutput("alphaGridKPIs"),
          div(class = "soft-card", DTOutput("tableAlphaSweep")),
          div(class = "soft-card", DTOutput("tableAlphaSweepBest")),
          div(class = "soft-card", plotOutput("plotAlphaSweepIndices", height = "420px")),
          div(class = "soft-card", plotOutput("plotAlphaSweepNMGI", height = "420px"))
        )
      )
    ),

    tabPanel(
      "Export / Validation",
      fluidRow(
        column(
          6,
          div(class = "soft-card",
              h3(class = "section-head", "Export workbook"),
              p(class = "help-note",
                "Choose the export mode. Validation workbook reproduces the compact baseline-plus-alpha deliverable used for case-study checking. Full studio workbook adds the current scenario trace, portfolio sheets, and alpha-grid best-point sheet."),
              radioButtons(
                "exportMode",
                "Export mode",
                choices = c(
                  "Validation workbook (recommended for model checking)" = "validation",
                  "Full studio workbook" = "full"
                ),
                selected = "validation"
              ),
              actionButton("exportExcel", "Export Excel results", class = "btn-success"),
              br(), br(),
              verbatimTextOutput("exportStatus")
          )
        ),
        column(
          6,
          div(class = "soft-card",
              h3(class = "section-head", "Validation checklist"),
              p(class = "help-note", "Recommended order for validation against the case-study workbook:"),
              tags$ol(
                tags$li("Check alpha-cut intervals"),
                tags$li("Check normalized intervals"),
                tags$li("Check conditional probability intervals and feasibility"),
                tags$li("Check entropy/diversification bounds"),
                tags$li("Check subjective, objective, and integrated weights"),
                tags$li("Check representative probabilities and integrated weights"),
                tags$li("Check final scores, ranking, ICI, and diagnostics")
              )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {

  rv <- reactiveValues(
    baseline_bundle = NULL,
    baseline_alpha = NULL,
    baseline_alpha_info = NULL,
    current_bundle = NULL,
    portfolio = list(),
    portfolio_summary = data.frame(),
    alpha_sweep = NULL,
    import_log = data.frame(),
    export_status = "No export has been produced yet."
  )

  baseline_results <- reactive({
    req(rv$baseline_bundle, rv$baseline_alpha)
    compute_fesmadmii_v2(
      dataMat = rv$baseline_bundle$dataMat,
      deltaMat = rv$baseline_bundle$deltaMat,
      sbjDF = rv$baseline_bundle$sbjDF,
      bcDF = rv$baseline_bundle$bcDF,
      alpha = rv$baseline_alpha
    )
  })

  current_results <- reactive({
    req(rv$current_bundle)
    compute_fesmadmii_v2(
      dataMat = rv$current_bundle$dataMat,
      deltaMat = rv$current_bundle$deltaMat,
      sbjDF = rv$current_bundle$sbjDF,
      bcDF = rv$current_bundle$bcDF,
      alpha = input$scenarioAlpha
    )
  })

  update_portfolio_choices <- function() {
    nm <- names(rv$portfolio)
    updateSelectInput(session, "portfolioRemoveName", choices = nm, selected = if (length(nm) > 0) nm[1] else character(0))
  }

  observeEvent(input$loadData, {
    req(input$fileExcel)

    path <- input$fileExcel$datapath
    sheets <- readxl::excel_sheets(path)

    missingSheets <- c(input$sheetMatrix, input$sheetDelta, input$sheetWeights, input$sheetBC)
    missingSheets <- missingSheets[!missingSheets %in% sheets]
    if (length(missingSheets) > 0) {
      showNotification(paste("Missing required sheets:", paste(missingSheets, collapse = ", ")), type = "error")
      return()
    }

    dfMatrix <- read_sheet_safe(path, input$sheetMatrix)
    dfDelta <- read_sheet_safe(path, input$sheetDelta)
    dfSbj <- read_sheet_safe(path, input$sheetWeights)
    dfBC <- read_sheet_safe(path, input$sheetBC)
    alpha_global_df <- if (input$sheetAlphaGlobal %in% sheets) read_sheet_safe(path, input$sheetAlphaGlobal) else NULL

    if (any(sapply(list(dfMatrix, dfDelta, dfSbj, dfBC), is.null))) {
      showNotification("Failed to read one or more required sheets.", type = "error")
      return()
    }

    alpha_global <- extract_single_alpha(alpha_global_df)
    alpha_final <- if (is.null(alpha_global)) input$fallbackAlpha else alpha_global
    alpha_final <- min(max(alpha_final, 0), 1)

    bundle <- build_model_bundle(dfMatrix, dfDelta, dfSbj, dfBC)

    tryCatch({
      validate_bundle(bundle$dataMat, bundle$deltaMat, bundle$sbjDF, bundle$bcDF, alpha_final)
    }, error = function(e) {
      showNotification(e$message, type = "error")
      stop(e)
    })

    rv$baseline_bundle <- bundle
    rv$current_bundle <- list(
      dataMat = bundle$dataMat,
      deltaMat = bundle$deltaMat,
      sbjDF = bundle$sbjDF,
      bcDF = bundle$bcDF
    )
    rv$baseline_alpha <- alpha_final
    rv$baseline_alpha_info <- data.frame(
      Source = if (is.null(alpha_global)) "Fallback UI alpha" else "GlobalAlpha sheet",
      Value = fmt6(alpha_final),
      stringsAsFactors = FALSE
    )
    rv$portfolio <- list()
    rv$portfolio_summary <- data.frame()
    rv$alpha_sweep <- NULL
    rv$import_log <- data.frame()
    rv$export_status <- "No export has been produced yet."

    critChoices <- rownames(bundle$dataMat)
    altChoices <- colnames(bundle$dataMat)

    updateSelectInput(session, "modWeightCrit", choices = critChoices, selected = critChoices[1])
    updateSelectInput(session, "modDataCrit", choices = critChoices, selected = critChoices[1])
    updateSelectInput(session, "modTypeCrit", choices = critChoices, selected = critChoices[1])
    updateSelectInput(session, "modDataAlt", choices = altChoices, selected = altChoices[1])
    updateNumericInput(session, "scenarioAlpha", value = alpha_final)
    updateTextInput(session, "scenarioName", value = "Scenario_1")
    update_portfolio_choices()

    showNotification(paste("Workbook loaded successfully. Baseline alpha =", fmt6(alpha_final)), type = "message")
  })

  observeEvent(input$modWeightCrit, {
    req(rv$current_bundle)
    row <- rv$current_bundle$sbjDF[rv$current_bundle$sbjDF$Criterion == input$modWeightCrit, , drop = FALSE]
    if (nrow(row) > 0) {
      updateNumericInput(session, "modWeightCentral", value = as.numeric(row$Central[1]))
      updateNumericInput(session, "modWeightDelta", value = as.numeric(row$Delta[1]))
    }
  }, ignoreInit = TRUE)

  observeEvent(input$modDataCrit, {
    req(rv$current_bundle)
    crit <- input$modDataCrit
    alt <- input$modDataAlt %||% colnames(rv$current_bundle$dataMat)[1]
    updateNumericInput(session, "modDataValue", value = as.numeric(rv$current_bundle$dataMat[crit, alt]))
    updateNumericInput(session, "modDataDelta", value = as.numeric(rv$current_bundle$deltaMat[crit, alt]))
  }, ignoreInit = TRUE)

  observeEvent(input$modDataAlt, {
    req(rv$current_bundle)
    crit <- input$modDataCrit %||% rownames(rv$current_bundle$dataMat)[1]
    alt <- input$modDataAlt
    updateNumericInput(session, "modDataValue", value = as.numeric(rv$current_bundle$dataMat[crit, alt]))
    updateNumericInput(session, "modDataDelta", value = as.numeric(rv$current_bundle$deltaMat[crit, alt]))
  }, ignoreInit = TRUE)

  observeEvent(input$modTypeCrit, {
    req(rv$current_bundle)
    row <- rv$current_bundle$bcDF[rv$current_bundle$bcDF$Criterion == input$modTypeCrit, , drop = FALSE]
    if (nrow(row) > 0) {
      updateRadioButtons(session, "modTypeValue", selected = tolower(as.character(row$Type[1])))
    }
  }, ignoreInit = TRUE)

  observeEvent(input$applyWeightChange, {
    req(rv$current_bundle)
    idx <- which(rv$current_bundle$sbjDF$Criterion == input$modWeightCrit)
    if (length(idx) == 1) {
      rv$current_bundle$sbjDF$Central[idx] <- input$modWeightCentral
      rv$current_bundle$sbjDF$Delta[idx] <- pmax(input$modWeightDelta, 0)
      showNotification("Current scenario subjective weight updated.", type = "message")
    }
  })

  observeEvent(input$applyDataChange, {
    req(rv$current_bundle)
    crit <- input$modDataCrit
    alt <- input$modDataAlt
    rv$current_bundle$dataMat[crit, alt] <- input$modDataValue
    rv$current_bundle$deltaMat[crit, alt] <- pmax(input$modDataDelta, 0)
    showNotification("Current scenario data and fuzzy deviation updated.", type = "message")
  })

  observeEvent(input$applyTypeChange, {
    req(rv$current_bundle)
    idx <- which(rv$current_bundle$bcDF$Criterion == input$modTypeCrit)
    if (length(idx) == 1) {
      rv$current_bundle$bcDF$Type[idx] <- input$modTypeValue
      showNotification("Current scenario criterion type updated.", type = "message")
    }
  })

  observeEvent(input$resetScenario, {
    req(rv$baseline_bundle)
    rv$current_bundle <- list(
      dataMat = rv$baseline_bundle$dataMat,
      deltaMat = rv$baseline_bundle$deltaMat,
      sbjDF = rv$baseline_bundle$sbjDF,
      bcDF = rv$baseline_bundle$bcDF
    )
    updateNumericInput(session, "scenarioAlpha", value = rv$baseline_alpha)
    showNotification("Current scenario reset to the baseline state.", type = "warning")
  })

  observeEvent(input$saveScenario, {
    req(rv$current_bundle)
    name <- trimws(input$scenarioName)
    if (!nzchar(name)) {
      showNotification("Please provide a scenario name.", type = "error")
      return()
    }

    res <- current_results()
    rv$portfolio[[name]] <- list(
      bundle = list(
        dataMat = rv$current_bundle$dataMat,
        deltaMat = rv$current_bundle$deltaMat,
        sbjDF = rv$current_bundle$sbjDF,
        bcDF = rv$current_bundle$bcDF
      ),
      alpha = input$scenarioAlpha,
      result = res,
      source = "Manual",
      origin = "Scenario Lab"
    )

    row <- scenario_summary_row(name, res, source = "Manual", origin = "Scenario Lab")
    if (nrow(rv$portfolio_summary) == 0) {
      rv$portfolio_summary <- row
    } else if (name %in% rv$portfolio_summary$Scenario) {
      idx <- which(rv$portfolio_summary$Scenario == name)[1]
      rv$portfolio_summary[idx, ] <- row
    } else {
      rv$portfolio_summary <- rbind(rv$portfolio_summary, row)
    }

    update_portfolio_choices()
    showNotification(paste("Scenario saved to portfolio:", name), type = "message")
  })

  observeEvent(input$importScenarios, {
    req(rv$baseline_bundle)
    req(input$scenarioImportFiles)

    files <- input$scenarioImportFiles
    if (is.null(files) || nrow(files) == 0) {
      showNotification("Please choose at least one scenario workbook.", type = "error")
      return()
    }

    if (isTRUE(input$importClearFirst)) {
      rv$portfolio <- list()
      rv$portfolio_summary <- data.frame()
    }

    import_rows <- list()

    for (i in seq_len(nrow(files))) {
      fpath <- files$datapath[i]
      fname <- files$name[i]
      base_name <- strip_extension(fname)
      target_name <- base_name
      status <- "Imported"
      detail <- ""

      tryCatch({
        loaded <- load_bundle_from_workbook(
          path = fpath,
          sheet_matrix = input$sheetMatrix,
          sheet_delta = input$sheetDelta,
          sheet_weights = input$sheetWeights,
          sheet_bc = input$sheetBC,
          sheet_alpha = input$sheetAlphaGlobal,
          fallback_alpha = input$fallbackAlpha
        )

        validate_scenario_compatibility(rv$baseline_bundle, loaded$bundle)

        existing_names <- names(rv$portfolio)
        if (target_name %in% existing_names) {
          if (identical(input$importDuplicateMode, "rename")) {
            target_name <- make_unique_scenario_name(target_name, existing_names)
            detail <- paste0("Imported with auto-renamed scenario name: ", target_name)
          } else {
            detail <- paste0("Existing scenario overwritten: ", target_name)
          }
        }

        res <- compute_fesmadmii_v2(
          dataMat = loaded$bundle$dataMat,
          deltaMat = loaded$bundle$deltaMat,
          sbjDF = loaded$bundle$sbjDF,
          bcDF = loaded$bundle$bcDF,
          alpha = loaded$alpha
        )

        rv$portfolio[[target_name]] <- list(
          bundle = list(
            dataMat = loaded$bundle$dataMat,
            deltaMat = loaded$bundle$deltaMat,
            sbjDF = loaded$bundle$sbjDF,
            bcDF = loaded$bundle$bcDF
          ),
          alpha = loaded$alpha,
          result = res,
          source = "Imported Excel",
          origin = fname
        )

        row <- scenario_summary_row(target_name, res, source = "Imported Excel", origin = fname)
        if (nrow(rv$portfolio_summary) == 0) {
          rv$portfolio_summary <- row
        } else if (target_name %in% rv$portfolio_summary$Scenario) {
          idx <- which(rv$portfolio_summary$Scenario == target_name)[1]
          rv$portfolio_summary[idx, ] <- row
        } else {
          rv$portfolio_summary <- rbind(rv$portfolio_summary, row)
        }
      }, error = function(e) {
        status <<- "Failed"
        detail <<- e$message
      })

      import_rows[[i]] <- data.frame(
        Workbook = fname,
        Scenario = target_name,
        Status = status,
        Detail = ifelse(nzchar(detail), detail, ifelse(status == "Imported", "Scenario imported successfully.", "")),
        stringsAsFactors = FALSE
      )
    }

    rv$import_log <- do.call(rbind, import_rows)
    update_portfolio_choices()

    n_ok <- sum(rv$import_log$Status == "Imported")
    n_fail <- sum(rv$import_log$Status == "Failed")
    showNotification(sprintf("Scenario import completed: %d imported, %d failed.", n_ok, n_fail), type = if (n_fail > 0) "warning" else "message")
  })

  observeEvent(input$removeScenario, {
    nm <- input$portfolioRemoveName %||% ""
    if (!nzchar(nm) || !(nm %in% names(rv$portfolio))) return()
    rv$portfolio[[nm]] <- NULL
    if (nrow(rv$portfolio_summary) > 0) {
      rv$portfolio_summary <- rv$portfolio_summary[rv$portfolio_summary$Scenario != nm, , drop = FALSE]
    }
    update_portfolio_choices()
    showNotification(paste("Scenario removed:", nm), type = "warning")
  })

  observeEvent(input$clearPortfolio, {
    rv$portfolio <- list()
    rv$portfolio_summary <- data.frame()
    update_portfolio_choices()
    showNotification("Scenario portfolio cleared.", type = "warning")
  })

  observeEvent(input$runAlphaSweep, {
    req(rv$baseline_bundle)
    rv$alpha_sweep <- tryCatch({
      compute_alpha_grid(rv$baseline_bundle, input$alphaFrom, input$alphaTo, input$alphaBy)
    }, error = function(e) {
      showNotification(e$message, type = "error")
      NULL
    })
    if (!is.null(rv$alpha_sweep)) showNotification("Alpha grid completed.", type = "message")
  })

  observeEvent(input$exportExcel, {
    req(rv$baseline_bundle, rv$baseline_alpha)
    mode <- input$exportMode %||% "validation"
    default_name <- if (identical(mode, "full")) "FESMADMIII_full_results.xlsx" else "FESMADMIII_validation_results.xlsx"
    path <- choose_save_path(default_name)
    if (!nzchar(path)) {
      rv$export_status <- "Export cancelled."
      return()
    }

    base <- baseline_results()
    curr <- current_results()

    export_list <- if (identical(mode, "full")) {
      build_full_export_list(
        base_res = base,
        current_res = curr,
        portfolio_summary = rv$portfolio_summary,
        portfolio = rv$portfolio,
        alpha_grid = rv$alpha_sweep,
        baseline_alpha = rv$baseline_alpha,
        current_alpha = input$scenarioAlpha
      )
    } else {
      build_validation_export_list(
        base_res = base,
        alpha_grid = rv$alpha_sweep
      )
    }

    ok <- tryCatch({
      writexl::write_xlsx(export_list, path = path)
      TRUE
    }, error = function(e) {
      showNotification(e$message, type = "error")
      FALSE
    })

    if (ok) {
      rv$export_status <- paste(
        if (identical(mode, "full")) "Full workbook saved to:" else "Validation workbook saved to:",
        normalizePath(path, winslash = "/", mustWork = FALSE)
      )
      showNotification("Excel workbook exported successfully.", type = "message")
    }
  })

  output$overviewSummary <- renderUI({
    if (is.null(rv$baseline_bundle)) {
      div(class = "soft-card",
          h3(class = "section-head", "Current state"),
          p(class = "help-note", "No workbook is loaded yet. Load the case-study workbook or your own workbook to activate the baseline, scenario, portfolio, alpha-grid, and export layers.")
      )
    } else {
      res <- baseline_results()
      div(class = "soft-card",
          h3(class = "section-head", "Baseline snapshot"),
          fluidRow(
            column(3, kpi_box("Global alpha", fmt6(res$alpha), paste(res$M, "criteria x", res$N, "alternatives"), "#1F4E79")),
            column(3, kpi_box("Top alternative", res$alternatives_table$Alternative[1], paste("Score", fmt6(res$alternatives_table$Score[1])), "#A16207")),
            column(3, kpi_box("NMI / CES", paste0(fmt6(res$NMI), " / ", fmt6(res$CES)), "Explanatory power and criteria effectiveness", "#0F766E")),
            column(3, kpi_box("ADI / NMGI", paste0(fmt6(res$ADI), " / ", fmt6(res$NMGI)), "Alternative distinction and global diagnostic performance", "#7C3AED"))
          )
      )
    }
  })

  output$tableMatrixPreview <- renderDT({
    req(rv$baseline_bundle)
    fes_dt(as.data.frame(rv$baseline_bundle$dataMat, stringsAsFactors = FALSE), pageLength = 8, rownames = TRUE) %>%
      formatRound(columns = seq_len(ncol(rv$baseline_bundle$dataMat)), digits = 4)
  })

  output$tableDeltaPreview <- renderDT({
    req(rv$baseline_bundle)
    fes_dt(as.data.frame(rv$baseline_bundle$deltaMat, stringsAsFactors = FALSE), pageLength = 8, rownames = TRUE) %>%
      formatRound(columns = seq_len(ncol(rv$baseline_bundle$deltaMat)), digits = 4)
  })

  output$tableWeightsPreview <- renderDT({
    req(rv$baseline_bundle)
    fes_dt(rv$baseline_bundle$sbjDF, pageLength = 8) %>%
      formatRound(columns = c("Central", "Delta"), digits = 4)
  })

  output$tableBCPreview <- renderDT({
    req(rv$baseline_bundle)
    fes_dt(rv$baseline_bundle$bcDF, pageLength = 8)
  })

  output$tableAlphaInfo <- renderDT({
    req(rv$baseline_alpha_info)
    fes_dt(rv$baseline_alpha_info, pageLength = 5, dom = "t")
  })

  output$baselineKPIs <- renderUI({
    if (is.null(rv$baseline_bundle)) return(NULL)
    res <- baseline_results()
    fluidRow(
      column(3, kpi_box("Global alpha", fmt6(res$alpha), paste(res$M, "criteria x", res$N, "alternatives"), "#1F4E79")),
      column(3, kpi_box("Winner", res$alternatives_table$Alternative[1], paste("Score", fmt6(res$alternatives_table$Score[1])), "#A16207")),
      column(3, kpi_box("NMI / CES", paste0(fmt6(res$NMI), " / ", fmt6(res$CES)), "Core diagnostic pair", "#0F766E")),
      column(3, kpi_box("ADI / NMGI", paste0(fmt6(res$ADI), " / ", fmt6(res$NMGI)), "Distinction and global diagnostic performance", "#7C3AED"))
    )
  })

  output$tableAltBaseline <- renderDT({
    req(rv$baseline_bundle)
    res <- baseline_results()
    fes_dt(res$alternatives_table, pageLength = 10, dom = "t") %>%
      formatRound(columns = "Score", digits = 6)
  })

  output$plotAltBaseline <- renderPlot({
    req(rv$baseline_bundle)
    res <- baseline_results()
    df <- res$alternatives_table
    df$Highlight <- ifelse(df$Rank == 1, "Winner", "Other")

    ggplot(df, aes(x = reorder(Alternative, Score), y = Score, fill = Highlight)) +
      geom_col(width = 0.72) +
      geom_text(aes(label = fmt4(Score)), hjust = -0.08, size = 4.1, fontface = "bold", colour = "#17324D") +
      coord_flip(clip = "off") +
      scale_fill_manual(values = c("Winner" = "#0F4C81", "Other" = "#BFD7EA")) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
      labs(title = "Baseline alternative ranking", subtitle = "Representative final scores", x = NULL, y = "Score") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })

  output$tableCritBaseline <- renderDT({
    req(rv$baseline_bundle)
    res <- baseline_results()
    fes_dt(res$criteria_table, pageLength = 10, scrollX = TRUE) %>%
      formatRound(columns = names(res$criteria_table)[sapply(res$criteria_table, is.numeric)], digits = 6)
  })

  output$plotWeightsBaseline <- renderPlot({
    req(rv$baseline_bundle)
    res <- baseline_results()
    df <- data.frame(
      Criterion = rep(res$critNames, 3),
      WeightType = factor(
        rep(c("Subjective midpoint", "Objective midpoint", "Integrated representative"), each = res$M),
        levels = c("Subjective midpoint", "Objective midpoint", "Integrated representative")
      ),
      Weight = c(
        mid_interval(res$wSBJ_lower, res$wSBJ_upper),
        mid_interval(res$wOBJ_lower, res$wOBJ_upper),
        res$wINT_star
      ),
      stringsAsFactors = FALSE
    )

    ggplot(df, aes(x = Criterion, y = Weight, fill = WeightType)) +
      geom_col(position = position_dodge(width = 0.72), width = 0.64) +
      scale_fill_manual(values = c("Subjective midpoint" = "#94A3B8", "Objective midpoint" = "#38B2AC", "Integrated representative" = "#1F4E79")) +
      labs(title = "Weight architecture", subtitle = "Subjective, objective, and final integrated criterion importance", x = NULL, y = "Weight") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })

  output$plotICIBaseline <- renderPlot({
    req(rv$baseline_bundle)
    res <- baseline_results()
    df <- res$criteria_table[order(res$criteria_table$ICI), c("Criterion", "ICI")]
    ggplot(df, aes(x = ICI, y = reorder(Criterion, ICI))) +
      geom_segment(aes(x = 0, xend = ICI, y = Criterion, yend = Criterion), colour = "#C8D6E5", linewidth = 1.2) +
      geom_point(size = 4.2, colour = "#7C3AED") +
      geom_text(aes(label = fmt4(ICI)), nudge_x = 0.01, hjust = 0, size = 4, colour = "#17324D", fontface = "bold") +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
      labs(title = "Integrated Criteria Importance (ICI)", subtitle = "Operational criterion contribution", x = "ICI", y = NULL) +
      theme_minimal(base_size = 13)
  })

  output$tableDiagBaseline <- renderDT({
    req(rv$baseline_bundle)
    res <- baseline_results()
    fes_dt(res$diagnostics_table, pageLength = 10, dom = "t") %>%
      formatRound(columns = "Value", digits = 6)
  })

  output$plotDiagBaseline <- renderPlot({
    req(rv$baseline_bundle)
    res <- baseline_results()
    df <- data.frame(
      Index = factor(c("NMI", "CES", "ADI", "NMGI"), levels = c("NMI", "CES", "ADI", "NMGI")),
      Value = c(res$NMI, res$CES, res$ADI, res$NMGI),
      stringsAsFactors = FALSE
    )

    ggplot(df, aes(x = Value, y = Index, fill = Index)) +
      geom_col(width = 0.64) +
      geom_text(aes(label = fmt4(Value)), hjust = -0.08, size = 4.1, fontface = "bold", colour = "#17324D") +
      scale_fill_manual(values = c("NMI" = "#1F4E79", "CES" = "#0F766E", "ADI" = "#A16207", "NMGI" = "#7C3AED")) +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
      labs(title = "Final diagnostic indices", subtitle = "NMI, CES, ADI, and NMGI", x = "Value", y = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")
  })

  output$tableAlphaCutsBaseline <- renderDT({
    req(rv$baseline_bundle)
    res <- baseline_results()
    df <- interval_matrix_to_df(res$xi_lower, res$xi_upper, row_label = "Criterion")
    fes_dt(df, pageLength = 10, scrollX = TRUE) %>%
      formatRound(columns = names(df)[-1], digits = 6)
  })

  output$tableNormBaseline <- renderDT({
    req(rv$baseline_bundle)
    res <- baseline_results()
    df <- interval_matrix_to_df(res$r_lower, res$r_upper, row_label = "Criterion")
    fes_dt(df, pageLength = 10, scrollX = TRUE) %>%
      formatRound(columns = names(df)[-1], digits = 6)
  })

  output$tableProbBaseline <- renderDT({
    req(rv$baseline_bundle)
    res <- baseline_results()
    df <- interval_matrix_to_df(res$P_lower, res$P_upper, row_label = "Criterion")
    fes_dt(df, pageLength = 10, scrollX = TRUE) %>%
      formatRound(columns = names(df)[-1], digits = 6)
  })

  output$tableEntropyBaseline <- renderDT({
    req(rv$baseline_bundle)
    res <- baseline_results()
    fes_dt(res$entropy_table, pageLength = 10, dom = "t") %>%
      formatRound(columns = c("h_L", "h_U", "d_L", "d_U"), digits = 6)
  })

  output$tableFeasibilityBaseline <- renderDT({
    req(rv$baseline_bundle)
    res <- baseline_results()
    fes_dt(res$feasibility_table, pageLength = 10, dom = "t") %>%
      formatRound(columns = c("Sum_P_Lower", "Sum_P_Upper"), digits = 6)
  })

  output$scenarioKPIs <- renderUI({
    if (is.null(rv$current_bundle) || is.null(rv$baseline_bundle)) return(NULL)
    base <- baseline_results()
    curr <- current_results()
    fluidRow(
      column(3, kpi_box("Scenario alpha", fmt6(curr$alpha), paste("Baseline", fmt6(base$alpha)), "#1F4E79")),
      column(3, kpi_box("Scenario winner", curr$alternatives_table$Alternative[1], paste("Score", fmt6(curr$alternatives_table$Score[1])), "#A16207")),
      column(3, kpi_box("Top-score shift", sprintf("%+.6f", curr$alternatives_table$Score[1] - base$alternatives_table$Score[1]), "Scenario minus baseline", "#0F766E")),
      column(3, kpi_box("NMGI shift", sprintf("%+.6f", curr$NMGI - base$NMGI), "Scenario minus baseline", "#7C3AED"))
    )
  })

  output$tableAltScenario <- renderDT({
    req(rv$current_bundle)
    res <- current_results()
    fes_dt(res$alternatives_table, pageLength = 10, dom = "t") %>%
      formatRound(columns = "Score", digits = 6)
  })

  output$tableDiagScenario <- renderDT({
    req(rv$current_bundle)
    res <- current_results()
    fes_dt(res$diagnostics_table, pageLength = 10, dom = "t") %>%
      formatRound(columns = "Value", digits = 6)
  })

  output$plotAltScenarioCompare <- renderPlot({
    req(rv$current_bundle, rv$baseline_bundle)
    base <- baseline_results()
    curr <- current_results()

    df <- merge(
      data.frame(Alternative = base$altNames, Baseline = base$scores, stringsAsFactors = FALSE),
      data.frame(Alternative = curr$altNames, Scenario = curr$scores, stringsAsFactors = FALSE),
      by = "Alternative",
      all = TRUE
    )
    df$Alternative <- factor(df$Alternative, levels = df$Alternative[order(df$Scenario)])

    ggplot(df, aes(y = Alternative)) +
      geom_segment(aes(x = Baseline, xend = Scenario, yend = Alternative), colour = "#CBD5E1", linewidth = 1.6) +
      geom_point(aes(x = Baseline), colour = "#94A3B8", size = 3.8) +
      geom_point(aes(x = Scenario), colour = "#0F4C81", size = 4.4) +
      labs(title = "Baseline vs current scenario scores", subtitle = "Grey = baseline, blue = current scenario", x = "Score", y = NULL) +
      theme_minimal(base_size = 13)
  })

  output$tableCompScenario <- renderDT({
    req(rv$current_bundle, rv$baseline_bundle)
    base <- baseline_results()
    curr <- current_results()

    df <- data.frame(
      Metric = c("NMI", "CES", "ADI", "NMGI", "Top score"),
      Baseline = c(base$NMI, base$CES, base$ADI, base$NMGI, base$alternatives_table$Score[1]),
      CurrentScenario = c(curr$NMI, curr$CES, curr$ADI, curr$NMGI, curr$alternatives_table$Score[1]),
      Delta = c(curr$NMI - base$NMI, curr$CES - base$CES, curr$ADI - base$ADI, curr$NMGI - base$NMGI, curr$alternatives_table$Score[1] - base$alternatives_table$Score[1]),
      stringsAsFactors = FALSE
    )

    fes_dt(df, pageLength = 10, dom = "t") %>%
      formatRound(columns = c("Baseline", "CurrentScenario", "Delta"), digits = 6)
  })

  output$plotCompScenario <- renderPlot({
    req(rv$current_bundle, rv$baseline_bundle)
    base <- baseline_results()
    curr <- current_results()

    df <- data.frame(
      Metric = factor(c("NMI", "CES", "ADI", "NMGI"), levels = c("NMI", "CES", "ADI", "NMGI")),
      Baseline = c(base$NMI, base$CES, base$ADI, base$NMGI),
      CurrentScenario = c(curr$NMI, curr$CES, curr$ADI, curr$NMGI),
      stringsAsFactors = FALSE
    )

    ggplot(df, aes(y = Metric)) +
      geom_segment(aes(x = Baseline, xend = CurrentScenario, yend = Metric), colour = "#D7E3EF", linewidth = 1.6) +
      geom_point(aes(x = Baseline), colour = "#94A3B8", size = 3.8) +
      geom_point(aes(x = CurrentScenario), colour = "#0E7490", size = 4.4) +
      labs(title = "Diagnostic movement under the current scenario", subtitle = "Baseline vs current scenario", x = "Value", y = NULL) +
      theme_minimal(base_size = 13)
  })

  output$tableTraceScenario <- renderDT({
    req(rv$current_bundle)
    curr <- current_results()
    trace_df <- merge(
      curr$criteria_table[, c("Criterion", "wINT_star", "ICI")],
      data.frame(Criterion = curr$critNames, ShiftApplied = curr$shifts, stringsAsFactors = FALSE),
      by = "Criterion",
      all.x = TRUE
    )
    fes_dt(trace_df, pageLength = 10, dom = "t") %>%
      formatRound(columns = c("wINT_star", "ICI", "ShiftApplied"), digits = 6)
  })

  output$scenarioImportKPIs <- renderUI({
    if (is.null(rv$import_log) || nrow(rv$import_log) == 0) {
      return(div(class = "soft-card",
                 h3(class = "section-head", "Import status"),
                 p(class = "help-note", "No scenario-import batch has been executed yet. Upload one or more scenario workbooks and import them into the shared portfolio.")))
    }
    df <- rv$import_log
    n_ok <- sum(df$Status == "Imported")
    n_fail <- sum(df$Status == "Failed")
    fluidRow(
      column(4, kpi_box("Imported in last batch", n_ok, "Successfully added to portfolio", "#1F4E79")),
      column(4, kpi_box("Failed in last batch", n_fail, "Validation or workbook incompatibility", "#A16207")),
      column(4, kpi_box("Portfolio size", length(rv$portfolio), "Manual + imported scenarios combined", "#0F766E"))
    )
  })

  output$tableScenarioImportLog <- renderDT({
    if (is.null(rv$import_log) || nrow(rv$import_log) == 0) {
      return(fes_dt(data.frame(Message = "No import log available yet."), pageLength = 5, dom = "t"))
    }
    fes_dt(rv$import_log, pageLength = 10, dom = "t")
  })

  output$portfolioKPIs <- renderUI({
    if (nrow(rv$portfolio_summary) == 0) {
      return(div(class = "soft-card",
                 h3(class = "section-head", "Portfolio status"),
                 p(class = "help-note", "No scenarios are saved yet. Use the Scenario Lab or the Scenario Import tab to add scenarios into the shared portfolio.")
      ))
    }

    df <- rv$portfolio_summary
    best_idx <- which.max(df$NMGI)[1]
    fluidRow(
      column(4, kpi_box("Saved scenarios", nrow(df), "Multi-scenario comparison is active", "#1F4E79")),
      column(4, kpi_box("Best NMGI scenario", df$Scenario[best_idx], paste("NMGI", fmt6(df$NMGI[best_idx])), "#7C3AED")),
      column(4, kpi_box("Best winner", df$Winner[best_idx], paste("Score", fmt6(df$WinnerScore[best_idx])), "#A16207"))
    )
  })

  output$tablePortfolio <- renderDT({
    req(nrow(rv$portfolio_summary) > 0)
    fes_dt(rv$portfolio_summary, pageLength = 10, scrollX = TRUE) %>%
      formatRound(columns = c("Alpha", "WinnerScore", "NMI", "CES", "ADI", "NMGI"), digits = 6)
  })

  output$plotPortfolioNMGI <- renderPlot({
    req(nrow(rv$portfolio_summary) > 0)
    df <- rv$portfolio_summary
    df$Scenario <- factor(df$Scenario, levels = df$Scenario[order(df$NMGI)])

    ggplot(df, aes(x = Scenario, y = NMGI)) +
      geom_col(fill = "#7C3AED", width = 0.68) +
      geom_text(aes(label = fmt4(NMGI)), vjust = -0.35, size = 4, fontface = "bold", colour = "#17324D") +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
      labs(title = "Scenario portfolio: NMGI comparison", subtitle = "Higher values indicate stronger overall diagnostic coherence", x = NULL, y = "NMGI") +
      theme_minimal(base_size = 13)
  })

  output$plotPortfolioTopScore <- renderPlot({
    req(nrow(rv$portfolio_summary) > 0)
    df <- rv$portfolio_summary
    df$Scenario <- factor(df$Scenario, levels = df$Scenario[order(df$WinnerScore)])

    ggplot(df, aes(x = WinnerScore, y = Scenario, fill = Winner)) +
      geom_col(width = 0.68) +
      geom_text(aes(label = paste0(Winner, " | ", fmt4(WinnerScore))), hjust = -0.04, size = 4, fontface = "bold", colour = "#17324D") +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
      labs(title = "Scenario portfolio: top-score comparison", subtitle = "Winning alternative and winning score by saved scenario", x = "Top score", y = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })

  output$plotPortfolioHeatmap <- renderPlot({
    req(nrow(rv$portfolio_summary) > 0)
    df <- rv$portfolio_summary
    long_df <- rbind(
      data.frame(Scenario = df$Scenario, Metric = "NMI", Value = df$NMI),
      data.frame(Scenario = df$Scenario, Metric = "CES", Value = df$CES),
      data.frame(Scenario = df$Scenario, Metric = "ADI", Value = df$ADI),
      data.frame(Scenario = df$Scenario, Metric = "NMGI", Value = df$NMGI)
    )

    ggplot(long_df, aes(x = Scenario, y = Metric, fill = Value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = fmt4(Value)), size = 4, fontface = "bold", colour = "#17324D") +
      scale_fill_gradient(low = "#E8F1FB", high = "#1F4E79") +
      labs(title = "Scenario portfolio diagnostic heatmap", subtitle = "Compact view of the four final indices across saved scenarios", x = NULL, y = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "right")
  })

  output$alphaGridKPIs <- renderUI({
    if (is.null(rv$alpha_sweep) || nrow(rv$alpha_sweep) == 0) return(NULL)
    df <- rv$alpha_sweep
    best_idx <- which.max(df$NMGI)[1]
    step_val <- ifelse(nrow(df) > 1, diff(df$Alpha)[1], 0)

    fluidRow(
      column(3, kpi_box("Grid points", nrow(df), paste("Step", fmt6(step_val)), "#1F4E79")),
      column(3, kpi_box("Peak NMGI alpha", fmt6(df$Alpha[best_idx]), paste("NMGI", fmt6(df$NMGI[best_idx])), "#7C3AED")),
      column(3, kpi_box("Winner at peak", df$Winner[best_idx], paste("Score", fmt6(df$WinnerScore[best_idx])), "#A16207")),
      column(3, kpi_box("Alpha range", paste0(fmt6(min(df$Alpha)), " to ", fmt6(max(df$Alpha))), "Sensitivity analysis across the confidence grid", "#0F766E"))
    )
  })

  output$tableAlphaSweep <- renderDT({
    req(!is.null(rv$alpha_sweep), nrow(rv$alpha_sweep) > 0)
    fes_dt(rv$alpha_sweep, pageLength = 10) %>%
      formatRound(columns = c("Alpha", "WinnerScore", "NMI", "CES", "ADI", "NMGI"), digits = 6)
  })

  output$tableAlphaSweepBest <- renderDT({
    req(!is.null(rv$alpha_sweep), nrow(rv$alpha_sweep) > 0)
    best <- rv$alpha_sweep[which.max(rv$alpha_sweep$NMGI), , drop = FALSE]
    fes_dt(best, pageLength = 5, dom = "t") %>%
      formatRound(columns = c("Alpha", "WinnerScore", "NMI", "CES", "ADI", "NMGI"), digits = 6)
  })

  output$plotAlphaSweepIndices <- renderPlot({
    req(!is.null(rv$alpha_sweep), nrow(rv$alpha_sweep) > 0)
    df <- rv$alpha_sweep
    plot_df <- rbind(
      data.frame(Alpha = df$Alpha, Index = "NMI", Value = df$NMI),
      data.frame(Alpha = df$Alpha, Index = "CES", Value = df$CES),
      data.frame(Alpha = df$Alpha, Index = "ADI", Value = df$ADI),
      data.frame(Alpha = df$Alpha, Index = "NMGI", Value = df$NMGI)
    )

    ggplot(plot_df, aes(x = Alpha, y = Value, color = Index)) +
      geom_line(linewidth = 1.15) +
      geom_point(size = 2.8) +
      scale_color_manual(values = c("NMI" = "#1F4E79", "CES" = "#0E7490", "ADI" = "#A16207", "NMGI" = "#7C3AED")) +
      labs(title = "Alpha-grid diagnostics", subtitle = "Evolution of the final indices across alpha", x = "Alpha", y = "Value") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })

  output$plotAlphaSweepNMGI <- renderPlot({
    req(!is.null(rv$alpha_sweep), nrow(rv$alpha_sweep) > 0)
    df <- rv$alpha_sweep
    best_idx <- which.max(df$NMGI)[1]

    ggplot(df, aes(x = Alpha, y = NMGI)) +
      geom_line(linewidth = 1.15, colour = "#7C3AED") +
      geom_point(size = 3, colour = "#7C3AED") +
      geom_point(data = df[best_idx, , drop = FALSE], colour = "#F59E0B", size = 4.8) +
      geom_text(data = df[best_idx, , drop = FALSE], aes(label = paste0("Peak ", fmt4(NMGI))), nudge_y = 0.0006, colour = "#17324D", fontface = "bold", size = 4) +
      labs(title = "Alpha-grid NMGI profile", subtitle = "Peak NMGI alpha highlighted", x = "Alpha", y = "NMGI") +
      theme_minimal(base_size = 13)
  })

  output$exportStatus <- renderText({
    rv$export_status
  })
}


if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2 && identical(args[1], "--validate")) {
    input_path <- args[2]
    output_path <- if (length(args) >= 3) args[3] else "FESMADMIII_validation_results.xlsx"
    alpha_from <- if (length(args) >= 4) as.numeric(args[4]) else NA_real_
    alpha_to   <- if (length(args) >= 5) as.numeric(args[5]) else NA_real_
    alpha_by   <- if (length(args) >= 6) as.numeric(args[6]) else NA_real_
    if (all(is.finite(c(alpha_from, alpha_to, alpha_by)))) {
      write_validation_workbook_from_input(input_path, output_path, alpha_from, alpha_to, alpha_by)
    } else {
      write_validation_workbook_from_input(input_path, output_path)
    }
    quit(save = "no")
  }
}

app <- shinyApp(ui = ui, server = server)

if (interactive()) {
  shiny::runApp(app, launch.browser = launch_standalone_window)
}
