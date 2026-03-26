########################################################################
# Analysis and Improvement of Space-Filling Design Expansion
# via Block Design
#
# Added Baseline Comparisons: Morris & Mitchell (lhs) and Sliced LHD (SLHD)
########################################################################

set.seed(2026)

# ==============================================================
# PART 1: Core Implementation of Algorithm 1
# ==============================================================

left_cyclic_latin_square <- function(k) {
  L <- matrix(0L, nrow = k, ncol = k)
  for (i in 1:k) {
    for (j in 1:k) {
      L[i, j] <- ((i - 1) + (j - 1)) %% k + 1  
    }
  }
  L
}

algorithm1 <- function(A, block_vectors, k) {
  N <- nrow(A)
  n <- ncol(A)
  s <- length(block_vectors[[1]])
  
  intermediates <- vector("list", k)
  for (i in 1:k) {
    Ci <- matrix(0L, nrow = N, ncol = n)
    bv <- block_vectors[[i]]
    for (r in 1:N) {
      for (cc in 1:n) {
        Ci[r, cc] <- bv[A[r, cc]]
      }
    }
    intermediates[[i]] <- Ci
  }
  
  L <- left_cyclic_latin_square(k)
  
  rows_list <- vector("list", k)
  for (l in 1:k) {
    perm <- L[l, ]
    row_block <- do.call(cbind, intermediates[perm])
    rows_list[[l]] <- row_block
  }
  
  DA <- do.call(rbind, rows_list)
  DA
}

# ==============================================================
# PART 2: Distance and Efficiency Metrics
# ==============================================================

L1_distance_design <- function(D) {
  N <- nrow(D)
  min_dist <- Inf
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      d <- sum(abs(D[i, ] - D[j, ]))
      if (d < min_dist) min_dist <- d
    }
  }
  as.integer(min_dist)
}

L1_upper_bound <- function(N, s, n) {
  denom <- 3 * N * s - 3 * s
  if (denom == 0) return(1L)
  as.integer(floor(N * (s^2 - 1) * n / denom))
}

column_correlation <- function(D) {
  n <- ncol(D)
  corrs <- c()
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      di <- as.numeric(D[, i])
      dj <- as.numeric(D[, j])
      di_c <- di - mean(di)
      dj_c <- dj - mean(dj)
      denom <- sqrt(sum(di_c^2) * sum(dj_c^2))
      if (denom < 1e-12) {
        corrs <- c(corrs, 0)
      } else {
        corrs <- c(corrs, abs(sum(di_c * dj_c) / denom))
      }
    }
  }
  if (length(corrs) == 0) return(c(rho_max = 0, rho_avg = 0))
  c(rho_max = max(corrs), rho_avg = mean(corrs))
}

# ==============================================================
# PART 3: Block Vector Construction Methods
# ==============================================================

block_vectors_row_order <- function(s, k) {
  lapply(0:(k - 1), function(i) ((i * s + 1):((i + 1) * s)))
}

block_vectors_column_order <- function(s, k) {
  vecs <- vector("list", k)
  for (i in 1:k) vecs[[i]] <- integer(0) 
  for (val in 1:(k * s)) {
    idx <- ((val - 1) %% k) + 1
    vecs[[idx]] <- c(vecs[[idx]], val)
  }
  vecs
}

block_vectors_zigzag <- function(s, k = 3) {
  b <- vector("list", k)
  for (i in 1:k) b[[i]] <- integer(0) 
  elements <- 1:(k * s)
  for (i in seq_along(elements)) {
    group <- (i - 1) %/% k  
    pos <- ((i - 1) %% k) + 1  
    if (group %% 2 == 0) {
      b[[pos]] <- c(b[[pos]], elements[i])
    } else {
      b[[k + 1 - pos]] <- c(b[[k + 1 - pos]], elements[i])
    }
  }
  b
}

# ==============================================================
# PART 4: Competitor Baselines (lhs & SLHD)
# ==============================================================

generate_maximin_lhs <- function(N_total, n_factors) {
  if (!requireNamespace("lhs", quietly = TRUE)) {
    return(NULL)
  }
  # Generate continuous LHS in [0,1] and convert to exact integer ranks 
  # to match the U-type design criteria of our expanded matrix
  D_norm <- lhs::maximinLHS(N_total, n_factors)
  D_int <- apply(D_norm, 2, rank) 
  return(D_int)
}

generate_maximin_slhd <- function(t_slices, m_runs, n_factors) {
  if (!requireNamespace("SLHD", quietly = TRUE)) {
    return(NULL)
  }
  res <- SLHD::maximinSLHD(t = t_slices, m = m_runs, k = n_factors)
  return(res$Design)
}

# ==============================================================
# PART 5: SA-TBO Optimization (Stages 1 & 2)
# ==============================================================

sa_optimize_blocks <- function(s, A, k = 3, max_iter = 10000) {
  best_blocks <- block_vectors_zigzag(s, k)
  best_DA <- algorithm1(A, best_blocks, k)
  best_score <- L1_distance_design(best_DA)
  
  current_blocks <- lapply(best_blocks, function(b) b)
  current_score <- best_score
  Temp <- 2.0
  alpha <- 0.9995
  
  for (it in 1:max_iter) {
    bi <- sample(k, 1)
    bj <- sample(setdiff(1:k, bi), 1)
    pi_idx <- sample(s, 1)
    pj_idx <- sample(s, 1)
    
    new_blocks <- lapply(current_blocks, function(b) b)
    tmp <- new_blocks[[bi]][pi_idx]
    new_blocks[[bi]][pi_idx] <- new_blocks[[bj]][pj_idx]
    new_blocks[[bj]][pj_idx] <- tmp
    
    new_DA <- algorithm1(A, new_blocks, k)
    new_score <- L1_distance_design(new_DA)
    
    delta <- new_score - current_score
    if (delta > 0 || runif(1) < exp(delta / max(Temp, 1e-10))) {
      current_blocks <- new_blocks
      current_score <- new_score
      if (current_score > best_score) {
        best_blocks <- lapply(current_blocks, function(b) b)
        best_score <- current_score
      }
    }
    Temp <- Temp * alpha
  }
  list(blocks = best_blocks, score = best_score)
}

two_stage_refine <- function(A, s, k = 3, max_iter_sa = 8000, max_iter_refine = 8000) {
  sa_result <- sa_optimize_blocks(s, A, k, max_iter = max_iter_sa)
  DA <- algorithm1(A, sa_result$blocks, k)
  
  kN <- nrow(DA)
  kn <- ncol(DA)
  best_score <- L1_distance_design(DA)
  best_DA <- DA
  current_DA <- DA
  current_score <- best_score
  Temp <- 1.0
  alpha <- 0.9993
  
  for (it in 1:max_iter_refine) {
    col <- sample(kn, 1)
    r1 <- sample(kN, 1)
    r2 <- sample(setdiff(1:kN, r1), 1)
    
    new_DA <- current_DA
    tmp <- new_DA[r1, col]
    new_DA[r1, col] <- new_DA[r2, col]
    new_DA[r2, col] <- tmp
    
    new_score <- L1_distance_design(new_DA)
    delta <- new_score - current_score
    if (delta > 0 || runif(1) < exp(delta / max(Temp, 1e-10))) {
      current_DA <- new_DA
      current_score <- new_score
      if (current_score > best_score) {
        best_DA <- current_DA
        best_score <- current_score
      }
    }
    Temp <- Temp * alpha
  }
  list(DA = best_DA, score = best_score, blocks = sa_result$blocks)
}

# ==============================================================
# PART 6: Initial Designs & Execution
# ==============================================================

get_initial_design <- function(s, n) {
  if (s == 3 && n == 3) {
    return(matrix(c(3,2,3, 2,1,1, 1,3,2), nrow = 3, byrow = TRUE))
  } else if (s == 4 && n == 3) {
    return(matrix(c(4,2,1, 1,3,2, 2,1,4, 3,4,3), nrow = 4, byrow = TRUE))
  } else if (s == 4 && n == 4) {
    return(matrix(c(1,4,2,3, 2,1,4,1, 3,3,1,4, 4,2,3,2), nrow = 4, byrow = TRUE))
  } else if (s == 5 && n == 5) {
    return(matrix(c(1,4,3,5,4, 2,2,1,1,3, 3,1,4,4,1, 4,3,5,2,5, 5,5,2,3,2),
                  nrow = 5, byrow = TRUE))
  } else {
    D <- matrix(0L, nrow = s, ncol = n)
    for (j in 1:n) D[, j] <- sample(s)
    return(D)
  }
}

cat(strrep("=", 70), "\n")
cat("SPACE-FILLING DESIGN EXPANSION: SA-TBO VS EXTERNAL BASELINES\n")
cat(strrep("=", 70), "\n")

test_cases <- list(
  list(s = 3, n = 3),
  list(s = 4, n = 3),
  list(s = 4, n = 4),
  list(s = 5, n = 5)
)
k <- 3

all_results <- data.frame(
  Case     = character(),
  Method   = character(),
  d1       = integer(),
  UB       = integer(),
  Eff_pct  = numeric(),
  rho_max  = numeric(),
  rho_avg  = numeric(),
  stringsAsFactors = FALSE
)

for (tc in test_cases) {
  s <- tc$s
  n <- tc$n
  label <- sprintf("s=%d,n=%d", s, n)
  
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("Case: %s, k=%d\n", label, k))
  
  A <- get_initial_design(s, n)
  ks <- k * s
  kn <- k * n
  ub <- L1_upper_bound(ks, ks, kn)
  
  add_result <- function(method_name, DA) {
    d1 <- L1_distance_design(DA)
    eff <- d1 / ub * 100
    cc <- column_correlation(DA)
    cat(sprintf("  %-25s d1=%5d, eff=%6.1f%%, rho_max=%.4f, rho_avg=%.4f\n",
                paste0(method_name, ":"), d1, eff, cc["rho_max"], cc["rho_avg"]))
    all_results <<- rbind(all_results, data.frame(
      Case = label, Method = method_name, d1 = d1, UB = ub,
      Eff_pct = round(eff, 1),
      rho_max = round(cc["rho_max"], 4),
      rho_avg = round(cc["rho_avg"], 4),
      stringsAsFactors = FALSE
    ))
  }
  
  # 1. Internal Baselines
  bv <- block_vectors_row_order(s, k)
  add_result("Row Order (Internal)", algorithm1(A, bv, k))
  
  bv <- block_vectors_zigzag(s, k)
  add_result("Zigzag (Internal)", algorithm1(A, bv, k))
  
  # 2. External Baselines (Morris & Mitchell & SLHD)
  cat("  [Running External Baselines...]\n")
  D_mm <- generate_maximin_lhs(N_total = ks, n_factors = kn)
  if (!is.null(D_mm)) add_result("Morris & Mitchell (lhs)", D_mm)
  
  D_slhd <- generate_maximin_slhd(t_slices = k, m_runs = s, n_factors = kn)
  if (!is.null(D_slhd)) add_result("Maximin SLHD (SLHD)", D_slhd)
  
  # 3. Proposed SA-TBO
  cat("  [Running SA-TBO Optimization...]\n")
  ts_res <- two_stage_refine(A, s, k, max_iter_sa = 10000, max_iter_refine = 10000)
  add_result("SA-TBO Framework", ts_res$DA)
}

cat("\nDone. Saving to CSV...\n")
write.csv(all_results, "Expanded_Design_Comparisons.csv", row.names = FALSE)