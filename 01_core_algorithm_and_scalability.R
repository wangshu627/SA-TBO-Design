########################################################################
# Analysis and Improvement of Space-Filling Design Expansion
# via Block Design
#
# This script:
# 1. Implements Algorithm 1 from the paper
# 2. Computes space-filling efficiency (L1 distance / upper bound)
# 3. Implements SA-optimized block vectors (Improvement 1)
# 4. Implements two-stage construction with refinement (Improvement 2)
# 5. Exhaustive search for small cases
# 6. Generates comparison tables
########################################################################

set.seed(42)

# ==============================================================
# PART 1: Core Implementation of Algorithm 1
# ==============================================================

left_cyclic_latin_square <- function(k) {
  L <- matrix(0L, nrow = k, ncol = k)
  for (i in 1:k) {
    for (j in 1:k) {
      L[i, j] <- ((i - 1) + (j - 1)) %% k + 1  # 1-indexed
    }
  }
  L
}

algorithm1 <- function(A, block_vectors, k) {
  # A: initial design (N x n), values in {1,...,s}
  # block_vectors: list of k vectors, each of length s
  # k: expansion factor
  # Returns: expanded design DA (kN x kn)
  N <- nrow(A)
  n <- ncol(A)
  s <- length(block_vectors[[1]])
  
  # Step 3: Generate intermediate designs C_i
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
  
  # Step 4: Left-cyclic arrangement
  L <- left_cyclic_latin_square(k)
  
  # Step 5: Construct DA
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

# L1_distance_design <- function(D) {
#   N <- nrow(D)
#   min_dist <- Inf
#   for (i in 1:(N - 1)) {
#     for (j in (i + 1):N) {
#       d <- sum(abs(D[i, ] - D[j, ]))
#       if (d < min_dist) min_dist <- d
#     }
#   }
#   as.integer(min_dist)
# }

# 改进后的距离计算（使用内置矩阵运算）
L1_distance_fast <- function(D) {
  # dist() 函数在 C 层面运行，远快于手写循环
  d_matrix <- as.matrix(dist(D, method = "manhattan"))
  # 排除对角线的 0，取最小值
  diag(d_matrix) <- Inf
  return(as.integer(min(d_matrix)))
}

L2_distance_design <- function(D) {
  N <- nrow(D)
  min_dist <- Inf
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      d <- sum((as.numeric(D[i, ]) - as.numeric(D[j, ]))^2)
      if (d < min_dist) min_dist <- d
    }
  }
  min_dist
}

L1_upper_bound <- function(N, s, n) {
  denom <- 3 * N * s - 3 * s
  if (denom == 0) return(1L)
  as.integer(floor(N * (s^2 - 1) * n / denom))
}

L2_upper_bound <- function(N, s, n) {
  denom <- 6 * N * s - 6 * s
  if (denom == 0) return(1L)
  as.integer(floor(N * (s^2 - 1) * n / denom))
}

column_correlation <- function(D) {
  N <- nrow(D)
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

projection_metric <- function(D) {
  N <- nrow(D)
  n <- ncol(D)
  total <- 0
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      prod_val <- 1
      for (l in 1:n) {
        diff2 <- (as.numeric(D[i, l]) - as.numeric(D[j, l]))^2
        if (diff2 == 0) return(Inf)
        prod_val <- prod_val * diff2
      }
      total <- total + 1 / prod_val
    }
  }
  C_N_2 <- N * (N - 1) / 2
  (total / C_N_2)^(1 / n)
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
    group <- (i - 1) %/% k  # 0-indexed group
    pos <- ((i - 1) %% k) + 1  # 1-indexed position
    if (group %% 2 == 0) {
      b[[pos]] <- c(b[[pos]], elements[i])
    } else {
      b[[k + 1 - pos]] <- c(b[[k + 1 - pos]], elements[i])
    }
  }
  b
}

# ==============================================================
# PART 4: IMPROVEMENT 1 - SA Block Vector Optimization
# ==============================================================

sa_optimize_blocks <- function(s, A, k = 3, criterion = "L1",
                               max_iter = 10000) {
  # Start from zigzag
  best_blocks <- block_vectors_zigzag(s, k)
  best_DA <- algorithm1(A, best_blocks, k)
  best_score <- L1_distance_design(best_DA)
  
  current_blocks <- lapply(best_blocks, function(b) b)
  current_score <- best_score
  Temp <- 2.0
  alpha <- 0.9995
  
  for (it in 1:max_iter) {
    # Propose swap between two blocks
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

# ==============================================================
# PART 5: IMPROVEMENT 2 - Two-Stage Construction
# ==============================================================

two_stage_refine <- function(A, s, k = 3,
                             max_iter_sa = 8000,
                             max_iter_refine = 8000) {
  # Stage 1: SA block optimization
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
  
  # Stage 2: within-column swaps
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
# PART 6: Exhaustive Search (small cases)
# ==============================================================

exhaustive_block_search <- function(s, A, k = 3) {
  ks <- k * s
  elements <- 1:ks
  best_d1 <- 0
  best_blocks <- NULL
  count <- 0
  
  b1_combos <- combn(elements, s, simplify = FALSE)
  
  for (b1 in b1_combos) {
    remaining <- setdiff(elements, b1)
    b2_combos <- combn(remaining, s, simplify = FALSE)
    
    for (b2 in b2_combos) {
      b3 <- setdiff(remaining, b2)
      blocks <- list(b1, b2, b3)
      DA <- algorithm1(A, blocks, k)
      d1 <- L1_distance_design(DA)
      
      if (d1 > best_d1) {
        best_d1 <- d1
        best_blocks <- blocks
      }
      count <- count + 1
    }
  }
  
  list(blocks = best_blocks, d1 = best_d1, count = count)
}

# ==============================================================
# PART 7: Initial Designs
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
  } else if (s == 6 && n == 6) {
    return(matrix(c(6,1,4,3,5,2, 1,5,2,6,3,4, 4,3,6,1,2,5,
                    3,6,1,5,4,3, 5,2,5,4,6,1, 2,4,3,2,1,6),
                  nrow = 6, byrow = TRUE))
  } else {
    D <- matrix(0L, nrow = s, ncol = n)
    for (j in 1:n) D[, j] <- sample(s)
    return(D)
  }
}

# # ==============================================================
# # PART 8: Run All Experiments
# # ==============================================================
# 
# cat(strrep("=", 70), "\n")
# cat("SPACE-FILLING DESIGN EXPANSION: ANALYSIS AND IMPROVEMENT\n")
# cat(strrep("=", 70), "\n")
# 
# test_cases <- list(
#   list(s = 3, n = 3),
#   list(s = 4, n = 3),
#   list(s = 4, n = 4),
#   list(s = 5, n = 5)
# )
# k <- 3
# 
# all_results <- data.frame(
#   Case     = character(),
#   Method   = character(),
#   d1       = integer(),
#   UB       = integer(),
#   Eff_pct  = numeric(),
#   rho_max  = numeric(),
#   rho_avg  = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# for (tc in test_cases) {
#   s <- tc$s
#   n <- tc$n
#   label <- sprintf("s=%d,n=%d", s, n)
#   
#   cat(sprintf("\n%s\n", strrep("=", 60)))
#   cat(sprintf("Case: %s, k=%d\n", label, k))
#   cat(sprintf("LH(%d,%d) -> LH(%d,%d)\n", s, n, k * s, k * n))
#   cat(sprintf("%s\n", strrep("=", 60)))
#   
#   A <- get_initial_design(s, n)
#   ks <- k * s
#   kn <- k * n
#   ub <- L1_upper_bound(ks, ks, kn)
#   d1_A <- L1_distance_design(A)
#   
#   cat(sprintf("Initial design L1 distance: %d\n", d1_A))
#   cat(sprintf("L1 upper bound for expanded: %d\n", ub))
#   
#   add_result <- function(method_name, DA) {
#     d1 <- L1_distance_design(DA)
#     eff <- d1 / ub * 100
#     cc <- column_correlation(DA)
#     cat(sprintf("  %-22s d1=%5d, eff=%6.1f%%, rho_max=%.4f, rho_avg=%.4f\n",
#                 paste0(method_name, ":"), d1, eff, cc["rho_max"], cc["rho_avg"]))
#     all_results <<- rbind(all_results, data.frame(
#       Case = label, Method = method_name, d1 = d1, UB = ub,
#       Eff_pct = round(eff, 1),
#       rho_max = round(cc["rho_max"], 4),
#       rho_avg = round(cc["rho_avg"], 4),
#       stringsAsFactors = FALSE
#     ))
#   }
#   
#   # Method 1: Row order
#   bv <- block_vectors_row_order(s, k)
#   DA <- algorithm1(A, bv, k)
#   add_result("Row Order (Paper)", DA)
#   
#   # Method 2: Column order
#   bv <- block_vectors_column_order(s, k)
#   DA <- algorithm1(A, bv, k)
#   add_result("Column Order (Paper)", DA)
#   
#   # Method 3: Zigzag
#   bv <- block_vectors_zigzag(s, k)
#   DA <- algorithm1(A, bv, k)
#   add_result("Zigzag Interleave", DA)
#   
#   # Method 4: SA-Optimized
#   cat("  [Running SA optimization...]\n")
#   t0 <- proc.time()
#   sa_res <- sa_optimize_blocks(s, A, k, max_iter = 12000)
#   DA <- algorithm1(A, sa_res$blocks, k)
#   t_sa <- (proc.time() - t0)["elapsed"]
#   add_result("SA-Optimized Blocks", DA)
#   cat(sprintf("    Optimal blocks: b1=(%s), b2=(%s), b3=(%s)\n",
#               paste(sa_res$blocks[[1]], collapse = ","),
#               paste(sa_res$blocks[[2]], collapse = ","),
#               paste(sa_res$blocks[[3]], collapse = ",")))
#   cat(sprintf("    Time: %.1f s\n", t_sa))
#   
#   # Method 5: Two-Stage
#   cat("  [Running two-stage refinement...]\n")
#   t0 <- proc.time()
#   ts_res <- two_stage_refine(A, s, k,
#                              max_iter_sa = 10000,
#                              max_iter_refine = 10000)
#   t_2s <- (proc.time() - t0)["elapsed"]
#   add_result("Two-Stage Refined", ts_res$DA)
#   cat(sprintf("    Time: %.1f s\n", t_2s))
# }
# 
# # Exhaustive search for s=3
# cat(sprintf("\n%s\n", strrep("=", 60)))
# cat("EXHAUSTIVE SEARCH: s=3, n=3\n")
# cat(sprintf("%s\n", strrep("=", 60)))
# A3 <- get_initial_design(3, 3)
# ex_res <- exhaustive_block_search(3, A3)
# ub3 <- L1_upper_bound(9, 9, 9)
# cat(sprintf("Checked %d partitions\n", ex_res$count))
# cat(sprintf("Optimal d1=%d, UB=%d, Eff=%.1f%%\n",
#             ex_res$d1, ub3, ex_res$d1 / ub3 * 100))
# cat(sprintf("Optimal blocks: (%s), (%s), (%s)\n",
#             paste(ex_res$blocks[[1]], collapse = ","),
#             paste(ex_res$blocks[[2]], collapse = ","),
#             paste(ex_res$blocks[[3]], collapse = ",")))
# 
# # Final summary
# cat(sprintf("\n\n%s\n", strrep("=", 70)))
# cat("FINAL SUMMARY TABLE\n")
# cat(sprintf("%s\n", strrep("=", 70)))
# print(all_results, row.names = FALSE)
# 
# cat("\nDone.\n")




# ==============================================================
# BATCH EXECUTION FOR SCALABILITY TABLE (Export to CSV)
# ==============================================================

# Define the large-scale test cases
large_test_cases <- list(
  list(s = 7, n = 5, k = 3),
  list(s = 7, n = 7, k = 3),
  list(s = 10, n = 5, k = 3),
  list(s = 10, n = 10, k = 3),
  list(s = 3, n = 3, k = 4),
  list(s = 3, n = 3, k = 5)
)

results_df <- data.frame(
  Case = character(),
  k = integer(),
  Baseline_d1 = integer(),
  SA_d1 = integer(),
  Theoretical_UB = integer(),
  SA_Efficiency_pct = numeric(),
  SA_Time_sec = numeric(),
  stringsAsFactors = FALSE
)

cat("Starting large-scale scalability tests...\n")

for (tc in large_test_cases) {
  s <- tc$s
  n <- tc$n
  k <- tc$k
  case_label <- sprintf("s=%d, n=%d", s, n)
  cat(sprintf("Processing %s with k=%d...\n", case_label, k))
  
  # Generate initial random design (or specific if available)
  A <- matrix(0L, nrow = s, ncol = n)
  for (j in 1:n) A[, j] <- sample(s)
  
  # Calculate upper bound
  ks <- k * s
  kn <- k * n
  ub <- L1_upper_bound(ks, ks, kn)
  
  # 1. Baseline (Row Order)
  bv_baseline <- block_vectors_row_order(s, k)
  DA_baseline <- algorithm1(A, bv_baseline, k)
  d1_baseline <- L1_distance_design(DA_baseline)
  
  # 2. SA-Optimized
  t0 <- proc.time()
  # Use fewer max_iter for s=10 to keep time reasonable, e.g., 5000
  iter_count <- ifelse(s >= 10, 10000, 20000) 
  sa_res <- sa_optimize_blocks(s, A, k, max_iter = iter_count)
  DA_sa <- algorithm1(A, sa_res$blocks, k)
  t_sa <- (proc.time() - t0)["elapsed"]
  
  d1_sa <- L1_distance_design(DA_sa)
  eff_pct <- round((d1_sa / ub) * 100, 2)
  
  # Append to dataframe
  results_df <- rbind(results_df, data.frame(
    Case = case_label,
    k = k,
    Baseline_d1 = d1_baseline,
    SA_d1 = d1_sa,
    Theoretical_UB = ub,
    SA_Efficiency_pct = eff_pct,
    SA_Time_sec = round(t_sa, 2)
  ))
}

# Export to CSV
write.csv(results_df, "scalability_results.csv", row.names = FALSE)
cat("Done! Results saved to 'scalability_results.csv'\n")
print(results_df)