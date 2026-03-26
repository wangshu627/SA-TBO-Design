########################################################################
# Chapter 7 Case Study: Surrogate Modeling of the Borehole Function
# Repeated Experiments (N trials) for Statistical Significance
########################################################################

library(DiceKriging)

# ==============================================================
# MODULE 1 & 2 & 3: Core SA-TBO Functions
# ==============================================================
left_cyclic_latin_square <- function(k) {
  L <- matrix(0L, nrow = k, ncol = k)
  for (i in 1:k) {
    for (j in 1:k) { L[i, j] <- ((i - 1) + (j - 1)) %% k + 1 }
  }
  L
}

algorithm1 <- function(A, block_vectors, k) {
  N <- nrow(A); n <- ncol(A); s <- length(block_vectors[[1]])
  intermediates <- vector("list", k)
  for (i in 1:k) {
    Ci <- matrix(0L, nrow = N, ncol = n)
    bv <- block_vectors[[i]]
    for (r in 1:N) {
      for (cc in 1:n) { Ci[r, cc] <- bv[A[r, cc]] }
    }
    intermediates[[i]] <- Ci
  }
  L <- left_cyclic_latin_square(k)
  rows_list <- vector("list", k)
  for (l in 1:k) {
    perm <- L[l, ]
    rows_list[[l]] <- do.call(cbind, intermediates[perm])
  }
  do.call(rbind, rows_list)
}

L1_distance_detailed <- function(D) {
  d_matrix <- as.matrix(dist(D, method = "manhattan"))
  diag(d_matrix) <- Inf
  min_d <- min(d_matrix)
  idx <- which(d_matrix == min_d, arr.ind = TRUE)
  list(score = as.integer(min_d), r1 = idx[1, 1], r2 = idx[1, 2])
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
    if (group %% 2 == 0) { b[[pos]] <- c(b[[pos]], elements[i])
    } else { b[[k + 1 - pos]] <- c(b[[k + 1 - pos]], elements[i]) }
  }
  b
}

optimize_initial_design <- function(s, n, max_iter = 5000) {
  A <- matrix(0L, nrow = s, ncol = n)
  for (j in 1:n) A[, j] <- sample(s)
  best_A <- A
  best_score <- L1_distance_detailed(A)$score
  current_A <- A; current_score <- best_score
  Temp <- 2.0; alpha <- 0.999
  
  for (it in 1:max_iter) {
    new_A <- current_A
    col <- sample(n, 1)
    r1 <- sample(s, 1); r2 <- sample(setdiff(1:s, r1), 1)
    tmp <- new_A[r1, col]; new_A[r1, col] <- new_A[r2, col]; new_A[r2, col] <- tmp
    new_score <- L1_distance_detailed(new_A)$score
    delta <- new_score - current_score
    if (delta > 0 || runif(1) < exp(delta / max(Temp, 1e-10))) {
      current_A <- new_A; current_score <- new_score
      if (current_score > best_score) { best_A <- current_A; best_score <- current_score }
    }
    Temp <- Temp * alpha
  }
  return(best_A)
}

targeted_deep_refine <- function(A, s, k = 3, max_iter = 15000) {
  blocks <- block_vectors_zigzag(s, k)
  DA <- algorithm1(A, blocks, k)
  kN <- nrow(DA); kn <- ncol(DA)
  info <- L1_distance_detailed(DA)
  best_score <- info$score; best_DA <- DA
  current_DA <- DA; current_score <- best_score
  Temp <- 1.5; alpha <- 0.9995
  
  for (it in 1:max_iter) {
    new_DA <- current_DA
    if (runif(1) < 0.5) {
      r_target <- ifelse(runif(1) < 0.5, info$r1, info$r2)
      r_swap <- sample(setdiff(1:kN, r_target), 1)
    } else {
      r_target <- sample(kN, 1)
      r_swap <- sample(setdiff(1:kN, r_target), 1)
    }
    col <- sample(kn, 1)
    tmp <- new_DA[r_target, col]; new_DA[r_target, col] <- new_DA[r_swap, col]; new_DA[r_swap, col] <- tmp
    
    new_info <- L1_distance_detailed(new_DA)
    delta <- new_info$score - current_score
    
    if (delta > 0 || (delta == 0 && runif(1) < 0.3) || runif(1) < exp(delta / max(Temp, 1e-10))) {
      current_DA <- new_DA; current_score <- new_info$score; info <- new_info 
      if (current_score > best_score) { best_DA <- current_DA; best_score <- current_score }
    }
    Temp <- Temp * alpha
  }
  list(DA = best_DA, score = best_score)
}

# ==============================================================
# MODULE 4: Borehole Function Definition
# ==============================================================
borehole_func <- function(x) {
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  
  rw <- x[, 1] * (0.15 - 0.05) + 0.05
  r  <- x[, 2] * (50000 - 100) + 100
  Tu <- x[, 3] * (115600 - 63070) + 63070
  Hu <- x[, 4] * (1110 - 990) + 990
  Tl <- x[, 5] * (116 - 63.1) + 63.1
  Hl <- x[, 6] * (820 - 700) + 700
  L  <- x[, 7] * (1580 - 1120) + 1120
  Kw <- x[, 8] * (12045 - 9855) + 9855
  
  numerator <- 2 * pi * Tu * (Hu - Hl)
  term1 <- log(r / rw)
  term2 <- 1 + (2 * L * Tu) / (term1 * (rw^2) * Kw) + (Tu / Tl)
  denominator <- term1 * term2
  
  return(numerator / denominator)
}

# ==============================================================
# MODULE 5 & 6: N Independent Trials Loop
# ==============================================================
N_trials <- 30  # 30独立试验保证统计显著性
s <- 15; n <- 3; k <- 3
levels_total <- s * k

results_all <- data.frame(Trial = integer(), Method = character(), 
                          RMSE = numeric(), MAE = numeric(), stringsAsFactors = FALSE)

# 用于存储第一次试验的数据，以便后续生成散点图和投影图
plot_data <- list()

cat(sprintf("Starting %d Independent Trials...\n", N_trials))

for(trial in 1:N_trials) {
  set.seed(2026 + trial) # 每次循环改变随机种子
  cat(sprintf("\n=== Running Trial %d / %d ===\n", trial, N_trials))
  
  # 1. Generate Initial Base Design
  A_init <- optimize_initial_design(s, n, max_iter = 8000)
  
  # SCENARIO 1: Initial 3D
  X_init_01 <- (A_init - 0.5) / s
  X_init_full <- cbind(X_init_01, matrix(0.5, nrow = s, ncol = 5)) 
  y_init <- borehole_func(X_init_full)
  model_init <- try(km(formula = ~1, design = X_init_01, response = y_init, 
                       covtype = "matern5_2", nugget.estim = TRUE, control = list(trace = FALSE)), silent=TRUE)
  
  # SCENARIO 2: Baseline Expanded Design
  D_base_raw <- algorithm1(A_init, block_vectors_column_order(s, k), k)
  X_base_01 <- (D_base_raw[, 1:8] - 0.5) / levels_total 
  y_base <- borehole_func(X_base_01)
  model_base <- try(km(formula = ~1, design = X_base_01, response = y_base, 
                       covtype = "matern5_2", nugget.estim = TRUE, control = list(trace = FALSE)), silent=TRUE)
  
  # SCENARIO 3: SA-TBO Optimized Expanded Design
  tbo_res <- targeted_deep_refine(A_init, s, k, max_iter = 25000)
  D_opt_raw <- tbo_res$DA
  X_opt_01 <- (D_opt_raw[, 1:8] - 0.5) / levels_total 
  y_opt <- borehole_func(X_opt_01)
  model_opt <- try(km(formula = ~1, design = X_opt_01, response = y_opt, 
                      covtype = "matern5_2", nugget.estim = TRUE, control = list(trace = FALSE)), silent=TRUE)
  
  # Testing
  n_test <- 10000
  X_test_01 <- matrix(runif(n_test * 8), nrow = n_test, ncol = 8)
  y_test_true <- borehole_func(X_test_01)
  
  # Predictions
  pred_init <- predict(model_init, newdata = X_test_01[, 1:3], type = "UK")$mean
  pred_base <- predict(model_base, newdata = X_test_01, type = "UK")$mean
  pred_opt  <- predict(model_opt, newdata = X_test_01, type = "UK")$mean
  
  # Calculate Errors
  rmse_init <- sqrt(mean((y_test_true - pred_init)^2))
  mae_init  <- mean(abs(y_test_true - pred_init))
  
  rmse_base <- sqrt(mean((y_test_true - pred_base)^2))
  mae_base  <- mean(abs(y_test_true - pred_base))
  
  rmse_opt <- sqrt(mean((y_test_true - pred_opt)^2))
  mae_opt  <- mean(abs(y_test_true - pred_opt))
  
  # 保存当前 Trial 的结果
  results_all <- rbind(results_all, data.frame(
    Trial = trial, Method = "Initial Small Design (3D)", RMSE = rmse_init, MAE = mae_init
  ))
  results_all <- rbind(results_all, data.frame(
    Trial = trial, Method = "Baseline Expansion (8D)", RMSE = rmse_base, MAE = mae_base
  ))
  results_all <- rbind(results_all, data.frame(
    Trial = trial, Method = "SA-TBO Expansion (8D)", RMSE = rmse_opt, MAE = mae_opt
  ))
  
  # 记录第一次循环的数据用于画图
  if(trial == 1) {
    plot_data$y_test_true <- y_test_true; plot_data$pred_init <- pred_init
    plot_data$pred_base <- pred_base; plot_data$pred_opt <- pred_opt
    plot_data$X_base_01 <- X_base_01; plot_data$X_opt_01 <- X_opt_01
  }
}

# ==============================================================
# MODULE 7: Statistical Aggregation & High-Quality Plotting
# ==============================================================
cat("\nAggregating Results over 30 Trials...\n")
summary_stats <- aggregate(cbind(RMSE, MAE) ~ Method, data = results_all, 
                           FUN = function(x) c(Mean = mean(x), SD = sd(x)))

summary_df <- data.frame(
  Method = summary_stats$Method,
  RMSE_Mean = round(summary_stats$RMSE[, "Mean"], 3),
  RMSE_SD   = round(summary_stats$RMSE[, "SD"], 3),
  MAE_Mean  = round(summary_stats$MAE[, "Mean"], 3),
  MAE_SD    = round(summary_stats$MAE[, "SD"], 3)
)
print(summary_df)
write.csv(summary_df, "Borehole_30Trials_Summary.csv", row.names = FALSE)

# --------- Plots ---------
cat("Generating PDF Figures for LaTeX...\n")

# Reorder factor levels for logical plotting order
results_all$Method <- factor(results_all$Method, 
                             levels = c("Initial Small Design (3D)", 
                                        "Baseline Expansion (8D)", 
                                        "SA-TBO Expansion (8D)"))

# 【新增】Figure: Boxplot of RMSE over 30 trials
pdf("fig_boxplot_rmse5.pdf", width = 8, height = 6)
par(mar = c(5, 5, 4, 2) + 0.1, cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.4)
boxplot(RMSE ~ Method, data = results_all, 
        col = c("#E0E0E0", "#FF9999", "#99CCFF"),
        pch = 16, outcol = rgb(0,0,0,0.5),
        ylab = "Root Mean Square Error (RMSE)", 
        xlab = "",
        main = "Distribution of Predictive RMSE over 30 Independent Trials")
grid(ny = NULL, nx = NA, lty = 2, col = "gray80")
boxplot(RMSE ~ Method, data = results_all, 
        col = c("#E0E0E0", "#FF9999", "#99CCFF"),
        pch = 16, outcol = rgb(0,0,0,0.5), add = TRUE)
dev.off()

# Original Fig 1: Initial Dilemma
pdf("fig_initial_dilemma5.pdf", width = 6, height = 5)
par(mar = c(4.5, 4.5, 3, 1), cex.main = 1.2, cex.lab = 1.1)
plot(plot_data$y_test_true, plot_data$pred_init, col = rgb(0.5, 0.5, 0.5, 0.2), pch = 16, 
     xlab = "True Borehole Yield", ylab = "Predicted Yield", 
     main = "Initial Dilemma: 3D Model Underfitting",
     xlim = range(plot_data$y_test_true), ylim = range(plot_data$y_test_true))
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
dev.off()

# Original Fig 2: True vs Predicted Scatter
pdf("fig_true_vs_predicted5.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), cex.main = 1.2, cex.lab = 1.1)
plot(plot_data$y_test_true, plot_data$pred_base, col = rgb(0.8, 0.2, 0.2, 0.3), pch = 16, 
     xlab = "True Borehole Yield", ylab = "Predicted Yield", 
     main = "Baseline Expansion Design (8D)",
     xlim = range(plot_data$y_test_true), ylim = range(plot_data$y_test_true))
abline(a = 0, b = 1, col = "black", lwd = 2, lty = 2)
plot(plot_data$y_test_true, plot_data$pred_opt, col = rgb(0.1, 0.5, 0.8, 0.3), pch = 16, 
     xlab = "True Borehole Yield", ylab = "Predicted Yield", 
     main = "SA-TBO Expansion Design (8D)",
     xlim = range(plot_data$y_test_true), ylim = range(plot_data$y_test_true))
abline(a = 0, b = 1, col = "black", lwd = 2, lty = 2)
dev.off()

# Original Fig 3: 2D Projection
pdf("fig_2d_projection5.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), cex.main = 1.2, cex.lab = 1.1)
plot(plot_data$X_base_01[, 1], plot_data$X_base_01[, 4], pch = 21, bg = "red", cex = 1.5,
     xlab = "Borehole Radius (rw)", ylab = "Upper Aquifer Head (Hu)",
     main = "Baseline: Severe Spatial Collapsing", xlim = c(0, 1), ylim = c(0, 1))
grid()
plot(plot_data$X_opt_01[, 1], plot_data$X_opt_01[, 4], pch = 21, bg = "blue", cex = 1.5,
     xlab = "Borehole Radius (rw)", ylab = "Upper Aquifer Head (Hu)",
     main = "SA-TBO: Comprehensive Coverage", xlim = c(0, 1), ylim = c(0, 1))
grid()
dev.off()

cat("All done! PDF images and summary CSV saved to current working directory.\n")