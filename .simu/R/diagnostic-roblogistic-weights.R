# Run from the repository root:
# Rscript .simu/R/diagnostic-roblogistic-weights.R
# No estimator, optimizer, or model-fitting function is called.
# Outputs: .simu/figures/roblogistic_weights.png,
#          .simu/figures/roblogistic_score_norm.png, and
#          .simu/data/roblogistic_weights.csv
# Rcpp and a C++ compiler are used only for the small validation below.
#
# Implementation correspondence (line numbers refer to the original source):
# src/roblogistic.cpp:815-816 adds an intercept to the supplied predictors.
# Lines 827-834: eta = x1 * start; mu = Sigmoid(eta);
# xp = ||x1_i||_2; s = |y_i - mu_i| * xp; z_i = (y_i-mu_i)*wc(s,c).
# Thus the individual log-likelihood score is U_i = x1_i*(y_i-mu_i),
# and Tukey receives s = ||U_i||_2, without variance/information scaling.
# src/common.h:55: wc(s,c) = 1 - 2*s^2/c^2 + s^4/c^4 for |s|<=c,
# and zero otherwise. Equivalently (1-(s/c)^2)^2 inside the cutoff.
# c = 4.685061 is the default in roblogisticWmle1 (line 797).
# The same multiplier occurs in roblogisticWmle1b, roblogisticWmle and
# wmle_logistic::f_grad. Wmle additionally subtracts Ewymu from z.
# The returned field "weights" is the Fisher-scoring diagonal w, NOT wc.
# Mqle is a DIFFERENT procedure: lines 595-597 pass the Pearson residual
# r=(y-mu)/sqrt(mu*(1-mu)) to psi(r,c). This script diagnoses score-based WMLE.

stopifnot(file.exists("src/roblogistic.cpp"), file.exists("src/common.h"))
set.seed(20260922)
n <- 200L
beta <- c(-0.5, 1, -0.5)
cutoff <- 4.685061
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n))
p <- plogis(drop(X %*% beta))
Y <- rbinom(n, size = 1, prob = p)
# Select an ordinary point using design alone: nearest to the predictor origin.
i0 <- which.min(rowSums(X[, 2:3]^2))
delta <- seq(-40, 40, by = 0.01)
path <- X[rep(i0, length(delta)), , drop = FALSE]
path[, 2] <- X[i0, 2] + delta
eta <- drop(path %*% beta)
design_norm <- sqrt(rowSums(path^2))

tukey_weight <- function(s, cutoff) {
  ans <- numeric(length(s))
  inside <- abs(s) <= cutoff
  ans[inside] <- (1 - (s[inside] / cutoff)^2)^2
  ans
}

evaluate_path <- function(y) {
  # Evaluate 1-p as plogis(-eta), avoiding cancellation when p rounds to 1.
  # This is algebraically y-plogis(eta), with stable agreeing-response tails.
  residual <- if (y == 1) plogis(-eta) else -plogis(eta)
  score <- path * residual
  s <- abs(residual) * design_norm
  data.frame(i0 = i0, simulated_y = Y[i0], y = y,
             actual_response = y == Y[i0], cutoff = cutoff,
             beta0 = beta[1], beta1 = beta[2], beta2 = beta[3],
             delta = delta, modified_covariate = path[, 2], x2 = path[, 3],
             eta = eta, probability = plogis(eta), residual = residual,
             score_intercept = score[, 1], score_x1 = score[, 2],
             score_x2 = score[, 3], design_norm = design_norm,
             score_norm = s, tukey_argument = s,
             robust_weight = tukey_weight(s, cutoff))
}
values <- rbind(evaluate_path(0L), evaluate_path(1L))

# wc has no R export. Compile its EXACT definition from common.h together
# with the original Sigmoid struct, using a temporary Rcpp wrapper only.
# The wrapper evaluates the score argument at fixed beta, reproducing
# roblogistic.cpp:827-834; it cannot fit or update coefficients.
header <- readLines("src/common.h", warn = FALSE)
wc_source <- grep("^inline double wc\\(", header, value = TRUE)
sigmoid_start <- grep("^struct Sigmoid \\{", header)
stopifnot(length(wc_source) == 1L, length(sigmoid_start) == 1L)
sigmoid_end <- sigmoid_start + which(header[seq.int(sigmoid_start + 1L,
                                                       length(header))] == "};")[1]
Rcpp::cppFunction(
  includes = paste(c("#include <cmath>", wc_source,
                     header[sigmoid_start:sigmoid_end]), collapse = "\n"),
  code = 'Rcpp::NumericMatrix fixed_weight_cpp(Rcpp::NumericMatrix x,
              Rcpp::NumericVector beta, Rcpp::NumericVector y, double c) {
    Rcpp::NumericMatrix out(x.nrow(), 3);
    Sigmoid sigmoid;
    for (int i = 0; i < x.nrow(); ++i) {
      double eta = 0.0, norm2 = 0.0;
      for (int j = 0; j < x.ncol(); ++j) {
        eta += x(i,j) * beta[j];
        norm2 += x(i,j) * x(i,j);
      }
      double mu = sigmoid(eta);
      double s = std::abs(y[i] - mu) * std::sqrt(norm2);
      out(i,0) = mu; out(i,1) = s; out(i,2) = wc(s,c);
    }
    return out;
  }',
  cacheDir = tempdir()
)
check_delta <- c(-40, -20, -8, -4, 0, 4, 8, 20, 40)
checked <- values[values$delta %in% check_delta, ]
cpp <- fixed_weight_cpp(cbind(1, checked$modified_covariate, checked$x2),
                        beta, checked$y, cutoff)
errors <- c(probability = max(abs(cpp[, 1] - checked$probability)),
            score_norm = max(abs(cpp[, 2] - checked$tukey_argument)),
            weight = max(abs(cpp[, 3] - checked$robust_weight)))
# Independently reconstruct ||x_i(y_i-p_i)||_2 from every saved score
# component. The intercept contribution is score_intercept.
reconstructed_score_norm <- sqrt(checked$score_intercept^2 +
                                 checked$score_x1^2 + checked$score_x2^2)
stopifnot(max(abs(reconstructed_score_norm - checked$score_norm)) < 1e-12,
          identical(checked$score_norm, checked$tukey_argument))
# Algebraic equivalence, not bitwise equality: polynomial ordering and the
# stable residual cause floating-point differences in the last few bits.
stopifnot(all(errors < 1e-12))
# Check the actual helper at and around its cutoff and at score zero.
boundary_s <- c(0, cutoff * (1 - 1e-8), cutoff, cutoff * (1 + 1e-8))
Rcpp::cppFunction(
  includes = paste("#include <cmath>", wc_source, sep = "\n"),
  code = 'double tukey_cpp(double s, double c) { return wc(s,c); }',
  cacheDir = tempdir()
)
stopifnot(max(abs(vapply(boundary_s, tukey_cpp, numeric(1), c = cutoff) -
                  tukey_weight(boundary_s, cutoff))) < 1e-12)
cat("C++ validation maximum absolute errors (both responses):\n")
print(errors)

dir.create(".simu/data", recursive = TRUE, showWarnings = FALSE)
dir.create(".simu/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(values, ".simu/data/roblogistic_weights.csv", row.names = FALSE)
png(".simu/figures/roblogistic_weights.png", width = 1500, height = 900,
    res = 150)
par(mar = c(4.5, 4.5, 4, 1))
colors <- c("#0072B2", "#D55E00")
plot(range(delta), c(0, 1), type = "n", ylim = c(-0.03, 1.04),
     xlab = expression(delta~~"(change in the first Gaussian predictor)"),
     ylab = "Robust observation weight", main = "Tukey weight at fixed beta")
abline(h = c(0, 1), col = "grey75", lty = 3)
abline(v = 0, col = "grey75", lty = 3)
for (y in 0:1) {
  v <- values[values$y == y, ]
  lines(v$delta, v$robust_weight, col = colors[y + 1],
        lwd = if (y == Y[i0]) 3 else 2,
        lty = if (y == Y[i0]) 1 else 2)
}
baseline <- values[values$actual_response & values$delta == 0, ]
points(0, baseline$robust_weight, pch = 19, col = colors[Y[i0] + 1])
legend("right", inset = 0.04,
       legend = paste0("y = ", 0:1,
                       ifelse(0:1 == Y[i0], " (simulated)", " (counterfactual)")),
       col = colors, lty = ifelse(0:1 == Y[i0], 1, 2),
       lwd = ifelse(0:1 == Y[i0], 3, 2), bg = "white", bty = "o")
mtext(sprintf("n = %d; i0 = %d; beta = (-0.5, 1, -0.5); c = %.6f; original x1 = %.3f",
              n, i0, cutoff, X[i0, 2]), side = 3, line = 0.4, cex = 0.8)
invisible(dev.off())

# Plot the exact input s = ||x_i(y_i-p_i)||_2 supplied to wc(). These values
# are the score_norm/tukey_argument columns already used by the weight plot.
png(".simu/figures/roblogistic_score_norm.png", width = 1500, height = 900,
    res = 150)
par(mar = c(4.5, 4.5, 4, 1))
plot(range(delta), range(values$score_norm), type = "n",
     xlab = expression(delta~~"(change in the first Gaussian predictor)"),
     ylab = "Euclidean norm of likelihood score ||U_i||_2",
     main = "Likelihood-score norm at fixed beta")
abline(h = c(0, cutoff), col = "grey75", lty = 3)
abline(v = 0, col = "grey75", lty = 3)
for (y in 0:1) {
  v <- values[values$y == y, ]
  lines(v$delta, v$score_norm, col = colors[y + 1],
        lwd = if (y == Y[i0]) 3 else 2,
        lty = if (y == Y[i0]) 1 else 2)
}
points(0, baseline$score_norm, pch = 19, col = colors[Y[i0] + 1])
legend("top", inset = 0.04,
       legend = paste0("y = ", 0:1,
                       ifelse(0:1 == Y[i0], " (simulated)", " (counterfactual)")),
       col = colors, lty = ifelse(0:1 == Y[i0], 1, 2),
       lwd = ifelse(0:1 == Y[i0], 3, 2), bg = "white", bty = "o")
mtext(sprintf("n = %d; i0 = %d; beta = (-0.5, 1, -0.5); c = %.6f; original x1 = %.3f",
              n, i0, cutoff, X[i0, 2]), side = 3, line = 0.4, cex = 0.8)
invisible(dev.off())
cat(sprintf("Selected i0 = %d, simulated y = %d, original x1 = %.6f, x2 = %.6f\n",
            i0, Y[i0], X[i0, 2], X[i0, 3]))
print(values[values$delta %in% c(-40, 0, 40),
             c("y", "delta", "eta", "tukey_argument", "robust_weight")],
      row.names = FALSE)
