# Companion calculations for "Multiperiod Attribution: Who Needs Linking Functions?"

CASE <- "base"  # "base" or "pure_selection"
START_WEALTH <- 100
TOL <- 1e-12
SECTORS <- list(c(1, 2), c(3, 4))
q <- function(x) round(x, 10)

W_BENCH <- matrix(c(
  0.20, 0.10, 0.10, 0.60,
  0.20, 0.20, 0.40, 0.20,
  0.15, 0.05, 0.10, 0.70,
  0.10, 0.10, 0.20, 0.60
), nrow = 4, byrow = TRUE)

W_PORT_BASE <- matrix(c(
  0.05, 0.05, 0.20, 0.70,
  0.05, 0.05, 0.40, 0.50,
  0.20, 0.10, 0.20, 0.50,
  0.20, 0.20, 0.20, 0.40
), nrow = 4, byrow = TRUE)

W_PORT_PURE_SEL <- matrix(c(
  0.10, 0.00, 0.70, 0.20,
  0.10, 0.00, 0.75, 0.15,
  0.30, 0.00, 0.30, 0.40,
  0.40, 0.00, 0.25, 0.35
), nrow = 4, byrow = TRUE)

W_BENCH_PURE_SEL <- W_PORT_BASE

RETURNS <- matrix(c(
   0.05,  0.05,  0.05, 0.10,
  -0.10, -0.05,  0.02, 0.08,
  -0.05,  0.06, -0.02, 0.03,
   0.05,  0.10,  0.04, 0.08
), nrow = 4, byrow = TRUE)

period_stats <- function(wp, wb, r) {
  ap <- sapply(SECTORS, function(s) sum(wp[s]))
  ab <- sapply(SECTORS, function(s) sum(wb[s]))
  rp <- sapply(seq_along(SECTORS), function(j) sum(wp[SECTORS[[j]]] * r[SECTORS[[j]]]) / ap[j])
  rb <- sapply(seq_along(SECTORS), function(j) sum(wb[SECTORS[[j]]] * r[SECTORS[[j]]]) / ab[j])

  RP <- sum(ap * rp)
  RB <- sum(ab * rb)
  A <- sum((ap - ab) * (rb - RB))
  S <- sum(ab * (rp - rb))
  I <- sum((ap - ab) * (rp - rb))
  RPB <- sum(ap * rb)
  RBP <- sum(ab * rp)
  list(RP = RP, RB = RB, g = c(A, S, I), RPB = RPB, RBP = RBP)
}

run <- function(w_port, w_bench) {
  n <- nrow(RETURNS)
  stats <- lapply(seq_len(n), function(t) period_stats(w_port[t, ], w_bench[t, ], RETURNS[t, ]))
  RP <- sapply(stats, `[[`, "RP")
  RB <- sapply(stats, `[[`, "RB")
  RA <- RP - RB
  G <- t(sapply(stats, `[[`, "g"))

  stopifnot(all(abs(rowSums(G) - RA) < TOL))

  WP <- WB <- START_WEALTH
  WA <- 0
  ewa <- c(0, 0, 0)
  carry <- 0

  cat("EWA account across time (per 100 initial):\n")
  cat(sprintf("%2s %9s %9s %9s %9s %9s %9s\n", "t", "A", "S", "I", "Carry", "sum", "dW^A"))
  for (t in seq_len(n)) {
    E <- WP * G[t, ]
    C <- WA * RB[t]
    dWA <- WP * RP[t] - WB * RB[t]
    stopifnot(abs(sum(E) + C - dWA) < TOL)

    ewa <- ewa + E
    carry <- carry + C
    cat(sprintf("%2d %+9.4f %+9.4f %+9.4f %+9.4f %+9.4f %+9.4f\n",
                t, q(E[1]), q(E[2]), q(E[3]), q(C), q(sum(E) + C), q(dWA)))

    WA <- WA + dWA
    WP <- WP * (1 + RP[t])
    WB <- WB * (1 + RB[t])
  }
  cat(strrep("-", 62), "\n")
  cat(sprintf("%2s %+9.4f %+9.4f %+9.4f %+9.4f %+9.4f %+9.4f\n",
              "", q(ewa[1]), q(ewa[2]), q(ewa[3]), q(carry), q(sum(ewa) + carry), q(WA)))

  F <- c(0, 0, 0)
  WP <- START_WEALTH
  for (t in seq_len(n)) {
    F <- (1 + RB[t]) * F + WP * G[t, ]
    WP <- WP * (1 + RP[t])
  }

  cPP <- cBB <- cPB <- cBP <- 0
  for (t in seq_len(n)) {
    cPP <- (1 + RP[t]) * (1 + cPP) - 1
    cBB <- (1 + RB[t]) * (1 + cBB) - 1
    cPB <- (1 + stats[[t]]$RPB) * (1 + cPB) - 1
    cBP <- (1 + stats[[t]]$RBP) * (1 + cBP) - 1
  }
  berg <- START_WEALTH * c(cPB - cBB, cBP - cBB, cPP - cPB - cBP + cBB)

  RPc <- prod(1 + RP) - 1
  RBc <- prod(1 + RB) - 1
  K <- (log1p(RPc) - log1p(RBc)) / (RPc - RBc)
  kt <- ifelse(abs(RA) < 1e-15, 1 / (1 + RP), (log1p(RP) - log1p(RB)) / RA)
  carino <- START_WEALTH * sapply(1:3, function(k) sum((kt / K) * G[, k]))

  mult <- sapply(seq_len(n), function(t) {
    left <- if (t == 1) 1 else prod(1 + RP[1:(t - 1)])
    right <- if (t == n) 1 else prod(1 + RB[(t + 1):n])
    left * right
  })
  rgrap <- START_WEALTH * sapply(1:3, function(k) sum(G[, k] * mult))

  rows <- list(
    list("EWA (explicit Carry)", ewa, carry),
    list("Frongello (2002)", F, NA),
    list("Berg (2014)", berg, NA),
    list("Carino (1999)", carino, NA),
    list("Reverse GRAP (1997)", rgrap, NA)
  )

  cat("\nFive methods on the same data (terminal, per 100 initial):\n")
  cat(sprintf("%-22s%9s%9s%9s%9s%9s\n", "Method", "A", "S", "I", "Carry", "Total"))
  for (row in rows) {
    name <- row[[1]]; effects <- row[[2]]; C <- row[[3]]
    total <- sum(effects) + ifelse(is.na(C), 0, C)
    stopifnot(abs(total - WA) < 1e-9)
    ctext <- ifelse(is.na(C), "-", sprintf("%+.4f", q(C)))
    cat(sprintf("%-22s%+9.4f%+9.4f%+9.4f%9s%+9.4f\n",
                name, q(effects[1]), q(effects[2]), q(effects[3]), ctext, q(total)))
  }

  D_path <- numeric(n)
  WP <- WB <- START_WEALTH
  for (t in seq_len(n)) {
    WP <- WP * (1 + RP[t]); WB <- WB * (1 + RB[t])
    D_path[t] <- WP - WB
  }

  ewa_c <- fro_c <- berg_c <- car_c <- rg_c <- numeric(n)
  WP <- START_WEALTH
  WAc <- Fc <- cPP <- cBB <- car_run <- rg_run <- 0
  for (t in seq_len(n)) {
    WAc <- WAc + WP * RA[t] + WAc * RB[t]
    Fc <- (1 + RB[t]) * Fc + WP * RA[t]
    cPP <- (1 + RP[t]) * (1 + cPP) - 1
    cBB <- (1 + RB[t]) * (1 + cBB) - 1
    car_run <- car_run + kt[t] / K * RA[t]
    rg_run <- rg_run + RA[t] * mult[t]

    ewa_c[t] <- WAc
    fro_c[t] <- Fc
    berg_c[t] <- START_WEALTH * (cPP - cBB)
    car_c[t] <- START_WEALTH * car_run
    rg_c[t] <- START_WEALTH * rg_run
    WP <- WP * (1 + RP[t])
  }

  cat("\nCumulative attribution through time vs realised active wealth\n")
  cat(sprintf("(horizon fixed at T = %d; per 100 initial):\n", n))
  cat(sprintf("%2s %9s %9s %10s %9s %9s %9s\n", "t", "W^A(t)", "EWA", "Frongello", "Berg", "Carino", "RevGRAP"))
  for (t in seq_len(n)) {
    cat(sprintf("%2d %+9.4f %+9.4f %+10.4f %+9.4f %+9.4f %+9.4f\n",
                t, q(D_path[t]), q(ewa_c[t]), q(fro_c[t]), q(berg_c[t]), q(car_c[t]), q(rg_c[t])))
  }
}

if (CASE == "base") {
  W_PORT <- W_PORT_BASE; W_B <- W_BENCH
} else if (CASE == "pure_selection") {
  W_PORT <- W_PORT_PURE_SEL; W_B <- W_BENCH_PURE_SEL
} else {
  stop("CASE must be 'base' or 'pure_selection'")
}

cat(sprintf("Explicit Wealth Attribution — case: %s\n\n", CASE))
run(W_PORT, W_B)
