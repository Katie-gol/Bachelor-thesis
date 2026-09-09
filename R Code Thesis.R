
##  Replication of Angrist & Lavy (1999), QJE 114(2), 533-575
##  "Using Maimonides' Rule to Estimate the Effect of Class Size on


## 0. Setup


library(haven)

OUT_DIR <- "output"
dir.create(OUT_DIR, showWarnings = FALSE)


grade4 <- read_dta("C:/Users/Katja/Documents/BA/final4.dta")
grade5 <- read_dta("C:/Users/Katja/Documents/BA/final5.dta")
grade4$grade <- 4L
grade5$grade <- 5L


## 1. Data preparation
##    Mirrors  AngristLavy_Table*.do files.


prepare <- function(dat) {
  
  d <- as.data.frame(zap_labels(dat))
  
  num <- c("schlcode", "c_size", "c_pik", "c_leom", "classize", "cohsize",
           "mathsize", "avgmath", "verbsize", "avgverb", "tipuach")
  d[num] <- lapply(d[num], as.numeric)
  
  ## Scores above 100 are coding errors; the authors subtract 100.
  i <- which(d$avgverb > 100); d$avgverb[i] <- d$avgverb[i] - 100
  i <- which(d$avgmath > 100); d$avgmath[i] <- d$avgmath[i] - 100
  
  ## Classes with no test takers get a missing score.
  d$avgverb[which(d$verbsize == 0)] <- NA
  d$avgmath[which(d$mathsize == 0)] <- NA
  
  ## Maimonides' rule, 
  d$func1 <- d$c_size / (floor((d$c_size - 1) / 40) + 1)
  
  ## Sample restrictions 
  ##   classes of 2-44 pupils, grade enrolment above 5,
  ##   c_leom == 1  Jewish public schools,
  ##   c_pik  <  3  secular or religious, excluding independent (haredi) schools
  keep <- with(d, classize > 1 & classize < 45 & c_size > 5 &
                 c_leom == 1 & c_pik < 3)
  keep[is.na(keep)] <- FALSE
  d <- d[keep, , drop = FALSE]
  
  d$c_size2 <- d$c_size^2 / 100
  

  e <- d$c_size
  d$trend <- ifelse(e <= 40,  e,
                    ifelse(e <= 80,  20      + e / 2,
                           ifelse(e <= 120, 100 / 3 + e / 3,
                                  ifelse(e <= 160, 130 / 3 + e / 4, NA))))
  
  ## +/-5 discontinuity sample
  d$disc <- as.numeric((e >= 36 & e <= 45) | (e >= 76 & e <= 85) |
                         (e >= 116 & e <= 125))
  
  ## Interaction terms for Table VII
  d$cs_pd <- d$classize * d$tipuach
  d$f_pd  <- d$func1    * d$tipuach
  
  rownames(d) <- NULL
  d
}

d5 <- prepare(grade5)   # 5th graders, tested 1991
d4 <- prepare(grade4)   # 4th graders, tested 1991


## 2. Moulton standard errors


moulton_factors <- function(resid, Xrho, cluster, kr) {
  
  n     <- length(resid)
  m     <- as.numeric(table(cluster))
  mbar  <- mean(m)
  varm  <- var(m)                    
  denom <- sum(m * (m - 1))
  gsize <- varm / mbar + mbar - 1     
  
  icc <- function(u, s2) {
    sj <- ave(u, cluster, FUN = sum)  
    sum(u * (sj - u)) / (s2 * denom)  
  }
  
  rho_e <- icc(resid, sum(resid^2) / (n - kr))
  
  fac <- numeric(kr)
  for (k in seq_len(ncol(Xrho))) {
    W  <- if (ncol(Xrho) == 1) matrix(1, n, 1)
    else cbind(1, Xrho[, -k, drop = FALSE])
    r  <- qr.resid(qr(W), Xrho[, k])         
    s2 <- sum(r^2) / (n - kr - 1)             
    fac[k] <- sqrt(1 + gsize * rho_e * icc(r, s2))
  }
  fac[kr] <- sqrt(1 + gsize * rho_e)        
  fac
}



## 3. Estimators.  The intercept is placed last, matching Stata


ols_moulton <- function(y, X, cluster) {
  
  n  <- length(y)
  Xm <- cbind(X, "(Intercept)" = 1)
  kr <- ncol(Xm)
  
  XtXinv <- chol2inv(chol(crossprod(Xm)))
  b      <- as.vector(XtXinv %*% crossprod(Xm, y))
  e      <- as.vector(y - Xm %*% b)
  s2     <- sum(e^2) / (n - kr)
  se_iid <- sqrt(s2 * diag(XtXinv))
  fac    <- moulton_factors(e, X, cluster, kr)
  
  list(term = colnames(Xm), coef = b, se = se_iid * fac, se_iid = se_iid,
       moulton = fac, rmse = sqrt(s2),
       r2 = 1 - sum(e^2) / sum((y - mean(y))^2), n = n, F1 = NA_real_)
}


tsls_moulton <- function(y, D, W, Zx, cluster) {
  
  n  <- length(y)
  D  <- as.matrix(D)
  Zx <- as.matrix(Zx)
  W  <- if (is.null(W)) NULL else as.matrix(W)
  
  X  <- cbind(D, W, "(Intercept)" = 1)
  Z  <- cbind(Zx, W, 1)
  kr <- ncol(X)
  
  Dhat <- qr.fitted(qr(Z), D)                
  Xhat <- cbind(Dhat, W, 1)
  
  b  <- as.vector(solve(crossprod(Xhat, X), crossprod(Xhat, y)))
  e  <- as.vector(y - X %*% b)
  s2 <- sum(e^2) / (n - kr)
  V  <- s2 * chol2inv(chol(crossprod(Xhat)))
  se_iid <- sqrt(diag(V))
  fac    <- moulton_factors(e, cbind(Dhat, W), cluster, kr)
  

  Wc   <- if (is.null(W)) matrix(1, n, 1) else cbind(W, 1)
  ssr1 <- sum(lm.fit(Z,  D[, 1])$residuals^2)
  ssr0 <- sum(lm.fit(Wc, D[, 1])$residuals^2)
  F1   <- ((ssr0 - ssr1) / ncol(Zx)) / (ssr1 / (n - ncol(Z)))
  
  list(term = colnames(X), coef = b, se = se_iid * fac, se_iid = se_iid,
       moulton = fac, rmse = sqrt(s2), r2 = NA_real_, n = n, F1 = F1)
}


## 4. Wrappers taking variable names, plus a tidy-output helper

est_ols <- function(dat, y, x, cluster = "schlcode") {
  s <- dat[complete.cases(dat[, c(y, x, cluster)]), ]
  ols_moulton(s[[y]], as.matrix(s[, x, drop = FALSE]), s[[cluster]])
}

est_iv <- function(dat, y, endog, inst, exog = character(0),
                   cluster = "schlcode") {
  v <- unique(c(y, endog, inst, exog, cluster))
  s <- dat[complete.cases(dat[, v]), ]
  W <- if (length(exog)) as.matrix(s[, exog, drop = FALSE]) else NULL
  tsls_moulton(s[[y]], as.matrix(s[, endog, drop = FALSE]), W,
               as.matrix(s[, inst, drop = FALSE]), s[[cluster]])
}

## one output row per coefficient
tidy_fit <- function(fit, info) {
  info <- info[rep(1L, length(fit$term)), , drop = FALSE]
  data.frame(info,
             term  = fit$term,
             coef  = round(fit$coef, 3),
             se    = round(fit$se, 3),
             tstat = round(fit$coef / fit$se, 2),
             rmse  = round(fit$rmse, 2),
             r2    = round(fit$r2, 3),
             F1    = round(fit$F1, 1),
             n     = fit$n,
             row.names = NULL)
}

rows <- function(lst) do.call(rbind, lst)



## 5. Table I 


vars_t1 <- c("classize", "c_size", "tipuach", "verbsize", "mathsize",
             "avgverb", "avgmath")

desc <- function(dat, label) {
  out <- list()
  for (v in vars_t1) {
    x <- dat[[v]][!is.na(dat[[v]])]
    q <- quantile(x, c(.10, .25, .50, .75, .90))
    out[[v]] <- data.frame(
      sample    = label,
      variable  = v,
      n_classes = nrow(dat),
      n_schools = length(unique(dat$schlcode)),
      mean = round(mean(x), 1), sd = round(sd(x), 1),
      p10 = round(q[1], 1), p25 = round(q[2], 1), p50 = round(q[3], 1),
      p75 = round(q[4], 1), p90 = round(q[5], 1), row.names = NULL)
  }
  rows(out)
}

tab1 <- list()
for (g in c("5th", "4th")) {
  dat <- if (g == "5th") d5 else d4
  dat <- dat[!is.na(dat$avgverb), ]
  tab1[[paste0(g, "_full")]] <- desc(dat, paste(g, "full sample"))
  tab1[[paste0(g, "_disc")]] <- desc(dat[dat$disc == 1, ],
                                     paste(g, "+/-5 discontinuity"))
}
tab1 <- rows(tab1)
print(tab1)
write.csv(tab1, file.path(OUT_DIR, "table1_descriptives.csv"),
          row.names = FALSE)



## 6. Table II 

specs_t2 <- list("no controls"      = c("classize"),
                 "+ PD"             = c("classize", "tipuach"),
                 "+ PD + enrolment" = c("classize", "tipuach", "c_size"))

tab2 <- list()
for (g in c("5th", "4th")) {
  dat <- if (g == "5th") d5 else d4
  for (y in c("avgverb", "avgmath")) {
    for (s in names(specs_t2)) {
      fit <- est_ols(dat, y, specs_t2[[s]])
      tab2[[paste(g, y, s)]] <- tidy_fit(
        fit, data.frame(grade = g, outcome = y, spec = s))
    }
  }
}
tab2 <- rows(tab2)
print(tab2)
write.csv(tab2, file.path(OUT_DIR, "table2_ols.csv"), row.names = FALSE)
mean_row <- function(dat, label) {
  do.call(rbind, lapply(c("avgverb", "avgmath"), function(v) {
    x <- dat[[v]][!is.na(dat[[v]])]
    data.frame(grade = label, outcome = v, n = length(x),
               mean = round(mean(x), 1), sd = round(sd(x), 1))
  }))
}
rbind(mean_row(d5, "5th"), mean_row(d4, "4th"))


## 7. Table III 


specs_t3 <- list("PD"             = c("func1", "tipuach"),
                 "PD + enrolment" = c("func1", "tipuach", "c_size"))

tab3 <- list()
for (g in c("5th", "4th")) {
  dat <- if (g == "5th") d5 else d4
  dat <- dat[!is.na(dat$avgverb), ]
  for (smp in c("full", "+/-5 disc")) {
    sub <- if (smp == "full") dat else dat[dat$disc == 1, ]
    for (y in c("classize", "avgverb", "avgmath")) {
      for (s in names(specs_t3)) {
        fit <- est_ols(sub, y, specs_t3[[s]])
        tab3[[paste(g, smp, y, s)]] <- tidy_fit(
          fit, data.frame(grade = g, sample = smp, outcome = y, spec = s))
      }
    }
  }
}
tab3 <- rows(tab3)
print(tab3)
write.csv(tab3, file.path(OUT_DIR, "table3_reduced_form.csv"),
          row.names = FALSE)



## 8. Tables IV and V 

specs_iv <- list(
  list(sample = "full",      controls = c("tipuach")),
  list(sample = "full",      controls = c("tipuach", "c_size")),
  list(sample = "full",      controls = c("tipuach", "c_size", "c_size2")),
  list(sample = "full",      controls = c("trend")),
  list(sample = "+/-5 disc", controls = c("tipuach")),
  list(sample = "+/-5 disc", controls = c("tipuach", "c_size"))
)

tab45 <- list()
for (g in c("5th", "4th")) {
  dat <- if (g == "5th") d5 else d4
  dat <- dat[!is.na(dat$avgverb), ]
  for (y in c("avgverb", "avgmath")) {
    for (i in seq_along(specs_iv)) {
      sp  <- specs_iv[[i]]
      sub <- if (sp$sample == "full") dat else dat[dat$disc == 1, ]
      fit <- est_iv(sub, y, "classize", "func1", sp$controls)
      tab45[[paste(g, y, i)]] <- tidy_fit(
        fit, data.frame(grade = g, outcome = y, sample = sp$sample,
                        controls = paste(sp$controls, collapse = " + ")))
    }
  }
}
tab45 <- rows(tab45)
print(tab45)
write.csv(tab45, file.path(OUT_DIR, "table4_5_2sls.csv"), row.names = FALSE)



## 9. Table VI 

disc_sample <- function(dat, w) {
  lo <- if (w == 5) c(36, 76, 116) else c(38, 78, 118)
  hi <- lo + 2 * w - 1
  e  <- dat$c_size
  s  <- dat[(e >= lo[1] & e <= hi[1]) | (e >= lo[2] & e <= hi[2]) |
              (e >= lo[3] & e <= hi[3]), ]
  e  <- s$c_size
  s$seg1 <- as.numeric(e >= lo[1] & e <= hi[1])
  s$seg2 <- as.numeric(e >= lo[2] & e <= hi[2])
  s$z1   <- as.numeric(e >= lo[1] + w & e <= hi[1])
  s$z2   <- as.numeric(e >= lo[2] + w & e <= hi[2])
  s$z3   <- as.numeric(e >= lo[3] + w & e <= hi[3])
  s
}

tab6 <- list()
for (g in c("5th", "4th")) {
  dat <- if (g == "5th") d5 else d4
  dat <- dat[!is.na(dat$avgverb), ]
  for (w in c(5, 3)) {
    sub   <- disc_sample(dat, w)
    ctrls <- if (w == 5) list(c("tipuach", "seg1", "seg2"))
    else        list(c("tipuach", "seg1", "seg2"), c("seg1", "seg2"))
    for (y in c("avgverb", "avgmath")) {
      for (j in seq_along(ctrls)) {
        fit <- est_iv(sub, y, "classize", c("z1", "z2", "z3"), ctrls[[j]])
        tab6[[paste(g, w, y, j)]] <- tidy_fit(
          fit, data.frame(grade = g, sample = paste0("+/-", w), outcome = y,
                          controls = paste(ctrls[[j]], collapse = " + ")))
      }
    }
  }
}
tab6 <- rows(tab6)
print(tab6)
write.csv(tab6, file.path(OUT_DIR, "table6_dummy_instruments.csv"),
          row.names = FALSE)



## 10. Table VII 


pool_vars <- c("schlcode", "grade", "classize", "c_size", "c_size2",
               "tipuach", "func1", "cs_pd", "f_pd", "avgverb", "avgmath",
               "disc", "trend")
pooled     <- rbind(d5[, pool_vars], d4[, pool_vars])
pooled$g4  <- as.numeric(pooled$grade == 4)   # the Grade 4 dummy of Table VII

tab7 <- list()
for (g in c("5th", "4th")) {
  dat <- if (g == "5th") d5 else d4
  dat <- dat[!is.na(dat$avgverb), ]
  for (y in c("avgverb", "avgmath")) {
    fit <- est_iv(dat, y, c("classize", "cs_pd"), c("func1", "f_pd"),
                  c("tipuach", "c_size"))
    tab7[[paste(g, y)]] <- tidy_fit(
      fit, data.frame(grade = g, outcome = y, spec = "interaction"))
  }
}

pdat <- pooled[!is.na(pooled$avgverb), ]
for (y in c("avgverb", "avgmath")) {
  fit <- est_iv(pdat, y, "classize", "func1",
                c("tipuach", "c_size", "g4"))
  tab7[[paste("pooled", y, "main")]] <- tidy_fit(
    fit, data.frame(grade = "pooled", outcome = y, spec = "no interaction"))
  
  fit <- est_iv(pdat, y, c("classize", "cs_pd"), c("func1", "f_pd"),
                c("tipuach", "c_size", "g4"))
  tab7[[paste("pooled", y, "int")]] <- tidy_fit(
    fit, data.frame(grade = "pooled", outcome = y, spec = "interaction"))
}
tab7 <- rows(tab7)
print(tab7)
write.csv(tab7, file.path(OUT_DIR, "table7_interactions.csv"),
          row.names = FALSE)


## 11. Figures


maimonides <- function(e) e / (floor((e - 1) / 40) + 1)

## Figure I
fig1 <- function(dat, main) {
  s  <- dat[dat$c_size <= 220, ]
  ag <- aggregate(list(classize = s$classize), by = list(e = s$c_size),
                  FUN = mean)
  plot(ag$e, ag$classize, type = "l", ylim = c(5, 42), xlim = c(0, 220),
       xlab = "Enrolment count", ylab = "Class size", main = main)
  ge <- 5:220
  lines(ge, maimonides(ge), lty = 2)
  abline(h = c(20.5, 26.67, 30, 32, 33.33, 40), col = "grey80", lty = 3)
  legend("bottomright", c("Actual class size", "Maimonides' rule"),
         lty = c(1, 2), bty = "n")
}

pdf(file.path(OUT_DIR, "figure1_class_size.pdf"), width = 8, height = 9)
par(mfrow = c(2, 1))
fig1(d5, "a. Fifth grade")
fig1(d4, "b. Fourth grade")
dev.off()

## Binning for Figures II and III: 

bin10 <- function(e) ifelse(e >= 160, 165, 10 * floor(e / 10) + 5)

binned <- function(dat, min_schools = 10) {
  s <- dat[!is.na(dat$avgverb) & dat$c_size >= 9 & dat$c_size <= 190, ]
  s$bin <- bin10(s$c_size)
  ag <- data.frame(bin = sort(unique(s$bin)))
  for (v in c("avgverb", "avgmath", "func1", "c_size", "tipuach")) {
    m <- tapply(s[[v]], s$bin, function(z) mean(z, na.rm = TRUE))
    ag[[v]] <- as.numeric(m[as.character(ag$bin)])
  }
  ns <- tapply(s$schlcode, s$bin, function(z) length(unique(z)))
  ag$n_schools <- as.numeric(ns[as.character(ag$bin)])
  ag <- ag[ag$n_schools >= min_schools, , drop = FALSE]
  
## Figure III 
  for (v in c("avgverb", "avgmath", "func1")) {
    ag[[paste0(v, "_res")]] <- residuals(lm(ag[[v]] ~ ag$c_size + ag$tipuach))
  }
  ag
}

a5 <- binned(d5)
a4 <- binned(d4)
print(a5[, c("bin", "n_schools", "avgverb", "func1")])

two_scale <- function(x, y1, y2, main, lab1, lab2,
                      leg1, leg2, ylim1 = NULL, ylim2 = NULL,
                      leg_pos = "bottomright") {
  par(mar = c(5, 4, 4, 5))
  plot(x, y1, type = "l", lty = 1, lwd = 1.6, xlim = c(0, 175), ylim = ylim1,
       xlab = "Enrolment count", ylab = lab1, main = main)
  par(new = TRUE)
  plot(x, y2, type = "l", lty = 2, lwd = 1.6, xlim = c(0, 175), ylim = ylim2,
       axes = FALSE, xlab = "", ylab = "")
  axis(4)
  mtext(lab2, side = 4, line = 3)
  legend(leg_pos, legend = c(leg1, leg2), lty = c(1, 2), lwd = 1.6,
         bty = "n", cex = 0.85)
}

pdf(file.path(OUT_DIR, "figure2_scores_and_rule.pdf"), width = 8, height = 9)
par(mfrow = c(2, 1))
two_scale(a5$bin, a5$avgverb, a5$func1, "a. Fifth grade",
          "Average reading score", "Average class-size function",
          "Average reading score", "Predicted class size",
          ylim2 = c(5, 40))
two_scale(a4$bin, a4$avgverb, a4$func1, "b. Fourth grade",
          "Average reading score", "Average class-size function",
          "Average reading score", "Predicted class size",
          ylim2 = c(5, 40))
dev.off()

pdf(file.path(OUT_DIR, "figure3_residuals.pdf"), width = 8, height = 11)
par(mfrow = c(3, 1))
two_scale(a5$bin, a5$avgverb_res, a5$func1_res, "a. Fifth grade (reading)",
          "Reading score residual", "Size-function residual",
          "Average test scores", "Predicted class size",
          ylim1 = c(-5, 5), ylim2 = c(-15, 15))
two_scale(a4$bin, a4$avgverb_res, a4$func1_res, "b. Fourth grade (reading)",
          "Reading score residual", "Size-function residual",
          "Average test scores", "Predicted class size",
          ylim1 = c(-5, 5), ylim2 = c(-15, 15))
two_scale(a5$bin, a5$avgmath_res, a5$func1_res, "c. Fifth grade (math)",
          "Math score residual", "Size-function residual",
          "Average test scores", "Predicted class size",
          ylim1 = c(-5, 5), ylim2 = c(-15, 15))
dev.off()


## Figure IV
fig4 <- function(dat, main) {
  s  <- dat[!is.na(dat$avgverb) & dat$disc == 1, ]
  lo <- s$classize[s$func1 <  32]
  hi <- s$classize[s$func1 >= 32]
  grid <- 15:41
  plot(grid, ecdf(lo)(grid), type = "s", lty = 2, lwd = 1.6,
       xlim = c(15, 41), ylim = c(0, 1),
       xlab = "Class size", ylab = "Cumulative share of classes", main = main)
  lines(grid, ecdf(hi)(grid), type = "s", lty = 1, lwd = 1.6)
  legend("topleft", bty = "n", lty = c(2, 1), lwd = 1.6, cex = 0.85,
         legend = c(sprintf("Maimonides' rule < 32  (n = %d)", length(lo)),
                    sprintf("Maimonides' rule >= 32 (n = %d)", length(hi))))
}

pdf(file.path(OUT_DIR, "figure4_cdfs.pdf"), width = 8, height = 9)
par(mfrow = c(2, 1))
fig4(d5, "a. Fifth grade")
fig4(d4, "b. Fourth grade")
dev.off()

################################################################################

# =====================================================================
# Replication Pop-Eleches & Urquiola (2013) 


library(haven)
library(fixest)
library(dplyr)
library(rdrobust)



data1 <- read_dta("C:/Users/Katja/Documents/BA/112645-V1/data/data-AER-1.dta")  
data2 <- read_dta("C:/Users/Katja/Documents/BA/112645-V1/data/data-AER-2.dta")   
data3 <- read_dta("C:/Users/Katja/Documents/BA/112645-V1/data/data-AER-3.dta")
data7 <- read_dta("C:/Users/Katja/Documents/BA/112645-V1/data/data-AER-7.dta")

#  Helpers 

get_bw <- function(data, outcome, running = "dzag", subset_expr = NULL) {
  d <- data %>% filter(.data[[running]] != 0, !is.na(.data[[outcome]]))
  if (!is.null(subset_expr)) d <- d %>% filter(!!rlang::parse_expr(subset_expr))
  rdbwselect(y = d[[outcome]], x = d[[running]], c = 0, bwselect = "mserd")$bws[1]
}

run_rd_spec <- function(data, outcome, bw = NULL, subset_expr = NULL) {
  d <- data %>% filter(dzag != 0)
  if (!is.null(subset_expr)) d <- d %>% filter(!!rlang::parse_expr(subset_expr))
  if (!is.null(bw) && !is.na(bw)) d <- d %>% filter(dzag >= -bw, dzag <= bw)
  fml <- as.formula(paste0(outcome, " ~ dga + dzag + dzag_after | uazY"))
  feols(fml, data = d, cluster = ~sid2)
}

run_rd_iv <- function(data, bw = NULL, subset_expr = NULL, outcome = "bcg") {
  d <- data %>% filter(dzag != 0)
  if (!is.null(subset_expr)) d <- d %>% filter(!!rlang::parse_expr(subset_expr))
  if (!is.null(bw) && !is.na(bw)) d <- d %>% filter(dzag >= -bw, dzag <= bw)
  fml <- as.formula(paste0(outcome, " ~ dzag + dzag_after | uazY | agus ~ dga"))
  feols(fml, data = d, cluster = ~sid2)
}

# Full 4-column block for panel Table 3 / Table 4
panel_4col <- function(data, outcome, label) {
  bw_all <- get_bw(data, outcome)
  bw_srv <- get_bw(data, outcome, subset_expr = "survey == 1")
  specs <- list(
    c1_all_1pt     = run_rd_spec(data, outcome, bw = 1),
    c2_all_ik      = run_rd_spec(data, outcome, bw = bw_all),
    c3_survey_1pt  = run_rd_spec(data, outcome, bw = 1,      subset_expr = "survey == 1"),
    c4_survey_ik   = run_rd_spec(data, outcome, bw = bw_srv, subset_expr = "survey == 1")
  )
  cat("\n########", label, " (MSE-opt bw: all =", round(bw_all, 3),
      ", survey =", round(bw_srv, 3), ") ########\n")
  print(etable(specs, keep = "dga", tex = FALSE))
  invisible(specs)
}


# TABLE 3 


t3_A <- panel_4col(data1, "agus",   "TABLE 3 Panel A: school-level score, between-school")
t3_B <- panel_4col(data3, "agus2B", "TABLE 3 Panel B: track-level score, between-TRACK")
t3_C <- panel_4col(data1, "agus2",  "TABLE 3 Panel C: track-level score, between-school")
t3_D <- panel_4col(data2, "agus",   "TABLE 3 Panel D: 2005-2007 cohorts, between-school")

# Panel D columns 5-6
bw_d7 <- get_bw(data7, "agus")
t3_D_survey <- list(
  c5_surveydata_1pt = run_rd_spec(data7, "agus", bw = 1),
  c6_surveydata_ik  = run_rd_spec(data7, "agus", bw = bw_d7)
)
cat("\n######## TABLE 3 Panel D, cols 5-6: survey data (data-AER-7) ########\n")
print(etable(t3_D_survey, keep = "dga", tex = FALSE))


# TABLE 4 


t4_A <- panel_4col(data1, "bct", "TABLE 4 Panel A: Bacc. taken, between-school")
t4_B <- panel_4col(data1, "bcg", "TABLE 4 Panel B: Bacc. grade, between-school")

# Panel C: IV

bw_bcg     <- get_bw(data1, "bcg")
bw_bcg_srv <- get_bw(data1, "bcg", subset_expr = "survey == 1")
t4_C <- list(
  c1_all_1pt    = run_rd_iv(data1, bw = 1),
  c2_all_ik     = run_rd_iv(data1, bw = bw_bcg),
  c3_survey_1pt = run_rd_iv(data1, bw = 1,          subset_expr = "survey == 1"),
  c4_survey_ik  = run_rd_iv(data1, bw = bw_bcg_srv, subset_expr = "survey == 1")
)
cat("\n######## TABLE 4 Panel C: Bacc. grade, IV specification ########\n")
print(etable(t4_C, keep = "agus", tex = FALSE))

t4_D <- panel_4col(data3, "bct", "TABLE 4 Panel D: Bacc. taken, between-TRACK")
t4_E <- panel_4col(data3, "bcg", "TABLE 4 Panel E: Bacc. grade, between-TRACK")



# TABLE 5 


t5_panel <- function(subset_expr, label) {
  specs <- list(
    c1_first_stage = run_rd_spec(data1, "agus", bw = 1, subset_expr = subset_expr),
    c2_bacc_taken  = run_rd_spec(data1, "bct",  bw = 1, subset_expr = subset_expr),
    c3_bacc_grade  = run_rd_spec(data1, "bcg",  bw = 1, subset_expr = subset_expr),
    c4_bacc_iv     = run_rd_iv(  data1,         bw = 1, subset_expr = subset_expr)
  )
  cat("\n########", label, "########\n")
  print(etable(specs, keep = c("dga", "agus"), tex = FALSE))
  invisible(specs)
}

t5_A <- t5_panel(NULL,             "TABLE 5 Panel A: full sample")
t5_B <- t5_panel("zga > 7.74",     "TABLE 5 Panel B: top tercile of cutoffs")
t5_C <- t5_panel("zga < 6.77",     "TABLE 5 Panel C: bottom tercile of cutoffs")
t5_D <- t5_panel("nusua >= 4",     "TABLE 5 Panel D: towns with 4+ schools")
t5_E <- t5_panel("nusua == 3",     "TABLE 5 Panel E: towns with 3 schools")
t5_F <- t5_panel("nusua == 2",     "TABLE 5 Panel F: towns with 2 schools")




library(haven)
library(fixest)
library(dplyr)
library(rdrobust)




data8 <-read_dta("C:/Users/Katja/Documents/BA/112645-V1/data/data-AER-8.dta")  

# 0. Helpers


# MSE-optimal bandwidth 
get_bw <- function(data, outcome, running = "dzag") {
  d <- data %>% filter(.data[[running]] != 0, !is.na(.data[[outcome]]))
  bw <- tryCatch(
    rdbwselect(y = d[[outcome]], x = d[[running]], c = 0, bwselect = "mserd")$bws[1],
    error = function(e) NA_real_
  )
  bw
}

# outcome ~ treat + running + running:treat | FE, clustered
rd_spec <- function(data, outcome, bw = NULL, subset_expr = NULL,
                    treat = "dga", running = "dzag", after = "dzag_after",
                    fe = "uazY", cluster = "sid2") {
  d <- data %>% filter(.data[[running]] != 0)
  if (!is.null(subset_expr)) d <- d %>% filter(!!rlang::parse_expr(subset_expr))
  if (!is.null(bw) && !is.na(bw)) {
    d <- d %>% filter(.data[[running]] >= -bw, .data[[running]] <= bw)
  }
  fml <- as.formula(paste0(outcome, " ~ ", treat, " + ", running, " + ", after,
                           " | ", fe))
  feols(fml, data = d, cluster = as.formula(paste0("~", cluster)))
}



# 1. TABLE 1 


d456 <- bind_rows(
  read_dta("C:/Users/Katja/Documents/BA/112645-V1/data/data-AER-4.dta"),
  read_dta("C:/Users/Katja/Documents/BA/112645-V1/data/data-AER-5.dta"),
  read_dta("C:/Users/Katja/Documents/BA/112645-V1/data/data-AER-6.dta")
)

table1_panel <- function(d, label) {
  ind <- d %>% group_by(year) %>%
    summarise(n            = n(),
              grade_mean   = mean(grade, na.rm = TRUE),
              grade_sd     = sd(grade,   na.rm = TRUE),
              bct_mean     = mean(bct,   na.rm = TRUE),
              bct_sd       = sd(bct,     na.rm = TRUE),
              bcg_mean     = mean(bcg,   na.rm = TRUE),
              bcg_sd       = sd(bcg,     na.rm = TRUE),
              bcg_n        = sum(!is.na(bcg)), .groups = "drop")
  
  trk <- d %>% count(year, us2) %>% group_by(year) %>%
    summarise(track_students_mean = mean(n), track_students_sd = sd(n),
              n_tracks = n(), .groups = "drop")
  
  sch <- d %>% count(year, us) %>% group_by(year) %>%
    summarise(school_students_mean = mean(n), school_students_sd = sd(n),
              n_schools = n(), .groups = "drop")
  
  twn <- d %>% count(year, ua) %>% group_by(year) %>%
    summarise(town_students_mean = mean(n), town_students_sd = sd(n),
              n_towns = n(), .groups = "drop")
  
  cat("\n==== TABLE 1 --", label, "====\n")
  print(as.data.frame(ind)); print(as.data.frame(trk))
  print(as.data.frame(sch)); print(as.data.frame(twn))
}

table1_panel(d456,                          "Panel A: all towns")
table1_panel(d456 %>% filter(survey == 1),  "Panel B: survey towns")


# 2. TABLE 2 


table2_vars <- list(
  `A. Socioeconomic` = c("head_sex", "head_age",
                         "head_nat_romanian", "head_nat_hungarian",
                         "head_nat_gypsy", "head_nat_other",
                         "head_educ_primary", "head_educ_sec",
                         "head_educ_tertiary", "ch_sex", "ch_age"),
  `B. Parental responses` = c("p_d_parent_volunteer", "p_tutoring",
                              "p_d_homework_help", "p_d_homework"),
  `C. Child responses` = c("ch_rank_peers", "ch_peers_index_bad",
                           "ch_d_homework", "ch_rank_homework_index"),
  `D. Teacher qualifications` = c("didactic_Romanian", "experience_Romanian",
                                  "novice_Romanian")
)

for (panel in names(table2_vars)) {
  cat("\n==== TABLE 2 --", panel, "====\n")
  out <- lapply(table2_vars[[panel]], function(v) {
    if (!v %in% names(data7)) return(NULL)
    data.frame(variable = v,
               mean = mean(data7[[v]], na.rm = TRUE),
               sd   = sd(data7[[v]],   na.rm = TRUE),
               n    = sum(!is.na(data7[[v]])))
  })
  print(do.call(rbind, out))
}


# 3. TABLES 6, 7, 8 

#   Panel A (school level): cluster school-cohort
#   Panel B (track  level): cluster school-track-cohort
#   Panel C (student level): cluster student 

survey_panels <- list(
  A_school  = list(suffix_school = TRUE,  cluster = "usY"),
  B_track   = list(cluster = "us2BY"),
  C_student = list(cluster = "sid2")
)

survey_block <- function(var_school, var_track, var_student, label) {
  cat("\n########", label, "########\n")
  specs <- list()
  
  if (!is.na(var_school) && var_school %in% names(data7)) {
    bw <- get_bw(data7, var_school)
    specs$A_1pt <- rd_spec(data7, var_school, bw = 1,  cluster = "usY")
    specs$A_ik  <- rd_spec(data7, var_school, bw = bw, cluster = "usY")
  }
  if (!is.na(var_track) && var_track %in% names(data7)) {
    bw <- get_bw(data7, var_track)
    specs$B_1pt <- rd_spec(data7, var_track, bw = 1,  cluster = "us2BY")
    specs$B_ik  <- rd_spec(data7, var_track, bw = bw, cluster = "us2BY")
  }
  if (!is.na(var_student) && var_student %in% names(data7)) {
    bw <- get_bw(data7, var_student)
    specs$C_1pt <- rd_spec(data7, var_student, bw = 1,  cluster = "sid2")
    specs$C_ik  <- rd_spec(data7, var_student, bw = bw, cluster = "sid2")
  }
  print(etable(specs, keep = "dga", tex = FALSE))
  invisible(specs)
}

## TABLE 6
t6_principal <- survey_block("sc_bestintown_teacherquality", NA, NA,
                             "T6 c1-2: principal says best teachers in town (school only)")
t6_certified <- survey_block("didactic_Romanian11", "didactic_Romanian2", "didactic_Romanian",
                             "T6 c3-4: Language teacher highest certification")
t6_experience <- survey_block("experience_Romanian11", "experience_Romanian2", "experience_Romanian",
                              "T6 c5-6: Language teacher experience (years)")
t6_novice <- survey_block("novice_Romanian11", "novice_Romanian2", "novice_Romanian",
                          "T6 c7-8: Language teacher is a novice (<2 yrs)")

## TABLE 7
t7_principal <- survey_block("sc_bestintown_parental", NA, NA,
                             "T7 c1-2: principal says best parental participation")
t7_volunteer <- survey_block("p_d_parent_volunteer1", "p_d_parent_volunteer2", "p_d_parent_volunteer",
                             "T7 c3-4: parents volunteered in past year")
t7_tutoring <- survey_block("p_tutoring1", "p_tutoring2", "p_tutoring",
                            "T7 c5-6: parents paid for tutoring")
t7_homework <- survey_block("p_d_homework_help1", "p_d_homework_help2", "p_d_homework_help",
                            "T7 c7-8: parents help child with homework often")

## TABLE 8
t8_principal <- survey_block("sc_bestintown_studentquality", NA, NA,
                             "T8 c1-2: principal says best student quality")
t8_rank <- survey_block("ch_rank_peers1", "ch_rank_peers2", "ch_rank_peers",
                        "T8 c3-4: child's perception of rank in track")
t8_negative <- survey_block("ch_peers_index_bad1", "ch_peers_index_bad2", "ch_peers_index_bad",
                            "T8 c5-6: child's negative interactions with peers")


## TABLE 9 
t9_spec <- function(outcome, cluster) {
  rd_spec(data8, outcome, bw = 1,
          treat = "dg_pr", running = "dzg_pr", after = "dzg_pr_after",
          fe = "usYprspB", cluster = cluster)
}

# Panel A 
t9_panelA <- list(
  c1_peer_quality = t9_spec("agus_class",             "usYclass"),
  c2_certified    = t9_spec("didactic_Romanian",      "usYclass"),
  c3_experience   = t9_spec("experience_Romanian",    "usYclass"),
  c4_novice       = t9_spec("novice_Romanian",        "usYclass"),
  c5_volunteer    = t9_spec("p_d_parent_volunteer3",  "usYclass"),
  c6_tutoring     = t9_spec("p_tutoring3",            "usYclass"),
  c7_hw_help      = t9_spec("p_d_homework_help3",     "usYclass"),
  c8_rank         = t9_spec("ch_rank_peers3",         "usYclass"),
  c9_negative     = t9_spec("ch_peers_index_bad3",    "usYclass"),
  c10_hw_child    = t9_spec("ch_d_homework3",         "usYclass"),
  c11_hw_parent   = t9_spec("p_d_homework3",          "usYclass"),
  c12_hw_easy     = t9_spec("ch_rank_homework_index3","usYclass")
)

print(etable(t9_panelA, keep = "dg_pr", tex = FALSE))

# Panel B 
t9_panelB <- list(
  c5_volunteer  = t9_spec("p_d_parent_volunteer",   "sid2"),
  c6_tutoring   = t9_spec("p_tutoring",             "sid2"),
  c7_hw_help    = t9_spec("p_d_homework_help",      "sid2"),
  c8_rank       = t9_spec("ch_rank_peers",          "sid2"),
  c9_negative   = t9_spec("ch_peers_index_bad",     "sid2"),
  c10_hw_child  = t9_spec("ch_d_homework",          "sid2"),
  c11_hw_parent = t9_spec("p_d_homework",           "sid2"),
  c12_hw_easy   = t9_spec("ch_rank_homework_index", "sid2")
)

print(etable(t9_panelB, keep = "dg_pr", tex = FALSE))


## Tabel 10

t10_vars <- c(c1_first_stage = "agus",
              c2_experience  = "experience_Romanian",
              c3_novice      = "novice_Romanian",
              c4_tutoring    = "p_tutoring",
              c5_hw_help     = "p_d_homework_help",
              c6_rank        = "ch_rank_peers",
              c7_negative    = "ch_peers_index_bad")

# Panel A
t10_panelA <- lapply(t10_vars, function(v)
  rd_spec(data7, v, bw = 1, cluster = "sid2"))

print(etable(t10_panelA, keep = "dga", tex = FALSE))

# Panel B
t10_panelB <- lapply(t10_vars, function(v) {
  d <- data7 %>% filter(dzag != 0, dzag >= -1, dzag <= 1)   # <= , matches Panel A
  fml <- as.formula(paste0(v, " ~ dga + dga_Y5 + dga_Y6 + dzag + dzag_after | uazY"))
  feols(fml, data = d, cluster = ~sid2)
})

print(etable(t10_panelB, keep = c("dga", "dga_Y5", "dga_Y6"), tex = FALSE))




###############################################################################
################################################################################
################################################################################
################################################################################

## Angrist Lavy MODERN

## ---------------------------------------------------------------------------
## 0. Setup and helpers
## ---------------------------------------------------------------------------



library(haven)
library(rdrobust)
library(rddensity)
library(cowplot)
library(RDHonest)

OUT_DIR <- "output_modern"
dir.create(OUT_DIR, showWarnings = FALSE)

CUTOFFS      <- c(40, 80, 120)          
TRUE_CUTOFFS <- c(40, 80, 120, 160)     


## Number of SCHOOLS 

school_counts <- function(dat, c0, h) {
  if (is.na(h)) return(c(schools_left = NA_integer_, schools_right = NA_integer_))
  s <- dat[!is.na(dat$c_size) & abs(dat$c_size - c0) <= h, , drop = FALSE]
  c(schools_left  = length(unique(s$schlcode[s$c_size <= c0])),
    schools_right = length(unique(s$schlcode[s$c_size >  c0])))
}


crosses_cutoff <- function(c0, h, cuts = TRUE_CUTOFFS) {
  if (is.na(h)) return(NA)
  other <- setdiff(cuts, c0)
  any(other >= c0 - h & other <= c0 + h)
}








## 1. Data preparation



prepare <- function(dat) {
  d <- as.data.frame(zap_labels(dat))
  num <- c("schlcode", "c_size", "c_pik", "c_leom", "classize", "cohsize",
           "mathsize", "avgmath", "verbsize", "avgverb", "tipuach")
  d[num] <- lapply(d[num], as.numeric)
  
  i <- which(d$avgverb > 100); d$avgverb[i] <- d$avgverb[i] - 100
  i <- which(d$avgmath > 100); d$avgmath[i] <- d$avgmath[i] - 100
  d$avgverb[which(d$verbsize == 0)] <- NA
  d$avgmath[which(d$mathsize == 0)] <- NA
  
  
  d$func1 <- d$c_size / (floor((d$c_size - 1) / 40) + 1)
  
  keep <- with(d, classize > 1 & classize < 45 & c_size > 5 &
                 c_leom == 1 & c_pik < 3)
  keep[is.na(keep)] <- FALSE
  d <- d[keep, , drop = FALSE]
  
  
  d <- d[!is.na(d$avgverb) | !is.na(d$avgmath), , drop = FALSE]
  rownames(d) <- NULL
  d
}

grade4 <- read_dta("C:/Users/Katja/Documents/BA/final4.dta")
grade5 <- read_dta("C:/Users/Katja/Documents/BA/final5.dta")

d5 <- prepare(grade5)
d4 <- prepare(grade4)
GRADES <- list("5th" = d5, "4th" = d4)


## 1b. Descriptive statistics 


desc <- do.call(rbind, lapply(names(GRADES), function(g) {
  dd <- GRADES[[g]]
  data.frame(
    grade          = g,
    classes        = nrow(dd),
    schools        = length(unique(dd$schlcode)),
    classes_verb   = sum(!is.na(dd$avgverb)),
    classes_math   = sum(!is.na(dd$avgmath)),
    ## how many maths classes v1 discarded by conditioning on reading
    math_lost_v1   = sum(!is.na(dd$avgmath) & is.na(dd$avgverb)),
    mean_avgverb   = round(mean(dd$avgverb, na.rm = TRUE), 2),
    sd_avgverb     = round(sd(dd$avgverb,   na.rm = TRUE), 2),
    mean_avgmath   = round(mean(dd$avgmath, na.rm = TRUE), 2),
    sd_avgmath     = round(sd(dd$avgmath,   na.rm = TRUE), 2),
    mean_classize  = round(mean(dd$classize, na.rm = TRUE), 2),
    sd_classize    = round(sd(dd$classize,   na.rm = TRUE), 2),
    mean_tipuach   = round(mean(dd$tipuach, na.rm = TRUE), 2),
    sd_tipuach     = round(sd(dd$tipuach,   na.rm = TRUE), 2))
}))
print(desc)
write.csv(desc, file.path(OUT_DIR, "descriptives.csv"), row.names = FALSE)



## 2. How much data sits near each cutoff



window_counts <- function(dat, label) {
  sch <- dat[!duplicated(dat$schlcode), ]
  out <- list()
  for (c0 in TRUE_CUTOFFS) {
    for (h in c(3, 5, 10, 15, 20)) {
      out[[paste(c0, h)]] <- data.frame(
        grade = label, cutoff = c0, h = h,
        schools       = sum(abs(sch$c_size - c0) <= h),
        schools_left  = sum(sch$c_size >= c0 - h & sch$c_size <= c0),
        schools_right = sum(sch$c_size >  c0     & sch$c_size <= c0 + h),
        classes       = sum(abs(dat$c_size - c0) <= h),
        classes_left  = sum(dat$c_size >= c0 - h & dat$c_size <= c0),
        classes_right = sum(dat$c_size >  c0     & dat$c_size <= c0 + h))
    }
  }
  do.call(rbind, out)
}

counts <- do.call(rbind, lapply(names(GRADES),
                                function(g) window_counts(GRADES[[g]], g)))
print(counts)
write.csv(counts, file.path(OUT_DIR, "window_counts.csv"), row.names = FALSE)




## 3. Manipulation 

dens <- list()
for (g in names(GRADES)) {
  sch <- GRADES[[g]][!duplicated(GRADES[[g]]$schlcode), ]
  for (c0 in CUTOFFS) {
    rdd <- tryCatch(rddensity(X = sch$c_size, c = c0 + 0.5), error = function(e) NULL)
    pv  <- if (is.null(rdd)) NA else rdd$test$p_jk
    
    for (w in c(3, 5)) {
      L <- sum(sch$c_size >= c0 - w + 1 & sch$c_size <= c0)
      R <- sum(sch$c_size >= c0 + 1     & sch$c_size <= c0 + w)
      bt <- binom.test(R, L + R, p = 0.5, alternative = "two.sided")
      dens[[paste(g, c0, w)]] <- data.frame(
        grade = g, cutoff = c0, window = w, left = L, right = R,
        binom_p = round(bt$p.value, 4), rddensity_p = round(pv, 4))
    }
  }
}
dens <- do.call(rbind, dens)
print(dens)
write.csv(dens, file.path(OUT_DIR, "manipulation_tests.csv"), row.names = FALSE)


## 4. First stage by cutoff


fs <- list()
for (g in names(GRADES)) {
  dd <- GRADES[[g]]
  for (c0 in CUTOFFS) {
    r <- tryCatch(
      rdrobust(y = dd$classize, x = dd$c_size, c = c0 + 0.5,
               cluster = dd$schlcode, kernel = "triangular", p = 1),
      error = function(e) NULL)
    if (is.null(r)) next
    h  <- r$bws[1, 1]
    sc <- school_counts(dd, c0, h)
    fs[[paste(g, c0)]] <- data.frame(
      grade = g, cutoff = c0,
      coef = round(r$coef[1], 3), coef_bc = round(r$coef[3], 3),
      se_conv = round(r$se[1], 3),
      se_robust = round(r$se[3], 3),
      p_conv   = round(r$pv[1], 4),
      p_robust = round(r$pv[3], 4),
      h = round(h, 2),
      classes_left = r$N_h[1], classes_right = r$N_h[2],
      schools_left = sc[["schools_left"]], schools_right = sc[["schools_right"]],
      window_hits_other_cutoff = crosses_cutoff(c0, h))
  }
}
fs <- do.call(rbind, fs)
print(fs)
write.csv(fs, file.path(OUT_DIR, "first_stage_by_cutoff.csv"), row.names = FALSE)



## 5. Fuzzy RD estimates



frd <- list()
for (g in names(GRADES)) {
  dd0 <- GRADES[[g]]
  for (y in c("avgverb", "avgmath")) {
    dd <- dd0[!is.na(dd0[[y]]), , drop = FALSE]
    for (c0 in CUTOFFS) {
      r <- tryCatch(
        rdrobust(y = dd[[y]], x = dd$c_size, c = c0 + 0.5, fuzzy = dd$classize,
                 cluster = dd$schlcode, kernel = "triangular", p = 1),
        error = function(e) NULL)
      if (is.null(r)) next
      h  <- r$bws[1, 1]
      sc <- school_counts(dd, c0, h)
      frd[[paste(g, y, c0)]] <- data.frame(
        grade = g, outcome = y, cutoff = c0,
        coef = round(r$coef[1], 3), coef_bc = round(r$coef[3], 3),
        se_conv = round(r$se[1], 3), se_robust = round(r$se[3], 3),
        ci_lo = round(r$ci[3, 1], 3), ci_hi = round(r$ci[3, 2], 3),
        h = round(h, 2),
        classes_left = r$N_h[1], classes_right = r$N_h[2],
        schools_left = sc[["schools_left"]], schools_right = sc[["schools_right"]],
        window_hits_other_cutoff = crosses_cutoff(c0, h))
    }
  }
}
frd <- do.call(rbind, frd)
print(frd)
write.csv(frd, file.path(OUT_DIR, "fuzzy_rd_by_cutoff.csv"), row.names = FALSE)



## 5b. Reduced form 


rf <- list()
for (g in names(GRADES)) {
  dd0 <- GRADES[[g]]
  for (y in c("avgverb", "avgmath")) {
    dd <- dd0[!is.na(dd0[[y]]), , drop = FALSE]
    for (c0 in CUTOFFS) {
      r <- tryCatch(
        rdrobust(y = dd[[y]], x = dd$c_size, c = c0 + 0.5,
                 cluster = dd$schlcode, kernel = "triangular", p = 1),
        error = function(e) NULL)
      if (is.null(r)) next
      h  <- r$bws[1, 1]
      sc <- school_counts(dd, c0, h)
      rf[[paste(g, y, c0)]] <- data.frame(
        grade = g, outcome = y, cutoff = c0,
        coef = round(r$coef[1], 3), coef_bc = round(r$coef[3], 3),
        se_robust = round(r$se[3], 3),
        ci_lo = round(r$ci[3, 1], 3), ci_hi = round(r$ci[3, 2], 3),
        p_robust = round(r$pv[3], 4), h = round(h, 2),
        schools_left = sc[["schools_left"]], schools_right = sc[["schools_right"]],
        window_hits_other_cutoff = crosses_cutoff(c0, h))
    }
  }
}
rf <- do.call(rbind, rf)
print(rf)
write.csv(rf, file.path(OUT_DIR, "reduced_form_by_cutoff.csv"), row.names = FALSE)



## 6. Bandwidth sensitivity at the first cutoff

bw_sens <- list()
for (g in names(GRADES)) {
  dd0 <- GRADES[[g]]
  for (y in c("avgverb", "avgmath")) {
    dd <- dd0[!is.na(dd0[[y]]), , drop = FALSE]
    for (h in c(5, 8, 10, 15, 20)) {
      r <- tryCatch(
        rdrobust(y = dd[[y]], x = dd$c_size, c = 40.5, fuzzy = dd$classize,
                 cluster = dd$schlcode, kernel = "triangular", p = 1, h = h),
        error = function(e) NULL)
      rfs <- tryCatch(
        rdrobust(y = dd$classize, x = dd$c_size, c = 40.5,
                 cluster = dd$schlcode, kernel = "triangular", p = 1, h = h),
        error = function(e) NULL)
      if (is.null(r)) next
      sc <- school_counts(dd, 40, h)
      bw_sens[[paste(g, y, h)]] <- data.frame(
        grade = g, outcome = y, h = h,
        coef = round(r$coef[1], 3), se_robust = round(r$se[3], 3),
        ci_lo = round(r$ci[3, 1], 3), ci_hi = round(r$ci[3, 2], 3),
        first_stage = if (is.null(rfs)) NA else round(rfs$coef[1], 3),
        classes = sum(r$N_h),
        schools_left = sc[["schools_left"]], schools_right = sc[["schools_right"]])
    }
  }
}
bw_sens <- do.call(rbind, bw_sens)
print(bw_sens)
write.csv(bw_sens, file.path(OUT_DIR, "bandwidth_sensitivity.csv"),
          row.names = FALSE)



## 7. Donut specifications


donut <- list()
for (g in names(GRADES)) {
  dd0 <- GRADES[[g]]
  for (y in c("avgverb", "avgmath")) {
    dda <- dd0[!is.na(dd0[[y]]), , drop = FALSE]
    for (k in 0:2) {
      sub <- dda[abs(dda$c_size - 40.5) > (k - 0.5), , drop = FALSE]
      r <- tryCatch(
        rdrobust(y = sub[[y]], x = sub$c_size, c = 40.5, fuzzy = sub$classize,
                 cluster = sub$schlcode, kernel = "triangular", p = 1, h = 10),
        error = function(e) NULL)
      rfs <- tryCatch(
        rdrobust(y = sub$classize, x = sub$c_size, c = 40.5,
                 cluster = sub$schlcode, kernel = "triangular", p = 1, h = 10),
        error = function(e) NULL)
      if (is.null(r)) next
      sc <- school_counts(sub, 40, 10)
      donut[[paste(g, y, k)]] <- data.frame(
        grade = g, outcome = y, donut = k,
        dropped = if (k == 0) "none"
        else paste0(41 - k, "-", 40 + k),
        coef = round(r$coef[1], 3), se_robust = round(r$se[3], 3),
        first_stage    = if (is.null(rfs)) NA else round(rfs$coef[1], 3),
        first_stage_se = if (is.null(rfs)) NA else round(rfs$se[3], 3),
        classes = sum(r$N_h),
        schools_left = sc[["schools_left"]], schools_right = sc[["schools_right"]])
    }
  }
}
donut <- do.call(rbind, donut)
print(donut)
write.csv(donut, file.path(OUT_DIR, "donut.csv"), row.names = FALSE)


## 8. Covariate balance and placebo cutoffs


## 8a. Balance on percent disadvantaged
bal <- list()
for (g in names(GRADES)) {
  dd <- GRADES[[g]]
  for (c0 in CUTOFFS) {
    r <- tryCatch(
      rdrobust(y = dd$tipuach, x = dd$c_size, c = c0 + 0.5,
               cluster = dd$schlcode, kernel = "triangular", p = 1),
      error = function(e) NULL)
    if (is.null(r)) next
    h  <- r$bws[1, 1]
    sc <- school_counts(dd, c0, h)
    
    mu <- mean(dd$tipuach, na.rm = TRUE)
    bal[[paste(g, c0)]] <- data.frame(
      grade = g, cutoff = c0, covariate = "percent disadvantaged",
      coef = round(r$coef[1], 3), se_robust = round(r$se[3], 3),
      p_robust = round(r$pv[3], 4),
      sample_mean = round(mu, 2),
      coef_as_pct_of_mean = round(100 * r$coef[1] / mu, 1),
      h = round(h, 2),
      schools_left = sc[["schools_left"]], schools_right = sc[["schools_right"]])
  }
}
bal <- do.call(rbind, bal)
print(bal)
write.csv(bal, file.path(OUT_DIR, "covariate_balance.csv"), row.names = FALSE)


## 8b. Placebo cutoffs.

PLACEBOS <- c(20, 25, 30, 55, 60, 65, 95, 100, 105)
placebo_fit <- function(dat, yvar, c0) {
  dd <- dat[!is.na(dat[[yvar]]), , drop = FALSE]
  r0 <- rdrobust(y = dd[[yvar]], x = dd$c_size, c = c0 + 0.5,
                 cluster = dd$schlcode, kernel = "triangular", p = 1)
  h  <- r0$bws[1, 1]
  sc <- school_counts(dd, c0, h)
  data.frame(
    placebo_cutoff = c0, outcome = yvar,
    h = round(h, 2),
    window_hits_cutoff = crosses_cutoff(c0, h),
    coef      = round(r0$coef[1], 3),
    se_robust = round(r0$se[3], 3),
    p_robust  = round(r0$pv[3], 4),
    schools_left  = sc[["schools_left"]],
    schools_right = sc[["schools_right"]])
}


placebo <- list()
for (g in names(GRADES)) {
  dd <- GRADES[[g]]
  for (c0 in PLACEBOS) {
    for (yv in c("classize", "avgverb")) {
      out <- tryCatch(placebo_fit(dd, yv, c0), error = function(e) NULL)
      if (is.null(out)) next
      placebo[[paste(g, c0, yv)]] <- cbind(grade = g, out)
    }
  }
}
placebo <- do.call(rbind, placebo)
print(placebo)
write.csv(placebo, file.path(OUT_DIR, "placebo_cutoffs.csv"), row.names = FALSE)

## Rejection counts
placebo_summary <- aggregate(
  cbind(rejects_uncapped = p_uncapped < 0.05,
        rejects_capped   = p_capped   < 0.05) ~ grade + outcome,
  data = placebo, FUN = sum)
print(placebo_summary)
write.csv(placebo_summary, file.path(OUT_DIR, "placebo_summary.csv"),
          row.names = FALSE)





## 9. Honest inference with a discrete running variable


pkg_M <- function(dd, yvar) {
  d <- dd[!is.na(dd[[yvar]]) & !is.na(dd$classize) & !is.na(dd$c_size), ,
          drop = FALSE]
  d$x  <- d$c_size - 40.5
  d$yy <- d[[yvar]]
  

  m <- tryCatch(as.numeric(RDHonest::MROT(yy | classize ~ x, data = d)),
                error = function(e) NULL)
  if (!is.null(m) && length(m) && all(is.finite(m)))
    return(list(M = m, how = "MROT() direct"))
  
 
  r <- tryCatch(RDHonest(yy | classize ~ x, data = d, kern = "triangular",
                         opt.criterion = "FLCI",
                         clusterid = d$schlcode, se.method = "EHW"),
                error = function(e) NULL)
  if (is.null(r)) return(list(M = NA_real_, how = "FAILED"))
  
  cf <- r$coefficients
  m  <- if (!is.null(cf) && "M" %in% names(cf)) as.numeric(cf$M) else
    if (!is.null(r$M))                     as.numeric(r$M)  else NA_real_
  
  ## If M comes back NA the field has been renamed; these two lines show where
  ## to look rather than leaving a silent NA in the table.
  if (all(is.na(m))) {
    message("Could not locate M. Fields available in $coefficients: ",
            paste(names(cf), collapse = ", "))
    message("Fields available at top level: ", paste(names(r), collapse = ", "))
  }
  list(M = m, how = "read back from fit")
}

mrot_tab <- list()
for (g in names(GRADES)) {
  for (y in c("avgverb", "avgmath")) {
    res <- pkg_M(GRADES[[g]], y)
    mrot_tab[[paste(g, y)]] <- data.frame(
      grade = g, outcome = y,
      M_outcome    = round(res$M[1], 4),
      M_firststage = round(if (length(res$M) > 1) res$M[2] else NA_real_, 4),
      source = res$how)
  }
}
mrot_tab <- do.call(rbind, mrot_tab)
print(mrot_tab)
write.csv(mrot_tab, file.path(OUT_DIR, "curvature_bounds.csv"), row.names = FALSE)



## clustered fit

fit_honest <- function(dat, My, Md) {
  base <- list(formula = as.formula("yy | classize ~ x"), data = dat,
               M = c(My, Md), kern = "triangular", opt.criterion = "FLCI")
  r <- tryCatch(do.call(RDHonest,
                        c(base, list(clusterid = dat$schlcode,
                                     se.method = "EHW"))),
                error = function(e) NULL)
  if (!is.null(r)) return(list(fit = r, clustered = TRUE))
  r <- tryCatch(do.call(RDHonest, base), error = function(e) NULL)
  list(fit = r, clustered = FALSE)
}

gv <- function(df, nm) {
  if (!is.null(df) && nm %in% names(df)) as.numeric(df[[nm]][1]) else NA_real_
}

honest <- list()
for (g in names(GRADES)) {
  for (y in c("avgverb", "avgmath")) {
    
    dd <- GRADES[[g]]
    dd <- dd[!is.na(dd[[y]]) & !is.na(dd$classize) & !is.na(dd$c_size), ,
             drop = FALSE]
    dd$x  <- dd$c_size - 40.5
    dd$yy <- dd[[y]]
    

    row    <- mrot_tab$grade == g & mrot_tab$outcome == y
    My_rot <- mrot_tab$M_outcome[row]
    Md_rot <- mrot_tab$M_firststage[row]
    if (is.na(My_rot)) next
    if (is.na(Md_rot)) Md_rot <- My_rot   
    
## Panel A
    for (mult in c(0.5, 1, 2)) {
      res <- fit_honest(dd, My_rot * mult, Md_rot * mult)
      if (is.null(res$fit)) next
      cf <- res$fit$coefficients
      honest[[paste(g, y, "rot", mult)]] <- data.frame(
        grade = g, outcome = y, panel = "rule of thumb",
        setting = paste0(mult, " x ROT"),
        M_outcome    = round(My_rot * mult, 4),
        M_firststage = round(Md_rot * mult, 4),
        estimate  = round(gv(cf, "estimate"), 3),
        ci_lo     = round(gv(cf, "conf.low"), 3),
        ci_hi     = round(gv(cf, "conf.high"), 3),
        bandwidth = round(gv(cf, "bandwidth"), 2),
        clustered = res$clustered)
    }
    
    
## Panel B
    for (My in c(0.01, 0.05, 0.10)) {
      res <- fit_honest(dd, My, Md_rot)
      if (is.null(res$fit)) next
      cf <- res$fit$coefficients
      honest[[paste(g, y, "fixed", My)]] <- data.frame(
        grade = g, outcome = y, panel = "fixed outcome M",
        setting = paste0("M_y = ", My),
        M_outcome = My, M_firststage = round(Md_rot, 4),
        estimate  = round(gv(cf, "estimate"), 3),
        ci_lo     = round(gv(cf, "conf.low"), 3),
        ci_hi     = round(gv(cf, "conf.high"), 3),
        bandwidth = round(gv(cf, "bandwidth"), 2),
        clustered = res$clustered)
    }
  }
}

honest <- do.call(rbind, honest)
print(honest)
write.csv(honest, file.path(OUT_DIR, "honest_ci.csv"), row.names = FALSE)

if (any(!honest$clustered)) {
  warning("RDHonest ran WITHOUT clustering for some rows. Check ?RDHonest for ",
          "the clustering argument in your version, and do not report those ",
          "intervals as clustered.")
} else {
  message("All honest intervals clustered on schools.")
}

## 12. Figures


PLOT_H <- 20

make_panel <- function(dat, yvar, c0, ylab) {
  sub <- abs(dat$c_size - c0) <= PLOT_H & !is.na(dat[[yvar]])
  p <- tryCatch(
    rdplot(y = dat[[yvar]], x = dat$c_size, c = c0 + 0.5, p = 1, h = PLOT_H,
           subset = sub, hide = TRUE,
           title = paste0(ylab, ", cutoff ", c0),
           x.label = "Grade enrolment", y.label = ylab),
    error = function(e) NULL)
  if (is.null(p)) NULL else p$rdplot
}

if (requireNamespace("cowplot", quietly = TRUE)) {
  panels <- list()
  for (c0 in CUTOFFS) {
    panels[[length(panels) + 1]] <- make_panel(d5, "classize", c0, "Class size")
    panels[[length(panels) + 1]] <- make_panel(d5, "avgverb",  c0,
                                               "Average reading score")
  }
  panels <- Filter(Negate(is.null), panels)
  if (length(panels)) {
    grid_plot <- cowplot::plot_grid(plotlist = panels, ncol = 2)
    ggplot2::ggsave(file.path(OUT_DIR, "rdplots.pdf"), grid_plot,
                    width = 9, height = 11)
  }
} else {
  message("cowplot not installed; writing rdplots one per page instead.")
  pdf(file.path(OUT_DIR, "rdplots.pdf"), width = 6, height = 4)
  for (c0 in CUTOFFS) {
    for (spec in list(c("classize", "Class size"),
                      c("avgverb", "Average reading score"))) {
      p <- make_panel(d5, spec[1], c0, spec[2])
      if (!is.null(p)) print(p)
    }
  }
  dev.off()
}

## School-level enrolment histogram around the first cutoff.

pdf(file.path(OUT_DIR, "enrolment_density.pdf"), width = 8, height = 9)
par(mfrow = c(2, 1), mar = c(4.5, 4.5, 3, 1))
LEVELS <- 20:60
for (g in names(GRADES)) {
  sch <- GRADES[[g]][!duplicated(GRADES[[g]]$schlcode), ]
  tb  <- table(factor(sch$c_size, levels = LEVELS))
  bp  <- barplot(tb, main = paste(g, "grade: schools by enrolment"),
                 xlab = "Grade enrolment", ylab = "Number of schools", las = 2)
  i40 <- which(LEVELS == 40)
  abline(v = mean(c(bp[i40], bp[i40 + 1])), lty = 2, lwd = 2)
}
dev.off()
################################################################################
################################################################################
################################################################################


# Modern  Pop-Eleches & Urquiola (2013) 

library(rdrobust)
library(rddensity)
library(haven)
library(dplyr)
library(ggplot2)
library(fixest)


data1 <- read_dta("C:/Users/Katja/Documents/BA/112645-V1/data/data-AER-1.dta")

# ---- 1. CCT robust RD estimate, main outcomes ---------------------------

rd_agus <- rdrobust(y = data1$agus, x = data1$dzag, c = 0,
                    cluster = data1$sid2)
rd_bct  <- rdrobust(y = data1$bct,  x = data1$dzag, c = 0,
                    cluster = data1$sid2)
rd_bcg  <- rdrobust(y = data1$bcg,  x = data1$dzag, c = 0,
                    cluster = data1$sid2)

summary(rd_agus)
summary(rd_bct)
summary(rd_bcg)



# ---- 2. Manipulation / sorting test at the cutoff ------------------------

dens_test <- rddensity(X = data1$dzag, c = 0)
summary(dens_test)
rdplotdensity(dens_test, data1$dzag)



dens_test_excl <- rddensity(X = data1$dzag[data1$dzag != 0], c = 0)
summary(dens_test_excl)
rdplotdensity(dens_test_excl, data1$dzag[data1$dzag != 0])


table(data1$dzag == 0)
sum(data1$dzag == 0)   


## 3.Bandwidth-sensitivity check

library(cowplot)


bw_grid <- seq(0.2, 2, by = 0.1)

bandwidth_sensitivity <- function(outcome, data = data1, grid = bw_grid) {
  sens <- lapply(grid, function(h) {
    d <- data %>% filter(dzag != 0, abs(dzag) <= h)
    fml <- as.formula(paste0(outcome, " ~ dga + dzag + dzag_after | uazY"))
    m <- feols(fml, data = d, cluster = ~sid2)
    data.frame(bw = h, coef = coef(m)["dga"], se = se(m)["dga"])
  })
  do.call(rbind, sens)
}

plot_sensitivity <- function(sens, outcome_label) {
  ggplot(sens, aes(bw, coef)) +
    geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se), alpha = .2) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Bandwidth (|dzag| <= h)", y = paste("RD estimate on", outcome_label),
         title = paste("Bandwidth sensitivity:", outcome_label))
}

sens_agus <- bandwidth_sensitivity("agus")
sens_bct  <- bandwidth_sensitivity("bct")
sens_bcg  <- bandwidth_sensitivity("bcg")

p_sens_agus <- plot_sensitivity(sens_agus, "transition score (agus)")
p_sens_bct  <- plot_sensitivity(sens_bct,  "baccalaureate taken (bct)")
p_sens_bcg  <- plot_sensitivity(sens_bcg,  "baccalaureate grade (bcg)")

figure_sensitivity <- cowplot::plot_grid(p_sens_agus, p_sens_bct, p_sens_bcg, ncol = 1)
figure_sensitivity
ggsave("figure_bandwidth_sensitivity.png", figure_sensitivity, width = 7, height = 12)


# 4. CUTOFF FIXED EFFECTS IN THE rdrobust ESTIMATES
#

resid_on_fe <- function(data, outcome, fe = "uazY") {
  d <- data %>% filter(dzag != 0, !is.na(.data[[outcome]]))
  m <- feols(as.formula(paste0(outcome, " ~ 1 | ", fe)), data = d)
  d$.resid <- resid(m, na.rm = FALSE)
  d %>% filter(!is.na(.resid))
}

rd_fe <- function(outcome, data = data1) {
  d <- resid_on_fe(data, outcome)
  rdrobust(y = d$.resid, x = d$dzag, c = 0, cluster = d$sid2)
}

rd_agus_fe <- rd_fe("agus")
rd_bct_fe  <- rd_fe("bct")
rd_bcg_fe  <- rd_fe("bcg")

summary(rd_agus_fe)
summary(rd_bct_fe)
summary(rd_bcg_fe)


# 5. THE MULTI-CUTOFF STRUCTURE


library(rdmulti)

# 5a. Cutoff-group-specific modern estimates.


rd_by_group <- function(outcome, group_expr, label, data = data1) {
  d <- data %>% filter(dzag != 0) %>% filter(!!rlang::parse_expr(group_expr))
  d <- d %>% filter(!is.na(.data[[outcome]]))
  r <- rdrobust(y = d[[outcome]], x = d$dzag, c = 0, cluster = d$sid2)
  data.frame(group = label,
             est   = r$coef["Robust", 1],
             ci_l  = r$ci["Robust", 1],
             ci_r  = r$ci["Robust", 2],
             h     = r$bws[1, 1],
             n     = r$N_h[1] + r$N_h[2])
}

groups <- list(
  c("zga > 7.74",  "Top tercile of cutoffs"),
  c("zga < 6.77",  "Bottom tercile of cutoffs"),
  c("nusua >= 4",  "Towns with 4+ schools"),
  c("nusua == 3",  "Towns with 3 schools"),
  c("nusua == 2",  "Towns with 2 schools")
)

multi_agus <- do.call(rbind, lapply(groups, function(g) rd_by_group("agus", g[1], g[2])))
multi_bcg  <- do.call(rbind, lapply(groups, function(g) rd_by_group("bcg",  g[1], g[2])))
print(multi_agus)
print(multi_bcg)









