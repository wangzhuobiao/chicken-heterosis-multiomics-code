classify_f1_pattern <- function(mean_YY, mean_YX, mean_XX,
                                parent_diff_min = PARENTAL_DIFF_MIN,
                                pos_tol = POSITION_TOLERANCE,
                                mean_tol_min = MEAN_TOLERANCE_MIN) {
  if (any(is.na(c(mean_YY, mean_YX, mean_XX)))) return("unresolved")
  lo <- min(mean_YY, mean_XX)
  hi <- max(mean_YY, mean_XX)
  range <- hi - lo
  if (!is.finite(range) || range < parent_diff_min) return("unresolved")
  tol <- max(range * pos_tol, mean_tol_min)
  if (mean_YX > hi + tol) return("YX_high_transgressive")
  if (mean_YX < lo - tol) return("YX_low_transgressive")
  if (abs(mean_YX - mean_YY) <= tol && abs(mean_YX - mean_XX) > tol) return("YY_dominant")
  if (abs(mean_YX - mean_XX) <= tol && abs(mean_YX - mean_YY) > tol) return("XX_dominant")
  if (mean_YX >= lo - tol && mean_YX <= hi + tol) return("YX_intermediate")
  "unresolved"
}

safe_numeric <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "-", "--")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

residualize <- function(y, data, covariates) {
  y <- safe_numeric(y)
  dat <- data
  dat$.response <- y
  keep_cols <- c(".response", covariates)
  complete <- complete.cases(dat[, keep_cols, drop = FALSE])
  out <- rep(NA_real_, nrow(dat))
  if (sum(complete) < 4) return(out)
  f <- as.formula(paste(".response ~", paste(covariates, collapse = " + ")))
  fit <- lm(f, data = dat[complete, , drop = FALSE])
  out[complete] <- residuals(fit)
  out
}

safe_cor_test <- function(x, y, method = "spearman") {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < MIN_COMPLETE_N || sd(x) == 0 || sd(y) == 0) {
    return(c(n = length(x), estimate = NA_real_, pvalue = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  c(n = length(x), estimate = unname(ct$estimate), pvalue = ct$p.value)
}

outcome_specs <- list(
  BW_wk15 = c("group"),
  ADG_0_15 = c("group"),
  Live_wt = c("group"),
  Carcass_wt = c("group"),
  Evis_wt = c("group"),
  AF_wt = c("Live_wt", "group"),
  AF_live_pct = c("group"),
  AF_resid_live = c("group")
)

spec <- list(BW_wk15=character(),ADG_0_15=character(),Live_wt=character(),Carcass_wt=character(),Evis_wt=character(),AF_wt='Live_wt',AF_live_pct=character(),AF_resid_live=character())

fit_one <- function(df,p,f,extra) {
 dat<-df[,unique(c('group',p,f,extra)),drop=FALSE];names(dat)[names(dat)==f]<-'factor_score'
 dat[[p]]<-safe_numeric(dat[[p]]);dat$factor_score<-safe_numeric(dat$factor_score)
 for(cc in extra) dat[[cc]]<-safe_numeric(dat[[cc]])
 dat<-dat[complete.cases(dat),,drop=FALSE]
 stopifnot(nrow(dat)>=10,sd(dat$factor_score)>0,sd(dat[[p]])>0)
 dat$phenotype_z<-as.numeric(scale(dat[[p]]));dat$factor_z<-as.numeric(scale(dat$factor_score));cov<-character()
 for(cc in extra){nm<-paste0(cc,'_z');dat[[nm]]<-as.numeric(scale(dat[[cc]]));cov<-c(cov,nm)}
 form<-as.formula(paste('phenotype_z ~',paste(c('group',cov,'factor_z'),collapse=' + ')))
 fit<-lm(form,data=dat);sm<-summary(fit);co<-sm$coefficients
 data.frame(phenotype=p,factor=f,n=nrow(dat),model=paste(deparse(form),collapse=''),standardized_beta=co['factor_z','Estimate'],SE=co['factor_z','Std. Error'],t_value=co['factor_z','t value'],pvalue=co['factor_z','Pr(>|t|)'],model_R2=sm$r.squared,adjusted_R2=sm$adj.r.squared)
}
