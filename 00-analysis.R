# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                                             #
#   Perils of Standardization: Divergences Between Change Score and Endpoint  # 
#   Effect Sizes in Meta-Analyses of Depression Psychotherapy                 #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

library(tidyverse)
library(data.table)
library(meta)
library(metafor)
library(metapsyTools)
library(foreach)
library(doParallel)
require(ellipse)
library(readxl)
library(xlsx)
library(sjPlot)
library(cowplot)

source("utils/utils.R")
source("utils/token.R") # Private PAT
source("utils/cm_conversion.R")


# 1. Univariate Meta-Analysis --------------------------------------------------

# Define rho values to test
rhos = c(.2, .4, .6, .8)

# Set up backend
cores = detectCores()
cl = makeCluster(cores[1]-1) 
registerDoParallel(cl)

sim.uni = foreach(i = 1:length(rhos), 
                  .verbose = TRUE) %dopar% {
                    
  require(tidyverse)
  require(data.table)
  require(meta)
  require(metafor)
  require(metapsyTools)
  require(foreach)
  require(doParallel)
  require(readxl)
 
  rho = rhos[i]
  
  # Import data
  dat = read_excel("data/data.xlsx")
  
  # Calculate SMD types
  calculateSMDs(dat, rho) -> dat
  
  # Define models to run
  mdls = c("overall", "combined", "lowest.highest",
           "influence", "threelevel.che")
  
  # Run models for different SMD types
  m.list = list(
    m.fo = 
      runMetaAnalysis(dat, which.run = mdls,
                      es.var = "d_fo", se.var = "d_fo_se",
                      which.influence = "combined") %>% 
      list(res = .$summary, varcomp = .$model.threelevel.che.var.comp),
    m.cs =
      runMetaAnalysis(dat, which.run = mdls,
                      es.var = "d_cs", se.var = "d_cs_se",
                      which.influence = "combined") %>% 
      list(res = .$summary, varcomp = .$model.threelevel.che.var.comp),
    m.cs.fo = 
      runMetaAnalysis(dat, which.run = mdls,
                      es.var = "d_cs_fo", se.var = "d_cs_fo_se",
                      which.influence = "combined") %>% 
      list(res = .$summary, varcomp = .$model.threelevel.che.var.comp),
    m.cs.cs = 
      runMetaAnalysis(dat, which.run = mdls,
                      es.var = "d_cs_cs", se.var = "d_cs_cs_se",
                      which.influence = "combined") %>% 
      list(res = .$summary, varcomp = .$model.threelevel.che.var.comp)
  ); m.list
}

# Stop cluster
stopCluster(cl)
save(sim.uni, file="results/sim.uni.rda")

sim.uni %>% 
  map2_dfr(as.list(rhos), function(x,y){
    x$m.fo$res$rho = y; x$m.fo$res$calc.type = "fo"
    x$m.cs$res$rho = y; x$m.cs$res$calc.type = "cs"
    x$m.cs.fo$res$rho = y; x$m.cs.fo$res$calc.type = "cs.fo"
    x$m.cs.cs$res$rho = y; x$m.cs.cs$res$calc.type = "cs.cs"
    x$m.fo$res$model = rownames(x$m.fo$res)
    x$m.cs$res$model = rownames(x$m.cs$res)
    x$m.cs.fo$res$model = rownames(x$m.cs.fo$res)
    x$m.cs.cs$res$model = rownames(x$m.cs.cs$res)
    rbind(x$m.fo$res, x$m.cs$res, x$m.cs.fo$res, x$m.cs.cs$res) %>% 
      as_tibble() %>% select(model, calc.type, rho, everything())
  }) %>% arrange(model, rho, desc(calc.type)) %>% 
  write.xlsx("results/univariate.xlsx", "effect")

sim.uni %>% 
  map2_dfr(as.list(rhos), function(x,y){
    x$m.cs$varcomp$model = "Three-Level Model (CHE)"
    x$m.cs$varcomp$calc.type = "cs"
    x$m.cs$varcomp$rho = y
    x$m.cs$varcomp$lvl = rownames(x$m.cs$varcomp)
    x$m.fo$varcomp$model = "Three-Level Model (CHE)"
    x$m.fo$varcomp$calc.type = "fo"
    x$m.fo$varcomp$rho = y
    x$m.fo$varcomp$lvl = rownames(x$m.fo$varcomp)
    x$m.cs.fo$varcomp$model = "Three-Level Model (CHE)"
    x$m.cs.fo$varcomp$calc.type = "cs.fo"
    x$m.cs.fo$varcomp$rho = y
    x$m.cs.fo$varcomp$lvl = rownames(x$m.cs.fo$varcomp)
    x$m.cs.cs$varcomp$model = "Three-Level Model (CHE)"
    x$m.cs.cs$varcomp$calc.type = "cs.cs"
    x$m.cs.cs$varcomp$rho = y
    x$m.cs.cs$varcomp$lvl = rownames(x$m.cs.cs$varcomp)
    rbind(x$m.cs$varcomp, x$m.fo$varcomp,
          x$m.cs.cs$varcomp, x$m.cs.fo$varcomp) %>% 
      as_tibble() %>% arrange(model, rho, desc(calc.type)) %>% 
      select(model, rho, calc.type, lvl, everything())
  }) %>% 
  write.xlsx("results/univariate.xlsx", "tau", append=T)



# 2. Bivariate Meta-Analysis ---------------------------------------------------

## 2.1 Run Simulation ----------------------------------------------------------

# Define rho values to test
rhos = c(.2, .4, .6, .8)

# Set up backend
cores = detectCores()
cl = makeCluster(cores[1]-1) 
registerDoParallel(cl)

sim.bv = foreach(i = 1:length(rhos), 
                  .verbose = TRUE) %dopar% {

    require(tidyverse)
    require(data.table)
    require(meta)
    require(metafor)
    require(metapsyTools)
    require(foreach)
    require(doParallel)
    require(readxl)
    
    rho = rhos[i]
    
    # Import data
    dat = read_excel("data/data.xlsx")
    
    # Calculate SMD types
    calculateSMDs(dat, rho) -> dat
    
    # Get combined SMD values
    m.fo = runMetaAnalysis(dat, "combined", es.var = "d_fo", se.var = "d_fo_se")
    m.cs = runMetaAnalysis(dat, "combined", es.var = "d_cs", se.var = "d_cs_se")
    m.cs.fo = runMetaAnalysis(dat, "combined", es.var = "d_cs_fo", se.var = "d_cs_fo_se") 
    m.cs.cs = runMetaAnalysis(dat, "combined",es.var = "d_cs_cs", se.var = "d_cs_cs_se")
    
    # Define ES indicator variable
    m.fo$model.combined$data$calc.type = "fo"
    m.cs$model.combined$data$calc.type = "cs"
    m.cs.fo$model.combined$data$calc.type = "cs"
    m.cs.cs$model.combined$data$calc.type = "cs"
    
    # Set up bivariate datasets
    as.data.frame(
      rbind(m.fo$model.combined$data,
            m.cs$model.combined$data)) %>% 
      arrange(study) %>% mutate(
        V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
        mean_age = parse_number(mean_age),
        percent_women = parse_number(percent_women) %>% 
          ifelse(.>1, ./100,.)) -> data.mv.pre
    
    as.data.frame(
      rbind(m.fo$model.combined$data,
            m.cs.fo$model.combined$data)) %>% 
      arrange(study) %>% mutate(
        V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
        mean_age = parse_number(mean_age),
        percent_women = parse_number(percent_women) %>% 
          ifelse(.>1, ./100,.)) -> data.mv.fo
    
    as.data.frame(
      rbind(m.fo$model.combined$data,
            m.cs.cs$model.combined$data)) %>% 
      arrange(study) %>% mutate(
        V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
        mean_age = parse_number(mean_age),
        percent_women = parse_number(percent_women) %>% 
          ifelse(.>1, ./100,.)) -> data.mv.cs
    
    # Run multivariate meta-analyses
    m.mv.pre = rma.mv(.TE ~ 0 + calc.type, V, random = ~ calc.type | study, 
                      struct = "UN", data = data.mv.pre)
    rg.mv.pre = matreg(y=1, x=2, R=m.mv.pre$G, cov=TRUE, 
                       means=coef(m.mv.pre), n=m.mv.pre$g.levels.comb.k)
    re.mv.pre = ranef(m.mv.pre)
    bl.mv.pre = data.frame(study = substr(re.mv.pre$`~calc.type | study` %>% 
                                            rownames, start = 6, stop = 1e5),
                           type = substr(re.mv.pre$`~calc.type | study` %>% 
                                           rownames, start = 1, stop = 2),
                           es = re.mv.pre$`~calc.type | study`$intrcpt,
                           se = re.mv.pre$`~calc.type | study`$se) %>% 
                           pivot_wider(id_cols = study, names_from = type, 
                                       values_from = c(es, se)) %>% 
                           mutate(es_cs = es_cs + coef(m.mv.pre)[1], 
                                  es_fo = es_fo + coef(m.mv.pre)[2])
    mv.pre = list(m = m.mv.pre, rg = rg.mv.pre, re = re.mv.pre, bl = bl.mv.pre)
     
    
    m.mv.fo = rma.mv(.TE ~ 0 + calc.type, V, random = ~ calc.type | study, 
                      struct = "UN", data = data.mv.fo)
    rg.mv.fo = matreg(y=1, x=2, R=m.mv.fo$G, cov=TRUE, 
                       means=coef(m.mv.fo), n=m.mv.fo$g.levels.comb.k)
    re.mv.fo = ranef(m.mv.fo)
    bl.mv.fo = data.frame(study = substr(re.mv.fo$`~calc.type | study` %>% 
                                            rownames, start = 6, stop = 1e5),
                           type = substr(re.mv.fo$`~calc.type | study` %>% 
                                           rownames, start = 1, stop = 2),
                           es = re.mv.fo$`~calc.type | study`$intrcpt,
                           se = re.mv.fo$`~calc.type | study`$se) %>% 
      pivot_wider(id_cols = study, names_from = type, 
                  values_from = c(es, se)) %>% 
      mutate(es_cs = es_cs + coef(m.mv.fo)[1], 
             es_fo = es_fo + coef(m.mv.fo)[2])
    mv.fo = list(m = m.mv.fo, rg = rg.mv.fo, re = re.mv.fo, bl = bl.mv.fo)
    
    
    m.mv.cs = rma.mv(.TE ~ 0 + calc.type, V, random = ~ calc.type | study, 
                      struct = "UN", data = data.mv.cs)
    rg.mv.cs = matreg(y=1, x=2, R=m.mv.cs$G, cov=TRUE, 
                       means=coef(m.mv.cs), n=m.mv.cs$g.levels.comb.k)
    re.mv.cs = ranef(m.mv.cs)
    bl.mv.cs = data.frame(study = substr(re.mv.cs$`~calc.type | study` %>% 
                                            rownames, start = 6, stop = 1e5),
                           type = substr(re.mv.cs$`~calc.type | study` %>% 
                                           rownames, start = 1, stop = 2),
                           es = re.mv.cs$`~calc.type | study`$intrcpt,
                           se = re.mv.cs$`~calc.type | study`$se) %>% 
      pivot_wider(id_cols = study, names_from = type, 
                  values_from = c(es, se)) %>% 
      mutate(es_cs = es_cs + coef(m.mv.cs)[1], 
             es_fo = es_fo + coef(m.mv.cs)[2])
    mv.cs = list(m = m.mv.cs, rg = rg.mv.cs, re = re.mv.cs, bl = bl.mv.cs)
    
    list(mv.pre = mv.pre, mv.fo = mv.fo, mv.cs = mv.cs)
}

# Stop cluster
stopCluster(cl)
save(sim.bv, file="results/sim.bv.rda")


## 2.2 Draw Bivariate Plot -----------------------------------------------------
pdf("results/plot_slopes.pdf", width=9.5, height = 3.6)

par(mfrow = c(1, 3))
cols = c("dodgerblue1", "dodgerblue2", 
         "dodgerblue3", "dodgerblue4")

plot(1, type="n", xlab=expression(SMD[EP/EP]), ylab=expression(SMD[CS/BL]), xlim=c(-1, 5), ylim=c(-1, 5))
title(expression(SMD[cs]~"standardized by baseline SD"))
polygon(c(-1.5, 5.5, 5.5), c(-1.5, 5.5, -1.5), col = "gray95", border=NA)
for (i in 1:length(rhos)) {
  x = sim.bv[[i]]; col = cols[i]
  points(x$mv.pre$bl[,3][[1]], x$mv.pre$bl[,2][[1]], cex=(200+(x$mv.pre$bl[,5][[1]]^2)^-1)/400, pch=19, col="lightgray")
  points(x$mv.pre$bl[,3][[1]], x$mv.pre$bl[,2][[1]], cex=(200+(x$mv.pre$bl[,5][[1]]^2)^-1)/400)
  abline(a=x$mv.pre$rg$tab$beta[1], b=x$mv.pre$rg$tab$beta[2], lwd=3, col=col)
  points(coef(x$mv.pre$m)[2], coef(x$mv.pre$m)[1], pch=19, cex=2, col = col)
}
abline(v=0, lwd=1, lty = 3)
abline(h=0, lwd=1, lty = 3)
abline(a=0, b=1, lwd=1, lty = 3)


plot(1, type="n", xlab=expression(SMD[EP/EP]), ylab=expression(SMD[CS/CS]), xlim=c(-1, 5), ylim=c(-1, 5))
title(expression(SMD[cs]~"standardized by change score SD"))
polygon(c(-1.5, 5.5, 5.5), c(-1.5, 5.5, -1.5), col = "gray95", border=NA)
for (i in 1:length(rhos)) {
  x = sim.bv[[i]]; col = cols[i]
  points(x$mv.cs$bl[,3][[1]], x$mv.cs$bl[,2][[1]], cex=(200+(x$mv.cs$bl[,5][[1]]^2)^-1)/400, pch=19, col="lightgray")
  points(x$mv.cs$bl[,3][[1]], x$mv.cs$bl[,2][[1]], cex=(200+(x$mv.cs$bl[,5][[1]]^2)^-1)/400)
  abline(a=x$mv.cs$rg$tab$beta[1], b=x$mv.cs$rg$tab$beta[2], lwd=3, col=col)
  points(coef(x$mv.cs$m)[2], coef(x$mv.cs$m)[1], pch=19, cex=2, col=col)
}
abline(v=0, lwd=1, lty = 3)
abline(h=0, lwd=1, lty = 3)
abline(a=0, b=1, lwd=1, lty = 3)
legend("bottomright", legend=c("r = 0.2", "r = 0.4", "r = 0.6", "r = 0.8"), 
       col=cols[1:length(rhos)], pch=15, cex=0.9, bty="n",
       x.intersp=1, y.intersp=0.5, horiz=TRUE)


plot(1, type="n", xlab=expression(SMD[EP/EP]), ylab=expression(SMD[CS/EP]), xlim=c(-1, 5), ylim=c(-1, 5))
title(expression(SMD[cs]~"standardized by endpoint SD"))
polygon(c(-1.5, 5.5, 5.5), c(-1.5, 5.5, -1.5), col = "gray95", border=NA)
for (i in 1:length(rhos)) {
  x = sim.bv[[i]]; col = cols[i]
  points(x$mv.fo$bl[,3][[1]], x$mv.fo$bl[,2][[1]], cex=(200+(x$mv.fo$bl[,5][[1]]^2)^-1)/400, pch=19, col="lightgray")
  points(x$mv.fo$bl[,3][[1]], x$mv.fo$bl[,2][[1]], cex=(200+(x$mv.fo$bl[,5][[1]]^2)^-1)/400)
  abline(a=x$mv.fo$rg$tab$beta[1], b=x$mv.fo$rg$tab$beta[2], lwd=3, col=col)
  points(coef(x$mv.fo$m)[2], coef(x$mv.fo$m)[1], pch=19, cex=2, col=col)
}
abline(v=0, lwd=1, lty = 3)
abline(h=0, lwd=1, lty = 3)
abline(a=0, b=1, lwd=1, lty = 3)

dev.off()



# 3. Meta-Regression -----------------------------------------------------------

## 3.1 Calculate Data ----------------------------------------------------------

# Import data & utils
dat = read_excel("data/data.xlsx")
source("utils/utils.R")

# Calculate p_dropout
within(dat, {
  p_dropout_arm1 = attr_arm1/rand_arm1
  p_dropout_arm2 = attr_arm2/rand_arm2
}) -> dat

# Calculate baseline difference effect size
within(dat, {
  smd.bl = (baseline_m_arm1 - baseline_m_arm2)/
               sqrt(((baseline_n_arm1-1)*baseline_sd_arm1^2 + 
                 (baseline_n_arm2-1)*baseline_sd_arm2^2)/
                (baseline_n_arm1 + baseline_n_arm2 - 2))
}) -> dat

# Calculate SMD types, assuming rho=0.2 to 0.8
calculateSMDs(dat, rho=.2) -> dat.rho.lo
calculateSMDs(dat, rho=.4) -> dat.rho.midlo
calculateSMDs(dat, rho=.6) -> dat.rho.midhi
calculateSMDs(dat, rho=.8) -> dat.rho.hi


## 3.2 Run Analyses ------------------------------------------------------------

### 3.2.1 Assuming Low Rho -----------------------------------------------------

# Get combined SMD values
m.fo = runMetaAnalysis(dat.rho.lo, "combined", es.var = "d_fo", se.var = "d_fo_se")
m.cs = runMetaAnalysis(dat.rho.lo, "combined", es.var = "d_cs", se.var = "d_cs_se")
m.cs.fo = runMetaAnalysis(dat.rho.lo, "combined", es.var = "d_cs_fo", se.var = "d_cs_fo_se") 
m.cs.cs = runMetaAnalysis(dat.rho.lo, "combined",es.var = "d_cs_cs", se.var = "d_cs_cs_se")

# Define ES indicator variable
m.fo$model.combined$data$calc.type = "fo"
m.cs$model.combined$data$calc.type = "cs"
m.cs.fo$model.combined$data$calc.type = "cs"
m.cs.cs$model.combined$data$calc.type = "cs"

# Set up bivariate datasets
as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.pre

as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs.fo$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.fo

as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs.cs$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.cs


# Run multivariate meta-analyses
rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractSubgroup()
) -> res.mv.pre

rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractSubgroup()
) -> res.mv.fo

rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractSubgroup()
) -> res.mv.cs

res.mv.pre$type = "pre"
res.mv.fo$type = "fo"
res.mv.cs$type = "cs"

rbind(res.mv.pre, res.mv.fo, res.mv.cs) %>% mutate(rho = 0.2) %>% 
  write.xlsx("results/regression/mr_low.xlsx", showNA = FALSE)


  

### 3.2.2 Assuming Mid-Low Rho -------------------------------------------------

# Get combined SMD values
m.fo = runMetaAnalysis(dat.rho.midlo, "combined", es.var = "d_fo", se.var = "d_fo_se")
m.cs = runMetaAnalysis(dat.rho.midlo, "combined", es.var = "d_cs", se.var = "d_cs_se")
m.cs.fo = runMetaAnalysis(dat.rho.midlo, "combined", es.var = "d_cs_fo", se.var = "d_cs_fo_se") 
m.cs.cs = runMetaAnalysis(dat.rho.midlo, "combined",es.var = "d_cs_cs", se.var = "d_cs_cs_se")

# Define ES indicator variable
m.fo$model.combined$data$calc.type = "fo"
m.cs$model.combined$data$calc.type = "cs"
m.cs.fo$model.combined$data$calc.type = "cs"
m.cs.cs$model.combined$data$calc.type = "cs"

# Set up bivariate datasets
as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.pre

as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs.fo$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.fo

as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs.cs$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.cs


# Run multivariate meta-analyses
rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractSubgroup()
) -> res.mv.pre

rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractSubgroup()
) -> res.mv.fo

rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractSubgroup()
) -> res.mv.cs

res.mv.pre$type = "pre"
res.mv.fo$type = "fo"
res.mv.cs$type = "cs"

rbind(res.mv.pre, res.mv.fo, res.mv.cs) %>% mutate(rho = 0.4) %>% 
  write.xlsx("results/regression/mr_midlow.xlsx", showNA = FALSE)


### 3.2.3 Assuming Mid-High Rho ------------------------------------------------

# Get combined SMD values
m.fo = runMetaAnalysis(dat.rho.midhi, "combined", es.var = "d_fo", se.var = "d_fo_se")
m.cs = runMetaAnalysis(dat.rho.midhi, "combined", es.var = "d_cs", se.var = "d_cs_se")
m.cs.fo = runMetaAnalysis(dat.rho.midhi, "combined", es.var = "d_cs_fo", se.var = "d_cs_fo_se") 
m.cs.cs = runMetaAnalysis(dat.rho.midhi, "combined",es.var = "d_cs_cs", se.var = "d_cs_cs_se")

# Define ES indicator variable
m.fo$model.combined$data$calc.type = "fo"
m.cs$model.combined$data$calc.type = "cs"
m.cs.fo$model.combined$data$calc.type = "cs"
m.cs.cs$model.combined$data$calc.type = "cs"

# Set up bivariate datasets
as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.pre

as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs.fo$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.fo

as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs.cs$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.cs


# Run multivariate meta-analyses
rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractSubgroup()
) -> res.mv.pre

rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractSubgroup()
) -> res.mv.fo

rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractSubgroup()
) -> res.mv.cs

res.mv.pre$type = "pre"
res.mv.fo$type = "fo"
res.mv.cs$type = "cs"

rbind(res.mv.pre, res.mv.fo, res.mv.cs) %>% mutate(rho = 0.6) %>% 
  write.xlsx("results/regression/mr_midhigh.xlsx", showNA = FALSE)



### 3.2.4 Assuming High Rho ----------------------------------------------------

# Get combined SMD values
m.fo = runMetaAnalysis(dat.rho.hi, "combined", es.var = "d_fo", se.var = "d_fo_se")
m.cs = runMetaAnalysis(dat.rho.hi, "combined", es.var = "d_cs", se.var = "d_cs_se")
m.cs.fo = runMetaAnalysis(dat.rho.hi, "combined", es.var = "d_cs_fo", se.var = "d_cs_fo_se") 
m.cs.cs = runMetaAnalysis(dat.rho.hi, "combined",es.var = "d_cs_cs", se.var = "d_cs_cs_se")

# Define ES indicator variable
m.fo$model.combined$data$calc.type = "fo"
m.cs$model.combined$data$calc.type = "cs"
m.cs.fo$model.combined$data$calc.type = "cs"
m.cs.cs$model.combined$data$calc.type = "cs"

# Set up bivariate datasets
as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.pre

as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs.fo$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.fo

as.data.frame(
  rbind(m.fo$model.combined$data,
        m.cs.cs$model.combined$data)) %>% 
  arrange(study) %>% mutate(
    V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
    mean_age = parse_number(mean_age),
    percent_women = parse_number(percent_women) %>% 
      ifelse(.>1, ./100,.)) -> data.mv.cs


# Run multivariate meta-analyses
rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.pre) %>% extractSubgroup()
) -> res.mv.pre

rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.fo) %>% extractSubgroup()
) -> res.mv.fo

rbind(
  rma.mv(.TE ~ 0 + calc.type*scale(rob), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("rob"),
  rma.mv(.TE ~ 0 + calc.type*scale(p_dropout_arm1), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("dropout_int"),
  rma.mv(.TE ~ 0 + calc.type*scale(abs(smd.bl)), V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractRegression("smd_bl"),
  rma.mv(.TE ~ 0 + calc.type:condition_arm1, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractSubgroup(),
  rma.mv(.TE ~ 0 + calc.type:condition_arm2, V, random = ~ calc.type | study, 
         struct = "UN", data = data.mv.cs) %>% extractSubgroup()
) -> res.mv.cs

res.mv.pre$type = "pre"
res.mv.fo$type = "fo"
res.mv.cs$type = "cs"

rbind(res.mv.pre, res.mv.fo, res.mv.cs) %>% mutate(rho = 0.8) %>% 
  write.xlsx("results/regression/mr_high.xlsx", showNA = FALSE)


## 3.3 Generate Plots ----------------------------------------------------------

rbind(
  read_excel("results/regression/mr_low.xlsx"),
  read_excel("results/regression/mr_midlow.xlsx"),
  read_excel("results/regression/mr_midhigh.xlsx"),
  read_excel("results/regression/mr_high.xlsx")
) %>% select(-1, -2) -> dat.plot

write.xlsx(dat.plot, file="results/regression/plot_dat.xlsx")
blues = c("dodgerblue1", "dodgerblue2", "dodgerblue3", "dodgerblue4")

dat.plot %>% 
  filter(fo %in% c("rob", "dropout_int", "smd_bl", "cbt")) %>% 
  mutate(type = recode(type, "cs" = "b_cs", "pre" = "a_pre", "fo"= "c_fo"),
         fo = recode(fo, "rob" = "Risk of Bias", "dropout_int" = "Attrition", 
                     "smd_bl" = "Baseline Imbalance", "cbt" = "zzz")) %>% 
  ggplot(aes(x = type, y = estimate, group = rho, color = as.factor(rho))) +
  geom_hline(yintercept = 0) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=estimate-qnorm(.975)*SE, 
                    ymax=estimate+qnorm(.975)*SE), 
                position=position_dodge(0.5), width=.2) +
  ylim(c(-.2,.2)) + ylab(expression(Delta[SMD])) + 
  xlab("") +
  scale_x_discrete(labels = c(expression(SMD[CS/BL]), 
                              expression(SMD[CS/CS]), expression(SMD[CS/EP]))) +
  theme_sjplot() + theme(legend.position = "none",
                         axis.text.x=element_text(size=8)) +
  scale_color_manual(values = blues) +
  facet_grid(cols = vars(fo)) -> p1


dat.plot %>% 
  filter(fo %in% c("cau", "other ctr", "wl", "cbt")) %>% 
  mutate(type = recode(type, "cs" = "b_cs", "pre" = "a_pre", "fo"= "c_fo"),
         fo = recode(fo, "cau" = "Care As Usual", "other ctr" = "Other Control", 
                     "wl" = "Waitlist", "cbt" = "zzz")) %>% 
  ggplot(aes(x = type, y = estimate, group = rho, color = as.factor(rho))) +
  geom_hline(yintercept = 0) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=estimate-qnorm(.975)*SE, 
                    ymax=estimate+qnorm(.975)*SE), 
                position=position_dodge(0.5), width=.2) +
  ylim(c(-.75,.75)) + ylab(expression(Delta[SMD])) + 
  xlab("") +
  scale_x_discrete(labels = c(expression(SMD[CS/BL]), 
                              expression(SMD[CS/CS]), expression(SMD[CS/EP]))) +
  theme_sjplot() + theme(legend.position = "none",
                         axis.text.x=element_text(size=8)) +
  scale_color_manual(values = blues) +
  facet_grid(cols = vars(fo)) -> p2

dat.plot %>% 
  filter(fo %in% c("3rd", "bat", "cbt", "dyn")) %>% 
  mutate(type = recode(type, "cs" = "b_cs", "pre" = "a_pre", "fo"= "c_fo"),
         fo = recode(fo, "cbt" = "CBT", "bat" = "Behavioral Activation", 
                     "3rd" = "3rd Wave CBT", "dyn" = "Psychodynamic Therapy")) %>% 
  ggplot(aes(x = type, y = estimate, group = rho, color = as.factor(rho))) +
    geom_hline(yintercept = 0) +
    geom_point(position=position_dodge(0.5)) +
    geom_errorbar(aes(ymin=estimate-qnorm(.975)*SE, 
                      ymax=estimate+qnorm(.975)*SE), 
                  position=position_dodge(0.5), width=.2) +
    ylim(c(-.75,.75)) + ylab(expression(Delta[SMD])) + 
    xlab("") +
    scale_x_discrete(labels = c(expression(SMD[CS/BL]), 
                                expression(SMD[CS/CS]), expression(SMD[CS/EP]))) +
    theme_sjplot() + theme(legend.position = "none",
                           axis.text.x=element_text(size=8)) +
  scale_color_manual(values = blues) +
  facet_grid(cols = vars(fo)) -> p3

dat.plot %>% 
  filter(fo %in% c("ipt", "lrt", "pst", "sup")) %>% 
  mutate(type = recode(type, "cs" = "b_cs", "pre" = "a_pre", "fo"= "c_fo"),
         fo = recode(fo, "ipt" = "Interpersonal Therapy", "pst" = "Problem-Solving", 
                     "sup" = "Supportive Counseling", "lrt" = "Life Review Therapy")) %>% 
  ggplot(aes(x = type, y = estimate, group = rho, color = as.factor(rho))) +
  geom_hline(yintercept = 0) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=estimate-qnorm(.975)*SE, 
                    ymax=estimate+qnorm(.975)*SE), 
                position=position_dodge(0.5), width=.2) +
  ylim(c(-.75,.75)) + ylab(expression(Delta[SMD])) + 
  xlab(expression("SD used to standardize"~SMD[CS])) +
  scale_x_discrete(labels = c(expression(SMD[CS/BL]), 
                              expression(SMD[CS/CS]), expression(SMD[CS/EP]))) +
  theme_sjplot() + theme(axis.text.x=element_text(size=8),
                         legend.position = "none") +
  facet_grid(cols = vars(fo)) +
  scale_color_manual(values = blues) -> p4
    

plot_grid(p1, p2, p3, p4, nrow=4, ncol=1) -> p
ggsave("results/regression/plot.png", bg = "white",
       width = 11, height = 11)




# 4. Effects compared to ADM ---------------------------------------------------

## 4.1 Univariate Meta-Analysis ------------------------------------------------

# Define rho values to test
rhos = c(.2, .4, .6, .8)

# Set up backend
cores = detectCores()
cl = makeCluster(cores[1]-1) 
registerDoParallel(cl)

sim.uni = foreach(i = 1:length(rhos), 
                  .verbose = TRUE) %dopar% {
                    
    require(tidyverse)
    require(data.table)
    require(meta)
    require(metafor)
    require(metapsyTools)
    require(foreach)
    require(doParallel)
    require(readxl)
    
    rho = rhos[i]
    
    # Import data
    dat = read_excel("data/data.pha.xlsx") %>% 
      filter(comparison %in% c("psy vs. Pha", "psy vs. pha"))
    
    # Calculate SMD types
    calculateSMDs(dat, rho) %>% 
      drop_na(d_fo, d_fo_se, d_cs, d_cs_se,
              d_cs_fo, d_cs_fo_se, d_cs_cs,
              d_cs_cs_se) -> dat
    
    # Define models to run
    mdls = c("overall", "combined", "lowest.highest",
             "threelevel.che", "influence")
    
    # Run models for different SMD types
    m.list = list(
      m.fo = 
        runMetaAnalysis(dat, which.run = mdls,
                        es.var = "d_fo", se.var = "d_fo_se",
                        which.influence = "combined") %>% 
        list(res = .$summary, varcomp = .$model.threelevel.che.var.comp),
      m.cs =
        runMetaAnalysis(dat, which.run = mdls,
                        es.var = "d_cs", se.var = "d_cs_se",
                        which.influence = "combined") %>% 
        list(res = .$summary, varcomp = .$model.threelevel.che.var.comp),
      m.cs.fo = 
        runMetaAnalysis(dat, which.run = mdls,
                        es.var = "d_cs_fo", se.var = "d_cs_fo_se",
                        which.influence = "combined") %>% 
        list(res = .$summary, varcomp = .$model.threelevel.che.var.comp),
      m.cs.cs = 
        runMetaAnalysis(dat, which.run = mdls,
                        es.var = "d_cs_cs", se.var = "d_cs_cs_se",
                        which.influence = "combined") %>% 
        list(res = .$summary, varcomp = .$model.threelevel.che.var.comp)
    ); m.list
}

# Stop cluster
stopCluster(cl)
save(sim.uni, file="results/sim.uni.pha.rda")

sim.uni %>% 
  map2_dfr(as.list(rhos), function(x,y){
    x$m.fo$res$rho = y; x$m.fo$res$calc.type = "fo"
    x$m.cs$res$rho = y; x$m.cs$res$calc.type = "cs"
    x$m.cs.fo$res$rho = y; x$m.cs.fo$res$calc.type = "cs.fo"
    x$m.cs.cs$res$rho = y; x$m.cs.cs$res$calc.type = "cs.cs"
    x$m.fo$res$model = rownames(x$m.fo$res)
    x$m.cs$res$model = rownames(x$m.cs$res)
    x$m.cs.fo$res$model = rownames(x$m.cs.fo$res)
    x$m.cs.cs$res$model = rownames(x$m.cs.cs$res)
    rbind(x$m.fo$res, x$m.cs$res, x$m.cs.fo$res, x$m.cs.cs$res) %>% 
      as_tibble() %>% select(model, calc.type, rho, everything())
  }) %>% arrange(model, rho, desc(calc.type)) %>% 
  write.xlsx("results/univariate.pha.xlsx", "effect")

sim.uni %>% 
  map2_dfr(as.list(rhos), function(x,y){
    x$m.cs$varcomp$model = "Three-Level Model (CHE)"
    x$m.cs$varcomp$calc.type = "cs"
    x$m.cs$varcomp$rho = y
    x$m.cs$varcomp$lvl = rownames(x$m.cs$varcomp)
    x$m.fo$varcomp$model = "Three-Level Model (CHE)"
    x$m.fo$varcomp$calc.type = "fo"
    x$m.fo$varcomp$rho = y
    x$m.fo$varcomp$lvl = rownames(x$m.fo$varcomp)
    x$m.cs.fo$varcomp$model = "Three-Level Model (CHE)"
    x$m.cs.fo$varcomp$calc.type = "cs.fo"
    x$m.cs.fo$varcomp$rho = y
    x$m.cs.fo$varcomp$lvl = rownames(x$m.cs.fo$varcomp)
    x$m.cs.cs$varcomp$model = "Three-Level Model (CHE)"
    x$m.cs.cs$varcomp$calc.type = "cs.cs"
    x$m.cs.cs$varcomp$rho = y
    x$m.cs.cs$varcomp$lvl = rownames(x$m.cs.cs$varcomp)
    rbind(x$m.cs$varcomp, x$m.fo$varcomp,
          x$m.cs.cs$varcomp, x$m.cs.fo$varcomp) %>% 
      as_tibble() %>% arrange(model, rho, desc(calc.type)) %>% 
      select(model, rho, calc.type, lvl, everything())
  }) %>% 
  write.xlsx("results/univariate.pha.xlsx", "tau", append=T)


## 4.2 Bivariate Meta-Analysis -------------------------------------------------

### 4.2.1 Run Simulation -------------------------------------------------------

# Define rho values to test
rhos = c(.2, .4, .6, .8)

# Set up backend
cores = detectCores()
cl = makeCluster(cores[1]-1) 
registerDoParallel(cl)

sim.bv = foreach(i = 1:length(rhos), 
                 .verbose = TRUE) %dopar% {
                   
  require(tidyverse)
  require(data.table)
  require(meta)
  require(metafor)
  require(metapsyTools)
  require(foreach)
  require(doParallel)
  require(readxl)
  
  rho = rhos[i]
  
  # Import data
  dat = read_excel("data/data.pha.xlsx") %>% 
    filter(comparison %in% c("psy vs. Pha", "psy vs. pha"))
  
  # Calculate SMD types
  calculateSMDs(dat, rho) %>% 
    drop_na(d_fo, d_fo_se, d_cs, d_cs_se,
            d_cs_fo, d_cs_fo_se,
            d_cs_cs, d_cs_cs_se) -> dat
  
  # Get combined SMD values
  m.fo = runMetaAnalysis(dat, "combined", es.var = "d_fo", se.var = "d_fo_se")
  m.cs = runMetaAnalysis(dat, "combined", es.var = "d_cs", se.var = "d_cs_se")
  m.cs.fo = runMetaAnalysis(dat, "combined", es.var = "d_cs_fo", se.var = "d_cs_fo_se") 
  m.cs.cs = runMetaAnalysis(dat, "combined",es.var = "d_cs_cs", se.var = "d_cs_cs_se")
  
  # Define ES indicator variable
  m.fo$model.combined$data$calc.type = "fo"
  m.cs$model.combined$data$calc.type = "cs"
  m.cs.fo$model.combined$data$calc.type = "cs"
  m.cs.cs$model.combined$data$calc.type = "cs"
  
  # Set up bivariate datasets
  as.data.frame(
    rbind(m.fo$model.combined$data,
          m.cs$model.combined$data)) %>% 
    arrange(study) %>% mutate(
      V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
      mean_age = parse_number(mean_age),
      percent_women = parse_number(percent_women) %>% 
        ifelse(.>1, ./100,.)) -> data.mv.pre
  
  as.data.frame(
    rbind(m.fo$model.combined$data,
          m.cs.fo$model.combined$data)) %>% 
    arrange(study) %>% mutate(
      V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
      mean_age = parse_number(mean_age),
      percent_women = parse_number(percent_women) %>% 
        ifelse(.>1, ./100,.)) -> data.mv.fo
  
  as.data.frame(
    rbind(m.fo$model.combined$data,
          m.cs.cs$model.combined$data)) %>% 
    arrange(study) %>% mutate(
      V = .seTE^2, condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
      mean_age = parse_number(mean_age),
      percent_women = parse_number(percent_women) %>% 
        ifelse(.>1, ./100,.)) -> data.mv.cs
  
  # Run multivariate meta-analyses
  m.mv.pre = rma.mv(.TE ~ 0 + calc.type, V, random = ~ calc.type | study, 
                    struct = "UN", data = data.mv.pre)
  rg.mv.pre = matreg(y=1, x=2, R=m.mv.pre$G, cov=TRUE, 
                     means=coef(m.mv.pre), n=m.mv.pre$g.levels.comb.k)
  re.mv.pre = ranef(m.mv.pre)
  bl.mv.pre = data.frame(study = substr(re.mv.pre$`~calc.type | study` %>% 
                                          rownames, start = 6, stop = 1e5),
                         type = substr(re.mv.pre$`~calc.type | study` %>% 
                                         rownames, start = 1, stop = 2),
                         es = re.mv.pre$`~calc.type | study`$intrcpt,
                         se = re.mv.pre$`~calc.type | study`$se) %>% 
    pivot_wider(id_cols = study, names_from = type, 
                values_from = c(es, se)) %>% 
    mutate(es_cs = es_cs + coef(m.mv.pre)[1], 
           es_fo = es_fo + coef(m.mv.pre)[2])
  mv.pre = list(m = m.mv.pre, rg = rg.mv.pre, re = re.mv.pre, bl = bl.mv.pre)
  
  
  m.mv.fo = rma.mv(.TE ~ 0 + calc.type, V, random = ~ calc.type | study, 
                   struct = "UN", data = data.mv.fo)
  rg.mv.fo = matreg(y=1, x=2, R=m.mv.fo$G, cov=TRUE, 
                    means=coef(m.mv.fo), n=m.mv.fo$g.levels.comb.k)
  re.mv.fo = ranef(m.mv.fo)
  bl.mv.fo = data.frame(study = substr(re.mv.fo$`~calc.type | study` %>% 
                                         rownames, start = 6, stop = 1e5),
                        type = substr(re.mv.fo$`~calc.type | study` %>% 
                                        rownames, start = 1, stop = 2),
                        es = re.mv.fo$`~calc.type | study`$intrcpt,
                        se = re.mv.fo$`~calc.type | study`$se) %>% 
    pivot_wider(id_cols = study, names_from = type, 
                values_from = c(es, se)) %>% 
    mutate(es_cs = es_cs + coef(m.mv.fo)[1], 
           es_fo = es_fo + coef(m.mv.fo)[2])
  mv.fo = list(m = m.mv.fo, rg = rg.mv.fo, re = re.mv.fo, bl = bl.mv.fo)
  
  
  m.mv.cs = rma.mv(.TE ~ 0 + calc.type, V, random = ~ calc.type | study, 
                   struct = "UN", data = data.mv.cs)
  rg.mv.cs = matreg(y=1, x=2, R=m.mv.cs$G, cov=TRUE, 
                    means=coef(m.mv.cs), n=m.mv.cs$g.levels.comb.k)
  re.mv.cs = ranef(m.mv.cs)
  bl.mv.cs = data.frame(study = substr(re.mv.cs$`~calc.type | study` %>% 
                                         rownames, start = 6, stop = 1e5),
                        type = substr(re.mv.cs$`~calc.type | study` %>% 
                                        rownames, start = 1, stop = 2),
                        es = re.mv.cs$`~calc.type | study`$intrcpt,
                        se = re.mv.cs$`~calc.type | study`$se) %>% 
    pivot_wider(id_cols = study, names_from = type, 
                values_from = c(es, se)) %>% 
    mutate(es_cs = es_cs + coef(m.mv.cs)[1], 
           es_fo = es_fo + coef(m.mv.cs)[2])
  mv.cs = list(m = m.mv.cs, rg = rg.mv.cs, re = re.mv.cs, bl = bl.mv.cs)
  
  list(mv.pre = mv.pre, mv.fo = mv.fo, mv.cs = mv.cs)
}

# Stop cluster
stopCluster(cl)
save(sim.bv, file="results/sim.bv.pha.rda")


### 4.2.2 Draw Bivariate Plot --------------------------------------------------

pdf("results/plot_slopes.pha.pdf", width=9.5, height = 3.6)

par(mfrow = c(1, 3))
cols = c("dodgerblue1", "dodgerblue2",  
         "dodgerblue3", "dodgerblue4")

plot(1, type="n", xlab=expression(SMD[EP/EP]), ylab=expression(SMD[CS/BL]), xlim=c(-1, 5), ylim=c(-1, 5))
title(expression(SMD[cs]~"standardized by baseline SD"))
polygon(c(-1.5, 5.5, 5.5), c(-1.5, 5.5, -1.5), col = "gray95", border=NA)
for (i in 1:length(sim.bv)) {
  x = sim.bv[[i]]; col = cols[i]
  points(x$mv.pre$bl[,3][[1]], x$mv.pre$bl[,2][[1]], cex=(200+(x$mv.pre$bl[,5][[1]]^2)^-1)/400, pch=19, col="lightgray")
  points(x$mv.pre$bl[,3][[1]], x$mv.pre$bl[,2][[1]], cex=(200+(x$mv.pre$bl[,5][[1]]^2)^-1)/400)
  abline(a=x$mv.pre$rg$tab$beta[1], b=x$mv.pre$rg$tab$beta[2], lwd=3, col=col)
  points(coef(x$mv.pre$m)[2], coef(x$mv.pre$m)[1], pch=19, cex=2, col = col)
}
abline(v=0, lwd=1, lty = 3)
abline(h=0, lwd=1, lty = 3)
abline(a=0, b=1, lwd=1, lty = 3)


plot(1, type="n", xlab=expression(SMD[EP/EP]), ylab=expression(SMD[CS/CS]), xlim=c(-1, 5), ylim=c(-1, 5))
title(expression(SMD[cs]~"standardized by change score SD"))
polygon(c(-1.5, 5.5, 5.5), c(-1.5, 5.5, -1.5), col = "gray95", border=NA)
for (i in 1:length(sim.bv)) {
  x = sim.bv[[i]]; col = cols[i]
  points(x$mv.cs$bl[,3][[1]], x$mv.cs$bl[,2][[1]], cex=(200+(x$mv.cs$bl[,5][[1]]^2)^-1)/400, pch=19, col="lightgray")
  points(x$mv.cs$bl[,3][[1]], x$mv.cs$bl[,2][[1]], cex=(200+(x$mv.cs$bl[,5][[1]]^2)^-1)/400)
  abline(a=x$mv.cs$rg$tab$beta[1], b=x$mv.cs$rg$tab$beta[2], lwd=3, col=col)
  points(coef(x$mv.cs$m)[2], coef(x$mv.cs$m)[1], pch=19, cex=2, col=col)
}
abline(v=0, lwd=1, lty = 3)
abline(h=0, lwd=1, lty = 3)
abline(a=0, b=1, lwd=1, lty = 3)
legend("bottomright", legend=c("r = 0.2", "r = 0.4", "r = 0.6", "r = 0.8"), 
       col=cols[1:length(rhos)], pch=15, cex=0.9, bty="n",
       x.intersp=1, y.intersp=0.5, horiz=TRUE)


plot(1, type="n", xlab=expression(SMD[EP/EP]), ylab=expression(SMD[CS/EP]), xlim=c(-1, 5), ylim=c(-1, 5))
title(expression(SMD[cs]~"standardized by endpoint SD"))
polygon(c(-1.5, 5.5, 5.5), c(-1.5, 5.5, -1.5), col = "gray95", border=NA)
for (i in 1:length(sim.bv)) {
  x = sim.bv[[i]]; col = cols[i]
  points(x$mv.fo$bl[,3][[1]], x$mv.fo$bl[,2][[1]], cex=(200+(x$mv.fo$bl[,5][[1]]^2)^-1)/400, pch=19, col="lightgray")
  points(x$mv.fo$bl[,3][[1]], x$mv.fo$bl[,2][[1]], cex=(200+(x$mv.fo$bl[,5][[1]]^2)^-1)/400)
  abline(a=x$mv.fo$rg$tab$beta[1], b=x$mv.fo$rg$tab$beta[2], lwd=3, col=col)
  points(coef(x$mv.fo$m)[2], coef(x$mv.fo$m)[1], pch=19, cex=2, col=col)
}
abline(v=0, lwd=1, lty = 3)
abline(h=0, lwd=1, lty = 3)
abline(a=0, b=1, lwd=1, lty = 3)

dev.off()



# 5. Save references -----------------------------------------------------------

## 5.1 Psychotherapy versus control --------------------------------------------

dat = read_excel("data/data.xlsx")

# Calculate SMD types
calculateSMDs(dat, rho) %>% 
  distinct(study, .keep_all = TRUE) %>% select(full_ref) %>% 
  write.xlsx("results/references.xlsx", "control")


## 5.2 Psychotherapy versus ADM ------------------------------------------------

# Import data
dat = read_excel("data/data.pha.xlsx") %>% 
  filter(comparison %in% c("psy vs. Pha", "psy vs. pha"))

# Calculate SMD types
calculateSMDs(dat, rho) %>% 
  drop_na(d_cs) %>% 
  distinct(study, .keep_all = TRUE) %>% select(full_ref) %>% 
  write.xlsx("results/references.xlsx", "adm", append=TRUE)
  



