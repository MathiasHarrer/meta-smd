# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  
#   Endpoint or change scores

# Perils of Standardization: Divergences between Change Score and Endpoint 
# Effect Sizes in Meta-Analyses of Psychotherapies for Depression



library(readxl)
library(tidyverse)
library(metapsyTools)
library(metafor)


dat = read_excel("data/data.xlsx")

rho = 0.35
m_t_pos = dat$mean_arm1
m_c_pos = dat$mean_arm2
m_t_pre = dat$baseline_m_arm1
m_c_pre = dat$baseline_m_arm2
n_t_pre = dat$baseline_n_arm1
n_c_pre = dat$baseline_n_arm2
s_t_pre = dat$baseline_sd_arm1
s_c_pre = dat$baseline_sd_arm2
s_t_pos = dat$sd_arm1
s_c_pos = dat$sd_arm2
s_pre = sqrt(((n_t_pre-1)*s_t_pre^2+(n_c_pre-1)*s_c_pre^2)/(n_t_pre + n_c_pre - 2))
s_pos = sqrt(((n_t_pre-1)*s_t_pos^2+(n_c_pre-1)*s_c_pos^2)/(n_t_pre + n_c_pre - 2))
s_c_cha = sqrt(s_c_pre^2 + s_c_pos^2 - (2 * rho * s_c_pre * s_c_pos))
s_t_cha = sqrt(s_t_pre^2 + s_t_pos^2 - (2 * rho * s_t_pre * s_t_pos))
s_cha = sqrt(((n_t_pre-1)*s_t_cha^2+(n_c_pre-1)*s_c_cha^2)/(n_t_pre + n_c_pre - 2))


df = n_t_pre + n_c_pre - 2
cf = 1 - (3/(4*df-1))
d_cs = (cf * (((m_t_pos-m_t_pre)-(m_c_pos-m_c_pre))/s_pre))*-1
d_cs_se = {cf^2 * 2 * (1-rho) * ((n_t_pre+n_c_pre)/(n_t_pre*n_c_pre)) * ((n_t_pre + n_c_pre - 2)/(n_t_pre + n_c_pre - 4)) *
  (1 + ((n_t_pre*n_c_pre)/(n_t_pre+n_c_pre)) * (d_cs^2/(2*(1-rho)))) - d_cs^2} %>% sqrt()


dat$d_cs = d_cs
dat$d_cs_se = d_cs_se
dat$d_fo = dat$.g
dat$d_fo_se = dat$.g_se


#es.var = ifelse(identical(es.measure[1], "g"), 
#                ".g", ".log_rr"),
#se.var = ifelse(identical(es.measure[1], "g"), 
#                ".g_se", ".log_rr_se"),


dat %>% drop_na(d_cs) -> data

runMetaAnalysis(data, which.run = c("overall", "combined"),
                es.var = "d_fo", se.var = "d_fo_se")

runMetaAnalysis(data, which.run = c("overall", "combined"),
                es.var = "d_cs", se.var = "d_cs_se")


runMetaAnalysis(data, which.run = c("overall", "combined",
                                    "lowest.highest", "outliers",
                                    "rob", "threelevel",
                                    "threelevel.che"),
                es.var = "d_fo", se.var = "d_fo_se") -> m.fo

runMetaAnalysis(data, which.run = c("overall", "combined",
                                    "lowest.highest", "outliers",
                                    "rob", "threelevel",
                                    "threelevel.che"),
                es.var = "d_cs", se.var = "d_cs_se") -> m.cs





rma.mv(.TE ~ 0 + calc.type, V, random = ~ calc.type | study, 
       struct = "UN", data=data.mv)


# Multivariate Meta-Analysis Model (k = 790; method: REML)
# 
# Variance Components:
#   
# outer factor: study     (nlvls = 395)
# inner factor: calc.type (nlvls = 2)
# 
#             estim    sqrt  k.lvl  fixed  level 
# tau^2.1    0.6745  0.8213    395     no     cs 
# tau^2.2    0.4751  0.6893    395     no     fo 
# 
#     rho.cs  rho.fo    cs   fo 
# cs       1             -  395 
# fo  0.9898       1    no    - 
#   
# Test for Residual Heterogeneity:
# QE(df = 788) = 9367.0165, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 434.8614, p-val < .0001
# 
# Model Results:
#   
#              estimate      se     zval    pval   ci.lb   ci.ub      
# calc.typecs    0.8782  0.0424  20.7186  <.0001  0.7951  0.9613  *** 
# calc.typefo    0.7395  0.0366  20.2021  <.0001  0.6677  0.8112  *** 
#   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


m.fo$model.combined$data$calc.type = "fo"
m.cs$model.combined$data$calc.type = "cs"
m.fo$model.combined$data$effsize = m.fo$model.combined$data$.TE
m.cs$model.combined$data$effsize = m.fo$model.combined$data$.TE

rbind(m.fo$model.combined$data,
      m.cs$model.combined$data) %>% 
  as.data.frame() %>% 
  arrange(study) %>% 
  mutate(V = .seTE^2,
         condition_arm2 = recode(condition_arm2, "wlc" = "wl"),
         mean_age = parse_number(mean_age),
         percent_women = parse_number(percent_women) %>% ifelse(.>1, ./100,.)) %>% 
  drop_na(mean_age, percent_women) -> data.mv


frml = .TE ~ calc.type*(scale(year) + scale(rob) + scale(percent_women) + scale(mean_age) +
                          country + condition_arm1 + condition_arm2 + recruitment)

rma.mv(.TE ~ calc.type*scale(effsize), V, random = ~ calc.type | study, 
       struct = "UN", data = data.mv)


rma.mv(.TE ~ 0 + calc.type, V, random = ~ calc.type | study, 
       struct = "UN", data = data.mv) -> tmp
reg <- matreg(y=1, x=2, R=tmp$G, cov=TRUE, means=coef(tmp), n=tmp$g.levels.comb.k)
reg

ranef(tmp) -> raneffs
data.frame(
study = substr(raneffs$`~calc.type | study` %>% rownames, start = 6, stop = 1e5),
type = substr(raneffs$`~calc.type | study` %>% rownames, start = 1, stop = 2),
es = raneffs$`~calc.type | study`$intrcpt,
se = raneffs$`~calc.type | study`$se) %>% 
  pivot_wider(id_cols = study, names_from = type, values_from = c(es, se)) %>% 
  mutate(es_cs = es_cs + coef(tmp)[1], es_fo = es_fo + coef(tmp)[2]) -> blups


require(ellipse)
plot(1, type="n", xlab="Post-Test (g)", ylab="Change Score (g)", xlim=c(-1, 5), ylim=c(-1, 5))
points(blups[,3][[1]], blups[,2][[1]], cex = (200+(blups[,5][[1]]^2)^-1)/400, pch=19, col="lightgray")
points(blups[,3][[1]], blups[,2][[1]], cex = (200+(blups[,5][[1]]^2)^-1)/400)
abline(a=reg$tab$beta[1], b=reg$tab$beta[2], lwd=3, col="blue")
points(coef(tmp)[2], coef(tmp)[1], pch=19, cex=2, col = "blue")
abline(v=0, lwd=1, lty = 3)
abline(h=0, lwd=1, lty = 3)
abline(v=coef(tmp)[2], lwd=1, lty = 3, col = "blue")
abline(h=coef(tmp)[1], lwd=1, lty = 3, col = "blue")
abline(a=0, b=1, lwd=1, lty = 3)
#xy <- ellipse(tmp$G, centre=coef(tmp), level=0.95)
#lines(xy[,2],xy[,1], col="gray50")

rma.mv(.TE ~ calc.type*scale(effsize), V, random = ~ 1| study/calc.type, 
       data = data.mv) -> tmp


rma.mv(frml, V, random = ~ calc.type | study, 
       struct = "UN", data = data.mv) -> M

emmprep(M, rg.limit = 12096) -> M.emm



predict(M, colMeans(mat)[-1])





