# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#     UTILITY FUNCTIONS                           #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

getIPD = function(study, token){
  
  # Requirements
  require(gh)
  require(base64enc)
  require(jsonlite)
  require(readr)
  require(magrittr)
  
  # Define path
  ghBase = "GET https://api.github.com/repos/prevdep-ipd/data/contents/tidy/"
  suffix = "/data.json"
  path = paste0(ghBase, study, suffix)
  
  # Send API request
  req = try({gh::gh(path, 
                    .token = token,
                    .accept = "application/vnd.github.raw")}, 
            silent = TRUE)
  
  if (identical(class(req)[1], "try-error"))
    stop(study, " could not be retrieved from Github. Either the entry does not exist, ",
         "or the personal access token is not valid.")
  
  # Make pretty
  req %>% 
    rawToChar() %>% jsonlite::fromJSON() %>% 
    apply(2, function(x){
      suppressWarnings({
        readr::parse_guess(x, locale = readr::locale(decimal_mark = "."))
      }) -> col
      ifelse(is.null(attr(col, "problems")), return(col), return(x))
    }, simplify = FALSE) %>%
    as.data.frame() -> dataClean
  
  return(dataClean)
  
}

extractSubgroup = function(x){
  require(emmeans)
  emmprep(x) -> tmp
  emmeans(tmp, "calc.type", type="response", weights="proportional")
  summary(pairs(tmp, adjust = "none")) %>% 
    separate(contrast, sep=" - ", into=c("cs", "fo")) %>% 
    filter(str_detect(cs, "^cs")) %>% 
    filter(str_detect(fo, "^fo")) %>% 
    separate(cs, sep="^([^\\s]*)\\s", into=c(NA, "cs")) %>% 
    separate(fo, sep="^([^\\s]*)\\s", into=c(NA, "fo")) %>% 
    filter(cs == fo)
}

extractRegression = function(x, predictor){
  summary(x) %>% 
    coefficients() %>% 
    mutate(cs = predictor, fo = predictor, SE = se, df = Inf, 
           z.ratio = zval, p.value = pval) %>% 
    select(cs, fo, estimate, SE, df, z.ratio, p.value) %>% 
    {.[4,]}
}

calculateSMDs = function(dat, rho) {
  
  require(tidyr)
  within(dat, {
    
    # Means
    m_t_pos = mean_arm1; m_c_pos = mean_arm2
    m_t_pre = baseline_m_arm1; m_c_pre = baseline_m_arm2
    
    # Sample sizes
    n_t_pos = n_arm1; n_c_pos = n_arm2
    
    # Standard deviations
    s_t_pre = baseline_sd_arm1; s_c_pre = baseline_sd_arm2
    s_t_pos = sd_arm1; s_c_pos = sd_arm2
    
    # Pooled pre-test, post-test, change sd
    s_pre = sqrt(((n_t_pos-1)*s_t_pre^2+(n_c_pos-1)*s_c_pre^2)/(n_t_pos + n_c_pos - 2))
    s_pos = sqrt(((n_t_pos-1)*s_t_pos^2+(n_c_pos-1)*s_c_pos^2)/(n_t_pos + n_c_pos - 2))
    .s_c_cha = sqrt(s_c_pre^2 + s_c_pos^2 - (2 * rho * s_c_pre * s_c_pos))
    .s_t_cha = sqrt(s_t_pre^2 + s_t_pos^2 - (2 * rho * s_t_pre * s_t_pos))
    s_cha = sqrt(((n_t_pos-1)*.s_t_cha^2+(n_c_pos-1)*.s_c_cha^2)/(n_t_pos + n_c_pos - 2))
    
    # Small-sample correction factor
    df = n_t_pos + n_c_pos - 2; cf = 1 - (3/(4*df-1))
    
    # Change SMD (pre-test)
    d_cs = (cf * (((m_t_pos-m_t_pre)-(m_c_pos-m_c_pre))/s_pre))*-1
    d_cs_se = {cf^2 * 2 * (1-rho) * ((n_t_pos+n_c_pos)/(n_t_pos*n_c_pos)) * 
        ((n_t_pos + n_c_pos - 2)/(n_t_pos + n_c_pos - 4)) *
        (1 + ((n_t_pos*n_c_pos)/(n_t_pos+n_c_pos)) * 
           (d_cs^2/(2*(1-rho)))) - d_cs^2} %>% sqrt()
    
    # Change SMD (post-test)
    d_cs_fo = (cf * (((m_t_pos-m_t_pre)-(m_c_pos-m_c_pre))/s_pos))*-1
    d_cs_fo_se = { cf^2 * (((2*(1-rho)*(n_t_pos+n_c_pos))/(n_t_pos*n_c_pos)) +
                             ((d_cs_fo^2)/(2*(n_t_pos+n_c_pos-2)))) } %>% sqrt()
    
    # Change SMD (change SD)
    d_cs_cs = (cf * (((m_t_pos-m_t_pre)-(m_c_pos-m_c_pre))/s_cha))*-1
    d_cs_cs_se = { cf^2 * (((2*(1-rho)*(n_t_pos+n_c_pos))/(n_t_pos*n_c_pos)) +
                             ((d_cs_cs^2)/(2*(n_t_pos+n_c_pos-2)))) } %>% sqrt() 
    
    # Post-test SMD (post-test)
    d_fo = (cf * ((m_t_pos - m_c_pos)/s_pos))*-1
    d_fo_se = { ((n_t_pos+n_c_pos)/(n_t_pos*n_c_pos)) +
        (1-((n_t_pos+n_c_pos-4)/((n_t_pos+n_c_pos-2)*cf^2)))*d_fo^2} %>% sqrt()
    
  }) %>% drop_na(d_cs) -> res
  
  return(res)
  
}









