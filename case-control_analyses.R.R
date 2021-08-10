### "CSF proteomics targeted for CNS processes in bipolar disorder"
### DOI: https://doi.org/10.1038/s41380-021-01236-5
### Correspondance: andreas.goteson@gu.se


# Introduction ------------------------------------------------------------

### The actual datasets are not published due to ethical concerns.

### This script is a generalized version of the actual scripts that 
### contain the real data, arranged in a list of preprocessed files 
### from the two cohorts (SBP-S and SBP-G)

library(tidyverse)
library(broom)

# Case-control analyses ---------------------------------------------------

scale2 = function(x) x = (x - mean(x, na.rm=T)) / sd(x, na.rm=T)

# New model from serum code (OR)
case_controlOR = function(datalist){
  # Logistic regressions
  logi_res = lapply(datalist, function(j){
    j %>% split(.$GENE) %>% 
      map(~ glm(PK~scale2(NPX) + age + sex + albR, family="binomial", data=.) %>% 
            broom::tidy(exponentiate=T, conf.int=T))
  })
  
  # Stack proteins
  logi_stacked = lapply(logi_res, function(j){
    out = data.frame()
    for(i in 1:length(j)){
      tmp = slice(j[[i]], 2)
      out = rbind(out, tmp)
    }
    out$term = names(j)
    return(out)
  })
  
  # And bind cohorts
  logi_rbind = bind_rows(
    logi_stacked$sthlm %>% mutate(site = "SBP-S"),
    logi_stacked$gbg %>% mutate(site="SBP-G")) %>% 
    rename("GENE" = term)
  
  # Get a replicated
  neuro_replicated = logi_rbind %>% 
    group_by(GENE) %>% 
    mutate(repli = all(p.value < 0.05)) %>% 
    filter(repli == T) %>% 
    pull(GENE) %>% unique
  
  out = list(
    res = logi_res,
    bind = logi_rbind,
    replicated = neuro_replicated)
  
  return(out)
}

# BD ~ Control
bdctrlOR = case_controlOR(datalist)

### and sim. for secondary analyses


# Association with clinical outcome ---------------------------------------

clinGlm2 = function(df, vars = c("episodes","BDduration", "CGI_life")){
  mods = list()
  
  for (i in vars){
    formula = as.formula(paste0(i, "~ NPX + age + sex + albR + diag"))
    mods[[i]] = glm(formula, family = "gaussian", 
                    data = df %>% filter(!is.na(i))) %>% tidy %>% slice(2)
  }
  
  out = cbind.data.frame(clinvar=vars,
                         term = sapply(mods, "[[", 1),
                         est = sapply(mods, "[[", 2),
                         p = sapply(mods, "[[", 5)) %>% 
    arrange(p)
  
  return(out)
}

prot4 = c("CLEC1B", "SPOCK1", "TNFRSF21", "DRAXIN")

# Combined cohorts
for (i in prot4) clinGlm2(neu$both %>% filter(GENE == i & PK == "BD")) %>% print


# Drug-protein analyses ---------------------------------------------------

glmDrug = datalist$both %>% filter(case_control=="BD") %>%
  split(.$GENE) %>% 
  map(~ lm(scale2(NPX) ~ age + sex + albR + diag + Li + AD + AP + AC, data = .) %>% 
        broom::tidy(conf.int=T))

glmDrug_out = data.frame()
for (i in 1:length(glmDrug)){
  tmp_df = glmDrug[[i]] %>% slice(7:10) 
  tmp_df$GENE = names(glmDrug)[i]
  glmDrug_out = rbind(glmDrug_out, tmp_df)
}

glmDrug_out$FDR = p.adjust(glmDrug_out$p.value, method="fdr")

