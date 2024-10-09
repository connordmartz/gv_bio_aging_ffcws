###### Proximity to Deadly Neighborhood Gun Violence Accelerates Biological Aging among Adolescents #######|
# created by Connor Martz
# last edit: 10/09/2024

# [1] load libraries ----
library(tidyverse)
library(magrittr)
library(dplyr)
library(plyr) 
library(ggplot2)
library(gridExtra)
library(grid)
library(fmsb)
library(sjmisc)
library(haven)
library(glue)
library(rhdf5)
library(here)
library(janitor)
library(gt)
library(reshape2)
library(filesstrings)
library(naniar)
library(readxl)
library(lmtest)
library(corrplot)
library(png)
library(lars)
library(knitr)
library(kableExtra)
library(lme4)
library(performance)
library(purrr)
library(lubridate)
library(car)
library(ggstatsplot)
library(cowplot)
library(patchwork)
library(MplusAutomation)
library(lars)
library(selectiveInference)
library(broom)
library(sandwich)
library(lmtest)
library(mgcv)
library(easystats)
library(DHARMa)
library(glmmTMB)
library(ggcorrplot)
library(MatchIt)
library(cobalt)
library(twang)
library(survey)
library(gbm)
library(Hmisc)
library(nnet)
library(tableone) 
library(DMwR2) 
library(weights) 
library(rbounds) 
library(randomForest)
library(osDesign)
library(quantreg)
library(forestplot)
library(marginaleffects)
library(DT)
library(GGally)
library(psych)
library(polycor)
library(table1)
# [2] DATA IMPORTS ####
# [2.1] core public data ----
core <- read_dta("~/Library/CloudStorage/Box-Box/FFCWS Data/FF_allwaves_2020v2.dta") 
# [2.2] raw DNAm clock data ---- 
dnam <- read_dta("~/Library/CloudStorage/Box-Box/FFCWS Data/biomarker_final_pub.dta")
# [2.3] tract-level data ----
census00 <- read_dta("/Volumes/FragileFamilies/Data/new_contract_data_sep23/ffgeo7_all_pub.dta") 
census00_reduced <- dplyr::select(census00, "idnum", 
                                  "tp6tract_cen00", # tract ID
                                  "tp6pblck_acs15", # % black
                                  "tp6phisp_acs15", # % hisp
                                  "tp6pwhte_acs15", # % white
                                  "tp6pfrgn_acs15", # % foreign born
                                  "tp6pfbpl_acs15", # % poverty
                                  "tp6mhinc_acs15", # med HH inc
                                  "tp6pfhhr_acs15", # % female-headed HH
                                  "tp6puemp_acs15", # % unemployed 
                                  "tp6ppuba_acs15", # % public assistance
                                  "tp6p25hs_acs15", # % HS degree +
                                  "tp6p25b_acs15",  # % bach degree +
                                  "tp6pfmga_acs15", # % managerial or professional occupation
                                  "tp6pvach_acs15", # % housing units vacant 
                                  "tp6pb10k_acs15", # % familes with income <10k
                                  "tp6p1014_acs15", # % familes with income 10k-14,999k
                                  "tp6p1524_acs15", # % familes with income 15k-24,999
                                  "tp6p2534_acs15", # % familes with income 25k-34,999
                                  "tp6p3549_acs15", # % familes with income 35k-49,999
                                  "tp6p5074_acs15", # % familes with income 50k-74,999
                                  "tp6p7599_acs15", # % familes with income 75k-99,999
                                  "tp6p100k_acs15", # % familes with income 100k-149,999
                                  "tp6p150k_acs15", # % familes with income >=150k
                                  "tm5pfbpl_cen00", 
                                  "tm5tract_cen00", 
                                  "cp6_moved_9_15")
# [2.4] geo identifiers ----
geoid <- read_dta("/Volumes/FragileFamilies/Data/new_contract_data_sep23/contractcity6pub.dta") 
# [2.5] gva data ----         
FF_gva <- read_dta("/Volumes/FragileFamilies/Data/new_contract_data_sep23/ff_gva_15y_pub1_20190426.dta")
# [2.6] county-level crime data (UCR) ----
FF_ucr <- read_dta("/Volumes/FragileFamilies/Data/new_contract_data_sep23/ff_UCR_pub1.dta")
# [3] DATA CLEANING ####
# [3.1] calculate residuals of clocks for full sample ----

# create separate datasets for W5/W6 epic and 450k
Epic.Dataset.W5 <- dnam %>%
  filter(!is.na(k5me_age))  %>%
  dplyr::select(k5me_age, k5me_pedbe, k5me_phenoage, 
                k5me_grim, k5me_poam45, k5me_epi, k5me_fib, k5me_ic, 
                k5me_epithelial, k5me_immune,
                k5me_plasmablast, k5me_cd8pcd28ncd45ran, k5me_cd8_naive, 
                k5me_gr, k5me_nk, k5me_b, k5me_cd4, k5me_cd8, k5me_mo, 
                idnum) 
Epic.Dataset.W6 <- dnam %>%
  filter(!is.na(k6me_age)) %>%
  dplyr::select(k6me_age, k6me_pedbe, k6me_phenoage, 
                k6me_grim, k6me_poam45, k6me_epi, k6me_fib, k6me_ic, 
                k6me_epithelial, k6me_immune,
                k6me_plasmablast, k6me_cd8pcd28ncd45ran,  k6me_cd8_naive, 
                k6me_gr, k6me_nk, k6me_b, k6me_cd4, k6me_cd8, k6me_mo, 
                idnum)
K450.Dataset.W5 <- dnam %>%
  filter(!is.na(k5mk_age)) %>%
  dplyr::select(k5mk_age,  k5mk_pedbe, k5mk_phenoage, 
                k5mk_grim, k5mk_poam45, k5mk_epi, k5mk_fib, k5mk_ic, 
                k5mk_epithelial, k5mk_immune,
                k5mk_plasmablast, k5mk_cd8pcd28ncd45ran,  k5mk_cd8_naive, 
                k5mk_gr, k5mk_nk, k5mk_b, k5mk_cd4, k5mk_cd8, k5mk_mo, 
                idnum) 
K450.Dataset.W6 <- dnam %>%
  filter(!is.na(k6mk_age)) %>%
  dplyr::select(k6mk_age, k6mk_pedbe, k6mk_phenoage, 
                k6mk_grim, k6mk_poam45, k6mk_epi, k6mk_fib, k6mk_ic,
                k6mk_epithelial, k6mk_immune,
                k6mk_plasmablast, k6mk_cd8pcd28ncd45ran,  k6mk_cd8_naive, 
                k6mk_gr, k6mk_nk, k6mk_b, k6mk_cd4, k6mk_cd8, k6mk_mo, 
                idnum) 

# calculate AgeAccelDNAm for overall sample by wave/chip
Epic.Dataset.W5 %<>% mutate(
  grim.raw.W5 = k5me_grim,
  phe.raw.W5 = k5me_phenoage, 
  pace.raw.W5 = k5me_poam45, 
  pbe.raw.W5 = k5me_pedbe, 
  pedbe.resid.W5 = residuals(lm(k5me_pedbe ~ k5me_age, data = .)), 
  phenoage.resid.W5 = residuals(lm(k5me_phenoage ~ k5me_age, data = .)), 
  grim.resid.W5 = residuals(lm(k5me_grim ~ k5me_age, data = .)),
  PoAm45.Z.W5 = scale(k5me_poam45) %>% as.vector, # no residuals for dunedinpace 
  chip = 1,
  age5 = 1*k5me_age,
  epi_t5 = 1*k5me_epi,
  fib_t5 = 1*k5me_fib,
  ic_t5 = 1*k5me_ic, 
  # round values of k5me_epithelial to account for high values of k5me_immune
  epithelial_t5 = case_when( 
    k5me_immune == 1 ~ 0, 
    k5me_immune < 1 ~ 1*k5me_epithelial,
    TRUE ~ NA),
  immune_t5 = 1*k5me_immune,
  plasmablast_t5 = 1*k5me_plasmablast, 
  cd8pcd28ncd45ran_t5 = 1*k5me_cd8pcd28ncd45ran,
  cd8naive_t5 = 1*k5me_cd8_naive,
  gr_t5 = 1*k5me_gr, 
  nk_t5 = 1*k5me_nk, 
  b_t5 = 1*k5me_b,
  cd4_t5 = 1*k5me_cd4,
  cd8_t5 = 1*k5me_cd8, 
  mo_t5 = 1*k5me_mo) 
Epic.Dataset.W6 %<>% mutate(
  grim.raw.W6 = k6me_grim,
  phe.raw.W6 = k6me_phenoage, 
  pace.raw.W6 = k6me_poam45, 
  pbe.raw.W6 = k6me_pedbe, 
  pedbe.resid.W6 = residuals(lm(k6me_pedbe ~ k6me_age, data = .)), 
  phenoage.resid.W6 = residuals(lm(k6me_phenoage ~ k6me_age, data = .)), 
  grim.resid.W6 = residuals(lm(k6me_grim ~ k6me_age, data = .)),
  PoAm45.Z.W6 = scale(k6me_poam45) %>% as.vector, # no residuals for dunedinpace 
  age6 = 1*k6me_age,
  epi_t6 = 1*k6me_epi,
  fib_t6 = 1*k6me_fib,
  ic_t6 = 1*k6me_ic,
  # round values of k5me_epithelial to account for high values of k5me_immune
  epithelial_t6 = case_when(
    k6me_immune == 1 ~ 0, 
    k6me_immune < 1 ~ 1*k6me_epithelial,
    TRUE ~ NA),
  immune_t6 = 1*k6me_immune,
  plasmablast_t6 = 1*k6me_plasmablast, 
  cd8pcd28ncd45ran_t6 = 1*k6me_cd8pcd28ncd45ran,
  cd8naive_t6 = 1*k6me_cd8_naive,
  gr_t6 = 1*k6me_gr, 
  nk_t6 = 1*k6me_nk, 
  b_t6 = 1*k6me_b,
  cd4_t6 = 1*k6me_cd4,
  cd8_t6 = 1*k6me_cd8, 
  mo_t6 = 1*k6me_mo) 
K450.Dataset.W5 %<>% mutate(
  grim.raw.W5 = k5mk_grim,
  phe.raw.W5 = k5mk_phenoage, 
  pace.raw.W5 = k5mk_poam45, 
  pbe.raw.W5 = k5mk_pedbe, 
  pedbe.resid.W5 = residuals(lm(k5mk_pedbe ~ k5mk_age, data = .)), 
  phenoage.resid.W5 = residuals(lm(k5mk_phenoage ~ k5mk_age, data = .)), 
  grim.resid.W5 = residuals(lm(k5mk_grim ~ k5mk_age, data = .)),
  PoAm45.Z.W5 = scale(k5mk_poam45) %>% as.vector, # no residuals for dunedinpace 
  chip = 0,
  age5 = 1*k5mk_age,
  epi_t5 = 1*k5mk_epi,
  fib_t5 = 1*k5mk_fib,
  ic_t5 = 1*k5mk_ic,
  # round values of k5me_epithelial to account for high values of k5me_immune
  epithelial_t5 = case_when(
    k5mk_immune == 1 ~ 0, 
    k5mk_immune < 1 ~ 1*k5mk_epithelial,
    TRUE ~ NA), 
  immune_t5 = 1*k5mk_immune,
  plasmablast_t5 = 1*k5mk_plasmablast, 
  cd8pcd28ncd45ran_t5 = 1*k5mk_cd8pcd28ncd45ran,
  cd8naive_t5 = 1*k5mk_cd8_naive,
  gr_t5 = 1*k5mk_gr, 
  nk_t5 = 1*k5mk_nk, 
  b_t5 = 1*k5mk_b,
  cd4_t5 = 1*k5mk_cd4,
  cd8_t5 = 1*k5mk_cd8, 
  mo_t5 = 1*k5mk_mo) 
K450.Dataset.W6 %<>% mutate(
  grim.raw.W6 = k6mk_grim,
  phe.raw.W6 = k6mk_phenoage, 
  pace.raw.W6 = k6mk_poam45, 
  pbe.raw.W6 = k6mk_pedbe, 
  pedbe.resid.W6 = residuals(lm(k6mk_pedbe ~ k6mk_age, data = .)), 
  phenoage.resid.W6 = residuals(lm(k6mk_phenoage ~ k6mk_age, data = .)), 
  grim.resid.W6 = residuals(lm(k6mk_grim ~ k6mk_age, data = .)),
  PoAm45.Z.W6 = scale(k6mk_poam45) %>% as.vector, # no residuals for dunedinpace 
  age6 = 1*k6mk_age,
  epi_t6 = 1*k6mk_epi,
  fib_t6 = 1*k6mk_fib,
  ic_t6 = 1*k6mk_ic,
  # round values of k5me_epithelial to account for high values of k5me_immune
  epithelial_t6 = case_when(
    k6mk_immune == 1 ~ 0, 
    k6mk_immune < 1 ~ 1*k6mk_epithelial,
    TRUE ~ NA),
  immune_t6 = 1*k6mk_immune,
  plasmablast_t6 = 1*k6mk_plasmablast, 
  cd8pcd28ncd45ran_t6 = 1*k6mk_cd8pcd28ncd45ran,
  cd8naive_t6 = 1*k6mk_cd8_naive,
  gr_t6 = 1*k6mk_gr, 
  nk_t6 = 1*k6mk_nk, 
  b_t6 = 1*k6mk_b,
  cd4_t6 = 1*k6mk_cd4,
  cd8_t6 = 1*k6mk_cd8, 
  mo_t6 = 1*k6mk_mo) 

# join epic and 450k for each wave
joined.comb.w5 <- bind_rows(Epic.Dataset.W5, K450.Dataset.W5)
joined.comb.w6 <- bind_rows(Epic.Dataset.W6, K450.Dataset.W6)

# combine y9 abd y15 into single df 
joined.comb_eaa <- joined.comb.w5 %>%
  left_join(joined.comb.w6, by = "idnum")

# [3.2] merge data frames and rename/code variables   ----
df_merged_clean <- core %>%
  inner_join(joined.comb_eaa, by = "idnum") %>%
  inner_join(geoid, by = "idnum") %>%
  inner_join(census00_reduced, by = "idnum")  %>%
  inner_join(FF_gva, by = "idnum") %>%
  inner_join(FF_ucr, by = "idnum") %>%
  mutate(
    # create date variables for Year 15 DNAm sample receipt 
    date_dnam15 = case_when(
      age6 > 0 ~ as.Date(paste(cm1intyr, cm1intmon, "15", sep = "-")) %m+% months(round(age6 * 12)),
      TRUE ~ NA),
    date_6mon_prior_dnam15 = date_dnam15 %m-% months(6),
    follow_up_month = month(date_dnam15),
    # rename epigenetic clock variables and create standardized variables 
    grmraw_t5 = grim.raw.W5,
    grmraw_t6 = grim.raw.W6,
    pheraw_t5 = phe.raw.W5,
    pheraw_t6 = phe.raw.W6,
    pceraw_t5 = pace.raw.W5,
    pceraw_t6 = pace.raw.W6,
    pberaw_t5 = pbe.raw.W5,
    pberaw_t6 = pbe.raw.W6,
    pbeeaa_t5 = pedbe.resid.W5,
    pbeeaa_t6 = pedbe.resid.W6,
    phe9 = phenoage.resid.W5,
    phe9_std = scale(phe9) %>% as.vector, 
    phe15 = phenoage.resid.W6,
    phe15_std = scale(phe15) %>% as.vector, 
    grm9 = grim.resid.W5,
    grm9_std = scale(grm9) %>% as.vector, 
    grm15 = grim.resid.W6,
    grm15_std = scale(grm15) %>% as.vector, 
    pce9 = PoAm45.Z.W5,
    pce15 = PoAm45.Z.W6,
    
# [3.3] create variables for years of dgv exposure data based on GVA start date (Jan 1, 2014) and date of DNAm sample receipt ----
    gva_start_date = as.Date("2014-01-01"),
    years_dgv_exposure_prior_to_dnam15 = (as.numeric(difftime(date_dnam15, gva_start_date, units = "days")) / 365),
    atleast_1year_exposure = case_when(years_dgv_exposure_prior_to_dnam15 < 1 ~ 0, 
                                      years_dgv_exposure_prior_to_dnam15 >= 1 ~ 1, 
                                      TRUE ~ NA),
# [3.4] create dichotomous variable for any DGV (gvh1600_any) ----
    gvh1600_any = case_when(
      rg6gva_totl_1600m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1600m_home_365d_k6 > 0 ~ 1,
      TRUE ~ NA),
# [3.5] create continuous (gvh1600count) and categorical (gvh1600cat) variables for incident DGV count ---- 
    gvh1600count = ifelse(rg6gva_totl_1600m_home_365d_k6 >= 0, rg6gva_totl_1600m_home_365d_k6, NA),
    gvh1600cat = case_when(
      gvh1600count == 0 ~ 0, # unexposed
      gvh1600count == 1 ~ 1, # low (n=330)
      gvh1600count > 1 & gvh1600count <=5 ~ 2, # moderate (n=396)
      gvh1600count > 5 ~ 3), # high (n=287)
# [3.6] create continuous variables for incident DGV count in 100m intervals ----
    gvh1500_1600count = case_when(
      rg6gva_totl_1600m_home_365d_k6 == 0 & rg6gva_totl_1500m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1600m_home_365d_k6 >= 0 & rg6gva_totl_1500m_home_365d_k6 >= 0 ~ rg6gva_totl_1600m_home_365d_k6 - rg6gva_totl_1500m_home_365d_k6,
      TRUE ~ NA), 
    gvh1400_1500count = case_when(
      rg6gva_totl_1500m_home_365d_k6 == 0 & rg6gva_totl_1400m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1500m_home_365d_k6 >= 0 & rg6gva_totl_1400m_home_365d_k6 >= 0 ~ rg6gva_totl_1500m_home_365d_k6 - rg6gva_totl_1400m_home_365d_k6,
      TRUE ~ NA), 
    gvh1300_1400count = case_when(
      rg6gva_totl_1400m_home_365d_k6 == 0 & rg6gva_totl_1300m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1400m_home_365d_k6 >= 0 & rg6gva_totl_1300m_home_365d_k6 >= 0 ~ rg6gva_totl_1400m_home_365d_k6 - rg6gva_totl_1300m_home_365d_k6,
      TRUE ~ NA), 
    gvh1250_1300count = case_when( 
      rg6gva_totl_1300m_home_365d_k6 == 0 & rg6gva_totl_1250m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1300m_home_365d_k6 >= 0 & rg6gva_totl_1250m_home_365d_k6 >= 0 ~ rg6gva_totl_1300m_home_365d_k6 - rg6gva_totl_1250m_home_365d_k6,
      TRUE ~ NA),
    gvh1200_1250count = case_when(
      rg6gva_totl_1250m_home_365d_k6 == 0 & rg6gva_totl_1200m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1250m_home_365d_k6 >= 0 & rg6gva_totl_1200m_home_365d_k6 >= 0 ~ rg6gva_totl_1250m_home_365d_k6 - rg6gva_totl_1200m_home_365d_k6,
      TRUE ~ NA),
    gvh1200_1300count = case_when(
      rg6gva_totl_1300m_home_365d_k6 == 0 & rg6gva_totl_1200m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1300m_home_365d_k6 >= 0 & rg6gva_totl_1200m_home_365d_k6 >= 0 ~ rg6gva_totl_1300m_home_365d_k6 - rg6gva_totl_1200m_home_365d_k6,
      TRUE ~ NA),
    gvh1100_1200count = case_when(
      rg6gva_totl_1200m_home_365d_k6 == 0 & rg6gva_totl_1100m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1200m_home_365d_k6 >= 0 & rg6gva_totl_1100m_home_365d_k6 >= 0 ~ rg6gva_totl_1200m_home_365d_k6 - rg6gva_totl_1100m_home_365d_k6,
      TRUE ~ NA),
    gvh1000_1100count = case_when(
      rg6gva_totl_1100m_home_365d_k6 == 0 & rg6gva_totl_1000m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1100m_home_365d_k6 >= 0 & rg6gva_totl_1000m_home_365d_k6 >= 0 ~ rg6gva_totl_1100m_home_365d_k6 - rg6gva_totl_1000m_home_365d_k6,
      TRUE ~ NA),
    gvh900_1000count = case_when(
      rg6gva_totl_1000m_home_365d_k6 == 0 & rg6gva_totl_900m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1000m_home_365d_k6 >= 0 & rg6gva_totl_900m_home_365d_k6 >= 0 ~ rg6gva_totl_1000m_home_365d_k6 - rg6gva_totl_900m_home_365d_k6,
      TRUE ~ NA ),
    gvh800_900count = case_when(
      rg6gva_totl_900m_home_365d_k6 == 0 & rg6gva_totl_800m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_900m_home_365d_k6 >= 0 & rg6gva_totl_800m_home_365d_k6 >= 0 ~ rg6gva_totl_900m_home_365d_k6 - rg6gva_totl_800m_home_365d_k6,
      TRUE ~ NA),
    gvh750_800count = case_when(
      rg6gva_totl_800m_home_365d_k6 == 0 & rg6gva_totl_750m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_800m_home_365d_k6 >= 0 & rg6gva_totl_750m_home_365d_k6 >= 0 ~ rg6gva_totl_800m_home_365d_k6 - rg6gva_totl_750m_home_365d_k6,
      TRUE ~ NA),
    gvh700_750count = case_when(
      rg6gva_totl_750m_home_365d_k6 == 0 & rg6gva_totl_700m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_750m_home_365d_k6 >= 0 & rg6gva_totl_700m_home_365d_k6 >= 0 ~ rg6gva_totl_750m_home_365d_k6 - rg6gva_totl_700m_home_365d_k6,
      TRUE ~ NA),
    gvh700_800count = case_when(
      rg6gva_totl_800m_home_365d_k6 == 0 & rg6gva_totl_700m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_800m_home_365d_k6 >= 0 & rg6gva_totl_700m_home_365d_k6 >= 0 ~ rg6gva_totl_800m_home_365d_k6 - rg6gva_totl_700m_home_365d_k6,
      TRUE ~ NA),
    gvh600_700count = case_when(
      rg6gva_totl_700m_home_365d_k6 == 0 & rg6gva_totl_600m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_700m_home_365d_k6 >= 0 & rg6gva_totl_600m_home_365d_k6 >= 0 ~ rg6gva_totl_700m_home_365d_k6 - rg6gva_totl_600m_home_365d_k6,
      TRUE ~ NA),
    gvh500_600count = case_when(
      rg6gva_totl_600m_home_365d_k6 == 0 & rg6gva_totl_500m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_600m_home_365d_k6 >= 0 & rg6gva_totl_500m_home_365d_k6 >= 0 ~ rg6gva_totl_600m_home_365d_k6 - rg6gva_totl_500m_home_365d_k6,
      TRUE ~ NA),
    gvh400_500count = case_when(
      rg6gva_totl_500m_home_365d_k6 == 0 & rg6gva_totl_400m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_500m_home_365d_k6 >= 0 & rg6gva_totl_400m_home_365d_k6 >= 0 ~ rg6gva_totl_500m_home_365d_k6 - rg6gva_totl_400m_home_365d_k6,
      TRUE ~ NA),
    gvh300_400count = case_when(
      rg6gva_totl_400m_home_365d_k6 == 0 & rg6gva_totl_300m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_400m_home_365d_k6 >= 0 & rg6gva_totl_300m_home_365d_k6 >= 0 ~ rg6gva_totl_400m_home_365d_k6 - rg6gva_totl_300m_home_365d_k6,
      TRUE ~ NA),
    gvh250_300count = case_when(
      rg6gva_totl_300m_home_365d_k6 == 0 & rg6gva_totl_250m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_300m_home_365d_k6 >= 0 & rg6gva_totl_250m_home_365d_k6 >= 0 ~ rg6gva_totl_300m_home_365d_k6 - rg6gva_totl_250m_home_365d_k6,
      TRUE ~ NA),
    gvh200_250count = case_when(
      rg6gva_totl_250m_home_365d_k6 == 0 & rg6gva_totl_200m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_250m_home_365d_k6 >= 0 & rg6gva_totl_200m_home_365d_k6 >= 0 ~ rg6gva_totl_250m_home_365d_k6 - rg6gva_totl_200m_home_365d_k6,
      TRUE ~ NA),
    gvh200_300count = case_when(
      rg6gva_totl_300m_home_365d_k6 == 0 & rg6gva_totl_200m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_300m_home_365d_k6 >= 0 & rg6gva_totl_200m_home_365d_k6 >= 0 ~ rg6gva_totl_300m_home_365d_k6 - rg6gva_totl_200m_home_365d_k6,
      TRUE ~ NA),
    gvh100_200count = case_when(
      rg6gva_totl_200m_home_365d_k6 == 0 & rg6gva_totl_100m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_200m_home_365d_k6 >= 0 & rg6gva_totl_100m_home_365d_k6 >= 0 ~ rg6gva_totl_200m_home_365d_k6 - rg6gva_totl_100m_home_365d_k6,
      TRUE ~ NA),
    gvh0_100count = case_when(
      rg6gva_totl_100m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_100m_home_365d_k6 > 0 ~ rg6gva_totl_100m_home_365d_k6,
      TRUE ~ NA), 
# [3.7] create distance-weighted incident DGV count using quadratic decay function (dist_weight_count & dist_weight_count_log) ----
    dist_weight_count = ((gvh0_100count * (17 - 1)^2) + 
                           (gvh100_200count * (17 - 2)^2) + 
                           (gvh200_300count * (17 - 3)^2) + 
                           (gvh300_400count * (17 - 4)^2) + 
                           (gvh400_500count * (17 - 5)^2) + 
                           (gvh500_600count * (17 - 6)^2) + 
                           (gvh600_700count * (17 - 7)^2) + 
                           (gvh700_800count * (17 - 8)^2) + 
                           (gvh800_900count * (17 - 9)^2) + 
                           (gvh900_1000count * (17 - 10)^2) + 
                           (gvh1000_1100count * (17 - 11)^2) + 
                           (gvh1100_1200count * (17 - 12)^2) + 
                           (gvh1200_1300count * (17 - 13)^2) + 
                           (gvh1300_1400count * (17 - 14)^2) + 
                           (gvh1400_1500count * (17 - 15)^2) + 
                           (gvh1500_1600count * (17 - 16)^2)),
    dist_weight_count_log = log(dist_weight_count + 1),
# [3.8] create standardize incident DGV rate by weighting count DGV incidents by years of available GVA data (gvh1600count_years_exposure & gvh1600count_years_exposure_log) ----
    gvh1600count_years_exposure = gvh1600count / years_dgv_exposure_prior_to_dnam15, 
    gvh1600count_years_exposure_log = log(gvh1600count_years_exposure + 1),
# [3.9] individual-level covariates ----
  # time-invariant 
  race = case_when(
      ck6ethrace == 1 ~ 1,         # white=1
      ck6ethrace == 2 ~ 2,         # black=2
      ck6ethrace == 3 ~ 3,         # hispanic=3
      ck6ethrace %in% c(4, 5) ~ 4, # multi-racial + other non-hispanic=4
      # if missing on youth self-identified race/ethnicity, use parent self-identified race/ethnicity
      ck6ethrace %in% c(-1, -2, -3) & cm1ethrace == 2 & (cf1ethrace == 2 | (m1i4 == 2 & m1i4a != 1)) ~ 2, 
      ck6ethrace %in% c(-1, -2, -3) & cm1ethrace == 1 & (cf1ethrace == 1 | (m1i4 == 1 & m1i4a != 1)) ~ 1,
      ck6ethrace %in% c(-1, -2, -3) & cm1ethrace == 3 & (cf1ethrace == 3 | (m1i4 == 5 & m1i4a == 1)) ~ 3,
      ck6ethrace %in% c(-1, -2, -3) & (cm1ethrace %in% c(4, 5) | cf1ethrace %in% c(4, 5)) ~ 4,
      ck6ethrace %in% c(-1, -2, -3) & cm1ethrace == 1 & cf1ethrace %in% c(2, 3) ~ 4,
      ck6ethrace %in% c(-1, -2, -3) & cf1ethrace == 1 & cm1ethrace %in% c(2, 3) ~ 4,
      ck6ethrace %in% c(-1, -2, -3) & cm1ethrace == 2 & cf1ethrace %in% c(1, 3) ~ 4,
      ck6ethrace %in% c(-1, -2, -3) & cf1ethrace == 2 & cm1ethrace %in% c(1, 3) ~ 4,
      ck6ethrace %in% c(-1, -2, -3) & cm1ethrace == 3 & cf1ethrace %in% c(1, 2) ~ 4,
      ck6ethrace %in% c(-1, -2, -3) & cf1ethrace == 3 & cm1ethrace %in% c(1, 2) ~ 4,
      ck6ethrace == -9 & cm1ethrace == 1 & cf1ethrace == 1 ~ 1, 
      ck6ethrace == -9 & cm1ethrace == 2 & cf1ethrace == 2 ~ 2, 
      ck6ethrace == -9 & cm1ethrace == 3 & cf1ethrace == 3 ~ 3, 
      ck6ethrace == -9 & cm1ethrace == 1 & cf1ethrace %in% c(2, 3, 4, 5) ~ 4, 
      ck6ethrace == -9 & cf1ethrace == 1 & cm1ethrace %in% c(2, 3, 4, 5) ~ 4, 
      TRUE ~ NA),
    black = case_when(race==2 ~ 1, 
                      race %in% c(1,3,4) ~ 0, 
                      TRUE ~ NA), 
    white = case_when(race==1 ~ 1, 
                      race %in% c(2,3,4) ~ 0, 
                      TRUE ~ NA), 
    hisp = case_when(race==3 ~ 1, 
                     race %in% c(1,2,4) ~ 0, 
                     TRUE ~ NA), 
    other = case_when(race==4 ~ 1, 
                      race %in% c(1,2,3) ~ 0, 
                      TRUE ~ NA), 
    female = case_when(
      cm1bsex == 1 ~ 0,
      cm1bsex == 2 ~ 1,
      TRUE ~ NA),
    momsmk = ifelse(is.na(m1g4), NA, 
                    ifelse(m1g4 == 4, 0, 1)),
  # time-varying 
    age_t5 = case_when(age5 > 0 ~ age5, TRUE ~ NA),
    age_t6 = case_when(age6 > 0 ~ age6, TRUE ~ NA), # age in months at DNAm sample receipt
    edu_t5 = case_when(cm5edu > 0 ~ cm5edu, TRUE ~ NA),
    nohs_t5 = case_when(
      edu_t5 == 1 ~ 1, 
      edu_t5 %in% c(2,3,4) ~ 0,
      TRUE ~ NA), 
    hs_t5 = case_when(
      edu_t5 == 2 ~ 1, 
      edu_t5 %in% c(1,3,4) ~ 0,
      TRUE ~ NA), 
    somcol_t5 = case_when(
      edu_t5 == 3 ~ 1, 
      edu_t5 %in% c(1,2,4) ~ 0,
      TRUE ~ NA), 
    col_t5 = case_when(
      edu_t5 == 4 ~ 1, 
      edu_t5 %in% c(1,2,3) ~ 0,
      TRUE ~ NA), 
    edu_t6 = case_when(cp6edu > 0 ~ cp6edu, TRUE ~ NA),
    nohs_t6 = case_when(
      edu_t6 == 1 ~ 1, 
      edu_t6 %in% c(2,3,4) ~ 0,
      TRUE ~ NA), 
    hs_t6 = case_when(
      edu_t6 == 2 ~ 1, 
      edu_t6 %in% c(1,3,4) ~ 0,
      TRUE ~ NA), 
    somcol_t6 = case_when(
      edu_t6 == 3 ~ 1, 
      edu_t6 %in% c(1,2,4) ~ 0,
      TRUE ~ NA), 
    col_t6 = case_when(
      edu_t6 == 4 ~ 1, 
      edu_t6 %in% c(1,2,3) ~ 0,
      TRUE ~ NA), 
    pov_t5 = case_when(cm5povco > 0 ~ cm5povco, TRUE ~ NA),
    pov_t6 = case_when(cp6povco > 0 ~ cp6povco, TRUE ~ NA),
    
# [3.10] tract-level covariates ---- 
    trct_id5 = case_when(tm5tract_cen00 >= 0 ~ tm5tract_cen00, TRUE ~ NA),
    trct_pct_pov_t5 = case_when(tm5pfbpl_cen00 >= 0 ~ tm5pfbpl_cen00*100, TRUE ~ NA),
    trct_id6 = case_when(tp6tract_cen00 >= 0 ~ tp6tract_cen00, TRUE ~ NA),
    trct_pct_pov_t6 = case_when(tp6pfbpl_acs15 >= 0 ~ tp6pfbpl_acs15*100, TRUE ~ NA),
    moved_tracts_9_15 = case_when(
      cp6_moved_9_15 == 1 ~ 1, 
      cp6_moved_9_15 == 2 ~ 0,
      TRUE ~ NA),

# [3.11] county-level crime rates at year 9 ----
    vc_rate9 = case_when(rg5ucr_mviort >= 0 ~ rg5ucr_mviort, TRUE ~ NA), 
    vc_rate9_1000 = vc_rate9/100) # convert rate to per 1,000 to reduce variance

# [3.12] calculate race-specific epigenetic age acceleration ----
# white 
Epic.Dataset.W5w <- df_merged_clean %>%
  filter(!is.na(k5me_age) & race == 1)  %>%
  dplyr::select(k5me_age, k5me_pedbe, k5me_phenoage, 
                k5me_grim, k5me_poam45, k5me_epi, k5me_fib, k5me_ic, idnum) 
Epic.Dataset.W6w <- df_merged_clean %>%
  filter(!is.na(k6me_age) & race == 1) %>%
  dplyr::select(k6me_age, k6me_pedbe, k6me_phenoage, 
                k6me_grim, 
                k6me_poam45, k6me_epi, k6me_fib, k6me_ic, idnum)
K450.Dataset.W5w <- df_merged_clean %>%
  filter(!is.na(k5mk_age) & race == 1) %>%
  dplyr::select(k5mk_age,  k5mk_pedbe, k5mk_phenoage, 
                k5mk_grim, 
                k5mk_poam45, k5mk_epi, k5mk_fib, k5mk_ic, idnum) 
K450.Dataset.W6w <- df_merged_clean %>%
  filter(!is.na(k6mk_age) & race == 1) %>%
  dplyr::select(k6mk_age, k6mk_pedbe, k6mk_phenoage, 
                k6mk_grim, 
                k6mk_poam45, k6mk_epi, k6mk_fib, k6mk_ic,idnum) 
Epic.Dataset.W5w %<>% mutate(
  grim.raw.W5w = k5me_grim,
  phe.raw.W5w = k5me_phenoage, 
  pace.raw.W5w = k5me_poam45, 
  pbe.raw.W5w = k5me_pedbe, 
  pedbe.resid.W5w = residuals(lm(k5me_pedbe ~ k5me_age, data = .)), 
  phenoage.resid.W5w = residuals(lm(k5me_phenoage ~ k5me_age, data = .)), 
  grim.resid.W5w = residuals(lm(k5me_grim ~ k5me_age, data = .)),
  phe9w_std = scale(phenoage.resid.W5w) %>% as.vector, 
  grm9w_std = scale(grim.resid.W5w) %>% as.vector, 
  PoAm45.Z.W5w = scale(k5me_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W5w:PoAm45.Z.W5w))
Epic.Dataset.W6w %<>% mutate(
  grim.raw.W6w = k6me_grim,
  phe.raw.W6w = k6me_phenoage, 
  pace.raw.W6w = k6me_poam45, 
  pbe.raw.W6w = k6me_pedbe, 
  pedbe.resid.W6w = residuals(lm(k6me_pedbe ~ k6me_age, data = .)), 
  phenoage.resid.W6w = residuals(lm(k6me_phenoage ~ k6me_age, data = .)), 
  grim.resid.W6w = residuals(lm(k6me_grim ~ k6me_age, data = .)),
  phe15w_std = scale(phenoage.resid.W6w) %>% as.vector, 
  grm15w_std = scale(grim.resid.W6w) %>% as.vector, 
  PoAm45.Z.W6w = scale(k6me_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W6w:PoAm45.Z.W6w))
K450.Dataset.W5w %<>% mutate(
  grim.raw.W5w = k5mk_grim,
  phe.raw.W5w = k5mk_phenoage, 
  pace.raw.W5w = k5mk_poam45, 
  pbe.raw.W5w = k5mk_pedbe, 
  pedbe.resid.W5w = residuals(lm(k5mk_pedbe ~ k5mk_age, data = .)), 
  phenoage.resid.W5w = residuals(lm(k5mk_phenoage ~ k5mk_age, data = .)), 
  grim.resid.W5w = residuals(lm(k5mk_grim ~ k5mk_age, data = .)),
  phe9w_std = scale(phenoage.resid.W5w) %>% as.vector, 
  grm9w_std = scale(grim.resid.W5w) %>% as.vector, 
  PoAm45.Z.W5w = scale(k5mk_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W5w:PoAm45.Z.W5w))
K450.Dataset.W6w %<>% mutate(
  grim.raw.W6w = k6mk_grim,
  phe.raw.W6w = k6mk_phenoage, 
  pace.raw.W6w = k6mk_poam45, 
  pbe.raw.W6w = k6mk_pedbe, 
  pedbe.resid.W6w = residuals(lm(k6mk_pedbe ~ k6mk_age, data = .)), 
  phenoage.resid.W6w = residuals(lm(k6mk_phenoage ~ k6mk_age, data = .)), 
  grim.resid.W6w = residuals(lm(k6mk_grim ~ k6mk_age, data = .)),
  phe15w_std = scale(phenoage.resid.W6w) %>% as.vector, 
  grm15w_std = scale(grim.resid.W6w) %>% as.vector, 
  PoAm45.Z.W6w = scale(k6mk_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W6w:PoAm45.Z.W6w))

joined.comb.W5w <- bind_rows(Epic.Dataset.W5w, K450.Dataset.W5w)
joined.comb.W6w <- bind_rows(Epic.Dataset.W6w, K450.Dataset.W6w)

joined.comb_eaaw <- joined.comb.W5w %>%
  left_join(joined.comb.W6w, by = "idnum")

# black 
Epic.Dataset.W5b <- df_merged_clean %>%
  filter(!is.na(k5me_age) & race == 2)  %>%
  dplyr::select(k5me_age, k5me_pedbe, k5me_phenoage, 
                k5me_grim, k5me_poam45, k5me_epi, k5me_fib, k5me_ic, idnum) 
Epic.Dataset.W6b <- df_merged_clean %>%
  filter(!is.na(k6me_age) & race == 2) %>%
  dplyr::select(k6me_age, k6me_pedbe, k6me_phenoage, 
                k6me_grim, 
                k6me_poam45, k6me_epi, k6me_fib, k6me_ic, idnum)
K450.Dataset.W5b <- df_merged_clean %>%
  filter(!is.na(k5mk_age) & race == 2) %>%
  dplyr::select(k5mk_age,  k5mk_pedbe, k5mk_phenoage, 
                k5mk_grim, 
                k5mk_poam45, k5mk_epi, k5mk_fib, k5mk_ic, idnum) 
K450.Dataset.W6b <- df_merged_clean %>%
  filter(!is.na(k6mk_age) & race == 2) %>%
  dplyr::select(k6mk_age, k6mk_pedbe, k6mk_phenoage, 
                k6mk_grim, 
                k6mk_poam45, k6mk_epi, k6mk_fib, k6mk_ic,idnum) 
Epic.Dataset.W5b %<>% mutate(
  grim.raw.W5b = k5me_grim,
  phe.raw.W5b = k5me_phenoage, 
  pace.raw.W5b = k5me_poam45, 
  pbe.raw.W5b = k5me_pedbe, 
  pedbe.resid.W5b = residuals(lm(k5me_pedbe ~ k5me_age, data = .)), 
  phenoage.resid.W5b = residuals(lm(k5me_phenoage ~ k5me_age, data = .)), 
  grim.resid.W5b = residuals(lm(k5me_grim ~ k5me_age, data = .)),
  phe9b_std = scale(phenoage.resid.W5b) %>% as.vector, 
  grm9b_std = scale(grim.resid.W5b) %>% as.vector, 
  PoAm45.Z.W5b = scale(k5me_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W5b:PoAm45.Z.W5b))
Epic.Dataset.W6b %<>% mutate(
  grim.raw.W6b = k6me_grim,
  phe.raw.W6b = k6me_phenoage, 
  pace.raw.W6b = k6me_poam45, 
  pbe.raw.W6b = k6me_pedbe, 
  pedbe.resid.W6b = residuals(lm(k6me_pedbe ~ k6me_age, data = .)), 
  phenoage.resid.W6b = residuals(lm(k6me_phenoage ~ k6me_age, data = .)), 
  grim.resid.W6b = residuals(lm(k6me_grim ~ k6me_age, data = .)),
  phe15b_std = scale(phenoage.resid.W6b) %>% as.vector, 
  grm15b_std = scale(grim.resid.W6b) %>% as.vector, 
  PoAm45.Z.W6b = scale(k6me_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W6b:PoAm45.Z.W6b))
K450.Dataset.W5b %<>% mutate(
  grim.raw.W5b = k5mk_grim,
  phe.raw.W5b = k5mk_phenoage, 
  pace.raw.W5b = k5mk_poam45, 
  pbe.raw.W5b = k5mk_pedbe, 
  pedbe.resid.W5b = residuals(lm(k5mk_pedbe ~ k5mk_age, data = .)), 
  phenoage.resid.W5b = residuals(lm(k5mk_phenoage ~ k5mk_age, data = .)), 
  grim.resid.W5b = residuals(lm(k5mk_grim ~ k5mk_age, data = .)),
  phe9b_std = scale(phenoage.resid.W5b) %>% as.vector, 
  grm9b_std = scale(grim.resid.W5b) %>% as.vector, 
  PoAm45.Z.W5b = scale(k5mk_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W5b:PoAm45.Z.W5b))
K450.Dataset.W6b %<>% mutate(
  grim.raw.W6b = k6mk_grim,
  phe.raw.W6b = k6mk_phenoage, 
  pace.raw.W6b = k6mk_poam45, 
  pbe.raw.W6b = k6mk_pedbe, 
  pedbe.resid.W6b = residuals(lm(k6mk_pedbe ~ k6mk_age, data = .)), 
  phenoage.resid.W6b = residuals(lm(k6mk_phenoage ~ k6mk_age, data = .)), 
  grim.resid.W6b = residuals(lm(k6mk_grim ~ k6mk_age, data = .)),
  phe15b_std = scale(phenoage.resid.W6b) %>% as.vector, 
  grm15b_std = scale(grim.resid.W6b) %>% as.vector, 
  PoAm45.Z.W6b = scale(k6mk_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W6b:PoAm45.Z.W6b))

# join epic and 450k for each wave
joined.comb.W5b <- bind_rows(Epic.Dataset.W5b, K450.Dataset.W5b)
joined.comb.W6b <- bind_rows(Epic.Dataset.W6b, K450.Dataset.W6b)

joined.comb_eaab <- joined.comb.W5b %>%
  left_join(joined.comb.W6b, by = "idnum")

# hispanic 
Epic.Dataset.W5h <- df_merged_clean %>%
  filter(!is.na(k5me_age) & race == 3)  %>%
  dplyr::select(k5me_age, k5me_pedbe, k5me_phenoage, 
                k5me_grim, k5me_poam45, k5me_epi, k5me_fib, k5me_ic, idnum) 
Epic.Dataset.W6h <- df_merged_clean %>%
  filter(!is.na(k6me_age) & race == 3) %>%
  dplyr::select(k6me_age, k6me_pedbe, k6me_phenoage, 
                k6me_grim, 
                k6me_poam45, k6me_epi, k6me_fib, k6me_ic, idnum)
K450.Dataset.W5h <- df_merged_clean %>%
  filter(!is.na(k5mk_age) & race == 3) %>%
  dplyr::select(k5mk_age,  k5mk_pedbe, k5mk_phenoage, 
                k5mk_grim, 
                k5mk_poam45, k5mk_epi, k5mk_fib, k5mk_ic, idnum) 
K450.Dataset.W6h <- df_merged_clean %>%
  filter(!is.na(k6mk_age) & race == 3) %>%
  dplyr::select(k6mk_age, k6mk_pedbe, k6mk_phenoage, 
                k6mk_grim, 
                k6mk_poam45, k6mk_epi, k6mk_fib, k6mk_ic,idnum) 
Epic.Dataset.W5h %<>% mutate(
  grim.raw.W5h = k5me_grim,
  phe.raw.W5h = k5me_phenoage, 
  pace.raw.W5h = k5me_poam45, 
  pbe.raw.W5h = k5me_pedbe, 
  pedbe.resid.W5h = residuals(lm(k5me_pedbe ~ k5me_age, data = .)), 
  phenoage.resid.W5h = residuals(lm(k5me_phenoage ~ k5me_age, data = .)), 
  grim.resid.W5h = residuals(lm(k5me_grim ~ k5me_age, data = .)),
  phe9h_std = scale(phenoage.resid.W5h) %>% as.vector, 
  grm9h_std = scale(grim.resid.W5h) %>% as.vector, 
  PoAm45.Z.W5h = scale(k5me_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W5h:PoAm45.Z.W5h))
Epic.Dataset.W6h %<>% mutate(
  grim.raw.W6h = k6me_grim,
  phe.raw.W6h = k6me_phenoage, 
  pace.raw.W6h = k6me_poam45, 
  pbe.raw.W6h = k6me_pedbe, 
  pedbe.resid.W6h = residuals(lm(k6me_pedbe ~ k6me_age, data = .)), 
  phenoage.resid.W6h = residuals(lm(k6me_phenoage ~ k6me_age, data = .)), 
  grim.resid.W6h = residuals(lm(k6me_grim ~ k6me_age, data = .)),
  phe15h_std = scale(phenoage.resid.W6h) %>% as.vector, 
  grm15h_std = scale(grim.resid.W6h) %>% as.vector, 
  PoAm45.Z.W6h = scale(k6me_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W6h:PoAm45.Z.W6h))
K450.Dataset.W5h %<>% mutate(
  grim.raw.W5h = k5mk_grim,
  phe.raw.W5h = k5mk_phenoage, 
  pace.raw.W5h = k5mk_poam45, 
  pbe.raw.W5h = k5mk_pedbe, 
  pedbe.resid.W5h = residuals(lm(k5mk_pedbe ~ k5mk_age, data = .)), 
  phenoage.resid.W5h = residuals(lm(k5mk_phenoage ~ k5mk_age, data = .)), 
  grim.resid.W5h = residuals(lm(k5mk_grim ~ k5mk_age, data = .)),
  phe9h_std = scale(phenoage.resid.W5h) %>% as.vector, 
  grm9h_std = scale(grim.resid.W5h) %>% as.vector, 
  PoAm45.Z.W5h = scale(k5mk_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W5h:PoAm45.Z.W5h))
K450.Dataset.W6h %<>% mutate(
  grim.raw.W6h = k6mk_grim,
  phe.raw.W6h = k6mk_phenoage, 
  pace.raw.W6h = k6mk_poam45, 
  pbe.raw.W6h = k6mk_pedbe, 
  pedbe.resid.W6h = residuals(lm(k6mk_pedbe ~ k6mk_age, data = .)), 
  phenoage.resid.W6h = residuals(lm(k6mk_phenoage ~ k6mk_age, data = .)), 
  grim.resid.W6h = residuals(lm(k6mk_grim ~ k6mk_age, data = .)),
  phe15h_std = scale(phenoage.resid.W6h) %>% as.vector, 
  grm15h_std = scale(grim.resid.W6h) %>% as.vector, 
  PoAm45.Z.W6h = scale(k6mk_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W6h:PoAm45.Z.W6h))

joined.comb.W5h <- bind_rows(Epic.Dataset.W5h, K450.Dataset.W5h)
joined.comb.W6h <- bind_rows(Epic.Dataset.W6h, K450.Dataset.W6h)

joined.comb_eaah <- joined.comb.W5h %>%
  left_join(joined.comb.W6h, by = "idnum")

# multiracial/some other race  
Epic.Dataset.W5o <- df_merged_clean %>%
  filter(!is.na(k5me_age) & race == 4)  %>%
  dplyr::select(k5me_age, k5me_pedbe, k5me_phenoage, 
                k5me_grim, k5me_poam45, k5me_epi, k5me_fib, k5me_ic, idnum) 
Epic.Dataset.W6o <- df_merged_clean %>%
  filter(!is.na(k6me_age) & race == 4) %>%
  dplyr::select(k6me_age, k6me_pedbe, k6me_phenoage, 
                k6me_grim, 
                k6me_poam45, k6me_epi, k6me_fib, k6me_ic, idnum)
K450.Dataset.W5o <- df_merged_clean %>%
  filter(!is.na(k5mk_age) & race == 4) %>%
  dplyr::select(k5mk_age,  k5mk_pedbe, k5mk_phenoage, 
                k5mk_grim, 
                k5mk_poam45, k5mk_epi, k5mk_fib, k5mk_ic, idnum) 
K450.Dataset.W6o <- df_merged_clean %>%
  filter(!is.na(k6mk_age) & race == 4) %>%
  dplyr::select(k6mk_age, k6mk_pedbe, k6mk_phenoage, 
                k6mk_grim, 
                k6mk_poam45, k6mk_epi, k6mk_fib, k6mk_ic,idnum) 
Epic.Dataset.W5o %<>% mutate(
  grim.raw.W5o = k5me_grim,
  phe.raw.W5o = k5me_phenoage, 
  pace.raw.W5o = k5me_poam45, 
  pbe.raw.W5o = k5me_pedbe, 
  pedbe.resid.W5o = residuals(lm(k5me_pedbe ~ k5me_age, data = .)), 
  phenoage.resid.W5o = residuals(lm(k5me_phenoage ~ k5me_age, data = .)), 
  grim.resid.W5o = residuals(lm(k5me_grim ~ k5me_age, data = .)),
  phe9o_std = scale(phenoage.resid.W5o) %>% as.vector, 
  grm9o_std = scale(grim.resid.W5o) %>% as.vector, 
  PoAm45.Z.W5o = scale(k5me_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W5o:PoAm45.Z.W5o))
Epic.Dataset.W6o %<>% mutate(
  grim.raw.W6o = k6me_grim,
  phe.raw.W6o = k6me_phenoage, 
  pace.raw.W6o = k6me_poam45, 
  pbe.raw.W6o = k6me_pedbe, 
  pedbe.resid.W6o = residuals(lm(k6me_pedbe ~ k6me_age, data = .)), 
  phenoage.resid.W6o = residuals(lm(k6me_phenoage ~ k6me_age, data = .)), 
  grim.resid.W6o = residuals(lm(k6me_grim ~ k6me_age, data = .)),
  phe15o_std = scale(phenoage.resid.W6o) %>% as.vector, 
  grm15o_std = scale(grim.resid.W6o) %>% as.vector, 
  PoAm45.Z.W6o = scale(k6me_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W6o:PoAm45.Z.W6o))
K450.Dataset.W5o %<>% mutate(
  grim.raw.W5o = k5mk_grim,
  phe.raw.W5o = k5mk_phenoage, 
  pace.raw.W5o = k5mk_poam45, 
  pbe.raw.W5o = k5mk_pedbe, 
  pedbe.resid.W5o = residuals(lm(k5mk_pedbe ~ k5mk_age, data = .)), 
  phenoage.resid.W5o = residuals(lm(k5mk_phenoage ~ k5mk_age, data = .)), 
  grim.resid.W5o = residuals(lm(k5mk_grim ~ k5mk_age, data = .)),
  phe9o_std = scale(phenoage.resid.W5o) %>% as.vector, 
  grm9o_std = scale(grim.resid.W5o) %>% as.vector, 
  PoAm45.Z.W5o = scale(k5mk_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W5o:PoAm45.Z.W5o))
K450.Dataset.W6o %<>% mutate(
  grim.raw.W6o = k6mk_grim,
  phe.raw.W6o = k6mk_phenoage, 
  pace.raw.W6o = k6mk_poam45, 
  pbe.raw.W6o = k6mk_pedbe, 
  pedbe.resid.W6o = residuals(lm(k6mk_pedbe ~ k6mk_age, data = .)), 
  phenoage.resid.W6o = residuals(lm(k6mk_phenoage ~ k6mk_age, data = .)), 
  grim.resid.W6o = residuals(lm(k6mk_grim ~ k6mk_age, data = .)),
  phe15o_std = scale(phenoage.resid.W6o) %>% as.vector, 
  grm15o_std = scale(grim.resid.W6o) %>% as.vector, 
  PoAm45.Z.W6o = scale(k6mk_poam45) %>% as.vector) %>%
  subset(select = c(idnum, phenoage.resid.W6o:PoAm45.Z.W6o))

joined.comb.W5o <- bind_rows(Epic.Dataset.W5o, K450.Dataset.W5o)
joined.comb.W6o <- bind_rows(Epic.Dataset.W6o, K450.Dataset.W6o)

joined.comb_eaao <- joined.comb.W5o %>%
  left_join(joined.comb.W6o, by = "idnum")

# merge race-specific clocks df with df_merged_clean
df_merged_clean <- df_merged_clean %>%
  left_join(joined.comb_eaaw, by = "idnum") %>%
  left_join(joined.comb_eaab, by = "idnum") %>%
  left_join(joined.comb_eaah, by = "idnum") %>%
  left_join(joined.comb_eaao, by = "idnum") %>%
  dplyr::rename(
    phe9w = phenoage.resid.W5w,
    phe15w = phenoage.resid.W6w,
    phe9b= phenoage.resid.W5b,
    phe15b = phenoage.resid.W6b,
    phe9h = phenoage.resid.W5h,
    phe15h = phenoage.resid.W6h,
    phe9o = phenoage.resid.W5o,
    phe15o = phenoage.resid.W6o,
    grm9w = grim.resid.W5w,
    grm15w = grim.resid.W6w,
    grm9b = grim.resid.W5b,
    grm15b = grim.resid.W6b,
    grm9h = grim.resid.W5h,
    grm15h = grim.resid.W6h,
    grm9o = grim.resid.W5o,
    grm15o = grim.resid.W6o,
    pce9w = PoAm45.Z.W5w,
    pce15w = PoAm45.Z.W6w,
    pce9b = PoAm45.Z.W5b,
    pce15b = PoAm45.Z.W6b,
    pce9h = PoAm45.Z.W5h,
    pce15h = PoAm45.Z.W6h,
    pce9o = PoAm45.Z.W5o,
    pce15o = PoAm45.Z.W6o)

# [4] final data prep ----
# classify categorical variables as factors
relevant_factor_vars <- c('edu_t6', 'edu_t5', 'gvh1600cat')
df_merged_clean[relevant_factor_vars] <- lapply(df_merged_clean[relevant_factor_vars], function(x) as.factor(x))
df_merged_clean$race <- factor(df_merged_clean$race, 
                               levels = c(1, 2, 3, 4),
                               labels = c("White", "Black", "Hispanic", "Other"))
df_merged_clean$race <- relevel(df_merged_clean$race, ref = "Black")


# create final df retaining only relevant variables 
df_final <- df_merged_clean %>% dplyr::select(idnum, 
                                        grm9, phe9, pce9, grm15, phe15, pce15, grm9_std, phe9_std, grm15_std, phe15_std,
                                        grm9b, phe9b, pce9b, grm15b, phe15b, pce15b, grm9b_std, phe9b_std, grm15b_std, phe15b_std,
                                        grm9w, phe9w, pce9w, grm15w, phe15w, pce15w, grm9w_std, phe9w_std, grm15w_std, phe15w_std,
                                        grm9h, phe9h, pce9h, grm15h, phe15h, pce15h, grm9h_std, phe9h_std, grm15h_std, phe15h_std,
                                        grm9o, phe9o, pce9o, grm15o, phe15o, pce15o, grm9o_std, phe9o_std, grm15o_std, phe15o_std,
                                        chip, plasmablast_t5, cd8pcd28ncd45ran_t5, cd8naive_t5, plasmablast_t6, cd8pcd28ncd45ran_t6, cd8naive_t6, 
                                        age_t5, age_t6, female, race, black, white, hisp, other, momsmk, 
                                        edu_t5, nohs_t5, hs_t5, somcol_t5, col_t5, edu_t6, nohs_t6, hs_t6, somcol_t6, col_t6, pov_t5, pov_t6, 
                                        trct_id5, trct_pct_pov_t5, trct_id6, trct_pct_pov_t6, moved_tracts_9_15, 
                                        vc_rate9, vc_rate9_1000, 
                                        gvh1600_any, gvh1600count, dist_weight_count, dist_weight_count_log, gvh1600cat, gvh1600count_years_exposure, gvh1600count_years_exposure_log,
                                        years_dgv_exposure_prior_to_dnam15, atleast_1year_exposure)



# [5] INSPECT CLEANED DATA -----
# [5.1] intraclass correlation coefficients ----
icc_grm <- lmer(grm15 ~ (1 | trct_id6), data = df_final)
icc_phe <- lmer(phe15 ~ (1 | trct_id6), data = df_final) 
icc_pce <- lmer(pce15 ~ (1 | trct_id6), data = df_final)
icc(icc_grm, by_group = TRUE) 
icc(icc_phe, by_group = TRUE) 
icc(icc_pce, by_group = TRUE) 

# [5.2] missing data ----
model_vars <- df_final %>% dplyr::select(grm15, pce15, phe15, chip, plasmablast_t6, cd8naive_t6, cd8pcd28ncd45ran_t6, 
                                         age_t6, female, race, momsmk, pov_t6, edu_t6, trct_pct_pov_t6, vc_rate9_1000, 
                                         gvh1600_any, gvh1600count, gvh1600cat, gvh1600count_years_exposure)
vis_miss(model_vars)
gg_miss_upset(model_vars, nsets = n_var_miss(model_vars))

# [6] DATA ANALYSIS ----   
# MAIN MODELS ----
# [Table 1] descriptive statistics -----
# generate table 
summary_table <- df_final %>%
  group_by(race) %>%
  summarise(
    Age_Mean = mean(age_t6, na.rm = TRUE),
    Age_SD = sd(age_t6, na.rm = TRUE),
    Female_Percent = mean(female) * 100,
    Income_Poverty_Mean = mean(pov_t6, na.rm = TRUE),
    Income_Poverty_SD = sd(pov_t6, na.rm = TRUE),
    Tract_Poverty_Mean = mean(trct_pct_pov_t6, na.rm = TRUE),
    Tract_Poverty_SD = sd(trct_pct_pov_t6, na.rm = TRUE),
    Any_DGV_Exposure_Percent = mean(gvh1600_any) * 100,
    Count_DGV_Exposures_Mean = mean(gvh1600count, na.rm = TRUE),
    Count_DGV_Exposures_SD = sd(gvh1600count, na.rm = TRUE),
    GrimAge_Mean = mean(grm15, na.rm = TRUE),
    GrimAge_SD = sd(grm15, na.rm = TRUE),
    PhenoAge_Mean = mean(phe15, na.rm = TRUE),
    PhenoAge_SD = sd(phe15, na.rm = TRUE),
    DunedinPACE_Mean = mean(pce15, na.rm = TRUE),
    DunedinPACE_SD = sd(pce15, na.rm = TRUE))
print(summary_table)
table1(~ age_t6 + factor(female) + pov_t6 + factor(edu_t6) + trct_pct_pov_t6 + factor(gvh1600_any) + gvh1600count + 
         grm15 + phe15 + pce15 | race, data=df_final)

# calculate percentages for maternal educational attainment by race
education_table <- df %>%
  group_by(race, maternal_education) %>%
  summarize(
    Percent = n() / sum(n()) * 100
  ) %>%
  pivot_wider(names_from = maternal_education, values_from = Percent)

# combine continuous/dichotomous variables with categorical educational attainment 
final_table <- summary_table %>%
  left_join(education_table, by = "race")

# print final table
print(final_table)

## test for sig differences in any DGV by race/ethnicity 
any_dgv_race_data <- table(df_final$race, df_final$gvh1600_any)
chi_sq_test <- chisq.test(any_dgv_race_data)
print(chi_sq_test)

## test for sig differences in incident DGV count by race/ethnicity 
anova_test <- aov(gvh1600count ~ race, data = df_final)
summary(anova_test)


    
# [In Text] unstandardized linear regression estimates for DGV count ----
out_1_list <- list(
  # m1: base covariates 
  m1_grm = lm(grm15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final),
  m1_phe = lm(phe15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final),
  m1_pce = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final),  
  # m2: + individual SES
  m2_grm = lm(grm15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final),
  m2_phe = lm(phe15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final),
  m2_pce = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final),
  # m3: + neighborhood SES
  m3_grm = lm(grm15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final),
  m3_phe = lm(phe15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final),
  m3_pce = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final),
  # m4: + race 
  m4_grm = lm(grm15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final),
  m4_phe = lm(phe15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final),
  m4_pce = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final), 
  # m5: Y9 variables 
  m5_grm = lm(grm15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final),
  m5_phe = lm(phe15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final),
  m5_pce = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final))

summarize_model <- function(model) {
  tidy_model <- tidy(model, conf.int = TRUE)
  tidy_model %>% mutate(summary = paste0(sprintf("%.3f", estimate), 
                                         " (", sprintf("%.3f", conf.low), 
                                         ", ", sprintf("%.3f", conf.high), ")")) %>% 
    dplyr::select(term, summary)}

# grimage 
models_grm <- list(
  m1_grm = out_1_list$m1_grm,
  m2_grm = out_1_list$m2_grm,
  m3_grm = out_1_list$m3_grm,
  m4_grm = out_1_list$m4_grm,
  m5_grm = out_1_list$m5_grm)
grm_summaries <- lapply(models_grm, summarize_model)
grm_summary_df <- Reduce(function(x, y) full_join(x, y, by = "term"), grm_summaries)
colnames(grm_summary_df) <- c("Variable", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
datatable(grm_summary_df, options = list(pageLength = 2, autoWidth = TRUE))

#dunedinpace
models_pce <- list(
  m1_pce = out_1_list$m1_pce,
  m2_pce = out_1_list$m2_pce,
  m3_pce = out_1_list$m3_pce,
  m4_pce = out_1_list$m4_pce,
  m5_pce = out_1_list$m5_pce)
pce_summaries <- lapply(models_pce, summarize_model)
pce_summary_df <- Reduce(function(x, y) full_join(x, y, by = "term"), pce_summaries)
colnames(pce_summary_df) <- c("Variable", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
datatable(pce_summary_df, options = list(pageLength = 2, autoWidth = TRUE))

# phenoage 
models_phe <- list(
  m1_phe = out_1_list$m1_phe,
  m2_phe = out_1_list$m2_phe,
  m3_phe = out_1_list$m3_phe,
  m4_phe = out_1_list$m4_phe,
  m5_phe = out_1_list$m5_phe)
phe_summaries <- lapply(models_phe, summarize_model)
phe_summary_df <- Reduce(function(x, y) full_join(x, y, by = "term"), phe_summaries)
colnames(phe_summary_df) <- c("Variable", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
datatable(phe_summary_df, options = list(pageLength = 2, autoWidth = TRUE))

# num obs 
num_observations <- sapply(out_1_list, function(model) {
  nobs(model)})
num_observations_df <- data.frame(
  model = names(num_observations),
  num_observations = num_observations)
print(num_observations_df)

# [In-Text] unstandardized ATT estimates for DGV count ----
# select model variables and omit obs missing on any 
df_psm <- df_final %>% dplyr::select(idnum, grm15, phe15, pce15, grm9, phe9, pce9, chip, 
                                     plasmablast_t6, cd8pcd28ncd45ran_t6, cd8naive_t6, plasmablast_t5, cd8pcd28ncd45ran_t5, cd8naive_t5, 
                                     age_t5, age_t6,  female, race, momsmk, 
                                     pov_t5, edu_t5, trct_pct_pov_t5, pov_t6, edu_t6, trct_pct_pov_t6, vc_rate9, vc_rate9_1000, 
                                     gvh1600_any, gvh1600count, dist_weight_count_log, gvh1600cat)
df_psm_nomiss <- na.omit(df_psm)
df_psm_nomiss <- as.data.frame(df_psm_nomiss)

# propensity score matching 
m.out <- matchit(gvh1600_any ~ race + pov_t6 + edu_t6 + trct_pct_pov_t6 + pov_t5 + edu_t5 + trct_pct_pov_t5 + vc_rate9_1000, 
                 data = df_psm_nomiss, 
                 method = "full",
                 distance = "glm", 
                 link = "logit", 
                 estimand = "att",
                 caliper = c(pov_t6 = .25, trct_pct_pov_t6 = .25, pov_t5 = .25, trct_pct_pov_t5 = .25, vc_rate9_1000 = 1), # larger vc rate caliper used to maximize sample size while maintaining covariate and PS balance 
                 std.caliper = c(TRUE, TRUE, TRUE, TRUE, TRUE),
                 exact = ~ race)

# check matching 
summary(m.out, addlvariables = ~ gvh1600count + grm15 + phe15 + pce15 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, 
        data = df_psm_nomiss, un = FALSE) 
bal.tab(m.out)

# check plots 
custom_labels <- c(
  "distance" = "Propensity Score",
  "race_1" = "White",
  "race_2" = "Black",
  "race_3" = "Hispanic",
  "race_4" = "Other/Multi-Racial",
  "pov_t6" = "Y15: Income-Poverty Ratio",
  "edu_t6_1" = "Y15: No High School",
  "edu_t6_2" = "Y15: High School",
  "edu_t6_3" = "Y15: Some College",
  "edu_t6_4" = "Y15: College",
  "trct_pct_pov_t6" = "Y15: Tract % Poverty", 
  "pov_t5" = "Y9: Income-Poverty Ratio",
  "edu_t5_1" = "Y9: No High School",
  "edu_t5_2" = "Y9: High School",
  "edu_t5_3" = "Y9: Some College",
  "edu_t5_4" = "Y9: College",
  "trct_pct_pov_t5" = "Y9: Tract % Poverty",
  "vc_rate9_1000," = "Y9: County-Level Violent Crime Rate")

love.plot(m.out, threshold = 0.10, var.names = custom_labels) # all covaraites balanced (all std. mean dif <.10)
plot(m.out, type = "jitter", interactive = FALSE)
plot(m.out, type='hist')
plot(m.out, type = "density", interactive = FALSE,
     which.xs = ~race + pov_t6 + trct_pct_pov_t6 + edu_t6)

# extract matched data and weights 
matched_data <- match.data(m.out, subclass = "subclass")

# outcome weighted regression models for gvh1600count among matched sample 
# grimage 
lm1 <- lm(grm15 ~ gvh1600count + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights)
coeftest(lm1, vcov. = vcovCL, cluster = ~subclass) # cluster robust estimates
coefci(lm1, vcov. = vcovCL, cluster = ~subclass, level = 0.95) # cluster robust 95% confidence intervals 
performance::check_model(lm1) # regression diagnostics 

# phenoage
lm2 <- lm(phe15 ~ gvh1600count + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights)
coeftest(lm2, vcov. = vcovCL, cluster = ~subclass) # cluster robust estimates
coefci(lm2, vcov. = vcovCL, cluster = ~subclass, level = 0.95) # cluster robust 95% confidence intervals 
performance::check_model(lm2) # regression diagnostics 

# dunedinpace 
lm3 <- lm(pce15 ~ gvh1600count + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights)
coeftest(lm3, vcov. = vcovCL, cluster = ~subclass) # cluster robust estimates
coefci(lm3, vcov. = vcovCL, cluster = ~subclass, level = 0.95) # cluster robust 95% confidence intervals 
performance::check_model(lm3) # regression diagnostics 

# [Figure 1] standardized linear regression estimates for DGV count ----
df_figure_1 <- df_final %>% dplyr::select(idnum, grm15, phe15, pce15, grm15_std, phe15_std, grm9, phe9, pce9, chip,
                                              plasmablast_t6, cd8pcd28ncd45ran_t6, cd8naive_t6, plasmablast_t5, cd8pcd28ncd45ran_t5, cd8naive_t5,
                                              age_t6, age_t5,  female, race, momsmk, 
                                              pov_t5, edu_t5, trct_pct_pov_t5, pov_t6, edu_t6, trct_pct_pov_t6, vc_rate9_1000,  
                                              gvh1600_any, gvh1600count, dist_weight_count, gvh1600cat, gvh1600count_years_exposure)

out_fig1_list <- list(
  # m1: base covariates 
  m1_grm = lm(grm15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  m1_pce = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),  
  m1_phe = lm(phe15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  # m2: + individual SES
  m2_grm = lm(grm15_std ~ gvh1600count + pov_t6 + edu_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  m2_pce = lm(pce15 ~ gvh1600count + pov_t6 + edu_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  m2_phe = lm(phe15_std ~ gvh1600count + pov_t6 + edu_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  # m3: + neighborhood SES
  m3_grm = lm(grm15_std ~ gvh1600count + pov_t6 + edu_t6 + trct_pct_pov_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  m3_pce = lm(pce15 ~ gvh1600count + pov_t6 + edu_t6 + trct_pct_pov_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  m3_phe = lm(phe15_std ~ gvh1600count + pov_t6 + edu_t6 + trct_pct_pov_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  # m4: + race 
  m4_grm = lm(grm15_std ~ gvh1600count + race + pov_t6 + edu_t6 + trct_pct_pov_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  m4_pce = lm(pce15 ~ gvh1600count + race + pov_t6 + edu_t6 + trct_pct_pov_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1), 
  m4_phe = lm(phe15_std ~ gvh1600count + race + pov_t6 + edu_t6 + trct_pct_pov_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  # m5: Y9 variables 
  m5_grm = lm(grm15_std ~ gvh1600count + vc_rate9_1000 + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + pov_t5 + edu_t5 + trct_pct_pov_t5 + race + pov_t6 + edu_t6 + trct_pct_pov_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  m5_pce = lm(pce15 ~ gvh1600count + vc_rate9_1000 + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + pov_t5 + edu_t5 + trct_pct_pov_t5 + race + pov_t6 + edu_t6 + trct_pct_pov_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1),
  m5_phe = lm(phe15_std ~ gvh1600count + vc_rate9_1000 + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + pov_t5 + edu_t5 + trct_pct_pov_t5 + race + pov_t6 + edu_t6 + trct_pct_pov_t6 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_figure_1))

# generate table with final estiamtes and confidence intervals for all models and all clocks 
results_fig1 <- lapply(names(out_fig1_list), function(model_name) {
  model <- out_fig1_list[[model_name]]
  tidy_model <- tidy(model)
  confint_model <- confint(model)
  tidy_model <- tidy_model %>%
    filter(term == "gvh1600count") %>%
    mutate(
      conf.low = confint_model[term, 1],
      conf.high = confint_model[term, 2])
  return(tidy_model)})
results_fig1_df <- do.call(rbind, results_fig1)
names(results_fig1) <- names(out_fig1_list)
df_results_fig1 <- bind_rows(results_fig1, .id = "model")
outcome_labels <- c("GrimAge", "DunedinPACE", "PhenoAge")
num_models <- length(results_fig1)
outcome_vector <- rep(outcome_labels, length.out = nrow(df_results_fig1))
df_results_fig1 <- df_results_fig1 %>% mutate(Outcome = outcome_vector)
df_results_fig1 <- df_results_fig1 %>%
  dplyr::rename(
    Beta = estimate,
    `CI_Lower` = conf.low,
    `CI_Upper` = conf.high,
    `P value` = p.value) %>%
  dplyr::mutate(
    Beta = round(Beta, 3),
    `CI_Lower` = round(`CI_Lower`, 3),
    `CI_Upper` = round(`CI_Upper`, 3),
    `P value` = round(`P value`, 6),
    Significance = case_when(
      `P value` < 0.05 ~ paste0(`P value`, "*"),
      `P value` < 0.1 ~ paste0(`P value`, "+"),
      TRUE ~ as.character(`P value`))) %>%
  dplyr::select(Outcome, model, term, Beta, `CI_Lower`, `CI_Upper`, Significance)
model_groups <- c("m1_", "m2_", "m3_", "m4_", "m5_")
df_results_fig1 <- df_results_fig1 %>%
  mutate(group = case_when(
    grepl("^m1_", model) ~ "Model 1",
    grepl("^m2_", model) ~ "Model 2",
    grepl("^m3_", model) ~ "Model 3",
    grepl("^m4_", model) ~ "Model 4",
    grepl("^m5_", model) ~ "Model 5",
    TRUE ~ model))
border_indices <- which(df_results_fig1$group != lag(df_results_fig1$group, default = first(df_results_fig1$group)))
results_fig1_kable <- df_results_fig1 %>% 
  arrange(group) %>%
  knitr::kable(
    caption = "Model Results for gvh1600count", 
    col.names = c("Outcome","model", "term", "Beta", "CI_Lower", "CI_Upper", "Significance", "group"),
    escape = FALSE) %>%
  kable_styling(full_width = FALSE, position = "center", bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  row_spec(border_indices, extra_css = "border-top: 2px solid black;")
num_observations <- sapply(out_fig1_list, function(model) {
  nobs(model)})
num_observations_df <- data.frame(
  model = names(num_observations),
  num_observations = num_observations)

# print estimates and n obs
print(num_observations_df)
print(results_fig1_kable)

## plot standardized estimates 
df_results_fig1$group <- factor(df_results_fig1$group, levels = c("Model 5", "Model 4", "Model 3", "Model 2", "Model 1"))
df_results_fig1$Outcome <- factor(df_results_fig1$Outcome, levels = c("GrimAge", "DunedinPACE", "PhenoAge"))
outcome_labs <- c("PhenoAge\nResiduals", "DunedinPACE", "GrimAge\nResiduals")
custom_labels <- c("Model 5 (n=1,781)", 
                   "Model 4 (n=1,868)", 
                   "Model 3 (n=1,868)", 
                   "Model 2 (n=1,869)", 
                   "Model 1 (n=1,888)")
lm_plot <- ggplot(df_results_fig1, aes(x = Beta, y = Outcome, shape = group)) +
  geom_point(size = 3, color = "black", position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, color = "black", position = position_dodge(width = 0.5)) + 
  geom_vline(xintercept = 0, linetype = "solid", color = "steelblue") +
  scale_y_discrete(limits = rev(levels(df_results_fig1$Outcome)), labels = outcome_labs) +
  scale_shape_manual(values = c(4, 18, 15, 17, 16), labels = custom_labels) +  
  labs(x = "Standardized Beta Coefficient and 95% CI's", y = NULL, shape = "Model") +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    plot.title.position = "plot",
    axis.text.y = element_text(size = 10, vjust = 0.5), 
    axis.ticks.length = unit(0.2, "cm"),  
    legend.position = "right") +
  guides(shape = guide_legend(reverse = TRUE)) + 
  ggtitle("") 
print(lm_plot)
# save plot 
ggsave(filename = "/Users/cdm/Library/CloudStorage/Box-Box/MANUSCRIPTS/GVA/plots/lm_plot.png", 
       plot = lm_plot, 
       width = 9, height = 4.5, units = "in", dpi = 300)

# [Figure 2] standardized PSM regression estimated for DGV count ----
# restrict to model variables and omit obs missing on any 
df_figure_2 <- df_final %>% dplyr::select(idnum, grm15, phe15, pce15, grm15_std, phe15_std, grm9, phe9, pce9, chip,
                                          plasmablast_t6, cd8pcd28ncd45ran_t6, cd8naive_t6, plasmablast_t5, cd8pcd28ncd45ran_t5, cd8naive_t5,
                                          age_t6, age_t5,  female, race, momsmk, 
                                          pov_t5, edu_t5, trct_pct_pov_t5, pov_t6, edu_t6, trct_pct_pov_t6, vc_rate9_1000,  
                                          gvh1600_any, gvh1600count, gvh1600cat, dist_weight_count, gvh1600count_years_exposure)
df_figure_2 <- na.omit(df_figure_2)
df_figure_2 <- as.data.frame(df_figure_2)

# propensity score matching 
m.out <- matchit(gvh1600_any ~ race + pov_t6 + edu_t6 + trct_pct_pov_t6 + pov_t5 + edu_t5 + trct_pct_pov_t5 + vc_rate9_1000, 
                 data = df_figure_2, 
                 method = "full",
                 distance = "glm", 
                 link = "logit", 
                 estimand = "att",
                 caliper = c(pov_t6 = .25, trct_pct_pov_t6 = .25, pov_t5 = .25, trct_pct_pov_t5 = .25, vc_rate9_1000 = 1), # larger vc rate caliper used to maximize sample size while maintaining covariate and PS balance 
                 std.caliper = c(TRUE, TRUE, TRUE, TRUE, TRUE),
                 exact = ~ race)

# check covariate balance and matching 
summary(m.out, addlvariables = ~ age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, 
        data = df_figure_2, un = FALSE) 
bal.tab(m.out)

# extract matched data and weights 
matched_data <- match.data(m.out, subclass = "subclass")

# outcome weighted regression models among matched sample (confirm correct estimates before plotting)
# grimage 
lm1 <- lm(grm15_std ~ gvh1600count + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights)
coef_lm1 <- coeftest(lm1, vcov. = vcovCL, cluster = ~subclass)["gvh1600count", ] # cluster robust estimates 
ci_lm1 <- coefci(lm1, vcov. = vcovCL, cluster = ~subclass, level = 0.95)["gvh1600count", ] # cluster robust 95% confidence intervals 

# phenoage
lm2 <- lm(phe15_std ~ gvh1600count + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights)
coef_lm2 <- coeftest(lm2, vcov. = vcovCL, cluster = ~subclass)["gvh1600count", ] # cluster robust estimates 
ci_lm2 <- coefci(lm2, vcov. = vcovCL, cluster = ~subclass, level = 0.95)["gvh1600count", ]  # cluster robust 95% confidence intervals 

# dunedinpace 
lm3 <- lm(pce15 ~ gvh1600count + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights)
coef_lm3 <- coeftest(lm3, vcov. = vcovCL, cluster = ~subclass)["gvh1600count", ] # cluster robust estimates 
ci_lm3 <- coefci(lm3, vcov. = vcovCL, cluster = ~subclass, level = 0.95)["gvh1600count", ]  # cluster robust 95% confidence intervals 

# plot 
plot_data <- data.frame(
  Outcome = rep(c("GrimAge\nResiduals", "PhenoAge\nResiduals", "DunedinPACE"), each = 1),
  Estimate = c(coef_lm1[1], coef_lm2[1], coef_lm3[1]),
  CI_Lower = c(ci_lm1[1], ci_lm2[1], ci_lm3[1]),
  CI_Upper = c(ci_lm1[2], ci_lm2[2], ci_lm3[2]))
plot_data$Outcome <- factor(plot_data$Outcome, levels = c("GrimAge\nResiduals", "DunedinPACE", "PhenoAge\nResiduals"))
psm_plot <- ggplot(plot_data, aes(x = Estimate, y = Outcome, color = Outcome)) +
  geom_point(size = 3, color = "black") +  # Set points to black
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, color = "black") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "steelblue") +
  scale_y_discrete(limits = rev(levels(factor(plot_data$Outcome)))) +
  labs(x = "Average Treatment Effect on the Treated (ATT) with Cluster-Robust 95% CI's", y = NULL) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    plot.title.position = "plot",
    axis.text.y = element_text(size = 10, vjust = 0.5),  # Adjust the size and vertical alignment
    axis.ticks.length = unit(0.2, "cm"),  # Adjust tick length to reduce spacing
    legend.position = "none") +
  ggtitle("") 
print(psm_plot)
# save plot 
ggsave(filename = "/Users/cdm/Library/CloudStorage/Box-Box/MANUSCRIPTS/GVA/plots/psm_plot.png", 
       plot = psm_plot, 
       width = 9, height = 4.5, units = "in", dpi = 300)

# SUPPLMENTAL MATERIALS ----
# [Tables S1-S3] unstandardized linear regression estimates for DGV count ----
summarize_model <- function(model) {
  tidy_model <- tidy(model, conf.int = TRUE)
  tidy_model %>% mutate(summary = paste0(sprintf("%.3f", estimate), 
                                         " (", sprintf("%.3f", conf.low), 
                                         ", ", sprintf("%.3f", conf.high), ")")) %>% 
    dplyr::select(term, summary)}

# grimage 
models_grm <- list(
  m1_grm = out_s1_list$m1_grm,
  m2_grm = out_s1_list$m2_grm,
  m3_grm = out_s1_list$m3_grm,
  m4_grm = out_s1_list$m4_grm,
  m5_grm = out_s1_list$m5_grm)
grm_summaries <- lapply(models_grm, summarize_model)
grm_summary_df <- Reduce(function(x, y) full_join(x, y, by = "term"), grm_summaries)
colnames(grm_summary_df) <- c("Variable", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
datatable(grm_summary_df, options = list(pageLength = 28, autoWidth = TRUE))

#dunedinpace
models_pce <- list(
  m1_pce = out_s1_list$m1_pce,
  m2_pce = out_s1_list$m2_pce,
  m3_pce = out_s1_list$m3_pce,
  m4_pce = out_s1_list$m4_pce,
  m5_pce = out_s1_list$m5_pce)
pce_summaries <- lapply(models_pce, summarize_model)
pce_summary_df <- Reduce(function(x, y) full_join(x, y, by = "term"), pce_summaries)
colnames(pce_summary_df) <- c("Variable", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
datatable(pce_summary_df, options = list(pageLength = 28, autoWidth = TRUE))

# phenoage
models_phe <- list(
  m1_phe = out_s1_list$m1_phe,
  m2_phe = out_s1_list$m2_phe,
  m3_phe = out_s1_list$m3_phe,
  m4_phe = out_s1_list$m4_phe,
  m5_phe = out_s1_list$m5_phe)
phe_summaries <- lapply(models_phe, summarize_model)
phe_summary_df <- Reduce(function(x, y) full_join(x, y, by = "term"), phe_summaries)
colnames(phe_summary_df) <- c("Variable", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
datatable(phe_summary_df, options = list(pageLength = 28, autoWidth = TRUE))

# num obs 
num_observations <- sapply(out_s1_list, function(model) {
  nobs(model)})
num_observations_df <- data.frame(
  model = names(num_observations),
  num_observations = num_observations)
print(num_observations_df)

# [Table S4] race-specific standardized linear regression estimates for DGV count  ----

# first prep race-specific datasets and standardize variables within race groups 
df_final_black <- df_final %>% filter(race=="Black") # black 
df_final_hisp <- df_final %>% filter(race=="Hispanic") # hispanic
df_final_other <- df_final %>% filter(race=="Other") # multi-racial/some other race 
df_final_white <- df_final %>% filter(race=="White") # white 

# race-specific models 1-4 
out_s4_list <- list(
  # m1 
  m1_grm_black = lm(grm15b_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_black),
  m1_grm_hisp = lm(grm15h_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_hisp),
  m1_grm_other = lm(grm15o_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_other),
  m1_grm_white = lm(grm15w_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_white),
  m1_pce_black = lm(pce15b ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_black),
  m1_pce_hisp = lm(pce15h ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_hisp),
  m1_pce_other = lm(pce15o ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_other),
  m1_pce_white = lm(pce15w ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_white),
  m1_phe_black = lm(phe15b_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_black),
  m1_phe_hisp = lm(phe15h_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_hisp),
  m1_phe_other = lm(phe15o_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_other),
  m1_phe_white = lm(phe15w_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_white),
  # m2
  m2_grm_black = lm(grm15b_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_black),
  m2_grm_hisp = lm(grm15h_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_hisp),
  m2_grm_other = lm(grm15o_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_other),
  m2_grm_white = lm(grm15w_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_white),
  m2_pce_black = lm(pce15b ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_black),
  m2_pce_hisp = lm(pce15h ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_hisp),
  m2_pce_other = lm(pce15o ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_other),
  m2_pce_white = lm(pce15w ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_white),
  m2_phe_black = lm(phe15b_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_black),
  m2_phe_hisp = lm(phe15h_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_hisp),
  m2_phe_other = lm(phe15o_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_other),
  m2_phe_white = lm(phe15w_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_white),
  # m3
  m3_grm_black = lm(grm15b_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_black),
  m3_grm_hisp = lm(grm15h_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_hisp),
  m3_grm_other = lm(grm15o_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_other),
  m3_grm_white = lm(grm15w_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_white),
  m3_pce_black = lm(pce15b ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_black),
  m3_pce_hisp = lm(pce15h ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_hisp),
  m3_pce_other = lm(pce15o ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_other),
  m3_pce_white = lm(pce15w ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_white),
  m3_phe_black = lm(phe15b_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_black),
  m3_phe_hisp = lm(phe15h_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_hisp),
  m3_phe_other = lm(phe15o_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_other),
  m3_phe_white = lm(phe15w_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_white),
  # m4
  m4_grm_black = lm(grm15b_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + grm9b + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_black),
  m4_grm_hisp = lm(grm15h_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + grm9h + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_hisp),
  m4_grm_other = lm(grm15o_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + grm9o + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_other),
  m4_grm_white = lm(grm15w_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + grm9w + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_white),
  m4_pce_black = lm(pce15b ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + pce9b + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_black),
  m4_pce_hisp = lm(pce15h ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + pce9h + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_hisp),
  m4_pce_other = lm(pce15o ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + pce9o + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_other),
  m4_pce_white = lm(pce15w ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + pce9w + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_white),
  m4_phe_black = lm(phe15b_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + phe9b + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_black),
  m4_phe_hisp = lm(phe15h_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + phe9h + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_hisp),
  m4_phe_other = lm(phe15o_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + phe9o + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_other),
  m4_phe_white = lm(phe15w_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + phe9w + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_white))

# generate results table for main predictor estimates 
extract_estimates <- function(model_list, predictor = "gvh1600count") {
  results <- lapply(model_list, function(model) {
    if (!is.null(model)) {
      tidy_model <- tidy(model, conf.int = TRUE)
      estimate_row <- tidy_model %>% filter(term == predictor)
      if (nrow(estimate_row) == 1) {
        return(paste0(round(estimate_row$estimate, 3), " (", 
                      round(estimate_row$conf.low, 3), ", ", 
                      round(estimate_row$conf.high, 3), ")"))
      }
    }
    return(NA)
  })
  return(results)
}
grm_black <- extract_estimates(out_s4_list[grep("^m\\d_grm_black", names(out_s4_list))])
grm_hisp <- extract_estimates(out_s4_list[grep("^m\\d_grm_hisp", names(out_s4_list))])
grm_other <- extract_estimates(out_s4_list[grep("^m\\d_grm_other", names(out_s4_list))])
grm_white <- extract_estimates(out_s4_list[grep("^m\\d_grm_white", names(out_s4_list))])

pce_black <- extract_estimates(out_s4_list[grep("^m\\d_pce_black", names(out_s4_list))])
pce_hisp <- extract_estimates(out_s4_list[grep("^m\\d_pce_hisp", names(out_s4_list))])
pce_other <- extract_estimates(out_s4_list[grep("^m\\d_pce_other", names(out_s4_list))])
pce_white <- extract_estimates(out_s4_list[grep("^m\\d_pce_white", names(out_s4_list))])

phe_black <- extract_estimates(out_s4_list[grep("^m\\d_phe_black", names(out_s4_list))])
phe_hisp <- extract_estimates(out_s4_list[grep("^m\\d_phe_hisp", names(out_s4_list))])
phe_other <- extract_estimates(out_s4_list[grep("^m\\d_phe_other", names(out_s4_list))])
phe_white <- extract_estimates(out_s4_list[grep("^m\\d_phe_white", names(out_s4_list))])

df_results <- data.frame(
  Outcome_Race_Ethnicity = c(
    "DGV Count: Black Sample", 
    "DGV Count: Hispanic Sample",
    "DGV Count: Multi/Other Race Sample",
    "DGV Count: White Sample",
    "DGV Count: Black Sample", 
    "DGV Count: Hispanic Sample",
    "DGV Count: Multi/Other Race Sample",
    "DGV Count: White Sample",
    "DGV Count: Black Sample", 
    "DGV Count: Hispanic Sample",
    "DGV Count: Multi/Other Race Sample",
    "DGV Count: White Sample"),
  Outcome = rep(c("GrimAge Residuals", "DunedinPACE", "PhenoAge Residuals"), each = 4),
  Model_1 = c(grm_black[[1]], grm_hisp[[1]], grm_other[[1]], grm_white[[1]], 
              pce_black[[1]], pce_hisp[[1]], pce_other[[1]], pce_white[[1]], 
              phe_black[[1]], phe_hisp[[1]], phe_other[[1]], phe_white[[1]]),
  Model_2 = c(grm_black[[2]], grm_hisp[[2]], grm_other[[2]], grm_white[[2]], 
              pce_black[[2]], pce_hisp[[2]], pce_other[[2]], pce_white[[2]], 
              phe_black[[2]], phe_hisp[[2]], phe_other[[2]], phe_white[[2]]),
  Model_3 = c(grm_black[[3]], grm_hisp[[3]], grm_other[[3]], grm_white[[3]], 
              pce_black[[3]], pce_hisp[[3]], pce_other[[3]], pce_white[[3]], 
              phe_black[[3]], phe_hisp[[3]], phe_other[[3]], phe_white[[3]]),
  Model_4 = c(grm_black[[4]], grm_hisp[[4]], grm_other[[4]], grm_white[[4]], 
              pce_black[[4]], pce_hisp[[4]], pce_other[[4]], pce_white[[4]], 
              phe_black[[4]], phe_hisp[[4]], phe_other[[4]], phe_white[[4]]))
datatable(df_results, options = list(pageLength = 12, autoWidth = TRUE))

# n obs 
num_observations <- sapply(out_s4_list, function(model) {
  nobs(model)})
num_observations_df <- data.frame(
  model = names(num_observations),
  num_observations = num_observations)
print(num_observations_df)

# [Figure S1] correlation matrix of DGV count, epigenetic clocks, race/ethnicity, and SES ---- 
# prep dataframe 
df_cor <- df_final %>%
  mutate(
    POC = factor(case_when(
      black == 1 ~ "2",
      white == 1 ~ "1",
      hisp == 1 ~ "2",
      other == 1 ~ "2"), 
      ordered = FALSE), 
    edu_cat = factor(case_when(
      nohs_t6 == 1 ~ "1", 
      hs_t6 == 1 ~ "2", 
      somcol_t6 == 1 ~ "3", 
      col_t6 == 1 ~ "4"), 
      ordered = TRUE)) %>%
  dplyr::select(gvh1600count, POC, edu_cat, pov_t6, trct_pct_pov_t6, grm15, phe15, pce15)

# ensure continuous variables are numeric 
df_cor <- df_cor %>%
  mutate(across(c(gvh1600count, pov_t6, trct_pct_pov_t6, grm15, phe15, pce15), as.numeric))

# remove missing data and convert to df 
df_cor_nomiss <- na.omit(df_cor) 
df_cor_nomiss <- as.data.frame(df_cor_nomiss)

# calculate correlations accounting for different variable types 
cor_matrix <- hetcor(df_cor_nomiss)$correlations

# apply labels 
labels <- c("gvh1600count" = "Incident DGV Count", 
            "POC" = "Non-White Race/Ethnicity",  
            "edu_cat" = "Educational Attainment", 
            "pov_t6" = "Income-Poverty Ratio", 
            "trct_pct_pov_t6" = "Tract % Poverty",
            "grm15" = "GrimAge Residuals", 
            "phe15" = "PhenoAge Residuals", 
            "pce15" = "DunedinPACE")

colnames(cor_matrix) <- labels
rownames(cor_matrix) <- labels

# plot and save correlation matrix 
png("/Users/cdm/Library/CloudStorage/Box-Box/MANUSCRIPTS/GVA/plots/cor_matrix.png", 
    width = 6, height = 6, units = "in", res = 300)

cor_matrix <- corrplot(cor_matrix, 
                       method = 'circle',
                       type = 'lower',
                       tl.col = 'black', 
                       #addCoef.col = 'black', 
                       col = COL2('RdBu', 10),
                       number.cex = 0.6, 
                       cl.ratio = 0.2, 
                       tl.srt = 45, 
                       diag = FALSE, 
                       tl.cex = 0.7)
dev.off()

# [Figure S2] PSM covariate love plot ----
custom_labels <- c(
  "distance" = "Propensity Score",
  "race_1" = "White",
  "race_2" = "Black",
  "race_3" = "Hispanic",
  "race_4" = "Other/Multi-Racial",
  "pov_t6" = "Y15: Income-Poverty Ratio",
  "edu_t6_1" = "Y15: No High School",
  "edu_t6_2" = "Y15: High School",
  "edu_t6_3" = "Y15: Some College",
  "edu_t6_4" = "Y15: College",
  "trct_pct_pov_t6" = "Y15: Tract % Poverty", 
  "pov_t5" = "Y9: Income-Poverty Ratio",
  "edu_t5_1" = "Y9: No High School",
  "edu_t5_2" = "Y9: High School",
  "edu_t5_3" = "Y9: Some College",
  "edu_t5_4" = "Y9: College",
  "trct_pct_pov_t5" = "Y9: Tract % Poverty",
  "vc_rate9" = "Y9: County-Level Violent Crime Rate")
love.plot(m.out, threshold = 0.10, var.names = custom_labels) # all covaraites balanced (all std. mean dif <.10)

# [Figure S3] distribution of propensity scores ----
plot(m.out, type = "jitter", interactive = FALSE)
# [Table S5] unstandardized logistic regression for PSM ----
ps_model <- m.out$model
summary(ps_model)
odds_ratios <- exp(coef(ps_model))
odds_ratios
coefs <- coef(ps_model)
se <- sqrt(diag(vcov(ps_model)))
ci_lower_log_odds <- coefs - 1.96 * se
ci_upper_log_odds <- coefs + 1.96 * se
odds_ratios <- exp(coefs)
ci_lower_odds <- exp(ci_lower_log_odds)
ci_upper_odds <- exp(ci_upper_log_odds)

# generate table 
summary_df <- data.frame(
  Variable = names(coefs),
  Odds_Ratio = odds_ratios,
  CI_Lower = ci_lower_odds,
  CI_Upper = ci_upper_odds)
summary_df

# n obs 
n_obs_manual <- nrow(na.omit(df_psm_nomiss))
n_obs_manual

# [Table S6] descriptive statistics for matched sample -----
# data prep 
vars_continuous <- c("age_t6", "pov_t6", "trct_pct_pov_t6", "gvh1600count", "grm15", "pce15", "phe15")
vars_categorical <- c("race", "female", "edu_t6", "momsmk")
overall_results_continuous <- list()
overall_results_categorical <- list()
results_continuous <- list()
results_categorical <- list()

# calculate overall statistics for continuous variables
for (var in vars_continuous) {
  variable <- matched_data[[var]]
  weights <- matched_data$weights
  weighted_mean <- sum(weights * variable) / sum(weights)
  weighted_variance <- sum(weights * (variable - weighted_mean) ^ 2) / sum(weights)
  weighted_sd <- sqrt(weighted_variance)
  overall_results_continuous[[var]] <- data.frame(
    var_name = var,
    category = NA,
    group = "Overall",
    weighted_mean = weighted_mean,
    weighted_sd = weighted_sd)}
overall_results_continuous_df <- do.call(rbind, overall_results_continuous)

# calculate overall statistics for categorical variables
for (var in vars_categorical) {
  variable <- matched_data[[var]]
  weights <- matched_data$weights
  levels_var <- levels(as.factor(variable))
  for (level in levels_var) {
    weighted_n <- sum(weights[variable == level])
    weighted_pct <- weighted_n / sum(weights) * 100
    overall_results_categorical[[paste(var, level)]] <- data.frame(
      var_name = var,
      category = level,
      group = "Overall",
      weighted_n = weighted_n,
      weighted_pct = weighted_pct)}}
overall_results_categorical_df <- do.call(rbind, overall_results_categorical)

# calculate statistics for continuous variables by group (control and treated)
for (var in vars_continuous) {
  # extract variable for treatment and control groups
  variable_treat <- matched_data[[var]][matched_data$gvh1600_any == 1]
  variable_control <- matched_data[[var]][matched_data$gvh1600_any == 0]
  
  # extract weights for treatment and control groups
  weights_treat <- matched_data$weights[matched_data$gvh1600_any == 1]
  weights_control <- matched_data$weights[matched_data$gvh1600_any == 0]
  
  # calculate weighted mean and SD for treatment group
  weighted_mean_treat <- sum(weights_treat * variable_treat) / sum(weights_treat)
  weighted_variance_treat <- sum(weights_treat * (variable_treat - weighted_mean_treat) ^ 2) / sum(weights_treat)
  weighted_sd_treat <- sqrt(weighted_variance_treat)
  
  # calculate weighted mean and SD for control group
  weighted_mean_control <- sum(weights_control * variable_control) / sum(weights_control)
  weighted_variance_control <- sum(weights_control * (variable_control - weighted_mean_control) ^ 2) / sum(weights_control)
  weighted_sd_control <- sqrt(weighted_variance_control)
  
  # store results for control and treated groups
  results_continuous[[var]] <- data.frame(
    var_name = var,
    category = NA,
    group = c("Control", "Treated"),
    weighted_mean = c(weighted_mean_control, weighted_mean_treat),
    weighted_sd = c(weighted_sd_control, weighted_sd_treat))}
results_continuous_df <- do.call(rbind, results_continuous)

# calculate statistics for categorical variables by group
for (var in vars_categorical) {
  # extract variable and weights for treatment and control groups
  variable_treat <- matched_data[[var]][matched_data$gvh1600_any == 1]
  variable_control <- matched_data[[var]][matched_data$gvh1600_any == 0]
  weights_treat <- matched_data$weights[matched_data$gvh1600_any == 1]
  weights_control <- matched_data$weights[matched_data$gvh1600_any == 0]
  
  # get levels of the categorical variable
  levels_var <- levels(as.factor(matched_data[[var]]))
  
  # calculate weighted counts and percentages for treatment and control groups
  for (level in levels_var) {
    # treatment 
    weighted_n_treat <- sum(weights_treat[variable_treat == level])
    weighted_pct_treat <- weighted_n_treat / sum(weights_treat) * 100
    
    # control 
    weighted_n_control <- sum(weights_control[variable_control == level])
    weighted_pct_control <- weighted_n_control / sum(weights_control) * 100
    
    # store results for control and treated groups
    results_categorical[[paste(var, level)]] <- data.frame(
      var_name = var,
      category = level,
      group = c("Control", "Treated"),
      weighted_n = c(weighted_n_control, weighted_n_treat),
      weighted_pct = c(weighted_pct_control, weighted_pct_treat))}}
results_categorical_df <- do.call(rbind, results_categorical)
combined_continuous <- rbind(
  overall_results_continuous_df %>% 
    mutate(weighted_stat = paste0(round(weighted_mean, 2), " (", round(weighted_sd, 2), ")")) %>% 
    select(var_name, category, group, weighted_stat),
    results_continuous_df %>% 
    mutate(weighted_stat = paste0(round(weighted_mean, 2), " (", round(weighted_sd, 2), ")")) %>% 
    select(var_name, category, group, weighted_stat))

# combine categorical results (overall, control, treated)
combined_categorical <- rbind(
  overall_results_categorical_df %>%
    mutate(weighted_stat = paste0(round(weighted_n), " (", round(weighted_pct, 2), "%)")) %>%
    select(var_name, category, group, weighted_stat),
  results_categorical_df %>%
    mutate(weighted_stat = paste0(round(weighted_n), " (", round(weighted_pct, 2), "%)")) %>%
    select(var_name, category, group, weighted_stat))
combined_results <- rbind(combined_continuous, combined_categorical)
formatted_table <- combined_results %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = weighted_stat,
    id_cols = c(var_name, category)) %>%
  select(var_name, category, Overall, Control, Treated)
colnames(formatted_table) <- c("Variable", "Category", "Overall", "Control", "Treated")

# generate table 
formatted_table %>%
  kbl(caption = "Weighted Descriptive Characteristics of the Matched Sample") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  add_header_above(c(" " = 2, "Overall (n=499)" = 1, "Control (n=265)" = 1, "Treated (n=265)" = 1))

# [Table S7] standardized ATT estimates for DGV count  ----
# grimage 
lm1 <- lm(grm15_std ~ gvh1600count + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights)
coeftest(lm1, vcov. = vcovCL, cluster = ~subclass) # cluster robust estimates 
coefci(lm1, vcov. = vcovCL, cluster = ~subclass, level = 0.95) # cluster robust 95% confidence intervals  

# dunedinpace 
lm2 <- lm(pce15 ~ gvh1600count + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights)
coeftest(lm2, vcov. = vcovCL, cluster = ~subclass) # cluster robust estimates 
coefci(lm2, vcov. = vcovCL, cluster = ~subclass, level = 0.95) # cluster robust 95% confidence intervals 

# phenoage
lm3 <- lm(phe15_std ~ gvh1600count + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights)
coeftest(lm3, vcov. = vcovCL, cluster = ~subclass) # cluster robust estimates 
coefci(lm3, vcov. = vcovCL, cluster = ~subclass, level = 0.95) # cluster robust 95% confidence intervals 

# generate table 
extract_estimates <- function(model) {
  coefs <- tidy(model, conf.int = TRUE) %>%
    mutate(
      estimate_CI = paste0(round(estimate, 3), " (", round(conf.low, 3), ", ", round(conf.high, 3), ")")
    ) %>%
    dplyr::select(term, estimate_CI)
}
extract_robust_estimates <- function(model, cluster_var) {
  robust_coefs <- coeftest(model, vcov. = vcovCL, cluster = cluster_var)
  robust_ci <- coefci(model, vcov. = vcovCL, cluster = cluster_var, level = 0.95)
  coefs_df <- data.frame(
    term = rownames(robust_coefs),
    estimate = robust_coefs[, 1],
    conf.low = robust_ci[, 1],
    conf.high = robust_ci[, 2])
  coefs_df %>%
    mutate(
      estimate_CI = paste0(round(estimate, 3), " (", round(conf.low, 3), ", ", round(conf.high, 3), ")")) %>%
    dplyr::select(term, estimate_CI)}
lm1_coefs <- extract_robust_estimates(lm1, ~subclass)
lm2_coefs <- extract_robust_estimates(lm2, ~subclass)
lm3_coefs <- extract_robust_estimates(lm3, ~subclass)
temp_results <- full_join(lm1_coefs, lm2_coefs, by = "term") %>%
  full_join(lm3_coefs, by = "term")
results <- temp_results %>%
  dplyr::rename(
    Variable = term,
    `GrimAge Residuals (ATT 95% CI)` = estimate_CI.x,
    `DunedinPACE (ATT 95% CI)` = estimate_CI.y,
    `PhenoAge Residuals (ATT 95% CI)` = estimate_CI)
datatable(results, options = list(pageLength = 18, autoWidth = TRUE))

# n obs 
nobs(lm1)
nobs(lm2)
nobs(lm3)

# [Table S8] summary table of standardized estimates from linear models of various exposure measurements ----
out_s8_list <- list(
  ## any DGV
  # m1: base covariates 
  m1_grm_any = lm(grm15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_phe_any = lm(phe15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_pce_any = lm(pce15 ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),  
  # m2: + individual SES
  m2_grm_any = lm(grm15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_phe_any = lm(phe15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_pce_any = lm(pce15 ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  # m3: + neighborhood SES
  m3_grm_any = lm(grm15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_phe_any = lm(phe15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_pce_any = lm(pce15 ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  # m4: + race 
  m4_grm_any = lm(grm15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_phe_any = lm(phe15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_pce_any = lm(pce15 ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std), 
  # m5: Y9 variables 
  m5_grm_any = lm(grm15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + grm9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_phe_any = lm(phe15_std ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + phe9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_pce_any = lm(pce15 ~ gvh1600_any + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  
  ## exposure count 
  # m1: base covariates 
  m1_grm_count = lm(grm15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_phe_count = lm(phe15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_pce_count = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),  
  # m2: + individual SES
  m2_grm_count = lm(grm15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_phe_count = lm(phe15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_pce_count = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  # m3: + neighborhood SES
  m3_grm_count = lm(grm15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_phe_count = lm(phe15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_pce_count = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  # m4: + race 
  m4_grm_count = lm(grm15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_phe_count = lm(phe15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_pce_count = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std), 
  # m5: Y9 variables 
  m5_grm_count = lm(grm15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + grm9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_phe_count = lm(phe15_std ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + phe9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_pce_count = lm(pce15 ~ gvh1600count + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std), 
  
  ## categorical exposure count 
  # m1: base covariates 
  m1_grm_cat = lm(grm15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_phe_cat = lm(phe15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_pce_cat = lm(pce15 ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),  
  # m2: + individual SES
  m2_grm_cat = lm(grm15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_phe_cat = lm(phe15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_pce_cat = lm(pce15 ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  # m3: + neighborhood SES
  m3_grm_cat = lm(grm15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_phe_cat = lm(phe15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_pce_cat = lm(pce15 ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  # m4: + race 
  m4_grm_cat = lm(grm15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_phe_cat = lm(phe15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_pce_cat = lm(pce15 ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std), 
  # m5: Y9 variables 
  m5_grm_cat = lm(grm15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + grm9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_phe_cat = lm(phe15_std ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + phe9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_pce_cat = lm(pce15 ~ gvh1600cat + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  
  ## intensity dgv exposure 
  # m1: base covariates 
  m1_grm_intensity = lm(grm15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_phe_intensity = lm(phe15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_pce_intensity = lm(pce15 ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),  
  # m2: + individual SES
  m2_grm_intensity = lm(grm15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_phe_intensity = lm(phe15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_pce_intensity = lm(pce15 ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  # m3: + neighborhood SES
  m3_grm_intensity = lm(grm15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_phe_intensity = lm(phe15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_pce_intensity = lm(pce15 ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  # m4: + race 
  m4_grm_intensity = lm(grm15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_phe_intensity = lm(phe15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_pce_intensity = lm(pce15 ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std), 
  # m5: Y9 variables 
  m5_grm_intensity = lm(grm15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + grm9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_phe_intensity = lm(phe15_std ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + phe9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_pce_intensity = lm(pce15 ~ dist_weight_count_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  
  ## std. rate count dgv  
  # m1: base covariates 
  m1_grm_rate = lm(grm15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_phe_rate = lm(phe15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),
  m1_pce_rate = lm(pce15 ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data = df_final_std),  
  # m2: + individual SES
  m2_grm_rate = lm(grm15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_phe_rate = lm(phe15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  m2_pce_rate = lm(pce15 ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6, data = df_final_std),
  # m3: + neighborhood SES
  m3_grm_rate = lm(grm15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_phe_rate = lm(phe15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  m3_pce_rate = lm(pce15 ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6, data = df_final_std),
  # m4: + race 
  m4_grm_rate = lm(grm15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_phe_rate = lm(phe15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std),
  m4_pce_rate = lm(pce15 ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race, data = df_final_std), 
  # m5: Y9 variables 
  m5_grm_rate = lm(grm15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + grm9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_phe_rate = lm(phe15_std ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + phe9_std + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std),
  m5_pce_rate = lm(pce15 ~ gvh1600count_years_exposure_std + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip + edu_t6 + pov_t6 + trct_pct_pov_t6 + race + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + edu_t5 + pov_t5 + trct_pct_pov_t5 + vc_rate9_1000, data = df_final_std))

# generate table 
summarize_model <- function(model) {
  tidy_model <- tidy(model, conf.int = TRUE)
  terms_to_include <- c("gvh1600_any", "gvh1600count", "gvh1600cat1", "gvh1600cat2", "gvh1600cat3",
                        "dist_weight_count_std", "gvh1600count_years_exposure_std")
  tidy_model %>%
    filter(term %in% terms_to_include) %>%
    mutate(
      summary = paste0(sprintf("%.3f", estimate), 
                       " (", sprintf("%.3f", conf.low), 
                       ", ", sprintf("%.3f", conf.high), ")")) %>%
    dplyr::select(term, summary)}
combine_results <- function(models) {
  model_summaries <- lapply(models, summarize_model)
  combined_summary <- Reduce(function(x, y) full_join(x, y, by = "term"), model_summaries)
  return(combined_summary)}
models_grm <- list(
  GrimAge_Residuals = list(
    m1_grm_any = out_s8_list$m1_grm_any, m2_grm_any = out_s8_list$m2_grm_any, m3_grm_any = out_s8_list$m3_grm_any, 
    m4_grm_any = out_s8_list$m4_grm_any, m5_grm_any = out_s8_list$m5_grm_any,
    m1_grm_count = out_s8_list$m1_grm_count, m2_grm_count = out_s8_list$m2_grm_count, m3_grm_count = out_s8_list$m3_grm_count, 
    m4_grm_count = out_s8_list$m4_grm_count, m5_grm_count = out_s8_list$m5_grm_count,
    m1_grm_cat = out_s8_list$m1_grm_cat, m2_grm_cat = out_s8_list$m2_grm_cat, m3_grm_cat = out_s8_list$m3_grm_cat, 
    m4_grm_cat = out_s8_list$m4_grm_cat, m5_grm_cat = out_s8_list$m5_grm_cat,
    m1_grm_intensity = out_s8_list$m1_grm_intensity, m2_grm_intensity = out_s8_list$m2_grm_intensity, 
    m3_grm_intensity = out_s8_list$m3_grm_intensity, m4_grm_intensity = out_s8_list$m4_grm_intensity, m5_grm_intensity = out_s8_list$m5_grm_intensity,
    m1_grm_rate = out_s8_list$m1_grm_rate, m2_grm_rate = out_s8_list$m2_grm_rate, m3_grm_rate = out_s8_list$m3_grm_rate, 
    m4_grm_rate = out_s8_list$m4_grm_rate, m5_grm_rate = out_s8_list$m5_grm_rate),
  DunedinPACE = list(
    m1_pce_any = out_s8_list$m1_pce_any, m2_pce_any = out_s8_list$m2_pce_any, m3_pce_any = out_s8_list$m3_pce_any, 
    m4_pce_any = out_s8_list$m4_pce_any, m5_pce_any = out_s8_list$m5_pce_any,
    m1_pce_count = out_s8_list$m1_pce_count, m2_pce_count = out_s8_list$m2_pce_count, m3_pce_count = out_s8_list$m3_pce_count, 
    m4_pce_count = out_s8_list$m4_pce_count, m5_pce_count = out_s8_list$m5_pce_count,
    m1_pce_cat = out_s8_list$m1_pce_cat, m2_pce_cat = out_s8_list$m2_pce_cat, m3_pce_cat = out_s8_list$m3_pce_cat, 
    m4_pce_cat = out_s8_list$m4_pce_cat, m5_pce_cat = out_s8_list$m5_pce_cat,
    m1_pce_intensity = out_s8_list$m1_pce_intensity, m2_pce_intensity = out_s8_list$m2_pce_intensity, 
    m3_pce_intensity = out_s8_list$m3_pce_intensity, m4_pce_intensity = out_s8_list$m4_pce_intensity, m5_pce_intensity = out_s8_list$m5_pce_intensity,
    m1_pce_rate = out_s8_list$m1_pce_rate, m2_pce_rate = out_s8_list$m2_pce_rate, m3_pce_rate = out_s8_list$m3_pce_rate, 
    m4_pce_rate = out_s8_list$m4_pce_rate, m5_pce_rate = out_s8_list$m5_pce_rate),
  PhenoAge_Residuals = list(
    m1_phe_any = out_s8_list$m1_phe_any, m2_phe_any = out_s8_list$m2_phe_any, m3_phe_any = out_s8_list$m3_phe_any, 
    m4_phe_any = out_s8_list$m4_phe_any, m5_phe_any = out_s8_list$m5_phe_any,
    m1_phe_count = out_s8_list$m1_phe_count, m2_phe_count = out_s8_list$m2_phe_count, m3_phe_count = out_s8_list$m3_phe_count, 
    m4_phe_count = out_s8_list$m4_phe_count, m5_phe_count = out_s8_list$m5_phe_count,
    m1_phe_cat = out_s8_list$m1_phe_cat, m2_phe_cat = out_s8_list$m2_phe_cat, m3_phe_cat = out_s8_list$m3_phe_cat, 
    m4_phe_cat = out_s8_list$m4_phe_cat, m5_phe_cat = out_s8_list$m5_phe_cat,
    m1_phe_intensity = out_s8_list$m1_phe_intensity, m2_phe_intensity = out_s8_list$m2_phe_intensity, 
    m3_phe_intensity = out_s8_list$m3_phe_intensity, m4_phe_intensity = out_s8_list$m4_phe_intensity, m5_phe_intensity = out_s8_list$m5_phe_intensity,
    m1_phe_rate = out_s8_list$m1_phe_rate, m2_phe_rate = out_s8_list$m2_phe_rate, m3_phe_rate = out_s8_list$m3_phe_rate, 
    m4_phe_rate = out_s8_list$m4_phe_rate, m5_phe_rate = out_s8_list$m5_phe_rate))
grm_summary <- combine_results(models_grm$GrimAge_Residuals)
dunedin_summary <- combine_results(models_grm$DunedinPACE)
phenoage_summary <- combine_results(models_grm$PhenoAge_Residuals)
grm_summary <- grm_summary %>% mutate(Outcome = "GrimAge Residuals")
dunedin_summary <- dunedin_summary %>% mutate(Outcome = "DunedinPACE")
phenoage_summary <- phenoage_summary %>% mutate(Outcome = "PhenoAge Residuals")
final_summary_df <- bind_rows(grm_summary, dunedin_summary, phenoage_summary)
final_summary_df <- final_summary_df %>%
  mutate(
    Variable = case_when(
      str_detect(term, "gvh1600_any") ~ "Any DGV",
      str_detect(term, "gvh1600count") ~ "Count DGV (log-transformed)",
      str_detect(term, "gvh1600cat1") ~ "Categorical Count DGV: Low",
      str_detect(term, "gvh1600cat2") ~ "Categorical Count DGV: Moderate",
      str_detect(term, "gvh1600cat3") ~ "Categorical Count DGV: High",
      str_detect(term, "dist_weight_count_std") ~ "Intensity DGV Exposure (log-transformed)",
      str_detect(term, "gvh1600count_years_exposure_std") ~ "Std. Rate DGV Count (log-transformed)")) %>%
  dplyr::select(Outcome, Variable, everything()) 
colnames(final_summary_df) <- c("Outcome", "Exposure Variable", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
final_summary_df <- final_summary_df %>%
  arrange(Outcome, `Exposure Variable`)
datatable(final_summary_df, options = list(pageLength = 28, autoWidth = TRUE))

# n obs 
num_observations <- sapply(out_s8_list, function(model) {
  nobs(model)})
num_observations_df <- data.frame(
  model = names(num_observations),
  num_observations = num_observations)
print(num_observations_df)

# [Table S9] summary table of standardized ATT estimates from PSM models of various exposure measurements ----
out_s9_list <- list(
  # grm
  grm_any = lm(grm15_std ~ gvh1600_any + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights), 
  grm_count = lm(grm15_std ~ gvh1600count + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  grm_cat = lm(grm15_std ~ gvh1600cat + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  grm_intensity = lm(grm15_std ~ dist_weight_count_std + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  grm_rate = lm(grm15_std ~ gvh1600count_years_exposure_std + grm9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  # pace
  pce_any = lm(pce15 ~ gvh1600_any + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights), 
  pce_count = lm(pce15 ~ gvh1600count + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  pce_cat = lm(pce15 ~ gvh1600cat + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  pce_intensity = lm(pce15 ~ dist_weight_count_std + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  pce_rate = lm(pce15 ~ gvh1600count_years_exposure_std + pce9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  # phe
  phe_any = lm(phe15_std ~ gvh1600_any + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights), 
  phe_count = lm(phe15_std ~ gvh1600count + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  phe_cat = lm(phe15_std ~ gvh1600cat + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  phe_intensity = lm(phe15_std ~ dist_weight_count_std + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights),
  phe_rate = lm(phe15_std ~ gvh1600count_years_exposure_std + phe9 + age_t5 + plasmablast_t5 + cd8pcd28ncd45ran_t5 + cd8naive_t5 + age_t6 + female + momsmk + plasmablast_t6 + cd8pcd28ncd45ran_t6 + cd8naive_t6 + chip, data=matched_data, weights = weights))

# generate table 
extract_att <- function(model, term) {
  if (is.null(model)) return(NA)  # Return NA if model is NULL
  
  # Get coefficient and cluster-robust standard error for the term
  att_est <- coeftest(model, vcov. = vcovCL, cluster = ~subclass)[term, ]
  
  # Get confidence interval for the term
  ci <- coefci(model, vcov. = vcovCL, cluster = ~subclass, level = 0.95)[term, ]
  
  # Return formatted ATT estimate and CI
  att_summary <- sprintf("%.3f (%.3f, %.3f)", att_est[1], ci[1], ci[2])
  return(att_summary)
}
terms_list <- list(
  grm = list(
    any = "gvh1600_any",
    count = "gvh1600count",
    cat = c("gvh1600cat1", "gvh1600cat2", "gvh1600cat3"),  
    intensity = "dist_weight_count_std",
    rate = "gvh1600count_years_exposure_std"),
  pce = list(
    any = "gvh1600_any",
    count = "gvh1600count",
    cat = c("gvh1600cat1", "gvh1600cat2", "gvh1600cat3"),
    intensity = "dist_weight_count_std",
    rate = "gvh1600count_years_exposure_std"),
  phe = list(
    any = "gvh1600_any",
    count = "gvh1600count",
    cat = c("gvh1600cat1", "gvh1600cat2", "gvh1600cat3"),
    intensity = "dist_weight_count_std",
    rate = "gvh1600count_years_exposure_std"))
results_grm <- list()
results_pce <- list()
results_phe <- list()
extract_categorical_terms <- function(model, terms) {
  results <- list()
  for (term in terms) {
    results[[term]] <- extract_att(model, term)
  }
  return(results)
}
for (term in names(terms_list$grm)) {
  if (term == "cat") {
    cat_model <- out_s9_list[["grm_cat"]]  
    results_grm <- c(results_grm, extract_categorical_terms(cat_model, terms_list$grm$cat))} 
  else {if (!is.null(out_s9_list[[paste0("grm_", term)]])) {
      results_grm[[term]] <- extract_att(out_s9_list[[paste0("grm_", term)]], terms_list$grm[[term]])} 
    else {results_grm[[term]] <- NA}}}
for (term in names(terms_list$pce)) {
  if (term == "cat") {    cat_model <- out_s9_list[["pce_cat"]]  
    results_pce <- c(results_pce, extract_categorical_terms(cat_model, terms_list$pce$cat))} 
  else {if (!is.null(out_s9_list[[paste0("pce_", term)]])) 
    {results_pce[[term]] <- extract_att(out_s9_list[[paste0("pce_", term)]], terms_list$pce[[term]])
    } else {results_pce[[term]] <- NA}}}
for (term in names(terms_list$phe)) {
  if (term == "cat") {cat_model <- out_s9_list[["phe_cat"]]  
    results_phe <- c(results_phe, extract_categorical_terms(cat_model, terms_list$phe$cat))}
  else {if (!is.null(out_s9_list[[paste0("phe_", term)]])) {
      results_phe[[term]] <- extract_att(out_s9_list[[paste0("phe_", term)]], terms_list$phe[[term]])}
    else {results_phe[[term]] <- NA}}}
results_df <- data.frame(
  Variable = c(
    "Any DGV", "Count DGV", "Categorical Count DGV: Low",
    "Categorical Count DGV: Moderate", "Categorical Count DGV: High",
    "Intensity DGV Exposure", "Rate DGV Count"),
  Y15_GrimAge_Residuals = c(
    results_grm$any, results_grm$count, 
    results_grm$gvh1600cat1, results_grm$gvh1600cat2, 
    results_grm$gvh1600cat3, results_grm$intensity, 
    results_grm$rate),
  Y15_DunedinPACE = c(
    results_pce$any, results_pce$count, 
    results_pce$gvh1600cat1, results_pce$gvh1600cat2, 
    results_pce$gvh1600cat3, results_pce$intensity, 
    results_pce$rate),
  Y15_PhenoAge_Residuals = c(
    results_phe$any, results_phe$count, 
    results_phe$gvh1600cat1, results_phe$gvh1600cat2, 
    results_phe$gvh1600cat3, results_phe$intensity, 
    results_phe$rate))
# print table 
datatable(results_df, options = list(pageLength = 10, autoWidth = TRUE))

# [Sens. Test] compare estimates accounting for movers/stayers from Y9-Y15 ----

# first examine differences in main study variables for movers/stayers 
table(df_final$moved_tracts_9_15)                             # n movers = 968 (54%), n missing on move/stay variable = 193
chisq.test(df_final$gvh1600_any, df_final$moved_tracts_9_15)  # no differences in DGV exposure at Y15 between movers and stayers from Y9 to Y15 (p=.32)
t.test(gvh1600count ~ moved_tracts_9_15, data = df_final)     # no differences in count DGV exposures at Y15 between movers and stayers (p=.37)
t.test(grm15 ~ moved_tracts_9_15, data = df_final)            # movers have higher GrimAge at Y15 (p=.064); residential mobility may influence outcome 
t.test(pce15 ~ moved_tracts_9_15, data = df_final)            # movers have higher DunedinPACE at Y15 (p=.079); residential mobility may influence outcome
t.test(pov_t6 ~ moved_tracts_9_15, data = df_final)           # movers have lower inc-pov ratios (higher poverty) at Y15
t.test(trct_pct_pov_t5 ~ moved_tracts_9_15, data = df_final)  # movers have higher tract poverty at Y15
t.test(vc_rate ~ moved_tracts_9_15, data = df_final)          # movers have higher county-level violent crime rates at Y15
t.test(vc_rate9 ~ moved_tracts_9_15, data = df_final)         # movers also had higher county-level violent crime rates at Y9
test <- df_final %>% mutate(                                  # [prep data to test for differences in change in SES and violent crime]
  t_pov_dif = trct_pct_pov_t6-trct_pct_pov_t5, 
  pov_dif = pov_t6-pov_t5, 
  vc_dif = vc_rate - vc_rate9)
t.test(pov_dif ~ moved_tracts_9_15, data=test)                # no differences in changes in inc-pov from Y9 to Y15 between movers and stayers 
t.test(t_pov_dif ~ moved_tracts_9_15, data=test)              # no differences in changes in tract poverty from Y9 to Y15 between movers and stayers 
t.test(vc_dif ~ moved_tracts_9_15, data=test)                 # movers experienced significantly greater reduction in violent crime rates from Y9 to Y15 than stayers 

# notes: 
# Given these tests, matching on propensity of Y15 DGV exposure based on Y15 and Y9 confounders as well as 
# Y9 violent crime rate (as a proxy for DGV exposure given lack of Y9 DGV data) may be an appropriate method to 
# account for selection bias and isolate incident exposure to Y15 DGV. 

# Restricting sample to those who did not move census tracts from Y9 to Y15 significantly reduces 
# sample size even before matching --> using alternative approach to account for selection bias 

# Matching and outcome regression estimates for the main model are consistent with those adjusting for mover/stayers 
# (moved_tracts_9_15) as covariate in outcome regression models and including matching covariate in PSM model. 
