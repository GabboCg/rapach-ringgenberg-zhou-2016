#'---
#' author: Gabriel Cabrera 
#' title: Asset Allocation Implementation
#' date: 30-12-2018  (last updated: 30-12-2018)
#'---

if(!require("pacman")) install.packages("pacman")
p_load("tidyverse", "xlsx")

GW_dataset <- read.xlsx("data/Returns_short_interest_data.xlsx", 
                        sheetName = "GW variables") 

GW_varset <- as.data.frame(apply(GW_dataset, 2, as.numeric)) %>% 
             filter(yyyymm >= 197212) %>% 
             select(yyyymm, Index, D12, E12, b.m, ntis, tbl, lty, ltr, BAA, 
               AAA, corpr, infl) %>%
             rename(SP = Index, D12 = D12, E12 = E12, BM = b.m, NTIS = ntis,
               TBL = tbl, LTY = lty, LTR = ltr, BAA = BAA, AAA = AAA,
               CORPR = corpr, INFL_lag = infl) %>%  
             mutate(log_DP = log(D12/SP), log_EP = log(E12/SP), 
               log_DE = log(D12/E12), log_DY = log(D12/lag(SP)),
               TMS = LTY - TBL, DFY = BAA - AAA, DFR = CORPR - LTR,
               INFL_lag = lag(INFL_lag)) %>% 
             select(log_DP, log_DY, log_EP, log_DE, BM, NTIS, TBL, LTY, 
               LTR, TMS, DFY, DFR, INFL_lag) %>% 
             na.omit()

# Stock excess return volatility (annualized) ----------------------------------
RVOL <- as.data.frame(apply(GW_dataset, 2, as.numeric)) %>% 
        filter(yyyymm >= 197201) %>% 
        select(yyyymm, CRSP_SPvw, Rfree) %>% 
        mutate(R_F_lag = lag(Rfree)) %>% 
        rename(SP_R = CRSP_SPvw) %>% 
        select(yyyymm, SP_R, R_F_lag) %>% 
        na.omit() 

RVOL_mat <- matrix(NaN, nrow = nrow(RVOL) - 11, ncol = 1)

for(i in seq_along(RVOL_mat)){
  
  RVOL_mat[i,] <- mean(abs(RVOL$SP_R[i:(i+11)] - RVOL$R_F_lag[i:(i+11)]))
  
}

RVOL_mat <- sqrt(pi/2)*sqrt(12)*RVOL_mat

# merge variables --------------------------------------------------------------
GW_predictor <- cbind(GW_varset, RVOL_mat) %>% 
                rename(RVOL = RVOL_mat) %>% 
                select(log_DP, log_DY, log_EP, log_DE, RVOL, BM, NTIS, TBL, LTY, 
                       LTR, TMS, DFY, DFR, INFL_lag)

rm(GW_varset, RVOL, RVOL_mat)

# Load equity risk premium data, 1973:01-2014:12 -------------------------------
equity_risk <- data.frame(apply(GW_dataset, 2, as.numeric)) %>% 
               filter(yyyymm >= 197212) %>% 
               select(yyyymm, Rfree, CRSP_SPvw) %>%
               rename(R_SP500 = CRSP_SPvw) %>% 
               mutate(ER = R_SP500 - lag(Rfree), Rfree_lag = lag(Rfree)) %>% 
               na.omit()

# short-interest rate ----------------------------------------------------------
short_interest <- read.xlsx("data/Returns_short_interest_data.xlsx", 
                             sheetName = "Short interest") %>% 
                  select(Date, EWSI) %>% 
                  mutate(log_EWSI = log(EWSI)) %>% 
                  select(Date, log_EWSI)

# Compute cumulative excess returns --------------------------------------------
h <- c(1, 3, 6, 12)

ER_h <- matrix(NaN, nrow = nrow(equity_risk), ncol = NROW(h))
R_f_h <- matrix(NaN, nrow = nrow(equity_risk), ncol = NROW(h))

for(i in seq_along(h)){
  for(j in 1:(nrow(equity_risk) - (h[i] - 1))){
    
    ER_h[j,i] <- prod(1 + equity_risk$R_SP500[j:(j+h[i]-1)]) - prod(1 + equity_risk$Rfree_lag[j:(j+h[i]-1)])
    R_f_h[j,i] <- prod(1 + equity_risk$Rfree_lag[j:(j+h[i]-1)]) - 1
    
  }
}


# Load function Aset Allocation ------------------------------------------------
source("Asset_Allocation_Function.R")

results_list <- list()
results <- matrix(0, nrow = 14, ncol = 9)

col_names <- c("mu_pred", "sd_pred", "Sr_pred",
               "mu_bench", "sd_bench", "Sr_bench",
               "CER_pred", "CER_bench", "CER_Var")

row_names <- c("DP", "DY", "EP", "DE", "RVOL", "BM", "NTIS", "TBL", "LTY", 
               "LTR", "TMS", "DFY", "DFR", "INFL") 

list_nanes <- c("Horizon_1", "Horizon_2", "Horizon_3", "Horizon_4")

h <- c(1, 3, 6, 12)

for(z in 1:4){
  for(i in 1:14){
    
    asset_allocation(equity_risk, GW_predictor[,i], 1989, 1972, 10, h, 3, -0.5, 1.5)
    
    metrics <- list(mean_results, Std_results, Sharpe_results, 
                    mean_results, Std_results, Sharpe_results, 
                    CER_PREDIC_results, CER_BENCH_results, CER_gain_results)
    
      for(j in 1:9){
        
      aux_metrics <- metrics[[j]]
        
      if(j <= 3){
          
        results[i,j] <- aux_metrics[2,z]
          
      }else{
          
        results[i,j] <- aux_metrics[1,z]
        
      }
    }
  }
  
  colnames(results) <- col_names
  row.names(results) <- row_names
  
  results_list[[z]] <- results   
  
}

names(results_list) <- list_nanes

lapply(results_list, function(x) round(x, 3))
