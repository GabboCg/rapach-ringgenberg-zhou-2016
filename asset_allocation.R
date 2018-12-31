#'---
#' author: Gabriel Cabrera 
#' title: Asset Allocation
#' date: 27-12-2018  (last updated: 27-12-2018)
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
   
     ER_h[j,i] <- prod(1 + equity_risk$R_SP500[j:(j+h[i]-1)] ) - prod(1 + equity_risk$Rfree_lag[j:(j+h[i]-1)])
     R_f_h[j,i] <- prod(1 + equity_risk$Rfree_lag[j:(j+h[i]-1)]) - 1
     
  }
}

# Take care of preliminaries ---------------------------------------------------
T <- nrow(equity_risk)
in_sample_end <- 1989
R <- (in_sample_end  - 1972)*12 # in-sample period
P <- T - R # out-of-sample period

FC_PM <- matrix(NaN, nrow = P, ncol = NROW(h))
w_PM <- matrix(NaN, nrow = P, ncol = NROW(h))
R_PM <- matrix(NaN, nrow = P, ncol = NROW(h))
ER_PM <- matrix(NaN, nrow = P, ncol = NROW(h))

FC_PR <- array(NaN,c(P, ncol(GW_predictor) + 1, NROW(h)))
w_PR <- array(NaN,c(P, ncol(GW_predictor) + 1, NROW(h)))
R_PR <- array(NaN,c(P, ncol(GW_predictor) + 1, NROW(h)))
ER_PR <- array(NaN,c(P, ncol(GW_predictor) + 1, NROW(h)))

R_BH <- matrix(NaN, nrow = P, ncol = NROW(h))
ER_BH <- matrix(NaN, nrow = P, ncol = NROW(h))
FC_vol <- matrix(NaN, nrow = P, ncol = NROW(h))

window <- 12*10
RRA <- 3
w_LB <- -0.5
w_UB <- 1.5

# Compute out-of-sample forecasts ----------------------------------------------
for(p in 1:P){
  print(p)
  for(j in 1:NROW(h)){
    
    # volatility
    if(R+p-h[j] <= window-1){
      
      FC_vol[p,j] <- sd(ER_h[1:(R+p-h[j]),j])
    
    }else{
      
      FC_vol[p,j] <- sd(ER_h[((R+p-h[j])-(window-1)):(R+p-h[j]),j])
      
    }
    
    # Prevailing mean benchmark
    FC_PM[p,j] <- mean(ER_h[1:(R+p-h[j]),j])
    
    # Predictive regressions
    for(i in 1:(ncol(GW_predictor)+1)){
      if(i <= ncol(GW_predictor)){
        
      X_i_j_p <- cbind(matrix(1, nrow = (R+(p-1)-h[j]), ncol = 1),
                       GW_predictor[1:(R+(p-1)-h[j]),i])
      results_i_j_p <- lm(ER_h[2:(R+p-h[j]),j] ~ X_i_j_p[,2])
      FC_PR[p,i,j]= cbind(1, GW_predictor[R+(p-1),i])%*%results_i_j_p$coefficients
        
      }else{
      
        trend_p <- as.data.frame(matrix(0, nrow = (R+(p-1)), ncol = 1)) %>% 
                   mutate(V1 = as.numeric(as.numeric(rownames(.)))) %>% 
                   as.matrix()
        
        X_linear_p <- cbind(matrix(1, nrow = (R+(p-1)), ncol = 1), trend_p)
        results_linear_p <- lm(short_interest$log_EWSI[1:(R+(p-1))] ~ X_linear_p[,2]) 
        SII_p <- scale(results_linear_p$residuals)
        X_SII_j_p <- cbind(matrix(1, nrow = (R+(p-1)-h[j]), ncol = 1), SII_p[1:(R+(p-1)-h[j])])
        results_SII_j_p <- lm(ER_h[2:(R+p-h[j]),j] ~ X_SII_j_p[,2])
        FC_PR[p,i,j] <- cbind(1, tail(SII_p,1))%*%results_SII_j_p$coefficients
        
      }
    }
  }
}

# Computing portfolio weights/returns ------------------------------------------
for(j in 1:NROW(h)){
  for(t in 1:(P/h[j])){
   
    FC_vol_j_t <- FC_vol[((t-1)*h[j]+1),j]
    FC_PM_j_t <- FC_PM[((t-1)*h[j]+1),j]
    
    w_PM_j_t <- (1/RRA)*FC_PM_j_t/(FC_vol_j_t)^2
    
    ifelse(w_PM_j_t > w_UB,
           w_PM[((t-1)*h[j]+1),j] <- w_UB,
           ifelse(w_PM_j_t < w_LB,
                  w_PM[((t-1)*h[j]+1),j] <- w_LB,
                  w_PM[((t-1)*h[j]+1),j] <- w_PM_j_t))
    
    R_PM[((t-1)*h[j]+1),j] <- R_f_h[(R+(t-1)*h[j]+1),j] + w_PM[((t-1)*h[j]+1),j]*ER_h[(R+(t-1)*h[j]+1),j]
    ER_PM[((t-1)*h[j]+1),j] <- R_PM[((t-1)*h[j]+1),j] - R_f_h[(R+(t-1)*h[j]+1),j]
    
    for(i in 1:ncol(FC_PR)){
      
      FC_PR_i_j_t <- FC_PR[((t-1)*h[j]+1),i,j]
      w_PR_i_j_t <- (1/RRA)*FC_PR_i_j_t/(FC_vol_j_t)^2
      
      ifelse(w_PR_i_j_t > w_UB,
             w_PR[((t-1)*h[j]+1),i,j] <- w_UB,
             ifelse(w_PR_i_j_t < w_LB,
                    w_PR[((t-1)*h[j]+1),i,j] <- w_LB,
                    w_PR[((t-1)*h[j]+1),i,j] <- w_PR_i_j_t))
      
     R_PR[((t-1)*h[j]+1),i,j] <- R_f_h[(R+(t-1)*h[j]+1),j] + w_PR[((t-1)*h[j]+1),i,j]*ER_h[(R+(t-1)*h[j]+1),j]
     ER_PR[((t-1)*h[j]+1),i,j] <- R_PR[((t-1)*h[j]+1),i,j] - R_f_h[(R+(t-1)*h[j]+1),j]
      
    }
    
    R_BH[((t-1)*h[j]+1),j] <- R_f_h[(R+(t-1)*h[j]+1),j] + ER_h[(R+(t-1)*h[j]+1),j]
    ER_BH[((t-1)*h[j]+1),j] <- ER_h[(R+(t-1)*h[j]+1),j]
    
  }
}


# Compute CER gains and Sharpe ratios for full OOS period ----------------------
CER_gain <- matrix(NaN, nrow = (ncol(FC_PR)+1), ncol = NROW(h))
Sharpe <- matrix(NaN, nrow = (ncol(FC_PR)+2), ncol = NROW(h))

for(j in 1:NROW(h)){
  
  R_PM_j <- R_PM[,j]
  R_PM_j <- na.omit(R_PM_j)
  ER_PM_j <- ER_PM[,j]
  ER_PM_j <- na.omit(ER_PM_j)
  CER_PM_j <- (12/h[j])*(mean(R_PM_j) - 0.5*RRA*sd(R_PM_j)^2)
  Sharpe[1,j] <- sqrt((12/h[j]))*mean(ER_PM_j)/sd(ER_PM_j)
  
  for(i in 1:NCOL(FC_PR)){
    
    R_PR_i_j <- R_PR[,i,j]
    R_PR_i_j <- na.omit(R_PR_i_j)
    ER_PR_i_j <- ER_PR[,i,j]
    ER_PR_i_j <- na.omit(ER_PR_i_j)
    CER_PR_i_j <- (12/h[j])*(mean(R_PR_i_j) - 0.5*RRA*sd(R_PR_i_j)^2)
    CER_gain[i,j] <- 100*(CER_PR_i_j - CER_PM_j)
    Sharpe[i+1,j] <- sqrt((12/h[j]))*mean(ER_PR_i_j)/sd(ER_PR_i_j)
    
  }
  
  R_BH_j <- R_BH[,j]
  R_BH_j <- na.omit(R_BH_j)
  ER_BH_j <- ER_BH[,j]
  ER_BH_j <- na.omit(ER_BH_j)
  CER_BH_j <- (12/h[j])*(mean(R_BH_j) - 0.5*RRA*sd(R_BH_j)^2)
  CER_gain[NROW(CER_gain),j] <- 100*(CER_BH_j - CER_PM_j)
  Sharpe[NROW(Sharpe),j] <- sqrt((12/h[j]))*mean(ER_BH_j)/sd(ER_BH_j)
  
}

CER_gain
Sharpe


# Compute CER gains and Sharpe ratios for Global Financial Crisis period -------
GFC_start <- (2006-1989)*12 + 1
CER_gain_GFC <- matrix(NaN, nrow = (ncol(FC_PR)+1), ncol = NROW(h))
Sharpe_GFC <- matrix(NaN, nrow = (ncol(FC_PR)+2), ncol = NROW(h))

for(j in 1:NROW(h)){
  
  R_PM_j <- R_PM[GFC_start:NROW(R_PM),j]
  R_PM_j <- na.omit(R_PM_j)
  ER_PM_j <- ER_PM[GFC_start:NROW(ER_PM),j]
  ER_PM_j <- na.omit(ER_PM_j)
  CER_PM_j <- (12/h[j])*(mean(R_PM_j) - 0.5*RRA*sd(R_PM_j)^2)
  Sharpe_GFC[1,j] <- sqrt((12/h[j]))*mean(ER_PM_j)/sd(ER_PM_j)
  
  for(i in 1:NCOL(FC_PR)){
    
    R_PR_i_j <- R_PR[GFC_start:NROW(R_PR),i,j]
    R_PR_i_j <- na.omit(R_PR_i_j)
    ER_PR_i_j <- ER_PR[GFC_start:NROW(ER_PR),i,j]
    ER_PR_i_j <- na.omit(ER_PR_i_j)
    CER_PR_i_j <- (12/h[j])*(mean(R_PR_i_j) - 0.5*RRA*sd(R_PR_i_j)^2)
    CER_gain_GFC[i,j] <- 100*(CER_PR_i_j - CER_PM_j)
    Sharpe_GFC[i+1,j] <- sqrt((12/h[j]))*mean(ER_PR_i_j)/sd(ER_PR_i_j)
    
  }
  
  R_BH_j <- R_BH[GFC_start:NROW(R_BH),j]
  R_BH_j <- na.omit(R_BH_j)
  ER_BH_j <- ER_BH[GFC_start:NROW(ER_BH),j]
  ER_BH_j <- na.omit(ER_BH_j)
  CER_BH_j <- (12/h[j])*(mean(R_BH_j) - 0.5*RRA*sd(R_BH_j)^2)
  CER_gain_GFC[NROW(CER_gain_GFC),j] <- 100*(CER_BH_j - CER_PM_j)
  Sharpe_GFC[NROW(Sharpe_GFC),j] <- sqrt((12/h[j]))*mean(ER_BH_j)/sd(ER_BH_j)
  
}

