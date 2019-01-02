#'---
#' author: Gabriel Cabrera 
#' title: Asset Allocation Function
#' date: 30-12-2018  (last updated: 02-91-2019)
#'---

asset_allocation <- function(forecast_variable, predictor, in_sample_end, 
                             start_sample, window, forecast_horizon, RRA, 
                             w_LB, w_UB){
  
  h <- forecast_horizon
  
  # Take care of preliminaries -------------------------------------------------
  T <- nrow(forecast_variable)
  in_sample_end <- in_sample_end
  R <- (in_sample_end  - start_sample)*12 # in-sample period 204
  P <- T - R # out-of-sample period 300
  
  FC_PM <- matrix(NaN, nrow = P, ncol = NROW(h))
  w_PM <- matrix(NaN, nrow = P, ncol = NROW(h))
  R_PM <- matrix(NaN, nrow = P, ncol = NROW(h))
  ER_PM <- matrix(NaN, nrow = P, ncol = NROW(h))
  
  FC_PR <- array(NaN,c(P, NCOL(predictor), NROW(h)))
  w_PR <- array(NaN,c(P, NCOL(predictor), NROW(h)))
  R_PR <- array(NaN,c(P, NCOL(predictor), NROW(h)))
  ER_PR <- array(NaN,c(P, NCOL(predictor), NROW(h)))
  
  R_BH <- matrix(NaN, nrow = P, ncol = NROW(h))
  ER_BH <- matrix(NaN, nrow = P, ncol = NROW(h))
  FC_vol <- matrix(NaN, nrow = P, ncol = NROW(h))
  
  window <- 12*window
  RRA <- 3
  w_LB <- w_LB
  w_UB <- w_UB
  
  # Compute out-of-sample forecasts --------------------------------------------
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
      for(i in 1:(NCOL(predictor)+1)){
        if(i <= NCOL(predictor)){
          
          X_i_j_p <- cbind(matrix(1, nrow = (R+(p-1)-h[j]), ncol = 1),
                           predictor[1:(R+(p-1)-h[j])])
          results_i_j_p <- lm(ER_h[2:(R+p-h[j]),j] ~ X_i_j_p[,2])
          FC_PR[p,i,j] <- cbind(1, predictor[R+(p-1)])%*%results_i_j_p$coefficients
          
        }
      }
    }
  }
  
  # Computing portfolio weights/returns ----------------------------------------
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
  
  
  # Compute CER gains and Sharpe ratios for full OOS period --------------------
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
  
  # results
  CER_gain_results <<- CER_gain
  Sharpe_results <<- Sharpe 
}
