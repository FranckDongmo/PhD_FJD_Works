setwd("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\Code Manuela2\\Plots")
rm(list = ls()) #remet à jour tout l'espace de travail
reset = function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="",axes=FALSE)
}

#Control implementation MDA & SMC
create_control= function(t_begin_Camp, Number_of_cycle, VectTime_between_cycles, 
                         Dur_cycle, Gap, time, Nt, nb_years = 6) {
  result= rep(0, Nt)
  if (Number_of_cycle < 1) return(list("control" = result))
  days_per_year= 360
  
  for (t in 1:Nt) {
    current_year= floor(time[t] / days_per_year)
    if (current_year >= nb_years) next
    time_in_year= time[t] - current_year * days_per_year
    
    for (k in 0:(Number_of_cycle - 1)) { 
      ta= if (k == 0) t_begin_Camp else t_begin_Camp + sum(Dur_cycle + VectTime_between_cycles[1:k])
      tb= ta + Dur_cycle - 1
      
      if (ta - Gap <= time_in_year && time_in_year < ta) {result[t]= (time_in_year - ta + Gap) / Gap}
      else if (ta <= time_in_year && time_in_year <= tb) { result[t]= 1 }
      else if (tb < time_in_year && time_in_year <= tb + Gap) { result[t]= (tb + Gap - time_in_year) / Gap}
    }
  }
  return(list("control" = result))
}

#eta et kappa pour MDA & SMC
Efficacy_fun= function(sigma, T50, shape) {
  eta=1- 0.95 / (1 + (sigma / T50)^shape)
  kappa= shape / T50 * ((sigma / T50)^(shape - 1)) / (1 + (sigma / T50)^shape)
  return(list(eta = eta, kappa = kappa))
}

initialize_age_parameters= function(Va, Na, p_f) {
  # Mortalité naturelle Ref Quentin Richard 2024
  muh_values= c(66.8, 7.6, 1.7, 0.9, 1.3, 1.9, 2.4, 2.8, 3.6, 4.7, 6.3, 8.9,
                13.2, 19.8, 31.1, 47.7, 71.3, 110.5, 186.7) / 1000
  age_labels0= c(0, 1, seq(5, 90, by = 5))
  muh= rep(0, Na)
  for (i in 1:Na) {
    age_group= findInterval(Va[i], age_labels0, rightmost.closed = TRUE)
    muh[i]= muh_values[min(age_group, length(muh_values))]
  }
  
  # Taux d'infection Ref Quentin Richard 2024
  alpha_1= 0.045; alpha_2= 0.181; #alpha_1 =0.071[0.023;0.175],  #alpha_2=0.302[0.16; 0.475]
  G= function(a) 22.7 * a * exp(-0.0934 * a)
  betah_s=  alpha_1 * (G(Va)^alpha_2)
  betah_r= betah_s
  # Progression vers cas clinique Ref chitnis Malaria 2008
  nu_0= 1/10; nu_1= 1/150; nu_2= 1/5
  nuh_m= nuh_f= rep(0, Na)
  for (i in 1:Na) {
    if (Va[i] <= 5) {  nuh_m[i]= nuh_f[i]= nu_2 }
    else if (Va[i] <= 15) { nuh_m[i]= nuh_f[i]= nu_0 }
    else { nuh_m[i]= nu_1; nuh_f[i]= if (Va[i] <= 45) p_f * nu_0 + (1 - p_f) * nu_1 else nu_1
    }
  }
  
  # Taux de guérison, gammah_su= Ref Chitnis Malaria 2008
  Recov=3 # les humains treated guerissent Recov fois plus vite
  #Ref pour Recov [2,3] ; Sumba PO, Wong SL, et al 2008, Recovery from Malaria Kenya
  bar_gamma= 0.7; gammah_su= rep(0.0035, Na); gammah_st= Recov*gammah_su 
  gammah_rt= bar_gamma * gammah_st; gammah_ru= bar_gamma * gammah_su;
  gammahM_s=Recov*gammah_su; gammahM_r= bar_gamma * gammah_st;
  gammahS_s=Recov*gammah_su; gammahS_r= bar_gamma * gammah_st;
  # Mortalité induite par la maladie Ref Quentin Richard 2024
  deltah_values= c(1.07e-3, 7.02e-4, 4.55e-4, 5.73e-5)
  bar_delta= 1.5; deltah_su= deltah_rt= rep(0, Na)
  
  age_breaks= c(0, 1, 5, 15, Inf)
  for (i in 1:Na) {
    age_group= findInterval(Va[i], age_breaks, rightmost.closed = TRUE)
    deltah_su[i]= deltah_values[min(age_group, length(deltah_values))]
    deltah_rt[i]= bar_delta * deltah_su[i]
  }
  return(list(muh = muh, betah_s = betah_s, betah_r=betah_r, nuh_m = nuh_m, nuh_f = nuh_f,gammah_su = gammah_su,
              gammah_ru=gammah_ru, gammah_rt=gammah_rt, gammah_st = gammah_st, deltah_su = deltah_su, deltah_rt = deltah_rt,
              gammahM_s=gammahM_s, gammahM_r=gammahM_r, gammahS_s=gammahS_s, gammahS_r=gammahS_r
  ))
}

Mainfunction= function(MDA_Dur_cycle, MDA_VectTime_between_cycles,
                       t_begin_Camp, Number_of_cycle, SMC_VectTime_between_cycles,
                       SMC_Dur_cycle, PropSMC, PropMDA, Gap, Nyr, LengYr,  PropTreatment, 
                       dt, p_f, a_lim, age_cible, hbr, epsi_s, bar_epss, InitPrevR,
                       initial_conditions = NULL) {
  # Efficacité MDA et SMC
  dsigmaSmc= 1; VsigmaSmc= seq(0, 70, by = dsigmaSmc); NsigmaSmc= length(VsigmaSmc)
  dsigmaMda= 1; VsigmaMda= seq(0, 70, by = dsigmaMda); NsigmaMda= length(VsigmaMda)
  Eff_MDA_50 = 35; shape_MDA = 8; Eff_SMC_50 = 28; shape_SMC = 6.5
  
  eff_MDA= Efficacy_fun(VsigmaMda, Eff_MDA_50, shape_MDA) 
  eff_SMC= Efficacy_fun(VsigmaSmc, Eff_SMC_50, shape_SMC)
  eta_M= eff_MDA$eta; k_M= eff_MDA$kappa; Eff_M_loss=1-eta_M
  eta_S= eff_SMC$eta; k_S= eff_SMC$kappa; Eff_S_loss=1-eta_S
  
  # Variables et Paramètres fixes
  {
    Tmax= Nyr * LengYr; time= seq(0, Tmax, by = dt);Nt= length(time)
    da= dt ; Va= seq(0, 90, by = da); Na= length(Va)
    dsigmaR= 1; VsigmaR= seq(0, 60, by = dsigmaR); NsigmaR= length(VsigmaR)
    
    bar_bm= 0.5; bar_deltam= 0.95; Lh= 30750; Lm= 2.1e5;mum= 0.13; betam_s= 0.25; deltam_s= 1/20; 
    deltah_st_base= 0; deltah_ru_base= 9e-5; rhoh_base= 5.5e-4; p= 0.52;
    
    deltam_r= bar_deltam * deltam_s; betam_r= bar_bm * betam_s
    epsi_r= bar_epss * epsi_s; betam= diag(c(betam_s, betam_r))
    deltam= diag(c(deltam_s, deltam_r))
    
    # Initialisation des paramètres dépendant de l'âge
    age_params= initialize_age_parameters(Va, Na, p_f); deltah_st= rep(deltah_st_base, Na)
    deltah_ru= rep(deltah_ru_base, Na); rhoh= rep(rhoh_base, NsigmaR)
    
    epsi= matrix(c(1-epsi_s, epsi_r, 1, 0, epsi_s,1-epsi_r,0,1), ncol = 4, byrow = TRUE)
    f= matrix(1, ncol = 4, byrow = TRUE); m= matrix(1, ncol = 2, byrow = TRUE)
    
    # Matrice theta pour traitement
    MatTheta= array(NA, dim = c(4, 2, Nt))
    for (t in 1:Nt) {
      MatTheta[,,t]= matrix(c(PropTreatment, 0, 0, PropTreatment,
                              1-PropTreatment, 0, 0, 1-PropTreatment), ncol = 2, byrow = TRUE)
    }
    
    barbetahS <- array(NA, dim = c(2, 2, Na, NsigmaSmc))
    barbetahM <- array(NA, dim = c(2, 2, Na, NsigmaMda))
    betah <- array(NA, dim = c(2, 2, Na))
    deltah <- array(NA, dim = c(4, 4, Na))
    gammah <- array(NA, dim = c(4, 4, Na))
    gammahS <- array(NA, dim = c(2, 2, Na, NsigmaSmc))
    gammahM <- array(NA, dim = c(2, 2, Na, NsigmaMda))
    
    for (a in 1:Na) {
      betah[,,a] <- diag(c(age_params$betah_s[a], age_params$betah_r[a]))
      deltah[,,a] <- diag(c(deltah_st[a], age_params$deltah_rt[a],
                            age_params$deltah_su[a], deltah_ru[a]))
      gammah[,,a] <- diag(c(age_params$gammah_st[a], age_params$gammah_rt[a],
                            age_params$gammah_su[a], age_params$gammah_ru[a]))
      
      ## Paramètres dépendant de l'âge ET du temps depuis traitement SMC et MDA
      for (sigmaSmc in 1:NsigmaSmc) {
        barbetahS[,,a,sigmaSmc] <- diag(c(age_params$betah_s[a] * (1 - Eff_S_loss[sigmaSmc]),
                                          age_params$betah_r[a] * (1 - Eff_S_loss[sigmaSmc])))
        gammahS[,,a,sigmaSmc]=diag(c(age_params$gammah_su[a]+(age_params$gammah_st[a]-age_params$gammah_su[a])*Eff_S_loss[sigmaSmc],
                                     age_params$gammah_ru[a]+(age_params$gammah_rt[a]-age_params$gammah_ru[a])*Eff_S_loss[sigmaSmc]))
      }
      
      for (sigmaMda in 1:NsigmaMda) {
        barbetahM[,,a,sigmaMda] <- diag(c(age_params$betah_s[a] * (1 - Eff_M_loss[sigmaMda]),
                                          age_params$betah_r[a] * (1 - Eff_M_loss[sigmaMda])))
        gammahM[,,a,sigmaMda]=diag(c(age_params$gammah_su[a]+(age_params$gammah_st[a]-age_params$gammah_su[a])*Eff_M_loss[sigmaMda],
                                     age_params$gammah_ru[a]+(age_params$gammah_rt[a]-age_params$gammah_ru[a])*Eff_M_loss[sigmaMda]))
      }
    }
  }
  
  # Stratégies SMC et MDA
  ModelSMC= create_control(t_begin_Camp, Number_of_cycle, SMC_VectTime_between_cycles, 
                           SMC_Dur_cycle, Gap, time, Nt, Nyr)
  ModelMDA= create_control(t_begin_Camp, Number_of_cycle, MDA_VectTime_between_cycles, 
                           MDA_Dur_cycle, Gap, time, Nt, Nyr)
  
  bar_phi=  ModelSMC$control
  bar_psi=  ModelMDA$control
  phi= psi_f= psi_m= matrix(0, nrow = Na, ncol = Nt)
  if (age_cible == 0) { ## dans tous le travail finalement age_cible=0
    for (t in 1:Nt) {
      for (a in 1:Na) {
        if (Va[a] <= a_lim) {  phi[a, t] = (PropSMC / SMC_Dur_cycle) * bar_phi[t]  }
        if (Va[a] > a_lim) {
          psi_m[a, t] = (PropMDA / MDA_Dur_cycle) * bar_psi[t]
          if (Va[a] <= 15 || Va[a] >= 45) { psi_f[a, t] = (PropMDA / MDA_Dur_cycle) * bar_psi[t] }
        }
      }
    }
  } else {
    Prop1 = 0.95; Prop2 = 0.8; Prop3 = 0.6
    
    for (t in 1:Nt) {
      for (a in 1:Na) {
        if (Va[a] <= a_lim) { phi[a, t] = (Prop1 / SMC_Dur_cycle) * bar_phi[t]}
        if (Va[a] > a_lim & Va[a] <= 15) {
          psi_m[a, t] = (Prop1 / MDA_Dur_cycle) * bar_psi[t]
          psi_f[a, t] = (Prop1 / MDA_Dur_cycle) * bar_psi[t]
        }
        if (Va[a] > 15 & Va[a] <= 40) { psi_m[a, t] = (Prop2 / MDA_Dur_cycle) * bar_psi[t] }
        if (Va[a] > 40) { psi_m[a, t] = (Prop3 / MDA_Dur_cycle) * bar_psi[t]}
        if (Va[a] > 45) { psi_f[a, t] = (Prop3 / MDA_Dur_cycle) * bar_psi[t]}
      }
    }
  }
  
  # Initialisation des variables d'état
  mSh= fSh= matrix(0, nrow = Na, ncol = Nt); mAh= fAh= array(0, dim = c(2, Na, Nt))
  mIh= fIh= array(0, dim = c(4, Na, Nt)); mRh= fRh= array(0, dim = c(NsigmaR, Na, Nt))
  
  mSh_smc= fSh_smc= array(0,dim=c(NsigmaSmc,Na,Nt));mAh_smc=fAh_smc= array(0, dim = c(2, NsigmaSmc, Na, Nt))
  mSh_mda= fSh_mda= array(0, dim = c(NsigmaMda, Na, Nt))
  mAh_mda= fAh_mda= array(0, dim = c(2, NsigmaMda, Na, Nt))
  
  Sm= Im= Dh= Nh= Nm= numeric(Nt)
  Im= matrix(0, nrow = 2, ncol = Nt); lambdam= lambdah= matrix(0, nrow = 2, ncol = Nt)
  
  # Conditions initiales
  if (is.null(initial_conditions)) {
    Int_values= 8136.10 * (c(12.9, 12.5, 11.5, 13.1, 11.9, 9.3, 7.3, 5.5, 4.4,
                             3.2, 2.6, 1.8, 1.4, 0.9, 0.7, 0.3, 0.2, 0.2) / 0.997)
    Prev2= 0.0015; Prev= 0.004; PrevI= 1/4
    
    for (i in seq_along(Int_values)) {
      if (i == 1) {
        idxs= which(Va <= 5)
      } else if (i == length(Int_values)) {
        idxs= which(Va > 85)
      } else {
        idxs= which(Va > (i-1)*5 & Va <= i*5)
      }
      
      if (length(idxs) > 0) {
        n_group= length(idxs)
        mSh[idxs, 1]= (1-Prev-Prev2) * (1-p) * Int_values[i] / n_group
        fSh[idxs, 1]= (1-Prev-Prev2) * p * Int_values[i] / n_group
        mAh[, idxs, 1]= (1-PrevI) * Prev * c((1-p) * Int_values[i] * (1-InitPrevR),
                                             (1-p) * Int_values[i] * InitPrevR) / n_group
        fAh[, idxs, 1]= (1-PrevI) * Prev * c(p * Int_values[i] * (1-InitPrevR),
                                             p * Int_values[i] * InitPrevR) / n_group
        mIh[, idxs, 1]= PrevI*Prev * c((1-p)*Int_values[i]*(1-InitPrevR), (1-p)*Int_values[i]*InitPrevR,
                                       (1-p)*Int_values[i]*(1-InitPrevR), (1-p)*Int_values[i]*InitPrevR)/n_group
        fIh[, idxs, 1]= PrevI*Prev * c(p*Int_values[i]*(1-InitPrevR), p*Int_values[i]*InitPrevR,
                                       p*Int_values[i]*(1-InitPrevR), p*Int_values[i]*InitPrevR)/n_group
        mRh[1, idxs, 1]= Prev2 * (1-p) * Int_values[i] / n_group
        fRh[1, idxs, 1]= Prev2 * p * Int_values[i] / n_group
      }
    }
    
    Dh[1]= 1; S0m= 1e6; Prev_m0= 0.1; Sm[1]= (1 - Prev_m0) * S0m
    Im[,1]= c((1 - InitPrevR) * Prev_m0 * S0m, InitPrevR * Prev_m0 * S0m)
  } else {
    # Charger les CI déjà fournies
    mSh[,1]  = initial_conditions$mSh_final; fSh[,1]  = initial_conditions$fSh_final
    mAh[,,1] = initial_conditions$mAh_final; fAh[,,1] = initial_conditions$fAh_final
    mIh[,,1] = initial_conditions$mIh_final; fIh[,,1] = initial_conditions$fIh_final
    mRh[,,1] = initial_conditions$mRh_final; fRh[,,1] = initial_conditions$fRh_final
    
    mSh_smc[,,1]= initial_conditions$mSh_smc_final; fSh_smc[,,1]= initial_conditions$fSh_smc_final
    mAh_smc[,,,1]= initial_conditions$mAh_smc_final; fAh_smc[,,,1]= initial_conditions$fAh_smc_final
    
    mSh_mda[,,1]= initial_conditions$mSh_mda_final; fSh_mda[,,1]= initial_conditions$fSh_mda_final
    mAh_mda[,,,1]= initial_conditions$mAh_mda_final; fAh_mda[,,,1]= initial_conditions$fAh_mda_final
    
    Sm[1] = initial_conditions$Sm_final; Im[,1]= initial_conditions$Im_final
    Dh[1] = initial_conditions$Dh_final
  }
  
  t=1
  Nh[t] =  sum(mSh[,t]+fSh[,t]) + sum(mAh[, , t]+fAh[ , , t]) + sum(mIh[, , t]+fIh[ , , t]) + 
    sum(mRh[, , t]+fRh[ , , t]) + sum(mSh_mda[,,t]+fSh_mda[,,t]) + sum(mSh_smc[,,t]+fSh_smc[,,t]) + 
    sum(mAh_mda[,,,t]+fAh_mda[,,,t])  + sum(mAh_smc[,,,t]+fAh_smc[,,,t])
  Nm[t] = Sm[t] + sum(Im[ , t])
  
  for (t in 1:(Nt - 1)) {
    lambdam[,t]= hbr[t] *betam %*% Im[,t] / Nh[t]
  
    lambdah[,t] = lambdah_sum = numeric(2)
    for(a in 1:Na) {
      term1 = betah[,,a] %*% (mAh[,a,t] + fAh[,a,t] + epsi %*% (mIh[,a,t]+fIh[,a,t]))
      term2 = term3 = numeric(2)
      for (sigmaSmc in 1:NsigmaSmc) {
        smc_sum = mAh_smc[,sigmaSmc,a,t] + fAh_smc[,sigmaSmc,a,t]
        term2 = term2 + barbetahS[,,a,sigmaSmc] %*% smc_sum 
      }
      for (sigmaMda in 1:NsigmaMda) {
        mda_sum = mAh_mda[,sigmaMda,a,t] + fAh_mda[,sigmaMda,a,t]
        term3 = term3 + barbetahM[,,a,sigmaMda] %*% mda_sum
      }
      lambdah_sum = lambdah_sum + term1 + term2 + term3
    }   
    lambdah[,t] = hbr[t]*lambdah_sum / Nh[t]
    
    for (a in 1:Na) {
      if (a == 1) {
        mSh[a,t+1]= ((mSh[a,t]/dt) + ((1-p)*Lh/da) + sum(rhoh * mRh[,a,t]) + sum(k_S * mSh_smc[,a,t]) + 
                       sum(k_M * mSh_mda[,a,t])) /((1/dt) + (1/da) + age_params$muh[a] +sum(lambdam[,t]) + phi[a,t] + psi_m[a,t])
        fSh[a,t+1]= ((fSh[a,t]/dt) + (p*Lh/da) + sum(rhoh * fRh[,a,t]) + sum(k_S * fSh_smc[,a,t]) + 
                       sum(k_M * fSh_mda[,a,t])) / ((1/dt) + (1/da) + age_params$muh[a] + sum(lambdam[,t]) + phi[a,t] + psi_f[a,t])
        
        mAh[,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+ age_params$muh[a] + age_params$nuh_m[a] +  phi[a,t] + psi_m[a,t]))%*%
          ((mAh[,a,t]/dt) + lambdam[,t] * mSh[a,t])
        fAh[,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+age_params$muh[a] +age_params$nuh_f[a] + phi[a,t] + psi_f[a,t])) %*% 
          ((fAh[,a,t]/dt) + lambdam[,t] * fSh[a,t])
        
        mIh[,a,t+1]= solve(diag(4)*((1/dt)+(1/da)+age_params$muh[a]) + deltah[,,a] + gammah[,,a]) %*%
          ((mIh[,a,t]/dt) + age_params$nuh_m[a] * MatTheta[,,t] %*% mAh[,a,t])
        fIh[,a,t+1]= solve(diag(4) * ((1/dt) + (1/da) + age_params$muh[a]) +  deltah[,,a] + gammah[,,a]) %*%
          ((fIh[,a,t]/dt) + age_params$nuh_f[a] * MatTheta[,,t] %*% fAh[,a,t])
      } else {
        mSh[a,t+1]= ((mSh[a,t]/dt) + (mSh[a-1,t+1]/da) + sum(rhoh * mRh[,a,t]) +sum(k_S * mSh_smc[,a,t]) +
                       sum(k_M * mSh_mda[,a,t])) / ((1/dt) + (1/da) + age_params$muh[a] + sum(lambdam[,t]) + phi[a,t] + psi_m[a,t])
        fSh[a,t+1]= ((fSh[a,t]/dt) + (fSh[a-1,t+1]/da) + sum(rhoh * fRh[,a,t]) +  sum(k_S * fSh_smc[,a,t]) + 
                       sum(k_M * fSh_mda[,a,t])) / ((1/dt) + (1/da) + age_params$muh[a] + sum(lambdam[,t]) +  phi[a,t] + psi_f[a,t])
        
        mAh[,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+age_params$muh[a]+age_params$nuh_m[a]+phi[a,t]+psi_m[a,t])) %*%
          ((mAh[,a,t]/dt) + (mAh[,a-1,t+1]/da) + lambdam[,t] * mSh[a,t])
        
        fAh[,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+age_params$muh[a] + age_params$nuh_f[a]+phi[a,t]+psi_f[a,t])) %*%
          ((fAh[,a,t]/dt) + (fAh[,a-1,t+1]/da) + lambdam[,t]*fSh[a,t])
        
        mIh[,a,t+1]= solve(diag(4)*((1/dt)+(1/da)+age_params$muh[a]) + deltah[,,a] + gammah[,,a]) %*%
          ((mIh[,a,t]/dt) + (mIh[,a-1,t+1]/da) +  age_params$nuh_m[a] * MatTheta[,,t] %*% mAh[,a,t])
        
        fIh[,a,t+1]= solve(diag(4)*((1/dt)+(1/da)+age_params$muh[a]) + deltah[,,a]+gammah[,,a]) %*%
          ((fIh[,a,t]/dt) + (fIh[,a-1,t+1]/da) + age_params$nuh_f[a] * MatTheta[,,t] %*% fAh[,a,t])
      }
      
      sigmaR= 1
      #Calculer la somme gammahM pour Ah_mda and Ah_smc
      gammahM_sum_m <- numeric(2); gammahM_sum_f <- numeric(2)
      gammahS_sum_m <- numeric(2); gammahS_sum_f <- numeric(2)
      for (sigmaMda in 1:NsigmaMda) {
        gammahM_sum_m <- gammahM_sum_m + gammahM[,,a,sigmaMda] %*% mAh_mda[,sigmaMda,a,t]
        gammahM_sum_f <- gammahM_sum_f + gammahM[,,a,sigmaMda] %*% fAh_mda[,sigmaMda,a,t]
      }
      for (sigmaSmc in 1:NsigmaSmc) {
        gammahS_sum_m <- gammahS_sum_m + gammahS[,,a,sigmaSmc] %*% mAh_mda[,sigmaSmc,a,t]
        gammahS_sum_f <- gammahS_sum_f + gammahS[,,a,sigmaSmc] %*% fAh_mda[,sigmaSmc,a,t]
      }
      
      if (a == 1) {
        mRh[sigmaR,a,t+1]= ((mRh[sigmaR,a,t]/dt)+(f %*% gammah[,,a] %*% mIh[,a,t] + m %*% gammahM_sum_m + m %*% gammahS_sum_m)/dsigmaR) /
          ((1/dt)+(1/da)+(1/dsigmaR) + age_params$muh[a] + rhoh[sigmaR] + phi[a,t] + psi_m[a,t])
        
        fRh[sigmaR,a,t+1]= ((fRh[sigmaR,a,t]/dt)+(f %*% gammah[,,a] %*% fIh[,a,t] + m %*% gammahM_sum_f + m %*% gammahS_sum_f)/dsigmaR) /
          ((1/dt)+(1/da) + (1/dsigmaR) + age_params$muh[a] + rhoh[sigmaR] + phi[a,t] + psi_f[a,t])
      } else {
        mRh[sigmaR,a,t+1]= ((mRh[sigmaR,a,t]/dt) + (mRh[sigmaR,a-1,t+1]/da) + 
                              (f %*% gammah[,,a] %*% mIh[,a,t] + m %*% gammahM_sum_m + m %*% gammahS_sum_m)/dsigmaR) /
          ((1/dt) + (1/da) + (1/dsigmaR) + age_params$muh[a] + rhoh[sigmaR] + phi[a,t] + psi_m[a,t])
        fRh[sigmaR,a,t+1]= ((fRh[sigmaR,a,t]/dt) + (fRh[sigmaR,a-1,t+1]/da) + 
                              (f %*% gammah[,,a] %*% fIh[,a,t] + m %*% gammahM_sum_f+ m %*% gammahS_sum_f)/dsigmaR) /
          ((1/dt) + (1/da) + (1/dsigmaR) + age_params$muh[a] + rhoh[sigmaR] + phi[a,t] + psi_f[a,t])
      }
      
      for (sigmaR in 2:NsigmaR) {
        if (a == 1) {
          mRh[sigmaR,a,t+1]= ((mRh[sigmaR,a,t]/dt)+(mRh[sigmaR-1,a,t+1]/dsigmaR))/((1/dt) + (1/da) + (1/dsigmaR)+age_params$muh[a] + 
                                                                                     rhoh[sigmaR] + phi[a,t] + psi_m[a,t])
          
          fRh[sigmaR,a,t+1]= ((fRh[sigmaR,a,t]/dt)+(fRh[sigmaR-1,a,t+1]/dsigmaR))/((1/dt) + (1/da) + (1/dsigmaR)+age_params$muh[a] + 
                                                                                       rhoh[sigmaR] + phi[a,t] + psi_f[a,t])
        } else {
          mRh[sigmaR,a,t+1]= ((mRh[sigmaR,a,t]/dt)+(mRh[sigmaR,a-1,t+1]/da) + (mRh[sigmaR-1,a,t+1]/dsigmaR))/
            ((1/dt) + (1/da) + (1/dsigmaR) + age_params$muh[a] +  rhoh[sigmaR] + phi[a,t] + psi_m[a,t])
          
          fRh[sigmaR,a,t+1]= ((fRh[sigmaR,a,t]/dt)+(fRh[sigmaR,a-1,t+1]/da) + (fRh[sigmaR-1,a,t+1]/dsigmaR))/
            ((1/dt) + (1/da) + (1/dsigmaR) + age_params$muh[a] +  rhoh[sigmaR] + phi[a,t] + psi_f[a,t])
        }
      }
    }
    
    # Compartiments avec SMC ET MDA
    for (a in 1:Na) {
      for (sigmaSmc in 1:NsigmaSmc) {
        if (a == 1) {
          if (sigmaSmc == 1) {
            mSh_smc[sigmaSmc,a,t+1]= ((mSh_smc[sigmaSmc,a,t]/dt)+(phi[a,t]*(mSh[a,t]+sum(mRh[,a,t]))/dsigmaSmc))/
              ((1/dt) + (1/da) + (1/dsigmaSmc) + age_params$muh[a] + eta_S[sigmaSmc]*sum(lambdam[,t]) + k_S[sigmaSmc])
            
            fSh_smc[sigmaSmc,a,t+1]= ((fSh_smc[sigmaSmc,a,t]/dt)+(phi[a,t]*(fSh[a,t]+sum(fRh[,a,t]))/dsigmaSmc)) /
              ((1/dt) + (1/da) + (1/dsigmaSmc) + age_params$muh[a] + eta_S[sigmaSmc]*sum(lambdam[,t]) + k_S[sigmaSmc])
            
            mAh_smc[,sigmaSmc,a,t+1]=solve(diag(2)*((1/dt) + (1/da) + (1/dsigmaSmc) + age_params$muh[a]) +gammahS[,,a,sigmaSmc]) %*%
              ((mAh_smc[,sigmaSmc,a,t]/dt)+(phi[a,t]*mAh[,a,t]/dsigmaSmc) + eta_S[sigmaSmc]*lambdam[,t]*mSh_smc[sigmaSmc,a,t])
            
            fAh_smc[,sigmaSmc,a,t+1]= solve(diag(2)*((1/dt)+(1/da) + (1/dsigmaSmc) + age_params$muh[a])+gammahS[,,a,sigmaSmc]) %*%
              ((fAh_smc[,sigmaSmc,a,t]/dt)+(phi[a,t]*fAh[,a,t]/dsigmaSmc) +eta_S[sigmaSmc]*lambdam[,t]*fSh_smc[sigmaSmc,a,t])
          } else {
            # Advection en sigma
            mSh_smc[sigmaSmc,a,t+1]= ((mSh_smc[sigmaSmc,a,t]/dt)+(mSh_smc[sigmaSmc-1,a,t+1]/dsigmaSmc)) /
              ((1/dt) + (1/da) + (1/dsigmaSmc) + age_params$muh[a] + eta_S[sigmaSmc]*sum(lambdam[,t]) + k_S[sigmaSmc])
            
            fSh_smc[sigmaSmc,a,t+1]= ((fSh_smc[sigmaSmc,a,t]/dt)+ (fSh_smc[sigmaSmc-1,a,t+1]/dsigmaSmc)) /
              ((1/dt) + (1/da) + (1/dsigmaSmc) + age_params$muh[a] +  eta_S[sigmaSmc]*sum(lambdam[,t]) + k_S[sigmaSmc])
            
            mAh_smc[,sigmaSmc,a,t+1]= solve(diag(2)*((1/dt) + (1/da) + (1/dsigmaSmc)+age_params$muh[a]) +gammahS[,,a,sigmaSmc]) %*%
              ((mAh_smc[,sigmaSmc,a,t]/dt)+ (mAh_smc[,sigmaSmc-1,a,t+1]/dsigmaSmc)+ eta_S[sigmaSmc]*lambdam[,t]*mSh_smc[sigmaSmc,a,t])
            
            fAh_smc[,sigmaSmc,a,t+1]= solve(diag(2)*((1/dt) + (1/da) + (1/dsigmaSmc)+age_params$muh[a]) +gammahS[,,a,sigmaSmc]) %*%
              ((fAh_smc[,sigmaSmc,a,t]/dt)+(fAh_smc[,sigmaSmc-1,a,t+1]/dsigmaSmc)+ eta_S[sigmaSmc]*lambdam[,t]*fSh_smc[sigmaSmc,a,t])
          }
        } else {
          if (sigmaSmc == 1) {
            mSh_smc[sigmaSmc,a,t+1]= ((mSh_smc[sigmaSmc,a,t]/dt) + (mSh_smc[sigmaSmc,a-1,t+1]/da) + 
                                        (phi[a,t] * (mSh[a,t] + sum(mRh[,a,t])) / dsigmaSmc)) /
              ((1/dt)+(1/da)+(1/dsigmaSmc)+age_params$muh[a]+ eta_S[sigmaSmc]*sum(lambdam[,t]) + k_S[sigmaSmc])
            
            fSh_smc[sigmaSmc,a,t+1]= ((fSh_smc[sigmaSmc,a,t]/dt) + (fSh_smc[sigmaSmc,a-1,t+1]/da) + 
                                        (phi[a,t] * (fSh[a,t] + sum(fRh[,a,t])) / dsigmaSmc)) /
              ((1/dt)+(1/da)+(1/dsigmaSmc)+age_params$muh[a]+ eta_S[sigmaSmc]*sum(lambdam[,t]) + k_S[sigmaSmc])
            
            mAh_smc[,sigmaSmc,a,t+1]= solve(diag(2)*((1/dt) + (1/da) + (1/dsigmaSmc) + age_params$muh[a])+gammahS[,,a,sigmaSmc]) %*%
              ((mAh_smc[,sigmaSmc,a,t]/dt)+(mAh_smc[,sigmaSmc,a-1,t+1]/da) + (phi[a,t]*mAh[,a,t]/dsigmaSmc) 
               + eta_S[sigmaSmc]*lambdam[,t]*mSh_smc[sigmaSmc,a,t])
            
            fAh_smc[,sigmaSmc,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaSmc) + age_params$muh[a]) +gammahS[,,a,sigmaSmc]) %*%
              ((fAh_smc[,sigmaSmc,a,t]/dt)+(fAh_smc[,sigmaSmc,a-1,t+1]/da)+(phi[a,t]*fAh[,a,t]/dsigmaSmc) + 
                 eta_S[sigmaSmc]*lambdam[,t]*fSh_smc[sigmaSmc,a,t])
          } else {
            mSh_smc[sigmaSmc,a,t+1]= ((mSh_smc[sigmaSmc,a,t]/dt) + (mSh_smc[sigmaSmc,a-1,t+1]/da) + (mSh_smc[sigmaSmc-1,a,t+1]/dsigmaSmc))/
              ((1/dt)+(1/da)+(1/dsigmaSmc) + age_params$muh[a]+eta_S[sigmaSmc]*sum(lambdam[,t])+k_S[sigmaSmc])
            
            fSh_smc[sigmaSmc,a,t+1]= ((fSh_smc[sigmaSmc,a,t]/dt) + (fSh_smc[sigmaSmc,a-1,t+1]/da) + (fSh_smc[sigmaSmc-1,a,t+1]/dsigmaSmc))/
              ((1/dt)+(1/da)+(1/dsigmaSmc) + age_params$muh[a]+eta_S[sigmaSmc]*sum(lambdam[,t]) + k_S[sigmaSmc])
            
            mAh_smc[,sigmaSmc,a,t+1]= solve(diag(2)*((1/dt) + (1/da) + (1/dsigmaSmc) + age_params$muh[a]) +gammahS[,,a,sigmaSmc]) %*%
              ((mAh_smc[,sigmaSmc,a,t]/dt)+(mAh_smc[,sigmaSmc,a-1,t+1]/da) + (mAh_smc[,sigmaSmc-1,a,t+1]/dsigmaSmc) 
               + eta_S[sigmaSmc]*lambdam[,t]*mSh_smc[sigmaSmc,a,t])
            
            fAh_smc[,sigmaSmc,a,t+1]= solve(diag(2)*((1/dt) + (1/da) + (1/dsigmaSmc) +age_params$muh[a]) +gammahS[,,a,sigmaSmc]) %*%
              ((fAh_smc[,sigmaSmc,a,t]/dt) + (fAh_smc[,sigmaSmc,a-1,t+1]/da) + 
                 (fAh_smc[,sigmaSmc-1,a,t+1]/dsigmaSmc) + eta_S[sigmaSmc]*lambdam[,t]*fSh_smc[sigmaSmc,a,t])
          }
        }
      }
      
      for (sigmaMda in 1:NsigmaMda) {
        if (a == 1) {
          if (sigmaMda == 1) {
            mSh_mda[sigmaMda,a,t+1]= ((mSh_mda[sigmaMda,a,t]/dt)+(psi_m[a,t]*(mSh[a,t]+sum(mRh[,a,t]))/dsigmaMda))/
              ((1/dt)+(1/da)+(1/dsigmaMda)+age_params$muh[a] + eta_M[sigmaMda]*sum(lambdam[,t]) + k_M[sigmaMda])
            
            fSh_mda[sigmaMda,a,t+1]= ((fSh_mda[sigmaMda,a,t]/dt) +(psi_f[a,t]*(fSh[a,t]+sum(fRh[,a,t]))/dsigmaMda))/
              ((1/dt) + (1/da) + (1/dsigmaMda) + age_params$muh[a] + eta_M[sigmaMda]*sum(lambdam[,t]) + k_M[sigmaMda])
            
            mAh_mda[,sigmaMda,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaMda)) + gammahM[,,a,sigmaMda] 
                                            + diag(2)*age_params$muh[a]) %*%((mAh_mda[,sigmaMda,a,t]/dt)+(psi_m[a,t]*mAh[,a,t]/dsigmaMda) + 
                 eta_M[sigmaMda]*lambdam[,t]*mSh_mda[sigmaMda,a,t])
            
            fAh_mda[,sigmaMda,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaMda))+gammahM[,,a,sigmaMda]
                                            +diag(2)*age_params$muh[a]) %*% ((fAh_mda[,sigmaMda,a,t]/dt)+(psi_f[a,t]*fAh[,a,t]/dsigmaMda) + 
                 eta_M[sigmaMda]*lambdam[,t]*fSh_mda[sigmaMda,a,t])
          } else {
            mSh_mda[sigmaMda,a,t+1]= ((mSh_mda[sigmaMda,a,t]/dt) +(mSh_mda[sigmaMda-1,a,t+1]/dsigmaMda))/
              ((1/dt)+(1/da)+(1/dsigmaMda)+age_params$muh[a] + 
                 eta_M[sigmaMda]*sum(lambdam[,t])+k_M[sigmaMda])
            
            fSh_mda[sigmaMda,a,t+1]= ((fSh_mda[sigmaMda,a,t]/dt) +(fSh_mda[sigmaMda-1,a,t+1]/dsigmaMda))/
              ((1/dt) + (1/da) + (1/dsigmaMda) + age_params$muh[a] + eta_M[sigmaMda] * sum(lambdam[,t]) + k_M[sigmaMda])
            
            mAh_mda[,sigmaMda,a,t+1]= solve(diag(2)*((1/dt) + (1/da) + (1/dsigmaMda)) +gammahM[,,a,sigmaMda] 
                                            + diag(2)*age_params$muh[a]) %*% ((mAh_mda[,sigmaMda,a,t]/dt)+(mAh_mda[,sigmaMda-1,a,t+1]/dsigmaMda)+ 
                 eta_M[sigmaMda]*lambdam[,t]*mSh_mda[sigmaMda,a,t])
            
            fAh_mda[,sigmaMda,a,t+1]= solve(diag(2)*((1/dt) + (1/da) + (1/dsigmaMda)) + gammahM[,,a,sigmaMda] + diag(2)*age_params$muh[a]) %*%
              ((fAh_mda[,sigmaMda,a,t]/dt)+(fAh_mda[,sigmaMda-1,a,t+1]/dsigmaMda)+eta_M[sigmaMda]*lambdam[,t]*fSh_mda[sigmaMda,a,t])
          }
        } else {
          if (sigmaMda == 1) {
            mSh_mda[sigmaMda,a,t+1]= ((mSh_mda[sigmaMda,a,t]/dt) + (mSh_mda[sigmaMda,a-1,t+1]/da) + 
                                        (psi_m[a,t] * (mSh[a,t] + sum(mRh[,a,t])) / dsigmaMda)) /
              ((1/dt) + (1/da) + (1/dsigmaMda) + age_params$muh[a] + eta_M[sigmaMda]*sum(lambdam[,t]) + k_M[sigmaMda])
            
            fSh_mda[sigmaMda,a,t+1]= ((fSh_mda[sigmaMda,a,t]/dt) + (fSh_mda[sigmaMda,a-1,t+1]/da) + 
                                        (psi_f[a,t] * (fSh[a,t] + sum(fRh[,a,t])) / dsigmaMda)) /
              ((1/dt) + (1/da) + (1/dsigmaMda) + age_params$muh[a] + eta_M[sigmaMda] * sum(lambdam[,t]) + k_M[sigmaMda])
            
            mAh_mda[,sigmaMda,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaMda)) + gammahM[,,a,sigmaMda] + 
                                              diag(2)*age_params$muh[a])%*%((mAh_mda[,sigmaMda,a,t]/dt)+(mAh_mda[,sigmaMda,a-1,t+1]/da)+ 
                 (psi_m[a,t]*mAh[,a,t]/dsigmaMda) + eta_M[sigmaMda]*lambdam[,t]*mSh_mda[sigmaMda,a,t])
            
            fAh_mda[,sigmaMda,a,t+1]= solve(diag(2) * ((1/dt) + (1/da) + (1/dsigmaMda)) + gammahM[,,a,sigmaMda] 
                                            +diag(2)*age_params$muh[a]) %*% ((fAh_mda[,sigmaMda,a,t]/dt)+(fAh_mda[,sigmaMda,a-1,t+1]/da)+ 
                                              (psi_f[a,t]*fAh[,a,t]/dsigmaMda) + eta_M[sigmaMda]*lambdam[,t]*fSh_mda[sigmaMda,a,t])
          } else {
            mSh_mda[sigmaMda,a,t+1]= ((mSh_mda[sigmaMda,a,t]/dt) +(mSh_mda[sigmaMda,a-1,t+1]/da) + 
                                        (mSh_mda[sigmaMda-1,a,t+1]/dsigmaMda)) /
              ((1/dt)+(1/da)+(1/dsigmaMda)+age_params$muh[a] + eta_M[sigmaMda]*sum(lambdam[,t]) + k_M[sigmaMda])
            
            fSh_mda[sigmaMda,a,t+1]= ((fSh_mda[sigmaMda,a,t]/dt) + (fSh_mda[sigmaMda,a-1,t+1]/da) + 
                                        (fSh_mda[sigmaMda-1,a,t+1]/dsigmaMda)) /
              ((1/dt)+(1/da)+(1/dsigmaMda) + age_params$muh[a] + eta_M[sigmaMda]*sum(lambdam[,t]) + k_M[sigmaMda])
            
            mAh_mda[,sigmaMda,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaMda)) + gammahM[,,a,sigmaMda] 
                                            + diag(2)*age_params$muh[a]) %*% ((mAh_mda[,sigmaMda,a,t]/dt) + (mAh_mda[,sigmaMda,a-1,t+1]/da) + 
                 (mAh_mda[,sigmaMda-1,a,t+1]/dsigmaMda) + eta_M[sigmaMda]*lambdam[,t]*mSh_mda[sigmaMda,a,t])
            
            fAh_mda[,sigmaMda,a,t+1]= solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaMda)) + gammahM[,,a,sigmaMda] + diag(2)*age_params$muh[a]) %*%
              ((fAh_mda[,sigmaMda,a,t]/dt) + (fAh_mda[,sigmaMda,a-1,t+1]/da) + 
                 (fAh_mda[,sigmaMda-1,a,t+1]/dsigmaMda) + eta_M[sigmaMda]*lambdam[,t]*fSh_mda[sigmaMda,a,t])
          }
        }
      }
    }
    
    Dh_increment= 0
    for (a in 1:Na) {Dh_increment= Dh_increment + sum(deltah[,,a]%*%(mIh[,a,t]+fIh[,a,t]))}
    Dh[t+1]= Dh[t] + dt * Dh_increment
    
    Sm[t+1]= (Sm[t]/dt + Lm) / (1/dt + mum + sum(lambdah[,t]))
    Im[,t+1]= solve(diag(2)/dt + deltam) %*% (Im[,t]/dt + lambdah[,t] * Sm[t])
    
    Nh[t+1]= sum(mSh[,t+1] + fSh[,t+1]) + sum(mAh[,,t+1] + fAh[,,t+1]) + sum(mIh[,,t+1] + fIh[,,t+1]) + 
      sum(mRh[,,t+1] + fRh[,,t+1]) +  sum(mSh_mda[,,t+1] + fSh_mda[,,t+1]) + sum(mSh_smc[,,t+1] + fSh_smc[,,t+1]) + 
      sum(mAh_mda[,,,t+1] + fAh_mda[,,,t+1]) + sum(mAh_smc[,,,t+1] + fAh_smc[,,,t+1])
    Nm[t+1]= Sm[t+1] + sum(Im[,t+1])
  }
  
  ShTot= AhTot= IhTot= RhTot= ResistantHumanCases= ImTot= numeric(Nt)
  TotPopSmc= TotPopMda= Totlambdah= numeric(Nt)
  
  for (t in 1:Nt) {
    ResistantHumanCases[t]= sum(mAh[2,,t] + fAh[2,,t] + mIh[2,,t] + mIh[4,,t] + fIh[2,,t] + fIh[4,,t]) + 
      sum(mAh_smc[2,,,t]) + sum(fAh_smc[2,,,t]) + sum(mAh_mda[2,,,t]) + sum(fAh_mda[2,,,t])
    
    ShTot[t]= sum(mSh[,t]+fSh[,t]) + sum(mSh_smc[,,t] + fSh_smc[,,t]) +sum(mSh_mda[,,t] + fSh_mda[,,t])
    AhTot[t]= sum(mAh[,,t]+fAh[,,t]) + sum(mAh_smc[,,,t]+fAh_smc[,,,t])+ sum(mAh_mda[,,,t] + fAh_mda[,,,t])
    IhTot[t]= sum(mIh[,,t] + fIh[,,t]); RhTot[t]= sum(mRh[,,t] + fRh[,,t])
    
    TotPopSmc[t]= sum(mAh_smc[,,,t] + fAh_smc[,,,t]) + sum(mSh_smc[,,t] + fSh_smc[,,t])
    TotPopMda[t]= sum(mAh_mda[,,,t] + fAh_mda[,,,t]) + sum(mSh_mda[,,t] + fSh_mda[,,t])
    
    Totlambdah[t]= sum(lambdah[,t]); ImTot[t]= sum(Im[,t])
  }
  
  Resis_Cases=sum(ResistantHumanCases[])
  PropSh= ShTot/Nh; PropAh= AhTot / Nh; PropIh= IhTot / Nh; Prophum_Inf= (IhTot+AhTot)/ Nh
  PropRh= RhTot / Nh; PropSm= Sm / Nm; PropIm= ImTot / Nm
  T_eq = ifelse(Nyr >= 1 & PropSMC == 0 & PropMDA == 0 & PropTreatment == 0, 
                t_begin_Camp / dt, Nt)
  
  return(list( "mSh" = mSh, "mAh" = mAh, "mIh" = mIh, "mRh" = mRh, "fSh" = fSh, "fAh" = fAh, "fIh" = fIh,
               "fRh" = fRh,"Sm" = Sm, "Im" = Im, "Dh" = Dh,"ShTot" = ShTot, "AhTot" = AhTot, "IhTot" = IhTot,
               "RhTot" = RhTot, "ImTot" = ImTot, "CumulDh" = Dh[Nt], "RhCases" = ResistantHumanCases, 
               "PrevRhCases" = ResistantHumanCases[Nt], "TotPopSmc" = TotPopSmc, "TotPopMda" = TotPopMda,
               "lambdam" = lambdam, "lambdah" = lambdah, "Totlambdah" = Totlambdah, "PropSh" = PropSh, 
               "PropAh" = PropAh, "PropIh" = PropIh,  "PropRh" = PropRh, "PropSm" = PropSm, "PropIm" = PropIm,
               "time" = time, "Nh" = Nh, "Nm" = Nm, "mSh_final" = mSh[, T_eq], "fSh_final" = fSh[, T_eq],
               "mAh_final" = mAh[,,T_eq], "fAh_final" = fAh[,,T_eq], "mIh_final" = mIh[,,T_eq], "fIh_final" = fIh[,,T_eq],
               "mRh_final" = mRh[,,T_eq], "fRh_final" = fRh[,,T_eq], "mSh_smc_final" = mSh_smc[,,T_eq], 
               "fSh_smc_final" = fSh_smc[,,T_eq], "mAh_smc_final" = mAh_smc[,,,T_eq], "fAh_smc_final" = fAh_smc[,,,T_eq],
               "mSh_mda_final" = mSh_mda[,,T_eq], "fSh_mda_final" = fSh_mda[,,T_eq], "mAh_mda_final" = mAh_mda[,,,T_eq],
               "fAh_mda_final" = fAh_mda[,,,T_eq], "Sm_final" = Sm[T_eq], "Im_final" = Im[,T_eq], 
               "Dh_final" = Dh[T_eq],"Prophum_Inf"=Prophum_Inf, "Resis_Cases"=Resis_Cases))
}

# Modelzero= No control
# ModelStrategy01= Treatement Alone
# ModelStrategy02= MDA+ standart SMC Alone
# ModelStrategy03= Treatment+ MDA+ standart SMC
# ModelStrategy04= MDA+ Expanded SMC Alone
# ModelStrategy05= Treatment+ MDA+ Expanded SMC

# Fonction de paramètres de base
params_base = function(t_begin_Camp = 120, epsi_s = 1e-1) {
  list(Gap = 1, dt = 1, p_f = 0.15, LengYr = 360, t_begin_Camp = t_begin_Camp, Time_plot0 = 0,
       epsi_s = epsi_s, MDA_Dur_cycle = 7, SMC_Dur_cycle = 7)
}

# Seasonality functions
{
  p_hbr.fun = function(time, hbr_high, hbr_low, T_high_Seas, t_begin_high_Seas, Gap, LengYr) {
    hbr = rep(hbr_low, length(time))
    
    for (i in 1:length(time)) {
      t_year = time[i] %% LengYr
      ta = t_begin_high_Seas; tb = ta + T_high_Seas #T_high_Seas Durée haute transmission
      
      if (t_year < ta - Gap) {
        hbr[i] = hbr_low
      } else if (t_year >= ta - Gap && t_year < ta) {
        hbr[i] = hbr_low + ((t_year - (ta - Gap)) / Gap) * (hbr_high - hbr_low)
      } else if (t_year >= ta && t_year < tb) {
        hbr[i] = hbr_high
      } else if (t_year >= tb && t_year <= tb + Gap) {
        hbr[i] = hbr_high - ((t_year - tb) / Gap) * (hbr_high - hbr_low)
      } else {
        hbr[i] = hbr_low
      }
    }
    return(hbr)
  }
  
  Model_hbr = function(Saisonality, Nyr, LengYr, Gap, dt) {
    hbr_high = 0.55; hbr_low = 0.35
    T_high_Seas = 120; t_begin_high_Seas = 120
    
    Tmax = Nyr * LengYr
    time = seq(0, Tmax, by = dt)
    Nt = length(time)
    
    if (Saisonality == 1) {
      hbr = p_hbr.fun(time, hbr_high, hbr_low, T_high_Seas, t_begin_high_Seas, Gap, LengYr)
    } else {
      hbr = rep(hbr_high, Nt)
    }
    return(hbr)
  }
}

calc_reduction = function(intervention, control, variable = "PropIh") {
  100 * (1 - intervention[[variable]] / (control[[variable]]))
}

# Main simulation functions
run_simulation_Eq = function(params, PropTreatment = 0, PropSMC = 0, PropMDA = 0, a_lim = 5,
                             age_cible=0, Number_of_cycle = 1, Saisonality = 1, Nyr = 1, 
                             bar_epss= 0.5,InitPrevR= 0.001,initial_conditions = NULL) {
  
  MDA_VectTime_between_cycles = if(Number_of_cycle > 1) rep(30, Number_of_cycle - 1) else c(0, 0, 0)
  SMC_VectTime_between_cycles = if(Number_of_cycle > 1) rep(30, Number_of_cycle - 1) else c(0, 0, 0)
  
  hbr = Model_hbr(Saisonality, Nyr, params$LengYr, params$Gap, params$dt)
  
  Mainfunction(params$MDA_Dur_cycle, MDA_VectTime_between_cycles, params$t_begin_Camp, 
               Number_of_cycle, SMC_VectTime_between_cycles, params$SMC_Dur_cycle, PropSMC,
               PropMDA, params$Gap, Nyr, params$LengYr, PropTreatment, params$dt, params$p_f,
               a_lim, age_cible, hbr, params$epsi_s, bar_epss, InitPrevR, initial_conditions)
}
#bar_epss= 0.5
run_simulation = function(params, PropTreatment = 0, PropSMC = 0, PropMDA = 0,  a_lim = 5, 
                          age_cible=0, Number_of_cycle = 1, Saisonality = 1,  Nyr = 1,
                          bar_epss= 0.5, InitPrevR= 0.001, initial_conditions = NULL) {
  
  MDA_VectTime_between_cycles = if(Number_of_cycle > 1) rep(30, Number_of_cycle - 1) else c(0, 0, 0)
  SMC_VectTime_between_cycles = if(Number_of_cycle > 1) rep(30, Number_of_cycle - 1) else c(0, 0, 0)
  
  hbr = Model_hbr(Saisonality, Nyr, params$LengYr, params$Gap, params$dt)
  
  Mainfunction(params$MDA_Dur_cycle, MDA_VectTime_between_cycles, params$t_begin_Camp, 
               Number_of_cycle, SMC_VectTime_between_cycles, params$SMC_Dur_cycle,  PropSMC, 
               PropMDA, params$Gap, Nyr, params$LengYr, PropTreatment,   params$dt, params$p_f, 
               a_lim, age_cible, hbr, params$epsi_s, bar_epss, InitPrevR, initial_conditions)
}

# Get equilibrium conditions
params_eq = params_base()
Equilibrium = run_simulation_Eq(params_eq, PropTreatment = 0, PropSMC = 0, PropMDA = 0)

initial_conditions_eq = list(
  mSh_final = Equilibrium$mSh_final, fSh_final = Equilibrium$fSh_final,
  mAh_final = Equilibrium$mAh_final, fAh_final = Equilibrium$fAh_final,
  mIh_final = Equilibrium$mIh_final, fIh_final = Equilibrium$fIh_final,
  mRh_final = Equilibrium$mRh_final, fRh_final = Equilibrium$fRh_final,
  mSh_smc_final = Equilibrium$mSh_smc_final, fSh_smc_final = Equilibrium$fSh_smc_final,
  mAh_smc_final = Equilibrium$mAh_smc_final, fAh_smc_final = Equilibrium$fAh_smc_final,
  mSh_mda_final = Equilibrium$mSh_mda_final, fSh_mda_final = Equilibrium$fSh_mda_final,
  mAh_mda_final = Equilibrium$mAh_mda_final, fAh_mda_final = Equilibrium$fAh_mda_final,
  Sm_final = Equilibrium$Sm_final, Im_final = Equilibrium$Im_final, Dh_final = Equilibrium$Dh_final
)

# Without control
{
  plot_results = function(test_result, params,Nyr = 1,  title_suffix = "") {
    time = test_result$time
    Tmax = Nyr * params$LengYr
    Tmaxmonths = (Tmax - params$Time_plot0) / 30
    
    pdf(paste0("Dynamic_Evolution_", title_suffix, ".pdf"), width = 7, height = 3)
    par(oma = c(1, .1, 1, .1), mar = c(2, 3.2, 1, .5), mfrow = c(1, 2))
    
    ColVect = c("#037153", "#8B4513", "#ff3355", "#205072")
    LineVect = c(1, 3, 2, 4)
    
    # Graphique populations humaines
    plot(-1, 1, type = "l", xlim = c(params$Time_plot0, Tmax), ylim = c(0, 0.42),
         xlab = "", ylab = "", yaxt = "n", xaxt = "n")
    axis(1, at = seq(params$Time_plot0, Tmax, by = 60), 
         labels = seq(0, Tmaxmonths, by = 2), cex.axis = 0.7)
    axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c(0, 10, 20, 30, 40), cex.axis = 0.8)
    
    # Tracer les courbes
    populations = list(test_result$ShTot, test_result$AhTot, test_result$IhTot, test_result$RhTot)
    labels = c(expression(S[h]), expression(A[h]), expression(I[h]), expression(R[h]))
    
    for(i in 1:4) {
      lines(time, populations[[i]]/test_result$Nh, lwd = 2, lty = LineVect[i], col = ColVect[i])
    }
    
    mtext("Human population (%)", side = 2, cex = 0.85, line = 2)
    mtext("Time (Months)", side = 1, cex = 0.9, line = 2)
    legend("bottomleft", legend = labels, lwd = 2, col = ColVect, lty = LineVect, 
           cex = 0.75, bty = "n", horiz = TRUE)
    text(params$Time_plot0 - 0.05, 0.42 * 1.15, "(A)", cex = 1, xpd = NA)
    
    # Graphique populations moustiques
    plot(-1, 1, type = "l", xlim = c(params$Time_plot0, Tmax), ylim = c(0, 0.82),
         xlab = "", ylab = "", yaxt = "n", xaxt = "n")
    axis(1, at = seq(params$Time_plot0, Tmax, by = 60), 
         labels = seq(0, Tmaxmonths, by = 2), cex.axis = 0.7)
    axis(2, at = seq(0, 0.75, 0.15), labels = seq(0, 75, 15), cex.axis = 0.8)
    
    lines(time, test_result$Sm/test_result$Nm, lwd = 2, lty = 1, col = "#037153")
    lines(time, test_result$ImTot/test_result$Nm, lwd = 2, lty = 2, col = "#ff3355")
    
    mtext("Mosquito population (%)", side = 2, cex = 0.85, line = 2)
    mtext("Time (Months)", side = 1, cex = 0.9, line = 2)
    legend("bottomleft", legend = c(expression(S[m]), expression(I[m])), 
           lwd = 2, col = c("#037153", "#ff3355"), lty = c(1, 2), cex = 0.75, bty = "n", horiz = TRUE)
    text(params$Time_plot0 - 0.05, 0.82 * 1.15, "(B)", cex = 1, xpd = NA)
    
    dev.off()
  }
  
  params=params_base()
  baseline = run_simulation(params,PropTreatment = 0, PropSMC = 0, 
                            PropMDA = 0, initial_conditions = initial_conditions_eq)
  plot_results(baseline, params, Nyr = 1, "baseline_No_Control")
}

# Scénario de visualisation test pour le modèle MDA/SMC 
params = params_base()
params$Nyr = 2  # 2 ans
params$p_f = 0.15  
{
  # Scénario 1: Contrôle complet (SMC + MDA)
  test_combined = run_simulation(params, PropTreatment = 0, PropSMC = 0.7,PropMDA = 0.7, a_lim = 5, 
                                 Number_of_cycle = 4, Saisonality = 1, Nyr = params$Nyr, 
                                 initial_conditions = initial_conditions_eq)
  
  # Scénario 2: SMC seul
  test_SMC_only = run_simulation(params, PropTreatment = 0, PropSMC = 0.7,PropMDA = 0, a_lim = 5,
                                 Number_of_cycle = 4, Saisonality = 1, Nyr = params$Nyr, 
                                 initial_conditions = initial_conditions_eq)
  
  # Scénario 3: MDA seul
  test_MDA_only = run_simulation(params, PropTreatment = 0, PropSMC = 0, PropMDA = 0.7, a_lim = 5,
                                 Number_of_cycle = 4, Saisonality = 1, Nyr = params$Nyr, 
                                 initial_conditions = initial_conditions_eq)
  
  # Scénario 4: Sans intervention (baseline)
  baseline_nocontrol = run_simulation(params, PropTreatment = 0, PropSMC = 0, PropMDA = 0, 
                                      Nyr = params$Nyr, initial_conditions = initial_conditions_eq)
  
  # Save results
  save(test_combined, test_SMC_only, test_MDA_only, baseline_nocontrol,
       file = "Scenario_Results_Dynamic_test_2years.RData")
  
  load("Scenario_Results_Dynamic_test_2years.RData")
  ls()
  # ========================================
  # PDF 1: PK/PD profile et Paramètres épidémiologiques
  # ========================================
  
  pdf("PK_PD_Epidemiological_Profile.pdf", width = 8, height = 3.2)
  layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE),
         widths = c(1,1), heights = c(1))
  par(oma = c(2, 1, 1, 0.5), mar = c(2, 2.2, 1, 1.5))
  LC = 2.5
  
  # FIGURE 1A: Paramètres épidémiologiques par âge
  {
    Va = seq(0, 90, by = 1)
    Na = length(Va)
    alpha_1 = 0.045; alpha_2 = 0.181
    G = function(a) 22.7 * a * exp(-0.0934 * a)
    betah_s = alpha_1 * (G(Va)^alpha_2)
    
    # Progression vers cas clinique par âge (valeurs réelles)
    nu_0 = 1/10; nu_1 = 1/150; nu_2 = 1/5
    nuh = rep(0, Na)
    for (i in 1:Na) {
      if (Va[i] <= 5) { nuh[i] = nu_2 }
      else if (Va[i] <= 15) { nuh[i] = nu_0 }
      else { nuh[i] = nu_1 }
    }
    
    plot(Va, betah_s, type = 'l', lwd = LC, xlim = c(0, 90), ylim = c(0, max(betah_s) * 1.1),
         xlab = "", ylab = "", axes = FALSE, frame.plot = TRUE)
    # SMC zone (0-5 ans)
    rect(0, 0, 5, max(betah_s) * 1.1, col = adjustcolor("#037153", alpha = 0.1), border = NA)
    
    # MDA zones pour tous (>5 ans)
    rect(5, 0, 90, max(betah_s) * 1.1, col = adjustcolor("#8B4513", alpha = 0.1), border = NA)
    
    # Zone d'exclusion pour les femmes en âge de procréer (15-45 ans) - hachuré
    rect(15, 0, 45, max(betah_s) * 0.55, col = adjustcolor("red", alpha = 0.1), border = NA)
    # Ajouter des hachures pour marquer l'exclusion
    for(i in seq(15, 45, by = 2)) {
      segments(i, 0, i, max(betah_s) * 0.55, col = adjustcolor("red", alpha = 0.3), lwd = 0.5)
    }
    
    # Lignes des paramètres (valeurs réelles)
    lines(Va, betah_s, lwd = LC, col = "darkblue", lty = 1)
    
    # Lignes verticales pour marquer les transitions d'âge importantes
    abline(v = 5, lty = 3, col = "gray50", lwd = 1)
    abline(v = 15, lty = 3, col = "gray50", lwd = 1)
    abline(v = 45, lty = 3, col = "gray50", lwd = 1)
    
    # Axes
    axis(1, at = seq(0, 90, by = 10), cex.axis = 0.7)
    axis(2, at = seq(0, max(betah_s), by = 0.02), 
         labels =seq(0, max(betah_s), by = 0.02),  cex.axis = 0.7, las = 1)
    
    mtext("Age−structured human population (years)", side = 1, line = 2.5, cex = 0.75)
    mtext(expression("Infection rate " * beta[h]), side = 2, line = 2, cex = 0.85)
    
    legend("topright", 
           legend = c(expression("Infection rate (" * beta[h] * ")"), 
                      "SMC target", 
                      "MDA target",
                      "MDA excluded females"),
           col = c("darkblue",  adjustcolor("#037153", alpha = 0.35),
                   adjustcolor("#8B4513", alpha = 0.35), adjustcolor( "red", alpha = 0.35)),
           lty = c(1,  NA, NA, NA), lwd = c(LC,  NA, NA, NA),
           pch = c(NA,  15, 15, 15), pt.cex = 1.5,
           bg = "white", box.lwd = 0.5, cex = 0.5)
    
    # Annotations des zones d'âge
    text(2.5, max(betah_s) * 1.04, "SMC\n0-5y", cex = 0.5, col = "#037153", font = 2)
    text(10, max(betah_s) * 0.65, "MDA\n5-15y", cex = 0.5, col = "#8B4513", font = 2)
    text(30, max(betah_s) * 0.35, "Females excluded", cex = 0.5, col = "red", font = 2)
    text(30, max(betah_s) * 0.65, "MDA males\n15-45y", cex = 0.5, col = "#8B4513", font = 2)
    text(67.5, max(betah_s) * 0.65, "MDA\n>45y", cex = 0.5, col = "#8B4513", font = 2)
    
    mtext("(A)", side = 3, adj = 0, cex = 0.8, font = 1, line = 0.5)
    mtext("Infection rate and intervention targets", side = 3, adj = 0.5, cex = 0.75, font = 1, line = 1)
  }
  
  
  # FIGURE 1B: Efficacité de MDA et SMC (PK/PD)
  {
    sigma_range = seq(0, 70, by = 0.5)
    T50_MDA = 35; shape_MDA = 8
    T50_SMC = 28; shape_SMC = 6.5
    
    Eff_MDA = 0.95 / (1 + (sigma_range / T50_MDA)^shape_MDA)
    Eff_SMC = 0.95 / (1 + (sigma_range / T50_SMC)^shape_SMC)
    
    plot(sigma_range, Eff_MDA, type = 'l', lwd = LC, col = "#8B4513", lty = 1, xlim = c(0, 64), 
         ylim = c(0, 1.05), xlab = "", ylab = "", axes = FALSE, frame.plot = TRUE)
    lines(sigma_range, Eff_SMC, lwd = LC, col = "#037153", lty = 2)
    
    axis(1, at = seq(0, 64, by = 7), labels = seq(0, 64, by = 7), cex.axis = 0.7)
    axis(2, at = seq(0, 1, by = 0.2), labels = seq(0, 100, by = 20), cex.axis = 0.7, las = 1)
    
    abline(v = T50_MDA, lty = 3, col = "#8B4513", lwd = 2)
    abline(v = T50_SMC, lty = 3, col = "#037153", lwd = 2)
    
    # Zone d'efficacité optimale avec gradient
    n_grad = 50
    y_seq = seq(0.5, 1, length.out = n_grad)
    for(i in 1:(n_grad-1)) {
      polygon(c(0, 64, 64, 0),
              c(y_seq[i], y_seq[i], y_seq[i+1], y_seq[i+1]),
              col = adjustcolor("#D3D3D3", alpha.f = 0.4 - 0.4*i/n_grad), 
              border = NA)
    }
    
    mtext("Time since administration (days)", side = 1, line = 2.5, cex = 0.85)
    mtext("Loss of drug efficacy (%)", side = 2, line = 2, cex = 0.85)
    
    legend("topright", legend = c("MDA", "SMC", expression(sigma[50]~"values"), "Optimal zone (>50%)"), 
           col = c("#8B4513", "#037153", "gray40", "#D3D3D3"), 
           lty = c(1, 2, 3, NA), lwd = c(LC, LC, 2, NA),
           pch = c(NA, NA, NA, 15), pt.cex = 1, bg = "white", box.lwd = 0.5, cex = 0.4)
    
    text(T50_MDA-1, 0.07, expression(sigma[50]^MDA), cex = 0.7, col = "#8B4513", font = 2)
    text(T50_SMC-1, 0.07, expression(sigma[50]^SMC), cex = 0.7, col = "#037153", font = 2)
    
    mtext("(B)", side = 3, adj = 0, cex = 0.8, font = 1, line = 0.5)
    mtext("PK/PD profile", side = 3, adj = 0.5, cex = 0.7, font = 1, line = 1)
  }
  
  dev.off()
  
  
  # ========================================
  # PDF 2: Scénarios principaux sur 2 ans (1 ligne x 3 colonnes)
  # ========================================
  
  pdf("Scenario_MDA_SMC_test.pdf", width = 11, height = 3.4)
  layout(matrix(c(1,2,3,4,4,4), nrow = 2, ncol = 3, byrow = TRUE),
         widths = c(1,1,1), heights = c(1, 0.15))
  par(oma = c(0.5, 1, 1, 1), mar = c(3.5, 3, 2, 2))
  LC = 2.5
  
  # FIGURE 2A: Population sous protection
  {
    time_plot = test_combined$time
    pop_smc = test_combined$TotPopSmc / test_combined$Nh * 100
    pop_mda = test_combined$TotPopMda / test_combined$Nh * 100
    pop_total_protected = pop_smc + pop_mda
    ymax = max(pop_total_protected) * 1.1
    
    plot(time_plot, pop_total_protected, type = 'n', xlim = c(0, 720), ylim = c(0, ymax),
         xlab = "", ylab = "", axes = FALSE, frame.plot = TRUE)
    
    # Zone d'intervention (de t_begin à Tmax)
    polygon(c(120, 720, 720, 120), c(0, 0, ymax, ymax),
            col = adjustcolor("#D3D3D3", alpha = 0.15), border = NA)
    
    # Zones remplies pour visualiser la couverture
    polygon(c(time_plot, rev(time_plot)), c(pop_smc, rep(0, length(time_plot))),
            col = rgb(3/255, 113/255, 83/255, alpha = 0.3), border = NA)
    polygon(c(time_plot, rev(time_plot)), c(pop_total_protected, rev(pop_smc)),
            col = rgb(139/255, 69/255, 19/255, alpha = 0.3), border = NA)
    
    # Marqueurs des interventions pour les deux années
    cycle_times = c(120)
    for (ct in cycle_times) {abline(v = ct, lty = 2, lwd = 1.5, col = "#205072") }
    
    # Lignes principales
    lines(time_plot, pop_total_protected, lwd = LC, col = "purple")
    lines(time_plot, pop_smc, lwd = LC, col = "#037153", lty = 4)
    lines(time_plot, pop_mda, lwd = LC, col = "#8B4513", lty = 2)
    
    axis(1, at = seq(0, 720, by = 120), labels = seq(0, 24, by = 4), cex.axis = 0.8)
    axis(2, at = seq(0, ceiling(ymax), by = 10), cex.axis = 0.7, las = 1)
    
    mtext("Time (months)", side = 1, line = 2.5, cex = 0.75)
    mtext("Protected human population (%)", side = 2, line = 2, cex = 0.75)
    
    mtext("(A)", side = 3, adj = 0, cex = 0.9, font = 1, line = 0.5)
    mtext("Coverage achieved by interventions", side = 3, adj = 0.5, cex = 0.8, font = 1, line = 0.5)
    # # Légende
    # legend("topright",  legend = c("MDA", "Standart SMC","Total protected" ),
    #        col = c("#8B4513","#037153","purple"), lty = c(2, 4, 1), lwd = c(LC, LC, LC),
    #        bg = "white", box.lwd = 0.5, cex = 0.7)
  }
  
  # FIGURE 2B: Force d'infection réduite 
  {
    time_plot = baseline_nocontrol$time
    # Forces d'infection pour chaque scénario
    lambda_h_baseline = baseline_nocontrol$Totlambdah
    lambda_h_SMC = test_SMC_only$Totlambdah
    lambda_h_MDA = test_MDA_only$Totlambdah
    lambda_h_combined = test_combined$Totlambdah
    
    # Réduction relative pour chaque intervention
    reduction_SMC = 100 * (1 - lambda_h_SMC / lambda_h_baseline)
    reduction_MDA = 100 * (1 - lambda_h_MDA / lambda_h_baseline)
    reduction_combined = 100 * (1 - lambda_h_combined / lambda_h_baseline)
    
    # Remplacer NaN et Inf par 0
    reduction_SMC[!is.finite(reduction_SMC)] = 0
    reduction_MDA[!is.finite(reduction_MDA)] = 0
    reduction_combined[!is.finite(reduction_combined)] = 0
    
    # Y max adaptatif
    ymax_reduction = max(c(reduction_SMC, reduction_MDA, reduction_combined), na.rm = TRUE) * 1.01
    
    plot(time_plot, reduction_combined, type = 'n',
         xlim = c(0, 720), ylim = c(0, ymax_reduction),
         xlab = "", ylab = "", axes = FALSE, frame.plot = TRUE)
    
    # Zone d'intervention (de t_begin à Tmax)
    polygon(c(120, 720, 720, 120),
            c(0, 0, ymax_reduction, ymax_reduction),
            col = adjustcolor("#D3D3D3", alpha = 0.15), border = NA)
    
    Time_plot_max=length(time_plot)-1
    lines(time_plot[1:Time_plot_max], reduction_SMC[1:Time_plot_max], lwd = LC, col = "#037153", lty = 4)
    lines(time_plot[1:Time_plot_max], reduction_MDA[1:Time_plot_max], lwd = LC, col = "#8B4513", lty = 2)
    lines(time_plot[1:Time_plot_max], reduction_combined[1:Time_plot_max], lwd = LC, col = "purple", lty = 1)
    
    # Marqueurs des interventions
    cycle_times = c(120)
    for (ct in cycle_times) {abline(v = ct, lty = 2, col = "#205072", lwd = 1.5)}
    
    axis(1, at = seq(0, 720, by = 120), labels = seq(0, 24, by = 4), cex.axis = 0.8)
    axis(2, at = seq(0, ceiling(ymax_reduction/10)*10, by = 10), cex.axis = 0.7, las = 1)
    
    # abline(h = 50, lty = 3, col = "gray50")
    
    mtext("Time (months)", side = 1, line = 2.5, cex = 0.75)
    mtext(expression("Reduction in force of infection " * lambda[h](t) * " (%)"), side = 2, line = 2.5, cex = 0.75)
    
    mtext("(B)", side = 3, adj = 0, cex = 0.9, font = 1, line = 0.5)
    mtext("Transmission reduction", side = 3, adj = 0.5, cex = 0.8, font = 1, line = 0.5)
  }
  
  # FIGURE 2C: Cas cliniques évités ET Proportion de cas résistants (double axe Y)
  {
    time_plot = baseline_nocontrol$time
    
    # Calcul des cas cliniques évités (en %)
    cases_baseline = baseline_nocontrol$IhTot
    cases_SMC = test_SMC_only$IhTot
    cases_MDA = test_MDA_only$IhTot
    cases_combined = test_combined$IhTot
    
    averted_SMC = 100 * (1 - cases_SMC / cases_baseline)
    averted_MDA = 100 * (1 - cases_MDA / cases_baseline)
    averted_combined = 100 * (1 - cases_combined / cases_baseline)
    
    # Proportion de cas résistants
    resist_SMC = 100 * test_SMC_only$RhCases / (test_SMC_only$IhTot )
    resist_MDA = 100 * test_MDA_only$RhCases / (test_MDA_only$IhTot )
    resist_combined = 100 * test_combined$RhCases / (test_combined$IhTot )
    resist_baseline = 100 * baseline_nocontrol$RhCases / (baseline_nocontrol$IhTot)
    
    # Remplacer NaN et Inf par 0
    averted_SMC[!is.finite(averted_SMC)] = 0
    averted_MDA[!is.finite(averted_MDA)] = 0
    averted_combined[!is.finite(averted_combined)] = 0
    resist_SMC[!is.finite(resist_SMC)] = 0
    resist_MDA[!is.finite(resist_MDA)] = 0
    resist_combined[!is.finite(resist_combined)] = 0
    resist_baseline[!is.finite(resist_baseline)] = 0
    
    # Y max adaptatifs
    ymax_averted = max(c(averted_SMC, averted_MDA, averted_combined), na.rm = TRUE) * 1.1
    ymax_resist = 5
    
    # Plot principal pour cas évités
    plot(time_plot, averted_combined, type = 'n',xlim = c(0, 720), ylim = c(0, ymax_averted),
         xlab = "", ylab = "", axes = FALSE, frame.plot = TRUE)
    axis(2, at = seq(0, ceiling(ymax_averted/10)*10, by = 10), cex.axis = 0.7, las = 1)
    
    # Zone d'intervention
    polygon(c(120, 720, 720, 120),
            c(0, 0, ymax_averted, ymax_averted),
            col = adjustcolor("#D3D3D3", alpha = 0.15), border = NA)
    
    # Lignes pour cas évités (traits pleins avec différents styles)
    lines(time_plot, averted_SMC, lwd = LC, col = "#037153", lty = 4)
    lines(time_plot, averted_MDA, lwd = LC, col = "#8B4513", lty = 2)
    lines(time_plot, averted_combined, lwd = LC, col = "purple", lty = 1)
    
    # Créer un second plot pour la résistance (axe Y4)
    par(new = TRUE)
    plot(time_plot, resist_combined, type = 'n',
         xlim = c(0, 720), ylim = c(0, 5),
         xlab = "", ylab = "", axes = FALSE)
    
    # # Lignes pour résistance (toutes en pointillés, lty=3)
    # lines(time_plot, resist_baseline, 
    #       lwd = 1.5, col = adjustcolor("black", alpha = 0.7), lty = 3)
    lines(time_plot, resist_SMC, 
          lwd =  LC, col = adjustcolor("#037153", alpha = 0.8), lty = 3)
    lines(time_plot, resist_MDA, 
          lwd =  LC, col = adjustcolor("#8B4513", alpha = 0.8), lty = 3)
    lines(time_plot, resist_combined, 
          lwd =  LC, col = adjustcolor("purple", alpha = 0.8), lty = 3)
    
    # Marqueurs des interventions
    cycle_times = c(120)
    for (ct in cycle_times) {
      abline(v = ct, lty = 2, col = "#205072", lwd = 1.5)
      # text(ct, ymax_averted * 0.95, paste0("Y", which(cycle_times == ct)), cex = 0.75, col = "#205072")
    }
    
    # Axes
    axis(1, at = seq(0, 720, by = 120), labels = seq(0, 24, by = 4), cex.axis = 0.8)
    axis(4, at = seq(0, ceiling(ymax_resist/10)*10, by = 1), cex.axis = 0.7, las = 1)
    # # Ligne de référence 50%
    # abline(h = 50, lty = 3, col = "gray50", lwd = 0.5)
    # 
    mtext("Time (months)", side = 1, line = 2, cex = 0.75)
    mtext("Clinical cases averted (%)", side = 2, line = 2, cex = 0.75)
    mtext("Clinical resistant cases (%)", side = 4, line = 2, cex = 0.75)
    
    mtext("(C)", side = 3, adj = 0, cex = 0.9, font = 1, line = 0.5)
    mtext("Clinical impact & resistance", side = 3, adj = 0.5, cex = 0.8, font = 1, line = 0.5)
  }
  
  # Légende commune horizontale en bas
  {
    reset = function() {
      par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
      plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
    }
    reset()
    
    legend("bottom", legend = c( "MDA alone", "Standart SMC alone", "MDA+standart SMC ", 
                                 "Resistance (%)"),
           xpd = NA, horiz = TRUE, inset = c(0, 0.01), bty = "n", lty = c( 2, 4, 1, 3), 
           lwd = c( 2, 2, 2, 2),  col = c( "#8B4513", "#037153", "purple", "black"), 
           cex = 0.8)
  }
  
  dev.off()
}

# Scenarios clinical cases with seasonality 
Vect_t_begin_Camp =c(90, 120, 150)
VectProp =c(0.7) #c(0.7, 0.9)
VectNumber_of_cycle =c(3, 4, 5)
CovTreatment = 0.7
Nyr = 2
{
  # Initialize result storage
  ModelZero = list()
  ModelStrategy01 = list() 
  ModelStrategy02 = list()
  ModelStrategy03 = list()
  # ModelStrategy04 = list()
  # ModelStrategy05 = list()
  
  # Run baseline scenarios for each t_begin_Camp
  for (p in 1:length(Vect_t_begin_Camp)) {
    params_p = params_base(t_begin_Camp = Vect_t_begin_Camp[p])
    ModelZero[[p]] = run_simulation(params_p, PropTreatment = 0, PropSMC = 0, 
                                    PropMDA = 0, Nyr = Nyr, 
                                    initial_conditions = initial_conditions_eq)
  }
  
  # Run intervention scenarios
  for (p in 1:length(Vect_t_begin_Camp)) {
    params_p = params_base(t_begin_Camp = Vect_t_begin_Camp[p])
    
    ModelStrategy01[[p]] = run_simulation(params_p, PropTreatment = CovTreatment,  PropSMC = 0, 
                                          PropMDA = 0, Nyr = Nyr, initial_conditions = initial_conditions_eq)
    
    ModelStrategy02[[p]] = list()
    ModelStrategy03[[p]] = list()
    # ModelStrategy04[[p]] = list()
    # ModelStrategy05[[p]] = list()
    
    for (i in 1:length(VectProp)) {
      ModelStrategy02[[p]][[i]] = list()
      ModelStrategy03[[p]][[i]] = list()
      # ModelStrategy04[[p]][[i]] = list()
      # ModelStrategy05[[p]][[i]] = list()
      
      for (k in 1:length(VectNumber_of_cycle)) {
        ModelStrategy02[[p]][[i]][[k]] = run_simulation( params_p, PropTreatment = 0,PropSMC = VectProp[i], a_lim = 5,
                                                         PropMDA = VectProp[i], Number_of_cycle = VectNumber_of_cycle[k], Nyr = Nyr,
                                                         initial_conditions = initial_conditions_eq)
        
        ModelStrategy03[[p]][[i]][[k]] = run_simulation( params_p, PropTreatment = CovTreatment, PropSMC = VectProp[i], a_lim = 5,
                                                         PropMDA = VectProp[i], Number_of_cycle = VectNumber_of_cycle[k], Nyr = Nyr,
                                                         initial_conditions = initial_conditions_eq )
        # ModelStrategy04[[p]][[i]][[k]] = run_simulation( params_p, PropTreatment = 0,PropSMC = VectProp[i], a_lim = 10,
        #                                                  PropMDA = VectProp[i], Number_of_cycle = VectNumber_of_cycle[k], Nyr = Nyr,
        #                                                  initial_conditions = initial_conditions_eq)
        # 
        # ModelStrategy05[[p]][[i]][[k]] = run_simulation( params_p, PropTreatment = CovTreatment, PropSMC = VectProp[i], a_lim = 10,
        #                                                  PropMDA = VectProp[i], Number_of_cycle = VectNumber_of_cycle[k], Nyr = Nyr,
        #                                                  initial_conditions = initial_conditions_eq )
      }
    }
  }
  
  # Save results
  save(ModelZero, ModelStrategy01, ModelStrategy02, ModelStrategy03, # ModelStrategy04,  ModelStrategy05,
       Vect_t_begin_Camp, VectProp, VectNumber_of_cycle, CovTreatment,
       params_base, Nyr, file = "Scenario_Results_Ih.RData")
  
  # PLOTTING FUNCTION
  plot_scenario_results = function() {
    load("Scenario_Results_Ih.RData")
    
    pdf("Scenario_cases_adverted30.pdf", width = 10.5, height = 8.2)
    par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.7, 2, .1))
    par(mfrow = c(length(Vect_t_begin_Camp), length(VectNumber_of_cycle)))
    
    ColVect = c("#ff3355","#037153", "purple", "#D98E04", "#8B4513", "#8B008B", "#205072")
    LineVect = c(3, 4, 1, 4, 2)
    LC = 2.5
    
    t_begin_high_Seas = 120; T_high_Seas = 120
    GrNumber = 0; i_prop = 1  # Using VectProp[1] = 0.7
    
    params_ref = params_base()
    Tmax = Nyr * params_ref$LengYr
    Tmaxmonths = (Tmax - params_ref$Time_plot0) / 30
    
    for (p in 1:length(Vect_t_begin_Camp)) {
      for (k in 1:length(VectNumber_of_cycle)) {
        GrNumber = GrNumber + 1
        time = ModelZero[[p]]$time
        
        # Calculer les réductions
        red_treatment = calc_reduction(ModelStrategy01[[p]], ModelZero[[p]], "PropIh")
        red_MDA_SMC_Standard = calc_reduction(ModelStrategy02[[p]][[i_prop]][[k]], ModelZero[[p]], "PropIh")
        red_combined_Standart = calc_reduction(ModelStrategy03[[p]][[i_prop]][[k]], ModelZero[[p]], "PropIh")
        # red_MDA_SMC_Extended = calc_reduction(ModelStrategy04[[p]][[i_prop]][[k]], ModelZero[[p]], "PropIh")
        # red_combined_Extended = calc_reduction(ModelStrategy05[[p]][[i_prop]][[k]], ModelZero[[p]], "PropIh")
        # 
        y_max=70
        # Plot
        plot(-1, 1, type = "l", xlab = "", xlim = c(params_ref$Time_plot0, Tmax),  ylab = "", 
             ylim = c(0, y_max), cex.lab = 1.2, yaxt = "n",xaxt = "n")
        
        # # Zone gradient
        # n_grad = 50
        # y_seq = seq(0, y_max, length.out = n_grad)
        # for(i in 1:(n_grad-1)) {
        #   polygon(c(t_begin_high_Seas, Tmax, Tmax, t_begin_high_Seas), c(y_seq[i], y_seq[i], y_seq[i+1], y_seq[i+1]),
        #           col = adjustcolor("#D3D3D3", alpha.f = 0.4 - 0.4*i/n_grad),  border = NA)
        # }
        # 
        # t_begin_Camp = 120; Tmax = Nyr * 360
        polygon(c(Vect_t_begin_Camp[p], Tmax, Tmax, Vect_t_begin_Camp[p]), 
                c(0, 0, y_max, y_max),
                col = adjustcolor("#D3D3D3", alpha.f = 0.1), border = NA)
        
        lines(time, red_treatment, lwd = LC, col = ColVect[1], lty = LineVect[1])
        lines(time, red_MDA_SMC_Standard, lwd = LC, col = ColVect[2], lty = LineVect[2])
        lines(time, red_combined_Standart, lwd = LC, col = ColVect[3], lty = LineVect[3])
        # lines(time, red_MDA_SMC_Extended, lwd = LC, col = ColVect[4], lty = LineVect[4])
        # lines(time, red_combined_Extended, lwd = LC, col = ColVect[5], lty = LineVect[5])
        
        # abline(h = 20, lty = 3, col = "gray50")
        # Ligne verticale début intervention
        abline(v = Vect_t_begin_Camp[p], lty = 2, col = "#205072", lwd = 1.5)
        text(Vect_t_begin_Camp[p], y_max* 0.9755, "Start", cex = 0.85, col = "#205072")
        
        axis(1, at = seq(params_ref$Time_plot0, Tmax, by = 180), 
             labels = seq(0, Tmaxmonths, by = 6), las = 1)
        axis(2, at = seq(0, y_max, by = 10), labels = seq(0, y_max, by = 10), cex.axis = 0.9)
        
        par(xpd = NA)
        text(params_ref$Time_plot0 - 0.05, y_max * 1.1, paste("(", LETTERS[GrNumber], ")",
                                                              sep = ""), cex = 1.6)
        par(xpd = FALSE)
        
        Peak_season = ifelse(p == 1, "MDA and SMC starts 1 month before peak season",
                             ifelse(p == 2, "MDA and SMC starts at peak season",
                                    "MDA and SMC starts 1 month after peak season"))
        
        if (GrNumber %in% c(1, 4, 7)) {
          mtext(Peak_season, side = 2, adj = 0.5, cex = 0.7, line = 3, font = 1)
          mtext("Clinical cases averted (%)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8)
        }
        
        if (GrNumber %in% c(1, 2, 3)) {
          mtext(bquote("#cycle for MDA and SMC: " ~ n[c] == .(VectNumber_of_cycle[k])),
                side = 3, cex = 0.95, line = 1)
        }
        
        if (GrNumber %in% c(7, 8, 9)) {
          mtext("Time (Months)", side = 1, adj = 0.5, cex = 0.9, line = 3, font = 0.8)
        }
      }
    }
    
    LegendTex = c( "Treatment", "MDA+Standart SMC", "Treatment+MDA+Standart SMC")
    # ,"MDA+Extended SMC", "Treatment+MDA+Extended SMC" )
    
    reset()
    legend("bottom", legend = LegendTex, xpd = NA, horiz = TRUE, inset = c(0, 0.01),
           bty = "n", lty = LineVect, lwd = 2, col =ColVect, cex = 0.85)
    
    dev.off()
  }
  
  # Run plotting
  plot_scenario_results()
}

#Scenario resistance et cases averted
{ 
  # SCENARIO PARAMETERS
  VecEsp = c(1/10, 1/100)
  VectProp = c(0.1, 0.5, 0.9)
  VectNumber_of_cycle =c(4,3,5) #c(3, 4, 5)
  VectCovTreatment = c(0.1, 0.5, 0.9)
  t_begin_Camp = 120
  Nyr=2
  
  # RUN SIMULATIONS AND SAVE
  for (idx_epsi in 1:length(VecEsp)) {
    epsi_s = VecEsp[idx_epsi]
    params_p = params_base(t_begin_Camp = t_begin_Camp, epsi_s = epsi_s)
    
    #Model zero resistance
    ModelZero_R = run_simulation(params_p, PropTreatment = 0, PropSMC = 0, 
                                 PropMDA = 0, Nyr = Nyr, initial_conditions = initial_conditions_eq)
    
    ModelStrategy03 = list()
    
    for (i in 1:length(VectProp)) {
      ModelStrategy03[[i]] = list()
      
      for (k in 1:length(VectNumber_of_cycle)) {
        ModelStrategy03[[i]][[k]] = list()
        
        for (j in 1:length(VectCovTreatment)) {
          ModelStrategy03[[i]][[k]][[j]] = run_simulation(params_p, PropTreatment = VectCovTreatment[j],
                                                          PropSMC = VectProp[i],  PropMDA = VectProp[i],  
                                                          Number_of_cycle = VectNumber_of_cycle[k], 
                                                          Nyr = Nyr,  initial_conditions = initial_conditions_eq )
        }
      }
    }
    
    filename = ifelse(idx_epsi == 1, 
                      "Scenario_Epsi_Fort.RData", 
                      "Scenario_Epsi_Faible.RData")
    
    save(ModelZero_R, ModelStrategy03, VecEsp, VectProp, VectNumber_of_cycle, 
         VectCovTreatment, t_begin_Camp, Nyr, epsi_s,  file = filename)
  }
  
  #Plot Functions
  plot_epsi_scenario = function(Nyr=2) {
    Tmax = Nyr * 365
    Tmaxmonths = Tmax/30
    
    # Load both data files
    load("Scenario_Epsi_Fort.RData")
    ModelStrategy03_Fort = ModelStrategy03
    ModelZero_R_Fort = ModelZero_R
    
    load("Scenario_Epsi_Faible.RData")
    ModelStrategy03_Faible = ModelStrategy03
    ModelZero_R_Faible = ModelZero_R
    
    ColVect = c("#8B4513", "#037153")
    VectPch = c(8, 16)
    LC = 2
    LineVect = c(1, 2)
    
    RowText = c("Coverage MDA and SMC: 10%", "Coverage MDA and SMC: 50%", "Coverage MDA and SMC: 90%")
    ColText = c("Treatment: 10%", "Treatment: 50%", "Treatment: 90%")
    
    # Create one figure for each VectNumber_of_cycle
    for (k in 1:length(VectNumber_of_cycle)) {
      
      pdf(paste0("Dyn_Cases_advert_Resist_Cycle30", VectNumber_of_cycle[k], ".pdf"), width = 10.5, height = 8.2)
      par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.7, 2, 3))
      par(mfrow = c(length(VectProp), length(VectCovTreatment)))
      
      GrNumber = 0; y_max = 82; y_min = 0
      
      # Find max for RhCases proportion
      MaxRhCases = 0; MinRhCases = 0
      for (i in 1:length(VectProp)) {
        for (j in 1:length(VectCovTreatment)) {
          MaxRhCases = max(c(MaxRhCases, 
                             ModelStrategy03_Fort[[i]][[k]][[j]]$RhCases / ModelStrategy03_Fort[[i]][[k]][[j]]$IhTot,
                             ModelStrategy03_Faible[[i]][[k]][[j]]$RhCases / ModelStrategy03_Faible[[i]][[k]][[j]]$IhTot))
        }
      }
      MaxRhCases=100*MaxRhCases
      
      for (i in 1:length(VectProp)) {
        for (j in 1:length(VectCovTreatment)) {
          GrNumber = GrNumber + 1
          
          time = ModelStrategy03_Fort[[i]][[k]][[j]]$time
          IdSeq = seq(1, length(time), by = 90)
          
          # Calculate reductions
          red_combined_Fort = calc_reduction(ModelStrategy03_Fort[[i]][[k]][[j]], 
                                             ModelZero_R_Fort, "PropIh")
          red_combined_Faible = calc_reduction(ModelStrategy03_Faible[[i]][[k]][[j]], 
                                               ModelZero_R_Faible, "PropIh")
          
          # Plot reduction (primary axis)
          plot(-1, 1, type = "l", xlab = "", xlim = c(0, max(time)), 
               ylab = "", ylim = c(y_min, y_max), cex.lab = 1.2, xaxt = "n")
          axis(1, at = seq(0, Tmax, by = 120), labels = seq(0, Tmaxmonths, by = 4), las = 1)
          
          # Highlight intervention period
          polygon(c(t_begin_Camp, Tmax, Tmax, t_begin_Camp), 
                  c(y_min, y_min, y_max, y_max),
                  col = adjustcolor("#D3D3D3", alpha.f = 0.1), border = NA)
          
          # Plot reduction curves
          lines(time, red_combined_Fort, lwd = LC, col = ColVect[1], lty = LineVect[1])
          points(time[IdSeq], red_combined_Fort[IdSeq], 
                 pch = VectPch[1], cex = 1.2, col = ColVect[1])
          
          lines(time, red_combined_Faible, lwd = LC, col = ColVect[1], lty = LineVect[1])
          points(time[IdSeq], red_combined_Faible[IdSeq], 
                 pch = VectPch[2], cex = 1.2, col = ColVect[1])
          
          # Vertical line at intervention start
          abline(v = t_begin_Camp, lty = 2, col = "#205072", lwd = 1.5)
          
          
          # Labels for primary axis
          if (GrNumber %in% c(1, 4, 7)) {
            mtext("Clinical cases averted (%)", side = 2, adj = 0.5, cex = .8, line = 2.5)
            mtext(RowText[i], side = 2, adj = 0.5, cex = 0.75, line = 3.5)
          }
          
          if (GrNumber %in% c(7, 8, 9)) {
            mtext("Time (months)", side = 1, adj = 0.5, cex = .85, line = 2.5)
          }
          
          if (GrNumber %in% c(1, 2, 3)) {
            mtext(ColText[j], side = 3, adj = 0.5, cex = 0.85, line = 0.5)
          }
          
          # Add secondary axis for RhCases proportion
          par(new = TRUE)
          
          plot(-1, 1, type = "l", xlab = "", xlim = c(0, max(time)), 
               ylab = "", axes = FALSE, ylim = c(MinRhCases, MaxRhCases))
          
          # Plot RhCases as proportion (secondary axis)
          lines(time, 100*ModelStrategy03_Fort[[i]][[k]][[j]]$RhCases / ModelStrategy03_Fort[[i]][[k]][[j]]$IhTot, 
                lwd = LC, col = ColVect[2], lty = LineVect[2])
          points(time[IdSeq], 100*ModelStrategy03_Fort[[i]][[k]][[j]]$RhCases[IdSeq] / ModelStrategy03_Fort[[i]][[k]][[j]]$IhTot[IdSeq], 
                 pch = VectPch[1], cex = 1.2, col = ColVect[2])
          
          lines(time, 100*ModelStrategy03_Faible[[i]][[k]][[j]]$RhCases / ModelStrategy03_Faible[[i]][[k]][[j]]$IhTot, 
                lwd = LC, col = ColVect[2], lty = LineVect[2])
          points(time[IdSeq], 100*ModelStrategy03_Faible[[i]][[k]][[j]]$RhCases[IdSeq] / ModelStrategy03_Faible[[i]][[k]][[j]]$IhTot[IdSeq], 
                 pch = VectPch[2], cex = 1.2, col = ColVect[2])
          
          axis(4, col = "black")
          if (GrNumber == 1) {
            text(t_begin_Camp, MaxRhCases * 0.95, "Start", cex = 0.75, col = "#205072")
          }
          
          if (GrNumber %in% c(3, 6, 9)) {
            mtext("Proportion of resistant cases (%)", side = 4, adj = 0.5, cex = .75, line = 2)
          }
          
          # Panel letter
          par(xpd = NA)
          text(min(time) - 0.05, MaxRhCases * 1.1, 
               paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.3, adj = 0)
          par(xpd = FALSE)
          
          if (GrNumber == 8) {
            legend("bottom", 
                   legend = c("Clinical cases averted", "Clinical resistant cases"),
                   xpd = NA, horiz = TRUE, lwd = 2, inset = c(-3, -0.5),
                   bty = "n", col = ColVect, lty = LineVect, cex = 1.3, text.col = 1)
          }
        }
      }
      
      # Legend for epsi values
      LEGEND = c(expression(epsilon[s] == 0.1 ~ "(Strong)"), 
                 expression(epsilon[s] == 0.01 ~ "(Weak)"))
      
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      
      # Title
      title(main = bquote("# cycle for MDA and SMC:" ~ n[c] == .(VectNumber_of_cycle[k])), 
            cex.main = 1.4, line = -1)
      
      # Legend
      legend("bottom", legend = LEGEND, xpd = NA, horiz = TRUE, 
             bty = "n", pch = VectPch, col = "black", cex = 1.4)
      
      dev.off()
    }
  }
  
  plot_epsi_scenario()
  
}

# SCÉNARIO DURABILITÉ
{
  scenario_sustainability = function(params_base, a_lim = 5, coverage_level = 0.7,
                                     treatment_levels = c(0.1, 0.5, 0.9),
                                     initial_conditions = initial_conditions_eq) {
    intervention_years = 1
    post_intervention_years = 1
    total_years = intervention_years + post_intervention_years
    cycle_options = c(3, 4, 5)
    
    results = list()
    total_sims = 1 + length(treatment_levels) * 3 * length(cycle_options)
    pb = txtProgressBar(min = 0, max = total_sims, style = 3)
    counter = 0
    
    # Baseline
    counter = counter + 1
    setTxtProgressBar(pb, counter)
    baseline = run_simulation(params_base, PropTreatment = 0, PropSMC = 0, PropMDA = 0, Nyr = total_years, 
                              initial_conditions = initial_conditions)
    results[["baseline"]] = baseline
    
    for (treat_level in treatment_levels) {
      for (n_cyc in cycle_options) {
        # S1: Treatment seul
        counter = counter + 1
        setTxtProgressBar(pb, counter)
        
        s1_intervention = run_simulation(params_base, PropTreatment = treat_level, PropSMC = 0, 
                                         PropMDA = 0, a_lim = a_lim, Nyr = intervention_years,
                                         initial_conditions = initial_conditions)
        
        s1_post_init = list(
          mSh_final = s1_intervention$mSh_final, fSh_final = s1_intervention$fSh_final,
          mAh_final = s1_intervention$mAh_final, fAh_final = s1_intervention$fAh_final,
          mIh_final = s1_intervention$mIh_final, fIh_final = s1_intervention$fIh_final,
          mRh_final = s1_intervention$mRh_final, fRh_final = s1_intervention$fRh_final,
          mSh_smc_final = s1_intervention$mSh_smc_final, fSh_smc_final = s1_intervention$fSh_smc_final,
          mAh_smc_final = s1_intervention$mAh_smc_final, fAh_smc_final = s1_intervention$fAh_smc_final,
          mSh_mda_final = s1_intervention$mSh_mda_final, fSh_mda_final = s1_intervention$fSh_mda_final,
          mAh_mda_final = s1_intervention$mAh_mda_final, fAh_mda_final = s1_intervention$fAh_mda_final,
          Sm_final = s1_intervention$Sm_final, Im_final = s1_intervention$Im_final,
          Dh_final = s1_intervention$Dh_final
        )
        
        s1_post = run_simulation(params_base, PropTreatment = treat_level, PropSMC = 0, PropMDA = 0,
                                 Nyr = post_intervention_years, initial_conditions = s1_post_init)
        
        results[[paste0("S1_treat", treat_level*100, "_cycle", n_cyc)]] = list(intervention = s1_intervention, post = s1_post)
        
        # S2: MDA+SMC
        counter = counter + 1
        setTxtProgressBar(pb, counter)
        
        s2_intervention = run_simulation(params_base, PropTreatment = 0, PropSMC = coverage_level, 
                                         PropMDA = coverage_level, a_lim = a_lim, Number_of_cycle = n_cyc,
                                         Nyr = intervention_years, initial_conditions = initial_conditions)
        
        s2_post_init = list(
          mSh_final = s2_intervention$mSh_final, fSh_final = s2_intervention$fSh_final,
          mAh_final = s2_intervention$mAh_final, fAh_final = s2_intervention$fAh_final,
          mIh_final = s2_intervention$mIh_final, fIh_final = s2_intervention$fIh_final,
          mRh_final = s2_intervention$mRh_final, fRh_final = s2_intervention$fRh_final,
          mSh_smc_final = s2_intervention$mSh_smc_final, fSh_smc_final = s2_intervention$fSh_smc_final,
          mAh_smc_final = s2_intervention$mAh_smc_final, fAh_smc_final = s2_intervention$fAh_smc_final,
          mSh_mda_final = s2_intervention$mSh_mda_final, fSh_mda_final = s2_intervention$fSh_mda_final,
          mAh_mda_final = s2_intervention$mAh_mda_final, fAh_mda_final = s2_intervention$fAh_mda_final,
          Sm_final = s2_intervention$Sm_final, Im_final = s2_intervention$Im_final,
          Dh_final = s2_intervention$Dh_final
        )
        
        s2_post = run_simulation(params_base, PropTreatment = 0, PropSMC = 0, PropMDA = 0,
                                 Nyr = post_intervention_years, initial_conditions = s2_post_init)
        
        results[[paste0("S2_treat", treat_level*100, "_cycle", n_cyc)]] = list(intervention = s2_intervention, post = s2_post)
        
        # S1+S2: Combinée
        counter = counter + 1
        setTxtProgressBar(pb, counter)
        
        s3_intervention = run_simulation(params_base, PropTreatment = treat_level, PropSMC = coverage_level,
                                         PropMDA = coverage_level, a_lim = a_lim, Number_of_cycle = n_cyc,
                                         Nyr = intervention_years, initial_conditions = initial_conditions)
        
        s3_post_init = list(
          mSh_final = s3_intervention$mSh_final, fSh_final = s3_intervention$fSh_final,
          mAh_final = s3_intervention$mAh_final, fAh_final = s3_intervention$fAh_final,
          mIh_final = s3_intervention$mIh_final, fIh_final = s3_intervention$fIh_final,
          mRh_final = s3_intervention$mRh_final, fRh_final = s3_intervention$fRh_final,
          mSh_smc_final = s3_intervention$mSh_smc_final, fSh_smc_final = s3_intervention$fSh_smc_final,
          mAh_smc_final = s3_intervention$mAh_smc_final, fAh_smc_final = s3_intervention$fAh_smc_final,
          mSh_mda_final = s3_intervention$mSh_mda_final, fSh_mda_final = s3_intervention$fSh_mda_final,
          mAh_mda_final = s3_intervention$mAh_mda_final, fAh_mda_final = s3_intervention$fAh_mda_final,
          Sm_final = s3_intervention$Sm_final, Im_final = s3_intervention$Im_final,
          Dh_final = s3_intervention$Dh_final
        )
        
        s3_post = run_simulation(params_base, PropTreatment = treat_level, PropSMC = 0, PropMDA = 0,
                                 Nyr = post_intervention_years, initial_conditions = s3_post_init)
        
        results[[paste0("S1+S2_treat", treat_level*100, "_cycle", n_cyc)]] = list(intervention = s3_intervention, post = s3_post)
      }
    }
    
    close(pb)
    
    # Calculer métriques
    metrics = list()
    baseline_prev_intervention = mean(baseline$PropIh[1:(intervention_years*360)])
    baseline_prev_post = mean(baseline$PropIh[(intervention_years*360+1):(total_years*360)])
    
    for (treat_level in treatment_levels) {
      for (n_cyc in cycle_options) {
        metrics_cycle = data.frame(
          strategy = character(), n_cycles = integer(), treatment_level = numeric(),
          reduction_intervention = numeric(), reduction_post = numeric(),
          resistance_intervention = numeric(), resistance_post = numeric(),
          rebound_speed = numeric(), sustained_benefit = numeric(),
          stringsAsFactors = FALSE
        )
        
        for (strat in c("S1", "S2", "S1+S2")) {
          strat_key = paste0(strat, "_treat", treat_level*100, "_cycle", n_cyc)
          
          # Réductions
          prev_interv = mean(results[[strat_key]]$intervention$PropIh[1:(intervention_years*360)])
          prev_post = mean(results[[strat_key]]$post$PropIh[1:(post_intervention_years*360)])
          
          reduction_interv = 100 * (1 - prev_interv / baseline_prev_intervention)
          reduction_post = 100 * (1 - prev_post / baseline_prev_post)
          
          # Résistance
          resist_interv = mean(results[[strat_key]]$intervention$RhCases[1:(intervention_years*360)] / 
                                 results[[strat_key]]$intervention$IhTot[1:(intervention_years*360)], na.rm = TRUE)
          resist_post = mean(results[[strat_key]]$post$RhCases[1:(post_intervention_years*360)] / 
                               results[[strat_key]]$post$IhTot[1:(post_intervention_years*360)], na.rm = TRUE)
          
          # Vitesse de rebond
          rebound_speed = 100 * (prev_post - prev_interv) / (baseline_prev_intervention - prev_interv)
          
          # Bénéfice soutenu
          sustained = reduction_post
          
          metrics_cycle = rbind(metrics_cycle,
                                data.frame(
                                  strategy = strat, n_cycles = n_cyc, treatment_level = treat_level,
                                  reduction_intervention = reduction_interv,
                                  reduction_post = reduction_post,
                                  resistance_intervention = resist_interv * 100,
                                  resistance_post = resist_post * 100,
                                  rebound_speed = rebound_speed,
                                  sustained_benefit = sustained,
                                  stringsAsFactors = FALSE
                                ))
        }
        
        metrics[[paste0("treat", treat_level*100, "_cycle", n_cyc)]] = metrics_cycle
      }
    }
    
    save(results, metrics, intervention_years, post_intervention_years, total_years, 
         coverage_level, treatment_levels, cycle_options, a_lim,
         file = paste0("scenario_sustainability_aLim", a_lim, ".RData"))
    return(list(results = results, metrics = metrics))
  }
  
  # FONCTION DE PLOT COMBINÉ
  plot_sustainability_combined = function(a_lim = 5) {
    load(paste0("scenario_sustainability_aLim", a_lim, ".RData"))
    
    
    pdf(paste0("Sustainability_aLim", a_lim, ".pdf"), width = 11, height = 6.8)
    par(mfrow = c(2, 3), mar = c(3, 2.5, 2, 2), oma = c(3, 2.5, 2, 2))
    
    colors = c("S1" = "#ff3355", "S2" = "#037153", "S1+S2" = "purple")
    strategies_ordered = c("S1", "S2", "S1+S2")
    lytvect = c("S1" = 2, "S2" = 4, "S1+S2" = 1)
    
    # PREMIÈRE LIGNE: Évolution temporelle combinée (nc=4)
    n_cyc = 4
    
    for (col in 1:3) {
      treat_level = treatment_levels[col]
      
      # Plot avec deux axes Y
      plot(NULL, xlim = c(0, total_years*360), ylim = c(0, 65), xlab = "", 
           ylab = "", main = "", frame.plot = FALSE, axes = FALSE)
      mtext(bquote("Treatment: " ~ .(treat_level * 100) * "%" ~ "(" ~ n[c] == .(n_cyc) ~ ")"), 
            side = 3, adj = 0.5, cex = 0.8, line = 0.5)
      
      axis(1, at = seq(0, total_years*360, by = 120), labels = seq(0, total_years*12, by = 4))
      axis(2, las = 1,  at = seq(0, 65, by = 15), labels = seq(0, 65, by = 15))
      if (col == 1) {
        mtext("Clinical cases averted (%)", side = 2, adj = 0.5, cex = 0.8, line = 2.5)
      }
      axis(4, las = 1,  at = seq(0, 65, by = 15), labels = seq(0, 65, by = 15))
      if (col == 3) {
        mtext("Clinical resistant cases (%)", side = 4, adj = 0.5, cex = 0.8, line = 2.5)
      }
      # Zone grisée et ligne verticale
      rect(0, 0, intervention_years*360, 65, col = adjustcolor("#D3D3D3", alpha = 0.25), border = NA)
      abline(v = intervention_years*360, lty = 2, lwd = 1.5, col = "black")
      if (col == 1) {
        text(intervention_years*360, 60, "Relaxing\nintervention", cex = 0.9, col = "black")
      }
      
      baseline_prev_int = results[["baseline"]]$PropIh[1:(intervention_years*360)]
      baseline_prev_post = results[["baseline"]]$PropIh[(intervention_years*360+1):(total_years*360)]
      
      # Tracer les courbes pour chaque stratégie
      for (strat in strategies_ordered) {
        strat_key = paste0(strat, "_treat", treat_level*100, "_cycle", n_cyc)
        
        # Phase intervention
        prev_interv = results[[strat_key]]$intervention$PropIh[1:(intervention_years*360)]
        time_interv = results[[strat_key]]$intervention$time[1:(intervention_years*360)]
        reduction_interv = 100 * (1 - prev_interv / baseline_prev_int)
        lines(time_interv, reduction_interv, col = colors[strat], lty = lytvect[strat], lwd = 2.5)
        
        # Phase post
        prev_post = results[[strat_key]]$post$PropIh[1:(post_intervention_years*360)]
        time_post = results[[strat_key]]$post$time[1:(post_intervention_years*360)] + intervention_years*360
        reduction_post = 100 * (1 - prev_post / baseline_prev_post)
        lines(time_post, reduction_post, col = colors[strat], lty = lytvect[strat], lwd = 2.5)
        
        # RESISTANCE
        # Phase intervention
        resist_interv_full = results[[strat_key]]$intervention$RhCases / results[[strat_key]]$intervention$IhTot
        resist_interv_full[!is.finite(resist_interv_full)] = NA
        resist_interv = resist_interv_full[1:(intervention_years*360)]
        # Mise à l'échelle pour l'axe Y gauche
        resist_interv_scaled = resist_interv * 100   
        lines(time_interv, resist_interv_scaled, col = colors[strat], lty = 3, lwd = 2.5)
        
        # Phase post
        resist_post_full = results[[strat_key]]$post$RhCases / results[[strat_key]]$post$IhTot
        resist_post_full[!is.finite(resist_post_full)] = NA
        resist_post = resist_post_full[1:(post_intervention_years*360)]
        resist_post_scaled = resist_post * 100
        lines(time_post, resist_post_scaled, col = colors[strat], lty = 3, lwd = 2.5)
      }
      mtext("Time (months)", side = 1, line = 2.5, cex = 0.8)
      
      # Légende uniquement pour le premier panel
      if (col == 1) {
        legend("topleft", legend = c(strategies_ordered, "Resistance"),
               col = c( colors, "black"), lty = c(lytvect,  3),
               lwd = c(rep(2, 3), 2), cex = 0.95, bty = "n")
      }
      
      # Étiquette du panel
      panel_labels = c("(A)", "(B)", "(C)")
      mtext(panel_labels[col], side = 3, adj = 0, cex = 0.9, line = 0.5)
      
      box()
    }
    
    # DEUXIÈME LIGNE: Barplots pour chaque niveau de traitement
    for (col in 1:3) {
      treat_level = treatment_levels[col]
      
      reduction_during = reduction_after = numeric(9)
      idx = 1
      
      for (cyc in cycle_options) {
        for (strat in strategies_ordered) {
          metric_data = metrics[[paste0("treat", treat_level*100, "_cycle", cyc)]]
          subset = metric_data[metric_data$strategy == strat,]
          if (nrow(subset) > 0) {
            reduction_during[idx] = subset$reduction_intervention[1]
            reduction_after[idx] = subset$reduction_post[1]
          }
          idx = idx + 1
        }
      }
      
      # Créer les données pour barplot côte à côte
      all_values = c(rbind(reduction_during, reduction_after))
      
      # Créer les couleurs alternées
      bar_colors = rep(c("gray30", "gray70"), 9)
      
      # Définir les espaces entre les barres
      spaces = numeric(18)
      for (i in 1:18) {
        if (i %% 2 == 1) {  # Première barre de chaque paire
          if (i %in% c(1, 7, 13)) {  # Premier groupe de chaque nc
            spaces[i] = 0.65
          } else {
            spaces[i] = 0.2
          }
        } else {  # Deuxième barre de chaque paire
          spaces[i] = 0.15
        }
      }
      
      bp = barplot(all_values, col = bar_colors, axes = FALSE, xlab = "", 
                   ylab = "",  main ="", space = spaces, border = NA, ylim = c(0, 65))
      
      axis(2, las = 1)
      if (col == 1) {
        mtext("Annual cases averted (%)", side = 2, adj = 0.5, cex = 0.8, line = 2.5)
      }
      
      # Ajouter les labels des stratégies
      label_positions = bp[seq(1, 17, 2)] + diff(bp[1:2])/2
      axis(1, at = label_positions, labels = rep(strategies_ordered, 3), cex.axis = 0.8, las = 1)
      
      # Ajouter les labels nc
      nc_positions = c(mean(bp[3:4]), mean(bp[9:10]), mean(bp[15:16]))
      text(x = nc_positions, y = -max(all_values)*0.19, 
           labels = c(expression(n[c]==3), expression(n[c]==4), expression(n[c]==5)),
           xpd = TRUE, cex = 1.1, font = 2)
      
      box()
      
      if (col == 1) {
        legend("topleft", legend = c("First year", "second year"),
               fill = c("gray30", "gray70"), bty = "n", cex = 1.2)
      }
      
      # Étiquette du panel
      panel_labels = c("(D)", "(E)", "(F)")
      mtext(panel_labels[col], side = 3, adj = 0, cex = 0.9, line = 0.5)
      
    }
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    
    legend("bottom", legend = c("S1: Treatment", "S2: MDA+SMC", "S1+S2: Treatment+MDA+SMC"),
           fill = colors, horiz = TRUE, bty = "n", cex = 1.1, xpd = NA, inset = c(0, 0.01))
    
    dev.off()
  }
  
  # FONCTION DE COMPARAISON STANDARD vs EXTENDED
  plot_comparison_combined = function() {
    # Charger les données Standard
    load("scenario_sustainability_aLim5.RData")
    results_std = results
    metrics_std = metrics
    treatment_levels_std = treatment_levels
    intervention_years_std = intervention_years
    post_intervention_years_std = post_intervention_years
    total_years_std = total_years
    cycle_options_std = cycle_options
    
    # Charger les données Extended
    load("scenario_sustainability_aLim10.RData")
    results_ext = results
    metrics_ext = metrics
    
    # Utiliser les variables de temps depuis le chargement
    intervention_years = intervention_years_std
    post_intervention_years = post_intervention_years_std
    total_years = total_years_std
    cycle_options = cycle_options_std
    
    pdf("Comparison_Standard_vs_Extended.pdf", width = 10.5, height = 6.8)
    par(mfrow = c(2, 3), mar = c(3, 2.5, 2, 2), oma = c(3, 2.5, 2, 2))
    
    colors = c("S2" = "#037153", "S1+S2" = "purple")
    strategies_ordered = c("S2", "S1+S2")
    lytvect = c("S2" = 4, "S1+S2" = 1)
    n_cyc = 4
    
    # PREMIÈRE LIGNE: Différences temporelles
    for (col in 1:3) {
      treat_level = treatment_levels_std[col]
      
      plot(NULL, xlim = c(0, total_years*360), ylim = c(-0.5, 3.5), xlab = "", ylab = "",
           main = "",  frame.plot = FALSE, axes = FALSE)
      
      mtext(bquote("Difference (Treatment: " ~ .(treat_level * 100) * "%" ~ "," ~ n[c] == .(n_cyc) ~ ")"), 
            side = 3, adj = 0.5, cex = 0.8, line = 0.5)
      mtext("Time (months)", side = 1, line = 2.5, cex = 0.8)
      
      axis(1, at = seq(0, total_years*360, by = 120), labels = seq(0, total_years*12, by = 4))
      axis(2, las = 1, at = seq(-0.5, 3.5, by = 1))
      if (col == 1) {
        mtext("Clinical cases averted (%)", side = 2, adj = 0.5, cex = 0.8, line = 2.5)
      }
      # Axe Y droite (Resistance diff)
      axis(4, las = 1, at = seq(-0.5, 3.5, by = 1))
      if (col == 3) {
        mtext("Clinical resistant cases(%)", side = 4, adj = 0.5, cex = 0.8, line = 2.5)
      }
      
      rect(0, -0.5, intervention_years*360, 3.5, col = adjustcolor("#D3D3D3", alpha = 0.25), border = NA)
      abline(h = 0, lty = 2, col = "gray50", lwd = 1)
      abline(v = intervention_years*360, lty = 2, lwd = 1.5, col = "black")
      if (col == 1) {
        text(intervention_years*360, 2.5, "Relaxing\nintervention", cex = 0.9, col = "black")
      }
      
      for (strat in strategies_ordered) {
        strat_key = paste0(strat, "_treat", treat_level*100, "_cycle", n_cyc)
        
        # Différence CASES AVERTED (ligne normale)
        baseline_prev = results_std[["baseline"]]$PropIh
        
        # Phase intervention
        prev_std = results_std[[strat_key]]$intervention$PropIh[1:(intervention_years*360)]
        prev_ext = results_ext[[strat_key]]$intervention$PropIh[1:(intervention_years*360)]
        reduction_std = 100 * (1 - prev_std / baseline_prev[1:(intervention_years*360)])
        reduction_ext = 100 * (1 - prev_ext / baseline_prev[1:(intervention_years*360)])
        diff_interv = reduction_std - reduction_ext
        
        # Phase post
        prev_std_post = results_std[[strat_key]]$post$PropIh[1:(post_intervention_years*360)]
        prev_ext_post = results_ext[[strat_key]]$post$PropIh[1:(post_intervention_years*360)]
        baseline_post = baseline_prev[(intervention_years*360+1):(total_years*360)]
        reduction_std_post = 100 * (1 - prev_std_post / baseline_post)
        reduction_ext_post = 100 * (1 - prev_ext_post / baseline_post)
        diff_post = reduction_std_post - reduction_ext_post
        
        lines(1:(intervention_years*360), diff_interv, 
              col = colors[strat], lty = lytvect[strat], lwd = 2.5)
        lines((intervention_years*360+1):(total_years*360), diff_post, 
              col = colors[strat], lty = lytvect[strat], lwd = 2.5)
        
        # Différence RESISTANCE (lty=3)
        resist_std_full = results_std[[strat_key]]$intervention$RhCases / results_std[[strat_key]]$intervention$IhTot
        resist_ext_full = results_ext[[strat_key]]$intervention$RhCases / results_ext[[strat_key]]$intervention$IhTot
        resist_std_full[!is.finite(resist_std_full)] = NA
        resist_ext_full[!is.finite(resist_ext_full)] = NA
        resist_std = resist_std_full[1:(intervention_years*360)]
        resist_ext = resist_ext_full[1:(intervention_years*360)]
        diff_resist_interv = (resist_std - resist_ext) * 100
        
        resist_std_post_full = results_std[[strat_key]]$post$RhCases / results_std[[strat_key]]$post$IhTot
        resist_ext_post_full = results_ext[[strat_key]]$post$RhCases / results_ext[[strat_key]]$post$IhTot
        resist_std_post_full[!is.finite(resist_std_post_full)] = NA
        resist_ext_post_full[!is.finite(resist_ext_post_full)] = NA
        resist_std_post = resist_std_post_full[1:(post_intervention_years*360)]
        resist_ext_post = resist_ext_post_full[1:(post_intervention_years*360)]
        diff_resist_post = (resist_std_post - resist_ext_post) * 100
        
        lines(1:(intervention_years*360), diff_resist_interv, 
              col = colors[strat], lty = 3, lwd = 2, na.rm = TRUE)
        lines((intervention_years*360+1):(total_years*360), diff_resist_post, 
              col = colors[strat], lty = 3, lwd = 2, na.rm = TRUE)
      }
      
      
      # Légende uniquement pour le premier panel
      if (col == 1) {
        legend("topleft", legend = c(strategies_ordered, "Resistance"),
               col = c( colors, "black"), lty = c(lytvect,  3),
               lwd = c(rep(2, 3), 2), cex = 0.95, bty = "n")
      }
      
      panel_labels = c("(A)", "(B)", "(C)")
      mtext(panel_labels[col], side = 3, adj = 0, cex = 0.9, line = 0.5)
      box()
    }
    
    # DEUXIÈME LIGNE: Barplots des différences
    for (col in 1:3) {
      treat_level = treatment_levels_std[col]
      
      # Préparer les données
      diff_during = numeric(6)  # 2 strategies × 3 cycles
      diff_post = numeric(6)
      idx = 1
      
      for (cyc in cycle_options) {
        for (strat in strategies_ordered) {
          metric_std = metrics_std[[paste0("treat", treat_level*100, "_cycle", cyc)]]
          metric_ext = metrics_ext[[paste0("treat", treat_level*100, "_cycle", cyc)]]
          
          diff_during[idx] = metric_std[metric_std$strategy == strat, "reduction_intervention"] - 
            metric_ext[metric_ext$strategy == strat, "reduction_intervention"]
          diff_post[idx] = metric_std[metric_std$strategy == strat, "reduction_post"] - 
            metric_ext[metric_ext$strategy == strat, "reduction_post"]
          idx = idx + 1
        }
      }
      
      # Créer les données pour barplot
      all_values = c(rbind(diff_during, diff_post))
      bar_colors = rep(c("gray30", "gray70"), 6)
      
      # Espaces entre barres
      spaces = numeric(12)
      for (i in 1:12) {
        if (i %% 2 == 1) {
          if (i %in% c(1, 5, 9)) {
            spaces[i] = 0.65
          } else {
            spaces[i] = 0.2
          }
        } else {
          spaces[i] = 0.15
        }
      }
      
      bp = barplot(all_values, col = bar_colors, axes = FALSE, xlab = "", 
                   ylab = "", main = "", space = spaces, border = NA, ylim = c(-0.2, 3.5))
      
      abline(h = 0, lty = 2, col = "gray50")
      axis(2, las = 1)
      if (col == 1) {
        mtext("Annual Cases averted(%)", side = 2, adj = 0.5, cex = 0.8, line = 2.5)
      }
      
      # Labels des stratégies
      label_positions = bp[seq(1, 11, 2)] + diff(bp[1:2])/2
      axis(1, at = label_positions, labels = rep(strategies_ordered, 3), cex.axis = 0.8, las = 1)
      
      # Labels nc
      nc_positions = c(mean(bp[2:3]), mean(bp[6:7]), mean(bp[10:11]))
      text(x = nc_positions, y = -0.5, 
           labels = c(expression(n[c]==3), expression(n[c]==4), expression(n[c]==5)),
           xpd = TRUE, cex = 1.1, font = 2)
      
      box()
      
      if (col == 1) {
        legend("topleft",legend = c("First year", "second year"),
               fill = c("gray30", "gray70"), bty = "n", cex = 1.2)
      }
      
      panel_labels = c("(D)", "(E)", "(F)")
      mtext(panel_labels[col], side = 3, adj = 0, cex = 0.9, line = 0.5)
      
    }
    dev.off()
    
    cat("Plots saved as: Sustainability_Combined_aLim*.pdf and Comparison_Standard_Extended_Combined.pdf\n")
  }
  
  # EXÉCUTION
  params_base = params_base(epsi_s = 1e-1)
  
  # Run Standard (a_lim = 5)
  result_std = scenario_sustainability(params_base, a_lim = 5, treatment_levels = c(0.1, 0.5, 0.9))
  plot_sustainability_combined(a_lim = 5)
  
  # Run Extended (a_lim = 10)
  params_base = function(t_begin_Camp = 120, epsi_s = 1e-1) {
    list(Gap = 1, dt = 1, p_f = 0.15, LengYr = 360, t_begin_Camp = t_begin_Camp, Time_plot0 = 0,
         epsi_s = epsi_s, MDA_Dur_cycle = 7, SMC_Dur_cycle = 7)
  }
  params_base = params_base(epsi_s = 1e-1)
  result_ext = scenario_sustainability(params_base, a_lim = 10, treatment_levels = c(0.1, 0.5, 0.9))
  plot_sustainability_combined(a_lim = 10)
  
  # Comparaison
  plot_comparison_combined()
}


# SCÉNARIO RECHERCHE DE COUVERTURE MINIMALE AVEC TRAITEMENT FIXÉ
scenario_Boxplot_min_coverage = function(params_base, a_lim = 5, treatment_levels = c(0.1, 0.5, 0.9), 
                                         thresholds = c(10, 20, 30),
                                         initial_conditions = initial_conditions_eq) {
  
  coverage_grid = c(0.05, seq(0.1, 0.9, by = 0.1))
  cycle_options = c(3, 4, 5)
  
  results = list()
  
  # Initialiser les résultats pour chaque traitement et seuil
  for (treat_val in treatment_levels) {
    for (thresh in thresholds) {
      key = paste0("treatment_", treat_val*100, "_threshold_", thresh)
      results[[key]] = data.frame(
        treatment_level = numeric(), strategy = character(), 
        coverage_smc = numeric(), coverage_mda = numeric(),
        n_cycles = integer(), a_lim = integer(), 
        mean_reduction = numeric(), mean_resistant_prop = numeric(), 
        achieves_threshold = logical(),  stringsAsFactors = FALSE
      )
    }
  }
  
  # Baseline sans intervention
  baseline = run_simulation(params_base, PropTreatment = 0, PropSMC = 0,PropMDA = 0, 
                            Nyr = 1, initial_conditions = initial_conditions)
  baseline_mean = mean(baseline$PropIh[1:360])
  
  # Compteur pour progression
  total_configs = length(treatment_levels) * (
    length(coverage_grid) * length(cycle_options) +  # SMC seul
      length(coverage_grid) * length(cycle_options) +  # MDA seul
      length(coverage_grid) * length(coverage_grid) * length(cycle_options)  # SMC+MDA (échantillonnage)
  )
  
  pb = txtProgressBar(min = 0, max = total_configs, style = 3)
  counter = 0
  
  # Boucle sur chaque niveau de traitement fixé
  for (treat_val in treatment_levels) {
    
    # Stratégie 1: Treatment + SMC seul
    for (n_cyc in cycle_options) {
      for (cov_smc in coverage_grid) {
        counter = counter + 1
        setTxtProgressBar(pb, counter)
        
        tryCatch({
          sim = run_simulation(params_base, PropTreatment = treat_val, PropSMC = cov_smc, 
                               PropMDA = 0, a_lim = a_lim,  Number_of_cycle = n_cyc, 
                               Nyr = 1,  initial_conditions = initial_conditions)
          
          sim_mean = mean(sim$PropIh[1:360])
          reduction = 100 * (1 - sim_mean / baseline_mean)
          resistant_prop = mean(sim$RhCases[1:360] / sim$IhTot[1:360], na.rm = TRUE)
          
          for (thresh in thresholds) {
            achieves = reduction >= thresh
            key = paste0("treatment_", treat_val*100, "_threshold_", thresh)
            new_row = data.frame(
              treatment_level = treat_val,  strategy = "S",  # SMC seul
              coverage_smc = cov_smc, coverage_mda = 0,
              n_cycles = n_cyc, a_lim = a_lim, 
              mean_reduction = reduction, mean_resistant_prop = resistant_prop,
              achieves_threshold = achieves, stringsAsFactors = FALSE
            )
            results[[key]] = rbind(results[[key]], new_row)
          }
        }, error = function(e) {})
      }
    }
    
    # Stratégie 2: Treatment + MDA seul
    for (n_cyc in cycle_options) {
      for (cov_mda in coverage_grid) {
        counter = counter + 1
        setTxtProgressBar(pb, counter)
        
        tryCatch({
          sim = run_simulation(params_base,  PropTreatment = treat_val, PropSMC = 0, 
                               PropMDA = cov_mda,  a_lim = a_lim, Number_of_cycle = n_cyc,  
                               Nyr = 1, initial_conditions = initial_conditions)
          
          sim_mean = mean(sim$PropIh[1:360])
          reduction = 100 * (1 - sim_mean / baseline_mean)
          resistant_prop = mean(sim$RhCases[1:360] / sim$IhTot[1:360], na.rm = TRUE)
          
          for (thresh in thresholds) {
            achieves = reduction >= thresh
            key = paste0("treatment_", treat_val*100, "_threshold_", thresh)
            new_row = data.frame(
              treatment_level = treat_val,  strategy = "M",  # MDA seul
              coverage_smc = 0,  coverage_mda = cov_mda, n_cycles = n_cyc, 
              a_lim = a_lim,   mean_reduction = reduction, 
              mean_resistant_prop = resistant_prop, achieves_threshold = achieves, 
              stringsAsFactors = FALSE
            )
            results[[key]] = rbind(results[[key]], new_row)
          }
        }, error = function(e) {})
      }
    }
    
    # Stratégie 3: Treatment + SMC + MDA
    for (n_cyc in cycle_options) {
      # Échantillonnage pour réduire le nombre de simulations
      smc_grid = c(0.05, seq(0.1, 0.9, by = 0.1))  
      mda_grid = c(0.05, seq(0.1, 0.9, by = 0.1))
      
      for (cov_smc in smc_grid) {
        for (cov_mda in mda_grid) {
          counter = counter + 1
          if (counter <= total_configs) setTxtProgressBar(pb, counter)
          
          tryCatch({
            sim = run_simulation(params_base, PropTreatment = treat_val, PropSMC = cov_smc, 
                                 PropMDA = cov_mda, a_lim = a_lim, Number_of_cycle = n_cyc, 
                                 Nyr = 1, initial_conditions = initial_conditions)
            
            sim_mean = mean(sim$PropIh[1:360])
            reduction = 100 * (1 - sim_mean / baseline_mean)
            resistant_prop = mean(sim$RhCases[1:360] / sim$IhTot[1:360], na.rm = TRUE)
            
            for (thresh in thresholds) {
              achieves = reduction >= thresh
              key = paste0("treatment_", treat_val*100, "_threshold_", thresh)
              new_row = data.frame(
                treatment_level = treat_val, strategy = "S+M",  # SMC + MDA
                coverage_smc = cov_smc, coverage_mda = cov_mda,
                n_cycles = n_cyc, a_lim = a_lim, 
                mean_reduction = reduction, mean_resistant_prop = resistant_prop,
                achieves_threshold = achieves, stringsAsFactors = FALSE
              )
              results[[key]] = rbind(results[[key]], new_row)
            }
          }, error = function(e) {})
        }
      }
    }
  }
  close(pb)
  
  # Déterminer les configurations minimales
  minimal_coverage = list()
  
  for (treat_val in treatment_levels) {
    for (thresh in thresholds) {
      key = paste0("treatment_", treat_val*100, "_threshold_", thresh)
      thresh_data = results[[key]]
      achieving = thresh_data[thresh_data$achieves_threshold,]
      
      if (nrow(achieving) > 0) {
        min_configs = data.frame()
        
        strategies = c("S", "M", "S+M")
        for (strat in strategies) {
          for (n_cyc in cycle_options) {
            strat_cyc_data = achieving[achieving$strategy == strat & 
                                         achieving$n_cycles == n_cyc,]
            
            if (nrow(strat_cyc_data) > 0) {
              # Calculer l'effort total (somme des couvertures)
              strat_cyc_data$total_effort = strat_cyc_data$coverage_smc +
                strat_cyc_data$coverage_mda
              
              # Trouver la configuration minimale
              min_idx = which.min(strat_cyc_data$total_effort)
              
              if (length(min_idx) > 0) {
                min_configs = rbind(min_configs, strat_cyc_data[min_idx[1],])
              }
            }
          }
        }
        minimal_coverage[[key]] = min_configs
      }
    }
  }
  
  save(results, minimal_coverage, treatment_levels, thresholds, a_lim,
       file = paste0("Scenario_FixedTreatment_aLim", a_lim, ".RData"))
  
  return(list(full_results = results, minimal = minimal_coverage))
}

# FONCTION DE PLOT Boxplot 
plot_scenario_fixed_treatment = function(a_lim = 5) {
  load(paste0("Scenario_FixedTreatment_aLim", a_lim, ".RData"))
  
  pdf(paste0("Scenario_FixedTreatment_aLim", a_lim, "_Boxplots.pdf"), width = 11, height = 6.6)
  par(mfrow = c(2, 3), oma = c(3, 2, 1.5, 1), mar = c(3.5, 4, 2, 1))
  
  strategy_colors = c("S" = "#037153", "M" = "#8B4513", "S+M" = "purple")
  strategies_ordered = c("S", "M", "S+M")
  cycle_options = c(3, 4, 5)
  panel_labels = c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)")
  
  panel_idx = 1
  
  # Ligne 1: Cases Averted
  for (treat_idx in 1:length(treatment_levels)) {
    treat_val = treatment_levels[treat_idx]
    
    # Combiner toutes les données pour ce niveau de traitement
    treat_data = do.call(rbind, results[grep(paste0("treatment_", treat_val*100), names(results))])
    
    # Préparer les données pour le boxplot
    boxplot_data = list()
    labels_vec = c()
    colors_vec = c()
    
    for (cyc in cycle_options) {
      for (strat in strategies_ordered) {
        data_subset = treat_data[treat_data$strategy == strat & treat_data$n_cycles == cyc,]
        if (nrow(data_subset) > 0) {
          valid_data = data_subset$mean_reduction[!is.na(data_subset$mean_reduction)]
          if (length(valid_data) > 0) {
            boxplot_data[[length(boxplot_data) + 1]] = valid_data
            labels_vec = c(labels_vec, "")
            colors_vec = c(colors_vec, strategy_colors[strat])
          }
        }
      }
    }
    
    boxplot(boxplot_data, col = colors_vec, xlab = "", ylab =  "", main = "",ylim = c(0, 65),
            names = labels_vec, las = 1, cex.axis = 0.8, axes = FALSE, frame.plot = TRUE)
    
    if (treat_idx == 1) {
      mtext("Annual cases averted (%)", side = 2, adj = 0.5, cex = 0.85, font = 1, line = 2)
    }
    
    # Titre du traitement (en haut de chaque colonne)
    mtext(paste0("Treatment = ", treat_val*100, "%"), side = 3, line = 1, cex = 0.9, font = 2)
    
    # Label du panel
    mtext(panel_labels[panel_idx], side = 3, adj = 0, cex = 1, font = 1, line = 0.5)
    panel_idx = panel_idx + 1
    
    # Titre de la ligne (première colonne seulement)
    if (treat_idx == 1) {
      mtext("Cases averted distribution", side = 2, adj = 0.5, cex = 0.85, font = 1, line = 3)
    }
    
    axis(2, las = 1)
    
    # Labels des stratégies
    strategy_positions = 1:9
    axis(1, at = strategy_positions, labels = rep(strategies_ordered, 3), cex.axis = 0.7, las = 1, line = 0)
    # 
    # # Labels des cycles
    # nc_positions = c(2, 5, 8)
    # mtext(expression(n[c] == 3), side = 1, at = nc_positions[1], line = 2.1, cex = 0.75, font = 2)
    # mtext(expression(n[c] == 4), side = 1, at = nc_positions[2], line = 2.1, cex = 0.75, font = 2)
    # mtext(expression(n[c] == 5), side = 1, at = nc_positions[3], line = 2.1, cex = 0.75, font = 2)
  }
  
  # Ligne 2: Resistant Cases
  for (treat_idx in 1:length(treatment_levels)) {
    treat_val = treatment_levels[treat_idx]
    
    # Combiner toutes les données pour ce niveau de traitement
    treat_data = do.call(rbind, results[grep(paste0("treatment_", treat_val*100), names(results))])
    
    # Préparer les données pour le boxplot
    boxplot_data_res = list()
    labels_vec_res = c()
    colors_vec_res = c()
    
    for (cyc in cycle_options) {
      for (strat in strategies_ordered) {
        data_subset = treat_data[treat_data$strategy == strat & treat_data$n_cycles == cyc,]
        if (nrow(data_subset) > 0) {
          valid_data = data_subset$mean_resistant_prop[!is.na(data_subset$mean_resistant_prop)]
          if (length(valid_data) > 0) {
            boxplot_data_res[[length(boxplot_data_res) + 1]] = valid_data * 100
            labels_vec_res = c(labels_vec_res, "")
            colors_vec_res = c(colors_vec_res, strategy_colors[strat])
          }
        }
      }
    }
    
    boxplot(boxplot_data_res, col = colors_vec_res, xlab = "",  ylab = "",main = "",ylim = c(0, 8),
            names = labels_vec_res, las = 1, cex.axis = 0.8, axes = FALSE, frame.plot = TRUE)
    
    if (treat_idx == 1) {
      mtext("Proportion of resistant cases (%)", side = 2, adj = 0.5, cex = 0.8, font = 1, line = 2)
    }
    
    # Label du panel
    mtext(panel_labels[panel_idx], side = 3, adj = 0, cex = 1, font = 1, line = 0.5)
    panel_idx = panel_idx + 1
    
    # Titre de la ligne (première colonne seulement)
    if (treat_idx == 1) {
      mtext("Resistant cases distribution", side = 2, adj = 0.5, cex = 0.85, font = 1, line = 3)
    }
    
    axis(2, las = 1)
    
    # Labels des stratégies
    strategy_positions = 1:9
    axis(1, at = strategy_positions, labels = rep(strategies_ordered, 3), 
         cex.axis = 0.7, las = 1, line = 0)
    
    # Labels des cycles
    nc_positions = c(2, 5, 8)
    mtext(expression(n[c] == 3), side = 1, at = nc_positions[1], line = 2.1, cex = 0.75, font = 2.1)
    mtext(expression(n[c] == 4), side = 1, at = nc_positions[2], line = 2.1, cex = 0.75, font = 2.1)
    mtext(expression(n[c] == 5), side = 1, at = nc_positions[3], line = 2.1, cex = 0.75, font = 2.1)
  }
  
  # Légende globale
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  
  legend("bottom", legend = c("S: SMC", "M: MDA", "S+M: SMC + MDA"),
         fill = strategy_colors, horiz = TRUE, bty = "n", 
         cex = 1, xpd = NA, inset = c(0, 0.02))
  
  dev.off()
}

# FONCTION DE PLOT COUVERTURES MINIMALES 
plot_minimal_coverage_fixed_treatment = function(a_lim = 5) {
  load(paste0("Scenario_FixedTreatment_aLim", a_lim, ".RData"))
  
  pdf(paste0("Minimal_Coverage_FixedTreatment_aLim", a_lim, ".pdf"), width = 12, height = 8.2)
  par(mfrow = c(3, 3), oma = c(3, 2, 2, 1), mar = c(3, 3.7, 1.5, .5))
  
  strategy_colors = c("S" = "#037153", "M" = "#8B4513", "S+M" = "purple")
  strategies_ordered = c("S", "M", "S+M")
  cycle_options = c(3, 4, 5)
  thresholds = c(10, 20, 30)
  panel_labels = c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)", "(G)", "(H)", "(I)")
  
  panel_idx = 1
  
  # Pour chaque seuil (lignes)
  for (thresh_idx in 1:length(thresholds)) {
    thresh = thresholds[thresh_idx]
    
    # Pour chaque niveau de traitement (colonnes)
    for (treat_idx in 1:length(treatment_levels)) {
      treat_val = treatment_levels[treat_idx]
      
      # Préparer les données pour le barplot
      bar_data_smc = numeric(9)   # 3 cycles × 3 stratégies
      bar_data_mda = numeric(9)
      bar_colors = character(9)
      
      # Récupérer les données de couverture minimale
      key = paste0("treatment_", treat_val*100, "_threshold_", thresh)
      min_data = minimal_coverage[[key]]
      
      # Construire les données pour le barplot
      idx = 1
      for (cyc in cycle_options) {
        for (strat in strategies_ordered) {
          bar_colors[idx] = strategy_colors[strat]
          
          if (!is.null(min_data) && nrow(min_data) > 0) {
            subset = min_data[min_data$strategy == strat & min_data$n_cycles == cyc,]
            
            if (nrow(subset) > 0) {
              bar_data_smc[idx] = subset$coverage_smc[1]
              bar_data_mda[idx] = subset$coverage_mda[1]
            } else {
              bar_data_smc[idx] = NA
              bar_data_mda[idx] = NA
            }
          } else {
            bar_data_smc[idx] = NA
            bar_data_mda[idx] = NA
          }
          
          idx = idx + 1
        }
      }
      
      # Créer la matrice pour le barplot empilé
      # Pour S+M, empiler SMC (gray30) et MDA (gray70)
      bar_matrix = rbind(bar_data_smc * 100, bar_data_mda * 100)
      
      # Identifier les positions NA
      na_positions = which(is.na(bar_data_smc) & is.na(bar_data_mda))
      
      # Remplacer NA par 0 pour l'affichage
      bar_matrix[is.na(bar_matrix)] = 0
      
      bp = barplot(bar_matrix, col = c(adjustcolor("gray30", alpha = 1),  # SMC
                                       adjustcolor("gray70", alpha = 1)),  # MDA
                   ylim = c(0, 160), ylab = "", xlab = "", border = NA, 
                   space = c(0.25, 0.25, 0.25, 0.6, 0.25, 0.25, 0.6, 0.25, 0.25), axes = FALSE)
      
      # Label pour le seuil (première colonne de chaque ligne)
      if (treat_idx == 1) {
        mtext(paste0(thresh, "% Reduction"), side = 2, line = 3.5, cex = 0.85, font = 1)
        mtext("Minimal coverage (%)", side = 2, adj = 0.5, cex = 0.75, line = 2)
      }
      
      # Titre du traitement (première ligne)
      if (thresh_idx == 1) {
        mtext(paste0("Treatment = ", treat_val*100, "%"), 
              side = 3, adj = 0.5, cex = 0.95, font = 1, line = 0.5)
      }
      
      # Label du panel
      mtext(panel_labels[panel_idx], side = 3, adj = 0, cex = 1, font = 1, line = 0.5)
      panel_idx = panel_idx + 1
      
      # Ajouter l'axe Y
      axis(2, las = 1)
      
      # Labels des stratégies
      axis(1, at = bp, labels = rep(strategies_ordered, 3), 
           cex.axis = 0.85, las = 1)
      
      # Marquer les cas non atteignables
      for (i in na_positions) {
        x_pos = bp[i]
        segments(x_pos - 0.3, 5, x_pos + 0.3, 15, col = "red", lwd = 1.5)
        segments(x_pos - 0.3, 15, x_pos + 0.3, 5, col = "red", lwd = 1.5)
        text(x_pos, 20, "NA", col = "red", cex = 0.75, font = 2)
      }
      
      # Labels des cycles
      cycle_positions = c(mean(bp[1:3]), mean(bp[4:6]), mean(bp[7:9]))
      mtext(expression(n[c] == 3), side = 1, at = cycle_positions[1], line = 2.5, cex = 0.7, font = 2)
      mtext(expression(n[c] == 4), side = 1, at = cycle_positions[2], line = 2.5, cex = 0.7, font = 2)
      mtext(expression(n[c] == 5), side = 1, at = cycle_positions[3], line = 2.5, cex = 0.7, font = 2)
      
      box(bty = "o")
      
      # Légende dans le premier panel
      if (panel_idx == 2) {
        legend("topleft", legend = c("SMC coverage", "MDA coverage", "Not achievable"),
               fill = c(adjustcolor("gray30", alpha = 1),  adjustcolor("gray70", alpha = 1), NA),
               pch = c(NA, NA, 4), 
               col = c(NA, NA, "red"),
               bty = "n", cex = 0.75, pt.cex = 1.5)
      }
    }
  }
  
  # Légende globale
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  
  legend("bottom", legend = c("S: SMC", "M: MDA", "S+M: SMC + MDA"), fill = "black",
         horiz = TRUE, bty = "n", cex = 1.1, xpd = NA, inset = c(0, 0.02))
  
  dev.off()
}

# EXÉCUTION
params_base = params_base(epsi_s=1e-1)
result = scenario_Boxplot_min_coverage(params_base, a_lim=5)
plot_scenario_fixed_treatment(a_lim=5)
plot_minimal_coverage_fixed_treatment(a_lim=5)


#Sensitivity analysis
{
  VectCovTreatment = c(0.1, 0.5, 0.95)
  VectNumber_of_cycle = c(3, 4, 5)
  VectPropMDA = c(0.1, 0.5, 0.95)
  VectPropSMC = c(0.1, 0.5, 0.95)
  VectInitPrevR=c(0.01,0.1)
  Vectp_f = c(0.1, 0.2)
  Vecta_lim = c(5, 10)
  
  
  SizePerCombination = length(VectCovTreatment) * length(VectPropMDA) * 
    length(VectPropSMC) * length(VectInitPrevR) * length(Vectp_f)
  print(paste("Simulations par combinaison (cycle, a_lim):", SizePerCombination))
  
  combinaisons_a_faire = list( c(4,5), c(4,10), c(5,5), c(5,10)) #list(c(3,5), c(3,10))
  
  # MAIN LOOP FOR EACH COMBINATION
  for (i in 1:length(combinaisons_a_faire)) {
    
    Number_of_cycle = combinaisons_a_faire[[i]][1]
    a_lim_fixe = combinaisons_a_faire[[i]][2]
    
    # PARAMETRIC EXPLORATION PHASE
    ColCovTreatment = rep(NA, SizePerCombination)
    ColPropMDA = rep(NA, SizePerCombination)
    ColPropSMC = rep(NA, SizePerCombination)
    ColInitPrevR = rep(NA,SizePerCombination)
    Colp_f = rep(NA, SizePerCombination)
    ColCases_advert = rep(NA, SizePerCombination)
    ColResis_Cases = rep(NA, SizePerCombination)
    
    
    RunNumber = 0
    
    for (p_f in Vectp_f) {
      # EQUILIBRIUM PHASE
      params_eq = params_base(epsi_s = 1e-1)
      params_eq$p_f = p_f
      Equilibrium = run_simulation_Eq(params_eq, PropTreatment = 0, PropSMC = 0, 
                                      PropMDA = 0, a_lim = a_lim_fixe)
      
      initial_conditions_eq = list(
        mSh_final = Equilibrium$mSh_final, fSh_final = Equilibrium$fSh_final,
        mAh_final = Equilibrium$mAh_final, fAh_final = Equilibrium$fAh_final,
        mIh_final = Equilibrium$mIh_final, fIh_final = Equilibrium$fIh_final,
        mRh_final = Equilibrium$mRh_final, fRh_final = Equilibrium$fRh_final,
        mSh_smc_final = Equilibrium$mSh_smc_final, fSh_smc_final = Equilibrium$fSh_smc_final,
        mAh_smc_final = Equilibrium$mAh_smc_final, fAh_smc_final = Equilibrium$fAh_smc_final,
        mSh_mda_final = Equilibrium$mSh_mda_final, fSh_mda_final = Equilibrium$fSh_mda_final,
        mAh_mda_final = Equilibrium$mAh_mda_final, fAh_mda_final = Equilibrium$fAh_mda_final,
        Sm_final = Equilibrium$Sm_final, Im_final = Equilibrium$Im_final, 
        Dh_final = Equilibrium$Dh_final
      )
      params_run = params_base(epsi_s = 1e-1)
      params_run$p_f = p_f
      Test_result_zero = run_simulation(params_run, PropTreatment = 0, PropSMC =0, PropMDA = 0, a_lim = a_lim_fixe,
                                        initial_conditions = initial_conditions_eq )
      
      for (CovTreatment in VectCovTreatment) {
        for (PropMDA in VectPropMDA) {
          for (PropSMC in VectPropSMC) {
            for (InitPrevR in VectInitPrevR) {
              
              
              RunNumber = RunNumber + 1
              print(paste("Cycle", Number_of_cycle, ", a_lim", a_lim_fixe,
                          "- Run", RunNumber, "/", SizePerCombination))
              
              params_run = params_base(epsi_s = 1e-1)
              params_run$p_f = p_f
              
              tryCatch({
                Test_result = run_simulation(params_run, PropTreatment = CovTreatment,
                                             PropSMC = PropSMC, PropMDA = PropMDA, InitPrevR=InitPrevR,
                                             a_lim = a_lim_fixe, Number_of_cycle = Number_of_cycle,
                                             Saisonality = 1, Nyr = 1, initial_conditions = initial_conditions_eq )
                
                Cases_advert_t = calc_reduction(Test_result, Test_result_zero, "PropIh")
                
                ColCovTreatment[RunNumber] = CovTreatment
                ColPropMDA[RunNumber] = PropMDA
                ColPropSMC[RunNumber] = PropSMC
                ColInitPrevR[RunNumber] = InitPrevR
                Colp_f[RunNumber] = p_f
                ColCases_advert[RunNumber] = sum(Cases_advert_t[])
                ColResis_Cases[RunNumber] = Test_result$Resis_Cases
                
              }, error = function(e) {
                ColCovTreatment[RunNumber] = CovTreatment
                ColPropMDA[RunNumber] = PropMDA
                ColPropSMC[RunNumber] = PropSMC
                ColInitPrevR[RunNumber] = InitPrevR
                Colp_f[RunNumber] = p_f
                ColCases_advert[RunNumber] = NA
                ColResis_Cases[RunNumber] = NA
              })
            }
          }
        }
      }
    }
    
    Tableau_combinaison = data.frame(CovTreatment = ColCovTreatment,PropMDA = ColPropMDA,
                                     PropSMC = ColPropSMC, InitPrevR = ColInitPrevR, 
                                     Number_of_cycle = Number_of_cycle, p_f = Colp_f,
                                     a_lim = a_lim_fixe,  Cases_advert = ColCases_advert,
                                     Resis_Cases = ColResis_Cases
    )
    
    nom_fichier_combinaison = paste0("TabMDA_Cycle_", Number_of_cycle, "_aLim_", a_lim_fixe, ".csv")
    write.table(Tableau_combinaison, file = nom_fichier_combinaison, sep = ",", row.names = FALSE)
  }
  
  # ASSEMBLY FUNCTION
  assembler_tableaux = function() {
    VectNumber_of_cycle = c(3, 4, 5)
    Vecta_lim = c(5, 10)
    tous_tableaux = list()
    
    for (cycle_val in VectNumber_of_cycle) {
      for (a_lim_val in Vecta_lim) {
        nom_fichier = paste0("TabMDA_Cycle_", cycle_val, "_aLim_", a_lim_val, ".csv")
        
        if (file.exists(nom_fichier)) {
          tableau_combinaison = read.csv(nom_fichier, stringsAsFactors = FALSE)
          tous_tableaux[[paste(cycle_val, a_lim_val, sep="_")]] = tableau_combinaison
        }
      }
    }
    
    if (length(tous_tableaux) > 0) {
      Tableau_complet = do.call(rbind, tous_tableaux)
      write.table(Tableau_complet, file = "Table_MDA_SMC_Treatment_complet10.csv", sep = ",", row.names = FALSE)
      return(Tableau_complet)
    }
  }
  
  # EXECUTION
  Tableau_final = assembler_tableaux()
  
  # SAVE WORKSPACE
  save.image(file = "ANOVA_Analysis_Results.RData")
}

getwd()
