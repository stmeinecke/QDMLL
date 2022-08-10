//no S - with p.DK[k]; interpolation selected at compile time by macros below
// #define XCat(tau, var) (Xhist->at_C(tau,var))
// #define XRat(tau, var) (Xhist->at_R(tau,var))

#define XCat(tau, var) (Xhist->at_IntPol_3rd_C(tau,var))
#define XRat(tau, var) (Xhist->at_IntPol_3rd_R(tau,var))

// #define XCat(tau, var) (Xhist->at_C_3rd_wSmallTau(tau,var))
// #define XRat(tau, var) (Xhist->at_R_3rd_wSmallTau(tau,var))

// #define XCat(tau, var) (Xhist->at_C_fast3rd(tau,var))
// #define XRat(tau, var) (Xhist->at_R_fast3rd(tau,var))



auto MLL_system_eqs = [](vars *X, vars_vec *Xhist, vars *d, parameters *p){
  

  
  //sum over subgroups
  for(unsigned int k = 0; k<NSEC;k++){
    X->P_plus[k] = 0.0;
    X->P_minus[k] = 0.0;
    for(unsigned int l = 0; l<NSAMPLES; l++){
      X->P_plus[k] += p->QD_inh[l] * (2.0 * X->rho_GS[k][l] - 1.0) * X->G_plus[k][l];
      X->P_minus[k] += p->QD_inh[l] * (2.0 * X->rho_GS[k][l] - 1.0) * X->G_minus[k][l];
    }
    X->P_plus[k] *= p->g[k];
    X->P_minus[k] *= p->g[k];
  }
  
  
  for(unsigned int k = 0; k<NSEC;k++){
    X->P_ES[k] = 0.0;
    for(unsigned int l = 0; l<NSAMPLES; l++){
      X->P_ES[k] += p->QD_inh[l] * p->idelta_omega_ES[k][l] * X->rho_ES[k][l];
    }
  }


  //left boundary
  X->A_plus[0] = (XCat(p->DT[0], IND::A_minus[0]) * sqrt(p->r_L) * ( p->alpha_int_A[0] - XRat(p->DT[0], IND::P_ES[0])*p->delta_z_alpha_int_S[0] )
                    + p->delta_z_alpha_int_S[0] * ( X->P_plus[0] + sqrt(p->r_L) * XCat(p->DT[0], IND::P_minus[0]) ) )
                    / ( 1.0 + X->P_ES[0]*p->delta_z_alpha_int_S[0]);
  
  for(int k = 1; k < NSEC ; k++){
    X->A_plus[k] = (XCat(p->DT[k], IND::A_plus[k-1]) * ( p->alpha_int_A[k-1] - XRat(p->DT[k], IND::P_ES[k-1])*p->delta_z_alpha_int_S[k-1] )
                    + (p->delta_z_alpha_int_S[k] * X->P_plus[k] + p->delta_z_alpha_int_S[k-1] * XCat(p->DT[k], IND::P_plus[k-1]) ) )
                    / ( 1.0 + X->P_ES[k]*p->delta_z_alpha_int_S[k]);
  }
  
  for(int k = 0; k < NSEC-1 ; k++){
    X->A_minus[k] = (XCat(p->DT[k+1], IND::A_minus[k+1]) * ( p->alpha_int_A[k] - XRat(p->DT[k+1], IND::P_ES[k+1])*p->delta_z_alpha_int_S[k+1] )
                    + (p->delta_z_alpha_int_S[k] * X->P_minus[k] + p->delta_z_alpha_int_S[k+1] * XCat(p->DT[k+1], IND::P_minus[k+1]) ) )
                    / ( 1.0 + X->P_ES[k]*p->delta_z_alpha_int_S[k]);
  }
  
  //right boundary
  X->A_minus[NSEC-1] = (XCat(p->DT[NSEC], IND::A_plus[NSEC-1]) * sqrt(p->r_R) * ( p->alpha_int_A[NSEC-1] - XRat(p->DT[NSEC], IND::P_ES[NSEC-1])*p->delta_z_alpha_int_S[NSEC-1] )
                    + p->delta_z_alpha_int_S[NSEC-1] * ( X->P_minus[NSEC-1] + sqrt(p->r_R) * XCat(p->DT[NSEC], IND::P_plus[NSEC-1]) ) )
                    / ( 1.0 + X->P_ES[NSEC-1]*p->delta_z_alpha_int_S[NSEC-1]);
                    
                    
  //feedback
  //simple
  X->A_out = sqrt(1.0-p->r_R) * XCat(p->DT[NSEC],IND::A_plus[NSEC-1]) + p->z[NSEC-1]*0.5*XCat(p->DT[NSEC],IND::P_plus[NSEC-1]);
  //multiple reflections in the external cavity
  X->A_out += sqrt(p->r_R*p->K)*exp(sm::img * p->C)* XCat(p->tau,IND::A_out);
  
  X->A_minus[NSEC-1] += (1.0 - p->r_R) * sqrt(p->K) * exp(sm::img * p->C) * XCat(p->DT[NSEC]+p->tau,IND::A_out);
  
  
  //dynamical equations
  for(int k = 0; k < NSEC ; k++){
    
    //carrier
    
    d->N[k] = p->J[k] - p->gamma_N[k] * X->N[k];
    
    for(unsigned int l = 0; l<NSAMPLES; l++){
      
      const double expGSES = exp((p->QD_eps_GS[l] - p->QD_eps_ES[l])*0.8/0.025);
      const double expESQW = exp((p->QD_eps_ES[l] - p->QD_confine)*0.8/0.025);
      
      const double rel = p->Rrel[k] * ( (1.0-X->rho_GS[k][l])*X->rho_ES[k][l] - X->rho_GS[k][l]*(1.0-X->rho_ES[k][l])*expGSES);
      const double cap = p->Rcap[k] * (1.0/(1.0 + expESQW/(exp(X->N[k]/p->D_N)-1.0)) - X->rho_ES[k][l]);
      
      d->rho_ES[k][l] = -p->gamma_ES[k] * X->rho_ES[k][l] + cap - 0.5 * rel;
      
      d->rho_GS[k][l] = -p->gamma_GS[k] * X->rho_GS[k][l] + rel - p->g[k] * p->invetasqrd[k] * (2.0 * X->rho_GS[k][l] - 1.0) * (real(X->A_plus[k] * conj(X->G_plus[k][l])) + real(X->A_minus[k] * conj(X->G_minus[k][l])));
      
      d->N[k] -= p->QD_ES_DOS[l]*cap;
      
    }
    
    //slow gain
    for(unsigned int l = 0; l < NSAMPLES; l++){
      d->G_plus[k][l] = p->gamma_2[k] * (X->A_plus[k] - X->G_plus[k][l]) + p->iDelta_omega[k][l] * X->G_plus[k][l];
      d->G_minus[k][l] = p->gamma_2[k] * (X->A_minus[k] - X->G_minus[k][l]) + p->iDelta_omega[k][l] * X->G_minus[k][l];
    }
   
  }
  
};



 
