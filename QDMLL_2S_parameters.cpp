struct parameters{
  //class to collect all parameters of the system.
  
  //define parameters here:
  
  double dt;
  
  double v_g; 
  double T0;
  //geometry of the laser
  unsigned int nSec;
  unsigned int nSecG;
  unsigned int nSecA;
  bool NoA;

  double device_length;
  double l_G;
  double l_A;
  
  double r_L;
  double r_R;
  
  double J_Q;
  double J_G;
  bool setJG;
  double P_G;
  double injEf;
  
  double gamma_N_Q;
  double gamma_N_G;
  double gamma_ES_Q;
  double gamma_ES_G;
  double gamma_GS_Q;
  double gamma_GS_G;
  
  double T_2_G;
  double T_2_Q;

  double g_Q;
  double g_G;
  
  double etasqrd;
  
  double Delta_omega_Q;
  double Delta_omega_G;
  
  
  double U;
  bool setU;
  
  double SARS;
  double SARSat10V;
  
  double gamma_ES_Q_0;
  double gamma_ES_Q_U;
  double Delta_omega_Q_U0;
  
  
  double Rcap_Q;
  double Rcap_G;
  double Rrel_Q;
  double Rrel_G;
  
  double N_QD;
  double D_N;
  double exp_ESGS;
  double exp_ESQW;
  double ESQW;
  
    
  double w_base;
  double A_G;
  double A;
  
  double alpha_int_base;
  
  double Gamma_r_base;
  double Gamma_base;
  
  double SqrtNoiseStr;
  double noiseStr;
  
  
  //gain paramters
  double Gamma_L;
  double Gamma;
  double e_ph;
  double omega;
  double a_L;
  double mu_GS;
  double mu_ES;
  double eps_r;
  double h_QW;
  double nu_GS;
  double nu_ES;
  
  
  //QD subgroup parametes
  double QD_GSES;
  double QD_ESQW;
  double QD_confine;
  double QD_GS_inh_FWHM;
  double QD_ES_inh_FWHM;
  double QD_eps_GS[NSAMPLES];
  double QD_eps_ES[NSAMPLES];
  double QD_omega_GS[NSAMPLES];
  double QD_omega_ES[NSAMPLES];
  double QD_inh[NSAMPLES];
  double QD_GS_DOS[NSAMPLES];
  double QD_ES_DOS[NSAMPLES];
  
  
  //feedback paras
  double K;
  double tau;
  double C;
  
  
  //section parameters
  double z[NSEC];
  double T[NSEC];
  double DT[NSEC+1];
  int SA[NSEC];
  
  double alpha_int[NSEC];
  double alpha_int_A[NSEC];
  double delta_z_alpha_int_S[NSEC];
  double g_A_Reff[NSEC];
  double gamma_2[NSEC];
  
  double gamma_N[NSEC];
  double gamma_ES[NSEC];
  double gamma_GS[NSEC];
  double Rcap[NSEC];
  double Rrel[NSEC];
  double J[NSEC];
  double g[NSEC];
  double invetasqrd[NSEC];
  std::complex<double> idelta_omega_ES[NSEC][NSAMPLES];
  
  std::complex<double> iDelta_omega[NSEC][NSAMPLES];

  double w[NSEC];
//   double Gamma[NSEC];
//   double Gamma_r[NSEC];
  
  
  double AvCounter;
  int secSEED;

  
  parameters(){}
}; 


void init_states(parameters &p){
    
    //transition energy sampling step -> considered range from -FWHM to FWHM
    double QD_eps_GS_step = 2.0 * p.QD_GS_inh_FWHM / (NSAMPLES-1.0);
    double QD_eps_ES_step = 2.0 * p.QD_ES_inh_FWHM / (NSAMPLES-1.0);
    
    //set up transition energy sampling points
    for(unsigned int k = 0; k<NSAMPLES; k++){
      p.QD_eps_GS[k] = -p.QD_GS_inh_FWHM + k*QD_eps_GS_step;
      p.QD_eps_ES[k] = p.QD_GSES - p.QD_ES_inh_FWHM + k*QD_eps_ES_step;
    }
    
    //transition frequencies
    for(unsigned int k = 0; k<NSAMPLES; k++){
     p.QD_omega_GS[k] = ( p.QD_eps_GS[k] * sm::e0 / sm::hbar ) * (1e-12/tsf) ; // in ps
     p.QD_omega_ES[k] = ( p.QD_eps_ES[k] * sm::e0 / sm::hbar ) * (1e-12/tsf) ; // in ps
    }
    
    //Gaussian distribution of the energy states
    double norm = 0.0;
    for(unsigned int k = 0; k<NSAMPLES; k++){
      p.QD_inh[k] = exp(-4.0*log(2)*((p.QD_eps_GS[k]*p.QD_eps_GS[k])/(p.QD_GS_inh_FWHM*p.QD_GS_inh_FWHM)));
      norm += p.QD_inh[k];
    }
    
    //normalization such that QD_inh[k]  in [0,1] and the sum = 1.0
    for(unsigned int k = 0; k<NSAMPLES; k++){
      p.QD_inh[k] /= norm;
    }
    
    //calc QD DOS
    for(unsigned int k = 0; k<NSAMPLES; k++){
      p.QD_GS_DOS[k] = p.nu_GS * p.N_QD * p.QD_inh[k];
      p.QD_ES_DOS[k] = p.nu_ES * p.N_QD * p.QD_inh[k];
    }
    
}


//function that computes parameters from other parameters
void compute_parameters(parameters &p){
    
  p.A_G = p.l_G * p.w_base;
  
  if(p.NoA == true) p.A_G = p.l_G * p.w_base + p.l_A * p.w_base;
  p.A = p.l_G * p.w_base + p.l_A * p.w_base;
  
  double PJconversionF = (p.a_L * sm::e0 / p.injEf) * 1e12;
  if(p.setJG == true) p.P_G = p.J_G * p.A_G * PJconversionF;
  else p.J_G = p.P_G / (p.A_G * PJconversionF);
  std::cout << "pump current density: " << p.J_G << std::endl;
  
  if(p.setU){
    p.gamma_ES_Q = p.gamma_ES_Q_0 * exp((p.U)/p.gamma_ES_Q_U);
    p.Delta_omega_Q = p.SARSat10V * (p.U/p.Delta_omega_Q_U0);
  }
  else p.U = 2.0 * log(p.gamma_ES_Q / p.gamma_ES_Q_0);
  
  
  
  p.a_L = (int)p.a_L;
  p.Gamma = p.Gamma_L * p.a_L;
  
  //calc gain
  const double T2G = p.T_2_G * 1e-12;
  const double T2Q = p.T_2_Q * 1e-12;
  const double NQD = p.N_QD * 1e4;
  const double vg = p.v_g * 1e10;
  
  const double gG = p.omega * p.Gamma * T2G * NQD * p.nu_GS * (p.mu_GS * p.mu_GS) / ( 2.0 * p.eps_r * sm::eps0 * sm::hbar * vg * p.h_QW );
  const double gQ = p.omega * p.Gamma * T2Q * NQD * p.nu_GS * (p.mu_GS * p.mu_GS) / ( 2.0 * p.eps_r * sm::eps0 * sm::hbar * vg * p.h_QW );
  
//   p.g_G = p.omega * p.Gamma * p.T_2_G * p.a_L * p.N_QD * p.nu_GS * (p.mu_GS * p.mu_GS) / ( 2.0 * p.eps_r * sm::eps0 * sm::hbar * p.v_g * p.h_QW );
//   p.g_Q = p.omega * p.Gamma * p.T_2_Q * p.a_L * p.N_QD * p.nu_GS * (p.mu_GS * p.mu_GS) / ( 2.0 * p.eps_r * sm::eps0 * sm::hbar * p.v_g * p.h_QW );
  std::cout << "g_G: " << gG << std::endl;
  std::cout << "g_Q: " << gQ << std::endl;
  
  p.g_G = gG * 1e-2;
  p.g_Q = gQ * 1e-2;
  
  
//   std::cout << "p.v_g: " << p.v_g << std::endl;
}


//function that computes parameters for the individual sections
void compute_section_parameters(parameters &p){
  
  double zmin = p.device_length;
  for(std::size_t k = 0; k<NSEC; k++){
    if(k < p.nSecA){
      p.z[k] = p.l_A / (double)(p.nSecA);
    }
    else{
      p.z[k] = p.l_G / (double)(p.nSecG);
    }
    if(p.z[k] < zmin) zmin = p.z[k];
  }
  
  // randomize section lengths to avoid spurious frequencies
  p.secSEED = 817;
  std::mt19937 secGENERATOR(p.secSEED);
  std::uniform_real_distribution<double> secDIST(-1.0, 1.0);
//   std::normal_distribution<double> secDIST{0.0,1.0};
  secDIST(secGENERATOR);
  
  for(unsigned int k = 0; k < NSEC-1; k++){
    if(k != p.nSecA-1){
      const double scale = std::min(p.z[k],p.z[k+1])/10.0;      
      const double shift = secDIST(secGENERATOR)*scale;
      p.z[k] += shift;
      p.z[k+1] -= shift;
    }
  }
  
  //section propagation times
  for(std::size_t k = 0; k<NSEC; k++){
    p.T[k] = p.z[k] / p.v_g;
  }
  
  //delay coupling times
  p.DT[0] = p.T[0];
  for(unsigned int k = 1; k < NSEC; k++){
    p.DT[k] = 0.5*(p.T[k-1] + p.T[k]);
  }
  p.DT[NSEC] = p.T[NSEC-1];
  
//   for(unsigned int k = 0; k < NSEC; k++){
//     std::cout << p.DT[k] << std::endl;
//   }
  
  
  //laser parameters
  for(std::size_t k = 0; k<NSEC; k++){
    
    if(k < p.nSecA){
      if(p.NoA) p.SA[k] = 0;
      else p.SA[k] = 1;
    }
    else{
      p.SA[k] = 0;
    }
    
    p.w[k] = p.w_base;
    p.alpha_int[k] = p.alpha_int_base;
//     p.Gamma_r[k] = p.Gamma_r_base;

    //equation parameters
    if(p.SA[k] == 0){
//       p.gamma_2[k] = p.gamma_Gain_G;
      p.gamma_2[k] = 1.0/p.T_2_G;
      p.gamma_N[k] = p.gamma_N_G;
      p.gamma_ES[k] = p.gamma_ES_G;
      p.gamma_GS[k] = p.gamma_GS_G;
      p.J[k] = p.J_G;
      p.g[k] = p.g_G;
      p.Rcap[k] = p.Rcap_G;
      p.Rrel[k] = p.Rrel_G;
//       p.idelta_omega_ES[k] = sm::img * p.g[k] * 2.0 * p.delta_omega_G;

      for(unsigned int l = 0; l<NSAMPLES;l++){
        p.iDelta_omega[k][l] = sm::img * (p.QD_omega_GS[l] + p.Delta_omega_G);
        p.idelta_omega_ES[k][k] = sm::img * (p.g[k] * 2.0 * (p.T_2_G * p.QD_omega_ES[l])/(1.0 + (p.T_2_G * p.QD_omega_ES[l]) * (p.T_2_G * p.QD_omega_ES[l])));
      }

    }
    if(p.SA[k] == 1){
//       p.gamma_2[k] = p.gamma_Gain_Q;
      p.gamma_2[k] = 1.0/p.T_2_Q;
      p.gamma_N[k] = p.gamma_N_Q;
      p.gamma_ES[k] = p.gamma_ES_Q;
      p.gamma_GS[k] = p.gamma_GS_Q;
      p.J[k] = p.J_Q;
      p.g[k] = p.g_Q;
      p.Rcap[k] = p.Rcap_Q;
      p.Rrel[k] = p.Rrel_Q;
//       p.idelta_omega_ES[k] = sm::img * p.g[k] * 2.0 * p.delta_omega_Q;

      for(unsigned int l = 0; l<NSAMPLES;l++){
        p.iDelta_omega[k][l] = sm::img * (p.QD_omega_GS[l] + p.Delta_omega_Q);
        p.idelta_omega_ES[k][k] = sm::img * (p.g[k] * 2.0 * (p.T_2_Q * p.QD_omega_ES[l])/(1.0 + (p.T_2_Q * p.QD_omega_ES[l]) * (p.T_2_Q * p.QD_omega_ES[l])));
      }
    }
    
    p.invetasqrd[k] = (1.0/p.etasqrd);
    p.alpha_int_A[k] = (1 - 0.25*p.z[k]*p.alpha_int[k])/(1 + 0.25*p.z[k]*p.alpha_int[k]);
    p.delta_z_alpha_int_S[k] = p.z[k] / ( 2 + 0.5*p.z[k]*p.alpha_int[k]);
    
  }
    
}



