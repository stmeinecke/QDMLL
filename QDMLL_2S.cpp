#define NSEC 10
#define NSAMPLES 41

const double tsf = 1.0; // time scale factor
const double lsf = 1.0; // length scale factor

#include "incl/global.cpp"
#include "incl/get_params_sm.cc"
#include "incl/get_specs_sm_v6.cc"
#include "incl/get_correlations_sm.cc"
#include "QDMLL_2S_parameters.cpp"
#include "QDMLL_2S_vars.cpp"
#include "incl/vars_vec_v9.cpp"
#include "incl/DDEintegrator_v5.cpp"
#include "incl/evalMLTS_v13.cpp"
#include "QDMLL_2S_eqs.cpp"


using namespace std;

int main(int argc, char* argv[]){
  
  bool bquiet = katana::getCmdOption_bool(argv, argv+argc, "-quiet" , false, true);
  
  //dynamical variables indices
  IND::initIndices();
  
  std::string str_suffix = katana::getCmdOption(argv, argv+argc, "-strSuffix" , "");
  
  //noise
  katana::Seed s;
  s.set_seed();
  
  //solver
  DDEintegrator DDEsolver;
  //solver parameters
  unsigned long long int tn = 0;
  double intTime = katana::getCmdOption(argv, argv+argc, "-intTime" , 30000.0 * tsf);
  double dt = katana::getCmdOption(argv, argv+argc, "-dt" , 0.01 * tsf);
  double outTime = katana::getCmdOption(argv, argv+argc, "-outTime" , 3000.0  * tsf);
  if(outTime < 0.0) outTime = intTime;
  unsigned long long int outTime_ntn = (unsigned long long int)(outTime/dt);
  
  //thresholds and triggers for time series analysis, i.e. for pulse detection
  double T_maxPW = katana::getCmdOption(argv, argv+argc, "-T_maxPW" , 0.3); //max pulse width search radius in units of T0
  double T_dCountTol = katana::getCmdOption(argv, argv+argc, "-T_dCountTol" , 0.1); // tolerance for detecting unique maxima
  double T_relMax = katana::getCmdOption(argv, argv+argc, "-T_relMax" , 0.05); //threshold for maxima = pulses
  double T_TVPD_trig = katana::getCmdOption(argv, argv+argc, "-T_TVPD_trig" , 0.05); //trigger value for pulse detection
  double T_TVPD_bound = katana::getCmdOption(argv, argv+argc, "-T_TVPD_bound" , 0.001); //pulse boundaries for weighted pulse position/pulse area
  
  //system parameters in cm in ps
  parameters p; 
  //set default parameters:
  p.dt = dt;
  
  //device geometry
  p.T0 = 25.0 * tsf; //in ps
  p.v_g = 0.1 * lsf / (12.5 * tsf); // in cm/ps
  
  p.device_length = 0.5 * p.T0 * p.v_g;
  p.l_A = katana::getCmdOption(argv, argv+argc, "-l_A" , 0.1) * p.device_length;
  p.l_G = p.device_length - p.l_A;
  
  p.nSec = NSEC;  
  p.nSecA =  std::max(1,(int)round(p.nSec * p.l_A/p.device_length)); //guarantee at least one absorber section
  p.nSecG = p.nSec - p.nSecA;
  
  p.w_base = katana::getCmdOption(argv, argv+argc, "-w_base" , 0.0004*lsf); //4microns
  
  p.NoA = katana::getCmdOption_bool(argv, argv+argc, "-NoA" , false);

  //field propagation
  p.alpha_int_base = katana::getCmdOption(argv, argv+argc, "-alpha_int_base" , 2.0 / lsf);
  p.Gamma_r_base = katana::getCmdOption(argv, argv+argc, "-Gamma_r_base" , 1.0);
  
  p.r_L = katana::getCmdOption(argv, argv+argc, "-r_L" , 0.95);
  p.r_R = katana::getCmdOption(argv, argv+argc, "-r_R" , 0.31);
  
  
  //QD parameters
  p.N_QD = katana::getCmdOption(argv, argv+argc, "-N_QD" , 0.3E11 / (lsf * lsf) ); 
  p.QD_GSES = katana::getCmdOption(argv, argv+argc, "-QD_GSES" , 0.080); // eV
  p.QD_ESQW = katana::getCmdOption(argv, argv+argc, "-QD_ESQW" , 0.030); // eV
  p.QD_confine = katana::getCmdOption(argv, argv+argc, "-QD_confine" , p.QD_GSES+p.QD_ESQW); // eV
  
  p.QD_GS_inh_FWHM = katana::getCmdOption(argv, argv+argc, "-QD_inh" , 0.050); // in eV
  p.QD_ES_inh_FWHM = katana::getCmdOption(argv, argv+argc, "-QD_ES_inh" , 0.050); // in eV
  
  p.T_2_G = katana::getCmdOption(argv, argv+argc, "-T_2_G" , 0.1 * tsf);
  p.T_2_Q = katana::getCmdOption(argv, argv+argc, "-T_2_Q" , 0.2 * tsf);
  
  p.nu_GS = katana::getCmdOption(argv, argv+argc, "-nu_GS" , 2.0);
  p.nu_ES = katana::getCmdOption(argv, argv+argc, "-nu_ES" , 4.0);
  
  
  //QD relaxation rates
  p.gamma_N_Q = katana::getCmdOption(argv, argv+argc, "-gamma_N_Q" , 1.0 / tsf); //absorber -> doesnt matter -> is empty anyways
  p.gamma_N_G = katana::getCmdOption(argv, argv+argc, "-gamma_N_G" , 0.001 / tsf);
  p.gamma_ES_Q = katana::getCmdOption(argv, argv+argc, "-gamma_ES_Q" , 1.0 / tsf); //absorber
  p.gamma_ES_G = katana::getCmdOption(argv, argv+argc, "-gamma_ES_G" , 0.001 / tsf);
  p.gamma_GS_Q = katana::getCmdOption(argv, argv+argc, "-gamma_GS_Q" , 0.001 / tsf);
  p.gamma_GS_G = katana::getCmdOption(argv, argv+argc, "-gamma_GS_G" , 0.001 / tsf);
  
  //QD scattering rates
  p.Rrel_Q = katana::getCmdOption(argv, argv+argc, "-Rrel_Q" , 5.0 / tsf);
  p.Rrel_G = katana::getCmdOption(argv, argv+argc, "-Rrel_G" , 5.0 / tsf);
  p.Rcap_Q = katana::getCmdOption(argv, argv+argc, "-Rcap_Q" , 0.0 / tsf);  // 'exclude' carrier reservoire
  p.Rcap_G = katana::getCmdOption(argv, argv+argc, "-Rcap_G" , 0.1 / tsf);
  p.D_N = 2.5E12 / (lsf * lsf) ; //(0.24 electron mass / (Pi * hbar^2) * (25meV))
  
  
  //gain coefficient - in SI units to be later converted
  p.a_L = katana::getCmdOption(argv, argv+argc, "-a_L" , 10);
  p.Gamma_L = katana::getCmdOption(argv, argv+argc, "-Gamma_L" , 0.005); //optical confinement factor per QD layer
  p.Gamma = katana::getCmdOption(argv, argv+argc, "-Gamma" , p.Gamma_L * p.a_L);
  p.e_ph = katana::getCmdOption(argv, argv+argc, "-e_ph" , 0.97) * sm::e0; //1280 nm
  p.omega = katana::getCmdOption(argv, argv+argc, "-omega" , p.e_ph / sm::hbar);
  p.mu_GS = katana::getCmdOption(argv, argv+argc, "-mu_GS" , 0.6 * 1e-9) * sm::e0; 
  p.mu_ES = katana::getCmdOption(argv, argv+argc, "-mu_ES" , 0.6 * 1e-9) * sm::e0;
  p.eps_r = katana::getCmdOption(argv, argv+argc, "-eps_r" , 14.2);
  p.h_QW = katana::getCmdOption(argv, argv+argc, "-h_QW" , 5.0 * 1e-9); 
  
  //default gain 
  p.g_Q = katana::getCmdOption(argv, argv+argc, "-g_Q" , 60.0 / lsf);
  p.g_G = katana::getCmdOption(argv, argv+argc, "-g_G" , 40.0 / lsf);
  

  p.etasqrd = katana::getCmdOption(argv, argv+argc, "-etasqrd" , 1);
  
  
  // reverse bias parameters
  p.setU = katana::getCmdOption_bool(argv, argv+argc, "-U" , false);
  p.U = katana::getCmdOption(argv, argv+argc, "-U" , 6.0);
  
  p.gamma_ES_Q_0 = katana::getCmdOption(argv, argv+argc, "-gamma_ES_Q_0" , 0.05 / tsf);
  p.gamma_ES_Q_U = katana::getCmdOption(argv, argv+argc, "-gamma_ES_Q_U" , 2.0);
  p.Delta_omega_Q_U0 = katana::getCmdOption(argv, argv+argc, "-Delta_omega_Q_U0" , 10.0);
  
  p.SARSat10V = katana::getCmdOption(argv, argv+argc, "-SARSat10V" , 0.015) * sm::e0 / (sm::hbar * 1e12); //15meV redshift @ 10V -> converted in 1/ps
  p.SARS = katana::getCmdOption(argv, argv+argc, "-SARS" , 0.015*0.6) * sm::e0 / (sm::hbar * 1e12);
  p.Delta_omega_Q = katana::getCmdOption(argv, argv+argc, "-Delta_omega_Q" , p.SARS);
  p.Delta_omega_G = katana::getCmdOption(argv, argv+argc, "-Delta_omega_G" , 0.0);

  // pump parameters
  p.J_Q = katana::getCmdOption(argv, argv+argc, "-J_Q" , 0.0);
  p.setJG = katana::getCmdOption_bool(argv, argv+argc, "-J_G" , false);
  p.J_G = katana::getCmdOption(argv, argv+argc, "-J_G" , 1.7E9 / tsf);
  p.P_G = katana::getCmdOption(argv, argv+argc, "-P_G" , 0.13); 
  p.injEf = katana::getCmdOption(argv, argv+argc, "-injEf" , 0.77);
  

  //feedback parameters
  p.K = katana::getCmdOption(argv, argv+argc, "-K" , 0.0);
  p.C = katana::getCmdOption(argv, argv+argc, "-C" , 0.0);
  p.tau = katana::getCmdOption(argv, argv+argc, "-tau" , 0.0);


  //compute inhomogenously broadened QD distribution
  init_states(p);
  
  outputfile.open("data/states");
  for(std::size_t k = 0; k < NSAMPLES; k++){
    outputfile << p.QD_eps_GS[k]  << "\t" << p.QD_omega_GS[k] << "\t" << p.QD_eps_ES[k]  << "\t" << p.QD_omega_ES[k] << "\t" << p.QD_inh[k] << "\t" << p.QD_GS_DOS[k] << "\t" << p.QD_ES_DOS[k] << endl;
  }
  outputfile.close();
  
  
  //compute parameters
  compute_parameters(p);
  compute_section_parameters(p);

  double maxT = 0.0;
  for(int k=0; k<=NSEC; k++) if(p.DT[k] > maxT) maxT = p.DT[k];
  double Z = 0;
  double halfT0 = 0;
  double halfDT0 = 0;
  double Zabs = 0;
  for(int k=0; k<NSEC; k++){
    Z += p.z[k];  
    halfT0 += p.T[k];
    if(p.SA[k]==1) Zabs += p.z[k];
  }
  halfDT0 += p.DT[0] + p.DT[NSEC];
  for(int k = 1; k < NSEC; k++) halfDT0 += 2.0*p.DT[k];
  
  cout << "device_length: " << Z << " cold-cavity rounttrip time: " << halfDT0 << " absorber length: " << Zabs << endl;
//   for(int k=0; k<NSEC; k++){
//     cout << k << "\t" << p.z[k] <<  "\t" << p.T[k] <<  "\t" << p.SA[k] << endl;
//   }
    

  
  //dynamical variables with delay
  double max_delay = katana::getCmdOption(argv, argv+argc, "-tau_max" , p.tau) + maxT + 20.0*dt;
  vars_vec Xhist(max_delay,dt);
  cout << "History initilized with maxT / dt: " << (maxT)/dt << endl;

  
  //////////////////////////////////////////
  //noise
  //////////////////////////////////////////
  
  //noise
  p.noiseStr = katana::getCmdOption(argv, argv+argc, "-noiseStr" , 1E-10);
  p.SqrtNoiseStr = sqrt( katana::getCmdOption(argv, argv+argc, "-SqrtNoiseStr" , p.noiseStr) * dt );
  
  
  std::function<void (vars*, vars_vec*, parameters*)> noise = [&](vars* X, vars_vec* Xhist, parameters* p){ 
    for(std::size_t k = 0;k<NSEC;k++){
      for(unsigned int l = 0;l<NSAMPLES;l++){
        X->G_plus[k][l] += p->SqrtNoiseStr*katana::gwNoise();
        X->G_minus[k][l] += p->SqrtNoiseStr*katana::gwNoise();
      }
    }
  };
  
  std::function<void (vars*, vars_vec*, parameters*)> uniNoise = [&](vars* X, vars_vec* Xhist, parameters* p){ 
    std::complex<double> uniformNoise = p->SqrtNoiseStr*katana::gwNoise();
    for(std::size_t k = 0;k<NSEC;k++){
      for(unsigned int l = 0;l<NSAMPLES;l++){
        X->G_plus[k][l] += uniformNoise; 
        X->G_minus[k][l] += uniformNoise; 
      }
    }
  };
    
  //gw noise for each section (same for + and -)
  if(katana::getCmdOption_bool(argv, argv+argc, "-uniNoise" , false)) noise = uniNoise;
  if(katana::getCmdOption_bool(argv, argv+argc, "-noNoise" , false)) noise = after_step_empty;
  
  
  //////////////////////////////////////////
  //outputfunctions
  //////////////////////////////////////////
  
  
  //conversion to Watts
  const double CtoWatts = 2.0 * p.e_ph * p.w_base * p.N_QD * 1e12;
//   cout << CtoWatts << endl;
  
  //for timeseries output
  outputfile.precision(10); //outputfile must be declared in global.cpp
  unsigned int out_k = 0;
  unsigned int out_noutk = (int)katana::getCmdOption(argv, argv+argc, "-nout" , 1);
  double out_dt = dt * (double)out_noutk;
  
  
  
  //output timeseries of time and all dynamical variables -> very slow and might take lots  of memory
  auto outputAllToFile = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn && out_k == 0){
      outputfile << tn*dt << "\t";
      for(std::size_t k = 0; k < NSEC; k++) outputfile << x->N[k] << "\t";
//       for(std::size_t k = 0; k < NSEC; k++) outputfile << x->rho_ES[k] << "\t";
//       for(std::size_t k = 0; k < NSEC; k++) outputfile << x->rho_GS[k] << "\t";
      for(std::size_t l = 0; l < NSAMPLES; l++) outputfile << x->rho_GS[NSEC-1][l] << "\t";
      for(std::size_t l = 0; l < NSAMPLES; l++) outputfile << x->rho_ES[NSEC-1][l] << "\t";
      for(std::size_t k = 0; k < NSEC; k++) outputfile << norm(x->A_plus[k])*CtoWatts << "\t";
      for(std::size_t k = 0; k < NSEC; k++) outputfile << norm(x->A_minus[k])*CtoWatts << "\t";
      outputfile << std::endl;
    }
    out_k = (out_k+1)%out_noutk;
  };
  
  
  //output to vectors for TS analysis
  std::vector<double> outDataVec;
  std::vector<std::complex<double>> outDataVecComplex;
//   outDataVec.reserve((int)(outTime/dt));
//   outDataVecComplex.reserve((int)(outTime/dt));
  auto outputToVector = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn && out_k == 0){
      outDataVec.push_back(norm(x->A_out)*CtoWatts);
      outDataVecComplex.push_back(x->A_out);
    }
    out_k = (out_k+1)%out_noutk;
  };
  
  //////////////////////////////////////////
  //for sweeps
  //////////////////////////////////////////
  
  std::ofstream out_sweep_ML_State, out_sweep_Ext, out_sweep_Ext_all, out_sweep_DeltaT, out_sweep_TS;
  out_sweep_ML_State.precision(10), out_sweep_Ext_all.precision(10), out_sweep_TS.precision(10), out_sweep_Ext.precision(10), out_sweep_DeltaT.precision(10);
  
  evalMLTS outEv_sweep;
  outEv_sweep.PW_rel_searchRadius = T_maxPW;
  outEv_sweep.doubleCountTol = T_dCountTol;
  outEv_sweep.pulse_rel_max_th = T_relMax;
  outEv_sweep.TVPD_rel_pulse_trig = T_TVPD_trig;
  outEv_sweep.TVPD_rel_pulse_bounds = T_TVPD_bound;
  outEv_sweep.TSSmoothing = 1;
  
  double* PtrToSweepPar;
  std::string str_sweep_parameter;
  double* PtrToSecPar;
  std::string str_sweep_sec_parameter;
  
  double sweep_IntTime = katana::getCmdOption(argv, argv+argc, "-sIntTime" , 15000.0);
  double sweep_OutTime = katana::getCmdOption(argv, argv+argc, "-sOutTime" , 5000.0);
  int nsweep_steps = katana::getCmdOption(argv, argv+argc, "-sSteps" , 100);
  double sweep_start = katana::getCmdOption(argv, argv+argc, "-sStart" , 0);
  double sweep_end = katana::getCmdOption(argv, argv+argc, "-sEnd" , 1);
  
  int sweep_nAllExt_out = 2000;
  int sweep_nDeltaT_out = 2000;
  
  std::string str_MLL_state;
  std::string str_extrema;
  std::string str_extrema_all;
  std::string str_DeltaT;
  std::string str_powerSpec;
  std::string str_opticalSpec;
  std::string str_AC;
  std::string str_IntAC;
  std::string str_TS;
  std::string str_bin;
  
  std::string str_sweep_data = "data/";
  std::string str_sweep_updown = "up";
  
  bool bool_powerSpec = katana::getCmdOption_bool(argv, argv+argc, "-wpowerSpec" , false);
  bool bool_opticalSpec = katana::getCmdOption_bool(argv, argv+argc, "-wopticalSpec" , false);
  bool bool_AC = katana::getCmdOption_bool(argv, argv+argc, "-wAC" , false);
  bool bool_IntAC = katana::getCmdOption_bool(argv, argv+argc, "-wIntAC" , false);
  bool bool_TS = katana::getCmdOption_bool(argv, argv+argc, "-wTS" , false);
  bool bool_bin = katana::getCmdOption_bool(argv, argv+argc, "-wBin" , false);
  bool bool_Extrema = katana::getCmdOption_bool(argv, argv+argc, "-wExtrema" , false);
  bool bool_Extrema_all = katana::getCmdOption_bool(argv, argv+argc, "-wExtrema_all" , false);
  bool bool_DeltaT = katana::getCmdOption_bool(argv, argv+argc, "-wDeltaT" , false);
  bool bool_NoSweep = katana::getCmdOption_bool(argv, argv+argc, "-wNoSweep" , false);
  bool bool_LogSweep = katana::getCmdOption_bool(argv, argv+argc, "-wLogSweep" , false);
  
  
  auto sweepMLL = [&](){
    
    
    std::vector<double> SweepVals(nsweep_steps+1);
    SweepVals[0] = sweep_start;

    if(bool_LogSweep){
      double SweepLogMult = pow( (sweep_end/sweep_start), 1.0/((double)(nsweep_steps)) );
      for(unsigned int k = 1; k < SweepVals.size(); k++){
        SweepVals[k] = SweepVals[k-1]*SweepLogMult;
//         cout << k << " " <<  SweepVals[k] << endl;
      }
    }
    else{
      double SweepLinAdd = (sweep_end-sweep_start)/(double)nsweep_steps;
      for(unsigned int k = 1; k < SweepVals.size(); k++){
        SweepVals[k] = SweepVals[k-1]+SweepLinAdd;
//         cout << k << " " <<  SweepVals[k] << endl;
      }
    }
    
    
    std::ostringstream ostringstream_DoubleToStr;
    ostringstream_DoubleToStr << std::fixed << std::setprecision(8);
//     ostringstream_DoubleToStr << std::scientific << std::setprecision(8);
    ostringstream_DoubleToStr << *PtrToSecPar;
    std::string StrSecParVal = ostringstream_DoubleToStr.str();
    ostringstream_DoubleToStr.str("");
    ostringstream_DoubleToStr.clear();

    
    str_MLL_state = str_sweep_data + str_sweep_updown + "_sweep_MLL_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_extrema = str_sweep_data + str_sweep_updown + "_sweep_extrema_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_powerSpec = str_sweep_data + "powerSpec/powerSpec_"+str_sweep_updown+"_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_opticalSpec = str_sweep_data + "opticalSpec/opticalSpec_"+str_sweep_updown+"_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_AC = str_sweep_data + "AC/AC_" + str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_IntAC = str_sweep_data + "IntAC/IntAC_"+str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_TS = str_sweep_data + "TS/TS_"+str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_bin = str_sweep_data + "bin/bin_"+str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_extrema_all = str_sweep_data + str_sweep_updown + "_sweep_all_extrema_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_DeltaT = str_sweep_data + str_sweep_updown + "_sweep_DeltaT_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    
    outTime = sweep_OutTime;
    outTime_ntn = (unsigned long long int)(outTime/dt);
    
    outDataVec.reserve(outTime_ntn);
    outDataVecComplex.reserve(outTime_ntn);
    
    out_sweep_ML_State.open(str_MLL_state);
    if(bool_Extrema) out_sweep_Ext.open(str_extrema);
    if(bool_Extrema_all) out_sweep_Ext_all.open(str_extrema_all);
    if(bool_DeltaT) out_sweep_DeltaT.open(str_DeltaT);

    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(histFile);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-6);
    }
    
    for(unsigned int k = 0; k < SweepVals.size(); k++){
    *PtrToSweepPar = SweepVals[k];  
    
    ostringstream_DoubleToStr << std::fixed << std::setprecision(8);
    ostringstream_DoubleToStr << *PtrToSweepPar;
    std::string StrSweepParVal = ostringstream_DoubleToStr.str();
    ostringstream_DoubleToStr.str("");
    ostringstream_DoubleToStr.clear();
    
      
      sm::clck = clock();
      if(bool_NoSweep) {
        //read history file
        if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
          std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
          Xhist.load(histFile);
          tn = Xhist.tn;
        }
        else{
          tn=0;
          Xhist.setToCnst(1E-6);
        }
      }
        
      tn=0;
      outDataVec.resize(0), outDataVecComplex.resize(0);
      compute_parameters(p);
      compute_section_parameters(p);
      DDEsolver.DDE_RK4(MLL_system_eqs, noise, &Xhist, &p, &tn, sweep_IntTime, dt, outputToVector);
      
      outEv_sweep.loadNewTS(&outDataVec, out_dt, p.T0);
      
      
      cout << str_sweep_parameter << ": " << *PtrToSweepPar;
      out_sweep_ML_State << *PtrToSweepPar << "\t" << *PtrToSecPar << "\t"; // 2
      //basics
      out_sweep_ML_State << outEv_sweep.TSAverage << "\t" << outEv_sweep.TSSmallestMin << "\t" << outEv_sweep.TSGreatestMax << "\t"; // 5
      //from max detection:
      out_sweep_ML_State << outEv_sweep.MLnMax() << "\t" << outEv_sweep.meanMax <<  "\t" << outEv_sweep.meanMaxTSep << "\t"; // 8
      out_sweep_ML_State << outEv_sweep.PW_mean_PFWHM << "\t" << outEv_sweep.PW_mean_PSD << "\t"; // 10
      out_sweep_ML_State << outEv_sweep.PtPJ_MaxPos() << "\t" << outEv_sweep.RelAmpJitter() << "\t"; // 12
      //from IntAC:
      out_sweep_ML_State << outEv_sweep.ML_Fund_RR_from_IntAC() << "\t" << outEv_sweep.ML_FWHM_IntAC() << "\t" << outEv_sweep.ML_StdDev_IntAC() << "\t"; //15
      //from TriggerValPulseDetection_Statistics
      out_sweep_ML_State << outEv_sweep.TVPD_meanPulseMax << "\t" << outEv_sweep.TVPD_meanPulseArea << "\t" << outEv_sweep.TVPD_meanPulseSD << "\t" << outEv_sweep.TVPD_meanPulseSep << "\t"; // 19
      out_sweep_ML_State << sqrt(outEv_sweep.Variance(outEv_sweep.TVPD_PM_vec))/outEv_sweep.TVPD_meanPulseMax << "\t" << sqrt(outEv_sweep.Variance(outEv_sweep.TVPD_PA_vec))/outEv_sweep.TVPD_meanPulseArea << "\t" << outEv_sweep.PtPJitter(outEv_sweep.TVPD_PS_vec) << "\t"; //22
      // changing triggers and thresholds for evaluation
      double ev_TSep, ev_TSsd, ev_PM, ev_PMsd, ev_PA, ev_PAsd;
      tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv_sweep.TriggerValPulseDetection_Statistics_forEval(0.05, 0.002);
      out_sweep_ML_State << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //28
      tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv_sweep.TriggerValPulseDetection_Statistics_forEval(0.05, 0.01);
      out_sweep_ML_State << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //34
      tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv_sweep.TriggerValPulseDetection_Statistics_forEval(0.10, 0.02);
      out_sweep_ML_State << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //40
      //IntAC Eval
      double ev_fMa, ev_fMap, ev_fMi, ev_fMip;
      tie(ev_fMa, ev_fMap, ev_fMi, ev_fMip) = outEv_sweep.EvalIntAC();
      out_sweep_ML_State << ev_fMa << "\t" << ev_fMap << "\t" << ev_fMi << "\t" << ev_fMip << "\t"; //44
      //new line
      out_sweep_ML_State << std::endl;
      
      
      if(bool_Extrema){
        out_sweep_Ext << *PtrToSweepPar << "\t" << outEv_sweep.TSGreatestMax << std::endl;
        for(std::size_t k = 0; k < outEv_sweep.uniqueMax.size(); k++) out_sweep_Ext << *PtrToSweepPar << "\t" << outEv_sweep.uniqueMax[k] << std::endl;      
        out_sweep_Ext << *PtrToSweepPar << "\t" <<  outEv_sweep.TSSmallestMin << std::endl;
        for(std::size_t k = 0; k < outEv_sweep.uniqueMin.size(); k++) out_sweep_Ext << *PtrToSweepPar << "\t" << outEv_sweep.uniqueMin[k] << std::endl;
      }
      
      if(bool_Extrema_all){
        out_sweep_Ext_all << *PtrToSweepPar << "\t" << outEv_sweep.TSGreatestMax << std::endl;
        for(int k = 0; k < std::min((int)outEv_sweep.maxima.size(),sweep_nAllExt_out); k++) out_sweep_Ext_all << *PtrToSweepPar << "\t" << outEv_sweep.maxima[k] << std::endl;      
        out_sweep_Ext_all << *PtrToSweepPar << "\t" <<  outEv_sweep.TSSmallestMin << std::endl;
        for(int k = 0; k < std::min((int)outEv_sweep.uniqueMin.size(),sweep_nAllExt_out); k++) out_sweep_Ext_all << *PtrToSweepPar << "\t" << outEv_sweep.uniqueMin[k] << std::endl;
      }
    
      if(bool_DeltaT){
        for(int k = 0; k < std::min((int)((int)outEv_sweep.maxima_tvec.size()-1),sweep_nDeltaT_out); k++) out_sweep_DeltaT << *PtrToSweepPar << "\t" << (outEv_sweep.maxima_tvec[k+1] - outEv_sweep.maxima_tvec[k]) << std::endl;
      }
      
      if(bool_TS){
        out_sweep_TS.open(str_TS+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix);
        for(std::size_t k = 0; k < outDataVec.size()/10; k++) out_sweep_TS << k*out_dt << "\t" << outDataVec[k] << std::endl;
        out_sweep_TS.close();
      }
    
      if(bool_IntAC){
        katana::sm_print_correlation(outEv_sweep.IntAC, out_dt, str_IntAC+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix, 100.0);
      }
      
      if(bool_AC){
        std::vector<double> AC;
        katana::get_autocorr_norm_out(outDataVecComplex, AC);
        katana::sm_print_correlation(AC, out_dt, str_AC+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix, 100);
      }
      
      if(bool_powerSpec){
        std::vector<double> powerSpec;
        std::vector<double> WindowedData;
        sm::window_Hann(outDataVec, WindowedData);
        sm::get_power_spec(WindowedData, powerSpec, out_dt);
        sm::dump_power_spec(powerSpec, out_dt*1E-3, 50.0, str_powerSpec+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix);
      }
      
      if(bool_opticalSpec){
        std::vector<double> opticalSpec;
        sm::get_optical_spec(outDataVecComplex, opticalSpec);
        sm::dump_optical_spec(opticalSpec, out_dt*1E-3, 1000.0, str_opticalSpec+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix);
      }
      
      if(bool_bin){
        Xhist.tn = tn;
        Xhist.save(str_bin+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix);
      }

      std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
      std::time_t now_p = std::chrono::system_clock::to_time_t(now);
      cout << ",\tintegration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << ",\tcurrent time and date: " << std::put_time(std::localtime(&now_p), "%F %T") << endl;
      
    }
    out_sweep_Ext.close();
    out_sweep_ML_State.close();
    if(bool_Extrema_all) out_sweep_Ext_all.close();
    if(bool_DeltaT) out_sweep_DeltaT.close();
    
  };
  
  
  
  //////////////////////////////////////////
  //simulations
  //////////////////////////////////////////
  
  

    
      
  //timeseries to file
  if(katana::getCmdOption_bool(argv, argv+argc, "-simpleTS" , false)){
    outputfile.open("data/out_simpleTS");
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(histFile);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-6);
    }
    sm::clck = clock();
    DDEsolver.DDE_RK4(MLL_system_eqs, noise, &Xhist, &p, &tn, intTime, dt, outputAllToFile);
//     DDEsolver.DDE_euler(MLL_system_eqs, noise, &Xhist, &p, &tn, intTime, dt, outputAllToFile);
    cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
    outputfile.close();
    
    //save history to binary file
    if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
      Xhist.tn = tn;
      Xhist.save("data/Xhist.bin");
    }
    
  }
  
  
  //compute and analize time series
  if(katana::getCmdOption_bool(argv, argv+argc, "-TS" , false)){
    evalMLTS outEv;
    outEv.PW_rel_searchRadius = T_maxPW;
    outEv.doubleCountTol = T_dCountTol;
    outEv.pulse_rel_max_th = T_relMax;
    outEv.TVPD_rel_pulse_trig = T_TVPD_trig;
    outEv.TVPD_rel_pulse_bounds = T_TVPD_bound;
    outEv.TSSmoothing = 1;
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(histFile);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-6);
    }
    
    outDataVec.reserve(outTime_ntn);
    outDataVecComplex.reserve(outTime_ntn);

    sm::clck = clock();
    DDEsolver.DDE_RK4(MLL_system_eqs, noise, &Xhist, &p, &tn, intTime, dt, outputToVector);
//     DDEsolver.DDE_euler(MLL_system_eqs, noise, &Xhist, &p, &tn, intTime, dt, outputToVector);
    cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
    
    outEv.loadNewTS(&outDataVec, out_dt, p.T0);
        
    outputfile.open("data/out_TS");
    for(std::size_t k = 0; k < outDataVec.size(); k++){
      outputfile << k*out_dt << "\t" << outDataVec[k] << "\t";
//       outputfile << real(outDataVecComplex[k]) << "\t" << imag(outDataVecComplex[k]) << "\t" << arg(outDataVecComplex[k]) << "\t";
//       outputfile << fb[k] << "\t";
      outputfile << endl;
    }
    outputfile.close();
    
    
//     std::vector<double> Phase,Freq;
//     double avOptFreq = outEv.opticalFreq(outDataVecComplex, Phase, Freq, dt);
//     cout << "avOptFreq: " << avOptFreq << endl;
//     
//     outputfile.open("data/out_TS_phase");
//     for(std::size_t k = 0; k < Freq.size(); k++) outputfile << k*dt << "\t" << Phase[k] << "\t" << Freq[k] << endl;
//     outputfile.close();
   
//     
//     cout << "----input----" << endl;
    cout << "----results----" << endl;
    cout << "average output power: " << outEv.TSAverage << " greatest peak power: " << outEv.TSGreatestMax << " TSSmallestMin: " << outEv.TSSmallestMin << endl;
    cout << "----from maxima detection:" << endl;
    cout << "number of unique ML max: " << outEv.MLnMax() << " mean peak power: " << outEv.meanMax <<  " mean max t seperation: " << outEv.meanMaxTSep << endl;
    cout << "PtP jitter from max pos: " << outEv.PtPJ_MaxPos() << " amplitude jitter from max pos: " << outEv.RelAmpJitter() << endl;
    cout << "PW from pulse FWHM: " << outEv.PW_mean_PFWHM << " PW from pulse StdDev: " << outEv.PW_mean_PSD << endl;
    cout << "----from the intensity AC:" << endl;
    cout << "ML fundamental RR from IntAC: " << outEv.ML_Fund_RR_from_IntAC() << endl;
    cout << "IntAC FWHM: " << outEv.ML_FWHM_IntAC() << " IntAC StdDev: " << outEv.ML_StdDev_IntAC() << endl;
    cout << "----from weighted pulse positions:" << endl;
    cout << "mean TVPD peak power: " << outEv.TVPD_meanPulseMax << " mean TVPD pulse area: " << outEv.TVPD_meanPulseArea << " TVPD PW from PStdDev: " << outEv.TVPD_meanPulseSD << " TVPD mean pulse sep: " << outEv.TVPD_meanPulseSep << endl;
    cout << "TVPD peak power rel jitter: " << sqrt(outEv.Variance(outEv.TVPD_PM_vec))/outEv.TVPD_meanPulseMax << " TVPD pulse area rel jitter: " << sqrt(outEv.Variance(outEv.TVPD_PA_vec))/outEv.TVPD_meanPulseArea << endl;
    cout << "TVPD pulse position jitter: " << outEv.PtPJitter(outEv.TVPD_PS_vec) << endl;
    double ev_TSep, ev_TSsd, ev_PM, ev_PMsd, ev_PA, ev_PAsd;
    tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv.TriggerValPulseDetection_Statistics_forEval(0.005, 0.001);
    cout << ev_TSep << " " << ev_TSsd << " " << ev_PM << " " << ev_PMsd << " " << ev_PA << " " << ev_PAsd << endl;
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-powerSpec" , false)){
      std::vector<double> powerSpec;
      std::vector<double> WindowedData;
      
      sm::get_power_spec(outDataVec, powerSpec, out_dt);
      sm::dump_power_spec(powerSpec, out_dt*1E-3, 100, "data/out_TS_powerSpec");
      
      sm::window_Hann(outDataVec, WindowedData);
      sm::get_power_spec(WindowedData, powerSpec, out_dt);
      sm::dump_power_spec(powerSpec, out_dt*1E-3, 100, "data/out_TS_powerSpec_Hann");
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-opticalSpec" , false)){
      std::vector<double> opticalSpec;
      std::vector<std::complex<double>> WindowedData;
      
      sm::get_optical_spec(outDataVecComplex, opticalSpec);
      sm::dump_optical_spec(opticalSpec, out_dt*1E-3, 4000, "data/out_TS_opticalSpec"+str_suffix);
      
      sm::window_Hann(outDataVecComplex, WindowedData);
      sm::get_optical_spec(WindowedData, opticalSpec);
      sm::dump_optical_spec(opticalSpec, out_dt*1E-3, 4000, "data/out_TS_opticalSpec_Hann"+str_suffix);
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-AC" , false)){
      double ACout = katana::getCmdOption(argv, argv+argc, "-ACout" , 100);
      std::vector<double> AC;
      katana::get_autocorr_norm_out(outDataVecComplex, AC);
      katana::sm_print_correlation(AC, out_dt, "data/out_TS_AC", ACout);
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-IntAC" , false)){
      katana::sm_print_correlation(outEv.IntAC, out_dt, "data/out_TS_IntAC", 100);
    }
    
    //save history to binary file
    if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
      Xhist.tn = tn;
      std::string histFile = katana::getCmdOption(argv, argv+argc, "-shistFile" , "data/Xhist.bin");
      Xhist.save(histFile);
    }
    
  }
  
  
  //calc fitness
  if(katana::getCmdOption_bool(argv, argv+argc, "-fitness" , false, bquiet)){
    evalMLTS outEv;
    outEv.PW_rel_searchRadius = T_maxPW;
    outEv.doubleCountTol = T_dCountTol;
    outEv.pulse_rel_max_th = T_relMax;
    outEv.TVPD_rel_pulse_trig = T_TVPD_trig;
    outEv.TVPD_rel_pulse_bounds = T_TVPD_bound;
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false, bquiet)){
      std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(histFile);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-6);
    }
    
    DDEsolver.DDE_RK4(MLL_system_eqs, noise, &Xhist, &p, &tn, intTime, dt, outputToVector);
    
    outEv.loadNewTS(&outDataVec, out_dt, p.T0);
    
    
    double meanPower = outEv.TSAverage;
    double meanPeakPower = outEv.TVPD_meanPulseMax;
    double pulseWidth_FWHM = outEv.ML_FWHM_IntAC();
    double pulseWidth_std = outEv.ML_StdDev_IntAC();
    double meanPulseSep = outEv.TVPD_meanPulseSep;
    
    double fitness = 0.0;
    
    if( meanPower > 0.0 && pulseWidth_std > 0.0 && pulseWidth_FWHM > 0.0 && meanPulseSep > 0.8*p.T0 ){
      fitness = meanPower / ( p.P_G * (pulseWidth_FWHM + pulseWidth_std) * (pulseWidth_FWHM + pulseWidth_std) );
    }
    
    cout << fitness << endl;
    
    
    std::string filename = katana::getCmdOption(argv, argv+argc, "-filename" , "data/fitness");
    
    outputfile.open(filename);
    outputfile << "fitness: " << fitness << endl;
    outputfile << "pump: " << p.P_G << endl;
    outputfile << "mean_power: " << meanPower << endl;
    outputfile << "mean_peak_power: " << meanPeakPower << endl;
    outputfile << "pulse_width_FWHM: " << pulseWidth_FWHM << endl;
    outputfile << "pulse_width_std: " << pulseWidth_std << endl;
    outputfile.close();
    
  }
  
  
  //   //high resolution power spectrum
//   if(katana::getCmdOption_bool(argv, argv+argc, "-HRPowerSpec" , false)){
//     
//     //read history file
//     if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
//       std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
//       Xhist.load(histFile);
//       tn = Xhist.tn;
//     }
//     else{
//       tn=0;
//       Xhist.setToCnst(1E-6);
//     }
//     
//     const double nsigmas = 10;
//     double SampleDt = katana::getCmdOption(argv, argv+argc, "-SampleDt" , 1.0);
//     unsigned int SampleDt_dtsteps = (int)(SampleDt/dt);
//     unsigned int ConvKernelSize = (int)(nsigmas*SampleDt/dt);
//     
//     intTime += nsigmas*SampleDt;
//     outTime_ntn += ConvKernelSize;
//     
//     cout << "SampleDt: " << SampleDt << ", SampleDt / dt: " << SampleDt_dtsteps << ", max frequency 1000/(2*SampleDt): " << 1000.0/(2.0*SampleDt) << endl;
//     double maxOutFreq = katana::getCmdOption(argv, argv+argc, "-maxOutFreq" , 280.0);
//     
//     std::vector<double> ConvKernel(ConvKernelSize);
//     if(katana::getCmdOption_bool(argv, argv+argc, "-SmplKrnl" , false)){
//       for(unsigned int k = 0; k<ConvKernelSize; k++){
//         if(k==0) ConvKernel[k] = 1.0;
//         else ConvKernel[k] = 0.0;
//       }
//     }
//     else{
//       double ConvNorm = 0.0;
//       for(unsigned int k = 0; k<ConvKernelSize; k++){
//           double v = exp(-0.5*(k-ConvKernelSize/2.0)*(k-ConvKernelSize/2.0)/((double)(0.5*SampleDt_dtsteps*0.5*SampleDt_dtsteps)));
//     //       cout << v << endl;
//           ConvKernel[k] = v;
//           ConvNorm += v;
//       }
//       for(unsigned int k = 0; k<ConvKernelSize; k++){
//         ConvKernel[k] = ConvKernel[k]/ConvNorm;
//       }
//     }
//     
// //     outputfile.open("data/test");
// //     for(std::size_t k = 0; k < ConvKernel.size(); k++){
// //       outputfile << k << "\t" << ConvKernel[k] << "\t";
// //       outputfile << endl;
// //     }
// //     outputfile.close();
//     
//     std::vector<double> ShortHist(ConvKernelSize);
//     unsigned int SHind = 0;
//     
//     //output intensity to vector 
//     std::vector<double> ConvDataVec;
//     ConvDataVec.reserve((int)(outTime/SampleDt));
//     cout << "reseverd ConvDataVec size = 8byte * " << (int)(outTime/SampleDt) << " = " << 8e-6*outTime/SampleDt << " Mbyte" << endl;
//     auto outputToVectorConv = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
//       if(tn >= tn_final - outTime_ntn){
//         ShortHist[SHind] = norm(x->A_out);
//         SHind = (SHind+1)%ConvKernelSize;
//         
//         if(tn >= tn_final - outTime_ntn + ConvKernelSize){
//           if( (SHind%SampleDt_dtsteps) == 0 ){
//             double ConvVal = 0.0;
//             for(unsigned int k = 0; k<ConvKernelSize; k++){
//               ConvVal += ShortHist[(k+SHind)%ConvKernelSize]*ConvKernel[k];
//             }
//             ConvDataVec.push_back(ConvVal);
//           }
//         }
//       }
//     };
//     
//     
//     if(katana::getCmdOption_bool(argv, argv+argc, "-HRS_buffer" , false)){
//       double HRS_bufferIntTime = katana::getCmdOption(argv, argv+argc, "-LTTJ_BIntTime" , 50000);
//       cout << "integrate buffer time " << HRS_bufferIntTime << endl;
//       DDEsolver.DDE_RK4(MLL_system_eqs, noise, &Xhist, &p, &tn, HRS_bufferIntTime, dt, outputfunc_empty);
//     }
//     
//     sm::clck = clock();
//     DDEsolver.DDE_RK4(MLL_system_eqs, noise, &Xhist, &p, &tn, intTime, dt, outputToVectorConv);
//     cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
// 
//     cout << "ConvDataVec capacity after computation: " << 8e-6*ConvDataVec.capacity() << " Mbyte" << endl;
//     cout << "ConvKernel capacity after computation: " << 8e-6*ConvKernel.capacity() << " Mbyte" << endl;
//     cout << "ShortHist capacity after computation: " << 8e-6*ShortHist.capacity() << " Mbyte" << endl;
//     
//     if(katana::getCmdOption_bool(argv, argv+argc, "-HRPowerSpec_Hann" , false)){
//       sm::window_Hann(ConvDataVec);
//     }
//     
// //     std::vector<double> powerSpec;
// //     sm::get_power_spec(ConvDataVec, powerSpec, SampleDt);
// //     sm::dump_power_spec(powerSpec, SampleDt*1E-3, maxOutFreq, "data/ConvTS_powerSpec"+str_suffix);
// //     cout << "powerSpec capacity after computation: " << 8e-6*powerSpec.capacity() << " Mbyte" << endl;
// 
//     sm::get_power_spec(ConvDataVec, SampleDt);
//     sm::dump_power_spec(ConvDataVec, SampleDt*1E-3, maxOutFreq, "data/ConvTS_powerSpec"+str_suffix);
//     
//     
// //     outputfile.open("data/out_ConvTS");
// //     for(std::size_t k = 0; k < ConvDataVec.size(); k++){
// //       outputfile << k*SampleDt << "\t" << ConvDataVec[k] << "\t";
// //       outputfile << endl;
// //     }
// //     outputfile.close();
//     
//     //save history to binary file
//     if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
//       Xhist.tn = tn;
//       Xhist.save("data/Xhist.bin");
//     }
//     
//   }
  
  
  
  
//   //position resolved time averaged dynamics
//   if(katana::getCmdOption_bool(argv, argv+argc, "-posTAvDyn" , false)){
//     evalMLTS outEv;
//     outEv.PW_rel_searchRadius = T_maxPW;
//     outEv.doubleCountTol = T_dCountTol;
//     outEv.pulse_rel_max_th = T_relMax;
//     outEv.TVPD_rel_pulse_trig = T_TVPD_trig;
//     outEv.TVPD_rel_pulse_bounds = T_TVPD_bound;
//     outEv.TSSmoothing = 1;
//     
//     
//     //output matrix
//     std::vector<std::vector<double>> outVecAllVars(2*NSEC,std::vector<double>((unsigned int)outTime/dt));
//     for(std::size_t k = 0; k < 2*NSEC; k++) outVecAllVars[k].resize(0);
//     
//     //output all intensities to matrix
//     auto outputToMatrix = [&](vars* X, unsigned long long int tn, unsigned long long int tn_final){
//       if(tn >= tn_final - outTime_ntn){
//         for(std::size_t k = 0; k < NSEC; k++) outVecAllVars[k].push_back(norm(X->A_plus[k]));
//         for(std::size_t k = NSEC; k < 2*NSEC; k++) outVecAllVars[k].push_back(norm(X->A_minus[k-NSEC]));
//       }
//     };
//     
//     //read history file
//     if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
//       std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
//       Xhist.load(histFile);
//       tn = Xhist.tn;
//     }
//     else{
//       tn=0;
//       Xhist.setToCnst(1E-6);
//     }
//     
//     sm::clck = clock();
//     DDEsolver.DDE_RK4(MLL_system_eqs, noise, &Xhist, &p, &tn, intTime, dt, outputToMatrix);
//     cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
// 
//     
// //     outputfile.open("data/posTAvDyn/out_TS"+str_suffix);
// //     for(std::size_t k = 0; k < outVecAllVars[0].size(); k++) outputfile << k*dt << "\t" << outVecAllVars[NSEC-1][k] << endl;
// //     outputfile.close();
//     
//     std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
//     std::time_t now_p = std::chrono::system_clock::to_time_t(now);
//     cout << "anayzing time series. Current time and date: " << std::put_time(std::localtime(&now_p), "%F %T") << endl;
//     
//     std::ofstream outExt,outPosTAvState,outPulsePositions;
//     outExt.precision(10),outPosTAvState.precision(10),outPulsePositions.precision(15);
//     
// //     outExt.open("data/posTAvDyn/posTAvDyn_extrema_R"+str_suffix);
//     outPosTAvState.open("data/posTAvDyn/posTAvDyn_State_R"+str_suffix);
//     
//     std::vector<std::vector<double>> PulsePositions (2*NSEC);
//     
//     
//     for(std::size_t k = 0; k <NSEC; k++){
//       outEv.loadNewTS(&(outVecAllVars[k]), dt, p.T0);
//       
//       outPosTAvState << k << "\t" << 1 << "\t"; // 2
//       //basics
//       outPosTAvState << outEv.TSAverage << "\t" << outEv.TSSmallestMin << "\t" << outEv.TSGreatestMax << "\t"; // 5
//       //from max detection:
//       outPosTAvState << outEv.MLnMax() << "\t" << outEv.meanMax <<  "\t" << outEv.meanMaxTSep << "\t"; // 8
//       outPosTAvState << outEv.PW_mean_PFWHM << "\t" << outEv.PW_mean_PSD << "\t"; // 10
//       outPosTAvState << outEv.PtPJ_MaxPos() << "\t" << outEv.RelAmpJitter() << "\t"; // 12
//       //from IntAC:
//       outPosTAvState << outEv.ML_Fund_RR_from_IntAC() << "\t" << outEv.ML_FWHM_IntAC() << "\t" << outEv.ML_StdDev_IntAC() << "\t"; //15
//       //from TriggerValPulseDetection_Statistics
//       outPosTAvState << outEv.TVPD_meanPulseMax << "\t" << outEv.TVPD_meanPulseArea << "\t" << outEv.TVPD_meanPulseSD << "\t" << outEv.TVPD_meanPulseSep << "\t"; // 19
//       outPosTAvState << sqrt(outEv.Variance(outEv.TVPD_PM_vec))/outEv.TVPD_meanPulseMax << "\t" << sqrt(outEv.Variance(outEv.TVPD_PA_vec))/outEv.TVPD_meanPulseArea << "\t" << outEv.PtPJitter(outEv.TVPD_PS_vec) << "\t"; //22
//       // changing triggers and thresholds for evaluation
//       double ev_TSep, ev_TSsd, ev_PM, ev_PMsd, ev_PA, ev_PAsd;
//       tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv.TriggerValPulseDetection_Statistics_forEval(0.05, 0.002);
//       outPosTAvState << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //28
//       tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv.TriggerValPulseDetection_Statistics_forEval(0.05, 0.01);
//       outPosTAvState << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //34
//       tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv.TriggerValPulseDetection_Statistics_forEval(0.10, 0.02);
//       outPosTAvState << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //40
//       //IntAC Eval
//       double ev_fMa, ev_fMap, ev_fMi, ev_fMip;
//       tie(ev_fMa, ev_fMap, ev_fMi, ev_fMip) = outEv.EvalIntAC();
//       outPosTAvState << ev_fMa << "\t" << ev_fMap << "\t" << ev_fMi << "\t" << ev_fMip << "\t"; //44
//       //new line
//       outPosTAvState << std::endl;
//       
//       PulsePositions[k] = outEv.TVPD_PS_vec;
// //       for(std::vector<double>::iterator it = outEv.maxima.begin(); it != outEv.maxima.end(); ++it) outExt << k-1 << "\t" << *it << endl;
//     }
// //     outExt.close();
//     outPosTAvState.close();
//     
// //     outExt.open("data/posTAvDyn/posTAvDyn_extrema_L"+str_suffix);
//     outPosTAvState.open("data/posTAvDyn/posTAvDyn_State_L"+str_suffix);
// 
//     for(std::size_t k = NSEC; k < 2*NSEC; k++){
//       outEv.loadNewTS(&(outVecAllVars[k]), dt, p.T0);
//       
//       outPosTAvState << k-NSEC << "\t" << -1 << "\t"; // 2
//       //basics
//       outPosTAvState << outEv.TSAverage << "\t" << outEv.TSSmallestMin << "\t" << outEv.TSGreatestMax << "\t"; // 5
//       //from max detection:
//       outPosTAvState << outEv.MLnMax() << "\t" << outEv.meanMax <<  "\t" << outEv.meanMaxTSep << "\t"; // 8
//       outPosTAvState << outEv.PW_mean_PFWHM << "\t" << outEv.PW_mean_PSD << "\t"; // 10
//       outPosTAvState << outEv.PtPJ_MaxPos() << "\t" << outEv.RelAmpJitter() << "\t"; // 12
//       //from IntAC:
//       outPosTAvState << outEv.ML_Fund_RR_from_IntAC() << "\t" << outEv.ML_FWHM_IntAC() << "\t" << outEv.ML_StdDev_IntAC() << "\t"; //15
//       //from TriggerValPulseDetection_Statistics
//       outPosTAvState << outEv.TVPD_meanPulseMax << "\t" << outEv.TVPD_meanPulseArea << "\t" << outEv.TVPD_meanPulseSD << "\t" << outEv.TVPD_meanPulseSep << "\t"; // 19
//       outPosTAvState << sqrt(outEv.Variance(outEv.TVPD_PM_vec))/outEv.TVPD_meanPulseMax << "\t" << sqrt(outEv.Variance(outEv.TVPD_PA_vec))/outEv.TVPD_meanPulseArea << "\t" << outEv.PtPJitter(outEv.TVPD_PS_vec) << "\t"; //22
//       // changing triggers and thresholds for evaluation
//       double ev_TSep, ev_TSsd, ev_PM, ev_PMsd, ev_PA, ev_PAsd;
//       tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv.TriggerValPulseDetection_Statistics_forEval(0.05, 0.002);
//       outPosTAvState << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //28
//       tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv.TriggerValPulseDetection_Statistics_forEval(0.05, 0.01);
//       outPosTAvState << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //34
//       tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv.TriggerValPulseDetection_Statistics_forEval(0.10, 0.02);
//       outPosTAvState << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //40
//       //IntAC Eval
//       double ev_fMa, ev_fMap, ev_fMi, ev_fMip;
//       tie(ev_fMa, ev_fMap, ev_fMi, ev_fMip) = outEv.EvalIntAC();
//       outPosTAvState << ev_fMa << "\t" << ev_fMap << "\t" << ev_fMi << "\t" << ev_fMip << "\t"; //44
//       //new line
//       outPosTAvState << std::endl;
//       
//       PulsePositions[k] = outEv.TVPD_PS_vec;
// //       for(std::vector<double>::iterator it = outEv.maxima.begin(); it != outEv.maxima.end(); ++it) outExt << k-1-NSEC << "\t" << *it << endl;
//     }
// //     outExt.close();
//     outPosTAvState.close();
//     
//     now = std::chrono::system_clock::now();
//     now_p = std::chrono::system_clock::to_time_t(now);
//     cout << "writing pulse interval times to output file. Current time and date: " << std::put_time(std::localtime(&now_p), "%F %T") << endl;
//     
//     //determine minimal number of pulses in the PulsePositions vectors (all section)
//     std::size_t PulsePositions_maxsize = PulsePositions[0].size();
//     for (std::size_t k = 1; k < 2*NSEC; k++) if(PulsePositions[k].size() < PulsePositions_maxsize) PulsePositions_maxsize = PulsePositions[k].size();
//     
// //     outPulsePositions.open("data/posTAvDyn/posTAvDyn_PulsePositions"+str_suffix);
// //     for (std::size_t m = 0; m < PulsePositions_maxsize; m++){
// //       outPulsePositions << m << "\t";
// //       for (std::size_t k = 0; k < 2*NSEC; k++) outPulsePositions << PulsePositions[k][m] << "\t";
// //       outPulsePositions << endl;
// //     }
// //     outPulsePositions.close();
//     
//     outPulsePositions.open("data/posTAvDyn/posTAvDyn_PulseIntervals"+str_suffix);
//     for (std::size_t m = 1; m < PulsePositions_maxsize; m++){
//       outPulsePositions << m << "\t";
//       for (std::size_t k = 0; k < 2*NSEC; k++) outPulsePositions << PulsePositions[k][m]-PulsePositions[k][m-1] << "\t";
//       outPulsePositions << endl;
//     }
//     outPulsePositions.close();
//     
//     now = std::chrono::system_clock::now();
//     now_p = std::chrono::system_clock::to_time_t(now);
//     cout << "run finished at current time and date: " << std::put_time(std::localtime(&now_p), "%F %T") << endl;
//    
//   }
  
  
//   //position and time resolved gain and phase dynamics
//   if(katana::getCmdOption_bool(argv, argv+argc, "-posTDyn" , false)){
//     
//     //dynamical variables with delay
// //     vars_vec XhistLong( ((int)((p.T0*1.1)/dt))*dt ,dt);
//     vars_vec_wdX XhistLong( ((int)((p.T0*1.1)/dt))*dt ,dt);
//     cout << "Long history initilized with T[0] / dt: " << ((int)((p.T0*1.1)/dt)) << endl;
//     XhistLong.dt_divides_tau = true;
//    
//     //read history file
//     if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
//       XhistLong.load("data/XhistLong.bin");
//       tn = XhistLong.tn;
//     }
//     else{
//       tn=0;
//       XhistLong.setToCnst(0.001);
//     }
//     
//     std::ofstream f_netGain ("data/posTDyn/netGain", std::ios::out | std::ios::trunc);
//     std::ofstream f_gain ("data/posTDyn/gain", std::ios::out | std::ios::trunc);
//     std::ofstream f_pulse ("data/posTDyn/pulse", std::ios::out | std::ios::trunc);
//     std::ofstream f_pulseAmp ("data/posTDyn/pulseAmp", std::ios::out | std::ios::trunc);
//     std::ofstream f_phase ("data/posTDyn/phase", std::ios::out | std::ios::trunc);
//     std::ofstream f_phaseChange ("data/posTDyn/phaseChange", std::ios::out | std::ios::trunc);
//     
//     bool onlyNetgain = katana::getCmdOption_bool(argv, argv+argc, "-onlyNetgain" , false);
//     
//     
//     f_netGain.precision(10), f_pulse.precision(10), f_pulseAmp.precision(10), f_phase.precision(10), f_phaseChange.precision(10);
//     
//     auto negatePhaseJump = [](double B, double A){
//       if((B-A) > M_PI) return B-A-2.0*M_PI;
//       if((B-A) < -M_PI) return B-A+2.0*M_PI;
//       else return B-A;
//     };
//     
//     
//     auto outputGainPhase = [&](vars* X, unsigned long long int tn, unsigned long long int tn_final){
//       if(tn >= tn_final - outTime_ntn){
// 	
//         double netGain = 0.0;
//         //including absorber detuning
// //         for(std::size_t k = 0; k < NSEC; k++) netGain += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5) / (1.0 + (imag(p.iDelta_omega[k])*imag(p.iDelta_omega[k]))/(p.gamma_2[k]*p.gamma_2[k])) - 0.5*p.alpha_int[k]);
// //         for(std::size_t k = 0; k < NSEC; k++) netGain += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5) / (1.0 + (imag(p.iDelta_omega[k])*imag(p.iDelta_omega[k]))/(p.gamma_2[k]*p.gamma_2[k])) - 0.5*p.alpha_int[k]);
// //         //excluding absorber detuning
//       	for(std::size_t k = 0; k < NSEC; k++) netGain += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5) - 0.5*p.alpha_int[k]);
//       	for(std::size_t k = 0; k < NSEC; k++) netGain += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5) - 0.5*p.alpha_int[k]);
//         
//         double GainL = 0.0;
//         for(std::size_t k = 0; k < p.nSecG; k++) GainL += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5));
//         for(std::size_t k = 0; k < p.nSecG; k++) GainL += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5));
//         double GainSA = 0.0;
// //         for(std::size_t k = p.nSecG; k < p.nSecG + p.nSecA; k++) GainSA += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5) / (1.0 + (imag(p.iDelta_omega[k])*imag(p.iDelta_omega[k]))/(p.gamma_2[k]*p.gamma_2[k])) );
// //         for(std::size_t k = p.nSecG; k < p.nSecG + p.nSecA; k++) GainSA += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5) / (1.0 + (imag(p.iDelta_omega[k])*imag(p.iDelta_omega[k]))/(p.gamma_2[k]*p.gamma_2[k])) );
//         //excluding absorber detuning
//         for(std::size_t k = p.nSecG; k < p.nSecG + p.nSecA; k++) GainSA += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5) );
//         for(std::size_t k = p.nSecG; k < p.nSecG + p.nSecA; k++) GainSA += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5) );
//         double GainT = 0.0;
//         for(std::size_t k = p.nSecG + p.nSecA; k < NSEC; k++) GainT += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5));
//         for(std::size_t k = p.nSecG + p.nSecA; k < NSEC; k++) GainT += p.z[k] * (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5));
//         double Losses = 0.0;
//       	for(std::size_t k = 0; k < NSEC; k++) Losses += p.z[k] * (-0.5*p.alpha_int[k]);
//       	for(std::size_t k = 0; k < NSEC; k++) Losses += p.z[k] * (-0.5*p.alpha_int[k]);
//         
//         double ImagGain = 0.0;
// //         for(std::size_t k = 0; k < NSEC; k++) ImagGain += p.z[k] * (imag(XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_ES[k] * -p.idelta_omega_ES[k]) + imag(p.iDelta_omega[k] * p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5)/p.gamma_2[k]));
// //         for(std::size_t k = 0; k < NSEC; k++) ImagGain += p.z[k] * (imag(XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_ES[k] * -p.idelta_omega_ES[k]) + imag(p.iDelta_omega[k] * p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5)/p.gamma_2[k]));
//         
//         for(std::size_t k = 0; k < NSEC; k++) ImagGain += p.z[k] * (imag(XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_ES[k] * -p.idelta_omega_ES[k]));
//         for(std::size_t k = 0; k < NSEC; k++) ImagGain += p.z[k] * (imag(XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_ES[k] * -p.idelta_omega_ES[k]));
//         
//         f_netGain << dt*tn << "\t" << norm(X->A_plus[NSEC-1]) << "\t";
//         f_netGain << netGain + log(sqrt(p.r_L*p.r_R)) << "\t" << GainL << "\t" << GainSA << "\t" << GainT << "\t" << Losses + log(sqrt(p.r_L*p.r_R)) << "\t" << ImagGain << "\t";
//         f_netGain << real(X->A_plus[NSEC-1]) << "\t" << imag(X->A_plus[NSEC-1]) << "\t";
//         
//         f_netGain << endl;
//         
//         if(onlyNetgain == false){
//         
//         
//           for(std::size_t k = 0; k < NSEC; k++) f_gain << (p.Gamma_r[k] * p.g[k] * (X->rho_GS[k] - 0.5) / (1.0 + (imag(p.iDelta_omega[k])*imag(p.iDelta_omega[k]))/(p.gamma_2[k]*p.gamma_2[k])) - 0.5*p.alpha_int[k]) << "\t";
//   //       	for(std::size_t k = 0; k < NSEC; k++) f_gain << p.Gamma_r[k] * p.g[k] * (X->rho_GS[k] - 0.5) - 0.5*p.alpha_int[k] << "\t";
//           for(std::size_t k = 0; k < NSEC; k++) f_gain << abs(p.idelta_omega_ES[k] * X->rho_ES[k]) << "\t";
// 
//           
//           //lab frame:::::::
//           
//           //pulse propagation
//           for(std::size_t k = 0; k < NSEC; k++) f_pulse << norm(X->A_plus[k]) << "\t";
//           for(std::size_t k = 0; k < NSEC; k++) f_pulse << norm(X->A_minus[k]) << "\t";
// 
//         // 	//pulse amplification/net gain 
//         // 	f_pulseAmp << norm(XhistLong.at_t0_rPtr()->A_plus[0]) - norm(XhistLong.at_Rnd_rPtr(p.T[0])->A_minus[0]) << "\t";
//         // 	for(std::size_t k = 1; k < NSEC; k++) f_pulseAmp << norm(XhistLong.at_t0_rPtr()->A_plus[k]) - norm(XhistLong.at_Rnd_rPtr(p.T[0])->A_plus[k-1]) << "\t";
//         // 	for(std::size_t k = 0; k < NSEC-1; k++) f_pulseAmp << norm(XhistLong.at_t0_rPtr()->A_minus[k]) - norm(XhistLong.at_Rnd_rPtr(p.T[0])->A_minus[k+1]) << "\t";
//         // 	f_pulseAmp << norm(XhistLong.at_t0_rPtr()->A_minus[NSEC-1]) - norm(XhistLong.at_Rnd_rPtr(p.T[0])->A_plus[NSEC-1]) << "\t";
//           
//           
//           
//           //co-moving frame::::::::::
//           
//           //net gain
//           //including absorber detuning
//           for(std::size_t k = 0; k < NSEC; k++) f_gain << (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5) / (1.0 + (imag(p.iDelta_omega[k])*imag(p.iDelta_omega[k]))/(p.gamma_2[k]*p.gamma_2[k])) - 0.5*p.alpha_int[k]) << "\t";
//           for(std::size_t k = 0; k < NSEC; k++) f_gain << (p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5) / (1.0 + (imag(p.iDelta_omega[k])*imag(p.iDelta_omega[k]))/(p.gamma_2[k]*p.gamma_2[k])) - 0.5*p.alpha_int[k]) << "\t";
//   //         //excluding absorber detuning
//   //       	for(std::size_t k = 0; k < NSEC; k++) f_gain << p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5) - 0.5*p.alpha_int[k] << "\t";
//   //       	for(std::size_t k = 0; k < NSEC; k++) f_gain << p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5) - 0.5*p.alpha_int[k] << "\t";
//           
//           //net igain
//   //         for(std::size_t k = 0; k < NSEC; k++) f_gain << imag(XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_ES[k] * -p.idelta_omega_ES[k]) + imag(p.iDelta_omega[k] * p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_GS[k] - 0.5)/p.gamma_2[k]) << "\t";
//   //         for(std::size_t k = 0; k < NSEC; k++) f_gain << imag(XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_ES[k] * -p.idelta_omega_ES[k]) + imag(p.iDelta_omega[k] * p.Gamma_r[k] * p.g[k] * (XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_GS[k] - 0.5)/p.gamma_2[k]) << "\t";
//           
//           for(std::size_t k = 0; k < NSEC; k++) f_gain << imag(XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->rho_ES[k] * -p.idelta_omega_ES[k]) << "\t";
//           for(std::size_t k = 0; k < NSEC; k++) f_gain << imag(XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->rho_ES[k] * -p.idelta_omega_ES[k]) << "\t";
//           
//           //pulse propagation
//           for(std::size_t k = 0; k < NSEC; k++) f_pulse << norm(XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->A_plus[k]) << "\t";
//           for(std::size_t k = 0; k < NSEC; k++) f_pulse << norm(XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->A_minus[k]) << "\t";
//           
//           //pulse amplification/net gain
//           f_pulseAmp << norm(XhistLong.at_Rnd_rPtr((NSEC)*p.T[0])->A_plus[0]) - norm(XhistLong.at_Rnd_rPtr((NSEC-1)*p.T[0])->A_minus[0]) << "\t";
//           for(std::size_t k = 1; k < NSEC; k++) f_pulseAmp << norm(XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->A_plus[k]) - norm(XhistLong.at_Rnd_rPtr((NSEC-k)*p.T[0])->A_plus[k-1]) << "\t";
//           //left moving
//           for(std::size_t k = 0; k < NSEC-1; k++) f_pulseAmp << norm(XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->A_minus[k]) - norm(XhistLong.at_Rnd_rPtr((NSEC+k+1)*p.T[0])->A_minus[k+1]) << "\t";
//           f_pulseAmp << norm(XhistLong.at_Rnd_rPtr((NSEC+NSEC-1)*p.T[0])->A_minus[NSEC-1]) - norm(XhistLong.at_Rnd_rPtr((NSEC+NSEC)*p.T[0])->A_plus[NSEC-1]) << "\t";
// 
//           
//           //pulse phase propagation
//           for(int k = 0; k < NSEC; k++) f_phase << arg(XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->A_plus[k]) << "\t";
//           for(std::size_t k = 0; k < NSEC; k++) f_phase << arg(XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->A_minus[k]) << "\t";
//           
//           //pulse phase change
//           
//           f_phaseChange << negatePhaseJump(arg(XhistLong.at_Rnd_rPtr((NSEC)*p.T[0])->A_plus[0]) , arg(XhistLong.at_Rnd_rPtr((NSEC-1)*p.T[0])->A_minus[0])) << "\t";
//           for(std::size_t k = 1; k < NSEC; k++) f_phaseChange << negatePhaseJump(arg(XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->A_plus[k]) , arg(XhistLong.at_Rnd_rPtr((NSEC-k)*p.T[0])->A_plus[k-1])) << "\t";
//           //left moving
//           for(std::size_t k = 0; k < NSEC-1; k++) f_phaseChange << negatePhaseJump(arg(XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->A_minus[k]) , arg(XhistLong.at_Rnd_rPtr((NSEC+k+1)*p.T[0])->A_minus[k+1])) << "\t";
//           f_phaseChange << negatePhaseJump(arg(XhistLong.at_Rnd_rPtr((NSEC+NSEC-1)*p.T[0])->A_minus[NSEC-1]) , arg(XhistLong.at_Rnd_rPtr((NSEC+NSEC)*p.T[0])->A_plus[NSEC-1])) << "\t";
// 
//         // 	f_phaseChange << arg(XhistLong.at_Rnd_rPtr((NSEC)*p.T[0])->A_plus[0]) - arg(XhistLong.at_Rnd_rPtr((NSEC-1)*p.T[0])->A_minus[0]) << "\t";
//         // 	for(int k = 1; k < NSEC; k++) f_phaseChange << arg(XhistLong.at_Rnd_rPtr((NSEC-k-1)*p.T[0])->A_plus[k]) - arg(XhistLong.at_Rnd_rPtr((NSEC-k)*p.T[0])->A_plus[k-1]) << "\t";
//         // 	//left moving
//         // 	for(int k = 0; k < NSEC-1; k++) f_phaseChange << arg(XhistLong.at_Rnd_rPtr((NSEC+k)*p.T[0])->A_minus[k]) - arg(XhistLong.at_Rnd_rPtr((NSEC+k+1)*p.T[0])->A_minus[k+1]) << "\t";
//         // 	f_phaseChange << arg(XhistLong.at_Rnd_rPtr((NSEC+NSEC-1)*p.T[0])->A_minus[NSEC-1]) - arg(XhistLong.at_Rnd_rPtr((NSEC+NSEC)*p.T[0])->A_plus[NSEC-1]) << "\t";
//           
//         f_gain << endl;
//         f_pulse << endl;
//         f_pulseAmp << endl;
//         f_phase << endl;
//         f_phaseChange << endl;
//         
//         }
//         
//       }
//     };
    
//     sm::clck = clock();
//     DDEsolver.DDE_RK4(MLL_system_eqs, noise, &XhistLong, &p, &tn, intTime, dt, outputGainPhase);
//     cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
// 
//     f_netGain.close();
//     f_pulse.close();
//     f_pulseAmp.close();
//     f_phase.close();
//     f_phaseChange.close();
//     
//     //save history to binary file
//     if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
//       XhistLong.tn = tn;
//       XhistLong.save("data/XhistLong.bin");
//     }
// 
//   }

  
  
  
  //P_G sweep with gamma_ES_Q
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_P_G_wgamma_ES_Q" , false)){
    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.P_G;
    str_sweep_parameter = "P_G";
    PtrToSecPar = &p.gamma_ES_Q;
    str_sweep_sec_parameter = "gamma_ES_Q";
    
    sweepMLL();
  }
  
  //P_G sweep with U
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_P_G_wU" , false)){

    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.P_G;
    str_sweep_parameter = "P_G";
    PtrToSecPar = &p.U;
    str_sweep_sec_parameter = "U";
    
    sweepMLL();
  }
  
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_P_G_wNoiseStr" , false)){

    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.P_G;
    str_sweep_parameter = "P_G";
    PtrToSecPar = &p.noiseStr;
    str_sweep_sec_parameter = "noiseStr";
    
    sweepMLL();
  }
  
  
  //P_G sweep with SARS
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_P_G_wSARS" , false)){

    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.P_G;
    str_sweep_parameter = "P_G";
    PtrToSecPar = &p.SARS;
    str_sweep_sec_parameter = "SARS";
    
    sweepMLL();
  }
  
  
  
  //tau sweep with C
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_tau_wC" , false)){
    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.tau;
    str_sweep_parameter = "tau";
    PtrToSecPar = &p.C;
    str_sweep_sec_parameter = "C";
      
    sweepMLL();
  }
  
  if(katana::getCmdOption_bool(argv, argv+argc, "-noiseTest" , false)){
    outputfile.open("data/noiseTest");
    unsigned int maxSteps = (unsigned int)katana::getCmdOption(argv, argv+argc, "-maxSteps" , 10000);
    
    double rw = 0.0;
    double n = 0.0;
    std::vector<double> noiseVec, rwVec;
    noiseVec.reserve(maxSteps), rwVec.reserve(maxSteps);
    for(unsigned int k = 0; k < maxSteps; k ++){
      n = real(katana::gwNoise());
      rw += n;
      noiseVec.push_back(n);
      rwVec.push_back(rw);
      outputfile << k << "\t" << n << "\t" << rw << endl;
    }
    outputfile.close();
    
    std::vector<double> PSpec;
    katana::get_power_spec(noiseVec, PSpec);
    katana::sm_dump_power_spec(PSpec, 1.0, 1.0, "data/noiseTest_PSn");
    
    katana::get_power_spec(rwVec, PSpec);
    katana::sm_dump_power_spec(PSpec, 1.0, 1.0, "data/noiseTest_PSrw");
  }
  
  
  //high resolution power spectrum
  if(katana::getCmdOption_bool(argv, argv+argc, "-testPS" , false)){
    
    double nu = katana::getCmdOption(argv, argv+argc, "-nu" , 1.0/75.16786333214512);
    double maxOutFreq = katana::getCmdOption(argv, argv+argc, "-maxOutFreq" , 280.0);
//     double SampleDt = katana::getCmdOption(argv, argv+argc, "-SampleDt" , 1.0);
    
    std::vector<double> TestDataVec;
    TestDataVec.reserve(outTime_ntn);
    TestDataVec.resize(outTime_ntn);
    
    for(unsigned int k = 0; k<TestDataVec.size(); k++){
      TestDataVec[k] = cos(2*M_PI*nu*dt*k);
    }
    
    
    std::vector<double> powerSpec;
    sm::get_power_spec(TestDataVec, powerSpec, dt);
    sm::dump_power_spec(powerSpec, dt*1E-3, maxOutFreq, "data/testPS");
    
    

  }
  
  
  return 0;
} 
