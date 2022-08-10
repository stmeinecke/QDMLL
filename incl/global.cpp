#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <functional>
#include <chrono>
#include <time.h>
#include <tuple>
#include <cstdlib>
#include <sys/time.h>
#include <climits>
#include <sstream>
#include <algorithm>
#include <iomanip> 
#include <random>


//for continous (slow) output 
std::ofstream outputfile;

namespace sm{
  
  //clock for timing
  clock_t clck;

  //constants
  const double PI =  3.1415926535897932384626433;

  //imagenary unit
  const std::complex<double> img = std::complex<double>(0.0,1.0);

  
  const double e0 = 1.60218E-19;
  const double m_e = 9.10938E-31;
  const double h = 6.62607E-34;
  const double hbar = 1.05457E-34;
  const double eps0 = 8.85419E-12;
  const double kB = 1.38065E-23;
  const double c0 = 2.99792E8;


  //input degrees into tan
  double tanDeg(double theta){
    return tan(theta * PI / 180.0);
  }


  long int DebugCounter = 0;



  ///////////////lookup table for exp (real valued) /////////////////
  #define EXP_SAMPLES 10000
  #define EXP_STEP 0.005
  #define EXP_STEP_invers 200.0 // == 1.0/EXP_STEP
  // exp_y = Schnittpunkt mit y-Achse
  // exp_slope = Steigung
  // Exponentialfunktionergibt sich aus exp(x) = exp_y[(int)xx] + exp_slope[(int)xx] * xx, mit xx := x / EXP_STEP
  double exp_y[EXP_SAMPLES], exp_slope[EXP_SAMPLES];
  void init_expf() {
      for(int i=0; i<EXP_SAMPLES; i++) {
	  exp_y[i] = exp(i*EXP_STEP);
      }
      for(int i=0; i<EXP_SAMPLES-1; i++) {
	  exp_slope[i] = exp_y[i+1]-exp_y[i];
	  exp_y[i] = exp_y[i]-i*(exp_y[i+1]-exp_y[i]);
      }
  }

  double expf( double ex ) {
      if(ex < 0) return 1.0/expf(-ex);
      double EE = ex*EXP_STEP_invers;
      int EEi = (int)EE;
      return exp_slope[EEi]*EE + exp_y[EEi];
  }


  inline double ExpR(double x){
    return exp(x);
  }

  inline std::complex<double> ExpC(std::complex<double> z){
    return exp(z);
  }
  
  
  void window_Hann(std::vector<double> &in, std::vector<double> &out){
    if(in.size() != out.size()) out.resize(in.size());
    int size = in.size();
    for(int k = 0; k < size; k++){
      out[k] = in[k] * (0.5 - 0.5*cos(2.0 * M_PI * k / (size-1)));
    }
  }
  
  void window_Hann(std::vector<double> &in){ //in place window
    window_Hann(in,in);
  }
  
  void window_Hann(std::vector<std::complex<double>> &in, std::vector<std::complex<double>> &out){
    if(in.size() != out.size()) out.resize(in.size());
    int size = in.size();
    for(int k = 0; k < size; k++){
      out[k] = in[k] * (0.5 - 0.5*cos(2.0 * M_PI * k / (size-1)));
    }
  }
  
  void window_Hamming(std::vector<double> &in, std::vector<double> &out){
    if(in.size() != out.size()) out.resize(in.size());
    int size = in.size();
    for(int k = 0; k < size; k++){
      out[k] = in[k] * (0.54 - 0.46*cos(2.0 * M_PI * k / (size-1)));
    }
  }
  
  void window_Hamming(std::vector<double> &in){ //in place window
    window_Hamming(in,in);
  }
  
  void window_Blackman(std::vector<double> &in, std::vector<double> &out){
    if(in.size() != out.size()) out.resize(in.size());
    int size = in.size();
    for(int k = 0; k < size; k++){
      out[k] = in[k] * (0.42 - 0.5*cos(2.0 * M_PI * k / (size-1)) + 0.08*cos(4.0 * M_PI * k / (size-1)));
    }
  }
  
  void window_Blackman(std::vector<double> &in){ //in place window
    window_Blackman(in,in);
  }
  
  void window_Blackman_exact(std::vector<double> &in, std::vector<double> &out){
    if(in.size() != out.size()) out.resize(in.size());
    int size = in.size();
    for(int k = 0; k < size; k++){
      out[k] = in[k] * (7938.0/18608.0 - 9240.0/18608.0*cos(2.0 * M_PI * k / (size-1)) + 1430.0/18608.0*cos(4.0 * M_PI * k / (size-1)));
    }
  }
  
  void window_Blackman_exact(std::vector<double> &in){ //in place window
    window_Blackman_exact(in,in);
  }
  
  void window_Nuttall(std::vector<double> &in, std::vector<double> &out){
    if(in.size() != out.size()) out.resize(in.size());
    int size = in.size();
    for(int k = 0; k < size; k++){
      out[k] = in[k] * (0.355768 - 0.487396*cos(2.0 * M_PI * k / (size-1)) + 0.144232*cos(4.0 * M_PI * k / (size-1)) - 0.012604*cos(6.0 * M_PI * k / (size-1)));
    }
  }
  
  void window_Nuttall(std::vector<double> &in){ //in place window
    window_Nuttall(in,in);
  }
  
  void window_Blackman_Nuttall(std::vector<double> &in, std::vector<double> &out){
    if(in.size() != out.size()) out.resize(in.size());
    int size = in.size();
    for(int k = 0; k < size; k++){
      out[k] = in[k] * (0.3635819 - 0.4891775*cos(2.0 * M_PI * k / (size-1)) + 0.1365995*cos(4.0 * M_PI * k / (size-1)) - 0.0106411*cos(6.0 * M_PI * k / (size-1)));
    }
  }
  
  void window_Blackman_Nuttall(std::vector<double> &in){ //in place window
    window_Blackman_Nuttall(in,in);
  }
  
}


namespace katana{
  //from Chris

  ////////////////////////////////// Noise ///////////////////
  std::complex<double> gwNoise(){
    double spont_phi= (rand()%RAND_MAX)/(double) RAND_MAX *2.0*3.14159265359;
    double spont_amp= sqrt(-2* log((rand()%RAND_MAX+1)/(double) RAND_MAX)); //gibt nur die Zufallszahl nach normierter Gaussverteilungi
    return  {spont_amp*cos(spont_phi), spont_amp*sin(spont_phi)};
  }

  class Seed{
  public:
    Seed();
    void set_seed();
    void set_seed(int );
  };

  Seed::Seed(){}
  void Seed::set_seed(){
    struct timeval t;
    gettimeofday(&t, NULL);

    unsigned int time_sec = t.tv_sec;     // Current seconds since epoc (midnight, jan 1, 1970)
    unsigned int time_musec = t.tv_usec; // Number of microseconds this second
    unsigned long long int rawseed;			// long long is sufficient in every case. for AMD64 same as u long int
    unsigned int seed;


    std::stringstream seedstream;
    std::string seedstring;
						    
    seedstream << time_sec << time_musec; 
    seedstream >> rawseed; // time in sec+musec into seedstring

    seed= (unsigned int) rawseed % UINT_MAX;
    std::srand(seed);									
  }
  void Seed::set_seed(int seed){ 
    srand(seed);
  }
  ////////////////////////////////// Noise ///////////////////
}
  
