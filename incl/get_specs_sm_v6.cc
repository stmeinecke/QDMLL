#include <cmath>
#include <iostream>
#include <vector>
#include <complex.h>
#include "fftw3.h"
// #include "get_specs.hh"
#include <fstream>
#include <string>



namespace sm{

  
  ////////////////////////
  /// power spectrum
  ////////////////////////
  
  
  void get_power_spec(std::vector<double> &timeseries, std::vector<double> &powerspec, double dt){
    fftw_complex *data;
    fftw_plan p;
    unsigned int N=timeseries.size();
    data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Write Data in  fftarrays
    for (unsigned int iter = 0; iter <N; iter++){
      data[iter][0]=timeseries[iter];
      data[iter][1]=0;
    };

    p = fftw_plan_dft_1d(N, data, data, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    // Write Output in Outvector

    if(powerspec.size()!=timeseries.size()){
      powerspec.resize(N);
    }

    double norm = dt / timeseries.size(); // = intTime / samples^2  = dt^2 / intTime

    for (unsigned int iter=0; iter<N; iter++){
      powerspec[iter] = norm * ( data[iter][0]*data[iter][0]+data[iter][1]*data[iter][1] );
    }

    fftw_destroy_plan(p);
    fftw_free(data);
  }
    
  void dump_power_spec(std::vector<double> &powerspec, double dt, double maxFreq, std::string fourieroutname){
    std::ofstream fourierout (fourieroutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
    fourierout <<"# (1)frequency (2)component; timestep dt: " << dt << " integration time: " << ((double)powerspec.size())*dt << std::endl;
    unsigned int maxN = (int)fabs(maxFreq*powerspec.size()*dt);
    if(maxN > powerspec.size()/2){
      maxN=powerspec.size()/2;
    }
    for(unsigned int iter=0; iter<=maxN; ++iter ){
      fourierout << iter/(powerspec.size()*dt) << "\t " << powerspec[iter] << std::endl;	
    }
    fourierout.close();
  }
  
  
  //in place powerspec  - powerspec is written into the input data vector
  void get_power_spec(std::vector<double> &data, double dt){
    get_power_spec(data, data, dt);
  }
  
  ////////////////////////
  /// optical spectrum
  ////////////////////////
  
  
  void get_optical_spec(std::vector<std::complex<double>> &timeseries, std::vector<double> &optspec){
    fftw_complex *in, *out;
    fftw_plan p;
    unsigned int N=timeseries.size();
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Write Data in  fftarrays
    for (unsigned int iter = 0; iter <N; iter++){
      in[iter][0]=timeseries[iter].real();
      in[iter][1]=timeseries[iter].imag();
    };

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    // Write Output in Outvector

    if(optspec.size()!=timeseries.size()){
      optspec.resize(N);
    }

    for (unsigned int iter=0; iter<N; iter++){
      optspec[iter] = out[iter][0]*out[iter][0]+out[iter][1]*out[iter][1];
    }

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);	
  }
  
  
  void dump_optical_spec(std::vector<double> &spec, double dt, double maxFreq, std::string fourieroutname){
    std::ofstream fourierout (fourieroutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
    fourierout <<"# (1)frequency (2)component; timestep dt: " << dt << "integration time: " << (double)spec.size()*dt << std::endl;
    unsigned int maxN = (int)fabs(maxFreq*spec.size()*dt);
    if(maxN > spec.size()/2){
      maxN=spec.size()/2;
    }
    for(int iter=(int) -maxN; iter<(int) maxN; ++iter ){
      fourierout << iter/(spec.size()*dt) << "\t " << spec[(iter+spec.size())%spec.size()] << std::endl;
    }
    fourierout.close();
  }
  
  
  ////////////////////////
  /// fourier spectrum
  ////////////////////////
  
  void get_fourier_spec(std::vector<std::complex<double>> &timeseries, std::vector<std::complex<double>> &fourierspec){
    fftw_complex *in, *out;
    fftw_plan p;
    unsigned int N=timeseries.size();
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Write Data in  fftarrays
    for (unsigned int iter = 0; iter <N; iter++){
      in[iter][0]=timeseries[iter].real();
      in[iter][1]=timeseries[iter].imag();
    };

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    // Write Output in Outvector

    if(fourierspec.size()!=timeseries.size()){
      fourierspec.resize(N);
    }

    for (unsigned int iter=0; iter<N; iter++){
      fourierspec[iter]=std::complex<double> (out[iter][0], out[iter][1]);
    }

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);	
  }
  
  void dump_fourier_spec(std::vector<std::complex<double>> &spec, double dt, double maxFreq, std::string fourieroutname){
    std::ofstream fourierout (fourieroutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
    fourierout <<"# (1)frequency (2)component; timestep dt: " << dt << "integration time: " << (double)spec.size()*dt << std::endl; 
    unsigned int maxN = (int)fabs(maxFreq*spec.size()*dt);
    if(maxN > spec.size()/2){
      maxN=spec.size()/2;
    }
    for(int iter=(int) -maxN; iter<(int) maxN; ++iter ){
      fourierout << iter/(spec.size()*dt) << "\t " << spec[(iter+spec.size())%spec.size()].real() << "\t" << spec[(iter+spec.size())%spec.size()].imag() << std::endl;
    }
    fourierout.close();
  }
  
  
}



namespace katana 
{

using Time=double;

void get_power_spec(std::vector<double> &timeseries, std::vector<double> &powerspec)
	{
	 fftw_complex *in, *out;
         fftw_plan p;
         unsigned int N=timeseries.size();
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	 // Write Data in  fftarrays
	 for (unsigned int iter = 0; iter <N; iter++)
		{
			in[iter][0]=timeseries[iter];
			in[iter][1]=0;
		};

         p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
         
	 fftw_execute(p); /* repeat as needed */

	 // Write Output in Outvector

	if(powerspec.size()!=timeseries.size())
		{powerspec.resize(N);}

       
	for (unsigned int iter=0; iter<N; iter++)
		{
// 			powerspec[iter]=std::sqrt(out[iter][0]*out[iter][0]+out[iter][1]*out[iter][1]);
      powerspec[iter]=out[iter][0]*out[iter][0]+out[iter][1]*out[iter][1]; //fix by SM 2020/02/10
		}

 
         fftw_destroy_plan(p);
         fftw_free(in); fftw_free(out);	

	}
	

void get_power_spec_normalized(std::vector<double> &timeseries, std::vector<double> &powerspec, double dt)
	{
	 fftw_complex *in, *out;
         fftw_plan p;
         unsigned int N=timeseries.size();
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	 // Write Data in  fftarrays
	 for (unsigned int iter = 0; iter <N; iter++)
		{
			in[iter][0]=timeseries[iter];
			in[iter][1]=0;
		};

         p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
         
	 fftw_execute(p); /* repeat as needed */

	 // Write Output in Outvector

	if(powerspec.size()!=timeseries.size())
		{powerspec.resize(N);}

  double norm = dt / timeseries.size(); // = intTime / samples^2  = dt^2 / intTime
       
	for (unsigned int iter=0; iter<N; iter++)
		{
      powerspec[iter] = norm * ( out[iter][0]*out[iter][0]+out[iter][1]*out[iter][1] );
		}

 
         fftw_destroy_plan(p);
         fftw_free(in); fftw_free(out);	

	}


void get_optical_spec(std::vector<std::complex<double>> &timeseries, std::vector<double> &optspec)
	{
	 fftw_complex *in, *out;
         fftw_plan p;
         unsigned int N=timeseries.size();
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	 // Write Data in  fftarrays
	 for (unsigned int iter = 0; iter <N; iter++)
		{
			in[iter][0]=timeseries[iter].real();
			in[iter][1]=timeseries[iter].imag();
		};

         p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
         
	 fftw_execute(p); /* repeat as needed */

	 // Write Output in Outvector

	if(optspec.size()!=timeseries.size())
		{optspec.resize(N);}

       
	for (unsigned int iter=0; iter<N; iter++)
		{
			optspec[iter]=std::sqrt(out[iter][0]*out[iter][0]+out[iter][1]*out[iter][1]);
		}

 
         fftw_destroy_plan(p);
         fftw_free(in); fftw_free(out);	
}

void get_fourier_spec(std::vector<std::complex<double>> &timeseries, std::vector<std::complex<double>> &fourierspec)
	{
	 fftw_complex *in, *out;
         fftw_plan p;
         unsigned int N=timeseries.size();
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	 // Write Data in  fftarrays
	 for (unsigned int iter = 0; iter <N; iter++)
		{
			in[iter][0]=timeseries[iter].real();
			in[iter][1]=timeseries[iter].imag();
		};

         p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
         
	 fftw_execute(p); /* repeat as needed */

	 // Write Output in Outvector

	if(fourierspec.size()!=timeseries.size())
		{fourierspec.resize(N);}

       
	for (unsigned int iter=0; iter<N; iter++)
		{
			fourierspec[iter]=std::complex<double> (out[iter][0], out[iter][1]);
		}

 
         fftw_destroy_plan(p);
         fftw_free(in); fftw_free(out);	
}




// MaxN deklariert die gesamtlänge des arrays. dt ist die interne Zeitskala und wird hier mit Einheit ns angenommen, sodass die 
// Frequenz Sekunden ergibt

void dump_power_spec(std::vector<double> &powerspec, Time dt, unsigned int maxN, std::string fourieroutname)
	{	std::ofstream fourierout (fourieroutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		fourierout <<"# Frequency [GHz] \t"  <<  "Fouriercomponent" << std::endl; 
		if(powerspec.size()<maxN) maxN=powerspec.size();
		for(int iter=(int) -maxN/2; iter<(int) maxN/2; ++iter )  // maxN auf int casten, for safety!
			{fourierout << iter/(powerspec.size()*dt) << "\t " << powerspec[(iter+powerspec.size())%powerspec.size()]/powerspec[0] << "\t" << (iter+powerspec.size())%powerspec.size() <<  std::endl;	}					

		fourierout.close();
	
	}
// eigentlich die gleiche Funktion. get_optical_spec() liefert den Betrag des optischen Spektrums.
void dump_optical_spec(std::vector<double> &spec, Time dt, unsigned int maxN, std::string fourieroutname)
	{	std::ofstream fourierout (fourieroutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		fourierout <<"# Frequency [GHz] \t"  <<  "Fouriercomponent" << std::endl; 
		if(spec.size()<maxN) maxN=spec.size();
		for(int iter=(int) -maxN/2; iter<(int) maxN/2; ++iter )
			{fourierout << iter/(spec.size()*dt) << "\t " << spec[(iter+spec.size())%spec.size()]/spec[0] << std::endl;	}					
		fourierout.close();
	
	}

void dump_fourier_spec(std::vector<std::complex<double>> &spec, Time dt, unsigned int maxN, std::string fourieroutname)
	{	std::ofstream fourierout (fourieroutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		fourierout <<"# Frequency [GHz] \t"  <<  "Fouriercomponent" << std::endl; 
		if(spec.size()<maxN) maxN=spec.size();
		for(int iter=(int) -maxN/2; iter<(int) maxN/2; ++iter )
			{fourierout << iter/(spec.size()*dt) << "\t " << spec[(iter+spec.size())%spec.size()].real()/spec[0] << std::endl;	}				//ANDRE, mach das mal zuende	
		fourierout.close();
	
	}

void dump_power_spec_pp(std::vector<double> &powerspec, Time dt, unsigned int maxN, std::string fourieroutname)
	{	std::ofstream fourierout (fourieroutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		fourierout <<"# Frequency [GHz] \t"  <<  "Fouriercomponent" << std::endl; 
		if(powerspec.size()<maxN) maxN=powerspec.size();
		for(int iter=0; iter<(int) maxN/2; ++iter )  // maxN auf int casten, for safety!
			{fourierout << iter/(powerspec.size()*dt) << "\t " << powerspec[(iter+powerspec.size())%powerspec.size()]/powerspec[0] << "\t" << (iter+powerspec.size())%powerspec.size() <<  std::endl;	}					

		fourierout.close();
	
	}
	
void sm_dump_power_spec(std::vector<double> &powerspec, Time dt, double maxFreq, std::string fourieroutname)
	{	std::ofstream fourierout (fourieroutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		fourierout <<"# Frequency [GHz] \t"  <<  "Fouriercomponent" << std::endl; 
		unsigned int maxN = (int)fabs(maxFreq*powerspec.size()*dt);
    if(powerspec.size()/2<maxN) maxN=powerspec.size()/2;
		for(unsigned int iter=0; iter<=maxN; ++iter )
// 			{fourierout << iter/(powerspec.size()*dt) << "\t " << powerspec[iter]/powerspec[0] << std::endl;	}
			{fourierout << iter/(powerspec.size()*dt) << "\t " << powerspec[iter] << std::endl;	}

		fourierout.close();
	
	}
	
void sm_dump_optical_spec(std::vector<double> &spec, Time dt, double maxFreq, std::string fourieroutname)
	{	std::ofstream fourierout (fourieroutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		fourierout <<"# Frequency [GHz] \t"  <<  "Fouriercomponent" << std::endl; 
    unsigned int maxN = (int)fabs(maxFreq*spec.size()*dt);
    if(spec.size()<maxN/2) maxN=spec.size()/2;
		for(int iter=(int) -maxN; iter<(int) maxN; ++iter )
			{fourierout << iter/(spec.size()*dt) << "\t " << spec[(iter+spec.size())%spec.size()] << std::endl;	}					
		fourierout.close();
	
	}



}
