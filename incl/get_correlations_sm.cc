#include <cmath>
#include <iostream>
#include <vector>
#include <complex.h>
#include <algorithm>
#include "fftw3.h"
// #include "get_correlations.hh"
#include <fstream>
#include <string>
// #include <armadillo>
// #include "arma_wrap.hh"

namespace katana 
{

using Time=double;

/*########################################################################################################
########                                                                                  ################
######## Using the Wiener–Khinchin theorem to compute auto and cross correlation funtions.################
######## This Module also provides an normalization function, and a routine to transform  ################
######## into amplitude and phase.                                                        ################
######## AUTHOR: Chris									  ################
########################################################################################################*/




double get_autocorr(const std::vector<std::complex<double>> &timeseries, std::vector<std::complex<double>> &autocorrelation)  // Das hier wird typischerweise die g2-funktion darstellen
	{
	 fftw_complex *in_ts, *out_spec, *power_spec, *autocorr; 
	 
         fftw_plan p, p2;
         unsigned int N=timeseries.size();
         in_ts = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // time-series vector, input
         out_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);      // intermediate variable, to calculate the initial fourier transform of in_ts
	
	 // Write Data in  fftarrays
	 for (unsigned int iter = 0; iter <N; iter++)
		{
			in_ts[iter][0]=timeseries[iter].real();
			in_ts[iter][1]=timeseries[iter].imag();
		};

         p = fftw_plan_dft_1d(N, in_ts, out_spec, FFTW_FORWARD, FFTW_ESTIMATE);
         
	 fftw_execute(p);
         fftw_destroy_plan(p);
	 fftw_free(in_ts);

	// Transform is done, now we need to compute the powerspectrum of it.

	 power_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // power_spec vec
	
	 for (unsigned int iter=0; iter<N; iter++)
		{
			power_spec[iter][0]=out_spec[iter][0]*out_spec[iter][0]+ out_spec[iter][1]*out_spec[iter][1];
			power_spec[iter][1]=0;
		}

	
	// Now, calculate inverse fft of power_spec 

	 fftw_free(out_spec);
	 autocorr=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         //  autocorr vector

         p2 = fftw_plan_dft_1d(N, power_spec, autocorr, FFTW_BACKWARD, FFTW_ESTIMATE);
	 fftw_execute(p2);
 	 fftw_destroy_plan(p2);

	 double returnval=power_spec[0][0];  // Store mean Value of the time series to calculate g2!

	// Write Output in Outvector
	 if(autocorrelation.size()!=timeseries.size())
		{autocorrelation.resize(N);}

	 for (unsigned int iter=0; iter<N; iter++)
		{
			autocorrelation[iter]={autocorr[iter][0],autocorr[iter][1]};
		}
	
	 
          fftw_free(power_spec); fftw_free(autocorr);	
	return returnval;
	}
	
double get_autocorr_norm_out(const std::vector<std::complex<double>> &timeseries, std::vector<double> &autocorrelation)  // Das hier wird typischerweise die g2-funktion darstellen
	{
	 fftw_complex *in_ts, *out_spec, *power_spec, *autocorr; 
	 
         fftw_plan p, p2;
         unsigned int N=timeseries.size();
         in_ts = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // time-series vector, input
         out_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);      // intermediate variable, to calculate the initial fourier transform of in_ts
	
	 // Write Data in  fftarrays
	 for (unsigned int iter = 0; iter <N; iter++)
		{
			in_ts[iter][0]=timeseries[iter].real();
			in_ts[iter][1]=timeseries[iter].imag();
		};

         p = fftw_plan_dft_1d(N, in_ts, out_spec, FFTW_FORWARD, FFTW_ESTIMATE);
         
	 fftw_execute(p);
         fftw_destroy_plan(p);
	 fftw_free(in_ts);

	// Transform is done, now we need to compute the powerspectrum of it.

	 power_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // power_spec vec
	
	 for (unsigned int iter=0; iter<N; iter++)
		{
			power_spec[iter][0]=out_spec[iter][0]*out_spec[iter][0]+ out_spec[iter][1]*out_spec[iter][1];
			power_spec[iter][1]=0;
		}

	
	// Now, calculate inverse fft of power_spec 

	 fftw_free(out_spec);
	 autocorr=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         //  autocorr vector

         p2 = fftw_plan_dft_1d(N, power_spec, autocorr, FFTW_BACKWARD, FFTW_ESTIMATE);
	 fftw_execute(p2);
 	 fftw_destroy_plan(p2);

	 double returnval=power_spec[0][0];  // Store mean Value of the time series to calculate g2!

	// Write Output in Outvector
	 if(autocorrelation.size()!=timeseries.size())
		{autocorrelation.resize(N);}

	 for (unsigned int iter=0; iter<N; iter++)
		{
			autocorrelation[iter]=autocorr[iter][0]*autocorr[iter][0] + autocorr[iter][1]*autocorr[iter][1];
		}
	
	 
          fftw_free(power_spec); fftw_free(autocorr);	
	return returnval;
	}

double get_autocorr_real(const std::vector<double> &realtimeseries, std::vector<double> &realautocorrelation)  // Das hier wird typischerweise die g2-funktion darstellen
	{
	 fftw_complex *in_ts, *out_spec, *power_spec, *autocorr; 
	 
         fftw_plan p, p2;
         unsigned int N=realtimeseries.size();
         in_ts = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // time-series vector, input
         out_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);      // intermediate variable, to calculate the initial fourier transform of in_ts
	
	 // Write Data in  fftarrays
	 for (unsigned int iter = 0; iter <N; iter++)
		{
			in_ts[iter][0]=realtimeseries[iter];
			in_ts[iter][1]=0;		};

         p = fftw_plan_dft_1d(N, in_ts, out_spec, FFTW_FORWARD, FFTW_ESTIMATE);
         
	 fftw_execute(p);
         fftw_destroy_plan(p);
	 fftw_free(in_ts);

	// Transform is done, now we need to compute the powerspectrum of it.

	 power_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // power_spec vec
	
	 for (unsigned int iter=0; iter<N; iter++)
		{
			power_spec[iter][0]=out_spec[iter][0]*out_spec[iter][0]+ out_spec[iter][1]*out_spec[iter][1];
			power_spec[iter][1]=0;
		}

	
	// Now, calculate inverse fft of power_spec 

	 double returnval=std::sqrt(power_spec[0][0]);  // Store mean Value of the time series to calculate g2!
	 fftw_free(out_spec);
	 autocorr=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         //  autocorr vector

         p2 = fftw_plan_dft_1d(N, power_spec, autocorr, FFTW_BACKWARD, FFTW_ESTIMATE);
	 fftw_execute(p2);
         fftw_destroy_plan(p2);


	// Write Output in Outvector
	 if(realautocorrelation.size()!=realtimeseries.size())
		{realautocorrelation.resize(N);}

	 for (unsigned int iter=0; iter<N; iter++)
		{
			realautocorrelation[iter]=autocorr[iter][0];
		}
	
	 
          fftw_free(power_spec); fftw_free(autocorr);	
	return returnval;
	}

double get_autocorr(const std::vector<double> &realtimeseries, std::vector<double> &autocorrelation)
	{	
		return get_autocorr_real( realtimeseries, autocorrelation);	
	}



// Calculation of thje cross correlation function using Wiener Khinchin Theorem! Here, the input is assumed to be real, yet the correlation function
// should in general be complex. Thus, we need a complex vector for the output!

std::vector<double> get_crosscorr(const std::vector<std::complex<double>> &timeseries1,const std::vector<std::complex<double>> &timeseries2, std::vector<std::complex<double>> &crosscorrelation)  // Das hier wird typischerweise die g2-funktion darstellen
	{
	 fftw_complex *in1_ts, *in2_ts, *out1_spec, *out2_spec, *corr_spec, *crosscorr; 
	 
         fftw_plan p_forward1, p_forward2, p_backward;  // Only one fourier transform for the back transformation needed.
         unsigned int N=timeseries1.size();
	 if(timeseries2.size()==N)
		{/*std::cout << "Sizeof Timeseries match! Continue calculating" << std::endl;*/}
	 else {std::cout << "Timeseries don't match! ABORTING!" << std::endl; return std::vector<double> {0,0};}

         in1_ts = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // time-series vector, input
         in2_ts = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // time-series vector, input
         out1_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);      // intermediate variable, to calculate the initial fourier transform of in_ts
         out2_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);      // intermediate variable, to calculate the initial fourier transform of in_ts
	
	 // Write Data in  fftarrays
	 for (unsigned int iter = 0; iter <N; iter++)
		{
			in1_ts[iter][0]=timeseries1[iter].real();
			in1_ts[iter][1]=timeseries1[iter].imag();
			in2_ts[iter][0]=timeseries2[iter].real();
			in2_ts[iter][1]=timeseries2[iter].imag();
		};

         p_forward1 = fftw_plan_dft_1d(N, in1_ts, out1_spec, FFTW_FORWARD, FFTW_ESTIMATE);
         p_forward2 = fftw_plan_dft_1d(N, in2_ts, out2_spec, FFTW_FORWARD, FFTW_ESTIMATE);
         
	 fftw_execute(p_forward1);
	 fftw_execute(p_forward2);
         fftw_destroy_plan(p_forward1);
         fftw_destroy_plan(p_forward2);
	 fftw_free(in1_ts);
	 fftw_free(in2_ts);

	// Transform is done, now we need to compute the powerspectrum of it.

	 corr_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // power_spec vec
	
	 for (unsigned int iter=0; iter<N; iter++)
		{
			corr_spec[iter][0]=out1_spec[iter][0]*out2_spec[iter][0]+ out1_spec[iter][1]*out2_spec[iter][1];  //This is the complex multiplation of
			corr_spec[iter][1]=out1_spec[iter][1]*out2_spec[iter][0]- out1_spec[iter][0]*out2_spec[iter][1];  //out_spec1 * complex_conj(out_spec2)
		}
	
	std::vector<double> returnvals={out1_spec[0][0], out2_spec[0][0]};  // Store mean Values of the time series' to calculate g2!
	
	// Now, calculate inverse fft of power_spec 

	 fftw_free(out1_spec);
	 fftw_free(out2_spec);
	 crosscorr=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         //  autocorr vector

         p_backward = fftw_plan_dft_1d(N, corr_spec, crosscorr, FFTW_BACKWARD, FFTW_ESTIMATE);
	 fftw_execute(p_backward);
         fftw_destroy_plan(p_backward);
	
	// Write Output in Outvector

	 if(crosscorrelation.size()!=N)
		{crosscorrelation.resize(N);}

	 for (unsigned int iter=0; iter<N; iter++)
		{
			crosscorrelation[iter]={crosscorr[iter][0], crosscorr[iter][1]};
		}
	
 
          fftw_free(corr_spec); fftw_free(crosscorr);	
	return returnvals;
	}

std::vector<double> get_crosscorr_real(const std::vector<double> &realtimeseries1,const std::vector<double> &realtimeseries2, std::vector<double> &crosscorrelation)  // Das hier wird typischerweise die g2-funktion darstellen
	{
	 fftw_complex *in1_ts, *in2_ts, *out1_spec, *out2_spec, *corr_spec, *crosscorr; 
	 
         fftw_plan p_forward1, p_forward2, p_backward;  // Only one fourier transform for the back transformation needed.
         unsigned int N=realtimeseries1.size();
	 if(realtimeseries2.size()==N)
		{/*std::cout << "Sizeof Timeseries match! Continue calculating" << std::endl;*/}
	 else {std::cout << "Timeseries don't match! ABORTING!" << std::endl; return std::vector<double> {0,0};}

         in1_ts = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // time-series vector, input
         in2_ts = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // time-series vector, input
         out1_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);      // intermediate variable, to calculate the initial fourier transform of in_ts
         out2_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);      // intermediate variable, to calculate the initial fourier transform of in_ts
	
	 // Write Data in  fftarrays
	 for (unsigned int iter = 0; iter <N; iter++)
		{
			in1_ts[iter][0]=realtimeseries1[iter];
			in2_ts[iter][0]=realtimeseries2[iter];
			in1_ts[iter][1]=in2_ts[iter][1]=0;
		};

         p_forward1 = fftw_plan_dft_1d(N, in1_ts, out1_spec, FFTW_FORWARD, FFTW_ESTIMATE);
         p_forward2 = fftw_plan_dft_1d(N, in2_ts, out2_spec, FFTW_FORWARD, FFTW_ESTIMATE);
         
	 fftw_execute(p_forward1);
	 fftw_execute(p_forward2);
         fftw_destroy_plan(p_forward1);
         fftw_destroy_plan(p_forward2);
	 fftw_free(in1_ts);
	 fftw_free(in2_ts);

	// Transform is done, now we need to compute the powerspectrum of it.

	 corr_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         // power_spec vec
	
	 for (unsigned int iter=0; iter<N; iter++)
		{
			corr_spec[iter][0]=out1_spec[iter][0]*out2_spec[iter][0]+ out1_spec[iter][1]*out2_spec[iter][1];  //This is the complex multiplation of
			corr_spec[iter][1]=out1_spec[iter][1]*out2_spec[iter][0]- out1_spec[iter][0]*out2_spec[iter][1];  //out_spec1 * complex_conj(out_spec2)
		}
	
	std::vector<double> returnvals={out1_spec[0][0], out2_spec[0][0]};  // Store mean Values of the time series' to calculate g2!
	
	// Now, calculate inverse fft of power_spec 

	 fftw_free(out1_spec);
	 fftw_free(out2_spec);
	 crosscorr=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);         //  autocorr vector

         p_backward = fftw_plan_dft_1d(N, corr_spec, crosscorr, FFTW_BACKWARD, FFTW_ESTIMATE);
	 fftw_execute(p_backward);
 	 fftw_destroy_plan(p_backward);
	
	// Write Output in Outvector

	 if(crosscorrelation.size()!=N)
		{crosscorrelation.resize(N);}

	 for (unsigned int iter=0; iter<N; iter++)
		{
			crosscorrelation[iter]=crosscorr[iter][0]; //Fuer reellen Input ist auch die Korrelationsfunktion reell!
		}
 
          fftw_free(corr_spec); fftw_free(crosscorr);	
	return returnvals;
	}

std::vector<double> get_crosscorr(const std::vector<double> &Realtimeseries1,const std::vector<double> &Realtimeseries2, std::vector<double> &g2)
	{
		return get_crosscorr_real(Realtimeseries1, Realtimeseries2, g2);

	}


void normalize_correlation(std::vector<double> &correlation, double normalizer)
	{	
		std::transform(correlation.begin(), correlation.end(), correlation.begin(), [&](double val) {return val/normalizer;});
	}

void normalize_correlation(std::vector<std::complex<double>> &correlation, double normalizer)
	{	
		std::transform(correlation.begin(), correlation.end(), correlation.begin(), [&](std::complex<double> val) {return val/normalizer;});
	}




void get_g2(const std::vector<double> &Timeseries, std::vector<double> &g2)
	{
	 	double normalizer=get_autocorr(Timeseries, g2);
		normalize_correlation(g2, normalizer*normalizer);	
	}



void get_g2(const std::vector<double> &Timeseries1,const std::vector<double> &Timeseries2, std::vector<double> &g2)
	{
		std::vector<double> normalizer=get_crosscorr(Timeseries1, Timeseries2, g2);
		double norm=normalizer[0]*normalizer[1];
		normalize_correlation(g2, norm);	
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copy functions but with little different handling, directly returns the correlation vector. For brevity in the main

std::vector<std::complex<double>> get_autocorr(const std::vector<std::complex<double>> &timeseries)
	{ std::vector<std::complex<double>> autocorrelation;
	  get_autocorr(timeseries, autocorrelation);
	  return autocorrelation;
	};
std::vector<double> get_autocorr_real(const std::vector<double> &timeseries)
	{ std::vector<double> autocorrelation;
	  get_autocorr(timeseries, autocorrelation);
	  return autocorrelation;
	};

std::vector<double> get_autocorr(const std::vector<double> &timeseries) //Overload of get_autocorr with exact definition of get_autocorr_real 
	{ std::vector<double> autocorrelation;
	  get_autocorr( timeseries, autocorrelation);
	  return autocorrelation;
	};



//Now for the g2 functions:

std::vector<double>  get_g2(const std::vector<double> &timeseries) 								//Auto g2 correlation
	{	std::vector<double> g2;
		double normalizer=get_autocorr(timeseries, g2);
		normalize_correlation(g2, normalizer*normalizer);	
		return g2;
	};
std::vector<double>  get_g2_cc(const std::vector<double> &timeseries1, const std::vector<double> &timeseries2) 		//Cross g2 correlation
	{	std::vector<double> g2;
		std::vector<double> normalizer=get_crosscorr(timeseries1, timeseries2, g2);
		double norm=normalizer[0]*normalizer[1];
		normalize_correlation(g2, norm);
		return g2;
	};






void print_correlation(const std::vector<double> &correlation, Time dt, std::string correlationoutname, std::size_t maxN , std::size_t every, std::string unit)
	{	std::ofstream corr_out (correlationoutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
// 		std::cout << "Write file: "<< correlationoutname << std::endl;
    corr_out.precision(10);
		size_t correlation_size=correlation.size();
		if(maxN>correlation_size || maxN==0) maxN=correlation_size;
		
		corr_out <<"# Delay" << "["<< unit << "]\t"  <<  "Correlation \t" << std::endl; 
		for(int iter= -int(maxN/2); iter<int(maxN/2); iter+=every )  
			{corr_out << iter*dt << "\t" << correlation[(iter+correlation_size)%correlation_size] << std::endl;
			}				
		corr_out.close();	
	
	}

void print_correlation(const std::vector<std::complex<double>> &correlation, Time dt, std::string correlationoutname, std::size_t maxN, std::size_t every, std::string unit)
	{	std::ofstream corr_out (correlationoutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
// 		std::cout << "Write file: "<< correlationoutname << std::endl;
    corr_out.precision(10);
		size_t correlation_size= correlation.size();
	
		if(maxN>correlation_size || maxN==0) maxN=correlation_size;
		corr_out <<"# Delay" << "["<< unit << "]\t"  <<  "Correlation \t" << std::endl; 
		for(int iter= -int(maxN/2); iter<int(maxN/2); iter+=every )  
			{corr_out << iter*dt << "\t" << correlation[(iter+correlation_size)%correlation_size] << std::endl;
			}				
		corr_out.close();	
	}

void print_correlations(const std::vector<std::vector<double>> &correlation, Time dt, std::string correlationoutname, std::string correlationheader, std::size_t maxN, int every, unsigned int precision)
	{	std::ofstream corrout (correlationoutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		corrout.precision(precision);
		size_t correlation_size= correlation[0].size();
// 		std::cout << "Write file: "<< correlationoutname << std::endl;
		corrout << correlationheader << std::endl; 

		if(maxN>correlation_size || maxN==0) maxN=correlation_size;
		for(int iter= -int(maxN/2); iter<int(maxN/2); iter+=every )  
			{corrout << iter*dt;
				for(std::size_t elem=0; elem<correlation[0].size(); elem++)
					{corrout << "\t" << correlation[(iter+correlation_size)%correlation_size][elem];}
				corrout << std::endl;
			}					
		corrout.close();
	
	}

////sm
void sm_print_correlation(const std::vector<double> &correlation, Time dt, std::string correlationoutname, double  maxT)
	{	std::ofstream corr_out (correlationoutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
    corr_out.precision(10);
		size_t correlation_size=correlation.size();
    unsigned int maxN = (unsigned int)(maxT/dt);
		if(maxN>correlation_size/2 || maxN==0) maxN=correlation_size/2;
		
		corr_out <<"# Delay" << "\t"  <<  "Correlation \t" << std::endl; 
		for(unsigned int iter= 0; iter<maxN; iter+=1 )  
			{corr_out << iter*dt << "\t" << correlation[iter] << std::endl;
			}				
		corr_out.close();	
	
	}
	
	


}
