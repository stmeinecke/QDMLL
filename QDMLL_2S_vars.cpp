class vars{ 
	 public:
    std::complex<double> A_out;
    std::complex<double> A_plus[NSEC];
		std::complex<double> A_minus[NSEC];
    std::complex<double> P_plus[NSEC];
    std::complex<double> P_minus[NSEC];
    std::complex<double> P_ES[NSEC];
    double N[NSEC];
		double rho_ES[NSEC][NSAMPLES];
		double rho_GS[NSEC][NSAMPLES];
    std::complex<double> G_plus[NSEC][NSAMPLES];
		std::complex<double> G_minus[NSEC][NSAMPLES];

		vars(){this->setTo(0.0);}

		void add(vars *v, vars *r){
			for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)r)[i] = ((double*)this)[i] + ((double*)v)[i];
			}
		} 

		void mult(vars *r, double m){
			for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)r)[i] = ((double*)this)[i]*m ;
			}
		}

		void setTo(double s){
			for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)this)[i] = s;
			}
		}

};

namespace IND{

  unsigned int A_out;
  unsigned int A_plus[NSEC]; 
  unsigned int A_minus[NSEC];
  unsigned int P_plus[NSEC];
  unsigned int P_minus[NSEC];
  unsigned int P_ES[NSEC];
//   unsigned int N[NSEC];
//   unsigned int rho_ES[NSEC];
//   unsigned int rho_GS[NSEC];
//   unsigned int G_plus[NSEC];
//   unsigned int G_minus[NSEC];

  void initIndices(){
    unsigned int k = 0;
    
    A_out = k;
    k += 2;
    for(int l = 0; l < NSEC; l++){
      A_plus[l] = k;
      k += 2;
    }
    for(int l = 0; l < NSEC; l++){
      A_minus[l] = k;
      k += 2;
    }
    for(int l = 0; l < NSEC; l++){
      P_plus[l] = k;
      k += 2;
    }
    for(int l = 0; l < NSEC; l++){
      P_minus[l] = k;
      k += 2;
    }
    for(int l = 0; l < NSEC; l++){
      P_ES[l] = k;
      k += 2;
    }
//     for(int l = 0; l < NSEC; l++){
//       N[l] = k;
//       k += 1;
//     }
//     for(int l = 0; l < NSEC; l++){
//       rho_ES[l] = k;
//       k += 1;
//     }
//     for(int l = 0; l < NSEC; l++){
//       rho_GS[l] = k;
//       k += 1;
//     }
//     for(int l = 0; l < NSEC; l++){
//       G_plus[l] = k;
//       k += 2;
//     }
//     for(int l = 0; l < NSEC; l++){
//       G_minus[l] = k;
//       k += 2;
//     }

    
  }
  
}
