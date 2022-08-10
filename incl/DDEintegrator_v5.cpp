class DDEintegrator{

  public:
    
    void DDE_euler(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> after_step_func, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    
    void DDE_RK2(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> after_step_func, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    
    void DDE_RK4(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> after_step_func, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    
    
    void DDE_euler(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> after_step_func, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    
    void DDE_RK2(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> after_step_func, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    
    void DDE_RK4(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> after_step_func, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    
    
    
    void DDE_adStpRK45(std::function<void (vars*, vars_vec_wdX_addt*, vars*, parameters*)> derivs, vars_vec_wdX_addt* Xhist, parameters* paras, const double T, const double dt_min, const double dt_max, const double dt_init, std::function<void (vars*, vars_vec_wdX_addt*, parameters*, double, double)> before_step_func, std::function<void (vars*, vars_vec_wdX_addt*, parameters*, double)> after_step_func);
    
    double S = 0.95; //safety factor to not increase the stepsize to much. Keep between 0.9-1.0
    double dt_scale_max_incr = 10; //maximal factor for increasing a step
    double dt_scale_max_decr = 0.2; //mininal factor for decreasing a step
    double beta = 0.08; //used for PI control for the next step; 0 = off; otherwise keep between 0.04 and 0.08
//     double alpha = 1.0/5.0 - beta*0.75; //used for PI control for the next step
    double atol = 0; //absolute error tolerance
    double rtol = 1E-6; //relative error tolerance
    
  private:
    
    
    //RK45 constants for Dormand-Prince method
    const double c2 = 1.0/5.0;
    const double c3 = 3.0/10.0;
    const double c4 = 4.0/5.0;
    const double c5 = 8.0/9.0;
    const double c6 = 1.0;
    const double c7 = 1.0;

    const double a21 = 0.2;
    const double a31 = 3.0/40.0;
    const double a32 = 9.0/40.0;
    const double a41 = 44.0/45.0;
    const double a42 = -56.0/15.0;
    const double a43 = 32.0/9.0;
    const double a51 = 19372.0/6561.0;
    const double a52 = -25360.0/2187.0;
    const double a53 = 64448.0/6561.0;
    const double a54 = -212.0/729.0;
    const double a61 = 9017.0/3168.0;
    const double a62 = -355.0/33.0;
    const double a63 = 46732.0/5247.0;
    const double a64 = 49.0/176.0;
    const double a65 = -5103.0/18656.0;
    const double a71 = 35.0/384.0;
    const double a72 = 0.0;
    const double a73 = 500.0/1113.0;
    const double a74 = 125.0/192.0;
    const double a75 = -2187.0/6784.0;
    const double a76 = 11.0/84.0;

    const double b1 = 35.0/384.0;
    const double b2 = 0.0;
    const double b3 = 500.0/1113.0;
    const double b4 = 125.0/192.0;
    const double b5 = -2187.0/6784.0;
    const double b6 = 11.0/84.0;
    const double b7 = 0.0;

    const double B1 = 5179.0/57600.0;
    const double B2 = 0.0;
    const double B3 = 7571.0/16695.0;
    const double B4 = 393.0/640.0;
    const double B5 = -92097.0/339200.0;
    const double B6 = 187.0/2100.0;
    const double B7 = 1.0/40.0;

    const double e1 = b1 - B1;
    const double e2 = b2 - B2;
    const double e3 = b3 - B3;
    const double e4 = b4 - B4;
    const double e5 = b5 - B5;
    const double e6 = b6 - B6;
    const double e7 = b7 - B7;
    
    //RK45 subroutines
    inline void adStpRK45_calcKi(const double dt, vars *x, vars_vec_wdX_addt* Xhist, parameters *paras, std::function<void (vars*, vars_vec_wdX_addt*, vars*, parameters*)> derivs, vars *k1, vars *k2, vars *k3, vars *k4, vars *k5, vars *k6, vars *k7, vars *xtmp);
    inline void adStpRK45_NextStep(vars *x_old, vars *x_new, vars *k1, vars *k2, vars *k3, vars *k4, vars *k5, vars *k6);
    inline void adStpRK45_calcDelta(vars *delta, vars *k1, vars *k2, vars *k3, vars *k4, vars *k5, vars *k6, vars *k7);
    inline double adStpRK45_calcErr(vars *delta, vars *scale);
    inline void adStpRK45_setScale(vars *x_old, vars *x_new, vars *scale, double atol, double rtol);
    
    
};


/////////////////////
//addaptive stepsize
/////////////////////


void DDEintegrator::DDE_adStpRK45(std::function<void (vars*, vars_vec_wdX_addt*, vars*, parameters*)> derivs, vars_vec_wdX_addt* Xhist, parameters* paras, const double T, const double dt_min, const double dt_max, const double dt_init, std::function<void (vars*, vars_vec_wdX_addt*, parameters*, double, double)> before_step_func, std::function<void (vars*, vars_vec_wdX_addt*, parameters*, double)> after_step_func){
  //output to a file
  std::ofstream outputfile;
  outputfile.open ("data/RK45");
  outputfile.precision(10);

  
  //initialize vars
  vars k1,k2,k3,k4,k5,k6,k7;
  vars xtmp;
  vars delta;
  vars scale;
  
  //initilize variables for stepper algorithm
  double err = 1.0;
  double errold = 1.0;
  double dt = dt_init;
  double dt_scale;
  double ttemp; 
  //set constants for stepper
  double alpha = 1.0/5.0 - beta*0.75; //used for PI control for the next step
  
//   const double S = 0.95; //safety factor to not increase the stepsize to much. Keep between 0.9-1.0
//   const double dt_scale_max_incr = 10; //maximal factor for increasing a step
//   const double dt_scale_max_decr = 0.2; //mininal factor for decreasing a step
//   const double beta = 0.4/5.0; //used for PI control for the next step; 0 = off; otherwise keep between 0.04 and 0.08
//   const double alpha = 1.0/5.0 - beta*0.75; //used for PI control for the next step
//   const double atol = 0; //absolute error tolerance
//   const double rtol = 1E-6; //relative error tolerance
  
  //integrate ODE using Dormand-Prince with adaptive step size
  double TMax = Xhist->t[Xhist->t0] + T;
  while(Xhist->t[Xhist->t0] < TMax){
    //output current state to file
//     outputfile << Xhist->t[Xhist->t0] << "\t" << norm(Xhist->at_t0_rPtr()->E) << "\t" << Xhist->at_t0_rPtr()->G << "\t" << Xhist->at_t0_rPtr()->Q << "\t" << dt << "\t" << err << "\t" << scale.E.real() << "\t" << delta.E.real() << std::endl;
//     outputfile << Xhist->t[Xhist->t0] << "\t" << dt << "\t" << err << std::endl;

    before_step_func(Xhist->at_t0_rPtr(), Xhist, paras, dt, err);
    
    errold = err; //remember previous error for PI control
    ttemp = Xhist->t[Xhist->t0];
    
    adStpRK45_calcKi(dt,Xhist->at_t0_rPtr(), Xhist, paras, derivs, &k1, &k2, &k3, &k4, &k5, &k6, &k7, &xtmp);
    adStpRK45_setScale(Xhist->at_t0_rPtr(), &xtmp, &scale, atol, rtol);
    adStpRK45_calcDelta(&delta, &k1, &k2, &k3, &k4, &k5, &k6, &k7);
    err = adStpRK45_calcErr(&delta, &scale);
    //accept step if err < 1: advance x and calculate dt for the next step
    if(err <= 1.0){
      adStpRK45_NextStep(Xhist->at_t0_rPtr(), Xhist->at_nxt_rPtr(), &k1, &k2, &k3, &k4, &k5, &k6);
//       dt = std::min(S*dt*pow(err,-alpha)*pow(errold,beta), dt_scale_max_incr*dt);
      if(err < 1e-4){
        dt_scale = dt_scale_max_incr;
      }
      else{
        dt_scale = S*pow(err,-alpha)*pow(errold,beta);
        if(dt_scale < dt_scale_max_decr) dt_scale = dt_scale_max_decr;
        if(dt_scale > dt_scale_max_incr) dt_scale = dt_scale_max_incr;
      }
      dt = dt * dt_scale;
      if(dt < dt_min) dt = dt_min;
      if(dt > dt_max) dt = dt_max;
    }
    //reject step if err > 1: decrease stepsize until err < 1 and then advance x and calculate dt for the next step
    else{
      while(err > 1.0 && dt > dt_min){
//         dt = std::max(S*dt*pow(err,-alpha)*pow(errold,beta),dt_scale_max_decr*dt);
//         std::cout << "error too large: " << err << " new dt: " << dt << std::endl;
        dt_scale = S*pow(err,-alpha);
        if(dt_scale < dt_scale_max_decr) dt_scale = dt_scale_max_decr;
        if(dt_scale > dt_scale_max_incr) dt_scale = dt_scale_max_incr;
        dt = dt * dt_scale;
        if(dt < dt_min) dt = dt_min;
        if(dt > dt_max) dt = dt_max;
        
        adStpRK45_calcKi(dt,Xhist->at_t0_rPtr(), Xhist, paras, derivs, &k1, &k2, &k3, &k4, &k5, &k6, &k7, &xtmp);
        adStpRK45_calcDelta(&delta, &k1, &k2, &k3, &k4, &k5, &k6, &k7);
        err = adStpRK45_calcErr(&delta, &scale);
      }
      adStpRK45_NextStep(Xhist->at_t0_rPtr(), Xhist->at_nxt_rPtr(), &k1, &k2, &k3, &k4, &k5, &k6);
//       dt = std::min(S*dt*pow(err,-alpha)*pow(errold,beta),dt_scale_max_incr*dt);
//       if(dt < dt_min) dt = dt_min;
//       if(dt > dt_max) dt = dt_max;
      dt_scale = S*pow(err,-alpha)*pow(errold,beta);
      if(dt_scale < dt_scale_max_decr) dt_scale = dt_scale_max_decr;
      if(dt_scale > dt_scale_max_incr) dt_scale = dt_scale_max_incr;
      dt = dt * dt_scale;
      if(dt < dt_min) dt = dt_min;
      if(dt > dt_max) dt = dt_max;
    }
    
    Xhist->incr_t0();
    Xhist->t[Xhist->t0] = ttemp + dt;
    after_step_func(Xhist->at_t0_rPtr(), Xhist, paras, dt);
  }
  outputfile.close();
}


inline void DDEintegrator::adStpRK45_calcKi(const double dt, vars *x, vars_vec_wdX_addt* Xhist, parameters *paras, std::function<void (vars*, vars_vec_wdX_addt*, vars*, parameters*)> derivs, vars *k1, vars *k2, vars *k3, vars *k4, vars *k5, vars *k6, vars *k7, vars *xtmp){
  //calculates the seven k_i used for the RK formula
//   vars *xtmp = new vars(); //position at which the derivative is evaluated at
  
  //k1
  Xhist->tau_offset = 0.0; //reset tau_offset
  derivs(x,Xhist,k1,paras);
  //write derivative to derivative history vector
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
    ((double*)Xhist->dX_at_t0_rPtr())[i] = ((double*)k1)[i];
  }
  k1->mult(k1,dt); //k1 = dt*f(x,t)
  
  //k2
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
	  ((double*)xtmp)[i] = ((double*)x)[i]
                        + ((double*)k1)[i] * a21;
  }
  Xhist->tau_offset = c2*dt;
  derivs(xtmp,Xhist,k2,paras);
  k2->mult(k2,dt); //k2 = dt*f(x + a21*k1,t + c2*dt) 
   
  //k3
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
	  ((double*)xtmp)[i] = ((double*)x)[i]
                        + ((double*)k1)[i] * a31
                        + ((double*)k2)[i] * a32;
  }
  Xhist->tau_offset = c3*dt;
  derivs(xtmp,Xhist,k3,paras);
  k3->mult(k3,dt); //k3 = dt*f(x + a31*k1 + a32*k2,t + c3*dt) 
  
  //k4
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
	  ((double*)xtmp)[i] = ((double*)x)[i]
                        + ((double*)k1)[i] * a41
                        + ((double*)k2)[i] * a42
                        + ((double*)k3)[i] * a43;
  }
  Xhist->tau_offset = c4*dt;
  derivs(xtmp,Xhist,k4,paras);
  k4->mult(k4,dt);
  
  //k5
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
	  ((double*)xtmp)[i] = ((double*)x)[i]
                        + ((double*)k1)[i] * a51
                        + ((double*)k2)[i] * a52
                        + ((double*)k3)[i] * a53
                        + ((double*)k4)[i] * a54;
  }
  Xhist->tau_offset = c5*dt;
  derivs(xtmp,Xhist,k5,paras);
  k5->mult(k5,dt);
  
  //k6
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
	  ((double*)xtmp)[i] = ((double*)x)[i]
                        + ((double*)k1)[i] * a61
                        + ((double*)k2)[i] * a62
                        + ((double*)k3)[i] * a63
                        + ((double*)k4)[i] * a64
                        + ((double*)k5)[i] * a65;
  }
  Xhist->tau_offset = c6*dt;
  derivs(xtmp,Xhist,k6,paras);
  k6->mult(k6,dt);
   
  //k7
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
	  ((double*)xtmp)[i] = ((double*)x)[i]
                        + ((double*)k1)[i] * a71
//                         + ((double*)k2)[i] * a72 // a72 = 0 anyways
                        + ((double*)k3)[i] * a73
                        + ((double*)k4)[i] * a74
                        + ((double*)k5)[i] * a75
                        + ((double*)k6)[i] * a76;
  }
  Xhist->tau_offset = c7*dt;
  derivs(xtmp,Xhist,k7,paras);
  k7->mult(k7,dt);
  
  
}


inline void DDEintegrator::adStpRK45_NextStep(vars *x_old, vars *x_new, vars *k1, vars *k2, vars *k3, vars *k4, vars *k5, vars *k6){
  //advance x by 5th order RK formula
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
	  ((double*)x_new)[i] = ((double*)x_old)[i]
                        + ((double*)k1)[i] * b1
                        + ((double*)k3)[i] * b3
                        + ((double*)k4)[i] * b4
                        + ((double*)k5)[i] * b5
                        + ((double*)k6)[i] * b6;
  }
  
}

inline void DDEintegrator::adStpRK45_calcDelta(vars *delta, vars *k1, vars *k2, vars *k3, vars *k4, vars *k5, vars *k6, vars *k7){
  //calculate difference vector delta using the k_i and e_i where the e_i are given by b_i - B_i (difference between the coefficients for the 5th order and embedded 4th order formula
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1){
    ((double*)delta)[i] = ((double*)k1)[i] * e1
                        + ((double*)k3)[i] * e3
                        + ((double*)k4)[i] * e4
                        + ((double*)k5)[i] * e5
                        + ((double*)k6)[i] * e6;
//     ((double*)delta)[i] += ((double*)k7)[i] * e7;
  }
}

inline double DDEintegrator::adStpRK45_calcErr(vars *delta, vars *scale){
  //calculates the scaled error err using the euclidian norm.
  double err = 0;
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1){
        err += (((double*)delta)[i] * ((double*)delta)[i]) / (((double*)scale)[i] * ((double*)scale)[i]);
  }
  err /= (sizeof(vars)/sizeof(double));
  err = sqrt(err);
  return err;
}

inline void DDEintegrator::adStpRK45_setScale(vars *x_old, vars *x_new, vars *scale, double atol, double rtol){
  //sets the error tolerances for each variable. scale[i] = atol + max(x_{n}[i], x_{n+1}[i])
  for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1){
    ((double*)scale)[i] = atol + std::max( fabs(((double*)x_old)[i]), fabs(((double*)x_new)[i]) )*rtol;
//     ((double*)scale)[i] = atol + fabs( ((double*)x_old)[i]) * rtol;
  }
}



/////////////////////
// fixed stepsize with cubic delay array interpolation (no dX)
/////////////////////


void DDEintegrator::DDE_euler(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> after_step_func, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for euler method
  vars dXdt;
  dXdt.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 0;
  
  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calculate derivatives at the current position
    derivs(Xhist->at_t0_rPtr(),Xhist, &dXdt, paras);
    dXdt.mult(&dXdt, dt);
    
    //calculate next step
    Xhist->at_t0_rPtr()->add(&dXdt,Xhist->at_nxt_rPtr());
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
    
    //calc noise/etc....
    after_step_func(Xhist->at_t0_rPtr(),Xhist,paras);
        
  }
}

void DDEintegrator::DDE_RK2(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> after_step_func, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for RK2 method
  vars dXdt, k1, k2;
  k1.setTo(0.0), k2.setTo(0.0), dXdt.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 1;
  
  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calculate trial step across dt and save it to k1
    derivs(Xhist->at_t0_rPtr(), Xhist, &dXdt, paras);
    dXdt.mult(&k1, dt/2.0);
    Xhist->at_t0_rPtr()->add(&k1, &k1);
    
    //calculate the step k2 using the trial step k1
    Xhist->tau_offset = 0.5*dt; //adjust tau_offset for the intermediate step
    derivs(&k1, Xhist, &dXdt, paras);
    Xhist->tau_offset = 0.0; //reset tau_offset
    dXdt.mult(&k2,dt);
    
    //calculate the next state using k2 and save it to the history vector
    Xhist->at_t0_rPtr()->add(&k2,Xhist->at_nxt_rPtr());
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
    
    //calc noise/etc....
    after_step_func(Xhist->at_t0_rPtr(),Xhist,paras);
        
  }
}

void DDEintegrator::DDE_RK4(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> after_step_func, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for RK4 method
  vars k1, k2, k3, k4, xtmp;
  k1.setTo(0.0), k2.setTo(0.0), k3.setTo(0.0), k4.setTo(0.0), xtmp.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 3;

  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calculate next step -> results are written to Xhist->at_next = one step in the future
    //calculate RK4 indvividual components ki that add to x as in x_{n+1} = x_{n} + dt*(k'1/6 + k'2/3 + k3/3 + k4/6) where
    
    //first step
    Xhist->tau_offset = 0.0; //reset tau_offset
    derivs(Xhist->at_t0_rPtr(),Xhist,&k1,paras);
    k1.mult(&k1,dt/2.0); //k1 = dt/2*f(x,t,paras) = k'1*0.5
    
    //second step
    Xhist->at_t0_rPtr()->add(&k1,&xtmp); //xtmp = x + k1
    Xhist->tau_offset = 0.5*dt;
    derivs(&xtmp,Xhist,&k2,paras); 
    k2.mult(&k2,dt/2.0); //k2 = dt/2*f(xtmp=x + k1,t + dt/2,paras) = k'2*0.5
    
    //third step
    Xhist->at_t0_rPtr()->add(&k2,&xtmp); //xtmp = x + k2
    derivs(&xtmp,Xhist,&k3,paras);
    k3.mult(&k3,dt); //k3 = dt*f(xtemp = x + k2, t + dt/2,paras)
    
    //fourth step
    Xhist->at_t0_rPtr()->add(&k3,&xtmp); //xtemp = x + k3
    Xhist->tau_offset = dt;
    derivs(&xtmp,Xhist,&k4,paras);
    k4.mult(&k4,dt); //k4 = dt*f(xtemp = x + k3, t + dt,paras)
    //advance x by  k1/3 + k2/1.5 + k3/3 + k4/6
    for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
	  ((double*)Xhist->at_nxt_rPtr())[i] = ((double*)Xhist->at_t0_rPtr())[i]
                                        + ((double*)&k1)[i] / 3.0
                                        + ((double*)&k2)[i] / 1.5
                                        + ((double*)&k3)[i] / 3.0
                                        + ((double*)&k4)[i] / 6.0;
    }
    
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
    
    //calc noise/etc....
    after_step_func(Xhist->at_t0_rPtr(),Xhist,paras);
    
  }
}

/////////////////////
//with vars_vec_wdX
/////////////////////

void DDEintegrator::DDE_euler(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> after_step_func, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for euler method
  vars dXdt;
  dXdt.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 0;
  
  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calculate derivatives at the current position
    derivs(Xhist->at_t0_rPtr(),Xhist, &dXdt, paras);
    dXdt.mult(&dXdt, dt);
    
    //calculate next step
    Xhist->at_t0_rPtr()->add(&dXdt,Xhist->at_nxt_rPtr());
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
    
    //calc noise/etc....
    after_step_func(Xhist->at_t0_rPtr(),Xhist,paras);

  }
}

void DDEintegrator::DDE_RK2(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> after_step_func, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for RK2 method
  vars dXdt, k1, k2;
  k1.setTo(0.0), k2.setTo(0.0), dXdt.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 1;
  
  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calculate trial step across dt and save it to k1
    derivs(Xhist->at_t0_rPtr(), Xhist, &dXdt, paras);
    dXdt.mult(&k1, dt/2.0);
    Xhist->at_t0_rPtr()->add(&k1, &k1);
    
    //calculate the step k2 using the trial step k1
    Xhist->tau_offset = 0.5*dt; //adjust tau_offset for the intermediate step
    derivs(&k1, Xhist, &dXdt, paras);
    Xhist->tau_offset = 0.0; //reset tau_offset
    dXdt.mult(&k2,dt);
    
    //calculate the next state using k2 and save it to the history vector
    Xhist->at_t0_rPtr()->add(&k2,Xhist->at_nxt_rPtr());
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
    
    //calc noise/etc....
    after_step_func(Xhist->at_t0_rPtr(),Xhist,paras);

  }
}


void DDEintegrator::DDE_RK4(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> after_step_func, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for RK4 method
  vars k1, k2, k3, k4, xtmp;
  k1.setTo(0.0), k2.setTo(0.0), k3.setTo(0.0), k4.setTo(0.0), xtmp.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 3;

  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calculate next step -> results are written to Xhist->at_next = one step in the future
    //calculate RK4 indvividual components ki that add to x as in x_{n+1} = x_{n} + dt*(k'1/6 + k'2/3 + k3/3 + k4/6) where
    
    //first step
    derivs(Xhist->at_t0_rPtr(),Xhist,&k1,paras);
    //write derivative to derivative history vector
    for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
      ((double*)Xhist->dX_at_t0_rPtr())[i] = ((double*)&k1)[i];
    }
    //mult k1 with dt/2
    k1.mult(&k1,dt/2.0); //k1 = dt/2*f(x,t,paras) = k'1*0.5
    
    //second step
    Xhist->at_t0_rPtr()->add(&k1,&xtmp); //xtmp = x + k1
    Xhist->tau_offset = 0.5*dt;
    derivs(&xtmp,Xhist,&k2,paras); 
    k2.mult(&k2,dt/2.0); //k2 = dt/2*f(xtmp=x + k1,t + dt/2,paras) = k'2*0.5
    
    //third step
    Xhist->at_t0_rPtr()->add(&k2,&xtmp); //xtmp = x + k2
    derivs(&xtmp,Xhist,&k3,paras);
    k3.mult(&k3,dt); //k3 = dt*f(xtemp = x + k2, t + dt/2,paras)
    
    //fourth step
    Xhist->at_t0_rPtr()->add(&k3,&xtmp); //xtemp = x + k3
    Xhist->tau_offset = dt;
    derivs(&xtmp,Xhist,&k4,paras);
    k4.mult(&k4,dt); //k4 = dt*f(xtemp = x + k3, t + dt,paras)
    //advance x by  k1/3 + k2/1.5 + k3/3 + k4/6
    for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
      ((double*)Xhist->at_nxt_rPtr())[i] = ((double*)Xhist->at_t0_rPtr())[i]
                                        + ((double*)&k1)[i] / 3.0
                                        + ((double*)&k2)[i] / 1.5
                                        + ((double*)&k3)[i] / 3.0
                                        + ((double*)&k4)[i] / 6.0;
    }
    
    Xhist->tau_offset = 0.0; //reset tau_offset
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn += 1;
    
    //calc noise/etc....
    after_step_func(Xhist->at_t0_rPtr(),Xhist,paras);
    
  }
}






auto outputfunc_empty = [](vars* X, unsigned long long int tn, unsigned long long int tn_final){};

auto after_step_empty = [](vars* X, vars_vec* Xhist, parameters* p){};
auto after_step_empty_wdX = [](vars* X, vars_vec_wdX* Xhist, parameters* p){};
auto after_step_empty_adSt = [](vars* X, vars_vec_wdX_addt* Xhist, parameters *p, double dt){};


