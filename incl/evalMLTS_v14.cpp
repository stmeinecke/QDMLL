class evalMLTS {
public:
  
  // basic parameters of the ML time series
  double T0, dt;
  
  // output of the basics () function
  double TSAverage, TSGreatestMax , TSSmallestMin;
  double meanMaxTSep, meanMax; //
  // output of the TriggerValPulseDetection_Statistics function called by the basics() function
  double TVPD_meanPulseSep, TVPD_meanPulseMax, TVPD_meanPulseArea, TVPD_meanPulseSD;
  std::vector<double> TVPD_PS_vec, TVPD_PA_vec, TVPD_PM_vec, TVPD_PSD_vec; //vectors
  // output of pulswidth measurment algorithms
  double PW_mean_PFWHM, PW_SD_PFWHM;
  double PW_mean_PSD, PW_SD_PSD;
  
  ///////
  // thresholds and triggers for detection and evaluation algorithms
  ///////
  int TSSmoothing = 0; // 0 = no smoothing; 1 = 'standard' smoothing; 2 = cosine smoothing
  double rel_cos_sm_width = 0.05; //smoothing width in units of T0
  ///basics
  double extremaTol = 1E-16; // relative threshhold to identify extrema by comparison to neighboring points
  double doubleCountTol = 0.01; // relative numerical tolerance to distinguish between unique extrema
  double max_min_dscr = 1E-5; //used as a cut-off for detecting cw emission by relativ difference between greatest max and smallest min
  double av_thresh = 1E-50; // threshhold to for the mean intensity
  double pulse_rel_max_th = 0.01; //minimal value for a pulse maximum to be accepted, i.e. included in maxima, in units of the greatest max
  //pulse width analysis
  double PW_rel_searchRadius = 0.3; //search radius for pulse width measurements in units of T0
  //weighted pulse position detection
  double TVPD_rel_vector_tBounds = 0.3; // boundaries at the beginning and the end of the array to avoid artifacts in case pulses are at the boundaries -> dont consider these part of the array. In units of T0
  double TVPD_rel_pulse_trig = 0.01; //triggers a pulse detection and subsequently calculates the center of mass of the pulse. In units of the greatest max.
  double TVPD_rel_pulse_bounds = 1E-3; //pulse center of mass withing theses bounds in units of the greatest max
  
  ///////
  // data and stuff
  ///////
  std::vector<double>* raw_DataVec_Ptr; //pointer to the input data
  std::vector<double>* sm_DataVec_Ptr = new std::vector<double>; //pointer to a smoothed input data vector
  std::vector<double>* DataVec_Ptr; // pointer to the vector that is to be evaluated
  std::vector<double> IntAC; // vector containing the intensity auto-correlation function of the input data
  
  std::vector<double> maxima; //maxima
  std::vector<double> maxima_tvec; //maxima positions
  std::vector<double> minima; //minima
  std::vector<double> minima_tvec; //minima positions
  std::vector<double> uniqueMax; // contains the unique maxima
  std::vector<double> uniqueMin; // contains the unique minima
  std::vector<double> uniqueMLMax; // contains the unique maxima
  
  
  ///////
  // constructor
  ///////
  evalMLTS();
  evalMLTS(std::vector<double>* DataVec_Ptr, double dt, double T0);
  

  ///////
  // methods
  ///////
  void loadNewTS(std::vector<double>* DataVec_Ptr, double dt, double T0); //resets vectors and variables, loads new time series and runs bascis() method
  
  
  std::tuple<double,double> findMaxima(std::vector<double> &max_vec, std::vector<double> &max_tvec, double rel_max_th);
  void findMinima(std::vector<double> &min_vec, std::vector<double> &min_tvec);
  int numberOfMax(void); //counts unique maxima
  int numberOfMin(void); //counts unique minima
  int MLnMax(void); //counts unique ML maxima
  
  std::tuple<double,double> PW_FWHM_fMaxP(void); //calculates mean FWHM pulse width from all pulses (and stdDev of that)
  std::tuple<double,double> PW_StdDev_fMaxP(void);  //calculates mean pulse stdDev as a measure of the pulse width from all pulses (and stdDev of that)
  
  double PtPJ_MaxPos();
  double AmpJitter();
  double RelAmpJitter();

  std::tuple<double,double,double,double> TriggerValPulseDetection_Statistics(std::vector<double> &PPos, std::vector<double> &PArea, std::vector<double> &PPeakPow, std::vector<double> &PWidSD);
  std::tuple<double,double,double,double> TriggerValPulseDetection_Statistics(std::vector<double> &PPos, std::vector<double> &PArea, std::vector<double> &PPeakPow, std::vector<double> &PWidSD, double trig_rel_th, double pulse_rel_bounds);
  std::tuple<double,double,double,double> TriggerValPulseDetection_Statistics(double trig_rel_th, double pulse_rel_bounds);
  std::tuple<double,double,double,double,double,double> TriggerValPulseDetection_Statistics_forEval(double trig_rel_th, double pulse_rel_bounds);
  
  double Weighted_PulsePositions(std::vector<double> &PP);
  double Weighted_PulsePositions(std::vector<double> &PP, double trig_rel_th, double pulse_rel_bounds);
  double PtPJ_WPP();
  double PtPJ_WPP(double trig_rel_th);
  double PtPJ_WPP(double trig_rel_th, double pulse_rel_bounds);
  
  double ML_FWHM_AC(std::vector<std::complex<double>> &outVC, double dt);
  double ML_FWHM_IntAC(void);
  double ML_StdDev_IntAC(void);
  int IntAC_TEP(double integrationRadius);
  std::tuple<double,double,double,double> EvalIntAC();
  int IntAC_MaxAboveTh(std::vector<double> &IACMaxPos, double max_thresh);
  void dump_IACmaxpos(double thresh);
  
  double ML_Fund_RR_from_AC(std::vector<std::complex<double>> &outVC, double dt);
  double ML_Fund_RR_from_IntAC(void);
  
  double opticalFreq(std::vector<std::complex<double>> &ComplexE, std::vector<double> &Phase, std::vector<double> &Frequency, double dt);
  void smoothTS(std::vector<double> &vec, std::vector<double> &out);
  void smoothTS_cos(std::vector<double> &in, std::vector<double> &out, double width);

  double Variance(std::vector<double> &vec);
  double PtPJitter(std::vector<double> &vec);
  
private:
  void basics(); //calculates basic figures of merit
  
  double PW_searchRadius; //absolut search radius for pulse width measurements
  double WPP_vector_tBounds; //absolut array bounds for weighted pulse position detection
  
  double firstMaxTPos, lastMaxTPos; //first an last max position for meanMaxTSep
  
  double tmpExtrm, tmpExtrm_t; //tmp Extremum and its position from quadratic interpolations
  std::tuple<double,double> interpolateQuadExtrm(const double X1, const double f1, const double X2, const double f2, const double X3, const double f3); //quadratic interpolation of an extremum and its position using its neighboring points
};
 
evalMLTS::evalMLTS(){};

evalMLTS::evalMLTS(std::vector<double>* in_DataVec_Ptr, double dt, double T0){
  if(in_DataVec_Ptr->size() < 2){
    std::cout << "error: input data vector too short" << std::endl;
  }
  else{
    this->raw_DataVec_Ptr = in_DataVec_Ptr;
    this->dt = dt;
    this->T0 = T0;
    
    switch(TSSmoothing){
      case 0 : {
        this->DataVec_Ptr = this->raw_DataVec_Ptr;
        break;
      }
      case 1 : {
        this->smoothTS(*raw_DataVec_Ptr, *sm_DataVec_Ptr);
        this->DataVec_Ptr = this->sm_DataVec_Ptr;
        break;
      }
      case 2 : {
        double SW = this->rel_cos_sm_width * this->T0;
        this->smoothTS_cos(*raw_DataVec_Ptr, *sm_DataVec_Ptr,SW);
        this->DataVec_Ptr = this->sm_DataVec_Ptr;
        break;
      }
    }
    
    basics();
  }
}

void evalMLTS::loadNewTS(std::vector<double>* in_DataVec_Ptr, double dt, double T0){
  if(in_DataVec_Ptr->size() < 2){
    std::cout << "error: input data vector too short" << std::endl;
  }
  else{
    this->raw_DataVec_Ptr = in_DataVec_Ptr;
    this->dt = dt;
    this->T0 = T0;
    
    switch(TSSmoothing){
      case 0 : {
        this->DataVec_Ptr = this->raw_DataVec_Ptr;
        break;
      }
      case 1 : {
        this->smoothTS(*raw_DataVec_Ptr, *sm_DataVec_Ptr);
        this->DataVec_Ptr = this->sm_DataVec_Ptr;
        break;
      }
      case 2 : {
        double SW = this->rel_cos_sm_width * this->T0;
        this->smoothTS_cos(*raw_DataVec_Ptr, *sm_DataVec_Ptr,SW);
        this->DataVec_Ptr = this->sm_DataVec_Ptr;
        break;
      }
    }
    
    basics();
  }
}

void evalMLTS::basics(){
  this->PW_searchRadius = this->PW_rel_searchRadius * this->T0;
  this->WPP_vector_tBounds = this->TVPD_rel_vector_tBounds * this->T0;
  
  maxima.resize(0);
  maxima_tvec.resize(0);
  minima.resize(0);
  minima_tvec.resize(0);
  
  uniqueMax.resize(0);
  uniqueMin.resize(0);
  uniqueMLMax.resize(0);

  IntAC.resize(0);
  
  TSAverage = 0;
  TSGreatestMax = 0;
  TSSmallestMin = 0;
  
  meanMaxTSep = 0;
  meanMax = 0;
  
  //find smalles minimum and greatest maximum
  double DataVec_Ptr_sum = 0;
  TSSmallestMin = DataVec_Ptr->at(0);
  TSGreatestMax = DataVec_Ptr->at(0);
  for(std::size_t k=0; k < DataVec_Ptr->size()-2; k++){
    //calc DataVec_Ptr_sum for TSAverage
    DataVec_Ptr_sum += DataVec_Ptr->at(k);
    //find greates maximum
    if(TSGreatestMax < DataVec_Ptr->at(k)){
      TSGreatestMax = DataVec_Ptr->at(k);
    }
    //find smallest mininum
    if(TSSmallestMin > DataVec_Ptr->at(k)){
      TSSmallestMin = DataVec_Ptr->at(k);
    }
  }
  //calculate TSAverage
  TSAverage = DataVec_Ptr_sum / (DataVec_Ptr->size()-2);

  //find maxima -> write to maxima and maxima_tvec; calc mean maximum and mean maximum time seperation
  std::tie(meanMax,meanMaxTSep) = findMaxima(maxima, maxima_tvec, pulse_rel_max_th);
  
  //calculate intensity auto-correlation
  katana::get_autocorr(*raw_DataVec_Ptr, IntAC);
  
  
  //calculate pulse train statisctics using weighted pulse positions and trigger detection
  TVPD_meanPulseSep=0, TVPD_meanPulseMax=0, TVPD_meanPulseArea=0, TVPD_meanPulseSD=0;
  TVPD_PS_vec.resize(0), TVPD_PA_vec.resize(0), TVPD_PM_vec.resize(0), TVPD_PSD_vec.resize(0);
  std::tie(TVPD_meanPulseSep, TVPD_meanPulseArea, TVPD_meanPulseMax, TVPD_meanPulseSD) = TriggerValPulseDetection_Statistics(TVPD_PS_vec,TVPD_PA_vec,TVPD_PM_vec,TVPD_PSD_vec);
  
  //calculate pulse pulse width
  PW_mean_PFWHM=0, PW_SD_PFWHM=0;
  PW_mean_PSD=0, PW_SD_PSD=0;
  std::tie(PW_mean_PFWHM, PW_SD_PFWHM) = PW_FWHM_fMaxP();
  std::tie(PW_mean_PSD, PW_SD_PSD) = PW_StdDev_fMaxP();
  
}

///////////////////////////////////////////////////////////////////
// extrema finding stuff
///////////////////////////////////////////////////////////////////

std::tuple<double,double> evalMLTS::findMaxima(std::vector<double> &max_vec, std::vector<double> &max_tvec, double rel_max_th){
  double max_th = rel_max_th * this->TSGreatestMax; //set ML pulse threshhold
  double mMax = 0; //mean maximum
  //find extrema
  for(std::size_t k=0; k < DataVec_Ptr->size()-2; k++){  
    if(DataVec_Ptr->at(k) < (1.0-extremaTol) * DataVec_Ptr->at(k+1) && (1.0-extremaTol)*DataVec_Ptr->at(k+1) > DataVec_Ptr->at(k+2) && DataVec_Ptr->at(k+1) > max_th){
      if(max_vec.size() == 0) firstMaxTPos = dt*(k+1);
      else lastMaxTPos = dt*(k+1);
      //interpolate true maximum from with next neighbors 
      std::tie(tmpExtrm, tmpExtrm_t) = interpolateQuadExtrm(dt*(k),DataVec_Ptr->at(k),dt*(k+1),DataVec_Ptr->at(k+1),dt*(k+2),DataVec_Ptr->at(k+2));
      max_vec.push_back(tmpExtrm);
      max_tvec.push_back(tmpExtrm_t);
      mMax += tmpExtrm;
    }
  }
  //calculate TSAverage distance between maxima
  double mMaxTSep = 0;
  if(max_vec.size() > 1 && TSGreatestMax/TSSmallestMin > (1.0 + max_min_dscr) && TSAverage > av_thresh) mMaxTSep = (lastMaxTPos - firstMaxTPos)/((int)max_vec.size()-1);
  else mMaxTSep = -1;
  //calculate TSAverage maximum
  if(max_vec.size() > 1 && TSGreatestMax/TSSmallestMin > (1.0 + max_min_dscr) && TSAverage > av_thresh) mMax = mMax/max_vec.size();
  else mMax = TSAverage;
  
  return std::make_tuple(mMax, mMaxTSep);
}

void evalMLTS::findMinima(std::vector<double> &min_vec, std::vector<double> &min_tvec){
  //find extrema
  for(std::size_t k=0; k < DataVec_Ptr->size()-2; k++){  
    //find minima
    if(DataVec_Ptr->at(k) > (1.0+extremaTol) * DataVec_Ptr->at(k+1) && (1.0+extremaTol)*DataVec_Ptr->at(k+1) < DataVec_Ptr->at(k+2)){
      //interpolate true maximum from with next neighbors
      std::tie(tmpExtrm, tmpExtrm_t) = interpolateQuadExtrm(dt*(k),DataVec_Ptr->at(k),dt*(k+1),DataVec_Ptr->at(k+1),dt*(k+2),DataVec_Ptr->at(k+2));
      minima.push_back(tmpExtrm);
      minima_tvec.push_back(tmpExtrm_t);
    }
  }
}

int evalMLTS::numberOfMax(){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return 0;
  else{
    int max_counter = 0;
    bool doublecount;
    for(std::size_t k = 0; k < maxima.size(); k++){
      doublecount = false;
      for(int l = 0; l < max_counter; l++){
        if(fabs(maxima[k]-TSAverage) > (1.0-doubleCountTol) * fabs(uniqueMax[l]-TSAverage) && fabs(maxima[k]-TSAverage) < (1.0+doubleCountTol) * fabs(uniqueMax[l]-TSAverage)){
          doublecount = true;
        }
      }
      if(doublecount == false){
        uniqueMax.push_back(maxima[k]);
        max_counter += 1;
      }
    }
    return max_counter;
  }
}

int evalMLTS::numberOfMin(void){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return 0;
  else{
    int min_counter = 0;
    bool doublecount;
    for(std::size_t k = 0; k < minima.size(); k++){
      doublecount = false;
      for(int l = 0; l < min_counter; l++){
        if(fabs(minima[k]-TSAverage) > (1.0-doubleCountTol) * fabs(uniqueMin[l]-TSAverage) && fabs(minima[k]-TSAverage) < (1.0+doubleCountTol) * fabs(uniqueMin[l]-TSAverage)){
          doublecount = true;
        }
      }
      if(doublecount == false){
        uniqueMin.push_back(minima[k]);
        min_counter += 1;
      }
    }
    return min_counter;
  }
}

int evalMLTS::MLnMax(void){
  //returns number of unique ML maxima
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return 0;
  else{
    int max_counter = 0;
    bool doublecount;
    for(std::size_t k = 0; k < maxima.size(); k++){
      doublecount = false;
      for(int l = 0; l < max_counter; l++){
        if(fabs(maxima[k]-TSAverage) > (1.0-doubleCountTol) * fabs(uniqueMLMax[l]-TSAverage) && fabs(maxima[k]-TSAverage) < (1.0+doubleCountTol) * fabs(uniqueMLMax[l]-TSAverage)){
          doublecount = true;
        }
      }
      if(doublecount == false){
        uniqueMLMax.push_back(maxima[k]);
        max_counter += 1;
      }
    }
    return max_counter;

  }
}

///////////////////////////////////////////////////////////////////
// pulse and pulse train properties and statistics
///////////////////////////////////////////////////////////////////

std::tuple<double,double> evalMLTS::PW_FWHM_fMaxP(void){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return std::make_tuple(-1.0,-1.0);
  else{
    std::vector<double> PW;
    for(std::size_t k = 0; k < maxima_tvec.size(); k++){
      if(maxima_tvec[k] > this->PW_searchRadius && maxima_tvec[k] < this->DataVec_Ptr->size()*this->dt - this->PW_searchRadius){
        int maxIndex = (int)(maxima_tvec[k]/this->dt);
        double halfMax = maxima[k]/2.0;
        int kk = maxIndex;
        while(kk < maxIndex + (int)(this->PW_searchRadius/this->dt) && (*DataVec_Ptr)[kk] > halfMax) kk++;
        double pulseEnd = kk*this->dt + this->dt * ((*DataVec_Ptr)[kk]-halfMax)/((*DataVec_Ptr)[kk]-(*DataVec_Ptr)[kk+1]);
        kk = maxIndex;
        while(kk > maxIndex - (int)(this->PW_searchRadius/this->dt) && (*DataVec_Ptr)[kk] > halfMax) kk--;
        double pulseStart = kk*this->dt - this->dt * ((*DataVec_Ptr)[kk]-halfMax)/((*DataVec_Ptr)[kk]-(*DataVec_Ptr)[kk-1]);
        
        PW.push_back(pulseEnd - pulseStart);
      }
    }
    double avPW = 0.0;
    for(std::size_t k = 0; k < PW.size(); k++) avPW += PW[k]/PW.size();
    
    double var = 0.0;
    for(std::size_t k = 0; k < PW.size(); k++) var += (avPW - PW[k]) * (avPW - PW[k]) / PW.size();
    
    return std::make_tuple(avPW, sqrt(var));
  }
}

std::tuple<double,double> evalMLTS::PW_StdDev_fMaxP(void){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return std::make_tuple(-1.0,-1.0);
  else{
    std::vector<double> Pmean, Pvar;
    for(std::size_t k = 0; k < maxima_tvec.size(); k++){
      if(maxima_tvec[k] > this->PW_searchRadius && maxima_tvec[k] < this->DataVec_Ptr->size()*this->dt -this->PW_searchRadius){
        int maxIndex = (int)(maxima_tvec[k]/this->dt);
        
        double mean = 0.0;
        double norm = 0.0;
        for(int kk = maxIndex-(int)(this->PW_searchRadius/this->dt); kk < maxIndex+(int)(this->PW_searchRadius/this->dt); kk++){
          norm += (*DataVec_Ptr)[kk];
          mean += (*DataVec_Ptr)[kk] * kk * this->dt;
        }
        mean = mean / norm;
        
        double var = 0.0;
        for(int kk = maxIndex-(int)(this->PW_searchRadius/dt); kk < maxIndex+(int)(this->PW_searchRadius/dt); kk++){
          var += (mean - kk*this->dt)*(mean - kk*this->dt) * (*DataVec_Ptr)[kk] / norm;
        }
        
        Pmean.push_back(mean);
        Pvar.push_back(var);
	
      }
    }
    double PstdDev_mean = 0.0;
    for(std::size_t k = 0; k < Pvar.size(); k++) PstdDev_mean += sqrt(Pvar[k]) / Pvar.size();
    
    double PstdDev_var = 0.0;
    for(std::size_t k = 0; k < Pvar.size(); k++) PstdDev_var += (sqrt(Pvar[k]) - PstdDev_mean) * (sqrt(Pvar[k]) - PstdDev_mean) / Pvar.size();
    
    return std::make_tuple(PstdDev_mean, sqrt(PstdDev_var));
  } 
}

std::tuple<double,double,double,double> evalMLTS::TriggerValPulseDetection_Statistics(std::vector<double> &PPos, std::vector<double> &PArea, std::vector<double> &PPeakPow, std::vector<double> &PWidSD){
  return TriggerValPulseDetection_Statistics(this->TVPD_PS_vec,this->TVPD_PA_vec,this->TVPD_PM_vec,this->TVPD_PSD_vec, this->TVPD_rel_pulse_trig, this->TVPD_rel_pulse_bounds);
}

std::tuple<double,double,double,double> evalMLTS::TriggerValPulseDetection_Statistics(std::vector<double> &PPos, std::vector<double> &PArea, std::vector<double> &PPeakPow, std::vector<double> &PWidSD, double trig_rel_th, double pulse_rel_bounds){
  //only use for FML or HML
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return std::make_tuple(-1,0,0,-1);
  else{
    PPos.reserve(1000),PArea.reserve(1000),PPeakPow.reserve(1000),PWidSD.reserve(1000);
    PPos.resize(0),PArea.resize(0),PPeakPow.resize(0),PWidSD.resize(0);
    
    double Inorm = 0.0;
    double Imax = 0.0;
    bool pulse = false;
    int lowerPbound = 0;
    int upperPbound = 0;
    double firstMoment = 0.0;
    double secondMoment = 0.0;
    
    for(std::size_t k = (int)(WPP_vector_tBounds/this->dt); k < (*DataVec_Ptr).size()-(int)(WPP_vector_tBounds/this->dt); k++){
      
      //detect pulses and calc norm/pulse area
      if(pulse == false && (*DataVec_Ptr)[k] > this->TSGreatestMax * trig_rel_th){
        pulse = true;
        int wk = k-1;
        while((*DataVec_Ptr)[wk] > this->TSGreatestMax * pulse_rel_bounds && (k-wk)*this->dt < WPP_vector_tBounds/2.0){
          Inorm += (*DataVec_Ptr)[wk];
          wk--;
        }
        lowerPbound = wk + 1;
      }
      if(pulse == true && (*DataVec_Ptr)[k] > this->TSGreatestMax * pulse_rel_bounds ){
        Inorm += (*DataVec_Ptr)[k];
        upperPbound = k;
      }
      
      //calc stuff
      if(pulse == true && (*DataVec_Ptr)[k] < this->TSGreatestMax * pulse_rel_bounds ){
        pulse = false;
        for(int pk = lowerPbound; pk <= upperPbound; pk++){
          firstMoment += ((*DataVec_Ptr)[pk] / Inorm) * pk * this->dt;
          secondMoment += ((*DataVec_Ptr)[pk] / Inorm) * (pk * this->dt) * (pk * this->dt);
          if(Imax < (*DataVec_Ptr)[pk]) Imax = (*DataVec_Ptr)[pk];
        }
        PPos.push_back(firstMoment);
        PArea.push_back(Inorm*this->dt);
        PPeakPow.push_back(Imax);
        PWidSD.push_back(sqrt(secondMoment - (firstMoment * firstMoment)));
        Inorm = 0.0;
        firstMoment = 0.0;
        secondMoment = 0.0;
        Imax = 0.0;
      }
    }
    
    double avPPosMaxDist = 0.0;
    double avPArea = 0.0;
    double avPPeakPow = 0.0;
    double avPWidSD = 0.0;
    if(PPos.size() > 1) for(int k = 0; k < (int)PPos.size()-1; k++) avPPosMaxDist += (PPos[k+1] - PPos[k])/(PPos.size()-1);
    else avPPosMaxDist = 0.0;
    
    for(std::size_t k = 0; k < PArea.size(); k++) avPArea += PArea[k]/PArea.size();
    for(std::size_t k = 0; k < PPeakPow.size(); k++) avPPeakPow += PPeakPow[k]/PPeakPow.size();
    for(std::size_t k = 0; k < PWidSD.size(); k++) avPWidSD += PWidSD[k]/PWidSD.size();
  

    return std::make_tuple(avPPosMaxDist, avPArea, avPPeakPow, avPWidSD);
  }
}

std::tuple<double,double,double,double> evalMLTS::TriggerValPulseDetection_Statistics(double trig_rel_th, double pulse_rel_bounds){
  std::vector<double> PS_vec, PA_vec, PM_vec, PSD_vec; //vectors
  return TriggerValPulseDetection_Statistics(PS_vec,PA_vec,PM_vec,PSD_vec, trig_rel_th, pulse_rel_bounds);
}

std::tuple<double,double,double,double,double,double> evalMLTS::TriggerValPulseDetection_Statistics_forEval(double trig_rel_th, double pulse_rel_bounds){
  std::vector<double> PS_vec, PA_vec, PM_vec, PSD_vec;
  double mS,mA,mM,mSD;
  std::tie(mS, mA, mM, mSD) = TriggerValPulseDetection_Statistics(PS_vec,PA_vec,PM_vec,PSD_vec, trig_rel_th, pulse_rel_bounds);
  return std::make_tuple(mS, this->PtPJitter(PS_vec), mM, sqrt(this->Variance(PM_vec))/mM, mA, sqrt(this->Variance(PA_vec))/mA);
}

double evalMLTS::PtPJitter(std::vector<double> &vec){
  if(vec.size() < 3) return -1;
  else{
    double mean = 0.0;
    double var = 0.0;
    for(int k = 0; k < (int)vec.size()-1; k++) mean += (vec[k+1] - vec[k])/(vec.size()-1);
    for(int k = 0; k < (int)vec.size()-1; k++) var += (mean - (vec[k+1] - vec[k])) * (mean - (vec[k+1] - vec[k])) / (vec.size()-1);
    return sqrt(var);
  }
}

double evalMLTS::PtPJ_MaxPos(){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    return this->PtPJitter(this->maxima_tvec);
  }
}


double evalMLTS::AmpJitter(){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    return sqrt(this->Variance(this->maxima));
  }
}

double evalMLTS::RelAmpJitter(){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    double mean = 0.0;
    for(std::size_t k = 0; k < maxima.size(); k++) mean += maxima[k] / maxima.size();
    
    double var = 0.0;
    for(std::size_t k = 0; k < maxima.size(); k++) var += (mean - maxima[k]) * (mean - maxima[k]) / maxima.size();
    
    return sqrt(var)/mean;
  }
}

///////////////////////////////////////////////////////////////////
// auto-correlation based stuff
///////////////////////////////////////////////////////////////////

double evalMLTS::ML_Fund_RR_from_AC(std::vector<std::complex<double>> &outVC, double dt){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    std::vector<double> AC;
    katana::get_autocorr_norm_out(outVC, AC);
    
    int T0_ind = (int)(this->T0/this->dt);
    double relT0sr = 0.1; //search radius
    
    if( AC.size() > (unsigned int)((1.0+relT0sr)*T0/this->dt) ){
      double tmpExt, tmpExt_t;
      double itmax = AC[T0_ind];
      int itk = T0_ind;
      ;
      for(int k = T0_ind - (int)(relT0sr*T0/this->dt); k < (T0_ind + (int)(relT0sr*T0/this->dt)); k++){
        if(AC[k] > itmax){
          itmax = AC[k];
          itk = k;
        }
      }
      std::tie(tmpExt, tmpExt_t) = interpolateQuadExtrm(dt*(itk-1),AC[itk-1],dt*(itk),AC[itk],dt*(itk+1),AC[itk+1]);
      return tmpExt_t;
      
    }
    return -1;
  }
}

double evalMLTS::ML_Fund_RR_from_IntAC(void){  
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    int T0_ind = (int)(this->T0/this->dt);
    double relT0sr = 0.1; //search radius
    
    if( this->IntAC.size() > (unsigned int)((1.0+relT0sr)*T0/this->dt) ){
      double tmpExt, tmpExt_t;
      double itmax = this->IntAC[T0_ind];
      int itk = T0_ind;
      ;
      for(int k = T0_ind - (int)(relT0sr*T0/this->dt); k < (T0_ind + (int)(relT0sr*T0/this->dt)); k++){
        if(this->IntAC[k] > itmax){
          itmax = this->IntAC[k];
          itk = k;
        }
      }
      std::tie(tmpExt, tmpExt_t) = interpolateQuadExtrm(dt*(itk-1),this->IntAC[itk-1],dt*(itk),this->IntAC[itk],dt*(itk+1),this->IntAC[itk+1]);
      return tmpExt_t;
      
    }
    return -1;
  }
}

double evalMLTS::ML_FWHM_AC(std::vector<std::complex<double>> &outVC, double dt){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    std::vector<double> AC;
    katana::get_autocorr_norm_out(outVC, AC);
    
    int k = 0;
    while(AC[k+1] > AC[0]/2.0 && k < (int)AC.size()-1) k++; //Half max (AC[0]/2.0) is between k and k+1
    
    if(dt*k > this->PW_searchRadius) return -1; //if pulses are that long, they are not mode-locked pulses or something has gone wrong
    return sqrt(2) * (dt*k + dt * (AC[k]-AC[0]/2.0)/(AC[k]-AC[k+1])); //calc half-width and multiply by 2/sqrt(2) = sqrt(2) to account for FW and AC 'broadening'; including linear interpolation of the 'true' (time-step independent) Half-Max position
  }
}

double evalMLTS::ML_FWHM_IntAC(void){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    
    int k = 0;
    while(this->IntAC[k+1] > this->IntAC[0]/2.0 && k < (int)this->IntAC.size()-1) k++; //Half max (AC[0]/2.0) is between k and k+1
    
    if(this->dt*k > this->PW_searchRadius) return -1; //if pulses are that long, they are not mode-locked pulses or something has gone wrong
    return sqrt(2) * (this->dt*k + this->dt * (this->IntAC[k]-this->IntAC[0]/2.0)/(this->IntAC[k]-this->IntAC[k+1])); //calc half-width and multiply by 2/sqrt(2) = sqrt(2) to account for FW and AC 'broadening'; including linear interpolation of the 'true' (time-step independent) Half-Max position
  }
}

double evalMLTS::ML_StdDev_IntAC(void){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{

    double IntgrAC = 0;
    for(int k = 0; k*this->dt < this->PW_searchRadius; k++) IntgrAC += this->IntAC[k];
    
    double Var = 0;
    for(int k = 0; k*this->dt < this->PW_searchRadius; k++) Var += this->IntAC[k] * (k*this->dt)*(k*this->dt) / IntgrAC;
    
    return sqrt(2) * sqrt(Var);
  }
}

int evalMLTS::IntAC_TEP(double integrationRadius){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{

    int maxcounter = 0;
    for(int k = 0; k*this->dt < integrationRadius; k++){
      if(this->IntAC[k] < (1.0-extremaTol) * this->IntAC[k+1] && (1.0-extremaTol)*this->IntAC[k+1] > this->IntAC[k+2] && this->IntAC[k+1] > this->IntAC[0]*1E-3) maxcounter++;
    }
    
    return maxcounter;
  }
}


std::tuple<double,double,double,double> evalMLTS::EvalIntAC(){
  unsigned int maxInd = (unsigned int)(1.1*this->T0/this->dt);
  double firstMax = -1;
  double firstMaxPos = -1;
  double firstMin = -1;
  double firstMinPos = -1;
  
  bool bMax = false;
  bool bMin = false;
  
  if(maxInd > this->IntAC.size()) return std::make_tuple(firstMax, firstMaxPos, firstMin, firstMinPos);
  
  
  for(unsigned int k = 2; k < maxInd; k++){
    if(this->IntAC[k-2] < this->IntAC[k-1] && this->IntAC[k-1] > this->IntAC[k] && this->IntAC[k-1] > this->IntAC[0] * 0.01 && bMax == false){
      firstMax = this->IntAC[k-1]/this->IntAC[0];
      firstMaxPos = this->dt * k;
      bMax = true;
    }
    if(this->IntAC[k-2] > this->IntAC[k-1] && this->IntAC[k-1] < this->IntAC[k] && bMin == false){
      firstMin = this->IntAC[k-1]/this->IntAC[0];
      firstMinPos = this->dt * k;
      bMin = true;
    }
    if(bMax == true && bMin == true) break;
  }
  
  return std::make_tuple(firstMax, firstMaxPos, firstMin, firstMinPos);
}
    
  
int evalMLTS::IntAC_MaxAboveTh(std::vector<double> &IACMaxPos, double max_thresh){
  IACMaxPos.resize(0);
  unsigned int maxInd = (unsigned int)(1.1*this->T0/this->dt);

  unsigned int max_counter = 0;
  if(maxInd > this->IntAC.size()) return max_counter;

  const unsigned int max_max_number = 100;
  bool bmax = true; //AC zero lag peak
  
  double tmp_maxp = 0.0;
  double tmp_max = 0.0;
  
    for(unsigned int k = 0; k < maxInd; k++){
      if(this->IntAC[k] > this->IntAC[0] * max_thresh && max_counter < max_max_number && bmax == false){
        bmax = true;
        tmp_max = 0.0;
      }
      if(bmax == true && this->IntAC[k] > this->IntAC[0] * max_thresh){
        if(this->IntAC[k] > tmp_max){
          tmp_max = this->IntAC[k];
          tmp_maxp = this->dt * k;
        }
      }
      if(bmax == true && this->IntAC[k] < this->IntAC[0] * max_thresh){
        bmax = false;
        IACMaxPos.push_back(tmp_maxp);
        max_counter += 1;
      }
    }
  
  return max_counter;
}

void evalMLTS::dump_IACmaxpos(double thresh){
  std::vector<double> IACmp;
  this->IntAC_MaxAboveTh(IACmp, thresh);
  std::cout << "IACmp size: " << IACmp.size() << std::endl;
  for(unsigned int k = 0; k < IACmp.size(); k++){
    std::cout << k << " " << IACmp[k] << std::endl;
  }
}


///////////////////////////////////////////////////////////////////
// misc
///////////////////////////////////////////////////////////////////

double evalMLTS::Variance(std::vector<double> &vec){
  if(vec.size() < 2) return -1;
  else{
    double mean = 0.0;
    double var = 0.0;
    for(std::size_t k = 0; k < vec.size(); k++) mean += vec[k]/vec.size();
    for(std::size_t k = 0; k < vec.size(); k++) var += (mean - vec[k]) * (mean - vec[k]) / vec.size();
    return var;
  }
}

void evalMLTS::smoothTS_cos(std::vector<double> &in, std::vector<double> &out, double width){
  int N = (int)(width / this->dt);
  N = N-(N%2)+1;
  out.resize(in.size()-N);
  std::vector<double> weight(N);
  double norm = 0.0;
  for(int k = 0; k < N; k++){
    weight[k] = (0.5 - 0.5*cos(2.0 * M_PI * k / (N-1.)));
    norm += weight[k];
  }
  for(std::size_t k = 0; k < out.size(); k++){
    out[k] = 0.0;
    for(int l = 0; l<N; l++) out[k] += weight[l] * in[k+l];
    out[k] /= norm;
  }
}

void evalMLTS::smoothTS(std::vector<double> &vec, std::vector<double> &out){
  out.resize(vec.size()-7);
  for(std::size_t k = 0; k < out.size(); k++) out[k] = (0.2*vec[k] + 0.5*vec[k+1] + 0.8*vec[k+2] + vec[k+3] + 0.8*vec[k+4] + 0.5*vec[k+5] + 0.2*vec[k+6])/4.0;
}

double evalMLTS::opticalFreq(std::vector<std::complex<double>> &ComplexE, std::vector<double> &Phase, std::vector<double> &Frequency, double dt){
  Phase.resize(ComplexE.size());
  Frequency.resize(ComplexE.size());
  int npi = 0;
  
  for(std::size_t k = 0; k < ComplexE.size() - 1; k++){
    Phase[k] = arg(ComplexE[k]) + npi * sm::PI;
    if (Phase[k] > sm::PI + arg(ComplexE[k+1]) + npi * sm::PI){
      npi += 2;
    }
    else if (Phase[k] < -sm::PI + arg(ComplexE[k+1]) + npi * sm::PI){
      npi -= 2;
    }
  }
  double OptFreqAv = 0;
  for(std::size_t k = 0; k < Phase.size()-2;k++){
    Frequency[k] = (Phase[k+1] - Phase[k])/dt;
    OptFreqAv += Frequency[k];
  }
  return OptFreqAv / (ComplexE.size()-2);
}

std::tuple<double,double> evalMLTS::interpolateQuadExtrm(const double X1, const double f1, const double X2, const double f2, const double X3, const double f3){
  //interpolates extremum where f2 is a local extremum at X2 in the DataVec_Ptr with neighboring points f1 and f3
  const double x1 = 0.0;
  const double x2 = X2 - X1;
  const double x3 = X3 - X1;
  
//   const double x1 = X1;
//   const double x2 = X2;
//   const double x3 = X3 - X1;
  
  const double k1 = f1/((x1-x2)*(x1-x3));
  const double k2 = f2/((x2-x1)*(x2-x3));
  const double k3 = f3/((x3-x1)*(x3-x2));
  
  const double a = k1 + k2 + k3;
  const double b = -(k1*(x2+x3) + k2*(x1+x3) + k3*(x1+x2));
  const double c = k1*x2*x3 + k2*x1*x3 + k3*x1*x2;
  
  const double xmax = -0.5*b/a+X1;
  const double fmax = 0.25*b*b/a - 0.5*b*b/a + c;
  
  return std::make_tuple(fmax,xmax);
}


///////////////////////////////////////////////////////////////////
// legacy stuff
///////////////////////////////////////////////////////////////////


double evalMLTS::Weighted_PulsePositions(std::vector<double> &PP){
  return evalMLTS::Weighted_PulsePositions(PP, this->TVPD_rel_pulse_trig, this->TVPD_rel_pulse_bounds);
}

double evalMLTS::Weighted_PulsePositions(std::vector<double> &PP, double trig_rel_th, double pulse_rel_bounds){
  //only use for FML or HML
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    PP.reserve(1000);
    PP.resize(0);
    
    double Inorm = 0.0;
    bool pulse = false;
    int lowerPbound = 0;
    int upperPbound = 0;
    double t_sum = 0.0;
    
    for(std::size_t k = (int)(WPP_vector_tBounds/this->dt); k < (*DataVec_Ptr).size()-(int)(WPP_vector_tBounds/this->dt); k++){
      if(pulse == false && (*DataVec_Ptr)[k] > this->TSGreatestMax * trig_rel_th){
        pulse = true;
        int wk = k-1;
        while((*DataVec_Ptr)[wk] > this->TSGreatestMax * pulse_rel_bounds && (k-wk)*this->dt < WPP_vector_tBounds/2.0){
          Inorm += (*DataVec_Ptr)[wk];
          wk--;
        }
        lowerPbound = wk + 1;
      }
      if(pulse == true && (*DataVec_Ptr)[k] > this->TSGreatestMax * pulse_rel_bounds ){
        Inorm += (*DataVec_Ptr)[k];
        upperPbound = k;
      }
      if(pulse == true && (*DataVec_Ptr)[k] < this->TSGreatestMax * pulse_rel_bounds ){
        pulse = false;
        for(int pk = lowerPbound; pk <= upperPbound; pk++){
          t_sum += ((*DataVec_Ptr)[pk] / Inorm) * pk * this->dt;
        }
        PP.push_back(t_sum);
        Inorm = 0.0;
        t_sum = 0.0;
      }
    }
    //return TSAverage maxima distance -> to be used to calc a clock time
    return (PP[PP.size()-1] - PP[0]) / (PP.size()-1);
  }
}

double evalMLTS::PtPJ_WPP(){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    std::vector<double> PP;
    this->Weighted_PulsePositions(PP);
    return this->PtPJitter(PP);
  }
}

double evalMLTS::PtPJ_WPP(double trig_rel_th){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    std::vector<double> PP;
    this->Weighted_PulsePositions(PP, trig_rel_th, this->TVPD_rel_pulse_bounds);
    return this->PtPJitter(PP);
  }
}

double evalMLTS::PtPJ_WPP(double trig_rel_th, double pulse_rel_bounds){
  if(TSGreatestMax/TSSmallestMin < (1.0 + max_min_dscr) || TSAverage < av_thresh) return -1;
  else{
    std::vector<double> PP;
    this->Weighted_PulsePositions(PP, trig_rel_th, pulse_rel_bounds);
    return this->PtPJitter(PP);
  }
}



