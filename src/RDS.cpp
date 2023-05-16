/*

#include "dy4.h"
#include "RDS.h"

//We weren’t able to finish RDS but did what we can and commented it out so it doesn’t interfere with the implementation of mono and stereo

void RRC(float Fs,int num_taps,std::vector<float> &h ){
  float Ts = 1/2375.0;
  float beta = 0.90;
  h.resize(num_taps,0);
  for (int k=0; k < num_taps; k++){
    float t;
    t = float((k - num_taps/2))/Fs;
    if(t == 0.0) {  h[k] = 1.0 + beta *((4/PI)-1); }
    else if ((t == -Ts/(4*beta)) || (t == Ts/(4*beta))) {
      h[k] = (beta/std::sqrt(2))*(((1+2/PI)* (std::sin(PI/(4*beta)))) + ((1-2/PI)*(std::cos(PI/(4*beta)))));
    } 
    else {
      h[k] = (std::sin(PI*t*(1-beta)/Ts) + 4*beta*(t/Ts)*std::cos(PI*t*(1+beta)/Ts))/(PI*t*(1-(4*beta*t/Ts)*(4*beta*t/Ts))/Ts);
    }
  }
}

void RDS_processing( const std::vector<float> &processed_data, const std::vector<float> &fm_demod,
		                     std::vector<float> &RDS_state, std::vector<float> &RDS_filt_state, pLL_state &pLLstate, int IFFs, int RDSfc1, int RDSfc2,
		                     int RDSfc3, int RDSfc4, unsigned short num_taps, int MonoFc, unsigned short num_taps){
				     
				     
    // These parameters would go into our modes in project.cpp. however, we did not want to add any RDS code into that file as it is incomplete
    float RDSfc1 = 54000;
    float RDSfc2 = 60000;
    float RDSfc3 = 113500;
    float RDSfc4 = 114500;
    float RDSfc5 = 3000;
    int SPS0 = 11;
    int SPS2 = 46;
    int GCD_RDS0 = std::__gcd(IFFs, SPS0*2375);	
    int GCD_RDS2 = std::__gcd(IFFs, SPS2*2375);	
	
    std::vector<float> RDS_carrier;
    std::vector<float> RDS_data;
    std::vector<float> RDS_coeff;
    std::vector<float> audio_coeff_RDS;
    std::vector<float> recovered_coeff_RDS;
    
    	// BPF for  recovery
	  impulseResponseBPF(IFFs, RDSfc3, RDSfc4, num_taps, recovered_coeff_RDS);
    	// BPF for channel extraction
	  impulseResponseBPF(IFFs, RDSfc1, RDSfc2, num_taps, audio_coeff_RDS);
	  

   	// RDS Channel Extraction
	resample(audio_coeff_RDS, fm_demod, RDS_data, RDS_state, 1, 1);

    	// RDS Carrier Recovery
	RDS_resized_data.resize(RDS_data.size());
	for (size_t i = 0; i < RDS_data.size(); ++i) {
   		RDS_resized_data[i] = RDS_data[i] * RDS_data[i];
	}
	resample(recovered_coeff_RDS, fm_demod, RDS_data, RDS_state, 1, 1);
	fmPll(recovered_coeff_RDS, RDS_carrier, pLL_freq, IFFs, pLLstate, 0.5, 0.0, 0.01);
	
	// Mixing
	std::vector<float> RDS_mixed(RDS_data.size(),0);
	for (int i = 0; i < RDS_mixed.size(); i++){
		RDS_mixed[i] = RDS_data[i] * RDS_carrier[i];
		RDS_mixed[i] *= 2;
	}
	
	impulseResponseLPF(IFFs , RDSFc5, num_taps, 1, RDS_mixed);
	
	// RDS Demodulation (Not complete)
	resample()..
*/
