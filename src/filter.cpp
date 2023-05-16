/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, int amplify, std::vector<float> &h)
{
	// allocate memory for the impulse response
	float h_i;
	float Norm_Cutoff = Fc/(Fs/2);
	h.clear(); h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (int i=0; i< num_taps; i++){
		if (i == (num_taps-1)/2){ h_i = Norm_Cutoff; }
		else { h_i = Norm_Cutoff * (sin(PI*Norm_Cutoff*(i-(num_taps-1)/2)))/(PI*Norm_Cutoff*(i-(num_taps-1)/2)); }
		h_i = amplify * (h_i * pow(sin((i*PI)/num_taps),2));
		h[i] = h_i;
	}
}

void impulseResponseBPF(float Fs, float fb, float fe, unsigned short int num_taps, std::vector<float>& h) {
    // Allocate memory for the impulse response
    h.clear(); h.resize(num_taps, 0.0);

    // Compute filter coefficients.
    float Fc = (fe + fb) / 2;
    float BW = fe - fb;
    float norm_pass = BW / (Fs/2);
    float center_freq =  Fc / (Fs/2);
    for (int i = 0; i < (num_taps - 1); i++) {
        if (i == (num_taps - 1) / 2) {
            h[i] = norm_pass;
        } else {
            h[i] = norm_pass*(std::sin(PI*(norm_pass/2)*(i-(num_taps-1)/2))/(PI*(norm_pass/2)*(i-(num_taps-1)/2)));
        }
        // Apply a frequency shift by the center frequency
        h[i] = h[i] * std::cos(i* PI * center_freq);
        // Apply the Hann window
        h[i] = h[i] * std::pow(std::sin(i * PI / num_taps), 2);
    }
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"

void resample(const std::vector<float> &coeff, const std::vector<float> &data, std::vector<float> &output, std::vector<float> &state, int down_sample, int up_sample)
{

		output.clear(); output.resize(data.size()*up_sample/down_sample,0.0);



	for(int i = 0;  i < (int)output.size(); i++){
			int phase = (int)((i * down_sample) % up_sample);
			output[i] = 0;
			for(int j=phase; j<(int)coeff.size();j+=up_sample){
					int A = (int)((i*down_sample-j)/(float)up_sample);

					if(A >=0 ){
							output[i] += coeff[j] * data[A];
					}
					else{
							output[i]+= coeff[j]*state[state.size() + A];
					}
			}
			output[i] = up_sample* output[i];
	}

	for(int k =0; k<(int)state.size();k++){
			state[k]=data[data.size()-state.size()+k];
	}
}

void split_iq_data(const std::vector<float> &bin_data, std::vector<float> &i_data, std::vector<float> &q_data){
	for (unsigned int i = 0; i < bin_data.size(); i+=2) {
			i_data[i/2] = bin_data[i];
			q_data[i/2] = bin_data[i+1];

	}
}


void fmDemod(std::vector<float>& fm_demod, const std::vector<float>& I, const std::vector<float>& Q, float& pre_I, float& pre_Q)
{

    float dQ, dI;

    for (int k = 0; k < (int)(I.size()); k++) {
        if (k == 0) {
            dQ = Q[k]-pre_Q;
            dI = I[k]-pre_I;
        }
        else {
            dQ = Q[k] - Q[k-1];
            dI = I[k] - I[k-1];
        }
        if ((pow(I[k], 2) + pow(Q[k], 2)) == 0) {
            fm_demod[k] = 0;
        }
        else {
            fm_demod[k] = (1 / (pow(I[k], 2) + pow(Q[k], 2))) * ((I[k] * dQ) - (Q[k] * dI));
        }
			}
			pre_I = I[I.size()-1];
			pre_Q = Q[Q.size()-1];
}

	// Changed ncoscale to 2 because that is what the stereo mode is
void fmPll(std::vector<float> &pllIn, std::vector<float> &ncoOut, float freq, float Fs, pLL_state &state, float ncoscale, float phaseAdjust, float normBandwidth){
	  float Cp = 2.666;
  	float Ci = 3.555;

  	// gain for the porportional term
  	float Kp = normBandwidth * Cp;
  	// gain for the intergrator term
  	float Ki = (normBandwidth * normBandwidth) * Ci;

  	// output array for the NCO
  	ncoOut.resize(pllIn.size()+1,0);

 	  // Initialize internal state
	  float integrator = state.integrator;
    float phaseEst = state.phaseEst;
    float feedbackI = state.feedbackI;
    float feedbackQ = state.feedbackQ;
	  ncoOut[0] = state.ncoOut1;
    float trigOffset = state.trigOffset;

  	for (int k = 0; k < pllIn.size(); k++) {
    		// Phase detector
    		float errorI = pllIn[k] * (+feedbackI);
    		float errorQ = pllIn[k] * (-feedbackQ);


   		  // Phase error detection
   		  float errorD = std::atan2(errorQ, errorI);

    		// Loop filter
    		integrator = integrator + Ki * errorD;

    		// Update phaseEst
    		phaseEst = phaseEst + Kp * errorD + integrator;

    		// NCO (internal oscillator)
    		trigOffset += 1;
    		float trigArg = 2 * PI * (freq/Fs) * trigOffset + phaseEst;
    		feedbackI = std::cos(trigArg);
    		feedbackQ = std::sin(trigArg);
    		ncoOut[k+1] = std::cos(trigArg * ncoscale + phaseAdjust);
   	}
		state.integrator = integrator;
	  state.phaseEst = phaseEst;
	  state.feedbackI = feedbackI;
	  state.feedbackQ = feedbackQ;
  	state.ncoOut1 = ncoOut[ncoOut.size()-1];
	  state.trigOffset = trigOffset;
}
