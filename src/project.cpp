/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"

int main(int argc, char* argv[])
{

	int mode = 0;
	int channel = 1;
	if (argc < 2) {
        std::cerr << "Operating in default mode 0" << std::endl;
    } else if (argc == 2) {
        mode = atoi(argv[1]);
        if (mode > 3) {
            std::cerr << "Wrong mode" << mode << std::endl;
            exit(1);
        }
    } else {
        std::cerr << "Usage: " << argv[0] << std::endl;
        std::cerr << "or " << std::endl;
        std::cerr << "Usage: " << argv[0] << " <mode> " << std::endl;
        std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;

        exit(1);
    }
		std::cerr << "Operating in mode " << mode << std::endl;



	float RFFs = 0;
	float IFFs = 0;
	int up1 = 0;
	int down1 = 0;
	int up2 = 0;
	int down2 = 0;

	if(mode == 0){
		RFFs = 2400000;
		IFFs = 240000;
		up1 = 1;
		down1 = 10;
		up2 = 1;
		down2 = 5;

	}
	else if(mode == 1){
		RFFs = 2880000;
		IFFs = 288000;
		up1 = 1;
		down1 = 10;
		up2 = 1;
		down2 = 6;
	}
	else if(mode == 2){
		RFFs = 2400000;
		IFFs = 240000;
		up1 = 1;
		down1 = 10;
		up2 = 147;
		down2 = 800;
	}
	else if(mode == 3){
		RFFs = 2304000;
		IFFs = 256000;
		up1 = 1;
		down1 = 9;
		up2 = 441;
		down2 = 2560;
	}


	float monoFc = 16000;
	float pilotFB = 18500;
	float pilotFE = 19500;
	float stereoFc1 = 23000;
	float stereoFc2 = 38000;
	float stereoFc3 = 53000;

	float rfFc = 100e3;

	unsigned int num_taps = 71;
  unsigned int num_taps_up = (num_taps-1)*up2+1;
	int block_size = int((1024*down1*down2/up2)*2);
	int i_data_size = int(block_size/2);

	pLL_state pLLstate;
  pLLstate.integrator = 0.0;
  pLLstate.phaseEst = 0.0;
  pLLstate.feedbackI = 1.0;
  pLLstate.feedbackQ = 0.0;
  pLLstate.ncoOut1 = 1.0;
  pLLstate.trigOffset = 0.0;

	std::vector<float> i_filt;
	std::vector<float> q_filt;

	std::vector<float> stereo_left_channel;
  std::vector<float> stereo_right_channel;

	std::vector<float> audio_coeff_mono;
	std::vector<float> audio_coeff_stereo_left;
	std::vector<float> audio_coeff_stereo_right;

	float pre_I;
	float pre_Q;

	std::vector<float> i_state(num_taps-1,0.0);
	// i_state.clear();
	std::vector<float> q_state(num_taps-1,0.0);
	std::vector<float> state(num_taps_up-1,0.0);

	std::vector<float> stereo_state(num_taps-1,0.0);
	std::vector<float> pilot_state(num_taps-1,0.0);
	std::vector<float> stereo_filt_state(num_taps_up-1);

	std::vector<float> rf_coeff;

	std::vector<float> processed_data;

	std::vector<float> stereo_carrier;
	std::vector<float> stereo_data;
	std::vector<float> pilot;
	std::vector<float> pilot_coeff;
	std::vector<float> stereo_coeff;
	std::vector<float> audio_coeff_stereo;

	std::vector<float> delay_state(up2*(num_taps-1)/(2*down2),0.0);
	std::vector<float> delay_data;

	impulseResponseLPF(RFFs, rfFc, num_taps, 1, rf_coeff);
	impulseResponseLPF(IFFs*up2, monoFc, num_taps_up, 1, audio_coeff_mono);

	impulseResponseBPF(IFFs, pilotFB, pilotFE, num_taps, pilot_coeff);
	impulseResponseBPF(IFFs, stereoFc1, stereoFc3, num_taps, audio_coeff_stereo);
	impulseResponseLPF(IFFs*up2 , monoFc, num_taps_up, 1, stereo_coeff);
//Block processing begins
	for(unsigned int block_id = 0;; block_id++){
		std::vector<float> block_data(block_size); //create array of block_size stored in block_data
		readStdinBlockData(block_size, block_id, block_data);
		//------------------------------------Mono path begins---------------------------------------------------------------
		std::vector<float> i_data(i_data_size);
		std::vector<float> q_data(i_data_size);
		split_iq_data(block_data,i_data,q_data);
		//---------------------------------------------------------------------------------------------------
		resample(rf_coeff, i_data, i_filt, i_state, down1, up1);
		resample(rf_coeff, q_data, q_filt, q_state, down1, up1);
		//---------------------------------------------------------------------------------------------------
		std::vector<float> fm_demod(i_filt.size());
		fmDemod(fm_demod, i_filt, q_filt, pre_I, pre_Q);
		//---------------------------------------------------------------------------------------------------
		resample(audio_coeff_mono, fm_demod, processed_data, state, down2, up2);
		//------------------------------Mono finished---------------------------------------------------------------------
		//------------------------------stereo begins--------------------------------------------------------------------
		//------------------------------pilot extraction and PLL---------------------------------------------------------------------
	  resample(pilot_coeff, fm_demod, pilot, pilot_state, 1, 1);
	  fmPll(pilot, stereo_carrier, pLL_freq, IFFs, pLLstate, 2.0, 0.0, 0.01);
//------------------------stereo channel extraction---------------------------------------------------------------------------
		resample(audio_coeff_stereo, fm_demod, stereo_data, stereo_state, 1, 1);
//------------------------stereo shifting---------------------------------------------------------------------------
		std::vector<float> mixed_data(stereo_data.size(), 0);
		for (int i = 0; i < mixed_data.size(); i++) {
		    mixed_data[i] = 2*stereo_data[i] * stereo_carrier[i];

		}

//------------------------------------------------------------------------------------------------
		std::vector<float> stereo_filt;
		resample(stereo_coeff, mixed_data, stereo_filt, stereo_filt_state, down2, up2);

	  stereo_left_channel.resize(stereo_filt.size(), 0);
		stereo_right_channel.resize(stereo_filt.size(), 0);
//------------------------------mono delay------------------------------------------------------------------
		delay_data.resize(processed_data.size(),0.0);
		for(int i =0; i < (int)delay_state.size();i++)
		{
			delay_data[i] = delay_state[i];
			delay_state[i] = processed_data[processed_data.size() - delay_state.size() + i];
		}
		for(int j =delay_state.size(); j < (int)delay_data.size();j++)
		{
			delay_data[j] = processed_data[j-delay_state.size()];
		}
//-------------------------------stereo mixer-----------------------------------------------------------------
		std::vector<float> both;
		for (int i=0;i<stereo_right_channel.size();i++){
			both.push_back((delay_data[i] + stereo_filt[i])/2);
			both.push_back((delay_data[i] - stereo_filt[i])/2);
		}
//----------------------------------Processing finished, write data--------------------------------------------------------------
		if (channel == 0){
			std::vector<short int> audio_data(processed_data.size());
			for(unsigned int k=0; k<(unsigned int)(processed_data.size());k++){

		    if(std::isnan(processed_data[k])){
					audio_data[k] = 0;

				}
		    else{
					audio_data[k] = static_cast<short int>(processed_data[k]*16384);
					//std::cout << "1";
				}
			}
			fwrite(&audio_data[0],sizeof(short int), audio_data.size(),stdout);
		}

		else if (channel == 1){
			std::vector<short int> audio_data(both.size());
			for(unsigned int k=0; k<(unsigned int)(both.size());k++){

		    if(std::isnan(both[k])){
					audio_data[k] = 0;

				}
		    else{
					audio_data[k] = static_cast<short int>(both[k]*16384);
					//std::cout << "1";
				}
			}
			fwrite(&audio_data[0],sizeof(short int), audio_data.size(),stdout);
		}


	}
	return 0;
}
