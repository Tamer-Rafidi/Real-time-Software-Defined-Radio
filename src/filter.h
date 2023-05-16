/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

#define pLL_freq 19e3

// add headers as needed
#include <iostream>
#include <vector>
#include <cmath>
//#include <numeric>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, int, std::vector<float> &);
void impulseResponseBPF(float, float , float, unsigned short int , std::vector<float>&);
void resample(const std::vector<float> &, const std::vector<float> &, std::vector<float> &, std::vector<float> &, int, int);
void split_iq_data(const std::vector<float> &, std::vector<float> &, std::vector<float> &);
void fmDemod(std::vector<float>&, const std::vector<float>&, const std::vector<float>&, float&,float&);

typedef struct pLL_state {
  float integrator;
  float phaseEst;
  float feedbackI;
  float feedbackQ;
  float	ncoOut1;
  float trigOffset;
} pLL_state;

void fmPll(std::vector<float> &, std::vector<float> &, float, float, pLL_state &, float, float, float);


#endif // DY4_FILTER_H
