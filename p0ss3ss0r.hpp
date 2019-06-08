/* p0ss3ss0r v0.1 by st3pan0va 2019 */
/* thanks to jeremysalwen for the lv2 template */

#ifndef __AGAIN_H
#define __AGAIN_H

#include <lv2plugin.hpp>
#include <fftw3.h>
#include "float.h"
#include "p0ss3ss0r.peg"

using namespace LV2;

class APossessor:public Plugin<APossessor> {
	public:
	APossessor(double srate);
	~APossessor();
	void run(uint32_t numsamples);
	void suspend();
	
protected:
	char programName[32];
	
private:
	double sampleRate;
	int lastSelected;
	int fftSize;
	static long rollOver;
	
	void processAudio(long numSampsToProcess, long fftFrameSize, long overlapFactor, double sampleRate, float *inDataL, float *inDataR, float *outData, double phaseMix, int loCut, double resonanceFactor);

	void makeLookup(int fftFrameSize);
	
	float* CarrierFIFO;
	float* ModulatorFIFO;
	float* OutBuffer;
	float* OutputAccum;
	float* FFTRealBuffer;
	float* PhaseL;
	float* MagnL;
	float* MagnR;
	double* Window;
	
	fftwf_complex * FFTworkspaceL;
	fftwf_complex * FFTworkspaceR;
	fftwf_plan forward_sp1;
	fftwf_plan forward_sp2;
	fftwf_plan backwards;

	fftwf_complex * FFTworkspaceL_alt;
	fftwf_complex * FFTworkspaceR_alt;
	fftwf_plan forward_sp1_alt;
	fftwf_plan forward_sp2_alt;
	fftwf_plan backwards_alt;

};

#endif
