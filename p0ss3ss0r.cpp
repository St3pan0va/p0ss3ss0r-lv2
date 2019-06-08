/* p0ss3ss0r v0.1 by st3pan0va 2019 */

#include "p0ss3ss0r.hpp"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PI 3.141592653589793
#define FFTWINDOW 2048 // max FFT frame length

//-----------------------------------------------------------------------------
APossessor::APossessor(double rate) : Plugin<APossessor>(p_n_ports)
{
	fftSize=FFTWINDOW; //initialise FFT size
	lastSelected=1; //used to check if FFT size has changed (need to clear buffer if so)
	sampleRate=rate;

	CarrierFIFO = new float [FFTWINDOW];
	ModulatorFIFO = new float [FFTWINDOW];
	OutBuffer = new float [FFTWINDOW];
	FFTRealBuffer=new float[FFTWINDOW];
	OutputAccum = new float [2*FFTWINDOW];
	PhaseL = new float [FFTWINDOW];
	MagnL = new float [FFTWINDOW];
	MagnR = new float [FFTWINDOW];

	Window = new double [FFTWINDOW];
	makeLookup(FFTWINDOW); // make lookup table for window function

	/* initialise plans and workspaces for FFTW - two different window sizes */
	/* FFTW needs an appropriate plan and a workspace for each FFT size */

	FFTworkspaceR = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * FFTWINDOW);
	FFTworkspaceL = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * FFTWINDOW);

	forward_sp1= fftwf_plan_dft_r2c_1d(FFTWINDOW, FFTRealBuffer, FFTworkspaceL,
                                    FFTW_ESTIMATE);
	forward_sp2= fftwf_plan_dft_r2c_1d(FFTWINDOW, FFTRealBuffer, FFTworkspaceR,
                                    FFTW_ESTIMATE);	
	backwards=fftwf_plan_dft_c2r_1d(FFTWINDOW, FFTworkspaceL, FFTRealBuffer,
                                    FFTW_ESTIMATE);

	FFTworkspaceL_alt = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * FFTWINDOW/4);
	FFTworkspaceR_alt = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * FFTWINDOW/4);

	forward_sp1_alt = fftwf_plan_dft_r2c_1d(FFTWINDOW/4, FFTRealBuffer, FFTworkspaceL_alt,
                                    FFTW_ESTIMATE);
	forward_sp2_alt = fftwf_plan_dft_r2c_1d(FFTWINDOW/4, FFTRealBuffer, FFTworkspaceR_alt,
                                    FFTW_ESTIMATE);
	backwards_alt = fftwf_plan_dft_c2r_1d(FFTWINDOW/4, FFTworkspaceL_alt, FFTRealBuffer,
                                    FFTW_ESTIMATE);

	suspend ();		// flush buffer

}

//-----------------------------------------------------------------------------------------
APossessor::~APossessor() // delete buffers in destructor
{
	delete[] CarrierFIFO;
	delete[] ModulatorFIFO;
	delete[] FFTRealBuffer;
	delete[] OutBuffer;
	delete[] OutputAccum;
	delete[] PhaseL;
	delete[] MagnL;
	delete[] MagnR;	
	delete[] Window;

	fftwf_free(FFTworkspaceL);
	fftwf_free(FFTworkspaceR);
	fftwf_free(FFTworkspaceL_alt);
	fftwf_free(FFTworkspaceR_alt);
}

//-----------------------------------------------------------------------------------------
void APossessor::suspend ()
{
	memset(CarrierFIFO, 0, FFTWINDOW*sizeof(float));
	memset(ModulatorFIFO, 0, FFTWINDOW*sizeof(float));
	memset(OutBuffer, 0, FFTWINDOW*sizeof(float));
	memset(FFTworkspaceL, 0, 2*FFTWINDOW*sizeof(float));
	memset(FFTworkspaceR, 0, 2*FFTWINDOW*sizeof(float)); 
	memset(FFTworkspaceL_alt, 0, FFTWINDOW/2*sizeof(float));
	memset(FFTworkspaceR_alt, 0, FFTWINDOW/2*sizeof(float));     
	memset(OutputAccum, 0, 2*FFTWINDOW*sizeof(float));
	memset(PhaseL, 0, FFTWINDOW*sizeof(float));
	memset(MagnR, 0, FFTWINDOW*sizeof(float));
	memset(MagnL, 0, FFTWINDOW*sizeof(float));
}

//-----------------------------------------------------------------------------------------
void APossessor::run(uint32_t sampleFrames)
{
	// get parameters from control ports 

	int lo = int (*p(p_lowcut)); // low cut threshold 0-128
	double resonance = *p(p_resonance); // should be 0.0-1.0
	resonance = resonance < 0 ? 0 : resonance; // range check
	double mix = *p(p_phasemix); // should be 0.0-1.0
	mix = mix < 0 ? 0 : mix; // range check
	double select = *p(p_fftselector);
	select = select < 0 ? 0 : select;
	int iselect = int (*p(p_fftselector) + 1.49); // fftselector input range is 0.0-1.0 but we want 1 or 2
	iselect *= iselect; // now int select == either 1 or 4 

	if (lastSelected != iselect) // clear output buffer on FFT size switch
	{
		memset(OutputAccum, 0, 2*FFTWINDOW*sizeof(float));
		memset(OutBuffer, 0, FFTWINDOW*sizeof(float));
	}
	lastSelected=iselect;
	
	lo = lo / iselect; // scale low cut threshold to FFT window size
	int overlap = 8; // must be a multiple of 4. 

	// call process

	processAudio(sampleFrames, iselect, overlap, sampleRate, p(p_left), p(p_right), p(p_out), mix, lo, resonance);
}

// -----------------------------------------------------------------------------------------------------------------


void APossessor::processAudio(long numSampsToProcess, long fftSelect, long overlapFactor, double sampleRate, float *inDataL, float *inDataR, float *outData, double phaseMix, int loCut, double resonanceFactor)

{
	static long rollOver=false;

	/* set FFT frame size from selector and max FFTWINDOW */
	
	long fftFrameSize = FFTWINDOW / fftSelect; // fftSelect should be int 1 or 4 from user input fftselector

	/* set up variables for overlap-add thank you Mr. Bernsee */

	long fftFrameSize2 = fftFrameSize/2;
	long stepSize = fftFrameSize/overlapFactor;
	long inFifoLatency = fftFrameSize-stepSize;
	double outputFactor = 1.0 / ((double)fftFrameSize2*(double)overlapFactor);

	if (rollOver == false) 
	{
		rollOver = inFifoLatency;
	}

	/* main processing loop */
	for (long i = 0; i < numSampsToProcess; i++)
	{
		/* As long as we have not yet collected enough data just read in */
		CarrierFIFO[rollOver] = inDataL[i]; //
		ModulatorFIFO[rollOver] = inDataR[i]; //  input gain factors could go here
		outData[i] = OutBuffer[rollOver-inFifoLatency]; // output, gain control could go here
		
		rollOver++;

		/* if we have enough data for processing */
		if (rollOver >= fftFrameSize) 
	       {
			rollOver = inFifoLatency;

			/* choose which FFT option to execute (or fall through if ffSelect is not 1 or 4) */

			if (fftSelect == 1) 
			{
				/* use longer FFT window (smoother result) */

				/* do windowing on left channel (carrier) input */
				for (long k = 0; k < fftFrameSize; k++) 
				{
					FFTRealBuffer[k]=CarrierFIFO[k] * Window[k * fftSelect]; 
				}

				/* FFT carrier input */
				fftwf_execute(forward_sp1);
			
				/* do windowing on right channel (modulator) input */
				for (long k = 0; k < fftFrameSize; k++) 
				{
					FFTRealBuffer[k]=ModulatorFIFO[k] * Window[k * fftSelect];
				}

				/* FFT modulator input */
				fftwf_execute(forward_sp2);

				/* processing in frequency domain from here */

				for (long k = loCut; k <= fftFrameSize2; k++) 
				{
					/* de-interlace FFT buffers */
					double realL = FFTworkspaceL[k][0];
					double imagL = FFTworkspaceL[k][1];
 					double realR = FFTworkspaceR[k][0];
					double imagR = FFTworkspaceR[k][1];	

					/* compute both magnitudes but only left (carrier) signal phases */
					MagnL[k] = 2.*sqrt(realL*realL + imagL*imagL);
					MagnR[k] = 2.*sqrt(realR*realR + imagR*imagR);
	
					PhaseL[k] = atan2(imagL, realL);
				}

				// lo cut
				for (long k = 0; k <= loCut; k++) 
				{  
					FFTworkspaceL[k][0]=0;
					FFTworkspaceL[k][1]=0;
				}

				for (long k = loCut; k <= fftFrameSize2; k++) 
				{
					/* something like convolution but not quite */
					// morph left and right magnitudes - with resonance
					double newMagn1 = resonanceFactor * (MagnL[k]*MagnR[k]) + (1.0-resonanceFactor) * sqrt(MagnL[k]*MagnR[k]); 
					// 'vocoder' uses right channel (modulator) only - with resonance
					double newMagn2 = MagnR[k] * (1.0-resonanceFactor) + MagnR[k]*MagnR[k] * resonanceFactor; 
					// mix between 'vocoder' and morph according to mix setting
					double newMagn = (newMagn1*phaseMix + newMagn2*(1.0-phaseMix)); 					
							
					/* back to complex coordinates (real, imag) from magnitude and phase */
					FFTworkspaceL[k][0] = newMagn*cos(PhaseL[k]); // use bin phase from left (carrier) input only
					FFTworkspaceL[k][1] = newMagn*sin(PhaseL[k]); //
				} 

				/* do inverse transform */
				fftwf_execute(backwards);

				/* do windowing and add to output accumulator */ 
				for(long k=0; k < fftFrameSize; k++) 
				{
					OutputAccum[k] += Window[k * fftSelect] * FFTRealBuffer[k] * outputFactor;
				}

			}

			/* if small fft window selected - execute the alt version of fftw instead */
			/* guess there's a "better" way to do this */
			/* but this was easy and costs nothing */

			if (fftSelect == 4) 
			{
				/* use shorter FFT window four times */

				/* do windowing on left channel (carrier) input */
				for (long k = 0; k < fftFrameSize; k++) 
				{
					FFTRealBuffer[k]=CarrierFIFO[k] * Window[k * fftSelect]; // fftSelect is also used to downsample the window function
				}
				/* transform left channel input alt version */
				fftwf_execute(forward_sp1_alt);
			
				/* do windowing on right channel (modulator) input */
				for (long k = 0; k < fftFrameSize; k++) 
				{
					FFTRealBuffer[k]=ModulatorFIFO[k] * Window[k * fftSelect];
				}
				/* transform right channel (modulator) input alt version */
				fftwf_execute(forward_sp2_alt);

				/* processing in frequency domain from here */

				for (long k = loCut; k <= fftFrameSize2; k++) 
				{
					/* de-interlace FFT buffers */
					double realL = FFTworkspaceL_alt[k][0];
					double imagL = FFTworkspaceL_alt[k][1];
	 				double realR = FFTworkspaceR_alt[k][0];
					double imagR = FFTworkspaceR_alt[k][1];	

					/* compute both magnitudes but only left (carrier) signal phases */
					MagnL[k] = 2.*sqrt(realL*realL + imagL*imagL);
					MagnR[k] = 2.*sqrt(realR*realR + imagR*imagR);

					PhaseL[k] = atan2(imagL, realL);
				}

				// lo cut just zeros fft bins - sinc filter considered harmful, we don't care here ;-)
				for (long k = 0; k <= loCut; k++) 
				{  
					FFTworkspaceL_alt[k][0]=0;
					FFTworkspaceL_alt[k][1]=0;
				}

				for (long k = loCut; k <= fftFrameSize2; k++) 
				{
					/* something like convolution but not quite */
					// morph left and right - with resonance
					double newMagn1 = resonanceFactor * (MagnL[k]*MagnR[k]) + (1.0-resonanceFactor) * sqrt(MagnL[k]*MagnR[k]); 
					// 'vocoder' using right channel (modulator) only - with resonance
					double newMagn2 = MagnR[k] * (1.0-resonanceFactor) + MagnR[k]*MagnR[k] * resonanceFactor; 
					// mix between 'vocoder' and morph according to mix setting
					double newMagn = (newMagn1*phaseMix + newMagn2*(1.0-phaseMix)); 					
							
					/* back to complex coordinates (real, imag) from magnitude and phase alt version */
					FFTworkspaceL_alt[k][0] = newMagn*cos(PhaseL[k]); // use bin phase from left (carrier) input only
					FFTworkspaceL_alt[k][1] = newMagn*sin(PhaseL[k]); //
				} 

				/* do inverse transform alt version */
				fftwf_execute(backwards_alt);

				/* do windowing and add to output accumulator */ 
				for(long k=0; k < fftFrameSize; k++) 
				{
					OutputAccum[k] += Window[k * fftSelect] * FFTRealBuffer[k] * outputFactor;
				}

			}

			/* transfer output accum to output buffer */
			memmove (OutBuffer, OutputAccum, stepSize*sizeof(float));

			/* shift accumulator */
			memmove (OutputAccum, OutputAccum+stepSize, fftFrameSize*sizeof(float));

			memmove (CarrierFIFO, CarrierFIFO+stepSize, inFifoLatency*sizeof(float));
			memmove (ModulatorFIFO, ModulatorFIFO+stepSize, inFifoLatency*sizeof(float));
					
		}
	}
}

// -----------------------------------------------------------------------------------------------------------------

/* make lookup table for window function */
/* the window implementation follows Stephan Bernsee's example */

void APossessor::makeLookup(int fftMax)
{
	for (int k=0; k<fftMax; k++) 
	{
		Window[k] = -.5*cos(2*PI*(double)k/(double)fftMax)+.5;
	}
}

/* register plugin - */

static int _ = APossessor::register_class(p_uri);
