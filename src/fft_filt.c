/*
 * fft_filt.c
 *
 *  Created on: Nov 24, 2016
 *      Author: Walter
 */

/*
 * FFTlength N = (L + M - 1)
 * first M - 1 samples are from previous block, or zero for first block
 * discard the first M - 1 samples from FFT output
 * example
 * 		N = 1024 = 8*128
 * 		M = 129
 * 		L = 7*128
 *
 * 		y = {y[L-M+1].... y[L-1], x[0],   , x[L-1]}
 *
 * 		process y -> z
 *
 * 		discard z[0].... z[M-2], keep L samples z[M-1] .... Z[N]
 */
#include "arm_math.h"
#include "arm_const_structs.h"

#define NN 1024
#define NF 128
#define ASTOP 90.0f

static float bb[2*NN];
static float yy[2*NN];
static float zz[2*NN];

//#define PI 3.14159265358979f // already defined in arm_math.h
#define TPI	 6.28318530717959f
#define PIH	 1.57079632679490f

float m_sinc(int m, float fc)
{	// fc is f_cut/(Fsamp/2)
	// m is between -M and M step 2
	//
	float x = m*PIH;
	if(m == 0)
		return 1.0f;
	else
		return sinf(x*fc)/(fc*x);
}

float32_t Izero (float32_t x)
{
    float32_t x2 = x / 2.0;
    float32_t summe = 1.0;
    float32_t ds = 1.0;
    float32_t di = 1.0;
    float32_t errorlimit = 1e-9;
    float32_t tmp;
    do
    {   tmp = x2 / di;
        tmp *= tmp;
        ds *= tmp;
        summe += ds;
        di += 1.0;
    }   while (ds >= errorlimit * summe);
    return (summe);
}  // END Izero

void calc_FIR_coeffs (float * coeffs, int numCoeffs, float32_t fc, float32_t Astop, int type, float dfc)
 {	// modified after
	// Wheatley, M. (2011): CuteSDR Technical Manual. www.metronix.com, pages 118 - 120, FIR with Kaiser-Bessel Window
	// assess required number of coefficients by
	//     numCoeffs = (Astop - 8.0) / (2.285 * TPI * normFtrans);
	// selecting high-pass, numCoeffs is adapter to become an even number as required for high-pass

	 int ii,jj;
     float32_t Beta;
     float32_t izb;
     float fcf;
     int nc;

     // calculate Kaiser-Bessel window shape factor beta from stop-band attenuation
     if (Astop < 20.96)
    	 Beta = 0.0;
     else if (Astop >= 50.0)
    	 Beta = 0.1102 * (Astop - 8.71);
     else
    	 Beta = 0.5842 * powf((Astop - 20.96), 0.4) + 0.7886 * (Astop - 20.96);

     izb = Izero (Beta);
     if(type == 0) // low pass filter
     {	fcf = fc;
     	nc =  numCoeffs;
     }
     else if(type == 1) // high-pass filter
     {	fcf = -fc;
     	nc =  2*(numCoeffs/2);
     }
     else if ((type == 2) || (type==3)) // band-pass filter
     {
    	 fcf = dfc;
     	 nc =  2*(numCoeffs/2); // maybe not needed
     }
     else if (type==4)	// Hilbert transform
	 {
    	 // clear coefficients
    	 for(ii=0; ii< 2*(nc-1); ii++) coeffs[ii]=0;
    	 // set real delay
    	 coeffs[nc]=1;

    	 // set imaginary Hilbert coefficients
    	 for(ii=1; ii< (nc+1); ii+=2)
    	 {
		    	 float x =(float)(2*ii - nc)/(float)nc;
		    	 float w = Izero(Beta*sqrtf(1.0f - x*x))/izb; // Kaiser window
				 coeffs[2*ii+1] = 1.0f/(PIH*(float)(ii-nc/2)) * w ;
    	 }
    	 return;
	 }

     for(ii= - nc, jj=0; ii< nc; ii+=2,jj++)
     {
    	 float x =(float)ii/(float)nc;
    	 float w = Izero(Beta*sqrtf(1.0f - x*x))/izb;	// Kaiser window
    	 coeffs[jj] = fcf * m_sinc(ii,fcf) * w;
     }

     if(type==1)
     {
    	 coeffs[nc/2] += 1;
     }
     else if (type==2)
     {
       	 for(jj=0; jj< nc+1; jj++) coeffs[jj] *= 2.0f*cosf(PIH*(2*jj-nc)*fc);
     }
     else if (type==3)
     {
       	 for(jj=0; jj< nc+1; jj++) coeffs[jj] *= -2.0f*cosf(PIH*(2*jj-nc)*fc);
    	 coeffs[nc/2] += 1;
     }

 } // END calc_lowpass_coeffs



void fft_filt_init(float fc, int type, float dfc)
{   int ii;
	for(ii=0; ii< 1024; ii++)
	{
		bb[2*ii]=0; bb[2*ii+1]=0;
		yy[2*ii]=0; yy[2*ii+1]=0;
	}
	float coeffs[NF+1];
	calc_FIR_coeffs(coeffs, NF, fc, ASTOP, type, dfc);

	for(ii=0; ii<(NF+1); ii++) bb[2*ii]=coeffs[ii];
	//
	arm_cfft_f32(&arm_cfft_sR_f32_len1024, bb, 0, 1);
}

 void fft_filt_exec(float *zr, float *zi, float *xr, float *xi, int nx, int MM)
{	int ii;
	int LL = NN - MM + 1;

	if(LL != nx) return;

	// copy data to FFT buffer
	float *yp = &yy[2*(MM-1)];
	for(ii = 0; ii < LL; ii++) {yp[2*ii] = xr[ii]; yp[2*ii+1] = xi[ii];}

	// perform complex FFT (in-place)
	arm_cfft_f32(&arm_cfft_sR_f32_len1024, yy, 0, 1);

	// complex multiply with filter spectrum
	arm_cmplx_mult_cmplx_f32 (yy, bb, zz, NN);

	// perform iFFT (in-place)
	arm_cfft_f32(&arm_cfft_sR_f32_len1024, zz, 1, 1);

	// copy overlap to beginning
	float *xrp = &xr[LL-MM+1];
	float *xip = &xi[LL-MM+1];
	for(ii = 0; ii < MM-1; ii++) {yy[2*ii] = xrp[ii]; yy[2*ii+1] = xip[ii];}

	// copy data to result buffer
	float *zp = &zz[2*(MM-1)];
	for(ii = 0; ii < LL; ii++) { zr[ii] = zp[2*ii]; zi[ii] = zp[2*ii+1];}
}

