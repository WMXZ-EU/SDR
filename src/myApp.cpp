/* myApp.cpp
 *
 *  Created on: Aug 28, 2016
 *      Author: Walter
 */
//
#include <math.h>

#include "kinetis.h"
#include "core_pins.h"
#include "usb_serial.h"

extern "C" void setup(void);
extern "C" void loop(void);

extern "C"
{
	void fft_filt_init(float fc, int type, float dfc);
	void fft_filt_exec(float *zr, float *zi, float *xr, float *xi, int nx, int MM);
}

///--------------------------- Software trigger -----------------------------------------
void doProcessing(void);

#define TRIGGER_SWI(irq, prio) {NVIC_CLEAR_PENDING(irq); NVIC_SET_PRIORITY(irq, prio); NVIC_TRIGGER_IRQ(irq);}
#define swiNum 71 // is unused on T3.2 and T3.6 others are possible
//
void SWI_init(void)
{	_VectorsRam[swiNum] = doProcessing;
	NVIC_ENABLE_IRQ(swiNum-16);
}
//
void SWI_trigger(int prio)
{
	TRIGGER_SWI(swiNum-16,prio*16) // do here the real triggering
}

///--------------------------- Periodic data simulation ----------------------------------------
void doSimulate(void);

void timer2_init(void)
{
	SIM_SCGC6 |= SIM_SCGC6_PIT;
	// turn on PIT
	PIT_MCR = 0x00;
}

void timer2_start(uint32_t frequency)
{//
	PIT_LDVAL2 = F_BUS/frequency; // setup timer 0 for (F_BUS/frequency) cycles
	PIT_TCTRL2 = 2; // enable Timer 2 interrupts
	NVIC_SET_PRIORITY(IRQ_PIT_CH2, 7*16); // 8 is normal priority (set in mk20dx128.c)
	NVIC_ENABLE_IRQ(IRQ_PIT_CH2);
	PIT_TCTRL2 |= 1; // start Timer 2
}

void timer2_stop(void)
{
	PIT_TCTRL2 &= ~1; // stop Timer 2

}

void pit2_isr(void)
{ //
	PIT_TFLG2=1;
	doSimulate();
}

///------------------------------------ Simulation ---------------------------------------------------
#define nf  129
#define mf (nf-1)
#define nx (1024 -mf)

short  sro[mf], sio[mf], nro[mf], nio[mf]; // for simulation
float sr[nx], si[nx];	// hold simulated data
float xr[nx], xi[nx], zr[nx], zi[nx];	// for processing

uint32_t t0,t1, t2,t3, t4,t5;
uint32_t procCount=0;
uint32_t dataCount=0;
void doSimulate(void)
{ 	int ii;

	t0 = micros();
	// copy data but change noise somewhat and scale signal and noise to reasonable values
	for(ii=0;ii<mf;ii++)
	{
#define shs 10
#define shn 10
		sr[dataCount*mf+ii] = (float) (sro[ii]>>shs) + (float)(nro[(dataCount+ii) % mf]>>shn);
		si[dataCount*mf+ii] = (float) (sio[ii]>>shs) + (float)(nio[(dataCount+ii) % mf]>>shn);
	}

	t1= micros();
	// if we have sufficient data, copy to working buffer and signal to loop() for (lower priority) processing
	if(++dataCount==(nx/mf))
	{
		for(int ii=0; ii < nx; ii++)
		{
			xr[ii]=sr[ii];
			xi[ii]=si[ii];
		}
		SWI_trigger(9);
		dataCount=0;
	}
}

void doProcessing(void)
{
	t2=micros();
	fft_filt_exec(zr,zi,xr,xi,nx, nf);
	procCount++;
    t3=micros();
}

void setup()
{
  // put your setup code here, to run once:
  //
	while(!Serial);
	while(Serial.available()); // clear input buffer
	Serial.printf("Starting SDR\n\r");

//	fft_filt_init(0.6, 0, 0.0);		// LP filter (3rf parameter irrelevant)
//	fft_filt_init(0.6, 1, 0.0);		// High pass filter (3rf parameter irrelevant)
//	fft_filt_init(0.75, 2, 0.05);	// band-pass filter (1st parameter is central, 3rd parameter is half bandwidth)
	fft_filt_init(0.75, 3, 0.05);	// stop-pass filter

	// initialize simulation
#define pi 3.14159f
#define no 48 // number of cycles in 128 samples, or frequency = no * (fsamp/128)
#define Sc 32767.0f

	for(int ii=0; ii < mf; ii++)
	{	float arg = 2.0f*pi*no*ii/128.0f;
		 sro[ii] = (short)(Sc*cosf(arg));
		 sio[ii] = (short)(Sc*sinf(arg));
		 nro[ii] = (short)((Sc*(float)rand()/(float) RAND_MAX) - 0.5);
		 nio[ii] = (short)((Sc*(float)rand()/(float) RAND_MAX) - 0.5);
	}
	SWI_init();
	//
	// following should be replaced with I2S driver
	timer2_init();
#define NINT 345	// interrupts / second (345 = 44100/128)
	timer2_start(NINT);

    t4=millis();
}

void loop()
{	int ii;

    if(procCount==50) // is about 1 sec of data
    {
    	t5=millis();
    	// stop processing
    	timer2_stop();
    	//
    	float ur[nx],ui[nx],vr[nx],vi[nx];

    	for(ii=0; ii<nx;ii++)
    	{	ur[ii]=xr[ii]; ui[ii]=xi[ii]; 	vr[ii]=zr[ii]; vi[ii]=zi[ii];
    	}

    	for(ii=0; ii<nx;ii++)
    	{	Serial.printf("%d, %d, %d, %d, %d\n\r", ii,	(int)(ur[ii]*1024.0f), (int)(ui[ii]*1024.0f),
    												(int)(vr[ii]*1024.0f), (int)(vi[ii]*1024.0f));
    		delay(10);
    	}
    	Serial.printf("\n\r%d: %d %d %d\n\r",procCount,t1-t0,t3-t2,t5-t4);
    	procCount=0;
    	while(1) delay(1000); // wait forever
    }
}

//#define DD4WH
#ifdef DD4WH
/*********************************************************************************************
 * (c) Frank DD4WH 2016_11_24
 *
 * FFT Fast Convolution = Digital Convolution
 * with overlap - save = overlap-discard
 * dynamically creates FIR coefficients
 *
 * in floating point 32bit
 * tested on Teensy 3.6
 *
 * audio queue optimized by Pete El Supremo 2016_10_27, thanks Pete!
 * final hint on the implementation came from Alberto I2PHD, thanks for that!
 * audio processing in float32_t with the NEW ARM CMSIS lib (see https://forum.pjrc.com/threads/38325-Excellent-results-with-Floating-Point-FFT-IFFT-Processing-and-Teensy-3-6?p=119797&viewfull=1#post119797),
 * and see here: https://forum.pjrc.com/threads/38325-Excellent-results-with-Floating-Point-FFT-IFFT-Processing-and-Teensy-3-6?p=120177&viewfull=1#post120177
 * MIT license
 */

#include <Audio.h>
#include <Wire.h>
#include <SPI.h>
#include <SD.h>
#include <Metro.h>
#include "font_Arial.h"
#include <ILI9341_t3.h>
#include <arm_math.h>
#include <arm_const_structs.h>

#include "Filter_coeffs.h"

#define BACKLIGHT_PIN 0
#define TFT_DC      20
#define TFT_CS      21
#define TFT_RST     32  // 255 = unused. connect to 3.3V
#define TFT_MOSI     7
#define TFT_SCLK    14
#define TFT_MISO    12

ILI9341_t3 tft = ILI9341_t3(TFT_CS, TFT_DC, TFT_RST, TFT_MOSI, TFT_SCLK, TFT_MISO);

Metro five_sec=Metro(2000); // Set up a Metro

// this audio comes from the codec by I2S2
AudioInputI2S            i2s_in;

AudioRecordQueue         Q_in_L;
AudioRecordQueue         Q_in_R;

AudioPlayQueue           Q_out_L;
AudioPlayQueue           Q_out_R;
AudioAnalyzeFFT256       myFFT;
AudioOutputI2S           i2s_out;

AudioConnection          patchCord1(i2s_in, 0, Q_in_L, 0);
AudioConnection          patchCord2(i2s_in, 1, Q_in_R, 0);

AudioConnection          patchCord5(Q_out_R,0,myFFT,0);

AudioConnection          patchCord3(Q_out_L, 0, i2s_out, 1);
AudioConnection          patchCord4(Q_out_R, 0, i2s_out, 0);

AudioControlSGTL5000     sgtl5000_1;     //xy=265.212

int idx_t = 0;
int idx = 0;
int64_t sum;
float32_t mean;
int n_L;
int n_R;
long int n_clear;

int peak[512];
int barm[512];

ulong samp_ptr = 0;
// bool FFT_state = false;

const int myInput = AUDIO_INPUT_LINEIN;

#define BUFFER_SIZE 128

// complex FFT with the new library CMSIS V4.5
const static arm_cfft_instance_f32 *S;

// create coefficients on the fly for custom filters
// and let the algorithm define which FFT size you need
// input variables by the user
//   first experiment: LOWPASS
// Fpass, Fstop, Astop (stopband attenuation)
// fixed: sample rate, scale = 1.0 ???
//
float32_t Fpass = 4860;
float32_t Fstop = 4980;
float32_t Astop = 90;
float32_t FIR_Coef[2049];
#define MAX_NUMCOEF 2049
#define K_2PI         2.0 * 3.14159265358979323846264338327950288419716939937510
#define K_PI          3.14159265358979323846264338327950288419716939937510
uint32_t m_NumTaps = 3;
const uint32_t FFT_L = 4096;
uint32_t FFT_length = FFT_L;
uint32_t N_BLOCKS;
const uint32_t N_B = FFT_L / 2 / BUFFER_SIZE;
float32_t float_buffer_L [BUFFER_SIZE * N_B];
float32_t float_buffer_R [BUFFER_SIZE * N_B];
float32_t FFT_buffer [FFT_L * 2] __attribute__ ((aligned (4)));
float32_t last_sample_buffer_L [BUFFER_SIZE * N_B];
float32_t last_sample_buffer_R [BUFFER_SIZE * N_B];

// complex iFFT with the new library CMSIS V4.5
const static arm_cfft_instance_f32 *iS;
float32_t iFFT_buffer [FFT_L * 2] __attribute__ ((aligned (4)));

// FFT instance for direct calculation of the filter mask
// from the impulse response of the FIR - the coefficients
const static arm_cfft_instance_f32 *maskS;
float32_t FIR_filter_mask [FFT_L * 2] __attribute__ ((aligned (4)));

int8_t first_block = 1;
uint32_t i = 0;

void setup() {
  Serial.begin(115200);
  delay(1000);

  // for the large FFT sizes we need a lot of buffers
  AudioMemory(60);

  // Enable the audio shield. select input. and enable output
  sgtl5000_1.enable();
  sgtl5000_1.inputSelect(myInput);
  sgtl5000_1.volume(0.65); // when I put this to 0.9, the digital noise of the Teensy is way too loud,
  // so 0.5 to 0.65 seems a good compromise
  sgtl5000_1.adcHighPassFilterDisable(); // does not help too much!

  pinMode( BACKLIGHT_PIN, OUTPUT );
  analogWrite( BACKLIGHT_PIN, 1023 );

  tft.begin();
  tft.setRotation( 3 );
  tft.fillScreen(ILI9341_BLACK);
  tft.setCursor(10, 1);
  tft.setTextSize(2);
  tft.setTextColor(ILI9341_WHITE);
  tft.setFont(Arial_14);
  tft.print("Fast convolution - FIR coeffs: ");

  // this routine does all the magic of calculating the FIR coeffs (Bessel-Kaiser window)
  calc_FIR_lowpass_coeffs (1.0, Astop, Fpass, Fstop, AUDIO_SAMPLE_RATE_EXACT);

  // adjust length of FFT according to number of coefficients calculated
      if (m_NumTaps > 1025)
      {
          FFT_length = 4096;
      }
      else if (m_NumTaps > 513)
      {
          FFT_length = 2048;
      }
      else if (m_NumTaps > 257)
      {
          FFT_length = 1024;
      }
      else if (m_NumTaps > 129)
      {
          FFT_length = 512;
      }
      else
      {
          FFT_length = 256;
      }
  // set no of BLOCKS
  N_BLOCKS = FFT_length / 2 / BUFFER_SIZE;

  // set to zero all other coefficients in coefficient array
  for(i = 0; i < MAX_NUMCOEF; i++)
  {
    Serial.print (FIR_Coef[i]); Serial.print(", ");
      if (i >= m_NumTaps) FIR_Coef[i] = 0.0;
  }
    // print number of FIR taps to impress others and show-off
    tft.print((int)m_NumTaps);
//    tft.print("  ");
//    tft.print(FFT_length);

 /****************************************************************************************
 *  init complex FFT
 ****************************************************************************************/
   switch (FFT_length) {
    case 256:
      S = &arm_cfft_sR_f32_len256;
      iS = &arm_cfft_sR_f32_len256;
      maskS = &arm_cfft_sR_f32_len256;
      break;
    case 512:
      S = &arm_cfft_sR_f32_len512;
      iS = &arm_cfft_sR_f32_len512;
      maskS = &arm_cfft_sR_f32_len512;
      break;
    case 1024:
      S = &arm_cfft_sR_f32_len1024;
      iS = &arm_cfft_sR_f32_len1024;
      maskS = &arm_cfft_sR_f32_len1024;
      break;
    case 2048:
      S = &arm_cfft_sR_f32_len2048;
      iS = &arm_cfft_sR_f32_len2048;
      maskS = &arm_cfft_sR_f32_len2048;
      break;
    case 4096:
      S = &arm_cfft_sR_f32_len4096;
      iS = &arm_cfft_sR_f32_len4096;
      maskS = &arm_cfft_sR_f32_len4096;
      break;
  }
 /****************************************************************************************
 *  Calculate the FFT of the FIR filter coefficients once to produce the FIR filter mask
 ****************************************************************************************/

// the FIR has exactly m_NumTaps and a maximum of (FFT_length / 2) + 1 taps = coefficients, so we have to add (FFT_length / 2) -1 zeros before the FFT
// in order to produce a FFT_length point input buffer for the FFT
    // copy coefficients into real values of first part of buffer, rest is zero
    for (i = 0; i < (FFT_length / 2) + 1; i+=2)
    {
        FIR_filter_mask[i * 2] = FIR_Coef [i];
        FIR_filter_mask[i * 2 + 1] =  0.0;
    }
// FFT of FIR_filter_mask
// perform FFT (in-place), needs only to be done once (or every time the filter coeffs change)
    arm_cfft_f32(maskS, FIR_filter_mask, 0, 1);
    Serial.println("Filter mask FFT done ");

///////////////////////////////////////////////////////////////77
// PASS-THRU only for TESTING
/////////////////////////////////////////////////////////////77
/*
// hmmm, unclear, whether [1,1] or [1,0] or [0,1] are pass-thru filter-mask-multipliers??
// empirically, [1,0] sounds most natural = pass-thru
  for(i = 0; i < FFT_length * 2; i+=2)
  {
        FIR_filter_mask [i] = 1.0; // real
        FIR_filter_mask [i + 1] = 0.0; // imaginary
  }
*/
 /****************************************************************************************
 *  begin to queue the audio from the audio library
 ****************************************************************************************/
    delay(100);
    Q_in_L.begin();
    Q_in_R.begin();
} // END SETUP

int16_t *sp_L;
int16_t *sp_R;

void loop() {

 elapsedMicros usec = 0;
/**********************************************************************************
 *  Get samples from queue buffers
 **********************************************************************************/
    // we have to ensure that we have enough audio samples: we need
    // N_BLOCKS = fft_length / 2 / BUFFER_SIZE
    // FFT1024 point --> n_blocks = 1024 / 2 / 128 = 4
    // when these buffers are available, read them in and perform
    // the FFT - cmplx-mult - iFFT
    //
    // are there at least N_BLOCKS buffers in each channel available ?
    if (Q_in_L.available() > N_BLOCKS && Q_in_R.available() > N_BLOCKS)
    {
// get audio samples from the audio  buffers and convert them to float
    for (i = 0; i < N_BLOCKS; i++)
    {
    sp_L = Q_in_L.readBuffer();
    sp_R = Q_in_R.readBuffer();

      // convert to float one buffer_size
     arm_q15_to_float (sp_L, &float_buffer_L[BUFFER_SIZE * i], BUFFER_SIZE); // convert int_buffer to float 32bit
     arm_q15_to_float (sp_R, &float_buffer_R[BUFFER_SIZE * i], BUFFER_SIZE); // convert int_buffer to float 32bit
     Q_in_L.freeBuffer();
     Q_in_R.freeBuffer();
    }

    // this is supposed to prevent overfilled queue buffers
    // rarely the Teensy audio queue gets a hickup
    // in that case this keeps the whole audio chain running smoothly
    if (Q_in_L.available() > N_BLOCKS + 4 || Q_in_R.available() > N_BLOCKS + 4)
    {
      Q_in_L.clear();
      Q_in_R.clear();
      n_clear ++; // just for debugging to check how often this occurs [about once in an hour of playing . . .]
    }

/**********************************************************************************
 *  Digital convolution
 **********************************************************************************/
//  basis for this was Lyons, R. (2011): Understanding Digital Processing.
//  "Fast FIR Filtering using the FFT", pages 688 - 694
//  numbers for the steps taken from that source
//  Method used here: overlap-and-save

// 4.) ONLY FOR the VERY FIRST FFT: fill first samples with zeros
      if(first_block) // fill real & imaginaries with zeros for the first BLOCKSIZE samples
      {
        for(i = 0; i < BUFFER_SIZE * 2 * N_BLOCKS; i++)
        {
          FFT_buffer[i] = 0.0;
        }
        first_block = 0;
      }
      else

// HERE IT STARTS for all other instances
// 6 a.) fill FFT_buffer with last events audio samples
      for(i = 0; i < BUFFER_SIZE * N_BLOCKS; i++)
      {
        FFT_buffer[i * 2] = last_sample_buffer_L[i]; // real
        FFT_buffer[i * 2 + 1] = last_sample_buffer_R[i]; // imaginary
      }

    // copy recent samples to last_sample_buffer for next time!
    arm_copy_f32(float_buffer_L, last_sample_buffer_L, BUFFER_SIZE * N_BLOCKS);
    arm_copy_f32(float_buffer_R, last_sample_buffer_R, BUFFER_SIZE * N_BLOCKS);

// 6. b) now fill audio samples into FFT_buffer (left channel: re, right channel: im)
      for(i = 0; i < BUFFER_SIZE * N_BLOCKS; i++)
      {
        FFT_buffer[FFT_length + i * 2] = float_buffer_L[i]; // real
        FFT_buffer[FFT_length + i * 2 + 1] = float_buffer_R[i]; // imaginary
      }

// perform complex FFT
// calculation is performed in-place the FFT_buffer [re, im, re, im, re, im . . .]
        arm_cfft_f32(S, FFT_buffer, 0, 1);

// complex multiply with filter mask
     arm_cmplx_mult_cmplx_f32 (FFT_buffer, FIR_filter_mask, iFFT_buffer, FFT_length);

/*  // just for testing: pass-thru !
     for(i = 0; i < FFT_length * 2; i++)
    {
      iFFT_buffer[i] = FFT_buffer[i];
    }
*/

// perform iFFT (in-place)
        arm_cfft_f32(iS, iFFT_buffer, 1, 1);

// our desired output is the real part (left channel) AND the imaginary part (right channel) of the second half of the FFT_buffer
      for(i = 0; i < BUFFER_SIZE * N_BLOCKS; i++)
      {
        float_buffer_L[i] = iFFT_buffer[FFT_length + (i * 2)];
        float_buffer_R[i] = iFFT_buffer[FFT_length + (i * 2) + 1];
//        float_buffer_L[i] = iFFT_buffer[(i * 2)];
//        float_buffer_R[i] = iFFT_buffer[(i * 2) + 1];
      }

    for (i = 0; i < N_BLOCKS; i++)
    {
      sp_L = Q_out_L.getBuffer();
      sp_R = Q_out_R.getBuffer();
      arm_float_to_q15 (&float_buffer_L[BUFFER_SIZE * i], sp_L, BUFFER_SIZE);
      arm_float_to_q15 (&float_buffer_R[BUFFER_SIZE * i], sp_R, BUFFER_SIZE);
      Q_out_L.playBuffer(); // play it !
      Q_out_R.playBuffer(); // play it !
    }

/**********************************************************************************
 *  PRINT ROUTINE FOR ELAPSED MICROSECONDS
 **********************************************************************************/

      sum = sum + usec;
      idx_t++;
      if (idx_t > 1000) {
          tft.fillRect(240,50,90,20,ILI9341_BLACK);
          tft.setCursor(240, 50);
          mean = sum / idx_t;
          tft.print (mean);
          Serial.print (mean);
          Serial.print (" microsec for ");
          Serial.print (N_BLOCKS);
          Serial.print ("  stereo blocks    ");
          Serial.println();
          Serial.print (" n_clear    ");
          Serial.println(n_clear);
          idx_t = 0;
          sum = 0;

      }

     }
/**********************************************************************************
 *  PRINT ROUTINE FOR AUDIO LIBRARY PROCESSOR AND MEMORY USAGE
 **********************************************************************************/
          if (five_sec.check() == 1)
    {
      Serial.print("Proc = ");
      Serial.print(AudioProcessorUsage());
      Serial.print(" (");
      Serial.print(AudioProcessorUsageMax());
      Serial.print("),  Mem = ");
      Serial.print(AudioMemoryUsage());
      Serial.print(" (");
      Serial.print(AudioMemoryUsageMax());
      Serial.println(")");
/*      tft.fillRect(100,120,200,80,ILI9341_BLACK);
      tft.setCursor(10, 120);
      tft.setTextSize(2);
      tft.setTextColor(ILI9341_WHITE);
      tft.setFont(Arial_14);
      tft.print ("Proc = ");
      tft.setCursor(100, 120);
      tft.print (AudioProcessorUsage());
      tft.setCursor(180, 120);
      tft.print (AudioProcessorUsageMax());
      tft.setCursor(10, 150);
      tft.print ("Mem  = ");
      tft.setCursor(100, 150);
      tft.print (AudioMemoryUsage());
      tft.setCursor(180, 150);
      tft.print (AudioMemoryUsageMax());
     */
      AudioProcessorUsageMaxReset();
      AudioMemoryUsageMaxReset();
    }
   spectrum();
}

 void spectrum() { // spectrum analyser code by rheslip - modified
     if (myFFT.available()) {
    int scale;
    scale = 2;
  for (int16_t x=2; x < 100; x+=2) {

     int bar = (abs(myFFT.output[x]) * scale);
     if (bar >180) bar=180;
     // this is a very simple IIR filter to smooth the reaction of the bars
     bar = 0.2 * bar + 0.8 * barm[x];
     if (bar > peak[x]) peak[x]=bar;
//     tft.drawFastVLine(x, 210-bar,bar, ILI9341_PURPLE);
     tft.drawFastVLine(x*2+10, 210-bar,bar, ILI9341_PINK);

     tft.drawFastVLine(x*2+10, 20, 210-bar-20, ILI9341_BLACK);

     tft.drawPixel(x*2+10,209-peak[x], ILI9341_YELLOW);

     if(peak[x]>0) peak[x]-=1;
     barm[x] = bar;
  }
  } //end if

   } // end void spectrum

void calc_FIR_lowpass_coeffs (float32_t Scale, float32_t Astop, float32_t Fpass, float32_t Fstop, float32_t Fsamprate)
{ // Wheatley, M. (2011): CuteSDR Technical Manual. www.metronix.com, pages 118 - 120, FIR with Kaiser-Bessel Window
    int n;
    float32_t Beta;
    float32_t izb;
    // normalized F parameters
    float32_t normFpass = Fpass / Fsamprate;
    float32_t normFstop = Fstop / Fsamprate;
    float32_t normFcut = (normFstop + normFpass) / 2.0; // lowpass filter 6dB cutoff
    // calculate Kaiser-Bessel window shape factor beta from stopband attenuation
    if (Astop < 20.96)
    {
      Beta = 0.0;
    }
    else if (Astop >= 50.0)
    {
      Beta = 0.1102 * (Astop - 8.71);
    }
    else
    {
      Beta = 0.5842 * powf((Astop - 20.96), 0.4) + 0.7886 * (Astop - 20.96);
    }
    // estimate number of taps
    m_NumTaps = (int32_t)((Astop - 8.0) / (2.285 * K_2PI * (normFstop - normFpass)) + 1.0);
    if (m_NumTaps > MAX_NUMCOEF)
    {
      m_NumTaps = MAX_NUMCOEF;
    }
    if (m_NumTaps < 3)
    {
      m_NumTaps = 3;
    }
    float32_t fCenter = 0.5 * (float32_t) (m_NumTaps - 1);
    izb = Izero (Beta);
    for (n = 0; n < m_NumTaps; n++)
    {
      float32_t x = (float32_t) n - fCenter;
      float32_t c;
      // create ideal Sinc() LP filter with normFcut
      if ((float32_t) n == fCenter)
      {
        c = 2.0 * normFcut;
      }
      else
      {
        c = (float32_t) sin(K_2PI * x * normFcut) / (K_PI * x);
      }
      x = ((float32_t) n - ((float32_t) m_NumTaps - 1.0) / 2.0) / (((float32_t) m_NumTaps - 1.0) / 2.0);
      FIR_Coef[n] = Scale * c * Izero (Beta * sqrt(1 - (x * x) )) / izb;
    }
} // END calc_lowpass_coeffs

float32_t Izero (float32_t x)
{
    float32_t x2 = x / 2.0;
    float32_t summe = 1.0;
    float32_t ds = 1.0;
    float32_t di = 1.0;
    float32_t errorlimit = 1e-9;
    float32_t tmp;
    do
    {
        tmp = x2 / di;
        tmp *= tmp;
        ds *= tmp;
        summe += ds;
        di += 1.0;
    }   while (ds >= errorlimit * summe);
    return (summe);
}  // END Izero
#endif
