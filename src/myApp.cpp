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
