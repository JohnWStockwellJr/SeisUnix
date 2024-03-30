/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */


/* SUATTRIBUTES:  $Revision: 1.35 $ ; $Date: 2016/05/09 16:43:24 $	*/


#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 									",
" SUATTRIBUTES - instantaneous trace ATTRIBUTES 			",
" 									",
" suattributes <stdin >stdout mode=amp					",
" 									",
" Required parameters:							",
" 	none								",
" 									",
" Optional parameter:							",
" 	mode=amp	output flag 					",
" 	       		=amp envelope traces				",
" 	       		=phase phase traces				",
" 	       		=uphase unwrapped phase traces			",
" 	       		=freq frequency traces				",
"			=freqw Frequency Weighted Envelope		",
"			=thin  Thin-Bed (inst. freq - average freq)	",
"			=bandwith Instantaneous bandwidth		",
"			=normamp Normalized Phase (Cosine Phase)	",
" 	       		=fdenv 1st envelope traces derivative		",
" 	       		=sdenv 2nd envelope traces derivative		",
" 	       		=q Ins. Q Factor				",
"	unwrap=		default unwrap=0 for mode=phase			",
" 			default unwrap=1 for freq, uphase, freqw, Q	",
" 			dphase_min=PI/unwrap				",
"	wint=		windowing for freqw				",
"			windowing for thin				",
"			default=1 					",
" 			o--------o--------o				",
" 			data-1	data	data+1				",
" 									",
" Notes:								",
" This program performs complex trace attribute analysis. The first three",
" attributes, amp,phase,freq are the classical Taner, Kohler, and	",
" Sheriff, 1979 attributes.						",
" 									",
" 									",
" The unwrap parameter is active only for mode=freq and mode=phase. The	",
" quantity dphase_min is the minimum change in the phase angle taken to be",
" the result of phase wrapping, rather than natural phase variation in the",
" data. Setting unwrap=0 turns off phase-unwrapping altogether. Choosing",
" unwrap > 1 makes the unwrapping function more sensitive to phase changes.",
" Setting unwrap > 1 may be necessary to resolve higher frequencies in	",
" data (or sample data more finely). The phase unwrapping is crude. The ",
" differentiation needed to compute the instantaneous frequency		",
" freq(t)= d(phase)/dt is a simple centered difference.			",
"	 					       			",
" The mode=uphase generates uwrapped phase traces by integrating the	",
" instantaneous amplitude traces.		       			",
"	 					       			",
" Examples:								",
" suvibro f1=10 f2=50 t1=0 t2=0 tv=1 | suattributes2 mode=amp | ...	",
" suvibro f1=10 f2=50 t1=0 t2=0 tv=1 | suattributes2 mode=phase | ...	",
" suvibro f1=10 f2=50 t1=0 t2=0 tv=1 | suattributes2 mode=freq | ...	",
" suplane | suattributes mode=... | supswigb |...       		",

NULL};

/* Credits:
 *	CWP: Jack Cohen
 *      CWP: John Stockwell (added freq and unwrap features)
 *	UGM (Geophysics Students): Agung Wiyono
 *           email:aakanjas@gmail.com (others)
 *	CSM: Kylee Brown and Steven Rennolet, Senior Design,
 *	     updates to instanteous phase, instantaneous frequency,
 *	     first time derivative of the envelope, second time derivative
 *	     of the envelope, instantaneous quality factor, and thin bed
 *	     indicator
 *					
 *
 * Algorithm:
 *	c(t) = hilbert_tranform_kernel(t) convolved with data(t)  
 *
 *  amp(t) = sqrt( c.re^2(t) + c.im^2(t))
 *  phase(t) = arctan( c.im(t)/c.re(t))
 *  freq(t) = d(phase)/dt
 *
 * Reference: 
 *  Taner, M. T., Koehler, A. F., and  Sheriff R. E.   "Complex seismic trace 
 *      analysis", Geophysics,  vol.44, p. 1041-1063, 1979
 *  Chopra, S. and K.  Marfurt, 2005, A historical perspective, Geophysics,
 *      vol. 70, no. 5, p.3SO-295SO, Society of Exploration Geophysicists.
 *  Barnes, A. E, 1992, The calculation of instantaneous frequency and 
 *      instantaneous bandwidth, Geophysics, vol. 57, no. 11, p. 1520-1524,
 *      Society of Exploration Geophysicists.
 *
 * Trace header fields accessed: ns, trid
 * Trace header fields modified: d1, trid

 */
/**************** end self doc ********************************/

#define	AMP		 1
#define	ARG		 2
#define	FREQ		 3
#define BANDWIDTH	 4
#define NORMAMP		 5
#define FREQW		 6
#define THIN		 7
#define FENV		 8
#define SENV		 9
#define Q		 10
#define UPHASE		 11

/* function prototype of functions used internally */
void unwrap_phase(int n, float unwrap, float *phase);
void differentate1d(int n, float h, float *f);
void twindow(int nt, int wtime, float *data);

segy tr;

int
main(int argc, char **argv)
{
	cwp_String mode;	/* display: real, imag, amp, arg	*/
	int imode=AMP;		/* integer abbrev. for mode in switch	*/
	register complex *ct;	/* complex trace			*/
	int nt;			/* number of points on input trace	*/
	float dt;		/* sample spacing			*/
	float *data;		/* array of data from each trace	*/
	float *hdata;		/* array of Hilbert transformed data	*/
	float unwrap;		/* PI/unwrap=min dphase assumed to by wrap*/
	int wint;		/* n time sampling to window */
	cwp_Bool seismic;	/* is this seismic data?		*/
	int ntout;
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
	

	/* Get info from first trace */
	if (!gettr(&tr)) err("can't get first trace");
	nt = tr.ns;
	dt = ((double) tr.dt)/1000000.0;
	ntout = nt + nt -1;

	/* check to see if data type is seismic */
	seismic = ISSEISMIC(tr.trid);

	if (!seismic)
		warn("input is not seismic data, trid=%d", tr.trid);

	/* Get mode; note that imode is initialized to AMP */
	if (!getparstring("mode", &mode))	mode = "amp";

	if      (STREQ(mode, "phase"))  imode = ARG;
	else if (STREQ(mode, "freq"))	imode = FREQ;
	else if (STREQ(mode, "uphase"))  imode = UPHASE;
	else if (STREQ(mode, "bandwidth")) imode = BANDWIDTH;
	else if (STREQ(mode, "normamp")) imode = NORMAMP;
	else if (STREQ(mode, "freqw")) imode = FREQW;
	else if (STREQ(mode, "thin")) imode = THIN;
	else if (STREQ(mode, "fdenv")) imode = FENV;
	else if (STREQ(mode, "sdenv")) imode = SENV;
	else if (STREQ(mode, "q")) imode = Q;
	else if (!STREQ(mode, "amp"))
		err("unknown mode=\"%s\", see self-doc", mode);

	/* getpar value of unwrap */
	switch(imode) {
	case FREQ:
		if (!getparfloat("unwrap", &unwrap))	unwrap=1;
	break;
	case UPHASE:
		if (!getparfloat("unwrap", &unwrap))	unwrap=1;
	break;
	case Q:
		if (!getparfloat("unwrap", &unwrap))	unwrap=1;
	break;
	case FREQW:
		if (!getparfloat("unwrap", &unwrap))	unwrap=1;
		if (!getparint("wint", &wint))	wint=3; 
	break;
	case THIN:
		if (!getparfloat("unwrap", &unwrap))	unwrap=1;
		if (!getparint("wint", &wint))	wint=3;
	break;
	case ARG:
		if (!getparfloat("unwrap", &unwrap))	unwrap=0;
	break;
	}

	/* allocate space for data and hilbert transformed data, cmplx trace */
	data = ealloc1float(nt);
	hdata = ealloc1float(nt);
	ct = ealloc1complex(nt);


	/* Loop over traces */
	do {
		register int i;

		/* Get data from trace */
		for (i = 0; i < nt; ++i)  data[i] = tr.data[i];

		
		/* construct quadrature trace with hilbert transform */
		hilbert(nt, data, hdata);

		/* build the complex trace */
		for (i = 0; i < nt; ++i)  ct[i] = cmplx(data[i],hdata[i]);

		/* Form absolute value, phase, or frequency */
		switch(imode) {
		case AMP:
			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;
				tr.data[i] = sqrt(re*re + im*im);
			}
			
			/* set trace id */
			tr.trid = ENVELOPE;
		break;
		case ARG:
		{
			float *phase = ealloc1float(nt);

			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;
				if (re*re+im*im)  phase[i] = atan2(im, re);
				else              phase[i] = 0.0;
			}

			/* phase unwrapping */
			/* default unwrap=0 for this mode */
			if (unwrap!=0) unwrap_phase(nt, unwrap, phase);
			
			/* write phase values to tr.data */
			for (i = 0; i < nt; ++i) tr.data[i] = phase[i];
			
			/* set trace id */
			tr.trid = INSTPHASE;
		}
		break;
		case FREQ:
		{
			float *phase = ealloc1float(nt);
			float	fnyq = 0.5 / dt;

			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;
				if (re*re+im*im) {
					phase[i] = atan2(im, re);
				} else {
					phase[i] = 0.0;
				}
				
			}

			/* unwrap the phase */
			if (unwrap!=0) unwrap_phase(nt, unwrap, phase);

			/* compute freq(t)=dphase/dt */
			differentate1d(nt, 2.0*PI*dt, phase);
			
			/* correct values greater nyquist frequency */
			for (i=0 ; i < nt; ++i)	{
				if (phase[i] > fnyq)
					phase[i] = 2 * fnyq - phase[i];
			}
                                        
			/* write freq(t) values to tr.data */
			for (i=0 ; i < nt; ++i) tr.data[i] = phase[i];

			/* set trace id */
			tr.trid = INSTFREQ;
		}
		break;
		case UPHASE:
		{
			float *phase = ealloc1float(nt);
			float	fnyq = 0.5 / dt;

			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;
				if (re*re+im*im) {
					phase[i] = atan2(im, re);
				} else {
					phase[i] = 0.0;
				}
				
			}

			/* unwrap the phase */
			if (unwrap!=0) unwrap_phase(nt, unwrap, phase);

			/* compute freq(t)=dphase/dt */
			differentate1d(nt, 2.0*PI*dt, phase);
			
			/* correct values greater nyquist frequency */
			for (i=0 ; i < nt; ++i)	{
				if (phase[i] > fnyq)
					phase[i] = 2 * fnyq - phase[i];
			}
			/* integrate instantaneous frequency values */
			/* and write unwrapped phase values to tr.data */
			for (i = 1; i < nt; ++i) {
				tr.data[0] = phase[0];
                                tr.data[i] += phase[i-1];
                        }

			/* set trace id */
			tr.trid = INSTPHASE;
		}
		break;
		case FREQW:
		{
			float	fnyq = 0.5 / dt;
			float *freqw = ealloc1float(nt);
			float *phase = ealloc1float(nt);
			float *envelop = ealloc1float(nt);
			float *envelop2 = ealloc1float(nt);
			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;
				if (re*re+im*im) {
					phase[i] = atan2(im, re);
					} else {
						phase[i] = 0.0;
						}
				envelop[i] = sqrt(re*re + im*im);
			}

			/* unwrap the phase */
			if (unwrap!=0) unwrap_phase(nt, unwrap, phase);

			/* compute freq(t)=dphase/dt */
			differentate1d(nt, 2.0*PI*dt, phase);
			
			/* correct values greater nyquist frequency */
			for (i=0 ; i < nt; ++i)	{
				if (phase[i] > fnyq)
					phase[i] = 2 * fnyq - phase[i];
			envelop2[i]=envelop[i]*phase[i];
			}
			twindow(nt, wint, envelop);
			twindow(nt, wint, envelop2);
			/* correct values greater nyquist frequency */
			for (i=0 ; i < nt; ++i) {
			freqw[i] = (envelop[i] == 0.0) ? 0.0 :envelop2[i]/envelop[i];
			}
			/* write freq(t) values to tr.data */
			for (i=0 ; i < nt; ++i) tr.data[i] = freqw[i];
			
			/* set trace id */
			tr.trid = INSTFREQ;
		}
		break;
		case THIN:
		{
			float	fnyq = 0.5 / dt;
			float *phase = ealloc1float(nt);
			float *freqw = ealloc1float(nt);
			float *phase2 = ealloc1float(nt);



			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;

				if (re*re+im*im) {
					phase[i] = atan2(im, re);
				} else {
					phase[i] = 0.0;
				}
			}

			/* unwrap the phase */
			if (unwrap!=0) unwrap_phase(nt, unwrap, phase);

			/* compute freq(t)=dphase/dt */
			differentate1d(nt, 2.0*PI*dt, phase);

			/* correct values greater nyquist frequency */
			for (i=0 ; i < nt; ++i)	{
				if (phase[i] > fnyq)
					phase[i] = 2 * fnyq - phase[i];
					phase2[i]= 2 * fnyq - phase[i];
			}
			/* Do windowing for Average Ins . Freq over wint*/
			twindow(nt, wint, phase2);

			for (i=0 ; i < nt; ++i)	{
				freqw[i] = phase[i] - phase2[i];
			/*	if (abs(freqw[i]) > fnyq)
				freqw[i] = 2 * fnyq - freqw[i];
			*/
			/* write Thin-Bed(t) values to tr.data */
				tr.data[i] = freqw[i];
				}
			/* set trace id */
			tr.trid = INSTFREQ;
		}
		break;
		case BANDWIDTH:
		{
			float *envelop = ealloc1float(nt);
			float *envelop2 = ealloc1float(nt);

		/* Bandwidth (Barnes 1992)

		          |d(envelope)/dt|
		band =abs |--------------|
		          |2 PI envelope |
	 	*/

			for (i = 0; i < nt; ++i) {
				float er = ct[i].r;
				float em = ct[i].i;
				envelop[i] = sqrt(er*er + em*em);
				envelop2[i]=sqrt(er*er + em*em);

			}
				differentate1d(nt, dt, envelop);

				for (i = 0; i < ntout; ++i) {
				   if (2.0*PI*envelop2[i]!=0.0) {
					tr.data[i] = ABS(envelop[i]/(2.0*PI*envelop2[i]));
				   } else {
				        tr.data[i]=0.0;
				   }
				}
				tr.trid = ENVELOPE;
		}
		break;
		case NORMAMP:
		{
			float phase;
			float *na = ealloc1float(nt);
			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;
				if (re*re+im*im)  phase = atan2(im, re);
				else              phase = 0.0;
				na[i] = cos(phase);
			}
			for (i=0 ; i < nt; ++i) tr.data[i] = na[i];
			
			/* set trace id */
			tr.trid = INSTPHASE;
			}
		break;
		case FENV:
		{
			float *amp = ealloc1float(nt);
			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;
				amp[i] = sqrt(re*re + im*im);
			}
		/*conv(nt, 0, envelop, nt, 0, time, ntout, 0, ouput);*/

		differentate1d(nt, 2.0*PI*dt, amp);
		for (i=0 ; i < nt; ++i) tr.data[i] = amp[i];
			/* set trace id */
			tr.trid = ENVELOPE;
		}
		break;
		case SENV:
		{
			float *amp = ealloc1float(nt);
			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;
				amp[i] = sqrt(re*re + im*im);
			}

		differentate1d(nt, 2.0*PI*dt, amp);
		differentate1d(nt, 2.0*PI*dt, amp);
		for (i=0 ; i < nt; ++i) tr.data[i] = amp[i];
			/* set trace id */
			tr.trid = ENVELOPE;
		}
		break;

		case Q:
		{
			float *envelop = ealloc1float(nt);
			float *envelop2 = ealloc1float(nt);
			float *phase = ealloc1float(nt);
			float	fnyq = 0.5 / dt;

		/* Banswith (Barnes 1992)

		        -PI Freq(t) d(envelope)/dt
		band =  --------------------------
		                 envelope(t)
	 	*/

			for (i = 0; i < nt; ++i) {
				float re = ct[i].r;
				float im = ct[i].i;
				envelop[i] = sqrt(re*re + im*im);
				envelop2[i]=sqrt(re*re + im*im);
				if (re*re+im*im) {
					phase[i] = atan2(im, re);
				} else {
					phase[i] = 0.0;
				}

			}
			/* get envelope diff */
			differentate1d(nt, dt, envelop);
			/* unwrap the phase */
			if (unwrap!=0) unwrap_phase(nt, unwrap, phase);
			/* compute freq(t)=dphase/dt */
			differentate1d(nt, 2.0*PI*dt, phase);

			for (i=0 ; i < nt; ++i)	{
				if (phase[i] > fnyq)
					phase[i] = 2 * fnyq - phase[i];
			}

			for (i = 0; i < ntout; ++i) {
				if (envelop[i]!=0.0)
				tr.data[i] = -1*PI*phase[i]*envelop2[i]/envelop[i];
				else
				tr.data[i]=0.0;
				}
				tr.trid = INSTFREQ;
		}
		break;
		default:
			err("%s: mysterious mode=\"%s\"", __LINE__, mode);
		}


		tr.d1 = dt;   /* for graphics */
		puttr(&tr);

	} while (gettr(&tr));


	return(CWP_Exit());
}


void differentate1d(int n, float h, float *f)
/************************************************************************
differentate1d - compute the 1st derivative of a function f[]
************************************************************************
Input:
n		number of samples
h		sample rate
f		array[n] of input values

Output:
f		array[n], the derivative of f
************************************************************************
Notes:
This is a simple 2 point centered-difference differentiator.
The derivatives at the endpoints are computed via 2 point leading and
lagging differences. 
************************************************************************
Author: John Stockwell, CWP, 1994
************************************************************************/
{
	int i;	
	float *temp;
	float h2=2*h;

	/* allocate space in temporary vector */
	temp = ealloc1float(n);

	/* do first as a leading difference */
	temp[0] = (f[1] - f[0])/h;

	/* do the middle values as a centered difference */
	for (i=1; i<n-1; ++i) temp[i] = (f[i+1] - f[i-1])/h2;

	/* do last value as a lagging difference */
	temp[n-1] = (f[n-1] - f[n-2])/h;

	for (i=0 ; i < n ; ++i) f[i] = temp[i];

	free1float(temp);
}

void unwrap_phase(int n, float w, float *phase)
/************************************************************************
unwrap_phase - unwrap the phase
*************************************************************************
Input:
n		number of samples
w		unwrapping flag; returns an error if w=0
phase		array[n] of input phase values

Output:
phase		array[n] of output phase values
*************************************************************************
Notes:
The phase is assumed to be continuously increasing. The strategy is
to look at the change in phase (dphase) with each time step. If it is larger
than PI/w, then use the previous value of dphase. No attempt is
made at smoothing the dphase curve.
*************************************************************************
Author: John Stockwell, CWP, 1994
************************************************************************/
{
	int i;
	float pibyw=0.0;
	float *dphase;
	float *temp;

	/* prevent division by zero in PI/w */
	if (w==0)  err("wrapping parameter is zero");
	else       pibyw = PI/w;

	/* allocate space */
	dphase = ealloc1float(n);
	temp = ealloc1float(n);

	/* initialize */
	temp[0]=phase[0];
	dphase[0]=0.0;

	/* compute unwrapped phase at each time step */
	for (i = 1; i < n; ++i) {

		/* compute jump in phase */
		dphase[i] = ABS(phase[i] - phase[i-1]);

		/* if dphase >= PI/w, use previous dphase value */
		if (ABS(dphase[i] - dphase[i-1]) >= pibyw )
			dphase[i] = dphase[i-1];

		/* sum up values in temporary vector */
		temp[i] = temp[i-1] + dphase[i];
	}

	/* assign values of temporary vector to phase[i] */
	for (i=0; i<n; ++i) phase[i] = temp[i];

	/* free space */
	free1float(temp);
	free1float(dphase);
}

void twindow(int nt, int wtime, float *data)
/************************************************************
twindow - simple time gating
*************************************************************
Input:
nt	number of time samples
wtime	= n*dt   where n are integer ex=1,2,3,4,5,...
          wtime=3 as default
 used for Frequency Weighted and Thin-bed attributes
*************************************************************
Author:	UGM (Geophysics Students): Agung Wiyono, 2005
************************************************************/
{
	float val;
	float *temp;
	int i;
	float sum;
	int nwin;
	
	nwin=2*wtime+1;
	temp = ealloc1float(nt);
	sum=0.0;
	for (i = 0; i< wtime+1; ++i) {
		val = data[i];
		sum +=val;
	}
	/* weighted */
	temp[0] = sum/nwin;
	
	/* dt<wtime */
	for (i = 1; i < wtime; ++i) {
		val = data[i+wtime];
		sum+=val;
		++nwin;
		temp[i] = sum/nwin;
		}
	/*wtime<dt<dt-wtime */
	for (i = wtime ; i < nt-wtime; ++i) {
		val = data[i+wtime];
		sum += val;
		val = data[i-wtime];
		sum -=val;
		temp[i] = sum/nwin;
	}

	/*dt-wtime<dt*/
	for (i = nt - wtime; i < nt; ++i) {
		val = data[i-wtime];
		sum -= val;
		--nwin;
		temp[i] = sum/nwin;
	}
	
	
	for (i=0; i<nt; ++i) data[i] = temp[i];

	/* Memori free */
	free1float(temp);
}



