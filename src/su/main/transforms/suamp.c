/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.			*/

/* SUAMP: $Revision: 1.24 $ ; $Date: 2015/08/07 22:31:17 $		*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 									",
" SUAMP - output amp, phase, real or imag trace from			",
" 	(frequency, x) domain data					",
" 									",
" suamp <stdin >stdout mode=amp						",
" 									",
" Required parameters:							",
" none									",
" Optional parameter:							",
" mode=amp	output flag		 				",
" 		=amp	output amplitude traces				",
" 		=logamp	output log(amplitude) traces			",
" 			=phase	output phase traces			",
" 			=ouphase output unwrapped phase traces (oppenheim)",
" 			=suphase output unwrapped phase traces (simple)	",
" 			=real	output real parts			",
" 	     	=imag	output imag parts	 			",
" jack=0	=1  divide value at zero frequency by 2   		",
"		(operative only for mode=amp)				",
" 									",
" .... phase unwrapping options	..... 					",
" unwrap=1	 |dphase| > pi/unwrap constitutes a phase wrapping	",
"			(operative only for mode=suphase)		",
" trend=1	remove linear trend from the unwrapped phase		",
" zeromean=0	assume phase(0)=0.0, else assume phase is zero mean	",
" smooth=0	apply damped least squares smoothing to unwrapped phase ",
" r=10.0	    ... damping coefficient, only active when smooth=1	",
"  									",
" Notes:								",
" 	The trace returned is half length from 0 to Nyquist. 		",
" 									",
" Example:								",
" 	sufft <data | suamp >amp_traces					",
" Example: 								",
"	sufft < data > complex_traces					",
" 	 suamp < complex_traces mode=real > real_traces			",
" 	 suamp < complex_traces mode=imag > imag_traces			",
"  									",
" Note: the inverse of the above operation is: 				",
"	suop2 real_traces imag_traces op=zipper > complex_traces	",
"  									",
" Note: Explanation of jack=1 						",
" The amplitude spectrum is the modulus of the complex output of	",
" the fft. f(0) is thus the average over the range of integration	",
" of the transform. For causal functions, or equivalently, half		",
" transforms, f(0) is 1/2 of the average over the full range.		",
" Most oscillatory functions encountered in wave applications are	",
" zero mean, so this is usually not an issue.				",
"  									",
" Note: Phase unwrapping: 						",
"  									",
" The mode=ouphase uses the phase unwrapping method of Oppenheim and	",
" Schaffer, 1975. 							",
" The mode=suphase generates unwrapped phase assuming that jumps	",
" in phase larger than pi/unwrap constitute a phase wrapping.		",
NULL};

/* Credits:
 *	CWP: Shuki Ronen, Jack K. Cohen c.1986
 *
 * Notes:
 *	If efficiency becomes important consider inverting main loop
 *      and repeating extraction code within the branches of the switch.
 *
 * Trace header fields accessed: ns, trid
 * Trace header fields modified: ns, trid
 */
/**************** end self doc ***********************************/

/* function prototype of a subroutine used internally */
void smooth_1(float *x, float *z, float r, int n);



#define	REAL	1
#define	IMAG	2
#define	AMP	3
#define	ARG	4
#define	SUPHASE	5
#define	OUPHASE	6
#define LOGAMP  7

segy tr;

int
main(int argc, char **argv)
{
	cwp_String mode;	/* display: real, imag, amp, arg	*/
	int imode=AMP;		/* integer abbrev. for mode in switch	*/
	int nfby2;		/* nf/2					*/
	register float *xr;	/* real part of trace			*/
	register float *xi;	/* imaginary part of trace		*/
	int jack=0;		/* flag for special treatment of zero omega */
	float unwrap;		/* PI/unwrap = min dphase assumed to wrap */

	int zeromean;		/* =0 assume phase(0)=0.0 ; =1  zero mean*/
	int verbose=0;		/* =1 chatty ; =0 silent		*/
	int trend ;		/* remove linear trend in phase unwrapping */

	float df;		/* frequency sampling interval		*/
	int nf;			/* number of samples on input trace	*/


	/* smoothing */
	int smooth=0;   /* =1 perform least squares damping */
	float r=0;      /* smoothing parameter, active only for smooth=1 */


	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);


	/* Get information from first trace */
	if (!gettr(&tr)) err("can't get first trace");
	if (tr.trid != FUNPACKNYQ)
		err("input not complex freq data, trid=%d", tr.trid);

	nf = tr.ns; /* always even--includes 1 or 2 "wasted" slots */
	nfby2 = nf/2;
	if (!getparint("verbose",&verbose))	verbose=0 ;

	df=tr.d1;
	if (!df) {
		float dt = ((double) tr.dt)/1000000.0;
		if (!dt) {
			dt = .004;
			if (verbose) warn("dt or d1 not set, assumed to be .004");
		}
		df = 1.0/nf*dt;
	}


	/* Get mode; note that imode is initialized to AMP */
	if (!getparstring("mode", &mode))	mode = "amp";
	if (!getparint("jack",&jack))	jack = 0;

	 /* get smoothing */
	if(!getparint("smooth",&smooth))	 smooth=0;
	if(!getparfloat("r",&r))		   r=10.0;



	if      (STREQ(mode, "phase")) imode = ARG;
	else if (STREQ(mode, "logamp"))  imode = LOGAMP;
	else if (STREQ(mode, "ouphase"))  imode = OUPHASE;
	else if (STREQ(mode, "suphase"))  imode = SUPHASE;
	else if (STREQ(mode, "real"))  imode = REAL;
	else if (STREQ(mode, "imag"))  imode = IMAG;
	else if (!STREQ(mode, "amp"))

		err("unknown mode=\"%s\", see self-doc", mode);


	if(imode==OUPHASE || imode==SUPHASE) {
		if (!getparint("trend",&trend))	trend = 1;
		if (!getparint("zeromean",&zeromean))	zeromean = 0;
		if (!getparfloat("unwrap",&unwrap))	unwrap = 1;
	}
		

	checkpars();

	/* Allocate arrays for real and imaginary parts */
	xr = ealloc1float(nfby2);
	xi = ealloc1float(nfby2);
	
		
	float *xout=NULL;
	float *phase=NULL;
	if(imode == OUPHASE || imode == SUPHASE) {
  	   /* allocate space for the phase */
	   phase = ealloc1float(nfby2);
	   xout = ealloc1float(nfby2);
	}

	/* Main loop over traces */
	do {
		register int i;

		/* Separate complex trace into real and imag parts */
		for (i = 0; i < nfby2; ++i) {
			xr[i] = tr.data[2*i];
			xi[i] = tr.data[2*i+1];
		}

		/* Compute the desired half-length trace */
		switch(imode) {
		case REAL:
			for (i = 0; i < nfby2; ++i) {
				tr.data[i] = xr[i];
			}
			tr.trid = REALPART;
		break;
		case IMAG:
			for (i = 0; i < nfby2; ++i) {
				tr.data[i] = xi[i];
			}
			tr.trid = IMAGPART;
		break;
		case AMP:
		{
	 		register float re, im;
	
			re = xr[0];
			im = xi[0];
			if (jack) {
				tr.data[0] = (float) sqrt (re * re + im * im) / 2.0;
			} else {
				tr.data[0] = (float) sqrt (re * re + im * im);
			}
				for (i = 1; i < nfby2; ++i) {
					re = xr[i];
					im = xi[i];
					tr.data[i] = (float) sqrt (re * re + im * im);
			}
			tr.trid = AMPLITUDE;
		}
		break;
		case LOGAMP:
		{
	 		register float re, im, tmpamp;
	
			re = xr[0];
			im = xi[0];
			if (jack) {
				tr.data[0] = (float) sqrt (re * re + im * im) / 2.0;
			} else {
				tr.data[0] = (float) sqrt (re * re + im * im);
			}
				for (i = 1; i < nfby2; ++i) {

					re = xr[i];
					im = xi[i];

					tmpamp=(float) sqrt (re * re + im * im);

					if (CLOSETO(tmpamp,0.0)){
						tr.data[i]=0.0;
					} else {
						tr.data[i] = log(tmpamp);
					}
			}
			tr.trid = LOGAMPLITUDE;
		}
		break;
		case ARG:
			for (i = 0; i < nfby2; ++i) {
				float re = xr[i];
				float im = xi[i];

				if (re*re+im*im) {
					tr.data[i] = atan2(im, re);
				} else {
					tr.data[i] = 0.0;
				}
			}
			tr.trid = PHASE;
		break;
		case SUPHASE:
		{
			for (i = 0; i < nfby2; ++i) {
				float re = xr[i];
				float im = xi[i];

				xout[i] = i;

				if (re*re+im*im) {
					phase[i] = atan2(im, re);
				} else {
					phase[i] = 0.0;
				}
			}

			/* unwrap the phase */
			if (unwrap)
				simple_unwrap_phase(nfby2, trend, zeromean, unwrap, phase);
 
			/* smooth */
			if (smooth) smooth_1(xout,phase,r,nfby2);

			
			/* write unwrapped phase */
			for ( i = 0; i < nfby2; ++i) tr.data[i] = phase[i];

			tr.trid = PHASE;
		}
		break;
		case OUPHASE:
		{
			memset((void *) phase,  0, nfby2*FSIZE);
			/* unwrap the phase */
			if (unwrap)
				oppenheim_unwrap_phase(nfby2, trend, zeromean, df, xr, xi, phase);

			for (i=0; i< nfby2; ++i) {
				 xout[i] = i;
			}

			/* smooth */
			if (smooth) smooth_1(xout,phase,r,nfby2);

			/* write unwrapped phase */
			for ( i = 0; i < nfby2; ++i) tr.data[i] = phase[i];

			tr.trid = PHASE;
		}
		break;
		default:
			err("mysterious mode=\"%s\"", mode);
		}

		/* Output the half-length trace */
		tr.ns = nfby2;
		puttr(&tr);

	} while (gettr(&tr));
	
	if(imode == OUPHASE || imode == SUPHASE) {
  	   /* free space for the phase */
	   free1float(phase);
	   free1float(xout);
	}
	free1float(xr);
	free1float(xi);

	return(CWP_Exit());
}
