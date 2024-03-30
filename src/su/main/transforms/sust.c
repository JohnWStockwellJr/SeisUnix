/* Copyright (c) Colorado School of Mines, 2016.*/
/* All rights reserved.                       */

/* SUST: $Revision: $ ; $Date: 2016/05/16 23:35:04 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>
#include <time.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"									",
" SUST -  Outputs a time-frequency representation of seismic data via",
"	   the Stockwell transform (S- transform)			",
"									",
"    sust <stdin >stdout [optional parameters]			",
"									",
" Required parameters:					 		",
"	if dt is not set in header, then dt is mandatory		",
"									",
" Optional parameters:							",
"	dt=(from header)	time sampling interval (sec)		",
"	fmin=0			minimum frequency of filter array (hz)	",
"	fmax=NYQUIST 		maximum frequency of filter array (hz)	",
"	verbose=0		=1 supply additional info		",
"									",
" Notes: The S transform provide a time dependent frequency distribution ",
" of the signal. It is similar to the Gabor transform which which utilizes" 
" a Gaussian window for for spectral location. In the S transform the    ",
" Gaussian window is scalable whith the frequency which provides a better",
"  time freqency resolution    			                ",
"									",
" Examples:								",
"    suvibro | sust | suximage					",
"    suvibro | sust | suxmovie n1= n2= n3= 				",
"     (because suxmovie scales it's amplitudes off of the first panel,  ",
"      may have to experiment with the wclip and bclip parameters	",
"    suvibro | sust | supsimage | ... ( your local PostScript utility)",
"									",
NULL};

/* Credits:
 *
 *	BRGM : Adnand BITRI 20/05/2016 
 *
 * References: 	Stockwell,R,G.,Mansinha, L., and Lowe, R, 1996, Localization
 *		of the comlex spectrum: The S transform,
 *		IEEE Trans. Signal Process, 44, No, 4. 2957-2962.
 *
 *              Brown, R. A., Lauzon, M. L., and Frayen, R., 2010
 *		A general description of linear time-frequency transforms
 *              and formulation of a fast, invertible transform that  
 *		samples the continuous S-transform spectrum non 
 *              nonredundantly
 *              IEEE Trans. Signal Process, 58, No. 1, 281-290
 *
 * Trace header fields accessed: ns, dt, trid, ntr
 * Trace header fields modified: tracl, tracr, d1, f2, d2, trid, ntr
 */
/**************** end self doc ***********************************/

/* Prototype of function used internally */
static void st(int len, int lo, int hi, float *data, float **result);

/* constants used internally */
#define FRAC0   0.0	/* Ratio of default fmin to Nyquist */
#define FRAC1   1.0	/* Ratio of default fmax to Nyquist */
#define LOOKFAC 2	/* Look ahead factor for npfao	*/
#define PFA_MAX 720720  /* Largest allowed nfft	   */

/* global SEGY declaration */
segy tr;

int
main(int argc, char **argv)
{
	register float *rt;	/* real trace				*/

	float **vgr;		/* S-transform data  array		*/
	float dt;		/* sample spacing			*/
	float nyq;		/* nyquist frequency			*/
	int nt;			/* number of points on input trace	*/
	int nfft;		/* number of points for fft trace	*/
	int nf;			/* number of frequencies (incl Nyq)	*/
	cwp_Bool seismic;	/* is this seismic data?		*/

	int ntr;		/* number of traces 			*/
	int nflow;		/* number of filters 			*/
	int nfupr;		/* number of filters 			*/
	int tracr=0;		/* tracr counter			*/
	int verbose=0;		/* verbose flag				*/
        float d1;               /* output sample interval in Hz         */

	float fmin;		/* minimum frequency to window		*/
	float fmax;		/* maximum frequency to window		*/

	/****************************************************************/
        
	
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);


	/* Get info from first trace */ 
	if (!gettr(&tr))  err("can't get first trace");
	seismic = ISSEISMIC(tr.trid);

	ntr = tr.ntr;
		
	if (!seismic)
		warn("input is not seismic data, trid=%d", tr.trid);
	nt = tr.ns;
	if (!getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
	if (!dt) err("dt field is zero and not getparred");

	nyq = 0.5/dt;
	
	/******************/
	
	
	/* Get parameters */
	if (!getparint("verbose", &verbose))	verbose=0;
	if (!getparfloat("fmin", &fmin))	fmin = FRAC0*nyq;
	if (!getparfloat("fmax", &fmax))	fmax = FRAC1*nyq;

	/* Set up FFT parameters */
	nfft = npfaro(nt, LOOKFAC * nt);
	if (nfft >= SU_NFLTS || nfft >= PFA_MAX)
		err("Padded nt=%d -- too big", nfft);
	
	d1 = 1.0/(nfft*dt);
	
	nflow = (int)(fmin/d1 +1);
	nfupr = (int)(fmax/d1 +1 +0.5);
	nf = nfupr - nflow +1;
	

        checkpars();


	/* Allocate arrays */
	rt = ealloc1float(nfft);
	vgr = ealloc2float(nt,nf);

	/* Main loop over traces */
	do {
		register int i,j;

		/* Load trace into rt (zero-padded) */
		memcpy((void *) rt, (const void *) tr.data, nt*FSIZE);
		memset((void *) (rt + nt),0, (nfft-nt)*FSIZE);

		/* Stockwell transform */
	
		st(nfft, nflow,nfupr,rt,vgr);
			
		tracr =0;
		for ( i=0; i< nf; i++){
			for(j=0; j< nt; j++){
			tr.data[j] = vgr[i][j];
                       }
		       tr.d1 = dt;
		       tr.tracr = ++tracr;
		       tr.gx = tr.tracr;
		       tr.f2=fmin*1.0;
		       tr.d2=d1*1.0;
		       tr.trid = AMPLITUDE;
		       puttr(&tr);
		}
	} while (gettr(&tr));

	return(CWP_Exit());
}

static double gaussa(int n, int m);

static void st(int len, int lo, int hi, float *data, float **result)
{

    int i, k, n, l2;
    double s,ab1,ab2;
    static double *g,*p;
    complex *h, *G;
   
    
	
	/* Check for frequency defaults. */
	
	if (lo == 0 && hi == 0) {
		 hi = len / 2;
	}

	/** Alloc array ***/


        p   = ealloc1double(len);
        g   = ealloc1double(len);

        h = ealloc1complex(len);
        G = ealloc1complex(len);
	
	
	/* Convert the input to complex. Also compute the mean. */
	
	s = 0.0;
	
	memset((void *) (h),0, (len)*CSIZE);
	
	for (i = 0; i < len; i++) {
		h[i].r = data[i];
		s += data[i];
	}
	
	s /= len;



	/***** Complex to complex FFT ***/

	pfacc(1, len, h);

	
	/* Hilbert transform. The upper half-circle gets multiplied by
    	two, and the lower half-circle gets set to zero.  The real axis
    	is left alone. */
	
	l2 = (len + 1) / 2;
	for (i = 1; i < l2; i++) {
		h[i].r *= 2.;
		h[i].i *= 2.;
	}
	
	l2 = len / 2 + 1;
	for (i = l2; i < len; i++) {
		 h[i].r = 0.0;
		 h[i].i = 0.0;

        }	
		
	
	
	n = lo;
	
	/*Rows contain the inverse FFT of the spectrum
    multiplied with the FFT of scaled gaussians. */

	
	while (n <= hi) {

	for( i=0; i< len; i++){
		p[i] = 0.0;
		g[i] = 0.0;
		}

	/* Scale the FFT of the gaussian. Negative frequencies
        wrap around. */

	g[0] = gaussa(n, 0);
	
	l2 = len / 2 + 1;
	for (i = 1; i < l2; i++) {
	     g[i] = g[len - i] = gaussa(n, i);
	}
	
	for (i = 0; i < len; i++) {
		s = g[i];
		k = n + i;
	 	if (k >= len) k -= len;

		G[i].r = h[k].r * s;
		G[i].i = h[k].i * s;
	}
	
	/* Inverse FFT the result to get the next row. */
	
	      pfacc(-1,len,G);
	
	for (i = 0; i < len; i++) {
	
		ab1 = G[i].r / len;
		ab1 =ab1*ab1;
		ab2 = G[i].i / len;
		ab2 = ab2*ab2;
		p[i] += (sqrt)(ab1+ab2);

	}

	/* Go to the next row. */

	for (i = 0; i < len; i++) 
	result[n-lo][i] = (float)(p[i]);
	
	n++;
	
   }

}
/*  Fourier Transform of a Gaussian. */
 
static double gaussa(int n, int m)
{ 
    return exp(-2. * PI * PI * m * m / (n * n));
} 

