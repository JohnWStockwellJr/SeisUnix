/* Copyright (c) Colorado School of Mines, 2015.*/
/* All rights reserved.			*/

/* SUOP2: $Revision: 1.25 $ ; $Date: 2015/05/22 18:31:12 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 									",
" SUOP2 - do a binary operation on two data sets			",
" 									",
" suop2 data1 data2 op=diff [trid=111] >stdout				",
" 									",
" Required parameters:							",
" 	none								",
" 									",
" Optional parameter:							",
" 	op=diff		difference of two panels of seismic data	",
" 			=sum  sum of two panels of seismic data		",
" 			=prod product of two panels of seismic data	",
" 			=quo quotient of two panels of seismic data	",
" 			=ptdiff differences of a panel and single trace	",
" 			=ptsum sum of a panel and single trace		",
" 			=ptprod product of a panel and single trace	",
" 			=ptquo quotient of a panel and single trace	",
" 			=zipper do \"zipper\" merge of two panels	",
"			=zippol convert polar to rectangular and then zip",
" 									",
"  trid=FUNPACKNYQ	output trace identification code. (This option  ",
" 			is active only for op=zipper)			",
"			For SU version 39-43 FUNPACNYQ=111		",
" 			(See: sukeyword trid     for current value)	",
" 									",
" 									",
" Note1: Output = data1 \"op\" data2 with the header of data1		",
" 									",
" Note2: For convenience and backward compatibility, this		",
" 	program may be called without an op code as:			",
" 									",
" For:  panel \"op\" panel  operations: 				",
" 	susum  file1 file2 == suop2 file1 file2 op=sum			",
" 	sudiff file1 file2 == suop2 file1 file2 op=diff			",
" 	suprod file1 file2 == suop2 file1 file2 op=prod			",
" 	suquo  file1 file2 == suop2 file1 file2 op=quo			",
" 									",
" For:  panel \"op\" trace  operations: 				",
" 	suptsum  file1 file2 == suop2 file1 file2 op=ptsum		",
" 	suptdiff file1 file2 == suop2 file1 file2 op=ptdiff		",
" 	suptprod file1 file2 == suop2 file1 file2 op=ptprod		",
" 	suptquo  file1 file2 == suop2 file1 file2 op=ptquo		",
" 									",
" Note3: If an explicit op code is used it must FOLLOW the		",
"	filenames.							",
" 									",
" Note4: With op=quo and op=ptquo, divide by 0 is trapped and 0 is returned.",
" 									",
" Note5: Weighted operations can be specified by setting weighting	",
"	coefficients for the two datasets:				",
"	w1=1.0								",
"	w2=1.0								",
" 									",
" Note6: With op=zipper, it is possible to set the output trace id code ",
" 		(See: sukeyword trid)					",
"  This option processes the traces from two files interleaving its samples.",
"  Both files must have the same trace length and must not longer than	",
"  SU_NFLTS/2  (as in SU 39-42  SU_NFLTS=32768).			",
"			  						",
"  Being \"tr1\" a trace of data1 and \"tr2\" the corresponding trace of",
"  data2, The merged trace will be :					",
"			  						",
"  tr[2*i]= tr1[i]							",
"  tr[2*i+1] = tr2[i]							",
" 									",
"  The default value of output tr.trid is that used by sufft and suifft,",
"  which is the trace id reserved for the complex traces obtained through",
"  the application of sufft. See also, suamp.				",
" 									",
" Note 7: op=zippol is like op=zipper, but the input samples are polar	",
"	(amplitude and phase) and are converted to cartesian (real, imag)",
"	before interleaving them.					",
" 									",
" For operations on non-SU binary files  use:farith 			",
NULL};

/* Credits:
 *	SEP: Shuki Ronen
 *	CWP: Jack K. Cohen
 *	CWP: John Stockwell, 1995, added panel op trace options.
 *	: Fernando M. Roxo da Motta <petro@roxo.org> - added zipper op
 *
 * Notes:
 *	If efficiency becomes important consider inverting main loop
 *	and repeating operation code within the branches of the switch.
 */
/**************** end self doc ***********************************/

#define	ADD	1
#define	SUB	2
#define	MUL	3
#define	DIV	4
#define	ZIPPER	5
#define	ZIPPOL	6
#define	TPANOP	6
#define	PTADD	7
#define	PTSUB	8
#define	PTMUL	9
#define	PTDIV	10

segy intrace1, intrace2;

int
main(int argc, char **argv)
{
	FILE *fp1=NULL;		/* file pointer for first file	*/
	FILE *fp2=NULL;		/* file pointer for second file	*/
	cwp_String op="diff";	/* operation: add, sub, mul, div*/
	int iop=SUB;		/* integer abbrev. for op in switch	*/
	int nt=0;		/* number of sample points on trace	*/
	int nbytes=0;		/* number of bytes on traces		*/
	int itr=0;		/* number of trace being processed	*/
	int trid=0;  		 /* number of trace being processed	*/

	float w1=0.0;	/* weighting function for first file*/
	float w2=0.0;	/* weighting function for second file  */

	/* Initialize */
	initargs(argc, argv);
	requestdoc(2); /* two file args required */


	/* Open two files given as arguments for reading */
	fp1 = efopen(argv[1], "r");
	fp2 = efopen(argv[2], "r");


	getparfloat( "w1" ,&w1 );
	getparfloat( "w2" ,&w2 );

	/* Get operation , recall iop initialized to the default SUB*/
	getparstring("op", &op);
	if	(STREQ(op, "sum"))	iop = ADD;
	else if (STREQ(op, "prod"))	iop = MUL;
	else if (STREQ(op, "quo"))	iop = DIV;
	else if (STREQ(op, "ptdiff"))	iop = PTSUB;
	else if (STREQ(op, "ptsum"))	iop = PTADD;
	else if (STREQ(op, "ptprod"))	iop = PTMUL;
	else if (STREQ(op, "ptquo"))	iop = PTDIV;
	else if (STREQ(op, "zipper"))	iop = ZIPPER;
	else if (STREQ(op, "zippol"))	iop = ZIPPOL;
	else if (!STREQ(op, "diff"))
		err("unknown operation=\"%s\", see self-doc", op);

	/* getpar trid if set */
	if ( iop == ZIPPER || iop == ZIPPOL ) {
		/* Get optional parameter trid */
		if (!getparint( "trid" ,&trid ))	trid = FUNPACKNYQ;

	}

	checkpars();

	/* panel operations */
	if (iop <= TPANOP) { /* do panel op panel operations */
		/* Loop over the traces */
		while (fgettr(fp1, &intrace1) &&
					(nbytes = fgettr(fp2, &intrace2))) {

			/* verify that ns is the same in both files */
			if ((nt = intrace1.ns) != intrace2.ns) {
				warn("trace %d:", itr);
				err("%s and %s have different ns (%d vs %d)",
				argv[1], argv[2], intrace1.ns, intrace2.ns);
			}
			if( intrace1.trid == FUNPACKNYQ && intrace2.trid == FUNPACKNYQ ){
			   warn( "performing complex arithmetic" );
			}

			/* apply scaling factors */
			if( w1 != 0.0 && w2 != 0.0 ){
				int i;
				for (i = 0; i < nt; ++i){
					intrace1.data[i] *= w1;
					intrace2.data[i] *= w2;
				}
	 		}
					
			/* Do the desired binary operation */
			switch(iop) { register int i;
			case SUB:
				for (i = 0; i < nt; ++i)
					intrace1.data[i] -= intrace2.data[i];
			break;
			case ADD:
				for (i = 0; i < nt; ++i)
					intrace1.data[i] += intrace2.data[i];
			break;
			case MUL:
				if( intrace1.trid == FUNPACKNYQ && intrace2.trid == FUNPACKNYQ  ){
				   float a,b ,c ,d;
   				   for (i = 2; i < nt-2; i+=2 ) {
					 a = intrace1.data[i];
					 b = intrace1.data[i+1];
					 c = intrace2.data[i];
					 d = intrace2.data[i+1];
					 intrace1.data[i] = a*c - b*d; 
					 intrace1.data[i+1] = a*d + b*c;
				   
				   }
				}else{
				   for (i = 0; i < nt; ++i){
					   intrace1.data[i] *= intrace2.data[i];
				   }
				}
			break;
			case DIV:

				if( intrace1.trid == FUNPACKNYQ && intrace2.trid == FUNPACKNYQ ){
				   float a,b ,c ,d;
   				   for (i = 2; i < nt; i+=2 ) {
				      a = intrace1.data[i];
				      b = intrace1.data[i+1];
				      c = intrace2.data[i];
				      d = intrace2.data[i+1];
				      if( c == 0.0 && d == 0.0 ){
					 intrace1.data[i] = 0.0;
					 intrace1.data[i+1] = 0.0;
				      }else if( c == 0.0 ){
					 intrace1.data[i] =  - b/d;
				      }else if( d == 0.0 ){
					 intrace1.data[i+1] = b/c;
				      }else{
					 intrace1.data[i] = a/c - b/d; 
					 intrace1.data[i+1] = a/d + b/c;
				   
				      }
				   }
				}else{
				   for (i = 0; i < nt; ++i) {
					   float denom = intrace2.data[i];
					   if (!denom) intrace1.data[i] = 0.0;
					   else	intrace1.data[i] /= denom;
				   }
				}
			break;
			case ZIPPER:
				for (i = nt-1; i >= 0; --i) {
					intrace1.data[2*i] = intrace1.data[i];
					intrace1.data[2*i+1] = intrace2.data[i];
				}

				/* update output trace header */
				intrace1.trid = trid;
				intrace1.ns= 2*nt;
			break;
			case ZIPPOL:
				for (i = nt-1; i >= 0; --i) {
					intrace1.data[2*i+1] = intrace1.data[i]
							*sin(intrace2.data[i]);
					intrace1.data[2*i] = intrace1.data[i]
							*cos(intrace2.data[i]);
				}
				/* update output trace header */
				intrace1.trid = trid;
				intrace1.ns= 2*nt;
			break;
			default:  /* defensive programming */
				err("mysterious operation=\"%s\"", op);
			}

			puttr(&intrace1);
			++itr;
		}

		/* See if both files exhausted; note if fd1 exhausted, then */
		/* we don't do an fgettr on fd2 on the final pass above */
		if (!nbytes) {
			warn("%s still had traces when %s was exhausted",
						argv[1], argv[2]);
			warn("processed %d pairs of traces before EOF", itr);
		} else if (fgettr(fp2, &intrace2)) {
			warn("%s still had traces when %s was exhausted",
						argv[2], argv[1]);
			warn("processed %d pairs of traces before EOF", itr);
		}
	} else { /* do panel op trace operations */

		/* get single seismic trace from file2 */
		fgettr(fp2,&intrace2);

		/* apply scaling factor for trace 2 */
		if(w2!=0.0){
 			int i;
	 		for (i = 0; i < nt; ++i){
				intrace2.data[i] *= w2;
	 		}
		}

		/* Loop over the traces */
		while (fgettr(fp1, &intrace1)) {

			/* apply scaling factor for trace 1 */
	 		if(w1!=0.0) {
				int i;
	  			for (i = 0; i < nt; ++i) {
					intrace1.data[i] *= w1;
	  			}
			}

			/* check that the number of samples is same on both */
			if ((nt = intrace1.ns) != intrace2.ns) {
				warn("trace %d:", itr);
				err("%s and %s have different ns (%d vs %d)",
				argv[1], argv[2], intrace1.ns, intrace2.ns);
			}

			/* Do the desired binary operation */
			switch(iop) { register int i;
			case PTSUB:
				for (i = 0; i < nt; ++i)
					intrace1.data[i] -= intrace2.data[i];
			break;
			case PTADD:
				for (i = 0; i < nt; ++i)
					intrace1.data[i] += intrace2.data[i];
			break;
			case PTMUL:
				if( intrace1.trid == FUNPACKNYQ && intrace2.trid == FUNPACKNYQ ){
				   float a,b ,c ,d;
   				   for (i = 0; i < nt; i+=2 ) {
					 a = intrace1.data[i];
					 b = intrace1.data[i+1];
					 c = intrace2.data[i];
					 d = intrace2.data[i+1];
					 intrace1.data[i] = a*c - b*d; 
					 intrace1.data[i+1] = a*d + b*c;
				   
				   }
				}else{
				   for (i = 0; i < nt; ++i)
					   intrace1.data[i] *= intrace2.data[i];
				}
			break;
			case PTDIV:
				if( intrace1.trid == FUNPACKNYQ  && intrace2.trid == FUNPACKNYQ ){
				   float a,b ,c ,d;
   				   for (i = 0; i < nt; i+=2 ) {
				      a = intrace1.data[i];
				      b = intrace1.data[i+1];
				      c = intrace2.data[i];
				      d = intrace2.data[i+1];
				      if( c == 0.0 && d == 0.0 ){
					 intrace1.data[i] = 0.0;
					 intrace1.data[i+1] = 0.0;
				      }else if( c == 0.0 ){
					 intrace1.data[i] =  - b/d;
				      }else if( d == 0.0 ){
					 intrace1.data[i+1] = b/c;
				      }else{
					 intrace1.data[i] = a/c - b/d; 
					 intrace1.data[i+1] = a/d + b/c;
				   
				      }
				   }
				}else{
				   for (i = 0; i < nt; ++i) {
					   float denom = intrace2.data[i];
					   if (!denom) intrace1.data[i] = 0.0;
					   else	intrace1.data[i] /= denom;
				   }
				}
			break;
			default:  /* defensive programming */
				err("mysterious operation=\"%s\"", op);
			}

			puttr(&intrace1);
			++itr;
		}

	}


	return(CWP_Exit());
}

