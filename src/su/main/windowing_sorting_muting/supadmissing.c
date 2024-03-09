/* Copyright (c) Colorado School of Mines, 2023.*/
/* All rights reserved.			*/

/* SUPADMISSING: $Revision: 1.0 $ ; $Date: 2024/02/15 11:00:01 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "headcase.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                       ",
" SUPADMISSING - Add Dead Traces For Missing Key Values.                ",
"                                                                       ",
" supadmissing <stdin >stdout                                           ",
"                                                                       ",
" keyloc=cdp  Input traces must be in increasing order of this key.     ",
" 									",
" incval=1.0  Increment. Must be greater than 0.0                       ",
"             Input traces should have keyloc values at whole number    ",
"             multiples of incval added to minval (minval + incval*N).  ",
"             So if minval=2 and incval=5 then 2,7,12,17...are expected.",
"             Any missing values will cause a padded trace to be output.",
"             The padded trace has keyloc set to the missing value.     ",
"             Padded traces have all amplitudes equal to 0.0            ",
"             Padded traces have nhs=0 set but all other key values     ",
"             are copied from the next higher input trace (or last).    ",
"        ***  This means the padded header values are NOT interpolated. ",
"                                                                       ",
" minval=     Minimum keyloc value to output. Input traces with keyloc  ",
"             less than this are not output. If first input trace is    ",
"             greater than this value, padding occurs up to this value. ",
"             Default is keyloc value from first input trace.           ",
"        ***  This default is often NOT what you want to do.            ",
" 									",
" maxval=     Maximum keyloc value to output. Input traces with keyloc  ",
"             greater than this are not output. If last input trace is  ",
"             less than this value, padding occurs up to this value.    ",
"             Default is keyloc value from last input trace.            ",
"        ***  This default is often NOT what you want to do.            ",
" 									",
" NOTE: Input trace range can be completely before or completely after  ",
"       the minval,maxval range in which case no input trace is output. ",
"       But padded traces are still output from minval to maxval using  ",
"       the header values from the first or last input trace.           ",
" 									",
" NOTE: This program does not honour scalco and scalel. So if you use   ",
"       any key associated with them as keyloc, then specify incval,    ",
"       minval, maxval based on the actual values in the trace headers. ",
" 									",
NULL};

/* Author:
 *	Andre Latour. Feb 2024   
 *	1. This program adds output traces for missing key values.
 *	   This is primarily intended to add missing cdps in a stack.     
 *	   This may need to be done before using certain programs.       
 *	   Notably, some old 2D migration programs assume that traces   
 *	   and velocity functions have a one-to-one sequential input
 *         relationship (and they do not check the cdp key at all).  
 */
/**************** end self doc ***********************************/

segy tr;
segy trz;

int main(int argc, char **argv) {

  cwp_String keyloc = NULL;
  int kaseloc=0;

  int i=0;
  double incval=1.0;
  int minhave=0;
  int maxhave=0;
  double minval=0.0;
  double maxval=0.0;
  double dval=0.0;
  double dpre=0.0;

/* hook up getpar */

  initargs(argc, argv);
  requestdoc(1);

/* Get key case numbers from GetCase.                                        */

  if(countparval("keyloc")>1) err("**** Error: keyloc cannot be multiple keys.");

  if(countparval("keyloc")>0) {
    getparstring("keyloc", &keyloc);
  }    
  else {
    keyloc = ealloc1(3,1);
    strcpy(keyloc,"cdp");
  }    

  kaseloc = GetCase(keyloc);
  if(kaseloc<1) err("**** Error: Specified keyloc name %s is not recognized.",keyloc);

  if (!getpardouble("incval",&incval)) incval = 1.0;
  if(incval<=0.0) err("**** Error: incval must be greater than 0.0 ");

  minhave = 0;
  if(countparval("minval")>0) {
    getpardouble("minval",&minval);
    minhave = 1;
  }
  maxhave = 0;
  if(countparval("maxval")>0) {
    getpardouble("maxval",&maxval);
    maxhave = 1;
  }

  checkpars();

  if (!gettr(&tr))  err("Error: cannot get first trace");

/* Copy header from first trace. Set samples to zero. Set stack fold to 0.    */

  memcpy(&trz,&tr,HDRBYTES);
  for (i=0; i<trz.ns; i++) trz.data[i] = 0.0;
  trz.nhs = 0; 

  dval = fromhead(tr,kaseloc);
  if(minhave==0) minval = dval; 
  dpre = dval; /* dpre is used to check input is increasing order             */

/* Set the maxval value so that pad-stopping is only triggered by last input. */

  if(maxhave==0) maxval = 1.0e30; 

  if (minval>maxval)  err("Error: minval cannot be greater than maxval");

/* To get started: bypass input traces less than minval parameter or          */
/*                 else output pad traces from first trace TO minval.         */

  if(minval>dval) {
    while(minval>dval) {
      if (!gettr(&tr)) {                                                            

/* False means all input traces are less than minval. For this specific       */
/* case tr gets output once by main loop (the lines after 162) so reset tr    */
/* to contain padded values. This is all twiddly, so extensive testing done.  */

        tr.nhs = 0;                 
        tohead(&tr,kaseloc,minval); 
        for (i=0; i<tr.ns; i++) tr.data[i] = 0.0; 
        dval = minval;                           
        break;                                  
      }
      dval = fromhead(tr,kaseloc);
      if(dval<=dpre) { 
        err("Error: keyloc=%s is not increasing. Next trace=%f  Previous=%f",keyloc,dpre,dval);
      }
      dpre = dval;
    }
    memcpy(&trz,&tr,HDRBYTES);
    trz.nhs = 0; 
  }
  else if(minval<dval) {
    while(minval<dval) {
      tohead(&trz,kaseloc,minval);
      puttr(&trz);
      minval+=incval;
      if(minval>maxval) break;
    }
  }

/* Next, main looping over the input traces. At this point dval contains the  */
/* value from the first INPUT trace that is going to be output.               */
/* But it still may be greater than minval here because the input traces may  */
/* have a gap surrounding the user specified minval value.                    */

  while(dval<=maxval) {

    if(minval<dval) {
      memcpy(&trz,&tr,HDRBYTES);
      trz.nhs = 0; 
      while(minval<dval) {
        tohead(&trz,kaseloc,minval);
        puttr(&trz);
        minval+=incval;
      }
    }

    puttr(&tr); /* output this input trace                                    */

    if (gettr(&tr)) { /* try to input another trace                           */
      dval = fromhead(tr,kaseloc);
      minval+=incval;
      if(dval<=dpre) {
        err("Error: keyloc=%s is not increasing. Next trace=%f  Previous=%f",keyloc,dpre,dval);
      }
      dpre = dval;
    }
    else {
      if(maxhave==1) {
        memcpy(&trz,&tr,HDRBYTES); /* note this assumes tr is still good copy */
        trz.nhs = 0;               /* of last trace after gettr==false return */
        while(dval<maxval) {
          dval+=incval;
          tohead(&trz,kaseloc,dval);
          puttr(&trz);
        }
      }
      break;
    }

  }

  return(CWP_Exit());
}

