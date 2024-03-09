/* Copyright (c) Colorado School of Mines, 2022.*/
/* All rights reserved.                       */

/* SUGETQCSV: $Revision: 1.01 $ ; $Date: 2023/11/26 00:00:01 $		*/
 
#include "su.h"
#include "segy.h" 
#include <stdbool.h>
#include "qdefine.h"
#include "headcase.h"


/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SUGETQCSV - Write Trace Header Key Values and Sample Values to a Q-file    ",
"									     ",
"  sugetqcsv [parameters].                                                   ",
"									     ",
" Parameters:	         						     ",
"                                                                            ",
" qout=qhed.csv  Output csv Q-file name (default is qhed.csv)                ",
"									     ",
" keys=     Key name list. Must be trace header key names, or must be        ",
"           one of the special names below:                                  ",
"             mx = midpoint X coordinate (sx+gx)/2.0                         ",
"             my = midpoint Y coordinate (sy+gy)/2.0                         ",
"             mg = midpoint station (grnlof+gaps)/2.0                        ",
"     Note: Names sx,sy,gx,gy,mx,my are scaled by scalco before output and   ",
"           selev,gelev,sdepth,sdel,gdel,swdep,gwdep are scaled by scalel.   ",
"									     ",
" formxy=%.20g  The C format code for printing all values to the q-records.  ",
"              Note that this default format prints up to 20 digits but      ",
"              not trailing zeroes to the right of the decimal point.        ",
"              (This format code is the default because some coordinates     ",
"               scaled by scalco MIGHT need it for full representation).     ",
"                                                                            ",
"                                                                            ",
" fsamp=0   First sample value to output. The first sample is number 0.      ",
"                                                                            ",
" lsamp=-1  Last sample value to output.                                     ",
"           The default (-1) means do not output any sample values.          ",
"           If this value is greater than the number of samples in the       ",
"           input traces, the output is restricted to the number of samples. ",
"     Note: In accordance with q-file standards, sample values are given the ",
"           names numb0, numb1, numb2... (numb and sample number appended).  ",
"      ***  To output sample values, all input traces must have the          ",
"           same amount of samples as the first input trace.                 ",
"                                                                            ",
" formsv=%.10g  The C format code for printing sample values to q-records.   ",
"               Note that this default format prints up to 10 digits but     ",
"               not trailing zeroes to the right of the decimal point.       ",
"               (This format code is the default because sample values are   ",
"               single precision float and need it for full representation). ",
"                                                                            ",
"                                                                            ",
NULL};

/* Created: Nov  2023: Andre Latour                                          */ 
/**************** end self doc *******************************************/

segy tr;

double scaledfromhead(segy tr, int k);

/*----------------------------------------------------------------------*/

int main(int argc, char **argv) {

  cwp_String qout=NULL;   /* name for output q-file                */
  FILE *fpQ=NULL;         /* file pointer for output q-file        */

  int i = 0;                                                                         

  cwp_String *keyn=NULL;  
  int *kase=NULL;
  double *dvals=NULL;
  int numd = 0;

  cwp_String formxyt=NULL;
  cwp_String formxy=NULL;
  cwp_String formxylong=NULL;
  int lenformxy = 0;

  int fsamp = 0;                                                                         
  int lsamp = -1;                                                                         
  cwp_String formsvt=NULL;
  cwp_String formsv=NULL;
  cwp_String formsvlong=NULL;
  int lenformsv = 0;

/* hook up getpar */

  initargs(argc, argv);
  requestdoc(1);

  if(countparval("qout")>0) {
    getparstring("qout", &qout);
    fpQ = fopen(qout,"w");
  }
  else {
    fpQ = fopen("qhed.csv","w");
  }
  if (fpQ == NULL) err("error: output q-file did not open correctly.");




  numd = countparval("keys");          
  if(numd<1) err("**** Error: andre some keys is needed (default?) ");

  keyn  = ealloc1(numd,sizeof(cwp_String *));
  kase  = ealloc1int(numd);
  dvals = ealloc1double(numd);

  getparstringarray("keys",keyn);

/* Resolve a few things. ---------------------------------------------------- */

  for(i=0; i<numd; i++) {

    kase[i] = -1;
    if(strcmp(keyn[i],"sx")==0) kase[i] = 101;
    else if(strcmp(keyn[i],"sy")==0) kase[i] = 102;
    else if(strcmp(keyn[i],"gx")==0) kase[i] = 103;
    else if(strcmp(keyn[i],"gy")==0) kase[i] = 104;
    else if(strcmp(keyn[i],"mx")==0) kase[i] = 105;
    else if(strcmp(keyn[i],"my")==0) kase[i] = 106;
    else if(strcmp(keyn[i],"mg")==0) kase[i] = 107;
    else if(strcmp(keyn[i],"selev")==0) kase[i] = 111;
    else if(strcmp(keyn[i],"gelev")==0) kase[i] = 112;
    else if(strcmp(keyn[i],"sdepth")==0) kase[i] = 113;
    else if(strcmp(keyn[i],"sdel")==0) kase[i] = 114;
    else if(strcmp(keyn[i],"gdel")==0) kase[i] = 115;
    else if(strcmp(keyn[i],"swdep")==0) kase[i] = 116;
    else if(strcmp(keyn[i],"gwdep")==0) kase[i] = 117;
    else {
      kase[i] = GetCase(keyn[i]);
      if(kase[i]<1) err("**** Error: Specified dimension name %s is not recognized.",keyn[i]);
    }

  } /* end of  for(i=0; i<numd; i++) { */

  getparstring("formxy",&formxyt);
  if(formxyt==NULL) {
    lenformxy = 5;
    formxy = ealloc1(lenformxy,1);
    strcpy(formxy,"%.20g");
  }
  else {
    lenformxy = strlen(formxyt);
    formxy = ealloc1(lenformxy,1);
    strcpy(formxy,formxyt);
  }

  formxylong = ealloc1(1+lenformxy,1);
  strcpy(formxylong,",");
  strcpy(formxylong+1,formxy);

  getparstring("formsv",&formsvt);
  if(formsvt==NULL) {
    lenformsv = 5;
    formsv = ealloc1(lenformsv,1);
    strcpy(formsv,"%.10g");
  }
  else {
    lenformsv = strlen(formsvt);
    formsv = ealloc1(lenformsv,1);
    strcpy(formsv,formsvt);
  }

  formsvlong = ealloc1(1+lenformsv,1);
  strcpy(formsvlong,",");
  strcpy(formsvlong+1,formsv);

  if(!getparint("fsamp",&fsamp)) fsamp = 0;
  if(fsamp<0) err("**** Error: fsamp cannot be less than 0"); 
  
  if(!getparint("lsamp",&lsamp)) lsamp = -1;
  if(lsamp<-1) err("**** Error: lsamp cannot be less than -1"); 
  
/*--------------------------------------------------------------------------  */

  if (!gettr(&tr))  err("Error: cannot get first trace");

  if(lsamp>-1) {
    if(lsamp>tr.ns) lsamp = tr.ns - 1;
    if(fsamp>lsamp) err("**** Error: fsamp must be less than or equal to lsamp");
  }

/*--------------------------------------------------------------------------  */
 
  fprintf(fpQ,"C_SU_MATCH,%s\nC_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID",keyn[0]);
  for(i=0; i<numd; i++) fprintf(fpQ,",%s",formxy);
  for(i=fsamp; i<=lsamp; i++) fprintf(fpQ,",%s",formsv);
  fprintf(fpQ,"\nC_SU_NAMES\nC_SU_ID");
  for(i=0; i<numd; i++) fprintf(fpQ,",%s",keyn[i]);
  for(i=fsamp; i<=lsamp; i++) fprintf(fpQ,",numb%d",i); 
  fprintf(fpQ,"\n");

/* loop over traces ------------------------------------------------------   */ 

  do {

    for(i=0; i<numd; i++) dvals[i] = scaledfromhead(tr,kase[i]);

/* Write the q-records.                                                       */

    fprintf(fpQ,"Q");
    for(i=0; i<numd; i++) fprintf(fpQ,formxylong,dvals[i]);
    for(i=fsamp; i<=lsamp; i++) fprintf(fpQ,formsvlong,tr.data[i]);
    fprintf(fpQ,"\n");

  } while (gettr(&tr));

  return 0;

} /* end of main for sugetqcsv */

double scaledfromhead(segy tr, int k) {

  double dval;

  switch (k) {
      
    case 101:
      dval = (double)(tr.sx);
      if(tr.scalco > 1) dval *= tr.scalco;
      else if(tr.scalco < 0) dval /= -tr.scalco;
    break;
        
    case 102:
      dval = (double)(tr.sy);
      if(tr.scalco > 1) dval *= tr.scalco;
      else if(tr.scalco < 0) dval /= -tr.scalco;
    break;
        
    case 103:
      dval = (double)(tr.gx);
      if(tr.scalco > 1) dval *= tr.scalco;
      else if(tr.scalco < 0) dval /= -tr.scalco;
    break;
        
    case 104:
      dval = (double)(tr.gy);
      if(tr.scalco > 1) dval *= tr.scalco;
      else if(tr.scalco < 0) dval /= -tr.scalco;
    break;
        
    case 105:
      dval = 0.5 * ((double)(tr.sx) + (double)(tr.gx));
      if(tr.scalco > 1) dval *= tr.scalco;
      else if(tr.scalco < 0) dval /= -tr.scalco;
    break;
        
    case 106:
      dval = 0.5 * ((double)(tr.sy) + (double)(tr.gy));
      if(tr.scalco > 1) dval *= tr.scalco;
      else if(tr.scalco < 0) dval /= -tr.scalco;
    break;

    case 107:
      dval = 0.5 * ((double)(tr.grnlof) + (double)(tr.gaps));
    break;

/* ---- */

    case 111:
      dval = (double)(tr.selev);
      if(tr.scalel > 1) dval *= tr.scalel;
      else if(tr.scalel < 0) dval /= -tr.scalel;
    break;

    case 112:
      dval = (double)(tr.gelev);
      if(tr.scalel > 1) dval *= tr.scalel;
      else if(tr.scalel < 0) dval /= -tr.scalel;
    break;
   
    case 113:
      dval = (double)(tr.sdepth);
      if(tr.scalel > 1) dval *= tr.scalel;
      else if(tr.scalel < 0) dval /= -tr.scalel;
    break;

    case 114:
      dval = (double)(tr.sdel);
      if(tr.scalel > 1) dval *= tr.scalel;
      else if(tr.scalel < 0) dval /= -tr.scalel;
    break;

    case 115:
      dval = (double)(tr.gdel);
      if(tr.scalel > 1) dval *= tr.scalel;
      else if(tr.scalel < 0) dval /= -tr.scalel;
    break;

    case 116:
      dval = (double)(tr.swdep);
      if(tr.scalel > 1) dval *= tr.scalel;
      else if(tr.scalel < 0) dval /= -tr.scalel;
    break;

    case 117:
      dval = (double)(tr.gwdep);
      if(tr.scalel > 1) dval *= tr.scalel;
      else if(tr.scalel < 0) dval /= -tr.scalel;
    break;

    default:
      dval = fromhead(tr, k);
    break;

  } /* end of switch */

  return (dval);
}
