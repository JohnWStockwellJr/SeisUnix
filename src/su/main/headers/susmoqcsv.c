/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUBINQCSV: $Revision: 1.00 $ ; $Date: 2024/02/25 00:00:01 $		*/
 
#include "su.h"
#include "segy.h" 
#include "qdefine.h"
#include "gridread.h"
#include "gridxy.h"
#include "linterpd.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SUSMOQCSV - Input Q-records, alter and spatially smooth, output Q-records. ",
"									     ",
"  susmoqcsv [parameters].   (No traces in or out).                          ",
"									     ",
" Parameters:	         						     ",
"                                                                            ",
" qin=       Input file containing q-records. Must be specified.             ",
"            For 3D there must be the same number of q-records in the        ", 
"            qin file as cdps (cells) defined by the 3D Grid.                ", 
"									     ",
" qout=      Output file for q-records. Must be specified.                   ",
"            Output records will be in cdp order.                            ",
"									     ",
" rfile=     If set, read a K-record containing 3D grid definition.          ", 
"            Assumed to be 2D if no file is specified here.                  ", 
"									     ",
" check=0   Do not print checking details. The options of this program       ",
"           cannot repair a corrupt input grid. Grid checking is always      ",
"           done on input grid definition and will error-halt if corrupt.    ",
"      =1   print checking details.                                          ",
"									     ",
" inloc=1   Location options. The option here determines which input names   ",
"           contain the locations of all other input values. The option here ",
"           must find the corresponding names within the non-tuple names in  ",
"           the input q-file.                                                ",
"      =1   cdp (default). Option 1 is required if no grid is input.         ",
"      =2   cdpt (grid transposed cdp numbers).                              ",
"      =3   igi,igc (grid index numbers).                                    ",
"      =4   gx,gy (grid coordinates).  Note option 6.                        ",
"      =5   sx,sy (world coordinates). Note option 7.                        ",
"      =6   gx,gy (world coordinates).                                       ",
"      =7   sx,sy (grid coordinates).                                        ",
"     Note: The q-records in qin file do not have to be in order but there   ",
"           must be exactly one record for each cdp defined by the 3D Grid.  ",
"									     ",
" keyloc=   Location key. Can only be specified if inloc=1 and 2D (no rfile).",
"           Default is cdp. Name must be in the input Q-file.                ",
"           This allows smoothing based on values such as point numbers,     ",
"           station numbers, reciever numbers, or shot numbers.              ",
"           Note that q-records are internally sorted by the values of       ",
"           this name. Only cdps are required to be contiguous with an       ",
"           increment of 1 (no gaps). For other names navrg, nsmth, nback    ",
"           refer to the sorted record order, not to their actual values.    ",
"									     ",
" extrapt=0        do not extrapolate at ends in independent values.         ",
"                  (values outside range of functions are held constant).    ",
"               =1 extrapolate outside range of functions                    ",
"               =2 extrapolate when less than range of functions             ",
"               =3 extrapolate when greater than range of functions          ",
"                                                                            ",
" outind=0   Output independent dimension range. One, two, or three numbers. ",
"            No default, but outlist can be specified instead.               ",
"            Cannot be specified if input has no tuples.                     ",
"       =f,l,i start at f, increment by i, end at last value <= l.           ",
"       =f,l   start at f, increment by 1000, end at last value <= l.        ",
"       =f     start at f, increment by 1000, end at last value <= 3000.     ",
"            The i value cannot be 0, but can be negative, which requires    ",
"            an l value less than or equal to f.                             ",
"									     ",
" outlist=   Output independent dimension value list. Cannot be specified    ",
"            if outind is used. Cannot be specified if input has no tuples.  ",
"									     ",
" ename=     Names to eliminate. These names will be ignored in the input    ",
"            q-records and therefore not be in the output q-records (note    ",
"            that standard 2d and 3d names will still be output even if you  ",
"            specify them here). A primary purpose of this option is to get  ",
"            rid of names in input tuples in order to reduce the complexity  ",
"            of output q-records (in case certain programs cannot handle it).",
"            But you can also get rid of non-tuple names here.               ",
"									     ",
" sname=     Names to smooth. By default all names with 1 element per record ",
"            are NOT averaged and smoothed. Listing names here overrides     ",   
"            that and averages and smooths their values.                     ",
"                                                                            ",
" cname=     Names to just copy. By default all names with multiple values   ",   
"            per record (tuples) are averaged and smoothed. Listing them here",
"            overrides that and does NOT average and smooth their values.    ",   
"                                                                            ",
"                                                                            ",
" navrg=10       CDP Averaging Operator Size.                                ",
"                The averaging is symmetrical and extends navrg on each side ",
"                of the location being averaged.                             ",
"                For 2D surveys, this is number of cin Q-records to average  ",
"                symmetrically with its neighbours (navrg on each side).     ",
"                For 3D surveys, two values can be listed. The first value   ",
"                is for the inline direction and the second is for the       ",
"                crossline direction. If no second value is specified, the   ",
"                same value is used for both directions.                     ",
"            =0  means no averaging is done.                                 ",
"            =-1 means average the values of all locations in this direction.",
"                You also MUST set nsmth to 0 for this direction.            ",
"                Typically this option is used for 3D surveys where the      ",
"                number of inline or crossline cells is less than 5 or one   ",
"                of the directions has little change (or you want to produce ",
"                a constant value in that direction).                        ",
"           ***  You must use 0 or -1 if number of locations is less than 5. ",
"									     ",
" nsmth=navrg    CDP Smoothing Operator Size.                                ",
"                This smoothing uses a technique which leaves the values     ",
"                approximately centred with their initial values.            ",
"                The smoothing affect is symmetrical and extends 2*nsmth on  ",
"                each side of the location being smoothed.                   ",
"                For 3D surveys, two values can be listed. The first value   ",
"                is for the inline direction and the second is for the       ",
"                crossline direction. If no second value is specified, the   ",
"                same value is used for both directions.                     ",
"             =0 means no smoothing is done.                                 ",
"           ***  You must use 0 if number of locations is <5 or navrg=-1     ",
"									     ",
" nback=3        Linear Extrapolation Location.                              ",
"                For 2D surveys, this is the number of cin Q-records inwards ",
"                from the ends of the survey used to linearly extrapolate    ",
"                some padding outside the ends of the survey. These          ",
"                extrapolated values affect the averaging and smoothing of   ",
"                the navrg and nsmth parameters (at both ends of the survey).",
"                This value determines which cin records inside the survey   ",
"                are used to determine extrapolation points.                 ",
"                Must be less than the number of Q-records in the cin file.  ",
"            =1  means use the first and first+1 records to extrapolate into ",
"                the padding at the beginning of the survey -and- use the    ",
"                last and last-1 records to extrapolate into the end padding.",
"            =n  means use the first and first+n records to extrapolate into ",
"                the padding at the beginning of the survey -and- use the    ",
"                last and last-n records to extrapolate into the end padding.",
"            =0  means duplicate the values of first record into the padding ",
"                at the beginning of the survey -and- duplicate the values   ",
"                of the last record into the padding at the end of survey.   ",
"           ***  For 3D surveys, two values can be listed. The first value   ",
"                indicates how many cells inwards in the inline direction    ",
"                and the second value indicates how many cells inwards in    ",
"                the crossline direction. If no second value is specified,   ",
"                the same value is used for both directions.                 ",
"                Must be less than the number of 3D Grid cells defined       ",
"                in the same direction.                                      ",
"         Note:  The AMOUNT of added padding is NOT set by this parameter,   ",
"                it is actually set as needed for navrg and nsmth parameters.",
"									     ",
" cdpt=0      Output in transposed cdp order.                                ",
"          =0 No. Output in cdp order.                                       ",
"          =1 Yes. Output in transposed cdp order. Can only specify for 3D.  ",
"             This has no effect on values within the q-records.             ",
"             This order may be convenient (notably 2D by 2D migrations).    ",
"									     ",
" formtv=%.2f The C format code for printing all values in the q-records     ",
"             except the standard names cdp,cdpt,igi,igc.                    ",
"             These names are specially handled within this program.         ",
"             For 3D, they are generated for output using grid definition.   ",
"             For 2D, igi,igc are not output.                                ",
"             The default is 2 digits below the decimal point.               ",
"									     ",
" formxy=%.2f The C format code for printing gx,gy,sx,sy in the q-records.   ",
"             These names are specially handled within this program.         ",
"             For 3D, they are generated for output using grid definition.   ",
"             For 2D, they are not output.                                   ",
"             The default is 2 digits below the decimal point.               ",
"									     ",
"   ------------------------------------------------------------------       ",
"   ------------------------------------------------------------------       ",
"									     ",
"  Notes:                                                       	     ",
"									     ",
"   All input values are stored in double-precision (8 byte float).          ",
"									     ",
"   This program does not care if the values are seconds or milliseconds     ",
"   or feet or metres. To this program they are all just numbers.            ",
"									     ",
NULL};

/*   Created: Feb 2024: Andre Latour                                         */ 
/**************** end self doc *******************************************/

void runsmo (double **dalla, int msize, int klast, int nback, int nmax, double *dtog, int *ierr) ;

void runav (double **dall, int ncdp, int klast, int nmin, int nmax, double *dtog, int *ierr) ;

segy tr;

struct QInfo *RecInfo; /* Storage for all function location value pointers */

int compSort (const void * q1, const void * q2) ; /* comparison function for qsort  */

int main(int argc, char **argv) {

  int ncdp = 0;		/* number of cdps specified */
  int icdp = 0;		/* index into arrays dimensioned by ncdp */
  int icdpt = 0;        /* 90 degree grid cdp number             */
  int jcdp = 0;         /* index into arrays dimensioned by ncdp */

  int is3d = 1;         /* flag for 3d */
  int inloc = 1;        /* flag for 3d input location names        */
  int jnloc1 = -1;      
  int jnloc2 = -1;      

  int mgtextr = 0;      /* for independent dimension extrapolation */
  int icheck = 0;       /* a command line parameter value          */

  int ifixd = 0;          /* flag for all tuples same size or vary   */
  int iztuple = 0;        /* element number where first tuple exists */
  int ktuple = 0;         /* type of tuples (2=pairs, 3=triplets)    */

  cwp_String Rname=NULL;  /* text file name for input grid        */
  FILE *fpR=NULL;         /* file pointer for Rname input file    */

  cwp_String Pname=NULL;  /* text file name for Q input file      */
  FILE *fpP=NULL;         /* file pointer for Q input file        */
  cwp_String Qname=NULL;  /* text file name for output Q values   */
  FILE *fpQ=NULL;         /* file pointer for Q output file       */

  double gvals[99];       /* more than enough to contain grid definition */

  cwp_String *pname = NULL;                                         
  cwp_String *ndims = NULL;                                                 
  cwp_String *qname = NULL;                                              
  cwp_String *ename = NULL;                                              
  cwp_String *sname = NULL;                                              
  cwp_String *cname = NULL;                                              
  cwp_String *qform = NULL;   
  int numpname = 0;
  int numsname = 0;
  int numcname = 0;
  int numstandard = 0;

  double *outind = NULL;  /* this set of variables is related to   */
  double *pindepa = NULL; /* input of independent dimension values */ 
  int numdind = 0;
  double **dalla = NULL;
  double *dinput = NULL;
  double *dind = NULL;
  double *dswap = NULL;
  double *dtog  = NULL;
	
  cwp_String formxyt=NULL; /* this set of variables is related to output format */
  cwp_String formxy=NULL;
  cwp_String formxylong=NULL;
  int lenformxy = 0;
  cwp_String formtvt=NULL;
  cwp_String formtv=NULL;
  cwp_String formtvlong=NULL;
  int lenformtv = 0;

  int i = 0;              /* this set of variables are either trivial or so   */
  int j = 0;              /* entangled that no short comment up here will help*/
  int k = 0;  
  int n = 0;  
  int ierr = 0;  
  int isize = 0;  
  int errwarn = 0; 
  int *kindx = NULL;
  int *ksize = NULL;
  int *ksmth = NULL;
  int klast = 0;
  int cdpt=0;   
  double xw = 0.;
  double yw = 0.;
  double xg = 0.;
  double yg = 0.;
 
  int nback[2];
  int navrg[2];
  int nsmth[2];
  int lastab=0;
  int lastac=0;

/* hook up getpar */
  initargs(argc, argv);
  requestdoc(1);

  if(isatty(STDIN_FILENO)!=1 || isatty(STDOUT_FILENO)!=1)  
    err("**** Error: this program does not input or output traces.");
          
/* ------------------------------------------------------------------- */

  i = 0;
  gridcommand(&i);
  if(i!=0) err("error: grid definition parameters on command line (input grid only via K-file).");

  getparstring("rfile", &Rname);
  if(Rname == NULL) is3d = 0;

  if (!getparint("extrapt", &mgtextr)) mgtextr = 0;
  if(mgtextr<0 || mgtextr>3) err ("error: extrapt= option not in range ");

  if (!getparint("inloc", &inloc)) inloc = 1;
  if(inloc<1 || inloc>7) err ("error: inloc= option not in range ");
  if(is3d==0 && inloc>1) err("**** Error: cannot specify inloc>1 with no input 3d grid.");

  cwp_String keyloc = NULL;  
  if(countparval("keyloc") > 0) {
    if(inloc!=1 || is3d==1) err ("Error: to use keyloc= must be inloc=1 and also not a 3d.");
    getparstring("keyloc", &keyloc);
  }
  else {
    keyloc = ealloc1(3,1);
    strcpy(keyloc,"cdp");
  }

  if (!getparint("check", &icheck)) icheck = 0;
  if(icheck<0 || icheck>1) err("**** Error: check= value out of range.");
  if(is3d==0 && icheck!=0) err("**** Error: cannot specify icheck= with no input 3d grid.");

  getparstring("qin", &Pname);

  getparstring("qout", &Qname);
  if(Pname != NULL && Qname != NULL && strcmp(Pname,Qname) == 0) 
    err("**** Error: qout= output Q-file name must be different than qin= input Q-file name.");
      
  fpQ = fopen(Qname, "w");
  if (fpQ == NULL) err("qfile error: output Q-file did not open correctly.");

  if(countparval("outind")>0) {
    if(countparval("outlist")>0) err("error: outlist= and outind= cannot both be specified.");
    if(countparval("outind")>3) err("error: outind= too long (maximum 3 values: first,last,inc).");
    outind = ealloc1double(3); 
    outind[0] = 0.;
    outind[1] = 3000.;
    outind[2] = 1000.;
    getpardouble("outind",outind);

    if(outind[2]==0. || (outind[1]-outind[0])/outind[2] < -0.000005) 
      err("error: outind= values are wrong.");

    numdind = (int) (1.00001 + (outind[1] - outind[0]) / outind[2]);

    dind = ealloc1double(numdind);
    for (i=0; i<numdind; ++i) dind[i] = outind[0] + i*outind[2];

  }
  else if(countparval("outlist")>0) {
    numdind = countparval("outlist");
    dind = ealloc1double(numdind); 
    getpardouble("outlist",dind);
  }

  if(countparval("navrg") > 2) err("**** Error: navrg= list can only have 1 or 2 values ");
  if(countparval("navrg") < 1) {
    navrg[0] = 10;
    navrg[1] = 10;
  }
  else {
    getparint("navrg",navrg);
    if(countparval("navrg") == 1) navrg[1] = navrg[0];
  }
  if(navrg[0]<-1) err("**** Error in first navrg value. It must be greater or equal to -1 ");
  if(navrg[1]<-1) err("**** Error in second navrg value. It must be greater or equal to -1 ");

  if(countparval("nsmth") > 2) err("**** Error: nsmth= list can only have 1 or 2 values ");
  if(countparval("nsmth") < 1) {
    nsmth[0] = navrg[0];
    nsmth[1] = navrg[1];
   if(nsmth[0] == -1) nsmth[0] = 0;
   if(nsmth[1] == -1) nsmth[1] = 0;
  }
  else {
    getparint("nsmth",nsmth);
    if(countparval("nsmth") == 1) nsmth[1] = nsmth[0];
  }
  if(nsmth[0]<0) err("**** Error in first nsmth value. It must be greater or equal to 0 ");
  if(nsmth[1]<0) err("**** Error in second nsmth value. It must be greater or equal to 0 ");

  if(countparval("nback") > 2) err("**** Error: nback= list can only have 1 or 2 values ");
  if(countparval("nback") < 1) {
    nback[0] = 3;
    nback[1] = 3;
  }
  else {
    getparint("nback",nback);
    if(countparval("nback") == 1) nback[1] = nback[0];
  }
  if(nback[0]<0) err("**** Error: first nback= must be greater than or equal to 0 ");
  if(nback[1]<0) err("**** Error: second nback= must be greater than or equal to 0 ");

  if(!getparint("cdpt",&cdpt)) cdpt = 0;   
  if(cdpt<0 || cdpt>1) err("**** Error: cdpt parameter out of range.");
  if(is3d==0 && cdpt==1) err("**** Error: cannot specify cdpt=1 with no input 3d grid.");

  getparstring("formxy",&formxyt);
  if(formxyt==NULL) {
    lenformxy = 4;
    formxy = ealloc1(lenformxy,1);
    strcpy(formxy,"%.2f");
  }
  else {
    if(is3d==0) err("**** Error: cannot specify formxy= with no input 3d grid.");
    lenformxy = strlen(formxyt);
    formxy = ealloc1(lenformxy,1);
    strcpy(formxy,formxyt);
  }

  formxylong = ealloc1(13+(lenformxy+1)*4,1);
  strcpy(formxylong,"Q,%d,%d,%d,%d,");
  j = 14;
  strcpy(formxylong+j,formxy);
  j += lenformxy;
  strcpy(formxylong+j,",");
  j += 1;
  strcpy(formxylong+j,formxy);
  j += lenformxy;
  strcpy(formxylong+j,",");
  j += 1;
  strcpy(formxylong+j,formxy);
  j += lenformxy;
  strcpy(formxylong+j,",");
  j += 1;
  strcpy(formxylong+j,formxy);

/*-------------------------------------------------------------------------  */
/* Specify some standard qnames whose values will be generated automatically */
/* for the output q-file rather than extracted from input parameters.        */ 
/* So we must get rid of them later after we know what names actually exist  */
/* in input data. Set the list up here to keep later code cleaner.           */
/* Also deal with qform format things which will be output to q-file header. */
/* Note that some of following formats are %.20g even though values actually */
/* get written using %d. The %.20g simply makes it easier to use sutoolcsv in*/
/* the rare case where it is desired. (The names on record after C_SU_NAMES  */
/* will still need manual editing by the users in order to use sutoolcsv).   */

  qname = ealloc1(20+ktuple,sizeof(cwp_String *)); /* more than enough memory*/
  qform = ealloc1(20+ktuple,sizeof(cwp_String *)); /* more than enough memory*/

  numstandard = 0;
  qname[numstandard] = ealloc1(strlen(keyloc),1);
  strcpy(qname[numstandard],keyloc); 
  qform[numstandard] = ealloc1(5,1);
  strcpy(qform[numstandard],"%.20g"); 
  numstandard++;
  qname[numstandard] = ealloc1(4,1);
  strcpy(qname[numstandard],"cdpt");
  qform[numstandard] = ealloc1(5,1);
  strcpy(qform[numstandard],"%.20g");  
  numstandard++;

  if(is3d>0) {
    qname[numstandard] = ealloc1(3,1);
    strcpy(qname[numstandard],"igi");
    qform[numstandard] = ealloc1(5,1);
    strcpy(qform[numstandard],"%.20g"); 
    numstandard++;
    qname[numstandard] = ealloc1(3,1);
    strcpy(qname[numstandard],"igc");
    qform[numstandard] = ealloc1(5,1);
    strcpy(qform[numstandard],"%.20g"); 
    numstandard++;
    qname[numstandard] = ealloc1(2,1);
    strcpy(qname[numstandard],"gx");
    qform[numstandard] = ealloc1(lenformxy,1);
    strcpy(qform[numstandard],formxy);
    numstandard++;
    qname[numstandard] = ealloc1(2,1);
    strcpy(qname[numstandard],"gy");
    qform[numstandard] = ealloc1(lenformxy,1);
    strcpy(qform[numstandard],formxy);
    numstandard++;
    qname[numstandard] = ealloc1(2,1);
    strcpy(qname[numstandard],"sx");
    qform[numstandard] = ealloc1(lenformxy,1);
    strcpy(qform[numstandard],formxy);
    numstandard++;
    qname[numstandard] = ealloc1(2,1);
    strcpy(qname[numstandard],"sy");
    qform[numstandard] = ealloc1(lenformxy,1);
    strcpy(qform[numstandard],formxy);
    numstandard++;
  }

  getparstring("formtv",&formtvt);
  if(formtvt==NULL) {
    lenformtv = 4;
    formtv = ealloc1(lenformtv,1);
    strcpy(formtv,"%.2f");
  }
  else {
    lenformtv = strlen(formtvt);
    formtv = ealloc1(lenformtv,1);
    strcpy(formtv,formtvt);
  }

  formtvlong = ealloc1(1+lenformtv,1);
  strcpy(formtvlong,",");
  strcpy(formtvlong+1,formtv);

  if(is3d>0) {

    fpR = fopen(Rname, "r");
    if(fpR==NULL) err("error: input K-file did not open correctly.");

    errwarn = 1;
    gridread(fpR,gvals,&errwarn);
    if(errwarn>0) err("error reading grid from K-file.");

    gridset(gvals,&errwarn);

    if(errwarn==1) err ("gridset error: grid_wb cell width must be positive.");
    else if(errwarn==2) err ("gridset error: grid_wc cell width must be positive.");
    else if(errwarn==3) err ("gridset error: corner B is within grid_wb cell width of corner A.");
    else if(errwarn>0) err ("gridset error: returned with some unrecognized error code.");
    else if(errwarn==-1) warn ("gridset warning: corner C is near A and is reset to A.");

    gridcheck(gvals,icheck,&errwarn);
    if(errwarn>0) err ("gridcheck error: returned with some unrecognized error code.");

  } /* end of  if(Rname != NULL) { */

/* ----------------------------------------------------------------------  */
/* At this point we have read-in most parameters.                          */
/* And we have read-in the 3d grid definition (if input).                  */
/* And we have set formats and various other things in preperation         */
/* for making the output q-file (as far as possible).                      */
/*     So, now start on actual read-in of values and interpolation.        */
/* ----------------------------------------------------------------------  */

  fpP = fopen(Pname, "r");
  if(fpP==NULL) err("error: input Q-file did not open correctly.");
    
/* User wants to get rid of some names from the input q-file? Put them on the pname list.*/

  if(countparval("ename")>0) {
    numpname = countparval("ename");
    ename = ealloc1(numpname,sizeof(cwp_String *));
    getparstringarray("ename", ename);
    pname = ealloc1(7+numpname,sizeof(cwp_String *));
    for (i=0; i<numpname; ++i) pname[i] = ename[i];
  }
  else {
    numpname = 0;
    pname = ealloc1(7,sizeof(cwp_String *));
  }

/* Overrdie smoothing default?                                                */

  numsname = 0;
  if(countparval("sname")>0) {
    numsname = countparval("sname");
    sname = ealloc1(numsname,sizeof(cwp_String *));
    getparstringarray("sname", sname);
  }

/* Overrdie just copy default?                                                */

  numcname = 0;
  if(countparval("cname")>0) {
    numcname = countparval("cname");
    cname = ealloc1(numcname,sizeof(cwp_String *));
    getparstringarray("cname", cname);
  }

/* Also get rid of standard names. For 3D, these are ignored later anyway (see nzfirst). */
/* The reason to get rid of them here is not just to save memory, but also to skip over  */
/* any badly formated values in these fields of the q-records.                           */

  if(is3d>0) { 
    if(inloc!=1) {
      pname[numpname] = ealloc1(3,1);
      strcpy(pname[numpname],"cdp");
      numpname++;
    }
    if(inloc!=2) {
      pname[numpname] = ealloc1(4,1);
      strcpy(pname[numpname],"cdpt");
      numpname++;
    }
    if(inloc!=3) {
      pname[numpname] = ealloc1(3,1);
      strcpy(pname[numpname],"igi");
      numpname++;
      pname[numpname] = ealloc1(3,1); 
      strcpy(pname[numpname],"igc"); 
      numpname++;
    }
    if(inloc!=4 && inloc!=6) {
      pname[numpname] = ealloc1(2,1);
      strcpy(pname[numpname],"gx");
      numpname++;
      pname[numpname] = ealloc1(2,1); 
      strcpy(pname[numpname],"gy");
      numpname++;
    }
    if(inloc!=5 && inloc!=7) {
      pname[numpname] = ealloc1(2,1);
      strcpy(pname[numpname],"sx");
      numpname++;
      pname[numpname] = ealloc1(2,1);
      strcpy(pname[numpname],"sy");
      numpname++;
    }
  }

/* 0 or negative numpname is a flag to NOT store values if they are on pname list. */
/* positive numpname is a flag to ONLY store values if they are on pname list.     */

  numpname = 0 - numpname; 

  getviaqfile(fpP, &pname, &numpname, &iztuple, numdind,   
              &ktuple, &ifixd, &RecInfo, &ncdp, 
              &pindepa,  &ndims, &errwarn) ;

  if(errwarn==1) err("getqinfo error: extra C_SU_NAMES record in q-file");
  else if(errwarn==2) err("getqinfo error: extra C_SU_NDIMS record in q-file");
  else if(errwarn==3) err("getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record.");
  else if(errwarn==11) 
    err("readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES");
  else if(errwarn==12) 
    err("readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record.");
  else if(errwarn==22) err("getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record.");
  else if(errwarn==23) err("getviaqfile error: C_SU_NAMES tupled names out-of-order, changed");
  else if(errwarn==24) err("getviaqfile error: C_SU_NDIMS blank where valid number expected");
  else if(errwarn==25) err("getviaqfile error: C_SU_NDIMS non-number where valid number expected");
  else if(errwarn==26) err("getviaqfile error: C_SU_NDIMS value must be same for all members of tuple");
  else if(errwarn==27) err("getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found.");
  else if(errwarn>100) 
    err("getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",errwarn-99);
  else if(errwarn>0) err("getviaqfile error: unrecognized error code %d",errwarn);

/* Note that the ifixd flag indicates what is going on with the INPUT      */
/* q-records or parameter sets. The output q-records from this program are */
/* always fixed (the values in dind are never put on output q-records,     */
/* they are put on the single output C_SU_NDIMS record).                   */
/*                                                                         */

  if(ifixd==0) { 
    err("error: input q-records have varying number of values.");
  }

  if(numdind>0 && ifixd==2) {
    err("error: input does not have tuples so outind= and outlist= cannot be specified.");
  }

  if(is3d>0) {
    lastab = (int) (0.1 + gvals[12]);
    lastac = (int) (0.1 + gvals[13]);
    if(lastab*lastac != ncdp) {
      err("error: input q-file does not have same amount of q-records as 3D Grid cdps (cells).");
    }
  }

/* Check that indepenent dimension values are input in increasing order.   */
/* Note that we could just sort the tuples into increasing order ourselves.*/
/* But that is far too likely to induce errors. One mis-typed value in the */
/* spreadsheet and the results will be subtly bad (and, if you don't know  */ 
/* by now, subtly bad is the worst-kind-of-bad in seismic data processing).*/
   
  if(ifixd==1) { 
    for(i=1; i<RecInfo[0].nto; i++) { /* ifixd==1 means nto same in all records */ 
      if(pindepa[i-1] >= pindepa[i]) err("error: independent dimension values not increasing.");
    }

/* Default the output indepenent dimension values to same as input?   */

    if(numdind==0) {
      numdind = RecInfo[0].nto;
      dind = ealloc1double(numdind); 
      for(i=0; i<numdind; i++) dind[i] = pindepa[i];
    }
  }

/* Honor inloc option (which determines what location values to use). */
/* Note that getviaqfile returned with whatever names it found        */
/* (as long as we did not tell it to eliminate a name).               */

  for (icdp=0; icdp<ncdp; ++icdp) RecInfo[icdp].kinf = ealloc1int(3);

  jnloc1 = -1;
  jnloc2 = -1;

  if(inloc==1) {     /* cdp (required if no grid is input) */
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],keyloc)==0) {
        jnloc1 = i;
        break;
      }
    }
    if(jnloc1<0) err("error: for inloc=1, input must have your keyloc=%s (among non-tuple names).",keyloc);

    if(is3d>0) {
      for (icdp=0; icdp<ncdp; ++icdp) {
        RecInfo[icdp].kinf[0] = lrint(RecInfo[icdp].dlots[jnloc1]); 
        gridcdpic(gvals,RecInfo[icdp].kinf[0],RecInfo[icdp].kinf+1,RecInfo[icdp].kinf+2);
        if(RecInfo[icdp].kinf[1] < -2147483644) 
          err("error: for inloc=1, input cdp %d is not in grid",RecInfo[icdp].kinf[0]);
      }
    }
    else {
      for (icdp=0; icdp<ncdp; ++icdp) {
        RecInfo[icdp].kinf[0] = lrint(RecInfo[icdp].dlots[jnloc1]);
        RecInfo[icdp].kinf[1] = RecInfo[icdp].kinf[0]; /* set these for 2d */ 
        RecInfo[icdp].kinf[2] = 1;                     /* set these for 2d */ 
      }
    }
  }

  else if(inloc==2) { /* cdpt (grid transposed cdps) */            
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"cdpt")==0) {
        jnloc1 = i;
        break;
      }
    }
    if(jnloc1<0) err("error: for inloc=2, input must have cdpt (among non-tuple names).");

    for (icdp=0; icdp<ncdp; ++icdp) {
      k = lrint(RecInfo[icdp].dlots[jnloc1]); 
      gridcdpic90(gvals,k,RecInfo[icdp].kinf+1,RecInfo[icdp].kinf+2); /* set igi,igc */
      if(RecInfo[icdp].kinf[1] < -2147483644) 
        err("error: for inloc=2, input 90 degree cdpt %d is not in grid",k);
      gridiccdp(gvals,RecInfo[icdp].kinf[1],RecInfo[icdp].kinf[2],RecInfo[icdp].kinf); /* set cdp */
    }

  }

  else if(inloc==3) { /* igi,igc (grid indexs)       */            
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"igi")==0) {
        jnloc1 = i;
        break;
      }
    }
    if(jnloc1<0) err("error: for inloc=3, input must have igi (among non-tuple names).");
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"igc")==0) {
        jnloc2 = i;
        break;
      }
    }
    if(jnloc2<0) err("error: for inloc=3, input must have igc (among non-tuple names).");

    for (icdp=0; icdp<ncdp; ++icdp) {
      RecInfo[icdp].kinf[1] = lrint(RecInfo[icdp].dlots[jnloc1]);
      RecInfo[icdp].kinf[2] = lrint(RecInfo[icdp].dlots[jnloc2]);
      gridiccdp(gvals,RecInfo[icdp].kinf[1],RecInfo[icdp].kinf[2],RecInfo[icdp].kinf); /* set cdp */
      if(RecInfo[icdp].kinf[0] < -2147483644) 
        err("error: for inloc=3 input igi,igc = %d,%d are not in grid",
        RecInfo[icdp].kinf[1],RecInfo[icdp].kinf[2]);
    }

  }
  else if(inloc==4) { /* gx,gy (grid coordinates)    */            
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"gx")==0) {
        jnloc1 = i;
        break;
      }
    }
    if(jnloc1<0) err("error: for inloc=4, input must have gx (among non-tuple names).");
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"gy")==0) {
        jnloc2 = i;
        break;
      }
    }
    if(jnloc2<0) err("error: for inloc=4, input must have gy (among non-tuple names).");

    for (icdp=0; icdp<ncdp; ++icdp) {
      xg = RecInfo[icdp].dlots[jnloc1];
      yg = RecInfo[icdp].dlots[jnloc2];
      gridgridxyrawxy(gvals,xg,yg,&xw,&yw); /* convert to world xys, then cdp,igi,igc */
      gridrawxycdpic(gvals,xw,yw,RecInfo[icdp].kinf,RecInfo[icdp].kinf+1,RecInfo[icdp].kinf+2);
      if(RecInfo[icdp].kinf[0] < -2147483644) 
        err("error: for inloc=4, input gx,gy = %g,%g are not in grid",xg,yg);
    }

  }
  else if(inloc==5) { /* sx,sy (world coordinates)   */            
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"sx")==0) {
        jnloc1 = i;
        break;
      }
    }
    if(jnloc1<0) err("error: for inloc=5, input must have sx (among non-tuple names).");
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"sy")==0) {
        jnloc2 = i;
        break;
      }
    }
    if(jnloc2<0) err("error: for inloc=5, input must have sy (among non-tuple names).");

    for (icdp=0; icdp<ncdp; ++icdp) {
      xw = RecInfo[icdp].dlots[jnloc1];
      yw = RecInfo[icdp].dlots[jnloc2];
      gridrawxycdpic(gvals,xw,yw,RecInfo[icdp].kinf,RecInfo[icdp].kinf+1,RecInfo[icdp].kinf+2);
      if(RecInfo[icdp].kinf[0] < -2147483644) 
        err("error: for inloc=5, input sx,sy = %g,%g are not in grid",xw,yw);
    }
  }
  else if(inloc==6) { /* gx,gy (world coordinates)   */            
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"gx")==0) {
        jnloc1 = i;
        break;
      }
    }
    if(jnloc1<0) err("error: for inloc=6, input must have gx (among non-tuple names).");
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"gy")==0) {
        jnloc2 = i;
        break;
      }
    }
    if(jnloc2<0) err("error: for inloc=6, input must have gy (among non-tuple names).");

    for (icdp=0; icdp<ncdp; ++icdp) {
      xw = RecInfo[icdp].dlots[jnloc1];
      yw = RecInfo[icdp].dlots[jnloc2];
      gridrawxycdpic(gvals,xw,yw,RecInfo[icdp].kinf,RecInfo[icdp].kinf+1,RecInfo[icdp].kinf+2);
      if(RecInfo[icdp].kinf[0] < -2147483644) 
        err("error: for inloc=6, input gx,gy = %g,%g are not in grid",xw,yw);
    }
  }
  else if(inloc==7) { /* sx,sy (grid coordinates)    */            
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"sx")==0) {
        jnloc1 = i;
        break;
      }
    }
    if(jnloc1<0) err("error: for inloc=7, input must have sx (among non-tuple names).");
    for (i=0; i<iztuple; ++i) {
      if(strcmp(pname[i],"sy")==0) {
        jnloc2 = i;
        break;
      }
    }
    if(jnloc2<0) err("error: for inloc=7, input must have sy (among non-tuple names).");

    for (icdp=0; icdp<ncdp; ++icdp) {
      xg = RecInfo[icdp].dlots[jnloc1];
      yg = RecInfo[icdp].dlots[jnloc2];
      gridgridxyrawxy(gvals,xg,yg,&xw,&yw); /* convert to world xys, then cdp,igi,igc */
      gridrawxycdpic(gvals,xw,yw,RecInfo[icdp].kinf,RecInfo[icdp].kinf+1,RecInfo[icdp].kinf+2);
      if(RecInfo[icdp].kinf[0] < -2147483644) 
        err("error: for inloc=7, input sx,sy = %g,%g are not in grid",xg,yg);
    }
  }

//} 

/*--------------------------------------------------------------*/

/* Sort into cdp order */

  qsort(RecInfo,ncdp,sizeof(struct QInfo),compSort);

  if(is3d>0) {
    for (jcdp=1; jcdp<ncdp; ++jcdp) {
      if(RecInfo[jcdp-1].kinf[0] == RecInfo[jcdp].kinf[0]) { 
        err("error: Two q-records input for cdp,igi,igc = %d %d %d",
        RecInfo[jcdp].kinf[0],RecInfo[jcdp].kinf[1],RecInfo[jcdp].kinf[2]);
      }
    }
  }
  else if(strcmp(keyloc,"cdp")==0) {
    for (jcdp=1; jcdp<ncdp; ++jcdp) {
      if(RecInfo[jcdp-1].kinf[0] + 1 != RecInfo[jcdp].kinf[0]) { 
        err("error: input q-records have missing cdp(s) after cdp= %d",
        RecInfo[jcdp-1].kinf[0]);
      }
    }
  }


  checkpars(); 

/* Change the stored dependent dimension(s) values in dlots to values at outind or outlist. */
/* For instance, use input tims,vels pairs and get vels corresponding to output tims.       */
/* First, use qelementout to compute where they are going to be, then put them there.       */
/*   If you look inside qelementout you will see that it is trivial and could easily be     */
/*   computed right here or on-the-fly in the code below. But the point here is to make it  */
/*   clear to future programmers where the dependent dimension values end up in dlots after */
/*   interpolating them to the independent dimension values input in outind or outlist.     */

  kindx = ealloc1int(iztuple+ktuple-1);   
  ksize = ealloc1int(iztuple+ktuple-1);   

  qelementout(iztuple,ktuple,numdind,kindx,ksize);

  klast = kindx[iztuple+ktuple-2] + ksize[iztuple+ktuple-2];

  dswap = ealloc1double(klast); 
  dtog  = ealloc1double(klast); /* could use dswap for dtop, but this clearer */
  ksmth = ealloc1int(klast);   

  for(jcdp=0; jcdp<ncdp; jcdp++) {

    for(j=0; j<iztuple; j++) dswap[kindx[j]] = RecInfo[jcdp].dlots[j];

    for(i=0; i<numdind; i++) {  
      for(k=0; k<ktuple-1; k++) { 
        linterpd(dind[i],pindepa,RecInfo[jcdp].dlots+iztuple+k*RecInfo[jcdp].nto,
                 RecInfo[jcdp].nto,mgtextr,dswap+kindx[k+iztuple]+i);
      } 
    } 

    for(j=0; j<klast; j++) RecInfo[jcdp].dlots[j] = dswap[j];

  } 

/* If their names are on the standard list, we will output them explicitly, */
/* so set their ksize negative as a flag so as not to output them twice.    */

  for(j=0; j<iztuple+ktuple-1; j++) { 
    k = 0;
    for(i=0; i<numstandard; i++) { 
      if(strcmp(qname[i],pname[j])==0) {
        k = 1;
        break;
      }
    }
    if(k==1) ksize[j] = 0 - ksize[j];
  }

/* Note that I chose to smooth all values BUT only copy SOME of them back     */
/* to RecInfo for output. ksmth contains flags that say whether to copy       */
/* the smoothed value from dswap to RecInfo for eventual output.              */
/* Reminder: ksize <= 0 means do not output the value at all.                 */

/* Process the smoothing override list.                                       */

  for(k=0; k<iztuple; k++) ksmth[k] = 0; 

  for(i=0; i<iztuple; i++) {
    if(ksize[i] > 0) {
      for(j=0; j<numsname; j++) {
        if(strcmp(sname[j],pname[i])==0) ksmth[i] = 1; /* not worth break out */
      }
    }
  }

/* Process the just copy override list. Note q-file standards mean all        */
/* multiple value names (tuples) are after all single value names.            */

  for(k=iztuple; k<klast; k++) ksmth[k] = 1; 

  for(k=0,i=ktuple-2; k<ktuple-1; k++,i--) {
    if(ksize[k+iztuple] > 0) {
      for(j=0; j<numcname; j++) {
        if(strcmp(cname[j],pname[i+iztuple])==0) {
          for(n=0; n<ksize[iztuple+k]; n++) ksmth[kindx[iztuple+k]+n] = 0;
        }
      }
    }
  }

/* Allocate some memory for averaging and smoothing.                         */

  isize = 2*navrg[0] + 3; /* plus 3 because navrg might be -1 */
  if(is3d==1 && 2*navrg[1] + 3 > isize) isize = 2*navrg[1] + 3;

  if(2*nsmth[0] + 1 > isize) isize = 2*nsmth[0] + 1;
  if(is3d==1 && 2*nsmth[1] + 1 > isize) isize = 2*nsmth[1] + 1;

  dalla = calloc(ncdp + isize,sizeof(double *));
  if(dalla == NULL) err("**** Unable to allocate dalla pp memory ");
  dinput = calloc((ncdp + isize)*klast,sizeof(double));
  if(dinput == NULL) err("**** Unable to allocate dalla p memory ");

  for (j=0; j<ncdp + isize; j++) {
    dalla[j] = dinput + j*klast;
    for(i=0; i<klast; i++) dalla[j][i] = 0.0;
  }

/*-------------------------------------------------------------------------- */
/* Write the q-file header records. The information is put onto those header */
/* records is what allows OTHER programs to know that tuples actually exist. */

  fprintf(fpQ,"C_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID");
  for(i=0; i<numstandard; i++) fprintf(fpQ,",%s",qform[i]); 
  for(i=0; i<iztuple; i++) { 
    if(ksize[i] > 0) fprintf(fpQ,",%s",formtv); 
  }
  for(i=0; i<numdind; i++) { /* note: if ifixd=2, numdind is 0 */ 
    for(k=0; k<ktuple-1; k++) {
      if(ksize[k+iztuple] > 0) fprintf(fpQ,",%s",formtv); 
    }
  } 

/* Note the additional commas added in order to align the dind values to the */
/* same (spreadsheet) column as their corresponding tuple value in q-records.*/
/* And note that the dind values are also duplicated so each tuple value gets*/
/* its corresponding dind value at the top of the same (spreadsheet) column. */

  if(ifixd!=2) {
    fprintf(fpQ,"\nC_SU_NDIMS,%s",ndims[0]); 
    for(i=1; i<numstandard; i++) fprintf(fpQ,","); 
    for(i=0; i<iztuple; i++) { 
      if(ksize[i] > 0) fprintf(fpQ,","); 
    }
    for(i=0; i<numdind; i++) { 
      for(k=0; k<ktuple-1; k++) {
        if(ksize[k+iztuple] > 0) fprintf(fpQ,",%.15g",dind[i]);
      }
    }
  }

  fprintf(fpQ,"\nC_SU_NAMES\nC_SU_ID");
  for(i=0; i<numstandard; i++) fprintf(fpQ,",%s",qname[i]);  
  for(i=0; i<iztuple; i++) { 
    if(ksize[i] > 0) fprintf(fpQ,",%s",pname[i]);
  }
  for(i=0; i<numdind; i++) { /* note: if ifixd=2, numdind is 0 */ 
    for(k=0; k<ktuple-1; k++) {
      if(ksize[k+iztuple] > 0) fprintf(fpQ,",%s",pname[k+iztuple]);
    }
  }  

  fprintf(fpQ,"\n");

/*--------------------------------------------------------------------------  */
/* perform smoothing? */

  if(is3d!=1) { /* is 2d */

    for (j=0; j<ncdp; j++) {
      for(i=0; i<klast; i++) dalla[j][i] = RecInfo[j].dlots[i];
    }

    if(navrg[0]>0) {
      runsmo (dalla,ncdp+2*navrg[0]+1,klast,nback[0],navrg[0],dtog,&ierr);
      if(ierr==1) err("**** Error in averaging. Less than 5 Q-records in cin file and navrg>0");
      else if(ierr>0) err("**** Error in navrg averaging for cin file. Parameter out-of-range.");
    }
    else if(navrg[0]==-1) {
      for(i=0; i<klast; i++) dtog[i] = 0.0;
      for (j=0; j<ncdp; j++) {
        for(i=0; i<klast; i++) dtog[i] += dalla[j][i];
      }
      for (j=0; j<ncdp; j++) {
        for(i=0; i<klast; i++) dalla[j][i] = dtog[i] / ncdp;
      }
    } 

    if(nsmth[0]>0) {
      runsmo (dalla,ncdp+2*nsmth[0]+1, klast,nback[0],0-nsmth[0],dtog,&ierr);
      if(ierr==1) err("**** Error in smoothing. Less than 5 Q-records in cin file and nsmth>0");
      else if(ierr>0) err("**** Error in nsmth smoothing for cin file. Parameter out-of-range.");
    }  

    for (j=0; j<ncdp; j++) {
      for(i=0; i<klast; i++) {
        if(ksmth[i] > 0) RecInfo[j].dlots[i] = dalla[j][i];
      }
    }

  }
  else { /* is 3d */

    if(navrg[0]!=0 || nsmth[0]>0) {

      for (n=0; n<lastac; n++) {

        for (j=0; j<lastab; j++) {
          for(i=0; i<klast; i++) dalla[j][i] = RecInfo[j+n*lastab].dlots[i];
        }

        if(navrg[0]>0) {
          runsmo(dalla,lastab+2*navrg[0]+1,klast,nback[0],navrg[0],dtog,&ierr);
          if(ierr==1) err("**** Error in averaging. Less than 5 inline cells defined in 3D Grid and navrg>0");
          else if(ierr>0) err("**** Error in navrg averaging for 3D Grid inlines. Parameter out-of-range.");
        }
        else if(navrg[0]==-1) {
          for(i=0; i<klast; i++) dtog[i] = 0.0;
          for (j=0; j<lastab; j++) {
            for(i=0; i<klast; i++) dtog[i] += dalla[j][i];
          }
          for (j=0; j<lastab; j++) {
            for(i=0; i<klast; i++) dalla[j][i] = dtog[i] / lastab;
          }
        }  

        if(nsmth[0]>0) {
          runsmo(dalla,lastab+2*nsmth[0]+1,klast,nback[0],0-nsmth[0],dtog,&ierr);
          if(ierr==1) err("**** Error in smoothing. Less than 5 inline cells defined in 3D Grid and nsmth>0.");
          else if(ierr>0) err("**** Error in nsmth smooth for 3D Grid inlines. Parameter out-of-range.");
        }

        for (j=0; j<lastab; j++) {
          for(i=0; i<klast; i++) {
            if(ksmth[i] > 0) RecInfo[j+n*lastab].dlots[i] = dalla[j][i];
          }
        }

      } /* end of  for (n=0; n<lastac; n++)  */

    } /* end of  if(navrg[0]!=0 || nsmth[0]>0) */

/* Repeat in the other direction. */

    if(navrg[1]!=0 || nsmth[1]>0) {

      for (n=0; n<lastab; n++) {

        for (j=0; j<lastac; j++) {
          for(i=0; i<klast; i++) dalla[j][i] = RecInfo[n+j*lastab].dlots[i];
        }

        if(navrg[1]>0) {
          runsmo(dalla,lastac+2*navrg[1]+1,klast,nback[1],navrg[1],dtog,&ierr);
          if(ierr==1) err("**** Error averaging. Less than 5 crossline cells defined in 3D Grid and navrg>0");
          else if(ierr>0) err("**** Error in navrg average for 3D Grid crosslines. Parameter out-of-range.");
        }

        if(nsmth[1]>0) {
          runsmo(dalla,lastac+2*nsmth[1]+1,klast,nback[1],0-nsmth[1],dtog,&ierr);
          if(ierr==1) err("**** Error smoothing. Less than 5 crossline cells defined in 3D Grid and nsmth>0");
          else if(ierr>0) err("**** Error in nsmth smooth for 3D Grid crosslines. Parameter out-of-range.");
        }
        else if(navrg[1]==-1) {
          for(i=0; i<klast; i++) dtog[i] = 0.0;
          for (j=0; j<lastac; j++) {
            for(i=0; i<klast; i++) dtog[i] += dalla[j][i];
          }
          for (j=0; j<lastac; j++) {
            for(i=0; i<klast; i++) dalla[j][i] = dtog[i] / lastac;
          }
        }  

        for (j=0; j<lastac; j++) {
          for(i=0; i<klast; i++) {
            if(ksmth[i] > 0) RecInfo[n+j*lastab].dlots[i] = dalla[j][i];
          }
        }

      } /* end of  for (n=0; n<lastab; n++)  */

    } /* end of  if(navrg[1]!=0 || nsmth[1]>0) */

  } /* end of  if(is3d!=1) ... else */ 

  if(cdpt==1) {
    for(jcdp=0; jcdp<ncdp; jcdp++) { /* set kinf[0] to transposed cdp number  */
      gridiccdp90(gvals,RecInfo[jcdp].kinf[1],RecInfo[jcdp].kinf[2],RecInfo[jcdp].kinf);
    }

    qsort(RecInfo,ncdp,sizeof(struct QInfo),compSort); /* sort to kinf[0]     */

    for(jcdp=0; jcdp<ncdp; jcdp++) { /* reset kinf[0] to cdp number */
      gridiccdp(gvals,RecInfo[jcdp].kinf[1],RecInfo[jcdp].kinf[2],RecInfo[jcdp].kinf); 
    }
  }

/* Output */

  for(jcdp=0; jcdp<ncdp; jcdp++) {

    if(is3d>0) {
      gridiccdp90(gvals,RecInfo[jcdp].kinf[1],RecInfo[jcdp].kinf[2],&icdpt);
      gridicgridxy(gvals,RecInfo[jcdp].kinf[1],RecInfo[jcdp].kinf[2],&xg,&yg);
      gridicrawxy(gvals,RecInfo[jcdp].kinf[1],RecInfo[jcdp].kinf[2],&xw,&yw);
      fprintf(fpQ,formxylong,RecInfo[jcdp].kinf[0],icdpt,
              RecInfo[jcdp].kinf[1],RecInfo[jcdp].kinf[2],xg,yg,xw,yw);
    }
    else {
      fprintf(fpQ,"Q,%d,1",RecInfo[jcdp].kinf[0]);
    }

    for(j=0; j<klast; j++) dswap[j] = RecInfo[jcdp].dlots[j];

    for(j=0; j<iztuple; j++) {
      if(ksize[j] > 0) fprintf(fpQ,formtvlong,dswap[kindx[j]]);
    }

    for(i=0; i<numdind; i++) { /* if ifixd=2, numdind is 0 */ 
      for(k=ktuple-2; k>=0; k--) { 
        if(ksize[k+iztuple] > 0) fprintf(fpQ,formtvlong,dswap[kindx[k+iztuple]+i]);
      } 
    } 

    fprintf(fpQ,"\n");

  }

  return(CWP_Exit());  

/*------------------------------------------------------------------------------*/
 
} /* end of susmoqcsv */

/*                                                                    */
/* ------------------------------------------------------------------ */
/*    Running average.  (copied from sufarldcsv Feb 2024).            */
/* Input arguments:                                                   */
/* dall[ncdp][klast] = the values.                                    */
/* ncdp  = length of first  dimension of dall.                        */
/* klast = length of second dimension of dall.                        */
/* nmin  = minimum length in operator. >=0 and <=nmax.                */
/* nmax  = maximum length in operator. !=0.                           */
/*         nmin and |nmax| refer to a number of points in the first   */
/*         dimension of dall. The average for point i is all points   */
/*         from i-nmax to i-nmin and i+nmin to i+nmax inclusive       */ 
/*         (i.e. a symmetric average). For nmin of 0 the average is   */
/*         all points from i-nmax to i+nmax.                          */
/*         If nmax is positive, output the average. If nmax negative, */
/*         output the difference between the average and the original */
/*         value, added back to the original value.                   */
/* dtog    work buffer of length klast                                */
/*                                                                    */
/* Note: Realistically, nmax should be considerably less than ncdp/2  */
/*       Returns ierr=1 if 2*nmax+6 > ncdp)                           */
/*                                                                    */
/* Output arguments:                                                  */
/* dall[ncdp][klast] = contains the averaged values from dall[nmax]   */
/*                     to dall[ncdp-nmax]. The values outside that    */
/*                     are nonsense. In other words, on input, at     */
/*                     least the first nmax and last nmax locations   */
/*                     must be extrapolated or extra points that you  */
/*                     do not mind having twimmed back. The ends are  */
/*                     used by this routine but then altered for its  */
/*                     its own purposes.                              */
/* ierr = 0 (good).                                                   */
/*      = 1 too few points (ncdp) for nmax (number to average).       */
/*      = 2 some input argument is out-of-range.                      */

void runav (double **dall, int ncdp, int klast, int nmin, int nmax, double *dtog, int *ierr) {

  int iavr = 1;
  if(nmax<0) {
    nmax = 0 - nmax;
    iavr = 0;
  }

  if(ncdp<1 || klast<1 || nmin<0 || nmax<nmin) {
    *ierr = 2;
    return;
  }
  if(2*(nmax+3) > ncdp) {
    *ierr = 1;
    return;
  }

  *ierr = 0;

  int i = 0;
  int k = 0;
  int n = 0;

/* Basically this routine proceeds by keeping a running total.   */
/* As we move to the next location, the operator adds the values */
/* for the location that enters the sum range, and subtracts     */
/* the values for the location that leaves the sum range.        */
/* When nmin is non-zero the sum operator also has an inside     */
/* range which is added/subtracted in the same way.              */

/* So, first initialize the sum for the starting location        */
/* (which is also nmax since summing operator is symmetrical).   */

  n = 0;
  for (k=0; k<klast; k++) dtog[k] = 0.0;

  for (i=0; i<nmax*2+1; i++) {
    n++;
    for (k=0; k<klast; k++) dtog[k] += dall[i][k];
  } 

  for (i=nmax-nmin+1; i<nmax+nmin; i++) {
    n--;
    for (k=0; k<klast; k++) dtog[k] -= dall[i][k];
  } 

/* Add next values (at i+1+nmax) and subtract trailing values (at i-nmax)     */ 
/* Note this means the resulting totals are actually for the i+1 location.    */ 

  for (i=nmax; i<ncdp-nmax-1; i++) {

    for (k=0; k<klast; k++) dtog[k] += dall[i+1+nmax][k] - dall[i-nmax][k];

    if(nmin>0) {
      for (k=0; k<klast; k++) dtog[k] += dall[i+1-nmin][k] - dall[i+nmin][k];
    }
  
/* Temporarily store the values into the location that we are just beyond.    */

    if(iavr==1) {
      for (k=0; k<klast; k++) dall[i-nmax][k] = dtog[k] / n;
    }
    else {
      for (k=0; k<klast; k++) dall[i-nmax][k] = dall[i+1][k] + dall[i+1][k] - dtog[k] / n;
    }

  }

/* Move values from where they were temporarily stored to where they belong.  */

  for (i=ncdp-nmax-1; i>=nmax; i--) {
    for (k=0; k<klast; k++) dall[i][k] = dall[i-nmax][k];
  }

  return; 

}

/*                                                                    */
/* ------------------------------------------------------------------ */
/*    Running Smoother (or Averager).                                 */
/* Input arguments:                                                   */
/* dalla[msize][klast] = the values plus 2*nmax+1 extra padding space */
/* msize = length of first dimension of dalla. The values must start  */
/*         dalla[0] and dalla must have 2*abs(nmaxi)+1 extra padding. */
/* klast = length of second dimension of dalla.                       */
/* nback = how to linearly extrapolate values for the nmaxi padding   */
/*       1 means use the values of first point and the first+1 point  */
/*               and the values of last point and the last-1 point.   */
/*       2 means use the values of first point and the first+2 point  */
/*               and the values of last point and the last-2 point.   */
/*       n means use the values of first point and the first+n point  */
/*               and the values of last point and the last-n point.   */
/*       0 means just hold values of first point constant to that end */
/*               and  hold values of last point constant to that end  */
/* nmaxi = maximum length in averging or smoothing operator.          */
/*         If >0 average all points around i from i-nmaxi to i+nmaxi. */
/*         If <0 smooth all points. Smoothing uses the averaging      */
/*         operator twice so affect is from i-2*nmaxi to i+2*nmaxi.   */
/*         (Smoothing does more than just apply the averaging twice). */
/* dtog    work buffer of length klast                                */
/*                                                                    */
/* Output arguments:                                                  */
/* dalla[msize][klast] = contains the output values from dalla[0]     */
/*                     to dalla[msize-2*abs(nmaxi)-1]. Values beyond  */
/*                     that are nonsense.                             */
/* ierr = 0 (good).                                                   */
/*      = 1 too few points (msize) for nmaxi.                         */
/*          (it is an error if 2*abs(nmaxi)+6 is greater than msize)  */
/*      = 2 some input argument is out-of-range.                      */

void runsmo (double **dalla, int msize, int klast, int nback, int nmaxi, double *dtog, int *ierr) {

  int i = 0;
  int j = 0;
  int nmin = 0;
  int nmax = 0;
  int lastcdp = 0;
  double dt = 0.;
  double dback = 0.;

  if(nmaxi<0) nmax = 0 - nmaxi;
  else nmax = nmaxi;

  if(msize<1 || klast<1 || nmax<1 || nback<0) {
    *ierr = 2;
    return;
  }
  if(2*nmax+6 > msize) {
    *ierr = 1;
    return;
  }

  if(nback>0) dback = 1. / (double)nback;

  *ierr = 0;

  lastcdp = msize - 2*nmax - 1;

/* For smoothing, runav is used twice, first with -nmaxi then with +nmaxi.    */
/* The runav function expects values to be padded by nmax before and after.   */
/* Those padded values should be extrapolated from the true values at ends.   */
/* The runav function shifts its output, which needs to be un-shifted.        */
/* The various shifts could probably be optimized (go ahead, if you dare).    */

  for (j=lastcdp-1; j>=0; j--) {
    for(i=0; i<klast; i++) dalla[j+1+nmax][i] = dalla[j][i]; /* shift it */
  }

  for (j=0; j<nmax+1; j++) {
    for(i=0; i<klast; i++) {
      dt = (dalla[nmax+1][i] - dalla[nmax+1+nback][i]) * dback;
      dalla[j][i] = dalla[nmax+1][i] + (nmax-j+1) * dt;
    }
  }   

  for (j=lastcdp+nmax+1; j<lastcdp+2*nmax+1; j++) {
    for(i=0; i<klast; i++) {
      dt = (dalla[lastcdp+nmax][i] - dalla[lastcdp+nmax-nback][i]) * dback;
      dalla[j][i] = dalla[lastcdp+nmax][i] + (j-lastcdp-nmax  ) * dt;
    }
  }   

  runav (dalla,lastcdp+2*nmax+1,klast,nmin,nmaxi,dtog,ierr); /* note nmaxi here */
  if(*ierr>0) return;

  for (j=0; j<=lastcdp-1; j++) {
    for(i=0; i<klast; i++) dalla[j][i] = dalla[j  +nmax][i]; /* unshift it */
  }

/*----------------------------------------------------------------------------*/
  if(nmaxi>0) return;
/*----------------------------------------------------------------------------*/

  for (j=lastcdp-1; j>=0; j--) {
    for(i=0; i<klast; i++) dalla[j+1+nmax][i] = dalla[j][i]; /* shift it */
  }

  for (j=0; j<nmax+1; j++) {
    for(i=0; i<klast; i++) {
      dt = (dalla[nmax+1][i] - dalla[nmax+1+nback][i]) * dback;
      dalla[j][i] = dalla[nmax+1][i] + (nmax-j+1) * dt;
    }
  }   

  for (j=lastcdp+nmax+1; j<lastcdp+2*nmax+1; j++) {
    for(i=0; i<klast; i++) {
      dt = (dalla[lastcdp+nmax][i] - dalla[lastcdp+nmax-nback][i]) * dback;
      dalla[j][i] = dalla[lastcdp+nmax][i] + (j-lastcdp-nmax  ) * dt;
    }
  }   

  runav (dalla,lastcdp+2*nmax+1,klast,nmin,0-nmaxi,dtog,ierr); /* note nmaxi here */
  if(*ierr>0) return;

  for (j=0; j<=lastcdp-1; j++) {
    for(i=0; i<klast; i++) dalla[j][i] = dalla[j  +nmax][i]; /* unshift it */
  }

  return; 

}
/* -----------------------------------------------------------         */
/* Specify compare function for qsort.                                 */

int compSort (const void * q1, const void * q2) {

  struct QInfo* p1 = (struct QInfo*) q1;
  struct QInfo* p2 = (struct QInfo*) q2;

/* cdp order                                                           */

  if(p1->kinf[0] < p2->kinf[0]) return (-1);
  if(p1->kinf[0] > p2->kinf[0]) return (1); 

  return (0); 

}
