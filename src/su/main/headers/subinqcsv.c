/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUBINQCSV: $Revision: 1.02 $ ; $Date: 2024/01/25 00:00:01 $		*/
 
#include "su.h"
#include "segy.h" 
#include "qdefine.h"
#include "gridread.h"
#include "gridxy.h"
#include "bilinear.h"
#include "linterpd.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SUBINQCSV - Input Q-records or parameters. interpolate, output Q-records.  ",
"									     ",
"  subinqcsv [parameters].   (No traces in or out).                          ",
"									     ",
" Parameters:	         						     ",
"                                                                            ",
" qin=       Input file containing q-records. This parameter is optional,    ",
"            but if you do not specify it, you must specify a set of         ",
"            parameters which you can see lower down in this document.       ",
"            Those parameters allow specification similar to sunmo.          ",
"									     ",
" qout=      Output file for q-records. Must be specified.                   ",
"									     ",
" rfile=     If set, read a K-file containing 3D grid definition. Assumed 2D ", 
"            if no file is specified here. A 3D forces restrictions on the   ", 
"            input locations (they have to form aligned rectangles).         ", 
"            A 2D means that parameters for the igi component of 3Ds will be ",
"            honored, but parameters for the igc component are not allowed.  ",
"									     ",
" check=0   Do not print checking details. The options of this program       ",
"           cannot repair a corrupt input grid. Grid checking is always      ",
"           done on input grid definition and will error-halt if corrupt.    ",
"      =1   print checking details.                                          ",
"									     ",
" inloc=1   Location options. The option here determines which input names   ",
"           contain the locations of all other input values. The option here ",
"           must find the corresponding names within the non-tuple names in  ",
"           the input q-file or in the parameter names you specify (later).  ",
"      =1   cdp (default). Option 1 is required if no grid is input.         ",
"      =2   cdpt (grid transposed cdp numbers).                              ",
"      =3   igi,igc (grid index numbers).                                    ",
"      =4   gx,gy (grid coordinates).  Note option 6.                        ",
"      =5   sx,sy (world coordinates). Note option 7.                        ",
"      =6   gx,gy (world coordinates).                                       ",
"      =7   sx,sy (grid coordinates).                                        ",
"           You may want to use iecho=1 initially with these options since   ",
"           the output q-records will contain the cell centre values         ",
"           corresponding to your choice here (giving you a chance to        ",
"           adjust your input values before going through the error-check    ",
"           that forces the input locations to be in aligned rectangles).    ",
"									     ",
" keyloc=   Location key. Can only be specified if inloc=1 and 2D (no rfile).",
"           Name must be in the input Q-file (or parameters listed below).   ",
"           This allows linear interpolation based on values such as point   ",
"           numbers, station numbers, reciever numbers, or shot numbers.     ",
"           Parameters extrapi and igiout refer to these values.             ",
"									     ",
" extrapi=0        do not extrapolate at ends in igi direction.              ",
"                  (values beyond outer functions are held constant).        ",
"               =1 extrapolate beyond outer functions at both ends           ",
"               =2 extrapolate beyond outer functions only at lower end      ",
"               =3 extrapolate beyond outer functions only at higher end     ",
"               For 2D, this parameter controls extrapolation at line ends.  ",
"                                                                            ",
" extrapc=0        do not extrapolate at ends in igc direction.              ",
"                  (values beyond outer functions are held constant).        ",
"               =1 extrapolate beyond outer functions at both ends           ",
"               =2 extrapolate beyond outer functions only at lower end      ",
"               =3 extrapolate beyond outer functions only at higher end     ",
"                                                                            ",
" extrapt=0        do not extrapolate at ends in independent values.         ",
"                  (values outside range of functions are held constant).    ",
"               =1 extrapolate outside range of functions                    ",
"               =2 extrapolate when less than range of functions             ",
"               =3 extrapolate when greater than range of functions          ",
"                                                                            ",
" igiout=1   Output igi range (default is entire igi range of grid).         ",
"            List with one, two, or three integers.                          ",
"       =f     start at f, increment by 1, end at maximum igi of grid.       ",
"       =f,l   start at f, increment by 1, end at l or maximum igi of grid.  ",
"       =f,l,i start at f, increment by i, end at l or nearest igi that      ",
"              is not greater than maximum igi of grid.                      ",
"             For 2D, f and l default to minimum and maximum numbers within  ",
"             the cdp= list or in the q-records input via qin=. For 2D, they ",
"             are allowed to be set outside that range (you can extrapolate  ",
"             to the true first and last cdp of line).                       ",
"									     ",
" igcout=1   Output igc range (default is entire igc range of grid).         ",
"            Same options as the previous igiout parameter.                  ",
"									     ",
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
" invp2=v    Square and inverse values before applying bilinear weighting.   ",
"            (This is how velocities are treated in sunmo and sunmocsv).     ",
"            The sequence is: linear interpolate to outind or outlist values,",
"            square and inverse, perform weighted-sum based on location,     ",
"            then inverse square-root that sum. The default is to do this if ",
"            the value name starts with v. You can also specify any other    ",
"            characters such as vq, tzm, qtal and so on.                     ",
"       =all use square and inverse for all values.                          ",
"      =none use square and inverse for no values.                           ",
"     =tuple use square and inverse for all values in the tuples.            ",
"  =nontuple use square and inverse for all values not in the tuples.        ",
"									     ",
" ename=     Names to eliminate. These names will be ignored in the input    ",
"            q-records and therefore not be in the output q-records (note    ",
"            that standard 2d and 3d names will still be output even if you  ",
"            specify them here). A primary purpose of this option is to get  ",
"            rid of names in input tuples in order to reduce the complexity  ",
"            of output q-records (in case certain programs cannot handle it).",
"            But you can also get rid of non-tuple names here.               ",
"									     ",
" iecho=0    Full spatial output (do not just output at the input locations).",
"      =1    Just output at the same locations as input.                     ",
"            This still performs interpolation in the independent dimension  ",
"            as specified by outind or outlist. (This option bypasses the    ",
"            error check that input locations must form aligned rectangles). ",
"              Parameters igiout= and igcout= cannot be specified.           ",
"                                                                            ",
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
"   The following parameters can only be specified if you                    ",
"   do not specify an input q-file using parameter qin=.                     ",
"									     ",
" cdp=             cdp is a likely parameter name to specify here but        ",
"                  all seismic unix key names are allowed as parameter names.",
"                  So, for instance, fldr= is also allowed to be specified.  ",
"                  These parameters are all optional but at least one of     ",
"                  them must be specified. These are the names which         ",
"                  contain single values at each location (also known as     ",
"            	   the non-tuple names). Each name specified here must have  ",
"            	   the same amount of values (separated by commas).          ",
"            	   For example: cdp=2,7,23,44 means that you must also       ",
"            	   list 4 values if you specify tstat=.                      ",
"              *** Note the names specified must include the names needed    ",
"                  for your inloc= option (by default inloc= requires cdp).  ",
" numa=        *** In addition to all seismic unix key names, the four names ",
" numb=            to the left are ALSO allowed.                             ",
" vuma=        *** See invp2= for why some names here start with v.          ",
" vumb=                                                                      ",
"									     ",
" tupa=            The parameters to the left are all optional. These are    ",
" offs=            the names which have multiple values at each location     ",
" tims=            (also known as the tuple names). Typically, tims= and     ",
" tnmo=            vnmo= are the most likely names to specify. Each name     ",
" dpth=      	   specified here must be repeated same number of times as   ",
" vels=            amount of values in non-tuple names above (except Note 1).",
" vnmo=            (See invp2= for why some names here start with v).        ",
" tupb=    Note 1: Usually, any of the previous parameter names must be      ",
" vupa=            repeated for each value in the non-tuple lists. But the   ",
"                  independent dimension name can also be specified just one ",
"                  time if all the other names here have the same amount of  ",
"                  values as that single set of independent dimension values.",
"          Note 2: The independent dimension defaults to first name that is  ",
"            	   found in order here (NOT their order on the command line).",
"									     ",
"	Example using 3 cpds with velocity semblance scan picks:             ",
"									     ",
"         cdp=7,27,57                                                        ",
"         tims=500,780,1240,2100,4000                                        ",
"         vels=1500,1930,2260,3100,4400                                      ",
"         tims=400,680,2400,3800                                             ",
"         vels=1300,1730,3700,4300                                           ",
"         tims=400,680,880,1140,2400,3800                                    ",
"         vels=1300,1500,1730,2300,3700,4300                                 ",
"									     ",
"       And, if you have a floating datum static at these 3 cdps, then also: ",
"         tstat=5.30,17.44,12.63                                             ",
"									     ",
"   ------------------------------------------------------------------       ",
"									     ",
"  Notes:                                                       	     ",
"									     ",
"   All input values are stored in double-precision (8 byte float).          ",
"									     ",
"   This program does not care if the values are seconds or milliseconds     ",
"   or feet or metres. To this program they are all just numbers. However,   ",
"   default of parameter invp2= treats names starting with v differently     ",
"   for weighting.                                                           ",
"									     ",
"									     ",
" For 3D, user needs to input locations which form aligned rectangles.       ",
" That is, how-ever-many grid inlines the user chooses to put locations on,  ",
" there must always be the same number of locations on each inline and those ",
" functions must be located at the same grid crosslines. For instance, if    ",
" the user inputs locations for inline 7 at crosslines 15,25,40 then the     ",
" the user must alos input locations at crosslines 15,25,40 for any other    ",
" inlines that the user wants to supply locations for. (If user is lucky     ",
" enough that the grid has 100 inlines, then the input locations could be at ",
" CDPs 715,725,740 and 1115,1125,1140 and 2015,2025,2040 and 2515,2525,2540. ",
" Note that the inloc= option allows different values to be used to specify  ",
" locations but they are all translated to inline and crossline numbers      ",
" of the cell centre using the input 3D grid definition.                     ",
"									     ",
" Bilinear interpolation is done if the output location is surrounded by 4   ",
" input locations. If the output location is not surrounded by 4 input       ",
" locations, the result depends on the extrapi and extrapc options.          ",
" The result can be any combination of linear interpolation, linear          ",
" extrapolation, and constant extrapolation. If input locations only exist   ",
" along 1 inline or 1 crossline the result is linear interpolation in that   ", 
" direction (and linear or constant extrapolation at locations beyond ends). ",
"									     ",
NULL};

/*   Created: Nov 2021: Andre Latour                                         */ 
/**************** end self doc *******************************************/

segy tr;

struct QInfo *RecInfo; /* Storage for all function location value pointers */

int compSort2 (const void * q1, const void * q2) ; /* comparison function for qsort  */

#define SQ(x) ((x))*((x))

int main(int argc, char **argv) {

  int ncdp = 0;		/* number of cdps specified */
  int icdp = 0;		/* index into arrays dimensioned by ncdp */
  int icdpt = 0;        /* 90 degree grid cdp number             */
  int jcdp = 0;         /* index into arrays dimensioned by ncdp */

  int is3d = 1;         /* flag for 3d */
  int inloc = 1;        /* flag for 3d input location names        */
  int jnloc1 = -1;      
  int jnloc2 = -1;      

  int mgiextr = 0;      /* for igi extrapolation option            */
  int mgcextr = 0;      /* for igc extrapolation option            */
  int mgtextr = 0;      /* for independent dimension extrapolation */
  int *igiout = NULL;   /* igi range for output to Q-records       */
  int *igcout = NULL;   /* igc range for output to Q-records       */
  int icheck = 0;       /* a command line parameter value          */
  int iecho = 0;        /* a command line parameter value          */
  int invp2 = 0;          /* a command line parameter value        */
  cwp_String nvp2 = NULL; /* related to invp2 parameter above      */


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

  int knownnames = 93;         /* this set of variables is related to   */
  int ltuple = 84;             /* handling of parameter and value names */
  cwp_String *zname = NULL;                                         
  cwp_String *pname = NULL;                                         
  cwp_String *ndims = NULL;                                                 
  cwp_String *qname = NULL;                                              
  cwp_String *ename = NULL;                                              
  cwp_String *qform = NULL;   
  int *ihere = NULL;
  int numpname = 0;
  int numstandard = 0;

  double *outind = NULL;  /* this set of variables is related to   */
  double *pindepa = NULL; /* input of independent dimension values */ 
  int numdind = 0;
  double *dind = NULL;
  double *dswap = NULL;
	
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
  int errwarn = 0; 
  int minigi = 0;                                                      
  int maxigi = 0;                                                      
  int minigc = 0;                                                          
  int maxigc = 0;                                                          
  int *kindx = NULL;
  int *ksize = NULL;
  int *knvp2 = NULL;
  int klast = 0;
  int kflag = 0;
  int mgi_tot = -1;
  int mgc_tot = 0;
  int *mgi = NULL;
  int *mgc = NULL;
  int kigi = 0;
  int kigc = 1;  
  int mgix = 0;
  int mgcx = 0;
  int ndxi = 0; 
  int ndxc = 0;
  int mdxi = 0; 
  int mdxc = 0;
  int mgi_totdeg = 0;
  double wi = 0.;
  double wc = 0.;
  double xw = 0.;
  double yw = 0.;
  double xg = 0.;
  double yg = 0.;

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

  if (!getparint("extrapi", &mgiextr)) mgiextr = 0;
  if(mgiextr<0 || mgiextr>3) err ("error: extrapi= option not in range ");

  if (!getparint("extrapc", &mgcextr)) mgcextr = 0;
  if(is3d==0 && mgcextr!=0) err("**** Error: cannot specify mgcextr= with no input 3d grid.");
  if(mgcextr<0 || mgcextr>3) err ("error: extrapc= option not in range ");

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

/* Read iecho and check igiout and igcout, (but do not read them until defaults available). */

  if (!getparint("iecho", &iecho)) iecho = 0;
  if(iecho<0 || iecho>1) err ("error: iecho= option not in range ");
  if(iecho>0 && (countparval("igiout")>0 || countparval("igcout")>0)) {
    err("error: iecho>0 cannot be specified with igiout= or igcout=");
  }
  if(countparval("igiout")>3) err("error: igiout= has too long a list.");
  if(countparval("igcout")>3) err("error: igcout= has too long a list.");

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

  getparstring("invp2",&nvp2);
  if(nvp2==NULL) { /* default to v */
    invp2 = 1;
    nvp2 = ealloc1(1,1);
    strcpy(nvp2,"v"); 
  }
  else if(strcmp(nvp2,"none")==0) invp2 = 0;
  else if(strcmp(nvp2,"all")==0) invp2 = 2;
  else if(strcmp(nvp2,"tuple")==0) invp2 = 3;
  else if(strcmp(nvp2,"nontuple")==0) invp2 = 4;
  else invp2 = 1; /* if a string is specified */

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
/*                                                                         */
/* Note: There is nothing special about using key names as parameter names.*/
/*       So, for instance, names numa,numb,vuma,vumb are also allowed.     */
/*       But I need to scan for the names, so I must have some idea what   */
/*       they are. This would be a lot easier if getparmnames existed.     */
/*   *** Yes, these values are often used to update trace headers so it    */
/*       is better for the USERs to use key names. But THIS program itself */
/*       does not care whether they are key names are not.                 */
/* ----------------------------------------------------------------------  */

  if(Pname == NULL) { /* no q-file, so readin via command line parameters */ 

    pname = ealloc1(knownnames,sizeof(cwp_String *)); 
    zname = ealloc1(knownnames,sizeof(cwp_String *));
    ihere = ealloc1int(knownnames);
    for(j=0; j<knownnames; j++) {
      pname[j] = ealloc1(8,1);
      zname[j] = ealloc1(8,1);
      strcpy(pname[j],"null");
      ihere[j] = 0;
    }
    strcpy(zname[0],"cdp");
    strcpy(zname[1],"cdpt");
    strcpy(zname[2],"fldr");
    strcpy(zname[3],"grnors");
    strcpy(zname[4],"grnofr");
    strcpy(zname[5],"grnlof");
    strcpy(zname[6],"gaps");
    strcpy(zname[7],"igi");
    strcpy(zname[8],"igc");
    strcpy(zname[9],"gx");
    strcpy(zname[10],"gy");
    strcpy(zname[11],"sx");
    strcpy(zname[12],"sy");
    strcpy(zname[13],"numa");
    strcpy(zname[14],"numb");
    strcpy(zname[15],"vuma");
    strcpy(zname[16],"vumb");
    strcpy(zname[17],"tracl");
    strcpy(zname[18],"tracr");
    strcpy(zname[19],"tracf");
    strcpy(zname[20],"ep");
    strcpy(zname[21],"trid");
    strcpy(zname[22],"nvs");
    strcpy(zname[23],"nhs");
    strcpy(zname[24],"duse");
    strcpy(zname[25],"offset"); 
    strcpy(zname[26],"gelev");
    strcpy(zname[27],"selev");
    strcpy(zname[28],"sdepth");
    strcpy(zname[29],"gdel");
    strcpy(zname[30],"sdel");
    strcpy(zname[31],"swdep");
    strcpy(zname[32],"gwdep");
    strcpy(zname[33],"scalel");
    strcpy(zname[34],"scalco");
    strcpy(zname[35],"counit");
    strcpy(zname[36],"wevel");
    strcpy(zname[37],"swevel");
    strcpy(zname[38],"sut");
    strcpy(zname[39],"gut");
    strcpy(zname[40],"sstat");
    strcpy(zname[41],"gstat");
    strcpy(zname[42],"tstat");
    strcpy(zname[43],"laga");
    strcpy(zname[44],"lagb");
    strcpy(zname[45],"delrt");
    strcpy(zname[46],"muts");
    strcpy(zname[47],"mute");
    strcpy(zname[48],"ns");
    strcpy(zname[49],"dt");
    strcpy(zname[50],"gain");
    strcpy(zname[51],"corr");
    strcpy(zname[52],"sfs");
    strcpy(zname[53],"sfe");
    strcpy(zname[54],"slen");
    strcpy(zname[55],"styp");
    strcpy(zname[56],"stas");
    strcpy(zname[57],"stae");
    strcpy(zname[58],"tatyp");
    strcpy(zname[59],"afilf");
    strcpy(zname[60],"afils");
    strcpy(zname[61],"nofilf");
    strcpy(zname[62],"nofils");
    strcpy(zname[63],"lcf");
    strcpy(zname[64],"hcf");
    strcpy(zname[65],"lcs");
    strcpy(zname[66],"hcs");
    strcpy(zname[67],"year");
    strcpy(zname[68],"day");
    strcpy(zname[69],"hour");
    strcpy(zname[70],"minute");
    strcpy(zname[71],"sec");
    strcpy(zname[72],"timbas");
    strcpy(zname[73],"trwf");
    strcpy(zname[74],"otrav");
    strcpy(zname[75],"d1");
    strcpy(zname[76],"f1");
    strcpy(zname[77],"d2");
    strcpy(zname[78],"f2");
    strcpy(zname[79],"ungpow");
    strcpy(zname[80],"unscale");
    strcpy(zname[81],"ntr");
    strcpy(zname[82],"mark");
    strcpy(zname[83],"shortpad");
    ltuple = 84; /* first tuple name is at 84 */
    strcpy(zname[84],"tupa"); 
    strcpy(zname[85],"offs");
    strcpy(zname[86],"tims");
    strcpy(zname[87],"tnmo");
    strcpy(zname[88],"dpth");
    strcpy(zname[89],"vels");
    strcpy(zname[90],"vnmo");
    strcpy(zname[91],"tupb");
    strcpy(zname[92],"vupa");

    ktuple = 0;
    numpname = 0;

    for(j=0; j<knownnames; j++) {
      if (countparname(zname[j])>0) ihere[j] += 1;
    }

    for(j=0; j<ltuple; j++) { 
      if(ihere[j] > 0) {
        pname[numpname] = zname[j];
        numpname++;
      }
    }
    iztuple = numpname;

    for(j=ltuple; j<knownnames; j++) { 
      if(ihere[j] > 0) {
        pname[numpname] = zname[j];
        numpname++;
      }
    }

/* Note that this program (subinqcsv) has a complicated default for pname */
/* because it allows various parameter set names for input values.        */
/* But for other programs, just set pname[0] to the independent values    */
/* parameter set name and pname[1] to the dependant values name.          */
/* For instance, sunmocsv uses pname[0]="tims" pname[1]="vels" ktuple=2   */
/* and  program sumutecsv uses pname[0]="offs" pname[1]="tims" ktuple=2   */

   getviacommand(&pname, &numpname, &iztuple, numdind,
                 &ktuple, &ifixd, &RecInfo, &ncdp,
                 &pindepa, &ndims, &errwarn) ;

    if(errwarn==1) err("getviacommand error: no non-tuple name passed in.");
    else if(errwarn==2) err("getviacommand error: non-tuple names have different amounts of values.");
    else if(errwarn==3) err("getviacommand error: independent dimension parameter is empty.");
    else if(errwarn==4) err("getviacommand error: an independent dimension parameter is empty.");
    else if(errwarn==5) err("getviacommand error: members of tuple have different amounts at same location."); 
    else if(errwarn>0) err("getviacommand error: returned unrecognized error number = %d",errwarn);

  }

  else { /* Read-in via q-file? */ 

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

  } /* end read-in via q-file */

/* A few error-checks now. Note here the conceptual difference between  */
/* these errors and errors returned from getviaqfile and getviacommand. */
/* These errors are things that subinqcsv decided to call errors even   */
/* though we could choose to ignore the tuples or outind= or outlist=   */

  if(numdind==0 && ifixd!=2) {
    err("error: input has tuples so outind= or outlist= must be specified.");
  }
  else if(numdind>0 && ifixd==2) {
    err("error: input does not have tuples so outind= and outlist= cannot be specified.");
  }

/* Note that the ifixd flag indicates what is going on with the INPUT      */
/* q-records or parameter sets. The output q-records from this program are */
/* always fixed (the values in dind are never put on output q-records,     */
/* they are put on the single output C_SU_NDIMS record).                   */
/*                                                                         */
/* The indepenent dimension name may be on pname list (if ifixd=0 the      */
/* independent values are actually in the input q-records or parameters).  */
/* For example: time,velocity pairs from velocity analysis (velans).       */
/* But, the output values are going to be at dind values, which are        */
/* not varying from record to record, so remove that name.                 */
/* In other words: on input the tims values are often part of tuples whose */
/* number varies from record to record. But on output tims is definitely   */
/* not going to be in the output q-records (the single set of tims values  */
/* specified via outind or outlist will end up on the C_SU_NDIMS record).  */

  if(ifixd==0) { 
    ndims[0] = pname[iztuple]; /* save it for print in C_SU_NDIMS.         */
    for(k=iztuple; k<iztuple+ktuple-1; k++) pname[k] = pname[k+1];
  }

/* Check that indepenent dimension values are input in increasing order.   */
/* Note that we could just sort the tuples into increasing order ourselves.*/
/* But that is far too likely to induce errors. One mis-typed value in the */
/* spreadsheet and the results will be subtly bad (and, if you don't know  */ 
/* by now, subtly bad is the worst-kind-of-bad in seismic data processing).*/
   
  if(ifixd==0) {
    for(jcdp=0; jcdp<ncdp; jcdp++) {
      pindepa = RecInfo[jcdp].dlots+iztuple+(ktuple-1)*RecInfo[jcdp].nto;
      for(i=1; i<RecInfo[jcdp].nto; i++) {
        if(pindepa[i-1] >= pindepa[i]) 
          err("error: independent dimension values not increasing (input location = %d)",jcdp+1);
      }
    }
  }
  else if(ifixd==1) {
    for(i=1; i<RecInfo[0].nto; i++) {
      if(pindepa[i-1] >= pindepa[i]) err("error: independent dimension values not increasing.");
    }
  }

/* Honor inloc option (which determines what location values to use). */
/* Note a small detail here: getviacommand has already error-halted   */
/* if it could not access the names needed - because we had to supply */
/* them to it. But getviaqfile returned with whatever names it found  */
/* (as long as we did not tell it to eliminate that name).            */

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

/* Get cdp range. */

  if(is3d>0) {
    minigi = 1;
    maxigi = (int) (0.1 + gvals[12]);
    minigc = 1;
    maxigc = (int) (0.1 + gvals[13]);
  }
  else { 
    minigi =  2147483645;
    maxigi = -2147483645;
    minigc = 1;
    maxigc = 1;
    for (icdp=0; icdp<ncdp; ++icdp) {
      if(RecInfo[icdp].kinf[0] < minigi) minigi = RecInfo[icdp].kinf[0];
      if(RecInfo[icdp].kinf[0] > maxigi) maxigi = RecInfo[icdp].kinf[0];
    }
  } 

/* Now that defaults are available, read igiout and igcout parameters. */

  igiout = ealloc1int(3); 
  igiout[0] = minigi;
  igiout[1] = maxigi;
  igiout[2] = 1;
  if(countparval("igiout")>0) getparint("igiout",igiout); 
  if(igiout[2]<1) err("error: igiout= increment is less than 1");
  if(is3d>0) {
    if(igiout[0]<minigi) err("error: igiout= first is less than minimum grid igi index %d",minigi);
    if(igiout[1]>maxigi) igiout[1] = maxigi;
  }
  if(igiout[0]>igiout[1]) err("error: igiout= first is greater than last. ");

  igcout = ealloc1int(3); 
  igcout[0] = minigc;
  igcout[1] = maxigc;
  igcout[2] = 1;
  if(countparval("igcout")>0) getparint("igcout",igcout);
  if(igcout[2]<1) err("error: igcout= increment is less than 1");
  if(is3d>0) {
    if(igcout[0]<minigc) err("error: igcout= first is less than minimum grid igc index %d",minigc);
    if(igcout[1]>maxigc) igcout[1] = maxigc;
  }
  if(igcout[0]>igcout[1]) err("error: igcout= first is greater than last. ");

/*--------------------------------------------------------------*/

  checkpars(); 

/* Does user just want to echo the input locations? If so, do not sort. And do not check to make */
/* sure user has not input the same location twice. And do not check to see if input locations   */
/* form aligned rectangles as required by bilinear interpolation.                                */

  if(iecho==0) { 

/* For bilinear interpolation, user must input function locations which form aligned rectangles.  */
/* That is, howevermany inlines the user chooses to put functions on, there must be the same      */
/* number of functions on each inline and those functions must be located at the same crosslines. */
/* For instance, if user inputs functions for inline 7 at crosslines 15,25,40 then the user must  */
/* input the functions at crosslines 15,25,40 for any other inlines that the user wants to supply */
/* functions for. The qsplit function enforces that restriction on the user input - and also      */
/* separates the igi and igc values into 2 simple arrays.                                         */
/*   The qsplit function expects the values to be sorted by igi,igc values.                       */
/*   The results of qsplit are then used by binterpfind and binterpapply.                         */

    qsort(RecInfo,ncdp,sizeof(struct QInfo),compSort2);

    for (jcdp=1; jcdp<ncdp; ++jcdp) {
      if(RecInfo[jcdp-1].kinf[1] == RecInfo[jcdp].kinf[1] &&
         RecInfo[jcdp-1].kinf[2] == RecInfo[jcdp].kinf[2]) {
           err("error: Two sets of values input for cdp,igi,igc = %d %d %d",
                RecInfo[jcdp].kinf[0],RecInfo[jcdp].kinf[1],RecInfo[jcdp].kinf[2]);
      }
    }

/* (Note: qsort is a standard c function, qsplit is just a function I wrote in su/lib/qdefine.c). */

    qsplit(RecInfo,ncdp,&mgi,&mgi_tot,&mgc,&mgc_tot,&errwarn);

    if(errwarn>0) err("qsplit error: locations do not form aligned rectangles (detected at cdp= %d)", 
                  RecInfo[errwarn-1].kinf[0]);

  } /* end of  if(iecho==0) { */ 

/* Change the stored dependent dimension(s) values in dlots to values at outind or outlist. */
/* For instance, use input tims,vels pairs and get vels corresponding to output tims.       */
/* First, use qelementout to compute where they are going to be, then put them there.       */
/*   If you look inside qelementout you will see that it is trivial and could easily be     */
/*   computed right here or on-the-fly in the code below. But the point here is to make it  */
/*   clear to future programmers where the dependent dimension values end up in dlots after */
/*   interpolating them to the independent dimension values input in outind or outlist.     */
/*   Follow knvp2 to see how kindx,ksize make things easier to understand and use.          */

  kindx = ealloc1int(iztuple+ktuple-1);   
  ksize = ealloc1int(iztuple+ktuple-1);   

  qelementout(iztuple,ktuple,numdind,kindx,ksize);

  klast = kindx[iztuple+ktuple-2] + ksize[iztuple+ktuple-2];

  dswap = ealloc1double(klast); 

  for(jcdp=0; jcdp<ncdp; jcdp++) {

    for(j=0; j<iztuple; j++) dswap[kindx[j]] = RecInfo[jcdp].dlots[j];

    if(ifixd==0) pindepa = RecInfo[jcdp].dlots+iztuple+(ktuple-1)*RecInfo[jcdp].nto;

    for(i=0; i<numdind; i++) {  
      for(k=0; k<ktuple-1; k++) { 
        linterpd(dind[i],pindepa,RecInfo[jcdp].dlots+iztuple+k*RecInfo[jcdp].nto,
                 RecInfo[jcdp].nto,mgtextr,dswap+kindx[k+iztuple]+i);
      } 
    } 

    for(j=0; j<klast; j++) RecInfo[jcdp].dlots[j] = dswap[j];

  } 

/* Set flags for values to square-and-inverse (usually velocity values). */
/* invp2=0=none  =2=all  =3=tuple  =4=nontuple  =1=string (default "v")  */

  knvp2 = ealloc1int(iztuple+ktuple-1);   
  for(k=0; k<iztuple; k++) { 
    knvp2[k] = 1; /* default is linear weight */
    if(invp2!=0) { 
      if(invp2==2 || invp2==4 || 
        (invp2==1 && strncmp(pname[k],nvp2,strlen(nvp2))) == 0) {
        knvp2[k] = -2; /* 1/(x*x) weight */
      }
    }
  } 
  for(k=0,i=ktuple-2; k<ktuple-1; k++,i--) { /* remember: stored in reverse   */ 
    knvp2[k+iztuple] = 1;                    /* relative to their pname order */
    if(invp2!=0) { 
      if(invp2==2 || invp2==3 || 
        (invp2==1 && strncmp(pname[i+iztuple],nvp2,strlen(nvp2))) == 0) { 
        knvp2[k+iztuple] = -2; /* 1/(x*x) weight */
      }
    }
  }

/* If their names are on the standard list, we will output them explicitly, */
/* so set their ksize negative as a flag so as not to output them twice.    */
/* Also set kflag if there are any that will need inverse square/root.      */

  kflag = 0;
  for(j=0; j<iztuple+ktuple-1; j++) { 
    k = 0;
    for(i=0; i<numstandard; i++) { 
      if(strcmp(qname[i],pname[j])==0) {
        k = 1;
        break;
      }
    }
    if(k==1) ksize[j] = 0 - ksize[j];
    else if(knvp2[j] == -2) kflag = 1;
  }

/* Apply forward invp2 option. */

  if(kflag>0) {
    for(jcdp=0; jcdp<ncdp; jcdp++) {
      for(k=0; k<iztuple+ktuple-1; k++) { 
        if(knvp2[k] == -2) {
          for(i=kindx[k]; i<kindx[k]+ksize[k]; i++) {  
            if(RecInfo[jcdp].dlots[i] > 1.e-98) RecInfo[jcdp].dlots[i] = 
              1./(RecInfo[jcdp].dlots[i] * RecInfo[jcdp].dlots[i]);
            else if(RecInfo[jcdp].dlots[i] < -1.e-98) RecInfo[jcdp].dlots[i] = 
             -1./(RecInfo[jcdp].dlots[i] * RecInfo[jcdp].dlots[i]);
          } 
        } 
      } 
    } 
  } /* end of  if(kflag>0) { */

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

/*-------------------------------------------------------------------------- */
/* OK, so the user just wants to echo the input locations without any        */
/* spatial interpolation. Note that this does not mean output values are the */
/* same as input since dind are the OUTPUT independent dimension values.     */
 
  if(iecho==1) {

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

/* To make iecho=1 as similar to code for iecho=0, copy to dswap. Then apply */
/* inverse of invp2. Result is "no difference" for this iecho=1 case since   */
/* it is not performing the bilinear summing. I could, of course, disable    */
/* forward and inverse of invp2 option for the iecho=1 case, reducing cpu.   */
/* But, it is better testing if the iecho=1 case replicates the full code as */
/* much as possible.                                                         */

      for(j=0; j<klast; j++) dswap[j] = RecInfo[jcdp].dlots[j];

      if(kflag>0) {
        for(k=0; k<iztuple+ktuple-1; k++) { 
          if(knvp2[k] == -2) {
            for(i=kindx[k]; i<kindx[k]+ksize[k]; i++) {  
              if(dswap[i] > 1.e-98) dswap[i] = 1./sqrt(dswap[i]);
              else if(dswap[i] < -1.e-98) dswap[i] = -1./sqrt(0.-dswap[i]);
            } 
          } 
        } 
      } /* end of  if(kflag>0) { */

/* Print out */

      for(j=0; j<iztuple; j++) {
        if(ksize[j] > 0) fprintf(fpQ,formtvlong,dswap[kindx[j]]);
      }

      for(i=0; i<numdind; i++) { /* if ifixd=2, numdind is 0 */ 
        for(k=ktuple-2; k>=0; k--) { 
          if(ksize[k+iztuple] > 0) fprintf(fpQ,formtvlong,dswap[kindx[k+iztuple]+i]);
        } 
      } 

      fprintf(fpQ,"\n");

    } /* end of  for(jcdp=0; jcdp<ncdp; jcdp++) {  */

    return(CWP_Exit()); 
  } /* end of  if(iecho==1) { */

/*------------------------------------------------------------------------------*/
/* OK, so the user wants full spatial interpolation, not just echo.             */

  mgi_totdeg = mgi_tot; /* read explanation later */
  if(mgi_tot==1 || mgc_tot==1) mgi_totdeg = 0;

  ndxi = 0; 
  ndxc = 0;
  mdxi = 0; 
  mdxc = 0;

/* So, cycle through the locations that we want to output.      */

  for(kigc=igcout[0]; kigc<=igcout[1]; kigc+=igcout[2]) {
    for(kigi=igiout[0]; kigi<=igiout[1]; kigi+=igiout[2]) {

/* Print the beginning of each q-record.                        */
/* If 3d, use grid indexes igi,igc to get cdp and coordinates.  */

      if(is3d>0) {
        gridiccdp(gvals,kigi,kigc,&icdp);
        gridiccdp90(gvals,kigi,kigc,&icdpt);
        gridicgridxy(gvals,kigi,kigc,&xg,&yg);
        gridicrawxy(gvals,kigi,kigc,&xw,&yw);
        fprintf(fpQ,formxylong,icdp,icdpt,kigi,kigc,xg,yg,xw,yw);
      }
      else {
        fprintf(fpQ,"Q,%d,1",kigi);
      }

/* find input cdp (higher) locations mgix,mgcx and weights wi,wc */

      if(ncdp>1) { 
        binterpfind(kigi,mgi,mgi_tot,mgiextr,kigc,mgc,mgc_tot,mgcextr,
                    &mgix,&mgcx,&wi,&wc);

/* mgix and mgcx are the locations computed for each direction seperately.     */
/* mgix and mgcx are always returned as the highest of the 2 (near) locations. */
/* (if mgi has 10 locations, mgix is only returned from 1 to 9, never 0).      */
/* So, the 4 locations are: mgix,mgcx  mgix-1,mgcx  mgix,mgcx-1  mgix-1,mgcx-1.*/
/* For those 4, compute the element numbers of the stored functions in RecInfo.*/
/* Note that for the degenerate cases of mgi_tot=1 or mgc_tot=1 the            */
/* mgi_totdeg=0, which results in ndxc=ndxi, which in turn means the second    */
/* two functions passed to binterpapply are the same as first two (which works */
/* because either weight wi or wc will be 0.0). Similarly, if ncdp<2 then      */
/* ndxi,ndxc,mdxi,mdxc remain 0 (binterpapply is passed RecInfo[0] four times).*/

        ndxi = mgix + mgi_tot * (mgcx-1);
        ndxc = ndxi + mgi_totdeg;
        mdxi = ndxi-1;
        mdxc = ndxc-1;
      }

      binterpapply(RecInfo[mdxi].dlots, RecInfo[ndxi].dlots, mgi_tot, wi,
                   RecInfo[mdxc].dlots, RecInfo[ndxc].dlots, mgc_tot, wc,
                   klast,dswap);     

/* Apply inverse of invp2 option. */

      if(kflag>0) {
        for(k=0; k<iztuple+ktuple-1; k++) { 
          if(knvp2[k] == -2) {
            for(i=kindx[k]; i<kindx[k]+ksize[k]; i++) {  
              if(dswap[i] > 1.e-98) dswap[i] = 1./sqrt(dswap[i]);
              else if(dswap[i] < -1.e-98) dswap[i] = -1./sqrt(0.-dswap[i]);
            } 
          } 
        } 
      } /* end of  if(kflag>0) { */

      for(j=0; j<iztuple; j++) {
        if(ksize[j] > 0) fprintf(fpQ,formtvlong,dswap[kindx[j]]);
      }

      for(i=0; i<numdind; i++) { /* if ifixd=2, numdind is 0 */ 
        for(k=ktuple-2; k>=0; k--) { 
          if(ksize[k+iztuple] > 0) fprintf(fpQ,formtvlong,dswap[kindx[k+iztuple]+i]);
        } 
      } 

      fprintf(fpQ,"\n");

    } /* end of   for(kigi=igiout[0]; kigi<igiout[1]; kigi+=igiout[2]) { */
  } /* end of   for(kigc=igcout[0]; kigc<igcout[1]; kigc+=igcout[2]) { */
 
} /* end of subinqcsv */

/* -----------------------------------------------------------         */
/* Specify compare function for qsort.                                 */

int compSort2 (const void * q1, const void * q2) {

  struct QInfo* p1 = (struct QInfo*) q1;
  struct QInfo* p2 = (struct QInfo*) q2;

/* Note I decided so sort to igc,igi order (kinf[2], then kinf[1])     */

  if(p1->kinf[2] < p2->kinf[2]) return (-1);
  if(p1->kinf[2] > p2->kinf[2]) return (1); 
  if(p1->kinf[1] < p2->kinf[1]) return (-1);
  if(p1->kinf[1] > p2->kinf[1]) return (1); 

  return (0); 

}
