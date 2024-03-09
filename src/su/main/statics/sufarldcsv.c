/* Copyright (c) Colorado School of Mines, 2023.*/
/* All rights reserved.                       */

/* SUFARLDCSV: $Revision: 1.01 $ ; $Date: 2023/08/01 00:00:01 $		*/
 
#include "su.h"
#include "segy.h" 
#include <stdbool.h>
#include "qdefine.h"
#include "gridread.h"
#include "gridxy.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SUFARLDCSV - Floating Approximate Robust Local Datum From Q-files.         ",
"									     ",
"              The Approximate Robust Local (ARL) algorithm is intended      ",
"              to produce a reasonable smooth floating datum for the         ",
"              typical spatial distributions of shots and receivers that     ",
"              exist in many 2D and 3D surveys.                              ",
"									     ",
"              This documentation uses the familiar words shot, receiver,    ",
"              cdp, coordinates, and statics. And this program defaults to   ",
"              Seismic Unix key names for various parameters. But nothing    ",
"              in this program behaves differently if other words or names   ",
"              are used. To this program all input values are just numbers.  ",
"              For example, the same values are output whether the input     ",
"              values are named statics, or elevations.                      ",
"									     ",
"     Caution: When using the output of this program you should use          ",
"              cdp or igi,igc in the SUSHIFTCSV program to find the          ",
"              cdp of the traces. Using nearest COORDINATES to find the      ",
"              cdp does not produce cdp-consistent results because cdps      ",
"              only contain one-half of their coordinate boundaries.         ",
"									     ",
"									     ",
"  sufarldcsv [parameters]. (No traces in or out).                           ",
"									     ",
" Parameters:	         						     ",
"                                                                            ",
" sin=sstat.csv  Input csv Q-file (usually containing Source information).   ",
"    =none       No sin file. Value 0 is used for all sin statics.           ",
"									     ",
" sloc=sx,sy     List of sin file location Names.                            ",
"                For 2D surveys, 1 or 2 names are allowed.                   ",
"                For 3D surveys, there must be 2 names and the first name    ",
"                must have the x value and the second name must have the y.  ",
"									     ",
" sstat=sstat    This is the default source static name.                     ",
"									     ",
" smult=1.0      Sin file Static Multiplier. Statics in Q-files should be    ",
"                milliseconds, if not, you can scale them before use here.   ",
"          Note: Also, to reverse static sign, specify negative here.        ",
"                                                                            ",
" rin=rstat.csv  Input csv Q-file (usually containing Receiver information). ",
"    =none       No rin file. Value 0 is used for all receiver statics.      ",
"          Note: Both sin and rin can be none since it might be usefull to   ",
"                get a dout file formatted with 2D values or 3D Grid values. ",
"									     ",
" rloc=gx,gy     List of rin file location Names.                            ",
"                For 2D surveys, 1 or 2 names are allowed.                   ",
"                For 3D surveys, there must be 2 names and the first name    ",
"                must have the x value and the second name must have the y.  ",
"									     ",
" rstat=gstat    This is the default rin file static name.                   ",
"									     ",
" rmult=1.0      Rin File Static Multiplier. Statics in Q-files should be    ",
"                milliseconds, if not, you can scale them before use here.   ",
"          Note: Also, to reverse static sign, specify negative here.        ",
"                                                                            ",
"  rfile=        If specified, read the 3D Grid defined in the K-record.     ", 
"           ***  MUST be specified for 3D and must NOT be specified for 2D.  ",
"           ***  3D Grid must have at least 5 cells in the corresponding     ",
"                direction if navrg>0,0 or nsmth>0,0                         ",
"                See subincsv for 3D Grid documentation.                     ",
"                If not a 3D survey, specify the next parameter (cin file).  ", 
"                                                                            ",
" cin=none       Input csv Q-file containing CDP Information.                ",
"    =none       Default. No cin file. An rfile with a 3D Grid must be input.",
"           ***  MUST be specified for 2D and must NOT be specified for 3D.  ",
"           ***  At least 5 Q-records must exist if nsmth>0 or navrg>0.      ",
"									     ",
" cloc=gx,gy     List of cin file Location Names. Only 1 or 2 names allowed. ",
"           ***  Can only be specified for 2D.                               ",
"          Note: The default is receiver coordinate key names because the    ",
"                cin file is often created from receiver locations.          ",
"									     ",
" cseq=cdp       This is the cin Q-file record sequencing name.              ",
"           ***  Can only be specified for 2D.                               ",
"     =asis      Means use cin Q-records in the order they are input.        ",
"          Note: The input cin Q-record are sorted by these values and that  ",
"                sequence is used by navrg, nsmth, and nback parameters.     ",
"                The averaging and smoothing controlled by those parameters  ",
"                just uses these sequential numbers. That is, those options  ",
"                do not care if cseq actually contains a cdp number, and     ",
"                they do not care whether gaps exist in input cseq numbers.  ",
"									     ",
" cstat=tstat    This is the default output static name. Typically it is the ",
"                sum of the (averaged, smoothed) output shot static and the  ",
"                (averaged, smoothed) output receiver static divided by 2,   ",
"                but if either sin or rin files are not input then it is     ",
"                just a copy of the other (averaged, smoothed) static.       ",
"									     ",
" avrad=1000.0   Radius to Average Shot statics and Average Receiver statics.",
"                Statics are averaged separately for shots and receivers.    ",
"          =-1.  Do not average any statics. Just use nearest shot static    ",
"                and nearest receiver static. You must also set nedge=0.     ",
"          Note: Specify this parameter in units of sloc, rloc, and cloc.    ",
"                By default those are coordinate names so this parameter     ",
"                is a true distance. But for 2D surveys sloc, rloc, and cloc ",
"                can be single values such as point numbers in which case    ",
"                this value is the difference in point numbers.              ",
"									     ",
" nedge=2        Minimum Number of Points To Perform Edge Compensation.      ",
"                Edge compensation is done by a relatively complicated       ",
"                method based on the difference between each shot location   ",
"                and the average location of all shots within avrad distance ",
"                of that shot (the centroid location). The larger the        ",
"                distance to the centroid, the larger the compensation.      ",
"                The same method is used independently for the receivers.    ",
"                Must be greater than 1.                                     ",
"           ***  Cannot specify if avrad is -1.                              ",
"          Note: For 2D surveys this compensation mostly affects points near ",
"                line ends, but will also be larger for gaps within the line.",
"                For 3D surveys this compensation mostly affects points near ",
"                the edges of the areas of coverage, but compensation will   ",
"                also be larger for gaps within the areas of coverage.       ",
"      =0        Do not perform Edge Compensation. Typically, for surveys    ",
"                with slants in the statics near edges or gaps, this results ",
"                in values curving in the opposite direction as the slant.   ",
"									     ",
" cedge=avrad/10 Minimum Distance to Centroid to Perform Edge Compensation.  ",
"                (See previous parameter). If the centroid location is       ",
"                within this distance of the shot or receiver location,      ",
"                do not perform edge compensation for this shot or receiver. ",
"           ***  Cannot specify if avrad is -1. Cannot be greater than avrad.",
"          Note: Specify this parameter in units of sloc, rloc, and cloc.    ",
"                By default those are coordinate names so this parameter     ",
"                is a true distance. But for 2D surveys sloc, rloc, and cloc ",
"                can be single values such as point numbers in which case    ",
"                this value is the difference in point numbers.              ",
"									     ",
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
" --------------                                                             ",
"                                                                            ",
" formxy=%.20g  The C format code for printing all values to the q-records.  ",
"              Note that the default format prints up to 20 digits           ",
"              (but not trailing zeroes to the right of the decimal point).  ",
"         Note: This default is used because all values in input q-files     ",
"               can have full double precision and this program computes in  ",
"               double precision. In particular, 2D cdp profile coordinates  ",
"               and 3D cdp cell centre coordinates often need to be output   ",
"               in full double precision (which is approximately 16 digits). ",
"               In some situations this makes it hard to easily check the    ",
"               output values so you might use %.10g here during testing.    ",
"									     ",
" dout=tstat.csv  Output csv Q-file name.                                    ",
"                 For 2D this file contains:                                 ",
"                    cseq,cloc1,cloc2,sstat,gstat,tstat                      ",
"                     - where cseq,cloc1,cloc2 are the names specified in    ",
"                       those parameters (with names set to numb1 and numb2  ",
"                       and values of 0.0 if they are not being used).       ",
"                     - and sstat,rstat,cstat are the statics associated     ",
"                       with the names specified in those parameters.        ",
"                 For 3D this file contains:                                 ",
"                    cdp,cdpt,igi,igc,sx,sy,gx,gy,sstat,rstat,cstat          ",
"                     - where cdp,cdpt,igi,igc are cell and transposed cell  ",
"                       number and inline and crossline index computed from  ",
"                       the 3D Grid defined in K-record of the input rfile.  ",
"                     - and sx,sy,gx,gy are real-world and grid-relative     ",
"                       cell centre coordinates computed from the 3D Grid.   ",
"                     - and sstat,rstat,cstat are the statics associated     ",
"                       with the names specified in those parameters.        ",
"									     ",
" --------------                                                             ",
"									     ",
" The following 2 parameters affect cpu time, but not results. The search    ",
" for the nearest shot and/or the nearest receiver is done using kdtrees.    ",
"                                                                            ",
" kdist=5   Initial search distance. If positive, this value is added to     ",
"           the distance between the previous 2D cin file sorted Q-record    ",
"           or 3D cell and its nearest shot and nearest receiver and used as ",
"           the initial search distance to find the nearest shot and nearest ",
"           receiver for the next 2D cin file sorted Q-record or 3D cell.    ",
"           If negative, the initial search distance for all locations       ",
"           is set to this absolute value.                                   ",
" kmult=2   Search multiplier.                                               ",
"           If the nearest shot or receiver is not found after searching with",
"           the initial search distance, the search distance is multiplied   ",
"           by this value, and the searches are performed again.             ",
"           This repeats until finding the nearest shot and receiver.        ",
"									     ",
NULL};

/* Created: Oct  2023: Andre Latour                                          */ 
/*  1. To compute a floating datum this program uses an algoithm             */ 
/*     which I call Approximate Robust Local, which just happens             */ 
/*     to also match my initials (Andre Roger Latour).                       */ 
/**************** end self doc ***********************************************/

segy tr;

struct QInfo *SInfo;
struct QInfo *RInfo;
struct QInfo *CInfo;

int locp = -1;

int compSort1 (const void * q1, const void * q2) ; /* comparison function for qsort  */

/* Note: I make no claim that this is a particularly good kd tree implementation.     */
/*       It is not explicitly balanced (it has an option to get approximate balance). */
/* Note: See sunearcsv for details and variations.                                    */

typedef struct node {
   unsigned long long elem; 
   struct node * l;
   struct node * r;
} node;

void connect_nodes (node *tree_nodes, unsigned long long tree_numc, double **tree_dl, int tree_numd,
                   int ihop);

void connect_all (node *tree_nodes, unsigned long long tree_numc, double **tree_dl, int tree_numd);

void find_in (node *tree_nodes, double **tree_dl, int tree_numd,
             double *extent_min, double *extent_max,
             unsigned long long *out_elem, unsigned long long *num_out);

void find_in_rcur (double **tree_dl, int tree_numd, node *now_node, int naxe,
                  double *extent_min, double *extent_max,
                  unsigned long long *out_elem, unsigned long long *num_out);

void find_near_rcur (double **tree_dl, int *tree_nt, int tree_numd, node *now_node, int naxe, 
                    double *extent_min, double *extent_max, double *target,
                    unsigned long long *near_elem, double *near_dist, unsigned long long *num_found);

void cycle_for_near (node *tree_nodes, double **tree_dl, int *tree_nt, int tree_numd,
                    double *extent_min, double *extent_max, double *target, 
                    double init_rad, double rad_scal,  
                    unsigned long long *near_elem, double *near_dist, 
                    unsigned long long *num_found, int *ncycles);

void runsmo (double **dalla, int msize, int klast, int nback, int nmax, double *dtog, int *ierr) ;

void runav (double **dall, int ncdp, int klast, int nmin, int nmax, double *dtog, int *ierr) ;

/*----------------------------------------------------------------------*/

int main(int argc, char **argv) {

  int npoints = 0;

  int ifixd = 0;          /* flag for all tuples same size or vary   */
  int iztuple = 0;        /* element number where first tuple exists */
  int ktuple = 0;         /* type of tuples (2=pairs, 3=triplets)    */

  cwp_String *pname = NULL;                                         
  cwp_String *ndims = NULL;                                                 
  int numpname = 0;

  double *pindepa = NULL;                                               
  int numdind = 0;
	
  cwp_String formxyt=NULL;
  cwp_String formxy=NULL;
  cwp_String formxylong=NULL;
  int lenformxy = 0;

  int i = 0;                                                                         
  int j = 0;                                                                       
  int k = 0;                                                                       
  int n = 0;                                                                       
  int errwarn = 0; 

  int iprint = 0;
  int maygrid = 0;
  int is3d    = 1;
  int icheck= 0;
  cwp_String Kname=NULL;  /* text file name for grid values       */
  FILE *fpK=NULL;         /* file pointer for grid input file     */
  double gvals[999];      /* to contain the grid definition       */
  int j3dcdp = 1;
  int igi = 1;
  int igc = 1;
  int lastab=0;
  int lastac=0;
  int lastcdp=0;

  double kdist = 5.0;
  double kmult = 2.;
  int kcycles = 0;
  int kdadd   = 1;

  cwp_String Sname=NULL;                                             
  FILE *fpS=NULL;                                                       
  cwp_String skey[9];  
  int slocn[9];
  int stree_numd = 0; 
  double smult=1.0;
  cwp_String sstat=NULL; 
  int osloc  = 0;  
  double *stree_dl[9];
  unsigned long long stree_npoints = 0;
  node *stree_nodes = NULL;
  double sdist = 100.0; 
  double sdist2 = 100.0; 
  unsigned long long snum_found = 0; 
  double snear_dist = 0.;

  cwp_String Rname=NULL;                                             
  FILE *fpR=NULL;                                                       
  cwp_String rkey[9];  
  int rlocn[9];
  int rtree_numd = 0; 
  double rmult=1.0;
  cwp_String rstat=NULL; 
  int orloc  = 0;  
  double *rtree_dl[9];
  unsigned long long rtree_npoints = 0;
  node *rtree_nodes = NULL;
  double rdist = 100.0;
  double rdist2 = 100.0;
  unsigned long long rnum_found = 0;
  double rnear_dist = 0.;

  cwp_String Cname=NULL;                                             
  FILE *fpC=NULL;                                                       
  cwp_String ckey[9];  
  int clocn[9];
  int ctree_numd = 0; 
  double cmult=1.0;
  cwp_String cstat=NULL; 
  cwp_String cseq = NULL;
  cwp_String ukey = NULL;

  unsigned long long *out_elem = NULL;
  unsigned long long num_out = 0;

  unsigned long long near_elem = 0;
  int ihop = 1;

  double extent_dmin[9];
  double extent_dmax[9];
  double target[9];
  int tree_nt[9];
  double extent_min[9];
  double extent_max[9];

  cwp_String dout=NULL;
  FILE *fpD=NULL; 

  double avrad=1000.0;
  int    nedge=2;
  double cedge=100.0;
  double xc=0.;
  double yc=0.;
  double dx=0.;
  double dy=0.;
  double dr=0.;
  double ds=0.;
  double dc=0.;
  double tx=0.;
  double d4print[11];

  double xl=0.;
  double yl=0.;
  double sl=0.;
  int    kl=0;
  double xh=0.;
  double yh=0.;
  double sh=0.;
  int    kh=0;

  double **dalla=NULL;
  double **dall3=NULL;
  double *dinput = NULL;
  double *dtog = NULL;
  int klast = 2;
  int nback[2];
  int navrg[2];
  int nsmth[2];
  int ierr  = 0;
  int isize = 0;

/* hook up getpar */
  initargs(argc, argv);
  requestdoc(1);

  if(isatty(STDIN_FILENO)!=1 || isatty(STDOUT_FILENO)!=1)
    err("**** Error: this program does not input or output traces.");

  if(getpardouble("avrad",&avrad)) {
    if(avrad<0.0 && avrad!=-1.0)  err("**** Error: avrad must be positive or -1."); 
  }

  if(getparint("nedge",&nedge)) {
    if(nedge==1 || nedge<0)  err("**** Error: nedge must greater than or equal to 2 (or 0)."); 
    if(avrad==-1.0 && nedge!=0) err("**** Error: nedge cannot be specified if avrad is -1");
  }

  if(getpardouble("cedge",&cedge)) {
    if(avrad==-1.0) err("**** Error: cedge cannot be specified if avrad is -1");
    if(cedge>avrad) err("**** Error: cedge cannot be greater than avrad");
    if(cedge<0.0)  err("**** Error: cedge must be greater than or equal to 0.0"); 
  }
  else cedge = avrad / 10.;

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

  if(!getpardouble("kdist",&kdist)) kdist = 5.;
  if(kdist==0.)  err("**** Error: kdist cannot be 0.");

  if(!getpardouble("kmult",&kmult)) kmult = 2.;
  if(kmult<0.)  err("**** Error: kmult must be non-negative.");

  if(kdist<0.) {
    kdist = 0. - kdist;
    kdadd = 0;
  }

  sdist = kdist;
  sdist2 = kdist;
  rdist = kdist;
  rdist2 = kdist;

/* Set infinit extents to find the nearests.                                */
/* And set tree_nt to use all dimensions in Pythagorean nearest.            */

  for(i=0; i<9; i++) {
    extent_dmin[i] = -DBL_MAX;
    extent_dmax[i] =  DBL_MAX;
    tree_nt[i] = 1;
  }

/* Process and set the grid definition values?                         */

  getparstring("rfile", &Kname);
      
  gridcommand(&maygrid);
      
  is3d = 1;
  if(maygrid==1  && Kname != NULL) err("error: input rfile not allowed when full grid on command line.");
  if(maygrid==-1 && Kname == NULL) err("error: input rfile required when partial grid on command line.");
  if(maygrid==0  && Kname == NULL) is3d = 0;

  if (!getparint("check", &icheck)) icheck = 0;
  if (!getparint("print", &iprint)) iprint = 0;

  if(is3d==1) {

    if(countparval("cin")>0) err("error: There is a 3D Grid in rfile. So cin= file cannot be input.");

    if(maygrid!=1) { /* open if not full grid on command line (else pass fpK still NULL) */
      fpK= fopen(Kname, "r");
      if(fpK==NULL) err("error: input rfile did not open correctly.");
    }

    errwarn = 1; /* print if error or unusual thing inside gridread */
    gridread(fpK,gvals,&errwarn); 
    if(errwarn>0) err("error reading grid (from rfile or command line parameters)");

    gridset(gvals,&errwarn); 

    if(errwarn==1) err ("gridset error: grid_wb cell width must be positive.");
    else if(errwarn==2) err ("gridset error: grid_wc cell width must be positive.");
    else if(errwarn==3) err ("gridset error: corner B is within grid_wb cell width of corner A.");
    else if(errwarn>0) err ("gridset error: returned with some unrecognized error code.");
    else if(errwarn==-1) warn ("gridset warning: corner C is near A and is reset to A.");
 
    gridcheck(gvals,icheck,&errwarn); 
    if(errwarn>0) err ("gridcheck error: returned with some unrecognized error code.");

    i = 1;
    j = 1;
    gridiccdp(gvals,i,j,&j3dcdp);       /* first cell number (3d cdp number)  */
    lastab  = gvals[12] + 0.001;        /* number of cells in A-->B direction */
    lastac  = gvals[13] + 0.001;        /* number of cells in A-->C direction */
    lastcdp = lastab * lastac;          /* total number of cells              */
    if(nback[0]>=lastab) err("error: 3D Grid inlines do not have enough cells for the nback parameter."); 
    if(nback[1]>=lastac) err("error: 3D Grid crosslines do not have enough cells for the nback parameter."); 

  }
  else {
    if(countparval("cin")<1) err("error: No 3D Grid in rfile. So cin= file MUST be input.");
  }

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


  if(countparval("dout")>0) {
    getparstring("dout", &dout);
    fpD = fopen(dout,"w");
  }
  else {
    fpD = fopen("tstat.csv","w");
  }
  if (fpD == NULL) err("error: output floating datum file did not open correctly.");

/*--------------------------------------------------------------------------  */
/*- Read the parameters related to sources ---------------------------------  */
/*--------------------------------------------------------------------------  */

  if(!getparstring("sstat", &sstat)) sstat = "sstat"; /* even when sin=none  */

  stree_numd = 1;
  if(countparval("sin")>0) {
    getparstring("sin", &Sname);
    if(strcmp(Sname,"none")==0) stree_numd = 0; /* just a flag here */
  }
  else {
    Sname = ealloc1(9,1);
    strcpy(Sname,"sstat.csv");
  }

  if(stree_numd>0) { /* flag not to read these parameters if no sin file */

    stree_numd = countparval("sloc");
    if(stree_numd<1) {
      stree_numd = 2;
      skey[0] = ealloc1(2,1);
      strcpy(skey[0],"sx");
      skey[1] = ealloc1(2,1);
      strcpy(skey[1],"sy");
    }
    else {
      if(stree_numd>2) err("**** Error for sin parameters: Only 1 or 2 sloc names can be specified.");
      getparstringarray("sloc",skey);
    }

    if(!getpardouble("smult",&smult)) smult = 1.0;
    if(smult==0.0)  err("**** Error for sin parameters: smult cannot be 0.0"); 

    fpS = fopen(Sname, "r");
    if(fpS==NULL) err("error sin file did not open correctly.");
    
/* Set input numpname,pname to just store what is going to be accessed.       */
/* numpname>0 is a flag to ONLY store values if they are on pname list.       */
/* numpname<1 is a flag to NOT store values if they are on pname list.        */

    pname = ealloc1(stree_numd + 1,sizeof(cwp_String *));
    numpname = 0;
  
    for (i=0; i<stree_numd; ++i) {
      pname[numpname] = ealloc1(strlen(skey[i]),1);
      strcpy(pname[numpname],skey[i]);
      numpname++;
    }
    pname[numpname] = ealloc1(strlen(sstat),1);
    strcpy(pname[numpname],sstat);
    numpname++;

    getviaqfile(fpS, &pname, &numpname, &iztuple, numdind,   
                &ktuple, &ifixd, &SInfo, &npoints, 
                &pindepa,  &ndims, &errwarn) ;

    if(errwarn==1) err("sin file getqinfo error: extra C_SU_NAMES record in q-file");
    else if(errwarn==2) err("sin file getqinfo error: extra C_SU_NDIMS record in q-file");
    else if(errwarn==3) err("sin file getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record.");
    else if(errwarn==11) 
      err("sin file readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES");
    else if(errwarn==12) 
      err("sin file readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record.");
    else if(errwarn==22) err("sin file getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record.");
    else if(errwarn==23) err("sin file getviaqfile error: C_SU_NAMES tupled names out-of-order, changed");
    else if(errwarn==24) err("sin file getviaqfile error: C_SU_NDIMS blank where valid number expected");
    else if(errwarn==25) err("sin file getviaqfile error: C_SU_NDIMS non-number where valid number expected");
    else if(errwarn==26) err("sin file getviaqfile error: C_SU_NDIMS value must be same for all members of tuple");
    else if(errwarn==27) err("sin file getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found.");
    else if(errwarn>100) 
      err("sin file getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",errwarn-99);
    else if(errwarn>0) err("sin file getviaqfile error: unrecognized error code %d",errwarn);

    if(ifixd==0) err("sin file error: input with varying number of tuples is not allowed.");

    for (i=0; i<stree_numd; ++i) {
      slocn[i] = -1;
      for (j=0; j<iztuple; ++j) {
        if(strcmp(pname[j],skey[i])==0) slocn[i] = j;
      }
      if(slocn[i] < 0) err("sin file error: sloc name %s not found in non-tuple part of sin Q-file.",skey[i]);
    }

    osloc = -1;
    for (j=0; j<iztuple; ++j) { 
      if(strcmp(pname[j],sstat)==0) osloc = j;
    }
    if(osloc<0) err("sin file error: Q-file must have sstat parameter value %s (among non-tuple names).",sstat);

    if(npoints<1) err("sin file error: no records starting with Q found in file."); 

    stree_npoints = npoints;
    stree_nodes = ealloc1(stree_npoints,sizeof(node));

    for(i=0; i<stree_numd; i++) {
      stree_dl[i] = ealloc1double(stree_npoints);
      for(j=0; j<stree_npoints; j++) stree_dl[i][j] = SInfo[j].dlots[slocn[i]];
    }

    ihop = 1;
    connect_nodes (stree_nodes, stree_npoints, stree_dl, stree_numd, ihop);

  } /* end of  if(stree_numd>0) */

/*--------------------------------------------------------------------------  */
/*- Read the parameters related to receivers -------------------------------  */
/*--------------------------------------------------------------------------  */

  if(!getparstring("rstat", &rstat)) rstat = "gstat"; /* even when rin=none   */

  rtree_numd = 1;
  if(countparval("rin")>0) {
    getparstring("rin", &Rname);
    if(strcmp(Rname,"none")==0) rtree_numd = 0; /* just a flag here */
  }
  else {
    Rname = ealloc1(9,1);
    strcpy(Rname,"rstat.csv");
  }

  if(rtree_numd>0) { /* flag not to read these parameters if no rin file */

    rtree_numd = countparval("rloc");
    if(rtree_numd<1) {
      rtree_numd = 2;
      rkey[0] = ealloc1(2,1);
      strcpy(rkey[0],"gx");
      rkey[1] = ealloc1(2,1);
      strcpy(rkey[1],"gy");
    }
    else {
      if(rtree_numd>2) err("**** Error for rin parameters: Only 1 or 2 rloc names can be specified.");
      getparstringarray("rloc",rkey);
    }

    if(!getpardouble("rmult",&rmult)) rmult = 1.0;
    if(rmult==0.0)  err("**** Error for rin parameters: rmult cannot be 0.0"); 

    fpR = fopen(Rname, "r");
    if(fpR==NULL) err("error: rin file did not open correctly.");
    
/* Set input numpname,pname to just store what is going to be accessed.       */
/* numpname>0 is a flag to ONLY store values if they are on pname list.       */
/* numpname<1 is a flag to NOT store values if they are on pname list.        */

    pname = ealloc1(rtree_numd + 1,sizeof(cwp_String *));
    numpname = 0;
  
    for (i=0; i<rtree_numd; ++i) {
      pname[numpname] = ealloc1(strlen(rkey[i]),1);
      strcpy(pname[numpname],rkey[i]);
      numpname++;
    }
    pname[numpname] = ealloc1(strlen(rstat),1);
    strcpy(pname[numpname],rstat);
    numpname++;

    getviaqfile(fpR, &pname, &numpname, &iztuple, numdind,   
                &ktuple, &ifixd, &RInfo, &npoints, 
                &pindepa,  &ndims, &errwarn) ;

    if(errwarn==1) err("rin file getqinfo error: extra C_SU_NAMES record in q-file");
    else if(errwarn==2) err("rin file getqinfo error: extra C_SU_NDIMS record in q-file");
    else if(errwarn==3) err("rin file getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record.");
    else if(errwarn==11) 
      err("rin file readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES");
    else if(errwarn==12) 
      err("rin file readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record.");
    else if(errwarn==22) err("rin file getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record.");
    else if(errwarn==23) err("rin file getviaqfile error: C_SU_NAMES tupled names out-of-order, changed");
    else if(errwarn==24) err("rin file getviaqfile error: C_SU_NDIMS blank where valid number expected");
    else if(errwarn==25) err("rin file getviaqfile error: C_SU_NDIMS non-number where valid number expected");
    else if(errwarn==26) err("rin file getviaqfile error: C_SU_NDIMS value must be same for all members of tuple");
    else if(errwarn==27) err("rin file getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found.");
    else if(errwarn>100) 
      err("rin file getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",errwarn-99);
    else if(errwarn>0) err("rin file getviaqfile error: unrecognized error code %d",errwarn);

    if(ifixd==0) err("rin file error: input with varying number of tuples is not allowed.");

    for (i=0; i<rtree_numd; ++i) {
      rlocn[i] = -1;
      for (j=0; j<iztuple; ++j) {
        if(strcmp(pname[j],rkey[i])==0) rlocn[i] = j;
      }
      if(rlocn[i] < 0) err("rin file error: rloc %s not found in non-tuple part of rin Q-file.",rkey[i]);
    }

    orloc = -1;
    for (j=0; j<iztuple; ++j) { 
      if(strcmp(pname[j],rstat)==0) orloc = j;
    }
    if(orloc<0) err("error: rin Q-file must have rstat parameter value %s (among non-tuple names).",rstat);

    if(npoints<1) err("rin file error: no records starting with Q found in file."); 

    rtree_npoints = npoints;
    rtree_nodes = ealloc1(rtree_npoints,sizeof(node));

    for(i=0; i<rtree_numd; i++) {
      rtree_dl[i] = ealloc1double(rtree_npoints);
      for(j=0; j<rtree_npoints; j++) rtree_dl[i][j] = RInfo[j].dlots[rlocn[i]];
    }

    ihop = 1;
    connect_nodes (rtree_nodes, rtree_npoints, rtree_dl, rtree_numd, ihop);

  } /* end of  if(rtree_numd>0) */

/*--------------------------------------------------------------------------  */
/*- Read the parameters related to cdps.           -------------------------  */
/*--------------------------------------------------------------------------  */

  ctree_numd = countparval("cloc");
  if(ctree_numd<1) {
    ctree_numd = 2;
    ckey[0] = ealloc1(2,1);
    strcpy(ckey[0],"gx");
    ckey[1] = ealloc1(2,1);
    strcpy(ckey[1],"gy");
  }
  else {
    if(ctree_numd>2) err("**** Error for cin parameters: Only 1 or 2 cloc names can be specified.");
    getparstringarray("cloc",ckey);
    if(ctree_numd==1) {
      ckey[1] = ealloc1(5,1);
      strcpy(ckey[1],"numb2");
    } 
  }

  if (!getparstring("cstat", &cstat)) cstat = "tstat";

  if(is3d!=1) {
    getparstring("cin", &Cname);

    fpC = fopen(Cname, "r");
    if(fpC==NULL) err("error: cin file did not open correctly.");
    
    if(countparval("cseq") > 0) {
      getparstring("cseq", &cseq);
    }
    else {
      cseq = ealloc1(3,1);
      strcpy(cseq,"cdp");
    }

/* Set input numpname,pname to just store what is going to be accessed.       */
/* numpname>0 is a flag to ONLY store values if they are on pname list.       */
/* numpname<1 is a flag to NOT store values if they are on pname list.        */

    pname = ealloc1(ctree_numd + 2,sizeof(cwp_String *));
    numpname = 0;
  
    for (i=0; i<ctree_numd; ++i) {
      pname[numpname] = ealloc1(strlen(ckey[i]),1);
      strcpy(pname[numpname],ckey[i]);
      numpname++;
    }
    pname[numpname] = ealloc1(strlen(cstat),1);
    strcpy(pname[numpname],cstat);
    numpname++;

    if(strcmp(cseq,"asis") != 0) {
      pname[numpname] = ealloc1(strlen(cseq),1);
      strcpy(pname[numpname],cseq);
      numpname++;
    }

    getviaqfile(fpC, &pname, &numpname, &iztuple, numdind,   
                &ktuple, &ifixd, &CInfo, &npoints, 
                &pindepa,  &ndims, &errwarn) ;

    if(errwarn==1) err("cin file getqinfo error: extra C_SU_NAMES record in q-file");
    else if(errwarn==2) err("cin file getqinfo error: extra C_SU_NDIMS record in q-file");
    else if(errwarn==3) err("cin file getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record.");
    else if(errwarn==11) 
      err("cin file readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES");
    else if(errwarn==12) 
      err("cin file readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record.");
    else if(errwarn==22) err("cin file getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record.");
    else if(errwarn==23) err("cin file getviaqfile error: C_SU_NAMES tupled names out-of-order, changed");
    else if(errwarn==24) err("cin file getviaqfile error: C_SU_NDIMS blank where valid number expected");
    else if(errwarn==25) err("cin file getviaqfile error: C_SU_NDIMS non-number where valid number expected");
    else if(errwarn==26) err("cin file getviaqfile error: C_SU_NDIMS value must be same for all members of tuple");
    else if(errwarn==27) err("cin file getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found.");
    else if(errwarn>100) 
      err("cin file getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",errwarn-99);
    else if(errwarn>0) err("cin file getviaqfile error: unrecognized error code %d",errwarn);

    if(ifixd==0) err("cin file error: input with varying number of tuples is not allowed.");

    for (i=0; i<ctree_numd; ++i) {
      clocn[i] = -1;
      for (j=0; j<iztuple; ++j) {
        if(strcmp(pname[j],ckey[i])==0) clocn[i] = j;
      }
      if(clocn[i] < 0) err("cin file error: cloc %s not found in non-tuple part of cin Q-file.",ckey[i]);
    }

    if(npoints<1) err("cin file error: no records starting with Q found in file."); 

    lastcdp = npoints; 

    if(nback[0]>=lastcdp) err("error: cin file does not have enough Q-records for the nback parameter."); 

    if(strcmp(cseq,"asis") == 0) {
      cseq = ealloc1(5,1);
      strcpy(cseq,"numb1");
    }
    else {
      locp = -1;
      for (i=0; i<iztuple; ++i) {
        if(strcmp(pname[i],cseq)==0) locp = i;
      }
      ukey = ealloc1(1,strlen(cseq)+2);
      for (j=0; j<strlen(cseq); ++j) ukey[j+1] = cseq[j];
      ukey[0] = '_';
      ukey[strlen(cseq)+1] = '_';
      for (i=0; i<iztuple; ++i) {
        if(strstr(pname[i],ukey)!=0) locp = i;
      }
      if(locp<0) err("error: input cin Q-file must have your cseq=%s (among non-tuple names).",cseq);

      qsort(CInfo,lastcdp,sizeof(struct QInfo),compSort1);
      for (j=1; j<lastcdp; j++) {
        if(CInfo[j-1].dlots[locp] == CInfo[j].dlots[locp]) {
             err("error: Two sets of values input for cin cseq = %f",CInfo[j].dlots[locp]);
        }
      }
    }

  } /* end of  if(is3d!=1) */

  if(stree_numd!=0 && stree_numd!=ctree_numd) err("**** Error: Must have the same number of sloc as cloc names.");
  if(rtree_numd!=0 && rtree_numd!=ctree_numd) err("**** Error: Must have the same number of rloc as cloc names.");

/* The easiest way to deal with just 1 location value is to pretend there     */
/* are 2 for some parts of the code. But note ctree_numd is not reset to 2    */

  if(ctree_numd==1) {
    clocn[1] = clocn[0];
    slocn[1] = slocn[0];
    rlocn[1] = rlocn[0];
    avrad *= sqrt(2.);
  }

/* Allocate various memory buffers.                                           */

  if(stree_npoints>rtree_npoints) 
       out_elem = (unsigned long long *)ealloc1(sizeof(unsigned long long),stree_npoints);
  else out_elem = (unsigned long long *)ealloc1(sizeof(unsigned long long),rtree_npoints);

  isize = 2*navrg[0] + 3; /* plus 3 because navrg might be -1 */
  if(is3d==1 && 2*navrg[1] + 3 > isize) isize = 2*navrg[1] + 3;

  if(2*nsmth[0] + 1 > isize) isize = 2*nsmth[0] + 1;
  if(is3d==1 && 2*nsmth[1] + 1 > isize) isize = 2*nsmth[1] + 1;

  dalla = calloc(lastcdp + isize,sizeof(double *));
  if(dalla == NULL) err("**** Unable to allocate dalla pp memory ");
  dinput = calloc((lastcdp + isize)*klast,sizeof(double));
  if(dinput == NULL) err("**** Unable to allocate dalla p memory ");

  dtog = ealloc1double(klast);

/* Note that klast is always 2 in initial version of program (Nov 2023).      */
/* This is why you find accesses such as dalla[j][0] and dalla[j][1] later.   */
/* Leaving the klast as a variable matches code in other programs...          */
/* ...and may even change in later versions of this program.                  */
/* If you are worried that the for-loops with klast are slower remember that  */
/* the compiler optimizer will un-roll these loops anyway.                    */

  for (j=0; j<lastcdp + isize; j++) {
    dalla[j] = dinput + j*klast;
    for(i=0; i<klast; i++) dalla[j][i] = 0.0;
  }

  if(is3d==1) {
    if(lastab>lastac) n = lastab;
    else n = lastac;
    dall3 = calloc(n + isize,sizeof(double *));
    if(dall3 == NULL) err("**** Unable to allocate dall3 pp memory ");
    dinput = calloc((n + isize)*klast,sizeof(double));
    if(dinput == NULL) err("**** Unable to allocate dalla p memory ");

    for (j=0; j<n + isize; j++) {
      dall3[j] = dinput + j*klast; 
      for(i=0; i<klast; i++) dall3[j][i] = 0.0;
    }
  }

/*--------------------------------------------------------------------------  */
  checkpars(); 
/*--------------------------------------------------------------------------  */

  for (j=0; j<lastcdp; j++) { 

    if(is3d==1) {
      gridcdpic(gvals,j3dcdp+j,&igi,&igc);  
      gridicrawxy(gvals,igi,igc,target,target+1);
    }
    else {
      target[0] = CInfo[j].dlots[clocn[0]];
      target[1] = CInfo[j].dlots[clocn[1]];
    }

    if(stree_numd>0) {
/* Find the distance to the nearest shot.  */
      if(kdadd==1) sdist2 = sdist + snear_dist;
      cycle_for_near(stree_nodes,stree_dl,tree_nt,stree_numd,
                     extent_dmin, extent_dmax, target, sdist2, kmult,
                     &near_elem, &snear_dist, &snum_found, &kcycles);

      snear_dist = sqrt(snear_dist);

      if(avrad<0.0) {
        dalla[j][0] = SInfo[near_elem].dlots[osloc];
      }
      else { 
/* Define a square centred on cdp using distance to nearest shot PLUS avrad.  */
        extent_min[0] = target[0] - snear_dist - avrad;
        extent_max[0] = target[0] + snear_dist + avrad;
        extent_min[1] = target[1] - snear_dist - avrad;
        extent_max[1] = target[1] + snear_dist + avrad;
/* Find all shots within the square. */
        find_in (stree_nodes,stree_dl,stree_numd,extent_min,extent_max,out_elem,&num_out);

        dalla[j][0] = 0.0; 
        xc = 0.;
        yc = 0.;
        k = 0; 
        for (i=0; i<num_out; i++) { /* Find all shots within enclosed circle.*/
          if((target[0]-SInfo[out_elem[i]].dlots[slocn[0]])*(target[0]-SInfo[out_elem[i]].dlots[slocn[0]]) +
             (target[1]-SInfo[out_elem[i]].dlots[slocn[1]])*(target[1]-SInfo[out_elem[i]].dlots[slocn[1]]) 
              <= (snear_dist+avrad)*(snear_dist+avrad)) {
            dalla[j][0] += SInfo[out_elem[i]].dlots[osloc];
            xc += SInfo[out_elem[i]].dlots[slocn[0]];
            yc += SInfo[out_elem[i]].dlots[slocn[1]];
            out_elem[k] = out_elem[i]; /* Reduce to just the shots in circle */
            k++;
          }
        }

        if(k>1) {
          dalla[j][0] /= k; /* Set to average incase of no edge compensation  */
          xc /= k;
          yc /= k;
        }

        dx = target[0] - xc; 
        dy = target[1] - yc;
        dr = sqrt(dx*dx+dy*dy);

        if(k>=nedge && dr>=cedge) { /* Edge compensation ? */

          ds = dy / dr; /* set sine   */
          dc = dx / dr; /* set cosine */

          xl = 0.;
          yl = 0.;
          sl = 0.;
          kl = 0;
          xh = 0.;
          yh = 0.;
          sh = 0.;
          kh = 0;

          for (i=0; i<k; i++) { 

            tx =  (xc - SInfo[out_elem[i]].dlots[slocn[0]])*dc 
               +  (yc - SInfo[out_elem[i]].dlots[slocn[1]])*ds;
/*          ty = -(xc - SInfo[out_elem[i]].dlots[slocn[0]])*ds   crossline    */
/*             +  (yc - SInfo[out_elem[i]].dlots[slocn[1]])*dc;  unneeded     */

            if(tx<0.) { /* which side of centroid is the point located?       */
              xl += SInfo[out_elem[i]].dlots[slocn[0]];
              yl += SInfo[out_elem[i]].dlots[slocn[1]];
              sl += SInfo[out_elem[i]].dlots[osloc];
              kl++;
            }
            else {
              xh += SInfo[out_elem[i]].dlots[slocn[0]];
              yh += SInfo[out_elem[i]].dlots[slocn[1]];
              sh += SInfo[out_elem[i]].dlots[osloc];
              kh++;
            }
   
          } /* end of  for (i=0; i<k; i++) { */

          if(kl>0 && kh>0 && (xh!=xl || yh!=yl)) {

            xl /= kl;
            yl /= kl;
            sl /= kl;
            xh /= kh;
            yh /= kh;
            sh /= kh;

            dx = xh - xl; 
            dy = yh - yl;
            dr = sqrt(dx*dx+dy*dy);
            ds = dy / dr; /* set sine   */
            dc = dx / dr; /* set cosine */
            tx =  (xh - target[0])*dc + (yh - target[1])*ds;
/*          ty = -(xh - target[0])*ds + (yh - target[1])*dc;   crossline unneeded */

            dalla[j][0] = tx/dr * sl - (tx/dr-1.0) * sh;

          } /* end of if(kl>0 && kh>0) { */
        } /* end of if(k>=nedge && dr>=cedge) {  */
      } /* end of  if(avrad<0.0) ... else { */

    } /* end of    if(stree_numd>0) */

    if(rtree_numd>0) {

      if(kdadd==1) rdist2 = rdist + rnear_dist;
      cycle_for_near(rtree_nodes,rtree_dl,tree_nt,rtree_numd,
                     extent_dmin, extent_dmax, target, rdist2, kmult,
                     &near_elem, &rnear_dist, &rnum_found, &kcycles);

      rnear_dist = sqrt(rnear_dist);

      if(avrad<0.0) {
        dalla[j][1] = RInfo[near_elem].dlots[orloc];
      }
      else {
        extent_min[0] = target[0] - rnear_dist - avrad;
        extent_max[0] = target[0] + rnear_dist + avrad;
        extent_min[1] = target[1] - rnear_dist - avrad;
        extent_max[1] = target[1] + rnear_dist + avrad;

        find_in (rtree_nodes,rtree_dl,rtree_numd,extent_min,extent_max,out_elem,&num_out);

        dalla[j][1] = 0.0; 
        xc = 0.;
        yc = 0.;
        k = 0; 
        for (i=0; i<num_out; i++) { 
          if((target[0]-RInfo[out_elem[i]].dlots[rlocn[0]])*(target[0]-RInfo[out_elem[i]].dlots[rlocn[0]]) +
             (target[1]-RInfo[out_elem[i]].dlots[rlocn[1]])*(target[1]-RInfo[out_elem[i]].dlots[rlocn[1]])
              <= (rnear_dist+avrad)*(rnear_dist+avrad)) {
            dalla[j][1] += RInfo[out_elem[i]].dlots[orloc];
            xc += RInfo[out_elem[i]].dlots[rlocn[0]];
            yc += RInfo[out_elem[i]].dlots[rlocn[1]];
            out_elem[k] = out_elem[i];
            k++;
          }
        }

        if(k>1) {
          dalla[j][1] /= k;
          xc /= k;
          yc /= k;
        }

        dx = target[0] - xc; 
        dy = target[1] - yc;
        dr = sqrt(dx*dx+dy*dy);

        if(k>=nedge && dr>=cedge) { /* Edge compensation ? */

          ds = dy / dr; 
          dc = dx / dr; 

          xl = 0.;
          yl = 0.;
          sl = 0.;
          kl = 0;
          xh = 0.;
          yh = 0.;
          sh = 0.;
          kh = 0;

          for (i=0; i<k; i++) { 

            tx =  (xc - RInfo[out_elem[i]].dlots[rlocn[0]])*dc 
               +  (yc - RInfo[out_elem[i]].dlots[rlocn[1]])*ds;

            if(tx<0.) { 
              xl += RInfo[out_elem[i]].dlots[rlocn[0]];
              yl += RInfo[out_elem[i]].dlots[rlocn[1]];
              sl += RInfo[out_elem[i]].dlots[orloc];
              kl++;
            }
            else {
              xh += RInfo[out_elem[i]].dlots[rlocn[0]];
              yh += RInfo[out_elem[i]].dlots[rlocn[1]];
              sh += RInfo[out_elem[i]].dlots[orloc];
              kh++;
            }
   
          } /* end of  for (i=0; i<k; i++) { */

          if(kl>0 && kh>0 && (xh!=xl || yh!=yl)) {
            xl /= kl;
            yl /= kl;
            sl /= kl;
            xh /= kh;
            yh /= kh;
            sh /= kh;
            dx = xh - xl; 
            dy = yh - yl;
            dr = sqrt(dx*dx+dy*dy);
            ds = dy / dr; 
            dc = dx / dr; 
            tx =  (xh - target[0])*dc + (yh - target[1])*ds;

            dalla[j][1] = tx/dr * sl - (tx/dr-1.0) * sh;

          } /* end of if(kl>0 && kh>0) { */
        } /* end of if(k>=nedge && dr>=cedge) {  */
      } /* end of  if(avrad<0.0) ... else { */

    } /* end of    if(rtree_numd>0) */

  } /* end of  for (j=0; j<lastcdp; j++) { */

/*--------------------------------------------------------------------------  */
/* perform smoothing? */

  if(is3d!=1) { /* is 2d */

    if(navrg[0]>0) {
      runsmo (dalla,lastcdp+2*navrg[0]+1,klast,nback[0],navrg[0],dtog,&ierr);
      if(ierr==1) err("**** Error in averaging. Less than 5 Q-records in cin file and navrg>0");
      else if(ierr>0) err("**** Error in navrg averaging for cin file. Parameter out-of-range.");
    }
    else if(navrg[0]==-1) {
      for(i=0; i<klast; i++) dtog[i] = 0.0;
      for (j=0; j<lastcdp; j++) {
        for(i=0; i<klast; i++) dtog[i] += dalla[j][i];
      }
      for (j=0; j<lastcdp; j++) {
        for(i=0; i<klast; i++) dalla[j][i] = dtog[i] / lastcdp;
      }
    } 

    if(nsmth[0]>0) {
      runsmo (dalla,lastcdp+2*nsmth[0]+1, klast,nback[0],0-nsmth[0],dtog,&ierr);
      if(ierr==1) err("**** Error in smoothing. Less than 5 Q-records in cin file and nsmth>0");
      else if(ierr>0) err("**** Error in nsmth smoothing for cin file. Parameter out-of-range.");
    }  

  }
  else { /* is 3d */

    if(navrg[0]!=0 || nsmth[0]>0) {

      for (n=0; n<lastac; n++) {

        for (j=0; j<lastab; j++) {
          for(i=0; i<klast; i++) dall3[j][i] = dalla[j+n*lastab][i];
        }

        if(navrg[0]>0) {
          runsmo(dall3,lastab+2*navrg[0]+1,klast,nback[0],navrg[0],dtog,&ierr);
          if(ierr==1) err("**** Error in averaging. Less than 5 inline cells defined in 3D Grid and navrg>0");
          else if(ierr>0) err("**** Error in navrg averaging for 3D Grid inlines. Parameter out-of-range.");
        }
        else if(navrg[0]==-1) {
          for(i=0; i<klast; i++) dtog[i] = 0.0;
          for (j=0; j<lastab; j++) {
            for(i=0; i<klast; i++) dtog[i] += dall3[j][i];
          }
          for (j=0; j<lastab; j++) {
            for(i=0; i<klast; i++) dall3[j][i] = dtog[i] / lastab;
          }
        }  

        if(nsmth[0]>0) {
          runsmo(dall3,lastab+2*nsmth[0]+1,klast,nback[0],0-nsmth[0],dtog,&ierr);
          if(ierr==1) err("**** Error in smoothing. Less than 5 inline cells defined in 3D Grid and nsmth>0.");
          else if(ierr>0) err("**** Error in nsmth smooth for 3D Grid inlines. Parameter out-of-range.");
        }

        for (j=0; j<lastab; j++) {
          for(i=0; i<klast; i++) dalla[j+n*lastab][i] = dall3[j][i];
        }

      } /* end of  for (n=0; n<lastac; n++)  */

    } /* end of  if(navrg[0]!=0 || nsmth[0]>0) */

/* Repeat in the other direction. */

    if(navrg[1]!=0 || nsmth[1]>0) {

      for (n=0; n<lastab; n++) {

        for (j=0; j<lastac; j++) {
          for(i=0; i<klast; i++) dall3[j][i] = dalla[n+j*lastab][i];
        }

        if(navrg[1]>0) {
          runsmo(dall3,lastac+2*navrg[1]+1,klast,nback[1],navrg[1],dtog,&ierr);
          if(ierr==1) err("**** Error averaging. Less than 5 crossline cells defined in 3D Grid and navrg>0");
          else if(ierr>0) err("**** Error in navrg average for 3D Grid crosslines. Parameter out-of-range.");
        }

        if(nsmth[1]>0) {
          runsmo(dall3,lastac+2*nsmth[1]+1,klast,nback[1],0-nsmth[1],dtog,&ierr);
          if(ierr==1) err("**** Error smoothing. Less than 5 crossline cells defined in 3D Grid and nsmth>0");
          else if(ierr>0) err("**** Error in nsmth smooth for 3D Grid crosslines. Parameter out-of-range.");
        }
        else if(navrg[1]==-1) {
          for(i=0; i<klast; i++) dtog[i] = 0.0;
          for (j=0; j<lastac; j++) {
            for(i=0; i<klast; i++) dtog[i] += dall3[j][i];
          }
          for (j=0; j<lastac; j++) {
            for(i=0; i<klast; i++) dall3[j][i] = dtog[i] / lastac;
          }
        }  

        for (j=0; j<lastac; j++) {
          for(i=0; i<klast; i++) dalla[n+j*lastab][i] = dall3[j][i];
        }

      } /* end of  for (n=0; n<lastab; n++)  */

    } /* end of  if(navrg[1]!=0 || nsmth[1]>0) */

  } /* end of  if(is3d!=1) ... else */ 

/*--------------------------------------------------------------------------  */
  if(stree_numd<1) smult = 0.0;
  if(rtree_numd<1) rmult = 0.0;
  if(stree_numd>0 && rtree_numd>0) cmult = 0.5;

  if(is3d!=1) { 

    fprintf(fpD,"C_SU_MATCH,%s,,,,,\nC_SU_SETID,Q,,,,,\nC_SU_FORMS\nC_SU_ID",cseq);
    for(i=0; i<6; i++) fprintf(fpD,",%s",formxy);
    fprintf(fpD,"\nC_SU_NAMES,,,,,,\nC_SU_ID,%s,%s,%s,%s,%s,%s\n",cseq,ckey[0],ckey[1],sstat,rstat,cstat);

    for (j=0; j<lastcdp; j++) { 

      if(locp>-1) d4print[0] = CInfo[j].dlots[locp];
      else d4print[0] = j+1;

      d4print[1] = CInfo[j].dlots[clocn[0]];

      if(ctree_numd==2) d4print[2] = CInfo[j].dlots[clocn[1]];
      else d4print[2] = 0.0;

      d4print[3] = dalla[j][0]*smult;
      d4print[4] = dalla[j][1]*rmult;
      d4print[5] = (dalla[j][0]*smult + dalla[j][1]*rmult) * cmult;

      fprintf(fpD,"Q");
      for(i=0; i<6; i++) fprintf(fpD,formxylong,d4print[i]);
      fprintf(fpD,"\n");

    }

  }
  else { /* is a 3D */

/* Writing the same number of commas in all records (just because I want to). */
/* Note that the names and values correspond to bintype=-30 in subincsv.      */
/* Specifically: sx,sy are cell-centre real-world coordinates.                */
/*               gx,gy are cell-centre grid-rotated-and-shifted coordinates.  */

    fprintf(fpD,"C_SU_MATCH,cdp,,,,,,\nC_SU_SETID,Q,,,,,,\nC_SU_FORMS\nC_SU_ID");
    for(i=0; i<10; i++) fprintf(fpD,",%s",formxy);
    fprintf(fpD,"\nC_SU_NAMES,,,,,,,,,,,\nC_SU_ID,cdp,cdpt,igi,igc,sx,sy,gx,gy,%s,%s,%s\n",sstat,rstat,cstat);

    for (j=0; j<lastcdp; j++) { 

      gridcdpic(gvals,j3dcdp+j,&igi,&igc);  
      gridiccdp90(gvals,igi,igc,&n);
      gridicrawxy(gvals,igi,igc,&xc,&yc);
      gridicgridxy(gvals,igi,igc,&dx,&dy);

      d4print[0] = j3dcdp+j;
      d4print[1] = n;
      d4print[2] = igi;
      d4print[3] = igc;
      d4print[4] = xc;
      d4print[5] = yc;
      d4print[6] = dx;
      d4print[7] = dy;
      d4print[8] = dalla[j][0]*smult;
      d4print[9] = dalla[j][1]*rmult;
      d4print[10]= (dalla[j][0]*smult + dalla[j][1]*rmult) * cmult;

      fprintf(fpD,"Q");
      for(i=0; i<11; i++) fprintf(fpD,formxylong,d4print[i]);
      fprintf(fpD,"\n");

    }

  } /* end of  if(is3d!=1) ...else */

  return 0;

} /* end of main for sushiftcsv */

/* --------------------------------------------------------------------------------------------------- */

void connect_nodes (node *tree_nodes, unsigned long long tree_numc, 
                    double **tree_dl, int tree_numd, int ihop) {

/*          This function connects the already allocated tree nodes.                                   */     
/*                                                                                                     */     
/* Input arguments:                                                                                    */     
/* tree_nodes    The already fully allocated tree nodes.                                               */     
/* tree_numc     Number of nodes in tree_nodes.                                                        */     
/* tree_dl       Pointers to first element of values of each dimension. That is, each coordinate is in */
/*               a separate contiguous array. These are the pointers to the start of those arrays.     */     
/*               Those arrays have size tree_numc.                                                     */     
/* tree_numd     Number of pointers in tree_dl (i.e. number of dimensions).                            */     
/* ihop      = 1 means process the points in dispersed order (approximately random).                   */     
/*               The worst search times for kd trees occur when there are just a few long branches.    */     
/*               This happens in simple kd tree code if the points are ordered in increasing or        */     
/*               decreasing coordinates. Unfortunately, that is often the case for survey files.       */     
/*               The literature about kd trees explains how to create trees with branches of the same  */     
/*               depth (perfectly balanced trees). That is not done here. Instead, this option hops    */     
/*               through the points and loads them in more-or-less random order. This creates a        */     
/*               reasonably-balanced tree (unless you are very unlucky).                               */     
/*               Note that seriously unbalanced trees still work, but they may be very slow.           */     
/*           = 0 process the points in the order they exist within their arrays.                       */     

  unsigned long long nrat = 0;
  unsigned long long nd = 0;
  unsigned long long n = 0;
  unsigned long long m = 0;

  for(m=0; m<tree_numc; m++) {
    tree_nodes[m].l    = 0; 
    tree_nodes[m].r    = 0;  
    tree_nodes[m].elem = m;
  }

/* nrat,nd and the for-loops are just heuristically set (I made them up).         */
/* The basic concept is: Survey points are usually in order. So, we want to load  */
/* them in big hops so that the top nodes will have far-apart coordinates. Then   */
/* reduce hop size so next deeper nodes will also be dispersed nicely. And so on. */

  if(ihop==1) {
    m = 0;
    for (nrat=tree_numc; nrat>0; nrat=nrat*0.6) {
      if(nrat>6) nd = sqrt((double) nrat); 
      else nd = 0;
      for (n=nd+nrat/2; n<tree_numc; n=n+nd+nrat) {
        if (tree_nodes[n].l == 0) {     /* pointer value just used here as a flag */
          tree_nodes[n].l = tree_nodes; /* (so we know we already did this point) */
          tree_nodes[m].elem = n;
          m++;
        }
      }  
    }  
    for(m=0; m<tree_numc; m++) tree_nodes[m].l = 0; /* un-flag it */  
  }

  connect_all(tree_nodes,tree_numc,tree_dl,tree_numd);

  return;

}

/* --------------------------------------------------------------------------------------------------- */

void connect_all(node *tree_nodes, unsigned long long tree_numc, double **tree_dl, int tree_numd) {

/*          This function connects the already allocated tree nodes in their order in elem.            */     
/*                                                                                                     */     
/* Input arguments:                                                                                    */     
/* tree_nodes    The already fully allocated tree nodes.                                               */     
/* tree_numc     Number of nodes in tree_nodes.                                                        */     
/* tree_dl       Pointers to first element of values of each dimension. That is, each coordinate is in */     
/*               a separate contiguous array. These are the pointers to the start of those arrays.     */     
/*               Those arrays have size tree_numc.                                                     */     
/* tree_numd     Number of pointers in tree_dl (i.e. number of dimensions).                            */     
/*                                                                                                     */     
/* Note that, technically, this function processes the nodes in their sequential order in tree_nodes.  */
/* But the coordinates of each node are pointed to via node->elem, which means the points              */
/* are actually loaded into the tree in a different order than the order in their arrays.              */

  unsigned long long m = 1; 
  int naxe = 0;

  for(m=1; m<tree_numc; m++) { /* First m is 1 since tree_nodes[0] is set correct.  */ 

    node *now_node, *next_node, *m_node;

    now_node  = NULL; /* just to prevent compiler warning of possible uninitialed use.*/
    m_node    = tree_nodes + m; /* the m-th node to connect                           */ 
    next_node = tree_nodes;     /* walk down tree from top for each node to connect   */

    naxe = 0;
    while(next_node!=0) {
      now_node = next_node;

      naxe++;
      if(naxe==tree_numd) naxe = 0;

      if(tree_dl[naxe][m_node->elem] < tree_dl[naxe][now_node->elem]) 
           next_node = now_node->l;
      else next_node = now_node->r;
    }

    if(tree_dl[naxe][m_node->elem] < tree_dl[naxe][now_node->elem]) 
         now_node->l = m_node; 
    else now_node->r = m_node;

  }

  return;
}

/* --------------------------------------------------------------------------------------------------- */

void find_in (node *tree_nodes, double **tree_dl, int tree_numd,
              double *extent_min, double *extent_max,  
              unsigned long long *out_elem, unsigned long long *num_out) {

/*          This function finds all points between specified extents of the dimensions.                */     
/*                                                                                                     */     
/* Input arguments:                                                                                    */     
/* tree_nodes    The already fully allocated and connected tree nodes.                                 */     
/* tree_dl       Pointers to first element of values of each dimension. That is, each coordinate is in */     
/*               a separate contiguous array. These are the pointers to the start of those arrays.     */     
/*               Those arrays have size tree_numc.                                                     */     
/* tree_numd     Number of pointers in tree_dl (i.e. number of dimensions).                            */     
/* extent_min    Minimum value of this dimension to find points. Size tree_numd.                       */     
/*               Greater OR EQUAL to this value is in range.                                           */     
/* extent_max    Maximum value of this dimension to find nearest point. Size tree_numd.                */     
/*               Strictly LESS than this value is in range.                                            */     
/*                                                                                                     */     
/* Output arguments:                                                                                   */     
/* out_elem      The element numbers of the points in the tree_dl arrays. With m meaning the           */     
/*               coordinates of the point are tree_dl[n][m] where n is the dimension number.           */     
/* num_out       Is the number of points found within the extent ranges.                               */
/*                                                                                                     */     
/* Return:false  means some kind of input argument error.                                              */     
/*               NOTE: inputting an impossible extent_min,extent_max range is not an error, you        */     
/*                     just get 0 for num_out.                                                         */     

  int naxe = 1;
  node *now_node;

  *num_out = 0;
  now_node = tree_nodes;
  find_in_rcur(tree_dl, tree_numd, now_node, naxe, extent_min, extent_max, out_elem, num_out);

  return;
}

/* --------------------------------------------------------------------------------------------------- */

void find_in_rcur (double **tree_dl, int tree_numd,
                   node * now_node, int naxe, 
                   double *extent_min, double *extent_max, 
                   unsigned long long *out_elem, unsigned long long *num_out) {

/*          This function finds all points between specified extents of the dimensions.                */     
/*                                                                                                     */     
/* Input arguments:                                                                                    */     
/* tree_dl       Pointers to first element of values of each dimension. That is, each coordinate is in */     
/*               a separate contiguous array. These are the pointers to the start of those arrays.     */     
/*               Those arrays have size tree_numc.                                                     */     
/* tree_numd     Number of pointers in tree_dl (i.e. number of dimensions).                            */     
/* extent_min    Minimum value of this dimension to find points. Size tree_numd.                       */     
/*               Greater OR EQUAL to this value is in range.                                           */     
/* extent_max    Maximum value of this dimension to find nearest point. Size tree_numd.                */     
/*               Strictly LESS than this value is in range.                                            */     
/*                                                                                                     */     
/* Input/Output arguments:                                                                             */     
/* now_node      The current node to start from.                                                       */     
/* naxe          The current splitting dimension.                                                      */     
/* out_elem      Accumulating element numbers of the points in the tree_dl arrays. With m meaning the  */     
/*               coordinates of the point are tree_dl[n][m] where n is the dimension number.           */     
/* num_out       Accumulating number of points found within the extent ranges.                         */

  if(now_node==0) return;

  bool in = true;                
  int i = 0;

  for(i=0; i<tree_numd; i++) {
    if(tree_dl[i][now_node->elem] < extent_min[i] || tree_dl[i][now_node->elem] >= extent_max[i]) { 
      in = false;
      break;
    }
  }
  if(in) { 
    out_elem[*num_out] = now_node->elem;  
    *num_out = *num_out + 1;
  }

/* Find more points. */

  if(naxe==tree_numd) naxe = 0;

  if(tree_dl[naxe][now_node->elem] >= extent_min[naxe] && now_node->l)
    find_in_rcur(tree_dl, tree_numd, now_node->l, naxe+1, extent_min, extent_max, out_elem, num_out);

  if(tree_dl[naxe][now_node->elem] <  extent_max[naxe] && now_node->r)
    find_in_rcur(tree_dl, tree_numd, now_node->r, naxe+1, extent_min, extent_max, out_elem, num_out);

  return;
}

/* --------------------------------------------------------------------------------------------------- */

void find_near_rcur (double **tree_dl, int *tree_nt, int tree_numd, node * now_node, int naxe, 
                     double *extent_min, double *extent_max, double *target, 
                     unsigned long long *near_elem, double *near_dist, unsigned long long *num_found) {

/*          This function finds a nearest point between specified extents of the dimensions.           */     
/*                                                                                                     */     
/* Input arguments:                                                                                    */     
/* tree_dl       Pointers to first element of values of each dimension. That is, each coordinate is in */     
/*               a separate contiguous array. These are the pointers to the start of those arrays.     */     
/*               Those arrays have size tree_numc.                                                     */     
/* tree_nt       Near type flags. 1 means standard pythagorean nearest. See note for what 0 means.     */     
/* tree_numd     Number of pointers in tree_dl and flags in tree_nt (i.e. number of dimensions).       */     
/* extent_min    Minimum value of this dimension to find nearest point. Size tree_numd.                */     
/*               Greater OR EQUAL to this value is in range.                                           */     
/* extent_max    Maximum value of this dimension to find nearest point. Size tree_numd.                */     
/*               Strictly LESS than this value is in range.                                            */     
/* target        Find the point nearest here considering extents and tree_nt. Size tree_numd.          */     
/*                                                                                                     */     
/* Input/Output arguments:                                                                             */     
/* now_node      The current node to start from.                                                       */     
/* naxe          The current splitting dimension.                                                      */     
/*                                                                                                     */     
/* Output arguments:                                                                                   */     
/* near_elem     The element number of a nearest point in the tree_dl arrays. For instance, a value of */     
/*               m means the coordinates of point are tree_dl[n][m] where n is the dimension number.   */     
/* near_dist     The SQUARED distance between the target and nearest point (if num_found>0).           */
/* num_found     0 if no point is found within the extent ranges.                                      */
/*              >0 is as many equally-near points as are found within the extent ranges.               */
/*                 The returned near_elem is the highest-numbered elem of all the nearest points.      */
/*                                                                                                     */     
/*    NOTE: A tree_nt value of 0 means the difference between the point coordinate and the target      */     
/*          coordinate are not added to the Pythagorean sum for the specified dimension(s).            */     
/*          So, the point distances are determined as if that dimension was not specified at all.      */     
/*          However, the extent_min,extent_max values for that dimension are still used.               */     
/*          The result is therefor the nearest point considering only the other dimensions, but        */     
/*          still restricted by the extent range of that dimension(s).                                 */     

  if(now_node==0) return; 

  bool in = true;          
  int i = 0;
  double rad = 0.;
 
  for(i=0; i<tree_numd; i++) {
    if(tree_dl[i][now_node->elem] < extent_min[i] || tree_dl[i][now_node->elem] >= extent_max[i]) { 
      in = false;
      break;
    }
  }

  if(in) { 
    rad = 0.;
    for(i=0; i<tree_numd; i++) {
      if(tree_nt[i] != 0) {
        rad += (target[i]-tree_dl[i][now_node->elem]) * (target[i]-tree_dl[i][now_node->elem]);
      }
    }

/* The following set of ifs could be done differently. But I think coding it  */
/* this way saves cpu time by immediately rejecting the > cases using 1 if.   */
/* But who knows how modern optimizers will alter this.                       */

    if(rad <= *near_dist) {
      if(rad < *near_dist) { /* if distance is smaller, reset count to 1 */
        *num_found = 1;
        *near_elem = now_node->elem;
        *near_dist = rad; 
      }
      else { /* so, distances are equal. Increment count, set output to higher elem.*/
        *num_found = *num_found + 1;
        if(now_node->elem > *near_elem) *near_elem = now_node->elem;
      }
    }

  }

/* Find more. */

  if(naxe==tree_numd) naxe = 0;

  if(tree_dl[naxe][now_node->elem] >= extent_min[naxe] && now_node->l)
    find_near_rcur(tree_dl, tree_nt, tree_numd, now_node->l, naxe+1,
                   extent_min, extent_max, target, 
                   near_elem, near_dist, num_found);

  if(tree_dl[naxe][now_node->elem] <  extent_max[naxe] && now_node->r)
    find_near_rcur(tree_dl, tree_nt, tree_numd, now_node->r, naxe+1, 
                   extent_min, extent_max, target, 
                   near_elem, near_dist, num_found);

  return;
}

/* --------------------------------------------------------------------------------------------------- */

void cycle_for_near(node *tree_nodes, double **tree_dl, int *tree_nt, int tree_numd,
                    double *extent_min, double *extent_max, double *target, 
                    double init_rad, double rad_scal,  
                    unsigned long long *near_elem, double *near_dist,
                    unsigned long long *num_found, int *ncycles) {

/*          This function finds a nearest point between specified extents of the dimensions.           */     
/*                                                                                                     */     
/* Input arguments:                                                                                    */     
/* tree_nodes    The already fully allocated and connected tree nodes.                                 */     
/* tree_dl       Pointers to first element of values of each dimension. That is, each coordinate is in */     
/*               a separate contiguous array. These are the pointers to the start of those arrays.     */     
/*               Those arrays have size tree_numc.                                                     */     
/* tree_nt       Near type flags. 1 means standard pythagorean nearest. See note for what 0 means.     */     
/* tree_numd     Number of pointers in tree_dl and flags in tree_nt (i.e. number of dimensions).       */     
/* extent_min    Minimum value of this dimension to find nearest point. Size tree_numd.                */     
/*               Greater OR EQUAL to this value is in range.                                           */     
/* extent_max    Maximum value of this dimension to find nearest point. Size tree_numd.                */     
/*               Strictly LESS than this value is in range.                                            */     
/* target        Find the point nearest here considering extents and tree_nt. Size tree_numd.          */     
/* init_rd       Initially search within radius init_rad. If no point within extents is found,         */     
/*               init_rad is multiplied by rad_scal and the seach is repeated. And so on.              */     
/* rad_scal      Scale value to apply to init_rd until a point is found within extents, or the         */     
/*               init_rad has been scaled to larger than the extents.                                  */     
/*                                                                                                     */     
/* Output arguments:                                                                                   */     
/* near_elem     The element number of a nearest point in the tree_dl arrays. For instance, a value of */     
/*               m means the coordinates of point are tree_dl[n][m] where n is the dimension number.   */     
/* near_dist     The SQUARED distance between the target and nearest point (if num_found>0).           */
/* num_found     0 if no point is found within the extent ranges.                                      */
/*              >0 is as many equally-near points as are found within the extent ranges.               */
/*                 The returned near_elem is the highest-numbered elem of all the nearest points.      */
/*                                                                                                     */     
/* ncycles       Number of cycles to find nearest point. This is just for informational purposes.      */
/*               This is the number of times that the rad_scal had to be applied before finding the    */
/*               nearest point. This could help set init_rad and rad_scal for faster searches.         */
/*                                                                                                     */     
/*    NOTE: A tree_nt value of 0 means the difference between the point coordinate and the target      */     
/*          coordinate are not added to the Pythagorean sum for the specified dimension(s).            */     
/*          So, the point distances are determined as if that dimension was not specified at all.      */     
/*          However, the extent_min,extent_max values for that dimension are still used.               */     
/*          The result is therefor the nearest point considering only the other dimensions, but        */     
/*          still restricted by the extent range of that dimension(s).                                 */     

  int naxe = 1;
  node *now_node;
  double now_rad = init_rad;
  double loc_min[9];
  double loc_max[9]; 
  int nset = 0;
  int i = 0;
  int nbig = 0;

  for(i=0; i<tree_numd; i++) {
    if(tree_nt[i] == 0) {
      loc_min[i] = extent_min[i];
      loc_max[i] = extent_max[i];
      nset+=2;
    }
  }

  *ncycles = 0;

CYCLE:

  *ncycles = *ncycles + 1;

/* Set extents using current radius. Also set nbig which tells us whether */
/* we are beyond the extents (therefor no reason to keep looking).        */

  nbig = nset;
  for(i=0; i<tree_numd; i++) {
    if(tree_nt[i] != 0) {
      loc_min[i] = target[i] - now_rad;
      loc_max[i] = target[i] + now_rad;
      if(loc_min[i] <= extent_min[i]) {
        loc_min[i] = extent_min[i];
        nbig++;
      }
      if(loc_max[i] >= extent_max[i]) { /* yes, >= is better than > here */ 
        loc_max[i] = extent_max[i];
        nbig++;
      }
    }
  }

  *num_found = 0;
  *near_dist = DBL_MAX;
  *near_elem = 0; 
  
  naxe = 1;
  now_node = tree_nodes; 
  find_near_rcur(tree_dl, tree_nt, tree_numd, now_node, naxe, 
                 loc_min, loc_max, target, 
                 near_elem, near_dist, num_found);

  if(*num_found<1) {
    if(nbig==2*tree_numd) return; /* None found. Are we outside the extents?  */
    now_rad = now_rad * rad_scal;
    goto CYCLE;
  }

/* Here, we need to consider the difference between a square and a circle.    */ 
/* The now_rad value is half the size of the square we just searched. So, if  */
/* current nearest point is in a circle with that radius, we are finished.    */
/* But, otherwise, there might be nearer points hiding in the area outside    */
/* the searched-square, but inside the circle. So increase the search size    */
/* to the CURRENT nearest point radius. On the next cycle, we will definitly  */        
/* get the nearest point because the CURRENT point is definitly within the    */
/* square with the now_rad that we are now setting. So the CURRENT point here */
/* will be among the points returned by the next cycle. And its radius will   */
/* definitly satisfy this condition because we explicitly made now_rad big    */
/* enough (but another point hiding in the square-circle area might sneak in).*/

  if(*near_dist > now_rad*now_rad && nbig!=2*tree_numd) {
    now_rad = sqrt(*near_dist) * 1.001;
    goto CYCLE;
  }

  return;
}
/*                                                                     */
/* -----------------------------------------------------------         */
/* Specify compare function for qsort.                                 */

int compSort1 (const void * q1, const void * q2) {

  struct QInfo* p1 = (struct QInfo*) q1; 
  struct QInfo* p2 = (struct QInfo*) q2; 

  if(p1->dlots[locp] < p2->dlots[locp]) return (-1);
  if(p1->dlots[locp] > p2->dlots[locp]) return (1); 

  return (0); 

}
/*                                                                    */
/* ------------------------------------------------------------------ */
/*    Running average.  (copied from suprofcsv Oct 2023, added dtog). */
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
