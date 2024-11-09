/* Copyright (c) Colorado School of Mines, 2022.*/
/* All rights reserved.                       */

/* SUNEARCSV: $Revision: 1.01 $ ; $Date: 2024/10/20 00:00:01 $		*/
 
#include "su.h"
#include "segy.h" 
#include <stdbool.h>
#include "qdefine.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SUNEARQCSV - Assign Values To Target Q-file From Nearest In Another Q-file.",
"									     ",
"  sunearqcsv [parameters].    (No traces in or out).                        ",
"									     ",
" Parameters:	         						     ",
"                                                                            ",
" qin=      Input q-file to copy nearest values FROM. The values to copy are ",
"           specified via the oreps parameter and the otuples parameter.     ",
"									     ",
" dimx=     Dimension name X. Name must exist in the qin q-file.             ",
" dimy=     Dimension name Y. Name must exist in the qin q-file.             ",
" dimr=     Dimension name R. Name must exist in the qin q-file.             ",
" dima=     Dimension name A. Name must exist in the qin q-file.             ",
" dimb=     Dimension name B. Name must exist in the qin q-file.             ",
" dimc=     Dimension name C. Name must exist in the qin q-file.             ",
" dimd=     Dimension name D. Name must exist in the qin q-file.             ",
" dime=     Dimension name E. Name must exist in the qin q-file.             ",
" dimf=     Dimension name F. Name must exist in the qin q-file.             ",
"     Note: No defaults. At least 1 dimension must be used. And for every    ",
"           qin dimension parameter name specified the corresponding         ",
"           tin target dimension parameter name must be specified.           ",
"									     ",
" The following 3 parameters can be used for any of the dimensions.	     ",
" Just substitute x,y,r,a,b,c,d,e,f as the ending character,                 ",
" for instance: typer,minr,maxr for dimension R.                             ",
" Typically these 3 are used just to specify an additional dimension that is ",
" used to restrict the search range (see typer options -1 and -2 below).     ",
"									     ",
" typer=1 Use in Pythagorean Nearest (squared difference to tin location),   ",
"         and minr,maxr are unchanging extent ranges for this dimension.     ",
"      =2 Use in Pythagorean Nearest and minr,maxr are the relative extent   ",
"         range to include for this dimension (the input tin value for       ",
"         this dimension is added to minr,maxr to form the extent range).    ",
"     =-1 Do not use in Pythagorean Nearest, and minr,maxr are unchanging    ",
"         extent ranges for this dimension.                                  ",
"     =-2 Do not use in Pythagorean Nearest and minr,maxr are the relative   ",
"         extent range to include for this dimension (the input tin value    ",
"         for this dimension is added to minr,maxr to form extent range).    ",
"    Note: Negative types means the difference between the points and        ",    
"       the tin is NOT added to Pythagorean sum for specified dimension(s).  ",
"       So, the nearest point is determined as if that dimension was NOT     ",
"       specified at all. However, the min,max for that dimension are still  ",     
"       used. This results in finding the nearest point considering only the ",
"       type>0 dimensions, but restricted by ranges of type<0 dimension(s).  ",     
"       Type=-2 can be used for situations where a profile approaches itself ",     
"       or intersects itself or overlaps itself. In those case, a third      ",     
"       dimension value (such as station number) can be used to restrict the ",     
"       search range for each tin to only the part of the profile with       ",     
"       approximately the same station numbers as the tin midpoint station.  ",     
" minr=     Extent Range Min. Typically negative. Default is no min limit.   ",
"           Greater or equal to this value is in range.                      ",
" maxr=     Extent Range Max. Typically positive. Default is no max limit.   ",
"           Strictly less than this value is in range. The min must be less  ",
"           than max, but the range does not have to be symmetric or centred ",
"           (that is, ranges such as min=-1000 and max=-200 are allowed).    ",
"									     ",
" ordp=asis Point Order Name. Sometimes there are two-or-more qin points     ",
"           which are equally near the tin target location. In those cases,  ",
"           the nearest point is set to the higher ordered point.            ",
"           The default (asis) means use the input order of qin records.     ",
"           Specifying a name here will cause the input qin records to be    ",
"           sorted internally by the values of the specified name.           ",
"     Note: The sunearcsv program has a similar parameter called keyp.       ",
"           This parameter matches its ordering functionality.               ",
"									     ",
" tin=      Target Input q-file to copy nearest values FROM qin q-file.      ",
"									     ",
" timx=     Target Dimension name X. Name must exist in the tin q-file.      ",
" timy=     Target Dimension name Y. Name must exist in the tin q-file.      ",
" timr=     Target Dimension name R. Name must exist in the tin q-file.      ",
" tima=     Target Dimension name A. Name must exist in the tin q-file.      ",
" timb=     Target Dimension name B. Name must exist in the tin q-file.      ",
" timc=     Target Dimension name C. Name must exist in the tin q-file.      ",
" timd=     Target Dimension name D. Name must exist in the tin q-file.      ",
" time=     Target Dimension name E. Name must exist in the tin q-file.      ",
" timf=     Target Dimension name F. Name must exist in the tin q-file.      ",
"     Note: A target dimension must be specified for each dim dimension used.",
"									     ",
" qout=qout.csv    Output q-file name.                                       ",
"									     ",
" oreps=    Output Replace Name List.                                        ",
"           Names listed here must be in the non-tuples part of qin q-file.  ",
"           If name here already exists in non-tuples part of tin file, its  ",
"           nearest value from the qin file are REPLACED in the output file. ",
"           If name does not exist in non-tuples part of tin file, its name  ",
"           and nearest value from the qin file are ADDED to the output file.",
"           Names and values NOT listed here are COPIED from tin to output.  ",
"     Note: This parameter does not have to be specified but the resulting   ",
"           output might just be a copy of input depending on the following  ",
"           otuples parameter and whether tuples exist in qin or tin.        ",
"									     ",
" otuples=1 Output the tuples asis from tin file (if they exist).            ",
"        =0 Do not output any tuples (even if they exist).                   ",
"        =2 Output the tuples of nearest point in qin file (if they exist).  ",
"           Tuples cannot be intermixed (either all from tin, or from qin).  ",
"           Tuples can have varying or fixed amounts (see Q File Standards). ",
"           Tuple names copied from the qin file are NOT error-checked to see",
"           if the same name exists in the non-tuple names in the tin file.  ",
"									     ",
" nopoint=1 If no qin point is found within the extents, error-halt (before  ",
"           outputting the first q-record that could not find a qin point).  ",
"        =0 Continue without any update of oreps values for q-records that   ",
"           could not find a qin point in extents (print their count at end).",
"       *** Since the default extent min,max are unlimited, this option has  ",
"           no effect unless you specify them.                               ",
"     Note: If you choose to use option 0, then values for oreps that exist  ",
"           in tin file are output with their tin values. But if the oreps   ",
"           name does not exist in the tin file, then output value is -9999. ",
"									     ",
" formxy=%.20g  The C format code for printing all values to qout records.   ",
"              Note that the default format prints up to 20 digits           ",
"              (but not trailing zeroes to the right of the decimal point).  ",
"                                                                            ",
"									     ",
" The following 3 parameters affect cpu time, but not results. The search    ",
" is done by building a kdtree from the dimension values in the qin file.    ",
" The code that builds the kdtree is reasonably standard and simplistic.     ",
" But the kdtree search code is slightly unusual due to some typer options.  ",
"									     ",
" sfunc=2   Search via the cycle_for_near function. This option is usually   ",
"           fastest. This option uses the sdist and smult parameters.        ",
"      =1   Search via the find_near function. This option may be faster     ",
"           if you are specifying small extent ranges.                       ",
"           This option does not use the sdist and smult parameters.         ",
"      =0   Search via the brute_near function. This option may be faster    ",
"           if there are only a small number of points in the qin file.      ",
"           This option does not use the sdist and smult parameters.         ",
"           (To keep the code simpler, this option still allocates and       ",
"            builds the kdtree. But the tree is not used for searching).     ",
"      =-1  Search via the find_near function and confirm the results using  ",
"           the brute_near function. This tests the find_near function and   ",
"           can also be used to determine the speed difference between them. ",
"           This option does not use the sdist and smult parameters.         ",
"      =-2  Search via the cycle_for_near function and confirm the results   ",
"           using brute_near function. This tests cycle_for_near function    ",
"           and can also be used to determine speed difference between them. ",
"           This option uses the sdist and smult parameters.                 ",
" sdist=100 Initial search distance. This parameter is only used             ",
"           if sfunc=2 or -2. If positive, this value is added to the        ",
"           distance between the previous tin and its nearest qin point      ",
"           and used as the initial search distance for the current tin.     ",
"           If negative, the initial search distance for all tin is set      ",
"           to this absolute value.                                          ",
" smult=2   Search multiplier. This parameter is only used if sfunc=2 or -2. ",
"           If the nearest qin point is not found after searching with the   ",
"           initial search distance, the search distance is multiplied by    ",
"           this value, and search is performed again. This repeats until    ",
"           finding the nearest point (or all min,max ranges are exceeded).  ",
"									     ",
"   ------------------------------------------------------------------       ",
"   ------------------------------------------------------------------       ",
"									     ",
NULL};

/* Created: Oct  2024: Andre Latour                                          */ 
/* This program started from sunearcsv.                                      */ 
/**************** end self doc *******************************************/

segy tr;

struct QInfo *RecInfo; /* Storage for all function location value pointers */
struct QInfo *RecInfot; /* Storage for all function location value pointers */
int locp = -1;      

int compSort1 (const void * q1, const void * q2) ; /* comparison function for qsort  */
 
/* Note: I make no claim that this is a particularly good kd tree implementation.     */
/*       It is not explicitly balanced (it has an option to get approximate balance). */
/* Note: Option tree_nt is unlikely to exist in other kd tree implementations.        */
/*       It exists because some crooked-profiles (land) or coil-profiles (marine)     */
/*       curve back-over-top-of-themselves. These self-intersections mean the profile */
/*       is not a function (in the mathematical sense). That is, just considering XYs */
/*       a tin midpoint can get confused as to which part of the profile it should    */
/*       belong to. Similar issues arise if the profile just curves back NEAR itself. */

typedef struct node {
   unsigned long long elem; 
   struct node * l;
   struct node * r;
} node;

/* Note: Yes, I could also have made a structure named tree, and then passed it into  */
/*       the functions instead of some of the individual arguments. That would have   */
/*       reduced the function arguments, but make it more confusing for new coders.   */

void connect_nodes (node *tree_nodes, unsigned long long tree_numc, double **tree_dl, int tree_numd,
                   int ihop);

void connect_nodes_your (node *tree_nodes, unsigned long long tree_numc, double **tree_dl, int tree_numd,
                        unsigned long long *your_order);

void connect_all (node *tree_nodes, unsigned long long tree_numc, double **tree_dl, int tree_numd);

void find_in (node *tree_nodes, double **tree_dl, int tree_numd,
             double *extent_min, double *extent_max,  
             unsigned long long *out_elem, unsigned long long *num_out);

void find_in_rcur (double **tree_dl, int tree_numd, node *now_node, int naxe,
                  double *extent_min, double *extent_max, 
                  unsigned long long *out_elem, unsigned long long *num_out);

void find_near (node *tree_nodes, double **tree_dl, int *tree_nt, int tree_numd,
               double *extent_min, double *extent_max, double *target, 
               unsigned long long *near_elem, double *near_dist, unsigned long long *num_found);

void find_near_rcur (double **tree_dl, int *tree_nt, int tree_numd, node *now_node, int naxe, 
                    double *extent_min, double *extent_max, double *target,
                    unsigned long long *near_elem, double *near_dist, unsigned long long *num_found);

void cycle_for_near (node *tree_nodes, double **tree_dl, int *tree_nt, int tree_numd,
                    double *extent_min, double *extent_max, double *target, 
                    double init_rad, double rad_scal,  
                    unsigned long long *near_elem, double *near_dist, 
                    unsigned long long *num_found, int *ncycles);

void brute_near (double **tree_dl, unsigned long long tree_numc, int *tree_nt, int tree_numd,
                 double *extent_min, double *extent_max, double *target, 
                 unsigned long long *near_elem, double *near_dist, unsigned long long *num_found);

/*----------------------------------------------------------------------*/

int main(int argc, char **argv) {

  int ncdp = 0;		
  int ifixd = 0;          /* flag for all tuples same size or vary   */
  int iztuple = 0;        /* element number where first tuple exists */
  int ktuple = 0;         /* type of tuples (2=pairs, 3=triplets)    */
  cwp_String Pname=NULL;  /* text file name for Q input file      */
  FILE *fpP=NULL;         /* file pointer for Q input file        */

  cwp_String *pname = NULL;                                         
  cwp_String *ndims = NULL;                                                 
  int numpname = 0;
  double *pindepa = NULL;                                               
  int numdind = 0;
	

  int ncdpt = 0;	                                     
  int ifixdt = 0;          /* flag for all tuples same size or vary   */
  int iztuplet = 0;        /* element number where first tuple exists */
  int ktuplet = 0;         /* type of tuples (2=pairs, 3=triplets)    */
  cwp_String Tname=NULL;  /* text file name for Q input file      */
  FILE *fpT=NULL;         /* file pointer for Q input file        */


  cwp_String *pnamet = NULL;                                         
  cwp_String *ndimst = NULL;                                                 
  int numpnamet = 0;
  double *pindepat = NULL;                                               
  int numdindt = 0;
  int moret = 0;

  cwp_String formxyt=NULL;
  cwp_String formxy=NULL;
  cwp_String formxylong=NULL;
  int lenformxy = 0;

  int jcdp = 0;      
  int i = 0;                                                                         
  int j = 0;                                                                       
  int k = 0;                                                                       
  int errwarn = 0; 

  cwp_String keyn[9];  
  int locn[9];
  int tree_numd = 0; 

  cwp_String keyt[9];  
  int locnt[9];

  int num_ordp = 0;

  int num_oreps = 0;
  cwp_String *oreps = NULL;  
  int *oloc  = NULL;  
  int *oloct = NULL;  
  int *olocm = NULL;  

  int otuples = 1;
  int nopoint = 1;
  double dmiss = -9999.;

  double sdist = 100.0;
  double smult = 2.;
  int ncycles = 0;
  int tcycles = 0;
  int sfunc   = 2;
  int scheck  = 0;
  int ncheck  = 0;
  int sdadd   = 1;
  double sdist2 = 100.0;

  unsigned long long near_elem = 0;
  unsigned long long num_found = 0;
  double near_dist = 0.;
  unsigned long long near_elemc = 0;
  unsigned long long num_foundc = 0;
  double near_distc = 0.;
  int ihop = 1;
  int nproct = 0;
  int noupdate = 0;

  int *in_rel_extent = NULL;
  double *extent_in_min = NULL;
  double *extent_in_max = NULL;
  double **tree_dl = NULL;
  int *tree_nt = NULL;
  double *extent_min = NULL;
  double *extent_max = NULL;
  double *target = NULL;
  unsigned long long tree_ncdp = 0;
  node *tree_nodes = NULL;

  cwp_String Oname=NULL;  /* text file name for O output file     */
  FILE *fpO=NULL;         /* file pointer for O output file       */

/* hook up getpar */
  initargs(argc, argv);
  requestdoc(1);

  if(isatty(STDIN_FILENO)!=1 || isatty(STDOUT_FILENO)!=1)
    err("**** Error: this program does not input or output traces.");

/* Maximum number of dimensions is 9.                                       */

  in_rel_extent = ealloc1int(9);
  extent_in_min = ealloc1double(9);
  extent_in_max = ealloc1double(9);
  tree_nt = ealloc1int(9);
  extent_min = ealloc1double(9);
  extent_max = ealloc1double(9);
  target = ealloc1double(9);
  tree_dl = ealloc1(9,sizeof(double *));

/* Set defaults for extent ranges.                                          */

  for(i=0; i<9; i++) {
    in_rel_extent[i] = 1; /* use in Pythagorean Nearest, has unchanging range */
    extent_in_min[i] = -DBL_MAX;
    extent_in_max[i] =  DBL_MAX;
  }

/* Note here that keyn is loaded with key names in the order in which the   */
/* dimension parameters are checked next. And keyn is never re-arranged.    */
/* This means that, eventually, the target array (named target) also must   */
/* have its values in the same order.                                       */

  if(countparval("dimx") > 0) {
    getparstring("dimx",keyn+tree_numd );
    if(countparval("typex")>0) getparint("typex",in_rel_extent+tree_numd);
    if(countparval("minx") >0) getpardouble("minx",extent_in_min+tree_numd);
    if(countparval("maxx") >0) getpardouble("maxx",extent_in_max+tree_numd);
    if(countparval("timx") > 0) getparstring("timx",keyt+tree_numd );
    else err("**** Error: Parameter dimx specified, but not timx.");
    tree_numd++;
  }
  else if(countparval("typex")>0 || countparval("minx")>0 || countparval("maxx")>0 || countparval("timx")>0) {
     err("**** Error: Parameter typex, minx, maxx, or timx specified, but not dimx.");
  }

  if(countparval("dimy") > 0) {
    getparstring("dimy",keyn+tree_numd );
    if(countparval("typey")>0) getparint("typey",in_rel_extent+tree_numd);
    if(countparval("miny") >0) getpardouble("miny",extent_in_min+tree_numd);
    if(countparval("maxy") >0) getpardouble("maxy",extent_in_max+tree_numd);
    if(countparval("timy") > 0) getparstring("timy",keyt+tree_numd );
    else err("**** Error: Parameter dimy specified, but not timy.");
    tree_numd++;
  }
  else if(countparval("typey")>0 || countparval("miny")>0 || countparval("maxy")>0 || countparval("timy")>0) {
     err("**** Error: Parameter typey, miny, maxy, timy specified, but not dimy.");
  }

  if(countparval("dimr") > 0) {
    getparstring("dimr",keyn+tree_numd );
    if(countparval("typer")>0) getparint("typer",in_rel_extent+tree_numd);
    if(countparval("minr") >0) getpardouble("minr",extent_in_min+tree_numd);
    if(countparval("maxr") >0) getpardouble("maxr",extent_in_max+tree_numd);
    if(countparval("timr") > 0) getparstring("timr",keyt+tree_numd );
    else err("**** Error: Parameter dimr specified, but not timr.");
    tree_numd++;
  }
  else if(countparval("typer")>0 || countparval("minr")>0 || countparval("maxr")>0 || countparval("timr")>0) {
     err("**** Error: Parameter typer, minr, maxr, timr specified, but not dimr.");
  }

  if(countparval("dima") > 0) {
    getparstring("dima",keyn+tree_numd );
    if(countparval("typea")>0) getparint("typea",in_rel_extent+tree_numd);
    if(countparval("mina") >0) getpardouble("mina",extent_in_min+tree_numd);
    if(countparval("maxa") >0) getpardouble("maxa",extent_in_max+tree_numd);
    if(countparval("tima") > 0) getparstring("tima",keyt+tree_numd );
    else err("**** Error: Parameter dima specified, but not tima.");
    tree_numd++;
  }
  else if(countparval("typea")>0 || countparval("mina")>0 || countparval("maxa")>0 || countparval("tima")>0) {
     err("**** Error: Parameter typea, mina, maxa, or tima specified, but not dima.");
  }

  if(countparval("dimb") > 0) {
    getparstring("dimb",keyn+tree_numd );
    if(countparval("typeb")>0) getparint("typeb",in_rel_extent+tree_numd);
    if(countparval("minb") >0) getpardouble("minb",extent_in_min+tree_numd);
    if(countparval("maxb") >0) getpardouble("maxb",extent_in_max+tree_numd);
    if(countparval("timb") > 0) getparstring("timb",keyt+tree_numd );
    else err("**** Error: Parameter dimb specified, but not timb.");
    tree_numd++;
  }
  else if(countparval("typeb")>0 || countparval("minb")>0 || countparval("maxb")>0 || countparval("timb")>0) {
     err("**** Error: Parameter typeb, minb, maxb, or timb specified, but not dimb.");
  }

  if(countparval("dimc") > 0) {
    getparstring("dimc",keyn+tree_numd );
    if(countparval("typec")>0) getparint("typec",in_rel_extent+tree_numd);
    if(countparval("minc") >0) getpardouble("minc",extent_in_min+tree_numd);
    if(countparval("maxc") >0) getpardouble("maxc",extent_in_max+tree_numd);
    if(countparval("timc") > 0) getparstring("timc",keyt+tree_numd );
    else err("**** Error: Parameter dimc specified, but not timc.");
    tree_numd++;
  }
  else if(countparval("typec")>0 || countparval("minc")>0 || countparval("maxc")>0 || countparval("timc")>0) {
     err("**** Error: Parameter typec, minc, maxc, or timc specified, but not dimc.");
  }

  if(countparval("dimd") > 0) {
    getparstring("dimd",keyn+tree_numd );
    if(countparval("typed")>0) getparint("typed",in_rel_extent+tree_numd);
    if(countparval("mind") >0) getpardouble("mind",extent_in_min+tree_numd);
    if(countparval("maxd") >0) getpardouble("maxd",extent_in_max+tree_numd);
    if(countparval("timd") > 0) getparstring("timd",keyt+tree_numd );
    else err("**** Error: Parameter dimd specified, but not timd.");
    tree_numd++;
  }
  else if(countparval("typed")>0 || countparval("mind")>0 || countparval("maxd")>0 || countparval("timd")>0) {
     err("**** Error: Parameter typed, mind, maxd, timd specified, but not dimd.");
  }

  if(countparval("dime") > 0) {
    getparstring("dime",keyn+tree_numd );
    if(countparval("typee")>0) getparint("typee",in_rel_extent+tree_numd);
    if(countparval("mine") >0) getpardouble("mine",extent_in_min+tree_numd);
    if(countparval("maxe") >0) getpardouble("maxe",extent_in_max+tree_numd);
    if(countparval("time") > 0) getparstring("time",keyt+tree_numd );
    else err("**** Error: Parameter dime specified, but not time.");
    tree_numd++;
  }
  else if(countparval("typee")>0 || countparval("mine")>0 || countparval("maxe")>0 || countparval("time")>0) {
     err("**** Error: Parameter typee, mine, maxe, or time specified, but not dime.");
  }

  if(countparval("dimf") > 0) {
    getparstring("dimf",keyn+tree_numd );
    if(countparval("typef")>0) getparint("typef",in_rel_extent+tree_numd);
    if(countparval("minf") >0) getpardouble("minf",extent_in_min+tree_numd);
    if(countparval("maxf") >0) getpardouble("maxf",extent_in_max+tree_numd);
    if(countparval("timf") > 0) getparstring("timf",keyt+tree_numd );
    else err("**** Error: Parameter dimf specified, but not timf.");
    tree_numd++;
  }
  else if(countparval("typef")>0 || countparval("minf")>0 || countparval("maxf")>0 || countparval("timf")>0) {
     err("**** Error: Parameter typef, minf, maxf, timf specified, but not dimf.");
  }

  if(tree_numd<1) err("**** Error: At least 1 dimension name must be specified.");

/* Resolve a few things. ---------------------------------------------------- */

  for(i=0; i<tree_numd; i++) {

    extent_min[i] = extent_in_min[i]; /* Set incase this is not relative */
    extent_max[i] = extent_in_max[i]; /* Set incase this is not relative */
    if(in_rel_extent[i]>0) {
      tree_nt[i] = 1; /* use as part of Pythagorean Nearest */
    }
    else {
      tree_nt[i] = 0; /* do not use as part of Pythagorean Nearest */
      in_rel_extent[i] = 0 - in_rel_extent[i]; /* reset to positive 1 or 2 */
    }

  } /* end of  for(i=0; i<tree_numd; i++) { */

  cwp_String ordp = NULL;  
  if(countparval("ordp") > 0) {
    getparstring("ordp", &ordp);
    if(strcmp(ordp,"asis")!=0) num_ordp = 1;   
  }

  if(countparval("oreps")>0) {
    num_oreps = countparval("oreps");
    oreps = ealloc1(num_oreps,sizeof(cwp_String *)); 
    getparstringarray("oreps", oreps);
    oloc  = ealloc1int(num_oreps); 
    oloct = ealloc1int(num_oreps); 
    olocm = ealloc1int(num_oreps); 
  }    

  if(!getparint("otuples",&otuples)) otuples = 1;
  if(otuples>2 || otuples<0) err("**** Error: otuples must be 0, 1, or 2."); 

  if(!getparint("nopoint",&nopoint)) nopoint = 1;
  if(nopoint>1 || nopoint<0) err("**** Error: nopoint must be 1 or 0."); 

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

  if(!getparint("sfunc",&sfunc)) sfunc = 2;
  if(sfunc<-2 || sfunc>2) err("**** Error: sfunc option is out-of-range."); 

  scheck = 0;
  if(sfunc<0) {
    sfunc = 0 - sfunc;
    scheck = 1;
  }

  if(!getpardouble("sdist",&sdist)) sdist = 100.;
  if(sdist==0.)  err("**** Error: sdist cannot be 0."); 

  if(!getpardouble("smult",&smult)) smult = 2.;
  if(smult<0.)  err("**** Error: smult must be non-negative."); 

  if(sdist<0.) {
    sdist = 0. - sdist;
    sdadd = 0;
  }
  sdist2 = sdist;

/*--------------------------------------------------------------------------  */

  getparstring("qin", &Pname);

  fpP = fopen(Pname, "r");
  if(fpP==NULL) err("error: qin input q-file did not open correctly.");
    
/* Set input numpname,pname to just store what is going to be accessed.       */
/* numpname>0 is a flag to ONLY store values if they are on pname list.       */
/* numpname<1 is a flag to NOT store values if they are on pname list.        */

/*otuples=1 Output the tuples asis from tin file.                             */
/*       =2 Output the tuples of nearest point in qin file.                   */
/*       =0 Do not output any tuples.                                         */

  pname = ealloc1(tree_numd + num_ordp + num_oreps,sizeof(cwp_String *));
  numpname = 0;
  
  if(otuples!=2) { /* Output tuples from qin. So input all values from qin.   */    
    for (i=0; i<tree_numd; ++i) {
      pname[numpname] = ealloc1(strlen(keyn[i]),1);
      strcpy(pname[numpname],keyn[i]);
      numpname++;
    }
    for (i=0; i<num_ordp; ++i) { /* currently just 1 sort name is allowed */
      pname[numpname] = ealloc1(strlen(ordp),1);
      strcpy(pname[numpname],ordp);
      numpname++;
    }
    for (i=0; i<num_oreps; ++i) { 
      pname[numpname] = ealloc1(strlen(oreps[i]),1);
      strcpy(pname[numpname],oreps[i]);
      numpname++;
    }
  }

  getviaqfile(fpP, &pname, &numpname, &iztuple, numdind,   
              &ktuple, &ifixd, &RecInfo, &ncdp, 
              &pindepa,  &ndims, &errwarn) ;

  if(errwarn==1) err("getqinfo error: qin file, extra C_SU_NAMES record");
  else if(errwarn==2) err("getqinfo error: qin file, extra C_SU_NDIMS record");
  else if(errwarn==3) err("getqinfo error: qin file, C_SU_ID record not found immediately after C_SU_NAMES record.");
  else if(errwarn==11) 
    err("readqhead error: qin file, if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES");
  else if(errwarn==12) 
    err("readqhead error: qin file, C_SU_ID record not found immediately after C_SU_NAMES record.");
  else if(errwarn==22) err("getviaqfile error: qin file, C_SU_NDIMS record not same length as C_SU_NAMES record.");
  else if(errwarn==23) err("getviaqfile error: qin file, C_SU_NAMES tupled names out-of-order, changed");
  else if(errwarn==24) err("getviaqfile error: qin file, C_SU_NDIMS blank where valid number expected");
  else if(errwarn==25) err("getviaqfile error: qin file, C_SU_NDIMS non-number where valid number expected");
  else if(errwarn==26) err("getviaqfile error: qin file, C_SU_NDIMS value must be same for all members of tuple");
  else if(errwarn==27) err("getviaqfile error: qin file, C_SU_NAMES record followed by C_SU_ID record not found.");
  else if(errwarn>100) 
    err("getviaqfile error: qin file, record %d (wrong comma count, damaged, non-numbers, ...)",errwarn-99);
  else if(errwarn>0) err("getviaqfile error: qin file, unrecognized error code %d",errwarn);

  getparstring("tin", &Tname);

  fpT = fopen(Tname, "r");
  if(fpT==NULL) err("error: tin input q-file did not open correctly.");

  errwarn = 0;  
  getviaqfile(fpT, &pnamet, &numpnamet, &iztuplet, numdindt,   
              &ktuplet, &ifixdt, &RecInfot, &ncdpt, 
              &pindepat,  &ndimst, &errwarn) ;

  if(errwarn==1) err("getqinfo error: tin file, extra C_SU_NAMES record");
  else if(errwarn==2) err("getqinfo error: tin file, extra C_SU_NDIMS record");
  else if(errwarn==3) err("getqinfo error: tin file, C_SU_ID record not found immediately after C_SU_NAMES record.");
  else if(errwarn==11) 
    err("readqhead error: tin file, if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES");
  else if(errwarn==12) 
    err("readqhead error: tin file, C_SU_ID record not found immediately after C_SU_NAMES record.");
  else if(errwarn==22) err("getviaqfile error: tin file, C_SU_NDIMS record not same length as C_SU_NAMES record.");
  else if(errwarn==23) err("getviaqfile error: tin file, C_SU_NAMES tupled names out-of-order, changed");
  else if(errwarn==24) err("getviaqfile error: tin file, C_SU_NDIMS blank where valid number expected");
  else if(errwarn==25) err("getviaqfile error: tin file, C_SU_NDIMS non-number where valid number expected");
  else if(errwarn==26) err("getviaqfile error: tin file, C_SU_NDIMS value must be same for all members of tuple");
  else if(errwarn==27) err("getviaqfile error: tin file, C_SU_NAMES record followed by C_SU_ID record not found.");
  else if(errwarn>100) 
    err("getviaqfile error: tin file, record %d (wrong comma count, damaged, non-numbers, ...)",errwarn-99);
  else if(errwarn>0) err("getviaqfile error: tin file, unrecognized error code %d",errwarn);

  getparstring("qout", &Oname);

  if(strcmp(Pname,Oname) == 0 || strcmp(Tname,Oname) == 0 )
    err("**** Error: qout output file name must be different than qin and tin input file names.");

  fpO = fopen(Oname, "w");
  if (fpO == NULL) err("qfile error: qout output file did not open correctly.");

  checkpars(); 

  if(otuples==1 && ifixdt==2) otuples = 0;
  if(otuples==2 && ifixd ==2) otuples = 0;

/*--------------------------------------------------------------------------  */

  for (i=0; i<tree_numd; ++i) {
    locn[i] = -1;
    for (j=0; j<iztuple; ++j) if(strcmp(pname[j],keyn[i])==0) locn[i] = j;
    if(locn[i] < 0) err("error: name %s not found in non-tuple part of qin q-file.",keyn[i]);
  }

  for (i=0; i<num_oreps; ++i) {
    oloc[i] = -1;
    for (j=0; j<iztuple; ++j) if(strcmp(pname[j],oreps[i])==0) oloc[i] = j;
    if(oloc[i]<0) err("error: input qin q-file must have your oreps=%s (among non-tuple names).",oreps[i]);
  }

  if(num_ordp>0) {
    locp = -1;
    for (i=0; i<iztuple; ++i) if(strcmp(pname[i],ordp)==0) locp = i;
    if(locp<0) err("error: input qin q-file must have your ordp=%s (among non-tuple names).",ordp);

/* Note that locp is used within compSort1. Note that this sort affects       */
/* cycle_for_near and find_near since they return the lowest-element point    */
/* when points are equally-near the tin (such as on boundary between 2 cdps). */

    qsort(RecInfo,ncdp,sizeof(struct QInfo),compSort1);

    for (jcdp=1; jcdp<ncdp; jcdp++) { 
      if(RecInfo[jcdp-1].dlots[locp] == RecInfo[jcdp].dlots[locp]) 
        warn("warning: Two input qin records same value %f  for ordp=%s",RecInfo[jcdp].dlots[locp],ordp);
    }
  }

/* similar things for the tin q-file */

  for (i=0; i<tree_numd; ++i) {
    locnt[i] = -1;
    for (j=0; j<iztuplet; ++j) if(strcmp(pnamet[j],keyt[i])==0) locnt[i] = j;
    if(locnt[i] < 0) err("error: tim dimension name %s not found in non-tuple part of tin q-file.",keyt[i]);
  }

  moret = 0;
  for (i=0; i<num_oreps; ++i) {
    oloct[i] = -1;
    olocm[i] = -1;
    for (j=0; j<iztuplet; ++j) if(strcmp(pnamet[j],oreps[i])==0) oloct[i] = j;
    if(oloct[i]<0) {
      olocm[moret] = oloc[i];
      moret++;
    }
  }

/*--------------------------------------------------------------------------  */

  tree_ncdp = ncdp;
  tree_nodes = ealloc1(tree_ncdp,sizeof(node));

  for(i=0; i<tree_numd; i++) {
    tree_dl[i] = ealloc1double(tree_ncdp);
    for(j=0; j<tree_ncdp; j++) tree_dl[i][j] = RecInfo[j].dlots[locn[i]];
  }

  ihop = 1;
  connect_nodes (tree_nodes, tree_ncdp, tree_dl, tree_numd, ihop);

/*--------------------------------------------------------------------------  */

  fprintf(fpO,"C_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID");
  for(i=0; i<iztuplet; i++) fprintf(fpO,",%s",formxy); 
  for (i=0; i<num_oreps; i++) if(oloct[i]<0) fprintf(fpO,",%s",formxy);
  if(otuples==1) { 
    if(ifixdt==0) { /* varying, just need to duplicate these formats once*/
      for(i=0; i<2; i++) {  
        for(k=0; k<ktuplet; k++) fprintf(fpO,",%s",formxy);
      }    
    }    
    else if(ifixdt==1) { 
      for(i=0; i<RecInfot[0].nto; i++) {  
        for(k=0; k<ktuplet-1; k++) fprintf(fpO,",%s",formxy);
      }    
    }
  }    
  else if(otuples==2) { 
    if(ifixd ==0) { /* varying, just need to duplicate tuple names once*/
      for(i=0; i<2; i++) { 
        for(k=0; k<ktuple; k++) fprintf(fpO,",%s",formxy);
      }    
    }    
    else if(ifixd ==1) { 
      for(i=0; i<RecInfo[0].nto; i++) { 
        for(k=0; k<ktuple-1; k++) fprintf(fpO,",%s",formxy);
      }    
    }    
  }    

  if(otuples==1 && ifixdt==1) {
    fprintf(fpO,"\nC_SU_NDIMS,%s",ndimst[0]);
    for(i=0; i<iztuplet-1+moret; i++) fprintf(fpO,",");
    for(i=0; i<RecInfot[0].nto; i++) { /* for fixed, use nto of first record  */ 
      for(k=0; k<ktuplet-1; k++) fprintf(fpO,",%.15g",pindepat[i]);
    }
  }
  if(otuples==2 && ifixd ==1) {
    fprintf(fpO,"\nC_SU_NDIMS,%s",ndims[0]);
    for(i=0; i<iztuplet-1+moret; i++) fprintf(fpO,","); /* yes, iztuplet,moret*/
    for(i=0; i<RecInfo[0].nto; i++) { /* for fixed, use nto of first record   */
      for(k=0; k<ktuple-1; k++) fprintf(fpO,",%.15g",pindepa[i]);
    }
  }

  fprintf(fpO,"\nC_SU_NAMES\nC_SU_ID");
  for(i=0; i<iztuplet; i++) fprintf(fpO,",%s",pnamet[i]);
  for (i=0; i<num_oreps; i++) if(oloct[i]<0) fprintf(fpO,",%s",oreps[i]); 
  
  if(otuples==1) { 
    if(ifixdt==0) { /* varying, just need to duplicate tuple names once*/
      for(i=0; i<2; i++) { 
        for(k=0; k<ktuplet; k++) fprintf(fpO,",%s",pnamet[k+iztuplet]);
      }    
    }
    else if(ifixdt==1) { 
      for(i=0; i<RecInfot[0].nto; i++) {  
        for(k=0; k<ktuplet-1; k++) fprintf(fpO,",%s",pnamet[k+iztuplet]);
      }    
    }
  }    
  else if(otuples==2) { 
    if(ifixd ==0) { /* varying, just need to duplicate tuple names once*/
      for(i=0; i<2; i++) { 
        for(k=0; k<ktuple; k++) fprintf(fpO,",%s",pname[k+iztuple]);
      }    
    }    
    else if(ifixd ==1) { 
      for(i=0; i<RecInfo[0].nto; i++) { 
        for(k=0; k<ktuple-1; k++) fprintf(fpO,",%s",pname[k+iztuple]);
      }    
    }    
  }    

  fprintf(fpO,"\n");

  if(otuples==1 && ifixdt==1) ktuplet--; 
  if(otuples==2 && ifixd ==1) ktuple--; 

/* loop over tin values   */ 

  for (jcdp=0; jcdp<ncdpt; jcdp++) {

    for(i=0; i<tree_numd; i++) {
      target[i] = RecInfot[jcdp].dlots[locnt[i]];  
      if(in_rel_extent[i]==2) { /* are extents relative to target location? */ 
        extent_min[i] = extent_in_min[i] + target[i];
        extent_max[i] = extent_in_max[i] + target[i];
      }
    }

/* Note that cycle_for_near and find_near return the higher-numbered element  */
/* when multiple points are equally-near the tin. This implies that each      */
/* point only contains one of the boundaries between itself and next points.  */
/* In particular, this is important for cdp-profiles since it explicitly      */
/* assigns tins that are mid-way between cdp centres to a specific cdp        */
/* no matter how the code is compilied/optimized (except, of course, double   */
/* precision optimization differences can still change distances so points    */
/* are not equally distant on all compilers/machines).                        */

    if(sfunc==2) {
      if(sdadd==1 && num_found>0) sdist2 = sdist + sqrt(near_dist);
      cycle_for_near(tree_nodes,tree_dl,tree_nt,tree_numd,
                     extent_min, extent_max, target,
                     sdist2, smult,
                     &near_elem, &near_dist, &num_found, &ncycles);
      tcycles += ncycles;
    }
    else if(sfunc==1) {
      find_near(tree_nodes,tree_dl,tree_nt,tree_numd,
                extent_min, extent_max, target,
                &near_elem, &near_dist, &num_found);
    }
    else {
      brute_near(tree_dl,tree_ncdp,tree_nt,tree_numd,
                 extent_min, extent_max, target,
                 &near_elem, &near_dist, &num_found);
    }

    if(scheck==1) {
      brute_near(tree_dl,tree_ncdp,tree_nt,tree_numd,
                 extent_min, extent_max, target,
                 &near_elemc, &near_distc, &num_foundc);
      if(near_elemc!=near_elem || near_distc!=near_dist || num_foundc!=num_found) {
        ncheck++;
        if(ncheck<10) {
          warn("near point: brute=%ld tree=%ld  dist diff (squared)=%g  equally near: brute=%ld tree=%ld  tin counter=%d ",
               near_elemc,near_elem,near_distc-near_dist,num_foundc,num_found,nproct+1);
        }
      }
    }

    if(num_found<1) {
      if(nopoint>0) err("error: No qin point within extents for tin q-record %d and nopoint option is 1.",nproct+1);
      noupdate++;
    }
    else {
      for (i=0; i<num_oreps; ++i) {
        if(oloct[i] > -1) RecInfot[jcdp].dlots[oloct[i]] = RecInfo[near_elem].dlots[oloc[i]];
      }
    }

    nproct++;

    fprintf(fpO,"Q");
    for (i=0; i<iztuplet; i++) fprintf(fpO,formxylong,RecInfot[jcdp].dlots[i]);
    if(num_found>0) {
      for (i=0; i<moret; i++) fprintf(fpO,formxylong,RecInfo[near_elem].dlots[olocm[i]]);
    }
    else {
      for (i=0; i<moret; i++) fprintf(fpO,formxylong,dmiss);
    }

    if(otuples==1) {
      for(i=0; i<RecInfot[jcdp].nto; i++) {  
        for(k=ktuplet-1; k>=0; k--) { 
          fprintf(fpO,formxylong,RecInfot[jcdp].dlots[iztuplet+k*RecInfot[jcdp].nto+i]);
        }
      }   
    }
    else if(otuples==2) {
      if(num_found>0) {
        for(i=0; i<RecInfo[near_elem].nto; i++) {  
          for(k=ktuple-1; k>=0; k--) { 
            fprintf(fpO,formxylong,RecInfo[near_elem].dlots[iztuple+k*RecInfo[near_elem].nto+i]);
          }
        }   
      }   
      else {
        for(k=ktuple-1; k>=0; k--) fprintf(fpO,formxylong,dmiss);
      }   
    }

    fprintf(fpO,"\n");

  } /* end of  for (jcdp=0; jcdp<ncdpt; jcdp++) { */

  if(sfunc==2) {
    warn("Number of tin q-records out=%d  Total search cycles=%d  Cycles per tin q-record=%d ",nproct,tcycles,tcycles/nproct);
  }
  else {
    warn("Number of q-records out =%d ",nproct);
  }
  if(noupdate>0) warn("Number of q-records not updated due to no qin point within extents = %d (permitted when nopoint=0).",noupdate);

  if(scheck==1) warn("There were %d total q-records out where brute_near disagreed with cycle_for_near or find_near.",ncheck);

  return 0;

} /* end of main for sunearqcsv */

/* -------------------------------------------------------------------------- */
/* Specify compare function for qsort.                                        */

int compSort1 (const void * q1, const void * q2) {

  struct QInfo* p1 = (struct QInfo*) q1;
  struct QInfo* p2 = (struct QInfo*) q2;

  if(p1->dlots[locp] < p2->dlots[locp]) return (-1);
  if(p1->dlots[locp] > p2->dlots[locp]) return (1); 

  return (0); 

}
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

void connect_nodes_your (node *tree_nodes, unsigned long long tree_numc, double **tree_dl, int tree_numd, 
                        unsigned long long *your_order) {

/*          This function connects the already allocated tree nodes in user specified sequence.        */     
/*                                                                                                     */     
/* Input arguments:                                                                                    */     
/* tree_nodes    The already fully allocated tree nodes.                                               */     
/* tree_numc     Number of nodes in tree_nodes.                                                        */     
/* tree_dl       Pointers to first element of values of each dimension. That is, each coordinate is in */     
/*               a separate contiguous array. These are the pointers to the start of those arrays.     */     
/*               Those arrays have size tree_numc.                                                     */     
/* tree_numd     Number of pointers in tree_dl (i.e. number of dimensions).                            */     
/* your_order    (length is tree_numc). Order in which points are to be processed.                     */         
/*               That is, this is the order in which the points should be added (loaded) into tree.    */         
/*               All values from 0 to tree_numc-1 should be somewhere in this array (none missing).    */         
/*               (To understand this, see the ihop argument of the connect_nodes function).            */         

  unsigned long long m = 0;

  for(m=0; m<tree_numc; m++) {
    tree_nodes[m].l    = 0; 
    tree_nodes[m].r    = 0;  
    tree_nodes[m].elem = your_order[m];
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

void find_near (node *tree_nodes, double **tree_dl, int *tree_nt, int tree_numd,
                double *extent_min, double *extent_max, double *target, 
                unsigned long long *near_elem, double *near_dist, unsigned long long *num_found) {

/*          This function finds a nearest point between specified extents of the dimensions.           */     
/*          This function is usually slower for wide extents compared to function cycle_for_near.      */
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

  int naxe = 1;
  node *now_node;

  *num_found = 0;
  *near_dist = DBL_MAX;
  *near_elem = 0; 

  now_node = tree_nodes;
  find_near_rcur(tree_dl, tree_nt, tree_numd, now_node, naxe,
                 extent_min, extent_max, target, 
                 near_elem, near_dist, num_found);

  return;
}

/* --------------------------------------------------------------------------------------------------- */

void find_near_rcur (double **tree_dl, int *tree_nt, int tree_numd, node * now_node, int naxe, 
                     double *extent_min, double *extent_max, double *target, 
                     unsigned long long *near_elem, double *near_dist, unsigned long long *num_found) {

/*          This function finds a nearest point between specified extents of the dimensions.           */     
/*          This function is usually slower for wide extents compared to function cycle_for_near.      */
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
/*          This function is usually faster for wide extents compared to function find_near.           */
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
/* init_rad      Initially search within radius init_rad. If no point within extents is found,         */     
/*               init_rad is multiplied by rad_scal and the seach is repeated. And so on.              */     
/* rad_scal      Scale value to apply to init_rad until a point is found within extents, or the        */     
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

/* --------------------------------------------------------------------------------------------------- */

void brute_near (double **tree_dl, unsigned long long tree_numc, int *tree_nt, int tree_numd,
                 double *extent_min, double *extent_max, double *target, 
                 unsigned long long *near_elem, double *near_dist, unsigned long long *num_found) {

/*          This function finds a nearest point between specified extents of the dimensions.           */     
/*          This function is usually much slower than cycle_for_near and find_near.                    */
/*          The original purpose of this function was to confirm that cycle_for_near and find_near     */
/*          worked as expected (modified kdtree searches are nothing to take for granted).             */
/*          However, for a small number of points (tree_numc) this function will also be faster.       */
/*                                                                                                     */     
/* Input arguments:                                                                                    */     
/* tree_dl       Pointers to first element of values of each dimension. That is, each coordinate is in */     
/*               a separate contiguous array. These are the pointers to the start of those arrays.     */     
/*               Those arrays have size tree_numc.                                                     */     
/* tree_numc     Number of points.                                                                     */     
/* tree_nt       Near type flags. 1 means standard pythagorean nearest. See note for what 0 means.     */     
/* tree_numd     Number of pointers in tree_dl and flags in tree_nt (i.e. number of dimensions).       */     
/* extent_min    Minimum value of this dimension to find nearest point. Size tree_numd.                */     
/*               Greater OR EQUAL to this value is in range.                                           */     
/* extent_max    Maximum value of this dimension to find nearest point. Size tree_numd.                */     
/*               Strictly LESS than this value is in range.                                            */     
/* target        Find the point nearest here considering extents and tree_nt. Size tree_numd.          */     
/*                                                                                                     */     
/* Input/Output arguments:                                                                             */     
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

  *num_found = 0;
  *near_dist = DBL_MAX;
  *near_elem = 0; 

  bool in = true;          
  int i = 0;
  double rad = 0.;
 
  unsigned long long ibrute = 0;

  for(ibrute=0; ibrute<tree_numc; ibrute++) {

    in = true;          
    for(i=0; i<tree_numd; i++) {
      if(tree_dl[i][ibrute] < extent_min[i] || tree_dl[i][ibrute] >= extent_max[i]) { 
        in = false;
        break;
      }
    }

    if(in) { 
      rad = 0.;
      for(i=0; i<tree_numd; i++) {
        if(tree_nt[i] != 0) {
          rad += (target[i]-tree_dl[i][ibrute]) * (target[i]-tree_dl[i][ibrute]);
        }
      }

      if(rad <= *near_dist) {
        if(rad < *near_dist) { /* if distance is smaller, reset count to 1 */
          *num_found = 1;
          *near_elem = ibrute;
          *near_dist = rad; 
        }
        else { /* so, distances are equal. Increment count, set to higher elem.      */
          *num_found = *num_found + 1;                 
          if(ibrute > *near_elem) *near_elem = ibrute;                                   
        }
      }

    }

  } /* end of  for(ibrute=0; ibrute<tree_numc; ibrute++) { */

  return;

}
