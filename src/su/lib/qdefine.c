/*  #include "segy.h" */

#include "su.h"
#include "qdefine.h"

/* - The subinqcsv program uses these routines and has extensive examples. - */
/* - It is easier to understand these routines and the related routines    - */
/* - for bilinear interpolation by running those examples (they produce    - */
/* - q-files, which are csv files that can be easily viewed).              - */
/*                                                                           */
/*                                                                           */
/* The getviacommand function reads values from command line parameters and  */
/* stores them into an array of QInfo structures. Typical stored values are  */
/* time.velocity pairs at cdps (velans). The command line parameters are     */
/* similar to those found in smnmo, sunmocsv, sumute, and sumutecsv programs.*/
/* But the concept is extended in several ways. Rather than just pairs of    */
/* values, this function can read-in and store triplets (and any other       */
/* number of values in a set; the sets are refered to herein as tuples).     */
/* And non-tuple values can also be read-in and stored along with the tuple  */
/* values. For instance, each cdp can also store a static, and an elevation, */
/* and so on.                                                                */
/*                                                                           */
/* Input Arguments:                                                          */
/* pnameA    - names of command line parameters                              */
/* numpnameA - number of names input in pnameA                               */
/* iztupleA  - first tuple name within pnameA                                */
/* numdind   - number of indepenent dimension values that will be output.    */
/*             This value results in better memory management since the      */
/*             stored input values can often be replaced later with          */
/*             intermediate values provided sufficient memory exists.        */
/*             If you have an unusual situation you might want to set this   */
/*             to a negative value (check the code herein).                  */
/*                                                                           */
/* Output Arguments:                                                         */
/* ktuple    - type of tuples (2=pairs, 3=triplets, etc..)                   */
/* ifixd     - =0 is varying (different number tuples per location)          */
/*             =1 is fixed (same number tuples per location, see pindepaA)   */
/*             =2 is none (zero tuples per location)                         */
/* RecInfoA  - array of structures containing all values for each location   */
/* ncdpA     - number of locations output in RecInfoA (array size)           */
/* pindepaA  - if ifixd=1, contains the one-and-only set of independent      */
/*             dimension values, otherwise the varying independent dimension */
/*             values are stored in RecInfoA along with everything else      */
/* ndimsA    - name of the independent dimension.                            */
/* errwarn   - >0 means nothing useful has been accomplished. Details:       */
/*           =1 getviacommand error: no non-tuple name passed in.            */
/*           =2 getviacommand error: non-tuple names have different amounts. */
/*           =3 getviacommand error: independent dimension parameter empty.  */
/*           =4 getviacommand error: an independent dimension paramter empty.*/
/*           =5 getviacommand error: elements of tuple have different amounts*/
/*              at same location                                             */
/*           >0 getviacommand error: returned unrecognized error (errwarn)   */

void getviacommand(cwp_String **pnameA, int *numpnameA, int *iztupleA, int numdind,
                   int *ktuple, int *ifixd, struct QInfo **RecInfoA, int *ncdpA,
                   double **pindepaA, cwp_String **ndimsA, int *errwarn) {

  cwp_String *pname = NULL;                                        
  cwp_String *andim = NULL;   
  int numpname = 0;
  int iztuple = 0;
  struct QInfo *RecInfo = NULL;
  int jtuple = 0;       
  int iqnto = 0;
  int i = 0;
  int k = 0;
  int ncdp = 0;
  int icdp = 0;
  double *pindepa = NULL;
  double *dvals = NULL;

  pname = *pnameA;                                              
  numpname = *numpnameA;
  iztuple = *iztupleA;

  *errwarn = 0;

  ncdp = countparval(pname[0]);

  if(iztuple<1 || numpname<iztuple || ncdp<1) { 
    *errwarn = 1;
    return;
  } 

  for (i=1; i<iztuple; ++i) {
    if(ncdp != countnparval(1,pname[i])) {
      *errwarn = 2;
      return;
    }
  }

  andim = ealloc1(1,sizeof(cwp_String *));

  dvals = ealloc1double(ncdp); 
  RecInfo = ealloc1(ncdp,sizeof(struct QInfo));

/* If no tupled values at all, set flag ifixd=2 and set other things */
/* logically/conveniently to get through the code and error checks.  */
/*   - iztuple was set to be 1 more than the last non-tuple value    */
/*     because it is defined as the first element of first tuple.    */
/*   - set ktuple same as jtuple just on general principle           */

  if(iztuple==numpname) {
    *ifixd = 2;    
    iqnto =  0;       
    jtuple = 1;     
    *ktuple = 1;        
    if(numdind>0) numdind = iztuple-1; 
    else          numdind = 0 - numdind;
  }
  else {

    *ifixd = 0;
    *ktuple = numpname - iztuple;
    jtuple = *ktuple; 

/* Here is why numdind is an input argument. It makes sure we have enough memory to  */
/* store values later AFTER interpolation to constant independent dimension values.  */
/* For example, for velocity analysis there will varying numbers of tims,vels pairs  */
/* for input locations (4 or 7 or 16 pairs, ...). But all of them are going to be    */
/* interpolated to same amount (say, a velocity every 100 milliseconds). The numdind */
/* value allows us to REPLACE the varying input values with the interpolated values. */
/* Note that it is ktuple-1 here because output never contains independent dimension */
/* in the q-records (it is always a fixed set of values on C_SU_NDIMS record).       */

    if(numdind>0) numdind = iztuple+numdind*(*ktuple-1);
    else          numdind = 0 - numdind;
 
/* Is there only 1 set of input independent dimension values?  */
/* If so, read the values, copy name to andim[0] for output    */
/* and then remove it from pname list. Note this is equivalent */
/* situation as the logic for the C_SU_NDIMS record in q-files.*/

    if(countparname(pname[iztuple])==1) {   
      *ifixd = 1;
      jtuple = *ktuple - 1; 
      iqnto = countnparval(1,pname[iztuple]);
      if (iqnto<1) {
        *errwarn = 3;
        return;
      }
      pindepa = ealloc1double(iqnto);
      getnpardouble(1,pname[iztuple],pindepa);  
      andim[0] = pname[iztuple];
      for (i=iztuple+1; i<numpname; ++i) pname[i-1] = pname[i];
      numpname--;
    }
    else {
      andim[0] = ealloc1(4,1);
      strcpy(andim[0],"vary");
    }

  } /* end of  if(iztuple==numpname) { */

/* ------------------------------------------------------------------  */

  for (icdp=0; icdp<ncdp; ++icdp) { 

    if(*ifixd==0) {
      iqnto = countnparval(icdp+1,pname[iztuple]);
      if (iqnto<1) {
        *errwarn = 4;
        return;
      }
    }
    RecInfo[icdp].nto = iqnto;


    if(numdind>iztuple+iqnto*jtuple) {
      RecInfo[icdp].dlots = ealloc1double(numdind);
    }
    else {
      RecInfo[icdp].dlots = ealloc1double(iztuple+iqnto*jtuple);
    }

    if(*ifixd!=2) {
      if(*ifixd==0) 
        getnpardouble(icdp+1,pname[iztuple],RecInfo[icdp].dlots+iztuple+(jtuple-1)*iqnto); 
      for (k=0,i=jtuple-1; k<jtuple; ++k,i--) {
        if (iqnto != countnparval(icdp+1,pname[iztuple+i])) {
          *errwarn = 5;
          return;
        }
        getnpardouble(icdp+1,pname[iztuple+i],RecInfo[icdp].dlots+iztuple+k*iqnto); 
      }
    }

  }

  for (i=0; i<iztuple; ++i) {
    getpardouble(pname[i],dvals);
    for (icdp=0; icdp<ncdp; ++icdp) RecInfo[icdp].dlots[i] = dvals[icdp];
  }

  *ndimsA   = andim;
  *iztupleA = iztuple; 
  *ncdpA    = ncdp;
  *pindepaA = pindepa;                                           
  *RecInfoA = RecInfo;

}

/* The getviaqfile function reads values from a standard q-file and stores   */
/* them into an array of QInfo structures. Some typical stored values are    */
/* time.velocity pairs at cdps (velans). The values in the q-file are        */
/* similar to those found in smnmo, sunmocsv, sumute, and sumutecsv programs.*/
/* But the concept is extended in several ways. Rather than just pairs of    */
/* values, this function can read-in and store triplets (and any other       */
/* number of values in a set; the sets are refered to herein as tuples).     */
/* And non-tuple values can also be read-in and stored along with the tuple  */
/* values. For instance, each cdp can also store a static, and an elevation, */
/* and so on.                                                                */
/*                                                                           */
/* Input Arguments:                                                          */
/* fpP       - an already open file pointer to a standard q-file.            */
/*             (Note: I defined q-file standard when I wrote this Dec 2021). */
/* pnameA    - names of values to determine what not to input from q-file.   */
/* numpnameA - absolute value is the number of names input in pnameA.        */
/*             If numpnameA is 0 or negative then only store a q-file value  */
/*             if it is NOT a name input in pnameA.                          */
/*             If numpnameA is positive then only store a q-file value if    */
/*             it IS a name input in pnameA.                                 */
/* numdind   - number of indepenent dimension values that will be output.    */
/*             This value results in better memory management since the      */
/*             stored input values can often be replaced later with          */
/*             intermediate values provided sufficient memory exists.        */
/*             If you have an unusual situation you might want to set this   */
/*             to a negative value (check the code herein).                  */
/*                                                                           */
/* Output Arguments:                                                         */
/* pnameA    - names of values from q-file (that were not eliminated         */
/*             because of the input pnameA,numpnameA options).               */
/* numpnameA - number of names output in pnameA                              */
/* iztupleA  - first tuple name within output pnameA                         */
/* ktuple    - type of tuples (2=pairs, 3=triplets, etc..)                   */
/* ifixd     - =0 is varying (different number tuples per location)          */
/*             =1 is fixed (same number tuples per location, see pindepaA)   */
/*             =2 is none (zero tuples per location)                         */
/* RecInfoA  - array of structures containing all values for each location   */
/* ncdpA     - number of locations output in RecInfoA (array size)           */
/* pindepaA  - if ifixd=1, contains the one-and-only set of independent      */
/*             dimension values, otherwise the varying independent dimension */
/*             values are stored in RecInfoA along with everything else      */
/* ndimsA    - parced fields from C_SU_NDIMS record. ndimsA[0] is the name   */
/*             of the independent dimension (or vary) and other fields are   */
/*             read as numbers and returned in pindepaA (if appropriate).    */
/* errwarn   - >0 means nothing useful has been accomplished. Details:       */
/*           =1 getqinfo error: extra C_SU_NAMES record in q-file            */
/*           =2 getqinfo error: extra C_SU_NDIMS record in q-file            */
/*           =3 getqinfo error: C_SU_ID record not found immediately         */
/*              after C_SU_NAMES record                                      */
/*          =11 readqhead error: if C_SU_NDIMS not vary, its numbers         */
/*              must align with C_SU_NAMES                                   */
/*          =12 readqhead error: C_SU_ID record not found immediately        */
/*              after C_SU_NAMES record.                                     */
/*          =22 getviaqfile error: C_SU_NDIMS record is not the same length  */
/*              as (the record after) the C_SU_NAMES record.                 */
/*          =23 getviaqfile error: C_SU_NAMES tupled names out-of-order, bad */
/*          =24 getviaqfile error: C_SU_NDIMS blank where valid number wanted*/  
/*          =25 getviaqfile error: C_SU_NDIMS non-number where valid number  */
/*              wanted                                                       */
/*          =26 getviaqfile error: C_SU_NDIMS value must be same for         */  
/*              all members of tuple                                         */  
/*          =27 getviaqfile error: C_SU_NAMES record followed                */  
/*              by C_SU_ID record not found                                  */  
/*         >100 getviaqfile error: at record (errwarn-100) wrong comma count,*/ 
/*              damaged, non-numbers where number wanted                     */ 
/*         >0   getviaqfile error: unrecognized error code (errwarn)         */

void getviaqfile(FILE *fpP, cwp_String **pnameA, int *numpnameA, int *iztupleA, 
                 int numdind,   int *ktupleA, int *ifixd, 
                 struct QInfo **RecInfoA, int *ncdpA, 
                 double **pindepaA,  cwp_String **ndimsA, int *errwarn) {

  struct QInfo *RecInfo;
  cwp_String *pname = NULL;                                           
  cwp_String *aname = NULL;                                              
  cwp_String *andim = NULL;   
  int numpname = 0;
  char *textraw = NULL;           
  char *textbeg = NULL;   
  int *nspot = NULL;
  double *dfield = NULL;
  int iztuple = -1;
  int ktuple = 0;
  int kspace = 0;
  int i = 0;   		                                             
  int j = 0;   	                                            
  int k = 0;   		                                             
  int n = 0;   		                                            
  int jtuple = 0;                                                     
  int iqnto = 0;
  int kqnto = 0;
  int numcases = 0;
  double *pindepa = NULL;
  double dtest = 0.;
  int izlast = 0;
  int notblank = 0;
  char rdel = ',';
  int icdp = 0;	                                              
  int irec = 0;
  int ianames = 0;
  int iandims = 0;
  int janames = 0;
  int jandims = 0;
  int iqflag = 1;
  int ibsize = 0;;
  int ivsize = 0;
  int irtotal = 0;
  int igot = 0;
  *ifixd = 0;
  *errwarn=0;
  *ncdpA = 0;

  pname = *pnameA;                                              
  numpname = *numpnameA;

  getqinfo(fpP,&ianames,&iandims,iqflag,&ibsize,&ivsize,&irtotal,errwarn);

  if(*errwarn>0) return;

  if(ianames<2) {
    *errwarn = 27; 
    return;
  }

  RecInfo = ealloc1(irtotal,sizeof(struct QInfo));

  textraw = ealloc1(ibsize,1);
  textbeg = ealloc1(ibsize,1);

  if(ianames>ivsize) ivsize = ianames;
  if(iandims>ivsize) ivsize = iandims;

  aname = ealloc1(ivsize,sizeof(cwp_String *));
  andim = ealloc1(ivsize,sizeof(cwp_String *));
  dfield = ealloc1double(ivsize);
  nspot = ealloc1int(ivsize); 

  readqhead(fpP,textraw,textbeg,ibsize,aname,andim,&janames,&jandims,errwarn); 

  if(*errwarn>0) { 
    *errwarn += 10;
    return;
  }

/* For ALL q-file records, there is always a record id before the first comma*/
/* Arrange to get rid of in aname,andim. Then honor name-eliminate requests. */

  if(strncmp(aname[0],"c_su_id",7) == 0) aname[0] = "null";

/* If names have lead/trail digits (like 01_gx_12), reduce them to root name.*/

  for(i=0; i<janames; i++) {  
    k = -1;
    n = -1;
    for(j=0; j<strlen(aname[i]); j++) {
      if(aname[i][j]=='_') {
        if(k==-1) k = j;
        else n = j;
      } 
    } 
    if(k!=-1 && n!=-1 && n-k>1) {

/* Make sure they are actually lead/trail digits (NOT like zw_no_hole)        */

      for(j=0; j<k; j++) {
        if(aname[i][j]!='0' && aname[i][j]!='1' && aname[i][j]!='2' && aname[i][j]!='3' &&
           aname[i][j]!='4' && aname[i][j]!='5' && aname[i][j]!='6' && aname[i][j]!='7' &&
           aname[i][j]!='8' && aname[i][j]!='9') {
          n = -1;
          break;
        }
      }
      if(n>-1) {
        for(j=n+1; j<strlen(aname[i]); j++) {
          if(aname[i][j]!='0' && aname[i][j]!='1' && aname[i][j]!='2' && aname[i][j]!='3' &&
             aname[i][j]!='4' && aname[i][j]!='5' && aname[i][j]!='6' && aname[i][j]!='7' &&
             aname[i][j]!='8' && aname[i][j]!='9') {
            k = -1;
            break;
          }
        }
      }

      if(k!=-1 && n!=-1) {
        for(j=0; j<n-k-1; j++) {
          aname[i][j] = aname[i][j+k+1];
        }
        aname[i][n-k-1] = '\0';
      } /* end of  if(k!=-1 && n!=-1) { */

    } /* end of  if(k!=-1 && n!=-1 && n-k>1) { */

  } 

  for(i=0; i<janames; i++) { 
    if(strncmp(aname[i],"null",4) != 0) { /*no store if name starts with null*/ 

      if(numpname<1) {
        k = 1;
        for(j=0; j<-numpname; j++) {        /*only store if name not in in-list*/ 
          if(strcmp(aname[i],pname[j]) == 0) {
            k = 0;
            break;
          }
        }
        if(k==1) {
          aname[numcases] = aname[i];
          if(jandims==janames) andim[numcases] = andim[i];
          nspot[numcases] = i;
          numcases++; 
        }
      }
      else {
        n = 0;
        for(j=0; j<numpname; j++) {         /*only store if name is in in-list*/ 
          if(strcmp(aname[i],pname[j]) == 0) {
            n = 1;
            break;
          }
        }
        if(n==1) {
          aname[numcases] = aname[i];
          if(jandims==janames) andim[numcases] = andim[i];
          nspot[numcases] = i;
          numcases++; 
        }
      }

    }
  }

  iztuple = -1; /* location of the first value of all tupled values */
  for(i=0; i<numcases; i++) { 
    for(j=i+1; j<numcases; j++) {
      if(strcmp(aname[j],aname[i])==0) {
        iztuple = i;         /* location of first element of first tupled value */
        ktuple = j-i;        /* how many in each tuple, except, see jtuple and pindepa */
        kspace = nspot[j]-nspot[i]; /* spacing between instances of first tupled value */
        break;
      }
    }
    if(iztuple>-1) break;
  }

/* If no tupled values at all, set flag ifixd=2 and set other things */
/* logically/conveniently to get through the code and error checks.  */
/*   - iztuple is set to be 1 more than the last non-tuple value     */
/*     because it is defined as the first element of first tuple.    */
/*   - set ktuple same as jtuple just on general principle           */

  if(iztuple<0) {
    iztuple = numcases; 
    *ifixd = 2;    
    iqnto =  0;       
    jtuple = 1;     
    ktuple = 1;        
    if(numdind>0) numdind = iztuple-1; 
    else          numdind = 0 - numdind;
  }

  else { /* tupled values exist */

    if(jandims<2 || strcmp(andim[0],"vary")==0) *ifixd = 0;
    else {
      *ifixd = 1;
      if(janames!=jandims) { 
        *errwarn = 22;
        return;
      }
    }

/* How many in each member of tuple? Note that for ifixd==0 this is not the  */
/* actual,used value of iqnto (it is not fixed). But the input names must be */
/* a representive set. That is, to determine what the tuples are, there must */
/* be at least two sets of names (the first set of names and then a duplicte */
/* set of the same names in same order).                                     */

    iqnto = 0;
    for(i=0; i<numcases; i++) { 
      if(strcmp(aname[i],aname[iztuple])==0) iqnto++; 
    }
    for(k=0; k<ktuple; k++) {
      kqnto = 0;
      for(i=iztuple+k; i<numcases; i+=ktuple) {
        if(strcmp(aname[iztuple+k],aname[i])==0) kqnto++;
      }
      if(kqnto != iqnto) {
        *errwarn = 23;
        return;
      }
    }

/* Using the values from the first representive set, continue nspot values up to*/
/* the maximum allocated size (ivsize, which is maximum number of fields found  */
/* in any record of the q-file). Note that this is done for both ifixd=1 and =0 */
/* even though ifixd=1 should not need it. Note also the use of kspace which is */
/* not the same as ktuple because some names in the input q-records may have    */
/* been in the list of names not to store - in other words maybe we want to     */
/* output 2-tuples, but input q-records have 3-tuples (need kspace=3 for that). */

    for(i=iztuple+ktuple; i<ivsize; i+=ktuple) { 
      for(k=0; k<ktuple; k++) {
        nspot[i+k] = nspot[i+k-ktuple] + kspace;
        numcases = i+k; /* set to last in nspot (just incase of trailing ivsize values) */
      }
    }

/* Note *ifixd==0 means independent dimension values are input on q-records. So values  */
/* can change, and so can how many exist on each individual q-record. For *ifixd==1,    */
/* their values here were read from the single C_SU_NDIMS record. Allocate, read, store */
/* them at pindepa. This means we will not be storing a bunch of copies in each RecInfo.*/
/* Note that when the values come from individual q-records, they are stored on end of  */
/* RecInfo[n].dlots, making it easier not to allocate memory in dlots for *ifixd==1.    */
 
    if(*ifixd==1) ktuple += 1;
    jtuple = ktuple;

/* Here is why numdind is an input argument. It makes sure we have enough memory to  */
/* store values later AFTER interpolation to constant independent dimension values.  */
/* For example, for velocity analysis there will varying numbers of tims,vels pairs  */
/* for input locations (4 or 7 or 16 pairs, ...). But all of them are going to be    */
/* interpolated to same amount (say, a velocity every 100 milliseconds). The numdind */
/* value allows us to REPLACE the varying input values with the interpolated values. */
/* Note that it is ktuple-1 here because output never contains independent dimension */
/* in the q-records (it is always a fixed set of values on C_SU_NDIMS record).       */

    if(numdind>0) numdind = iztuple+numdind*(ktuple-1);
    else          numdind = 0 - numdind;

/* Use and enforce C_SU_NDIMS standards here. For ifixd=1, each tupled value in q-record */
/* gets its own value in C_SU_NDIMS. For instance, first tims value will be at iztuple,  */
/* which is the same comma location as the first vels comma location in q-records.       */
/* When the q-records contain something like vels,numb pairs then the C_SU_NDIMS must    */
/* contain pairs of tims. The second tims value in the pair is allowed to be null or     */
/* blank or equal to the first times value. Similar if there are 3 members in tuples.    */

    if(*ifixd==1) { 
      jtuple = ktuple - 1; 
      pindepa = ealloc1double(iqnto);

      for(j=0; j<iqnto; j++) {
        if(strcmp(andim[iztuple+j*jtuple],"") == 0) { 
          *errwarn = 24;
          return;
        }

        igot = sscanf(andim[iztuple+j*jtuple],"%lf",pindepa+j);  
        if(igot<1) {
          *errwarn = 25;
          return;
        }

        for (k=1; k<jtuple; ++k) { 
          if(strncmp(andim[iztuple+j*jtuple+k],"null",4) != 0 && 
             strcmp(andim[iztuple+j*jtuple+k],"") != 0) { 

            igot = sscanf(andim[iztuple+j*jtuple+k],"%lf",&dtest);  
            if(igot<1) {
              *errwarn = 25;
              return;
            }
            if(pindepa[j] != dtest) {
              *errwarn = 26;
              return;
            }
          }
        }
      }
    }

  } /* end of   if(iztuple<0) { */

/* -------------------------------------------------------  */

  izlast = 0;
  rdel = ',';
  irec = 99; /* start above maximum errwarn value above */

  icdp = 0;

  fseek(fpP, 0L, SEEK_SET); /* reposition file to beginning record  */ 

  while (fgets(textraw, ibsize, fpP) != NULL) {   
    irec += 1;
    if(strncmp(textraw,"Q,",2) == 0 || strncmp(textraw,"q,",2) == 0) {

      getqvalscsv(textraw,textbeg,ibsize,rdel,
                  nspot,numcases,dfield,&izlast,&notblank,errwarn);

      if(*errwarn>0) {
         *errwarn = irec;
         return;
      }

      if(*ifixd>0) {                           /* for non-vary or no tuple cases,   */
        if(izlast != iqnto*jtuple+iztuple) {   /* we are assuming blanks are = 0.0  */
          *errwarn = irec;                     /* as long as enough commas exist to */
          return;                              /* make sure record was not cut off  */
        }
      }
      else {                                   /* for vary, blanks are also = 0.0   */
        iqnto = (izlast-iztuple)/jtuple;       /* but last non-blank must be in a   */
        if(iqnto*jtuple+iztuple < notblank) {  /* complete tuple (there must be     */
          *errwarn = irec;                     /* enough commas remaining to make   */
          return;                              /* sure the tuple was not cut off)   */
        }
        iqnto = (notblank-iztuple)/jtuple;
      }

      RecInfo[icdp].nto = iqnto;

      if(numdind>iztuple+iqnto*jtuple) {
        RecInfo[icdp].dlots = ealloc1double(numdind); 
      }
      else {
        RecInfo[icdp].dlots = ealloc1double(iztuple+iqnto*jtuple);
      }

      for (j=0; j<iztuple; ++j) RecInfo[icdp].dlots[j] = dfield[j];

/* The i in the inner loop results in storing values in reverse name order.     */ 
/* This puts independent dimension on end of dlots memory. Find pindepa for why.*/ 

      for (j=0; j<iqnto; ++j) {
        for (k=0,i=jtuple-1; k<jtuple; ++k,i--) { 
          RecInfo[icdp].dlots[iztuple+i*iqnto+j] = dfield[iztuple+j*jtuple+k];
        }
      }

      icdp++;

    }
  } /* end of  while (fgets(textraw, maxtext, fpP) != NULL) { */

  *pnameA    = aname;     
  *numpnameA = iztuple+ktuple;
  *ndimsA    = andim;     
  *iztupleA  = iztuple;
  *ktupleA   = ktuple;  
  *ncdpA     = icdp;
  *pindepaA  = pindepa;                                           
  *RecInfoA  = RecInfo;

}

/* RecInfo needs to be sorted on 2 values (inline and crossline) which are   */
/* stored in kinf[1] and kinf[2] within the QInfo structure. This function   */
/* checks/enforces the restriction that there always be the same specified   */
/* crossline locations ON every specified inline. That is, this function     */
/* makes sure that RecInfo contains 5inlines by 3crosslines or 17inlines by  */
/* 12crosslines or whatever (this is also known as aligned rectangles).      */
/* The unique inline and crossline numbers are returned in separate arrays.  */
/* (Later, binterpfind uses this to binary-search each direction separately  */
/* to find the 4 surrounding locations needed for bilinear interpolation).   */
/*                                                                           */
/* Input Arguments:                                                          */
/* RecInfo - array containing all values stored for each location            */
/*           RecInfo[n].kinf[1] must contain igi numbers.                    */
/*           RecInfo[n].kinf[2] must contain igc numbers.                    */
/*           RecInfo must be sorted by kinf[1] and kinf[2].                  */
/* ncdp    - number of locations in RecInfo (array size)                     */
/*                                                                           */
/* Output Arguments:                                                         */
/* ngi     - array of unique igi values in low-to-high order                 */
/* ngi_tot - number of values in ngi                                         */
/* ngc     - array of unique igc values in low-to-high order                 */
/* ngc_tot - number of values in ngc                                         */
/* errwarn - >0 means input RecInfo does not contain aligned rectangles.     */
/*           The actual returned value is index in RecInfo where issue was   */
/*           found (but that does not necessarily localize the issue much).  */

void qsplit(struct QInfo *RecInfo, int ncdp,
            int **ngi, int *ngi_tot, int **ngc, int *ngc_tot, int *errwarn) {

  int mgi_tot = -1;
  int mgc_tot = 0;
  int jcdp = 0;
  int igc_set = 0;
  int iset_tot= 0;
  int k = 0;
  int *mgi = NULL;
  int *mgc = NULL;

  *errwarn = 0;

  for (jcdp=0; jcdp<ncdp; ++jcdp) {
    if(RecInfo[jcdp].kinf[2] != RecInfo[0].kinf[2]) {
      mgi_tot = jcdp; /* since jcdp starts at 0 */
      break;
    }
  }

  if(mgi_tot==-1) mgi_tot = ncdp; /* just incase all kinf[2] are the same value */

  igc_set = RecInfo[0].kinf[2]; /* igc value of first set, so cannot match next set */
  iset_tot = mgi_tot;
  for (jcdp=mgi_tot; jcdp<ncdp; ++jcdp) {
    if(RecInfo[jcdp-mgi_tot].kinf[1] != RecInfo[jcdp].kinf[1]) {
      *errwarn = jcdp + 1;
      return;
    }
    if(igc_set == RecInfo[jcdp].kinf[2]) iset_tot++;
    else {
      if(iset_tot != mgi_tot) {
        *errwarn = jcdp + 1;
        return;
      }
      igc_set = RecInfo[jcdp].kinf[2];
      iset_tot = 1;
    }
  }

  mgc_tot = ncdp/mgi_tot;

  mgi = ealloc1int(mgi_tot);
  mgc = ealloc1int(mgc_tot);

  for (k=0; k<mgi_tot; ++k) mgi[k] = RecInfo[k].kinf[1];
  for (k=0; k<mgc_tot; ++k) mgc[k] = RecInfo[k*mgi_tot].kinf[2];

  *ngi_tot = mgi_tot;
  *ngc_tot = mgc_tot;
  *ngi = mgi;
  *ngc = mgc;

  return;

}

/* The qelementout function computes the (starting) number and number of     */
/* values for each non-tuple and tuple value in the dlots array and returns  */
/* them in kindx and ksize arrays. The kindx and ksize arrays return their   */
/* values in the same sequence as the OUTPUT pname sequence.                 */
/* Note: The code herein is trivial. Basically, everything up to iztuple     */
/*       just increments by 1, everything greater than or equal to iztuple   */
/*       simply increments by numdind. Why make a routine just for that?     */
/*       See subinqcsv for an example of how this makes it much easier       */
/*       to understand where to access the values in the dlots buffers.      */
/*                                                                           */
/* Input Arguments:                                                          */
/* iztuple - location of first tuple value in output dlots.                  */
/*           (This is an output argument of getviacommand and getviaqfile).  */
/* ktuple  - type of tuple (2=pairs, 3=triplets).                            */
/*           (This is an output argument of getviacommand and getviaqfile).  */
/* numdind - number of output values within each member of a tuple.          */
/*           (This is an input argument of getviacommand and getviaqfile).   */
/*                                                                           */
/* Output Arguments:                                                         */
/* kindx   - pass in any array of size iztuple+ktuple-1 or more.             */
/*           On output, contains the element number within dlots associated  */
/*           with the pnames output by getviacommand and getviaqfile.        */
/* ksize   - pass in any array of size iztuple+ktuple-1 or more.             */
/*           On output, contains number of elements within dlots associated  */
/*           with the pnames output by getviacommand and getviaqfile.        */
/*                                                                           */
void qelementout(int iztuple,int ktuple,int numdind, int *kindx,int *ksize) {
  int k = 0; 
  for(k=0; k<iztuple; k++) { 
    kindx[k] = k; 
    ksize[k] = 1; 
  }
  for(k=0; k<ktuple-1; k++) { 
    kindx[k+iztuple] = iztuple+k*numdind;
    ksize[k+iztuple] = numdind;
  }
}

/*                                                                                  */
/* Get numbers of values from q-file.                                               */
/*                                                                                  */
/* Note this routine uses fseek to reposition to beginning of the file on entry,    */
/* multiple times during analysis, and also just before return. This routine also   */
/* only looks for the C_SU_NAMES, and C_SU_NDIMS records before the                 */
/* first record that starts with Q, (q and a comma).                                */
/*                                                                                  */
/* Input Arguments:                                                                 */
/* fpR         = open FILE* for a Q-file.                                           */
/* iqflag      = 1 means get ibsize, ivsize, and irtotal. This means reading thru   */
/*               all the q-records.                                                 */
/*             = 0 means do not get ibsize, ivsize, and irtotal values.             */
/*               They are returned as 0.                                            */
/*                                                                                  */
/* Output Arguments:                                                                */
/* num_names   = count of fields on (the record after) the C_SU_NAMES record        */
/*               (actually, number of commas plus 1).                               */
/* num_ndims   = count of fields on the C_SU_NDIMS record                           */
/*               (actually, number of commas plus 1).                               */
/*                                                                                  */
/* ibsize      = maximum size (bytes) of any record (ANY record).                   */
/* ivsize      = maximum count of fields on any record (ANY record).                */
/* irtotal     = number of records (ALL records).                                   */
/*                                                                                  */
/* errwarn    >0 means some kind of error                                           */
/*            =1 more than one C_SU_NAMES parameter record.                         */
/*            =2 more than one C_SU_NDIMS parameter record.                         */
/*   Note that having no C_SU_NAMES or C_SU_NDIMS records does not set errwarn>0.   */
/*   Check the num_names, num_ndims return values.                                  */

void getqinfo(FILE *fpR, int *num_names, int *num_ndims, 
              int iqflag, int *ibsize, int *ivsize, int *irtotal, int *errwarn) {

  *num_names = 1;
  *num_ndims = 1;
  *ibsize    = 0;
  *ivsize    = 0;
  *irtotal   = 0;
  *errwarn = 0;

  int maxtext = 10001;
  char textraw[10001]; 
  char textbeg[21];   

  int num_c_su_names = 0;
  int num_c_su_ndims = 0;
  int idone = 0;
  int nextone = 0;
  int n = 0;
  int isiz = 0;
  int icom = 0;

/* Q-records can contain many values. In order to allocate for them  */
/* we need to know how many by examining the records associated with */
/* C_SU_NDIMS, and/or C_SU_NAMES records. Those records              */
/* may be larger than maxtext above. When they are larger, fgets     */
/* will just return with characters up to maxtext-1. So we need to   */
/* call fgets until we get an end-of-line.                           */
/* In other words it may take several fgets to get a complete record.*/

  fseek(fpR, 0L, SEEK_SET);  
  
  idone = 0;
  nextone = 0;
  while (fgets(textraw, maxtext, fpR) != NULL) { 
    if(nextone>0) {
      for(n=0; n<7; n++) textbeg[n] = tolower(textraw[n]);
      if(strncmp(textbeg,"c_su_id",7) != 0) {
        *errwarn = 3;
        return;
      }
      for(n=0; n<maxtext; n++) { 
        if(textraw[n] == '\n' || textraw[n] == '\r') {
          idone = 1;
          break;
        }
        if(textraw[n] == ',') *num_names += 1;
      }
    }

    nextone = 0; 
    for(n=0; n<10; n++) textbeg[n] = tolower(textraw[n]);
    if(strncmp(textbeg,"c_su_names",10) == 0) {
      num_c_su_names++;
      nextone = 1;
    }
    if(strncmp(textbeg,"q,",2) == 0) break;
  }

  fseek(fpR, 0L, SEEK_SET);  

  idone = 0;
  nextone = 0;
  while (fgets(textraw, maxtext, fpR) != NULL) { 
    for(n=0; n<10; n++) textbeg[n] = tolower(textraw[n]);
    if(strncmp(textbeg,"c_su_ndims",10) == 0) {
      num_c_su_ndims++;
      nextone = 1;
    }
    if(nextone>0) {
      for(n=0; n<maxtext; n++) { 
        if(textraw[n] == '\n' || textraw[n] == '\r') {
          idone = 1;
          break;
        }
        if(textraw[n] == ',') *num_ndims += 1;
      }
    }
    if(strncmp(textbeg,"q,",2) == 0) break;
    nextone = 0;
  }

  fseek(fpR, 0L, SEEK_SET);  

  if(iqflag==1) { 
    isiz = 0;
    icom = 0;
    idone = -1;
    while (fgets(textraw, maxtext, fpR) != NULL) { 
      for(n=0; n<maxtext; n++) { 
        if(textraw[n] == ',') icom += 1;
        if(textraw[n] == '\n' || textraw[n] == '\r') {
          *irtotal += 1;
          idone = n + 1;
          break;
        }
      }
      if(idone == -1) { 
        isiz += maxtext;
      }
      else {
        if(isiz+idone > *ibsize) *ibsize = isiz+idone;
        isiz = 0;
        if(icom > *ivsize) *ivsize = icom;
        icom = 0;
        idone = -1;
      }
    }
    *ivsize += 1;
    fseek(fpR, 0L, SEEK_SET);  
  }

  if(num_c_su_names>1) *errwarn = 1; 
  else if(num_c_su_ndims>1) *errwarn = 2;

  return;

}    

/* Input Argument:                                                                  */
/* fpR         = open FILE*                                                         */
/* textraw     = allocated large enough to contain a record (ibsize from getqinfo). */
/* textbeg     = allocated large enough to contain a record (ibsize from getqinfo). */
/* maxtext     = size of textraw and textbeg                                        */
/*                                                                                  */
/* Output Arguments:                                                                */
/* names       = names read from (the record after) the C_SU_NAMES record           */
/* ndims       = values read from C_SU_NDIMS record                                 */
/* errwarn    >0 means some kind of error                                           */
/*            =1 more than one C_SU_NAMES parameter record.                         */
/*            =2 does not have a C_SU_ID record right after C_SU_NAMES record.      */

void readqhead(FILE *fpR, char *textraw, char *textbeg, int maxtext, 
               cwp_String *names, cwp_String *ndims, 
               int *janames,int *jandims,int *errwarn) {

  *errwarn = 0;
  int num_names = 0;
  int num_ndims = 0;
  int num_c_su_names = 0;
  int num_c_su_ndims = 0;
  int read_names = 0;
  int read_ndims = 0;
  int n = 0;
  int tsize = 0;

  while (fgets(textraw, maxtext, fpR) != NULL) { /* read a line */

/* Stop this looping? */
                   
    if(read_names==-1 && read_ndims==-1) break;

/* Remove all blanks and tabs because tparseq is not designed to handle them.*/

    tsize = 0;
    for(n=0; n<maxtext; n++) { /*   linux \n            windows \r */
      if(textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') break;
      if(textraw[n] != ' ' && textraw[n] != '\t') {
        textbeg[tsize] = textraw[n];
        tsize++;
      }
    }

/* Careful here, sizeof(textbeg) is pointer size, not textbeg string size  */ 

    for(n=0; n<10; n++) textbeg[n] = tolower(textbeg[n]); 
    if(strncmp(textbeg,"c_su_names",10) == 0) num_c_su_names++;
    if(strncmp(textbeg,"c_su_ndims",10) == 0) num_c_su_ndims++;

    if(read_names>0) {
      if(strncmp(textbeg,"c_su_id",7) != 0) {
        *errwarn = 2;
        return;
      }
      textbeg[tsize] = '\0'; 
      tparseq(textbeg, ',', names, &num_names) ; 
      read_names = -1;
    }
    if(strncmp(textbeg,"c_su_names",10) == 0) read_names = 1;

    if(strncmp(textbeg,"c_su_ndims",10) == 0) {
      textbeg[tsize] = '\0';  
      tparseq(textbeg, ',', ndims, &num_ndims) ; 
      read_ndims = -1;
    }

    if(strncmp(textbeg,"q,",2) == 0) break; /* also stop at first q-record */

  } /* end of while (fgets(textraw,..... */
  
  fseek(fpR, 0L, SEEK_SET); /* reposition file to beginning record */ 

  *janames = num_names;
  *jandims = num_ndims;

  if(num_c_su_names>1) *errwarn = 1;
  
  return;
}    



void getqvalscsv(char *textraw, char *textbeg, int maxtext, char rdel, 
                 int *nspot, int numcases, double *dfield, int *izlast, int *notblank, int *errwarn) {

  *errwarn = 0;
  int nbeg = -1;
  int nfield = 0;
  int ineed  = 0;
  int igot = 0;
  int n = 0;
  int m = 0;
  int nb = 0;
  int ib = -1;
  double dval;
  *notblank = -1;

  for(n=0; n<maxtext; n++) {                            /* linux \n            windows \r */
    if(textraw[n] == rdel || textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') {
      if(nfield == nspot[ineed]) {
        dval = 1.1e308; 
        nb = -1;
        if(n-nbeg-1 > 0) {
          strncpy(textbeg,textraw+nbeg+1,n-nbeg-1);
          textbeg[n-nbeg-1] = '\0'; /* so sscanf knows where to stop */
          ib = -1;
          for (m=0; m<n-nbeg-1; m++) {
            if(textbeg[m] != ' ') {
              nb = m;
              if(ib>-1) {
                *errwarn = 1;  
                return;
              }
            }
            if(textbeg[m] == ' ' && nb>-1) ib = m;
          }  

          if(nb>-1) {
            igot = sscanf(textbeg,"%lf",&dval);  
            if(igot<1) {
              *errwarn = 2;  
              return;
            }
          }
        } /* end of  if(n-nbeg-1 > 0) { */
        if(nb<0) dval = 0.;
        else *notblank = ineed+1;
        dfield[ineed] = dval;
        ineed++;
        if(ineed>=numcases) break; 
      }
      if(textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') {
        if(ineed<numcases) break;
      }
      nbeg = n;
      nfield++;
    }
  } 

  *izlast = ineed; 

  return;
}

/* --------------------------- */
/* expects a string with no blanks and no tabs \t   */
void tparseq(char *tbuf, char d, char **fields, int *numfields) { 

  int nbeg = -1;
  int n = 0;
  *numfields = 0;

  for(n=0; ; n++) {
    if(tbuf[n] == d || tbuf[n] == '\0') {
      if(n-nbeg-1 > 0) {
        fields[*numfields] = ealloc1(n-nbeg-1,1);
        strncpy(fields[*numfields],tbuf+nbeg+1,n-nbeg-1);
      }
      else {
        fields[*numfields] = ealloc1(4,1);
        fields[*numfields][0] = 'n';
        fields[*numfields][1] = 'u';
        fields[*numfields][2] = 'l';
        fields[*numfields][3] = 'l';
/*      strncpy(fields[*numfields],"null",4);  makes compilor unhappy */
      }
      nbeg = n;
      *numfields = *numfields + 1;
    }
    if(tbuf[n] == '\0') break;
  }
 
}

