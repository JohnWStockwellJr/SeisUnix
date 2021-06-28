/* Copyright (c) Colorado School of Mines, 2002.*/
/* All rights reserved.                       */

/*
	FROM THE PUBLISHED VERSION IN GEOPHYSICS, SEPT '90
	Original Header:
	NAME:            SEG2SEGY.C
	VERSION:         V1.0
	CREATION DATE:   6/15/90
	AUTHOR:          Brett Bennett
	PURPOSE:         Conversion of the SEG2-PC standard data set 
		 to SEG-Y 32-BIT format.
	LIMITATIONS:     In this implementation of the conversion only
	the data generated in the IBM byte ordering method can be
	converted.  i.e.  byte ordering must be LSB-MSB Data generated
	by a 68000 pro- cessor, for instance, will not convert.  A flag
	is however in the code for future addition of this capability.
	The SEG-Y file also differs slightly from Barry et al standard.
	Word 88 is defined as "last-trace-flag" rather than "Geophone
	group number of last trace within original file."  This
	simplification greatly reduces the amount of code needed for
	conversion.  Header word 92 is defined here as the
	Shot-Sequence-Number.  This is used to assign each trace group a
	unique number.

	EXTENSIVE REWRITE/HACK AND PORT: 11/95 
	AUTHOR: Ian Kay.  
	-Added code to handle running on big-endian machines.
	-Replaced all non-unix portable code and library calls.
	-Typed all the ints so they are the right types
	    (on DOS int=short, on Unix int=long)
	-removed calls to "ieeeibm" since workstations are ieee floating point
	    format as well.
	-partial addition of structures to handle the keywords and headers.
	but I eventually gave up.

	$Id: seg2segy.c,v 1.6 2019/12/10 18:40:34 john Exp $

	Changes:
	Feb 96.  Trace header was not being written out with swap on bytes above
	89 although info is stored in 92,93,94.  Error in original, fixed.

	Oct 97.  Changed the input to command line, inspired by the seg2segy
	version in CWP/SU package.
	
	Jan 98.  
	fixed malloc problem; not enough space was malloced to hold the 
	string defpath. 
	freed the malloced space when finished. 
	closed the segykeyword file when finished.
	changed curf and lastf to int from short; this should allow digits
	in the file name prefix (but what will happen in DOS?).

	Sep 98.
	Ok.  So no one can figure out the segykeyword file mechanism, 
	so I've changed it.  Now, if the environment variable SEGKEYPATH
	is set , then that file is used for the segykeywords. If that
	isn't set then a system default is used /usr/local/lib/segykeyword
	or "segykeyw.ord" in the local directory in DOS.  Finally, neither
	of those exist, then the default values are built in.  

	NOTE that this means that you can do away with the segykeyword
	file altogether!

	$Log: seg2segy.c,v $
	Revision 1.6  2019/12/10 18:40:34  john
	changed tmpnam to mkstemp

	Revision 1.5  2019/12/10 17:57:25  john
	changed long to int and ld to d

	Revision 1.2  2003/02/06 19:05:45  pm
	fixed sawtooth on 20bit packed data conversion case 3
	
	Revision 1.1  2003/02/06 18:44:59  pm
	working on 3rdparty seismic unix
	
	Revision 1.2  2002/11/12 16:01:43  john
	type changed on i,j,k to in
	t

	Revision 1.1  2001/04/12 17:55:11  john
	Initial revision

	Revision 1.7  1998/10/08 20:17:15  kay
	A little more cleanup...^[

	Revision 1.6  1998/10/07 18:54:40  kay
	Built in the defaults for the segykeyword file.  This file is no longer
	needed.

	Revision 1.5  1998/01/07 13:36:15  kay
	added new changes and a Log for the future

 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef __MSDOS__   
#include <unistd.h>
#define MAXSAMP 16384
#define MAXTRACES 16384
#else
#include <io.h>
#define MAXSAMP 4096
#define MAXTRACES 512
#endif

#ifndef NULL
#define NULL    ((void *)0)
#endif 

#define STRINGWIDTH 100
#define MAXKEYWORDS 100
#define MAXPARMS 10

/* global data... */
short int segyreelheader[1800];
int totalkeys;
char string1[200];
short int outhead[120];


typedef struct _SegyKeyTable
	{
	char segykeyword[STRINGWIDTH];
	int segyfunction;
	int segyheader;
	double segyparms[MAXPARMS];
	} SegyKeyTable;


void
swapshort (short *si)
{
/*      void swab(char *from, char *to, size_t nbytes); */
	short tempc;
	tempc = *si;
	swab ((char *) &tempc, (char *) si, sizeof (short));
}

void
swapint (int *li)
{
	short sl, sh;
	int temp, temp1;
	temp = *li;
	sl = temp & 0x0000ffff;
	sh = (temp >> (8 * sizeof (short))) & 0x0000ffff;
	swapshort (&sl);
	swapshort (&sh);
	temp = (sl << (8 * sizeof (short)));
	temp1 = (sh & 0x0000ffff);
	temp = temp | temp1;

	*li = temp;
}

void
float_to_ibm(int from[], int to[], int n)
{
	register int fconv, fmant, ii, t;

	for (ii=0;ii<n;++ii) {
	fconv = from[ii];
	if (fconv) {
	fmant = (0x007fffff & fconv) | 0x00800000;
	t = (int) ((0x7f800000 & fconv) >> 23) - 126;
	while (t & 0x3) { ++t; fmant >>= 1; }
	fconv = (0x80000000 & fconv) | (((t>>2) + 64) << 24) | fmant; 
	} 
	to[ii] = fconv; 
	} 
	return; 
}

char *
strrev (char *s)
{
	char *tmpstr;
	int i, len;
	tmpstr = (char *) malloc (sizeof (char) * strlen (s));
	len = 0;
	while (s[len] != '\0')
	 len++;
	for (i = 0; i < len; i++) {
	tmpstr[i] = s[len - 1 - i];
	  tmpstr[i + 1] = '\0';
	}
	tmpstr[++i] = '\0';
	strcpy (s, tmpstr);
	free (tmpstr);
	tmpstr = (char *) NULL;
	return s;
}


void
parseSegyKeys(FILE *keyfptr, SegyKeyTable *keywords)
{
	int i, j;
	char input[STRINGWIDTH], inputbuf[STRINGWIDTH];
	char *token;
	i = 0;
	while (fgets (input, STRINGWIDTH, keyfptr)) {
	j = 0;
	if ( strlen(input) > (size_t)STRINGWIDTH) {
	fprintf (stderr, "String too int!\n (line:%d)\n", __LINE__);
	exit (-12);
	}
	/*  
	 * if left most character = "*" then 
	 * this line is a comment and should be ignored. 
	 */
	if (input[0] == '*' || input[0] == '\n')  
	continue;

	strcpy (inputbuf, input);  /* make a working copy of input */
	token = strtok (inputbuf, " ");  /* search out space. */

	/* copy keyword into array */
	strncpy (keywords[i].segykeyword, token, 1 + strlen (token));

	/* get function value. */
	token = strtok (NULL, " ");
	keywords[i].segyfunction = atoi (token);  /* convert to value. */
	token = strtok (NULL, " ");

	/* get segy header pointer */
	keywords[i].segyheader = atoi (token);

	/* now start getting any parameters on the line. */
	token = strtok (NULL, " ");
	j = 0;
	while (token != NULL) {
	/* get next token and start putting them in their array location */
	keywords[i].segyparms[j] = atof (token);
	token = strtok (NULL, " ");
	j++;
	if (j > MAXPARMS) {
	printf ("Too many parameters in %s keyword\n", keywords[i].segykeyword);
	printf ("No more than %d allowed per function\n", j - 1);
	exit (-13);
	}
	}      /* end parameter extraction while loop */
	i++;      /* inc counter to keyword number */
	}        /* end keyword string while loop */
	totalkeys = i;
	
}        /* end parseSegyKeys() */

void
loadSegyKeyTableDefaults(SegyKeyTable * keywords)
{
char *defaults[] = {
"ACQUISITION_DATE 3 6 ",
"ACQUISITION_TIME 3 7 ",
"ALIAS_FILTER 5 0 1 71 72 ",
"CDP_NUMBER 5 1 1 12 ",
"CDP_TRACE 5 1 1 14 ",
"CHANNEL_NUMBER 5 1 1 8 ",
"CLIENT 3 1 ", 
"COMPANY 3 2 ",
"DATUM 5 1 1 28 ",
"DELAY 1 55 1000 ",
"END_OF_GROUP 1 88 1 ",
"FIXED_GAIN 1 61 1 ",
"FIXED_GAIN 1 62 1 ",
"HIGH_CUT_FILTER 5 0 1 76 78 ",
"INSTRUMENT 3 4 ",
"JOB_ID 3 9 ",
"LINE_ID 3 3 ",
"LOW_CUT_FILTER 5 0 1 75 77 ",
"NOTCH_FREQUENCY 5 0 1 73 ",
"OBSERVER 3 5 ",
"RAW_RECORD 5 1 1 6 ",
"RECEIVER_LOCATION 5 1 1 42 44 22 ",
"RECEIVER_STATION_NUMBER 1 93 1 ",
"SAMPLE_INTERVAL 1 59 1000000 ",
"SOURCE_LOCATION 5 1 1 38 40 26 ",
"SHOT_SEQUENCE_NUMBER 1 92 1 ",
"SOURCE_STATION_NUMBER 1 94 1 ",
"STACK 1 16 1 ",
"STATIC_CORRECTIONS 5 0 1 50 51 52 ",
"TRACE_SORT 2 1615 ",
"TRACE_TYPE 4 15 ",
"UNITS 3 8 ",
"AMPLITUDE_RECOVERY 0 0 ",
"BAND_REJECT_FILTER 0 0 ",
"DESCALING_FACTOR 0 0 ",
"DIGITAL_BAND_REJECT_FILTER 0 0 ",
"DIGITAL_HIGH_CUT_FILTER 0 0 ",
"DIGITAL_LOW_CUT_FILTER 0 0 ",
"GENERAL_CONSTANT 0 0 ",
"NOTE 0 0  ",
"POLARITY 0 0 ",
"PROCESSING_DATE 0 0 ",
"PROCESSING_TIME 0 0 ",
"RECEIVER 0 0  ",
"RECEIVER_GEOMETRY 0 0 ",
"RECEIVER_SPECS 0 0 ",
"SKEW 0 0  ",
"SOURCE 0 0 ",
"SOURCE_GEOMETRY 0 0 ",
NULL};
	char **defaultsptr=defaults;

/*
	char tmpfilename[L_tmpnam];
	FILE *tmpfileptr;
#ifdef __MSDOS__
	/* MSDOS doesn't seem to have a tmpnam_s() so continue to use the tmpnam() 
         * function.  tmpnam_s() is available on WIN32 with possible code:
         * errno_t err;
         * err = tmpnam_s(tmpfilename, L_tmpnam);
         * if (err)
         * {
         * printf("Error occurred creating temp file.\n");
         * exit(1);
         * }
        */
        tmpname(tmpfilename);
#else
        mkstemp(tmpfilename);
#endif
	tmpfileptr=fopen(tmpfilename, "w");
*/

	FILE *tmpfileptr;
	char template_name[]="/tmp/cmguiXXXXXX";
	int temp_fd;

	temp_fd=mkstemp(template_name);
	tmpfileptr=fdopen(temp_fd,"w");

	/*fwrite(defaults, sizeof(char), strlen(defaults), tmpfileptr); */

	while(*defaultsptr) fprintf(tmpfileptr,"%s\n", *defaultsptr++);

/*	 tmpfileptr=freopen(tmpfilename,"r",tmpfileptr); */

	tmpfileptr=freopen(template_name,"r",tmpfileptr);

	parseSegyKeys(tmpfileptr, keywords);
	fclose(tmpfileptr);
/*	remove(tmpfilename); */
	
}        /* end loadSegyKeyTableDefaults() */


/*-----------------------------------------------------------------------
	READSEGYKEYS
	This routine reads in the valid keywords in file SEGKEYW.ORD 
	and stores the results in global arrays
----------------------------------------------------------------------*/
void
readSegyKeys (SegyKeyTable *keywords)
{
	char *envkeypath;
	char *syskeypath;
#ifdef __MSDOS__
	char *defpath = "segykeyw.ord";
#else
	char *defpath = "/usr/local/lib/segykeyword";
#endif

	FILE *keyfile;

	syskeypath = (char *) malloc ((strlen(defpath)+1) * sizeof(char));
	sprintf (syskeypath, "%s", defpath);
	envkeypath = getenv ("SEGKEYPATH");

	if (envkeypath != (char *) NULL)
	{
	if((keyfile=fopen(envkeypath,"rb")) != (FILE *)NULL)
	fprintf (stderr, "%s used for header word mapping.\n", envkeypath); 
	parseSegyKeys(keyfile, keywords);
	fclose (keyfile);
	} 
	else if ((keyfile = fopen (syskeypath, "rb")) != (FILE *)NULL)
	{
	fprintf (stderr, "Using header word mappings in %s\n",defpath);
	parseSegyKeys(keyfile, keywords);
	fclose (keyfile);
	} 
	else 
	{ 
	fprintf (stderr, "Using default header word mappings. \n");
	loadSegyKeyTableDefaults(keywords);
	}
	free(syskeypath);
	return;
}        /* end readsegyKeys() */



/*------------------------------------------------------------
 * KEYCHECK
 * This routine compares list of keywords with the global string 1 and inserts
 * values found in the input string into specified header locations (see
 * segykeyw.ord).  KEYCHECK checks EVERY keyword for a match.  Repeated 
 * keywords are valid.
 * GLOBALS USED:
 * totalkeys int
 *------------------------------------------------------------*/
void
keycheck (SegyKeyTable *keywords)
{
	int i;
	int matchfound;
	char string2[STRINGWIDTH];
	char *token;

	/* start a loop for total keys and compare the input with the list */
	/* make a copy of string1. strtok destroys it's input during token search*/
	/* therefore need to keep a copy for each keyword checked */
	strcpy (string2, string1);

	/* clear the match found flag */
	matchfound = 0;
	for (i = 0; i < totalkeys; i++) 
	{
	/* 
	 * restore string1.  It will have been destroyed if an 
	 * earlier match has occured. 
	 */
	strcpy (string1, string2);  

	if (0 == strncmp (string1, keywords[i].segykeyword, 
		strlen (keywords[i].segykeyword) ) )
	{
	/* have found a match! Set the matchfound flag */
	matchfound = 1;
	/* look up the function and implement it. */
	switch (keywords[i].segyfunction)
	{
	case 0:
		break;    /* null case nothing to do */
	case 1:
	{
	/* function 1 keywords have a single parameter 
	 * in the input;  assumed to be a number.  */
	/* find the parameter. */
	/* find first token which will be a keyword */
	token = strtok (string1, " ");  
	/* should be a number. */
	token = strtok (NULL, " ");  	
	/* parameter found. pointed to by token. */
	/* function 1 calls for the value on the input line to be mulitplied */
	/* by segykeyword parm 1. and then inserted into header. */
	outhead[keywords[i].segyheader - 1] = 
		atof (token) * keywords[i].segyparms[0];
	break;
	}
	case 2:
	{
	/* function 2 is a special function.  It has to deal with the special 
	 * case of trace sorting. In this case the parameter on the input 
	 * line is not a number  but a char string.  
	 * Notice that spaces in the keywords (i.e. AS ACQUIRED) cause 
	 * the parsing to fail.  Words must be 'spaceless' */
	token = strtok (string1, " ");  	/* pointing to keyword. */
	token = strtok (NULL, " ");  	/* token points to input char. */
	if (0 == strcmp ("AS_ACQUIRED", token))
		segyreelheader[keywords[i].segyheader - 1] = 1;
	else if (0 == strcmp ("CDP_GATHER", token))
		segyreelheader[keywords[i].segyheader - 1] = 2;
	else if (0 == strcmp ("CDP_STACK", token))
		segyreelheader[keywords[i].segyheader - 1] = 4;
	else if (0 == strcmp ("COMMON_OFFSET", token))
		segyreelheader[keywords[i].segyheader - 1] = 3;
	else if (0 == strcmp ("COMMON_RECEIVER", token))
		segyreelheader[keywords[i].segyheader - 1] = 1;
	else if (0 == strcmp ("COMMON_SOURCE", token))
		segyreelheader[keywords[i].segyheader - 1] = 1;
	else if (0 == strcmp ("METERS",token))
		segyreelheader[keywords[i].segyheader -1] = 1;
	else if (0 == strcmp ("FEET",token))
		segyreelheader[keywords[i].segyheader -1] = 2;

	break;
	}      /* end case 2 */
	case 3:
	{
	/* 
	 * this case requires the text string found on the input
	 * line to be copied to the SEGY-REEL header at the line
	 * indicated in the segyheader index. to compute this
	 * location it is (line#-1)*80
	 */
	strncpy ((char *) &segyreelheader[80 * (keywords[i].segyheader - 1)], 
		string1, 80);
	break;
	}      /* end case 3 */
	case 4:
	{
	/* this case, like 2, calls for special string parsing. */
	/* for TRACE_TYPE, see code for allowable inputs. */
	/* pointing to keyword. */
	token = strtok (string1, " "); 
	/* assume its seismic */
	outhead[keywords[i].segyheader - 1] = 1;  
	/* token points to input char. */
	token = strtok (NULL, " "); 
	if (0 == strcmp ("SEISMIC_DATA", token))
		outhead[keywords[i].segyheader - 1] = 1;
	else if (0 == strcmp ("DEAD", token))
		outhead[keywords[i].segyheader - 1] = 2;
	else if (0 == strcmp ("TEST_DATA", token))
		outhead[keywords[i].segyheader - 1] = 3;
	else if (0 == strcmp ("UPHOLE", token))
		outhead[keywords[i].segyheader - 1] = 5;
	/* RADAR_DATA not defined in SEG-Y assume it to be normal seismic */
	else if (0 == strcmp ("RADAR_DATA", token))
		outhead[keywords[i].segyheader - 1] = 1;
	break;
	}      /* end case 4 */
	case 5:
	{
	/* this case calls for the input data to be inserted into the 
	 * sepcified header location of parms 2-10. The normal 
	 * 'header' location indicates the type of data, 0=int, 
	 * 1=int int, 2=floating point. Parm 1 is the multiplier. */
	token = strtok (string1, " ");
	/* do the integer case first */
	if (keywords[i].segyheader == 0) {
		int paramcount = 1;
		int headindex;
		token = strtok (NULL, " ");    	/* token points to input char. */
		while (token != NULL && paramcount < 10) {
		headindex = keywords[i].segyparms[paramcount] - 1;
		outhead[headindex] = atof (token) * keywords[i].segyparms[0];
		paramcount++;
		token = strtok (NULL, " ");  /* token points next input char. */
		}
	}
		/* do the int integer case next */
	else if (keywords[i].segyheader == 1) {
		int paramcount = 1;
		int headindex;
		int *outpoint;
		token = strtok (NULL, " ");    /* token points to input char. */
		while (token != NULL && paramcount < 10) {
		/* set outpoint to beginning of WORD of interest */
		/* need to subtract 2 from the location value. */
		headindex = keywords[i].segyparms[paramcount] - 2;
		outpoint = (int *) &outhead[headindex];
		/* outpoint now points to the area of interest and is 
		 * typed correctly. Compute the value, cast it and send it in. */
		outpoint[0] = (int) (atof (token) * keywords[i].segyparms[0]);
		paramcount++;
		token = strtok (NULL, " ");  /* token points next input char. */
		}
	}
	/* finally do floating point case */
	else if (keywords[i].segyheader == 2) {
		int paramcount = 1;
		int headindex;
		float *outpoint;
		token = strtok (NULL, " ");    /* token points to input char. */
		while ((token != NULL) && (paramcount < 10))
		{
		/* set outpoint to point to beginning of WORD of interest */
		/* need to subtract 2 from the location value. */
		headindex = keywords[i].segyparms[paramcount] - 2;
		outpoint = (float *) &outhead[headindex];
		/* outpoint now points to the area of interest and is typed 
		 * correctly. Compute the value, cast it and put it in 
		 * outhead (using outpoint). */
		outpoint[0] = (float) (atof (token) * keywords[i].segyparms[0]);
		paramcount++;
		/*   ieee2ibm(outpoint,0); *//* convert to ibm format if necessary */
		token = strtok (NULL, " ");  /* token points next input char. */
		}
	}
	break;
	}      /* end case 5 */
	case 6: 
	{
	short day, year;
	token=strtok(string1," "); 
	token=strtok(NULL,"/"); 
	day=atoi(token);
	token=strtok(NULL,"/");
	if (0==strcmp("FEB",token) || 0==strcmp("02",token)) day+=31;
	else if(0 == strcmp("MAR",token) || 0== strcmp("03",token)) day+=59;
	else if(0 == strcmp("APR",token) || 0== strcmp("04",token)) day+=90;
	else if(0 == strcmp("MAY",token) || 0== strcmp("05",token)) day+=120;
	else if(0 == strcmp("JUN",token) || 0== strcmp("06",token)) day+=151;
	else if(0 == strcmp("JUL",token) || 0== strcmp("07",token)) day+=181;
	else if(0 == strcmp("AUG",token) || 0== strcmp("08",token)) day+=212;
	else if(0 == strcmp("SEP",token) || 0== strcmp("09",token)) day+=243;
	else if(0 == strcmp("OCT",token) || 0== strcmp("10",token)) day+=273;
	else if(0 == strcmp("NOV",token) || 0== strcmp("11",token)) day+=304;
	else if(0 == strcmp("DEC",token) || 0== strcmp("12",token)) day+=334;
		token=strtok(NULL," ");/* token points to input char. */
	year=atoi(token);
	if(!year%4 && day>59) day+=1;  /* Yikes.  This may break! */
	outhead[keywords[i].segyheader - 2] = year;
	outhead[keywords[i].segyheader - 1] = day;
	break; 
	} 
	default:    /* case where function not found.. should never happen */
	{
	printf ("Function %d not defined.\n", keywords[i].segyfunction);
	break;
	}
	}            /* end switch */
	break;        /* don't go through rest of for() loop, go to next string */
	}              /* end if */
		      	/* loop through and see if it can be found again */
	}                /* end of i loop */
	if (!matchfound)      /* did a match occur */
	 	printf ("No match found for %s\n", string1);
}             /* end of keysegy */
 



int
main (int argc, char *argv[])
{

#define NFILECHAR 32
	char  prefix[NFILECHAR], seg2file[NFILECHAR], 
	segyfile[NFILECHAR], suffix[NFILECHAR], str[NFILECHAR];
	char *digits = "1234567890";
	char stringtermcount, stringterm1, stringterm2;
	char linetermcount, lineterm1, lineterm2;
	char reserved[19];
	unsigned char datatype;
	int i, j, k; 
	short int ssn; 
	short int reversed = 0;
	int curf, lastf;
	short int blockid, revnum, pointerbytecount, numtrace, stringlength;
	short int first = 1;
	size_t l,ln;
	unsigned short int blockleng;
	int tracepointers[MAXTRACES];
	int outbuf[MAXSAMP];
	int *outheadl=NULL;
	unsigned int numsamples, datalength;
	FILE *f1=NULL, *f2=NULL;
	SegyKeyTable keywords[MAXKEYWORDS];

	double dinbuf[MAXSAMP];
	float *finbuf = (float *)dinbuf;
	int *linbuf = (int *)dinbuf;
	short int *iinbuf = (short int *)dinbuf;
	unsigned char *cinbuf = (unsigned char *)dinbuf;

	float scale=1;

	outheadl = (int *) &outhead[0];
	for (i = 0; i < 1800; i++) segyreelheader[i] = 0;

	if (argc < 3 || argc > 4) { 
	printf("Usage: seg2segy first-seg2file number-of-files [shot-number]\n");
	exit(-1);
	}
	sprintf(seg2file,"%s", argv[1]);
	if (strchr(seg2file,'.')==NULL)
	strcat(seg2file, ".dat");
	l = strcspn(seg2file,".");
	strncpy(segyfile, seg2file, l); 
	segyfile[l]='\0';
	strcpy(suffix,seg2file+l);
	l=strcspn(segyfile,digits);
	if (l==strlen(segyfile) 
	|| strspn(segyfile+l,digits)!=strlen(segyfile+l)){ 
	printf("file name seg2 %s invalid\n", seg2file);
	exit(-2);
	}
	strncpy(prefix, segyfile, l); prefix[l]='\0';
	curf=atoi(segyfile+l);
	ln=strlen(segyfile+l);
	strcat(segyfile, ".sgy");
	lastf=curf+atoi(argv[2])-1;
	if (argc==4) ssn=atoi(argv[3]); else ssn=1;

	f2 = fopen(segyfile, "wb");
	if (f2 == (FILE *) NULL) {
	fprintf (stderr, "**OUTPUT FILE OPEN FAILURE**\n **ABORTING**\n");
	  exit (-3);
	}

	readSegyKeys (keywords);

	/* start the big loop ... */
	for(; curf<=lastf; curf++) {
	strcpy(seg2file,prefix);
	sprintf(str,"%d",curf);
	l=strlen(str);
	while(l<ln) { strcat(seg2file,"0"); l++; }
	strcat(seg2file,str);
	strcat(seg2file,suffix);
	if ((f1=fopen(seg2file,"rb")) == (FILE *)NULL) {
	fprintf (stderr, "\n***ERROR OPENING FILE %s***\n", seg2file);
	fprintf (stderr, "Skipping to next file number.\n");
	continue;    /* go to end of loop, try next file */
	}
	fread (&blockid, 2, 1, f1);
	if (blockid == 0x553a) 
	reversed = 1;
   
	if (blockid != 0x3a55) {
	if (!reversed) {
	fprintf (stderr, "Not SEG-2 data can not continue\n");
	exit (-4);
	}
	}

	fread (&revnum, 2, 1, f1);
	fread (&pointerbytecount, 2, 1, f1);
	fread (&numtrace, 2, 1, f1);
	if (reversed) {
	swapshort (&revnum);
	swapshort (&pointerbytecount);
	swapshort (&numtrace);
	}

	printf ("File %s, Data Format Revision: %d, Number of traces: %d\n", seg2file, revnum, numtrace);
	fread (&stringtermcount, 1, 1, f1);
	fread (&stringterm1, 1, 1, f1);
	fread (&stringterm2, 1, 1, f1);
	fread (&linetermcount, 1, 1, f1);
	fread (&lineterm1, 1, 1, f1);
	fread (&lineterm2, 1, 1, f1);
	fread (reserved, 1, 18, f1);  /* reserved block, not used */

	if (numtrace > (pointerbytecount / 4)) {
	fprintf (stderr, "Number of traces greater than number of trace pointers\n");
	fprintf (stderr, "Number of pointers = %d\n", pointerbytecount / 4);
	fprintf (stderr, "Due to this inconsistency processing must stop\n");
	exit (-5);
	}
	fread (tracepointers, 4, numtrace, f1);
	if (reversed) {
	for (i = 0; i < numtrace; i++)
	swapint (&tracepointers[i]);
	}

	/* now read file descriptor block.  */
	fread (&stringlength, 2, 1, f1);
	if (reversed)
	swapshort (&stringlength);

	while (0 != stringlength) {
	  fread (string1, 1, stringlength - 2, f1);
	  keycheck (keywords);
	  fread (&stringlength, 2, 1, f1);
	  if (reversed) swapshort (&stringlength);
	}
	for (j = 0; j < numtrace; j++) {
	  for (i = 0; i < 120; i++)
	    outhead[i] = 0;
	  printf ("trace-%d-\r", j + 1); fflush(stdout);
	  fseek (f1, tracepointers[j], SEEK_SET);
	  fread (&blockid, 2, 1, f1);
	  if (reversed)
	    swapshort (&blockid);
	  if (blockid == 0x2244) {
	    /* reversed=1; should already know this */
	    fprintf (stderr, "Opps, I've blown it here.... (line:%d)\n", __LINE__);
	    exit (-6);
	  }
	  if (blockid != 0x4422) {
	    fprintf (stderr, "Not a SEG-2 trace header.  Can not process %x (line %d)\n", blockid, __LINE__);
	    exit (-7);
	  }
	  fread (&blockleng, 2, 1, f1);
	  if (reversed)
	    swapshort ((short *) &blockleng);
	  fread (&datalength, 4, 1, f1);
	  if (reversed)
	    swapint ((int *) &datalength);
	  fread (&numsamples, 4, 1, f1);
	  if (reversed)
	   swapint ((int *) &numsamples);
	if (numsamples >= MAXSAMP){ 
	fprintf(stderr, "Your data contains more samples than I can handle\n");
	exit(-8);
	}
	  fread (&datatype, 1, 1, f1);
	  fprintf(stderr,"Data type = %d \n", (int) datatype);
	  if (datatype > 5 || datatype < 1) {
	    fprintf (stderr,"Data type %d not available/valid\n", (int) datatype);
	    break;
	  }
	  outhead[57] = numsamples;
	  fread (reserved, 1, 19, f1);
	  fread (&stringlength, 2, 1, f1);
	  if (reversed)
	    swapshort (&stringlength);
	  while (0 != stringlength) {
	    fread (string1, 1, stringlength - 2, f1);
	    keycheck (keywords);
	    fread (&stringlength, 2, 1, f1);
	    if (reversed)
	      swapshort (&stringlength);
	  }
	
	  fseek (f1, blockleng + tracepointers[j], SEEK_SET);
	  switch ((int) datatype)
	  {
	case 1:
		fread (iinbuf, 2, (int) numsamples, f1);
		if (reversed) {
		for (i = 0; i < (int) numsamples; i++)
		swapshort (&iinbuf[i]);
		}
		for (k = 0; k < numsamples; k++)
		outbuf[k] = iinbuf[k];
		break;
	case 2:
		fread (linbuf, 4, (int) numsamples, f1);
		if (reversed) {
		for (i = 0; i < (int) numsamples; i++)
		swapint (&linbuf[i]);
		}
		for (k = 0; k < numsamples; k++)
		  outbuf[k] = linbuf[k];
		break;
	case 3:
	{
		unsigned int totalbytes, subpointer;
		unsigned int expo;
		int intdat;
		short int sdata; /*modified version by PM */
		totalbytes = (numsamples * 5) / 2;
		fread (cinbuf, 1, (size_t) totalbytes, f1);
		/*  Original Code SU-36 3rd Party
		for (k = 0; k < (numsamples);) {
		subpointer = (k / 4) * 5;
		expo = (unsigned) iinbuf[subpointer++];
		for (i = 0; i < 4; i++) {
		if (0x8000 & iinbuf[subpointer]) 
		intdat = 0xffff8000;
		else 
		intdat = 0;
		intdat = intdat | ((int) iinbuf[subpointer++] << (0x000f & expo));
		expo >>= 4;
		outbuf[k++] = intdat;
		}
		}
		*/
	/*modified version by P.Michaels <pm@cgiss.boisestate.edu> */
	/*fixes sawtooth conversion error on large negative values */
		for (k = 0; k < (numsamples);) {
		subpointer = (k / 4) * 5;
		expo = (unsigned) iinbuf[subpointer++];
		for (i = 0; i < 4; i++) {
	                sdata=iinbuf[subpointer++];
                     		 if (sdata>0) 
                         	  	{       
                                	intdat=(int) sdata;
                                	intdat=intdat << (0x000f & expo);
                        	   	}       
                	        	else    
                         	{       
                               	 	intdat=(int) -sdata; 
                                	intdat=intdat <<  (0x000f & expo);
                              	  	intdat = -intdat; 
                         	 	} /* endif */ 
                		expo >>= 4;
		outbuf[k++] = intdat;
                	} /* next i */
        	} /* next k */
		/* end of modifications */
		if (reversed) {
		for (i = 0; i < numsamples; i++)
		swapint (&outbuf[i]);
		}
	}
	break;
	case 4:
	fprintf(stderr,"Reading type 4/n");
	fread (finbuf, 4, (int) numsamples, f1);
		if (reversed) {
		int *buf = (int *) dinbuf;
		for (i = 0; i < numsamples; i++)
		swapint (&buf[i]);
		}

		/* scale values */
		for (k = 0; k < numsamples; k++) 
		outbuf[k] = finbuf[k]/scale;

		for (k = 0; k < numsamples; k++) 
		fprintf(stderr,"dataval = %d \n", outbuf[k]); 

		fprintf(stderr," if dataval = 0 then increase the value of scale, recompile\n");
		break;
	case 5:
	fread (dinbuf, 8, (int) numsamples, f1);
		if (reversed) {
		int *buf = (int *) dinbuf;
		for (i = 0; i < numsamples * 2; i++)
		swapint (&buf[i]);
		}
		for (k = 0; k < numsamples; k++)
		outbuf[k] = dinbuf[k];

		break;
	}      /* end switch */
	/* set vertical stack traces=1 */

/* 
	if(outhead[15]==0) outhead[15]=1; 
	if(outheadf[59]==0.0) gain=1.; else gain=outheadf[59];
	for(i=0;i<numsamples;i++)
			outbuf[i]*=gain/outhead[15];
	*/
		/* assign the original field record number as the current file number */
		if (outheadl[2] == 0) outheadl[2] = (int) curf;
		/* set trace type  as  sesmic data */
		if (outhead[14] == 0) outhead[14] = 1;

		/* set last trace flag (modified segy) */
		if (j == numtrace - 1 && outhead[87] == 0) {
		outhead[87] = 1;
		ssn = ssn + 1;
		}

		/* from rec-station-number and source-station-number (93 and 94) */
		/* distance from source to receiver */
		/* outheadl[9]=(int)(abs(outhead[93]-outhead[92])); */
		/* set group for trace one and roll switch position  */
		outhead[85] = outhead[86] = 
		(int) (1 + labs ((int) outhead[92] - outheadl[3]));

		/* special case, execute on first pass only... */
		if (first == 1) {
		first = 0;
		segyreelheader[1606] = numtrace;
		segyreelheader[1608] = outhead[58];
		segyreelheader[1609] = outhead[58];
		segyreelheader[1610] = numsamples;
		segyreelheader[1611] = numsamples;
		segyreelheader[1612] = 2;

		if (!reversed) /* swap only if we are on a little endian  machine */
		{  
			for (k = 1600; k < 1606; k += 2)
			swapint ((int *)&segyreelheader[k]);
		  for (k = 1606; k < 1630; k++)
		    swapshort((short *)&segyreelheader[k]);
		}
		fwrite (segyreelheader, 1, 3600, f2);  /*  create the segy headers */
		}

		if (!reversed) { /* swap only if we are on a little endian  machine */
		/* first swap ints */
		for (k = 0; k < 7; k++)
			swapint((int *)&outheadl[k]);
		for (k = 9; k < 17; k++)
		  swapint((int *)&outheadl[k]);
		for (k = 18; k < 22; k++)
		  swapint((int *)&outheadl[k]);
		/* now swap the shorts */ 
		for (k = 14; k < 18; k++)
		  swapshort((short *)&outhead[k]);
		for (k = 34; k < 36; k++)
		  swapshort((short *)&outhead[k]);
		/* for(k=44;k<90;k++)  *//* error: should have gone beyond 95 word */
		for (k = 44; k < 95; k++)
		  swapshort((short *)&outhead[k]);
		}
		if (120 != (k = fwrite (outhead, 2, 120, f2))) {      /* write header */
		fprintf (stderr,"\nWrite failure during header write\n");
		exit (-9);
		}
		if (!reversed) {
	for (k = 0; k < numsamples; k++)
		swapint ((int *)&outbuf[k]);
	}
	
	if ((int) numsamples != (k = fwrite (outbuf, 4, (int) numsamples, f2))) {
	fprintf (stderr,"Write failure during trace write\n");
	exit (-10);
	}

	}      /* end trace loop */
	fclose (f1);
	outhead[87] = 0;    /* reset last trace flag. */
	}        /* end kk loop */
	fclose(f2);
	return 0;
}        /* end main */
