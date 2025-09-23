/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* PSEPSI: $Revision: 1.10 $ ; $Date: 2011/11/17 00:10:53 $	*/

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include "psplot.h"
#include "par.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 									",
" PSEPSI - add an EPSI formatted preview bitmap to an EPS file		",
" 									",
" psepsi <epsfile >epsifile						",
" 									",
" Note:									",
" This application requires						",
" (1) that gs (the Ghostscript interpreter) exist, and			",
" (2) that the input EPS file contain a BoundingBox and EndComments.	",
" Ghostscript is used to build the preview bitmap, which is then		",
" merged with the input EPS file to make the output EPSI file.		",
" 									",
NULL};
/**************** end self doc ********************************/

/*
 * AUTHOR:  Dave Hale, Colorado School of Mines, 03/01/92
 */

#define MAXLINE 2048

/* variables defined at the end of this file */
extern char *preview[];
extern int lpreview;

int main (int argc, char **argv)
{
	int llx,lly,urx,ury,i,bboxfound=0,pbwritten=0;
	char line[MAXLINE],line2[MAXLINE],pbname[BUFSIZ],*pbnamefp=NULL;
	FILE *gsfp,*infp,*pbfp;

	/* initialize getpar (unnecessary now, since no parameters to get) */
	initargs(argc,argv);
	requestdoc(1);

	/* open a pipe to the Ghostscript interpreter */
	if ((gsfp=popen("gs -dNODISPLAY - >/dev/null","w"))==NULL)
		err("cannot find gs - the Ghostscript interpreter!");
	setbuf(gsfp,NULL);

	/* create temporary file used to hold input */
	if ((infp=tmpfile())==NULL)
		err("cannot create temporary file for input EPS!");
	setbuf(infp,NULL);

	/* copy input to a temporary file, while searching for bbox */
	while (fgets(line,MAXLINE,stdin)!=NULL) {
		if (!bboxfound && line[0]=='%') {
			if (!strncmp(line,"%%BoundingBox",13)) {
				if (strstr(line,"atend")!=NULL) continue;
				if (sscanf(line,"%*s %d %d %d %d",
					&llx,&lly,&urx,&ury)<4)
					err("bad BoundingBox in input EPS!");
				bboxfound=1;
			}
		}
		if (fputs(line,infp)==EOF)
			err("cannot copy input EPS to temporary file!");
	}

	/* ensure bbox was found */
	if (!bboxfound) err("input EPS(?) does not have a BoundingBox!");

	/* create temporary file to hold preview bitmap */
	pbnamefp = temporary_filename(pbname);	
	if ((pbfp=fopen(pbnamefp,"w+"))==NULL)
		err("cannot create temporary file for preview bitmap!");
	setbuf(pbfp,NULL);

	/* tell gs to make the preview bitmap in temporary file */
	for (i=0; i<lpreview; ++i)
		fprintf(gsfp,"%s\n",preview[i]);
	fprintf(gsfp,"%d %d (%s) makepreviewbitmap\n",urx-llx,ury-lly,pbname);

	/* send input EPS to Ghostscript interpreter */
	fprintf(gsfp,"gsave\n");
	fprintf(gsfp,"%d %d translate\n",-llx,-lly);
	efseeko(infp,(off_t) 0,SEEK_SET);
	while (fgets(line,MAXLINE,infp)!=NULL)
		fputs(line,gsfp);
	fprintf(gsfp,"grestore\n");
	fflush(gsfp);
	pclose(gsfp);

	/* insert preview bitmap before EndComments */
	efseeko(infp,(off_t) 0,SEEK_SET);
	efseeko(pbfp,(off_t) 0,SEEK_SET);
	while (fgets(line,MAXLINE,infp)!=NULL) {
		if (!pbwritten && line[0]=='%') {
			if (!strncmp(line,"%%EndComments",13)) {
				while (fgets(line2,MAXLINE,pbfp)!=NULL) {
					if (fputs(line2,stdout)==EOF) {
					err("cannot write EPSI to stdout!");
					}
				}
				pbwritten = 1;
			}
		}
		if (fputs(line,stdout)==EOF)
			err("cannot write EPSI to stdout!");
	}

	/* remove temporary file containing preview bitmap */
	if (remove(pbname)!=0)
		warn("cannot remove temporary file %s!",pbname);

	return 0;
}

/* Ghostscript code that makes the preview bitmap */
char *preview[] = {
"% <width> <height> <filename> makepreviewbitmap",
"%",
"% Replaces the current device with a memory device with the",
"% specified width and height.  When showpage is executed, this",
"% memory device will write an EPSI preview bitmap to a file with",
"% the specified filename.  The file will contain the following:",
"%",
"% %%BeginPreview: width height 1 nlines",
"% <nlines of hex data representing the bitmap>",
"% %%EndPreview",
"%",
"% where nlines is the number of lines of hexadecimal characters",
"% representing the bitmap.  The bitmap may then be merged into",
"% the prolog of an EPSI file.",
"%",
"% Author:  Dave Hale, Colorado School of Mines, 02/29/92",
"",
"% make dictionary",
"/pdict 25 dict def",
"pdict begin",
"",
"	% save the showpage operator",
"	/realshowpage /showpage load def",
"",
"	% define a monochrome palette with black and white reversed",
"	/monopalette [ 1 1 1 rgbcolor 0 0 0 rgbcolor ] def",
"",
"	% the main procedure",
"	/makepreviewbitmap {",
"",
"		% get file, height, and width",
"		(w) file /pfile exch def",
"		/pheight exch def",
"		/pwidth exch def",
"",
"		% make memory device",
"		[1 0 0 -1 0 pheight] pwidth pheight monopalette",
"		makeimagedevice /pdevice exch def",
"",
"		% set the device (which does an erasepage)",
"		pdevice setdevice",
"",
"		% replace definition of showpage",
"		userdict /showpage {pdict begin writepreview end} bind put",
"	} def",
"",
"	% procedure to write the preview bitmap to file",
"	/writepreview {",
"",
"		% length of row in bytes an number of rows in device",
"		/lrow pwidth 7 add 8 idiv def",
"		/nrow pheight def",
"",
"		% length of line and number of lines in preview bitmap",
"		/lline 38 def",
"		/nlineperrow lrow lline 1 sub add lline idiv def",
"		/nline nrow nlineperrow mul def",
"",
"		% write BeginPreview comment",
"		pfile (%%BeginPreview: ) writestring",
"		pfile pwidth write==only pfile ( ) writestring",
"		pfile pheight write==only pfile ( ) writestring",
"		pfile (1 ) writestring",
"		pfile nline write==only pfile (\\n) writestring",
"",
"		% string to contain one row from device",
"		/row lrow string def",
"",
"		% zero last byte in row (which may contain unused bits)",
"		row lrow 1 sub 0 put",
"",
"		% loop over rows in device",
"		0 1 nrow 1 sub {",
"",
"			% get row from device",
"			pdevice exch row copyscanlines pop",
"",
"			% loop over lines within row",
"			0 1 nlineperrow 1 sub {",
"",
"				% each line is a comment",
"				pfile (%) writestring",
"",
"				% index of first byte in line",
"				lline mul",
"",
"				% number of bytes in line; may be less",
"				% than lline if this is last line in row",
"				dup lrow exch sub dup lline lt",
"				exch lline ifelse",
"",
"				% row (before index and count)",
"				row 3 1 roll",
"",
"				% write line",
"				getinterval pfile exch writehexstring",
"				pfile (\\n) writestring",
"			} for",
"		} for",
"",
"		% write EndPreview comment",
"		pfile (%%EndPreview\\n) writestring",
"",
"		% flush file and prepare for next page",
"		pfile flushfile",
"		erasepage initgraphics",
"	} bind def",
"end",
"",
"/makepreviewbitmap { pdict begin makepreviewbitmap end } bind def",
};
int lpreview = sizeof(preview)/sizeof(preview[0]);
