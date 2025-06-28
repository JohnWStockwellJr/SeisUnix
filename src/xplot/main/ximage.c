/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* XIMAGE: $Revision: 1.47 $ ; $Date: 2011/11/21 17:03:51 $	*/

#include "par.h"
#include "xplot.h"
#include <X11/Xatom.h>
#include <X11/keysym.h>

/* ZM: interpolate the amplitude from dataset */
float getamp(float *zz,float f1,float d1,int n1,
             float f2,float d2,int n2,float x1,float x2,int verbose);

/*********************** self documentation **********************/
char *sdoc[] = {
"									",
" XIMAGE - X IMAGE plot of a uniformly-sampled function f(x1,x2)     	",
"									",
" ximage n1= [optional parameters] <binaryfile			        ",
"									",
" X Functionality:							",
" Button 1	Zoom with rubberband box				",
" Button 2	Show mouse (x1,x2) coordinates while pressed		",
" q or Q key	Quit							",
" s key		Save current mouse (x1,x2) location to file		",
" p or P key	Plot current window with pswigb (only from disk files)	",
" a or page up keys		enhance clipping by 10%			",
" c or page down keys		reduce clipping by 10%			",
" up,down,left,right keys	move zoom window by half width/height	",
" i or +(keypad) 		zoom in by factor 2 			",
" o or -(keypad) 		zoom out by factor 2 			",
"									",
" ... change colormap interactively					",
" r	     install next RGB - colormap				",
" R	     install previous RGB - colormap				",
" h	     install next HSV - colormap				",
" H	     install previous HSV - colormap				",
" H	     install previous HSV - colormap				",
" (Move mouse cursor out and back into window for r,R,h,H to take effect)",
" 									",
" Required Parameters:							",
" n1			 number of samples in 1st (fast) dimension	",
"									",
" Optional Parameters:							",
" d1=1.0		 sampling interval in 1st dimension		",
" f1=0.0		 first sample in 1st dimension			",
" n2=all		 number of samples in 2nd (slow) dimension	",
" d2=1.0		 sampling interval in 2nd dimension		",
" f2=0.0		 first sample in 2nd dimension			",
" mpicks=/dev/tty	 file to save mouse picks in			",
" perc=100.0		 percentile used to determine clip		",
" clip=(perc percentile) clip used to determine bclip and wclip		",
" bperc=perc		 percentile for determining black clip value	",
" wperc=100.0-perc	 percentile for determining white clip value	",
" bclip=clip		 data values outside of [bclip,wclip] are clipped",
" wclip=-clip		 data values outside of [bclip,wclip] are clipped",
" balance=0		 bclip & wclip individually			",
"			 =1 set them to the same abs value		",
"			   if specified via perc (avoids colorbar skew)	",
" cmap=hsv\'n\' or rgb\'m\'	\'n\' is a number from 0 to 13		",
"				\'m\' is a number from 0 to 11		",
"				cmap=rgb0 is equal to cmap=gray		",
"				cmap=hsv1 is equal to cmap=hue		",
"				(compatibility to older versions)	",
" legend=0	        =1 display the color scale			",
" units=		unit label for legend				",
" legendfont=times_roman10    font name for title			",
" verbose=1		=1 for info printed on stderr (0 for no info)	",
" xbox=50		x in pixels of upper left corner of window	",
" ybox=50		y in pixels of upper left corner of window	",
" wbox=550		width in pixels of window			",
" hbox=700		height in pixels of window			",
" lwidth=16		colorscale (legend) width in pixels		",
" lheight=hbox/3	colorscale (legend) height in pixels		",
" lx=3			colorscale (legend) x-position in pixels	",
" ly=(hbox-lheight)/3   colorscale (legend) y-position in pixels	",
" x1beg=x1min		value at which axis 1 begins			",
" x1end=x1max		value at which axis 1 ends			",
" d1num=0.0		numbered tic interval on axis 1 (0.0 for automatic)",
" f1num=x1min		first numbered tic on axis 1 (used if d1num not 0.0)",
" n1tic=1		number of tics per numbered tic on axis 1	",
" grid1=none		grid lines on axis 1 - none, dot, dash, or solid",
" label1=		label on axis 1					",
" x2beg=x2min		value at which axis 2 begins			",
" x2end=x2max		value at which axis 2 ends			",
" d2num=0.0		numbered tic interval on axis 2 (0.0 for automatic)",
" f2num=x2min		first numbered tic on axis 2 (used if d2num not 0.0)",
" n2tic=1		number of tics per numbered tic on axis 2	",
" grid2=none		grid lines on axis 2 - none, dot, dash, or solid",
" label2=		label on axis 2					",
" labelfont=Erg14	font name for axes labels			",
" title=		title of plot					",
" titlefont=Rom22	font name for title				",
" windowtitle=ximage	title on window					",
" labelcolor=blue	color for axes labels				",
" titlecolor=red	color for title					",
" gridcolor=blue	color for grid lines				",
" style=seismic	        normal (axis 1 horizontal, axis 2 vertical) or  ",
"			seismic (axis 1 vertical, axis 2 horizontal)	",
" blank=0		This indicates what portion of the lower range  ",
"			to blank out (make the background color).  The  ",
"			value should range from 0 to 1.			",
" plotfile=plotfile.ps  filename for interactive ploting (P)  		",
" curve=curve1,curve2,...  file(s) containing points to draw curve(s)   ",
" npair=n1,n2,n2,...            number(s) of pairs in each file         ",
" curvecolor=color1,color2,...  color(s) for curve(s)                   ",
" curvewidth=w1,w2,...		line width(s) for curve(s)		",
" blockinterp=0       whether to use block interpolation (0=no, 1=yes)  ",
"									",
"									",
" NOTES:								",
" The curve file is an ascii file with the points  specified as x1 x2	",
" pairs separated by a space, one pair to a line.  A \"vector\" of curve",
" files and curve colors may be specified as curvefile=file1,file2,etc. ",
" and curvecolor=color1,color2,etc, and the number of pairs of values   ",
" in each file as npair=npair1,npair2,... .                             ",
"									",
NULL};
/*
 * AUTHOR:  Dave Hale, Colorado School of Mines, 08/09/90
 *
 * Stewart A. Levin, Mobil - Added ps print option
 *
 * Brian Zook, Southwest Research Institute, 6/27/96, added blank option
 *
 * Toralf Foerster, Baltic Sea Research Institute, 9/15/96, new colormaps
 *
 * Berend Scheffers, Delft, colorbar (legend)
 *
 * Brian K. Macy, Phillips Petroleum, 11/27/98, added curve plotting option
 *
 * G.Klein, GEOMAR Kiel, 2004-03-12, added cursor scrolling and
 *                                   interactive change of zoom and clipping.
 *
 * Zhaobo Meng, ConocoPhillips, 12/02/04, added amplitude display
 *
 * Garry Perratt, Geocon, 08/04/05, modified perc handling to center colorbar if balance==1.
 */
/**************** end self doc ********************************/

/* functions defined and used internally */
static void zoomBox (int x, int y, int w, int h,
	int xb, int yb, int wb, int hb,
	int nx, int ix, float x1, float x2,
	int ny, int iy, float y1, float y2,
	int *nxb, int *ixb, float *x1b, float *x2b,
	int *nyb, int *iby, float *y1b, float *y2b);
static unsigned char *newInterpBytes (int n1in, int n2in, unsigned char *bin,
	int n1out, int n2out, int newInterpBytes);
void xMouseLoc(Display *dpy, Window win, XEvent event, int style, Bool show,
	int x, int y, int width, int height,
	float x1begb, float x1endb, float x2begb, float x2endb,
        float *z, float f1, float d1, int n1,float f2, float d2,
        int n2, int verbose);
void xMousePrint(XEvent event, int style, FILE *mpicksfp,
	int x, int y, int width, int height,
	float x1begb, float x1endb, float x2begb, float x2endb);
/* JG... */
void intl2b_block(int nxin, float dxin, float fxin,
				  int nyin, float dyin, float fyin, unsigned char *zin,
				  int nxout, float dxout, float fxout,
				  int nyout, float dyout, float fyout, unsigned char *zout);
/* .... JG */


int
main (int argc,char **argv)
{
	int n1,n2,n1tic,n2tic,
		i1,i2,grid1,grid2,style,
		n1c,n2c,i1beg,i1end,i2beg,i2end,i1c,i2c,
		nz,iz,i1step,i2step,verbose,
		xbox,ybox,wbox,hbox,
		xb,yb,wb,hb,
		x,y,width,height,
		i,j,nx,ny,nxb,nyb,ixb,iyb,
		imageOutOfDate,winwidth=-1,winheight=-1,
		showloc=0,
		balance,
		legend,lwidth,lheight,lx,ly; /* BEREND */

	int blockinterp=0; /* JG */
        int base;
        float fact;
        unsigned char* ptr;

	unsigned long nfloats;
	float labelsize,titlesize,perc,clip,bperc,wperc,bclip,wclip,
		d1,f1,d2,f2,*z,*temp,zscale,zoffset,zi,dx,dy,
		x1beg,x1end,x2beg,x2end,
		x1min,x1max,x2min,x2max,
		d1num,f1num,d2num,f2num,
		x1begb,x1endb,x2begb,x2endb,blank; /* dx,dy added GK */

	unsigned char *cz,*czp,*czb,*czbp,*czbi=NULL;
	char *label1="",*label2="",*title="",*windowtitle="ximage",
		*units="", *legendfont="times_roman10",
		*labelfont="Erg14",*titlefont="Rom22",
		*styles="seismic",*grid1s="none",*grid2s="none",
		*labelcolor="blue",*titlecolor="red",
		*gridcolor="blue",*cmap="",keybuf[256],*mpicks;
	FILE *infp=stdin, *mpicksfp;
	Display *dpy;
	Window win;
	XEvent event;
	KeySym keysym;
	XComposeStatus keystat;
	XImage *image=NULL;
	XImage *image_legend=NULL; /* BEREND */
	unsigned char *data_legend; /* BEREND */
	GC gci;
	int scr;
	unsigned long black,white,pmin,pmax;

        float **x1curve,**x2curve;
        int curve,*npair=NULL,ncurvecolor=0,ncurvewidth=0,*curvewidth=NULL;
        char **curvefile,**curvecolor=NULL;
        FILE *curvefp;

	char *plotfile;         /* filename of plotfile GK */
	int lock=0;		/* lock/unlock zoom while scrolling */
	float mve;		/* distance for scrolling */
	float mvefac=8.;	/* window factor for scrolldistance
	                         * 2=half window size;
				 * 8=one eighths of the window size */
	char  *msg="";		/* message on screen */

	/* initialize getpar */
	initargs(argc,argv);
	requestdoc(1);

	/* get parameters describing 1st dimension sampling */
	if (!getparint("n1", &n1))
		err("Must specify number of samples in 1st dimension!\n");
	d1 = 1.0;  getparfloat("d1",&d1);
	f1 = 0.0;  getparfloat("f1",&f1);
	x1min = (d1>0.0)?f1:f1+(n1-1)*d1;
	x1max = (d1<0.0)?f1:f1+(n1-1)*d1;

	/* get parameters describing 2nd dimension sampling */
	if (!getparint("n2",&n2)) {
		if (efseeko(infp, (off_t) 0, SEEK_END)!=0)
			err("must specify n2 if in a pipe!");

		nfloats = (int) (eftello(infp)/( (off_t) sizeof(float)));
		efseeko(infp, (off_t) 0,SEEK_SET);
		n2 = (int) (nfloats/n1);
	}
	d2 = 1.0;  getparfloat("d2",&d2);
	f2 = 0.0;  getparfloat("f2",&f2);
	x2min = (d2>0.0)?f2:f2+(n2-1)*d2;
	x2max = (d2<0.0)?f2:f2+(n2-1)*d2;

	/* set up file to save postscript plot * GK */
	if (!getparstring("plotfile", &plotfile))  plotfile = "plotfile.ps" ;

	/* set up file to save mouse picks */
	if (!getparstring("mpicks", &mpicks)) mpicks = "/dev/tty";
	mpicksfp = efopen(mpicks, "w");

	/* set up curve plotting */
	if ((curve=countparval("curve"))!=0) {
		curvefile=(char**)ealloc1(curve,sizeof(void*));
		getparstringarray("curve",curvefile);
		if ((x1curve=(float**)malloc(curve*sizeof(void*)))==NULL)
			err("Could not allocate x1curve pointers\n");
		if ((x2curve=(float**)malloc(curve*sizeof(void*)))==NULL)
			err("Could not allocate x2curve pointers\n");
                npair=ealloc1int(curve);
                getparint("npair",npair);
	} else {
		npair=(int *)NULL;
		curvefile=(char **)NULL;
		x1curve=(float **)NULL;
		x2curve=(float **)NULL;
	}
	if ((ncurvecolor=countparval("curvecolor"))<curve) {
		curvecolor=(char**)ealloc1(curve,sizeof(void*));
		if (!getparstringarray("curvecolor",curvecolor)) {
			curvecolor[0]=(char *)cwp_strdup("blue\0");
			ncurvecolor=1;
		}
		for (i=ncurvecolor; i<curve; i++)
			curvecolor[i]=(char *)cwp_strdup(curvecolor[ncurvecolor-1]);
	} else if( ncurvecolor ) {
		curvecolor=(char**)ealloc1(ncurvecolor,sizeof(void*));
		getparstringarray("curvecolor",curvecolor);
	}
	ncurvewidth=countparval("curvewidth");
	curvewidth=ealloc1int(ncurvewidth);
	getparint("curvewidth",curvewidth);
	for (j=0; j<curve; j++) {
		curvefp=efopen(curvefile[j],"r");
		x1curve[j]=ealloc1float(npair[j]);
		x2curve[j]=ealloc1float(npair[j]);
		for (i=0; i<npair[j]; i++) {
			fscanf(curvefp,"%f",&x1curve[j][i]);
			fscanf(curvefp,"%f",&x2curve[j][i]);
		}
		efclose(curvefp);
	}

        if (!getparfloat("d2",&d2)) d2 = 1.0;

	f2 = 0.0;  getparfloat("f2",&f2);
	x2min = (d2>0.0)?f2:f2+(n2-1)*d2;
	x2max = (d2<0.0)?f2:f2+(n2-1)*d2;

	/* read binary data to be plotted */
	nz = n1*n2;
	z = ealloc1float(nz);

	if (fread(z,sizeof(float),nz,infp)!=nz)
		err("error reading input file");

	/* if necessary, determine clips from percentiles */
	if (getparfloat("clip",&clip)) {
		bclip = clip;
		wclip = -clip;
	}
	if ((!getparfloat("bclip",&bclip) || !getparfloat("wclip",&wclip)) &&
		!getparfloat("clip",&clip)) {
		perc = 100.0;  getparfloat("perc",&perc);
		balance=0 ; getparint("balance",&balance);
		temp = ealloc1float(nz);
		/* Modded by GCP to balance bclip & wclip */

		if (balance==0)
			for (iz=0; iz<nz; ++iz) {
				temp[iz] = z[iz];
			} else {
				for (iz=0; iz<nz; ++iz) temp[iz] = abs(z[iz]);
				perc=100.0;
			}

		/* End of modded code */
		if (!getparfloat("bclip",&bclip)) {
			bperc = perc;	getparfloat("bperc",&bperc);
			iz = (nz*bperc/100.0);
			if (iz<0) iz = 0;
			if (iz>nz-1) iz = nz-1;
			qkfind(iz,nz,temp);
			bclip = temp[iz];
		}
		if (!getparfloat("wclip",&wclip)) {
			wperc = 100.0-perc;  getparfloat("wperc",&wperc);
			iz = (nz*wperc/100.0);
			if (iz<0) iz = 0;
			if (iz>nz-1) iz = nz-1;
			qkfind(iz,nz,temp);
			/* Modded by GCP to balance bclip & wclip */
			if (balance==0) wclip = temp[iz];
			else wclip = -1*bclip;
			/* End of modded code */
		}
		free1float(temp);
	}
	verbose = 1;  getparint("verbose",&verbose);
	if (verbose) warn("bclip=%g wclip=%g",bclip,wclip);

	/* get colormap specification */
	if (!(getparstring("cmap",&cmap))) {
		cmap = (char *)alloc1(5,1);
		sprintf(cmap,"%s","gray");
	}

	/* get interpolation style JG */
	if (!(getparint("blockinterp", &blockinterp))) blockinterp=0;

	/* get legend specs BEREND */
	legend = 0; getparint("legend", &legend); /* BEREND */
	getparstring("units", &units); /* BEREND */
	getparstring("legendfont", &legendfont);     /* BEREND */

	blank = 0; getparfloat("blank",&blank);

	/* get axes parameters */
	xbox = 50; getparint("xbox",&xbox);
	ybox = 50; getparint("ybox",&ybox);
	wbox = 550; getparint("wbox",&wbox);
	hbox = 700; getparint("hbox",&hbox);

	/* legend dimensions */
	if (!(getparint("lwidth",&lwidth)))	lwidth = 16;
	if (!(getparint("lheight",&lheight)))	lheight = hbox/3;
	if (!(getparint("lx",&lx)))	lx = 3;
	if (!(getparint("ly",&ly)))	ly = (hbox-lheight)/3;

	x1beg = x1min; getparfloat("x1beg",&x1beg);
	x1end = x1max; getparfloat("x1end",&x1end);
	d1num = 0.0; getparfloat("d1num",&d1num);
	f1num = x1min; getparfloat("f1num",&f1num);
	n1tic = 1; getparint("n1tic",&n1tic);
	getparstring("grid1",&grid1s);
	if (STREQ("dot",grid1s)) grid1 = DOT;
	else if (STREQ("dash",grid1s)) grid1 = DASH;
	else if (STREQ("solid",grid1s)) grid1 = SOLID;
	else grid1 = NONE;
	getparstring("label1",&label1);
	x2beg = x2min; getparfloat("x2beg",&x2beg);
	x2end = x2max; getparfloat("x2end",&x2end);
	d2num = 0.0; getparfloat("d2num",&d2num);
	f2num = 0.0; getparfloat("f2num",&f2num);
	n2tic = 1; getparint("n2tic",&n2tic);
	getparstring("grid2",&grid2s);
	if (STREQ("dot",grid2s)) grid2 = DOT;
	else if (STREQ("dash",grid2s)) grid2 = DASH;
	else if (STREQ("solid",grid2s)) grid2 = SOLID;
	else grid2 = NONE;
	getparstring("label2",&label2);
	getparstring("labelfont",&labelfont);
	labelsize = 18.0; getparfloat("labelsize",&labelsize);
	getparstring("title",&title);
	getparstring("titlefont",&titlefont);
	titlesize = 24.0; getparfloat("titlesize",&titlesize);
	getparstring("style",&styles);
	if (STREQ("normal",styles)) style = NORMAL;
	else style = SEISMIC;
	getparstring("titlecolor",&titlecolor);
	getparstring("labelcolor",&labelcolor);
	getparstring("gridcolor",&gridcolor);
	getparstring("windowtitle",&windowtitle);
        checkpars();

	/* adjust x1beg and x1end to fall on sampled values */
	i1beg = NINT((x1beg-f1)/d1);
	i1beg = MAX(0,MIN(n1-1,i1beg));
	x1beg = f1+i1beg*d1;
	i1end = NINT((x1end-f1)/d1);
	i1end = MAX(0,MIN(n1-1,i1end));
	x1end = f1+i1end*d1;

	/* adjust x2beg and x2end to fall on sampled values */
	i2beg = NINT((x2beg-f2)/d2);
	i2beg = MAX(0,MIN(n2-1,i2beg));
	x2beg = f2+i2beg*d2;
	i2end = NINT((x2end-f2)/d2);
	i2end = MAX(0,MIN(n2-1,i2end));
	x2end = f2+i2end*d2;

	/* allocate space for image bytes */
	n1c = 1+abs(i1end-i1beg);
	n2c = 1+abs(i2end-i2beg);
	cz = ealloc1(n1c*n2c,sizeof(unsigned char));

	/* convert data to be imaged into signed characters */
	zscale = (wclip!=bclip)?255.0/(wclip-bclip):1.0e10;
	zoffset = -bclip*zscale;
	i1step = (i1end>i1beg)?1:-1;
	i2step = (i2end>i2beg)?1:-1;
	if (style==NORMAL) {
		for (i2c=0,i2=i2beg; i2c<n2c; i2c++,i2+=i2step) {
			czp = cz+n1c*n2c-(i2c+1)*n1c;
			for (i1c=0,i1=i1beg; i1c<n1c; i1c++,i1+=i1step) {
				zi = zoffset+z[i1+i2*n1]*zscale;
				if (zi<0.0) zi = 0.0;
				if (zi>255.0) zi = 255.0;
				*czp++ = (unsigned char)zi;
			}
		}
	} else {
		czp = cz;
		for (i1c=0,i1=i1beg; i1c<n1c; i1c++,i1+=i1step) {
			for (i2c=0,i2=i2beg; i2c<n2c; i2c++,i2+=i2step) {
				zi = zoffset+z[i1+i2*n1]*zscale;
				if (zi<0.0) zi = 0.0;
				if (zi>255.0) zi = 255.0;
				*czp++ = (unsigned char)zi;
			}
		}
	}
/*	free1float(z);      keep data for plotting GK */

	/* initialize zoom box parameters */
	dx = (style==NORMAL ? d1 : d2);
	dy = (style==NORMAL ? d2 : d1);
	nxb = nx = (style==NORMAL ? n1c : n2c);
	nyb = ny = (style==NORMAL ? n2c : n1c);
	ixb = iyb = 0;
	czb = cz;
	x1begb = x1beg;	 x1endb = x1end;
	x2begb = x2beg;	 x2endb = x2end;

	/* connect to X server */
	dpy = XOpenDisplay(NULL);
	if (dpy==NULL)
		err("Cannot connect to display %s!\n",XDisplayName(NULL));
	scr = DefaultScreen(dpy);
	black = BlackPixel(dpy,scr);
	white = WhitePixel(dpy,scr);

	/* create window */
	win = xNewWindow(dpy,xbox,ybox,wbox,hbox,(int) black,(int) white,windowtitle);

	/* exit event protocol */
	Atom wm_delete_window = XInternAtom(dpy, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(dpy, win, &wm_delete_window, 1);

	/* backwards compatibility */
	if (STREQ(cmap,"gray")) {
		sprintf(cmap,"%s","rgb0");

	} else if (STREQ(cmap,"hue")) {
		/* free1(cmap); */
		cmap = (char *)alloc1(5,1);
		sprintf(cmap,"%s","hsv1");

	} else  if ((strncmp(cmap,"hsv",3)) && (strncmp(cmap,"rgb",3))){
			if (verbose) warn ("cmap=%s using cmap=gray", cmap);

			/* free1(cmap); */
			cmap = (char *)alloc1(5,1);
       			sprintf (cmap, "%s", "rgb0");
	}

	/* here are the new colormaps				*/
	if (strncmp(cmap, "rgb", 3) == 0)
		XSetWindowColormap(dpy,win,
			xCreateRGBColormap(dpy,win, cmap, verbose));
	else if (strncmp (cmap, "hsv", 3) == 0)
		XSetWindowColormap(dpy,win,
			xCreateHSVColormap(dpy,win, cmap, verbose));

	/* determine min and max pixels from standard colormap */
	pmin = xGetFirstPixel(dpy);
	pmax = xGetLastPixel(dpy);
	if (verbose) warn("pmin=%x,pmax=%x\n",pmin,pmax);
	if(pmax==0L)pmax=255L;

	if (verbose) warn("pmin=%x,pmax=%x\n",pmin,pmax);
	data_legend = (unsigned char *) malloc(lwidth * lheight);

        if( bclip < wclip ){
           base=256;
           fact=-256.0;
        }else{
           base=0;
           fact=256.0;
        }
        ptr = data_legend;
	for (i=0; i<lheight; i++){
           for( j=0; j<lwidth; j++ ){
	      *ptr++ = (unsigned char)
                            (base + (fact*i)/lheight);
           }
           /* fprintf(stderr," %d ",*(ptr-1) ); */
	}

	/* make GC for image */
	gci = XCreateGC(dpy,win,0,NULL);

	/* set normal event mask */
	XSelectInput(dpy,win,
		StructureNotifyMask |
		ExposureMask |
		KeyPressMask |
		PointerMotionMask |
		ButtonPressMask |
		ButtonReleaseMask |
		Button1MotionMask |
		Button2MotionMask);

	/* map window */
	XMapWindow(dpy,win);

	/* determine good size for axes box */
	xSizeAxesBox(dpy,win,
		labelfont,titlefont,style,
		&x,&y,&width,&height);

	/* clear the window */
	XClearWindow(dpy,win);

	/* note that image is out of date */
	imageOutOfDate = 1;

	/* main event loop */
	while(imageOutOfDate|(~imageOutOfDate)/*True*/) {
		XNextEvent(dpy,&event);

		/* if window was resized */
		if (event.type==ConfigureNotify &&
			(event.xconfigure.width!=winwidth ||
			 event.xconfigure.height!=winheight)) {
			winwidth = event.xconfigure.width;
			winheight = event.xconfigure.height;

			/* determine good size for axes box */
			xSizeAxesBox(dpy,win,
				labelfont,titlefont,style,
				&x,&y,&width,&height);

			/* clear the window */
			XClearWindow(dpy,win);

			/* note that image is out of date */
			imageOutOfDate = 1;

		/* else if window exposed */
		} else if (event.type==Expose) {

			/* clear all expose events from queue */
			while (XCheckTypedEvent(dpy,Expose,&event));

			/* if necessary, make new image */
			if (imageOutOfDate) {
				 czbi = newInterpBytes(nxb,nyb,czb,
							width,height,blockinterp);

				if (image!=NULL) XDestroyImage(image);
				image = xNewImage(dpy,pmin,pmax,
					width,height,blank,czbi);

				/* BEREND create image */
				if (legend) {
					if (image_legend!=NULL) XDestroyImage(image_legend);
					image_legend = xNewImage(dpy,pmin,pmax,lwidth,lheight,0,data_legend);
				}

				imageOutOfDate = 0;
			}

			/* draw image (before axes so grid lines visible) */
			XPutImage(dpy,win,gci,image,0,0,x,y,
				image->width,image->height);

			/* BEREND display image */
			if (legend)
				XPutImage(dpy,win,gci,image_legend,
					0,0,lx,y+ly,lwidth,lheight);

			/* BEREND draw legend axes on top of image */
			if (legend)
				xDrawLegendBox(dpy,win,
					lx,y+ly,lwidth,lheight,
					bclip,wclip,units,legendfont,
					labelfont,title,titlefont,
					labelcolor,titlecolor,gridcolor,
					style);

                        /* draw curve on top of image */
			for (i=0; i<curve; i++)
				xDrawCurve(dpy,win,
					   x,y,width,height,
					   x1begb,x1endb,0.0,0.0,
					   x2begb,x2endb,0.0,0.0,
					   x1curve[i],x2curve[i],npair[i],
					   curvecolor[i],curvewidth[i],style);

	                /* draw axes on top of image */
			xDrawAxesBox(dpy,win,
				x,y,width,height,
				x1begb,x1endb,0.0,0.0,
				d1num,f1num,n1tic,grid1,label1,
				x2begb,x2endb,0.0,0.0,
				d2num,f2num,n2tic,grid2,label2,
				labelfont,title,titlefont,
				labelcolor,titlecolor,gridcolor,
				style);

		/* else if key down */
		} else if (event.type==KeyPress) {

			XLookupString(&(event.xkey),keybuf,0,&keysym,&keystat);

                     /*  added moving, clipping and zooming GK */
			if (keysym==XK_s) {
				xMousePrint(event,style, mpicksfp,
					    x,y,width,height,
					    x1begb,x1endb,x2begb,x2endb);

			} else if (keysym==XK_l ) {
				/* set lock */
			     lock = 1 ;
			  if (verbose) warn("zoom lock set  %d\n",lock);

 			} else if (keysym==XK_u ) {
				/* unset lock */
			     lock = 0 ;
			  if (verbose) warn("zoom lock released %d\n",lock);

 			} else if (keysym==XK_Shift_L ) {
			     /* if (verbose)
			     fprintf(stderr,"Shift Left pressed \n");*/
			} else if (keysym==XK_KP_1 || keysym==XK_1 ) {
			     mvefac=1.;
			     fprintf(stderr,"Zoom/Move factor = 1 \n");
			} else if (keysym==XK_KP_2 || keysym==XK_2 ) {
			     mvefac=2.;
			     fprintf(stderr,"Zoom/Move factor = 2 \n");
			} else if (keysym==XK_KP_3 || keysym==XK_3 ) {
			     mvefac=3.;
			     if (verbose)
			     fprintf(stderr,"Zoom/Move factor = 3 \n");
			} else if (keysym==XK_KP_4 || keysym==XK_4 ) {
			     mvefac=4.;
			     if (verbose)
			     fprintf(stderr,"Zoom/Move factor = 4 \n");
			} else if (keysym==XK_KP_5 || keysym==XK_5 ) {
			     mvefac=5.;
			     if (verbose)
			     fprintf(stderr,"Zoom/Move factor = 5 \n");
			} else if (keysym==XK_KP_6 || keysym==XK_6 ) {
			     mvefac=6.;
			     if (verbose)
			     fprintf(stderr,"Zoom/Move factor = 6 \n");
			} else if (keysym==XK_KP_7 || keysym==XK_7 ) {
			     mvefac=7.;
			     if (verbose)
			     fprintf(stderr,"Zoom/Move factor = 7 \n");
			} else if (keysym==XK_KP_8 || keysym==XK_8 ) {
			     mvefac=8.;
			     if (verbose)
			     fprintf(stderr,"Zoom/Move factor = 8\n");
			} else if (keysym==XK_KP_9 || keysym==XK_9 ) {
			     mvefac=9.;
			     if (verbose)
			     fprintf(stderr,"Zoom/Move factor = 9\n");
			} else if (keysym==XK_Left ) {
 			  /* move zoom box to left by half window width */
			  mve = (x2endb - x2begb)/mvefac ;
			  x2begb = x2begb - mve ;
			  x2endb = x2endb - mve ;
			  msg="move ";
			  /* check for bounds of full window */
			  if (x2begb < x2beg){
			    if ( lock ) { x2begb = x2begb + mve ;
			                  x2endb = x2endb + mve ;
					  msg="limit ";
					  mve=0;
			    } else {  x2begb = x2beg ;
			              nxb=(int)((x2endb-x2begb)/dx);
				   }
			  }

			  if (verbose) fprintf(stderr,"%s %g\n",msg,mve);

			   ixb+=-(int)(mve/dx);
			   if ( (ixb<0) ||
			        ((ixb+nxb)>nx) ||
			        (nxb<2) ||
			        (nxb>nx)) {ixb=0;nxb=nx;
				           x2begb=x2beg;
					   x2endb=x2end;}

			   if (czb!=cz) free1(czb);
			   czb = ealloc1(nxb*nyb,
				     sizeof(signed char));
			   for (i=0,czbp=czb; i<nyb; i++) {
			       czp = cz+(iyb+i)*nx+ixb;
			       for (j=0; j<nxb; j++)
				    *czbp++ = *czp++;
			   }

			  /* clear area and force an expose event */
			  XClearArea(dpy,win,0,0,0,0,True);
			  /* note that image is out of date */
			  imageOutOfDate = 1;

			} else if (keysym==XK_Right ) {
			  /* move zoom box to right by half window width*/
			  mve = (x2endb - x2begb)/mvefac ;
			  x2begb = x2begb + mve ;
			  x2endb = x2endb + mve ;
			  msg="move ";
			  /* check for bounds of full window */
			  if (x2endb > x2end){
			    if ( lock ) { x2begb = x2begb - mve ;
			                  x2endb = x2endb - mve ;
					  msg="limit ";
					  mve=0;
			    } else { x2endb = x2end;
			             nxb=(int)((x2endb-x2begb)/dx);
				   }
			  }
			  if (verbose) fprintf(stderr,"%s %g\n",msg,mve);

			  /* for replot require
			   * ixb,iyb   start samples of image
			   * nxb,nyb   number of samples of image */

			   ixb+=(int)(mve/dx);
			   if ( (ixb<0) ||
			        ((ixb+nxb)>nx) ||
			        (nxb<2) ||
			        (nxb>nx)) {ixb=0;nxb=nx;
				           x2begb=x2beg;
					   x2endb=x2end;}


     			   if (czb!=cz) free1(czb);
			   czb = ealloc1(nxb*nyb,
				     sizeof(signed char));
			   for (i=0,czbp=czb; i<nyb; i++) {
			       czp = cz+(iyb+i)*nx+ixb;
			       for (j=0; j<nxb; j++)
				    *czbp++ = *czp++;
			   }


			  /* clear area and force an expose event */
			  XClearArea(dpy,win,0,0,0,0,True);
			  /* note that image is out of date */
			  imageOutOfDate = 1;

			} else if (keysym==XK_Down ) {
			  /* move zoom box down by half window height */
			  mve = (x1endb - x1begb)/mvefac ;
			  x1begb = x1begb + mve ;
			  x1endb = x1endb + mve ;
			  msg="move ";
			  /* check for bounds of full window */
			  if (x1endb > x1end){
			    if ( lock ) { x1begb = x1begb - mve ;
			                  x1endb = x1endb - mve ;
					  msg="limit ";
					  mve=0;
			    } else { x1endb = x1end;
			             nyb=(int)((x1endb-x1begb)/dy);
				   }
			  }
			  if (verbose) fprintf(stderr,"%s %g\n",msg,mve);

			   iyb+=(int)(mve/dy);

			   /* reset to original if out of range */
			   if ( (iyb<0) ||
			        ((iyb+nyb)>ny) ||
			        (nyb<2) ||
			        (nyb>ny)) {iyb=0;nyb=ny;
				           x1begb=x1beg;
					   x1endb=x1end;}


			  /* clear area and force an expose event */
			  XClearArea(dpy,win,0,0,0,0,True);
			  /* note that image is out of date */
			  imageOutOfDate = 1;

     			   if (czb!=cz) free1(czb);
			   czb = ealloc1(nxb*nyb,
				     sizeof(signed char));
			   for (i=0,czbp=czb; i<nyb; i++) {
			       czp = cz+(iyb+i)*nx+ixb;
			       for (j=0; j<nxb; j++)
				    *czbp++ = *czp++;
			   }

			} else if (keysym==XK_Up || keysym==XK_KP_Up ) {
			  /***********
			   * move zoom box up in .... vertical*
			   ***********                          */
			  mve = (x1endb - x1begb)/mvefac ;
			  x1begb = x1begb - mve ;
			  x1endb = x1endb - mve ;
			  msg="move ";
			  /* check for bounds of full window */
			  if (x1begb < x1beg){
			    if ( lock ) { x1begb = x1begb + mve ;
			                  x1endb = x1endb + mve ;
					  msg="limit ";
					  mve=0;
			    } else { x1begb = x1beg ;
			             nyb=(int)((x1endb-x1begb)/dy);
				   }
			  }
			  if (verbose) fprintf(stderr,"%s %g\n",msg,mve);

			  iyb+=-(int)(mve/dy);

			  /* reset to original if out of range */
			   if ( (iyb<0) ||
			        (nyb<2) ||
			        (nyb>ny)) {iyb=0;nyb=ny;
				           x1begb=x1beg;
					   x1endb=x1end;}

     			   if (czb!=cz) free1(czb);
			   czb = ealloc1(nxb*nyb,
				     sizeof(signed char));
			   for (i=0,czbp=czb; i<nyb; i++) {
			       czp = cz+(iyb+i)*nx+ixb;
			       for (j=0; j<nxb; j++)
				    *czbp++ = *czp++;
			   }

			/* clear area and force an expose event */
			XClearArea(dpy,win,0,0,0,0,True);

			/* note that image is out of date */
			imageOutOfDate = 1;

			} else if (keysym==XK_o || keysym==XK_KP_Subtract ) {
			  /***********
			   *zoom out .... vertical*
			   ***********            */
			  mve = (x1endb - x1begb)/mvefac ;
			  x1begb = x1begb - mve ;
			  x1endb = x1endb + mve ;
			  /* check for bounds of full window */
			  if (x1begb < x1beg){
			    if ( lock ) { x1begb = x1begb + mve ;
					  msg="limit ";
					  mve=0;
			    } else { x1begb = x1beg ;}
			  }
			  if (x1endb > x1end){
			    if ( lock ) { x1endb = x1endb - mve ;
					  msg="limit ";
					  mve=0;
			    } else { x1endb = x1end ;}
			  }
 		           nyb=(int)((x1endb-x1begb)/dy);
			   iyb+=-(int)(mve/dy);
			   if ( (iyb<0) || (nyb>ny)) {iyb=0;nyb=ny;}

			  /*   .... and horizontal */
			  mve = (x2endb - x2begb)/mvefac ;
			  x2begb = x2begb - mve ;
			  x2endb = x2endb + mve ;
			  /* check bounds of original image */
			  if (x2begb < x2beg){
			    if ( lock ) { x2begb = x2begb + mve ;
					  msg="limit ";
					  mve=0;
			    } else { x2begb = x2beg ;}
			  }
			  if (x2endb > x2end){
			    if ( lock ) { x2endb = x2endb - mve ;
					  msg="limit ";
					  mve=0;
			    } else { x2endb = x2end ;}
			  }
			   nxb=(int)((x2endb-x2begb)/dx);
			   ixb+=-(int)(mve/dx);
 			   if ( (ixb<0)        ||
			        ((ixb+nxb)>nx) ||
			        (nxb<0)        ||
			        (nxb>nx))  { ixb=0;nxb=nx;
				             x2begb=x2beg;
					     x2endb=x2end;}

			   if (czb!=cz) free1(czb);
			   czb = ealloc1(nxb*nyb,
				     sizeof(signed char));
			   for (i=0,czbp=czb; i<nyb; i++) {
			       czp = cz+(iyb+i)*nx+ixb;
			       for (j=0; j<nxb; j++)
				    *czbp++ = *czp++;
			   }
			  /* clear area and force an expose event */
		   	  XClearArea(dpy,win,0,0,0,0,True);

			  /* note that image is out of date */
			  imageOutOfDate = 1;

			} else if (keysym==XK_i || keysym==XK_KP_Add ) {
			  /***********
			   *zoom in .... vertical*
			   ***********           */
			  mve = (x1endb - x1begb)/(2.*mvefac) ;
			  x1begb = x1begb + mve ;
			  x1endb = x1endb - mve ;
			   iyb+=(int)(mve/dy);

			  /*   .... and horizontal */
			  mve = (x2endb - x2begb)/(2.*mvefac) ;
			  x2begb = x2begb + mve ;
			  x2endb = x2endb - mve ;
			   ixb+=(int)(mve/dx);

			   nxb=(int)((x2endb-x2begb)/dx);
			   nyb=(int)((x1endb-x1begb)/dy);
			   if ( (ixb<0) ||
			        (nxb>nx)||
				(ixb>nx)||
				(nxb<0) ) {ixb=0;nxb=nx;
				             x2begb=x2beg;
					     x2endb=x2end;}
			   if ( (iyb<0) ||
			        (nyb>ny)||
				(iyb>ny)||
				(nyb<0) ) {iyb=0;nyb=ny;
				             x1begb=x1beg;
					     x1endb=x1end;}

			  /* clear area and force an expose event */
			  XClearArea(dpy,win,0,0,0,0,True);

			  /* note that image is out of date */
			  imageOutOfDate = 1;
			if (czb!=cz) free1(czb);
					czb = ealloc1(nxb*nyb,
						sizeof(signed char));
					for (i=0,czbp=czb; i<nyb; i++) {
					    czp = cz+(iyb+i)*nx+ixb;
					    for (j=0; j<nxb; j++)
						    *czbp++ = *czp++;
					}
			} else if (keysym==XK_c || keysym==XK_Page_Down) {

				/* Change clip for image */
 		       		clip += clip/10. ;
				if (verbose) warn("clip=%g\n",clip);
 				/* note that image is out of date */
				 imageOutOfDate = 1;

			} else if (keysym==XK_a || keysym==XK_Page_Up) {

				/* Change clip for image */
			        clip -= clip/10. ;
				if (verbose) warn("clip=%g\n",clip);
				/* note that image is out of date */
				imageOutOfDate = 1;

				if (czb!=cz) free1(czb);
					czb = ealloc1(nxb*nyb,
						sizeof(signed char));
					for (i=0,czbp=czb; i<nyb; i++) {
					    czp = cz+(iyb+i)*nx+ixb;
					    for (j=0; j<nxb; j++)
						    *czbp++ = *czp++;
					}
			/* end of section for moving clipping and zooming GK */

			/* This is the exit from the event loop */
			} else if (keysym==XK_q || keysym==XK_Q) {
				break;

			} else if (keysym==XK_p || keysym==XK_P) {
			/* invoke pswigb with appropriate data */
				char *cmdline, cmdtemp[256];
				float cmdfloat;
				int iargc;
				FILE *plotfp;	/*fp for plot data*/

				cmdline = (char *) emalloc(BUFSIZ);
				strcpy(cmdline,"psimage");
				for(iargc = 1; iargc < argc; iargc++) {
					strcat(cmdline," ");
					strcat(cmdline,argv[iargc]);
					}
				/* override incompatible arguments */
				sprintf(cmdtemp," axescolor=%s",labelcolor);
				strcat(cmdline,cmdtemp);
				cmdfloat = DisplayWidthMM(dpy,scr)/25.4;
				cmdfloat /= DisplayWidth(dpy,scr);
				sprintf(cmdtemp," wbox=%g", cmdfloat*width);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," xbox=%g", 0.5+cmdfloat*xbox);
				strcat(cmdline,cmdtemp);
				cmdfloat = DisplayHeightMM(dpy,scr)/25.4;
				cmdfloat /= DisplayHeight(dpy,scr);
				sprintf(cmdtemp," hbox=%g", cmdfloat*height);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," ybox=%g", 0.5+cmdfloat*ybox);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x1beg=%g", x1begb);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x1end=%g",x1endb);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x2beg=%g", x2begb);
				strcat(cmdline,cmdtemp);
				sprintf(cmdtemp," x2end=%g",x2endb);
				strcat(cmdline,cmdtemp);
				if (STREQ(cmap,"gray")) {
					strcat(cmdline," brgb=0.0,0.0,0.0");
					strcat(cmdline," wrgb=1.0,1.0,1.0");
					}
				else if (STREQ(cmap,"hue")) {
					strcat(cmdline," bhls=0.75,0.5,1.0");
					strcat(cmdline," whls=0.0,0.5,1.0");
					}
				strcat(cmdline," title=\"");
				strcat(cmdline,title);
				strcat(cmdline,"\"");
				strcat(cmdline," label1=\"");
				strcat(cmdline,label1);
				strcat(cmdline,"\"");
				strcat(cmdline," label2=\"");
				strcat(cmdline,label2);
				strcat(cmdline,"\"");
				sprintf(cmdtemp," > %s ", plotfile);
				strcat(cmdline,cmdtemp);
				fprintf(stderr,"%i * %i = %i\n",n1,n2,nz);
				fprintf(stderr,"%s\n",cmdline);

                                /* pipe data to psimage and write to plotfile *GK*/
				plotfp = epopen(cmdline, "w");
				free(cmdline);
				efwrite(z,sizeof(float),nz,plotfp);
				epclose(plotfp);

			} else if (keysym==XK_r) {
				Colormap mycp=xCreateRGBColormap(dpy,win,"rgb_up",verbose);

				XSetWindowColormap(dpy,win,mycp);
				XInstallColormap(dpy,mycp);
                                /* clear area and force an expose event */
                                XClearArea(dpy,win,0,0,0,0,True);
                                /* note that image is out of date */
                                imageOutOfDate = 1;


			} else if (keysym==XK_R) {
                                Colormap mycp=xCreateRGBColormap(dpy,win,"rgb_down",verbose);

                                XSetWindowColormap(dpy,win,mycp);
				XInstallColormap(dpy,mycp);

                                /* clear area and force an expose event */
                                XClearArea(dpy,win,0,0,0,0,True);
                                /* note that image is out of date */
                                imageOutOfDate = 1;

			} else if (keysym==XK_h) {
                                Colormap mycp=xCreateHSVColormap(dpy,win,"hsv_up",verbose);

                                XSetWindowColormap(dpy,win,mycp);
				XInstallColormap(dpy,mycp);

                                /* clear area and force an expose event */
                                XClearArea(dpy,win,0,0,0,0,True);
                                /* note that image is out of date */
                                imageOutOfDate = 1;

			} else if (keysym==XK_H) {

                                Colormap mycp=xCreateHSVColormap(dpy,win,"hsv_down",verbose);

                                XSetWindowColormap(dpy,win,mycp);
				XInstallColormap(dpy,mycp);

                                /* clear area and force an expose event */
                                XClearArea(dpy,win,0,0,0,0,True);
                                /* note that image is out of date */
                                imageOutOfDate = 1;

			} else {
				continue;
			}

		/* else if button down (1 == zoom, 2 == mouse tracking */
		} else if (event.type==ButtonPress) {
			/* if 1st button: zoom */
			if (event.xbutton.button==Button1) {

				/* track pointer and get new box */
				xRubberBox(dpy,win,event,&xb,&yb,&wb,&hb);

				/* if new box has tiny width or height */
				if (wb<4 || hb<4) {

					/* reset box to initial values */
					x1begb = x1beg;
					x1endb = x1end;
					x2begb = x2beg;
					x2endb = x2end;
					nxb = nx;
					nyb = ny;
					ixb = iyb = 0;
					if (czb!=cz) free1(czb);
					czb = cz;

				/* else, if new box has non-zero width */
				/* and height */
				} else {

					/* calculate new box parameters */
					if (style==NORMAL) {
					    zoomBox(x,y,width,height,
						    xb,yb,wb,hb,
						    nxb,ixb,x1begb,x1endb,
						    nyb,iyb,x2endb,x2begb,
						    &nxb,&ixb,&x1begb,&x1endb,
						    &nyb,&iyb,&x2endb,&x2begb);
					} else {
					    zoomBox(x,y,width,height,
						    xb,yb,wb,hb,
						    nxb,ixb,x2begb,x2endb,
						    nyb,iyb,x1begb,x1endb,
						    &nxb,&ixb,&x2begb,&x2endb,
						    &nyb,&iyb,&x1begb,&x1endb);
					}

					/* make new bytes in zoombox */
					if (czb!=cz) free1(czb);
					czb = ealloc1(nxb*nyb,
						sizeof(signed char));
					for (i=0,czbp=czb; i<nyb; i++) {
					    czp = cz+(iyb+i)*nx+ixb;
					    for (j=0; j<nxb; j++)
						    *czbp++ = *czp++;
					}
				}

				/* clear area and force an expose event */
				XClearArea(dpy,win,0,0,0,0,True);

				/* note that image is out of date */
				imageOutOfDate = 1;

			/* else if 2nd button down: display mouse coords */
			} else if (event.xbutton.button==Button2) {

				showloc = 1;
				xMouseLoc(dpy,win,event,style,showloc,
					  x,y,width,height,x1begb,x1endb,
					  x2begb,x2endb,z,f1,d1,n1,f2,d2,
                                          n2,verbose);
			} else if (event.xbutton.button==Button3) {
                                /* ZM: Mouse-Button-3 shows the amplitude constantly */

				showloc = 1;
				xMouseLoc(dpy,win,event,style,showloc,
					  x,y,width,height,x1begb,x1endb,
					  x2begb,x2endb,z,f1,d1,n1,f2,d2,
                                          n2,verbose);

			} else {
				continue;
			}

		/* else if pointer has moved */
		} else if (event.type==MotionNotify) {

			/* if button2 down, show mouse location */
			if (showloc)
				xMouseLoc(dpy,win,event,style,True,
					x,y,width,height,x1begb,x1endb,
                                        x2begb,x2endb,z,f1,d1,n1,f2,d2,
                                        n2,verbose);

		/* else if button2 released, stop tracking */
		} else if (event.type==ButtonRelease &&
			   	event.xbutton.button==Button2) {
			showloc = 0;

		/* exit */
		} else if (event.type==ClientMessage &&
				event.xclient.data.l[0] == wm_delete_window) {
			break;
		}

	} /* end of event loop */

	/* close connection to X server */
	XCloseDisplay(dpy);
	if (curve) {
		free1int(npair);
		free1int(curvewidth);
		for (i=0; i<curve; i++) {
			free1float(x1curve[i]);
			free1float(x2curve[i]);
		}
		free((void**)x1curve);
		free((void**)x2curve);
		free((void**)curvefile);
		free((void**)curvecolor);
	}
    	free1float(z);

	return EXIT_SUCCESS;
}

/* update parameters associated with zoom box */
static void zoomBox (int x, int y, int w, int h,
	int xb, int yb, int wb, int hb,
	int nx, int ix, float x1, float x2,
	int ny, int iy, float y1, float y2,
	int *nxb, int *ixb, float *x1b, float *x2b,
	int *nyb, int *iyb, float *y1b, float *y2b)
{
	/* if width and/or height of box are zero, just copy values */
	if (wb==0 || hb==0) {
		*nxb = nx; *ixb = ix; *x1b = x1; *x2b = x2;
		*nyb = ny; *iyb = iy; *y1b = y1; *y2b = y2;
		return;
	}

	/* clip box */
	if (xb<x) {
		wb -= x-xb;
		xb = x;
	}
	if (yb<y) {
		hb -= y-yb;
		yb = y;
	}
	if (xb+wb>x+w) wb = x-xb+w;
	if (yb+hb>y+h) hb = y-yb+h;

	/* determine number of samples in rubber box (at least 2) */
	*nxb = MAX(nx*wb/w,2);
	*nyb = MAX(ny*hb/h,2);

	/* determine indices of first samples in box */
	*ixb = ix+(xb-x)*(nx-1)/w;
	*ixb = MIN(*ixb,ix+nx-*nxb);
	*iyb = iy+(yb-y)*(ny-1)/h;
	*iyb = MIN(*iyb,iy+ny-*nyb);


	/* determine box limits to nearest samples */
	*x1b = x1+(*ixb-ix)*(x2-x1)/(nx-1);
	*x2b = x1+(*ixb+*nxb-1-ix)*(x2-x1)/(nx-1);
	*y1b = y1+(*iyb-iy)*(y2-y1)/(ny-1);
	*y2b = y1+(*iyb+*nyb-1-iy)*(y2-y1)/(ny-1);
}

/* return pointer to new interpolated array of bytes */
static unsigned char *newInterpBytes (int n1in, int n2in, unsigned char *bin,
	int n1out, int n2out, int useBlockInterp) /* JG */
{
	unsigned char *bout;
	float d1in,d2in,d1out,d2out,f1in,f2in,f1out,f2out;

	f1in = f2in = f1out = f2out = 0.0;
	d1in = d2in = 1.0;
	d1out = d1in*(float)(n1in-1)/(float)(n1out-1);
	d2out = d2in*(float)(n2in-1)/(float)(n2out-1);
	bout = ealloc1(n1out*n2out,sizeof(unsigned char));
	 /* JG .... */
	if (!useBlockInterp)
	  {
		intl2b(n1in,d1in,f1in,n2in,d2in,f2in,bin,
			   n1out,d1out,f1out,n2out,d2out,f2out,bout);
	  }
	else
	  {
		intl2b_block(n1in,d1in,f1in,n2in,d2in,f2in,bin,
		n1out,d1out,f1out,n2out,d2out,f2out,bout);
	  }
	/* .... JG */
	return bout;
}

/*********************** self documentation **********************/
/*****************************************************************************
INTL2B_block - blocky interpolation of a 2-D array of bytes

intl2b_block		blocky interpolation of a 2-D array of bytes

******************************************************************************
Function Prototype:
void intl2b_block(int nxin, float dxin, float fxin,
	int nyin, float dyin, float fyin, unsigned char *zin,
	int nxout, float dxout, float fxout,
	int nyout, float dyout, float fyout, unsigned char *zout);

******************************************************************************
Input:
nxin		number of x samples input (fast dimension of zin)
dxin		x sampling interval input
fxin		first x sample input
nyin		number of y samples input (slow dimension of zin)
dyin		y sampling interval input
fyin		first y sample input
zin		array[nyin][nxin] of input samples (see notes)
nxout		number of x samples output (fast dimension of zout)
dxout		x sampling interval output
fxout		first x sample output
nyout		number of y samples output (slow dimension of zout)
dyout		y sampling interval output
fyout		first y sample output

Output:
zout		array[nyout][nxout] of output samples (see notes)

******************************************************************************
Notes:
The arrays zin and zout must passed as pointers to the first element of
a two-dimensional contiguous array of unsigned char values.

Constant extrapolation of zin is used to compute zout for
output x and y outside the range of input x and y.

******************************************************************************
Author:  James Gunning, CSIRO Petroleum 1999. Hacked from
intl2b() by Dave Hale, Colorado School of Mines, c. 1989-1991
*****************************************************************************/
/**************** end self doc ********************************/

void intl2b_block(int nxin, float dxin, float fxin,
	int nyin, float dyin, float fyin, unsigned char *zin,
	int nxout, float dxout, float fxout,
	int nyout, float dyout, float fyout, unsigned char *zout)
/*****************************************************************************
blocky interpolation of a 2-D array of bytes: gridblock effect
******************************************************************************
Input:
nxin		number of x samples input (fast dimension of zin)
dxin		x sampling interval input
fxin		first x sample input
nyin		number of y samples input (slow dimension of zin)
dyin		y sampling interval input
fyin		first y sample input
zin		    array[nyin][nxin] of input samples (see notes)
nxout		number of x samples output (fast dimension of zout)
dxout		x sampling interval output
fxout		first x sample output
nyout		number of y samples output (slow dimension of zout)
dyout		y sampling interval output
fyout		first y sample output

Output:
zout		array[nyout][nxout] of output samples (see notes)
******************************************************************************
Notes:
The arrays zin and zout must passed as pointers to the first element of
a two-dimensional contiguous array of unsigned char values.

Constant extrapolation of zin is used to compute zout for
output x and y outside the range of input x and y.

Mapping of bytes between arrays is done to preserve appearance of `gridblocks':
no smooth interpolation is performed.

*****************************************************************************/
{
	int ixout,iyout,iin,jin;
	float xoff,yoff;

	xoff=fxout+0.5*dxin-fxin;
	yoff=fyout+0.5*dyin-fyin;
	for (iyout=0;iyout<nyout;iyout++) {
		jin=(int)((iyout*dyout+yoff)/dyin);
		jin=MIN(nyin-1,MAX(jin,0));
		for (ixout=0;ixout<nxout;ixout++) {
			iin=(int)((ixout*dxout+xoff)/dxin);
			iin=MIN(nxin-1,MAX(iin,0));
			zout[nxout*iyout+ixout]=zin[nxin*jin+iin];
		}
	}
}

void xMouseLoc(Display *dpy, Window win, XEvent event, int style, Bool show,
	int x, int y, int width, int height,
	float x1begb, float x1endb, float x2begb, float x2endb,
        float *z, float f1, float d1, int n1,float f2, float d2,
        int n2, int verbose)
{
	static XFontStruct *fs=NULL;
	static XCharStruct overall;
	static GC gc;
	int dummy,xoffset=5,yoffset=5;
	float x1,x2,amp;
	char string[256];

	/* if first time, get font attributes and make gc */
	if (fs==NULL) {
		fs = XLoadQueryFont(dpy,"fixed");
		gc = XCreateGC(dpy,win,0,NULL);

		/* make sure foreground/background are black/white */
		XSetForeground(dpy,gc,BlackPixel(dpy,DefaultScreen(dpy)));
		XSetBackground(dpy,gc,WhitePixel(dpy,DefaultScreen(dpy)));

		XSetFont(dpy,gc,fs->fid);
		overall.width = 1;
		overall.ascent = 1;
		overall.descent = 1;
	}

	/* erase previous string */
	XClearArea(dpy,win,xoffset,yoffset,
		overall.width,overall.ascent+overall.descent,False);

	/* if not showing, then return */
	if (!show) return;

	/* convert mouse location to (x1,x2) coordinates */
	if (style==NORMAL) {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.x-x)/width;
		x2 = x2endb+(x2begb-x2endb)*(event.xmotion.y-y)/height;
	} else {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
		x2 = x2begb+(x2endb-x2begb)*(event.xmotion.x-x)/width;
	}

	/* draw string indicating mouse location */

        /* ZM: computing amplitude at the poked point from dataset */
        amp = getamp(z,f1,d1,n1,f2,d2,n2,x1,x2,verbose);
        sprintf(string,"(%0.6g,%0.6g,%0.6g)",x2,x1,amp); /* ZM */

	XTextExtents(fs,string,(int)strlen(string),&dummy,&dummy,&dummy,&overall);
	XDrawString(dpy,win,gc,xoffset,yoffset+overall.ascent,
		string,(int) strlen(string));
}

void xMousePrint(XEvent event, int style, FILE *mpicksfp,
		 int x, int y, int width, int height,
		 float x1begb, float x1endb, float x2begb, float x2endb)
{
	float x1,x2;

	/* convert mouse location to (x1,x2) coordinates */
	if (style==NORMAL) {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.x-x)/width;
		x2 = x2endb+(x2begb-x2endb)*(event.xmotion.y-y)/height;
	} else {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
		x2 = x2begb+(x2endb-x2begb)*(event.xmotion.x-x)/width;
	}

	/* write string indicating mouse location */
	fprintf(mpicksfp, "%0.6g  %0.6g\n", x1, x2);
}

float getamp(float *zz,float f1,float d1,int n1,
             float f2,float d2,int n2,float x1,float x2,int verbose)

/*****************************************************************************
return the amplitude value at x1,x2
******************************************************************************
Input:
zz		zz(n1*nz) is the data
f1		coordinate in first x
d1		x coordinate increment
n1		number samples in x
f2		coordinate in first y
d2		y coordinate increment
n2		number samples in y
x1              x1 coordinate of the probed point
x2              x2 coordinate of the probed point
******************************************************************************
Author: Zhaobo Meng, ConocoPhillips, Feb. 03,2004
*****************************************************************************/
{
        float x1last,x2last,x1i,x2i,xfrac,zfrac,xfrac0,zfrac0,temp;
        int ix1,ix2,i00;

        if (d1==0.0) err("d1 can not be 0.0");
        if (d2==0.0) err("d2 can not be 0.0");

	x1last = f1 + (n1-1)*d1;
	if (x1<f1 || x1>x1last) return -999.0;
	x2last = f2 + (n2-1)*d2;
	if (x2<f2 || x2>x2last) return -999.0;

        x1i = (x1-f1)/d1;
        ix1 = MIN(x1i,n1-1);
        xfrac = x1i-ix1;
        xfrac0 = 1.0-xfrac;

        x2i = (x2-f2)/d2;
        ix2 = MIN(x2i,n2-1);
        zfrac = x2i-ix2;
        zfrac0 = 1.0-zfrac;

        i00 = ix1 + n1*ix2;
        temp = zfrac *( xfrac*zz[i00+1+n1] + xfrac0*zz[i00+n1])
             + zfrac0*( xfrac*zz[i00+1   ] + xfrac0*zz[i00   ]);

        /*
        if (verbose) warn("x1=%g x2=%g,value=%g,ix1=%d,ix2=%d f1=%g d1=%g n1=%d f2=%g d2=%g n2=%d",
          x1,x2,temp,ix1,ix2,f1,d1,n1,f2,d2,n2);
        */

        return(temp);
}
