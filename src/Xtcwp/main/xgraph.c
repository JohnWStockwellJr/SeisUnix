/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* XGRAPH: $Revision: 1.18 $ ; $Date: 2011/11/30 21:15:49 $	*/

#include "par.h"
#include "Xtcwp/Xtcwp.h"
#include "Xtcwp/Axes.h"
#include <X11/keysym.h>
#include <X11/Shell.h>

/*********************** self documentation **********************/
char *sdoc[] = {
	" 									",
	" XGRAPH - X GRAPHer							",
	" Graphs n[i] pairs of (x,y) coordinates, for i = 1 to nplot.		",
	" 									",
	" xgraph n= [optional parameters] <binaryfile 				",
	" 									",
	" X Functionality:                                                      ",
	" q or Q key    Quit                                                    ",
	" 									",
	" Required Parameters:							",
	" n                      array containing number of points per plot	",
	" 									",
	" Optional Parameters:							",
	" nplot=number of n's    number of plots				",
	" d1=0.0,...             x sampling intervals (0.0 if x coordinates input)",
	" f1=0.0,...             first x values (not used if x coordinates input)",
	" d2=0.0,...             y sampling intervals (0.0 if y coordinates input)",
	" f2=0.0,...             first y values (not used if y coordinates input)",
	" pairs=1,...            =1 for data pairs in format 1.a, =0 for format 1.b",
	" linewidth=1,1,...      line widths in pixels (0 for no lines)		",
	" linecolor=2,3,...      line colors (black=0, white=1, 2,3,4 = RGB, ...)",
	" mark=0,1,2,3,...       indices of marks used to represent plotted points",
	" marksize=0,0,...       size of marks in pixels (0 for no marks)	",
	" x1beg=x1min            value at which axis 1 begins			",
	" x1end=x1max            value at which axis 1 ends			",
	" x2beg=x2min            value at which axis 2 begins			",
	" x2end=x2max            value at which axis 2 ends			",
	" reverse=0              =1 to reverse sequence of plotting curves	",
	" 									",
	" Optional resource parameters (defaults taken from resource database):	",
	" windowtitle=      	 title on window				",
	" width=                 width in pixels of window			",
	" height=                height in pixels of window			",
	" nTic1=                 number of tics per numbered tic on axis 1	",
	" grid1=                 grid lines on axis 1 - none, dot, dash, or solid",
	" label1=                label on axis 1				",
	" nTic2=                 number of tics per numbered tic on axis 2	",
	" grid2=                 grid lines on axis 2 - none, dot, dash, or solid",
	" label2=                label on axis 2				",
	" labelFont=             font name for axes labels			",
	" title=                 title of plot					",
	" titleFont=             font name for title				",
	" titleColor=            color for title				",
	" axesColor=             color for axes					",
	" gridColor=             color for grid lines				",
	" style=                 normal (axis 1 horizontal, axis 2 vertical) or	",
	"                        seismic (axis 1 vertical, axis 2 horizontal)	",
	" 									",
	" Data formats supported:						",
	" 	1.a. x1,y1,x2,y2,...,xn,yn					",
	" 	  b. x1,x2,...,xn,y1,y2,...,yn					",
	" 	2. y1,y2,...,yn (must give non-zero d1[]=)			",
	" 	3. x1,x2,...,xn (must give non-zero d2[]=)			",
	" 	4. nil (must give non-zero d1[]= and non-zero d2[]=)		",
	"   The formats may be repeated and mixed in any order, but if		",
	"   formats 2-4 are used, the d1 and d2 arrays must be specified including",
	"   d1[]=0.0 d2[]=0.0 entries for any internal occurences of format 1.	",
	"   Similarly, the pairs array must contain place-keeping entries for	",
	"   plots of formats 2-4 if they are mixed with both formats 1.a and 1.b.",
	"   Also, if formats 2-4 are used with non-zero f1[] or f2[] entries, then",
	"   the corresponding array(s) must be fully specified including f1[]=0.0",
	"   and/or f2[]=0.0 entries for any internal occurences of format 1 or	",
	"   formats 2-4 where the zero entries are desired.			",
	" mark index:                                                           ",
	" 1. asterisk                                                           ",
	" 2. x-cross                                                            ",
	" 3. open triangle                                                      ",
	" 4. open square                                                        ",
	" 5. open circle                                                        ",
	" 6. solid triangle                                                     ",
	" 7. solid square                                                       ",
	" 8. solid circle                                                       ",
	" 									",
	" Note:	n1 and n2 are acceptable aliases for n and nplot, respectively.	",
	" 									",
	" Example:								",
	" xgraph n=50,100,20 d1=2.5,1,0.33 <datafile				",
	"   plots three curves with equally spaced x values in one plot frame	",
	"   x1-coordinates are x1(i) = f1+i*d1 for i = 1 to n (f1=0 by default)	",
	"   number of x2's and then x2-coordinates for each curve are read	",
	"   sequentially from datafile.						",
	" 									",
	NULL};
/**************** end self doc ********************************/

/* AUTHOR:    Dave Hale and Lydia Deng, Colorado School of Mines, 01/17/91
 * MODIFIED:  Dave Hale, Colorado School of Mines, 04/18/91
 *            fixed conversion of command line resource parameters
 * MODIFIED:  Craig Artley, Colorado School of Mines
 *            Added optional paired data sub-format
 * MODIFIED:  Stewart A. Levin, Mobil
 *            Added keyQuit() action and translation for keyboard "Q"
 * MODIFIED:  James Gunning (CSIRO) 
 * 			  Added reverse= option.
 * MODIFIED:  Bruno Vermeulen
 *            Added XQuit() action exiting app by pressing 'X'
 */

/* client data structures for callbacks */
typedef struct ExposeCDStruct
{
	int nplot;
	int *n;
	float **data;
	int *linewidth;
	int *linecolor;
	int *mark;
	int *marksize;
	int plot_direction;
} ExposeCD;

/* Bill Wingle's typedefs */
typedef struct
{
	int x1, y1, x2, y2;
	int (*func)();
	GC gc;
} GBUFFER;

typedef struct
{
	int start_x, start_y, last_x, last_y;
	GC xorgc;
	GC gc;
	int (*current_func)();
	int foreground, background;
	int id;
	GBUFFER buffer[10];
	int next_pos;
} graphics_data;

graphics_data key_data;
Atom wm_delete_window;

/* callback functions */
void resizeCB(Widget w, XtPointer clientdata, XtcwpAxesCallbackStruct *calldata);
void exposeCB(Widget w, ExposeCD *clientdata, XtcwpAxesCallbackStruct *calldata);
void inputCB(Widget w, XtPointer clientdata, XtcwpAxesCallbackStruct *calldata);

/* Program exit action */
void keyQuit(Widget w, graphics_data *data, XKeyEvent *event);
void XQuit(Widget w, XEvent *ev, String *vector, Cardinal *count);

/* functions defined and used internally */
static void xDrawMark(Display *dpy, Drawable d, GC gc,
					  int x, int y, int index, int size);

XtActionsRec Actions[] = {
	{"XQuit", XQuit}};
int ActionsCount = XtNumber(Actions);

#define NPMAX 4096 /* Arbitrary maximum number of plots allowed	*/
#define NMARKS 9
#define MPLUS 0
#define MASTERISK 1
#define MCROSS 2
#define MTRIANGLE 3
#define MSQUARE 4
#define MCIRCLE 5
#define MFILLEDTRIANGLE 6
#define MFILLEDSQUARE 7
#define MFILLEDCIRCLE 8

int main(int argc, char **argv)
{
	int nplot, n[NPMAX], nn, iplot, npoint, nTic1, nTic2,
		i, j, npar, nd1, nf1, nd2, nf2, npairs,
		linewidth[NPMAX], linecolor[NPMAX],
		mark[NPMAX], pairs[NPMAX], marksize[NPMAX], width, height;
	float d1[NPMAX], f1[NPMAX], d2[NPMAX], f2[NPMAX],
		x1beg, x2beg, x1end, x2end,
		x1min, x1max, x2min, x2max;
	char *label1 = "", *label2 = "", *title = "", *windowtitle = "",
		 *labelFont = "", *titleFont = "",
		 *axesColor = "", *gridColor = "", *titleColor = "",
		 *style = "normal", *grid1 = "none", *grid2 = "none";
	int plot_direction = 1; /* JGHACK */
	float **data;
	XrmValue from, to = {.addr = NULL};
	ExposeCD exposeCD;
	Widget toplevel, axes;
	Arg args[100];
	int nargs;

	/* initialize getpar */
	initargs(argc, argv);
	requestdoc(1);

	/* get parameters needed to interpret datafile */
	for (i = 0; i < NPMAX; i++)
	{
		pairs[i] = 1;
		d1[i] = 0.0;
		f1[i] = 0.0;
		d2[i] = 0.0;
		f2[i] = 0.0;
	}
	npairs = getparint("pairs", pairs);
	nd1 = getparfloat("d1", d1);
	nf1 = getparfloat("f1", f1);
	nd2 = getparfloat("d2", d2);
	nf2 = getparfloat("f2", f2);

	if (!(nn = getparint("n", n)))
		nn = getparint("n1", n);
	if (nn == 0)
		err("Must specify n, the number of points per plot!");
	nplot = nn;
	getparint("n2", &nplot);
	getparint("nplot", &nplot);
	if (nplot > NPMAX)
		err("too many plots");
	for (i = nn; i < nplot; ++i)
		n[i] = n[nn - 1];
	if (npairs > 0)
		for (i = npairs; i < nplot; ++i)
			pairs[i] = pairs[npairs - 1];
	if (nd1 > 0)
		for (i = nd1; i < nplot; ++i)
			d1[i] = d1[nd1 - 1];
	if (nf1 > 0)
		for (i = nf1; i < nplot; ++i)
			f1[i] = f1[nf1 - 1];
	if (nd2 > 0)
		for (i = nd2; i < nplot; ++i)
			d2[i] = d2[nd2 - 1];
	if (nf2 > 0)
		for (i = nf2; i < nplot; ++i)
			f2[i] = f2[nf2 - 1];

	/* read, regularize and compute extreme values of data */
	data = (float **)malloc(nplot * sizeof(float *));
	x2max = x1max = -FLT_MAX;
	x2min = x1min = FLT_MAX;
	for (iplot = 0; iplot < nplot; ++iplot)
	{

		npoint = n[iplot];
		data[iplot] = (float *)malloc(npoint * 2 * sizeof(float));

		/* read data for this plot */
		if (d1[iplot] && d2[iplot])
		{ /* straight line */
			float x, y;
			register int i;
			float *pdata = data[iplot];

			for (i = 0; i < npoint; ++i)
			{
				x = f1[iplot] + i * d1[iplot];
				y = f2[iplot] + i * d2[iplot];
				x1max = MAX(x, x1max);
				x1min = MIN(x, x1min);
				x2max = MAX(y, x2max);
				x2min = MIN(y, x2min);
				*pdata++ = x;
				*pdata++ = y;
			}
		}
		else if (d1[iplot])
		{ /* equally spaced x's */
			float x, y;
			register int i;
			float *pdata = data[iplot];

			for (i = 0; i < npoint; ++i)
			{
				efread(&y, FSIZE, 1, stdin);
				x = f1[iplot] + i * d1[iplot];
				x1max = MAX(x, x1max);
				x1min = MIN(x, x1min);
				x2max = MAX(y, x2max);
				x2min = MIN(y, x2min);
				*pdata++ = x;
				*pdata++ = y;
			}
		}
		else if (d2[iplot])
		{ /* equally spaced y's */
			float x, y;
			register int i;
			float *pdata = data[iplot];

			for (i = 0; i < npoint; ++i)
			{
				efread(&x, FSIZE, 1, stdin);
				y = f2[iplot] + i * d2[iplot];
				x1max = MAX(x, x1max);
				x1min = MIN(x, x1min);
				x2max = MAX(y, x2max);
				x2min = MIN(y, x2min);
				*pdata++ = x;
				*pdata++ = y;
			}
		}
		else
		{ /* pairs */
			float x, y;
			register int i;
			float *pdata = data[iplot];

			if (pairs[iplot])
			{ /* x1,y1,x2,y2,...,xn,yn */
				for (i = 0; i < npoint; ++i)
				{
					efread(&x, FSIZE, 1, stdin);
					efread(&y, FSIZE, 1, stdin);
					x1min = MIN(x, x1min);
					x2max = MAX(y, x2max);
					x2min = MIN(y, x2min);
					*pdata++ = x;
					*pdata++ = y;
				}
			}
			else
			{ /* x1,x2,...,xn,y1,y2,...,yn */
				pdata = data[iplot];
				for (i = 0; i < npoint; ++i)
				{
					if (efread(&x, FSIZE, 1, stdin) != 1)
						err("Error reading input!\n");
					x1max = MAX(x, x1max);
					x1min = MIN(x, x1min);
					*pdata = x;
					pdata += 2;
				}
				pdata = data[iplot] + 1;
				for (i = 0; i < npoint; ++i)
				{
					if (efread(&y, FSIZE, 1, stdin) != 1)
						err("Error reading input!\n");
					x2max = MAX(y, x2max);
					x2min = MIN(y, x2min);
					*pdata = y;
					pdata += 2;
				}
			}
		}
	}

	/* cope with special cases */
	if (x1min == FLT_MAX)
		x1min = x1max = 0.0;
	if (x2min == FLT_MAX)
		x2min = x2max = 0.0;
	if (x1min == x1max)
	{
		x1min -= 1.0;
		x1max += 1.0;
	}
	if (x2min == x2max)
	{
		x2min -= 1.0;
		x2max += 1.0;
	}

	/* initialize toolkit intrinsics and set toplevel parameters */
	toplevel = XtInitialize(argv[0], "XGraph_XAPP_DEF", NULL, 0, &xargc, xargv);
	nargs = 0;
	if (getparstring("windowtitle", &windowtitle))
	{
		XtSetArg(args[nargs], XtNtitle, windowtitle);
		nargs++;
	}
	if (getparstring("windowtitle", &windowtitle))
	{
		XtSetArg(args[nargs], XtNiconName, windowtitle);
		nargs++;
	}
	if (getparint("width", &width))
	{
		XtSetArg(args[nargs], XtNwidth, width);
		nargs++;
	}
	if (getparint("height", &height))
	{
		XtSetArg(args[nargs], XtNheight, height);
		nargs++;
	}
	XtSetArg(args[nargs], XtNinput, TRUE);
	nargs++;
	XtSetValues(toplevel, args, nargs);
	XtAddActions(Actions, ActionsCount);

	/* create axes and set axes parameters */
	axes = XtCreateManagedWidget("axes", xtcwpAxesWidgetClass,
								 toplevel, NULL, 0);
	nargs = 0;
	if (getparstring("axesColor", &axesColor))
	{
		from.addr = (char *)axesColor;
		XtConvertAndStore(axes, XtRString, &from, XtRPixel, &to);
		if (!(to.addr == NULL))
			XtSetArg(args[nargs], XtNaxesColor,
					 *((unsigned long *)to.addr));
		nargs++;
	}
	if (getparstring("gridColor", &gridColor))
	{
		from.addr = (char *)gridColor;
		XtConvertAndStore(axes, XtRString, &from, XtRPixel, &to);
		if (to.addr)
			XtSetArg(args[nargs], XtNgridColor,
					 *((unsigned long *)to.addr));
		nargs++;
	}
	if (getparstring("titleColor", &titleColor))
	{
		from.addr = (char *)titleColor;
		XtConvertAndStore(axes, XtRString, &from, XtRPixel, &to);
		if (to.addr)
			XtSetArg(args[nargs], XtNtitleColor,
					 *((unsigned long *)to.addr));
		nargs++;
	}
	if (getparstring("labelFont", &labelFont))
	{
		from.addr = (char *)labelFont;
		XtConvertAndStore(axes, XtRString, &from, XtRFont, &to);
		if (to.addr)
			XtSetArg(args[nargs], XtNlabelFont, *((Font *)to.addr));
		nargs++;
	}
	if (getparstring("titleFont", &titleFont))
	{
		from.addr = (char *)titleFont;
		XtConvertAndStore(axes, XtRString, &from, XtRFont, &to);
		if (to.addr)
			XtSetArg(args[nargs], XtNtitleFont, *((Font *)to.addr));
		nargs++;
	}
	if (getparstring("grid1", &grid1))
	{
		from.addr = (char *)grid1;
		XtConvertAndStore(axes, XtRString, &from, XtcwpRAxesGrid, &to);
		if (!(to.addr == NULL))
			XtSetArg(args[nargs], XtNgrid1, *((int *)to.addr));
		nargs++;
	}
	if (getparstring("grid2", &grid2))
	{
		from.addr = (char *)grid2;
		XtConvertAndStore(axes, XtRString, &from, XtcwpRAxesGrid, &to);
		if (!(to.addr == NULL))
			XtSetArg(args[nargs], XtNgrid2, *((int *)to.addr));
		nargs++;
	}
	if (getparstring("style", &style))
	{
		from.size = (unsigned int)strlen(style);
		from.addr = (char *)style;
		XtConvertAndStore(axes, XtRString, &from, XtcwpRAxesStyle, &to);
		if (!(to.addr == NULL))
			XtSetArg(args[nargs], XtNstyle, *((int *)to.addr));
		nargs++;
	}
	if (getparint("nTic1", &nTic1))
	{
		XtSetArg(args[nargs], XtNnTic1, nTic1);
		nargs++;
	}
	if (getparint("nTic2", &nTic2))
	{
		XtSetArg(args[nargs], XtNnTic2, nTic2);
		nargs++;
	}
	if (getparstring("label1", &label1))
	{
		XtSetArg(args[nargs], XtNlabel1, label1);
		nargs++;
	}
	if (getparstring("label2", &label2))
	{
		XtSetArg(args[nargs], XtNlabel2, label2);
		nargs++;
	}
	if (getparstring("title", &title))
	{
		XtSetArg(args[nargs], XtNtitle, title);
		nargs++;
	}

	XtSetValues(axes, args, nargs);
	x1beg = x1min;
	getparfloat("x1beg", &x1beg);
	x1end = x1max;
	getparfloat("x1end", &x1end);
	x2beg = x2min;
	getparfloat("x2beg", &x2beg);
	x2end = x2max;
	getparfloat("x2end", &x2end);
	XtcwpSetAxesValues(axes, x1beg, x1end, x2beg, x2end);

	/* get graph parameters */
	for (i = 0; i < nplot; i++)
	{
		mark[i] = i % NMARKS;
		marksize[i] = 0;
		linewidth[i] = 1;
		linecolor[i] = 2 + i % 6;
	}
	if ((npar = getparint("mark", mark)))
		for (j = npar; j < nplot; j++)
			mark[j] = mark[npar - 1];
	if ((npar = getparint("marksize", marksize)))
		for (j = npar; j < nplot; j++)
			marksize[j] = marksize[npar - 1];
	if ((npar = getparint("linewidth", linewidth)))
		for (j = npar; j < nplot; j++)
			linewidth[j] = linewidth[npar - 1];
	if ((npar = getparint("linecolor", linecolor)))
		for (j = npar; j < nplot; j++)
			linecolor[j] = linecolor[npar - 1];

	/* get plotting order */
	if (getparint("reverse", &plot_direction))
		plot_direction = -1;

	checkpars();

	/* add callbacks to axes widget */
	XtAddCallback(axes, XtNresizeCallback, (XtCallbackProc)resizeCB, NULL);
	exposeCD.nplot = nplot;
	exposeCD.n = &n[0];
	exposeCD.data = data;
	exposeCD.linewidth = &linewidth[0];
	exposeCD.linecolor = &linecolor[0];
	exposeCD.mark = &mark[0];
	exposeCD.marksize = &marksize[0];
	exposeCD.plot_direction = plot_direction;
	XtAddCallback(axes, XtNexposeCallback, (XtCallbackProc)exposeCB, &exposeCD);
	XtAddCallback(axes, XtNinputCallback, (XtCallbackProc)inputCB, NULL);
	
	/* set up keyQuit */
	XtAddEventHandler(axes, KeyPress, FALSE,
		(XtEventHandler)keyQuit, &key_data);
		
	/* setup X exit event __ */
	XtOverrideTranslations(
		toplevel, XtParseTranslationTable("<Message>WM_PROTOCOLS: XQuit()\n"));
			
	XtRealizeWidget(toplevel);
	Display *dpy = XtDisplay(toplevel);
	wm_delete_window = XInternAtom(dpy, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(dpy, XtWindow(toplevel), &wm_delete_window, 1);

	XtMainLoop();
	return EXIT_SUCCESS;
}

void resizeCB(Widget w, XtPointer clientdata, XtcwpAxesCallbackStruct *ca)
{
	if (((char *)clientdata) - ((char *)clientdata))
		resizeCB(w, clientdata, ca);
}

#define NCOLOR 8
static char *color[NCOLOR] = {
	"black",
	"white",
	"red",
	"green",
	"blue",
	"cyan",
	"magenta",
	"yellow",
};

void exposeCB(Widget w, ExposeCD *cd, XtcwpAxesCallbackStruct *ca)
{
	int nplot = cd->nplot;
	int *n = cd->n;
	float **data = cd->data;
	int *linewidth = cd->linewidth;
	int *linecolor = cd->linecolor;
	int *mark = cd->mark;
	int *marksize = cd->marksize;
	int nplot_start, nplot_end;
	Position x = ca->x, y = ca->y;
	Dimension width = ca->width, height = ca->height;
	float x1beg = ca->x1beg, x1end = ca->x1end,
		  x2beg = ca->x2beg, x2end = ca->x2end;
	int style = ca->style;
	Display *dpy = XtDisplay(w);
	Window win = XtWindow(w);
	XWindowAttributes wa;
	Colormap cmap;
	XColor scolor, ecolor;
	XRectangle rect;
	unsigned long black = BlackPixelOfScreen(XtScreen(w));
	GC gc;
	int iplot, icolor;
	float xmin, xmax, ymin, ymax, xbase, xscale, ybase, yscale;
	static int firstexpose = 1;

	/* if first expose, check style and, if necessary, swap x,y coords */
	if (firstexpose)
	{
		if (style == XtcwpSEISMIC)
		{
			if (cd->plot_direction == 1)
			{
				nplot_start = 0;
				nplot_end = nplot - 1;
			}
			else
			{
				nplot_start = nplot - 1;
				nplot_end = 0;
			}
			for (iplot = nplot_start; (cd->plot_direction == 1) ? (iplot <= nplot_end) : (iplot >= nplot_end);
				 iplot += cd->plot_direction)
			{
				int ni = n[iplot], i;
				for (i = 0; i < ni; ++i)
				{
					float temp = data[iplot][2 * i];
					data[iplot][2 * i] = data[iplot][2 * i + 1];
					data[iplot][2 * i + 1] = temp;
				}
			}
		}
		firstexpose = 0;
	}

	/* determine current colormap */
	XGetWindowAttributes(dpy, win, &wa);
	cmap = wa.colormap;

	/* create GC */
	gc = XCreateGC(dpy, win, 0, NULL);

	/* set clip to axes box */
	rect.x = x;
	rect.y = y;
	rect.width = width;
	rect.height = height;
	XSetClipRectangles(dpy, gc, 0, 0, &rect, 1, Unsorted);

	/* determine min and max coordinates */
	xmin = x;
	xmax = x + ((int)width);
	ymin = y;
	ymax = y + ((int)height);

	/* determine base and scale factors */
	if (style == XtcwpNORMAL)
	{
		xscale = ((float)width) / (x1end - x1beg);
		xbase = x - x1beg * xscale;
		yscale = ((float)height) / (x2beg - x2end);
		ybase = y - x2end * yscale;
	}
	else
	{
		xscale = ((float)width) / (x2end - x2beg);
		xbase = x - x2beg * xscale;
		yscale = ((float)height) / (x1end - x1beg);
		ybase = y - x1beg * yscale;
	}

	if (cd->plot_direction == 1)
	{
		nplot_start = 0;
		nplot_end = nplot - 1;
	}
	else
	{
		nplot_start = nplot - 1;
		nplot_end = 0;
	}

	/* loop over plots */
	for (iplot = nplot_start; (cd->plot_direction == 1) ? (iplot <= nplot_end) : (iplot >= nplot_end);
		 iplot += cd->plot_direction)
	{
		/* set line width */
		XSetLineAttributes(dpy, gc,
						   linewidth[iplot],
						   LineSolid,
						   CapButt,
						   JoinMiter);

		/* set line color */
		icolor = linecolor[iplot];
		if (icolor < 0)
			icolor = 0;
		else if (icolor >= NCOLOR)
			icolor = NCOLOR - 1;
		if (XAllocNamedColor(dpy, cmap, color[icolor], &scolor, &ecolor))
			XSetForeground(dpy, gc, ecolor.pixel);
		else
			XSetForeground(dpy, gc, black);

		/* draw lines between points */
		if (linewidth[iplot] != 0)
		{
			int xi, yi, in, xilast, yilast, inlast, i, ni = n[iplot];
			float *pdata;

			pdata = data[iplot];
			xi = xbase + xscale * (*pdata++);
			yi = ybase + yscale * (*pdata++);
			in = (xi >= xmin && xi <= xmax && yi >= ymin && yi <= ymax);
			for (i = 1; i < ni; ++i)
			{
				xilast = xi;
				yilast = yi;
				inlast = in;
				xi = xbase + xscale * (*pdata++);
				yi = ybase + yscale * (*pdata++);
				in = (xi >= xmin && xi <= xmax &&
					  yi >= ymin && yi <= ymax);
				if (in || inlast)
					XDrawLine(dpy, win, gc,
							  xilast, yilast, xi, yi);
			}
		}

		/* draw marks at points */
		if (marksize[iplot] != 0)
		{
			int xi, yi, in, i, ni = n[iplot];
			float *pdata;

			pdata = data[iplot];
			for (i = 0; i < ni; ++i)
			{
				xi = xbase + xscale * (*pdata++);
				yi = ybase + yscale * (*pdata++);
				in = (xi >= xmin && xi <= xmax &&
					  yi >= ymin && yi <= ymax);
				if (in)
					xDrawMark(dpy, win, gc, xi, yi,
							  mark[iplot], marksize[iplot]);
			}
		}
	}

	/* free GC */
	XFreeGC(dpy, gc);
}

void inputCB(Widget w, XtPointer clientdata, XtcwpAxesCallbackStruct *ca)
{
	int x = ca->x, y = ca->y, width = ca->width, height = ca->height;
	float x1beg = ca->x1beg, x1end = ca->x1end, x2beg = ca->x2beg, x2end = ca->x2end;
	int style = ca->style;
	XEvent *event = ca->event;
	int xb, yb, wb, hb;
	float x1begn, x1endn, x2begn, x2endn;
	static int firstinput = 1;
	static float x1begs, x1ends, x2begs, x2ends;

	/* if first input, save initial axes limits */
	if (firstinput)
	{
		x1begs = x1beg;
		x1ends = x1end;
		x2begs = x2beg;
		x2ends = x2end;
		firstinput = (int)(((char *)clientdata) - ((char *)clientdata));
	}

	/* track pointer and get rubber box */
	XtcwpRubberBox(XtDisplay(w), XtWindow(w), *event, &xb, &yb, &wb, &hb);

	/* if new box has zero width or height */
	if (wb == 0 || hb == 0)
	{
		/* restore initial limits */
		XtcwpSetAxesValues(w, x1begs, x1ends, x2begs, x2ends);

		/* else if non-zero zoom box */
	}
	else
	{
		/* clip box */
		if (xb < x)
		{
			wb -= x - xb;
			xb = x;
		}
		if (yb < y)
		{
			hb -= y - yb;
			yb = y;
		}
		if (xb + wb > x + width)
			wb = x - xb + width;
		if (yb + hb > y + height)
			hb = y - yb + height;

		/* determine axes limits */
		if (style == XtcwpNORMAL)
		{
			x1begn = x1beg + (xb - x) * (x1end - x1beg) / width;
			x1endn = x1beg + (xb + wb - x) * (x1end - x1beg) / width;
			x2begn = x2end + (yb + hb - y) * (x2beg - x2end) / height;
			x2endn = x2end + (yb - y) * (x2beg - x2end) / height;
		}
		else
		{
			x1endn = x1beg + (yb + hb - y) * (x1end - x1beg) / height;
			x1begn = x1beg + (yb - y) * (x1end - x1beg) / height;
			x2begn = x2beg + (xb - x) * (x2end - x2beg) / width;
			x2endn = x2beg + (xb + wb - x) * (x2end - x2beg) / width;
		}

		/* set axes limits */
		XtcwpSetAxesValues(w, x1begn, x1endn, x2begn, x2endn);
	}

	/* force an expose event */
	XClearArea(XtDisplay(w), XtWindow(w), 0, 0, 0, 0, True);
}

static void xDrawMark(Display *dpy, Drawable d, GC gc, int x, int y, int index, int size)
{
	XPoint points[4];
	switch (index % NMARKS)
	{
	case MPLUS: /* plus */
		XDrawLine(dpy, d, gc, x - (int)(0.5 * size), y, x + (int)(0.5 * size), y);
		XDrawLine(dpy, d, gc, x, y - (int)(0.5 * size), x, y + (int)(0.5 * size));
		break;
	case MASTERISK: /* asterisk */
		XDrawLine(dpy, d, gc, x - (int)(0.5 * size), y, x + (int)(0.5 * size), y);
		XDrawLine(dpy, d, gc, x - (int)(0.25 * size), y - (int)(0.433 * size),
				  x + (int)(0.25 * size), y + (int)(0.433 * size));
		XDrawLine(dpy, d, gc, x + (int)(0.25 * size), y - (int)(0.433 * size),
				  x - (int)(0.25 * size), y + (int)(0.433 * size));
		break;
	case MCROSS: /* X */
		XDrawLine(dpy, d, gc, x - (int)(0.5 * size), y - (int)(0.5 * size),
				  x + (int)(0.5 * size), y + (int)(0.5 * size));
		XDrawLine(dpy, d, gc, x + (int)(0.5 * size), y - (int)(0.5 * size),
				  x - (int)(0.5 * size), y + (int)(0.5 * size));
		break;
	case MTRIANGLE: /* triangle */
		points[0].x = x - (int)(0.5 * size);
		points[0].y = y + (int)(0.25 * size);
		points[1].x = x + (int)(0.5 * size);
		points[1].y = y + (int)(0.25 * size);
		points[2].x = x;
		points[2].y = y - (int)(0.559 * size);
		points[3].x = x - (int)(0.5 * size);
		points[3].y = y + (int)(0.25 * size);
		XDrawLines(dpy, d, gc, points, 4, CoordModeOrigin);
		break;
	case MSQUARE: /* square */
		XDrawRectangle(dpy, d, gc, x - (int)(0.5 * size), y - (int)(0.5 * size),
					   size, size);
		break;
	case MCIRCLE: /* circle */
		XDrawArc(dpy, d, gc, x - (int)(0.5 * size), y - (int)(0.5 * size),
				 size, size, 0, 360 * 64);
		break;
	case MFILLEDTRIANGLE: /* filled triangle */
		points[0].x = x - (int)(0.75 * size);
		points[0].y = y + (int)(0.375 * size);
		points[1].x = x + (int)(0.75 * size);
		points[1].y = y + (int)(0.375 * size);
		points[2].x = x;
		points[2].y = y - (int)(0.838 * size);
		points[3].x = x - (int)(0.75 * size);
		points[3].y = y + (int)(0.375 * size);
		XFillPolygon(dpy, d, gc, points, 4, Convex, CoordModeOrigin);
		break;
	case MFILLEDSQUARE: /* filled square */
		XFillRectangle(dpy, d, gc, x - (int)(0.5 * size), y - (int)(0.5 * size),
					   size, size);
		break;
	case MFILLEDCIRCLE: /* filled circle */
		XFillArc(dpy, d, gc, x - (int)(0.5 * size), y - (int)(0.5 * size),
				 size, size, 0, 360 * 64);
		break;
	}
}

/*****************************************************************************
Functions to exit by key press (q or Q) and pointer press of top right "X"
*****************************************************************************/
void keyQuit(Widget w, graphics_data *data, XKeyEvent *event)
{
	char buffer[2];
	int bufsize = 2;
	KeySym key;
	XComposeStatus compose;

	XLookupString(event, buffer, bufsize, &key, &compose);
	if (key == XK_q || key == XK_Q)
	{
		exit(0);
	}
}
void XQuit(Widget w, XEvent *ev, String *vector, Cardinal *count)
{
	if (ev->type == ClientMessage &&
		((Atom)ev->xclient.data.l[0]) == wm_delete_window)
		exit(0);
}
