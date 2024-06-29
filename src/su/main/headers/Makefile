# Makefile for ...su/main

include $(CWPROOT)/src/Makefile.config


D = $L/libcwp.a $L/libpar.a $L/libsu.a


LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm $(POSTLFLAGS)


PROGS =			\
	$B/bhedtopar	\
	$B/su3dchart	\
	$B/suazimuth	\
	$B/suabshw	\
	$B/suahw	\
	$B/suaddhead	\
	$B/subincsv	\
	$B/subinqcsv	\
	$B/subinxxcsv	\
	$B/subinbigcsv	\
	$B/sucdpbin	\
	$B/suchart	\
	$B/suchw	\
	$B/sucliphead	\
	$B/sucountkey	\
	$B/sudumptrace	\
	$B/suedit	\
	$B/suflex2d     \
	$B/sugethw	\
	$B/sugeom	\
	$B/sugeomcsv	\
	$B/suhtmath	\
	$B/sukeycount	\
	$B/sulcthw	\
	$B/sulhead	\
	$B/sunearcsv    \
	$B/sunomstreamers \
	$B/supaste	\
	$B/suprofcsv    \
	$B/surange	\
	$B/surandhw	\
	$B/susehw	\
	$B/sushw	\
	$B/susmoqcsv	\
	$B/sustrip	\
	$B/sutoolcsv	\
	$B/sutrcount	\
	$B/suutm	\
	$B/suxedit	


INSTALL	:	$(PROGS)
	@-rm -f INSTALL
	@touch $@


$(PROGS):	$(CTARGET) $D 
	-$(CC) $(CFLAGS) $(@F).c $(LFLAGS) -o $@
	@$(MCHMODLINE)
	@echo $(@F) installed in $B

remake	:
	-rm -f $(PROGS) INSTALL
	$(MAKE) 
	
clean::
	rm -f a.out junk* JUNK* core
