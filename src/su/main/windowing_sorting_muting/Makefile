# Makefile for ...su/main

include $(CWPROOT)/src/Makefile.config


D = $L/libcwp.a $L/libpar.a $L/libsu.a


LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm $(POSTLFLAGS)

PROGS =			\
	$B/subraid	\
	$B/sucleave 	\
	$B/sucommand	\
	$B/sugetgthr	\
	$B/sugprfb	\
	$B/sukill	\
	$B/sumixgathers	\
	$B/sumute	\
	$B/sumutecsv   	\
	$B/supad	\
	$B/supadmissing \
	$B/suputgthr	\
	$B/susort	\
	$B/susorty	\
	$B/susplit	\
	$B/suwind	\
	$B/suwindpoly	\


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
