# Makefile for ...su/main

include $(CWPROOT)/src/Makefile.config


D = $L/libcwp.a $L/libpar.a $L/libsu.a


LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm $(POSTLFLAGS)


PROGS =			\
	$B/suaddstatics	\
	$B/sufarldcsv	\
	$B/sugetqcsv	\
	$B/surandstat	\
	$B/surescsv 	\
	$B/suresstat	\
	$B/sushiftcsv   \
	$B/sustatic	\
	$B/sustaticB	\
	$B/sustaticrrs


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
