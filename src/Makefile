#
# Makefile for the CWP/SU free software distribution
#
# Please read Makefile.config and make the necessary changes there
# then type:      make install    to install the standard CWP/SU codes
#                 make xtinstall  to install X-toolkit codes
#                 make finstall   to install Fortran codes
#                 make mglinstall to install Mesa/Open GL codes
#                 make xminstall  to install X-Motif codes (optional)
#		  make sfinstall  to install SFIO materials and SEGDREAD
#
# or if remaking:
# type:           make remake     to remake the standard CWP/SU codes
#                 make xtremake   to remake X-toolkit codes
#                 make fremake    to remake Fortran codes
#                 make mglremake  to remake Mesa/Open GL codes
#                 make xmremake   to remake X-Motif codes
#		  make sfremake	  to install SFIO materials and SEGDREAD
#
# Note: SEGDREAD is the program for reading SEG-D format tapes 
#
# Additional items are included in: cwputils 
# to compile:
# 		make utils		
#
# Do not try to install all of the codes at one time via:
# % make install ; make xtinstall ; make xminstall
# If you get an error message about not finding /src/Makefile.config
# the explicitly define your CWPROOT path on the next line
#CWPROOT =

include $(CWPROOT)/src/Makefile.config


donothing: # To protect against an idle "make" to "see what happens"
	@echo ""
	@echo "This is a dangerous makefile--so the default is do_nothing"
	@echo "Please read and edit Makefile.config appropriately"
	@echo "then type: make install   (to install the basic set of codes)"
	@echo "           make xtinstall (to install the X-toolkit applications)"
	@echo "           make finstall  (to install the Fortran codes)"
	@echo "           make mglinstall (to install the Mesa/ Open GL items)"
	@echo "           make utils     (to install libcwputils) (nonessential) "
	@echo "           make xminstall (to install the Motif application)"
	@echo "           make sfinstall (to install the SFIO version of SEGDREAD)"
	@echo ""
	@echo "or if remaking:"
	@echo "type:      make remake   (to remake the basic set of codes)"
	@echo "           make xtremake (to remake the X-toolkit applications)"
	@echo "           make fremake  (to install the Fortran codes)"
	@echo "           make mglremake (to install the Mesa/ Open GL items)"
	@echo "           make uremake  (to remake libcwputils.a)(nonessential)"
	@echo "           make xmremake (to remake the Motif-based applications)"
	@echo "           make sfremake (to remake the SFIO version of SEGDREAD)"
	@echo " "
	@echo "See the README_ files in ./Portability for more information."

install: checkroot LICENSE_44R23_ACCEPTED MAILHOME_44R23 makedirs cwpstuff plot sustuff tristuff tetrastuff compstuff reflstuff

xtinstall: xtcwp_
xminstall: xmcwp_
mglinstall: mglstuff

# automatic mail message (ask once)
LICENSE_44R23_ACCEPTED	:
	@./license.sh
	@touch $@

# automatic mail message (ask once)
MAILHOME_44R23	:
	@./mailhome.sh
	@touch $@

# check to see if the CWPROOT path is set
checkroot	:
	@./chkroot.sh
	@echo $(CWPROOT)

# Make the bin/include/lib directories if not already there
makedirs: 
	@echo "Making necessary directories"
	@./mkdirectories.sh

cwpstuff:
	cd ./cwp; $(MAKE)
	cd ./par; $(MAKE)

plot	:
	cd ./psplot; $(MAKE)

sustuff	:
	cd ./su ; $(MAKE)

xtcwp_	:
	cd ./Xtcwp; $(MAKE)
	cd ./xplot; $(MAKE)
	cd ./xtri; $(MAKE)

xmcwp_	:
	cd ./Xmcwp; $(MAKE)

complex	:
	cd ./Complex; $(MAKE)

utils	:
	cd ./cwputils; $(MAKE)

tristuff:
	cd ./tri; $(MAKE)
	cd ./Trielas; $(MAKE)

tetrastuff:
	cd ./tetra; $(MAKE)

compstuff:
	cd ./comp; $(MAKE)


reflstuff:
	cd ./Refl ; $(MAKE)

mglstuff:
	cd ./Mesa; $(MAKE)

finstall:
	cd ./Fortran; $(MAKE)

sfinstall:
	cd ./Sfio; $(MAKE)

remake	:
	@./chkroot.sh
	cd ./cwp; $(MAKE) remake
	cd ./par; $(MAKE) remake
	cd ./psplot; $(MAKE) remake
	cd ./su/include; $(MAKE) remake
	cd ./su; $(MAKE) remake
	cd ./tri; $(MAKE) remake
	cd ./tetra; $(MAKE) remake
	cd ./comp; $(MAKE) remake
	cd ./Refl; $(MAKE) remake

xtremake	:
	cd ./Xtcwp; $(MAKE) remake
	cd ./xplot; $(MAKE) remake
	cd ./xtri; $(MAKE) remake

xmremake	:
	cd ./Xmcwp; $(MAKE) remake

mglremake:
	cd ./Mesa; $(MAKE) remake

uremake:
	cd ./cwputils; $(MAKE) remake

tremake:
	cd ./tri; $(MAKE) remake

compmake:
	cd ./comp; $(MAKE) remake

fremake:
	cd ./Fortran; $(MAKE) remake

sfremake:
	cd ./Sfio; $(MAKE) remake

