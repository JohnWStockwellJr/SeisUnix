#! /bin/sh
# /*********************** self documentation **********************/
# CWPFIND - look for files with patterns in CWPROOT/src/cwp/lib
#
# Usage: cwpfind  pattern_fragment
#        cwpfind -e exact_pattern
#
# /**************** end self doc ********************************/

# Copyright 1985 by Jack K. Cohen


##################$############################################################
# test for CWPROOT
###############################################################################
if test "${CWPROOT}" = ""
then
	echo "The environment variable \"CWPROOT\" "
	echo "is not set in the user's working shell environment."
	echo "To set this variable in C-shell, use the command: "
	echo "  setenv  CWPROOT  /your/cwp/root/path"
	echo "To set this variable in Bourne or Korn-shell, use the command:"
	echo "  export  CWPROOT=/your/cwp/root/path" ; exit 1

fi

################################################################################
# test for CWPSRC, use value if set, define as $CWPROOT if not set
# (CWPSRC allows one set of source code and documentation for multiple machines)
################################################################################
if test "${CWPSRC}" = ""
then
CWPSRC=$CWPROOT
fi

ROOT=${CWPROOT}
BIN=$ROOT/bin
CWP=$CWPSRC/src/cwp

PATH=/bin:/usr/bin:/usr/ucb:$BIN
CWPLIB=$CWP/lib
PAGE_PROGRAM=more

cmd=`basename $0`
WHITE='[ 	]'

case $# in
1|2)	# OK
;;
*)
    	echo "CWPFIND - search for files with patterns in ../src/cwp/lib"
	echo
    	echo "Usage: $cmd keyword_fragment"
	echo "or:    $cmd -e exact_keyword"
	1>&2; exit 1
esac

option=fragment
word=
for i
do
	case $i in
	-e)
		option=wholeword
	;;
	-*)
    		echo "Illegal option: $i ($cmd [-e] word)"
		1>&2; exit 1
	;;
	*)
		word="$i"
	;;
	esac
done

echo "Scanning $CWPLIB ..."
case $option in
fragment)
	grep -li "$word" $CWPLIB/*.c
;;
*)
	egrep -l "$WHITE$word$WHITE|^$word$WHITE|$WHITE$word$|^$word$"\
		 $CWPLIB/*.c
;;
esac |
sed 's:.*/::'


exit 0
