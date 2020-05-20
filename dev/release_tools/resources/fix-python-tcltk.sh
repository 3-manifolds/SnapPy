#!/bin/bash
FIX_TK=false
FIX_TCL=false
PY=3.8
TCLTK=8.6 # Changing this probably requires a new tkinter.

options=$(getopt --long python,tk,tcl -- "$@")
eval set -- "$options"
shift
shift
while true; do
    case "$1" in
	--python) PY=$2 ; shift ; shift ;;
	--tk) FIX_TK=true ; shift ;;
	--tcl) FIX_TCL=true ; shift ;;
	--) shift ;;
	*) break ;;
    esac
done
if [ ${FIX_TK} == false ] && [ ${FIX_TCL} == false ] ; then
    echo "Specify --tk or --tcl or both." >&2 ; exit 1
fi

TK_FRAMEWORK=/Library/Frameworks/Tk.framework/Versions/${TCLTK}/
TCL_FRAMEWORK=/Library/Frameworks/Tcl.framework/Versions/${TCLTK}/
PYTHON_LIB=/Library/Frameworks/Python.framework/Versions/${PY}/lib/
echo
if [ $FIX_TK == true ] ; then 
    echo Installing Tk ${TCLTK} in Python ${PY}.
    sudo cp  ${TK_FRAMEWORK}Tk ${PYTHON_LIB}libtk${TCLTK}.dylib
    sudo cp ${TK_FRAMEWORK}libtkstub${TCLTK}.a ${PYTHON_LIB}
    sudo cp -R ${TK_FRAMEWORK}Resources/Scripts/ ${PYTHON_LIB}tk${TCLTK}/
fi
if [ $FIX_TCL == true ] ; then 
    echo Installing Tcl ${TCLTK} in Python ${PY}.
    sudo cp  ${TCL_FRAMEWORK}Tcl ${PYTHON_LIB}libtcl${TCLTK}.dylib
    sudo cp ${TCL_FRAMEWORK}libtclstub${TCLTK}.a ${PYTHON_LIB}
    sudo cp -R ${TCL_FRAMEWORK}Resources/Scripts/ ${PYTHON_LIB}tcl${TCLTK}/
fi
