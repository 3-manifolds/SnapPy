#! /bin/bash
if [ "$1" == "AMD64" ] ; then
    SOURCE_DIR="Release_VC1900_AMD64/Togl"
    TARGET_DIR="../../../python/togl/win32VC-x86_64-tk8.5"
elif [ "$1" == "X86" ] ; then
    SOURCE_DIR="Release_VC1900_X86/Togl"
    TARGET_DIR="../../../python/togl/win32VC-tk8.5"
else
    echo "usage: update_SnapPy.sh AMD64|X86"
    exit 1
fi
echo cp $SOURCE_DIR/Toglstub21.lib $TARGET_DIR
cp $SOURCE_DIR/Toglstub21.lib $TARGET_DIR
for file in pkgIndex.tcl Togl21.dll Togl21.lib ;
do
    echo cp $SOURCE_DIR/$file $TARGET_DIR/Togl2.1
    cp $SOURCE_DIR/$file $TARGET_DIR/Togl2.1
done
	    
