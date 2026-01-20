# Build SnapPy.app in ./dist when provided all the needed universal2
# wheels are available in ./wheelhouse or on PyPI.

set -e  # Fail if any step fails.

# Delete current app and using bundle_app to create the framework
# scaffolding.

rm -rf SnapPy.app dist
python3.14 -m bundle_app.build
mkdir dist
mv SnapPy.app dist/

# Now install all the needed wheels, and we're done!

TARGET="dist/SnapPy.app/Contents/Frameworks/Python.framework/Versions/\
Current/lib/python3.14/site-packages/"

PYTHON="dist/SnapPy.app/Contents/MacOS/Python"

UNIVERSAL_BINARY=" --platform=macosx_10_15_universal2 --only-binary :all: "

PIP_INSTALL="$PYTHON -m pip install --target $TARGET --upgrade \
 --find-links ./wheelhouse $UNIVERSAL_BINARY"

echo $PIP_INSTALL
$PIP_INSTALL pip
$PIP_INSTALL low_index
$PIP_INSTALL cypari
$PIP_INSTALL knot_floer_homology
$PIP_INSTALL plink==2.4.7
$PIP_INSTALL spherogram==2.4
$PIP_INSTALL snappy==3.3

