#!/bin/sh
#Routine for building a setup.py file
COMPONENT="paramgui"
AIPSROOT=`echo $CASAPATH | awk '{print $1}'`
ARCHLIB=`echo $CASAPATH | awk '{printf "%s/%s/lib", $1,$2}'`
ARCH=`echo $CASAPATH | awk '{print $2}'`
SITE=`echo $CASAPATH | awk '{print $3}'`
MAKEDEFS=$AIPSROOT/$ARCH/$SITE/makedefs
VARS="CPPSTD PYTHONLIBD PYTHONVER CCMTOOLSLIBD CCMTOOLSLIB CCMTOOLSINCD CORELIB CORELIBD COREINCD QT4LIBD QT4LIB QT4INCD CFITSIOLIBD CFITSIOINCD XTRNLIBS_rpath"
eval `gmake -f $AIPSROOT/$ARCH/makedefs VARS="$VARS" eval_vars`
DEFINES2=`for i in $CPPSTD; do echo $i | grep "\-D"; done`
DEFINES2="$DEFINES2 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE_64_SOURCE"
COMMA=""
SETUPDOTPY="setup.py"

CFITSIO_LIBPATH=""
if [ "$CFITSIOLIBD" != "" ]
then
    CFITSIO_LIBPATH=", '$CFITSIOLIBD'"
fi

# Load makedefs variables into the environment if necessary.
  echo "#$SETUPLIBS#$SETUPPYTHONLIB#SETUPEXTRALINK" | grep '##' > /dev/null 2>&1
  if [ "$?" = 0 ]
  then
     VARS="SETUPLIBS SETUPPYTHONLIB SETUPEXTRALINK"
     eval `gmake -f $AIPSROOT/$ARCH/makedefs VARS="$VARS" eval_vars`
  fi
#
# OK check again if they setup variables are not set then set them to a sensible default
#
  echo "#$SETUPLIBS#$SETUPEXTRALINK#" | grep '##' > /dev/null 2>&1
  if [ "$?" = 0 ]
  then
     if [ "$SETUPEXTRALINK" = "" ]; then
        if [ "$XTRNLIBS_rpath" = "" ]; then
            SETUPEXTRALINK="'-Xlinker', '-rpath', '-Xlinker', '$PYTHONLIBD'"
        else
            SETUPEXTRALINK="'$XTRNLIBS_rpath'"
        fi
     fi
  fi

# cd $AIPSROOT/gcode_$ARCH/Python_Converter
theFiles=`ls *.cc`
echo "from distutils.core import setup, Extension" > $SETUPDOTPY
echo "setup( name=\"$COMPONENT\", version=\"1.0\"," >> $SETUPDOTPY
echo "      ext_modules=[Extension(\"$COMPONENT\"," >> $SETUPDOTPY
echo "               [" >>$SETUPDOTPY
for afile in $theFiles
do
echo "            $COMMA'$afile'" >> $SETUPDOTPY
COMMA=", "
done
echo "               ]," >> $SETUPDOTPY
if [ "$SETUPLIBS" != "" ]; then
   if [ "$QT4LIBD" != "" ]; then
      echo "               library_dirs=['$ARCHLIB', $SETUPLIBS, '$PYTHONLIBD' $CFITSIO_LIBPATH, '$CORELIBD', '$QT4LIBD']," >> $SETUPDOTPY
   else
      echo "               library_dirs=['$ARCHLIB', $SETUPLIBS, '$PYTHONLIBD' $CFITSIO_LIBPATH, '$CORELIBD']," >> $SETUPDOTPY
   fi
else
   if [ "$QT4LIBD" != "" ]; then
      echo "               library_dirs=['$ARCHLIB', '$PYTHONLIBD' $CFITSIO_LIBPATH, '$QT4LIBD', '$CORELIBD']," >> $SETUPDOTPY
   else
      echo "               library_dirs=['$ARCHLIB', '$PYTHONLIBD' $CFITSIO_LIBPATH, '$CORELIBD']," >> $SETUPDOTPY
   fi
fi

core_libraries=`echo "$CORELIB" | perl -pe "s/-l//g; s/(\S+)/'\\\$1'/g; \\\$_ = join(', ', split(/\s+/))"`
echo "               libraries=[ 'nrao', 'xmlcasa', 'display', 'flagging', 'calibration', 'msvis'," >> $SETUPDOTPY
echo "                           'synthesis', 'graphics', 'casaqt', 'qwt'," >> $SETUPDOTPY
echo "                           $core_libraries," >> $SETUPDOTPY
#
#If it's a flavour of darwin the skip lapack and blas
tjunk=`echo "$ARCH" | grep "^darwin"`
if [ "$tjunk" ]; then
   echo "                           'cfitsio', " >> $SETUPDOTPY
else
   echo "                           'lapack', 'blas', 'cfitsio', " >> $SETUPDOTPY
   if [ "$QT4LIB" != "" ]; then
   qtlibs=''
   for i in `echo $QT4LIB | sed 's/-l//g'`; do
       qtlibs="$qtlibs '$i',";
   done
   echo "                          $qtlibs" >> $SETUPDOTPY
   fi
fi
if [ "$SETUPPYTHONLIB" != "" ]; then
echo "$SETUPPYTHONLIB," >> $SETUPDOTPY
fi
ccmtools=''
for i in `echo $CCMTOOLSLIB | sed 's/-l//g'`; do
    ccmtools="$ccmtools '$i',";
done
echo "                          $ccmtools" >> $SETUPDOTPY
echo "                           'c', 'm' ]," >> $SETUPDOTPY
echo "               extra_compile_args = [" >> $SETUPDOTPY
for i in $DEFINES2
do
echo                                       \'$i\', >> $SETUPDOTPY
done
echo "                                     ]," >> $SETUPDOTPY
if [ "$SETUPEXTRALINK" != "" ]; then
echo "               extra_link_args = [$SETUPEXTRALINK])], " >> $SETUPDOTPY
fi

echo "       include_dirs = ['$AIPSROOT/code/include', '$AIPSROOT/code/casa', '$COREINCD', '..', '$CCMTOOLSINCD', '$CFITSIOINCD', '../impl'," >> $SETUPDOTPY
if [ "$QT4INCD" != "" ]; then
   for i in $QT4INCD
   do
      echo "                     '$i'," >> $SETUPDOTPY
   done
fi
echo "                       '.' ])" >> $SETUPDOTPY
