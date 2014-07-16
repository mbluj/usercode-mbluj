#!/bin/bash

if [ -f "src/MyRootDict.cpp" ]; then 
    echo "The src/MyRootDict.cpp exists -> remove it ";
    rm src/MyRootDict.cpp;
fi
if [ -f "src/MyRootDict.h" ]; then 
    echo "The src/MyRootDict.h exists -> remove it ";
    rm src/MyRootDict.h;
fi

echo "Recreation of the src/MyRootDict.{h,cpp} dictionary files";

if [[ -z $ROOTSYS ]] ; then
    eval `scram runtime -sh`
fi
cd $CMSSW_BASE/src; 
$ROOTSYS/bin/rootcint -f Test/MiniAnalyzer/src/MyRootDict.cpp  -c -p \
    -I${CMSSW_BASE}/src/  \
    -I${CMSSW_RELEASE_BASE}/src/  \
    -I`scram tool info clhep | grep ^INCLUDE | sed 's/INCLUDE=//'` \
    Test/MiniAnalyzer/interface/*Obj.h \
    Test/MiniAnalyzer/interface/MyRoot_Linkdef.h \
    ;
cd -

ls -l src/MyRootDict.*
echo Done

exit