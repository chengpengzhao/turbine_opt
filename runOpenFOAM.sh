#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory
#------------------------------------------------------------------------------
python3 PARSEC_turbine.py;

caseName="turbine_run";
runpath="./cases/${caseName}";
if [ -d ${runpath} ];then
echo "Delete files......\n\n"
rm -rf ${runpath}
fi

foamCloneCase turbine_template ${runpath};
cp blockMeshDict ${runpath}/system;
cd cases/turbine_run;
blockMesh;
checkMesh;
#renumberMesh -overwrite;
paraFoam 
