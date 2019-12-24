#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory
#------------------------------------------------------------------------------
# generate mesh
echo "Generation mesh for turbineRegion and outRegion..."

if [ -d  outRegionMesh ];then
rm -rf  outRegionMesh
fi
if [ -d turbineRegionMesh ];then
rm -rf turbineRegionMesh
fi

mkdir outRegionMesh;
mkdir turbineRegionMesh;

cp -r turbine_template/* outRegionMesh;
cp -r turbine_template/* turbineRegionMesh;
python3 PARSEC_turbine.py;
python3 out_region.py;
cp blockMeshDict  turbineRegionMesh/system;
cp outRegion_blockMeshDict outRegionMesh/system/blockMeshDict;

cd outRegionMesh;
blockMesh > log.blockMesh;
renumberMesh -overwrite >log.renumberMesh;
checkMesh > log.checkMesh;
mv system/outRegion_createPatchDict system/createPatchDict;
createPatch -overwrite > log.createPatch;
echo "outRegion mesh, accomplished..."
paraFoam &
cd ..;
cd turbineRegionMesh;
blockMesh > log.blockMesh;
renumberMesh -overwrite >log.renumberMesh;
checkMesh > log.checkMesh;
mv system/turbineRegion_createPatchDict system/createPatchDict;
createPatch -overwrite > log.createPatch;
echo "turbineRegionMesh, accomplished..."
#caseName="turbine_run";
#runpath="./${caseName}";
#if [ -d ${runpath} ];then
#echo "Delete files......\n\n"
#rm -rf ${runpath}
#fi
#echo "Clone Case..."
