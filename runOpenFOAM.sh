#!/bin/bash
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

# adjust patches and cellZones
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
mv system/outRegion_topoSetDict system/topoSetDict;
topoSet > log.topoSet;
echo "outRegion mesh, accomplished..."

cd ..;
cd turbineRegionMesh;
blockMesh > log.blockMesh;
renumberMesh -overwrite >log.renumberMesh;
checkMesh > log.checkMesh;
mv system/turbineRegion_createPatchDict system/createPatchDict;
createPatch -overwrite > log.createPatch;
mv system/turbineRegion_topoSetDict system/topoSetDict;
topoSet > log.topoSet;
echo "turbineRegionMesh, accomplished..."
cd ..;

#------------------------------------------------------------------------------
# build final mesh
caseName="turbine_run";
runpath="./${caseName}";
if [ -d ${runpath} ];then
echo "Delete files......"
rm -rf ${runpath}
fi

cp -r turbineRegionMesh/ ${runpath}
cd ${runpath}
rm log.*
echo "Merge meshes..."
mergeMeshes . ../outRegionMesh -overwrite > log.mergeMeshes
renumberMesh -overwrite >log.renumberMesh;
checkMesh > log.checkMesh;
mv system/final_createPatchDict system/createPatchDict;
createPatch -overwrite > log.createPatch;
cp 0.orig/* 0
#decomposePar > log.decomposePar;
#mpirun --allow-run-as-root -np 12 simpleFoam -parallel | tee solve.log;
#reconstructPar -constant > log.reconstructPar;
echo "All done! ..."
