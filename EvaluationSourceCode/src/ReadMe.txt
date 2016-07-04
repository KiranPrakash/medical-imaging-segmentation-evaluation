This program was written for use in the SKI10 Challenge.

Website: http://www.ski10.org

------------------------------------------------------------------------------

EvaluateSegmentationResult is a program to compare a 3D knee segmentation to 
a reference. Both files have to be in ITK-readable format (e.g. MHD) with 
8-bit pixel values. A value of 0 represents background, 1 is femur bone, 
2 is femoral cartilage, 3 is tibia bone, and 4 is tibial cartilage. Cartilage 
is evaluated only in certain regions of interest, which are defined by 8-bit 
ROI images. Segmentations of the femoral cartilage will be evaluated in 
regions where bit 1 is set (i.e. values 1 and 3). Segmentations of the tibial 
cartilage will be evaluated in regions where bit 2 is set (i.e. values 2 and 
3).

To compile, you need a recent ITK (downloadable from http://www.itk.org), 
CMake (http://www.cmake.org) and the ANN-library 
(http://www.cs.umd.edu/~mount/ANN/).

First compile ANN (or download a precompiled version).
If you have not done so yet, download CMake and use it to compile ITK.
When configuring EvaluateSegmentationResult with CMake, set ANN_PATH to the 
ANN include directory (for Windows precompile it's the main directory) where 
the ANN subdirectory in located.

After successful compilation, you can run the program with segmentation, 
reference, and ROI filename as parameters. It will output the results to the 
screen and append them to evaluation.txt (or another file specified with the 
-o option). The output format is:

SegmentationFileName; ReferenceFileName; 
Femur-AvgSurfaceDist[mm]; Femur-RMSSurfaceDist[mm]; Femur-MaxSurfaceDist[mm];
Tibia-AvgSurfaceDist[mm]; Tibia-RMSSurfaceDist[mm]; Tibia-MaxSurfaceDist[mm];
FemoralCartilage-OverlapError[%]; FemoralCartilage-VolumeError[%]; 
TibialCartilage-OverlapError[%]; TibialCartilage-VolumeError[%]; 
BoneScore; CartilageScore; OverallScore;

The bone score is calculated from average and RMS surface distances for femur 
and tibia. The cartilage score is calculated from overlap and volume errors 
for femoral and tibial cartilage.

With small changes in the main routine, you can make the program evaluate a 
number of segmentations in one run (batch mode). Please read the respective
comments in the source code. 

For any questions or comments contact Tobias Heimann (t.heimann@dkfz.de)

