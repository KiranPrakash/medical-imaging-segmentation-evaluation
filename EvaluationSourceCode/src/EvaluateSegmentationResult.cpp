#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkPoint.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h> 
#include <itkSubtractImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <ANN/ANN.h>




typedef itk::Image<unsigned char, 3>                      SegmentationType;
typedef SegmentationType::Pointer                         SegmentationPointer;
typedef itk::ImageFileReader< SegmentationType >          SegmentationReaderType;
typedef itk::ImageFileWriter< SegmentationType >          SegmentationWriterType;
typedef itk::ImageRegionIterator<SegmentationType>        IteratorType;
typedef itk::BinaryBallStructuringElement<char, 3>        StructuringType;
typedef itk::BinaryErodeImageFilter<SegmentationType, 
  SegmentationType, StructuringType>                      ErodeFilterType;
typedef itk::BinaryThresholdImageFilter< SegmentationType,
  SegmentationType>                                       ThreshFilterType;
typedef itk::SubtractImageFilter<SegmentationType, 
  SegmentationType, SegmentationType>                     SubFilterType;
typedef itk::ResampleImageFilter<SegmentationType,
  SegmentationType>                                       ResampleFilterType;
typedef itk::NearestNeighborInterpolateImageFunction<
  SegmentationType, double>                               InterpolatorType;
typedef itk::GrayscaleFillholeImageFilter< SegmentationType,
  SegmentationType>                                       FillHoleFilterType;
typedef itk::Point<float,3>                               PointType;



void getSurfaceDistance( SegmentationPointer resultImage, SegmentationPointer validationImage, unsigned char startLabel, 
                         double &avgDist, double &rmsDist, double &maxDist )
{
  StructuringType structuringBall;
  structuringBall.SetRadius( 1 );
  structuringBall.CreateStructuringElement();
  PointType pnt;
 
  // compute border pixels and init kd-tree for image 1
  ThreshFilterType::Pointer threshFilter1 = ThreshFilterType::New();
  threshFilter1->SetInput( resultImage );
  threshFilter1->SetLowerThreshold(startLabel);
  threshFilter1->SetUpperThreshold(startLabel+1);
  threshFilter1->SetOutsideValue (0);
  threshFilter1->SetInsideValue (1);

  FillHoleFilterType::Pointer fillFilter1 = FillHoleFilterType::New();
  fillFilter1->SetInput( threshFilter1->GetOutput() );
  
  ResampleFilterType::Pointer resampleFilter1 = ResampleFilterType::New();
  resampleFilter1->SetInput( fillFilter1->GetOutput() );
  InterpolatorType::Pointer interpolator1 = InterpolatorType::New();
  resampleFilter1->SetInterpolator( interpolator1 );
  SegmentationType::SpacingType spacing1 = resultImage->GetSpacing();
  SegmentationType::SizeType size1 = resultImage->GetLargestPossibleRegion().GetSize();
  if (0.5 * spacing1[2] > spacing1[0]) { spacing1[2] *= 0.5; size1[2]*=2; }
  else printf( "WARNING: Strange spacing!\n" );
  resampleFilter1->SetOutputSpacing( spacing1 );
  resampleFilter1->SetOutputOrigin( resultImage->GetOrigin() );
  resampleFilter1->SetSize( size1 );
  
  ErodeFilterType::Pointer erodeFilter1 = ErodeFilterType::New();
  erodeFilter1->SetInput( resampleFilter1->GetOutput() );
  erodeFilter1->SetKernel( structuringBall );
  erodeFilter1->SetErodeValue( 1 );
  SubFilterType::Pointer subFilter1 = SubFilterType::New();
  subFilter1->SetInput2( erodeFilter1->GetOutput() );
  subFilter1->SetInput1( resampleFilter1->GetOutput() );
  subFilter1->UpdateLargestPossibleRegion();
  SegmentationPointer borderImg1 = subFilter1->GetOutput();

  IteratorType it1( borderImg1, borderImg1->GetLargestPossibleRegion() );
  unsigned int numBorderPts1 = 0;
  for ( it1.GoToBegin(); !it1.IsAtEnd(); ++it1 ) {
    if (it1.Get() != 0) numBorderPts1++;
  }
  ANNpointArray borderPts1 = annAllocPts( numBorderPts1, 3 );
  numBorderPts1 = 0;
  for ( it1.GoToBegin(); !it1.IsAtEnd(); ++it1 ) {
    if (it1.Get() != 0) {
      borderImg1->TransformIndexToPhysicalPoint( it1.GetIndex(), pnt );
      for (int d=0; d<3; d++) borderPts1[numBorderPts1][d] = pnt[d];
      numBorderPts1++;
    }
  }
  ANNkd_tree *borderTree1 = new ANNkd_tree( borderPts1, numBorderPts1, 3 );

  // compute border pixels and init kd-tree for image 2
  ThreshFilterType::Pointer threshFilter2 = ThreshFilterType::New();
  threshFilter2->SetInput( validationImage );
  threshFilter2->SetLowerThreshold(startLabel);
  threshFilter2->SetUpperThreshold(startLabel+1);
  threshFilter2->SetOutsideValue (0);
  threshFilter2->SetInsideValue (1);

  FillHoleFilterType::Pointer fillFilter2 = FillHoleFilterType::New();
  fillFilter2->SetInput( threshFilter2->GetOutput() );
  
  ResampleFilterType::Pointer resampleFilter2 = ResampleFilterType::New();
  resampleFilter2->SetInput( fillFilter2->GetOutput() );
  InterpolatorType::Pointer interpolator2 = InterpolatorType::New();
  resampleFilter2->SetInterpolator( interpolator2 );
  SegmentationType::SpacingType spacing2 = validationImage->GetSpacing();
  SegmentationType::SizeType size2 = validationImage->GetLargestPossibleRegion().GetSize();
  if (0.5 * spacing2[2] > spacing2[0]) { spacing2[2] *= 0.5;  size2[2] *= 2; }
  else printf( "WARNING: Strange spacing!\n" );
  resampleFilter2->SetOutputSpacing( spacing2 );
  resampleFilter2->SetOutputOrigin( validationImage->GetOrigin() );
  resampleFilter2->SetSize( size2 );

  ErodeFilterType::Pointer erodeFilter2 = ErodeFilterType::New();
  erodeFilter2->SetInput( resampleFilter2->GetOutput() );
  erodeFilter2->SetKernel( structuringBall );
  erodeFilter2->SetErodeValue( 1 );
  SubFilterType::Pointer subFilter2 = SubFilterType::New();
  subFilter2->SetInput2( erodeFilter2->GetOutput() );
  subFilter2->SetInput1( resampleFilter2->GetOutput() );
  subFilter2->UpdateLargestPossibleRegion();
  SegmentationPointer borderImg2 = subFilter2->GetOutput();

  IteratorType it2( borderImg2, borderImg2->GetLargestPossibleRegion() );
  unsigned int numBorderPts2 = 0;
  for ( it2.GoToBegin(); !it2.IsAtEnd(); ++it2 ) {
    if (it2.Get() != 0) numBorderPts2++;
  }
  ANNpointArray borderPts2 = annAllocPts( numBorderPts2, 3 );
  numBorderPts2 = 0;
  for ( it2.GoToBegin(); !it2.IsAtEnd(); ++it2 ) {
    if (it2.Get() != 0) {
      borderImg2->TransformIndexToPhysicalPoint( it2.GetIndex(), pnt );
      for (int d=0; d<3; d++) borderPts2[numBorderPts2][d] = pnt[d];
      numBorderPts2++;
    }
  }
  ANNkd_tree *borderTree2 = new ANNkd_tree( borderPts2, numBorderPts2, 3 );
  
  // calculate surface distance measures
  avgDist = 0;
  rmsDist = 0;
  maxDist = 0;
  ANNidxArray  nnIdx = new ANNidx[1];
  ANNdistArray dists = new ANNdist[1];
  
  for(unsigned int idx1=0; idx1<numBorderPts1; idx1++) {
    borderTree2->annkSearch( borderPts1[idx1], 1, nnIdx, dists);
    rmsDist += dists[0];
    double d = sqrt( dists[0] );
    avgDist += d;
    if (d>maxDist) maxDist = d;
  }

  for(unsigned int idx2=0; idx2<numBorderPts2; idx2++) {
    borderTree1->annkSearch( borderPts2[idx2], 1, nnIdx, dists);
    rmsDist += dists[0];
    double d = sqrt( dists[0] );
    avgDist += d;
    if (d>maxDist) maxDist = d;
  }

  double numBorderPts = numBorderPts1 + numBorderPts2;
  avgDist /= numBorderPts;
  rmsDist /= numBorderPts;
  rmsDist = sqrt( rmsDist );
  
  // clean up
  annDeallocPts( borderPts1 );
  annDeallocPts( borderPts2 );
  delete borderTree1;
  delete borderTree2;
  delete[] nnIdx;
  delete[] dists;
}



bool evaluateImage( const char *resultFilename, const char *referenceFilename, const char *roiFilename, const char *txtFilename, double &score )
{
  // initialisation
  SegmentationReaderType::Pointer resultReader = SegmentationReaderType::New();
  SegmentationReaderType::Pointer roiReader = SegmentationReaderType::New();
  SegmentationReaderType::Pointer validationReader = SegmentationReaderType::New();
  SegmentationPointer resultImage = resultReader->GetOutput();
  SegmentationPointer roiImage = roiReader->GetOutput();
  SegmentationPointer validationImage = validationReader->GetOutput();

  // read result image:
  try {
    resultReader->SetFileName( resultFilename );
    resultReader->Update();
  }
  catch(...) {
    printf( "Could not load image %s\n", resultFilename );
    return false;
  }
  // read roi image:
  try {
    roiReader->SetFileName( roiFilename );
    roiReader->Update();
  }
  catch(...) {
    printf( "Could not load image %s\n", roiFilename );
    return false;
  }
  // read validation image:
  try {
    validationReader->SetFileName( referenceFilename );
    validationReader->Update();
  }
  catch(...) {
    printf( "Could not load image %s\n", referenceFilename );
    return false;
  }
  // check if images have the same size
  SegmentationType::RegionType resRegion = resultImage->GetLargestPossibleRegion();
  SegmentationType::RegionType roiRegion = roiImage->GetLargestPossibleRegion();
  SegmentationType::RegionType valRegion = validationImage->GetLargestPossibleRegion();
  if (resRegion.GetSize() != valRegion.GetSize() || resRegion.GetSize() != roiRegion.GetSize())
  {
    printf( "Image sizes are different!\n" );
    return false;
  }

  SegmentationType::SpacingType valSpacing = validationImage->GetSpacing();
  double volumeFactor = 0.001*valSpacing[0]*valSpacing[1]*valSpacing[2];
  // adjust spacing and origin of result:
  resultImage->SetSpacing( valSpacing );
  resultImage->SetOrigin( validationImage->GetOrigin() );  
  
  double avgDistFemur, rmsDistFemur, maxDistFemur;
  getSurfaceDistance( resultImage, validationImage, 1, avgDistFemur, rmsDistFemur, maxDistFemur );
  double avgDistTibia, rmsDistTibia, maxDistTiba;
  getSurfaceDistance( resultImage, validationImage, 3, avgDistTibia, rmsDistTibia, maxDistTiba );
  
  char resultBuffer[1024];

  // Tanimoto overlap metric
  unsigned long vol1Fc=0, vol2Fc=0, isFc=0, vol1Tc=0, vol2Tc=0, isTc=0;
  IteratorType resIt( resultImage, resRegion ), valIt( validationImage, valRegion ), roiIt( roiImage, roiRegion );
  for ( resIt.GoToBegin(), valIt.GoToBegin(), roiIt.GoToBegin(); !resIt.IsAtEnd(); ++resIt, ++valIt, ++roiIt ) {
    if (roiIt.Get()==0) continue;

    // femoral ROI
    if (roiIt.Get() & 1) {
      if (resIt.Get()==2) {
        vol1Fc++;
        if (valIt.Get()==2) {
          vol2Fc++;
          isFc++;
        }
      }
      else {
        if (valIt.Get()==2) vol2Fc++;
      }
    }

    // tiboral ROI
    if (roiIt.Get() & 2) {
      if (resIt.Get()==4) {
        vol1Tc++;
        if (valIt.Get()==4) {
          vol2Tc++;
          isTc++;
        }
      }
      else {
        if (valIt.Get()==4) vol2Tc++;
      }
    }

  }

  double tanimotoValFc = 100.0 * (double)(isFc) / ((double)(vol1Fc+vol2Fc-isFc));
  double tanimotoErrorFc = 100.0 - tanimotoValFc;
  double volumeSegFc = (double)vol1Fc * volumeFactor;
  double volumeRefFc = (double)vol2Fc * volumeFactor;
  double volumeDifFc = volumeSegFc - volumeRefFc;
  double volumeDifPercFc = 100.0 * volumeDifFc / volumeRefFc;

  double tanimotoValTc = 100.0 * (double)(isTc) / ((double)(vol1Tc+vol2Tc-isTc));
  double tanimotoErrorTc = 100.0 - tanimotoValTc;
  double volumeSegTc = (double)vol1Tc * volumeFactor;
  double volumeRefTc = (double)vol2Tc * volumeFactor;
  double volumeDifTc = volumeSegTc - volumeRefTc;
  double volumeDifPercTc = 100.0 * volumeDifTc / volumeRefTc;

  // calculate score:
  double scrFemurAvg = 100.0 - 25.0 * avgDistFemur / 0.58;     if (scrFemurAvg < 0) scrFemurAvg=0;
  double scrFemurRMS = 100.0 - 25.0 * rmsDistFemur / 0.88;     if (scrFemurRMS < 0) scrFemurRMS=0;
  double scrTibiaAvg = 100.0 - 25.0 * avgDistTibia / 0.39;     if (scrTibiaAvg < 0) scrTibiaAvg=0;
  double scrTibiaRMS = 100.0 - 25.0 * rmsDistTibia / 0.60;     if (scrTibiaRMS < 0) scrTibiaRMS=0;
  double scrFCartVOE = 100.0 - 25.0 * tanimotoErrorFc / 26.6;  if (scrFCartVOE < 0) scrFCartVOE=0;
  double scrFCartVD  = 100.0 - 25.0 * fabs(volumeDifPercFc) / 5.4;   if (scrFCartVD  < 0) scrFCartVD=0;
  double scrTCartVOE = 100.0 - 25.0 * tanimotoErrorTc / 28.4;  if (scrTCartVOE < 0) scrTCartVOE=0;
  double scrTCartVD  = 100.0 - 25.0 * fabs(volumeDifPercTc) / 6.3;   if (scrTCartVD  < 0) scrTCartVD=0;
  
  double boneScore    = (scrFemurAvg+scrFemurRMS+scrTibiaAvg+scrTibiaRMS) / 4.0;
  double cartScore    = (scrFCartVOE+scrFCartVD+scrTCartVOE+scrTCartVD) / 4.0;
  double overallScore = (boneScore+cartScore) / 2.0;
  score = overallScore;

  // generate output:
  sprintf( resultBuffer, "%s; %s; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.1f; %.1f; %.1f\n", 
	   resultFilename, referenceFilename,
     avgDistFemur, rmsDistFemur, maxDistFemur, avgDistTibia, rmsDistTibia, maxDistTiba,
	   tanimotoErrorFc, volumeDifPercFc, tanimotoErrorTc, volumeDifPercTc, boneScore, cartScore, overallScore );

  printf( resultBuffer );

  // append info to specified text file:
  FILE *file = fopen( txtFilename, "a" );
  fputs( resultBuffer, file );
  fclose( file );
  return true;
}




int main(int argc, char* argv[])
{
  if (argc < 4)
  {
    std::cout << "Usage: " << argv[0] << " segmentation reference roi [-o OutputName]" << std::endl;
    return -1;
  }

/*
  // use this for batch evaluation of multiple consecutively enumerated images
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " imgBase [-o OutputName]" << std::endl;
    return -1;
  }
*/

  // save evaluation data with this name
  std::string evaluationName = "";
  for (int i=4; i<argc; i++) {
    if ( strcmp( argv[i], "-o" )==0 && i<(argc-1)) {
      evaluationName = argv[i+1];
      break;
    }
  }
  if (evaluationName == "") {
    evaluationName = "evaluation.txt";
  }
  
  double score = 0;

  // single evaluation
  bool ok = evaluateImage( argv[1], argv[2], argv[3], evaluationName.c_str(), score );

/*
  // use this for batch evaluation of multiple consecutively enumerated images
  double avgScore = 0;
  for (int imgIdx=61; imgIdx<=100; imgIdx++)
  {
    char resultFilename[512];
    sprintf( resultFilename, "%s-%03i.mhd", argv[1], imgIdx );
    char referenceFilename[512];
    sprintf( referenceFilename, "Insert-path-to-reference-images-here/labels-%03i.mhd", imgIdx );
    char roiFilename[512];
    sprintf( roiFilename, "Insert-path-to-reference-images-here/roi-image-%03i.mhd", imgIdx );
    bool ok = evaluateImage( resultFilename, referenceFilename, roiFilename, evaluationName.c_str(), score );
    avgScore += score;
  }
  avgScore /= 10.0;
  printf ("Average score = %.1f\n", avgScore );
*/  

  return 0;
}
