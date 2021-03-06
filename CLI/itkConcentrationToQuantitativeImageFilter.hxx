#ifndef _itkConcentrationToQuantitativeImageFilter_hxx
#define _itkConcentrationToQuantitativeImageFilter_hxx
#endif

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_diag_matrix.h" 
//"stores a diagonal matrix as a single vector"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/algo/vnl_qr.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_inverse.h"
#include "vnl/algo/vnl_matrix_inverse.h"

// work around compile error on Windows
#define M_PI 3.1415926535897932384626433832795

#include "itkConcentrationToQuantitativeImageFilter.h"

#include <iostream>

#include <math.h>

#include "nlopt.hpp"


//last ditch try 
#include <itkMatrix.h>
#include <itkVector.h>

namespace itk
{







template <class TInputImage, class TMaskImage, class TOutputImage>
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>::ConcentrationToQuantitativeImageFilter()
{

  //this is the class constructor isn't it? as evidenced by the fact that the above is outputted first 
  //ok try this 

  this->DynamicMultiThreadingOff();

  m_TE = 0.0f;
  m_FA = 0.0f;
  m_RGD_relaxivity = 4.9E-3f;
  m_S0GradThresh = 15.0f;
  m_fTol = 1e-4f;
  m_gTol = 1e-4f;
  m_xTol = 1e-5f;
  m_epsilon = 1e-9f;
  m_maxIter = 200;
  m_aifAUC = 0.0f;
  m_AIFBATIndex = 0;
  m_UsePopulationAIF = false;
  m_UsePrescribedAIF = false;
  m_MaskByRSquared = true;
  m_ModelType = itk::LMCostFunction::TOFTS_2_PARAMETER;
  m_constantBAT = 0;
  m_BATCalculationMode = "PeakGradient";
  this->Superclass::SetNumberOfRequiredInputs(1);
  this->Superclass::SetNthOutput(1, static_cast<TOutputImage*>(this->MakeOutput(1).GetPointer()));  // K2
  this->Superclass::SetNthOutput(2, static_cast<TOutputImage*>(this->MakeOutput(2).GetPointer()));  // K1
  this->Superclass::SetNthOutput(3, static_cast<TOutputImage*>(this->MakeOutput(3).GetPointer()));  // Max slope
  this->Superclass::SetNthOutput(4, static_cast<TOutputImage*>(this->MakeOutput(4).GetPointer()));  // AUC
  this->Superclass::SetNthOutput(5, static_cast<TOutputImage*>(this->MakeOutput(5).GetPointer()));  // R^2
  this->Superclass::SetNthOutput(6, static_cast<TOutputImage*>(this->MakeOutput(6).GetPointer()));  // BAT
  this->Superclass::SetNthOutput(7, static_cast<TOutputImage*>(this->MakeOutput(7).GetPointer()));  // MTT
  this->Superclass::SetNthOutput(8, static_cast<TOutputImage*>(this->MakeOutput(8).GetPointer()));  // CBF
  this->Superclass::SetNthOutput(9, static_cast<VectorVolumeType*>(this->MakeOutput(9).GetPointer())); // fitted

  //why is this 1 above the GetOutput number? 
  //GetOutput is 0-indexed but this isn't? 
}


bool check_tissue(const float* concentration, int timeSize) {

  int x = 0; 

  for (x = 0; x < timeSize; x++) {
    if (concentration[x] < 0 || concentration[x] > 100) {

      return false; 

    }
  }


  return true; 
} 


// start with time to peak --> should be really red on the left hemisphere 
// fit to some sort of curve 

//apply 3 x 3 Guassian filter (voxel, 3 x 3 box around that, look at the distribution/mean, make them all the same )


//"Integral of gamma variate fit"
//fitting of a gamma-variate function required to correct for 
//tracer recirculation 
//remove any elevated post-bolus baseline 

//let's try again but just 
//correct for baseline, not fit it to a gamma-variate function
float integrate(float* yValues, float* xValues, int size) {
// do we have to integrate with constant time points? 
// is that why this isn't working? 
// find out! 

//"baseline corrected"
//"gamma variate function fit "
//"correct curve by gamma variate fit"


  float area= 0.0f;

  //std::cout << "timing!" << xValues[0] << " " << xValues[1]; 

  int x; 
  
  for (x = 0; x < (size - 1); x++) {
   
    area += ((xValues[x + 1] - xValues[x]) * (yValues[x] + yValues[x + 1]))/2;

  }

  if (!check_tissue(yValues, size)) {
    return 0.0; 
  } else {
    return area; 
  }

}


float time_to_peak(float* yValues, std::vector<float>& xValues)
{

  float baseline; 

  baseline = yValues[0]; 


  float peak, time_of_peak;

  int x; 

  peak = 0.0; 

  if (!check_tissue(yValues, xValues.size())) {
    return 0.0; 
  }
  
  for (x = 0; x < (xValues.size() - 1); x++) {
    
    if (yValues[x] > peak) {
        peak = yValues[x]; 
        time_of_peak = xValues[x]; 
    }
    
  }
  
  return time_of_peak;

}



bool my_solver(int signalSize, const float* timeAxis, 
               const float* PixelConcentrationCurve, 
               const float* BloodConcentrationCurve, 
               float& alpha, float& beta, 
               float fTol, float gTol, float xTol,
               float epsilon, int maxIter,
               itk::LevenbergMarquardtOptimizer* optimizer,
               LMCostFunction2* costFunction,
               int modelType
               )
{
  //std::cout << "in pk solver" << std::endl;
  // probe.Start("pk_solver");
  // Note the unit: timeAxis should be in minutes!! This could be related to the following parameters!!
  // fTol      =  1e-4;  // Function value tolerance
  // gTol      =  1e-4;  // Gradient magnitude tolerance 
  // xTol      =  1e-5;  // Search space tolerance
  // epsilon   =  1e-9;    // Step
  // maxIter   =   200;  // Maximum number of iterations
  //std::cerr << "In pkSolver!" << std::endl;

  /* 
  m_BATCalculationMode = BATCalculationMode;
  m_ConstantBAT = constantBAT;
  */

  // Levenberg Marquardt optimizer  
        
  //////////////
  LMCostFunction2::ParametersType initialValue;
  if(modelType == itk::LMCostFunction2::TOFTS_2_PARAMETER)
    {
    initialValue = LMCostFunction2::ParametersType(2); ///...
    }
  else
    {
    initialValue = LMCostFunction2::ParametersType(3);
    initialValue[2] = 0.1;     //f_pv //...
    }

    
  initialValue[0] = 0.1;     //alpha //...
  initialValue[1] = 0.5;     //beta //...
        
  costFunction->SetNumberOfValues (signalSize);
  

  costFunction->SetCb (BloodConcentrationCurve, signalSize); //BloodConcentrationCurve
  costFunction->SetCv (PixelConcentrationCurve, signalSize); //Signal Y
  costFunction->SetTime (timeAxis, signalSize); //Signal X
  costFunction->SetIntCb (BloodConcentrationCurve, timeAxis, signalSize); // integral concentration == leakage
  costFunction->GetValue (initialValue); //...
  costFunction->SetModelType(modelType);

  optimizer->UseCostFunctionGradientOff();   

  try {
     optimizer->SetCostFunction( costFunction ); 
  }
  catch ( itk::ExceptionObject & e ) {
  std::cout << "Exception thrown ! " << std::endl;
  std::cout << "An error ocurred during Optimization" << std::endl;
  std::cout << e << std::endl;
  return false;
  }   
        
  itk::LevenbergMarquardtOptimizer::InternalOptimizerType * vnlOptimizer = optimizer->GetOptimizer();//...

  vnlOptimizer->set_f_tolerance( fTol ); //...
  vnlOptimizer->set_g_tolerance( gTol ); //...
  vnlOptimizer->set_x_tolerance( xTol ); //...
  vnlOptimizer->set_epsilon_function( epsilon ); //...
  vnlOptimizer->set_max_function_evals( maxIter ); //...
        
  // We start not so far from the solution 
        
  optimizer->SetInitialPosition( initialValue ); //...       
  
  try {
  //  probe.Start("optimizer");
  optimizer->StartOptimization();
  //   probe.Stop("optimizer");
  }
  catch( itk::ExceptionObject & e ) {
  std::cerr << "Exception thrown ! " << std::endl;
  std::cerr << "An error ocurred during Optimization" << std::endl;
  std::cerr << "Location    = " << e.GetLocation()    << std::endl;
  std::cerr << "Description = " << e.GetDescription() << std::endl;
  return false;
  }
  //vnlOptimizer->diagnose_outcome();
  //std::cerr << "after optimizer!" << std::endl;
  itk::LevenbergMarquardtOptimizer::ParametersType finalPosition;
  finalPosition = optimizer->GetCurrentPosition();
  //std::cerr << finalPosition[0] << ", " << finalPosition[1] << ", " << finalPosition[2] << std::endl;

        
  //Solution: remove the scale of 100  
  alpha = finalPosition[0];
  beta = finalPosition[1];

  // "Project" back onto the feasible set.  Should really be done as a
  // constraint in the optimization.
  if(beta<0) beta = 0;
  //if(K1>1) K1 = 1;
  if(alpha<0) alpha = 0;
  //if(K2>5) K2 = 5;
		
  //if((Fpv>1)||(Fpv<0)) Fpv = 0;
  //  probe.Stop("pk_solver");
  return true;
}




template< class TInputImage, class TMaskImage, class TOutputImage >
typename ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>::DataObjectPointer
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::MakeOutput(DataObjectPointerArraySizeType idx)
{
  if(idx<8)
  {
    return TOutputImage::New().GetPointer();
  }
  else if (idx==8)
  {
    return VectorVolumeType::New().GetPointer();
  } 
  else 
  {
    //return VectorVolumeType::New().GetPointer(); 
    //return; 
  }
}

// Set a prescribed AIF.  This is not currrently in the input vector,
// though it could be if we used a Decorator.
template< class TInputImage, class TMaskImage, class TOutputImage >
void
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::SetPrescribedAIF(const std::vector<float>& timing, const std::vector<float>& aif)
{
  if (aif.size() < 2)
    {
    itkExceptionMacro(<< "Prescribed AIF must contain at least two time points");
    }
  if (aif.size() != timing.size())
    {
    itkExceptionMacro("Timing vector and concentration vector for AIF must be the same size.");
    }

  m_PrescribedAIF = aif;
  m_PrescribedAIFTiming = timing;
}

// Set 3D AIF mask as second input
template< class TInputImage, class TMaskImage, class TOutputImage >
void
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::SetAIFMask(const TMaskImage* volume)
{
  this->SetNthInput(1, const_cast<TMaskImage*>(volume) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
const TMaskImage*
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::GetAIFMask() const
{
  return dynamic_cast< const TMaskImage * >( this->ProcessObject::GetInput(1) );
}

// Set 3D ROI mask as third input
template< class TInputImage, class TMaskImage, class TOutputImage >
void
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::SetROIMask(const TMaskImage* volume)
{
  this->SetNthInput(2, const_cast<TMaskImage*>(volume) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
const TMaskImage*
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::GetROIMask() const
{
  return dynamic_cast< const TMaskImage * >( this->ProcessObject::GetInput(2) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetK2Output()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(0) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetK1Output()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(1) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetMaxSlopeOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(2) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetAUCOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(3) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetRSquaredOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(4) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetBATOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(5) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetMTTOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(6) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetCBFOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(7) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TInputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetFittedDataOutput()
{
  return dynamic_cast< TInputImage * >( this->ProcessObject::GetOutput(8) );
}

template <class TInputImage, class TMaskImage, class TOutputImage>
void
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::BeforeThreadedGenerateData()
{
  const VectorVolumeType* inputVectorVolume = this->GetInput();
  const MaskVolumeType* maskVolume = this->GetAIFMask();

  //std::cout << "Model type: " << m_ModelType << std::endl;

  const int timeSize = (int)inputVectorVolume->GetNumberOfComponentsPerPixel();

  int   aif_FirstPeakIndex = 0;
  float aif_MaxSlope = 0.0f;

  // calculate AIF
  if (m_UsePrescribedAIF)
    {
    // resample the prescribed AIF vector to be at the specifed
    // m_Timing points and then assign to m_AIF
    m_AIF = std::vector<float>(timeSize);

    std::vector<float>::iterator ait = m_AIF.begin();
    std::vector<float>::iterator tit = m_Timing.begin();

    std::vector<float>::iterator pait = m_PrescribedAIF.begin();
    std::vector<float>::iterator ptit = m_PrescribedAIFTiming.begin();

    std::vector<float>::iterator paitnext = pait;
    paitnext++;
    std::vector<float>::iterator ptitnext = ptit;
    ptitnext++;

    for (; tit != m_Timing.end(); ++tit, ++ait)
      {
      // Three cases
      // (1) extrapolate the aif on the low end of the range of prescribed timings
      // (2) interpolate the aif
      // (3) extrapolate the aif on the high end of the range of prescribed timings
      //
      // Case (1) is handled implictly by the initialization and conditionals.
      if (*ptit <= *tit)
        {
        // Case (2) from above)
        // find the prescribed times that straddle the current time to interpolate
        while (*ptitnext < *tit && ptitnext != m_PrescribedAIFTiming.end())
          {
          ++ptit;
          ++ptitnext;
          ++pait;
          ++paitnext;
          }
        }
      if (ptitnext == m_PrescribedAIFTiming.end())
        {
        // we'll need to extrapolate (Case (3) from above)
        ptitnext = ptit;
        --ptit;
        paitnext = pait;
        --pait;
        }

      // interpolate aif;
      float a;
      a = *pait + ((*tit-*ptit) / (*ptitnext - *ptit)) * (*paitnext - *pait);
      *ait = a;
      }
    }
  else if (maskVolume && ! m_UsePopulationAIF)
    {
    // calculate the AIF from the image using the data under the
    // specified mask

    //std::cout << "Does this work? Test 2 \n"; 

    m_AIF = this->CalculateAverageAIF(inputVectorVolume, maskVolume); //not the problem!

    }
  else if (m_UsePopulationAIF)
    {
    m_AIF = this->CalculatePopulationAIF(m_Timing, 0.1);
    }
  else
    {
    itkExceptionMacro("A mask image over which to establish the AIF or a prescribed AIF must be assigned. If prescribing an AIF, then UsePrescribedAIF must be set to true.");
    }
  // Compute the bolus arrival time
    compute_bolus_arrival_time (m_AIF.size(), &m_AIF[0], m_AIFBATIndex, aif_FirstPeakIndex, aif_MaxSlope);


  // Compute the area under the curve for the AIF
  // probs doesn't work 
  m_aifAUC = area_under_curve(timeSize, &m_Timing[0], &m_AIF[0], m_AIFBATIndex, m_AUCTimeInterval);
  //printf("m_aifAUC = %f\n", m_aifAUC);
  

}

template <class TInputImage, class TMaskImage, class TOutputImage>
void
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::ThreadedGenerateData( const OutputVolumeRegionType& outputRegionForThread, ThreadIdType threadId )
{


  test_func(); 



  VectorVoxelType vectorVoxel, fittedVectorVoxel;

  float tempK2 = 0.0f;
  float tempTTP = 0.0f; 
  float tempK1 = 0.0f;
  float tempBeta = 0.0f; 
  float tempAlpha = 0.0f; 
  float tempMaxSlope = 0.0f;
  float tempAUC = 0.0f;
  float tempMTT = 0.0f;
  float tempCBF = 0.0f;
  int   BATIndex = 0;
  int   FirstPeakIndex = 0;

  const VectorVolumeType* inputVectorVolume = this->GetInput();

  VectorVolumeConstIterType inputVectorVolumeIter(inputVectorVolume, outputRegionForThread);
  //OutputVolumeIterType k2VolumeIter(this->GetK2Output(), outputRegionForThread);
  OutputVolumeIterType ttpVolumeIter(this->GetK2Output(), outputRegionForThread);
  OutputVolumeIterType k1VolumeIter(this->GetK1Output(), outputRegionForThread);

  typename VectorVolumeType::Pointer fitted = this->GetFittedDataOutput();
  VectorVolumeIterType fittedVolumeIter(fitted, outputRegionForThread);

  MaskVolumeConstIterType roiMaskVolumeIter;
  if(this->GetROIMask())
    {
    roiMaskVolumeIter = MaskVolumeConstIterType(this->GetROIMask(), outputRegionForThread);
    }

  OutputVolumeIterType maxSlopeVolumeIter(this->GetMaxSlopeOutput(), outputRegionForThread);
  OutputVolumeIterType aucVolumeIter(this->GetAUCOutput(), outputRegionForThread);
  OutputVolumeIterType mttVolumeIter(this->GetMTTOutput(), outputRegionForThread);
  OutputVolumeIterType cbfVolumeIter(this->GetCBFOutput(), outputRegionForThread);
  OutputVolumeIterType rsqVolumeIter(this->GetRSquaredOutput(), outputRegionForThread);
  OutputVolumeIterType batVolumeIter(this->GetBATOutput(), outputRegionForThread);

  //set up optimizer and cost function
  itk::LevenbergMarquardtOptimizer::Pointer optimizer = itk::LevenbergMarquardtOptimizer::New();
  LMCostFunction::Pointer                   costFunction = LMCostFunction::New();
  
  const int timeSize = (int)inputVectorVolume->GetNumberOfComponentsPerPixel();
  // each pixel has 14 'components' i.e. acquisition times. 

  /* ok here we go: 
  Ca array!
  */ 


 //TODO: 
 //1. check that the m_AIF is actually correct (correct values, is an actual array etc)
 //2. "standard SVD" --> improve by changing to bSVD or oSVD

  vnl_matrix<float> camat(timeSize, timeSize); 

  int i, j; 

  for (i = 0; i < timeSize; i++)  {
    for (j = 0; j < timeSize; j++) {
      if (j > i) {

        camat(i, j) = 0; 

      } else if (i == j) {

        camat(i, j) = m_AIF[0]; //might have to be 0-indexed 
        //yes...
        //"R arrays start at 1"
        //do vnl_matrix's have 0-indexed co-ordinates
        //better check this!
        //yes vnl_matrices are 0-indexed 
        
      } else if (i > j) {

        camat(i, j) = m_AIF[(i-j)]; 

      }
    }
  }

  /* 
  vnl_matrix<float> B(timeSize, timeSize, 0.0); 

  for (i = 0; i < timeSize; i++) {
    for (j = 0; j < timeSize; j++) {

        if (j > i) {
          B(i, j) = m_AIF[timeSize - (j - i)]; //0 - indexed, unlike in R (pretty sure!)
        } else {
          B(i, j) = 0.0; //test 
        }
    }
  }
  */ 

  /* 
  //can I even do a column bind? 
  //shonky cbind/rbind idea: 2 * timeSize, fill out manually.... 
  //actually should look up convolution matrix you dingus! 

  // cmat = rbind(cbind(camat,B),cbind(B,camat))
  vnl_matrix<float> cmat(2 * timeSize, 2 * timeSize); 

  for (i = 0; i < (2 * timeSize); i++) {
    for (j = 0; j < (2 * timeSize); j++) {
        if (i < timeSize && j < timeSize) {
          cmat(i, j) = camat(i, j); 
        } else if (i > timeSize && j > timeSize) {
          cmat(i, j) = camat((i - timeSize), (j - timeSize)); 
        } else if (i < timeSize && j > timeSize) {
          cmat(i, j) = B(i, (j - timeSize)); 
        } else if (i > timeSize && j < timeSize) {
          cmat(i, j) = B((i - timeSize), j); 
        }
    }
  }

  */ 

  //might be a good idea to use different SVD source code
  //one with better documentation!

  vnl_svd<float> svd (camat); 

  vnl_matrix<float> u, v, d; 


  u = svd.U(); 
  v = svd.V(); 
  d = svd.W(); // 'd' is the diagonal matrix, i.e. S in the paper!
  //diagonal matrix with the singular values of Ca 
  //also sorted decreasingly, no problem there! 
  // in vnl, this is not a proper vnl_matrix 
  // this is a "diag_matrix", i.e. a compressed form/ vector!!!!!!!!!!!!!!!!

  //cut off set to some percentage of the maximum singular value of Ca 
  float psvd = 0.05 * d.max_value(); 

  std::cout << "This is the max value: " << d.max_value() << " " << psvd << "finish!\n"; 
  //1383.16 --> seems fine 

  //can use "vnl_diagonal_matrix" as well 
  vnl_vector<float> vector_sr (timeSize, 0); //values initialised as 0!
  vnl_vector<float> vector_d (timeSize); 

  vector_d = d.get_diagonal(); 
  //diagonal is longer than timeSize... 
  //isn't it 

  //this is where the error is I THINK! 


  for (i = 0; i < (timeSize); i++) {

    if (vector_d.get(i) > psvd) {
      vector_sr.put(i, 1/vector_d.get(i)); 
    }

  }

  vnl_diag_matrix<float> sr(vector_sr); //pretty sure I can do this... 
  //double check ! 


  vnl_matrix<float> ca_inv (timeSize, timeSize); 


  //check that this is actually how you multiply! 
  //c'mon!
  
  ca_inv = v * sr * u.transpose(); //fixed !  

  

  std::vector<float> timeMinute;
  timeMinute = m_Timing;
  for(unsigned int i = 0; i < timeMinute.size(); i++)
    {
    timeMinute[i] = m_Timing[i]/60.0;
    }

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());


  VectorVoxelType shiftedVectorVoxel(timeSize); //as in also 14 components 
  int shift;
  unsigned int shiftStart = 0, shiftEnd = 0;
  bool success = true;

  
  vnl_matrix<float> Ct (timeSize, 1); 
  //m x n, m rows, n columns 
  //Ct is just the concentration values in the voxel at each time point ! 
  //vnl_vector<float> Ct(vectorVoxel, timeSize); 

  /* 
  vnl_matrix<float> Ctb ((2 * timeSize), 1, 0); 
  */ 

  vnl_matrix<float> rt_hat((timeSize), 1); //rt_hat is only a single column, for whatever reason


 

  int index = 0; 
  // int k = 0; 
 // int v = 0; can't call it v! v is part of your decomposition !  
  // float sumabsf; 




struct Point {
  double _x, _y;
};

struct Line {
  double _slope, _yInt;
  double getYforX(double x) {
    return _slope*x + _yInt;
  }

  // Construct line from points
  bool fitPoints(const std::vector<Point> &pts) {
    int nPoints = pts.size();
    if( nPoints < 2 ) {
      // Fail: infinitely many lines passing through this single point
      return false;
    }
    double sumX=0, sumY=0, sumXY=0, sumX2=0;
    for(int i=0; i<nPoints; i++) {
      sumX += pts[i]._x;
      sumY += pts[i]._y;
      sumXY += pts[i]._x * pts[i]._y;
      sumX2 += pts[i]._x * pts[i]._x;
    }
    double xMean = sumX / nPoints;
    double yMean = sumY / nPoints;
    double denominator = sumX2 - sumX * xMean;
    // You can tune the eps (1e-7) below for your specific task
    if( std::fabs(denominator) < 1e-7 ) {
      // Fail: it seems a vertical line
      return false;
    }
    _slope = (sumXY - sumX * yMean) / denominator;
    _yInt = yMean - _slope * xMean;
    return true;
  }
};

  
  aucVolumeIter.GoToBegin();
  mttVolumeIter.GoToBegin();
  cbfVolumeIter.GoToBegin(); 
  ttpVolumeIter.GoToBegin(); 
  fittedVolumeIter.GoToBegin();

  while (!aucVolumeIter.IsAtEnd()) {

    /* PLAN: 
      1. try with just shifting for BAT and see what happens 
      2. try to modify existing curve fitting 
        2. a. try to get entire PK solver locally, and modify it and see if it changes anything - worked 
      3. that failing, try and fit using linear regression instead
        3. a. using ITK 
        3. b. using my own implementation 
      4. you can do linear regression using matrixes, faster, and last shot if nothing else works 
      5. also try the different algorithms is Seth Lirette's code and see if they produce a better result - maybe do this first 
      6. try the MATLAB code - see how fast it is - maybe this will be enough anyway 
    */ 

      vectorVoxel = inputVectorVolumeIter.Get();
      fittedVectorVoxel = inputVectorVolumeIter.Get();
      
      tempAUC = tempMaxSlope = tempMTT = tempCBF = tempTTP = tempAlpha = tempBeta = tempK1 = tempK2 = 0.0;

      success = true; 

      tempK2 = tempK1 = tempMaxSlope = tempAUC = tempMTT = tempCBF = 0.0;
      BATIndex = FirstPeakIndex = 0;



      tempTTP = time_to_peak(const_cast<float *>(vectorVoxel.GetDataPointer() ), m_Timing); 
      ttpVolumeIter.Set(static_cast<OutputVolumePixelType>(tempTTP) );
      ++ttpVolumeIter; 

      // linear regression 

      /*
      vnl_matrix<double> ind_var(timeSize, 2);  
      //using MatrixTypeInd = itk::Matrix<float, 14, 2>; 
      //MatrixTypeInd ind_var; 


      vnl_matrix<double> dep_var(timeSize, 1); 
      //using MatrixTypeDep = itk::Matrix<float, 14, 1>; 
      //MatrixTypeDep dep_var; 


      vnl_matrix<double> weights(2, 1); 
      /* 
      using MatrixTypeWeights = itk::Matrix<float, 2, 1>; 
      MatrixTypeWeights weights; 
      */
 

      //starters, assuming t0 = 0, i.e. the bolus arrival is instantaneous 
      //ind_var.fill(0.0); 

      /* 
      for (index = 0; index < timeSize; index++) {
        ind_var(index, 0) =  1;
        if (m_Timing[index] != 0) {
            ind_var(index, 1) = (1 + log(m_Timing[index]/tempTTP) - (m_Timing[index]/tempTTP));
        } else {
            ind_var(index, 1) = 0.0; 
        }
         
        //row, column - right? 
        //std::cout<< "This is the x value here: " << log(m_Timing[index]/tempTTP) << "\n"; 
        // It may possible be working, just that... 
        // when you cout you're only seeing the initial values, with TTP = 0 because 
        // they are not tissue, just the air around the pt's head 
        // so continue forward and see if you do output something 
        if (vectorVoxel[index] != 0) {
             //dep_var.put(index, log(vectorVoxel[index]));
             dep_var(index, 0) = log(vectorVoxel[index]); 
        } else {
             //dep_var.put(index, 0.0); 
            //dep_var(index, 0) = 0.0; 
            dep_var(index, 0) = 0.0; 
        }
        

        if (tempTTP != 0) {
         // std::cout << " " << ind_var(index, 0) << " " << "\n"; 

        }
      }
      */ 

      std::vector<Point> my_points(timeSize); 

      for (index = 0; index < timeSize; index++) {

        struct Point temp; 
        
        if (m_Timing[index] != 0) {
          
          temp._x = (1 + log(m_Timing[index]/tempTTP) - (m_Timing[index]/tempTTP));

        } else {

          temp._x = 0.0; 

        }

        if (vectorVoxel[index] != 0) {

          temp._y = log(vectorVoxel[index]); 

        } else {

          temp._y = 0.0; 

        }

        my_points[index] = temp; 


      }

      struct Line my_line; 

      //std::cout << my_line.fitPoints(my_points) << "\n"; 

      /*
      if (my_line.fitPoints(my_points)) {
        std::cout << "this is the slope: " << my_line._slope << "\n"; 
        std::cout << "this is the yint: " << my_line._yInt << "\n"; 
      }
      */



      //std::cout << "number of rows of ind_var " << (vnl_inverse(ind_var.transpose() * ind_var) * ind_var.transpose() ).rows() << " " << (vnl_inverse(ind_var.transpose() * ind_var) * ind_var.transpose() ).columns() << " " ;
      //is this transposing this in place or returning a new copy -> probs important! 
      //ind_var transpose has 1 row, timeSize columns 
      //ind var has timeSize rows, 1 column 
      
      //product should have timeSize rows, timeSize columns => instead it's got 1 row x 1 column -> but this is correct

      //i think it's literally just slow... :( 

      //weights = (ind_var.transpose() * ind_var) * ind_var.transpose() * dep_var; 

      //maybe vnl_inverse is the problem 
      //you can tell just by taking it out 

      //all that hair pulling and it turned out just to be 
      //because i was trying to find the log of 0 
      
      /*
      float alpha = weights.get(0); 
      float beta = tempTTP/alpha; 
      */

      //"fitted curve" - lets see if this worked
      //doubt it worked, i think there's a problem with my maths, also just setting log(0) = 0 - not correct 
      //keep posting on stack overflow... it's helping you!

      float step = m_Timing[timeSize - 1]/timeSize; 
      for (index = 0; index < timeSize; index++) {
       // fittedVectorVoxel[index] = pow((step*index), alpha) * exp(-(step*index)/beta);
      
       if (index = 0) {
         //fittedVectorVoxel[index] = weights.get(0); 
       } else {
         //fittedVectorVoxel[index] = weights.get(1) * (step*index);  
       }

        //fittedVectorVoxel[index] = weights(0, 0) + weights(1, 0) * step * index; 
       

      }


      if (tempTTP != 0) {
        
           // std::cout << "Doubt this worked: " << holder(0,0) << holder(0, 1) << holder(1, 0) << holder(1, 1) << "\n";        
      }

     
      //Ct = inimage[i,j,]
      for (index = 0; index < timeSize; index++) {
        Ct(index, 0) = vectorVoxel[index]; 
      } 
      
      /* 
      Ctb.fill(0.0); 

      for (k = 0; k < timeSize; k++) { //did they make a mistake? 
        Ctb(k, 0) = Ct.get(k, 0);  // padding zeros at the end
        //shouldn't this be a copy! rather than, like, isn't this a pointer 
        // think! 
        //about it! 
      }
      */ 

      //MATRICES IN VNL ARE 0-INDEXED!!! 

      rt_hat = (ca_inv * Ct); 

      //this is the oscillating bit, will work on it! 
      
      /* 
      sumabsf = 0.0; 

      //std::cout << "the timesize is" << timeSize << " \n"; 

      //the first thing that needs fixing ! 
      for (index = 2; index < (2 * timeSize) - 1; index++) {
        sumabsf = sumabsf + abs(rt_hat.get(index, 0) - (2 * rt_hat.get(index - 1, 0)) + rt_hat.get(index - 2, 0)); 
        //std::cout << "this is the value " << rt_hat.get(index, 0); 
        //std::cout << "this is the absoulate value" << " " << sumabsf << "\n"; 
      }
      */ 
      
      
      if (check_tissue(const_cast<float *>(vectorVoxel.GetDataPointer()), timeSize)) {
           
          //tempCBF = (sumabsf) * (0.178); 
         
          
          tempCBF = rt_hat.max_value() * 80; 


      } else {

          tempCBF = 0; 

      }


        

      cbfVolumeIter.Set(static_cast<OutputVolumePixelType>(tempCBF) );
      ++cbfVolumeIter;


      //tempMTT = (mean_transit_time(timeSize, &m_Timing[0], const_cast<float *>(vectorVoxel.GetDataPointer() ), BATIndex,  m_AUCTimeInterval) );
      //tempMTT = (time_to_peak(const_cast<float *>(shiftedVectorVoxel.GetDataPointer() ), m_Timing) );
    
      //tempMTT = tempAUC/tempCBF; 

      
      //will this work? 

      float area = 0.0;

      //std::cout << "timing!" << xValues[0] << " " << xValues[1]; 

      int x; 
      
      for (x = 0; x < (timeSize - 1); x++) {
      
        area += ((m_Timing[x + 1] - m_Timing[x]) * (rt_hat.get(x, 0) + rt_hat.get(x+1, 0)))/2;

      }

      if (check_tissue(const_cast<float *>(vectorVoxel.GetDataPointer()), timeSize)) {

          tempMTT = area/rt_hat.max_value(); 
          //needs to be multiplied by 'constant' interval between scans, ▵T
          tempAUC = tempCBF/tempMTT; 

          //test TTP, putting it into 
          //the MTT because I am lazy 
          // tempMTT = m_Timing[rt_hat.arg_max()]; 
          //tempMTT = rt_hat.arg_max(); 
         // std::cout << "These are all the values of rt_hat: "; 
          for (index = 0; index < timeSize; index++) {
            if (rt_hat.arg_max() == rt_hat(index, 0)) {
              tempMTT = m_Timing[index]; 
            }
          }
          /*
          std::cout << "this is max value: " << rt_hat.max_value() << "\n"; 
          std::cout << "this should be the index of the max value " << rt_hat.arg_max() << "\n"; 
          */

      } else {

          tempMTT = 0; 
          tempAUC = 0; 

      }

      
      
      mttVolumeIter.Set(static_cast<OutputVolumePixelType>(tempMTT) );
      ++mttVolumeIter;

      aucVolumeIter.Set(static_cast<OutputVolumePixelType>(tempAUC) ); 
      ++aucVolumeIter;
      
      //actually time to peak, just didn't want to create a new variable! 

      ++inputVectorVolumeIter;
      ++batVolumeIter;

      fittedVolumeIter.Set(fittedVectorVoxel); 
      ++fittedVolumeIter;

      progress.CompletedPixel(); 
      
  } 


}




// Calculate a population AIF.
//
// See "Experimentally-Derived Functional Form for a Population-Averaged High-
// Temporal-Resolution Arterial Input Function for Dynamic Contrast-Enhanced
// MRI" - Parker, Robers, Macdonald, Buonaccorsi, Cheung, Buckley, Jackson,
// Watson, Davies, Jayson.  Magnetic Resonance in Medicine 56:993-1000 (2006)
template <class TInputImage, class TMaskImage, class TOutputImage>
std::vector<float>
ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
::CalculatePopulationAIF(std::vector<float> signalTime, const float bolusArrivalTimeFraction)
{

    // Inputs
    // ------
    // signalTime : sequence time, presumed in units of seconds.
    // bolusArrivalTimeFraction : fractional point between 0 and 1 when the bolus is
    //     desired to arrive.  Choose 0.0 to have it at the very beginning,
    //     1.0 to have it at the end.
    //
    // Outputs
    // -------
    // AIF : arterial input function as a function of time

    std::vector<float> AIF;

    // Make a high resolution timing vector as input to the AIF construction.
    std::vector<float> aif_time(signalTime.size() * 10);
    float final_time_point = signalTime[signalTime.size()-1];
    float resolution = final_time_point / (aif_time.size() - 1);
    for (size_t j = 0; j < aif_time.size(); ++j) {
	    aif_time[j] = resolution * j;
    }

    size_t bolus_arrival_time_idx = (float)aif_time.size() * bolusArrivalTimeFraction;

    size_t n = aif_time.size();
    AIF.resize(n);

    size_t numTimePoints = n - bolus_arrival_time_idx;
    std::vector<float> timeSinceBolus(numTimePoints);


    // t=FR*[0:numTimePoints-1]/60;
    // These time points "start" when the bolus arrives.
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        //timeSinceBolus[j] = FR * j / 60.0;
        timeSinceBolus[j] = aif_time[bolus_arrival_time_idx + j] - aif_time[bolus_arrival_time_idx];
    }

    // Parker
    // defining parameters
    double a1(0.809);
    double a2(0.330);
    double T1(0.17406);
    double T2(0.365);
    double sigma1(0.0563);
    double sigma2(0.132);
    double alpha(1.050);
    double beta(0.1685);
    double s(38.078);
    double tau(0.483);


    // term0=alpha*exp(-beta*t)./(1+exp(-s*(t-tau)));
    // Here the assumption is that time is in minutes, so must scale accordingly.
    // see Parker.
    std::vector<double> term0(numTimePoints);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        term0[j] = alpha * exp(-beta*timeSinceBolus[j]/60.0) 
		 / (1 + exp( -s * (timeSinceBolus[j]/60.0 - tau)));
    }


    // term1=[];
    // term2=[];
    double A1 = a1 / (sigma1 * pow((2*M_PI), 0.5));

    // B1=exp(-(t-T1).^2./(2.*sigma1^2));
    double numerator, denominator;
    std::vector<double> B1(numTimePoints);
    denominator = 2.0 * pow(sigma1, 2.0);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        numerator = -1 * pow(-(timeSinceBolus[j]/60.0 - T1), 2.0);
        B1[j] = exp( numerator / denominator );
    }

    // term1=A1.*B1;
    std::vector<double> term1(numTimePoints);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        term1[j] = A1 * B1[j];
    }


    // A2=a2/(sigma2*((2*pi)^0.5));
    double A2 = a2 / (sigma2 * pow(2*M_PI, 0.5));

    //B2=exp(-(t-T2).^2./(2.*sigma2^2));
    std::vector<double> B2(numTimePoints);
    denominator = 2.0 * pow(sigma2, 2.0);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        numerator = -1 * pow(-(timeSinceBolus[j]/60.0 - T2), 2.0);
        B2[j] = exp(numerator / denominator);
    }

    // term2=A2.*B2;
    std::vector<double> term2(numTimePoints);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        term2[j] = A2 * B2[j];
    }

    // aifPost=term0+term1+term2;
    std::vector<double> aifPost(numTimePoints);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        aifPost[j] = term0[j] + term1[j] + term2[j];
    }

    // Initialize values before bolus arrival.
    for ( size_t j = 0; j < bolus_arrival_time_idx; ++j ) {
        AIF[j] = 0;
    }

    // Shift the data to take into account the bolus arrival time.
    // sp=timeOfBolus+1;
    // AIF(sp:end)=aifPost;
    for ( size_t j = bolus_arrival_time_idx; j < AIF.size(); ++j ) {
        AIF[j] = aifPost[j - bolus_arrival_time_idx];
    }

    // Resample back to signal (sequence) time.
    std::vector<float> rAIF = this->ResampleAIF(aif_time, AIF, signalTime);

    return rAIF;

}


template <class TInputImage, class TMaskImage, class TOutputImage>
std::vector<float>
ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
::ResampleAIF(std::vector<float> t1, std::vector<float> y1, std::vector<float> t2)
{
    // Resample time1, y1 to time2
    size_t timeSize = t2.size();
    std::vector<float> y2(timeSize);

    std::vector<float>::iterator y2it = y2.begin();
    std::vector<float>::iterator t2it = t2.begin();

    std::vector<float>::iterator y1it = y1.begin();
    std::vector<float>::iterator t1it = t1.begin();

    std::vector<float>::iterator y1itnext = y1it;
    y1itnext++;
    std::vector<float>::iterator t1itnext = t1it;
    t1itnext++;

    for (; t2it != t2.end(); ++t2it, ++y2it)
      {
      // Three cases
      // (1) extrapolate the aif on the low end of the range of prescribed timings
      // (2) interpolate the aif
      // (3) extrapolate the aif on the high end of the range of prescribed timings
      //
      // Case (1) is handled implictly by the initialization and conditionals.
      if (*t1it <= *t2it)
        {
        // Case (2) from above)
        // find the prescribed times that straddle the current time to interpolate
        while (*t1itnext < *t2it && t1itnext != t1.end())
          {
          ++t1it;
          ++t1itnext;
          ++y1it;
          ++y1itnext;
          }
        }
      if (t1itnext == t1.end())
        {
        // we'll need to extrapolate (Case (3) from above)
        t1itnext = t1it;
        --t1it;
        y1itnext = y1it;
        --y1it;
        }

      // interpolate aif;
      float a;
      a = *y1it + ((*t2it-*t1it) / (*t1itnext - *t1it)) * (*y1itnext - *y1it);
      *y2it = a;
      }

    return y2;
}


// Calculate average AIF according to the AIF mask 
template <class TInputImage, class TMaskImage, class TOutputImage>
std::vector<float>
ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
::CalculateAverageAIF(const VectorVolumeType*  inputVectorVolume, const MaskVolumeType* maskVolume)
{
  std::vector<float> averageAIF;

  //std::cout << "Does this work? Test 3 \n"; 

  VectorVolumeConstIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetRequestedRegion() );
  MaskVolumeConstIterType  maskVolumeIter(maskVolume, maskVolume->GetRequestedRegion() );

  inputVectorVolumeIter.GoToBegin();
  maskVolumeIter.GoToBegin();

  VectorVoxelType vectorVoxel;
  long            numberVoxels = 0;
  long            numberOfSamples = inputVectorVolume->GetNumberOfComponentsPerPixel();
  averageAIF = std::vector<float>(numberOfSamples, 0.0);

  //std::cout << "Does this work? Test 4 \n"; 

  while (!inputVectorVolumeIter.IsAtEnd() )
    {
    //std::cout << "Does this work? Test 5 \n"; 
    if (maskVolumeIter.Get()!=0) // Mask pixel with value !0 will is part of AIF
      {
      numberVoxels++;
      vectorVoxel = inputVectorVolumeIter.Get();

      for(long i = 0; i < numberOfSamples; i++)
        {
        averageAIF[i] += vectorVoxel[i];
        }
      }
    ++maskVolumeIter;
    ++inputVectorVolumeIter;
    }


  //std::cout << "Does this work? Test 6 \n"; 

  for(long i = 0; i < numberOfSamples; i++)
    {
    averageAIF[i] /= (double)numberVoxels;
    }

  //std::cout << "Does this work? Test 7 \n"; 


  for (std::vector<float>::const_iterator i = averageAIF.begin(); i != averageAIF.end(); ++i)
    //std::cout << *i << ' ';

  return averageAIF;
}


template <class TInputImage, class TMaskImage, class TOutputImage>
void ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::SetTiming(const std::vector<float>& inputTiming)
{
  m_Timing = inputTiming;
}

template <class TInputImage, class TMaskImage, class TOutputImage>
const std::vector<float>& ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::GetTiming()
{
  return m_Timing;
}

template <class TInputImage, class TMaskImage, class TOutputImage>
void ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Function tolerance: " << m_fTol << std::endl;
  os << indent << "Gradient tolerance: " << m_gTol << std::endl;
  os << indent << "Parameter tolerance: " << m_xTol << std::endl;
  os << indent << "Epsilon: " << m_epsilon << std::endl;
  os << indent << "Maximum number of iterations: " << m_maxIter << std::endl;
}

} // end namespace itk

