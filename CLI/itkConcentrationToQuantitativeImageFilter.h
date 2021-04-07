/*=========================================================================
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConcentrationToQuantitativeImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2012/05/01 14:28:51 $
  Version:   $Revision: 0.0 $
=========================================================================*/
#ifndef __itkConcentrationToQuantitativeImageFilter_h
#define __itkConcentrationToQuantitativeImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkArray.h"
#include "itkArray2D.h"
#include "itkVectorImage.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"
#include "../PkSolver/PkSolver.h"
#include <string>
#include <math.h>

namespace itk
{



/* to make a correction to this estimate for contrast agent leakage, 

CONTRAST AGENT LEAKAGE!

GAMMA VARIATE FUNCTION: (not the same as the gamma function)

y(t) = C0(t - t0)^a x exp^-(t - t0)/b

y(t)= A(t- t0)^alpha . exp(-(t - t0)/beta)

C(t) = K(t - AT)^alpha . e^(-(t - AT)/beta)

valid for t > t0

t0 IS THE DELAY TIME FROM T = 0 TO THE TIME THE FUNCTION BEGINS, A = constant scale factor, alpha and beta are free parameters 
t0 is also called AT i.e. arrival time

i.e. t0 is the bolus arrival time 

can be simplified to (assuming t0 = 0): 

y(t) = A.t^alpha.exp(-t/beta) -> so I guess three parameters, A, alpha, and beta 
LET'S START BY DOING THIS --> MEANS WE DON'T HAVE TO FIND THE BOLUS ARRIVAL TIME AND 
"SHIFT OUR VOXEL VECTOR"

CAN output the fitted function as well so we can see if it works!


AND --> tmax (the time when y(t) is maximum, is just)

y'(tmax) = 0
tmax = alpha.beta 


model fitting, which has the additional advantage
of obtained interpolated values between time points, which is a prerequisite for the
current perfusion analysis methods when dealing with non-equidistantly time sampling
intervals.

The current standard is the gamma-variate model (Thompson et al., 1964). This
model was proposed because it visually resembles a measured tissue attenuation curve





*/ 




class LMCostFunction2: public itk::MultipleValuedCostFunction
{
public:
  typedef LMCostFunction2                    Self;
  typedef itk::MultipleValuedCostFunction   Superclass;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;
  itkNewMacro( Self );
        
  enum { SpaceDimension =  2 };
  unsigned int RangeDimension; 

  enum ModelType { TOFTS_2_PARAMETER = 1 };
        
  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;
  typedef Superclass::MeasureType                 MeasureType, ArrayType;
  typedef Superclass::ParametersValueType         ValueType;
		      
        
  int m_ModelType;

  LMCostFunction2()
  {
  }
        
  void SetModelType (int model) {
    m_ModelType = model;
  }
        
  void SetNumberOfValues(unsigned int NumberOfValues)
  {
    RangeDimension = NumberOfValues;
  }

  //blood concentration curve...     
  //THIS IS THE ARTERIAL INPUT FUNCTION!!!! 
  //I guess we don't need this... 
  void SetCb (const float* cb, int sz) //BloodConcentrationCurve.
  {
    Cb.set_size(sz);
    for( int i = 0; i < sz; ++i )
      Cb[i] = cb[i];
    //std::cout << "Cb: " << Cb << std::endl;
  }
        
  //what is cv? THE SHIFTED VECTOR VOXEL 
  //I.E. IF WE DON'T SHIFT, THE VECTOR VOXEL     
  void SetCv (const float* cv, int sz) //Self signal Y
  {    
    Cv.set_size (sz);
    for (int i = 0; i < sz; ++i)
      Cv[i] = cv[i];
    //std::cout << "Cv: " << Cv << std::endl;
  }
        
  //time !
  //THEY PASS TIME IN MINUTES SO... 
  //MAYBE DO THIS ALSO 
  void SetTime (const float* cx, int sz) //Self signal X
  {
    Time.set_size (sz);
    for( int i = 0; i < sz; ++i )
      Time[i] = cx[i];
    //std::cout << "Time: " << Time << std::endl;
  }
  


  //integral of the AIF from time j - 1 to time j 
  void SetIntCb (const float* cb, const float* cx, int sz) //leakage curve.
  {
    IntCb.set_size (sz);
    for( int i = 0; i < sz; ++i ) {
        IntCb[i] = 0;
        for (int j = 1; j < i; j++ ) 
            IntCb[i]+= (cx[j]-cx[j-1])*(cb[j]+cb[j-1])/2; //integrating using the trapezoidal rule 
        }
    //std::cout << "IntCb: " << IntCb << std::endl;
  }
        
  //for each pixel, we determined K1 and K2 by using a linear least squares fit to 
  //Equation (1)

  /*

    R2(t) = K1.R2(t)bar - K2.int(0 to t) of R2(t)bar 

    whole brain average of nonenhancing voxels and its time integral 

  */
  MeasureType GetValue( const ParametersType & parameters) const
  {
    MeasureType measure(RangeDimension);

    ValueType alpha = parameters[0];
    ValueType beta = parameters[1];
            
    ArrayType betaTerm;
    betaTerm = -alpha/beta*Time;
    ValueType deltaT = Time(1) - Time(0);
    
    if(m_ModelType == TOFTS_2_PARAMETER)
      {
      //probs won't work!
     // measure = Cv - (Power(Time, alpha) * PowerR(Time, beta));
     measure = Cv - 5; 
      //A * t^alpha * exp(-t/beta)
      //i.e. 
      //AND HERE'S THE KICKER 
      /* 
        they've substituted "whole brain average of nonenhancing tissues" with the AIF 
        and then the 'cost' is the vector Voxel (the TDCs) MINUS the function you're trying to fit 
        i.e. the difference between the real values and the fitted function 
        and the optimizer iteratively tries to minimise this cost 
      */ 

      }
            
    return measure; 
  }
  
  MeasureType GetFittedFunction( const ParametersType & parameters) const
  {
    MeasureType measure(RangeDimension);

    ValueType alpha = parameters[0];
    ValueType beta = parameters[1];
            
    ArrayType betaTerm;
    betaTerm = -alpha/beta*Time;
    ValueType deltaT = Time(1) - Time(0);
    
    if(m_ModelType == TOFTS_2_PARAMETER)
      {
      //measure = (Power(Time, alpha) * PowerR(Time, beta));
      measure = Time; 

      }
            
    return measure; 
  }
    
  //Not going to be used
  void GetDerivative( const ParametersType & /* parameters*/,
                      DerivativeType  & /*derivative*/ ) const
  {   
  }
        
  unsigned int GetNumberOfParameters(void) const
  {
    if(m_ModelType == TOFTS_2_PARAMETER)
      {
      return 2;
      }
  }
        
  unsigned int GetNumberOfValues(void) const
  {
    return RangeDimension;
  }
        
protected:
  virtual ~LMCostFunction2(){}
private:
        
  ArrayType Cv, Cb, Time, IntCb, Extime;
        
  ArrayType Convolution(ArrayType X, ArrayType Y) const
  {
    ArrayType Z;
    Z = vnl_convolve(X,Y).extract(X.size(),0);
    return Z;
  };

  ArrayType Power(ArrayType X, ValueType y) const
  {
    ArrayType Z;
    Z.set_size(X.size());
    for (unsigned int i=0; i<X.size(); i++)
      {
      Z[i] = pow(X[i], y);
      }
    return Z;
  };

  ArrayType PowerR(ArrayType X, ValueType y) const
  {
    double e_const = 2.71828; 

    ArrayType Z;
    ArrayType Holder; 
    Z.set_size(X.size());
    for (unsigned int i=0; i<X.size(); i++)
      {
      
      Z[i] = pow(X[i], y);
      }
    return Z;
  };

  /*
  vnl_vector<double> Exponential(vnl_vector<double> X, double y) const
  {
    vnl_vector<double> Z;
    Z.set_size(X.size());
    for (unsigned int i=0; i<X.size(); i++)
      {
      Z.put(i, exp(-X.get(i)));
      }
    return Z;
  };
  */
  

 ArrayType Exponential(ArrayType X, ValueType n) const
  {
    ArrayType Z;
    Z.set_size(X.size());
    
    for (unsigned int i=0; i<X.size(); i++)
      {
      Z[i] = exp(-X[i]);
      }
    return Z;
  };
        
  int constraintFunc(ValueType x) const
  {
    if (x<0||x>1)
      return 1;
    else
      return 0;
            
  };
        
        
};





/** \class ConcentrationToQuantitativeImageFilter
 * \brief Calculates quantitative imaging parameters from concentration curves.
 *
 * This filter computes Pk modelling quantitative images from
 * concentration curves. The input volume is a vector image of
 * concentration curves represented in floating point. (may not be in floating point because I removed the conversion filter) 
 * The output is
 * a series of floating point images of quantitative parameters.
 *
 * An second input, specifying the location of the arterial input
 * function, allows for the calculation to be adjusted for blood
 * verses tissue.
 *
 * \note
 * This work is part of the National Alliance for Medical Image Computing
 * (NAMIC), funded by the National Institutes of Health through the NIH Roadmap
 * for Medical Research, Grant U54 EB005149.
 *
 */

template <class TInputImage, class TMaskImage, class TOutputImage>
class ITK_EXPORT ConcentrationToQuantitativeImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef ConcentrationToQuantitativeImageFilter Self;

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                             VectorVolumeType;
  typedef typename VectorVolumeType::Pointer      VectorVolumePointerType;
  typedef typename VectorVolumeType::ConstPointer VectorVolumeConstPointerType;
  typedef typename VectorVolumeType::PixelType    VectorVolumePixelType;
  typedef typename VectorVolumeType::RegionType   VectorVolumeRegionType;
  typedef typename VectorVolumeType::SizeType     VectorVolumeSizeType;
  typedef itk::ImageRegionConstIterator<VectorVolumeType> VectorVolumeConstIterType;
 // A multi-dimensional iterator templated over image type that walks a region of pixels.


  typedef itk::ImageRegionIterator<VectorVolumeType> VectorVolumeIterType;

  typedef TMaskImage                            MaskVolumeType;
  typedef typename MaskVolumeType::Pointer      MaskVolumePointerType;
  typedef typename MaskVolumeType::ConstPointer MaskVolumeConstPointerType;
  typedef typename MaskVolumeType::PixelType    MaskVolumePixelType;
  typedef typename MaskVolumeType::RegionType   MaskVolumeRegionType;
  typedef typename MaskVolumeType::SizeType     MaskVolumeSizeType;
  typedef itk::ImageRegionConstIterator<MaskVolumeType> MaskVolumeConstIterType;

  typedef TOutputImage                                    OutputVolumeType;
  typedef typename OutputVolumeType::Pointer              OutputVolumePointerType;
  typedef typename OutputVolumeType::ConstPointer         OutputVolumeConstPointerType;
  typedef typename OutputVolumeType::PixelType            OutputVolumePixelType;
  typedef typename OutputVolumeType::RegionType           OutputVolumeRegionType;
  typedef typename OutputVolumeType::IndexType            OutputVolumeIndexType;
  typedef itk::ImageRegionIterator<OutputVolumeType>      OutputVolumeIterType;
  typedef itk::ImageRegionConstIterator<OutputVolumeType> OutputVolumeConstIterType;

  /** Methods to implement smart pointers and work with the itk object factory
   **/
  typedef ProcessObject              ProcessObjectSuperclass;
  typedef SmartPointer< Self >       Pointer;
  //typedef SmartPointer< const Self > ConstPointer;
  typedef itk::DataObject::Pointer   DataObjectPointer;

  typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
  using ProcessObjectSuperclass::MakeOutput;
  virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx);

  typedef itk::VariableLengthVector<float>          VectorVoxelType;

  /** Standard class typedefs. */
  typedef ImageToImageFilter<VectorVolumeType,OutputVolumeType> Superclass;
  typedef SmartPointer<const Self>                  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ConcentrationToQuantitativeImageFilter, ImageToImageFilter );

  /** ImageDimension enumeration */
  itkStaticConstMacro(VectorVolumeDimension, unsigned int,
                      VectorVolumeType::ImageDimension);
  itkStaticConstMacro(MaskVolumeDimension, unsigned int,
                      MaskVolumeType::ImageDimension);
  itkStaticConstMacro(OutputVolumeDimension, unsigned int,
                      OutputVolumeType::ImageDimension);

  /** Set and get the parameters to control the calculation of
  quantified valued */
  itkGetMacro( TE, float);
  itkSetMacro( TE, float);
  itkGetMacro( FA, float);
  itkSetMacro( FA, float);
  itkGetMacro( RGD_relaxivity, float);
  itkSetMacro( RGD_relaxivity, float);
  itkGetMacro( S0GradThresh, float);
  itkSetMacro( S0GradThresh, float);
  itkGetMacro( fTol, float);
  itkSetMacro( fTol, float);
  itkGetMacro( gTol, float);
  itkSetMacro( gTol, float);
  itkGetMacro( xTol, float);
  itkSetMacro( xTol, float);
  itkGetMacro( epsilon, float);
  itkSetMacro( epsilon, float);
  itkGetMacro( maxIter, int);
  itkSetMacro( maxIter, int);
  itkGetMacro( AUCTimeInterval, float);
  itkSetMacro( AUCTimeInterval, float);
  itkGetMacro( ModelType, int);
  itkSetMacro( ModelType, int);
  itkGetMacro( constantBAT, int);
  itkSetMacro( constantBAT, int);
  itkGetMacro( BATCalculationMode, std::string);
  itkSetMacro( BATCalculationMode, std::string);

  void SetTiming(const std::vector<float>& inputTiming);
  const std::vector<float>& GetTiming();

  /// Control whether a prescribed AIF vector is used or whether the
  /// AIF is specified by a mask. If UsePrescribedAIF is true, then
  /// an AIF supplied as a vector is used rather than being derived
  /// from a mask applied to the input concentration values. Default
  /// is off.
  itkSetMacro( UsePrescribedAIF, bool );
  itkGetMacro( UsePrescribedAIF, bool );
  itkBooleanMacro( UsePrescribedAIF );

  /// Control whether a population AIF vector is used.
  itkSetMacro( UsePopulationAIF, bool );
  itkGetMacro( UsePopulationAIF, bool );
  itkBooleanMacro( UsePopulationAIF );

  /// Set a mask to specify where the AIF is be calculated from the
  /// input concentration image.
  void SetAIFMask(const MaskVolumeType* volume);

  /// Get the mask that specifies from where the AIF is calculated
  const TMaskImage* GetAIFMask() const;

  /// Set a mask to specify where the model fit is be calculated from the
  /// input concentration image.
  void SetROIMask(const MaskVolumeType* volume);

  /// Get the mask that specifies from where the model fit is calculated
  const TMaskImage* GetROIMask() const;


  /// Set the AIF as a vector of timing and concentration
  /// values. Timing specified in seconds.
  void SetPrescribedAIF(const std::vector<float>& timing,
                        const std::vector<float>& aif);

  /// Get the prescribed AIF
  const std::vector<float>& GetPrescribedAIF()
  { return m_PrescribedAIF; }

  /// Get the timing of the prescribed AIF (ms)
  const std::vector<float>& GetPrescribedAIFTiming()
  { return m_PrescribedAIFTiming; }

  /// Control whether the output volumes are masked by a threshold on
  /// the R-squared goodness of fit
  itkSetMacro( MaskByRSquared, bool);
  itkGetMacro( MaskByRSquared, bool);
  itkBooleanMacro( MaskByRSquared );

  /// Get the quantitative output images
  TOutputImage* GetK2Output();
  TOutputImage* GetK1Output();
  TOutputImage* GetMaxSlopeOutput();
  TOutputImage* GetAUCOutput();
  TOutputImage* GetMTTOutput();
  TOutputImage* GetCBFOutput();
  TOutputImage* GetRSquaredOutput();
  TOutputImage* GetBATOutput();

  VectorVolumeType* GetFittedDataOutput();

protected:
  ConcentrationToQuantitativeImageFilter();
  ~ConcentrationToQuantitativeImageFilter(){
  }
  void PrintSelf(std::ostream& os, Indent indent) const;

  void BeforeThreadedGenerateData();

/* 
#if ITK_VERSION_MAJOR < 4
  void ThreadedGenerateData( const OutputVolumeRegionType& outputRegionForThread, int threadId );

#else
  void ThreadedGenerateData( const OutputVolumeRegionType& outputRegionForThread,
                             ThreadIdType threadId );

#endif
*/

void ThreadedGenerateData( const OutputVolumeRegionType& outputRegionForThread,
                             ThreadIdType threadId );

  //std::vector<float> CalculatePopulationAIF( const size_t time_of_bolus, std::vector<float> timing );
  std::vector<float> CalculatePopulationAIF( std::vector<float> timing, float bolus_arrival_fraction );
  std::vector<float> CalculateAverageAIF(const VectorVolumeType* inputVectorVolume, const MaskVolumeType* maskVolume);
  std::vector<float> ResampleAIF(std::vector<float> t1, std::vector<float> y1, std::vector<float> t2);

private:
  ConcentrationToQuantitativeImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  float  m_TE;
  float  m_FA;
  float  m_RGD_relaxivity;
  float  m_S0GradThresh;
  float  m_fTol;
  float  m_gTol;
  float  m_xTol;
  float  m_epsilon;
  int    m_maxIter;
  float  m_AUCTimeInterval;
  int    m_AIFBATIndex;
  int    m_ModelType;
  bool   m_MaskByRSquared;
  int m_constantBAT;
  std::string m_BATCalculationMode;

  std::vector<float> m_Timing;

  bool m_UsePopulationAIF;
  bool m_UsePrescribedAIF;
  std::vector<float> m_PrescribedAIF;
  std::vector<float> m_PrescribedAIFTiming;

  // variables to cache information to share between threads
  std::vector<float> m_AIF;
  float  m_aifAUC;
};

}; // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConcentrationToQuantitativeImageFilter.hxx"
#endif

#endif
