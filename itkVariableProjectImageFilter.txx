#ifndef __itkVariableProjectImageFilter_txx
#define __itkVariableProjectImageFilter_txx

#include "itkVariableProjectImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkStatisticsImageFilter.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage, template< typename , typename > class TProjectorTemplate >
VariableProjectImageFilter<TInputImage,TOutputImage,TProjectorTemplate >
::VariableProjectImageFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  m_AccumulateDimension=InputImageDimension-1;
}


template <class TInputImage, class TOutputImage, template< typename , typename > class TProjectorTemplate >
void
VariableProjectImageFilter<TInputImage,TOutputImage,TProjectorTemplate>
::GenerateOutputInformation()
{
  itkDebugMacro("GenerateOutputInformation Start");

  typename TOutputImage::RegionType outputRegion;
  typename TInputImage::IndexType inputIndex;
  typename TInputImage::SizeType  inputSize;
  typename TOutputImage::SizeType  outputSize;
  typename TOutputImage::IndexType outputIndex;
  typename TInputImage::SpacingType inSpacing;
  typename TInputImage::PointType inOrigin;
  typename TOutputImage::SpacingType outSpacing;
  typename TOutputImage::PointType outOrigin;

  // Get pointers to the input and output
  typename Superclass::OutputImagePointer output = this->GetOutput();
  typename Superclass::InputImagePointer input = const_cast< TInputImage * >( this->GetInput() );

  inputIndex = input->GetLargestPossibleRegion().GetIndex();
  inputSize = input->GetLargestPossibleRegion().GetSize();
  inSpacing = input->GetSpacing();
  inOrigin = input->GetOrigin();

  // Set the LargestPossibleRegion of the output.
  // Reduce the size of the accumulated dimension.
  for(unsigned int i = 0; i<InputImageDimension; i++)
    {
    if (i != m_AccumulateDimension)
      {
      outputSize[i]  = inputSize[i];
      outputIndex[i] = inputIndex[i];
      outSpacing[i] = inSpacing[i];
      outOrigin[i]  = inOrigin[i];
      }
    else
      {
      outputSize[i]  = 1;
      outputIndex[i] = 0;
      outSpacing[i] = inSpacing[i]*inputSize[i];
      outOrigin[i]  = inOrigin[i] + (i-1)*inSpacing[i]/2;
      }
    }

  outputRegion.SetSize(outputSize);
  outputRegion.SetIndex(outputIndex);
  output->SetOrigin(outOrigin);
  output->SetSpacing(outSpacing);
  output->SetLargestPossibleRegion(outputRegion);

  itkDebugMacro("GenerateOutputInformation End");
}


template <class TInputImage, class  TOutputImage, template< typename , typename > class TProjectorTemplate >
void
VariableProjectImageFilter<TInputImage,TOutputImage,TProjectorTemplate>
::GenerateInputRequestedRegion()
{
  itkDebugMacro("GenerateInputRequestedRegion Start");
  Superclass::GenerateInputRequestedRegion();

  if ( this->GetInput() )
    {
    typename TInputImage::RegionType RequestedRegion;
    typename TInputImage::SizeType  inputSize;
    typename TInputImage::IndexType inputIndex;
    typename TInputImage::SizeType  inputLargSize;
    typename TInputImage::IndexType inputLargIndex;
    typename TOutputImage::SizeType  outputSize;
    typename TOutputImage::IndexType outputIndex;

    outputIndex = this->GetOutput()->GetRequestedRegion().GetIndex();
    outputSize = this->GetOutput()->GetRequestedRegion().GetSize();
    inputLargSize = this->GetInput()->GetLargestPossibleRegion().GetSize();
    inputLargIndex = this->GetInput()->GetLargestPossibleRegion().GetIndex();

    for(unsigned int i=0; i<TInputImage::ImageDimension; i++)
      {
      if(i!=m_AccumulateDimension)
        {
        inputSize[i] = outputSize[i];
        inputIndex[i] = outputIndex[i];
        }
      else
        {
        inputSize[i]=inputLargSize[i];
        inputIndex[i]=inputLargIndex[i];
        }
      }

    RequestedRegion.SetSize(inputSize);
    RequestedRegion.SetIndex(inputIndex);
    InputImagePointer input = const_cast< TInputImage * > ( this->GetInput() );
    input->SetRequestedRegion (RequestedRegion);
    }

  itkDebugMacro("GenerateInputRequestedRegion End");
}


/**
 * GenerateData Performs the accumulation
 */
template <class TInputImage, class TOutputImage, template< typename , typename > class TProjectorTemplate >
void
VariableProjectImageFilter<TInputImage,TOutputImage,TProjectorTemplate>
::GenerateData( void )
{
  if(m_AccumulateDimension>=TInputImage::ImageDimension)
    {
    itkExceptionMacro(<<"VariableProjectImageFilter: invalid dimension to accumulate. AccumulateDimension = " << m_AccumulateDimension);
    }

  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef typename NumericTraits<OutputPixelType>::AccumulateType AccumulateType;

  typename Superclass::InputImageConstPointer  inputImage = this->GetInput();
  typename TOutputImage::Pointer outputImage = this->GetOutput();
  outputImage->SetBufferedRegion( outputImage->GetRequestedRegion() );
  outputImage->Allocate();
  

// Accumulate over the Nth dimension ( = m_AccumulateDimension)
// and divide by the size of the accumulated dimension.
  typedef ImageRegionIterator<TOutputImage> outputIterType;
  outputIterType outputIter(outputImage, outputImage->GetBufferedRegion());
  typedef ImageRegionConstIterator<TInputImage> inputIterType;
  
  typename TInputImage::RegionType AccumulatedRegion;
  typename TInputImage::SizeType AccumulatedSize = inputImage->GetLargestPossibleRegion().GetSize();
  typename TInputImage::IndexType AccumulatedIndex = inputImage->GetLargestPossibleRegion().GetIndex();

//  unsigned long SizeAccumulateDimension = AccumulatedSize[m_AccumulateDimension];
  
  long IndexAccumulateDimension = AccumulatedIndex[m_AccumulateDimension];
  for(unsigned int i=0; i< InputImageDimension; i++)
    {
    if (i != m_AccumulateDimension )
      {
      AccumulatedSize[i] = 1;
      }
    }
  AccumulatedRegion.SetSize(AccumulatedSize);
  outputIter.GoToBegin();
  while(!outputIter.IsAtEnd())
    {
    typename TOutputImage::IndexType OutputIndex = outputIter.GetIndex();
    for(unsigned int i=0; i<InputImageDimension; i++)
      {
      if (i != m_AccumulateDimension)
        {
        AccumulatedIndex[i] = OutputIndex[i];
        }
      else
        {
        AccumulatedIndex[i] = IndexAccumulateDimension;
        }
      }
    AccumulatedRegion.SetIndex(AccumulatedIndex);
    inputIterType inputIter(inputImage, AccumulatedRegion);
    inputIter.GoToBegin();
    m_Projector.Clear();
    while(!inputIter.IsAtEnd()) {
	    m_Projector.AddValue( inputIter.Get() );
	    ++inputIter;
	}
    outputIter.Set( static_cast<OutputPixelType>( m_Projector.GetProjectedValue() ) );
    ++outputIter;
    }
}


template <class TInputImage, class TOutputImage, template< typename , typename > class TProjectorTemplate >
void
VariableProjectImageFilter<TInputImage,TOutputImage,TProjectorTemplate>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "AccumulateDimension: " << m_AccumulateDimension << std::endl;
}


} // end namespace itk


#endif
