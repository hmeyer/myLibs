#ifndef __itkVariableProjectImageFilter_h
#define __itkVariableProjectImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkConceptChecking.h"
#include <numeric>


#define __itkVariableProjectImageFilter_MULTITHREADED
//#undef __itkVariableProjectImageFilter_MULTITHREADED

namespace itk
{
  
/** \class VariableProjectImageFilter
 * \brief Implements an accumulation of an image along a selected direction.
 *
 *    This class accumulates an image along a dimension and reduce the size 
 * of this dimension to 1. The dimension being accumulated is set by 
 * AccumulateDimension. 
 *
 *   Each pixel is the cumulative sum of the pixels along the collapsed
 * dimension and reduce the size of the accumulated dimension to 1 (only 
 * on the accumulated). 
 *
 *   The dimensions of the InputImage and the OutputImage must be the same.
 *
 *
 *
 * This class is parameterized over the type of the input image and
 * the type of the output image.
 *
 *
 * \author Emiliano Beronich
 *
 * This filter was contributed by Emiliano Beronich
 *
 * \ingroup   IntensityImageFilters     Singlethreaded
 */
template <class TInputImage, class TOutputImage, template< typename , typename > class TProjectorTemplate >
class ITK_EXPORT VariableProjectImageFilter : public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef VariableProjectImageFilter  Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VariableProjectImageFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef typename    InputImageType::Pointer    InputImagePointer;
  typedef typename    InputImageType::RegionType InputImageRegionType;
  typedef typename    InputImageType::PixelType  InputImagePixelType;
  typedef TOutputImage OutputImageType;
  typedef typename     OutputImageType::Pointer    OutputImagePointer;
  typedef typename     OutputImageType::RegionType OutputImageRegionType;
  typedef typename     OutputImageType::PixelType  OutputImagePixelType;
  
  typedef TProjectorTemplate< typename TInputImage::PixelType, typename TOutputImage::PixelType > ProjectorType;


  /** ImageDimension enumeration */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Input and output images must be the same dimension. */
  itkConceptMacro(ImageDimensionCheck,
      (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),
                              itkGetStaticConstMacro(OutputImageDimension)>));

  /** Set the direction in which to accumulate the data.  It must be
   * set before the update of the filter. Defaults to the last
   * dimension. */
  itkGetMacro( AccumulateDimension, unsigned int );
  itkSetMacro( AccumulateDimension, unsigned int );

  void SetProjector( const ProjectorType& proj) { m_Projector = proj; }
  const ProjectorType& GetProjector( void ) const { return m_Projector; }

protected:
  VariableProjectImageFilter();
  virtual ~VariableProjectImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Apply changes to the output image information. */
  virtual void GenerateOutputInformation();

  /** Apply changes to the input image requested region. */
  virtual void GenerateInputRequestedRegion();

  /** This method implements the actual accumulation of the image.
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */
  
#ifndef __itkVariableProjectImageFilter_MULTITHREADED
  void GenerateData(void);
#else  
  void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int threadId);
#endif  

private:
  VariableProjectImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  unsigned int m_AccumulateDimension;
  ProjectorType m_Projector;
};

} // end namespace itk

/*
template <class TInputPixel, class TOutputPixel, const int ptype>
class MinMaxProjector {
	public:
	static const int MaxType = 1;
	static const int MinType = 2;
	MinMaxProjector() { Clear(); }
	void Clear() { value = (ptype == MaxType)?itk::NumericTraits< TOutputPixel >::min() : itk::NumericTraits< TOutputPixel >::max(); }
	void AddValue( const TInputPixel& val) { value = (ptype == MaxType)?
		 std::max( value, static_cast<TOutputPixel>( val ) ) :
		 std::min( value, static_cast<TOutputPixel>( val ) )
		; }
	TOutputPixel GetProjectedValue(void) { return value; }
	bool operator!=(MinMaxProjector<TInputPixel, TOutputPixel, ptype> const& obj)
	{  return obj.value != value; }	
    private:
	TOutputPixel value;
};
template <class TInputPixel, class TOutputPixel, const int ptype>
std::ostream& operator<<(std::ostream& stream, MinMaxProjector<TInputPixel, TOutputPixel, ptype> const& obj)
	{  stream << toString( (ptype == MinMaxProjector<TInputPixel, TOutputPixel, ptype>::MaxType)?"Maximum":"Minimum") << "Projector" << std::endl; return stream; }	


template <class TInputPixel, class TOutputPixel>
class MaximumProjector: public MinMaxProjector< TInputPixel, TOutputPixel, MinMaxProjector<TInputPixel,TOutputPixel,0>::MaxType > {
};

template <class TInputPixel, class TOutputPixel>
class MinimumProjector: public MinMaxProjector< TInputPixel, TOutputPixel, MinMaxProjector<TInputPixel,TOutputPixel,0>::MinType > {
};


template <class TInputPixel, class TOutputPixel>
class softMipProjector {
	public:

    void Print(std::ostream& stream) const {
		stream << "softMipProjector:" << std::endl;
		for( softMipVertexContainer::const_iterator it = softMipFunction.begin(); it != softMipFunction.end(); ++it) {
			stream << "[" << it->first << "," << it->second * NormalizeFactor << "] ";
		}
		stream << std::endl;
		for( double t = 0.0; t <= 1.0; t += 0.1) {
			stream << "[" << t << "," << GetValue( t ) << "] ";
		}
		stream << std::endl;
    }

	softMipProjector():NormalizeFactor(0.0), MaxPosition(0.0) {}
	explicit softMipProjector(double f) {
	  AddVertex( 0.0, 1.0 );
	  if (f != 0.0) AddVertex( f, f * f );
	  else AddVertex( 0.000000001, 0);
	  AddVertex( 1.0, 0.0 );
	}
	void Clear() { pixels.clear(); }
	void AddValue( const TInputPixel& val) { pixels.push_back( val ); }
	TOutputPixel GetProjectedValue(void) {
		typename PixelContainer::size_type quantil_index = static_cast<typename PixelContainer::size_type>(pixels.size() * MaxPosition);
//		std::nth_element( pixels.begin(), pixels.begin() + quantil_index, pixels.end(), std::greater< TInputPixel >() );
		std::partial_sort( pixels.begin(), pixels.begin() + quantil_index, pixels.end(), std::greater< TInputPixel >() );
		
		double val = 0.0;
		typename PixelContainer::const_iterator pixit = pixels.begin();

		for(int i = 0; i < quantil_index; ++i) 
			val += *pixit++ * GetIntegral( i );
		return static_cast<TOutputPixel>(val); 
	}
	void AddVertex( double position, double value) {
		Integrals.clear();
		if (position >= 0.0 && position <= 1.0 ) {
			softMipFunction[ position ] = value;
			MaxPosition = std::max( position, MaxPosition );
			softMipVertexContainer::const_iterator it = softMipFunction.begin();
			NormalizeFactor = it->first * it->second;
			softMipVertexContainer::const_iterator old = it;
			++it;
			while (it != softMipFunction.end()) {
				NormalizeFactor += ( it->first - old->first) * ( it->second + old->second ) / 2;
				old = it;
				++it;
			} 
			NormalizeFactor += ( 1.0 - old->first) * old->second;
			if (NormalizeFactor != 0.0) NormalizeFactor = 1.0 / NormalizeFactor;
		}
	}
	bool operator!=(softMipProjector<TInputPixel, TOutputPixel> const& obj)
	{  return false; }	
    private:
    double GetValue( double pos ) const {
		if (softMipFunction.empty()) return 0.0;
		softMipVertexContainer::const_iterator it = softMipFunction.lower_bound( pos );
		if ( it == softMipFunction.end() ) return NormalizeFactor * (--it)->second;
		if ( it == softMipFunction.begin() ) return NormalizeFactor * it->second;
		softMipVertexContainer::const_iterator pre = it; --pre;
		return NormalizeFactor * 
			( pre->second + ( it->second - pre->second ) * (pos - pre->first) / ( it->first - pre->first ) );
    }
    double GetIntegral( int i) {
		if (Integrals.size() != pixels.size() ) {
			Integrals.resize( pixels.size() );
			double last;
			double current = GetValue( 0.0 );
			double step = 1.0 /pixels.size();
			double pos = step;
			for(int k = 0; k < pixels.size(); ++k) {
				last = current;
				current = GetValue( pos );
				Integrals[k] = step * ( last + current ) / 2.0;
				pos += step;
			}
		}
		return Integrals[ i ];
    }


    double NormalizeFactor;
    double MaxPosition;
    typedef std::vector< TOutputPixel > PixelContainer;
    typedef std::vector< double > IntegralContainer;
    IntegralContainer Integrals;
    typedef std::map< double, double  > softMipVertexContainer;
    softMipVertexContainer softMipFunction;
    PixelContainer pixels;
};
template <class TInputPixel, class TOutputPixel>
std::ostream& operator<<(std::ostream& stream, softMipProjector<TInputPixel, TOutputPixel> const& obj)
	{  stream << "softMipProjector" << std::endl; return stream; }	
	
	
template <class TInputPixel, class TOutputPixel>
class AverageProjector {
	public:
	AverageProjector() { Clear(); }
	void Clear() { sum = .0; num = 0; }
	void AddValue( const TInputPixel& val) { sum += val; num++; }
	TOutputPixel GetProjectedValue(void) {
	    return static_cast<TOutputPixel>( sum / num );
	}
	bool operator!=(AverageProjector<TInputPixel, TOutputPixel> const& obj)
	{  return false; }	
    private:
    double sum;
    unsigned int num;
};

template <class TInputPixel, class TOutputPixel>
std::ostream& operator<<(std::ostream& stream, AverageProjector<TInputPixel, TOutputPixel> const& obj)
	{  stream << "AverageProjector" << std::endl; return stream; }	
*/

template<class TImagePointerType, template < typename, typename > class Projector>
TImagePointerType ImageProjector(TImagePointerType image, 
				     Projector<typename TImagePointerType::ObjectType::PixelType, typename TImagePointerType::ObjectType::PixelType> proj = 
					Projector<typename TImagePointerType::ObjectType::PixelType, typename TImagePointerType::ObjectType::PixelType>(),
					int AccumulateDimension = -1) {
        typedef typename TImagePointerType::ObjectType TImageType;
	typedef itk::VariableProjectImageFilter< TImageType, TImageType, Projector > ProjectFilterType;
        typename ProjectFilterType::Pointer pf = ProjectFilterType::New();
	pf->SetProjector( proj );
	if (AccumulateDimension != -1) pf->SetAccumulateDimension( AccumulateDimension );
	pf->SetInput( image );
	pf->Update();
        TImagePointerType PImage;
	PImage = pf->GetOutput();
	return PImage;
}



#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVariableProjectImageFilter.txx"
#endif

#endif


