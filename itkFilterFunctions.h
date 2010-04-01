#ifndef __itkFilterFunctions_h
#define __itkFilterFunctions_h

#include "itkConceptChecking.h"
#include <numeric>

namespace itk
{
  

	template <class TInputPixel, class TOutputPixel>
	class MaximumProjector {
		public:
		MaximumProjector():value( itk::NumericTraits< TOutputPixel >::min() ) {}
		void Clear() { value = itk::NumericTraits< TOutputPixel >::min(); }
		void AddValue( const TInputPixel& val) { value = std::max( value, static_cast<TOutputPixel>( val ) ); }
		TOutputPixel GetProjectedValue(void) { return value; }
		bool operator!=(MaximumProjector<TInputPixel, TOutputPixel> const& obj)
		{  return false; }	
		private:
		TOutputPixel value;
	};
	template <class TInputPixel, class TOutputPixel>
	std::ostream& operator<<(std::ostream& stream, MaximumProjector<TInputPixel, TOutputPixel> const& obj)
		{  stream << "MaximumProjector" << std::endl; return stream; }	


	template <class TInputPixel, class TOutputPixel>
	class softMipProjector {
		public:
		typedef std::map< double, double  > softMipVertexContainer;
		typedef typename NumericTraits<TInputPixel>::RealType InputRealType;
		softMipProjector():MaxPosition(0.0), MinPosition(1.0) {}
		void Clear() { pixels.clear(); }
		void AddValue( const TInputPixel& val) { pixels.push_back( val ); }
		TOutputPixel GetProjectedValue(void) {
			typename PixelContainer::size_type quantil_begin = static_cast<typename PixelContainer::size_type>(pixels.size() * MinPosition);
			typename PixelContainer::size_type quantil_end = static_cast<typename PixelContainer::size_type>(pixels.size() * MaxPosition);
			std::nth_element( pixels.begin(), pixels.begin() + quantil_begin, pixels.end(), std::greater< TInputPixel >() );
			std::partial_sort( pixels.begin()+quantil_begin+1, pixels.begin() + quantil_end, pixels.end(), std::greater< TInputPixel >() );
			
			InputRealType val = 0.0;
			typename PixelContainer::const_iterator pixit = pixels.begin();

			for(unsigned int i = 0; i < pixels.size(); ++i) 
				val += *pixit++ * GetIntegral( i );
			return static_cast<TOutputPixel>(val); 
		}
		void AddVertex( double position, double value) {
			Integrals.clear();
			if (position >= 0.0 && position <= 1.0 ) {
				softMipFunction[ position ] = value;
				MaxPosition = std::max( position, MaxPosition );
				MinPosition = std::min( position, MinPosition );
			}
		}
		bool operator!=(softMipProjector<TInputPixel, TOutputPixel> const& obj)
		{  return false; }	
		private:
		double GetValue( double pos, softMipVertexContainer::const_iterator &next ) const {
			if (softMipFunction.empty()) {
				next = softMipFunction.end();
				return 0.0;
			}
			softMipVertexContainer::const_iterator it = softMipFunction.lower_bound( pos );
			if ( it == softMipFunction.end() ) {
				next = softMipFunction.end();
				return (--it)->second;
			}
			next = it;
			if ( it == softMipFunction.begin() ) {
				++next;
				return it->second;
			}
			--it;
			return  
				( it->second + ( next->second - it->second ) * (pos - it->first) / ( next->first - it->first ) );
		}
		void GetValues( double pos_start, double pos_end, softMipVertexContainer &result ) const {
			result.clear();
			softMipVertexContainer::const_iterator start_next;
			softMipVertexContainer::const_iterator end_next;
			result[ pos_start ] = GetValue( pos_start, start_next );
			result[ pos_end ] = GetValue( pos_end, end_next );
			while( start_next != end_next ) {
				result[ start_next->first ] = start_next->second;
				++start_next;
			}
		}
		double GetIntegral( int i) {
			if (Integrals.size() != pixels.size() ) {
				Integrals.resize( pixels.size() );
				double step = 1.0 /pixels.size();
				double pos = 0;
				double sum = 0;
				softMipVertexContainer vcont;
				for(unsigned int k = 0; k < pixels.size(); ++k) {
					GetValues( pos, pos+step, vcont );
					softMipVertexContainer::const_iterator vcontIt = vcont.begin();
					softMipVertexContainer::const_iterator preIt = vcontIt++;
					double val = 0.0;
					while( vcontIt!=vcont.end() ) {
						val += (vcontIt->first - preIt->first) * (vcontIt->second + preIt->second);
						preIt = vcontIt++;
					}
					Integrals[k] = val / 2.0;
					sum += Integrals[k];
					pos+=step;
				}
				double norm = (sum!=0.0)?1.0/sum:0.0;
				for(unsigned int k = 0; k < pixels.size(); ++k)
					Integrals[k] *= norm;
			}
			return Integrals[ i ];
		}
		double MaxPosition;
		double MinPosition;
		typedef std::vector< TOutputPixel > PixelContainer;
		typedef std::vector< double > IntegralContainer;
		IntegralContainer Integrals;
		softMipVertexContainer softMipFunction;
		PixelContainer pixels;
	};
	template <class TInputPixel, class TOutputPixel>
	std::ostream& operator<<(std::ostream& stream, softMipProjector<TInputPixel, TOutputPixel> const& obj)
		{  stream << "softMipProjector" << std::endl; return stream; }	
		
		
	template <class TInputPixel, class TOutputPixel>
	class AverageProjector {
		public:
		typedef typename NumericTraits<TInputPixel>::RealType InputRealType;
		AverageProjector() {}
		void Clear() { pixSum = NumericTraits<InputRealType>::Zero; }
		void AddValue( const TInputPixel& val) { pixSum += val; ++num; }
		TOutputPixel GetProjectedValue(void) {
			return static_cast<TOutputPixel>(pixSum / num); 
		}
		bool operator!=(AverageProjector<TInputPixel, TOutputPixel> const& obj)
		{  return false; }	
		private:
		InputRealType pixSum;
		unsigned long num;
	};
	template <class TInputPixel, class TOutputPixel>
	std::ostream& operator<<(std::ostream& stream, AverageProjector<TInputPixel, TOutputPixel> const& obj)
		{  stream << "AverageProjector" << std::endl; return stream; }	

}



#endif


