#include <iostream>
#include <exception>

#include "itkbasics.h"


#include "itkImportImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesWriter.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkResampleImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkScalarToRGBPixelFunctor.h"
#include "itkUnaryFunctorImageFilter.h"


#ifdef useFastMemoryHungryGaussian
#include "itkRecursiveGaussianImageFilter.h"
#else
#include "itkDiscreteGaussianImageFilter.h"
#endif


namespace itkBasic {
  
typedef itk::IdentityTransform< double, Dimension > IdentityTransformType;
typedef itk::LinearInterpolateImageFunction< FilterImageType, double > InterpolatorType;
typedef itk::RecursiveGaussianImageFilter< FilterImageType, FilterImageType > GaussFilterType;
typedef itk::GDCMSeriesFileNames     NamesGeneratorType;
typedef itk::ResampleImageFilter< FilterImageType, FilterImageType > ResampleFilterType;
typedef itk::SubtractImageFilter< FilterImageType, FilterImageType, FilterImageType  > SubtractFilter;

typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef itk::Image<RGBPixelType, Dimension> RGBImageType;
typedef itk::Functor::ScalarToRGBPixelFunctor<unsigned long> ColorMapFunctorType;
typedef itk::UnaryFunctorImageFilter<LabeledImageType, RGBImageType, ColorMapFunctorType> ColorMapFilterType;

typedef itk::PermuteAxesImageFilter< FilterImageType > PermuteAxesFilter;
typedef itk::FlipImageFilter< FilterImageType > FlipFilter;
typedef itk::GDCMImageIO ImageIOType;


void memtest(void) {
	long unsigned size = 0;
	long unsigned asize = 1024*1024*1024;
	long unsigned max = asize;
	long unsigned min = 0;
	std::list< void * > plist;
	bool working = true;
	while( working ) {
		void *t = malloc(asize);
		if (t!=NULL) {
			if (asize == 1024*1024*1024) {
				size += asize;
				plist.push_back( t );
			} else {
				free( t );
				if ((max - min)> 1024) {
					min = asize;
					asize = (asize + max) / 2;
				} else {
					size += asize;
					working = false;
				}
			}
		} else {
			max = asize;
			asize = (asize + min) / 2;
		}
	}
	std::cerr << "could allocate " << size/ (1024.0*1024.0) << " MB" << std::endl;
	for( std::list< void * >::iterator i = plist.begin(); i!=plist.end(); ++i) free( *i );
}	



DicomInputImageType::Pointer getDicomSerie(const FileNamesContainer &filenames, ReaderType *reader, unsigned int scaleValue) {
	ReaderType::Pointer smart_reader;
	if (reader == NULL) {
		smart_reader = ReaderType::New();
		reader = smart_reader;
	}
	ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
	reader->SetImageIO( gdcmImageIO );
	reader->SetFileNames(filenames);
	try {
		reader->Update();
	}catch( itk::ExceptionObject & excep ) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	DicomInputImageType::Pointer result = reader->GetOutput();
	if (scaleValue > 1) result = ImageShrink( result, scaleValue );

	return result;
}



void SeriesReader::getSeriesFileNames(unsigned int num, FileNamesContainer &fc) {
	if (slist.size() > num) fc = slist[num];
	else fc.clear();
};

void SeriesReader::readSeriesData(unsigned int minSlices) {
	nameGenerator->SetUseSeriesDetails( true );
	nameGenerator->SetDirectory( inputDir );
	typedef std::vector< std::string > SeriesIdContainer;
	const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
	SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
	SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
	while( seriesItr != seriesEnd )	{
		FileNamesContainer fileNames = nameGenerator->GetFileNames(seriesItr->c_str());
		if (fileNames.size() >= minSlices) 
			slist.push_back( fileNames );
		seriesItr++;
	}
}






/*
void saveCoroImage(const std::string &fname, FilterImageType::ConstPointer image, int windowWidth, int windowLevel) {
	PermuteAxesFilter::Pointer permutator = PermuteAxesFilter::New();
	permutator->SetInput( image );
	PermuteAxesFilter::PermuteOrderArrayType order;
	order[0] = 0; order[1] = 2; order[2] = 1;
	permutator->SetOrder( order );

	FlipFilter::Pointer flipper = FlipFilter::New();

	flipper->SetInput( permutator->GetOutput() );
	FlipFilter::FlipAxesArrayType faxes;
	faxes[0] = false; faxes[1] = true; faxes[2] = false;
	flipper->SetFlipAxes( faxes );
	try {
		flipper->Update();
	}catch( itk::ExceptionObject & excep ) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	saveImage(fname, flipper->GetOutput(), windowWidth, windowLevel);
}

*/

void saveCoroSubtractionImage(const std::string &fname, FilterImageType::ConstPointer minuent, FilterImageType::ConstPointer subtrahent) {

std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
	SubtractFilter::Pointer subtractor = SubtractFilter::New();
	subtractor->SetInput1( minuent );
	subtractor->SetInput2( subtrahent );
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
	try {
		subtractor->Update();
	}catch( itk::ExceptionObject & excep ) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	PermuteAxesFilter::Pointer permutator = PermuteAxesFilter::New();
	permutator->SetInput( subtractor->GetOutput() );
	PermuteAxesFilter::PermuteOrderArrayType order;
	order[0] = 0; order[1] = 2; order[2] = 1;
	permutator->SetOrder( order );

	FlipFilter::Pointer flipper = FlipFilter::New();

	flipper->SetInput( permutator->GetOutput() );
	FlipFilter::FlipAxesArrayType faxes;
	faxes[0] = false; faxes[1] = true; faxes[2] = false;
	flipper->SetFlipAxes( faxes );
	try {
		flipper->Update();
	}catch( itk::ExceptionObject & excep ) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	FilterImageType::Pointer flipImage = flipper->GetOutput();
	ImageSave(flipImage, fname);
}


/*
FilterImageType::Pointer gaussianScale(const FilterImageType::Pointer sourceImage, double scaleValue, float gaussSigma) {
	FilterImageType::SizeType destSize;
	for(unsigned int i = 0; i < Dimension; ++i)
		destSize[i] = static_cast<FilterImageType::SizeType::SizeValueType>(sourceImage->GetBufferedRegion().GetSize(i) * scaleValue);
	return gaussTransform(sourceImage, destSize, gaussSigma);
}
*/

FilterImageType::Pointer gaussTransform(const FilterImageType::Pointer sourceImage, FilterImageType::SizeType destSize, 
										float gaussSigma, const BaseTransformType *transform) {
std::cerr << "Filter, Transformation and Output start" << std::endl;
	itk::TimeProbesCollectorBase collector;
	collector.Start( "Filter, Transformation and Output" );
	collector.Start( "Filter" );
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;

	FilterImageType::SizeType sourceSize = sourceImage->GetBufferedRegion().GetSize();
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;

	FilterImageType::ConstPointer temp;
#ifdef useFastMemoryHungryGaussian
	{
		GaussFilterType::Pointer filter = GaussFilterType::New();
		filter->SetDirection( 0 );   // 0 --> X direction
		filter->SetOrder( GaussFilterType::ZeroOrder );
		filter->SetSigma( gaussSigma );
		filter->SetInput( sourceImage );
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
		try { filter->Update(); }
		catch( itk::ExceptionObject & excep ) { std::cerr << "Exception caught !" << std::endl << excep << std::endl; }
		temp = filter->GetOutput();
	}
	{
		GaussFilterType::Pointer filter = GaussFilterType::New();
		filter->SetDirection( 1 );   // 1 --> Y direction
		filter->SetOrder( GaussFilterType::ZeroOrder );
		filter->SetSigma( gaussSigma );
		filter->SetInput( temp );
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
		try { filter->Update(); }
		catch( itk::ExceptionObject & excep ) { std::cerr << "Exception caught !" << std::endl << excep << std::endl; }
		temp = filter->GetOutput();
	}
	{
		GaussFilterType::Pointer filter = GaussFilterType::New();
		filter->SetDirection( 2 );   // 2 --> Z direction
		filter->SetOrder( GaussFilterType::ZeroOrder );
		filter->SetSigma( gaussSigma );
		filter->SetInput( temp );
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
		try { filter->Update(); }
		catch( itk::ExceptionObject & excep ) { std::cerr << "Exception caught !" << std::endl << excep << std::endl; }
		temp = filter->GetOutput();
	}
#else
	{
		GaussFilterType::Pointer filter = GaussFilterType::New();
		filter->SetUseImageSpacingOn();
		filter->SetVariance( std::sqrt(gaussSigma) );
		filter->SetInput( sourceImage );
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
		try { filter->Update(); }
		catch( itk::ExceptionObject & excep ) { std::cerr << "Exception caught !" << std::endl << excep << std::endl; }
		temp = filter->GetOutput();
	}
#endif

	collector.Stop( "Filter" );
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetInterpolator( interpolator );

	IdentityTransformType::Pointer identPointer;
	if (transform == NULL) {
			identPointer  = IdentityTransformType::New();
			transform = identPointer;
	}

	resampler->SetTransform( transform );
	FilterImageType::SpacingType destSpacing;
	for(unsigned int i = 0; i < Dimension; ++i) {
		destSpacing[i] = sourceImage->GetSpacing()[i] * float(sourceSize[i]) / float(destSize[i]);
	}
	resampler->SetInput( temp );
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
	resampler->SetSize( destSize );
	resampler->SetOutputSpacing( destSpacing );
	resampler->SetOutputOrigin( sourceImage->GetOrigin() );
	resampler->SetDefaultPixelValue( 0 );
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;

	FilterImageType::Pointer result;
	try { resampler->Update(); result = resampler->GetOutput(); }
	catch( itk::ExceptionObject & excep ) { std::cerr << "Exception caught !" << std::endl << excep << std::endl; }
std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
	collector.Stop( "Filter, Transformation and Output" );
	collector.Report();
std::cerr << "Filter, Transformation and Output end" << std::endl;

	return result;
}



SeriesFilelist getSeriesFileNames(const std::string &inputDir, unsigned int minSlices) {
	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	nameGenerator->SetUseSeriesDetails( true );
	nameGenerator->SetDirectory( inputDir );
	SeriesFilelist slist;

	typedef std::vector< std::string > SeriesIdContainer;
	const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
	SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
	SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
	while( seriesItr != seriesEnd )	{
		FileNamesContainer fileNames = nameGenerator->GetFileNames(seriesItr->c_str());
		if (fileNames.size() >= minSlices) 
			slist.push_back( fileNames );
		seriesItr++;
	}

	return slist;
}





FilterImageType::Pointer getImage(const FileNamesContainer &filenames, DictionaryArray &dictArray) {
	typedef itk::ImageSeriesReader< FilterImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
	reader->SetImageIO( gdcmImageIO );
	reader->SetFileNames(filenames);
	reader->Update();
	ReaderType::DictionaryArrayRawPointer dict;

	dict = reader->GetMetaDataDictionaryArray();
	for( ReaderType::DictionaryArrayType::const_iterator it = dict->begin(); it != dict->end(); ++it) {
		dictArray.push_back( **it );
	}
	
	return reader->GetOutput();
}



void writeLabelImage( LabeledImageType::Pointer image, const std::string &fname) {

	try {
		ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();
		colormapper->SetInput( image );
		colormapper->Update();
		RGBImageType::ConstPointer colorImage = colormapper->GetOutput();
		colormapper = NULL;

		typedef itk::Image< RGBPixelType, 2 > OutputImage2DType;
		typedef itk::ImageSeriesWriter< RGBImageType, OutputImage2DType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		typedef itk::NumericSeriesFileNames NamesGeneratorType;
		NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

		nameGenerator->SetSeriesFormat( fname.c_str() );
		nameGenerator->SetStartIndex( 1 );
		nameGenerator->SetEndIndex( colorImage->GetBufferedRegion().GetSize()[2] );
		writer->SetFileNames( nameGenerator->GetFileNames() );
		writer->SetInput( colorImage );
		writer->Update();

	} catch( itk::ExceptionObject & exp ) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << exp << std::endl;
	}

}

}