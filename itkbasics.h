#ifndef ITKBASICS_H
#define ITKBASICS_H

#include "itkImage.h"
#include "itkMetaDataDictionary.h"
#include "itkTransform.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkShrinkImageFilter.h"
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLaplacianSharpeningImageFilter.h>
#include <itkWeightedAddImageFilter.h>
#include <itkExpandImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkImageSeriesWriter.h>
#include <itkImageRegionIterator.h>
#include <itkNumericSeriesFileNames.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkCastImageFilter.h>
#include <itkPermuteAxesImageFilter.h>
#include <itkFlipImageFilter.h>
#include <itkOrImageFilter.h>
#include <itkAndImageFilter.h>
#include <itkNotImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkGDCMImageIO.h>
#include "gdcmUIDGenerator.h"

#define useFastMemoryHungryGaussian

#ifdef useFastMemoryHungryGaussian
#include "itkRecursiveGaussianImageFilter.h"
#else
#include "itkDiscreteGaussianImageFilter.h"
#endif
#include "itkLinearInterpolateImageFunction.h"
#include "itkIdentityTransform.h"


namespace itkBasic {


const unsigned int Dimension = 3;



typedef float FilterPixelType;

typedef signed short DicomInputPixelType;
typedef itk::Image< DicomInputPixelType, Dimension > DicomInputImageType;


typedef itk::Image< FilterPixelType, Dimension > FilterImageType;

typedef unsigned long LabelPixelType;
typedef itk::Image<LabelPixelType, Dimension> LabeledImageType;

typedef bool BinaryPixelType;
typedef itk::Image< BinaryPixelType, Dimension > BinaryImageType;


typedef std::vector<std::string>    FileNamesContainer;
typedef std::vector< FileNamesContainer > SeriesFilelist;
typedef itk::MetaDataDictionary DictionaryType;
typedef std::vector< DictionaryType > DictionaryArray;
typedef itk::Transform< double, Dimension, Dimension> BaseTransformType;
typedef itk::GDCMSeriesFileNames NamesGeneratorType;
typedef itk::ImageSeriesReader< DicomInputImageType > ReaderType;



enum Orientation {
	NoChange,
	Coronal,
	Saggital,
	RotXY,
};


template <class T1, class T2>
bool equal(T1 input1, T2 input2) {
	return false;
}

template <class T>
bool equal(T input1, T input2) {
	return input1 == input2;
}



FilterImageType::Pointer itkRegister(FilterImageType::Pointer fixedImage, FilterImageType::Pointer movingImage,
									 unsigned int destDim[], float gaussSigma);

void saveCoroSubtractionImage(const std::string &fname, FilterImageType::Pointer minuent, FilterImageType::Pointer subtrahent);

void writeLabelImage( LabeledImageType::Pointer image, const std::string &fname);

DicomInputImageType::Pointer getDicomSerie(const FileNamesContainer &filenames, ReaderType *reader = NULL, unsigned int scaleValue = 1);



void memtest(void);

template<class TImagePointerType>
TImagePointerType ImagePermute(TImagePointerType in_image, Orientation orient) {
	TImagePointerType image;
	typedef typename TImagePointerType::ObjectType TImageType;
	typedef itk::PermuteAxesImageFilter< TImageType > PermuteAxesFilter;
	typedef itk::FlipImageFilter< TImageType > FlipFilter;

	if (orient == NoChange) image = in_image;
	else {
		typename PermuteAxesFilter::PermuteOrderArrayType order;
		typename FlipFilter::FlipAxesArrayType faxes;
		faxes[0] = false; faxes[1] = true; faxes[2] = false;
		bool flipY = true;
		bool flipX = false;
		if (orient == Coronal) {
			order[0] = 0; order[1] = 2; order[2] = 1;
		}
		if (orient == Saggital) {
			order[0] = 1; order[1] = 2; order[2] = 0;
		}
		if (orient == RotXY) {
			order[0] = 1; order[1] = 0; order[2] = 2;
			faxes[0] = true; faxes[1] = false;
		}
		typename PermuteAxesFilter::Pointer permutator = PermuteAxesFilter::New();
		permutator->SetInput( in_image );
		permutator->SetOrder( order );
		typename FlipFilter::Pointer flipper = FlipFilter::New();
		flipper->SetInput( permutator->GetOutput() );
		flipper->SetFlipAxes( faxes );
		try {
			flipper->Update();
		}catch( itk::ExceptionObject & excep ) {
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
			return NULL;
		}
		image = flipper->GetOutput();
	}
	return image;
}

template<class TImagePointerType>
TImagePointerType ImageFlip(TImagePointerType in_image, bool *flipAxes) {
	TImagePointerType image;
	typedef typename TImagePointerType::ObjectType TImageType;
	typedef itk::FlipImageFilter< TImageType > FlipFilter;
//	unsigned int scalefact[ ImageType::ImageDimension ];

	typename FlipFilter::FlipAxesArrayType faxes;
	faxes[0] = false; faxes[1] = true; faxes[2] = false;
	typename FlipFilter::Pointer flipper = FlipFilter::New();
	flipper->SetInput( in_image );
	flipper->SetFlipAxes( flipAxes );
	flipper->Update();
	image = flipper->GetOutput();
	return image;
}

template<class TSourceImageType, class TDestImageType>
typename TDestImageType::Pointer ImageCast(typename TSourceImageType::Pointer in_image) {
	typename TDestImageType::Pointer image;
	typedef itk::CastImageFilter< TSourceImageType, TDestImageType > CastFilter;
	typename CastFilter::Pointer caster = CastFilter::New();
	caster->SetInput( in_image );
	caster->Update();
	image = caster->GetOutput();
	return image;
}

void replaceMetaData( ReaderType::DictionaryArrayRawPointer dictArrayPointer, const std::string tag, const std::string newvalue, bool additive = false);

template<class TImagePointerType>
void writeDicomSeries(TImagePointerType image, const std::string &fileName, typename itk::ImageSeriesReader< typename TImagePointerType::ObjectType >::DictionaryArrayRawPointer dictArrayPointer, int &index) {
	typedef typename TImagePointerType::ObjectType TImageType;
	typedef itk::GDCMImageIO ImageIOType;
	typedef signed short OutputPixelType;
	const unsigned int OutputDimension = 2;
	typedef itk::Image< OutputPixelType, OutputDimension > Image2DType;
	typedef itk::ImageSeriesWriter<	TImageType, Image2DType > SeriesWriterType;
		itk::NumericSeriesFileNames::Pointer namesGenerator = itk::NumericSeriesFileNames::New();
	namesGenerator->SetSeriesFormat(fileName.c_str());
	namesGenerator->SetStartIndex( index );
	if (image->GetImageDimension() == 3)
	  index += image->GetBufferedRegion().GetSize()[2];
	else index += 1;

	namesGenerator->SetEndIndex( index - 1 );

	typename SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
	seriesWriter->SetInput( image );
	ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
	gdcm::UIDGenerator uid;
	uid.SetRoot( gdcmImageIO->GetUIDPrefix() );
	replaceMetaData( dictArrayPointer, "0020|000e", uid.Generate() );
	seriesWriter->SetImageIO( gdcmImageIO );
	seriesWriter->SetFileNames( namesGenerator->GetFileNames() );
	seriesWriter->SetMetaDataDictionaryArray( dictArrayPointer );
	seriesWriter->Update();
}

template<class TImagePointerType>
void writeDicomSeries(TImagePointerType image, const std::string &fileName, ReaderType::Pointer reader, int &index) {
  writeDicomSeries( image, fileName, reader->GetMetaDataDictionaryArray(), index);
}




template<class  TImagePointerType>
void ImageSave_core(TImagePointerType image, const std::string &fname) {
	typedef typename TImagePointerType::ObjectType OutputImageType;
	typedef typename OutputImageType::PixelType OutputPixelType;
	typedef itk::Image< OutputPixelType, 2 > OutputImage2DType;
	typedef itk::ImageSeriesWriter< OutputImageType, OutputImage2DType > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	typedef itk::NumericSeriesFileNames NamesGeneratorType;
	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	nameGenerator->SetSeriesFormat(fname.c_str());
	nameGenerator->SetStartIndex( 1 );
	if (image->GetImageDimension() == 3)
	  nameGenerator->SetEndIndex( image->GetBufferedRegion().GetSize()[2] );
	else nameGenerator->SetEndIndex(1);
	writer->SetFileNames( nameGenerator->GetFileNames() );
	writer->SetInput( image );
	try {
		writer->Update();
	} catch( std::exception & exp ) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << exp.what() << std::endl;
	}
}



template<class TImagePointerType, class TOutImageType >
typename TOutImageType::Pointer ImageRescaleFlexOutput(TImagePointerType image, typename TImagePointerType::ObjectType::PixelType min, typename TImagePointerType::ObjectType::PixelType max) {
	typedef typename TImagePointerType::ObjectType TImageType;
	typename TOutImageType::Pointer rescaledImage;
	typedef itk::RescaleIntensityImageFilter< TImageType, TOutImageType > RescaleFilterType;
	typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	rescaler->SetOutputMinimum(   min );
	rescaler->SetOutputMaximum( max );
	rescaler->SetInput( image );
	try {
		rescaler->Update();
	}catch( itk::ExceptionObject & excep ) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	rescaledImage = rescaler->GetOutput();
	return rescaledImage;
}

template<class TImagePointerType>
TImagePointerType ImageRescale(TImagePointerType image, typename TImagePointerType::ObjectType::PixelType min, typename TImagePointerType::ObjectType::PixelType max) {
    return ImageRescaleFlexOutput<TImagePointerType, typename TImagePointerType::ObjectType>(image, min, max);
}



template<class TImagePointerType>
void ImageSave(TImagePointerType image, const std::string &fname) {
	typedef typename TImagePointerType::ObjectType TImageType;
	typedef unsigned char OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
	OutputImageType::Pointer rescaledImage;
	typedef itk::RescaleIntensityImageFilter< TImageType, OutputImageType > RescaleFilterType;
	typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	rescaler->SetOutputMinimum(   0 );
	rescaler->SetOutputMaximum( 255 );
	rescaler->SetInput( image );
	try {
		rescaler->Update();
	}catch( itk::ExceptionObject & excep ) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	rescaledImage = rescaler->GetOutput();
	ImageSave_core( rescaledImage, fname );
	ImageSave_core( ImageRescaleFlexOutput<TImagePointerType,OutputImageType>(image, 0, 255 ), fname );
}


template<class TImagePointerType>
void ImageSave(TImagePointerType image, const std::string &fname, 
	typename TImagePointerType::ObjectType::PixelType wlevel,  typename TImagePointerType::ObjectType::PixelType wwidth) {
	typedef typename TImagePointerType::ObjectType TImageType;
	typedef typename TImageType::PixelType InputPixelType;
	typedef unsigned char OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
	typedef itk::IntensityWindowingImageFilter< TImageType, OutputImageType > RescaleFilterType;
	typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	typename OutputImageType::Pointer rescaledImage;
	
	InputPixelType minPixel = wlevel-wwidth/2;
	InputPixelType maxPixel = wlevel+wwidth/2;
	
	if ( wwidth==0 ) { // Auto-Window
	  typedef itk::ImageRegionIterator< TImageType > FilterImageIterator;
	  FilterImageIterator it( image, image->GetLargestPossibleRegion() );
	  maxPixel = itk::NumericTraits< InputPixelType >::min();
	  minPixel = itk::NumericTraits< InputPixelType >::max();
	  for(it.GoToBegin(); !it.IsAtEnd(); ++it) {
	    maxPixel = std::max( it.Get(), maxPixel );
	    minPixel = std::min( it.Get(), minPixel );
	  }	  
	}
	rescaler->SetOutputMinimum(   0 );
	rescaler->SetOutputMaximum( 255 );
	rescaler->SetWindowMinimum( minPixel );
	rescaler->SetWindowMaximum( maxPixel );
	rescaler->SetInput( image );
	try {
		rescaler->Update();
	}catch( itk::ExceptionObject & excep ) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	rescaledImage = rescaler->GetOutput();
	ImageSave_core( rescaledImage, fname );
}

template<class TImagePointerType>
TImagePointerType ImageROI(TImagePointerType image, const typename TImagePointerType::ObjectType::RegionType &roi) {
	typedef typename TImagePointerType::ObjectType TImageType;
	typedef itk::RegionOfInterestImageFilter< TImageType, TImageType > RegionOfInterestFilterType;
	typename RegionOfInterestFilterType::Pointer roif = RegionOfInterestFilterType::New();
	roif->SetInput( image );
	roif->SetRegionOfInterest( roi );
	roif->Update();
	TImagePointerType result = roif->GetOutput();
	return result;
}



template<class TImagePointerType, class ScalesType>
TImagePointerType ImageShrink(TImagePointerType input, ScalesType scales) {
	if (equal( scales, 1)) return input;
	typedef typename TImagePointerType::ObjectType ImageType;
	typedef itk::ShrinkImageFilter< ImageType, ImageType > ShrinkFilter;
	typename ShrinkFilter::Pointer shrinker = ShrinkFilter::New();
	shrinker->SetShrinkFactors( scales );
	shrinker->SetInput( input );
	shrinker->Update();
	typename ImageType::Pointer result = shrinker->GetOutput();
	return result;
}

template<class TImagePointerType>
TImagePointerType ImageSharp(TImagePointerType input, float strength) {
	typedef typename TImagePointerType::ObjectType ImageType;
	typedef itk::LaplacianSharpeningImageFilter< ImageType, ImageType > LaplacianFilter;
	typename LaplacianFilter::Pointer lap = LaplacianFilter::New();
	lap->SetInput( input );
	typedef itk::WeightedAddImageFilter< ImageType, ImageType, ImageType> AddFilter;
	typename AddFilter::Pointer adder = AddFilter::New();
	adder->SetInput1( lap->GetOutput() );
	adder->SetInput2( input );
	adder->SetAlpha( strength );
	adder->Update();
	typename ImageType::Pointer result = adder->GetOutput();
	return result;
}



template<class ImagePointerType, class ScalesType, template<typename, typename> class Interpolator>
ImagePointerType ImageExpand(ImagePointerType input, ScalesType scales) {
	if (equal( scales, 1)) return input;
	typedef typename ImagePointerType::ObjectType ImageType;
	typedef itk::ExpandImageFilter< ImageType, ImageType > ExpandFilter;
	typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double > InterpolatorType;

	typename ExpandFilter::Pointer expander = ExpandFilter::New();
	expander->SetInterpolator( InterpolatorType::New() );
	expander->SetInput( input );
	expander->SetExpandFactors( scales );
	expander->Update();
	typename ImageType::Pointer result = expander->GetOutput();
	return result;
}


template<class ImagePointerType, template<typename, typename> class Interpolator>
ImagePointerType ImageExpand( ImagePointerType input, 
		typename ImagePointerType::ObjectType::SpacingType spacing, 
		typename ImagePointerType::ObjectType::SizeType size, 
		typename ImagePointerType::ObjectType::PixelType outsideValue ) {
	typedef typename ImagePointerType::ObjectType ImageType;
	if (spacing == input->GetSpacing() && size == input->GetBufferedRegion().GetSize() ) return input;
	unsigned int scalefact[ ImageType::ImageDimension ];
	bool bloat = false;
	for(int i = 0; i < ImageType::ImageDimension; ++i) {
		double scale = input->GetSpacing()[i] / spacing[i];
		scalefact[i] = static_cast<unsigned int>(scale + 0.5);
		if (std::abs( scalefact[i] - scale) > 0.001) throw std::runtime_error("no integer scaling");
		if (scalefact[i]) bloat = true;
	}
	typename ImageType::Pointer bloated;
	if (bloat) {
		bloated = ImageExpand<ImagePointerType, unsigned int[], Interpolator>( input, scalefact);
	} else {
		bloated = input;
	}
	bool resize = false;
	if (bloated->GetBufferedRegion().GetSize() !=  size) resize = true;
	typename ImageType::Pointer result = ImageType::New();
	typename ImageType::RegionType region;
	region.SetSize( size );
	result->SetRegions( region );
	result->SetSpacing( spacing );
	result->SetOrigin( input->GetOrigin() );
	result->Allocate();
	typedef itk::ImageRegionIterator< ImageType > irItType;
	irItType iter( result, result->GetBufferedRegion() );
	for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
		iter.Set( outsideValue );
	}
	
	typename ImageType::RegionType iterateregion;
	for(int i = 0; i < ImageType::ImageDimension; ++i) {
		iterateregion.SetIndex(i,0);
		iterateregion.SetSize(i, std::min( size[i], bloated->GetBufferedRegion().GetSize(i)) );
	}
	irItType result_it( result, iterateregion );
	irItType bloat_it( bloated, iterateregion );
	result_it.GoToBegin();
	bloat_it.GoToBegin();
	while( !result_it.IsAtEnd() ) {
		result_it.Set( bloat_it.Get() );
		++result_it;
		++bloat_it;
	}
	return result;
}


template<class ImagePointerType, template<typename, typename> class Interpolator>
ImagePointerType ImageScale(ImagePointerType input, float scale) {
	typedef typename ImagePointerType::ObjectType ImageType;
	const float phi = 0.0001;
	assert(scale > phi && "scale is negative or zero");
	if ((scale >= (1.0-phi)) && (scale <= (1.0+phi))) return input;
	if (scale < 1.0) {
		scale = 1.0 / scale;
		assert( (scale - int(scale + phi)) < phi && "scale is not integer or recriproc integer");
		return ImageShrink( input, int( scale ) );
	}
	assert ( (scale - int(scale + phi)) < phi && "scale is not integer or recriproc integer");
	return ImageExpand<ImagePointerType, int, Interpolator>( input, int( scale ) );
}

template<class ImagePointerType, typename ScalesType> 
ImagePointerType ImageExpandNN(ImagePointerType input, ScalesType scales) {
	return ImageExpand<ImagePointerType, ScalesType, itk::NearestNeighborInterpolateImageFunction>( input, scales );
}
template<class ImagePointerType> 
ImagePointerType ImageExpandNN( ImagePointerType input, 
		typename ImagePointerType::ObjectType::SpacingType spacing, 
		typename ImagePointerType::ObjectType::SizeType size, 
		typename ImagePointerType::ObjectType::PixelType outsideValue ) {
	return ImageExpand<ImagePointerType, itk::NearestNeighborInterpolateImageFunction>( input, spacing, size, outsideValue);
}
template<class ImagePointerType> 
ImagePointerType ImageScaleNN(ImagePointerType input, float scale ) {
	return ImageScale<ImagePointerType, itk::NearestNeighborInterpolateImageFunction>( input, scale );
}


template<class ImagePointerType, template <typename,typename,typename> class MorphologicFilter>
ImagePointerType ImageMorphology(ImagePointerType input, int radius, typename ImagePointerType::ObjectType::PixelType foreground) {
	typedef typename ImagePointerType::ObjectType ImageType;
	typedef itk::BinaryBallStructuringElement< typename ImageType::PixelType, ImageType::ImageDimension > StructuringElementType;
	typedef MorphologicFilter< ImageType, ImageType, StructuringElementType > MorphologicFilterType;		

	typename MorphologicFilterType::Pointer morphFilter = MorphologicFilterType::New();			
	StructuringElementType structuringElement;
	structuringElement.SetRadius( radius ); 
	structuringElement.CreateStructuringElement();
	morphFilter->SetKernel( structuringElement );		
	morphFilter->SetInput( input );		
	morphFilter->SetForegroundValue( foreground );
	morphFilter->Update();
	typename ImageType::Pointer result = morphFilter->GetOutput();
	return result;
}

template<class ImagePointerType>
ImagePointerType ImageErosion(ImagePointerType input, int radius, typename ImagePointerType::ObjectType::PixelType foreground) {
	return ImageMorphology<ImagePointerType, itk::BinaryErodeImageFilter>(input, radius, foreground);
}
template<class ImagePointerType>
ImagePointerType ImageDilation(ImagePointerType input, int radius, typename ImagePointerType::ObjectType::PixelType foreground) {
	return ImageMorphology<ImagePointerType, itk::BinaryDilateImageFilter>(input, radius, foreground);
}

template<class InputImageType, class OutputImageType>
typename OutputImageType::Pointer ImageThreshold_TmplResult(typename InputImageType::Pointer input, 
	typename InputImageType::PixelType lowerThreshold, typename InputImageType::PixelType upperThreshold, 
	typename OutputImageType::PixelType insideValue, typename OutputImageType::PixelType outsideValue) {
	typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType> BinaryThresholdFilter;
	typename BinaryThresholdFilter::Pointer binThreshold = BinaryThresholdFilter::New();
	binThreshold->SetUpperThreshold( upperThreshold );
	binThreshold->SetLowerThreshold( lowerThreshold  );
	binThreshold->SetInput( input );
	binThreshold->SetInsideValue( insideValue );
	binThreshold->SetOutsideValue( outsideValue );
	binThreshold->Update();
	typename OutputImageType::Pointer result = binThreshold->GetOutput();
	return result;
}


template<class ImagePointerType>
typename BinaryImageType::Pointer ImageThreshold(ImagePointerType input, 
	typename ImagePointerType::ObjectType::PixelType lowerThreshold, typename ImagePointerType::ObjectType::PixelType upperThreshold, 
	typename BinaryImageType::PixelType insideValue, typename BinaryImageType::PixelType outsideValue) {
	return ImageThreshold_TmplResult<typename ImagePointerType::ObjectType, BinaryImageType>(input, lowerThreshold, upperThreshold, 
		insideValue, outsideValue);
}


template<class ImagePointerType>
typename BinaryImageType::Pointer ImageConnectedThreshold(ImagePointerType input, 
	typename ImagePointerType::ObjectType::PixelType lowerThreshold, typename ImagePointerType::ObjectType::PixelType upperThreshold, 
	typename BinaryImageType::PixelType insideValue, std::list< typename ImagePointerType::ObjectType::PointType > seeds) {
	typedef typename ImagePointerType::ObjectType ImageType;
	typedef itk::ConnectedThresholdImageFilter< ImageType, BinaryImageType> ConnectedThresholdFilter;
	typename ConnectedThresholdFilter::Pointer regionGrow = ConnectedThresholdFilter::New();
	regionGrow->SetUpper( upperThreshold );
	regionGrow->SetLower( lowerThreshold );
	regionGrow->SetInput( input );
	for(typename std::list< typename ImagePointerType::ObjectType::PointType >::const_iterator seed = seeds.begin();
		seed != seeds.end(); ++seed ){
		typename ImageType::IndexType seedIndex;
		input->TransformPhysicalPointToIndex(*seed, seedIndex);	
		regionGrow->AddSeed( seedIndex );
	}
	regionGrow->SetReplaceValue( insideValue );
	regionGrow->Update();		
	typename BinaryImageType::Pointer result = regionGrow->GetOutput();
	return result;
}


template<class ImagePointerType>
typename BinaryImageType::Pointer ImageConnectedThreshold(ImagePointerType input, 
	typename ImagePointerType::ObjectType::PixelType lowerThreshold, typename ImagePointerType::ObjectType::PixelType upperThreshold, 
	typename BinaryImageType::PixelType insideValue, typename ImagePointerType::ObjectType::IndexType seedIndex) {
	typedef typename ImagePointerType::ObjectType ImageType;
	typename ImageType::PointType seed;
	input->TransformIndexToPhysicalPoint( seedIndex, seed );
	return ImageConnectedThreshold( input, lowerThreshold, upperThreshold,
		insideValue, std::list< typename ImageType::PointType >( 1, seed ));
}






template<template <typename,typename,typename> class BinaryFilter>
class ImageBinaryFunctor {
	public:
	template<class Image1PointerType, class Image2PointerType>
	Image1PointerType operator()(Image1PointerType input1, Image2PointerType input2) {
		typedef typename Image1PointerType::ObjectType Image1Type;
		typedef typename Image2PointerType::ObjectType Image2Type;
		typedef BinaryFilter< Image1Type, Image2Type, Image1Type > BinFilter;
		typename BinFilter::Pointer bf = BinFilter::New();
		bf->SetInput1( input1 );
		bf->SetInput2( input2 );
		bf->Update();
		Image1PointerType result = bf->GetOutput();
		return result;
	}
};

template<class Image1PointerType, class Image2PointerType>
Image1PointerType ImageAnd(Image1PointerType input1, Image2PointerType input2) {
	return ImageBinaryFunctor<itk::AndImageFilter>()( input1, input2 );
}
template<class Image1PointerType, class Image2PointerType>
Image1PointerType ImageOr(Image1PointerType input1, Image2PointerType input2) {
	return ImageBinaryFunctor<itk::OrImageFilter>()( input1, input2 );
}


template<class ImagePointerType>
ImagePointerType ImageNot(ImagePointerType input) {
	typedef typename ImagePointerType::ObjectType ImageType;
	typedef itk::NotImageFilter< ImageType, ImageType > NotFilter;
	typename NotFilter::Pointer notf = NotFilter::New();
	notf->SetInput( input );
	notf->Update();
	ImagePointerType result = notf->GetOutput();
	return result;
}




class SeriesReader {
public:
	SeriesReader(const std::string &dir): inputDir(dir), nameGenerator( NamesGeneratorType::New())  { }
	void readSeriesData(unsigned int minSlices=20);
	int numSeries(void) { return slist.size(); }
	void getSeriesFileNames(unsigned int num, FileNamesContainer &fc);
private:
	std::string inputDir;
	NamesGeneratorType::Pointer nameGenerator;
	std::vector< FileNamesContainer > slist;
};

/*
template<class InputImagePointerType, class OutputImagePointerType>
OutputImagePointerType gaussTransform(const InputImagePointerType sourceImage, typename OutputImagePointerType::ObjectType::SizeType destSize, 
										float gaussSigma, const BaseTransformType *transform) {
	typedef typename InputImagePointerType::ObjectType InputImageType;
	typedef typename OutputImagePointerType::ObjectType OutputImageType;
	typedef itk::RecursiveGaussianImageFilter< InputImageType, OutputImageType > GaussFilterType;
	typedef itk::LinearInterpolateImageFunction< OutputImageType, double > InterpolatorType;
	typedef itk::ResampleImageFilter< OutputImageType, OutputImageType > ResampleFilterType;
	typedef itk::IdentityTransform< double, Dimension > IdentityTransformType;
	
#undef DEBUGMSG
#ifdef DEBUGMSG
    std::cerr << "Filter, Transformation and Output start" << std::endl;
	    itk::TimeProbesCollectorBase collector;
	    collector.Start( "Filter, Transformation and Output" );
	    collector.Start( "Filter" );
    std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
#endif
    OutputImageType::SizeType sourceSize = sourceImage->GetBufferedRegion().GetSize();

#ifdef DEBUGMSG
    std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
#endif

	OutputImageType::ConstPointer temp;
#ifdef useFastMemoryHungryGaussian
	{
		GaussFilterType::Pointer filter = GaussFilterType::New();
		filter->SetDirection( 0 );   // 0 --> X direction
		filter->SetOrder( GaussFilterType::ZeroOrder );
		filter->SetSigma( gaussSigma );
		filter->SetInput( sourceImage );
#ifdef DEBUGMSG
    std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
#endif
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
#ifdef DEBUGMSG
    std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
#endif
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
#ifdef DEBUGMSG
    std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
#endif
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
#ifdef DEBUGMSG
    std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
#endif
		try { filter->Update(); }
		catch( itk::ExceptionObject & excep ) { std::cerr << "Exception caught !" << std::endl << excep << std::endl; }
		temp = filter->GetOutput();
	}
#endif

#ifdef DEBUGMSG
collector.Stop( "Filter" );
#endif
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetInterpolator( interpolator );

	IdentityTransformType::Pointer identPointer;
	if (transform == NULL) {
			identPointer  = IdentityTransformType::New();
			transform = identPointer;
	}

	resampler->SetTransform( transform );
	OutputImageType::SpacingType destSpacing;
	for(unsigned int i = 0; i < Dimension; ++i) {
		destSpacing[i] = sourceImage->GetSpacing()[i] * float(sourceSize[i]) / float(destSize[i]);
	}
	resampler->SetInput( temp );
#ifdef DEBUGMSG
    std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
#endif
	resampler->SetSize( destSize );
	resampler->SetOutputSpacing( destSpacing );
	resampler->SetOutputOrigin( sourceImage->GetOrigin() );
	resampler->SetDefaultPixelValue( 0 );
#ifdef DEBUGMSG
    std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
#endif

	OutputImageType::Pointer result;
	try { resampler->Update(); result = resampler->GetOutput(); }
	catch( itk::ExceptionObject & excep ) { std::cerr << "Exception caught !" << std::endl << excep << std::endl; }
#ifdef DEBUGMSG
  std::cerr << "dbg:" << __FILE__ << " line:" << __LINE__ << std::endl;
	  collector.Stop( "Filter, Transformation and Output" );
	  collector.Report();
  std::cerr << "Filter, Transformation and Output end" << std::endl;
#endif
	return result;
}
*/

/*
template<class InputImagePointerType, class OutputImagePointerType>
OutputImagePointerType gaussianScale(const InputImagePointerType sourceImage, double scaleValue, float gaussSigma) {
	typedef typename OutputImagePointerType::ObjectType ImageType;
	OutputImagePointerType::SizeType destSize;
	for(unsigned int i = 0; i < Dimension; ++i)
		destSize[i] = static_cast<ImageType::SizeType::SizeValueType>(sourceImage->GetBufferedRegion().GetSize(i) * scaleValue);
	return gaussTransform<InputImagePointerType,OutputImagePointerType>(sourceImage, destSize, gaussSigma);
}
*/



}

#endif
