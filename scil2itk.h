#ifndef SCIL2ITK_H
#define SCIL2ITK_H

#include <scilPrototypes.h>
#include <cvstructs.h>
#include "itkImportImageFilter.h"

const unsigned int scilDimension = 3;


typedef signed short scilPixelType;
typedef itk::ImportImageFilter< scilPixelType, scilDimension > scilImportFilterType;
typedef itk::Image< scilPixelType, scilDimension > scilImageType;



scilImportFilterType::Pointer scil2ITKImportFilter(IMAGE *img, bool convertValues = true);
void convertImageValues(IMAGE *img);
void convertBackImageValues(IMAGE *img);


template<class ImageType>
void ScilGeoFromImage(typename ImageType::Pointer image, CV_GEOMETRY &g) {
		g.axisX.x = 1; g.axisX.y = 0; g.axisX.z = 0;
		g.axisY.x = 0; g.axisY.y = 1; g.axisY.z = 0;
		g.axisZ.x = 0; g.axisZ.y = 0; g.axisZ.z = 1;
		g.extent.x = image->GetBufferedRegion().GetSize()[0] * image->GetSpacing()[0];
		g.extent.y = image->GetBufferedRegion().GetSize()[1] * image->GetSpacing()[1];
		g.extent.z = image->GetBufferedRegion().GetSize()[2] * image->GetSpacing()[2];
		g.origin.x = image->GetOrigin()[0] + g.extent.x/2;
		g.origin.y = image->GetOrigin()[1] + g.extent.y/2;
		g.origin.z = image->GetOrigin()[2] + g.extent.z/2;
}


template<class ImageType>
void ITK2scilImage(typename ImageType::Pointer image, IMAGE *dest) {
	image->Update();
	ImageType::SizeType size;
	size[0] = ImageWidth(dest); // size along X
	size[1] = ImageHeight(dest); // size along Y
	size[2] = ImageDepth(dest); // size along Z
	const unsigned int numberOfPixels = size[0] * size[1] * size[2];
	if (numberOfPixels == image->GetBufferedRegion().GetNumberOfPixels()) {
		signed short *destData = (signed short *)ImageOutData(dest);
		itk::ImageRegionConstIterator< ImageType > iter( image, image->GetBufferedRegion() );
		unsigned long i=0;
		for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
			destData[i++] = iter.Get();
		}
		CV_GEOMETRY g;
		ScilGeoFromImage( image, g);
		cvGeoSetGeometry(dest, &g);
	} else {
		std::cerr << "ERROR: different number of pixels in scil and itk images:"
		 << " ITK:[" 
			<< image->GetBufferedRegion().GetSize(0) << ", "
			<< image->GetBufferedRegion().GetSize(1) << ", "
			<< image->GetBufferedRegion().GetSize(2) << "] "
		 << " scil:[" 
			<< ImageWidth(dest) << ", "
			<< ImageHeight(dest) << ", "
			<< ImageDepth(dest) << "] "
		 << std::endl;
	}
}





#endif