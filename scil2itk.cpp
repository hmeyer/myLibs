

#include <cvFunctions.h>

#include <iostream>
#include <exception>

#include "scil2itk.h"

void convertImageValues(IMAGE *img) {
	if (!img) throw std::runtime_error("image is null!");
	const unsigned long numberOfPixels = ImageWidth(img) * ImageHeight(img) * ImageDepth(img);
	signed short *data = (signed short *)ImageOutData(img);
	if (data[0] != 0)
		for(unsigned long i = 0; i < numberOfPixels; ++i) data[i] ^= 0x8000;
	else 
		for(unsigned long i = 0; i < numberOfPixels; ++i) data[i] -= 1000;
}

void convertBackImageValues(IMAGE *img) {
	if (!img) throw std::runtime_error("image is null!");
	const unsigned long numberOfPixels = ImageWidth(img) * ImageHeight(img) * ImageDepth(img);
	signed short *data = (signed short *)ImageOutData(img);
	if (data[0] != -1000)
		for(unsigned long i = 0; i < numberOfPixels; ++i) data[i] ^= 0x8000;
	else 
		for(unsigned long i = 0; i < numberOfPixels; ++i) data[i] += 1000;
}

scilImportFilterType::Pointer scil2ITKImportFilter(IMAGE *img, bool convertValues) {
	if (!img) throw std::runtime_error("image is null!");
	scilImportFilterType::Pointer importFilter = scilImportFilterType::New();
	
	if (convertValues) convertImageValues(img);
	scilImportFilterType::SizeType size;
	size[0] = ImageWidth(img); // size along X
	size[1] = ImageHeight(img); // size along Y
	size[2] = ImageDepth(img); // size along Z
	const unsigned long numberOfPixels = size[0] * size[1] * size[2];

	scilImportFilterType::IndexType start;
	start.Fill( 0 );
	scilImportFilterType::RegionType region;
	region.SetIndex( start );
	region.SetSize( size );
	importFilter->SetRegion( region );

	CV_GEOMETRY g;
	cvGeoGetGeometry(img, &g);
	double origin[ scilDimension ];
	origin[0] = g.origin.x - 0.5 * (g.extent.x * g.axisX.x + g.extent.y * g.axisY.x + g.extent.z * g.axisZ.x); // X coordinate
	origin[1] = g.origin.y - 0.5 * (g.extent.x * g.axisX.y + g.extent.y * g.axisY.y + g.extent.z * g.axisZ.y); // Y coordinate
	origin[2] = g.origin.z - 0.5 * (g.extent.x * g.axisX.z + g.extent.y * g.axisY.z + g.extent.z * g.axisZ.z); // Z coordinate
	importFilter->SetOrigin( origin );
	double spacing[ scilDimension ];
	spacing[0] = g.extent.x / ImageWidth(img); // along X direction
	spacing[1] = g.extent.y / ImageHeight(img); // along Y direction
	spacing[2] = g.extent.z / ImageDepth(img); // along Z direction
	importFilter->SetSpacing( spacing );

	const bool applicationWillOwnTheBuffer = true;
	importFilter->SetImportPointer( (scilPixelType*)ImageOutData(img), numberOfPixels, applicationWillOwnTheBuffer );

	return importFilter;
}


