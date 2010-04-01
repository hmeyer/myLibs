#include "scil_helper.h"

#include <string>
#include <list>

std::list< std::string > imagelist;



IMAGE *myMakeImage(const char *name, int type, int lenx, int leny, int lenz) {
	IMAGE *im;
	im = make_image((char *)name, type, lenx, leny, lenz, 0 ,0);
	if (im != NULL)	imagelist.push_front( name );
	return im;
}

int registerImageName(const char *name) {
	imagelist.push_front( name );
	return OK;
}



int destroyImage(const char *name) {
	IMAGE *image = get_image_by_name((char *)name, 0);
    if (is_image(image)) {
		scil_printf("destroying image %s\n", name);
        return destroy_image(image);
    }
	imagelist.remove( name );
	return OK;
}


int destroyAllImages(void) {
	while(imagelist.size()) {
		destroyImage( imagelist.begin()->c_str() );
	}
	return OK;
}