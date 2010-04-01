#ifndef SCIL_HELPER_H
#define SCIL_HELPER_H

#include <cvFunctions.h>
#include <scilPrototypes.h>
#include <cvSubscribe.h>

#ifdef __cplusplus
extern "C" {
#endif

IMAGE *myMakeImage(const char *name, int type, int lenx, int leny, int lenz);
int registerImageName(const char *name);
int destroyImage(const char *name);
int destroyAllImages(void);


#ifdef __cplusplus
}
#endif

#endif

