//routines.h

#ifndef __ROUTINES_H__
#define __ROUTINES_H__

#include "required_libs.h"
#include "defs.h"

cholmod_common * routines_create_cholmod_common();

void routines_destroy_cholmod_common( cholmod_common *c );

#endif 
