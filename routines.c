//routines.c

#include "routines.h"

cholmod_common * routines_create_cholmod_common()
{
	//declaration of cholmod common for the computations
	
	cholmod_common *c = (cholmod_common *)malloc( sizeof(cholmod_common) );
	cholmod_start( c );
	c->supernodal = CHOLMOD_SIMPLICIAL;
	c->print = 5;
	c->nmethods=1;
	c->method[0].ordering = 4; 
	c->postorder=FALSE;
	c->final_ll = TRUE;
	return( c );
}

void routines_destroy_cholmod_common( cholmod_common *c )
{
	cholmod_finish( c );
	free( c );
	return;
}
