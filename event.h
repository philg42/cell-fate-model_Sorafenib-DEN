/*
 *  event.h
 *  stem_cells_lesions
 *
 *  Created by Philip Greulich on 08.05.14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  global.h
 *  popdyn_1D_discr
 *
 *  Created by Philip Greulich on 07.10.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

// -----------------------------
// Datei: global.h
// Deklaration globaler Variablen
// -----------------------------


#ifndef _EVENT_
#define _EVENT_


#include "global.h"

class cell;
class run;

class event                   
{
public:

	event();
	event(cell* cell0, int eventtype0, double time0, run* run0);
	~event();
	
	void execute();
	double get_time();
	
	friend class cell;
	
	
private:
	
	cell* cell_ex;
	int eventtype;
	double time0;
	run* this_run;
};

#endif