/*
 *  event.cpp
 *  stem_cells_lesions
 *
 *  Created by Philip Greulich on 08.05.14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "event.h"

event::event()
{
	cell_ex = NULL;
	eventtype=0;
	time0 = 0;
	
};

event::event(cell* cell0, int eventtype0, double time1, run* run0)
{

	cell_ex = cell0;
	eventtype = eventtype0;
	time0 = time1;
	this_run = run0;
	
}

event::~event()
{
};

void event::execute()                  // Attention !!! Depends on details of model
{
	this_run->time_sim = time0;
//	cerr << "event executed\n";
	
	double time_next=0;
	
	if (eventtype == 0)              // Cell Division
	{
		*cell_ex->divide();
	}
	else if (eventtype == 1)
	{
		cell_ex->stratify();
	}
	else if (eventtype == 2)
	{
		cell_ex->kill();
	}
	else 
	{
		cerr << "error: event type not defined\n";
		exit(1);
	};
	
	this_run->eventlist.erase(this_run->eventlist.begin());

};

double event::get_time()
{
	return time0;
};
	