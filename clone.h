/*
 *  clone.h
 *  popdyn_1D_discr
 *
 *  Created by Philip Greulich on 07.10.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */




#ifndef _CLONE_
#define _CLONE_

#include "global.h"

//#include <vector>
//#include <list>
//#include <string>

using namespace std;

class cell;

class Clone                      
{
public:
	
	Clone();
	Clone(vector<cell*>& cells0);
	Clone(const Clone& clone0);
	
	void set_config(vector<int> Cnumber0);
	
	int get_Cnumber(int m);
	int get_size();
	int get_size_prog();
	list<cell*>* get_cells();
	
	cell* addcell(cell* cell0);
	void removecell(cell* cell0);
	
	double get_av_divctr();
	
	
	
private:
	
	list<cell*> cells;
};


#endif
