/*
 *  clone.cpp
 *  popdyn_1D_discr
 *
 *  Created by Philip Greulich on 07.10.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "clone.h"

class cell;

Clone::Clone()
{
	cells.clear();
};

Clone::Clone(vector<cell*>& cells0)
{
	for (int i=0; i<cells0.size(); i++) 
	{
		addcell(cells0[i]);
	};
};

Clone::Clone(const Clone& clone0)
{
	cells=clone0.cells;
}

cell* Clone::addcell(cell* cell0)
{
	cells.push_back(cell0);
	
	return cell0;
};

void Clone::removecell(cell* cell0)
{
	if (cell0->get_clone() != this)                                     // Debug !!!
	{
		cerr << "error: cell0, address " << cell0 << ", misassigned to clone\n";
		exit(1);
	};
	cells.remove(cell0);
};


int Clone::get_size()
{
	return cells.size();
};

int Clone::get_size_prog()
{
	int size=0;
	for (list<cell*>::iterator itr=cells.begin(); itr != cells.end(); itr++)
	{
		if ((*itr)->get_celltype() == 0) 
		{
			size++;
		};
	};
	
	return size;
};


list<cell*>* Clone::get_cells()
{
	return &cells;
};

double Clone::get_av_divctr()
{
	double av=0.0;
	double av_ctr=0;
	if (cells.size() > 0)
	{
		for (list<cell*>::iterator itr = cells.begin(); itr != cells.end(); itr++) 
		{
			if ((*itr)->get_celltype() == 0) 
				{
					av = av + (*itr)->get_division_ctr();
					av_ctr++;
				};
		};
		av = (double)av/av_ctr;
	};

	
	return av;
}	






