/*
 *  cell.h
 *  rev_diff_qt
 *
 *  Created by Philip Greulich on 20.11.12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */




#ifndef _CELL_
#define _CELL_

#include "global.h"

//#include <vector>
//#include <list>
//#include <string>

//using namespace std;

class Clone;
class run;
class event;

class cell                      
{
public:
	 
	cell(Clone* clone0, run* run0, int type0);
	cell();
	cell(const cell& cell0);
	
	cell* divide();
	cell* divide(int m1,int m2);
	void stratify();
	void turninto(int m);
	void switch_state(int s);
	void kill();
	void set_index(int i);
	void set_reprodstep(int reprodstep0);
	Clone* get_clone();
	int get_celltype();
	int get_index();
	run* get_thisrun();
	int get_division_ctr();
    int get_mut_ctr();
	
	void update();
	
	friend class run;

private:

	vector<vector<double> > parameters;
	Clone* its_clone;
	int celltype;
	int state;
	double time_last;
	int index;
	run* thisrun;
	int reprod_step;
	int division_ctr;						// measures predicted H2B concentration (concentration = 1/(2^division_ctr);
    int mut_ctr;
	
	event* create_new_event(int eventtype, double eventtime = -1);
	
};


#endif
