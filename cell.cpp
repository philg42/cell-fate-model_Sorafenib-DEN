/*
 *  cell.cpp
 *  rev_diff_qt
 *
 *  Created by Philip Greulich on 20.11.12.
 *  Copyright 2012 __MyCompanyName__Cnumber_glob. All rights reserved.
 *
 */

#include "cell.h"

//////////// Constructors /////////////////////////

cell::cell(Clone* clone0, run* run0, int type0)
{
	its_clone=clone0;
	celltype=type0;
	state=0;
	its_clone->addcell(this);
	time_last=-1;
	index=-1;
	thisrun=run0; 
	reprod_step=0;
	division_ctr=0;
    mut_ctr=0;
};

cell::cell()
{
	cerr << "error: non-initialized cell generated\n";
	exit(1);
	its_clone=NULL;
	celltype=-1;
	state=-1;
	time_last=-1;
	index=-1;
	thisrun=NULL;
	reprod_step=0;
	division_ctr=0;
    mut_ctr=0;
};

cell::cell(const cell& cell0)
{
	its_clone = cell0.its_clone;
	celltype = cell0.celltype;
	state = cell0.state;
	time_last = cell0.time_last;
	thisrun = cell0.thisrun;
	its_clone->addcell(this);
	index = -1;
	reprod_step=cell0.reprod_step;
	division_ctr=cell0.division_ctr;
    mut_ctr=cell0.mut_ctr;
};

/////////// End Constructors //////////////////////



cell* cell::divide()
{
	division_ctr++;
    
    drand48();
    if (drand48() < mut_rate0)
    {
        mut_ctr++;
    };
	cell* new_cell_pt = new cell(*this);
	
	double r1;
	double delta1;
	double direct_strat_frac1;
	
	if (state == 1) 
	{
		r1 = thisrun->parameters_switch[1][0];
		delta1 = thisrun->parameters_switch[3][0];
		direct_strat_frac1 = thisrun->parameters_switch[6][0];
	}
	else 
	{
		r1 = thisrun->parameters[1][0];
		delta1 = thisrun->parameters[3][0];
		direct_strat_frac1 = thisrun->parameters[6][0];
	};
	
	thisrun->cells.push_back(new_cell_pt);
	new_cell_pt->set_index(thisrun->cells.size()-1);
	new_cell_pt->set_reprodstep(0);
	
	Cnumber_glob[celltype]++;
	
	mrand48();
	mrand48();
	double rand1 = (norm * mrand48() + 0.5);
	
	if (rand1 < r1*(1.0+delta1))
	{
	}
	else if (rand1 < 2*r1)
	{
		new_cell_pt->turninto(celltype+1);
		turninto(celltype+1);
	}
	else 
	{
		new_cell_pt->turninto(celltype+1);
		if (drand48() < direct_strat_frac1)
		{
			if (celltype == 1) turninto(celltype+1);
			if (new_cell_pt->celltype == 1) new_cell_pt->turninto(new_cell_pt->celltype+1);
		};
	};
	
	
#if instant_strat == 1	
	if (celltype == 1) turninto(celltype+1);
	if (new_cell_pt->celltype == 1) new_cell_pt->turninto(new_cell_pt->celltype+1);
#else
	if (strat_infinity_flag == 1) 
	{
		if (celltype == 1) turninto(celltype+1);
		if (new_cell_pt->celltype == 1) new_cell_pt->turninto(new_cell_pt->celltype+1);
	}
#endif
	
#if event_queue == 1
	create_new_event(celltype);
	new_cell_pt->create_new_event(celltype);
#endif
	
	
	return new_cell_pt;
};



cell* cell::divide(int m1,int m2)
{
	cell* cell2=divide();
	turninto(m1);
	cell2->turninto(m2);
	cell2->set_reprodstep(0);
	
	return cell2;
}

void cell::stratify()
{
	turninto(celltype+1);
#if event_queue == 1
	create_new_event(1); 
#endif
};

void cell::turninto(int m)
{
	Cnumber_glob[celltype]--;
	celltype=m;
	Cnumber_glob[m]++;
};

void cell::switch_state(int s)
{
	state = s;	
}

void cell::kill()
{
//	cerr << "kill! \n";
	if (its_clone == NULL || index < 0 || index > thisrun->cells.size() || its_clone->get_cells()->size()==0)         // Debug !!!
	{
		cerr << "error: cell not assigned to clone or cell list \n";
		exit(1);
	}

	its_clone->removecell(this);
	Cnumber_glob[celltype]--;
	thisrun->cells[index] = thisrun->cells.back();
	thisrun->cells[index]->set_index(index);
	thisrun->cells.pop_back();
	
	delete this;
};

void cell::set_index(int i)
{
	index=i;
};

Clone* cell::get_clone()
{
	return its_clone;
};

int cell::get_celltype()
{
	return celltype;
};

int cell::get_index()
{
	return index;
};

run* cell::get_thisrun()
{
	return thisrun;
};

int cell::get_division_ctr()
{
	return division_ctr;
};

int cell::get_mut_ctr()
{
    return mut_ctr;
};

void cell::set_reprodstep(int reprodstep0)
{
	reprod_step = reprodstep0;
};

event* cell::create_new_event(int eventtype, double eventtime)
{	
	double timestep=0;
	double event_rate=0;
	
	if (eventtype == 0) 
	{
		event_rate = thisrun->reprod_rate[0];
	}
	else if (eventtype == 1)
	{
		event_rate = thisrun->strat_rate[1];
	}
	else if (eventtype == 2)
	{
		event_rate = thisrun->shed_rate[2];
	}
	else { cerr << "error: non-defined cell type\n";};
	
	if (eventtime == -1) 
	{
		drand48();
#if sync_reprod == 1
		timestep = 1.0/event_rate;
#else
		timestep = log(1.0/drand48())/event_rate;
#endif
		eventtime = thisrun->time_sim + timestep; 
	};

	if (thisrun->eventlist.empty()) 
	{
		event event_new(this,celltype,eventtime,thisrun);
		thisrun->eventlist.push_back(event_new);
	}
	else 
	{
		list<event>::reverse_iterator itr = thisrun->eventlist.rbegin();
		for (itr = thisrun->eventlist.rbegin(); itr != thisrun->eventlist.rend() && itr->time0 > eventtime; ++itr) 
		{
			
		};
		
		event event_new(this,celltype,eventtime,thisrun);
		thisrun->eventlist.insert(itr.base(),event_new);
	};
}

//--------------------------------------------------------
//       Update Dynamics: Contains Model Implemention
//--------------------------------------------------------
	
void cell::update()                              // Attention !!! Model-dependent
{
	drand48();
	rand1 = drand48();
	int m = celltype;
	
	double reprod_rate1;
    double reprod_rate2;
	double r1;
	double delta1;
	double strat_rate1;
	double shed_rate1;
	double direct_strat_frac1;
	double stateswitch_rate_forw1;
	double stateswitch_rate_backw1;
	double maxrate1 = thisrun->maxrate;
	
	
	if (state == 1) 
	{
		reprod_rate1 = thisrun->parameters_switch[0][0];
        reprod_rate2 = thisrun->parameters_switch[0][1];
		r1 = thisrun->parameters_switch[1][0];
		delta1 = thisrun->parameters_switch[3][0];
		strat_rate1 = thisrun->parameters_switch[4][1];
		shed_rate1 = thisrun->parameters_switch[5][2];
		direct_strat_frac1 = thisrun->parameters_switch[6][0];
		stateswitch_rate_forw1 = thisrun->parameters_switch[7][0];
		stateswitch_rate_backw1 = thisrun->parameters_switch[8][1];
	}
	else 
	{
		reprod_rate1 = thisrun->parameters[0][0];
        reprod_rate2 = thisrun->parameters[0][1];
		r1 = thisrun->parameters[1][0];
		delta1 = thisrun->parameters[3][0];
		strat_rate1 = thisrun->parameters[4][1];
		shed_rate1 = thisrun->parameters[5][2];
		direct_strat_frac1 = thisrun->parameters[6][0];
		stateswitch_rate_forw1 = thisrun->parameters[7][0];
		stateswitch_rate_backw1 = thisrun->parameters[8][1];
	};
	
//	cerr << "stateswitch_rate_forw1 = " << stateswitch_rate_forw1 << endl;
//	cerr << "stateswitch_rate_backw1 = " << stateswitch_rate_backw1 << endl;
//	cerr << "maxrate1 in cell() = " << maxrate1 << endl;
//	cerr << "reprod_rate1=" << reprod_rate1 << endl;
//	cerr << "r = " << r1 << endl;
//	cerr << "delta1 = " << delta1 << endl;
//	cerr << "strat_rate = " << strat_rate1 << endl;
//	cerr << "shed_rate = " << shed_rate1 << endl << endl;


	
//	cerr << maxrate1 << '\t' << reprod_rate1 << '\t' << r1 << '\t' << strat_rate1 << '\t' << shed_rate1 << endl;
	
    if (m == 0)
    {
        
#if many_step_process == 1
        if (reprod_step > reprod_step_number)
        {
            reprod_step = 0;
        };
        if (rand1 < reprod_rate1)
        {
            reprod_step++;
        };
#endif
        
#if many_step_process == 1
        if (reprod_step == reprod_step_number)
        {
#endif
            if (rand1 < (double)reprod_rate1/maxrate1)  // Division
            {
                divide();
                
#if many_step_process == 1
                reprod_step = 0;
#endif
            }
#if many_step_process == 1
            reprod_step = 0;
#endif
            
#if many_step_process == 1
            
        };
#endif
        //		else if (state == 0)
        //		{
        //			if (rand1 < (double)(reprod_rate1 + stateswitch_rate_forw1)/maxrate1)
        //			{
        //				switch_state(state+1);
        //			};
        //		}
        //		else if (state == 1)
        //		{
        //			if (rand1 < (double)(reprod_rate1 + stateswitch_rate_backw1)/maxrate1)
        //			{
        //				switch_state(state-1);
        //			};
        //		}
        else if (rand1 < (double)(reprod_rate1 + stateswitch_rate_forw1)/maxrate1)
        {
            turninto(m+1);
        };
        
        
    }
    else if (m == 1)
    {
        if (rand1 < (double)strat_rate1/maxrate1)
        {
            stratify();
        }
        else if (rand1 < (double)(strat_rate1 + stateswitch_rate_backw1)/maxrate1)
        {
            turninto(m-1);
        }
        else if (rand1 < (double)(strat_rate1 + stateswitch_rate_backw1 + reprod_rate2)/maxrate1)  // Division
        {
            divide();
        }
        
    }
    else
    {
        //		if (m != 2)                                   // Debug !!!
        //		{
        //			cerr << "error: not a suprabasal differentiated cell\n";
        //			exit(1);
        //		}
        if (rand1 < (double)shed_rate1/maxrate1)
        {
            kill();
        };
    }

	
};
	
	//       End Update Dynamics
	// -------------------------------------------------
	

	
