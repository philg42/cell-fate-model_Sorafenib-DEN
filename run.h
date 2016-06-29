/*
 *  run.h
 *  stem_cells_noqt
 *
 *  Created by Philip Greulich on 08.03.13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _RUN_
#define _RUN_

#include "global.h"

class event;

class run                    
{
public:

	run();
	run(const run& run0);
	run(const vector<vector<double> >& params, int clonenumber0 = -1);
	~run();
	
	run* replicate(int clonenumber0 = -1);
	
	void setup();
	void set_cell_config();
	void set_cell_config(const vector<vector<vector<int> > >& joint_CSD)  ;     // Attention !!! only for 3 celltypes    // not finished !!!
	void set_cell_config(vector<vector<int> > joint_CSD);     // Attention !!! only for 3 celltypes    // not finished !!!
	vector<vector<double> > prepare_initclone();
	
	int go();
	void variate_params(vector<int> var_vector0);
	void switch_params();      
	void switch_single_param(int param_index);
	
	int log_CSD();
	
	double loglikely(const vector<vector<vector<double> > >& data0);
	double loglikely_basal(vector<vector<vector<double> > > data_basal);
	double loglikely_total(vector<vector<vector<double> > > data_total);
	void print_averages(int runctr, int esmblctr);

	double av_clonesize();
	double av_clonesize_basal();
	
	int catch_global_events();
	
	int iterate();
	
	void link_params();
	double max_rate();
	
	vector<double> reprod_rate; 
	vector<double> r;
	vector<double> dediff_rate;
	vector<double> delta;
	vector<double> strat_rate;
	vector<double> shed_rate;
	vector<double> direct_strat_frac;
	vector<double> switchrate_forw;
	vector<double> switchrate_backw;
	vector<double> init2prog_frac;
	
	vector<double> reprod_rate_switch; 
	vector<double> r_switch;
	vector<double> dediff_rate_switch;
	vector<double> delta_switch;
	vector<double> strat_rate_switch;
	vector<double> shed_rate_switch;
	vector<double> direct_strat_frac_switch;
	vector<double> switchrate_forw_switch;
	vector<double> switchrate_backw_switch;
	vector<double> init2prog_frac_switch;
	
	vector<vector<double> > parameters;
	vector<vector<double> > parameters_switch;
	
	int clonenumber;
	
	vector<double> clonesizedistr;
	vector<double> clonesizedistr_norm;
	vector<double> clonesizedistr_basal;
	vector<double> clonesizedistr_basal_norm;
	vector<double> clonesizedistr_supra;
	vector<double> cumuldistr_basal;
	vector<double> cumuldistr_total;
	
	vector<vector<double> > joint_distr;
	vector<vector<double> > joint_distr_norm;
	
	vector<vector<vector<double> > > joint_distr_intm;
	vector<vector<vector<double> > > joint_distr_intm_norm;
	vector<vector<double> > CSD_intm_basal;
	vector<vector<double> > CSD_intm_basal_norm;
	vector<vector<double> > CSD_intm_total;
	vector<vector<double> > CSD_intm_total_norm;
	
	vector<int> clonenumber_data;
	
	vector<cell*> cells;
	vector<Clone> clones;
	
	vector<vector<double> > init_clone;
	
	double rho;
	
	double maxrate;
	double time_sim;
	
	double basal_av;
	double total_av;
	
	int flag_switch;
	int flag_switch_2;
	int logflag;
	int log_av_flag;
	int clonetobig_flag;
	int logtime_flag1;
	int logtime_flag2;
	double floating_clone_frac;
	
	int logctr;
	int evalctr;
	int completed;
	
	list<event> eventlist;
	
private:
	

	
};

#endif
