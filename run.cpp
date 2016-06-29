/*
 *  run.cpp
 *  stem_cells_noqt
 *
 *  Created by Philip Greulich on 08.03.13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "run.h"


//-----------------------------------------------------------------------
// Constructors
//-----------------------------------------------------------------------

run::run()
{
	cerr << "error: non-initialized run generated\n";
	exit(1);
	reprod_rate.resize(celltypenumber,0);	
	r.resize(celltypenumber,0);	
	dediff_rate.resize(celltypenumber,0);	
	delta.resize(celltypenumber,0);	
	strat_rate.resize(celltypenumber,0);	
	shed_rate.resize(celltypenumber,0);
	direct_strat_frac.resize(celltypenumber,0);	
	switchrate_forw.resize(celltypenumber,0);	
	switchrate_backw.resize(celltypenumber,0);	
	
	reprod_rate_switch.resize(celltypenumber,0);	
	r_switch.resize(celltypenumber,0);	
	dediff_rate_switch.resize(celltypenumber,0);	
	delta_switch.resize(celltypenumber,0);	
	strat_rate_switch.resize(celltypenumber,0);	
	shed_rate_switch.resize(celltypenumber,0);
	direct_strat_frac_switch.resize(celltypenumber,0);	
	switchrate_forw_switch.resize(celltypenumber,0);	
	switchrate_backw_switch.resize(celltypenumber,0);	

	link_params();
	
	time_sim=0;
	flag_switch=0;
	flag_switch_2=0;
	logflag=0;
	log_av_flag=0;
	rho=0;
	clonetobig_flag=0;
	floating_clone_frac=0;
	maxrate=0;
	evalctr=0;
	
	prepare_initclone();
	
	clonenumber = clonenumber_def;
	
	clones.clear();
	clones.resize(clonenumber);
	cells.clear();
	eventlist.clear();
	
	
//	setup();
	
}

run::run(const run& run0)
{
	cells.clear();

	for (int i=0; i<run0.cells.size(); i++) 
	{
		cell* cell_new = new cell(*(run0.cells[i]));
		cells.push_back(cell_new);
	};
	
	clones = run0.clones;
	
	reprod_rate = run0.reprod_rate; 
	r = run0.r;
	dediff_rate = run0.dediff_rate;
	delta = run0.delta;
	strat_rate = run0.strat_rate;
	shed_rate = run0.shed_rate;
	direct_strat_frac = run0.direct_strat_frac;
	switchrate_forw = run0.switchrate_forw;
	switchrate_backw = run0.switchrate_backw;
	
	reprod_rate_switch = run0.reprod_rate_switch; 
	r_switch = run0.r_switch;
	dediff_rate_switch = run0.dediff_rate_switch;
	delta_switch = run0.delta_switch;
	strat_rate_switch = run0.strat_rate_switch;
	shed_rate_switch = run0.shed_rate_switch;
	direct_strat_frac_switch = run0.direct_strat_frac_switch;
	switchrate_forw_switch = run0.switchrate_forw_switch;
	switchrate_backw_switch = run0.switchrate_backw_switch;
	
	link_params();
	
	clonesizedistr = run0.clonesizedistr;
	clonesizedistr_norm = run0.clonesizedistr_norm;
	clonesizedistr_basal = run0.clonesizedistr_basal;
	clonesizedistr_supra = run0.clonesizedistr_supra;
	joint_distr = run0.joint_distr;
	joint_distr_norm = run0.joint_distr_norm;
	
	clonenumber = run0.clonenumber;
	
	joint_distr_intm = run0.joint_distr_intm;
	joint_distr_intm_norm = run0.joint_distr_intm_norm;
	CSD_intm_basal = run0.CSD_intm_basal;
	CSD_intm_basal_norm = run0.CSD_intm_basal_norm;
	
	init_clone = run0.init_clone;
	
	clonenumber_data = run0.clonenumber_data;
	
	floating_clone_frac = run0.floating_clone_frac;
	
	maxrate = run0.maxrate;
	evalctr = run0.evalctr;
	time_sim = run0.time_sim;
	flag_switch = run0.flag_switch;
	flag_switch_2 = run0.flag_switch_2;
	logflag=run0.logflag;
	log_av_flag=run0.log_av_flag;
	clonetobig_flag = run0.clonetobig_flag;
	rho = run0.rho;
	completed = run0.completed;
	eventlist = run0.eventlist;
	
};

run::run(const vector<vector<double> >& params, int clonenumber0)
{
//	reprod_rate.resize(celltypenumber,0);	
//	r.resize(celltypenumber,0);	
//	dediff_rate.resize(celltypenumber,0);	
//	delta.resize(celltypenumber,0);	
//	strat_rate.resize(celltypenumber,0);	
//	shed_rate.resize(celltypenumber,0);
//	direct_strat_frac.resize(celltypenumber,0);
	
	reprod_rate = params[0]; 
	r = params[1];
	dediff_rate = params[2];
	delta = params[3];
	strat_rate = params[4];
	shed_rate = params[5];
	direct_strat_frac = params[6];
	switchrate_forw = params[7];
	switchrate_backw = params[8];
	init2prog_frac = params[9];
#if raw_data == 8
	rho = params[7][0];
#endif
	
	reprod_rate_switch = reprod_rate; 
	r_switch = r;
	dediff_rate_switch = dediff_rate;
	delta_switch = delta;
	strat_rate_switch = strat_rate;
	shed_rate_switch = shed_rate;
	direct_strat_frac_switch = direct_strat_frac;
	switchrate_forw_switch = switchrate_forw;
	switchrate_backw_switch = switchrate_backw;
	init2prog_frac_switch = init2prog_frac;
	
	reprod_rate_switch[0] = switch_reprodrate;
	r_switch[0] = switch_r;
	delta_switch[0] = switch_delta;
	
	link_params();
	
//	cerr << "rho in run(params) = " << rho << endl;       // Debug !!!
	
//	cerr << "new parameters in run(parameters): \nreprod rate = " << reprod_rate[0] << endl << "strat rate = " << strat_rate[1] << endl << "shed rate = " << shed_rate[2] << endl;   // Debug !!!
//	exit(2);	
	
	if (clonenumber0 == -1) clonenumber0 = clonenumber_def;
	
	prepare_initclone();
	
	clonenumber = clonenumber0;
	cerr << "clone number = " << clonenumber << endl;
	
//	if (clonenumber != clonenumber_def)					// Debug !!!
//	{
//		cerr << "error: run with incorrect clone number\n";
//		cerr << "clone number = " << clonenumber << endl;
//		exit(1);
//	}; 
	
	clonesizedistr.clear();
	clonesizedistr_norm.clear();
	clonesizedistr_basal.clear();
	clonesizedistr_supra.clear();
	joint_distr.clear();
	joint_distr_norm.clear();
	
	
	maxrate = max_rate();            // Attention !!! Depends on model details 
//	rho=0;

	cells.clear();
	clones.clear();
	clones.resize(clonenumber);
	eventlist.clear();
	
//	setup();
}

//-----------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------

run::~run()
{
	for (int i=0; i<cells.size(); i++) 
	{
		delete cells[i];
	};	
};


//------------- Link parameters in common pointer vector ----------------

void run::link_params()
{
	parameters.clear();
	parameters.push_back(reprod_rate);
	parameters.push_back(r);
	parameters.push_back(dediff_rate);
	parameters.push_back(delta);
	parameters.push_back(strat_rate);
	parameters.push_back(shed_rate);
	parameters.push_back(direct_strat_frac);
	parameters.push_back(switchrate_forw);
	parameters.push_back(switchrate_backw);
	parameters.push_back(init2prog_frac);
	
	parameters_switch.clear();
	parameters_switch.push_back(reprod_rate_switch);
	parameters_switch.push_back(r_switch);
	parameters_switch.push_back(dediff_rate_switch);
	parameters_switch.push_back(delta_switch);
	parameters_switch.push_back(strat_rate_switch);
	parameters_switch.push_back(shed_rate_switch);
	parameters_switch.push_back(direct_strat_frac_switch);
	parameters_switch.push_back(switchrate_forw_switch);
	parameters_switch.push_back(switchrate_backw_switch);
	parameters_switch.push_back(init2prog_frac_switch);
};

double run::max_rate()
{
	vector<double> maxrate_vect;
	maxrate_vect.clear();
	maxrate_vect.push_back(reprod_rate[0]+strat_rate[0]+switchrate_forw[0]);
	maxrate_vect.push_back(strat_rate[1]+reprod_rate[1]+switchrate_backw[1]);
	maxrate_vect.push_back(shed_rate[2]);
	
//#if raw_data == 10 || raw_data == 18 || raw_data == 24
//	maxrate_vect.push_back(reprod_rate_switch[0]+switchrate_forw_switch[0]);
//	maxrate_vect.push_back(reprod_rate_switch[0]+switchrate_backw_switch[0]);
//	maxrate_vect.push_back(strat_rate_switch[1]);
//	maxrate_vect.push_back(shed_rate_switch[2]);
//#endif
	
	cerr << "maxrate in function = " << *max_element(maxrate_vect.begin(),maxrate_vect.end());
	return *max_element(maxrate_vect.begin(),maxrate_vect.end());
	
};
	

//-----------------------------------------------------------------------
// Replicate run parameters
//-----------------------------------------------------------------------


run* run::replicate(int clonenumber0)
{
	if (clonenumber0 == -1) clonenumber0 = clonenumber;
	
	run* run1 = new run(parameters,clonenumber0);
	
#if raw_data == 2 || raw_data == 3 || raw_data == 4 || raw_data == 22
	
	reprod_rate0_soraf = run1->reprod_rate[0];
	r0_soraf = run1->r[0];
	delta0_soraf = run1->delta[0];
	strat_rate1_soraf=run1->strat_rate[1];
	shedrate2_soraf=run1->shed_rate[2];
	direct_strat_frac0_soraf=run1->direct_strat_frac[0];
	
	run1->reprod_rate[0] = reprod_rate0;
	run1->r[0] = r0;
	run1->delta[0] = delta0;
	run1->shed_rate[2]=shedrate2;
	run1->strat_rate[1]=strat_rate1;
#endif
	
	run1->init_clone=init_clone;
	return run1;
	
};

//-----------------------------------------------------------------------
// Setup Run
//-----------------------------------------------------------------------

vector<vector<double> > run::prepare_initclone()          // init_clone: 1. index i = cell index, 2.index j = cell type, value = probability for cell i to have type j
{
	prog_cell_frac_2 = init2prog_frac[0];
	init_clone.resize(init_clone_number);
	init_clone[0].resize(basalcelltype_number);
	init_clone[0][0]=prog_cell_frac_1;
	init_clone[0][1]=1.0-prog_cell_frac_1;
	if (init_clone.size() > 1)
	{
		init_clone[1].resize(basalcelltype_number);
		init_clone[1][0]=prog_cell_frac_2;
		init_clone[1][1]=1.0-prog_cell_frac_2;		
	};
	
	return init_clone;
};

void run::setup()
{
	time_sim = 0;
	flag_switch = 0;
	flag_switch_2 = 1;                // If == 0 all parameters are switched, else only 1 paramter is switched
	logflag = 0;
	clonetobig_flag = 0;
	logtime_flag1 = 0;
	logtime_flag2 = 0;
	logctr=0;
	evalctr=0;
	clonenumber_data.resize(timepoints);
	
	file_av_clonesize_basal.close();
	file_av_clonesize_basal.open("av_clonesize_basal.dat");
	
	for (int t=0; t<timepoints; t++)
	{
		clonenumber_data[t] = clonenumber_exp[t];
	};

	
	completed=0;
	floating_clone_frac=0;
	
	cells.clear();
	cells.reserve(1000000);                   // Attention !!! Much less for simulating experiments
	clones.clear();
	clones.resize(clonenumber);
	clonesizedistr.clear();
	joint_distr.clear();
	joint_distr_norm.clear();
	joint_distr_intm.clear();
	joint_distr_intm_norm.clear();
	CSD_intm_basal.clear();
	CSD_intm_basal_norm.clear();
	
#if reversible_diff == 1
//	r[0]=0.0;
//	parameters[1][0]=0.0;
//	maxrate=max_rate();
#endif
	
	
#if model == 1
	maxrate=max_rate();      // Attention !!! Model dependent
#else
	cerr "error: Model not yet fully defined.\nPlease implement choice of maxrate in run::setup()\n"; exit(1);
#endif
	
	cerr << "maxrate = " << maxrate << endl;
	
	cerr << "parameters:\n";
	cerr << "reprodrate = " << reprod_rate[0] << '\t' << "r = " << r[0] << endl;
	cerr << "stratrate = " << strat_rate[1] << '\t' << "shed_rate = " << shed_rate[2] << endl;
	cerr << "delta = " << delta[0] << endl;

	
	cellnumber_tot=0;
	Cnumber_glob.resize(celltypenumber);
	for (int i=0; i<celltypenumber; i++) 
	{
		Cnumber_glob[i]=0;
	};
    
    file_clonenumber_t.close();
    file_clonenumber_t.open("clonenumber_t.dat");
	
#if raw_data == 8
	set_cell_config(data_init);    
#else
	set_cell_config();
#endif
	
#if event_queue == 1
	// Create initial events
	for (int i=0; i<cells.size(); i++) 
	{
		if (i % 100 == 0) cerr << "i = " << i << endl;
		cell* curr_cell = cells[i];
		int eventtype = curr_cell->get_celltype();
		double rate = 0;
		
		if (eventtype == 0) 
		{
			rate = reprod_rate[0];
		}
		else if (eventtype == 1)
		{
			rate = strat_rate[1];
		}
		else if (eventtype == 2)
		{
			rate = shed_rate[2];
		};
		
#if sync_reprod == 1
		double eventtime = drand48()*1.0/rate;
#else
		double eventtime = log(1.0/drand48())/rate;
#endif
		curr_cell->create_new_event(eventtype,eventtime);
//		cerr << "initial event added\n";                   // Debug !!!
	};
	
	for (list<event>::iterator itr=eventlist.begin(); itr!=eventlist.end(); ++itr) 
	{
		cerr << "event time = " << itr->get_time() << endl;
	}
	
#endif
	
	cerr << "reprod_rate[0] = " << reprod_rate[0] << endl << "r[0] = " << r[0] << endl << "delta[0] = " << delta[0] << endl << "strat_rate[0] = " << strat_rate[0] << endl << "shed_rate[0] = " << shed_rate[0] << endl << endl; 
	
	vector<double> testparams;
	testparams.push_back(reprod_rate[0]);
	testparams.push_back(strat_rate[1]);
	
	cerr << "test prior function: " << reprod_rate[0] << '\t' << strat_rate[1] << '\t' << prior(reprod_rate[0],r[0],strat_rate[1],shed_rate[2],delta[0]) << endl;
	
//	for (int i=0; i<cells.size(); i++) 
//	{
//		if (cells[i]->celltype == 0) 
//		{
//			
//		}
//	}
	
}

//-----------------------------------------------------------------------
// Setup initial cell configurations
//-----------------------------------------------------------------------

void run::set_cell_config()
{
	if (clones.size()==0) 
	{
		cerr << "error: no clones\n";
		exit(1);
	}
	for (int i=0; i<clones.size(); i++) 
	{
		
		Clone* curr_clone=&clones[i];
		for (int j=0; j<init_clone.size(); j++) 
		{
			for (int k=0; k<init_clone[j].size(); k++) 
			{
				if (drand48() < init_clone[j][k]) 
				{
					cell* cell0 = new cell(curr_clone,this,k);
					cells.push_back(cell0);
					cell0->set_index(cells.size()-1);
					Cnumber_glob[k]++;
				};
			};
		};	
		if (curr_clone->get_cells()->size()==0) 
		{
			cell* cell0 = new cell(curr_clone,this,0);
			cells.push_back(cell0);
			cell0->set_index(cells.size()-1);
			Cnumber_glob[0]++;
		}
	};	
	
	for (int i=0; i<clones.size() ; i++)                            // Debug !!!
	{
		list<cell*>* cell_list = clones[i].get_cells();
		if ((*(cell_list->begin()))->get_celltype() != 0)
		{
	//		cerr << "error: celltype at beginning not 0\n";
	//		exit(1);
		};
		if (cell_list->size() != 1) 
		{
			if (init_clone_number < 2)
			{
	//			cerr << "error: more than one cell in clone at beginning\n";
				exit(1);
			};
		};
	};
};

void run::set_cell_config(vector<vector<int> > joint_CSD)       // Attention !!! only for 3 celltypes    // not finished !!!
{
//	cerr << "rho at run initialisation = " << rho << endl;    // Debug !!!
	clones.clear();
	int joint_CSD_norm=0;
	for (int i=0; i < joint_CSD.size(); i++) 
	{
		for (int j=0; j < joint_CSD[i].size(); j++) 
		{
			joint_CSD_norm = joint_CSD_norm + joint_CSD[i][j];
		};
	};
	for (int i=0; i < joint_CSD.size(); i++) 
	{
		for (int j=0; j < joint_CSD[i].size(); j++) 
		{
			joint_CSD[i][j] = joint_CSD[i][j]*(double)clonenumber/joint_CSD_norm;
		};
	};
	
	for (int i=0; i<joint_CSD.size(); i++)
	{
		for (int j=0; j<joint_CSD[i].size(); j++) 
		{
			for (int l=0; l<joint_CSD[i][j]; l++) 
			{
				
				Clone clone0;
				clones.push_back(clone0);
				Clone* curr_clone = &clones.back();
				for (int k=0; k<i; k++) 
				{
					if (drand48()<rho) 
					{
						cell* cell0 = new cell(curr_clone,this,0);	
						cells.push_back(cell0);
						cell0->set_index(cells.size()-1);
						Cnumber_glob[0]++;
//						cerr << "cell added!\n";    // Debug !!!
					}
					else 
					{
						cell* cell0 = new cell(curr_clone,this,1);
						cells.push_back(cell0);
						cell0->set_index(cells.size()-1);
						Cnumber_glob[1]++;
//						cerr << "cell added!\n";     // Debug !!!
					};
				};
				
				for (int k=0; k<j; k++) 
				{
					cell* cell0 = new cell(curr_clone,this,2); 		
					cells.push_back(cell0);
					cell0->set_index(cells.size()-1);
					Cnumber_glob[1]++;
//					cerr << "cell added!\n";     // Debug !!!
				};
			}
		};
	};
	
	cerr << "cell number in setup = " << cells.size() << endl;
	
};

void run::set_cell_config(const vector<vector<vector<int> > >& joint_CSD)       // Attention !!! only for 3 celltypes    // not finished !!!
{
	clones.clear();
	for (int i=0; i<joint_CSD.size(); i++)
	{
		for (int j=0; j<joint_CSD[i].size(); j++) 
		{
			for (int k=0; k<joint_CSD[i][j].size(); k++) 
			{
				
				Clone* clone0 = new Clone();
				for (int l=0; l<i; l++) 
				{
					cell* cell0 = new cell(clone0,this,0); 
					clone0->addcell(cell0);

				};
				for (int l=0; l<j; l++) 
				{
					cell* cell0 = new cell(clone0,this,1); 
					clone0->addcell(cell0);
				};
				for (int l=0; l<j; l++) 
				{
					cell* cell0 = new cell(clone0,this,2); 
					clone0->addcell(cell0);
				};
					clones.push_back(*clone0);
				delete clone0;
			}
		};
	};
	
};

//----------------- end setup cell configurations ----------------------------

//-----------------------------------------------------------------------
// Start Run
//-----------------------------------------------------------------------

int run::go()
{
	setup();
	cerr << "go\n";
	cerr << "time sim = " << time_sim << endl;             // Debug !!!
	while (time_sim < runtime)// && extinct_flag==0) 
	{
		iterate();
	};
	
	cerr << "log data at time " << time_sim << endl;
	if (logflag == 0) 
	{
		log_CSD();
        cerr << "logged\n";
	};
	
	
	
	int sum_model=0;
	double sum_exp_basal=0.0;
	double sum_exp_supra=0.0;
	double sum_norm=0.0;
	for (int i=0; i<clones.size(); i++) 
	{
		if (clones[i].get_cells()->size() > 0) sum_model = sum_model + clones[i].get_cells()->size();
	};
	for (int i=0; i<data0[0].size() ; i++) 
	{
		for (int j=0; j<data0[0][i].size(); j++) 
		{
			sum_norm = sum_norm + data0[0][i][j]; 
			sum_exp_basal = sum_exp_basal + i*data0[0][i][j];
		};
	};
	for (int i=0; i<data0[0].size(); i++) 
	{
		for (int j=0; j<data0[0][i].size(); j++) 
		{
			sum_exp_supra = sum_exp_supra + j*data0[0][i][j];
		};
	};
	
	
	cerr << "average clone size = " << (double)sum_model/clones.size() << endl;
	cerr << "average basal clones size data = " << (double)sum_exp_basal/sum_norm << endl;
	cerr << "average suprabasal clones size data  = " << (double)sum_exp_supra/sum_norm << endl;
	cerr << "suprabasal/basal ratio data = " << ((double)sum_exp_supra/sum_norm)/((double)sum_exp_basal/sum_norm) << endl;
	
	return 0;
	
};

//-----------------------------------------------------------------------
// Variate Parameters
//-----------------------------------------------------------------------

void run::variate_params(vector<int> var_vector0)                     // Attention !!! Model-dependent
{
#if raw_data == 2 || raw_data == 3 || raw_data == 4
	reprod_rate0_soraf =  min(reprodratelimit_high,max(reprodratelimit_low,reprod_rate0_soraf + stepsize_reprod*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[0]));
	r0_soraf = min(0.5,max(0.0,r0_soraf + stepsize_r*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[1]));
//	delta[0] = min(1.0, max(0.0,delta[0] + stepsize_delta*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[2]));
	strat_rate1_soraf = min(stratratelimit_high,max(stratratelimit_low,strat_rate1_soraf + stepsize_strat*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[3]));
	shedrate2_soraf = min(shedratelimit_high,max(shedratelimit_low,shedrate2_soraf + stepsize_shed*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[4]));
#elif raw_data == 0 && init_cells == 2
	init_clone[1][0] = min(1.0,max(0.0,init_clone[1][0] + stepsize_progfrac*T*log(1.0/drand48())*sgn(mrand48())));
	init_clone[1][1]=1-init_clone[1][0];
	init_clone[0][0] = min(1.0,max(0.0,init_clone[0][0] + stepsize_progfrac*T*log(1.0/drand48())*sgn(mrand48())));
	init_clone[0][1]=1-init_clone[0][0];
#else
	reprod_rate[0] = min(reprodratelimit_high,max(reprodratelimit_low,reprod_rate[0] + stepsize_reprod*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[0]));
	r[0] = min(0.5,max(0.0,r[0] + stepsize_r*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[1]));
//	delta[0] = min(1.0, max(0.0,delta[0] + stepsize_delta*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[2]));
//	strat_rate[1] = min(stratratelimit_high,max(stratratelimit_low,strat_rate[1] + stepsize_strat*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[3]));
	shed_rate[2] = min(shedratelimit_high,max(shedratelimit_low,shed_rate[2] + stepsize_shed*T*log(1.0/drand48())*sgn(mrand48())*var_vector0[4]));
#endif
//	direct_strat_frac[0] = direct_strat_frac[0] + stepsize_directstrat*log(1/drand48());
//	delta[0] = min(1.0,delta[0] + stepsize_delta*T*log(1.0/drand48())*sgn(mrand48()));
	
#if raw_data == 8
	rho = min(1.0,max(0.0,rho + stepsize_r*T*log(1.0/drand48())*sgn(mrand48())));
#endif
	maxrate = max_rate();           // Attention !!! Model-dependent
};

void run::switch_params()           // Attention !!! Model-dependent
{
#if raw_data != 22							// Attention !!! Maybe other assay runs will also be switched
	reprod_rate[0] = reprod_rate0_soraf;
	r[0] = r0_soraf;
	delta[0] = delta0_soraf;
	shed_rate[2]=shedrate2_soraf;
	strat_rate[1]=strat_rate1_soraf;
#else
	reprod_rate[0] = reprod_rate0_hom;
	r[0] = r0_hom;
	delta[0] = delta0_hom;
	shed_rate[2]=shedrate2_hom;
	strat_rate[1]=strat_rate1_hom;
#endif
	
//	reprod_rate[0] = reprod_rate0_MAML;
//	r[0] = r0_MAML;
//	delta[0] = delta0_MAML;
//	shed_rate[2]=shedrate2_MAML;
//	strat_rate[1]=strat_rate1_MAML;
	
	maxrate = max_rate();                // Attention !!! Model-dependent
	cerr << "maxrate in run() = " << maxrate << endl;
};

void run::switch_single_param(int param_index)           // Attention !!! Model-dependent
{
#if raw_data != 22
	if (param_index==0) 
	{
	reprod_rate[0] = reprod_rate0_hom;
	}
	else if (param_index==1)
	{
	r[0] = r0_hom;
	}
	else if (param_index==2)
	{
	delta[0] = delta0_hom;
	}
	else if (param_index==3)
	{
	shed_rate[2]=shedrate2_hom;
	}
	else if (param_index==4)
	{
	strat_rate[1]=strat_rate1_hom;
	}
	else 
	{
		cerr << "error: invalid parameter index\n"; exit(1);
	};
#else
	if (param_index==0) 
	{
		reprod_rate[0] = reprod_rate0_soraf;
	}
	else if (param_index==1)
	{
		r[0] = r0_soraf;
	}
	else if (param_index==2)
	{
		delta[0] = delta0_soraf;
	}
	else if (param_index==3)
	{
		shed_rate[2]=shedrate2_soraf;
	}
	else if (param_index==4)
	{
		strat_rate[1]=strat_rate1_soraf;
	}
	else 
	{
		cerr << "error: invalid parameter index\n"; exit(1);
	};
#endif

	
	maxrate = max_rate();                // Attention !!! Model-dependent
	cerr << "maxrate = " << maxrate << endl;
};

//-----------------------------------------------------------------------
// Log Clone Size Distributions
//-----------------------------------------------------------------------

int run::log_CSD()
{
	logflag=1;

//	cout << "program manually finished\n";
	
	cerr << "test: total cell number at begin of end of run = " << cells.size() << endl;
	
	clonesizedistr.clear();
	clonesizedistr_norm.clear();
	clonesizedistr_basal.clear();
	clonesizedistr_supra.clear();
	clonesizedistr_basal_norm.clear();
	cumuldistr_basal.clear();
	cumuldistr_total.clear();
	
	
	clonesizedistr.reserve(5000);
	clonesizedistr_norm.reserve(5000);
	clonesizedistr_basal.reserve(5000);
	clonesizedistr_supra.reserve(5000);
	clonesizedistr_basal_norm.reserve(5000);
	cumuldistr_basal.reserve(5000);
	cumuldistr_total.reserve(5000);
	
	joint_distr.clear();
	joint_distr_norm.clear();
	
	int maxsize=0;
	
	vector<vector<int> > clone_Cnumbers;
	clone_Cnumbers.clear();
	
	cerr << "test = " << *((&clones[5])->get_cells()->begin()) << endl;
	
	cerr << "clone_Cnumbers last added: ";
	for (int i=0; i<clones.size(); i++) 
	{
		if (clones[i].get_cells()->size() > 50000)
		{
			cerr << "clonesize too large\n"; 
			clonetobig_flag=1;
			return 1;
	//		time_sim = 10000;
		};
		
		vector<int> Cnumbers(celltypenumber,0);
		
		Clone* curr_clone = &clones[i];
		if (curr_clone->get_cells()->size()>0) 
		{
			for (list<cell*>::iterator itr = curr_clone->get_cells()->begin(); itr !=curr_clone->get_cells()->end(); itr++) 
			{
				Cnumbers[(*itr)->get_celltype()]++;
			};
		};
		
		clone_Cnumbers.push_back(Cnumbers);
//		cerr << "clone_Cnumbers last added = " << clone_Cnumbers.back()[1] << endl;
	}
	cerr << endl;                   // Debug !!!
	
#if raw_data == 9
	if (logctr==0)
	{
		limit_basal=min(limit_basal_1,(int)data0[logctr].size());
		limit_supra=min(limit_supra_1,(int)data0[logctr][0].size());
	}
	else if (logctr==1)
	{
		limit_basal=min(limit_basal_2,(int)data0[logctr].size());
		limit_supra=min(limit_supra_2,(int)data0[logctr][0].size());
	}
	else 
	{
		cerr << "error: invalid logctr\n";
		exit(1);
	};

#else
	limit_basal=min(limit_basal,(int)data0[logctr].size());
	limit_supra=min(limit_supra,(int)data0[logctr][0].size());
#endif
	
	validclonecount=0;
	for (int i = 0; i < clones.size(); i++) 
	{
	  if ((clone_Cnumbers[i][0]+clone_Cnumbers[i][1] >= start_valid_b || (clone_Cnumbers[i][0]+clone_Cnumbers[i][1] == start_valid_b - 1 && clone_Cnumbers[i][2] >= start_valid_s))/* && clone_Cnumbers[i][0]+clone_Cnumbers[i][1] < limit_basal && clone_Cnumbers[i][2] < limit_supra */)  
		{
			validclonecount++;
		};
	};
		
	cerr << "validclonecount = " << validclonecount << endl;
	
	for (int i=0; i<clones.size(); i++) 
	{
		int clonesize = clones[i].get_size();
		if (clonesize > maxsize) 
		{
			maxsize = clonesize;
			clonesizedistr.resize(maxsize+1,0);
		};
#if partial_basal_count == 1
		if (clone_Cnumbers[i][0]+clone_Cnumbers[i][1] >= start_valid_b || (clone_Cnumbers[i][0]+clone_Cnumbers[i][1] == start_valid_b - 1 && clone_Cnumbers[i][2] >= start_valid_s)) // Attention !!! Only if basal clones are counted only partially
		{
#endif
			clonesizedistr[clonesize]++;
#if partial_basal_count == 1
		};
#endif
		
	};
	maxsize=0;
	for (int i=0; i<clones.size(); i++)                                // Attention !!! Model-dependent , celltypenumber-dependent
	{
		int basalclone_size = clone_Cnumbers[i][0]+clone_Cnumbers[i][1];
		if (basalclone_size > maxsize) 
		{
			maxsize = basalclone_size;
			clonesizedistr_basal.resize(maxsize+1,0);
		};
#if partial_basal_count == 1
		if (basalclone_size >= start_valid_b || (basalclone_size == start_valid_b - 1 && clone_Cnumbers[i][2] >= start_valid_s)) // Attention !!! Only if basal clones are counted only partially
		{
#endif
			clonesizedistr_basal[basalclone_size]++;
#if partial_basal_count == 1
		};
#endif
	};
	maxsize=0;
	for (int i=0; i<clone_Cnumbers.size(); i++) 
	{
		int supraclone_size = clone_Cnumbers[i][2];
		if (supraclone_size > maxsize) 
		{
			maxsize = supraclone_size;
			clonesizedistr_supra.resize(maxsize+1,0);
		};
		clonesizedistr_supra[supraclone_size]++;
	};
	
	
	int maxsize_basal=data0[logctr].size();               
	int maxsize_supra=data0[logctr][0].size();
	
	joint_distr.reserve(600);
	joint_distr.resize(maxsize_basal+1);
	for (int i=0; i < joint_distr.size(); i++) 
	{
		joint_distr[i].clear();
		joint_distr[i].reserve(600);
		joint_distr[i].resize(maxsize_supra+1,0);
	}
	
	if (clone_Cnumbers.size() != clones.size())                // Debug !!!
	{
		cerr << "error: clones_Cnumbers not well-sized\n";
		exit(1);
	};
	for (int i=0; i<clones.size(); i++) 
	{
		int basalclone_size = clone_Cnumbers[i][0]+clone_Cnumbers[i][1];
		int supraclone_size = clone_Cnumbers[i][2];
		if (basalclone_size > maxsize_basal) 
		{
			maxsize_basal = basalclone_size;
			joint_distr.resize(maxsize_basal+1);
			for (int j=0; j<joint_distr.size(); j++) 
			{
				joint_distr[j].resize(maxsize_supra+1,0);
			};
		};
		if (supraclone_size > maxsize_supra) 
		{
			maxsize_supra = supraclone_size;
			for (int j=0; j<joint_distr.size(); j++) 
			{
				joint_distr[j].resize(maxsize_supra+1,0);
			};
		};
		
		//	cerr << "clone sizes: " << basalclone_size << '\t' << supraclone_size << endl;
		//	cerr << "dim. of joint distr: " << joint_distr.size() << '\t' << joint_distr[1].size() << endl; 
#if partial_basal_count == 1
		if (basalclone_size >= start_valid_b || (basalclone_size == start_valid_b - 1 && clone_Cnumbers[i][2] >= start_valid_s))
		{
#endif
		joint_distr[basalclone_size][supraclone_size]++;
#if partial_basal_count == 1
		};
#endif
		
	};
	
	clonesizedistr_norm = clonesizedistr;
	clonesizedistr_basal_norm.resize(clonesizedistr_basal.size());
	cellnumber_tot = cells.size();
	double clone_counter=0;
	cumuldistr_total.resize(clonesizedistr.size());
	for (int i=0; i < clonesizedistr.size(); i++) 
	{
		clonesizedistr_norm[i] = (double)clonesizedistr[i]/validclonecount;
		clonesizedistr[i] = (double)clonesizedistr[i]/validclonecount*clone_norm;
		clonesizedistr_supra[i] = (double)clonesizedistr_supra[i]/validclonecount*clone_norm;
		
		if (i >= start_valid_b - 1)
		{
			clone_counter = clone_counter + clonesizedistr[i];
			cumuldistr_total[i] = clone_norm - clone_counter;
		}
		else cumuldistr_total[i] = clone_norm;
	};
	clone_counter=0;
	cumuldistr_basal.resize(clonesizedistr_basal.size());
	for (int i=0; i < clonesizedistr_basal.size(); i++) 
	{
		clonesizedistr_basal_norm[i] = (double)clonesizedistr_basal[i]/validclonecount;
		clonesizedistr_basal[i] = clonesizedistr_basal_norm[i]*clone_norm;// *clonenumber_data[t];  // Attention !!! Depends on data set to be compared with
		
		if (i >= start_valid_b-1)
		{
			clone_counter = clone_counter + clonesizedistr_basal[i];
			cumuldistr_basal[i] = clone_norm - clone_counter;
		}
		else cumuldistr_basal[i] = clone_norm;
	};
	

	double norm_joint=0;
	for (int i=start_valid_b-1; i<clonesizedistr_basal_norm.size(); i++) 
	{
		norm_joint = norm_joint + clonesizedistr_basal_norm[i];
	};
	
	
	cerr << "norm of joint_distr = " << norm_joint << endl;
	if (norm_joint>1.01 || norm_joint<0.99)
	{	
		cerr << "error: CSD_basal_norm not normalized (in logCSD); norm = " << norm_joint << endl;
	//	exit(1);
	};
	
	
	joint_distr_norm.clear();
//	cerr << "joint_distr.size() = " << joint_distr.size() << endl;
	joint_distr_norm.reserve(600);
	joint_distr_norm.resize(joint_distr.size());
//	cerr << "joint_distr_norm.size() = " << '\t' << joint_distr_norm.size() << endl << "joint_distr.size() = " << '\t' << joint_distr.size() << endl;        // Debug !!!
	if (logctr>1) {cerr << "error: logctr to large, logctr = " << logctr << endl; exit(1);}          // Debug !!!
	for (int i=0; i < joint_distr.size(); i++) 
	{
		joint_distr_norm[i].reserve(600);
		joint_distr_norm[i].resize(joint_distr[i].size());
		for (int j=0; j < joint_distr[i].size(); j++) 
		{
//			cerr << "joint_distr_norm["<< i<< "].size() = " << endl << joint_distr_norm[i].size() << endl;        // Debug !!!
			joint_distr_norm[i][j] = (double)joint_distr[i][j]/validclonecount;
			joint_distr[i][j] = (double)joint_distr[i][j]/validclonecount*clone_norm;//clonenumber_data[logctr];
		};
	};
	
	for (int j=1; j<joint_distr[0].size(); j++) 
	{
		floating_clone_frac = floating_clone_frac + joint_distr[0][j];
	};
	cerr << "floating_clone_frac = " << floating_clone_frac << endl;
	
	Cnumber_glob[0] = 0;
	Cnumber_glob[1] = 0;
	Cnumber_glob[2] = 0;
	for (int i=0; i<clone_Cnumbers.size(); i++) 
	{
		Cnumber_glob[0] = Cnumber_glob[0] + clone_Cnumbers[i][0];
		Cnumber_glob[1] = Cnumber_glob[1] + clone_Cnumbers[i][1];
		Cnumber_glob[2] = Cnumber_glob[2] + clone_Cnumbers[i][2];
	};
	
	rho = (double)Cnumber_glob[0]/(Cnumber_glob[0] + Cnumber_glob[1]);
	
	cerr << "CSD.size() before = " << clonesizedistr_basal_norm.size() << endl;         // Debug !!!
	
	cerr << "CSD basal model after log:" << endl;                 // Debug !!!
	
	for (int i=0; i<clonesizedistr_basal.size(); i++)
	{
		cerr << clonesizedistr_basal[i] << endl;
	};
	
	cerr << "CSD basal model after log:" << endl;                 // Debug !!!
	
	for (int i=0; i<clonesizedistr_basal.size(); i++)
	{
		cerr << clonesizedistr_basal_norm[i] << endl;
	};
	
	CSD_intm_basal.push_back(clonesizedistr_basal);
	CSD_intm_basal_norm.push_back(clonesizedistr_basal_norm);
	CSD_intm_total.push_back(clonesizedistr);
	CSD_intm_total_norm.push_back(clonesizedistr_norm);
	cerr << "CSD.size() after = " << CSD_intm_basal_norm.back().size() << endl;         // Debug !!!
	joint_distr_intm.push_back(joint_distr);
	joint_distr_intm_norm.push_back(joint_distr_norm);
	
	cerr << "rho = " << rho << endl; 
	file_av << time_sim << '\t' <<  av_clonesize_basal() << endl;
	
	logctr++;
	cerr << "log data, logctr = " << logctr << endl;
	
	return 0;

};


//-----------------------------------------------------------------------
// Evaluate Goodness of parameter choice
//-----------------------------------------------------------------------




double run::loglikely(const vector<vector<vector<double> > >& data0)       // Joint data distr with basal and suprabasal celltypes + time points
{

	
	vector<vector<double> > joint_cumul_distr(joint_distr_norm.size());
	
	double norm_joint=0;
	for (int i=0; i<joint_distr_norm.size(); i++) 
	{
		joint_cumul_distr[i].resize(joint_distr_norm[i].size());
		for (int j=start_val(i); j<joint_distr_norm[i].size(); j++) 
		{
			norm_joint = norm_joint + joint_distr_norm[i][j];
			joint_cumul_distr[i][j] = norm_joint;
		};
	};
	
	
	cerr << "norm of joint_distr = " << norm_joint << endl;

	if (norm_joint>1.01 || norm_joint<0.99)
	{	
		cerr << "error: joint_distr_norm not normalized (in loglikely()); norm = " << norm_joint << endl;
		exit(1);
	};
	
	if (timepoints != joint_distr_intm_norm.size())                         // Debug !!!
	{
		cerr << "error: number of timepoints doesn't fit\n";
		cerr << "CSD_intm_basal_norm.size() = " <<  CSD_intm_basal_norm.size() << endl;
		cerr << "data0.size() = " <<  data0.size() << endl;
		cerr << "logctr = " << logctr << endl;
		cerr << "time = " << time_sim << endl;
		cerr << "cell number = " << cells.size() << endl;
//		exit(1);
	};
	
	double sum_t=0;
	for (int t=0; t < data0.size(); t++) 
	{
		int limit_basal0=0;
		int limit_supra0=0;
		
		norm_joint=0;
#if raw_data == 9 || raw_data == 24
		if (t==0)
		{
			limit_basal0=min(limit_basal_1,(int)data0[t].size());
			limit_supra0=min(limit_supra_1,(int)data0[t][0].size());
		}
		else if (t==1)
		{
			limit_basal0=min(limit_basal_2,(int)data0[t].size());
			limit_supra0=min(limit_supra_2,(int)data0[t][0].size());
		};
#else
		limit_basal0=min(limit_basal,(int)data0[t].size());
		limit_supra0=min(limit_supra,(int)data0[t][0].size());
#endif

		int data_cumul_ctr=0;
		for (int i=0; i<limit_basal0; i++) 
		{
			for (int j=start_val(i); j<limit_supra0; j++) 
			{
				norm_joint = norm_joint + joint_distr_intm_norm[t][i][j];
				data_cumul_ctr = data_cumul_ctr + data0[t][i][j];
			};
		};
//		if (norm_joint>1.001 || norm_joint<0.999)					// Debug !!!
//		{	
//			cerr << "error: not normalized; norm = " << norm_joint << endl;
//			exit(1);
//		};
      
		for (int i=0; i<limit_basal0; i++) 
		{
			for (int j=start_val(i); j<limit_supra0; j++) 
			{
				if (joint_distr_intm_norm[t][i][j]>0) 
				{
				  sum_t = sum_t + clonenumber_exp[t]*(double)data0[t][i][j]/clone_norm*log((double)joint_distr_intm_norm[t][i][j]/norm_joint);
				}
				else if (data0[t][i][j]>0)
				{
				//	sum_t = sum_t + data0[t][i][j]*log(0.3/validclonecount);
					sum_t = sum_t - 100000;
				};

			};
		};
		if (norm_joint < 1)
		{
//			sum_t = sum_t + (clonenumber_data[t]-data_cumul_ctr) * log(1 - norm_joint);
		};
		
	};
	
	return sum_t+log(prior(reprod_rate[0],r[0],strat_rate[1],shed_rate[2],delta[0]));
};

double run::loglikely_basal(vector<vector<vector<double> > > data0)       // Joint data distr with basal and suprabasal celltypes + time points
{
	vector<vector<double> > CSD_intm_basal_norm0(data0.size());
	
	for (int t=0; t < data0.size(); t++) 
	{
		data0[t]=bin(data0[t],binmode,binsize);
		CSD_intm_basal_norm0[t]=bin(CSD_intm_basal_norm[t],binmode,binsize);
		if (CSD_intm_basal_norm0[t].size() < data0[t].size()) CSD_intm_basal_norm0[t].resize(data0[t].size(),0);
		double norm_joint=0;
		for (int i=(start_valid_b-1)/binsize; i<CSD_intm_basal_norm0[t].size(); i++) 
		{
			norm_joint = norm_joint + CSD_intm_basal_norm0[t][i];
		};
		cerr << "norm of joint_distr = " << norm_joint << endl;
		if (norm_joint>1.01 || norm_joint<0.99)
		{	
			cerr << "error: CSD_basal_norm not normalized; norm = " << norm_joint << endl;
			exit(1);
		};
		
	};
	
	if (timepoints != CSD_intm_basal_norm0.size())                         // Debug !!!
	{
		cerr << "error: number of timepoints doesn't fit\n";
		cerr << "CSD_intm_basal_norm.size() = " <<  CSD_intm_basal_norm0.size() << endl;
		cerr << "data0.size() = " <<  data0.size() << endl;
		cerr << "logctr = " << logctr << endl;
		cerr << "time = " << time_sim << endl;
		cerr << "cell number = " << cells.size() << endl;
		exit(1);
	};
	
	double sum_t=0;
	for (int t=0; t < data0.size(); t++) 
	{
		int limit_basal0=0;
//		int limit_supra0=0;
		
		double norm_joint=0;
#if raw_data == 9 || raw_data == 24
		if (t==0)
		{
			limit_basal0=min(limit_basal_1,(int)data0[t].size())/binsize;
		}
		else if (t==1)
		{
			limit_basal0=min(limit_basal_2,(int)data0[t].size())/binsize;
		};
#else
		limit_basal0=min(limit_basal,(int)data0[t].size())/binsize;
//		limit_supra0=min(limit_supra,(int)data0[t][0].size())/binsize;
#endif
		int data_cumul_ctr=0;
		cerr << "fitted basal CSD (model + data)\n";
		for (int i=(start_valid_b-1)/binsize; i<limit_basal0; i++) 
		{
			data_cumul_ctr = data_cumul_ctr + data0[t][i][1];
			norm_joint = norm_joint + CSD_intm_basal_norm0[t][i];
			
		};
		cerr << "limit_basal0 = " << limit_basal0 << endl << "norm_joint = " << norm_joint << endl;
		norm_joint = norm_joint*clone_norm;
		if (norm_joint == 0) {cerr << "error: norm_joint=0" << endl; exit(1);}
		for (int i=(start_valid_b-1)/binsize; i<limit_basal0; i++) 
		{
			cerr << CSD_intm_basal_norm0[t][i]/norm_joint*data_cumul_ctr*data_cumul_ctr << '\t' << data0[t][i][1] << endl;
			//norm_joint = norm_joint + CSD_intm_basal_norm0[t][i];
			if (CSD_intm_basal_norm0[t][i] > 0) 
			{
				sum_t = sum_t + data0[t][i][1]*log((double)CSD_intm_basal_norm0[t][i]/norm_joint*data_cumul_ctr);
			}
			else if (data0[t][i][1]>0)
			{
				//	sum_t = sum_t + data0[t][i][1]*log(0.3/validclonecount);
				sum_t = sum_t - 100000;
			};
			
		};
		if (norm_joint < 1) 
		{
			//		sum_t = sum_t + (clonenumber_data[t]-data_cumul_ctr) * log(1 - norm_joint);
		};
	};
	cerr << "limit_basal0 = " << limit_basal0 << endl;
	
#if raw_data == 10 || raw_data == 18 || raw_data == 1 || raw_data == 30 || raw_data == 24
	
	sum_t = sum_t+log(prior(reprod_rate[0],r[0],strat_rate[1],shed_rate[2],delta[0]));
	
#endif
	
	return sum_t;
};


double run::loglikely_total(vector<vector<vector<double> > > data0)       // Joint data distr with total and supratotal celltypes + time points
{
	vector<vector<double> > CSD_intm_total_norm0(data0.size());
	
	for (int t=0; t < timepoints; t++) 
	{
		data0[t]=bin(data0[t],binmode,binsize);
		cerr << "CSD.size() in likely_total = " << CSD_intm_total_norm[t].size() << endl;         // Debug !!!
		CSD_intm_total_norm0[t]=bin(CSD_intm_total_norm[t],binmode,binsize);
		if (CSD_intm_total_norm0[t].size() < data0[t].size()) CSD_intm_total_norm0[t].resize(data0[t].size(),0);
		
		double norm_joint=0;
		for (int i=(start_valid_b-1)/binsize; i<CSD_intm_total_norm0[t].size(); i++) 
		{
			norm_joint = norm_joint + CSD_intm_total_norm0[t][i];
		};
		cerr << "norm of joint_distr = " << norm_joint << endl;
		if (norm_joint>1.01 || norm_joint<0.99)
		{	
			cerr << "error: CSD_total_norm not normalized; norm = " << norm_joint << endl;
			exit(1);
		};
		
	};
	
	if (timepoints != CSD_intm_total_norm0.size())                         // Debug !!!
	{
		cerr << "error: number of timepoints doesn't fit\n";
		cerr << "CSD_intm_total_norm.size() = " <<  CSD_intm_total_norm0.size() << endl;
		cerr << "data0.size() = " <<  data0.size() << endl;
		cerr << "logctr = " << logctr << endl;
		cerr << "time = " << time_sim << endl;
		cerr << "cell number = " << cells.size() << endl;
		exit(1);
	};
	
	double sum_t=0;
	for (int t=0; t < data0.size(); t++) 
	{
		int limit_total0 = min(limit_basal0+limit_supra0,(int)data0[t].size());
		
		double norm_joint=0;
		int data_cumul_ctr=0;
		for (int i=(start_valid_b-1)/binsize; i<limit_total0; i++) 
		{
			data_cumul_ctr = data_cumul_ctr + data0[t][i][1];
			norm_joint = norm_joint + CSD_intm_total_norm0[t][i];
			if (CSD_intm_total_norm0[t][i]<0) 
			{
				cerr << "CSD_intm_total_norm[t][i] = " << CSD_intm_total_norm0[t][i] << endl;
				exit(1);
			}
		};
		cerr << "limit_basal0 = " << limit_basal0 << endl << "norm_joint = " << norm_joint << endl;
		norm_joint = norm_joint*clone_norm;
		if (norm_joint == 0) {cerr << "error: norm_joint=0" << endl; exit(1);}
		for (int i=((start_valid_b-1)+start_valid_s)/binsize; i<limit_total0; i++) 
		{
			//norm_joint = norm_joint + CSD_intm_total_norm0[t][i];
			if (CSD_intm_total_norm0[t][i] > 0) 
			{
				sum_t = sum_t + data0[t][i][1]*log((double)CSD_intm_total_norm0[t][i]/norm_joint*data_cumul_ctr);
			}
			else if (data0[t][i][1]>0)
			{
				//	sum_t = sum_t + data0[t][i][1]*log(0.3/validclonecount);
				sum_t = sum_t - 100000;
			};
			
		};
		if (norm_joint < 1) 
		{
			//		sum_t = sum_t + (clonenumber_data[t]-data_cumul_ctr) * log(1 - norm_joint);
		};
	};
	
#if raw_data == 10 || raw_data == 18 || raw_data == 1 || raw_data == 30 || raw_data == 24 || raw_data == 25
	
	sum_t = sum_t+log(prior(reprod_rate[0],r[0],strat_rate[1],shed_rate[2],delta[0]));
	
#endif
	
	return sum_t;
};

//-------------------------------------------------------------
// Log Averages
//-------------------------------------------------------------


void run::print_averages(int runctr, int esmblctr)
{
	double av_ctr=0;
	basal_av=0;
	for (int i=1; i < clonesizedistr_basal.size(); i++) 
	{
		av_ctr = av_ctr + clonesizedistr_basal[i];
		basal_av = basal_av + i*clonesizedistr_basal[i];
	};
	basal_av = (double)basal_av/av_ctr;
	
	av_ctr=0;
	total_av=0;
	for (int i=1; i < clonesizedistr.size(); i++) 
	{
		av_ctr = av_ctr + clonesizedistr[i];
		total_av = total_av + i*clonesizedistr[i];
	};
	total_av = (double)total_av/av_ctr;
	cerr << "average basal clone size = " << basal_av << endl;
	file_av << esmblctr << '\t' << runctr << '\t' << basal_av << '\t' << total_av << '\t' << rho << endl;
};								

double run::av_clonesize()
{
	vector<vector<int> > clone_Cnumbers;
	clone_Cnumbers.clear();
	
	cerr << "test = " << *((&clones[5])->get_cells()->begin()) << endl;
	
	cerr << "clone_Cnumbers last added: ";
	for (int i=0; i<clones.size(); i++) 
	{
		if (clones[i].get_cells()->size() > 50000)
		{
			cerr << "clonesize too large\n"; 
			clonetobig_flag=1;
			return 1;
			//		time_sim = 10000;
		};
		
		vector<int> Cnumbers(celltypenumber,0);
		
		Clone* curr_clone = &clones[i];
		if (curr_clone->get_cells()->size()>0) 
		{
			for (list<cell*>::iterator itr = curr_clone->get_cells()->begin(); itr !=curr_clone->get_cells()->end(); itr++) 
			{
				Cnumbers[(*itr)->get_celltype()]++;
			};
		};
		
		clone_Cnumbers.push_back(Cnumbers);
		//		cerr << "clone_Cnumbers last added = " << clone_Cnumbers.back()[1] << endl;
	}
	
	
	double cellsum=0;
	int ctr=0;
	for (int i=0; i<clones.size(); i++)
    {
		int basalclone_size_curr = clones[i].get_size(); 
		if (basalclone_size_curr > 0) {cellsum = cellsum + basalclone_size_curr; ctr++;}
		
    };
	
	cellsum = (double)cellsum/ctr;
	
	return cellsum;
}

double run::av_clonesize_basal()
{
	vector<vector<int> > clone_Cnumbers;
	clone_Cnumbers.clear();

	for (int i=0; i<clones.size(); i++) 
	{
		vector<int> Cnumbers(celltypenumber,0);
		
		Clone* curr_clone = &clones[i];
		if (curr_clone->get_cells()->size()>=start_valid_b-1) 
		{
			for (list<cell*>::iterator itr = curr_clone->get_cells()->begin(); itr !=curr_clone->get_cells()->end(); itr++) 
			{
				Cnumbers[(*itr)->get_celltype()]++;
			};
		};
		
		clone_Cnumbers.push_back(Cnumbers);
		//		cerr << "clone_Cnumbers last added = " << clone_Cnumbers.back()[1] << endl;
	};
	
	double cellsum=0;
	int ctr=0;
	for (int i=0; i<clone_Cnumbers.size(); i++) 
	{
		if (clone_Cnumbers[i][0] + clone_Cnumbers[i][1] > 0) 
		{
			cellsum = cellsum + clone_Cnumbers[i][0] + clone_Cnumbers[i][1];
			ctr++;
		};
    };
	
	cellsum = (double)cellsum/ctr;
	return cellsum;

}


int run::catch_global_events()
{
	if (time_sim > t_switch && flag_switch==0)
	{
		flag_switch = 1;
		switch_single_param(2);
		cerr << endl << "!!!!!!!!!!!!!!!!!! switched single parameter at time " << time_sim << " !!!!!!!!!!!!!!!!!!!\n\n";                  // Debug!
	};
	if (time_sim > t_switch && flag_switch_2 == 0) 
	{
		flag_switch_2 = 1;
		switch_params();
		
		cerr << endl << "!!!!!!!!!!!!!!!!!! switched other parameters at time " << time_sim << " !!!!!!!!!!!!!!!!!!!\n\n";                  // Debug!
	};
	
#if raw_data == 9 || raw_data == 24
	if (time_sim > logtime1 && logctr==0)              // Data-dependent !!! depends on number of time points used simultanuously
	{
		cerr << "logtime1 = " << logtime1 << endl;     // Debug !!!
		log_CSD();
		cerr << "logged data at time " << time_sim << endl;
	};
#endif
#if raw_data == 20 || raw_data == 21
	if (log_av_flag == 0)
	{	
		//	    cerr << "log_av_flag == 0 !!!\n";	
		if (time_sim > timepoints_av[0] && logctr==0)              // Data-dependent !!! depends on number of time points used simultanuously
		{
			log_CSD();
		};
	}
	else if (log_av_flag == 1)
	{
		//  cerr << "log_av_flag == 1 !!!\n";
		if (time_sim > timepoints_av[0] && logctr==0)              // Data-dependent !!! depends on number of time points used simultanuously
		{
			file_av_clonesize_basal << time_sim << '\t' <<  av_clonesize_basal() << endl;
			log_CSD();
			
		};
		
		for (int j=1; j < timepoints_av.size(); j++)
		{
			if (time_sim > timepoints_av[j] && logctr==j)              // Data-dependent !!! depends on number of time points used simultanuously
			{
				file_av_clonesize_basal << time_sim << '\t' <<  av_clonesize_basal() << endl;
				logctr++;	
			};
		};
	};
#endif
	
	if (time_sim > runtime) 
	{
		cerr << "log data at time " << time_sim << endl;
		log_CSD();
		completed=1;
		
		return 0;
	};

	return 0;
};



//-----------------------------------------------------------------------
// Iteration Routine
//-----------------------------------------------------------------------

int run::iterate()
{
	cerr << "time = " << time_sim << endl << "cell number = " << cells.size() << endl << "maxrate = " << maxrate << endl; 
	
    file_clonenumber_t << time_sim << '\t' << clones.size() << endl;
    
	int N_curr = cells.size();

#if discrete_time == 1
	vector<cell*> cells_aux = cells;
#endif
	
	for (int i=0; i<N_curr; i++) 
	{
#if discrete_time == 1
        if (cells_aux[i]->get_celltype()==0)
        {
            cells_aux[i]->divide();
        }
		
//		time_sim = time_sim + 1.0/reprod_rate[0]/N_curr;
		
#elif event_queue == 1
	eventlist.begin()->execute();
		catch_global_events();
	
	if (completed==1) return 0;
#else

		//		cerr << "site of cell" << i << " = " << cells[i][0];            // Debug !!!
		//	cerr << "population of site 0 = " << lattice[0].get_population()[0] << endl;
		
		
		int k = (int)floor(((norm * mrand48() + 0.5) * cells.size()));         // Choose cell randomly
		
		while (k > cells.size()-1)                                          // Check if chosen index is not too high
		{
			k = (int)floor(((norm * mrand48() + 0.5) * cells.size()));
		};
		
		
		curr_cell = cells[k];
		
		if (cells.size() < 1) 
		{
			cerr << "error: no cells\n";
			exit(1);
		};
		
		///////////////// Update Cell ///////////////////////////////
		
		curr_cell->update();
		
		drand48();
		
//		cerr << "parameters: \n reprod rate = " << reprod_rate[0] << endl << "strat rate = " << strat_rate[1] << endl << "shed rate = " << shed_rate[2] << endl;
		
	//	cerr << "time before = " << time_sim << endl << "cell number = " << cells.size() << endl << "maxrate = " << maxrate << endl;                    // Debug !!!
		
		time_sim = time_sim + log(1/drand48())/(cells.size()*maxrate);
//		time_sim = time_sim - log(1-drand48())/(cells.size()*maxrate);
//		cerr << "time after = " << time_sim << endl;
		
		if (completed == 1) return 0;
		
		catch_global_events();
		

		
		
	
#endif

    };
    
#if discrete_time == 1
    for (int i=0; i<cells.size(); i++)
    {
        drand48();
        if (cells[i]->get_celltype() == 1 && drand48() < strat_rate[1]/(reprod_rate[0]+strat_rate[1]))
        {
            cells[i]->kill();
            i--;
        };
    };
    if (completed==1) return 0;
    time_sim = time_sim + 1.0/reprod_rate[0];
    catch_global_events();

#endif
    cerr << "time = " << time_sim << endl << endl;
	
	
#if Qt_flag==1	
	int test_sum = 0;
	if ((int)floor(time_sim)%10 == 0) 
	{
		for (int i=1; i<displays.size(); i++) 
		{
			displays[i]->display_val(Cnumber_glob[i]);
		};
		display_counter->display_val(time_sim);
		
		for (int i=0; i<clonesizedistr.size(); i++) 
		{
			clonesizedistr[i]=0;
		};
		
		maxsize=0;
		for (int i=0; i<clones.size(); i++) 
		{
			int clonesize = clones[i].get_size();
			//			if (clonesize>0) 
			//			{
			if (clonesize > maxsize) 
			{
				maxsize = clonesize;
				clonesizedistr.resize(maxsize+1,0);
			};
			clonesizedistr[clonesize]++;
			//			}
			
		};
		
	};
	
	
#endif	
	
	if (time_sim > evalctr*evalintv) 
	{
		evalctr++;
		file_av_clonesize_basal << time_sim << '\t' << av_clonesize_basal() << endl;
	}
	
	if (cells.size() > N_max || cells.size() < N_min) 
	{
		cerr << "cell number too large\n" << endl << endl;
		time_sim = 100000;
		completed=0;
		return 1;
	//	exit(2);
	}
    
    cerr << "all fine at time = ";
    cerr << time_sim << "\n\n";
	
	return 0;
};

//----------------------------------------------------------------------------
//         End Iteration Routine
//----------------------------------------------------------------------------

