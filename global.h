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


#ifndef _GLOBAL_
#define _GLOBAL_

#define parameter reprod_rate0
#define parameter_str "reprod_rate0"
#define par_type double
#define par_type_str "double"
#define raw_data 1
// 0 = control, 10d, 1 = Sorafenib pencil (10d), 2 = Sorafenib 22+10d, 3 = Sorafenib 4+10d 4 = Sorafenib 4+28d
					// 5 = MAML 3d, 6 = MAML 7d, 7 = MAML 10d, 8 = 4+22d with 4+10d initialization, 9 = MAML 7d+10d simulatanuously, 10 = Soraf+DEN
					// in lesions 10d, 11 = Soraf+DEN out of lesions 10d, 13 = EdU MAML 24h, 14 = EdU MAML 48h, 15 = EdU MAML 72h. 16 = YFP 3d (MAML control)
					// 17 = YFP 7d (MAML control), 18 = Soraf+DEN in lesions 20d, 19 = MAML 15d, 20 = Log averages Control 10+30+90+180+120d,
					// 21 = Log averages MAML 10+30+90+120d, 22 = 15d switch to normal, 23 = plain 3m, 24 = Soraf+DEN 10d+20d YFP, 25 = Soraf+DEN 10d confetti out of lesions,
					// 26 = Log averages normal 10+30+90+120d, 27 =  no data, free runtime, 28 = DEN+Ras 42d all, 29 = DEN+Ras 56d,
                    // 30 = Sorafenib pencil (21d)

#define reversible_diff 0     // Attention !!! r is set zero in run::setup() for reversible_diff == 1
#define Dcell_progenitors 0

#define plain 0

#define confetti 0

#define sim_mode 0     // 0 = single run, 1 = evolution, 2 = fixed parameters, get confidence interval

#define fitbasal 0   // 0 = fit joint CSD, 1 = fit basal CSD, 2 = fit basal CSD model to total CSD data, 3 = fit total CSD to  total CSD
#define nodata 0

#define model 1			// 1 = single prog. + DB + DSB cells, 2 = same as '1', but only basal, 3 = 2 types of progenitors.

#define many_step_process 0
#define init_cells 1

#define param_scan 0
#define partial_basal_count 1

#define get_bestparam 0

#define discrete_time 0
#define instant_strat 0

#define event_queue 0
#define sync_reprod 0
#define scan_rho 0

#define shedrate_fixed 0            // Attention !!! Ticking this may lead to wrong results !!!

#define fit_directstrat 0
#define fit_init2prog_frac 0


#define errorbars 1			// definition of errorbars: 0 = confidence interval, 1 = standard deviation, 2 = distinguished upper/lower std dev

#define dim 0         // System Dimension

#define param_fixed 0   // Fix a parameter at parameter variation, 0 = none, 1 = reprod, 2 = r, 3 = delta, 4 = strat, 5 = shed

#define Qt_flag 0     // 1, if Qt is used, else 0 

using namespace std;

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <strstream>
#include <algorithm>
#include <map>

#include "MainWindow.h"
#include "funcs.h"
#include "cell.h"
#include "clone.h"
#include "run.h"
#include "event.h"

#if Qt_flag==1            // include Qt library classes

#include <QtGui>
#include <QTimer>
#include <QSlider>
#include "LabeledSlider.h"
#include "LabeledDisplay.h"
#include "PaintWindow.h"

#endif


class Clone;
class cell;
class MainWindow;

extern MainWindow* window;

extern char fileoption;
extern int index_run;

extern int itr_intv;
extern int output_intv;
extern int image_intv;

//extern int logflag;     // Debug !!!

extern int timepoints;

extern double t_switch;
extern double runtime;

extern int start_valid_b;
extern int start_valid_s;


extern int runnumber;
extern int runnumber2;
extern int esmblnumber;
extern int run_index;
extern int extinct_flag;

extern double runtime0;

extern int excluded_clones;

extern cell* curr_cell;

extern int init_clone_number;
extern int basalcelltype_number;
extern double prog_cell_frac_1;
extern double prog_cell_frac_2;
extern double prog_cell_frac_3;

extern int image_ctr;
//extern double time_sim;

extern int celltypenumber;
extern int prolif_celltype_limit;
extern int shed_celltype_limit;
extern const int clonenumber_def;
extern int clonenumber_bestfit;
extern int clone_norm;
extern vector<int> clonenumber_exp;
//extern vector<int> clonenumber_exp_time;
extern vector<vector<vector<double> > > data_cumul;

extern const int param_number;
extern vector<int> var_vector;

extern int N_max;
extern int N_min;

extern int scatterpt_number;

extern double logtime1;
#if raw_data == 20 || raw_data == 21
extern double logtime2;
extern double logtime3;
extern double logtime4;
#endif

extern double siglevel;
extern double siglevel2;

extern double reprodrate_bayes_av;
extern double var_bayes;
extern double r_bayes_av;
extern double var_r_bayes;

extern double SCfrac;

extern map<vector<double>,double,bool(*)(vector<double>,vector<double>)> prior_map;

extern int validclonecount;

extern int limit_basal0;
extern int limit_supra0;

extern int limit_basal;
extern int limit_supra;

#if raw_data == 9 || raw_data == 24
extern int limit_basal_1;
extern int limit_supra_1;
extern int limit_basal_2;
extern int limit_supra_2;
#endif

extern vector<run*> run_old;


//extern vector<cell*> cells;           // cell[i][0]=clone of cell[i]; microbes[i][1] = type of microbe[i]
//extern vector<Clone> clones;

extern int reprod_step_number;
extern double SB_B_ratio;

extern vector<double> reprod_rate_init; 
extern vector<double> r_init;
extern vector<double> dediff_rate_init;
extern vector<double> delta_init;
extern vector<double> strat_rate_init;
extern vector<double> shed_rate_init;
extern vector<double> direct_strat_frac_init;
extern vector<double> switchrate_forw_init;
extern vector<double> switchrate_backw_init;
extern vector<double> init2prog_frac_init;

extern double delay;

extern double strat_rate_eff;

extern double reprod_rate0;
#if Dcell_progenitors == 1
extern double reprod_rate1;
#endif
extern double r0;
extern double delta0;
extern double strat_rate1;
extern double shedrate2;
extern double direct_strat_frac0;
extern double init2prog_frac0;
extern double mut_rate0;

extern double reprod_rate0_soraf;
extern double r0_soraf;
extern double delta0_soraf;
extern double strat_rate1_soraf;
extern double shedrate2_soraf;
extern double direct_strat_frac0_soraf;

extern double reprod_rate0_hom;
extern double r0_hom;
extern double delta0_hom;
extern double strat_rate1_hom;
extern double shedrate2_hom;
extern double direct_strat_frac0_hom;

extern double reprod_rate0_MAML;
extern double r0_MAML;
extern double delta0_MAML;
extern double strat_rate1_MAML;
extern double shedrate2_MAML;
extern double direct_strat_frac0_MAML;

extern double switchrate_forw0;
extern double switchrate_backw0;

extern double switch_reprodrate;
extern double switch_delta;
extern double switch_r;

extern double reprodratelimit_low;
extern double reprodratelimit_high;
extern double stratratelimit_low;
extern double stratratelimit_high;
extern double shedratelimit_low;
extern double shedratelimit_high;
extern double directstratlimit_low;
extern double directstratlimit_high;
extern double init2proglimit_low;
extern double init2proglimit_high;

extern double stratrate_intv;
extern double shedrate_intv;
extern double reprodrate_intv;
extern double init2prog_intv;

extern int strat_infinity_flag;
extern double rholimit_low;
extern double rholimit_high;
extern double rho_intv;
extern double directstrat_intv;

extern double stepsize_reprod;
extern double stepsize_r;
extern double stepsize_strat;
extern double stepsize_shed;
extern double stepsize_directstrat;
extern double stepsize_delta;
extern double stepsize_progfrac;

extern double T_max;  
extern double T;
extern double Tdropruns;
extern double Tdrop1;
extern double T0;

extern int evalintv;

extern int binsize;
extern double binsize_lesion;

#if raw_data == 20 || raw_data == 21 ||raw_data == 26
extern vector<double> timepoints_av;
#endif

extern int maxCS_av;

//extern double maxrate;

extern int cellnumber_tot;
extern vector<int> Cnumber_glob;

//extern vector<double> clonesizedistr;
//extern vector<double> clonesizedistr_norm;
//extern vector<double> clonesizedistr_basal;
//extern vector<double> clonesizedistr_supra;
//extern vector<vector<double> > joint_distr;

//extern int maxsize;

extern const double norm;              // Normalization factor for random numbers
extern int ibm;                               // random number seed. Attention: fixed in current version !!!
extern int seed;
extern double rand1;
extern double rand2;

extern vector<int> bins;
extern char binmode; 

extern vector<vector<double> > joint_distr_av;
extern vector<double> CSD_basal_av;
extern vector<double> CSD_total_av;
extern vector<vector<double> > lesion_SD;
extern vector<double> cumuldistr_basal_av;
extern vector<double> cumuldistr_total_av;
extern unsigned int i_curr;
extern unsigned int j_curr;

extern ofstream file_CSD;
extern ofstream file_CSD_norm;
extern ofstream file_CSD_basal;
extern ofstream file_CSD_supra;
extern ofstream file_CSD_basal_resc;
extern ofstream file_CSD_total_resc;
extern ofstream file_cumuldistr_total;
extern ofstream file_cumuldistr_basal;
extern ofstream file_cumuldistr_supra;
extern ofstream file_Jointdistr;
extern ofstream file_Jointdistr_norm;
extern ofstream file_Jointdistr_matrix;
extern ofstream file_Jointdistr_transp;
extern ofstream file_parameters;
extern ofstream file_loglikely;
extern ofstream file_av;
extern ofstream file_avsupra_basal;
extern ofstream file_celltypes_basal;
extern ofstream file_clonenumber_t;

extern ofstream file_data_basal_out;
extern ofstream file_data_total_out;
extern ofstream file_data_basal_resc_out;
extern ofstream file_data_total_resc_out;
extern ofstream file_data_avsupra_out;
extern ofstream file_data_cumul_basal_out;
extern ofstream file_data_cumul_total_out;
extern ofstream file_lesion_SD_out;
extern ofstream file_pval;

extern ofstream file_CSD_basal_intm;
extern ofstream file_Jointdistr_intm;
extern ofstream file_Jointdistr_matrix_intm;

extern ofstream file_H2B_divctr_scatter;
extern ofstream file_H2B_av_clonesize;

extern ofstream file_temp;         // Debug !!!

extern ifstream file_data;
extern ifstream file_data_total;
extern ifstream file_data_init;
extern ifstream file_data_basal;

#if raw_data == 9 || raw_data == 20 || raw_data == 21 || raw_data == 24
extern ifstream file_data2;
extern ifstream file_data_basal2;
#endif

#if raw_data == 24
extern ifstream file_data_total2;
#endif

#if raw_data == 10 || raw_data == 18
extern ifstream file_data_lesions;
#endif

extern vector<vector<int> > data_init;
extern vector<vector<vector<double> > > data0;
extern vector<vector<vector<double> > > data_basal;
extern vector<vector<vector<double> > > data_total;

extern int fileerrorbar_nr;
extern vector<ofstream*> file_errorbars;
extern ofstream file_basal_errorbars;
extern ofstream file_basal_errorbars_top;
extern ofstream file_basal_errorbars_bottom;
extern ofstream file_total_errorbars;
extern ofstream file_basal_errorbars_resc;
extern ofstream file_basal_cumul_errorbars;
extern ofstream file_total_cumul_errorbars;
extern vector<ofstream*> file_errorbars_xmgr;
extern ofstream file_basal_errorbars_xmgr;

extern ofstream file_param_error;
extern ofstream file_param_stddev;

extern ifstream file_prior;

extern ofstream file_av_clonesize;
extern ofstream file_av_clonesize_basal;
extern ofstream file_mut_distr;

extern double floatingclone_frac_av;
extern vector<double> floatingclone_errorbar;

extern ofstream* file_loglikelies;

extern int fail_ctr;
extern int toolarge_ctr;

extern int discrete;

#endif
