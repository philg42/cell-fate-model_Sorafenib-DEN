/*
 *  global.cpp
 *  popdyn_1D_discr
 *
 *  Created by Philip Greulich on 07.10.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "global.h"

MainWindow* window;

//------------- Variable Simulation Parameters

#if param_scan == 1
int runnumber=0;
int runnumber2=1;
int esmblnumber=1;
#else
int runnumber=0;
int runnumber2=500;
int esmblnumber=1;
#endif

double siglevel = 0.68;
double siglevel2 = 0.68;

#if raw_data == 1 || raw_data == 30 || raw_data == 30
double reprodrate_bayes_av = 0.997*7;
double var_bayes = 0.031*7*0.031*7;
double r_bayes_av = 0.1;
double var_r_bayes = 0.01/2*0.01/2;
#elif raw_data == 11
double reprodrate_bayes_av = 0.970*7;
double var_bayes = 0.058*7*0.058*7;
#elif raw_data == 10 || raw_data == 24
double reprodrate_bayes_av = 0.994*7;
double var_bayes = 0.005*7*0.005*7;
#else
double reprodrate_bayes_av = 10;
double var_bayes = 100000;
double r_bayes_av = 0.1;
double var_r_bayes = 100000;
#endif

double runtime0=20.0/(4.0/9*7.0);

int celltypenumber=3;                    // Model-dependent !!!
#if Dcell_progenitor == 0
int prolif_celltype_limit=1;
#elif Dcell_progenitor == 1
int prolif_celltype_limit=2;
#endif
int shed_celltype_limit=2;
const int clonenumber_def=150000;
int clonenumber_bestfit=150000;

double T_max=1000;                       // Quasi-temperatur
double T=0;
double Tdropruns=200;
double Tdrop1=0.2;
double T0 = 100;

int evalintv=2.0/7;

int binsize=1;
char binmode='l';

double binsize_lesion=0.5;


#if init_cells == 1
int init_clone_number=1;
#elif init_cells == 2
int init_clone_number=2;
#endif
int basalcelltype_number=celltypenumber-1;
double prog_cell_frac_1=1;
double prog_cell_frac_2=0;       // Attention !!! This is overwritten by init2cell_frac
double prog_cell_frac_3=0;

#if raw_data == 27
int start_valid_b=3;
int start_valid_s=0;
#else
#if fitbasal == 3
int start_valid_b=2;//binsize+1;
int start_valid_s=1;
#elif raw_data == 10
int start_valid_b=3;//binsize+1;
int start_valid_s=0;
#elif fitbasal == 2 || raw_data == 18 || raw_data == 24 || raw_data == 26 || raw_data == 28 || raw_data == 29
int start_valid_b=3;//binsize+1;
int start_valid_s=0;
#elif raw_data == 0 || raw_data == 16 || raw_data == 17  || raw_data == 23
int start_valid_b=2;
int start_valid_s=0;            // If no clones below start_valid_b are counted type 100000, otherwise type lowest suprabasal number (i.e. 1)
#else
int start_valid_b=2;            // Attention !!! Adjust when counting floating clones
int start_valid_s=1;			// If no clones below start_valid_b are counted type 100000, otherwise type lowest suprabasal number (i.e. 1)
#endif
#endif

#if raw_data == 1
double SB_B_ratio = 1.25;
#elif raw_data == 11
double SB_B_ratio = 1.16;
#elif raw_data == 25
double SB_B_ratio = 1.37;
#elif raw_data == 0
#if plain == 0
double SB_B_ratio = 0.64;
#else
double SB_B_ratio = 0.64;
#endif
#else
double SB_B_ratio = 1.25;
#endif

// double reprodratelimit_low=7.3;          // Limits for parameter variation: fine
// double reprodratelimit_high=7.4;
// double stratratelimit_low=0;
// double stratratelimit_high=14;
// double shedratelimit_low=10.0;
// double shedratelimit_high=10.1;
//
// double stratrate_intv=1.0;
// double shedrate_intv=1.0;
// double reprodrate_intv=2.0;

double reprodratelimit_low=6.8;          // Limits for parameter variation, Soraf
double reprodratelimit_high=7.15;
double stratratelimit_low=1.2;//0.2;
double stratratelimit_high=3.2;//3.5;
double shedratelimit_low=7;
double shedratelimit_high=7.01;
double directstratlimit_low=0;
double directstratlimit_high=directstratlimit_low+0.01;
double init2proglimit_low=0.0;
double init2proglimit_high=1.05;

double reprodrate_intv=0.05;
double stratrate_intv=0.2;
double shedrate_intv=0.4;
double directstrat_intv=1;
double init2prog_intv=2;

// Alternative scanning with rho; stratrate not ued in that case

int strat_infinity_flag=0;

double rholimit_low=0.7;
double rholimit_high=1.01;
double rho_intv=0.05;

//double reprodratelimit_low=0.5;          // 2 param fit
//double reprodratelimit_high=10.0;
//double stratratelimit_low=0;
//double stratratelimit_high=0.1;
//double shedratelimit_low=0;
//double shedratelimit_high=0.1;
//
//double stratrate_intv=0.2;
//double shedrate_intv=0.3;
//double reprodrate_intv=0.1;

double stepsize_reprod = 0.3;       // Step sizes for parameter variation
double stepsize_r = 0.01;
double stepsize_strat = 0.2;
double stepsize_shed = 0.1;
double stepsize_directstrat = 0.01;
double stepsize_delta = 0.05;
double stepsize_progfrac = 0.05;

int scatterpt_number = 10000;

#if raw_data == 10
int excluded_clones = 0;
#elif raw_data == 18
int excluded_clones = 0;
#else
int excluded_clones = 0;
#endif


//-----------------------------------
// Initial parameter sets : Model dependent !!! Adjust to 
//-------------------------------------


#if many_step_process == 1
int reprod_step_number = 20;
#else
int reprod_step_number = 0;
#endif

double strat_rate_eff = 6.0;
double init2prog_frac0=0.9;

double mut_rate0 = 0.1;

#if raw_data == 0 || raw_data == 2 || raw_data == 3 || raw_data == 4 || raw_data == 16 || raw_data == 17  || raw_data == 20 || raw_data == 23 || raw_data == 26


double reprod_rate0 = 1.9;                         // Doupe Science paper parameters
double reprod_rate1 = 0;
double r0 = 0.0;
double delta0 = 0;
double strat_rate1=3.5;
double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0+strat_rate1)/SB_B_ratio;
double direct_strat_frac0=0;

//double reprod_rate0 = 1.9;                         // Doupe Science paper parameters direct shed
//double reprod_rate1 = 0;
//double r0 = 0.1;
//double delta0 = 0;
//double strat_rate1=0.8;
//double shedrate2=strat_rate1;
//double direct_strat_frac0=0;

//double reprod_rate0 = 1.5;                         // Soraf control, manual best fit all combined, r=0.09 fixed, init2prog_frac=0.75  
//double reprod_rate1 = 0;
//double r0 = 0.09;
//double delta0 = 0;
//double strat_rate1=7.5;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0+strat_rate1)/SB_B_ratio;
//double direct_strat_frac0=0;

//double reprod_rate0 = 2.5;                         // Soraf control, manual best fit all combined, r=0.09 fixed  
//double reprod_rate1 = 0;
//double r0 = 0.09;
//double delta0 = 0;
//double strat_rate1=12.5;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0+strat_rate1)/SB_B_ratio;
//double direct_strat_frac0=0.0;

//                        // Play
//#if Dcell_progenitors == 1
//double reprod_rate0 = 0;
//double reprod_rate1 = 5;
//#else
//double reprod_rate0 = 1.9;
//double reprod_rate1 = 0;
//#endif
//double r0 = 0.1;
//double delta0 = 0;
//double strat_rate1=3.5;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0+strat_rate1)/SB_B_ratio;
//double direct_strat_frac0=0;


//double reprod_rate0 = 3;       // Soraf control best fit manual basal
//double reprod_rate1 = 0;
//double r0 = 0.0;
//double delta0 = 0;
//double strat_rate1=2.6;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;
//double direct_strat_frac0=0;

//double reprod_rate0 = 2.6;       // Soraf control best fit joint
//double reprod_rate1 = 0;
//double r0 = 0.005;
//double delta0 = 0;
//double strat_rate1=2.2;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;
//double direct_strat_frac0=0;


//double reprod_rate0 = 10;       // Play
//double reprod_rate1 = 0;
//double r0 = 0.5;
//double delta0 = 0;
//double strat_rate1=10;
//double shedrate2=10;//reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;
//double direct_strat_frac0=0;


#elif raw_data == 1 || raw_data == 30 || raw_data == 25

//double reprod_rate0 = 7.3;  	   // Sorafenib best fit total
//double reprod_rate1 = 0;
//double r0 = 0.015;
//double delta0 = 0.0;
//double strat_rate1=2;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;//6.67109;//6.67109;
//double direct_strat_frac0=0.0;

//double reprod_rate0 = 7.3;  	   // Sorafenib best fit no distinct layers
//double reprod_rate1 = 0;
//double r0 = 0.015;
//double delta0 = 0.0;
//double strat_rate1=1.87;//6.67109;//14.438;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;//6.67109;
//double direct_strat_frac0=0;

double reprod_rate0 = 7.0;       // Sorafenib best fit joint
double reprod_rate1 = 0;
double r0 = 0.02;
double delta0 = 0;
double strat_rate1=6.8;
double shedrate2=2.6;
double direct_strat_frac0=0;

//double reprod_rate0 = 7.3;       // Sorafenib best fit basal
//double reprod_rate1 = 0;
//double r0 = 0.015;
//double delta0 = 0;
//double strat_rate1=2.2;
//double shedrate2=reprod_rate0*strat_rate_eff/(reprod_rate0 + strat_rate_eff)/SB_B_ratio;
//double direct_strat_frac0=0.63;

//double reprod_rate0 = 7.3;  	   // Sorafenib best fit direct strat
//double reprod_rate1 = 0;
//double r0 = 0.015;
//double delta0 = 0.0;
//double strat_rate1=6;
//double shedrate2=reprod_rate0*strat_rate_eff/(reprod_rate0 + strat_rate_eff)/SB_B_ratio;//6.67109;
//double direct_strat_frac0=0;//1-(double)strat_rate1/strat_rate_eff;

//double reprod_rate0 = 7.3;  	   // Sorafenib+DEN out 20d confetti best fit 
//double reprod_rate1 = 0;
//double r0 = 0.04;
//double delta0 = 0.0;
//double strat_rate1=2;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;//6.67109;
//double direct_strat_frac0=0;

#elif raw_data == 11

double reprod_rate0 = 6.8;       // DEN+Sorafenib out best fit
double reprod_rate1 = 0;
double r0 = 0.04;
double delta0 = 0;
double strat_rate1=6.4;
double shedrate2=2.4;//reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;
double direct_strat_frac0=0;

//double reprod_rate0 = 4.9;       // Play
//double reprod_rate1 = 0;
//double r0 = 0.07;
//double delta0 = 0;
//double strat_rate1=5.9;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;
//double direct_strat_frac0=0;

#elif raw_data == 5 || raw_data == 6 || raw_data == 7 || raw_data == 9 || raw_data == 13 || raw_data == 14 || raw_data == 15 || raw_data == 19 || raw_data == 21 || raw_data == 22



//double reprod_rate0 = 7.76015;        // MAML optimal start 3d
//double reprod_rate1 = 0;
//double r0 = 0;
//double delta0 = 1.0;
//double strat_rate1=3.24195;
//double shedrate2=8.7472;
//double direct_strat_frac0=0;

//double reprod_rate0=6.4;        // MAML best fit 7d
//double reprod_rate1 = 0;
//double r0=0.065;
//double delta0=0.9;
//double strat_rate1=1.1;
//double shedrate2=1.8;
//double direct_strat_frac0=0;


//double reprod_rate0 = 6.8;       // MAML best fit 10d
//double reprod_rate1 = 0;
//double r0 = 0.045;
//double delta0 = 1.0;
//double strat_rate1=0.8;
//double shedrate2=1.2;
//double direct_strat_frac0=0;

double reprod_rate0 = 6;       // MAML best fit 7d10d
double reprod_rate1 = 0;
double r0 = 0.055;
double delta0 = 1;
double strat_rate1=0.8;
double shedrate2=0.6;
double direct_strat_frac0=0;

//double reprod_rate0=5.5;        // MAML best fit 7d
//double reprod_rate1 = 0;
//double r0=0.065;
//double delta0=1.0;
//double strat_rate1=0.8;
//double shedrate2=0.4;
//double direct_strat_frac0=0;

//double reprod_rate0 = 5.5;       // MAML best fit 7d10d, no cut-off
//double reprod_rate1 = 0;
//double r0 = 0.065;
//double delta0 = 1;
//double strat_rate1=0.8;
//double shedrate2=5;
//double direct_strat_frac0=0;

//double reprod_rate0 = 5.5;       // MAML alt
//double reprod_rate1 = 0;
//double r0 = 0.065;
//double delta0 = 1;
//double strat_rate1=0;
//double shedrate2=5;
//double direct_strat_frac0=0;

//double reprod_rate0 = 5.4;       // MAML best fit 2 parameters
//double reprod_rate1 = 0;
//double r0 = 0.065;
//double delta0 = 1;
//double strat_rate1=0;
//double shedrate2=5;
//double direct_strat_frac0=0;


#elif raw_data == 10 || raw_data == 18 || raw_data == 24
#if confetti == 0


//double reprod_rate0 = 7.3;  	   // Sorafenib in lesion, alternative
//double reprod_rate1 = 0;
//double r0 = 0.03;
//double delta0 = 0.4;
//double strat_rate1=0.8;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;//6.67109;
//double direct_strat_frac0=0;

double reprod_rate0 = 6.95;  	   // Sorafenib in lesion, best fit 10d20d YFP
double reprod_rate1 = 0;
double r0 = 0.05;
double delta0 = 0.4;
double strat_rate1=0.9;
double shedrate2=10.0;//6.67109;
double direct_strat_frac0=0;

#else

//double reprod_rate0 = 6.95;  	   // Sorafenib in lesion, confetti 12d+22d
//double reprod_rate1 = 0;
//double r0 = 0.045;
//double delta0 = 0.3;
//double strat_rate1=0.8;
//double shedrate2=7.0;//6.67109;
//double direct_strat_frac0=0;

double reprod_rate0 = 6.95;  	   // Sorafenib in lesion, confetti alternative
double reprod_rate1 = 0;
double r0 = 0.03;
double delta0 = 0.7;
double strat_rate1=2;
double shedrate2=15.0;//6.67109;
double direct_strat_frac0=0;

//double reprod_rate0 = 7.3;  	   // Sorafenib+DEN bestfit 21d confetti; in lesion
//double reprod_rate1 = 0;
//double r0 = 0.04;
//double delta0 = 0.0;
//double strat_rate1=1;
//double shedrate2=15.0;//6.67109;
//double direct_strat_frac0=0;

#endif


#elif raw_data == 27
#if Dcell_progenitors == 1
double reprod_rate0 = 7.0;
double reprod_rate1 = 0;
#else
double reprod_rate0 = 7.0;
double reprod_rate1 = 0;
#endif
double r0 = 0;//2.0/9;
double delta0 = 0.0;
#if discrete_time == 0
double strat_rate1=7.0;//0.7;
#else
double strat_rate1=7.0;//0.7/(reprod_rate0+0.7)*reprod_rate0;
#endif
double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0+strat_rate1)/SB_B_ratio;
double direct_strat_frac0=0;

#elif raw_data == 28 || raw_data == 29

//double reprod_rate0 = 7.3;  	   // Sorafenib in lesion best fit 12d22d
//double reprod_rate1 = 0;
//double r0 = 0.03;
//double delta0 = 0.4;
//double strat_rate1=0.8;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;//6.67109;
//double direct_strat_frac0=0;

double reprod_rate0 = 6.95;  	   // Sorafenib in lesion play
double reprod_rate1 = 0;
double r0 = 0.045;
double delta0 = 0.3;
double strat_rate1=0.8;
double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;//6.67109;
double direct_strat_frac0=0;

//double reprod_rate0 = 7.3;  	   // Sorafenib in lesion no delta
//double reprod_rate1 = 0;
//double r0 = 0.05;
//double delta0 = 0.0;
//double strat_rate1=0.7;
//double shedrate2=reprod_rate0*strat_rate1/(reprod_rate0 + strat_rate1)/SB_B_ratio;//6.67109;
//double direct_strat_frac0=0;


#endif

///////////////////////////////////////////////////////////////////

#if reversible_diff == 1
double switchrate_forw0=reprod_rate0;
double switchrate_backw0=strat_rate1;
#else
double switchrate_forw0=0;
double switchrate_backw0=0;
#endif

double reprod_rate0_soraf = 15;
double r0_soraf = 0.016;
double delta0_soraf = 0;
double strat_rate1_soraf=13.49;
double shedrate2_soraf=6.67;
double direct_strat_frac0_soraf=0;

double reprod_rate0_hom = 1.9;
double r0_hom = 0.1;
double delta0_hom = 0;
double strat_rate1_hom = 3.5;
double shedrate2_hom = reprod_rate0*strat_rate1/(reprod_rate0+strat_rate1)/SB_B_ratio;
double direct_strat_frac0_hom=0;

double switch_reprodrate = 2;//reprod_rate0;
double switch_delta = 1;
double switch_r = 0.5;

//---------------------------------------

#if raw_data==2 || raw_data == 4 || raw_data == 6 || raw_data == 7 || raw_data == 9 || raw_data == 10 || raw_data == 18 || raw_data == 24 || raw_data == 25 || raw_data == 27
int maxCS_av=3000;
#elif raw_data == 20 || raw_data == 21 || raw_data == 26  || raw_data == 28 || raw_data == 29
int maxCS_av=100000;
#elif raw_data == 22
int maxCS_av=100;
#else
int maxCS_av=5000;
#endif
int N_max=maxCS_av*clonenumber_def;
int N_min=clonenumber_def/50;

int limit_basal0;
int limit_supra0;

#if raw_data == 0 || raw_data == 1 || raw_data == 30 || raw_data == 13 || raw_data == 14 || raw_data == 15 || raw_data == 16  || raw_data == 17    // Limits for evaluation of Likelihood
#if fitbasal == 0 || fitbasal == 1
int limit_basal=10;
int limit_supra=20;
#else
int limit_basal=20;
int limit_supra=30;
#endif
#elif raw_data == 2
int limit_basal=20;
int limit_supra=20;
#elif raw_data == 3
int limit_basal=10;
int limit_supra=10;
#elif raw_data == 4 || raw_data == 8
int limit_basal=20;
int limit_supra=20;
#elif raw_data == 11
int limit_basal=13;
int limit_supra=20;
#elif raw_data == 5
int limit_basal=25;
int limit_supra=15;
#elif raw_data == 6
int limit_basal=17;
int limit_supra=8;
#elif raw_data == 7
int limit_basal=26;
int limit_supra=13;
#elif raw_data == 9 ||  raw_data == 20 || raw_data == 21 
int limit_basal_1=26;//16;
int limit_supra_1=12;//8;
int limit_basal_2=26;
int limit_supra_2=12;
int limit_basal=10000;
int limit_supra=10000;
#elif raw_data == 10
int limit_basal=45;
int limit_supra=100;
#elif raw_data == 18
int limit_basal=80;
int limit_supra=80;
#elif raw_data == 23 || raw_data == 26
int limit_basal=15;
int limit_supra=20;
#elif raw_data == 24
#if confetti == 1
int limit_basal_1=20;//16;
int limit_supra_1=30;//8;
int limit_basal_2=30;
int limit_supra_2=80;
int limit_basal=10000;
int limit_supra=1000;
#else
int limit_basal_1=45;//16;
int limit_supra_1=30;//8;
int limit_basal_2=80;
int limit_supra_2=80;
int limit_basal=10000;
int limit_supra=1000;
#endif
#elif raw_data == 19 || raw_data == 22
int limit_basal=40;
int limit_supra=40;
#elif raw_data == 25
int limit_basal=40;
int limit_supra=40;
#elif raw_data == 27
int limit_basal=40;
int limit_supra=40;
#elif raw_data == 28
int limit_basal = 70;
int limit_supra = 100;
#elif raw_data == 29
int limit_basal = 280;
int limit_supra = 100;
#endif

//------------- End Simulation Paramters ----------------------------------------------------

char fileoption='d';
int index_run=20;

int itr_intv=2000;
int output_intv=500;
int image_intv=1;

#if raw_data == 9 || raw_data == 20 || raw_data == 24
int timepoints = 2;   
#else
int timepoints = 1;
#endif

#if raw_data == 20 || raw_data == 21 || raw_data == 26
vector<double> timepoints_av;
#endif

#if raw_data == 1
double t_switch = 400000.0/7 + 0.0/7;
#elif raw_data == 2
double t_switch = 22.0/7 + 0.0/7;
#elif raw_data == 3 || raw_data == 4
double t_switch = 4.0/7;
#elif raw_data == 22
double t_switch = 11.0/7;
#else
double t_switch = 400000.0/7 + 0.0/7;
#endif

#if raw_data == 27
double runtime = runtime0;
#elif (raw_data == 1 || raw_data == 0) && plain == 0
#if init_cells > 1
double runtime = (12.0-7.0/reprod_rate0)/7;
#else
double runtime = 12.0/7;
#endif
#elif (raw_data == 1 || raw_data == 30 || raw_data == 0) && plain == 1
#if init_cells > 1
double runtime = (10.0-7.0/reprod_rate0)/7;
#else
double runtime = 10.0/7;
#endif
#elif raw_data == 10
#if confetti == 1
double runtime = 10.0/7;
#else
double runtime = 13.0/7;
#endif
#elif raw_data == 11
#if init_cells > 1
double runtime = (13.0-7.0/reprod_rate0)/7;
#else
double runtime = 13.0/7;
#endif
#elif raw_data == 2 || raw_data == 4
double runtime = 32.0/7;
#elif raw_data == 3
double runtime = 14.0/7;
#elif raw_data == 5 || raw_data == 16
#if init_cells > 1
double runtime = (3.0-7.0/reprod_rate0*log(2))/7;
#else
double runtime = 3.0/7;
#endif
#elif raw_data == 6 || raw_data == 17
#if init_cells > 1
double runtime = (7.0-7.0/reprod_rate0*log(2))/7;
#else
double runtime = 7.0/7;
#endif
#elif raw_data == 7 || raw_data == 9
#if init_cells > 1
double runtime = (10.0-7.0/reprod_rate0)/7;
#else
double runtime = 10.0/7;
#endif
#elif raw_data == 8
double runtime = 18.0/7;
#elif raw_data == 13 
#if init_cells > 1
double runtime = (1.0-7.0/reprod_rate0/2)/7;
#else
double runtime = 1.0/7;
#endif
#elif raw_data == 14 
#if init_cells > 1
double runtime = (2.0-7.0/reprod_rate0/2)/7;
#else
double runtime = 2.0/7;
#endif
#elif raw_data == 15 
#if init_cells > 1
double runtime = (3.0-7.0/reprod_rate0*1.4)/7;
#else
double runtime = 3.0/7;
#endif
#elif raw_data == 18 || raw_data == 24 || raw_data == 25
#if confetti == 1
double runtime = 21.0/7;
#else
double runtime = 22.0/7;
#endif
#elif raw_data == 20 || raw_data == 21	|| raw_data == 26
#if init_cells > 1
double runtime = (51-7.0/reprod_rate0*1.4)/7;
#else
double runtime = 51.0/7;
#endif
#elif raw_data == 19 || raw_data == 22
#if init_cells > 1
double runtime = (15.0-7.0/reprod_rate0*1.4)/7;
#else
double runtime = 15.0/7;
#endif
#elif raw_data == 23 || raw_data == 26
#if init_cells > 1
double runtime = (12*7-7.0/reprod_rate0*1.4)/7;
#else
double runtime = 12*7.0/7;
#endif
#elif raw_data == 28
double runtime = 42.0/7;
#elif raw_data == 29
double runtime = 56.0/7;
#elif raw_data == 30
double runtime = 21.0/7;
#endif



//int logflag=0;     // Debug !!!

int run_index=10;
int extinct_flag=0;

cell* curr_cell=NULL;

int image_ctr=0;
//double time_sim = 0.0;


const int param_number=5;                       // Model-dependent !!!
vector<int> var_vector(param_number);       // position 0 = reprod, 1 = r, 2 = delta, 3 = strat, 4 = shed

int clone_norm=100;

#if raw_data == 0
vector<int> clonenumber_exp(1,253);    // for Science paper data
//int clonenumber_exp=296;      // for MAML YFP 10d 
//int clonenumber_exp=232;		// for Sorafenib control
#else
vector<int> clonenumber_exp(timepoints);
#endif
////#if raw_data == 0
////int clonenumber_exp=296;      // for MAML YFP 10d    
//#elif raw_data == 1
//int clonenumber_exp=233;    // pencil: 233 // 4d+28d: 304       // Attention !!! Data-dependent 
//#elif raw_data == 2
//int clonenumber_exp=288;
//#elif raw_data == 3
//int clonenumber_exp=270;
//#elif raw_data == 4 || raw_data == 8
//int clonenumber_exp=304;
//#elif raw_data == 5 || raw_data == 6
//int clonenumber_exp=300;
//#elif raw_data == 7
//int clonenumber_exp=250;
//#elif raw_data == 9
//int clonenumber_exp=250;             // choose clonenumber=250 when testing 10d and 300 otherwise 
//#elif raw_data == 10
//int clonenumber_exp=100;
//#elif raw_data == 11
//int clonenumber_exp=245;
//#else
//int clonenumber_exp=0;
//#endif

//vector<int> clonenumber_exp_time(timepoints);
vector<vector<vector<double> > > data_cumul;



#if raw_data == 9
double logtime1 = 7.0/7;
#elif raw_data == 20 || raw_data == 21 || raw_data == 26
double logtime1 = 7.0/7;
double logtime2 = 10.0/7;
double logtime3 = 30.0/7;
int logtime4 = 50.0/7;
#elif raw_data == 24
#if confetti == 1
double logtime1 = 10.0/7;
#else
double logtime1 = 13.0/7;
#endif
#else
double logtime1 = 200000;
#endif

int validclonecount=0;

vector<run*> run_old;

double SCfrac=1;

bool(*fn_pt)(vector<double>,vector<double>) = compare_pairs;
map<vector<double>,double,bool(*)(vector<double>,vector<double>)> prior_map(fn_pt);

//vector<cell*> cells;           // cell[i][0]=clone of cell[i]; microbes[i][1] = type of microbe[i]
//vector<Clone> clones;

vector<double> reprod_rate_init(celltypenumber,0); 
vector<double> r_init(celltypenumber,0);
vector<double> dediff_rate_init(celltypenumber,0);
vector<double> delta_init(celltypenumber,0);
vector<double> strat_rate_init(celltypenumber,0);
vector<double> shed_rate_init(celltypenumber,0);
vector<double> direct_strat_frac_init(celltypenumber,0);
vector<double> switchrate_forw_init(celltypenumber,0);
vector<double> switchrate_backw_init(celltypenumber,0);
vector<double> init2prog_frac_init(celltypenumber,0);

double delay=0;

//double maxrate=0;

int cellnumber_tot=0;
vector<int> Cnumber_glob(celltypenumber);

//vector<double> clonesizedistr;
//vector<double> clonesizedistr_norm;
//vector<double> clonesizedistr_basal;
//vector<double> clonesizedistr_supra;
//vector<vector<double> > joint_distr;



//int maxsize=0;

const double norm=0.5/2147483647;              // Normalization factor for random numbers
int ibm=252;                               // random number seed. Attention: fixed in current version !!!
int seed=ibm*2+1;
double rand1;
double rand2;

vector<int> bins;

vector<vector<double> > joint_distr_av;
vector<double> CSD_basal_av;
vector<double> CSD_total_av;
vector<vector<double> > lesion_SD;
vector<double> cumuldistr_basal_av;
vector<double> cumuldistr_total_av;
unsigned int i_curr;
unsigned int j_curr;

ofstream file_CSD;
ofstream file_CSD_norm;
ofstream file_CSD_basal;
ofstream file_CSD_supra;
ofstream file_CSD_basal_resc;
ofstream file_CSD_total_resc;
ofstream file_cumuldistr_total;
ofstream file_cumuldistr_basal;
ofstream file_cumuldistr_supra;
ofstream file_Jointdistr;
ofstream file_Jointdistr_norm;
ofstream file_Jointdistr_matrix;
ofstream file_Jointdistr_transp;
ofstream file_parameters;
ofstream file_pval;
ofstream file_loglikely;
ofstream file_av;
ofstream file_avsupra_basal;
ofstream file_celltypes_basal;
ofstream file_clonenumber_t;

ofstream file_data_basal_out;
ofstream file_data_total_out;
ofstream file_data_basal_resc_out;
ofstream file_data_total_resc_out;
ofstream file_data_avsupra_out;
ofstream file_data_cumul_basal_out;
ofstream file_data_cumul_total_out;
ofstream file_lesion_SD_out;

ofstream file_CSD_basal_intm;
ofstream file_Jointdistr_intm;
ofstream file_Jointdistr_matrix_intm;

ofstream file_H2B_divctr_scatter;
ofstream file_H2B_av_clonesize;

ofstream file_temp;         // Debug !!!

ofstream file_av_clonesize;
ofstream file_av_clonesize_basal;
ofstream file_mut_distr;

ifstream file_data;
ifstream file_data_total;
ifstream file_data_init;
ifstream file_data_basal;

#if raw_data == 9 || raw_data == 20 || raw_data == 21 || raw_data == 24 || raw_data == 26
ifstream file_data2;
ifstream file_data_basal2;
#endif

#if raw_data == 24
ifstream file_data_total2;
#endif

#if raw_data == 10 || raw_data == 18
ifstream file_data_lesions;
#endif

int fileerrorbar_nr=13;
vector<ofstream*> file_errorbars;
ofstream file_basal_errorbars;
ofstream file_basal_errorbars_top;
ofstream file_basal_errorbars_bottom;
ofstream file_basal_errorbars_resc;
ofstream file_total_errorbars;
ofstream file_basal_cumul_errorbars;
ofstream file_total_cumul_errorbars;
vector<ofstream*> file_errorbars_xmgr;
ofstream file_basal_errorbars_xmgr;

ofstream file_param_error;
ofstream file_param_stddev;

ifstream file_prior;

double floatingclone_frac_av;
vector<double> floatingclone_errorbar(3);

ofstream* file_loglikelies;

vector<vector<int> > data_init;
vector<vector<vector<double> > > data0(timepoints);
vector<vector<vector<double> > > data_basal(timepoints);
vector<vector<vector<double> > > data_total(timepoints);

int fail_ctr=0;
int toolarge_ctr=0;

#if discrete_time == 1
int discrete = 1;
#else
int discrete = 0;
#endif

//------ to do for parameter scan -----------
//	switch param_scan = 1 in global.h
//	switch get_bestparams = 0 in global.h
//	adjust scan range and intervals of outer parameters in script_simu_param (on cluster)
//	check order of parameters to prog (on cluster)
//	adjust scan range of inner most parameter in global.cpp, top
//	
//
//


