/*
 *  MainWindow.cpp
 *  test2
 *
 *  Created by Philip Greulich on 22.07.10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
// Class MainWindow
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
/* This class is the main interface for simulations. It starts and stops the 
 simulation or provides a graphical interface to do things by hand. The Widgets are set up
 by giving the window layout and connecting signals with executive functions. MainWindow also 
 initiates all other objects and manages the update process by calling other functions (see func.cpp).
 */


//#include <iostream>
//
//#if Qt_flag==1
//#include <QtGui>
//#include <QTimer>
//#include <QSlider>
//#endif
//
//#include <vector>
//#include <list>
//#include <set>
//#include <cmath>
//#include <fstream>
//#include <string>
//
using namespace std;
//
#include "global.h"

#include "MainWindow.h"



//#if Qt_flag==1
//#include "PaintWindow.h"
//#include "LabeledSlider.h"
//#include "LabeledDisplay.h"
//#endif
//
//#include "grid.h"


#define Pi 3.141593


//-------------------------------------------------
// Constructor 
//-------------------------------------------------

#if Qt_flag==1 

MainWindow::MainWindow(QWidget *parent) : QDialog(parent)
{
	setup_sim();
	
	setup_widgets();
};
#else
MainWindow::MainWindow()
{
//setup_sim();
	
};
#endif

/////////////// Perform simulation without Qt /////////////////////////////
void MainWindow::start_sim()
{	
	start();

	
#if raw_data == 20 || raw_data == 21 || raw_data == 26
	timepoints_av.resize(10);
	timepoints_av[0]=7.0/7;
	timepoints_av[1]=10.0/7;
	timepoints_av[2]=15.0/7;
	timepoints_av[3]=20.0/7;
	timepoints_av[4]=25.0/7;
	timepoints_av[5]=30.0/7;
	timepoints_av[6]=40.0/7;
	timepoints_av[7]=50.0/7;
	timepoints_av[8]=60.0/7;
	timepoints_av[9]=70.0/7;
#endif

#if Qt_flag == 0
	read_files_to_data();
	
	prepare_data();
	vector<vector<double> > parameters;
	prepare_parameters(parameters);
	
#if get_bestparam == 1
	parameters=get_best_params(parameters);	
#endif
	
	fail_ctr=0;
#if param_scan == 1
	scan_likelihood(parameters,data0);
#elif fitbasal == 0
#if shedrate_fixed == 1
	shed_rate_init[2] = reprod_rate_init[0]*strat_rate_init[1]/(reprod_rate_init[0]+strat_rate_init[1])/SB_B_ratio;
#endif
	run* best_run = fit_evol(parameters,data0);
#if raw_data != 20 && raw_data != 21 && raw_data != 26
	find_errorbars(best_run,data0);
#endif
#elif fitbasal == 1 || fitbasal == 2 || fitbasal == 3
#if shedrate_fixed == 1
	shed_rate_init[2] = reprod_rate_init[0]*strat_rate_init[1]/(reprod_rate_init[0]+strat_rate_init[1])/SB_B_ratio;
#endif
	run* best_run = fit_evol(parameters,data_basal);
#if raw_data != 20 && raw_data != 21 && raw_data != 26
	find_errorbars(best_run,data_basal);
#endif

#endif
	
#endif

};

//-----------------------------------------------------------
//-----------------------------------------------------------
//     Functions for setting up simulation
//-----------------------------------------------------------
//-----------------------------------------------------------

//----------------- Read files and prepare them for internal handling

void MainWindow::read_files_to_data()
{
	
	file_CSD.open("CSD.dat");
	file_CSD_norm.open("CSD_norm.dat");
	file_CSD_basal.open("CSD_basal.dat");
	file_CSD_supra.open("CSD_supra.dat");
	file_CSD_basal_resc.open("CSD_basal_resc.dat");
	file_CSD_total_resc.open("CSD_total_resc.dat");
	file_Jointdistr.open("Jointdistr.dat");
	file_cumuldistr_basal.open("cumuldistr_basal.dat");
	file_cumuldistr_total.open("cumuldistr_total.dat");
	file_Jointdistr_norm.open("Jointdistr_norm.dat");
	file_Jointdistr_matrix.open("Jointdistr_matrix.dat");
	file_Jointdistr_transp.open("Jointdistr_data_transp.dat");
	file_parameters.open("parameters.dat");
	file_pval.open("pvals.dat");
	file_loglikely.open("loglikely.dat");
	file_av.open("averages.dat");
	file_avsupra_basal.open("average_supra_basal.dat");
    file_celltypes_basal.open("celltypes_basal.dat");
    file_clonenumber_t.open("clonenumber_t.dat");
	
	file_CSD_basal_intm.open("CSD_basal_intm.dat");
	file_Jointdistr_intm.open("Jointdistr_intm.dat");
	file_Jointdistr_matrix_intm.open("Jointdistr_matrix_intm.dat");
	
	file_lesion_SD_out.open("lesions_SD_data_curr.dat");
	
	file_H2B_divctr_scatter.open("H2B_divctr_scatter.dat");
	file_H2B_av_clonesize.open("H2B_av_clonesize.dat");
	
	file_temp.open("temp.dat");

	file_av_clonesize.open("av_clonesize.dat");
	file_av_clonesize_basal.open("av_clonesize_basal.dat");
    
    file_mut_distr.open("mut_distr.dat");
	
	string str="errorbars";
	file_errorbars=make_filelist(str,fileerrorbar_nr);
#if discrete_time == 1
    file_basal_errorbars.open("basal_errorbars_discr.dat");
#else
	file_basal_errorbars.open("basal_errorbars.dat");
#endif
	file_basal_errorbars_top.open("basal_errorbars_top.dat");
	file_basal_errorbars_bottom.open("basal_errorbars_bottom.dat");
	file_basal_errorbars_resc.open("basal_errorbars_resc.dat");
	file_total_errorbars.open("total_errorbars.dat");
	file_basal_cumul_errorbars.open("basal_cumul_errorbars.dat");
	file_total_cumul_errorbars.open("total_cumul_errorbars.dat");
	str="errorbars_xmgr";
	file_errorbars_xmgr=make_filelist(str,fileerrorbar_nr);
	file_basal_errorbars_xmgr.open("basal_errorbars_xmgr.dat");
	
#if get_bestparam == 1
	file_param_error.open("param_error.dat");
	file_param_stddev.open("param_stddev.dat");
#endif
	
#if param_scan == 1
	str = "loglikelies_";
	file_loglikelies = make_file(str,(int)(delta0*100+0.001),(int)(r0*1000+0.001));
#endif
	
	
#if raw_data == 0 || raw_data == 27									    // Attention !!! raw_data==0 only temporary since no suprabasal data for plain available yet
	//file_data.open("jointdistr_data_YFP_10d.dat");
	file_data.open("jointdistr_data_pencil_control.dat");
	//file_data_basal.open("CSD_basal_data_plain.dat");
	//file_data_basal.open("CSD_basal_data_YFP_10d.dat");
	file_data_basal.open("CSD_basal_data_pencil_control.dat");
	file_data_total.open("CSD_total_data_pencil_control.dat");
#elif raw_data == 1      
	file_data.open("jointdistr_data_soraf_10d.dat");
	file_data_basal.open("CSD_basal_data_soraf_pencil.dat");
	file_data_total.open("CSD_data_soraf_pencil.dat");
#elif raw_data == 2
	file_data.open("jointdistr_data_22d+10d.dat");
	file_data_basal.open("CSD_basal_data_22d+10d.dat");
#elif raw_data == 3
	file_data.open("jointdistr_data_soraf_4+10d.dat");
	file_data_basal.open("CSD_basal_data_4d+10d.dat");
#elif raw_data == 4
	file_data.open("jointdistr_data_soraf_4+28d.dat");
	file_data_basal.open("CSD_basal_data_4d+28d.dat");
#elif raw_data == 5
	file_data.open("jointdistr_data_MAML_3d.dat");
	file_data_total.open("CSD_total_data_MAML_3d.dat");
	file_data_basal.open("CSD_basal_data_MAML_3d.dat");
	
#elif raw_data == 6
	file_data.open("jointdistr_data_MAML_7d.dat");
	file_data_basal.open("CSD_basal_data_MAML_7d.dat");
#elif raw_data == 7
	file_data.open("jointdistr_data_MAML_10d.dat");
	file_data_basal.open("CSD_basal_data_MAML_10d.dat");
#elif raw_data == 8
	file_data_init.open("jointdistr_data_soraf_4+10d.dat");
	file_data.open("jointdistr_data_soraf_4+28d.dat");
	file_data_basal.open("CSD_basal_data_soraf_4+28d.dat");
#elif raw_data == 9 || raw_data == 20 || raw_data == 21
	file_data.open("jointdistr_data_MAML_7d.dat");
	file_data2.open("jointdistr_data_MAML_10d.dat");
	file_data_basal.open("CSD_basal_data_MAML_7d.dat");
	file_data_basal2.open("CSD_basal_data_MAML_10d.dat");
#elif raw_data == 10
#if confetti == 1
	file_data.open("jointdistr_data_DEN+Soraf_in_confetti_10d.dat");
	file_data_basal.open("CSD_basal_data_DEN_in_10d_confetti.dat");
	file_data_total.open("CSD_total_data_DEN+Soraf_in_10d_confetti.dat");  
#else
	file_data.open("jointdistr_data_DEN+Soraf_in_10d.dat");
	file_data_basal.open("CSD_basal_data_DEN+Soraf_in_10d.dat");
	file_data_total.open("CSD_total_data_DEN+Soraf_in_10d.dat");
#endif
	
// file_data_lesions.open("lesions_sizes_control_28d.dat");
	file_data_lesions.open("lesions_sizes_28d.dat");
//	file_data_lesions.open("lesions_sizes_28+10d.dat");
//	file_data_lesions.open("lesions_sizes_28+20d.dat");
#elif raw_data == 11
	file_data.open("jointdistr_data_DEN+Soraf_out_10d.dat");
	file_data_basal.open("CSD_basal_data_DEN+Soraf_out_10d.dat");
#elif raw_data == 12
	file_data.open("CSD_basal.dat");
	file_data_basal.open("Jointdistr_matrix.dat");
#elif raw_data == 13
	file_data.open("jointdistr_EdU_MAML_24h.dat");
	file_data_total.open("CSD_total_EdU_MAML_24h.dat");
	file_data_basal.open("CSD_EdU_MAML_24h.dat");
#elif raw_data == 14
	file_data.open("jointdistr_EdU_MAML_48h.dat");
	file_data_total.open("CSD_total_EdU_MAML_48h.dat");
	file_data_basal.open("CSD_EdU_MAML_48h.dat");
#elif raw_data == 15
	file_data.open("jointdistr_EdU_MAML_72h.dat");
	file_data_total.open("CSD_total_EdU_MAML_72h.dat");
	file_data_basal.open("CSD_EdU_MAML_72h.dat");
#elif raw_data == 16
	file_data.open("jointdistr_YFP_3d.dat");
	file_data_total.open("CSD_total_YFP_7d.dat");
	file_data_basal.open("CSD_YFP_3d.dat");
#elif raw_data == 17
	file_data.open("jointdistr_YFP_7d.dat");
	file_data_total.open("CSD_total_YFP_7d.dat");
	file_data_basal.open("CSD_YFP_7d.dat");
#elif raw_data == 18
#if confetti == 1
	file_data.open("jointdistr_data_DEN+Soraf_in_confetti_20d.dat");             // Attention !!! No detailed data available yet, take 10d as dummy
	file_data_basal.open("CSD_basal_data_DEN+Soraf_in_20d_confetti.dat");
	file_data_total.open("CSD_total_data_DEN_in_20d_confetti.dat");
#else
	file_data.open("jointdistr_data_DEN+Soraf_in_10d.dat");             // Attention !!! No detailed data available yet, take 10d as dummy
	file_data_basal.open("CSD_basal_data_DEN+Soraf_in_20d.dat");
	file_data_total.open("CSD_total_data_DEN_in_20d.dat");
#endif
	file_data_lesions.open("lesions_sizes_28d.dat");
#elif raw_data == 19  || raw_data == 22
	file_data.open("jointdistr_data_MAML_15d.dat");
	file_data_basal.open("CSD_basal_data_MAML_15d.dat");
#elif raw_data == 23 || raw_data == 26  
	file_data.open("jointdistr_data_YFP_10d.dat");                // Attention !!! No data for joint distributions available; YFP 10d as dummy
	file_data_basal.open("CSD_basal_data_plain_3m.dat");
	file_data_total.open("CSD_basal_data_plain_3m.dat");               // Attention !!! No data for total distributions available; basl distr as dummy
#elif raw_data == 24
#if confetti == 1
    file_data.open("jointdistr_data_DEN+Soraf_in_confetti_10d.dat");
    file_data_basal.open("CSD_basal_data_DEN_in_10d_confetti.dat");
    file_data_total.open("CSD_total_data_DEN+Soraf_in_10d_confetti.dat");
    file_data2.open("jointdistr_data_DEN+Soraf_in_confetti_20d.dat");
    file_data_basal2.open("CSD_basal_data_DEN+Soraf_in_20d_confetti.dat");
    file_data_total2.open("CSD_total_data_DEN_in_20d_confetti.dat");
#else
	file_data.open("jointdistr_data_DEN+Soraf_in_10d.dat");           // Attention !!! No detailed data available yet, take 10d as dummy 
	file_data2.open("jointdistr_data_DEN+Soraf_in_10d.dat");
	file_data_basal.open("CSD_total_data_DEN+Soraf_in_10d.dat");
	file_data_basal2.open("CSD_total_data_DEN_in_20d.dat");
	file_data_total.open("CSD_total_data_DEN+Soraf_in_10d.dat");
	file_data_total2.open("CSD_total_data_DEN_in_20d.dat");
#endif
#elif raw_data == 25
	file_data.open("jointdistr_data_DEN+Soraf_out_confetti_20d.dat");
	file_data_basal.open("CSD_total_data_DEN+Soraf_out_confetti_20d.dat");
	file_data_total.open("CSD_total_data_DEN+Soraf_out_confetti_20d.dat");
#elif raw_data == 28
    file_data.open("jointdistr_data_DEN+Soraf_in_10d.dat");             // Attention !!! No detailed data available yet, take 10d as dummy
    file_data_basal.open("CSD_basal_data_DEN+Soraf_in_20d.dat");
    file_data_total.open("CSD_DEN+Ras_42d.dat");
#elif raw_data == 29
    file_data.open("jointdistr_data_DEN+Soraf_in_10d.dat");             // Attention !!! No detailed data available yet, take 10d as dummy
    file_data_basal.open("CSD_basal_data_DEN+Soraf_in_20d.dat");
    file_data_total.open("CSD_DEN+Ras_56d.dat");
#elif raw_data == 30
    file_data.open("jointdistr_data_soraf_21d.dat");
    file_data_basal.open("CSD_basal_data_soraf_pencil_21d.dat");
    file_data_total.open("CSD_data_soraf_pencil_21d.dat");
#endif
	
	file_data_basal_out.open("CSD_basal_data_curr.dat");
	file_data_total_out.open("CSD_total_data_curr.dat");
	file_data_basal_resc_out.open("CSD_basal_data_resc_curr.dat");
	file_data_total_resc_out.open("CSD_total_data_resc_curr.dat");
	file_data_avsupra_out.open("avsupra_basal_data.dat");
	file_data_cumul_basal_out.open("cumuldistr_basal_data_curr.dat");
	file_data_cumul_total_out.open("cumuldistr_total_data_curr.dat");
	
	file_prior.open("prior.dat");
	
	
	///////// Create Prior Table /////////////////////////
	cerr << "start reading prior\n";
	
	vector<vector<double> > prior_matrix=read_file(&file_prior);
	
	for (int i=0;i<prior_matrix.size(); i++)
	{
		for (int j=0; j<prior_matrix[i].size(); j++) 
		{
			cerr << prior_matrix[i][j] << '\t';
		};
		cerr << endl;
	};
	
	cerr << endl;
	cerr << "create prior map:\n";
	
	for (int i=0; i<prior_matrix.size(); i++)
	{
		vector<double> params;
		params.push_back(prior_matrix[i][0]);
		params.push_back(prior_matrix[i][1]);
		
		pair<vector<double>,double> pair0(params,prior_matrix[i][2]);
		prior_map.insert(pair0);
		
		
//		cerr << params[0] << '\t' << params [1] << '\t' << prior_map.at(params) << endl;
	};
	/////////////////////////////////////////////
	

	//exit(0);
	
#if raw_data == 8
	
	data_init.clear();
	data_init = read_file_int(&file_data_init);
	
	for (int i=0; i<data_init.size(); i++)             
	{
		for (int j=0; j<data_init[i].size(); j++)
		{
			cerr << data_init[i][j] << '\t';
		};
		cerr << endl;
	};
	
	
	
#endif
	
	///////// Read data files ////////////////////////////////////////////
	
	data0[0] = read_file(&file_data);
	
	cerr << "Check fit data file:\n";                 // Debug !!!
	for (int i=0; i<data0[0].size(); i++)                  
	{
		for (int j=0; j<data0[0][i].size(); j++) 
		{
			cerr << data0[0][i][j] << '\t';
		};
		cerr << endl;
	};
	
	data_basal[0] = read_file(&file_data_basal);
	
	cerr << "Check fit basal data file:\n";                 // Debug !!!
	for (int i=0; i<data_basal[0].size(); i++)                  
	{
		for (int j=0; j<data_basal[0][i].size(); j++) 
		{
			cerr << data_basal[0][i][j] << '\t';
		};
		cerr << endl;
	};
	
	data_total[0] = read_file(&file_data_total);
    
    cerr << "Check fit total data file:\n";                 // Debug !!!
    for (int i=0; i<data_total[0].size(); i++)
    {
        for (int j=0; j<data_total[0][i].size(); j++)
        {
            cerr << data_total[0][i][j] << '\t';
        };
        cerr << endl;
    };
	
#if raw_data == 9 || raw_data == 20 || raw_data == 21 || raw_data == 24
	
	data0[1] = read_file(&file_data2);
	
	cerr << "Check fit data file:\n";                 // Debug !!!
	for (int i=0; i<data0[1].size(); i++)                  
	{
		for (int j=0; j<data0[1][i].size(); j++) 
		{
			cerr << data0[1][i][j] << '\t';
		};
		cerr << endl;
	};
	
	data_basal[1] = read_file(&file_data_basal2);
	
	cerr << "Check fit basal data file:\n";                 // Debug !!!
	for (int i=0; i<data_basal[1].size(); i++)                  
	{
		for (int j=0; j<data_basal[1][i].size(); j++) 
		{
			cerr << data_basal[1][i][j] << '\t';
		};
		cerr << endl;
	};
	
#if raw_data == 24
	
	data_total[1] = read_file(&file_data_total2);
	
	cerr << "Check fit basal data file:\n";                 // Debug !!!
	for (int i=0; i<data_total[1].size(); i++)                  
	{
		for (int j=0; j<data_total[1][i].size(); j++) 
		{
			cerr << data_total[1][i][j] << '\t';
		};
		cerr << endl;
	};
#endif
	
#endif
	
};

//-----------------------------------------------------------

//-------------- Normalize data and prepare basal data from full distributions -----------------

void MainWindow::prepare_data()
{
    
	//------------------ Normalize Data -----------------------------------------
	data_cumul.clear();
	data_cumul.resize(timepoints);
	vector<double> data_counter(timepoints,0);
	
	
	
#if plain == 1
	for (int t=0; t<timepoints; t++) 
	{
//#if raw_data == 9  || raw_data == 24
//		if (t==0)
//		{
//			limit_basal=min(limit_basal_1,(int)data0[t].size());
//			limit_supra=min(limit_supra_1,(int)data0[t][0].size());
//		}
//		else if (t==1)
//		{
//			limit_basal=min(limit_basal_2,(int)data0[t].size());
//			limit_supra=min(limit_supra_2,(int)data0[t][0].size());
//		};
//#else
		limit_basal0=data_basal[t].size();
		limit_supra0=data_basal[t][0].size();
//#endif
		
		for (int i=start_valid_b-1; i<limit_basal0; i++) 
		{
			data_counter[t] = data_counter[t] + data_basal[t][i][1];	
		};
	};
#elif raw_data == 18 || raw_data == 24 || raw_data == 10 || raw_data == 28 || raw_data == 29 || fitbasal == 2
	for (int t=0; t<timepoints; t++) 
	{

		limit_basal0=data_total[0].size();
		limit_supra0=data0[t][0].size();
		
		for (int i=(start_valid_b-1); i<limit_basal0; i++)
		{
			data_counter[t] = data_counter[t] + data_total[t][i][1];	
		};
	};
#else
	for (int t=0; t<timepoints; t++) 
	{

        limit_basal0=min(limit_basal,(int)data0[t].size());
        limit_supra0=min(limit_supra,(int)data0[t][0].size());
		
		data_cumul[t].resize(data0[t].size());
		for (int i=0; i<limit_basal0; i++) 
		{
			data_cumul[t][i].resize(data0[t][i].size());
			for (int j=start_val(i); j<limit_supra0; j++) 
			{
				data_counter[t] = data_counter[t] + data0[t][i][j];	
				data_cumul[t][i][j] = data_counter[t];
			};
		};
	};
#endif
	
#if plain != 1 && raw_data != 18  && raw_data != 24  && raw_data != 10 && raw_data != 28 && raw_data != 29 && fitbasal != 2       // Exclude data which does not have detailed CSDs (tumour data)
	
	for (int t = 0; t < timepoints ; t++) 
	{
		clonenumber_exp[t] = data_counter[t];
	};
	for (int t=0; t<timepoints; t++) 
	{
		for (int i=0; i<data0[t].size(); i++) 
		{
			for (int j=start_val(i); j<data0[t][i].size(); j++) 
			{
				data0[t][i][j] = (double)data0[t][i][j]/data_counter[t]*clone_norm;
			};
			if (i<data_basal[t].size()) data_basal[t][i][1] = (double)data_basal[t][i][1]/data_counter[t]*clone_norm;
		};
	};
#elif raw_data == 18 || raw_data == 24 || raw_data == 10 || raw_data == 28 || raw_data == 29 || fitbasal == 2         // Normalize total distributions
	
	for (int t = 0; t < timepoints ; t++) 
	{
		clonenumber_exp[t] = data_counter[t];
	};
	for (int t=0; t<timepoints; t++) 
	{
		for (int i=0; i<data_total[t].size(); i++) 
		{
			data_total[t][i][1] = (double)data_total[t][i][1]/data_counter[t]*clone_norm;
		};
		data_basal[t].clear();
		data_basal[t] = data_total[t];         // Attention !!! identification of basal and total cell number, since no discrimination of basal and suprabasal layer possible
		cerr << "data_basal.size() = " << data_basal[t].size() << endl;
		if (data_basal[t].size() != data_total[t].size())		// Debug !!!
		{
			cerr << "error: basal and total don't have same size\n" << endl; 
			exit(1);
		};
	};
#elif plain == 1
	for (int t = 0; t < timepoints ; t++) 
	{
		clonenumber_exp[t] = data_counter[t];
	};
	for (int t=0; t<timepoints; t++) 
	{
		for (int i=start_valid_b-1; i<data_basal[t].size(); i++) 
		{
			data_basal[t][i][1] = (double)data_basal[t][i][1]/data_counter[t]*clone_norm;
		};
	};

#endif
	
	if (data_counter[0] != clonenumber_exp[0]) 
	{
		cerr << "error: clonenumber_exp does not match\n";
		cerr << "clonenumber_exp = " << clonenumber_exp[0] << endl << "clonenumber data = " << data_counter[0] << endl;
#if raw_data != 0 && raw_data !=23
		exit(1);
#endif
	};
	
	//------------------------------------------------------------------------------
	
	//------------------ Make basal CSD -----------------------------------------
	
#if plain != 1 && raw_data != 18 && raw_data != 24 && raw_data != 10 && raw_data != 28 && raw_data != 29
	for (int t=0; t<timepoints; t++) 
	{
		
		data_basal[t].resize(data0[t].size());
		for (int i=0; i<data0[t].size(); i++) 
		{
			data_basal[t][i].resize(2);
			data_basal[t][i][0]=i;
			data_basal[t][i][1]=0;
			for (int j=start_val(i); j<data0[t][i].size(); j++) 
			{
				data_basal[t][i][1] = data_basal[t][i][1] + data0[t][i][j];
			};
		};
	};
	
#endif
	
	//------------------ Make total CSD -----------------------------------------
	
#if plain != 1 && raw_data != 18 && raw_data != 24 && raw_data != 10 && raw_data != 28 && raw_data != 29 && fitbasal != 2
	for (int t=0; t<timepoints; t++) 
	{
		data_total[t].clear();
		data_total[t].resize(data0[t].size()+data0[t][1].size());
		for (int k=0; k<data_total[t].size(); k++) 
		{
			data_total[t][k].resize(2);
			data_total[t][k][0]=k;
			data_total[t][k][1]=0;
		}
		
		for (int i=1; i<data0[t].size(); i++) 
		{
			for (int j=start_val(i); j<data0[t][i].size(); j++) 
			{
				data_total[t][i+j][1] = data_total[t][i+j][1] + data0[t][i][j];
			};
		};
		
		double total_sum=0;
		for (int k=2; k<data_total[t].size(); k++)
		{
			total_sum = total_sum + data_total[t][k][1];
		};
		for (int k=0; k<data_total[t].size(); k++)
		{
			data_total[t][k][1]=(double)data_total[t][k][1]/total_sum*clone_norm;
		};
	};
#endif
	
	vector<vector<vector<double> > > data_basal_print(timepoints);
	vector<vector<vector<double> > > data_total_print(timepoints);
	for (int t=0; t<timepoints; t++)
	  {
          data_basal_print[t] = bin(data_basal[t],binmode,binsize);   // Remove !!! No binning !!!
          data_total_print[t] = bin(data_total[t],binmode,binsize);
	  };

	//------------------ Print data in auxiliary file -----------------------------
	
#if raw_data == 10 || raw_data == 18
	
	vector<vector<double> > lesions_sizes_help = read_file(&file_data_lesions);
	vector<double> lesions_sizes(lesions_sizes_help.size());
	for (int i=0; i<lesions_sizes_help.size(); i++)
	{
		lesions_sizes[i] = lesions_sizes_help[i][0];
	};

	lesion_SD = bin_cont(lesions_sizes,binsize_lesion);
#endif
	
	//--------------- Rescale data ------------------------------------------------
	
	double av_clonesize_basal_data=0; 
	double av_clonesize_total_data=0;
	int av_ctr_basal=0;
	int av_ctr_total=0;
	for (int i=start_valid_b-1; i<data_basal[0].size(); i++) 
	{
		av_clonesize_basal_data = av_clonesize_basal_data + data_basal[0][i][0]*data_basal[0][i][1];
		av_ctr_basal = av_ctr_basal + data_basal[0][i][1];
	};
	for (int i=start_valid_b-1; i<data_total[0].size(); i++) 
	{
		av_clonesize_total_data = av_clonesize_total_data + data_total[0][i][0]*data_total[0][i][1];
		av_ctr_total = av_ctr_total + data_total[0][i][1];
	};
	
	av_clonesize_basal_data = (double)av_clonesize_basal_data/av_ctr_basal;
	av_clonesize_total_data = (double)av_clonesize_total_data/av_ctr_total;
	
//	cerr << "av_clonesize_basal_data = " << av_clonesize_basal_data << endl;
//	exit(2);
	
	for (int i=start_valid_b-1; i < data_basal[0].size(); i++) 
	{
		file_data_basal_resc_out << (double)data_basal[0][i][0]/av_clonesize_basal_data << '\t' << data_basal[0][i][1]*av_clonesize_basal_data << endl;      // Attention !!! Not correct for binning
	};
	for (int i=0; i < data_total[0].size(); i++) 
	{
		file_data_total_resc_out << (double)data_total[0][i][0]/av_clonesize_total_data << '\t' << data_total[0][i][1]*av_clonesize_total_data << endl;      // Attention !!! Not correct for binning
	};
	
	//--------------- End rescale data ------------------------------------------------
	
	if (binmode == 'e')
	{
		file_data_basal_out << data_basal_print[0][0][0] << '\t' << data_basal_print[0][0][1] << endl;
	};
	
	data_counter[0]=0;
	for (int i=(start_valid_b-1)/binsize; i<data_basal_print[0].size(); i++)			   // Attention !!! Adjust starting point, pertinent to binsize
	{
		file_data_basal_out << data_basal_print[0][i][0] << '\t' << data_basal_print[0][i][1] << endl;
		data_counter[0] = data_counter[0] + data_basal[0][i][1];
		file_data_cumul_basal_out << data_basal[0][i][0] << '\t' << ((double)clone_norm - data_counter[0]-excluded_clones)*(clone_norm)/(clone_norm-excluded_clones) << endl;
	}; 
	data_counter[0]=0;
	for (int i=(start_valid_b-1)/binsize; i<data_total_print[0].size(); i++)		// Attention !!! Adjust starting point, pertinent to binsize
	{
		file_data_total_out << data_total_print[0][i][0] << '\t' << data_total_print[0][i][1] << endl;
        cerr << data_total_print[0][i][0] << '\t' << data_total_print[0][i][1] << endl;     //Debug !!!
		data_counter[0] = data_counter[0] + data_total[0][i][1];
		file_data_cumul_total_out << data_total[0][i][0] << '\t' << (clone_norm - data_counter[0]-excluded_clones)*(clone_norm)/(clone_norm-excluded_clones) << endl;
	};

	for (int i=0; i<data0[0].size(); i++) 
	{
		double sum=0;
		double av_ctr=0;
		for (int j=0; j<data0[0][i].size(); j++) 
		{
			av_ctr = av_ctr + data0[0][i][j];
			sum = sum + j*data0[0][i][j];
		};
		if (av_ctr>0) file_data_avsupra_out << i << '\t' << (double)sum/av_ctr << endl;
	};
	
	for (int i=0; i<lesion_SD.size(); i++) 
	{
		file_lesion_SD_out << lesion_SD[i][0] << '\t' << lesion_SD[i][1] << endl;
	};
	//------------------------------------------------------------------------------
};
//-----------------------------------------------------------

//----------------- Prepare initial parameters --------------------------------------
vector<vector<double> > MainWindow::prepare_parameters(vector<vector<double> >& parameters0)
{
#if many_step_process == 1
	reprod_rate0 = (double)reprod_rate0*(reprod_step_number+1);
#endif
	parameters0.clear();
	
	reprod_rate_init[0] = reprod_rate0;
#if Dcell_progenitors == 1
    reprod_rate_init[1] = reprod_rate1;
#endif
	r_init[0] = r0;
	delta_init[0] = delta0;
	strat_rate_init[1] = strat_rate1;
	shed_rate_init[2] = shedrate2;
	direct_strat_frac_init[0] = direct_strat_frac0;
	switchrate_forw_init[0] = switchrate_forw0;
	switchrate_backw_init[1] = switchrate_backw0;
	init2prog_frac_init[0] = init2prog_frac0;
	
	parameters0.push_back(reprod_rate_init);         // 0       // Attention !!! Check order 
	parameters0.push_back(r_init);					// 1
	parameters0.push_back(dediff_rate_init);        // 2
	parameters0.push_back(delta_init);              // 3
	parameters0.push_back(strat_rate_init);         // 4
	parameters0.push_back(shed_rate_init);          // 5
	parameters0.push_back(direct_strat_frac_init);  // 6
	parameters0.push_back(switchrate_forw_init);    // 7
	parameters0.push_back(switchrate_backw_init);   // 8
	parameters0.push_back(init2prog_frac_init);     // 9
	
#if raw_data == 8
	cerr << "error: order of parameters not well-defined"; exit(1);
	vector<double> rho_vector(celltypenumber,0);
	rho_vector[0]=0.5;
	parameters0.push_back(rho_vector);
#endif
	
	for (int i=0; i<parameters0.size(); i++) 
	{
		for (int j=0; j<parameters0[i].size(); j++) 
		{
			cerr << parameters0[i][j] << '\t';
		};
		cerr << endl;
	};

	return parameters0;
};

//-----------------------------------------------------------

//--------------------- Find parameter set with maximum likelihood from files of completed parameter scan ----------------

vector<vector<double> > MainWindow::get_best_params(vector<vector<double> > parameters0)
{
	ifstream file_parlimits("parlimits.dat");
	ofstream file_bestparams("bestparams.dat");
	
	int min_index1 = 0;        // Attention !!! Depends on the fitting intervals. Adjust accordingly
	int max_index1 = 0;
	int index1_intv = 0;
	int min_index2 = 0;
	int max_index2 = 0;
	int index2_intv = 0;
//	int min_index3 = 0;
//	int max_index3 = 0;
//	int index3_intv = 0;
	
	file_parlimits >> min_index1;
	file_parlimits >> max_index1;
	file_parlimits >> index1_intv;
	file_parlimits >> min_index2;
	file_parlimits >> max_index2;
	file_parlimits >> index2_intv;
//	file_parlimits >> min_index3;
//	file_parlimits >> max_index3;
//	file_parlimits >> index3_intv;
	
	cerr <<  min_index1 << '\t' << max_index1 << '\t' <<  index1_intv << '\t' << min_index2 << '\t' << max_index2 << '\t' << index2_intv << endl; // '\t' << min_index3 << '\t' << max_index3 << '\t' << index3_intv << endl;		// Debug !!!
	vector<vector<int> > indices;
	for (int i=min_index1; i<=max_index1; i=i+index1_intv) 
	{
		for (int j=min_index2; j<=max_index2; j=j+index2_intv) 
		{
//			for (int k=min_index3; k<=max_index3; k=k+index3_intv) 
//			{
				vector<int> aux(0);
				aux.push_back(i);
				aux.push_back(j);
//				aux.push_back(k);
				indices.push_back(aux);
//			};
		};
	};
	
	string aux_str = "loglikelies";
	vector<ifstream*> loglikelies_in = make_filelist(aux_str,indices);
	vector<double> best_parameters = eval_loglikelyfiles(loglikelies_in,(double)index1_intv/100,reprodrate_intv,(double)index2_intv/1000,stratrate_intv,shedrate_intv);
	//vector<double> best_parameters = eval_loglikelyfiles(loglikelies_in,0,reprodrate_intv,0.01,stratrate_intv,shedrate_intv);
	cerr << "best parameters = " << '\t' << best_parameters[0] << '\t' << best_parameters[1] << '\t' << best_parameters[2] << '\t' << best_parameters[3] << '\t' << best_parameters[4] << endl;
	
	parameters0[0][0]=best_parameters[1];   // reprod rate
	parameters0[1][0]=best_parameters[2];   // r
	parameters0[3][0]=best_parameters[0];   // delta
	parameters0[4][1]=best_parameters[3];   // strat rate
	parameters0[5][2]=best_parameters[4];   // shed rate
	
	for (int i=0; i<parameters0.size(); i++) 
	{
		for (int j=0; j<parameters0[i].size(); j++) 
		{
			file_bestparams << parameters0[i][j] << '\t';
		};
		file_bestparams << endl;
	};
	
	file_bestparams.close();
	exit(0);
	return parameters0;
};


//-------------------------------------------------------------------		
//-------------------------------------------------------------------
// Functions for fitting the data
//-------------------------------------------------------------------
//-------------------------------------------------------------------

//------------------ Scan parameters for likelihood -------------------------------------
void MainWindow::scan_likelihood(vector<vector<double> > param_init, const vector<vector<vector<double> > >& data0)
{	
	*file_loglikelies << delta0 << '\t' << reprod_rate0 << '\t' << r0 << '\t' << 0.0 <<  endl;
	for (int i=0; i < (reprodratelimit_high - reprodratelimit_low)/reprodrate_intv; i++) 
	{
#if scan_rho == 1
		for (int j=0; j<(rholimit_high - rholimit_low)/rho_intv; j++)
		{
#else
		for (int j=0; j<(stratratelimit_high - stratratelimit_low)/stratrate_intv; j++)
		{
#endif
			for (int k=0; k<(shedratelimit_high - shedratelimit_low)/shedrate_intv; k++) 
			{
#if fit_init2prog_frac == 1
				for (int l=0; l<(init2proglimit_high - init2proglimit_low)/init2prog_intv; l++) 
				{
#endif
					strat_infinity_flag = 0;
					
					reprod_rate_init[0] = reprodratelimit_low + i*reprodrate_intv;
					param_init[0][0] = reprod_rate_init[0];			// Check this !!!
					
#if scan_rho == 1
					if (rholimit_low + j*rho_intv <= 0.99)
					{
						strat_rate_init[1] = (rholimit_low + j*rho_intv)/(1-(rholimit_low + j*rho_intv))*reprod_rate_init[0];
						param_init[4][1] = strat_rate_init[1];			// Check this !!!
					}
					else 
					{
						strat_infinity_flag = 1;
						strat_rate_init[1] = 10.0;
						param_init[4][1] = strat_rate_init[1];	
						
					}
					
					
#else
					strat_rate_init[1] = stratratelimit_low + j*stratrate_intv;
					param_init[4][1] = strat_rate_init[1];			// Check this !!!
#if fit_directstrat == 1
					direct_strat_frac_init[0] = 1 - (double)strat_rate_init[1]/strat_rate_eff; 	
					param_init[6][0] = direct_strat_frac_init[0];
#endif
#endif
					
					shed_rate_init[2] = shedratelimit_low + k*shedrate_intv;
					param_init[5][2] = shed_rate_init[2];
					
#if shedrate_fixed == 1
#if fit_directstrat == 1
					shed_rate_init[2] = reprod_rate_init[0]*strat_rate_eff/(reprod_rate_init[0]+strat_rate_eff)/SB_B_ratio;
					param_init[5][2] = shed_rate_init[2];
#else
					shed_rate_init[2] = reprod_rate_init[0]*strat_rate_init[1]/(reprod_rate_init[0]+strat_rate_init[1])/SB_B_ratio;
					param_init[5][2] = shed_rate_init[2];
#endif
#endif
					
#if fit_init2prog_frac == 1
					init2prog_frac_init[0] = init2proglimit_low + l*init2prog_intv;
					param_init[9][0] = init2prog_frac_init[0];
#endif
					
					run* run_temp = new run(param_init);
					
					run_temp->go();
#if fitbasal == 1
#if raw_data == 10 || raw_data == 18 || raw_data == 24 || raw_data == 28 || raw_data == 29
					double loglikely_curr = run_temp->loglikely_basal(data_total);
#else
					double loglikely_curr = run_temp->loglikely_basal(data_basal);
#endif
#elif fitbasal == 2
					double loglikely_curr = run_temp->loglikely_basal(data_total);
#elif fitbasal == 3
					double loglikely_curr = run_temp->loglikely_total(data_total);
#else 
					double loglikely_curr = run_temp->loglikely(data0);
#endif
					
#if scan_rho == 1
					*file_loglikelies << run_temp->reprod_rate[0] << '\t' << run_temp->strat_rate[1]/(run_temp->strat_rate[1] + run_temp->reprod_rate[0]) << '\t' << run_temp->shed_rate[2] << '\t' << loglikely_curr << endl;
#elif fit_directstrat == 1
					*file_loglikelies << run_temp->reprod_rate[0] << '\t' << run_temp->strat_rate[1] << '\t' << run_temp->direct_strat_frac[0] << '\t' << loglikely_curr << endl;
#elif fit_init2prog_frac== 1
					*file_loglikelies << run_temp->reprod_rate[0] << '\t' << run_temp->strat_rate[1] << '\t' << run_temp->init2prog_frac[0] << '\t' << loglikely_curr << endl;
#else
					*file_loglikelies << run_temp->reprod_rate[0] << '\t' << run_temp->strat_rate[1] << '\t' << run_temp->shed_rate[2] << '\t' << loglikely_curr << endl;
#endif
					
					delete run_temp;
#if fit_init2prog_frac == 1
				};
#endif
			};
		};
	};
	
};

//------------------ Fit with genetic algorithm -------------------------------------
run* MainWindow::fit_evol(vector<vector<double> > param_init, const vector<vector<vector<double> > >& data0)
{
	
	run_old.clear();
	
	toolarge_ctr=0;
	for (int esmblctr=0; esmblctr<esmblnumber; esmblctr++) 
	{
		cerr << "esmblctr = " << esmblctr << endl;

		run* run_temp = new run(param_init);
		run_old.push_back(run_temp);
		
		run_temp->log_av_flag = 1;
		run_temp->go();
		
#if raw_data !=21 && raw_data != 26
		while (run_temp->clonetobig_flag == 1) 
		{
			fail_ctr++;
//			run_temp->setup();
//			run_temp->strat_rate[1]= min(30.0,run_temp->strat_rate[1]+10.0);
			run_temp->go();
		};
#endif
		if (run_temp->logflag==0)                   // Debug !!!
		{
			cerr << "error: data not logged, since too many failed attempts\n";
			exit(1);
		};
		
//		for (int i=0; i<run_temp->joint_distr_norm.size(); i++)                // Debug !!!   
//		{
//			for (int j=0; j<run_temp->joint_distr_norm[i].size(); j++) 
//			{
//				cerr << run_temp->joint_distr_norm[i][j] << '\t';
//			};
//			cerr << endl;
//		};
		
#if raw_data != 20 && raw_data != 21 && raw_data != 26		
#if fitbasal == 1
		file_loglikely << -1 << '\t' << esmblctr << '\t' << run_temp->loglikely_basal(data_basal) << endl;
#else 
		file_loglikely << -1 << '\t' << esmblctr << '\t' << run_temp->loglikely(data0) << endl;
#endif
#endif
		
		for (int runctr=0; runctr<runnumber; runctr++)
		{
			run_temp = run_old[esmblctr];
			cerr << "runctr at start = " << runctr << endl << endl;           // Debug !!!
			T = (double)T0/(Tdrop1*runctr+1)*(1+T_max/T0*exp(-(double)runctr/Tdropruns)); //*exp(-1*runctr/anneal);
			run* run_new = run_temp->replicate();
			
			cerr << "cell number after run_new initialisation = " << run_new->cells.size() << endl;
			cerr << "run_temp->rho = " << run_temp->rho << endl;
				cerr << "run_old[i]->rho = " << run_old[esmblctr]->rho << endl;
//			cerr << "parameters of old run: \n reprod rate = " << run_temp->reprod_rate[0] << endl << "strat rate = " << run_temp->strat_rate[1] << endl << "shed rate = " << run_temp->shed_rate[2] << endl;      // Debug !!!
//			cerr << "parameters of new run: \n reprod rate = " << run_new->reprod_rate[0] << endl << "strat rate = " << run_new->strat_rate[1] << endl << "shed rate = " << run_new->shed_rate[2] << endl;
//			cerr << "cellnumber of new run = " << run_new->cells.size() << endl;
//			cerr << "clonenumber of new run = " << run_new->clones.size() << endl;
//			cerr << "cell type of last cell = " << run_new->cells.back()->get_celltype() << endl << endl;
			run_new->variate_params(var_vector);                
			
			cerr << "reprodrate after variation = " << run_new->reprod_rate[0] << endl;      // Debug !!!
			cerr << "cellnumber = " << run_new->cells.size() << endl;            // Debug !!!
			cerr << "clonenumber = " << run_new->clones.size() << endl << endl;
			
			run_new->go();
			
			cerr << "runctr = " << runctr << endl << "esmblctr = " << esmblctr << endl;
			while (run_new->clonetobig_flag == 1) 
			{
				fail_ctr++;
				run_new->setup();
				run_new->variate_params(var_vector);             // Attention !!! Don't use when parameters should be fixed
				run_new->go();
			};
			cerr << "joint_distr_norm : \n";
//			for (int i=0; i<run_new->joint_distr_norm.size(); i++)                // Debug !!!   
//			{
//				for (int j=0; j<run_new->joint_distr_norm[i].size(); j++) 
//				{
//					cerr << run_new->joint_distr_norm[i][j] << '\t';
//				};
//				cerr << endl;
//			};
			cerr << endl;
			cerr << "number of clones = " << run_temp->clones.size() << endl;
			cerr << "sizes of clones: \n";
			cerr << endl;

			cerr << "runctr at end = " << runctr << endl << endl;           // Debug !!!
			
			//------- new run completed, now compare likelihood --------------------
			
#if raw_data != 20 && raw_data !=21 && raw_data != 26
			if (run_new->completed == 1 && run_temp->completed == 1) 
			{
				
#if fitbasal == 1
				double loglikely_old = run_temp->loglikely_basal(data0);
				double loglikely_new = run_new->loglikely_basal(data0);
#else
				double loglikely_old = run_temp->loglikely(data0);
				double loglikely_new = run_new->loglikely(data0);
#endif
			
				file_temp << esmblctr << '\t' << runctr << '\t' << T << '\t' << T*stepsize_strat << '\t' << exp(1.0/T*(loglikely_new - loglikely_old)) << endl;
				
				cerr << "loglikely_old = "<< loglikely_old << endl << "loglikely_new = " << loglikely_new << endl;
				
				double rand1 = drand48();
				if (rand1 < exp(1.0/T*(loglikely_new - loglikely_old))) 
				{
					cerr << "new run survived\n" << endl;        // Debug !!!
					delete run_temp;
					run_old[esmblctr] = run_new;
					file_loglikely << esmblctr << '\t' << loglikely_new << endl;
				}
				else 
				{
					cerr << "old run survived\n" << endl;
					run_old[esmblctr] = run_temp;
					file_loglikely << esmblctr << '\t' << loglikely_old << endl;
					delete run_new;
				};
				
			}
			else if (run_temp->completed==0)
			{
				cerr << "new run survived\n" << endl;        // Debug !!!
				run_old[esmblctr] = run_new;
				delete run_temp;

			}
			else
			{
				toolarge_ctr++;
				run_old[esmblctr] = run_temp;
				delete run_new;
			};

#endif
		

			
			cerr << "failctr = " << fail_ctr << endl;
			file_parameters << run_old[esmblctr]->reprod_rate[0] << '\t' << run_old[esmblctr]->r[0] << '\t' << run_old[esmblctr]->delta[0] << '\t' << run_old[esmblctr]->strat_rate[1] << '\t' << run_old[esmblctr]->shed_rate[2] << '\t' << run_old[esmblctr]->delta[0];
#if init_cells != 1
			file_parameters << '\t' << run_old[esmblctr]->init_clone[0][0] << '\t' << run_old[esmblctr]->init_clone[1][0] << endl;
#else
			file_parameters << endl;
#endif
			run_old[esmblctr]->print_averages(runctr,esmblctr);
		};
		file_temp << endl;
		file_av << endl;
		file_parameters << endl << run_old[esmblctr]->reprod_rate[0] << '\t' << run_old[esmblctr]->r[0] << '\t' << run_old[esmblctr]->delta[0] << '\t' << run_old[esmblctr]->strat_rate[1] << '\t' << run_old[esmblctr]->shed_rate[2] << '\t' << run_old[esmblctr]->delta[0] << endl << endl;
		run_old[esmblctr]->print_averages(runnumber-1,esmblctr);
		file_av << endl;

	};
	
	//---------- End Ensemble Loop -------------------------------------------------------------------
	
		//---------- Select optimal fit of all run sequences -------------------------------------------------------------------
	
	run* maxrun = run_old[0];
	for (int i=0; i<esmblnumber-1 && run_old[i]->completed==0; i++) 
	{
		maxrun = run_old[i+1];
		cerr << "maxrun " << i+1 << " completed? " << maxrun->completed << endl << endl;		// Debug !!!
	};
	if (maxrun->completed==0) {cerr << "error: no valid run found\n"; exit(1);}
	
#if raw_data != 20 && raw_data != 21 && raw_data != 26
	for (int i=0; i<esmblnumber; i++) 
	{
		if (run_old[i]->completed==1 && maxrun->completed==1) 
		{
		
#if fitbasal == 1
			if (run_old[i]->loglikely_basal(data0) > maxrun->loglikely_basal(data0)) 
			{
				maxrun = run_old[i];
			};
#else
			if (run_old[i]->loglikely(data0) > maxrun->loglikely(data0))
			{
				maxrun = run_old[i];
			};
#endif
			file_parameters << run_old[i]->reprod_rate[0] << '\t' << run_old[i]->r[0] << '\t' << run_old[i]->delta[0] << '\t' << run_old[i]->strat_rate[1] << endl;
		};
	};
#endif	

	file_parameters << endl << endl;
	file_parameters << maxrun->reprod_rate[0] << '\t' << maxrun->r[0] << '\t' << maxrun->delta[0] << '\t' << maxrun->strat_rate[1] << endl;
	
	return maxrun;
};

//-------------------------------------------------------------------
// End Fitting
//-------------------------------------------------------------------



//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Run best fit parameter set and find Errorbars
//-------------------------------------------------------------------
//-------------------------------------------------------------------
	
void MainWindow::find_errorbars(run* maxrun, const vector<vector<vector<double > > >& data0)
{

	//-------------- Get reference CSD --------------------------
	N_max=maxCS_av*clonenumber_bestfit;
	N_min=clonenumber_bestfit/3;
	run* run_bestfit= maxrun->replicate(clonenumber_bestfit);
	cerr << run_bestfit->reprod_rate[0] << '\t' << run_bestfit->r[0]<< '\t' << run_bestfit->strat_rate[1] << endl;
//	exit(2);

	run_bestfit->go();
//    for (int i=0; i<run_bestfit->cells.size(); i++)                          // Debug !!!
//    {
//        drand48();
//        if (run_bestfit->cells[i]->get_celltype() == 1 /*&& drand48() <= strat_rate[1]/(reprod_rate[0]+strat_rate[1])+0.15*/)
//        {
//            cerr << " error: cell with celltype 1\n";
//            exit(1);
//        };
//    };
	joint_distr_av = run_bestfit->joint_distr;
	CSD_basal_av = run_bestfit->clonesizedistr_basal;
	CSD_total_av = run_bestfit->clonesizedistr;
	cumuldistr_basal_av = run_bestfit->cumuldistr_basal;
	cumuldistr_total_av = run_bestfit->cumuldistr_total;
	
	cerr << "check cumuldistr_basal_av (in best fit)): first element = " << run_bestfit->cumuldistr_basal.back() << endl;               // Debug !!!
//	cerr << "check cumuldistr_basal_av: first element = " << cumuldistr_basal_av.back() << endl;               // Debug !!!
	
	floatingclone_frac_av = run_bestfit->floating_clone_frac;
	
	vector<vector<double> > joint_cumul_distr_av(joint_distr_av.size());
	
	print_data(run_bestfit);
	
	CSD_basal_av = bin(CSD_basal_av,binmode,binsize);
	CSD_total_av = bin(CSD_total_av,binmode,binsize);
	//-------------- Collect ensembles of sims, size correspionding to data, for errorbars and p-values -----------------
	
	int dataset=timepoints-1;
	
	int clonenumber_simexp = (1/*+(double)run_bestfit->clonesizedistr_basal[0]/clone_norm*/)*clonenumber_exp[0];
#if raw_data == 0 || raw_data == 23 || raw_data == 26
	clonenumber_simexp = clonenumber_simexp + run_bestfit->clonesizedistr_basal[1];
	clonenumber_simexp = clonenumber_simexp + run_bestfit->clonesizedistr_basal[2];
#endif
	N_max=maxCS_av*clonenumber_simexp;
	N_min=clonenumber_simexp/3;
	
	vector<vector<vector<double> > > joint_distr_esmbl(runnumber2);
	vector<vector<vector<double> > > joint_cumul_distr_esmbl(runnumber2);
	vector<vector<double> > CSD_basal_esmbl(runnumber2);
	vector<vector<double> > CSD_total_esmbl(runnumber2);
	vector<vector<double> > cumuldistr_basal_esmbl(runnumber2);
	vector<vector<double> > cumuldistr_total_esmbl(runnumber2);
	
	vector<double> floatclone_freq(runnumber2);

	
	for (int runctr=0; runctr<runnumber2; runctr++) 
	{
		cerr << "clonenumber_simex = " << clonenumber_simexp << endl;
		cerr << "best fit clonenumber = " << run_bestfit->clonenumber << endl;
		run* run_indv = run_bestfit->replicate(clonenumber_simexp);
		run_indv->setup();
		run_indv->go();
		joint_distr_esmbl[runctr] = run_indv->joint_distr;
		CSD_basal_esmbl[runctr] = run_indv->clonesizedistr_basal;
//		CSD_basal_esmbl[runctr].resize(data0[dataset].size()+1,0);
		CSD_total_esmbl[runctr] = run_indv->clonesizedistr;
		CSD_basal_esmbl[runctr] = bin(CSD_basal_esmbl[runctr],binmode,binsize);

		cumuldistr_basal_esmbl[runctr] = run_indv->cumuldistr_basal; 
		cumuldistr_total_esmbl[runctr] = run_indv->cumuldistr_total; 
		
		double norm_joint=0;
		joint_cumul_distr_esmbl[runctr].resize(joint_distr_esmbl[runctr].size());
		for (int i=1; i<joint_distr_esmbl[runctr].size(); i++) 
		{
			joint_cumul_distr_esmbl[runctr][i].resize(joint_distr_esmbl[runctr][i].size());
			for (int j=0; j<joint_distr_esmbl[runctr][i].size(); j++) 
			{
				norm_joint = norm_joint + joint_distr_esmbl[runctr][i][j];
				joint_cumul_distr_esmbl[runctr][i][j] = norm_joint;
			};
		};
		
		floatclone_freq[runctr]=run_indv->floating_clone_frac;
		
		delete run_indv;
		
		
	};
	
	//-------------- Ensemble created; now evaluate error bars and p-values -------------------
	
	// reorder
	
	vector<vector<vector<double> > > clonefreq_plus(data0[dataset].size());
	vector<vector<vector<double> > > clonefreq_minus(data0[dataset].size());
	
	vector<vector<double> > errorbar_plus(data0[dataset].size());
	vector<vector<double> > errorbar_minus(data0[dataset].size());
	
	vector<vector<double> > clonefreq_basal_plus(data0[dataset].size());
	vector<vector<double> > clonefreq_basal_minus(data0[dataset].size());
	vector<double> errorbar_basal_plus(data0[dataset].size());
	vector<double> errorbar_basal_minus(data0[dataset].size());
	
	vector<vector<double> > p_val(data0[dataset].size());
	
	vector<vector<vector<double> > > joint_distr_esmbl_sort(data0[dataset].size()); 
	vector<vector<double> > CSD_basal_esmbl_sort(data0[dataset].size());
	
#if raw_data != 10 && raw_data != 18 && raw_data != 24 && raw_data != 23 && raw_data != 26 && raw_data != 28 && raw_data != 29 && !(raw_data == 0 && plain == 1)
	
	for (int i=0; i<data0[dataset].size(); i++) 
	{
		clonefreq_plus[i].resize(data0[dataset][i].size());
		clonefreq_minus[i].resize(data0[dataset][i].size());
		errorbar_plus[i].resize(data0[dataset][i].size());
		errorbar_minus[i].resize(data0[dataset][i].size());
		p_val[i].resize(data0[dataset][i].size());
		joint_distr_esmbl_sort[i].resize(data0[dataset][i].size());
		
		for (int j=start_val(i); j<data0[dataset][i].size(); j++) 
		{
			clonefreq_plus[i][j].clear();
			clonefreq_plus[i][j].reserve(runnumber2+1);
			clonefreq_minus[i][j].clear();
			clonefreq_minus[i][j].reserve(runnumber2+1);
			joint_distr_esmbl_sort[i][j].clear();
			joint_distr_esmbl_sort[i][j].reserve(runnumber2+1);
			p_val[i][j]=0;
			
			
			for (int runctr=0; runctr<runnumber2; runctr++) 
			{
				joint_distr_esmbl_sort[i][j].push_back(joint_distr_esmbl[runctr][i][j]);
				
				if (joint_distr_esmbl[runctr][i][j] > joint_distr_av[i][j]) 
				{
					clonefreq_plus[i][j].push_back(joint_distr_esmbl[runctr][i][j]);
				}
				else if (joint_distr_esmbl[runctr][i][j] < joint_distr_av[i][j])
				{
					clonefreq_minus[i][j].push_back(joint_distr_esmbl[runctr][i][j]);
				}
				else 
				{
					if (drand48() < 0.5) clonefreq_plus[i][j].push_back(joint_distr_esmbl[runctr][i][j]);
					else clonefreq_minus[i][j].push_back(joint_distr_esmbl[runctr][i][j]);
				};

				if (data0[dataset][i][j] > joint_distr_av[i][j] && joint_distr_esmbl[runctr][i][j] >= data0[dataset][i][j]) 
				{
					p_val[i][j] = p_val[i][j] + 1.0;
				}
				else if (data0[dataset][i][j] < joint_distr_av[i][j] && joint_distr_esmbl[runctr][i][j] <= data0[dataset][i][j])
				{
					p_val[i][j] = p_val[i][j] + 1.0;
				};
				
			};
			

			
			p_val[i][j] = (double)p_val[i][j]/runnumber2;
			
			i_curr=i;
			j_curr=j;
			sort(clonefreq_plus[i][j].begin(),clonefreq_plus[i][j].end());
			sort(clonefreq_minus[i][j].begin(),clonefreq_minus[i][j].end());
			sort(joint_distr_esmbl_sort[i][j].begin(),joint_distr_esmbl_sort[i][j].end(),compare_joint);
			

			joint_distr_esmbl_sort[i][j].resize(floor(siglevel*(double)joint_distr_esmbl_sort[i][j].size()));
			
			errorbar_plus[i][j] = *max_element(joint_distr_esmbl_sort[i][j].begin(),joint_distr_esmbl_sort[i][j].end());
			errorbar_minus[i][j] = *min_element(joint_distr_esmbl_sort[i][j].begin(),joint_distr_esmbl_sort[i][j].end());
			
			if (siglevel>1 || siglevel<0) 
			{
				cerr << "error: siglevel out of bounds\n";
				exit(1);
			};
			
		};
		

	};
	
#endif
	
	clonefreq_basal_plus.resize(run_bestfit->clonesizedistr_basal.size());
	clonefreq_basal_minus.resize(run_bestfit->clonesizedistr_basal.size());
	
		//------- Compute Basal errorbars ---------------------------------------
	int maxclone=0;
	int clonesize=0;
	for (int runctr=0; runctr<runnumber2; runctr++) 
	{
		clonesize=CSD_basal_esmbl[runctr].size();
		if (clonesize > maxclone) maxclone = clonesize;
	};
	clonesize = CSD_basal_av.size();
	if (clonesize > maxclone) maxclone = clonesize;
	CSD_basal_av.resize(maxclone+1,0);
	
	double last_element=cumuldistr_basal_av.back();
	cumuldistr_basal_av.resize(maxclone+1,last_element);
	
	for (int runctr=0; runctr<runnumber2; runctr++) 
	{
		CSD_basal_esmbl[runctr].resize(maxclone+1,0);
#if raw_data != 23 && raw_data != 26
		last_element=cumuldistr_basal_esmbl[runctr].back();
#endif
		cumuldistr_basal_esmbl[runctr].resize(maxclone+1,last_element);
	};
	
	CSD_basal_esmbl_sort.resize(maxclone+1);
	errorbar_basal_plus.resize(maxclone+1);
	errorbar_basal_minus.resize(maxclone+1);
	vector<vector<double> > cumuldistr_basal_esmbl_sort(maxclone+1);
	vector<double> errorbar_cumul_basal_plus(maxclone+1);
	vector<double> errorbar_cumul_basal_minus(maxclone+1);
	
	for (int i=0; i<CSD_basal_esmbl_sort.size(); i++) 
	{
		CSD_basal_esmbl_sort[i].clear();
		CSD_basal_esmbl_sort[i].reserve(runnumber2);
		cumuldistr_basal_esmbl_sort[i].clear();
		cumuldistr_basal_esmbl_sort[i].reserve(runnumber2);
		
		for (int runctr=0; runctr<runnumber2; runctr++) 
		{
			CSD_basal_esmbl_sort[i].push_back(CSD_basal_esmbl[runctr][i]);
			cumuldistr_basal_esmbl_sort[i].push_back(cumuldistr_basal_esmbl[runctr][i]);
		};
		i_curr=i;
		
#if errorbars == 1 || errorbars == 2
		vector<double> stddev_CSD_basal(maxclone+1,0);
		vector<int> error_plus_ctr(maxclone+1,0);
		vector<int> error_minus_ctr(maxclone+1,0);
        
		for (int runctr=0; runctr<runnumber2; runctr++)
		{
			if (CSD_basal_esmbl_sort[i][runctr] > CSD_basal_av[i]) 
			{
				errorbar_basal_plus[i] = errorbar_basal_plus[i] + (CSD_basal_esmbl_sort[i][runctr] - CSD_basal_av[i]);
				error_plus_ctr[i]++;
			}
			else if (CSD_basal_esmbl_sort[i][runctr] < CSD_basal_av[i])
			{
				errorbar_basal_minus[i] = errorbar_basal_minus[i] - (CSD_basal_esmbl_sort[i][runctr] - CSD_basal_av[i]);
				error_minus_ctr[i]++;
			}
			else 
			{
				error_plus_ctr[i]++;
				error_minus_ctr[i]++;
			}

			stddev_CSD_basal[i] = stddev_CSD_basal[i] + (CSD_basal_esmbl_sort[i][runctr] - CSD_basal_av[i])*(CSD_basal_esmbl_sort[i][runctr] - CSD_basal_av[i]);
		};
		stddev_CSD_basal[i] = sqrt((double)stddev_CSD_basal[i]/runnumber2);
#if errorbars == 1
        errorbar_basal_plus[i] = CSD_basal_av[i] + stddev_CSD_basal[i];
        errorbar_basal_minus[i] = CSD_basal_av[i] - stddev_CSD_basal[i];
        errorbar_cumul_basal_plus[i] = *max_element(cumuldistr_basal_esmbl_sort[i].begin(),cumuldistr_basal_esmbl_sort[i].end());
        errorbar_cumul_basal_minus[i] = *min_element(cumuldistr_basal_esmbl_sort[i].begin(),cumuldistr_basal_esmbl_sort[i].end());
#elif errorbars == 2
		errorbar_basal_plus[i] = CSD_basal_av[i] + (double)errorbar_basal_plus[i]/error_plus_ctr[i];
		errorbar_basal_minus[i] = CSD_basal_av[i] - (double)errorbar_basal_minus[i]/error_minus_ctr[i];
#endif
        
#else
		
		sort(CSD_basal_esmbl_sort[i].begin(),CSD_basal_esmbl_sort[i].end(),compare_basal);
		sort(cumuldistr_basal_esmbl_sort[i].begin(),cumuldistr_basal_esmbl_sort[i].end(),compare_cumul_basal);
		
		CSD_basal_esmbl_sort[i].resize(floor(siglevel*(double)CSD_basal_esmbl_sort[i].size()));
		cumuldistr_basal_esmbl_sort[i].resize(floor(siglevel*(double)cumuldistr_basal_esmbl_sort[i].size()));

		errorbar_basal_plus[i] = *max_element(CSD_basal_esmbl_sort[i].begin(),CSD_basal_esmbl_sort[i].end());
		errorbar_basal_minus[i] = *min_element(CSD_basal_esmbl_sort[i].begin(),CSD_basal_esmbl_sort[i].end());
		errorbar_cumul_basal_plus[i] = *max_element(cumuldistr_basal_esmbl_sort[i].begin(),cumuldistr_basal_esmbl_sort[i].end());
		errorbar_cumul_basal_minus[i] = *min_element(cumuldistr_basal_esmbl_sort[i].begin(),cumuldistr_basal_esmbl_sort[i].end());
		cerr << "errorbar_basal_minus[" << i << "] = " << errorbar_basal_minus[i] << endl;
		cerr << "errorbar_basal_plus[" << i << "] = " << errorbar_basal_plus[i] << endl;
		cerr << "CSD_basal_av[" << i << "] = " << CSD_basal_av[i] << endl;
        
        
#endif
	};
	
	
	
	maxclone=0;
	clonesize=0;
	for (int runctr=0; runctr<runnumber2; runctr++) 
	{
		clonesize=CSD_total_esmbl[runctr].size();
		if (clonesize > maxclone) maxclone = clonesize;
	};
	clonesize = CSD_total_av.size();
	if (clonesize > maxclone) maxclone = clonesize;
	CSD_total_av.resize(maxclone+1,0);
	
	for (int runctr=0; runctr<runnumber2; runctr++) 
	{
		CSD_total_esmbl[runctr].resize(maxclone+1,0);
		cerr << "cumuldistr_total_esmbl.size() = " << cumuldistr_total_esmbl.size() << endl;           // Debug
		last_element=cumuldistr_total_esmbl[runctr].back();
		cumuldistr_total_esmbl[runctr].resize(maxclone+1,last_element);
	};
	vector<vector<double> > CSD_total_esmbl_sort(maxclone+1);
	vector<vector<double> > cumuldistr_total_esmbl_sort(maxclone+1);
	
	vector<double> errorbar_total_plus(maxclone+1);
	vector<double> errorbar_total_minus(maxclone+1);
	vector<double> errorbar_cumul_total_plus(maxclone+1);
	vector<double> errorbar_cumul_total_minus(maxclone+1);
	
	for (int i=0; i<CSD_total_esmbl_sort.size(); i++) 
	{
		for (int runctr=0; runctr<runnumber2; runctr++) 
		{
			CSD_total_esmbl_sort[i].push_back(CSD_total_esmbl[runctr][i]);
			cumuldistr_total_esmbl_sort[i].push_back(cumuldistr_total_esmbl[runctr][i]);
		};
		i_curr=i;
        
#if errorbars == 1 || errorbars == 2
        
        vector<double> stddev_CSD_total(maxclone+1,0);
        
        vector<int> error_minus_ctr(maxclone+1,0);
        vector<int> error_plus_ctr(maxclone+1,0);
        
        for (int runctr=0; runctr<runnumber2; runctr++)
        {
            if (CSD_total_esmbl_sort[i][runctr] > CSD_total_av[i])
            {
                errorbar_total_plus[i] = errorbar_total_plus[i] + (CSD_total_esmbl_sort[i][runctr] - CSD_total_av[i]);
                error_plus_ctr[i]++;
            }
            else if (CSD_total_esmbl_sort[i][runctr] < CSD_total_av[i])
            {
                errorbar_total_minus[i] = errorbar_total_minus[i] - (CSD_total_esmbl_sort[i][runctr] - CSD_total_av[i]);
                error_minus_ctr[i]++;
            }
            else
            {
                error_plus_ctr[i]++;
                error_minus_ctr[i]++;
            }
            
            stddev_CSD_total[i] = stddev_CSD_total[i] + (CSD_total_esmbl_sort[i][runctr] - CSD_total_av[i])*(CSD_total_esmbl_sort[i][runctr] - CSD_total_av[i]);
        };
        stddev_CSD_total[i] = sqrt((double)stddev_CSD_total[i]/runnumber2);
#if errorbars == 1
        errorbar_total_plus[i] = CSD_total_av[i] + stddev_CSD_total[i];
        errorbar_total_minus[i] = CSD_total_av[i] - stddev_CSD_total[i];
        errorbar_cumul_total_plus[i] = *max_element(cumuldistr_total_esmbl_sort[i].begin(),cumuldistr_total_esmbl_sort[i].end());
        errorbar_cumul_total_minus[i] = *min_element(cumuldistr_total_esmbl_sort[i].begin(),cumuldistr_total_esmbl_sort[i].end());
#elif errorbars == 2
        errorbar_total_plus[i] = CSD_total_av[i] + (double)errorbar_total_plus[i]/error_plus_ctr[i];
        errorbar_total_minus[i] = CSD_total_av[i] - (double)errorbar_total_minus[i]/error_minus_ctr[i];
#endif
        
#else
        
		sort(CSD_total_esmbl_sort[i].begin(),CSD_total_esmbl_sort[i].end(),compare_total);
		sort(cumuldistr_total_esmbl_sort[i].begin(),cumuldistr_total_esmbl_sort[i].end(),compare_cumul_total);

		CSD_total_esmbl_sort[i].resize(floor(siglevel*(double)CSD_total_esmbl_sort[i].size()));
		cumuldistr_total_esmbl_sort[i].resize(floor(siglevel*(double)cumuldistr_total_esmbl_sort[i].size()));

		errorbar_total_plus[i] = *max_element(CSD_total_esmbl_sort[i].begin(),CSD_total_esmbl_sort[i].end());
		errorbar_total_minus[i] = *min_element(CSD_total_esmbl_sort[i].begin(),CSD_total_esmbl_sort[i].end());
		errorbar_cumul_total_plus[i] = *max_element(cumuldistr_total_esmbl_sort[i].begin(),cumuldistr_total_esmbl_sort[i].end());
		errorbar_cumul_total_minus[i] = *min_element(cumuldistr_total_esmbl_sort[i].begin(),cumuldistr_total_esmbl_sort[i].end());
		cerr << "errorbar_total_minus[" << i << "] = " << errorbar_total_minus[i] << endl;
		cerr << "errorbar_total_plus[" << i << "] = " << errorbar_total_plus[i] << endl;
		cerr << "CSD_total_av[" << i << "] = " << CSD_total_av[i] << endl;
        
#endif
        
	};
	
	// Print new data
	
#if raw_data != 10 && raw_data != 18 && raw_data != 24 && raw_data != 23 && raw_data != 26 && raw_data != 28 && raw_data != 29 && !(raw_data == 0 && plain == 1)
	
	for (int i=0; i<data0[dataset].size(); i++) 
	{
		for (int j=start_val(i); j<data0[dataset][i].size(); j++) 
		{
//			file_pval << i << '\t' << j << '\t' << p_val[i][j] << endl;
			if (i<fileerrorbar_nr) 
			{
				(*file_errorbars[i]) << j << '\t' << joint_distr_av[i][j] << '\t' << errorbar_minus[i][j] << '\t' << errorbar_plus[i][j] << endl;
				(*file_errorbars_xmgr[i]) << j << '\t' << joint_distr_av[i][j] << '\t' << errorbar_plus[i][j] - joint_distr_av[i][j] << '\t' <<  joint_distr_av[i][j] - errorbar_minus[i][j] << endl;
			};
		};

	};
#endif
	
	bins.clear();
	int maxsize = max(CSD_basal_av.size(),CSD_total_av.size()); 
	if (binmode == 'l')
	{
	    for (int i=0; binsize*i<maxsize; i++)
		{
			bins.push_back(binsize*i);
		};
	}
	else if (binmode == 'e')
	{
	    for (int i=0; powl(binsize,i)<maxsize; i++)
		{
			bins.push_back(powl(binsize,i));
		};
	}
	
	double data_counter_total=0;
	double data_counter_basal=0;
	double av_clonesize_basal=0; 
	double av_clonesize_total=0;
	int av_ctr_basal=0;
	int av_ctr_total=0;
	
	for (int i=(start_valid_b-1)/binsize; i<CSD_basal_av.size(); i++) 
	{
		av_clonesize_basal = av_clonesize_basal + i*CSD_basal_av[i];
		av_ctr_basal = av_ctr_basal + CSD_basal_av[i];
	};
	for (int i=1; i<CSD_basal_av.size(); i++) 
	{
		av_clonesize_total = av_clonesize_total + i*CSD_basal_av[i];
		av_ctr_total = av_ctr_total + CSD_basal_av[i];
	};
	
	av_clonesize_basal = (double)av_clonesize_basal/av_ctr_basal;
	av_clonesize_total = (double)av_clonesize_total/av_ctr_total;
	
	for (int i=0; i<CSD_basal_av.size(); i++) // !!!
	{
		//		file_errorbars[i] << endl;
		file_pval << endl;
		file_basal_errorbars << bins[i] << '\t' << CSD_basal_av[i] << '\t' << errorbar_basal_minus[i] << '\t' << errorbar_basal_plus[i] << endl;
		file_basal_cumul_errorbars << bins[i] << '\t' << cumuldistr_basal_av[i] << '\t' << errorbar_cumul_basal_minus[i] << '\t' << errorbar_cumul_basal_plus[i] << endl;
		file_basal_errorbars_xmgr << bins[i] << '\t' << CSD_basal_av[i] << '\t' << errorbar_basal_plus[i] - CSD_basal_av[i] << '\t' << CSD_basal_av[i] - errorbar_basal_minus[i] << endl;
		file_basal_errorbars_top << bins[i] << '\t' << errorbar_basal_plus[i] << endl;
		file_basal_errorbars_bottom << bins[i] << '\t' << errorbar_basal_minus[i] << endl;
        file_basal_errorbars_resc << (double)bins[i]/av_clonesize_basal << '\t' << CSD_basal_av[i]*av_clonesize_basal << '\t' << errorbar_basal_minus[i]*av_clonesize_basal << '\t' << errorbar_basal_plus[i]*av_clonesize_basal << endl;

	}
	
	for (int i=0; i < CSD_total_av.size(); i++) // !!!
	{
		file_total_errorbars << bins[i] << '\t' << CSD_total_av[i] << '\t' << errorbar_total_minus[i] << '\t' << errorbar_total_plus[i] << endl;
		file_total_cumul_errorbars << bins[i] << '\t' << cumuldistr_total_av[i] << '\t' << errorbar_cumul_total_minus[i] << '\t' << errorbar_cumul_total_plus[i] << endl;
	};
	
	sort(floatclone_freq.begin(),floatclone_freq.end(),compare_floating);
	floatclone_freq.resize(floor(siglevel*(double)floatclone_freq.size()));
	floatingclone_errorbar[0] = floatingclone_frac_av;					   
	floatingclone_errorbar[1] = *min_element(floatclone_freq.begin(),floatclone_freq.end());					   
	floatingclone_errorbar[2] = *max_element(floatclone_freq.begin(),floatclone_freq.end());
	cerr << "floatingclone + errorbars: " << endl << floatingclone_errorbar[0] << '\t' << floatingclone_errorbar[1] << '\t' << floatingclone_errorbar[2] << endl; 
	
		double p_val_av=0;
		double p_val_msq=0;
	
#if raw_data == 6 || raw_data == 7
		int i_limit = limit_basal;
		int j_limit = limit_supra;
#else
	int i_limit = 7;
	int j_limit = 7;
#endif
	
//	for (int i=0; i<i_limit; i++) 
//	{
//		for (int j=start_val(i); j<i_limit; j++) 
//		{
//			p_val_av = p_val_av + p_val[i][j];
//			p_val_msq = p_val_msq + p_val[i][j]*p_val[i][j];
//		};
//	};
//	
//	p_val_av = (double)p_val_av/(i_limit*j_limit);
//		p_val_msq = sqrt((double)p_val_msq/(i_limit*j_limit));
//	
//	cerr << "clonenumber = " << run_bestfit->clonenumber << endl;
//	cerr << "average p-value = " << p_val_av << endl;
//		cerr << "mean square p-value = " << p_val_msq << endl;
	
	for (int i=0; i<esmblnumber; i++) delete run_old[i];

	
//	for (int i=0; i<fileerrorbar_nr; i++) delete file_errorbars[i];

//#endif
	
}

//-----------------------------------------------------------------
// Parameter Tuning Functions
//-----------------------------------------------------------------
// These functions are called by tuning the sliders and adjust parameter values

void MainWindow::set_reprodrate0(int rate0_int)
{
	reprod_rate_init[0]=(double)rate0_int*minstep_reprod;
};
void MainWindow::set_reprodrate1(int rate0_int)
{
	reprod_rate_init[1]=(double)rate0_int*minstep_reprod;
};
void MainWindow::set_r0(int rate0_int)
{
	r_init[0]=(double)rate0_int*minstep_r;
};
void MainWindow::set_r1(int rate0_int)
{
	r_init[1]=(double)rate0_int*minstep_r;
};
void MainWindow::set_dediffrate(int rate0_int)
{
	dediff_rate_init[1]=(double)rate0_int*minstep_dediff;
};
void MainWindow::set_delta0(int rate0_int)
{
	delta_init[0]=(double)rate0_int*minstep_delta;
};
void MainWindow::set_delta1(int rate0_int)
{
	delta_init[1]=(double)rate0_int*minstep_delta;
};
void MainWindow::set_shedrate(int rate0_int)
{
	shed_rate_init[2]=(double)rate0_int*minstep_shed;
};
void MainWindow::set_stratrate(int rate0_int)
{
	strat_rate_init[2]=(double)rate0_int*minstep_strat;
};
void MainWindow::set_directstratfrac(int rate0_int)
{
	direct_strat_frac_init[2]=(double)rate0_int*minstep_direct_strat;
};




//--------------------------------------------------------------
// Set up Widgets
//--------------------------------------------------------------

#if Qt_flag==1

void MainWindow::setup_widgets()
{
	
	timer = new QTimer;
	output_timer = new QTimer;
	image_timer = new QTimer;
	
	
	double minstep_reprod = 0.001;                        // Attention: new ones might need to be added for more processes
	double minstep_symmdiv = 0.001;
	double minstep_dediff = 0.001;
	double minstep_delta = 0.001;
	double minstep_shed = 0.001;
	double minstep_strat = 0.001;
	double minstep_direct_strat = 0.001;

	
	//------------------ Create Widgets -----------------------------------
	
	label = new QLabel(tr("Count:"));
	
	startbutton = new QPushButton("Start");
	stopbutton = new QPushButton("Stop");
	continuebutton = new QPushButton("Continue");
	resetbutton = new QPushButton("Reset");
	quitbutton = new QPushButton("Quit");
	snapshotbutton = new QPushButton("SnapSh");
	startvideobutton = new QPushButton("Start Video");
	stopvideobutton = new QPushButton("Stop Video");
	
	
	LabeledDisplay* display_cellnumbertot= new LabeledDisplay(this,"total cell number",0);
	displays.push_back(display_cellnumbertot);
	
	for (int i=0; i<celltypenumber; i++) 
	{
		string displayname;
		displayname.append("display");
		displayname.append("_");
		displayname.append(itoa(i));

		LabeledDisplay* display_Cnumber = new LabeledDisplay(this,displayname.c_str(),0);
		displays.push_back(display_Cnumber);

	};
	
	display_counter = new LabeledDisplay(this,"run time",0);
	
	

	LabeledSlider* slider_reprod0 = new LabeledSlider(this,"reprod. rate 0",minstep_reprod,1,reprod_rate_init[0],0); // Attention: Model-dependent
	sliders.push_back(slider_reprod0);
	LabeledSlider* slider_r = new LabeledSlider(this,"r0",minstep_r,0.3,r[0],0);
	sliders.push_back(slider_r);
	LabeledSlider* slider_delta0 = new LabeledSlider(this,"delta0",minstep_delta,1,delta[0],0);
	sliders.push_back(slider_delta0);
	LabeledSlider* slider_strat = new LabeledSlider(this,"strat",minstep_strat,10,strat_rate_init[1],0);
	sliders.push_back(slider_strat);
	LabeledSlider* slider_shed = new LabeledSlider(this,"shed",minstep_shed,10,shed_rate_init[2],0);
	sliders.push_back(slider_shed);
	LabeledSlider* slider_direct_strat0 = new LabeledSlider(this,"direct strat",minstep_direct_strat,1,direct_strat_frac[0],0);
	sliders.push_back(slider_direct_strat0);
	
	
	output = new PaintWindow(this);
	
	

	//spinbox = new QSpinBox;
	
	//slider->setRange(0, 100);
	//spinbox->setRange(0, 100);
	//------------------------------------------------------------------------------
	
	
	//--------------------------- Connect Widgets with Functions ---------------------------------------------
	
	connect(startbutton, SIGNAL(clicked()), this, SLOT(start()));
	connect(stopbutton, SIGNAL(clicked()), this, SLOT(stop()));
	connect(continuebutton, SIGNAL(clicked()), this, SLOT(continue_sim()));
	connect(quitbutton,SIGNAL(clicked()),this,SLOT(finish()));
	connect(resetbutton, SIGNAL(clicked()), this, SLOT(reset()));
	connect(snapshotbutton, SIGNAL(clicked()), output, SLOT(paint_file()));
	connect(startvideobutton, SIGNAL(clicked()), this, SLOT(start_video()));
	connect(stopvideobutton, SIGNAL(clicked()), this, SLOT(stop_video()));
	
	connect(sliders[0]->get_slider(), SIGNAL(valueChanged(int)), this, SLOT(set_reprodrate0(int)));     // Debug !!! Attention: Change slider indices
	connect(sliders[1]->get_slider(), SIGNAL(valueChanged(int)), this, SLOT(set_r0(int)));     // when altering applied sliders
	connect(sliders[2]->get_slider(), SIGNAL(valueChanged(int)), this, SLOT(set_stratrate(int)));
	connect(sliders[3]->get_slider(), SIGNAL(valueChanged(int)), this, SLOT(set_shedrate(int)));
	connect(sliders[4]->get_slider(), SIGNAL(valueChanged(int)), this, SLOT(set_delta0(int)));
	connect(sliders[5]->get_slider(), SIGNAL(valueChanged(int)), this, SLOT(set_directstratfrac(int)));
	
	connect(timer, SIGNAL(timeout()), this, SLOT(iterate()));
	connect(timer, SIGNAL(timeout()), displays[0], SLOT(display_val()));
	connect(output_timer, SIGNAL(timeout()), output, SLOT(update()));
	connect(image_timer, SIGNAL(timeout()), output, SLOT(paint_file()));
	
	
	//	slider->setValue(int(w_create1*100));
	//spinbox->setValue(int(w_create1*100));	
	
	
	
	//	QVBoxLayout *sliderlayout = new QVBoxLayout; 
	//	sliderlayout->addWidget(sliderlabel); 
	//	sliderlayout->addWidget(slider_create1); 
	
	//---------------- Arrange Widgets in Window -----------------------
	
	QHBoxLayout *toptoplayout = new QHBoxLayout; 
	toptoplayout->addWidget(startbutton); 
	toptoplayout->addWidget(resetbutton);
	toptoplayout->addWidget(quitbutton);
		toptoplayout->addWidget(stopbutton);
		toptoplayout->addWidget(continuebutton);
	
	QHBoxLayout *topbottomlayout = new QHBoxLayout; 
	topbottomlayout->addWidget(snapshotbutton);
	topbottomlayout->addWidget(startvideobutton);
	topbottomlayout->addWidget(stopvideobutton);
	
	
	QVBoxLayout *toplayout = new QVBoxLayout; 
	toplayout->addLayout(toptoplayout); 
	toplayout->addLayout(topbottomlayout);
	
	QVBoxLayout *middlelayout = new QVBoxLayout; 
	for (int i=0; i<sliders.size(); i++) 
	{
	middlelayout->addWidget(sliders[i]); 
	};

	
	QVBoxLayout *toprightlayout = new QVBoxLayout;
	for (int i=0; i<displays.size(); i++) 
	{
		toprightlayout->addWidget(displays[i]);
	};
	
	
	toprightlayout->addWidget(display_counter);
	
	QVBoxLayout *rightlayout = new QVBoxLayout; 
	rightlayout->addLayout(toprightlayout);
	rightlayout->addLayout(middlelayout);
	
	QVBoxLayout *leftlayout = new QVBoxLayout; 
	leftlayout->addLayout(toplayout); 
	leftlayout->addWidget(output); 

	
	QHBoxLayout *mainlayout = new QHBoxLayout;
	mainlayout->addLayout(leftlayout);
	
	mainlayout->addLayout(rightlayout);
	
	
	setLayout(mainlayout);
	//-----------------------------------------------------------------
	
};
//---------------------------------------------- End Set Up of Simulation

#endif
	
	
//-----------------------------------------------------------------
// Start Routine
//-----------------------------------------------------------------

void MainWindow::start()
{
	ibm = mrand48();
	
#if Qt_flag==1                   // Start timers for iteration and display configuration
	timer->start(itr_intv);
	output_timer->start(output_intv);
#endif

}
//--------------------------------------------------------------------- End Start Routine	
	

//---------------------------------------------------------------------
// Other Window Buttons
//---------------------------------------------------------------------


#if Qt_flag==1
void MainWindow::stop()
{
	cout << "program manually stopped/n" ;
	timer->stop();
	output_timer->stop();
	//image_timer->stop();
};

void MainWindow::continue_sim()
{
	cout << "continue\n";
	timer->start(itr_intv);
	output_timer->start(output_intv);
	//image_timer->start(image_intv);
};

void MainWindow::start_video()
{
	image_timer->start(image_intv);
};

void MainWindow::stop_video()
{
	image_timer->stop();
};


//------------------------ Reset and Restart Simulation -------------------------------------
void MainWindow::reset()
{
	stop();
	reset_sim();
	setup_sim();
	
//	lattice.clear_grid();
	start();
	
}

#endif

//---------------------------------------------------------------------
// Finishing Routine
//---------------------------------------------------------------------
	

void MainWindow::finish()
{
	cerr << "program ended normally\n";
	exit(0);
};

void MainWindow::print_data(run* run0)
{
	cerr << "print data!" << endl;
	double data_counter_total=0;
	double data_counter_basal=0;
	double av_clonesize_basal=0; 
	double av_clonesize_total=0;
	int av_ctr_basal=0;
	int av_ctr_total=0;
	for (int i=(start_valid_b-1)/binsize; i<run0->clonesizedistr_basal.size(); i++) 
	{
		av_clonesize_basal = av_clonesize_basal + i*run0->clonesizedistr_basal[i];
		av_ctr_basal = av_ctr_basal + run0->clonesizedistr_basal[i];
	};
	for (int i=1; i<run0->clonesizedistr.size(); i++) 
	{
		av_clonesize_total = av_clonesize_total + i*run0->clonesizedistr[i];
		av_ctr_total = av_ctr_total + run0->clonesizedistr[i];
	};
	
	av_clonesize_basal = (double)av_clonesize_basal/av_ctr_basal;
	av_clonesize_total = (double)av_clonesize_total/av_ctr_total;
	
	
	for (int i=0; i<run0->clones.size(); i++) 
	{
		file_H2B_divctr_scatter << run0->clones[i].get_av_divctr() << '\t' << run0->clones[i].get_size() << endl;
	};
	
	vector<double> H2B_av_clonesize(2,0);
	vector<double> H2B_av_clonesize_ctr(2,0);
	int intv=1;
	int max_divctr=2;
	for (int i=0; i<run0->clones.size(); i++) 
	{
		double curr_divctr=floor(run0->clones[i].get_av_divctr());
		int curr_clonesize=run0->clones[i].get_size_prog();
		
		if (curr_divctr > max_divctr) 
		{
			max_divctr = (int)curr_divctr;
			H2B_av_clonesize.resize(max_divctr+1,0);
			H2B_av_clonesize_ctr.resize(max_divctr+1,0);
		};
		
		if (curr_clonesize>0) 
		{
			H2B_av_clonesize[(int)curr_divctr/intv] = H2B_av_clonesize[(int)curr_divctr/intv] + curr_clonesize;
			H2B_av_clonesize_ctr[(int)curr_divctr/intv]++;
		};
	};
	
	for (int i=0; i<max_divctr; i++) 
	{
		if(H2B_av_clonesize_ctr[i] > 0)
		{
			H2B_av_clonesize[i] = (double)H2B_av_clonesize[i]/H2B_av_clonesize_ctr[i];
		};
	};
	
    Cnumber_glob.clear();
    Cnumber_glob.resize(celltypenumber,0);
    vector<cell*> cells0 = run0->cells;
    
    for (int i=0; i<cells0.size(); i++)
    {
        Cnumber_glob[cells0[i]->get_celltype()]++;
    };
    
    for (int i=0; i<Cnumber_glob.size(); i++)
    {
        file_celltypes_basal << i << '\t' << Cnumber_glob[i] << endl;
    };
    
    
	for (int i=0; i<H2B_av_clonesize.size(); i++) 
	{
		file_H2B_av_clonesize << i << '\t' << H2B_av_clonesize[i] << endl;
	};
	
	//cerr << "av_clonesize_basal = " << av_clonesize_basal << endl;                  // Debug !!!
//	exit(2);
	
	for (int i=start_valid_b-1; i<run0->clonesizedistr.size(); i++) 
	{
		file_CSD << i << '\t' << run0->clonesizedistr[i] << endl; 
		file_CSD_norm << i << '\t' << run0->clonesizedistr_norm[i] << endl;
		
		file_CSD_basal << i << '\t' << run0->clonesizedistr_basal[i] << endl; 
		file_CSD_supra << i << '\t' << run0->clonesizedistr_supra[i] << endl;
		
		file_CSD_basal_resc << (double)i/av_clonesize_basal << '\t' << run0->clonesizedistr_basal[i]*av_clonesize_basal << endl;
		file_CSD_total_resc << (double)i/av_clonesize_total << '\t' << run0->clonesizedistr[i]*av_clonesize_total << endl;
				
		data_counter_total = data_counter_total + run0->clonesizedistr[i];
		data_counter_basal = data_counter_basal + run0->clonesizedistr_basal[i];
		
		file_cumuldistr_basal << i << '\t' << (double)clone_norm - data_counter_basal << endl;
		file_cumuldistr_total << i << '\t' << (double)clone_norm - data_counter_total << endl;
		
	};
	
	for (int i=0; i<run0->CSD_intm_basal[0].size(); i++) 
	{
		file_CSD_basal_intm << i << '\t' << run0->CSD_intm_basal[0][i] << endl;
	};
	
	for (int i=0; i<run0->joint_distr.size(); i++) 
	{
		for (int j=0; j < run0->joint_distr[i].size(); j++) 
		{
			file_Jointdistr << i << '\t' << j << '\t' << run0->joint_distr[i][j] << endl;
			file_Jointdistr_norm << i << '\t' << j << '\t' << run0->joint_distr_norm[i][j] << endl;
			file_Jointdistr_matrix << run0->joint_distr[i][j] << '\t';

		};
		file_Jointdistr << endl;
		file_Jointdistr_norm << endl;
		file_Jointdistr_matrix << endl;
		
		double sum=0;
		double av_ctr=0;
		for (int j=0; j < run0->joint_distr[i].size(); j++) 
		{
			av_ctr = av_ctr + run0->joint_distr[i][j];
			sum = sum + j*run0->joint_distr[i][j];
		};
		if (av_ctr > 0) file_avsupra_basal << i << '\t' << (double)sum/av_ctr << endl;
	};
	
	for (int i=0; i<run0->joint_distr_intm[0].size(); i++) 
	{
		for (int j=0; j < run0->joint_distr_intm[0][i].size(); j++) 
		{
			file_Jointdistr_intm << i << '\t' << j << '\t' << run0->joint_distr_intm[0][i][j] << endl;
			file_Jointdistr_matrix_intm << run0->joint_distr_intm[0][i][j] << '\t';
			
		};
		file_Jointdistr_intm << endl;
		file_Jointdistr_matrix_intm << endl;
	};
	
	int dataset=run0->logctr-1;
    
 //////////// Determine mutation distribution
    int max_mutnr=2;
    vector<int> mut_distr(max_mutnr+1,0);
    mut_distr.reserve(100);
    cells0 = run0->cells;
    int ctr=0;
    
    for (int i=0; i<cells0.size(); i++)
    {
        if (i%1000 == 0) cerr << i << endl;              // Debug !!!
        if (cells0[i]->get_celltype() == 0)
        {
            
            int curr_mutctr = cells0[i]->get_mut_ctr();
            if (curr_mutctr > mut_distr.size())
            {
                max_mutnr = curr_mutctr;
                cerr << "maxmutnr = " << max_mutnr << endl;             // Debug !!!
                mut_distr.resize(max_mutnr+2);
            };
            
            mut_distr[curr_mutctr]++;
            ctr++;
        };
    };
    
    for (int i=0; i<mut_distr.size(); i++)
    {
        file_mut_distr << i << '\t' << (double)mut_distr[i]/ctr*100 << endl;
    };

    //////// End mutation distribution
    
	for (int j=0; j<data0[dataset][0].size(); j++)
	{
		for (int i=0; i < data0[dataset].size(); i++) 
		{
			file_Jointdistr_transp << j << '\t' << data0[dataset][i][j] << '\t';
		};
		file_Jointdistr_transp << endl;
	
	}

	file_parameters << run0->reprod_rate[0] << endl << run0->r[0] << endl << run0->delta[0] << endl << run0->strat_rate[1] << endl << run0->shed_rate[2] << endl;  // Attention !!! Model-dependent
#if raw_data == 8
	file_parameters << run0->rho << endl;
#endif
	
	file_CSD.close();
	file_CSD_norm.close();
	file_CSD_basal.close();
	file_CSD_supra.close();
	file_Jointdistr.close();
		
};		
		
//---------------------------------------- End Finishing Routine

