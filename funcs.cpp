/*
 *  funcs.cpp
 *  popdyn_1D_discr
 *
 *  Created by Philip Greulich on 08.10.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "funcs.h"


string itoa(int n)
{
	string str;
	strstream help;
	help << n;
	help >> str;
	return str;
	
};

double sgn (double x)
{
	if (x>=0) 
	{
		return 1.0;
	}
	else 
	{
		return -1.0;
	};
};



vector<ofstream*> make_filelist(string& prefix, int filenumber)
{
	vector<ofstream*> files;
	
	for (int i=0; i<filenumber; i++)  // Attention: start at 0 or at 1 ??? !!! Bugsource
	{
		string filename;
		
		filename.append(prefix);
		filename.append("_");
		filename.append(itoa(i));
		filename.append(".dat");
		
		ofstream* curr_file;
		curr_file = new ofstream;
		
		cerr << "filename="<< filename << endl;
		
		(*curr_file).open(filename.c_str());
		
		if (!(*curr_file).is_open())
		{
			cerr << "error: file could not been opened\n";
			exit(1);
		};
		
		files.push_back(curr_file);
		
	};
	
	return files;
	
};

ofstream* make_file(string& prefix, int number1, int number2, int number3, char fileoption0)
{
	ofstream* file;
	
	string filename;
	
	filename.append(prefix);
//	filename.append("_");
	filename.append(itoa(number1));
	if (number2 >= 0)
	{
		filename.append("_");
		filename.append(itoa(number2));
	};
	if (number3 >= 0)
	{
		filename.append("_");
		filename.append(itoa(number3));
	};
	filename.append(".dat");
	file = new ofstream;
	
	cerr << "filename="<< filename << endl;
	
	if (fileoption0 != 'a') 
	{
		(*file).open(filename.c_str());
	}
	else if (fileoption0 == 'a')                     // Debug !!! Only else-if for checking 
	{
		(*file).open(filename.c_str(),ios::app);
	}
	else 
	{
		cerr << "error: invalid file option\n";
		exit(1);
	}


	
	if (!(*file).is_open())
	{
		cerr << "error: file could not been opened\n";
		exit(1);
	};
	
	
	
	return file;
	
};

vector<vector<double> > read_file(ifstream* file)
{
	vector<vector<double> > matrix;
	int comment_flag=0;
	
	while (file->good()) 
	{
		comment_flag=0;
		vector<double> line;
		
		double x;
		
		if (file->peek() == '#') comment_flag=1;       // '#' character at beginning of line is interpreted as comment and line ignored
		
		(*file) >> x;
		line.push_back(x);
		
		while (file->peek() != '\n' && file->peek() != '\r' && file->peek() != '\r\n' && file->good()) 
		{
			
			//	double x;
			
			(*file) >> x;
			
			
			line.push_back(x);
			
		};
		
		if(file->good() && comment_flag == 0)                            // Bug Source !!! file must have '\n' at end, else last line not read. 
		{
			matrix.push_back(line);
		};
	};
	
	
	return matrix;
}

vector<vector<int> > read_file_int(ifstream* file)
{
	vector<vector<int> > matrix;
	int comment_flag=0;
	
	while (file->good()) 
	{
		comment_flag=0;
		vector<int> line;
		
		int x;
		
		if (file->peek() == '#') comment_flag=1;       // '#' character at beginning of line is interpreted as comment and line ignored
		
		(*file) >> x;
		line.push_back(x);
		
		while (file->peek() != '\n' && file->peek() != '\r' && file->peek() != '\r\n' && file->good()) 
		{
			
			//	double x;
			
			(*file) >> x;
			
			
			line.push_back(x);
			
		};
		
		if(file->good() && comment_flag == 0)                            // Bug Source: file must have '\n' at end, else last line not read. 
		{
			matrix.push_back(line);
		};
	};
	
	
	return matrix;
}

vector<vector<double> > random_parameters()
{
	vector<vector<double> > parameters;
	parameters.clear();
	
	vector<double> parameter0;
	
	parameter0.clear();
	parameter0.resize(celltypenumber,0);
	
	parameter0[0] = reprodratelimit_low + drand48()*(reprodratelimit_high - reprodratelimit_low);
	parameters.push_back(parameter0);          // reprod_rate
	
	parameter0.clear();
	parameter0.resize(celltypenumber,0);
	
	parameter0[0] = drand48()*0.5;
	parameters.push_back(parameter0);        // r
	
	parameter0.clear();
	parameter0.resize(celltypenumber,0);
	parameters.push_back(parameter0);         // dediff_rate
	

	parameter0.clear();
	parameter0.resize(celltypenumber,0);       
#if raw_data == 5 || raw_data == 6 || raw_data == 7
	parameter0[0] = 1.0;
#endif
	parameters.push_back(parameter0);           // delta	
	
	parameter0.clear();
	parameter0.resize(celltypenumber,0);
	parameter0[1] = stratratelimit_low + drand48()*(stratratelimit_high - stratratelimit_low);
	parameters.push_back(parameter0);         // strat_rate
	
	parameter0.clear();
	parameter0.resize(celltypenumber,0);
	parameter0[2] = shedratelimit_low + drand48()*(shedratelimit_high - shedratelimit_low);
	parameters.push_back(parameter0);         // shed_rate
	
	parameter0.clear();
	parameter0.resize(celltypenumber,0);
	parameters.push_back(parameter0);         // direct_strat_frac
	
#if raw_data == 8
	parameter0.clear();
	parameter0.resize(celltypenumber,0);
	parameter0[0]=drand48()*0.99+0.01;
	parameters.push_back(parameter0);
	
	cerr << "parameter[7] = " << parameters[7][0] << endl;         // Debug !!!
#endif
	return parameters;
}

// ----------------------------------------------------------------------
// Functions to compare whether smaller or larger than reference value
// ----------------------------------------------------------------------

bool compare_joint(double x, double y)
{
	return(abs(x - joint_distr_av[i_curr][j_curr]) < abs(y - joint_distr_av[i_curr][j_curr]));
};

bool compare_basal(double x, double y)
{
	return(abs(x - CSD_basal_av[i_curr]) < abs(y - CSD_basal_av[i_curr]));
};

bool compare_total(double x, double y)
{
	return(abs(x - CSD_total_av[i_curr]) < abs(y - CSD_total_av[i_curr]));
};

bool compare_cumul_basal(double x, double y)
{
	return(abs(x - cumuldistr_basal_av[i_curr]) < abs(y - cumuldistr_basal_av[i_curr]));
};

bool compare_cumul_total(double x, double y)
{
	return(abs(x - cumuldistr_total_av[i_curr]) < abs(y - cumuldistr_total_av[i_curr]));
};

bool compare_floating(double x, double y)
{
	return(abs(x - floatingclone_frac_av) < abs(y - floatingclone_frac_av));
};

bool compare_likely(vector<double> x, vector<double> y)
{
	return(x.back() > y.back());
};

bool compare_pairs(vector<double> x, vector<double> y)
{
	if (x[0] == y[0]) 
	{
		return(x[1] > y[1]);
	}
	else 
	{
		return(x[0]>y[0]);
	};
};

//----------------------------------------------------------------------

int start_val(int i)
{
	if (i >= start_valid_b) return 0;
	else if (i == start_valid_b - 1) return start_valid_s;
	else return 100000;
};


vector<vector<double> > make_switchparams(const vector<vector<double> >& parameters)
{
	vector<vector<double> > parameters_switch = parameters;
	
	parameters_switch[0][0] = switch_reprodrate;
	parameters_switch[1][0] = switch_r;
	parameters_switch[3][0] = switch_delta;
	
	return parameters_switch;
}

//-----------------------------------------------------------------------------------
//		Functions for binning
//-----------------------------------------------------------------------------------

vector<vector<double> > bin(const vector<vector<double> >& CSD, char mode, int binsize)    // for joint CSD
{
  vector<vector<double> > CSD_binned;
  CSD_binned.reserve(((int)CSD.size()+2*binsize)/binsize);
  bins.clear();

	if (mode == 'l')
	  {
	    for (int i=0; binsize*i<CSD.size(); i++)
	      {
			  bins.push_back(binsize*i);
			  vector<double> aux(2);
			  aux[0]=bins.back();
			  aux[1]=0;
			  CSD_binned.push_back(aux);
	      };
	  }
	else if (mode == 'e')
	  {
	    for (int i=0; powl(binsize,i)<CSD.size(); i++)
	      {
		bins.push_back(powl(binsize,i));
		vector<double> aux(2);
		aux[0]=powl(binsize,i);//bins.back();
		aux[1]=0;
		CSD_binned.push_back(aux);
	      };
	  }
	else {cerr << "error: wrong binning mode\n"; exit(1);}

	if (CSD_binned.size() != bins.size()) {cerr << "error: CSD_binned not large enough\n"; exit(1);}

	for (int i=0; i<bins.size()-1; i++)
	  {
	    for (int j=bins[i]; j<bins[i+1]; j++)
	      {
			  if (j >= start_valid_b-1) CSD_binned[i][1] = CSD_binned[i][1] + CSD[j][1];
	      };
	  };
	for (int j=bins.back(); j<CSD.size(); j++)
	{
		CSD_binned.back()[1] = CSD_binned.back()[1] + CSD[j][1];
	};

	return CSD_binned;

}

vector<double> bin(const vector<double>& CSD, char mode, int binsize)
{
	vector<double> CSD_binned;
	CSD_binned.reserve(((int)CSD.size()+2*binsize)/binsize);
	bins.clear();

	if (mode == 'l')
	  {
	    for (int i=0; binsize*i<CSD.size()+1; i++)
	      {
		bins.push_back(binsize*i);
		//		vector<double> aux(2);
		//	aux[0]=bins.back();
		//	aux[1]=0;
		//	CSD_binned.push_back(aux);
		CSD_binned.push_back(0);
	      };
	  }
	else if (mode == 'e')
	  {
	    for (int i=0; powl(binsize,i)<CSD.size()+1; i++)
	      {
		bins.push_back(powl(binsize,i));
		//		vector<double> aux(2);
		//	aux[0]=bins.back();
		//	aux[1]=0;
		//	CSD_binned.push_back(aux);
		CSD_binned.push_back(0);
	      };
	  }
	else {cerr << "error: wrong binning mode\n"; exit(1);}

	if (CSD_binned.size() != bins.size()) {cerr << "error: CSD_binned not large enough\n"; exit(1);}

	for (int i=0; i<bins.size()-1; i++)
	  {
	    for (int j=bins[i]; j<bins[i+1]; j++)
	      {
			  if (j >= start_valid_b-1) CSD_binned[i] = CSD_binned[i] + CSD[j];
	      };
	  };
	for (int j=bins.back(); j<CSD.size(); j++)
	{
		CSD_binned.back() = CSD_binned.back() + CSD[j];
	};

	return CSD_binned;

}

vector<vector<double> > bin_cont(vector<double> sizes, double binsize_cont)
{
	double maxsize = *max_element(sizes.begin(),sizes.end());
	
	vector<vector<double> > size_distr(maxsize/binsize_cont+2);  
	for (int i=0; i<size_distr.size(); i++) 
	{
		size_distr[i].resize(2);
		size_distr[i][0] = binsize_cont*i;
	};
	
	for (int i=0; i<sizes.size(); i++) 
	{
		int bin_flag=0;
		for (int j=0; j<size_distr.size() && bin_flag==0; j++) 
		{
			if (sizes[i] < binsize_cont*j)
			{
				size_distr[j][1] = size_distr[j][1] + 1.0;
				bin_flag=1;
			};
		};
	};
	
	for (int i=0; i<size_distr.size(); i++) 
	{
		size_distr[i][1] = (double)size_distr[i][1]/sizes.size()*100;
	};
	
	return size_distr;
}

//-----------------------------------------------------------------------------------

vector<vector<vector<double> > > read_filesmbl(vector<ifstream*>& files)
{
	vector<vector<vector<double > > > filelist;
	
	for (int i=0; i<files.size(); i++) 
	{
		vector<vector<double> > data=read_file(files[i]);
		
		filelist.push_back(data);
		cerr << endl;
	};
	
	return filelist;
};

vector<ifstream*> make_filelist(string& prefix, vector<vector<int> >& filenumbers)
{
	vector<ifstream*> files;
	
	for (int i=0; i<filenumbers.size(); i++)
	{
		string filename;
		
		filename.append(prefix);
		
		for (int j=0; j<filenumbers[i].size(); j++) 
		{
			
			filename.append("_");
			filename.append(itoa(filenumbers[i][j]));
		};
		
		filename.append(".dat");
		
		ifstream* curr_file;
		curr_file = new ifstream;
		
		cerr << "filename="<< filename << endl;
		
		(*curr_file).open(filename.c_str());
		
		if (!(*curr_file).is_open())
		{
			cerr << "error: file could not been opened\n";
			cerr << "in make_filelist()\n";
			exit(1);
		}
		else 
		{
			cerr << "file " << filenumbers[i][0] << '\t' << filenumbers[i][1] << " opened\n";
			files.push_back(curr_file);
		}
		
		
		
	};
	
	return files;
	
};

double prior(double reprodrate0, double r0, double stratrate0, double shedrate, double delta)         // Attention !!! Adjust for each individual BAyesian prior function
{
  return sqrt(2*(4*atan(1))*var_bayes)*exp(-(reprodrate0 - reprodrate_bayes_av)*(reprodrate0 - reprodrate_bayes_av)/2/var_bayes);
};


//------------------------------------------------------------------------------
// Evaluate best fit 
//------------------------------------------------------------------------------

vector<double> eval_loglikelyfiles(vector<ifstream*>& file_loglikelies_in, double delta_intv, double reprod_intv, double r_intv, double strat_intv, double shed_intv)
{
	vector<vector<vector<double> > > loglikelies = read_filesmbl(file_loglikelies_in);
	vector<vector<double> > likelies_prep;
	vector<vector<double> > sig_params;
	vector<vector<double> > likely_scatter_sort;
	vector<double> error_plus(param_number,0);
	vector<double> error_minus(param_number,0);
	vector<double> standard_dev(param_number,0);
	vector<int> error_plus_ctr(param_number,0);
	vector<int> error_minus_ctr(param_number,0);
	vector<int> standard_dev_ctr(param_number,0);
	list<vector<double> > likely_list;
	
	ofstream file_likely_scatter("likely_scatter.dat");
	ofstream file_sigparams_limits("sigparams_limits.dat");
	ofstream file_likely_countour("likely_countour.dat");
	double max_llikely = -10000; 
	vector<double> best_params(param_number);
	
	for (int i=0; i<param_number; i++) 
	{
		best_params[i] = -1;
	};
	
	for (int index = 0; index < loglikelies.size(); index++) 
	{
		for (int i=1; i<loglikelies[index].size(); i++) 
		{
			vector<double> aux(0);
			aux.push_back(loglikelies[index][0][0]);        // delta
			aux.push_back(loglikelies[index][i][0]);		// reprodrate
			aux.push_back(loglikelies[index][0][2]);		// r
			aux.push_back(loglikelies[index][i][1]);		// stratrate
			aux.push_back(loglikelies[index][i][2]);		// shedrate
			aux.push_back(loglikelies[index][i][3]);		// likelihood
			likelies_prep.push_back(aux);
			
			double curr_llikely = loglikelies[index][i][3];
			if (curr_llikely > max_llikely)
			{
				max_llikely = curr_llikely;
				best_params[0] = loglikelies[index][0][0];
				best_params[1] = loglikelies[index][i][0];
				best_params[2] = loglikelies[index][0][2];
				best_params[3] = loglikelies[index][i][1];
				best_params[4] = loglikelies[index][i][2];
			};
		};
	};
	
	double likely_sum=0;
	
	likely_list.clear();
	for (int i = 0; i < likelies_prep.size(); i++) 
	{
		likelies_prep[i].back() = exp(likelies_prep[i].back()-max_llikely);
		likely_sum = likely_sum + likelies_prep[i].back();
		likely_list.push_back(likelies_prep[i]);
	};
	
	int col1=2;
	int col2=3;

    vector<vector<double> > likelies_projected(0);
	while (!likely_list.empty())
	{
		double x = likely_list.begin()->at(col1);
		double y = likely_list.begin()->at(col2);
		
		double likely_curr=likely_list.begin()->back();
		
		likely_list.erase(likely_list.begin());
		
		for (list<vector<double> >::iterator itr = likely_list.begin(); itr != likely_list.end(); itr++)
		{
			if (itr->at(col1)==x && itr->at(col2)==y) 
			{
				likely_curr = likely_curr + itr->back();
				list<vector<double> >::iterator help_itr = itr;
				itr--;
				likely_list.erase(help_itr);
//				cerr << "erase" << endl;                    // Debug !!!
			};
		};
		
        vector<double> aux(0);
        
		aux.push_back(x);
		aux.push_back(y);
		aux.push_back(likely_curr);
        
        likelies_projected.push_back(aux);
        
		
	};
	
    for(int i=0; i< likelies_projected.size(); i++)
    {
        if (likelies_projected[i][2] > max_llikely) max_llikely = likelies_projected[i][2];
    };
    
    for (int i=0; i<likelies_projected.size(); i++)
    {
        file_likely_countour << likelies_projected[i][0] << '\t' << likelies_projected[i][1] << '\t' << (double)likelies_projected[i][2]/max_llikely << endl;
    };

	
	for (int ctr=0; ctr < scatterpt_number; ctr++)
	{
		double rand0 = drand48()*(double)likely_sum;
		double rand_sum=0;
		double rand_sum_help=0;
		
		int k=0;
		for (k=0; k<likelies_prep.size() && rand_sum < rand0; k++) 
		{
			rand_sum_help = rand_sum;
			rand_sum = rand_sum + likelies_prep[k].back();
		};
	
		vector<double> aux;
		aux.push_back(likelies_prep[k-1][0] + (drand48()-0.5)*delta_intv);
		aux.push_back(likelies_prep[k-1][1]  + (drand48()-0.5)*reprod_intv);
		aux.push_back(likelies_prep[k-1][2] + (drand48()-0.5)*r_intv);
		aux.push_back(likelies_prep[k-1][3] + (drand48()-0.5)*strat_intv);
		aux.push_back(likelies_prep[k-1][4] + (drand48()-0.5)*shed_intv);
		aux.push_back(likelies_prep[k-1].back());
		
		for (int j = 0; j < param_number; j++) 
		{
			if (aux[j] > best_params[j]) 
			{
				error_plus[j] = error_plus[j] + (aux[j] - best_params[j]);
				error_plus_ctr[j]++;
			}
			else 
			{
				error_minus[j] = error_minus[j] - (aux[j] - best_params[j]);
				error_minus_ctr[j]++;
			};
			
			standard_dev[j] = standard_dev[j] + (aux[j] - best_params[j])*(aux[j] - best_params[j]);
			standard_dev_ctr[j]++;
		};
		
		likely_scatter_sort.push_back(aux);
		file_likely_scatter << aux[0] << '\t' << aux[1] << '\t' << aux[2] << '\t' << aux[3] << '\t' << aux[4] << endl; 
	};
	
	for (int j=0; j<param_number; j++) 
	{
		error_plus[j] = (double)error_plus[j]/error_plus_ctr[j];
		error_minus[j] = (double)error_minus[j]/error_minus_ctr[j];
		standard_dev[j] = sqrt((double)standard_dev[j]/standard_dev_ctr[j]);
		
		file_param_error << best_params[j] << '\t' << error_minus[j] << '\t' << error_plus[j] << endl;
		file_param_stddev << best_params[j] << '\t' << standard_dev[j] << endl;
	};
	
	
//	cerr << "likelysum = " << likely_sum << endl;
//	cerr << "likelies:" << endl;
//	for (int i=0; i<likelies_prep.size(); i++)                       // Debug !!!
//	{
//		cerr << likelies_prep[i][0] << '\t' << likelies_prep[i][1] << '\t' << likelies_prep[i][2] << '\t' << likelies_prep[i][3] << '\t' << likelies_prep[i][4] << '\t' << likelies_prep[i][5] << endl;
//	};
//	sort(likely_scatter_sort.begin(), likely_scatter_sort.end(),compare_likely);
//	likely_scatter_sort.resize(floor(siglevel*(double)likely_scatter_sort.size()));
//	vector<vector<double> > likelies_prep_sig = likely_scatter_sort;
	
	vector<vector<double> > likelies_prep_sig;
	sort(likelies_prep.begin(),likelies_prep.end(),compare_likely);
	//----------- Non-Bayesian ---------------------------
	//	for (int i=0; i<likelies_prep.size(); i++) 
	//{
	//	if (likelies_prep[i].back() > (1.0-siglevel2)) 
	//	{
	//		likelies_prep_sig.push_back(likelies_prep[i]);
	//	};
	//};
	//---------------------------------------------------
	

	
	int itr;
	double likely_count = 0;
//	//----------- Bayesian ---------------------------
	for (itr = 0; itr < likelies_prep.size() && likely_count < siglevel*likely_sum; itr++)
	{
		likely_count = likely_count + likelies_prep[itr].back();
		likelies_prep_sig.push_back(likelies_prep[itr]);
	};
	//---------------------------------------------------

//	cerr << "likelysum = " << likely_sum << endl;
//	cerr << "likely_count = " << likely_count << endl;
//	cerr << "siglevel*likely_sum = " << siglevel*likely_sum << endl;
//	cerr << "itr = " << itr << endl;
//	likelies_prep.resize(itr-1);
	
	cerr << "likelies after resizing:" << endl;
	for (int i=0; i<likelies_prep.size(); i++)                       // Debug !!!
	{
		cerr << likelies_prep[i][0] << '\t' << likelies_prep[i][1] << '\t' << likelies_prep[i][2] << '\t' << likelies_prep[i][3] << '\t' << likelies_prep[i][4] << '\t' << likelies_prep[i][5] << endl;
	};
	
	vector<double> max_sigparams(param_number,-1);
	vector<double> min_sigparams(param_number,10000);
//	for (int i=0; i < likelies_prep.size(); i++) 
//	{
//		for (int j=0; j<param_number; j++) 
//		{
//			if (likelies_prep[i][j] > max_sigparams[j]) max_sigparams[j] = likelies_prep[i][j];
//			if (likelies_prep[i][j] < min_sigparams[j]) min_sigparams[j] = likelies_prep[i][j];
//		};
//		
//	};
	
	for (int i=0; i < likelies_prep_sig.size(); i++) 
	{
		for (int j=0; j<param_number; j++) 
		{
			if (likelies_prep_sig[i][j] > max_sigparams[j]) max_sigparams[j] = likelies_prep_sig[i][j];
			if (likelies_prep_sig[i][j] < min_sigparams[j]) min_sigparams[j] = likelies_prep_sig[i][j];
		};
		
	};
	
	for (int j=0; j<param_number; j++) 
	{
		file_sigparams_limits << min_sigparams[j] << '\t' << max_sigparams[j] << endl;
	};
	
	file_likely_scatter.close();
	file_sigparams_limits.close();
	
	return best_params;
};











