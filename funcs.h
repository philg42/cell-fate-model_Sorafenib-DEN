/*
 *  funcs.h
 *  popdyn_1D_discr
 *
 *  Created by Philip Greulich on 08.10.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _FUNCS_
#define _FUNCS_

#include "global.h"
#include <string>

using namespace std;

string itoa(int n);

double sgn (double x);

vector<ofstream*> make_filelist(string& prefix, int filenumber);

ofstream* make_file(string& prefix, int number1, int number2=-1, int number3=-1, char fileoption0='n');

vector<vector<double> > read_file(ifstream* file);
vector<vector<int> > read_file_int(ifstream* file);

vector<vector<double> > random_parameters();

bool compare_joint(double x, double y);
bool compare_basal(double x, double y);
bool compare_total(double x, double y);
bool compare_cumul_basal(double x, double y);
bool compare_cumul_total(double x, double y);
bool compare_floating(double x, double y);
bool compare_likely(vector<double> x, vector<double> y);
bool compare_pairs(vector<double> x, vector<double> y);

int start_val(int i);

vector<vector<double> > make_switchparams(const vector<vector<double> >& parameters);

vector<vector<double> > bin(const vector<vector<double> >& CSD, char mode='l', int binsize=1);
vector<double> bin(const vector<double>& CSD, char mode='l', int binsize=1);
vector<vector<double> > bin_cont(vector<double> sizes, double binsize_cont);

vector<vector<vector<double> > > read_filesmbl(vector<ifstream*>& files);
vector<ifstream*> make_filelist(string& prefix, vector<vector<int> >& filenumbers);

double prior(double reprodrate0, double r0, double stratrate0, double shedrate, double delta);    // Attention !!! Adjust for each individual BAyesian prior function

vector<double> eval_loglikelyfiles(vector<ifstream*>& file_loglikelies_in, double delta_intv, double reprod_intv, double r_intv, double strat_intv, double shed_intv);









#endif
