/*
 *  MainWindow.h
 *  popdyn_1D_discr
 *
 *  Created by Philip Greulich on 07.10.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MAINWINDOW_H_
#define MAINWINDOW_H_

#include "global.h"

#include <iostream>

#if Qt_flag==1
#include <QtGui>
#include <QTimer>
#include <QSlider>
#include "PaintWindow.h"
#include "LabeledSlider.h"
#include "LabeledDisplay.h"
#endif

#include <vector>
#include <list>
#include <set>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;



//#include "func_per.h"


#define Pi 3.141593

using namespace std;

#if Qt_flag==1
class PaintWindow;
class LabeledSlider;
class LabeledDisplay;
#endif

class run;

#if Qt_flag==1
class MainWindow : public QDialog         // Main interface for simulation; manages all objects and functions
{
	Q_OBJECT
	
#else
	
	class MainWindow
	{
		
#endif
		
	public:
		
#if Qt_flag==1
		MainWindow(QWidget *parent = 0);
#else
		MainWindow();
#endif
		
		friend class PaintWindow;
		
		void start_sim();                // Start simulation manually
		void read_files_to_data();
		void prepare_data();
		vector<vector<double> > prepare_parameters(vector<vector<double> >& parameters0);
		vector<vector<double> > get_best_params(vector<vector<double> > parameters0);
		
		run* fit_evol(vector<vector<double> > param_init, const vector<vector<vector<double> > >& data0);   // Fit "data0" starting with parameters "param_init"
		void scan_likelihood(vector<vector<double> > param_init, const vector<vector<vector<double> > >& data0); // Scan parameter space for likelihoods
		void find_errorbars(run* maxrun, const vector<vector<vector<double > > >& data0);     // Find error bars and p-valus respective data, taking "maxrun" as test function
		
#if Qt_flag==1
		private slots:                
#endif
		
		void start();                  // Start simulation
		void finish();                 // Finish simulation
		
#if Qt_flag==1
		void stop();                 // Pause simulation
		void continue_sim();         // dito
		void start_video();          // Start making a video
		void stop_video();           // Stop making a video
#endif
		
		void print_data(run* run0);               // dito
		
		
		
		////////// Set parameters by sliders ////////////////
		void set_reprodrate0(int rate0_int);
		void set_reprodrate1(int rate0_int);
		void set_r0(int rate0_int);
		void set_r1(int rate0_int);
		void set_dediffrate(int rate0_int);
		void set_delta0(int rate0_int);
		void set_delta1(int rate0_int);
		void set_shedrate(int rate0_int);
		void set_stratrate(int rate0_int);
		void set_directstratfrac(int rate0_int);

		///////////////////////////////////////////
		
		
		
	private:
		
		
//		void setup_sim();            // Setup and initiate objects and files
//		void reset_sim();
		
#if Qt_flag==1
		
		void setup_widgets();         // Setup graphical interface
		
		QTimer *timer;                // Sends periodic signal for system update
		QTimer *output_timer;         // Sends periodic signal for animation update
		QTimer *image_timer;          // Sends periodic signal for video update (add image)
		QLabel *label;                
		QPushButton *startbutton;     // Button to start simulation
		QPushButton *stopbutton;     // Button to pause simulation
		QPushButton *continuebutton;     // Button to continue simulation
		QPushButton *resetbutton;     // Button to reset simulation
		QPushButton *quitbutton;     // Button to finish simulation
		QPushButton *snapshotbutton;     // Button to make snapshot image
		QPushButton *startvideobutton;     // Button to start making video (image sequence) 
		QPushButton *stopvideobutton;     // Button to stop making video
		

		
		
		//QSlider *slider;
		QSpinBox *spinbox;
		PaintWindow *output;            // panel to display system configurations
		QLabel *sliderlabel;
		
		vector<LabeledSlider*> sliders;		// Sliders to tune parameters

		
		vector<LabeledDisplay*> displays;	// Displays
		
		LabeledDisplay* display_counter;
	


		
#endif
		/////////// Discretization steps for sliders ////////////////////
		double minstep_reprod;                        // Attention: new ones might need to be added for more processes
		double minstep_r;
		double minstep_dediff;
		double minstep_delta;
		double minstep_shed;
		double minstep_strat;
		double minstep_direct_strat;
		
	    /////////////// Files ////////////////////////////////////////////	




	};

#endif

