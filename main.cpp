//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// Main Programm: reads parameters and calls MainWindow, which organizes the simulation
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
 

#include <iostream>

#include "global.h"

//#if Qt_flag==1
//#include <QApplication>
//#include <QtGui>
//#include <QLabel>
//#include <QPushButton>
//#include <QHBoxLayout>	
//#include <QSlider>
//#include <QDialog>
//#endif

//#include "MainWindow.h"

class MainWindow;


using namespace std;



int main (int argc, char *argv[]) 
{
	int seed;
	srand48(time(0));
	/////////// Assign parameters ///////////////////////////////
	// 1. Parameter (arv[0]) = program name
	// 2. Parameter (argv[1]) = reprod rate
	// 3. Parameter (argv[2]) = r
	// 4. Parameter (argv[3]) = delta
	// 5. Parameter (argv[4]) = strat rate
	// 6. Parameter (argv[5]) = random number seed
	// 7. Parameter (argv[6]) = number of runs
	// 8. Parameter (argv[7]) = option for file opening. 'a'=append, else delete
	// 9. Parameter (argv[8]) = run index for filenames

	if (argc==1)
	{
		seed = time(0);
		index_run = 0;
	//	parameter = 0;
	}
	else if (argc > 6) 
    {
		fprintf(stderr,"wrong number of arguments\n");
		exit(1);
    }
	else if (argc==2)
	{
		seed = time(0);
		reprod_rate0 = atof(argv[1]);
	}
	else if (argc==3)
	{
		seed = time(0);
		reprod_rate0 = atof(argv[1]);
		r0 = atof(argv[2]);
	}
	else if (argc==4)
	{
		seed = time(0);
		reprod_rate0 = atof(argv[1]);
		r0 = atof(argv[2]);
		delta0 = atof(argv[3]);
	}
	else if (argc==5)
	{
		
		seed = atoi(argv[5]);
		
		reprod_rate0 = atof(argv[1]);
		r0 = atof(argv[2]);
		delta0 = atof(argv[3]);
		strat_rate1 = atof(argv[4]);
		
	}
	else if (argc==6) 
	{
		seed = atoi(argv[5]);
		
		reprod_rate0 = atof(argv[1]);
		r0 = atof(argv[2]);
		delta0 = atof(argv[3]);
		strat_rate1 = atof(argv[4]);
		
		runnumber=atoi(argv[6]);
	}
	else if (argc==7) 
	{
		seed = atoi(argv[5]);
		
		reprod_rate0 = atof(argv[1]);
		r0 = atof(argv[2]);
		delta0 = atof(argv[3]);
		strat_rate1 = atof(argv[4]);
		
		runnumber=atoi(argv[6]);
		
		fileoption = *(argv[7]);
	}
	else if (argc==8) 
	{
		seed = atoi(argv[5]);
		
		reprod_rate0 = atof(argv[1]);
		r0 = atof(argv[2]);
		delta0 = atof(argv[3]);
		strat_rate1 = atof(argv[4]);
		
		runnumber=atoi(argv[6]);
		
		fileoption = *(argv[7]);
		
		index_run = atoi(argv[8]);
	}

	//////////////////////////////////////////////////////////////
	
	printf("seed=%d\n",seed);
	//printf("actindens0=%f\n",actindens);
	cout << parameter_str << "=" << parameter << endl;
	
	
	/////////// Setup Random Number Generators ///////////////////////
	srand48(seed);
	ibm=2*seed+1;
	ibm=ibm*16807;
	ibm=ibm*16807;
	///////////////////////////////

#if fitbasal == 1
	if (param_number==5) 
	{
		var_vector[0]=1; var_vector[1]=1; var_vector[2]=0; var_vector[3]=1; var_vector[4]=0;
	}
	else {cerr << "error: wrong number of parameters\n", exit(1);};

#else
	var_vector[0]=1; var_vector[1]=1; var_vector[2]=0; var_vector[3]=1; var_vector[4]=1;
#endif
		
#if param_fixed == 1
	var_vector[0]=0;
#elif param_fixed == 2
	var_vector[1]=0;
#elif param_fixed == 3
	var_vector[2]=0;
#elif param_fixed == 4
	var_vector[3]=0;
#elif param_fixed == 5
	var_vector[4]=0;
#endif
			
	//////////// Call MainWindow to make Graphical Interface and/or start Simulation /////////////////////////////
#if Qt_flag==1
	QApplication app(argc, argv); 
#endif
	
	window = new MainWindow;    
	
#if Qt_flag==1
	window->show();     // Graphical Interface; simulation has to be started manually
	return app.exec();
#else
	window->start_sim();  // If Qt off, simulation has to be started automatically
#endif
			
			
			
}
