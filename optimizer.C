#include <cmath>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <complex>
#include <unistd.h>
#include "TF1.h"
#include "TGraphErrors.h"

typedef struct{
	int EpicsN;
	float val;
	float orig;
	float max;
	float min;
	float stp;
} Element;

typedef struct{
	float val;
	int Nions;
} scan;

typedef struct{
	float ToF;
	float Sig;
	float Nions;
} meas;

#define MAX_elements 300

Element elements[MAX_elements];
int Nelements;
int Nsteps;
int Run;
int ResOnOff;
float scalestep;
float RunWaitTime;// in seconds
int RunDate;
ofstream LogFile("OptimizerLog.txt");

void wait_seconds(int time_s);
void wait_minutes(int time_m);
void define_element(int EpicsN, float orig, float max, float min, float stp);
void set_element(int el, float val);
void scan_element(int el, char mode);
meas run_and_read();
void wait_for_run();
void restore_elements();
void TurnQuadAFG(int onoff);
void CleanEvaFile(int RunToClean);
float OptimizerCriterea(float Counts, float ToFeff);
void SetTune(int tune);
void PearlRCstart();

void optimizer(){
	
	int el, round, naux;
	meas Measurement;
	Nelements = 0;
	

	//  Settings:
	////////////////////////
	int RunStart = 411826;
	RunWaitTime = 20.0;// seconds
	Nsteps = 3; // total number of steps will be (2*Nsteps+1)
	ResOnOff = 0; // 2D analysis On/Off resonance? (1) yes, (0) no
	scalestep = 1.0; // globally manipulate step size for all devices
	////////////////////////	


	//  List of variables:
 	//  Syntax: define_element(int EpicsN, float orig, float max, float min, float stp);
	////////////////////////
	define_element(56, 17812, 20000, 17000, 3); //				TRFC:PB5 (V)	
	define_element(57, -2820, -200.0, -5000.0, 150.0); //			TRFC:EL5 (V)	
	define_element(49, 308, 400, 50.0, 0.6); //				TRFCBL:B1:POS
 	define_element(52, 13, 100.0, -100.0, 3.0); //				TRFCBL:Q1:A (V)
	define_element(51, 230, 500.0, 50.0, 10.0); //				TRFCBL:Q1:C
	define_element(53, 82.5, 115, 80, 3.0); //				TRFCBL:XCB4 (V)
	define_element(54, 200, 225, 195, 3.0); //				TRFCBL:YCB4 (V)
	define_element(39, 2600, 4000.0, 500.0, 50.0); //			TSYBL:EL1 (V)
	define_element(45, 840, 4000.0, 500.0, 100.0); //			TSYBL:EL3 (V)
	define_element(67, 260, 260.0, 240.0, 4.0); //				TSYBL:XCB9 (V)
	define_element(68, 260, 260.0, 240.0, 4.0); //				TSYBL:YCB9 (V)
	define_element(80, 1130, 4000.0, 500.0, 50.0); // 			TSYBL:EL4 (V)
	define_element(34, 1750, 4000.0, 500.0, 50.0); // 			MPETBL:EL2 (V)
	define_element(36, 235, 265.0, 235.0, 2.5); //				MPETBL:XCB3 (V)
	define_element(37, 254, 265.0, 235.0, 2.5); //				MPETBL:YCB3 (V)
	define_element(25, 2100, 3000, 0.0, 3.0); // 				MPET:PLTP (V)
	define_element(24, -290, 0.0, -500, 3.0); // 				MPET:PLTM (V)	
	////////////////////////


	restore_elements();
	Run = RunStart;


	//  Scan Strategy:
	////////////////////////

	naux  = run_and_read().Nions;
	LogFile << "  Initial counts: " << naux <<endl<<endl;

	for(round=0;round<100;round++){

		LogFile << "----------------------------------------------------------------------------------------------------" <<endl;
		LogFile << "  ROUND " << round+1 <<endl;
		LogFile << "  Initial counts: " << naux <<endl<<endl;
		
		for(el=0;el<Nelements;el++) 
			scan_element(el,'o');

		naux  = run_and_read().Nions;
		LogFile <<endl<< "  Final counts: " << naux <<endl;
	}
}


void scan_element(int el, char mode){
	
	float Val[50], Nion[50], NionErr[50], ErrZero[50], ToF[50], ToFSig[50], Nion2[50], NionErr2[50], ToF2[50], ToFSig2[50], ToFeff[50], ToFeffErr[50], ToFErr[50], ToFErr2[50], Opt[50], OptErr[50];
	meas Measurement;
	int i;
	float A,B,C,Chisqparabolafit,OptVal, OptCount;
	TF1 *func;
	TGraphErrors *gr;

	LogFile << "====================================================" << endl;
	LogFile << "== Device " << elements[el].EpicsN << ": \t ";
	LogFile << "start at run " << Run << endl;

	
	for(i=0;i<(Nsteps*2 + 1);i++){
		Val[i] = elements[el].val + ((float)i-Nsteps)*(elements[el].stp);
		if(Val[i] > elements[el].max) Val[i] = elements[el].max;
		else if(Val[i] < elements[el].min) Val[i] = elements[el].min;
		
		if(ResOnOff) TurnQuadAFG(0);
		set_element(el,Val[i]);
		Measurement = run_and_read();
		Nion[i] = Measurement.Nions;
		NionErr[i] = sqrt(Nion[i]);
		ToF[i] = Measurement.ToF;
		ToFSig[i] = Measurement.Sig;
		if(NionErr[i]>0.0) ToFErr[i] = ToFSig[i]/NionErr[i];
		ErrZero[i] = 0.0;
		if(ResOnOff){ 
			TurnQuadAFG(1);		
			Measurement = run_and_read();
			Nion2[i] = Measurement.Nions;
			NionErr2[i] = sqrt(Nion2[i]);
			ToF2[i] = Measurement.ToF;
			ToFSig2[i] = Measurement.Sig;
			if(NionErr2[i]>0.0) ToFErr2[i] = ToFSig2[i]/NionErr2[i];
			if(Nion[i]>1.0 && Nion2[i]>1.0 && ToF[i]>0.0 && ToF2[i]>0.0 && ToFSig[i]>0.0 && ToFSig2[i]>0.0){
				ToFeff[i] = 1.0- ToF2[i]/ToF[i];
				ToFeffErr[i] = ToFeff[i]*sqrt( (ToFSig[i]/(NionErr[i]*ToF[i]))*(ToFSig[i]/(NionErr[i]*ToF[i])) + (ToFSig2[i]/(NionErr2[i]*ToF2[i]))*(ToFSig2[i]/(NionErr2[i]*ToF2[i]))     ) ;
				if(mode=='o'){
					Opt[i] = OptimizerCriterea((Nion2[i]+Nion[i]),ToFeff[i]);
					OptErr[i] = abs(Opt[i]*( 1.0/sqrt(Nion2[i]+Nion[i])   ));
				}			

			} else {ToFeff[i] = ToFeffErr[i] = 0.0;}


		} else if(mode=='o'){
					Opt[i] = OptimizerCriterea(Nion[i],1.0);
					OptErr[i] = abs(Opt[i]*( 1.0/sqrt(Nion[i])   ));
		}

		if(Nion[i]==0.0) Nion[i]=NionErr[i] = 1.0; // ugly, but that's how Ill handle zeros for now
		
		LogFile << Val[i] << "\t" <<  Nion[i]  << "\t" << NionErr[i]  << "\t" <<  ToF[i] << "\t" << ToFErr[i] << "\t" <<  ToFSig[i] << "\t";
		if(ResOnOff) LogFile <<   Nion2[i]  << "\t" << NionErr2[i]  << "\t" <<  ToF2[i] << "\t"  << ToFErr2[i]  << "\t" <<  ToFSig2[i] << "\t" << ToFeff[i] << "\t" <<  ToFeffErr[i] << "\t" ;
		LogFile << endl;
	}
		
	if(mode=='o'){
		func = new TF1("FitParabola","pol2",Val[0],Val[Nsteps*2]);
		gr = new TGraphErrors((Nsteps*2 + 1),Val,Opt,ErrZero,OptErr);
		gr->Fit("FitParabola", "Q");
			
		Chisqparabolafit = func->GetChisquare();
		A = func->GetParameter(2); 
		B = func->GetParameter(1); 
		C = func->GetParameter(0); 
		OptVal = -0.5*B/A;
		OptCount = A*OptVal*OptVal + B*OptVal +C;
		
		LogFile << " fit results {A,B,C,Chi2,OptVal,OptCount} = {" << A << " , " << B << " , "<< C << " , "<< Chisqparabolafit << " , "<< OptVal << " , "<< OptCount <<  " } \t ";
		
		if(A<0 && Chisqparabolafit < 25.0){ // remains unchanged if criteria not satisfied. Ensure we look to a peak (A<0) and fit was more or less reasonable (chi2<20)
			if(OptVal > Val[Nsteps*2]) elements[el].val = Val[Nsteps*2];
			else if(OptVal < Val[0]) elements[el].val = Val[0];
			else elements[el].val = OptVal;
		}
		else if(A>0){
			for(i=0;i<(Nsteps*2 + 1);i++)
			if(Opt[i] > OptCount){  OptCount = Opt[i]; OptVal = Val[i]; elements[el].val = OptVal;}
		}
		
		LogFile << " taking " << elements[el].val << " V" << endl;
	}

	set_element(el, elements[el].val);
	
}

void set_element(int el, float val){
	char command[500];

	sprintf(command, "/home/mpet/packages/midas/linux/bin/odbedit -e mpet -c \"set /Equipment/Beamline/Variables/Demand[%d] %.1f \";", elements[el].EpicsN, val);

	gSystem->Exec(command);
}


void define_element(int EpicsN, float orig, float max, float min, float stp){
	
	elements[Nelements].EpicsN = EpicsN;
	elements[Nelements].val = elements[Nelements].orig = orig;
	elements[Nelements].max = max;
	elements[Nelements].min = min;
	elements[Nelements].stp = stp*scalestep;
	
	Nelements++;
}

void wait_for_run(){
	usleep(1000000*RunWaitTime); // RunWaitTime in seconds
}
void wait_seconds(int time_s){
	int i;
	for(i=0; i<time_s; i++)
		usleep(1000000);
}
void wait_minutes(int time_m){
	int i;
	for(i=0; i<time_m; i++)
		wait_seconds(60);
}

meas run_and_read(){
	
	meas Measurement;
	char command[500];
	char line[255];
	float auxN, auxToF, auxSig;
	
	sprintf(command, "/home/mpet/packages/midas/linux/bin/odbedit -e mpet -c \"start %d now\";", Run);
	gSystem->Exec(command);
	wait_for_run();
	//sprintf(command, "/home/mpet/local/scripts/m2e.py /titan/data1/mpet/%d/run%d.mid ;", RunDate,Run);
	//gSystem->Exec(command);
	//sprintf(command, "/home/mpet/Rene/eva_read_counts.py /titan/data1/mpet/evafiles/run%d_eva.dat > auxout.dat ;", Run);
	sprintf(command, "/home/mpet/Rene/eva_read_counts_all.py /titan/data1/mpet/evafiles/run%d_eva.dat > auxout.dat ;", Run);
	gSystem->Exec(command);	
	ifstream result("auxout.dat");
    result.getline(line,255);
    sscanf(line, "%f %f %f", &auxN, &auxToF, &auxSig);
	result.close(); 
	 CleanEvaFile(Run);
	Run++;

	Measurement.Nions = auxN;
	Measurement.ToF = auxToF;
	Measurement.Sig = auxSig;

	
	return Measurement;
}


void restore_elements(){
	int el;
	
	for(el=0;el<Nelements;el++)	
		set_element(el,elements[el].orig);
}


void SetTune(int tune){
	char command[500];
	
	if(tune==39) sprintf(command, "/home/mpet/local/scripts/perlrc.sh tune K39");
	else if(tune==41) sprintf(command, "/home/mpet/local/scripts/perlrc.sh tune K41"); 
	else if(tune==87) sprintf(command, "/home/mpet/local/scripts/perlrc.sh tune Rb87"); 
	else if(tune==85) sprintf(command, "/home/mpet/local/scripts/perlrc.sh tune Rb85"); 
	else if(tune==1) sprintf(command, "/home/mpet/local/scripts/perlrc.sh tune ResMode"); 
	
	gSystem->Exec(command);
	wait_seconds(10);
	

}

void PearlRCstart(){
	char command[500];
	sprintf(command, "/home/mpet/local/scripts/perlrc.sh start &");
	gSystem->Exec(command);
}


void TurnQuadAFG(int onoff){
	char command[500];
	sprintf(command, "/home/mpet/packages/midas/linux/bin/odbedit -e mpet -c \"set '/Experiment/Variables/AFG On' '%d' \";", onoff);
	gSystem->Exec(command);
}

void CleanEvaFile(int RunToClean){
	char command[500];
	sprintf(command, "rm /titan/data1/mpet/evafiles/run%d_eva.dat ;", RunToClean);
	gSystem->Exec(command);	
}

float OptimizerCriterea(float Counts, float ToFeff){

	if(ToFeff<0.01) ToFeff=0.01;	//0.01 is too small to be relevant and may introduce too much weight to ToFeff
	
	return (Counts/(ToFeff));
}
