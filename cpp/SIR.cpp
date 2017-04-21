#include<stdio.h>
#include<iostream>
#include<fstream>
#include<conio.h>

double beta = 0.001;
double gamma = 0.1;
double tstep, nsteps;

	void rk4(double *ss, double *ii, double *rr, double *tt);
	double derivS(double tdummy, double sdummy, double idummy, double rdummy);
	double derivI(double tdummy, double sdummy, double idummy, double rdummy);
	double derivR(double tdummy, double sdummy, double idummy, double rdummy);

using namespace std;


 void main(){
	int j;
	double ss = 190;
	double ii = 10;
	double rr = 0;
	double tt = 0;
  
	  tstep = 0.1;
	  nsteps = 1000;
  
		
	//write file
	ofstream myfile ("data_SIR.csv");
	  if (myfile.is_open())
	  {
		//header
		myfile << "time, S, I, R\n";
		
		//Loop #####
			//call rk4 function
			for(j = 1; j <= nsteps; j++){
				 tt = j*tstep;
				
				//call rk4
				rk4(&ss, &ii, &rr, &tt);
				//write
				myfile << tt <<"," << ss << "," << ii << "," << rr <<"\n";
			}
		
		//done
			myfile.close();
			cout << "Output csv completed..." << endl;
			
	  }
	  else
	  {	  
		  cout << "Unable to open file";
	  }
	  _getch();
	  
 }
 
void rk4(double *ss, double *ii, double *rr, double *tt){
  
  double s, i, r, t;
  double h;
  double s1, s2, s3, s4;
  double i1, i2, i3, i4;
  double r1, r2, r3, r4;

  s = *ss;
  i = *ii;
  r = *rr;
  t = *tt;
  
  h = tstep/2.0;
  
  s1 = tstep * derivS(t, s, i, r);
  i1 = tstep * derivI(t, s, i, r);
  r1 = tstep * derivR(t, s, i, r);
  
  s2 = tstep * derivS(t+h, s+s1/2.0, i+i1/2.0, r+r1/2.0);
  i2 = tstep * derivI(t+h, s+s1/2.0, i+i1/2.0, r+r1/2.0);
  r2 = tstep * derivR(t+h, s+s1/2.0, i+i1/2.0, r+r1/2.0);
  
  s3 = tstep * derivS(t+h, s+s2/2.0, i+i2/2.0, r+r2/2.0);
  i3 = tstep * derivI(t+h, s+s2/2.0, i+i2/2.0, r+r2/2.0);
  r3 = tstep * derivR(t+h, s+s2/2.0, i+i2/2.0, r+r2/2.0);
  
  s4 = tstep * derivS(t+tstep, s+s3, i+i3, r+r3);
  i4 = tstep * derivI(t+tstep, s+s3, i+i3, r+r3);
  r4 = tstep * derivR(t+tstep, s+s3, i+i3, r+r3);

  s = s + (s1 + (2.0*(s2 + s3)) + s4)/6.0;
  i = i + (i1 + (2.0*(i2 + i3)) + i4)/6.0;
  r = r + (r1 + (2.0*(r2 + r3)) + r4)/6.0;
  
  *ss = s;
  *ii = i;
  *rr = r;
  *tt = t;
  
}

double derivS(double tdummy, double sdummy, double idummy, double rdummy){
	return(-1.0*beta*sdummy*idummy);
  
}


double derivI(double tdummy, double sdummy, double idummy, double rdummy){
	return(beta*sdummy*idummy - gamma*idummy);
}

double derivR(double tdummy, double sdummy, double idummy, double rdummy){
	return(gamma*idummy);
}
