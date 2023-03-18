/*:
COMPUTATION OF THE MASS-RADIUS DIAGRAM AND THE LOVE NUMBER K2 FOR A NEUTRON STAR WITH BL EQUATION OF STATE
*/

//Libraries
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//def pi
#define PI 3.141592653589

//Dimensionless quantity
double r; //radius (integration variable)
double p; //pressure at distance r from the centre of the star
double *pPtr=&p; //pressur poynter
double m; //mass enclosed in a shell of radius r
double *mPtr=&m; // mass poynter
double y; //function to integrate for the calculation of Love number k2
double *yPtr=&y; //y function poynter
double eps(double p); // function that returns the energy density at a given pressure
double depsdp(double p); //function that returns the derivative of energy density respect to pressure at a given pressure


/*
Integration method: Runge-Kutta IV order 
dp/dr=f(r,p,m,y)
dm/dr=g(r,p,m,y)
dy/dr=u(r,p,m,y)
NB: f and g are given by TOV equations, u given by Love coefficient's theory
*/

//declaration of function's parameters
double fargs[]={};  
double *fargsPtr=&fargs[0]; 
double gargs[]={};  
double *gargsPtr=&gargs[0]; 
double uargs[]={};  
double *uargsPtr=&uargs[0]; 

//declaration of functions
double f(double r, double p, double m, double y, double *fargPtr);
double g(double r, double p, double m, double y, double *gargsPtr);
double u(double r, double p, double m, double y, double *uargsPtr);


//declaration of Runge-Kutta IV order method
int RK4(double r, double *pPtr, double *mPtr, double *yPtr, double h, double(*f)(double, double, double, double, double*), double(*g)(double, double, double, double, double*), double(*u)(double, double, double, double, double*));

//declaration of functions that enter in the calculation of y
double F(double r, double mcorr, double pcorr, double ycorr); 
double Q(double r, double mcorr, double pcorr, double ycorr); 

//declaration of function that calculates the love coefficient of a star
double love(double C, double YR);

//declaration of functions and arrays for the input of the EOS
#define LENGTH 200 //lenght of arrays pressure and energy (from EOS)
double pressure_array[LENGTH];
double log_pressure_array[LENGTH];
double energy_array[LENGTH];
double log_energy_array[LENGTH];

int read_pressure();
int read_energy();


int main(){

	//read of file in input (pressure and energy density)
	read_pressure();
	read_energy();
	
	//dimensional parameters
	double c=299792458; // [m/s]
	double G=6.6725e-11; // [Nm^2/kg]
	double Ms=1.989e30; // [kg] Solar mass
	double r0=10000.0; // [m]
	double m0=c*c*r0/G; //[kg]
	double eps0=m0*c*c/(r0*r0*r0); //[J/m^3]
	double p0=eps0; //[J/m^3]
	
	//Integration parameters
	double h=1e-3; //integration step
	double minimun=1e-17; // minimum to controll when p\simeq 0
	
	//File of output
	FILE * NSptr;
	NSptr=fopen("NS_output.txt", "w");

	//Print
	fprintf(NSptr, "#Pc[J/m^3] \t \t R[km] \t \t M[Ms] \t \t \t C \t \t \t Y(R) \t \t \t k2 \n");
	
	for(int n=0; n<=120; n++){  // cycle that vary the boundary conditions
		//boundary condition (at the center of the star)
		r=0.0; // radius at the center of the star
		m=0.0; // mass at the center of the star
		p=1e-1/(pow(1.05, n)); //variation of the central pressure at every cycle 
		y=2.0; //boundary condition y function
		
		fprintf(NSptr, "%10.8e \t", p*p0); //print
		
		//variables to save final value of parameters
		double rf=0;
		double pf=0;
		double mf=0;
		double yf=0;
	
		//Integration
		while(p>minimun){
			r=r+h; //step
			rf=r; 
			pf=p;
			mf=m;
			yf=y;
			RK4(r, pPtr, mPtr, yPtr, h, f, g, u); //integration step on p, m, y		
		}
		
		//Results extrapolation
		double R=rf*r0*1e-3; //Radius in km
		double M=mf*m0/Ms; //Mass in Solar masses
		double YR=yf; //Y value at the radius star R
		double C=mf/rf; //Compactness
		double k2=love(C, YR); //Love number k2
		
		//Print of results
		fprintf(NSptr, "%10.8e \t %10.8e \t %10.8e \t %10.8e \t %10.8e \n", R, M, C, YR, k2);	
	}
	
	//file closure
	fclose(NSptr); 

return 0;	
}

//def Runge-Kutta IV: makes a step h of integration
int RK4(double r, double *pPtr, double *mPtr, double *yPtr, double h, double(*f)(double, double, double, double, double*), double(*g)(double, double, double, double, double*), double(*u)(double, double, double, double, double*)){
	//corrent variables
	double pcorr=*pPtr; 
	double mcorr=*mPtr;
	double ycorr=*yPtr;
	//first step of the method
	double k1=h*f(r, pcorr, mcorr, ycorr, fargsPtr);
	double l1=h*g(r, pcorr, mcorr, ycorr, gargsPtr);
	double d1=h*u(r, pcorr, mcorr, ycorr, uargsPtr);
	//second step of the method
	double k2=h*f(r+h/2, pcorr+k1/2, mcorr+l1/2, ycorr+d1/2, fargsPtr);
	double l2=h*g(r+h/2, pcorr+k1/2, mcorr+l1/2, ycorr+d1/2, gargsPtr);
	double d2=h*u(r+h/2, pcorr+k1/2, mcorr+l1/2, ycorr+d1/2, uargsPtr);
	//third step of the method
	double k3=h*f(r+h/2, pcorr+k2/2, mcorr+l2/2, ycorr+d2/2, fargsPtr);
	double l3=h*g(r+h/2, pcorr+k2/2, mcorr+l2/2, ycorr+d2/2, gargsPtr);
	double d3=h*u(r+h/2, pcorr+k2/2, mcorr+l2/2, ycorr+d2/2, uargsPtr);
	//fourth step of the method
	double k4=h*f(r+h, pcorr+k3, mcorr+l3, ycorr+d3, fargsPtr);
	double l4=h*g(r+h, pcorr+k3, mcorr+l3, ycorr+d3, gargsPtr);
	double d4=h*u(r+h, pcorr+k3, mcorr+l3, ycorr+d3, uargsPtr);
	
	//Calculation of p(r+h), m(r+h), y(r+h)
	*pPtr=pcorr+(k1+2*k2+2*k3+k4)/6.0;
	*mPtr=mcorr+(l1+2*l2+2*l3+l4)/6.0;
	*yPtr=ycorr+(d1+2*d2+2*d3+d4)/6.0;
return 0;
}


//def of functions that enter in RK4: from TOV and Love number's theory

double f(double rl, double pl, double ml, double yl, double *fargsPtr){
	double out=0;
	out=-(pl+eps(pl))*(ml+4.0*PI*rl*rl*rl*pl)/(rl*rl-2.0*ml*rl);
return out;
}

double g(double rl, double pl, double ml, double yl, double *gargsPtr){
	double out=0;
	out=4.0*PI*rl*rl*eps(pl);
return out;
}

double u(double rl, double pl, double ml, double yl, double *uargsPtr){
	double out=0;
	out=-yl*yl/(rl)-F(rl, pl, ml, yl)*yl/(rl)-rl*Q(rl, pl, ml, yl);
return out;
}


//def of function eps: gives the energy density at a given pressure (from linear interpolation of EOS)
double eps(double pl){
	double out=0;

//Linear interpolation
	double logpl=log10(pl);
	int nn=0;
	if(logpl<log_pressure_array[0]){
		printf("WARNING: pressure value to low %10.8e \t %10.8e\n", logpl, log_pressure_array[0]);
		nn=0;
	}
	if(logpl>log_pressure_array[LENGTH-1]){
		printf("WARNING: pressure value to high %10.8e \t %10.8e\n", logpl, log_pressure_array[LENGTH-1]);
		nn=LENGTH-2;
	}
	
	while(((logpl<log_pressure_array[nn])||(logpl>log_pressure_array[nn+1]))){
		nn++;
	}
	double alpha=(logpl-log_pressure_array[nn])/(log_pressure_array[nn+1]-log_pressure_array[nn]);
	out=log_energy_array[nn]*(1-alpha)+alpha*log_energy_array[nn+1];
	out=pow(10.0, out);	
return out;
}

//def of the derivative of energy density repsect to pressure at a given pressure (from linear interpolation of EOS)
double depsdp(double pl){
	double out=0;
	
//Linear interpolation
	double logpl=log10(pl);
	int nn=0;
	if(logpl<log_pressure_array[0]){
		printf("WARNING: pressure value to low \n");
		nn=0;
	}
	if(logpl>log_pressure_array[LENGTH-1]){
		printf("WARNING: pressure value to high \n");
		nn=LENGTH-2;
	}
	while(((logpl<log_pressure_array[nn])||(logpl>log_pressure_array[nn+1]))){
		nn++;
	}
	out=(energy_array[nn+1]-energy_array[nn])/(pressure_array[nn+1]-pressure_array[nn]);
return out;
}


//def of functions F and Q that enters in the calculation of y

double F(double rl, double pcorr, double mcorr, double ycorr){
	double out=0;
	out=(rl-4.0*PI*rl*rl*rl*(eps(pcorr)-pcorr))/(rl-2.0*mcorr);
return out;	
}

double Q(double rl, double pcorr, double mcorr, double ycorr){
	double out=0;
	out=4.0*PI*rl*(5.0*eps(pcorr)+9.0*pcorr+(eps(pcorr)+pcorr)*depsdp(pcorr)-6.0/(4.0*PI*rl*rl))/(rl-2.0*mcorr)-4.0*(mcorr+4.0*PI*pcorr*rl*rl*rl)*(mcorr+4.0*PI*pcorr*rl*rl*rl)/((rl*rl-2*mcorr*rl)*(rl*rl-2*mcorr*rl));
	return out;	
}


//Def of Love function: give the Love number K2 of a star
double love(double C, double YR){
	double out=0;
	double out1=8.0*C*C*C*C*C/5.0*(1.0-2.0*C)*(1.0-2.0*C);
	double out2=2.0+2.0*C*(YR-1.0)-YR;
	double out3=2.0*C*(6.0-3.0*YR+3.0*C*(5.0*YR-8.0));
	double out4=4.0*C*C*C*(13.0-11.0*YR+C*(3.0*YR-2.0)+2*C*C*(1+YR));
	double out5=3.0*(1.0-2.0*C)*(1.0-2.0*C)*(2.0-YR+2.0*C*(YR-1.0))*log(1.0-2.0*C);
	out=out1*out2/(out3+out4+out5);
return out;
}

//def of the function that reads in input the pressure from EOS
int read_pressure(){
	FILE * pressurePtr;
	pressurePtr=fopen("pressure.txt", "r");
	for(int ll=0; ll<LENGTH; ll++){
		fscanf(pressurePtr, "%lf", &pressure_array[ll]);
		log_pressure_array[ll]=log10(pressure_array[ll]);
	}
	fclose(pressurePtr);
return 0;
}

//def of the function that reads in input the energy density from EOS
int read_energy(){
	FILE * energyPtr;
	energyPtr=fopen("energy.txt", "r");
	for(int ll=0; ll<LENGTH; ll++){
		fscanf(energyPtr, "%lf", &energy_array[ll]);
		log_energy_array[ll]=log10(energy_array[ll]);
	}
	fclose(energyPtr);
}

