
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "particle.h"
#include "parameters.h"
#include "solve_growth.h"
#include "msg.h"

//double gb_om = 0.3071, gb_oml=0.6929, gb_w=-1.5, gb_H0=70;
//double gb_om = gb_om, gb_w = gb_w, gb_oml = gb_oml, gb_H0=70.;
double gb_om,gb_w, gb_oml, gb_H0;
int gb_use_solve_growth;
double D_final, E_final, Dp_final, Ep_final;
//double gb_H0 = 70.;


int solvegrowth_flag=0;


#define GB_NUM_INTPL  3000
double const gb_basenumber = 1. + 1./512.;
double const gb_logbasenumber = log(1. + 1./512.);
//#define GB_NUM_INTPL  6000
//double const gb_basenumber = 1. + 1./1024.;
//double const gb_logbasenumber = log(1. + 1./1024.);

double gb_maxintplz, gb_minintpla;

double gb_zdata[GB_NUM_INTPL], gb_Ddata[GB_NUM_INTPL], gb_Edata[GB_NUM_INTPL], gb_QDdata[GB_NUM_INTPL], gb_QEdata[GB_NUM_INTPL];;

double gb_zi(int i){
	return pow(gb_basenumber,(double)(i-1)) -1.;
}

// compute T^2 D1 = Q (d/da) ( Q (d/da) D1 ) = Q dQ/da dD1/da + Q^2  (d/da)^2 D1
//
//
double Ha(double a){
	return gb_H0 * pow( gb_om*pow(a,-3.0) + gb_oml*pow(a,-3.0*(1.0+gb_w)) + (1.-gb_om-gb_oml)/(a*a), 0.5);
	//return gb_H0 * pow( gb_om*pow(a,-3.0) + gb_oml*pow(a,-3.0*(1.0+gb_w) + (1.-gb_om-gb_oml)/(a*a)), 0.5);
}
double dHada(double a){
	double y, y1, y2;
//	y = pow(a,3.*gb_w);
//	y1 = 3. * pow(a,-4.) / y * gb_H0 * (1.+gb_om*(-1.+y-gb_w)+gb_w);
//	y2 = pow(a,-3.)/y*(1.+(-1.+y)*gb_om);
//	y2 = 2. * pow(y2,0.5);
//	return - y1 / y2;

//	y = a*a*a;
//	y1 = -3.*gb_H0*gb_om + 2.*a*gb_H0*(-1.+gb_oml+gb_om);
//	y2 = sqrt(gb_oml+(gb_om-a*(-1.+gb_oml+gb_om))/y);
//	y2 = y2*2.*a*y;
//	return y1/y2;

	y = pow(a, 3.*gb_w);
	y1 = -3.*y * gb_om + 2.*a * y * (gb_oml+gb_om-1.) - 3.*gb_oml*(1.+gb_w);
	y1 = y1 * gb_H0 / pow(a,4.) / y;
	y2 = gb_oml / y + gb_om - a*(gb_oml+gb_om-1.);
	y2 = y2 / pow(a,3.);
	y2 = 2. * pow(y2, 0.5);
	return y1/y2;
}
void gb_Tsquare_DE(double a, double D, double E, double QdDda, double QdEda, float *q1, float *q2){
	double H = Ha(a), Hp=dHada(a), Q=a*a*a*Ha(a)/gb_H0, y1, Dpp, Epp;
	double Dprime=QdDda/Q, Eprime=QdEda/Q;
	double dQda = 3./a*Q + Q*Hp/H;

	D*= D_final; E*=E_final; Dprime*=D_final; Eprime*=E_final;

//	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * D - 2*H*a*H*Dprime;
//	y1 = y1 / (a*H) - H*Dprime - a*Hp*Dprime;
//	y1 = y1 / (a*H);
	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * D - 2*H*a*H*Dprime;
	y1 = y1 / (a*H) - H*Dprime - a*Hp*Dprime;
	y1 = y1 / (a*H); Dpp = y1; ///D_final/gb_H0;


	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * (E-D*D);
        y1 = y1 - 2*H*a*H*Eprime;
	y1 = y1 / (a*H) - H*Eprime - a*Hp*Eprime;
	y1 = y1 / (a*H); Epp=y1; ///E_final/gb_H0;


	Dprime/= D_final; Dpp/= D_final; Eprime/=E_final; Epp/=E_final;
	*q1 = Q*dQda*Dprime + Q*Q*Dpp;
	*q2 = Q*dQda*Eprime + Q*Q*Epp;
}



int gb_iz(double z){
	int iz;
	iz = ceil(log(1.0+z)/gb_logbasenumber + 1.0);
	if(iz<2) iz=2;
	return iz;
}

double gb_intpl_vl(double x, double x1, double f1, double x2, double f2, double x3, double f3){
		double  d1, d2, d3, d12, d13, d23;
		d1  = x - x1;
		d2  = x - x2;
		d3  = x - x3;
		d12 = x1 - x2;
		d13 = x1 - x3;
		d23 = x2 - x3;
		return f1*d2*d3/(d12*d13) - f2*d1*d3/(d12*d23) + f3*d1*d2/(d13*d23);
}

void gb_solve_D_E_Qs();
//void intpl_D_E_Qs(double z, double *Dp, double *Ep, double *Q1p, double *Q2p);

/*int gb_iz(double z){
//	int iz;
//	iz = ceil(log(1.0+z)/gb_logbasenumber + 1.0);
//	if(iz<2) iz=2;
//	return iz;
//}

double gb_intpl_vl(double x, double x1, double f1, double x2, double f2, double x3, double f3){
		double  d1, d2, d3, d12, d13, d23;
		d1  = x - x1;
		d2  = x - x2;
		d3  = x - x3;
		d12 = x1 - x2;
		d13 = x1 - x3;
		d23 = x2 - x3;
		return f1*d2*d3/(d12*d13) - f2*d1*d3/(d12*d23) + f3*d1*d2/(d13*d23);
}


void gb_solve_D_E_Qs();*/
void intpl_D_E_Qs(double z, double *Dp, double *Ep, double *Q1p, double *Q2p){
	int i, i1, i2, i3;
	double z1, z2, z3, f1, f2, f3;
	double *p;

	if(!solvegrowth_flag) gb_solve_D_E_Qs();

	i = gb_iz(z);

	if(i < 2) i=2;
	if(i > GB_NUM_INTPL - 1) i = GB_NUM_INTPL -1;
	
	i-=1;

	i1 = i-1; i2 = i; i3 = i+1;
	//printf("i1,i2,i3= %d %d %d \n",i1,i2,i3);

	z1 = gb_zdata[i1];
	z2 = gb_zdata[i2];
	z3 = gb_zdata[i3];
	//printf("z,z1,z2,z3 = %f %f %f %f \n",z,z1,z2,z3);

	p = gb_Ddata;
	f1 = p[i1]; f2 = p[i2]; f3 = p[i3];
	//printf("f1,f2,f3= %f %f %f \n",f1,f2,f3);
	*Dp = gb_intpl_vl(z, z1, f1, z2, f2, z3, f3);

	p = gb_Edata;
	f1 = p[i1]; f2 = p[i2]; f3 = p[i3];
	//printf("f1,f2,f3= %f %f %f \n",f1,f2,f3);
	*Ep = gb_intpl_vl(z, z1, f1, z2, f2, z3, f3);

	p = gb_QDdata;
	f1 = p[i1]; f2 = p[i2]; f3 = p[i3];
	//printf("f1,f2,f3= %f %f %f \n",f1,f2,f3);
	*Q1p = gb_intpl_vl(z, z1, f1, z2, f2, z3, f3);

	p = gb_QEdata;
	f1 = p[i1]; f2 = p[i2]; f3 = p[i3];
	//printf("f1,f2,f3= %f %f %f \n",f1,f2,f3);
	*Q2p = gb_intpl_vl(z, z1, f1, z2, f2, z3, f3);
	//printf(" D, E, Q1, Q2 = %f %f %f %f \n",*Dp,*Ep,*Q1p,*Q2p);
}







// Now define equations
// D' = Dpda
// DpDa' = 
//
double dDda(double a, double ys[]){
//	double D = ys[0], Dprime = ys[1];
	double Dprime = ys[1];
	return Dprime; }
double dDprimeda(double a, double ys[]){
	double D = ys[0], Dprime = ys[1];
	double y1, H, Hp;
	H = Ha(a); Hp = dHada(a);
	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * D - 2*H*a*H*Dprime;
	y1 = y1 / (a*H) - H*Dprime - a*Hp*Dprime;
	y1 = y1 / (a*H);
	return y1; }
double dEda(double a, double ys[]){
//	double D = ys[0], Dprime = ys[1], E = ys[2], Eprime = ys[3];
	double Eprime = ys[3];
	return Eprime; }
double dEprimeda(double a, double ys[]){
//	double D = ys[0], Dprime = ys[1], E = ys[2], Eprime = ys[3];
	double D = ys[0], E = ys[2], Eprime = ys[3];
	double y1, H, Hp;
	H = Ha(a); Hp = dHada(a);
	y1 = (3.*gb_H0*gb_H0*gb_om)/(2.*a*a*a) * (E-D*D);
        y1 = y1 - 2*H*a*H*Eprime;
	y1 = y1 / (a*H) - H*Eprime - a*Hp*Eprime;
	y1 = y1 / (a*H);
	return y1; }


void runge_kutta( double (*fs[])(), int numf, double x0, double f0[], double x1, double f1[], int nstep  )
{
	int ifun;
//	double h, k1,k2,k3,k4, x, *ys, *k1s, *k2s, *k3s, *k4s, *ks;
	double h, x, *ys, *k1s, *k2s, *k3s, *k4s, *ks;
	for(ifun=0;ifun<numf;ifun++){f1[ifun]=f0[ifun];}
	h = (x1-x0) / (double) nstep;

	int istep, iy; 
	x = x0;
	ys = malloc(sizeof(double)*numf);
	k1s = malloc(sizeof(double)*numf);
	k2s = malloc(sizeof(double)*numf);
	k3s = malloc(sizeof(double)*numf);
	k4s = malloc(sizeof(double)*numf);
	ks = malloc(sizeof(double)*numf);
	for(istep=1;istep<=nstep;istep++){
		for(iy=0;iy<numf;iy++) ys[iy]=f1[iy];
		for(iy=0;iy<numf;iy++) k1s[iy]=h*(*fs[iy])(x,ys);
		for(iy=0;iy<numf;iy++) ys[iy]=f1[iy]+k1s[iy]*0.5;
		for(iy=0;iy<numf;iy++) k2s[iy]=h*(*fs[iy])(x+h*0.5,ys);
		for(iy=0;iy<numf;iy++) ys[iy]=f1[iy]+k2s[iy]*0.5;
		for(iy=0;iy<numf;iy++) k3s[iy]=h*(*fs[iy])(x+h*0.5,ys);
		for(iy=0;iy<numf;iy++) ys[iy]=f1[iy]+k3s[iy];
		for(iy=0;iy<numf;iy++) k4s[iy]=h*(*fs[iy])(x+h,ys);
		for(iy=0;iy<numf;iy++) ks[iy]=(k1s[iy]+2.*k2s[iy]+2.*k3s[iy]+k4s[iy])/6.;
		for(iy=0;iy<numf;iy++) f1[iy]+=ks[iy];
		x += h;
	}
	free(ys);
	//printf("Stop at x= %lf\n",x);
}

double dy1dx(double x, double ys[]){return 6.*x*x*x*x*x;}
//double dy1dx(double x, double ys[]){return 6.*ys[1]*x;}
//double dy2dx(double x, double ys[]){return 4.*x*x*x;}
double dy2dx(double x, double ys[]){return 4.*pow(ys[0],0.5);}

void check_funs(){
	double (*fs[2])() = {dy1dx,dy2dx};
	int numf = 2;
	double x0 = 0.0, x1=3.0;
	double f0[2] = {0.0,0.0};
	double f1[2];
	int nstep = 128;

	double a = 0.5;
	printf("(check_funs check our result)   a = %lf; Ha(a) =%lf dHada(a) =%lf\n", a, Ha(a), dHada(a));
	printf("(check_funs Mathematica result) a = 0.5; Ha(a) =124.631  dHada = -292.51\n");

	runge_kutta(fs, numf, x0, f0, x1, f1, nstep);
	int ifun;
	printf("(check_funs) value of x^6, x^4 at x = %lf :\n\t",x1);
	for(ifun=0;ifun<numf;ifun++)
		printf("%lf ",f1[ifun]);
	putchar('\n');
}

void get_D_E(double nowa, double *Da, double *Ea, double *QdDda, double *QdEda){

	//gb_om = gb_om; gb_w = gb_w; gb_oml = gb_oml;
	//printf(" (get_D_E) gb_om, gb_oml, gb_w = %f %f %f \n",gb_om,gb_oml,gb_w);

	// functions of dy/da;
	double (*fs[4])() = {dDda, dDprimeda, dEda, dEprimeda};
	int numf = 4;
	// initial and final value of scale factor
	double ainit = 0.00001, afinal = 1.0;

	// initial values of D, Dprime, E, Eprime
	double D_init = 0.00001, Dprime_init = 0.0, E_init = -3./7.*D_init*D_init, Eprime_init = 0.0;
	//double f_init[4] = {D_init, Dprime_init, E_init, Eprime_init};
	double f0[4] = {D_init, Dprime_init, E_init, Eprime_init};

	// final values of D, Dprime, E, Eprime
	double f1[4]; 
	double f_final[4];

	// time steps of RG4
	int nstep = 4096*4;

	//printf("  (get_D_E) called at a = %.4f \n",nowa);

	// solve the growth factors
	// 1. get value of growth factor at a = 1.0
	runge_kutta(fs, numf, ainit, f0, afinal, f_final, nstep*10);
	// 2. value of growth factor at current scale factor
	runge_kutta(fs, numf, ainit, f0, nowa, f1, (int)(nstep*10*nowa));
	*Da = f1[0]/f_final[0]; *Ea = f1[2]/f_final[2];

	// Q related varaibles
	double Q = nowa*nowa*nowa*Ha(nowa);
	*QdDda = Q*f1[1]/f_final[0]/gb_H0;
	*QdEda = Q*f1[3]/f_final[2]/gb_H0;
}

void test(){
	check_funs();

	// functions of dy/da;
	double (*fs[4])() = {dDda, dDprimeda, dEda, dEprimeda};
	int numf = 4;

	// initial and final value of scale factor
	double ainit = 0.00001;//, afinal = 1.0;

	// initial values of D, Dprime, E, Eprime
	double D_init = 0.00001, Dprime_init = 0.0, E_init = -3./7.*D_init*D_init, Eprime_init = 0.0;
	//double f_init[4] = {D_init, Dprime_init, E_init, Eprime_init};
	double f0[4] = {D_init, Dprime_init, E_init, Eprime_init};

	// final values of D, Dprime, E, Eprime
	double f1[4]; 
	double f_final[4];

	// time steps of RG4
	int nstep = 4096*4;

	// solve the growth factors
	double a; int ifun;
	// 1. get value of growth factor at a = 1.0
	runge_kutta(fs, numf, ainit, f0, 1.0, f_final, nstep*10);
	// 2. value of growth factor at many as
	for(a=0.1; a<=1.0; a+= 0.1){
		runge_kutta(fs, numf, ainit, f0, a, f1, nstep);
		ainit = a;
		printf("At a = %lf, D = %lf, E = %lf \n", a, f1[0]/f_final[0], f1[2]/f_final[2]);
		for(ifun=0;ifun<numf;ifun++) f0[ifun]=f1[ifun];
	}
	
	
}

void gb_solve_D_E_Qs(){

	//gb_om = gb_om; gb_w = gb_w; gb_oml = gb_oml;
	msg_printf(normal, " (gb_solve_D_E_Qs) gb_om, gb_oml, gb_w = %f %f %f \n",gb_om,gb_oml,gb_w);

	// functions of dy/da;
	double (*fs[4])() = {dDda, dDprimeda, dEda, dEprimeda};
	int numf = 4;
	// initial and final value of scale factor
	double ainit = 1. / (1.+gb_zi(GB_NUM_INTPL));
	double afinal = 1.0;

	// initial values of D, Dprime, E, Eprime
	double D_init = ainit, Dprime_init = 0.0, E_init = -3./7.*D_init*D_init, Eprime_init = 0.0;
	//double f_init[4] = {D_init, Dprime_init, E_init, Eprime_init};
	double f0[4] = {D_init, Dprime_init, E_init, Eprime_init};
	//double fa[4], fb[4];

	// final values of D, Dprime, E, Eprime
	double f1[4],f2[4]; 
	double f_final[4];

	// time steps of RG4
	int nstep = 4096*4;

	//printf("  (gb_solve_D_E_Qs) called!\n ");


	runge_kutta(fs, numf, ainit, f0, afinal, f_final, nstep*10);
	// 2. value of growth factor at current scale factor
	f1[0] = f0[0]; f1[1] = f0[1]; f1[2] = f0[2]; f1[3] = f0[3];

	double a1,a2,Q; //= nowa*nowa*nowa*Ha(nowa);
	int i;
	for(i=GB_NUM_INTPL; i>=1; i--){
		a1 =1. /(1.+ gb_zi(i));
		a2 =1. /(1.+ gb_zi(i-1));
		Q = a2*a2*a2*Ha(a2);
		runge_kutta(fs, numf, a1, f1, a2, f2, 4);
		gb_Ddata[i-1] = f2[0]/f_final[0]; 
		gb_Edata[i-1] = f2[2]/f_final[2];
		gb_QDdata[i-1] = Q*f2[1]/f_final[0]/gb_H0;
		gb_QEdata[i-1] = Q*f2[3]/f_final[2]/gb_H0;
		gb_zdata[i-1] = gb_zi(i-1);
		//if(i%100==0){printf("   i,a,D,E,QD,QE = %i %f %f %f %f %f \n",i,a2,gb_Ddata[i-1],gb_Edata[i-1],gb_QDdata[i-1],gb_QEdata[i-1]);}
		f1[0]=f2[0]; f1[1]=f2[1]; f1[2]=f2[2]; f1[3]=f2[3];
	}

	D_final=f_final[0]; Dp_final=f_final[1]; E_final=f_final[2]; Ep_final=f_final[3];
	solvegrowth_flag = 1;
	
	//runge_kutta(fs, numf, ainit, f0, nowa, f1, (int)(nstep*10*nowa));
	//*Da = f1[0]/f_final[0]; *Ea = f1[2]/f_final[2];

	// Q related varaibles
	//*QdDda = Q*f1[1]/f_final[0]/gb_H0;
	//*QdEda = Q*f1[3]/f_final[2]/gb_H0;
}
