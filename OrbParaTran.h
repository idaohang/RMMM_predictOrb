#pragma once 
#ifndef ORBIT_PARA_TRANSLATION_H_
#define ORBIT_PARA_TRANSLATION_H_


//wkliu  该变量是暂时存放在这里
const double utcf = 0.;

void OsculToMean(int ymdhm[5], double sec, double osc[6], double mean[6], 
				 double x1_osc[6], double x1_mean[6]);


int OrbParaCon(double utc, double* xyz, double* xin, double* xout, double* xit, 
				double* xot, double& da, double& c1, double& c2, double& c3);

void FuncMixIne(int flag, int Day, double SecOfDay, double* osc, double* ggs);

void FuncTrKfe(int flag, double* Norm_A, double* Norm_B);

void RadianNorm(double& ang);

void FuncMast1(int flag, int Day, double SecOfDay, double* b1, double* sb);

void BMoon(int Day, double SecOfDay, double* bm);

void BSun(int Day, double SecOfDay, double* bs);

void FuncGkepl(double* s, double* sb);

void FuncGecu(double* s, double* sb);

void FuncSHO1(double* s, double* sb, double* ss);

void FuncSHO2(double* s,double* sb,double as1,double* ss);

void FuncQfNmp(int n,double six[17],double cix[17],double nmp[5][5][5],double nmpb[5][5][5]);

void FuncGtSco(double dt, double s0, int n, double* s);

void FuncGbese(int n, double z, double* bi);

void FuncGbesq(double z, double* bi);

void FuncGtSho(double dt, double s0, int n, double nmp[5][5][5], 
			   double nmpb[5][5][5], double* s, double* sb, double* ss);

void FuncGdrto(double* GS,double* GSB);

void FuncGdrsh(double dt, double bl0, double* s, double* sb, double* ss);

void FuncKepl(double e, double m, double& x, double& y, int& status);

int anint(double val);

inline double dmax1(double x,double y);

void CalDaySec(double utc, int* ITime, double& stime);

void FuncLdlt(int n, int m, double* Arra, double* Arrp, 
		  double* Arrb, double* Arrx, int& isw);

void FuncAbMul(double* Arra, double* Arrb, double* Arrr, int l, int m, int n);

void NomialWsh(int num, double* tim, double* xtp, double& da);

int VTVMul(int n, int m, double* a, double* el, double* p, double* x, double* px, 
			double* atp, double* an, double* atpl, double* at, int& isw);

#endif