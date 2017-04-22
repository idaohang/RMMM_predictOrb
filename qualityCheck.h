#ifndef _QC_
#define _QC_

#include "DataProcess.h"
#include "RTOD_Structs.h"
double LSRel_PC(const int satNum, int * PRN,int * banPRN,const double * satPos,const double * StateInECEFA,
	const double * obs ,double * AprioriPosRel,double * AprioriClk, PPRESULT * Result);
double LEGE(GPSTIME Time, const int satNum, int * PRN,int * banPRN,const double * satPos,const double * StateInECEFA,const double * obs ,double * AprioriPosRel,double * AprioriClk, PPRESULT * Result);
int valResult(const double sigmaPos,const double sigma0,int nv,int nx);
//Ö»¿¼ÂÇÒ»¸ö´Ö²î
double Raim_fde1(int * satNum, int * PRN,int * banPRN,const double * satPos,const double * StateInECEFA,const double * obs ,double * AprioriPosRel,double * AprioriClk, PPRESULT * Result);
#endif
