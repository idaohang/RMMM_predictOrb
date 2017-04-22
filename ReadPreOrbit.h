#pragma once
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "stdlib.h"
#include <cstring>
#include "gpstime.h"
#include "RTOD_Structs.h"

struct ONESATEHP //一个历元一颗卫星的精密星历（卫星号、位置、速度、钟差） 
{
	int PRN;     //卫星的PRN号

	bool   ISExist; //判断该卫星的精密星历是否存在

	double  satposx;            //卫星位置-x
	double  satposy;            //卫星位置-y
	double  satposz;            //卫星位置-z

	double  satvelx;            //卫星速度-vx
	double  satvely;            //卫星速度-vy
	double  satvelz;            //卫星速度-vz
	
	double  satclk;   //卫星钟差
	double  satvelc;  //卫星钟速 
};

struct ONEPOCHOFSP3 //一个历元所有卫星的精密星历
{
	int     WeekNumber;         //GPS周的周号

	double  dWeeksec;           //GPS周的周秒

	double  totalWeeksec;       //总的GPS周秒，用于拉格朗日内插

	COMMONTIME   pT;                 //通用时

	int     Satnum;             //单个历元的卫星的数目

	struct  ONESATEHP obs[32];  //保存非差观测数据的数组

};
struct ONESATCLK
{
	int PRN;     //卫星的PRN号

	bool   ISExist; //判断该卫星的精密星历是否存在
	double  satclk;   //卫星钟差
	double  satvelc;  //卫星钟速 
};
struct ONEPOCHPRECLK
{
	int     WeekNumber;         //GPS周的周号

	double  dWeeksec;           //GPS周的周秒

	double  totalWeeksec;       //总的GPS周秒，用于拉格朗日内插

	COMMONTIME   pT;                 //通用时
	
	struct ONESATCLK obs[32];//保存一颗卫星的精密钟差


};
void ReadSp3File(char* filename,struct ONEPOCHOFSP3 *OBSP);
void ReadPreClkFile(char* filename,struct ONEPOCHPRECLK *CLKP);
int FindFitPoint_SP3(const int PRN,GPSTIME OneEpoch,int n,struct ONEPOCHOFSP3 *OBSP, double t[], 
	double  satPos[3][30], double clk[]);//2n为拟合阶数
int FindFitPoint_CLK(GPSTIME OneEpoch,int n,struct ONEPOCHPRECLK *CLKP,struct ONEPOCHPRECLK *FITPOINTCLK);
void Lagrange_SP3(GPSTIME OneEpoch,int PRN,int n,struct ONEPOCHOFSP3 *OBSP,ONESATEHP *PREORBIT_SP3);
int Lagrange_CLK(GPSTIME OneEpoch,int PRN,int n,struct ONEPOCHPRECLK *CLKP,ONESATCLK *PREORBIT_CLK);
void Rotation(int n,int PRN,double dt,ONEPOCHOFSP3 *FITPOINTSP3);
