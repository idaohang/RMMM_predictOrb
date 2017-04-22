/****************************************************************************
目的：    对观测数据进行处理，生成共同观测值，探测周跳
编写时间：2016.07.27
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _DATAPROCESS_
#define _DATAPROCESS_

#include "RTOD.h"
#include "RTOD_Structs.h"
#include "EKFilter.h"
#include "CommonFuncs.h"
#include "Trace.h"
#include "qualityCheck.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <afx.h>	
#define  Me_L		0.005								//载波相位
//#define  Me_P1		0.20
//#define  Me_P2		0.25
#define	 Me_dPC		0.57								/*无电离层伪距观测噪声*/	
#define  Me_dLC		4.0*Me_L							/*单差无电离层载波观测噪声*/
#define  C_Light	299792458.0							/* Speed of light  [m/s]; IAU 1976  */
#define  SIGMW		10000								/*初始宽巷模糊度方差,1.37是通过初始化公式MW推导而得*/
#define  SIGL1		10000								/*初始L1模糊度方差*/
#define  GF_THRES	0.05								/*LC探测周跳限值,单位m*/
#define	 LEMTA_NL	0.106953378142147					/*窄巷波长，单位为米*/
#define	 LEMTA_WL	0.861918400322005					/*宽巷波长，单位为米*/
#define	 RATIO		0.779220779220779					/*ratio=(f2/f1)*/
#define	 RATIO2		0.607185022769438					/*ratio2=(f2/f1)*(f2/f1)*/
#define	 RATIOK		0.124087591240876					/*ratiok=(f1-f2)/(f1+f2)*/
/*无电离层组合宽巷模糊度前观测系数Hwl*/
#define  Hwl		0.37748251108920522					/*LEMTA_1*RATIO/(1.0-RATIO2)*/
#define  NC_SIG     0.5									/*无电离层模糊度中误差设为0.5m*/
#define  MW_SIG		0.5									/*中误差为0.5周*/
#define	 INTERVAL	600.0								/*数据空缺时段大于10分钟则重新全部初始化载波*/
//生成A/B共同的观测卫星数据
int GenCommomObs(EpochObsData * Epoch_A,EpochObsData * Epoch_B,Common11Obs * tempComObs);
/*初始化模糊度的协方差，通过协方差传播定律*/
void initialAmbsCov(const double B[48],const double Dl[12],EKFSTATE * EkfState,const double sigma);
void getBpos(double * StateA,double * relState,double * StateB);
void getRelPos(double * StateA,double * StateB,double * relState);
//单差周跳探测，若发生周跳则标识为-1，不发生则标识为0，新观测值标识为1
int CySlipDetection_SD(double * Mat, Common11Obs * preComObs,Common11Obs * curComObs,EKFSTATE * EkfState);
//非差周跳探测
int CySlipDetection_ND(double * Mat, Common11Obs * preComObs,Common11Obs * curComObs,EKFSTATE * EkfState);
//周跳探测方法
/*电离层探测周跳法-单差观测值*/
int GFDetection_SD(OneSat11Obs * PreObs,OneSat11Obs * CurObs, double * dGF);
/*MW探测周跳法-单差观测值*/
int MWDetection_SD(OneSat11Obs * PreObs,OneSat11Obs * CurObs, double * tempMW);
/*电离层探测周跳法-非差观测值*/
int GFDetection_ND(OneSat11Obs * PreObs,OneSat11Obs * CurObs, double * dGF);
/*MW探测周跳法-非差观测值*/
int MWDetection_ND(OneSat11Obs * PreObs,OneSat11Obs * CurObs,double * tempMW);
//通过后处理的周跳探测结果产生的LLI来进行周跳判断
int CySlipDetection_LLI(double *Mat,Common11Obs * preComObs,Common11Obs * curComObs,EKFSTATE * EkfState);
int PPRTOD_LC(Common11Obs * CurObs,EKFSTATE * KFState, PPRESULT * Result, double * Mat);	/*最小二乘探测周跳法*/
//伪距几何法相对定轨
int PPRTOD_PC(Common11Obs * CurObs,EKFSTATE * KFState, PPRESULT * Result, double * Mat);
void EKFTimeUpdateRel( double OIState[108], EKFSTATE* KFState );
void outRes(int satNumWhole,FILE *fLC,double *OC_LC0,double *OC_LC1,EKFSTATE * KFState);
//相对测量更新
int relEKFilter_PC(double Mat[9],EKFSTATE * KFState);
int relEKFilter_LC(double Mat[9],EKFSTATE * KFState);
//
double getMap(double radEle);									/*得到投影函数值*/
/*输出伪距滤波验前后残差均值*/
void outResEKFPC(int satNumWhole,double *OC_PC0,double *OC_PC1,EKFSTATE * KFState);
/*输出载波滤波验前后残差均值*/
void outResEKFLC(int satNumWhole,double *OC_LC0,double *OC_LC1,EKFSTATE * KFState);
/*输出LC-PC的残差均值，检查是否为0，以判断周跳和模糊度处理正确*/
void outResAmbSD(const double sow,const double LC_PC0,const double LC_PC1);
#endif