/****************************************************************************
目的：    考虑动力学补偿算法, 用动力学轨道积分预报轨道作为先验轨道, 用伪距观测
          数据, 推广卡尔曼滤波更新卫星的运行轨道参数.

编写时间：2008.12.19
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _EXTEND_KALMAN_FILTER_
#define _EXTEND_KALMAN_FILTER_

#include "RTOD_Structs.h"
#include "DataProcess.h"
#include <stdio.h>
#define	 TAO_VTEC	10
//根据伪距计算出的倾斜电离层延迟最小二乘得到垂直电离层延迟
void initialVTEC(const EpochObsData* Epoch,
	const SATMIDRESULT* SatMidInfo,const PPRESULT Result,double * const VTEC,double *const sigma);
// 目的: 初始化AB的卫星卡尔曼滤波的状态
int init2Sat(FILE * Fobs_A,OBSTYPELIST * ObsTypeList_A,EpochObsData * Epoch_A,GPSEPHREC* GPSEph_A,SATMIDRESULT * SatMidInfo_A,
	FILE * Fobs_B,OBSTYPELIST * ObsTypeList_B,EpochObsData * Epoch_B,GPSEPHREC* GPSEph_B,SATMIDRESULT * SatMidInfo_B,
	double OiState[108],EKFSTATE * EkfState);

/****************************************************************************
initASat

目的: 初始化参考星A的状态
参数: FILE * Fobs_A
参数: OBSTYPELIST * ObsTypeList_A
参数: EpochObsData * Epoch_A
参数: GPSEPHREC * GPSEph_A
参数: SATMIDRESULT * SatMidInfo_A
参数: double OiState[54]
参数: EKFSTATE * EkfState
****************************************************************************/
int initASat(FILE * Fobs_A,OBSTYPELIST * ObsTypeList_A,EpochObsData * Epoch_A,GPSEPHREC* GPSEph_A,SATMIDRESULT * SatMidInfo_A,
	double OiState[54],EKFSTATE * EkfState);

/****************************************************************************
initRel

目的: 初始化相对状态
参数: FILE * Fobs_B
参数: OBSTYPELIST * ObsTypeList_B
参数: EpochObsData * Epoch_B
参数: GPSEPHREC * GPSEph_B
参数: SATMIDRESULT * SatMidInfo_B
参数: double OiState[54]
参数: EKFSTATE * EkfState
****************************************************************************/
int initRel(FILE * Fobs_A, OBSTYPELIST * ObsTypeList_A, EpochObsData * Epoch_A, GPSEPHREC* GPSEph_A, SATMIDRESULT * SatMidInfo_A, FILE * Fobs_B,OBSTYPELIST * ObsTypeList_B,EpochObsData * Epoch_B,GPSEPHREC* GPSEph_B,SATMIDRESULT * SatMidInfo_B, double OiState[54],EKFSTATE * EkfState);
// 目的: 对准观测值与卡尔曼滤波的时间 //如果确定一定是0 30s 60s这样，则不需要该函数
void alignObsEKF(int graceType,FILE * Fobs_A,OBSTYPELIST * ObsTypeList_A,EpochObsData * Epoch_A,GPSEPHREC* GPSEph_A,SATMIDRESULT * SatMidInfo_A,
	FILE * Fobs_B,OBSTYPELIST * ObsTypeList_B,EpochObsData * Epoch_B,GPSEPHREC* GPSEph_B,SATMIDRESULT * SatMidInfo_B,
	EKFSTATE * EkfState);

// 目的: 对准两卫星观测值的时间 
void alignEpoch( FILE * Fobs_A,OBSTYPELIST * ObsTypeList_A,EpochObsData * Epoch_A,FILE * Fobs_B,OBSTYPELIST * ObsTypeList_B,EpochObsData * Epoch_B);

//因长时间中断的更新
int intUpdate(FILE * Fout_A,EpochObsData * Epoch_A,
	FILE * Fout_B,
	double   OiState[108],EKFSTATE * EkfState);

// 目的: 清空EKFSTATE结构体
void EmptyEKFSTATEStruct( EKFSTATE* EkfState );
//初始化双差						
void EmptyddObsStruct(DDPSEUDOOBS *ddObs);
//初始化共同卫星
void EmptyCurComObsStruct(Common11Obs* ComObs);
//初始化之前的共同卫星结构体
void EmptyComObssStruct(Common11Obs * ComObss);
//初始化单颗共同卫星结构体
void EmptyOneSat11Obs(OneSat11Obs * oneSat11Obs);

void InitAStateAndCov(const GPSTIME* Time, double Cd, double Cr, double Tau, 
                        PPRESULT* Result, EKFSTATE* EkfState );



void InitRelStateAndCov(const GPSTIME* Time, PPRESULT Result, EKFSTATE* EkfState );
/***************************************************************************
//
// FormStateTransMatrixFromDyn
//
// 目的: 根据动力学积分生成的状态转移矩阵(6*6), 补充接收机钟差/动力学补偿参数, 
         形成最终的推广卡尔曼滤波的状态转移矩阵(EKFDIMENSION*EKFDIMENSION].
//
// 输入参数:

    Tao      补偿加速度的一阶高斯马尔可夫模型的相关时间
    Step     推广卡尔曼滤波的时间更新间隔[s]
    w[3]     补偿加速度[3维]  m/s^2
    STMDyn   动力学积分生成的状态转移矩阵, 只包含位置和速度[6*8]
    
   输出参数
   
    STM     推广卡尔曼的状态转移矩阵[EKFDIMENSION*EKFDIMENSION]

***************************************************************************/

void FormStateTransMatrixFromDynRTN4A( const double Tao, const double Step, double State[ADIMENSION], double STMDyn[48], double STM[] );

void FormStateTransMatrixFromDynRTN4rel( const double Tao, const double Step, double State[RELDIMENSION], double STMDyn[48], double STM[] );
void DMCrtn4A( const double Tao, const double Step, double State[]);
//补偿加速度是在rtn方向，上面是在xyz方向
void UpdateDeterministicComponentFromDMCrtn( const double Tao, const double Step,
	double State[]);

/***************************************************************************
//
// UpdateGNSSReceiverClockOffset
//
// 目的: 根据一阶高斯马尔可夫模型, 时间更新GNSS接收机的GPS钟差和GLONASS钟差.
         看技术文档, 仔细优化.

// 输入参数:

    Tao      补偿加速度的一阶高斯马尔可夫模型的相关时间
    Step     推广卡尔曼滤波的时间更新间隔[s]
    w[3]     补偿加速度[3维]  m/s^2
    State    动力学积分预报的卫星轨道, 只包含位置和速度[6]
    

***************************************************************************/

void UpdateGNSSReceiverClockOffset( double Tao, double Step,
                                 double w[3], double State[6] );
void DMCrtn4Rel( const double Tao, const double Step, double State[]);
void UdProcessNoiseCovRTN4A(double Tao, double Step, double Sigma,
	double State[6], double Qcov[] );
                   
void UdProcessNoiseCovRTN4rel( double Tao, double Step, double Sigma, double State[6], double Qcov[] );
int intUpdate(FILE * Fout_A,EpochObsData * Epoch_A,
	FILE * Fout_B,
	double   OiState[108],EKFSTATE * EkfState);
void EKFTimeUpdateA(double OIState[54], EKFSTATE * KFState );
bool EKFMeasureUpdateA( double O_C, double R, double H[], EKFSTATE* KFState );
/***************************************************************************
//
// EKFMeasureUpdate
//
// 目的: 使用UD滤波的原理, 将每个历元的多颗卫星的观测值组成序列, 依次对接收机的
         状态进行测量更新.

// 输入参数:

    O_C      观测值减计算值[m]
    R        观测值的精度指标[m]
    H        观测值线性化的系数向量
    
    输出参数
        
    KFState  输入为经过时间更新后的卡尔曼滤波状态与协方差矩阵信息, 经过一颗GNSS
             卫星的观测数据更新后, 输出更新后的状态和协方差矩阵信息

***************************************************************************/
bool EKFMeasureUpdateRel( double O_C, double R, double H[], EKFSTATE* KFState );
//总体更新
int WholeMeasureUpdate(FILE * Fout_A,EpochObsData * Epoch_A,GPSEPHREC* GPSEph_A,SATMIDRESULT * SatMidInfo_A,
	FILE * Fout_B,EpochObsData * Epoch_B,GPSEPHREC* GPSEph_B,SATMIDRESULT * SatMidInfo_B,
	double  OiState[108],EKFSTATE* KFState);
void ABTimeUpdate(double OiState[108],EKFSTATE* KFState);
#endif




















