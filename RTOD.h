/****************************************************************************
目的：    单点实时定轨, 结合动力学模型滤波实时定轨

编写时间：2008.12.5
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _REALTIME_ORBIT_DETERMINATION_
#define _REALTIME_ORBIT_DETERMINATION_

#include "ReadObs.h"
#include "EphProc.h"
#include "EKFilter.h"
#include "ReadPreOrbit.h"
#include "PhaseCent_Cor.h"
extern ONEPOCHOFSP3     *PRE_EPH;
extern pcvs_t pcvs;
/*
  计算信号发射时刻卫星的位置和速度, GLONASS卫星位置转换为WGS84.
  计算信号发射时刻卫星的钟差,

  计算相对论效应改正
  计算电离层延迟改正

  单点定位子程序
  单点测速子程序  用于初始化滤波定轨

  单点定位残差检查
  单点定位数据预处理

  滤波定轨子程序


 */

/***************************************************************************
//
// EmptyAPRIORBITStruct
//
// 目的: 初始化APRIORBIT结构体

   输入参数

    ApriOrbit      待初始化的变量

***************************************************************************/

void EmptyAPRIORBITStruct( APRIORBIT* ApriOrbit );


/***************************************************************************
//
// EmptySATLISTStruct
//
// 目的: 初始化SATLIST结构体

   输入参数

    SatList      待初始化的变量

***************************************************************************/

void EmptySATLISTStruct( SATLIST* SatList );


/***************************************************************************
//
// EmptyPPRESULTStruct
//
// 目的: 初始化PPRESULT结构体

   输入参数

    PPResult      待初始化的变量

***************************************************************************/

void EmptyPPRESULTStruct( PPRESULT* PPResult );
                                  

/***************************************************************************
//
// EmptySATMIDRESULTStruct
//
// 目的: 初始化SATMIDRESULT结构体

   输入参数

    Num             SATMIDRESULT数组的维数
    SatMidInfo      待初始化的变量

***************************************************************************/

void EmptySATMIDRESULTStruct( int Num, SATMIDRESULT* SatMidInfo );


/***************************************************************************
//
// ComputeGPSSatOrbitAtSignalTrans
//
// 目的: 计算信号发射时刻的GPS卫星轨道,卫星钟差,卫星高度角
         考虑地球自转改正(顾及卫星健康标记).

  URA指标( 导航电文中的SV accuracy以米为单位 )
  If the value of N is 6 or less, X=pow(2.0, 1+N/2 );
  If the value of N is 6 or more, but less than 15, X=pow(2,N-2);
  N=15 shall indicate the absence of an accuracy prediction and shall advise the
  unauthorized user to use that SV at his own risk.

  SV Healthy (GPS)
  The MSB shall indicate a summary of the health of the NAV data, where
  0 = all NAV data are OK
  1 = some or all NAV data are bad.
  The five LSBs shall indicate the health of the signal components in accordance
  with the codes given as following:
  
  Estimated Group Delay Differential (TGD)
  L1-L2 correction term, for the benifit of "L1 only" or "L2 only" users.
  instruct:  (deltaTsv)L1 = deltaTsv - Tgd
 
  // 输入参数:
  //
  Prn             GPS卫星号
  Time            观测时刻
  PreLeoOrb[3]    预报的星载接收机位置[m]
  PreLeoClk       预报的星载接收机GPS钟差[m]
  Height          星载接收机的高程[m]
  GPSEph          某颗GPS卫星的星历
  IonPara         电离层参数
  
	输出参数
    
	  SatMidInfo      卫星轨道, 速度等中间计算过程信息
	  
		返回值
		
		  计算成功返回true, 否则返回false.
		  
***************************************************************************/

bool ComputeGPSSatOrbitAtSignalTrans( const short Prn, const GPSTIME* Time,
									 double PreLeoOrb[3], double PreLeoClk, double &Height,
                                     GPSEPHREC* GPSEph, IONOPARA* IonPara,
                                     SATMIDRESULT* SatMidInfo );
                                  
/***************************************************************************
//
// ComputeGLONASSSatOrbitAtSignalTrans
//
// 目的: 计算卫星发射时刻的GLONASS卫星轨道,卫星钟差,相对论效应,卫星高度角
         时间以GPS时间系统为参考, 卫星轨道转换到WGS84坐标系.
         考虑地球自转改正(顾及卫星健康标记).

  // 输入参数:
  //
  // 输入参数:
  //
  Slot            GLONASS卫星号
  Time            观测时刻的GPS时间
  PreLeoOrb[3]    预报的星载接收机位置[m]
  PreLeoClk       预报的星载接收机GLONASS钟差[m]
  Height          星载接收机的高程[m]
  GLOEph          某颗GLONASS卫星的星历
  GloTmCorr       GLONASS时间系统改正参数
  IonPara         电离层参数
  
	输出参数
    
	  SatMidInfo      卫星轨道, 速度等中间计算过程信息
	  
		返回值
		
		  计算成功返回true, 否则返回false.
		  
***************************************************************************/

bool ComputeGLONASSSatOrbitAtSignalTrans( const short Slot, const GPSTIME* Time,
                                  double PreLeoOrb[3], double PreLeoClk, double& Height,
                                  GLONASSEPHREC* GLOEph, GLOTIMECORR* GloTmCorr,
                                  IONOPARA* IonPara, SATMIDRESULT* SatMidInfo );
                                  
/***************************************************************************
//
// ComputeQvv
//
// 目的: 单点定位后, 计算观测值改正数的协因数矩阵, 用于探测验后观测值粗差数据.

//
// 输入参数:
//
    m               观测值个数, 相对于使用的卫星数
	n               估计的参数个数, 对于单系统定位为4, 双系统组合定位为5
	A               观测系数矩阵 [m*n]
	InvN            Inv(ATPA)
	AT              A矩阵的转置
	Pll             观测值的权阵,为对角线元素

  输出参数
    
	Qvv             改正数的协因数矩阵, 为对角线元素

***************************************************************************/
void ComputeQvv( int m, int n, double A[], double InvN[], double AT[], 
				double Pll[], double Qvv[] );

/***************************************************************************
//
// Klobuchar
//
// 目的: 使用Klobuchar模型, 计算单频接收机的电离层延迟改正量

//
// 输入参数:
//
    Time      观测时刻,以GPS系统时间
    SatPos    卫星位置
    RcvPos    接收机位置
    IonoPara  电离层参数

  返回值
    
    信号路径方向电离层延迟改正值.

***************************************************************************/

double Klobuchar( const GPSTIME* Time, double SatPos[3],
                                     double RcvPos[3], IONOPARA* IonoPara);

/***************************************************************************
//
// KlobucharWithScale
//
// 目的: 使用Klobuchar模型, 计算单频接收机的电离层延迟改正量. 根据LEO和GPS
         卫星的位置, 计算比例因子. 该模型适用于低轨卫星的单频电离层改正

   模型参考文献:
   THE JOURNAL OF NAVIGATION (2002), 55, 293±304. 
   The Royal Institute of Navigation 
   DOI: 10.1017}S0373463302001789 Printed in the United Kingdom
   Ionospheric Correction for GPS Tracking of LEO Satellites
   Oliver Montenbruck and Eberhard Gill (German Aerospace Center DLR)

//
// 输入参数:
//
    Time      观测时刻,以GPS系统时间
    SatPos    卫星位置
    RcvPos    接收机位置
    IonoPara  电离层参数

  返回值
    
    信号路径方向电离层延迟改正值.

***************************************************************************/

double KlobucharForLEO( const GPSTIME* Time, double SatPos[3],
                                     double RcvPos[3], IONOPARA* IonoPara);


/***************************************************************************
//
// PointPositionRTOD
//
// 目的: 根据卫星星历和一个历元的观测数据C1, 进行单点定位(包括仅GPS, 
         仅GLONASS, GPS/GLONASS组合定位), 为滤波提供初始坐标信息.

//
// 输入参数:
//
    GPSEph          GPS卫星星历[32]
    IonPara         电离层参数
    GLOEph          GLONASS卫星星历[32]
    GloTmCorr       GLONASS时间与UTC(SU)系统间的改正值
    EpochObs        历元的观测数据, 以此顺序计算卫星轨道
    PreOrb          预报的星载接收机位置[m],预报的星载接收机钟差[m](GPS, GLONASS)

  输出参数
  
    SatMidInfo      卫星轨道, 速度等中间计算过程信息, 用于滤波定轨
    Result          单点定位结果
  
  返回值
    
    参与定位的卫星数

***************************************************************************/
int PointPositionRTOD( int graceType, GPSEPHREC* GPSEph, IONOPARA* IonPara,
                      EpochObsData* Epoch, APRIORBIT* PreOrb, 
                      SATMIDRESULT* SatMidInfo, PPRESULT* Result );

/***************************************************************************
//
// PointPositionRTODInECI
//
// 目的: 根据卫星星历和一个历元的观测数据C1, 进行单点定位(包括仅GPS, 
         仅GLONASS, GPS/GLONASS组合定位), 为滤波提供初始坐标信息.

//
// 输入参数:
//
    GPSEph          GPS卫星星历[32]
    IonPara         电离层参数
    GLOEph          GLONASS卫星星历[32]
    GloTmCorr       GLONASS时间与UTC(SU)系统间的改正值
    EpochObs        历元的观测数据, 以此顺序计算卫星轨道
    PreOrb          预报的星载接收机位置[m],预报的星载接收机钟差[m](GPS, GLONASS)

  输出参数
  
    SatMidInfo      卫星轨道, 速度等中间计算过程信息, 用于滤波定轨
    Result          单点定位结果
  
  返回值
    
    参与定位的卫星数

**************************************************************************
int PointPositionRTODInECI( GPSEPHREC* GPSEph, IONOPARA* IonPara,
                      GLONASSEPHREC* GLOEph, GLOTIMECORR* GloTmCorr,
                      EpochObsData* Epoch, APRIORBIT* PreOrb, 
                      SATMIDRESULT* SatMidInfo, PPRESULT* Result );
*/

/***************************************************************************
//
// PointPositionVelocityDetermination
//
// 目的: 根据卫星星历和一个历元的观测数据的D1, 进行单点测速(包括仅GPS, 
         仅GLONASS, GPS/GLONASS组合定位), 为滤波提供初始速度信息. 调试情况下,
         CHAMP卫星等数据没有多普勒观测数据, 不能用该子程序.

        GPS和GLONASS时间系统之间的钟差变化率在此被认为是相同的,
        所以在参数估计时只设置一个频差参数, 估计的参数是4维的.

  使用方法: 该子程序只能在PointPositionRTOD程序后调用, 因为它用到SatMidInfo
            和单点定位的位置信息.

//
// 输入参数:
//
    EpochObs        历元的观测数据, 
    GLOEph          GLONASS卫星星历, 从中获取每颗卫星的频率序数
    SatMidInfo      卫星轨道, 速度等中间计算过程信息
    Result          单点定位结果, 其中包含单点定位的坐标信息
  
  返回值
    
    参与测速的卫星数

***************************************************************************/
int PointPositionVelocityDetermination( EpochObsData* Epoch, 
                                       SATMIDRESULT* SatMidInfo, PPRESULT* Result );


/***************************************************************************
//
// DetectPseudoRangeOutliers
//
// 目的: 根据先验坐标或最小二乘估计的近似坐标, 计算出所有卫星的O-C
         (观测值减去计算值), 比较每颗卫星计算出来的接收机钟差是否一致.
         如果该历元包含GPS和GLONASS数据, 则分开探测, 因为两者计算的接收机钟差不一致.
         
//
// 输入参数:
//
    Accuracy        先验轨道数据的精度指标, 用来确定粗差剔除的阈值.
    O_C             伪距观测值减去观测值的计算值, 即每颗卫星的接收机钟差.
    Epoch           历元观测数据, 探测到粗差后将该卫星设置为粗差观测值,不可用
    SatList         与O_C相对应, 以确知粗差观测值来自于那颗卫星

  返回参数

  当探测到粗差时, 返回1, 否则返回0.

***************************************************************************/
bool DetectPseudoRangeOutliers( int graceType, double Accuracy, double O_C[], 
                              EpochObsData* Epoch, SATLIST* SatList );
//LC探测粗差
bool DetectOutliers_rel( double Accuracy, double O_C[], Common11Obs * CurComObs, SATLIST* SatList, double sigma );
void DeleteOutlierFrom_rel( int* GPSSatNum, double Resid[], double A[], double P[], SATLIST* SatList );
/***************************************************************************
//
// DeleteOutlierFromObsList
//
// 目的: 当探测出有粗差观测值后, 此函数将删除粗差卫星的观测矩阵/加权系数/残差数据等,
         避免粗差数据参与最小二乘估计
         
//
// 输入参数:
//
    GPSSatNum        GPS卫星数.
    GLOSatNum        GLONASS卫星数 
    Resid            观测值-计算值
    A                观测系数矩阵
    P                观测值权矩阵
    SatList          与O_C相对应, 以确知粗差观测值来自于那颗卫星

***************************************************************************/
void DeleteOutlierFromObsList( int* GPSSatNum, double Resid[],
                                double A[], double P[], SATLIST* SatList );



/***************************************************************************
//
// CheckOutPPPostResidual
//
// 目的: 在单点定位之后, 对验后残差进行假设检验. 检查结果的正确性和每颗卫星的残差
         是否合理. 如果残差被确认为粗差, 将不参与滤波定轨.
         
//
// 输入参数:
//
    n                历元卫星总数
	Residual         单点定位观测值残差
	Qvv              观测值残差改正数
    Result           单点定位的结果, 探测到粗差数据, 设置Result中SatList参数
***************************************************************************/
void CheckPostResidual_W( int n, double sigma0, double Residual[], double Qvv[], PPRESULT* Result );


/***************************************************************************
//
// EKFilterRTOD
//
// 目的: 在实时单点定位的基础上，用推广卡尔曼滤波对每颗导航卫星的观测数据进行
         测量更新。

  滤波定轨流程：

  初始化卫星位置和速度，转换到地心惯性系

  动力学轨道积分预报
  
  加密轨道内插

  滤波的时间更新

  单点定轨，地固系与惯性系的转换矩阵
  滤波测量更新（直接在惯性系中更新）

//
// 输入参数:
//

    Mat[9]          地固系转换到惯性系的转换矩阵
	SatMidInfo      卫星轨道, 速度等中间计算过程信息, 用于滤波定轨
    EpochObs        历元的观测数据, 以此顺序计算卫星轨道
	PPResult        单点定位结果

  输出参数
  
    KFState         滤波定轨结果
  
  返回值
    
    参与定位的卫星数

***************************************************************************/
int EKFilterRTODA( double Mat[9], SATMIDRESULT* SatMidInfo_A, EpochObsData* Epoch_A,
				 PPRESULT* PPResult, EKFSTATE* KFState );

/***************************************************************************
//
// EKFConvergencyCheck
//
// 目的: 检查卡尔曼滤波是否发散,用来决定滤波是否重置.

   方法: 
   
	 1. 根据滤波更新前的残差平方和与滤波更新后的残差平方和的特性是否满足假设检验;
	 2. 单点定位的轨道(在较好的状态下)与滤波更新的轨道比较, 看两者是否一致.

  更新的卫星数, 更新前的sigma, 更新后的sigma
	 
//
// 输入参数:
//

	PPResult        单点定位结果
	KFState         滤波估计结果

  输出参数

    设置KFState中的KFConvergency参数.

***************************************************************************/
void EKFConvergencyCheck( int graceType, EKFSTATE* KFState );
/*计算平差后观测值协因数阵*/
void ComputeQLL( int m, int n, const double A[],const double InvN[],const double AT[], const double Pll[], double QLL[] );
void initKFState(EKFSTATE *KFState);
void freeKFState(EKFSTATE KFState);
/*通过精密星历计算卫星位置和速度*/
void computeGPSfromSP3(const GPSTIME T_Tm,const int Prn,const int order, double * Pos0,double * Vel0,double * ClkOft,double * ClkSft);
#endif