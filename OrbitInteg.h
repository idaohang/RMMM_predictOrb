/****************************************************************************
目的：    动力学轨道积分,使用RK4, RKF4,以及轨道差值等函数

编写时间：2008.11.26
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _ORBIT_INTEGRATION_H_
#define _ORBIT_INTEGRATION_H_

#include "DynaModel.h"
#include "GPSTime.h"



/***************************************************************************
//
// RK4Step
//
// Purpose:
//
//   龙格-库塔4阶单步积分器函数
//
// Input/Output:
//
//   Mjd_GPS    初始时刻MJD时间(GPS Time),积分后为mjd_gps+step
//   Step       积分时间间隔[s]
//   Y0         积分初值,Mjd_GPS时刻的卫星位置和速度[m, m/s],
//                 积分后为下一时刻的卫星位置和速度
//   Para       动力学模型参数
//
***************************************************************************/
//void RK4Step( MJDTIME* Mjd_GPS, double step, double Y0[6], DYNMODELPARA* Para );

/***************************************************************************
//
// RKF4Step
//
// Purpose:
//
//   龙格-库塔-Felberg4阶单步积分器函数
//
// Input/Output:
//
//   Mjd_GPS    初始时刻MJD时间(GPS Time),积分后为mjd_gps+step
//   Step       积分时间间隔[s]
//   Y0         积分初值,Mjd_GPS时刻的卫星位置和速度[m, m/s],
//                 积分后为下一时刻的卫星位置和速度
//   Para       动力学模型参数
//
***************************************************************************/
//void RKF4Step( MJDTIME* Mjd_GPS, double step, double Y0[6], DYNMODELPARA* Para );

/***************************************************************************
//
// RKF4OrbitSTM
//
// Purpose:
//
//   龙格-库塔-Felberg4阶单步积分器函数, 进行轨道和状态转移矩阵积分
//
// Input/Output:
//
//   Mjd_GPS    初始时刻MJD时间(GPS Time),积分后为mjd_gps+step
//   Step       积分时间间隔[s]
//   Y0         轨道和转移矩阵初值,Mjd_GPS时刻的卫星位置和速度[m, m/s],
//                 积分后为下一时刻的卫星位置和速度[48维]
//   IntState   用Hermit5方法内插所需要的状态参数，在此处初始化[位置/速度/加速度]
//   Para       动力学模型参数
//
***************************************************************************/
void RKF4OrbitSTM( int graceType, MJDTIME* Mjd_GPS, double step, double Y0[54], SCState Stat[2], DYNMODELPARA* Para );

/***************************************************************************
//
// Hermite5
//
// Purpose:
//
//   Hermite 5阶多项式内插卫星的轨道
//
// Input/Output:
//
//   S0         起始时刻卫星的状态参数[位置/速度/加速度]
//   S1         终止时刻卫星的状态参数
//   CurrState  内插时刻的卫星状态
****************************************************************************/
void Hermite5( const SCState* S0, const SCState* S1, SCState* CurrState );

/***************************************************************************
//
// OrbitIntegToGivenTime
//
// Purpose:
//
//   使用RKF4单步积分器函数，进行轨道积分，预报任意时刻卫星的轨道
//
// Input/Output:
//
//   Mjd_GPS             初始时刻是以MJD时间表示的GPS时
//   Mjd_GivenTime       需要预报轨道的观测时刻，与Mjd_GPS的意义相同
//   Y0                  积分初值,Mjd_GPS时刻的卫星位置和速度[m, m/s],
//                              积分后为下一时刻的卫星位置和速度
//   Para                动力学模型参数
//
***************************************************************************/
//void OrbitIntegToGivenTime( MJDTIME* Mjd_GPS, MJDTIME* Mjd_GivenTime,
//						   double Y0[6], DYNMODELPARA* Para );


/***************************************************************************
//
// InitStateTranMatrix
//
// Purpose:
//
//   初始化轨道积分向量中的状态转移矩阵
//
// Input/Output:
//
//   row    状态转移矩阵的行数
     col    状态转移矩阵的列数
	 STM    状态转移矩阵向量

****************************************************************************/

void InitStateTranMatrix( int row, int col, double STM[] );


#endif
