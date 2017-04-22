/****************************************************************************
目的：    从Rinex格式的数据中读取GPS和GLONASS星历数据, 计算导航卫星
          的轨道和钟差

编写时间：2008.12.4
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _EPH_READ_PROCESS_H_
#define _EPH_READ_PROCESS_H_

#include <stdio.h>
#include "GPSTime.h"
#include "RTOD_Structs.h"


/***************************************************************************
//
// ReadGPSEph
//
// 目的: 读取观测时刻最接近的某颗GPS卫星星历.
//
// 输入参数:
//
//  Prn   GPS卫星号
//  Time  观测时刻
//  IonPara  电离层延迟改正参数
//
//  输出参数:
//  GPSEphRecord   GPS卫星星历结构体, 存储最接近时刻的卫星数据
//
//  返回值:   
// 
     如果文件读取有问题, 返回0, 否则返回1
***************************************************************************/

bool ReadGPSEph( const short Prn, const GPSTIME* Time,
               IONOPARA* IonPara, GPSEPHREC* GPSEphRecord );

/***************************************************************************
//
// CheckGPSEph
//
// 目的: 检查星历记录中星历数据是否为当前观测时刻的GPS卫星星历, 
         如果不是最新, 从星历数据文件读取
//
// 输入参数:
//
//  prn   GPS卫星号
//  Time  观测时刻

  //  输出参数:

//  IonPara  电离层延迟改正参数
//  GPSEphRecord   GPS卫星星历, 存储最接近时刻的卫星数据
//
//  返回值:   

// 如果星历是当前时刻, 返回true, 否则返回false
***************************************************************************/

bool CheckGPSEph( const short prn, const GPSTIME* Time,
               IONOPARA* IonPara, GPSEPHREC* GPSEphRecord );

/***************************************************************************
//
// ComputeGPSOrbitAndClockFullInfo
//
// 目的: 根据卫星星历, 计算卫星位置/速度/钟差/钟速全部信息
//
// 输入参数:
//
//  Prn   待计算的GPS卫星号
//  Time  信号发射时刻的GPS时间
//  Eph   GPS卫星星历
//  IonPara  电离层延迟改正参数
//
//  输出参数:
//  Pos, Vel 信号发射时刻的GPS卫星位置和速度[m, m/s]
//  clkoft, ClkSft  信号发射时刻的GPS卫星钟差和钟速[s, s/s]
//
返回值

  计算成功返回true, 否则返回false.
***************************************************************************/
bool ComputeGPSOrbitAndClockFullInfo( const int Prn, const GPSTIME* t,
                                     GPSEPHREC* Eph, IONOPARA* IonPara, 
                                     double Pos[3], double Vel[3],
                                     double* ClkOft, double* ClkSft );

/***************************************************************************
//
// ComputeGPSOrbitAndClockPartInfo
//
// 目的: 根据卫星星历, 计算卫星位置/钟差信息
//
// 输入参数:
//
//  Prn   待计算的GPS卫星号
//  Time  信号发射时刻的GPS时间
//  Eph   GPS卫星星历
//  IonPara  电离层延迟改正参数
//
//  输出参数:
//  Pos     信号发射时刻的GPS卫星位置和速度[m ]
//  clkoft  信号发射时刻的GPS卫星钟差和钟速[m ]
//
   返回值

  计算成功返回true, 否则返回false.
  ***************************************************************************/

bool ComputeGPSOrbitAndClockPartInfo( const int PRN, const GPSTIME* t,
                                     GPSEPHREC* Eph, IONOPARA* IonPara, 
                                     double Pos[3], double* ClkOft );

/***************************************************************************
//
// ReadGLONASSEph
//
// 目的: 读取观测时刻最接近的GLONASS卫星星历, 保存到星历记录数组.
//
// 输入参数:
//
//  slot  GLONASS卫星号
//  Time  观测时刻
//  GloTmCorr GLONASS时间与UTC(SU)系统间的改正值
//  
//  输出参数:
//  GLOEphRecord   GLONASS卫星星历, 存储最接近时刻的卫星数据
//
//  返回值:   
// 如果文件读取有问题, 返回0, 否则返回1
***************************************************************************/

bool ReadGLONASSEph( const short Slot, const GPSTIME* Time,
                    GLOTIMECORR* GloTmCorr, GLONASSEPHREC* GLOEphRecord );

/***************************************************************************
//
// CheckGLOEph
//
// 目的: 检查星历记录数据是否为当前观测时刻的GLONASS卫星星历, 如果不是最新,
         从星历数据文件读取.
//
// 输入参数:
//
//  slot  GLONASS卫星号
//  Time  观测时刻(UTC)
//
//  输出参数:

  //  GloTmCorr      GLONASS时间与UTC(SU)系统间的改正值
  //  GLOEphRecord   GLONASS卫星星历, 存储最接近时刻的卫星数据
  //
//  返回值:   
// 如果文件读取有问题, 返回0, 否则返回1
***************************************************************************/

bool CheckGLOEph( int Slot, const GPSTIME* Time, GLOTIMECORR* GloTmCorr, 
                 GLONASSEPHREC* GLOEphRecord );

/***************************************************************************
//
// ComputeRightFunc
//
// 目的: GLONASS卫星轨道的积分右函数, 用于计算卫星加速度和速度
//
// 输入参数:
//
//  Pos   当前时刻的卫星位置
//  Vel   当前时刻的卫星速度
//  Acc   当前时刻的卫星加速度

   输出参数:
   dY[6]  下一时刻的卫星速度和加速度

***************************************************************************/
void ComputeRightFunc( double Pos[3], double Vel[3], double Acc[3],
                      double dY[6] );


/***************************************************************************
//
// ComputeGLONASSSOrbit
//
// 目的: 根据卫星星历, 计算卫星位置/速度
//
// 输入参数:
//
//  Slot  待计算的GLONASS卫星号
//  Time  信号发射时刻的GLONASS系统时间
//  Eph   Slot卫星的GLONASS卫星星历
//  GloTmCorr GLONASS时间与UTC(SU)系统间的改正值
//  输出参数:
//  Pos, Vel 信号发射时刻的GPS卫星位置和速度[m, m/s]
//
   返回值

  计算成功返回true, 否则返回false.
***************************************************************************/
bool ComputeGLONASSSOrbit( const int Slot, const GPSTIME* t,
                             GLONASSEPHREC* Eph, GLOTIMECORR* GloTmCorr,
                             double Pos[3], double Vel[3] );


/***************************************************************************
//
// ComputeGLONASSSClockInfo
//
// 目的: 根据卫星星历, 计算钟差/钟速信息
//
// 输入参数:
//
//  Slot  待计算的GLONASS卫星号
//  Time  信号发射时刻的GLONASS系统时间
//  Eph   Slot卫星的GLONASS卫星星历
//  GloTmCorr GLONASS时间与UTC(SU)系统间的改正值

//  输出参数:
//  clkoft, ClkSft  信号发射时刻的GPS卫星钟差和钟速[m, m/s]

  返回值

  计算成功返回true, 否则返回false.
***************************************************************************/
bool ComputeGLONASSSClockInfo( const int Prn, const GPSTIME* t,
                             GLONASSEPHREC* Eph, GLOTIMECORR* GloTmCorr,
                             double* ClkOft, double* ClkSft );

#endif