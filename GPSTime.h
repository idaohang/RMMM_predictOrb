/****************************************************************************
目的：    定义时间结构体及其相互转换函数,
          空间直角坐标和大地坐标相互转换
编写时间：2008.11.22
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _GPS_TIME_H_
#define _GPS_TIME_H_

#include "RTOD_Structs.h"

/* 通用时,GPS时和简化儒略日之间的相互转换函数*/

void CommonTimeToMJDTime( const COMMONTIME* CT, MJDTIME* MJDT);
void MJDTimeToCommonTime( const MJDTIME* MJDT, COMMONTIME* CT );
void GPSTimeToMJDTime( const GPSTIME* GT, MJDTIME* MJDT );
void MJDTimeToGPSTime ( const MJDTIME* MJDT, GPSTIME* GT );
void CommonTimeToGPSTime ( const COMMONTIME* CT, GPSTIME* GT );
void GPSTimeToCommonTime ( const GPSTIME* GT, COMMONTIME* CT );

/* 检查GPS时间格式是否正确, 如周秒是否大于604800或小于0  */
void CheckGPSTime( GPSTIME* GT );

/* 空间直角坐标,大地坐标的相互转换函数 */

void XYZToBLH( const double xyz[], double blh[], double R, double F );
void BLHToXYZ( const double BLH[], double XYZ[], double R, double F );

/***************************************************************************
//
// BLHToNEUMatrix
//
// 目的: 地面测站点的NEU旋转矩阵计算

//
// 输入参数:
//
    BLH      地面测站点坐标[Rad, m]

  输出参数
    
    H        旋转矩阵

***************************************************************************/

void BLHToNEUMatrix( double BLH[3], double H[9]);


/* PZ90 与WGS84坐标转换 */

void CoorTranFromPZ90ToWGS84( double PZ90[3], double WGS84[3] );
void VeloTranFromPZ90ToWGS84( double PZ90[3], double WGS84[3] );

/***************************************************************************
//
// PhaseCentToMassCent
//
// 目的: 将GNSS接收机天线相位中心坐标转换到质量中心

//
// 输入参数:
//
   Flag      1为相位中心转换到质量中心, 0为质量中心转换相位中心
   Bias      星固系中的GNSS接收机天线相位中心与质量中心的偏差参数[m]

   输出参数
    
   State[6]  输入为星载接收机的位置和速度, 输出为经过改正后的位置

***************************************************************************/

void PhaseCentToMassCent( bool Flag, const double Bias[3], double State[6] );


/***************************************************************************
//
// XYZToRTN
//
// 目的: 将空间直角坐标系下的坐标分量转换为RTN轨道坐标系的分量

   说明: RTN指的是径向、切向和法向

//
// 输入参数:
//
   State     卫星的运动状态，包括位置和速度，用于定义RTN坐标系
   dXYZ      空间直角坐标系下的指标分量
   dRTN      RTN坐标系下的坐标分量

***************************************************************************/

void XYZToRTN( double State[6], double dXYZ[3], double dRTN[3] );
//生成矩阵
void MatrixXYZ2RTN(double State[6],double transMat[9],int flag);
/*从satods结构体GPSTIME转到rtklib通用时ep,为了利用rtklib里的精密星历相位中心改正模块*/
void GPSTimeToep ( const GPSTIME* GT, double ep[6]);
void SecTimeToCT(double sec,COMMONTIME * CT);
#endif