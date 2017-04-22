/****************************************************************************
目的：    定义定轨系统需要的时间系统和坐标系统及其相互转换函数,
          包括世界时/恒星时/GPS时/动力学时等时间系统;J2000惯性系
          /地心地固系等坐标系统的转换

编写时间：2008.11.22
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _REFSYS_H_
#define _REFSYS_H_

#include "RTOD_Structs.h"


double TT_UTC( void );
void   SetGPST_UTC( double val );
double GetGPST_UTC();
int    InitEOPPara( const MJDTIME* CurrTime );
void   InterposeEOP(const MJDTIME* time, EOPPARA* CurrEop);

/****************************************************************************
CheckLeapSecInEOP

  目的：检查EOP[3]中的LeapSec是否连续，如果不连续，进行调整。
  
****************************************************************************/

void CheckLeapSecInEOP();


/***************************************************************************
//
// MeanObliquity
//
// Purpose:
//
//   Computes the mean obliquity of the ecliptic
//
// Input/Output:
//
//   Mjd_TT    Modified Julian Date (Terrestrial Time)
//   <return>  Mean obliquity of the ecliptic,[Rad]
//
***************************************************************************/
double MeanObliquity ( const MJDTIME* Mjd_TT); 


/***************************************************************************
//
// EclMatrix
//
// Purpose:
//
//   Transformation of to ecliptical coordinates
// equatorial
// Input/Output:
//
//   Mjd_TT    Modified Julian Date (Terrestrial Time)
//   <return>  Transformation matrix
//
***************************************************************************/

void EclMatrix ( const MJDTIME* Mjd_TT, double Mat[]);  


/***************************************************************************
//
// PrecMatrix
//
// Purpose:
//
//   Precession transformation of equatorial coordinates
//
// Input/Output:
//
//   Mjd_1     Epoch given (Modified Julian Date TT)
//   MjD_2     Epoch to precess to (Modified Julian Date TT)
//   <return>  Precession transformation matrix
//
***************************************************************************/

void PrecMatrix (double Mjd_1, const MJDTIME* Mjd_2, double Mat[]);

//------------------------------------------------------------------------------
//
// NutAngles 
//
// Purpose:
//
//   Nutation in longitude and obliquity
//
// Input/Output:
//
//   Mjd_TT    Modified Julian Date (Terrestrial Time)
//   <return>  Nutation matrix
//
//------------------------------------------------------------------------------

void NutAngles (const MJDTIME* Mjd_TT, double* dpsi, double* deps);

/***************************************************************************
//
// NutMatrix 
//
// Purpose:
//
//   Transformation from mean to true equator and equinox
//
// Input/Output:
//
//   Mjd_TT    Modified Julian Date (Terrestrial Time)
//   <return>  Nutation matrix
//
***************************************************************************/

void NutMatrix ( const MJDTIME* Mjd_TT, double Mat[]);


/***************************************************************************
//
// NutMatrixSimple 
//
// Purpose:
//
//   Transformation from mean to true equator and equinox (low precision)
//
// Input/Output:
//
//   Mjd_TT    Modified Julian Date (Terrestrial Time)
//   <return>  Nutation matrix
//
***************************************************************************/
void NutMatrixSimple ( const MJDTIME* Mjd_TT, double Mat[]);

/***************************************************************************
//
// EqnEquinox 
//
// Purpose:
//
//   Computation of the equation of the equinoxes
//
// Input/Output:
//
//   Mjd_TT    Modified Julian Date (Terrestrial Time)
//   <return>  Equation of the equinoxes
//
// Notes:
//
//   The equation of the equinoxes dpsi*cos(eps) is the right ascension of the 
//   mean equinox referred to the true equator and equinox and is equal to the 
//   difference between apparent and mean sidereal time.
//
***************************************************************************/

double EqnEquinox (const MJDTIME* Mjd_TT);


/***************************************************************************
//
// GMST
//
// Purpose:
//
//   Greenwich Mean Sidereal Time
//
// Input/Output:
//
//   Mjd_UT1   Modified Julian Date UT1
//   <return>  GMST in [rad]
//
***************************************************************************/

double GMST ( const MJDTIME* Mjd_UT1 );


/***************************************************************************
//
// GAST
//
// Purpose:
//
//   Greenwich Apparent Sidereal Time
//
// Input/Output:
//
//   Mjd_UT1   Modified Julian Date UT1
//   <return>  GMST in [rad]
//
***************************************************************************/

double GAST ( const MJDTIME* Mjd_UT1);


/***************************************************************************
//
// GHAMatrix
//
// Purpose:
//
//   Transformation from true equator and equinox to Earth equator and 
//   Greenwich meridian system 
//
// Input/Output:
//
//   Mjd_UT1   Modified Julian Date UT1
//   <return>  Greenwich Hour Angle matrix
//
***************************************************************************/

void GHAMatrix ( const MJDTIME* Mjd_UT1, double Mat[] );


/***************************************************************************
//
// PoleMatrix
//
// Purpose:
//
//   Transformation from pseudo Earth-fixed to Earth-fixed coordinates
//   for a given date
//
// Input/Output:
//
//   x, y      Pole Motion para[Rad]
//   <return>  Pole matrix
//
***************************************************************************/

void PoleMatrix ( double x, double y, double Mat[]);

/***************************************************************************
//
// ICRF_ITRF_MJD
//
// 目的: 天球系与地固系坐标间的相互转换,其中包含位置和速度分量, 输入时间
         为MJD表示法的TT时.
//
// 参数:
//
//   Mjd1     Epoch given (Modified Julian Date TT)
//   Mjd2     Epoch to precess to (Modified Julian Date TT)
//   flag     1: 从天球系到地固系 
              0: 从地固系到天球系
//   ICRF     天球系坐标，其中包含位置和速度
//   ITRF     地固系坐标
//
***************************************************************************/
void ICRF_ITRF_MJD( const double Mjd1, const MJDTIME* Mjd2,
               int flag, double ICRF[6], double ITRF[6] );
			   
/***************************************************************************
//
// ICRF_ITRF_GPST
//
// 目的: 天球系与地固系坐标间的相互转换,其中包含位置和速度分量, 输入时间
         为MJD表示法的TT时.
//
// 参数:
//
//   Mjd1     Epoch given (Modified Julian Date TT)
//   CT       GPS时间
//   flag     1: 从天球系到地固系 
              0: 从地固系到天球系
//   ICRF     天球系坐标，其中包含位置和速度
//   ITRF     地固系坐标
//
***************************************************************************/
void ICRF_ITRF_GPST( const double Mjd1, const GPSTIME* GT,
               int flag, double ICRF[6], double ITRF[6] );


/***************************************************************************
//
// ICRF_ITRF_Matrix
//
// 目的: 天球系与地固系坐标间的相互转换矩阵,其中包含位置和速度分量.
         时间以GPS时.
//
// 参数:
//
//   Mjd1     Epoch given (Modified Julian Date TT, MJD_J2000)
//   CT       GPS时间
//   flag     1: 从天球系到地固系 
              0: 从地固系到天球系
//   Mat      转换矩阵
//
***************************************************************************/
void ICRF_ITRF_Matrix( const double Mjd1, const GPSTIME* CT,
               int flag, double Mat[9] );



#endif