/****************************************************************************
目的：    低轨卫星动力学模型, 包括重力场/日月质点力模型/大气阻力/固体潮
          摄动等

编写时间：2008.11.25
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _DYNAMODEL_H_
#define _DYNAMODEL_H_

#include "GPSTime.h"

/***************************************************************************
//
// Sun
//
// Purpose:
//
//   Computes the Sun's geocentric position using a low precision 
//   analytical series
//
// Input/Output:
//
//   Mjd_TT    Terrestrial Time (Modified Julian Date)
//   SunPos    Solar position [m] with respect to the 
//             mean equator and equinox of J2000 (EME2000, ICRF)
//
****************************************************************************/

void Sun ( const MJDTIME* Mjd_TT, double SunPos[] );

/***************************************************************************
//
// Moon
//
// Purpose:
//
//   Computes the Moon's geocentric position using a low precision
//   analytical series
//
// Input/Output:
//
//   Mjd_TT    Terrestrial Time (Modified Julian Date)
//   MoonPos   Lunar position [m] with respect to the 
//             mean equator and equinox of J2000 (EME2000, ICRF)
//
****************************************************************************/

void Moon ( const MJDTIME* Mjd_TT, double MoonPos[] );


/***************************************************************************
//
// AccelHarmonic
//
// Purpose:
//
//   Computes the acceleration due to the harmonic gravity field of the 
//   central body
//
// Input/Output:
//
//   Pos         Satellite position vector in the inertial system
//   E           Transformation matrix to body-fixed system
//   Para        Parameter of dynamic model 
//   Accel       Acceleration (a=d^2r/dt^2)
//
****************************************************************************/

void AccelHarmonic ( const double Pos[], const double E[9], 
                    DYNMODELPARA* Para, double Accel[3] );

/***************************************************************************
//
// AccelPointMass
//
// Purpose:
//
//   Computes the perturbational acceleration due to a point mass such as Sun 
     and Moon
//
// Input/Output:
//
//   Pos         Satellite position vector (r)
//   S           Point mass position vector (s)
//   GM          Gravitational coefficient of point mass
//   Accel       Acceleration (a=d^2r/dt^2)
//
****************************************************************************/

void AccelPointMass (const double Pos[3], const double S[], 
                     double GM, double Accel[3] );

/***************************************************************************
//
// AccelEarthTides
// 
//   Computer the perturbational acceleration due to solid earth tide.
//
//   Input/Output
//
//   Pos         Spacecraft position[m] 
//   SunPos      Sun position vector[m] 
//   MoonPos     Moon position vector[m] 
//   Accel       Acceleration (a=d^2r/dt^2)
****************************************************************************/

void AccelEarthTides( const double Pos[3], const double SunPos[3],
                     const double MoonPos[3], double Accel[3]);

/***************************************************************************
//
// AccelDrag
//
// Purpose:
//
//   Computes the acceleration due to the atmospheric drag.
//
// Input/Output:
//
//   Mjd_TT      Terrestrial Time (Modified Julian Date)
//   T           Transformation matrix to true-of-date inertial system
//   Pos         Satellite position vector in the inertial system [m]
//   Vel         Satellite velocity vector in the inertial system [m/s]
//   Para        Parameter of dynamic model 
//   Accel       Acceleration (a=d^2r/dt^2) [m/s^2]
//
****************************************************************************/

void AccelDrag ( const MJDTIME* Mjd_TT, const double T[9], 
                const double Pos[3], const double Vel[3], 
                DYNMODELPARA* Para, double Accel[3] );

/***************************************************************************
//
// Density_HP
//
// Purpose:
//
//   Computes the atmospheric density for the modified Harris-Priester model.
//
// Input/Output:
//
//   Mjd_TT      Terrestrial Time (Modified Julian Date)
//   r_tod       Satellite position vector in the inertial system [m]
//   <return>    Density [kg/m^3]
//
****************************************************************************/

double Density_HP ( const MJDTIME* Mjd_TT, const double Pos[3] );

/***************************************************************************
//
// AccelSolrad
//
// Purpose:
//
//   Computes the acceleration due to solar radiation pressure assuming 
//   the spacecraft surface normal to the Sun direction
//
// Input/Output:
//
//   SatPos         Satellite position vector in the inertial system [m]
     SunPos         Sun position vector in the inertial system[m]
//   Para           Parameter of dynamic model 
//   Acc            Acceleration (a=d^2r/dt^2) [m/s^2]
//
****************************************************************************/

void AccelSolrad ( const double SatPos[3], const double SunPos[3],
				  DYNMODELPARA* Para, double Acc[3]);

/***************************************************************************
//
// AccelMain
//
// Purpose:
//
//   Computes the acceleration of an Earth orbiting satellite due to 
//    - the Earth's harmonic gravity field, 
//    - the gravitational perturbations of the Sun and Moon
//    - the solid earth tide and
//    - the atmospheric drag
//
// Input/Output:
//
//   Mjd_TT      Terrestrial Time (Modified Julian Date)
//   Pos         Satellite position vector in the ICRF/EME2000 system
//   Vel         Satellite velocity vector in the ICRF/EME2000 system
//   E           Transformation matrix to body-fixed system
//   T           Transformation matrix to true-of-date inertial system
//   Para        Parameter of dynamic model 
//   Accel       Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
//
****************************************************************************/

void AccelMain( int graceType, const MJDTIME* Mjd_TT, const double Pos[3], const double Vel[3],
                const double E[9], const double T[9],
                DYNMODELPARA* Para, double Accel[3] );



/***************************************************************************
//
// Gradient
//
// Purpose:
//
//   Computes the gradient of the Earth's harmonic gravity field 
//
// Input/Output:
//
//   Mjd_UT      Modified Julian Date (Universal Time)
//   E           Transformation matrix to body-fixed system
//   Pos         Satellite position vector in the ICRF/EME2000 system
//   Para        Parameter of dynamic model 
//   Grad        da/dr
//
****************************************************************************/

void Gradient ( const MJDTIME* Mjd_TT, const double E[9], const double Pos[3],
               DYNMODELPARA* Para, double Grad[9] );


/***************************************************************************
//
// GradientCd
//
// Purpose:
//
//   Computes the gradient of the Earth's atmospheric drag
//
// Input/Output:
//
//   Mjd_UT      Modified Julian Date (Universal Time)
//   T           Transformation matrix to true-of-date inertial system
//   Pos         Satellite position vector in the ICRF/EME2000 system
//   Para        Parameter of dynamic model 
//   Grad        da/dCd
//
****************************************************************************/
void GradientCd (int graceType, const MJDTIME* Mjd_TT, const double T[9], const double Pos[9],
                 DYNMODELPARA* Para, double GradCd[3] );

/***************************************************************************
//
// GradientCr
//
// Purpose:
//
//   Computes the gradient of solar radiation pressure.
//
// Input/Output:
//
//   SatPos         Satellite position vector in the ICRF/EME2000 system
//   SunPos         Satellite velocity vector in the ICRF/EME2000 system
//   Para           Parameter of dynamic model 
//   Grad           da/dCd
//
****************************************************************************/
void GradientCr( int graceType, const MJDTIME* Mjd_TT, const double SatPos[3], const double SunPos[3],
				 DYNMODELPARA* Para, double GradCr[3] );

/***************************************************************************
//
// VarEquation
//
// Purpose:
//
//   计算卫星运动加速度和变分方程, 积分计算轨道和状态转移矩阵
//   
//   Phi_dot = F(t) * Phi
//
// Input/Output:
//
//   Mjd_GPS     Modified Julian Date (GPS Time)
//   Phi         包括位置/速度/状态转移矩阵[6+6*8维]
//   dPhi        Phi的导数
//   Para        Parameter of dynamic model 
****************************************************************************/
void VarEquation(int graceType, const MJDTIME* Mjd_GPS, const double Phi[54], 
                 double dPhi[54], DYNMODELPARA* Para );

#endif
