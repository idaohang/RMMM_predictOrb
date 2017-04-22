/****************************************************************************
目的：    定义LEOSAOD软件需要的常量参数
编写时间：2008.11.22
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _RTOD_CONST_H_
#define _RTOD_CONST_H_

/* Mathematical constants  */

#define IonFree1 2.54572778016316
#define IonFree2 1.54572778016316
#define pi 3.1415926535897932384626433832795
#define pi2 (2.0*pi)                    /* 2pi */
#define Rad (pi/180.0)                  /* Radians per degree */
#define Deg (180.0/pi)                  /* Degrees per radian */
#define Arcs (3600.0*180.0/pi)          /* Arcseconds per radian */

/* GPS time constants  */
#define JAN61980  44244                 /* MJD of 1980.1.6 */
#define JAN11901  15385                 /* MJD of 1901.1.1 */
#define SECPERHOUR 3600.0;                /* Seconds per hour */
#define SECPERDAY  86400.0                /* Seconds per day */
#define SECPERWEEK 604800.0               /* Seconds per week */

/* General constants */

#define MJD_J2000  51544.5       /* Modif. Julian Date of J2000.0 */
#define AU 149597870000.0        /* Astronomical unit [m]; IAU 1976  */
#define C_Light 299792458.0      /* Speed of light  [m/s]; IAU 1976  */


/* Physical parameters of the Earth, Sun and Moon  */

#define R_WGS84  6378137.0          /* Radius Earth [m]; WGS-84  */
#define F_WGS84  1.0/298.257223563  /* Flattening; WGS-84   */
#define R_Sun    696000.0e3         /* Radius Sun [m]  */
#define R_Moon   1738.0e3           /* Radius Moon [m] */
#define Omega_WGS 7.2921151467e-5   /*[rad/s], the earth rotation rate */

/* Geodetic constants and parameters of PZ-90.02 common terrestrial ellipsoid*/

#define R_PZ90 6378136.0          /* Radius Earth [m]; PZ-90.02 from ICD  */
#define F_PZ90 1/298.257839303    /* Flattening; PZ-90.02 from ICD   */
#define Omega_PZ90 7.292115e-5    /* Earth rotation rate[Rad/s]   */
#define GM_PZ90    398600.44e9    /* Gravitational constant */
#define Geopo_J02  1082625.7e-9   // 1082625.7E-9  /*Second zonal harmonic of the geopotential */



/* Gravitational coefficient  */

#define GM_Earth   398600.4415e+9     /* [m^3/s^2]; JGM3  */
#define R_Earth    6378136.3          /* Radius Earth [m]; JGM3  */
#define GM_Sun     1.32712438e+20     /* [m^3/s^2]; IAU 1976 */
#define GM_Moon  GM_Earth/81.300587   /* [m^3/s^2]; DE200  */


/* Solar radiation pressure at 1 AU  */

#define P_Sol    4.560E-6           /* [N/m^2] (~1367 W/m^2); IERS 96 */

/* some constants about GPS satellite signal */
#define  OSC_Freq  10.23E6
#define  FG1_Freq  (154.0 * OSC_Freq)    /* L1信号频率 */
#define  FG2_Freq  (120.0 * OSC_Freq)    /* L2信号频率 */
#define  FG12Ratio ( 77/60.0 )           /* FG1_Freq/FG2_Freq */
#define  WAVELENGTHL1 (C_Light/FG1_Freq)
#define  WAVELENGTHL2 (C_Light/FG2_Freq)

/* some constants about GLONASS satellite signal */
#define  FG01_GLO  1602E6               /* L1信号的基准频率 */
#define  DFG1_GLO  562.5E3              /* L1信号的频率增量 */
#define  FG02_GLO  1246E6               /* L2信号的基准频率 */
#define  DFG2_GLO  437.5E3              /* L2信号的频率增量 */
#define  FR12Ratio ( 9/7.0 )            /* FG01_GLO/FG02_GLO */

/* Receiver clock stochastic process parameters  */

/*参考《Precise Relative Positioning of Formation Flying Spacecraft using GPS》*/
#define  t_nsd	10E16					 /* clock noise spectral density*/

#define  Sf		0.0089875517873681764    /* clock offset PSD */
#define  Sg		0.0354814322702509941    /* clock shift PSD */
#define  Svtec	0.00000001				/*vtec方差*/


/* Orbit parameter translation consts */

#define AE  6.378140e6
#define CKV 0.12649638204184154e-3


#endif  

