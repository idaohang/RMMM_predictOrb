/****************************************************************************
RefSys.c
目的：    定义定轨系统需要的时间系统和坐标系统及其相互转换函数,
包括世界时/恒星时/GPS时/动力学时等时间系统;J2000惯性系
/地心地固系等坐标系统的转换

编写时间：2008.11.22
版本:     V1.1
版权：    武汉大学
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RefSys.h"
#include "GPSTime.h"
#include "RTOD_Const.h"
#include "commonfuncs.h"



static double GPST_UTC = 10.0;        /* 从1980年1月1日开始, 跳秒数[s],在上星前要设置为正确的跳秒 */
static EOPPARA EOP[3] = { 0.0 };


/****************************************************************************
InitEOPPara
目的：从EOP参数文件中读取当前时刻的EOP参数,
要求EOP.txt文件必须在执行程序目录下,且文件名固定.

参数：
CurrTime   当前时刻的MJD
EOP[2]     为全局静态变量, 当前时刻前后的EOP参数

返回值
正确读取相应时刻的EOP参数,返回1;
文件中没有相应时刻EOP参数,返回0;

****************************************************************************/
int InitEOPPara( const MJDTIME* CurrTime )
{
	char line[256];
	int index;
	double CurrMjd;
	static bool Initialized = false;
	FILE *in;
	double LeapSec = 0.0;

	EOPPARA TempEOP;

	CurrMjd = CurrTime->Days + CurrTime->FracDay;

	printf( "Initialize EOP again.\n" );

	if( (in  = fopen( "EOP.txt", "rt" )) == NULL )
	{
		printf( "The file 'EOP.txt' was not opened\n" );
		exit(0);
	}

	index = 0;
	if( Initialized == false )   /* EOP参数为空 */
	{
		while( !feof( in ) )
		{
			if( fgets( line, 256, in ) != NULL )
			{
				sscanf( line, "%lf  %lf  %lf  %lf", 
					&TempEOP.Mjd, &TempEOP.x, &TempEOP.y, &TempEOP.dUT1 );

				TempEOP.LeapSec = GPST_UTC;  /* 先用星历文件初始化，或者在上星工作前设置正确的初始值 */
				TempEOP.Status  = 1;

				if( ( (CurrMjd - TempEOP.Mjd) > 1.0E-5 ) 
					&&( (CurrMjd - TempEOP.Mjd) < 1.1 ) )
				{
					EOP[0] = TempEOP;   /* Set the first EOP struct */

					if( fgets( line, 256, in ) != NULL )    /* Set the second EOP struct */
					{
						sscanf( line, "%lf  %lf  %lf  %lf", 
							&EOP[1].Mjd, &EOP[1].x, &EOP[1].y, &EOP[1].dUT1 );

						EOP[1].LeapSec = EOP[0].LeapSec;
						EOP[1].Status  = 1;

					}
					else
					{
						printf( "No new EOP data.\n" );
						fclose( in );
						return 0;
					}

					if( fgets( line, 256, in ) != NULL )   /* Set the third EOP struct */
					{
						sscanf( line, "%lf  %lf  %lf  %lf", 
							&EOP[2].Mjd, &EOP[2].x, &EOP[2].y, &EOP[2].dUT1 );

						EOP[2].LeapSec =  EOP[0].LeapSec;;
						EOP[2].Status  = 1;

						Initialized = true;
						CheckLeapSecInEOP();
						fclose( in );
						return 1;
					}
					else
					{
						printf( "No new EOP data.\n" );
						fclose( in );
						return 0;
					}
				}
			}
			else
			{
				EOP[0].Status = 0;
				printf( "No new EOP data.\n" );
				fclose( in );
				return 0;
			}
		}
	}
	else
	{
		EOP[0] = EOP[1];
		EOP[1] = EOP[2];

		while( !feof( in ) )
		{
			if( fgets( line, 256, in ) != NULL )
			{
				sscanf( line, "%lf  %lf  %lf  %lf", 
					&TempEOP.Mjd, &TempEOP.x, &TempEOP.y, &TempEOP.dUT1 );

				TempEOP.LeapSec = GPST_UTC;
				TempEOP.Status  = 1;

				if( fabs(EOP[2].Mjd - TempEOP.Mjd) < 0.1 )
				{
					if( fgets( line, 256, in ) != NULL )
					{
						sscanf( line, "%lf  %lf  %lf  %lf", 
							&EOP[2].Mjd, &EOP[2].x, &EOP[2].y, &EOP[2].dUT1 );

						EOP[2].Status  = 1;
						EOP[2].LeapSec = EOP[1].LeapSec;

						CheckLeapSecInEOP();

						fclose( in );
						return 1;

					}
					else
					{
						printf( "No new EOP data.\n" );
						fclose( in );
						return 0;
					}
				}
			}
			else
			{
				printf( "No new EOP data.\n" );
				fclose( in );
				return 0;
			}
		}
	}

	fclose( in );
	return 0;
}

/****************************************************************************
CheckLeapSecInEOP

目的：检查EOP[3]中的LeapSec是否连续，如果不连续，进行调整。

****************************************************************************/

void CheckLeapSecInEOP()
{
	double dSec;

	if( (EOP[0].Status==1) && (EOP[1].Status==1) && (EOP[2].Status==1) )
	{
		dSec = EOP[1].dUT1 - EOP[0].dUT1;

		if( fabs(dSec)>= 0.8 )                  /* EOP文件检查出跳秒 */
		{
			EOP[0].LeapSec = EOP[0].LeapSec + 1.0;
			EOP[0].dUT1    = EOP[0].dUT1 + 1.0;
			EOP[1].LeapSec = EOP[0].LeapSec;
		}

		dSec = EOP[2].dUT1 - EOP[1].dUT1;

		if( fabs(dSec)>= 0.8 )                  /* EOP文件检查出跳秒 */
		{
			EOP[0].LeapSec = EOP[0].LeapSec + 1.0;
			EOP[1].LeapSec = EOP[1].LeapSec + 1.0;
			EOP[0].dUT1    = EOP[0].dUT1 + 1.0;
			EOP[1].dUT1    = EOP[1].dUT1 + 1.0;
			EOP[2].LeapSec = EOP[1].LeapSec;

			printf( "Leap Second happened at Mjd: %12.1f\n", EOP[2].Mjd );

			GPST_UTC = EOP[0].LeapSec;
		}

	}

}

/****************************************************************************
InterposeEOP
目的：根据所读取的EOP参数,线性内插当前时刻的EOP参数

参数：
time       当前时刻的MJD
EOP[2]     为全局静态变量, 当前时刻前后的EOP参数
CurrEOP    插值结果
****************************************************************************/

void InterposeEOP(const MJDTIME* time, EOPPARA* CurrEop)
{
	int i, j;
	double t;
	double temp;

	t = time->Days + time->FracDay;

	if( fabs( t - EOP[1].Mjd ) > 0.8 )  /* time不在EOP[2]的时间段内 */
	{
		if( InitEOPPara( time ) == 0 )
		{
			printf( " No new EOP data.\n " );
			exit(0);
		}
	}

	CurrEop->x = 0.0;
	CurrEop->y = 0.0;
	CurrEop->dUT1 = 0.0;

	if( (EOP[0].Status==1) && (EOP[1].Status==1) && (EOP[2].Status==1) )
	{
		for( i=0; i<3; i++ )
		{
			temp = 1.0;

			for( j=0; j<3; j++ )
			{
				if( i==j )
				{
					continue;
				}

				temp = temp * (t-EOP[j].Mjd) / (EOP[i].Mjd-EOP[j].Mjd);
			}

			CurrEop->x = CurrEop->x + temp*EOP[i].x;
			CurrEop->y = CurrEop->y + temp*EOP[i].y;
			CurrEop->dUT1 = CurrEop->dUT1 + temp*EOP[i].dUT1;
		}
		CurrEop->Mjd = t;
		CurrEop->LeapSec = EOP[1].LeapSec;
		CurrEop->Status = 1;
	}
	else
	{
		printf( "No new EOP data.\n" );
		exit(1);
	}

}

/****************************************************************************
InterposeEOP
目的：根据所读取的EOP参数,线性内插当前时刻的EOP参数

参数：
time       当前时刻的MJD
EOP[2]     为全局静态变量, 当前时刻前后的EOP参数
CurrEOP    插值结果
****************************************************************************/

// void InterposeEOP(const MJDTIME* time, EOPPARA* CurrEop)
// {
//     double t;
// 	
//     t = time->Days + time->FracDay;
// 	
//    if( fabs( t - EOP[1].Mjd ) > 0.8 )  /* time不在EOP[2]的时间段内 */
//    {
// 		printf( "EOP initialized again!\n" );
//         if( InitEOPPara( time ) == 0 )
// 		{
// 			printf( " No new EOP data.\n " );
// 			exit(0);
// 		}
// 	}
// 	
// 	CurrEop->x = 0.0;
// 	CurrEop->y = 0.0;
// 	CurrEop->dUT1 = 0.0;
// 	
// 	if( (EOP[0].Status==1) && (EOP[1].Status==1) && (EOP[2].Status==1) )
// 	{
// 		if( t>=EOP[0].Mjd && t<EOP[1].Mjd )
// 		{
// 			CurrEop->Mjd = t;
// 			CurrEop->x   = EOP[0].x + (t-EOP[0].Mjd) * (EOP[1].x-EOP[0].x)
// 				/ (EOP[1].Mjd-EOP[0].Mjd);
// 			CurrEop->y   = EOP[0].y + (t-EOP[0].Mjd) * (EOP[1].y-EOP[0].y)
// 				/ (EOP[1].Mjd-EOP[0].Mjd);
// 			CurrEop->dUT1= EOP[0].dUT1 + (t-EOP[0].Mjd) * (EOP[1].dUT1-EOP[0].dUT1)
// 				/ (EOP[1].Mjd-EOP[0].Mjd);
// 			CurrEop->LeapSec = EOP[1].LeapSec;
// 			CurrEop->Status = 1;
// 		}
// 		if( t>=EOP[1].Mjd && t<EOP[2].Mjd )
// 		{
// 			CurrEop->Mjd = t;
// 			CurrEop->x   = EOP[1].x + (t-EOP[1].Mjd) * (EOP[2].x-EOP[1].x)
// 				/ (EOP[2].Mjd-EOP[1].Mjd);
// 			CurrEop->y   = EOP[1].y + (t-EOP[1].Mjd) * (EOP[2].y-EOP[1].y)
// 				/ (EOP[2].Mjd-EOP[1].Mjd);
// 			CurrEop->dUT1= EOP[1].dUT1 + (t-EOP[1].Mjd) * (EOP[2].dUT1-EOP[1].dUT1)
// 				/ (EOP[2].Mjd-EOP[1].Mjd);
// 			CurrEop->LeapSec = EOP[1].LeapSec;
// 			CurrEop->Status = 1;
// 		}
// 	}
// 	else
// 	{
// 		printf( "No new EOP data.\n" );
// 		exit(1);
// 	}
// 
// }

/****************************************************************************
TT_UTC
目的：计算某一时刻的TT时与UTC时的差值

返回值
TT-UTC    差值[s]
****************************************************************************/
double TT_UTC()
{
	return ( TT_TAI - GPST_TAI + GPST_UTC );
}

void SetGPST_UTC( double val )
{
	if( fabs( val-GPST_UTC ) > 2.0 )
		GPST_UTC = val;                /* 尽量不用观测文件更新，使用EOP文件中跳秒来更新 */
}

double GetGPST_UTC()
{
	return EOP[1].LeapSec;
}

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
double MeanObliquity ( const MJDTIME* Mjd_TT)
{
	double T, Val;

	T = ( Mjd_TT->Days + Mjd_TT->FracDay - MJD_J2000 ) / 36525.0;

	Val = Rad * ( 23.43929111 - 
		( 46.8150 + ( 0.00059 - 0.001813 * T ) *T ) *T/3600.0 );

	return ( Val );

}


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

void EclMatrix (MJDTIME* Mjd_TT, double Mat[])
{
	double Ang = MeanObliquity( Mjd_TT );

	Rotation_x( Ang, Mat );

}


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

void PrecMatrix (double Mjd_1, const MJDTIME* Mjd_2, double Mat[])
{
	double T, dT;
	double zeta,z,theta;
	double M1[9], M2[9], M3[9];

	T  = ( Mjd_1 - MJD_J2000 ) / 36525.0;
	dT = ( Mjd_2->Days -Mjd_1 + Mjd_2->FracDay ) / 36525.0;

	/* Precession angles */

	zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+
		((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/Arcs;
	z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/Arcs;
	theta =  ( (2004.3109-(0.85330+0.000217*T)*T)-
		((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/Arcs;


	Rotation_z( -z, M1 );
	Rotation_y( theta, M2 );
	MatrixMultiply( 3, 3, 3, 3, M1, M2, M3 );

	Rotation_z( -zeta, M1 );
	MatrixMultiply( 3, 3, 3, 3, M3, M1, Mat );

}


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

void NutAngles (const MJDTIME* Mjd_TT, double* dpsi, double* deps)
{
	int i;
	double T, T2, T3, rev;

	double  l, lp, F, D, Om;
	double  arg;

	static int C[106][9] = {

		/* l  l' F  D Om    dpsi    *T     deps     *T     */

		{  0, 0, 0, 0, 1,-1719960,-1742,  920250,   89 },   //   1
		{  0, 0, 0, 0, 2,   20620,    2,   -8950,    5 },   //   2
		{ -2, 0, 2, 0, 1,     460,    0,    -240,    0 },   //   3
		{  2, 0,-2, 0, 0,     110,    0,       0,    0 },   //   4
		{ -2, 0, 2, 0, 2,     -30,    0,      10,    0 },   //   5
		{  1,-1, 0,-1, 0,     -30,    0,       0,    0 },   //   6
		{  0,-2, 2,-2, 1,     -20,    0,      10,    0 },   //   7
		{  2, 0,-2, 0, 1,      10,    0,       0,    0 },   //   8
		{  0, 0, 2,-2, 2, -131870,  -16,   57360,  -31 },   //   9
		{  0, 1, 0, 0, 0,   14260,  -34,     540,   -1 },   //  10
		{  0, 1, 2,-2, 2,   -5170,   12,    2240,   -6 },   //  11
		{  0,-1, 2,-2, 2,    2170,   -5,    -950,    3 },   //  12
		{  0, 0, 2,-2, 1,    1290,    1,    -700,    0 },   //  13
		{  2, 0, 0,-2, 0,     480,    0,      10,    0 },   //  14
		{  0, 0, 2,-2, 0,    -220,    0,       0,    0 },   //  15
		{  0, 2, 0, 0, 0,     170,   -1,       0,    0 },   //  16
		{  0, 1, 0, 0, 1,    -150,    0,      90,    0 },   //  17
		{  0, 2, 2,-2, 2,    -160,    1,      70,    0 },   //  18
		{  0,-1, 0, 0, 1,    -120,    0,      60,    0 },   //  19
		{ -2, 0, 0, 2, 1,     -60,    0,      30,    0 },   //  20
		{  0,-1, 2,-2, 1,     -50,    0,      30,    0 },   //  21
		{  2, 0, 0,-2, 1,      40,    0,     -20,    0 },   //  22
		{  0, 1, 2,-2, 1,      40,    0,     -20,    0 },   //  23
		{  1, 0, 0,-1, 0,     -40,    0,       0,    0 },   //  24
		{  2, 1, 0,-2, 0,      10,    0,       0,    0 },   //  25
		{  0, 0,-2, 2, 1,      10,    0,       0,    0 },   //  26
		{  0, 1,-2, 2, 0,     -10,    0,       0,    0 },   //  27
		{  0, 1, 0, 0, 2,      10,    0,       0,    0 },   //  28
		{ -1, 0, 0, 1, 1,      10,    0,       0,    0 },   //  29
		{  0, 1, 2,-2, 0,     -10,    0,       0,    0 },   //  30
		{  0, 0, 2, 0, 2,  -22740,   -2,    9770,   -5 },   //  31
		{  1, 0, 0, 0, 0,    7120,    1,     -70,    0 },   //  32
		{  0, 0, 2, 0, 1,   -3860,   -4,    2000,    0 },   //  33
		{  1, 0, 2, 0, 2,   -3010,    0,    1290,   -1 },   //  34
		{  1, 0, 0,-2, 0,   -1580,    0,     -10,    0 },   //  35
		{ -1, 0, 2, 0, 2,    1230,    0,    -530,    0 },   //  36
		{  0, 0, 0, 2, 0,     630,    0,     -20,    0 },   //  37
		{  1, 0, 0, 0, 1,     630,    1,    -330,    0 },   //  38
		{ -1, 0, 0, 0, 1,    -580,   -1,     320,    0 },   //  39
		{ -1, 0, 2, 2, 2,    -590,    0,     260,    0 },   //  40
		{  1, 0, 2, 0, 1,    -510,    0,     270,    0 },   //  41
		{  0, 0, 2, 2, 2,    -380,    0,     160,    0 },   //  42
		{  2, 0, 0, 0, 0,     290,    0,     -10,    0 },   //  43
		{  1, 0, 2,-2, 2,     290,    0,    -120,    0 },   //  44
		{  2, 0, 2, 0, 2,    -310,    0,     130,    0 },   //  45
		{  0, 0, 2, 0, 0,     260,    0,     -10,    0 },   //  46
		{ -1, 0, 2, 0, 1,     210,    0,    -100,    0 },   //  47
		{ -1, 0, 0, 2, 1,     160,    0,     -80,    0 },   //  48
		{  1, 0, 0,-2, 1,    -130,    0,      70,    0 },   //  49
		{ -1, 0, 2, 2, 1,    -100,    0,      50,    0 },   //  50
		{  1, 1, 0,-2, 0,     -70,    0,       0,    0 },   //  51
		{  0, 1, 2, 0, 2,      70,    0,     -30,    0 },   //  52
		{  0,-1, 2, 0, 2,     -70,    0,      30,    0 },   //  53
		{  1, 0, 2, 2, 2,     -80,    0,      30,    0 },   //  54
		{  1, 0, 0, 2, 0,      60,    0,       0,    0 },   //  55
		{  2, 0, 2,-2, 2,      60,    0,     -30,    0 },   //  56
		{  0, 0, 0, 2, 1,     -60,    0,      30,    0 },   //  57
		{  0, 0, 2, 2, 1,     -70,    0,      30,    0 },   //  58
		{  1, 0, 2,-2, 1,      60,    0,     -30,    0 },   //  59
		{  0, 0, 0,-2, 1,     -50,    0,      30,    0 },   //  60
		{  1,-1, 0, 0, 0,      50,    0,       0,    0 },   //  61
		{  2, 0, 2, 0, 1,     -50,    0,      30,    0 },   //  62
		{  0, 1, 0,-2, 0,     -40,    0,       0,    0 },   //  63
		{  1, 0,-2, 0, 0,      40,    0,       0,    0 },   //  64
		{  0, 0, 0, 1, 0,     -40,    0,       0,    0 },   //  65
		{  1, 1, 0, 0, 0,     -30,    0,       0,    0 },   //  66
		{  1, 0, 2, 0, 0,      30,    0,       0,    0 },   //  67
		{  1,-1, 2, 0, 2,     -30,    0,      10,    0 },   //  68
		{ -1,-1, 2, 2, 2,     -30,    0,      10,    0 },   //  69
		{ -2, 0, 0, 0, 1,     -20,    0,      10,    0 },   //  70
		{  3, 0, 2, 0, 2,     -30,    0,      10,    0 },   //  71
		{  0,-1, 2, 2, 2,     -30,    0,      10,    0 },   //  72
		{  1, 1, 2, 0, 2,      20,    0,     -10,    0 },   //  73
		{ -1, 0, 2,-2, 1,     -20,    0,      10,    0 },   //  74
		{  2, 0, 0, 0, 1,      20,    0,     -10,    0 },   //  75
		{  1, 0, 0, 0, 2,     -20,    0,      10,    0 },   //  76
		{  3, 0, 0, 0, 0,      20,    0,       0,    0 },   //  77
		{  0, 0, 2, 1, 2,      20,    0,     -10,    0 },   //  78
		{ -1, 0, 0, 0, 2,      10,    0,     -10,    0 },   //  79
		{  1, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  80
		{ -2, 0, 2, 2, 2,      10,    0,     -10,    0 },   //  81
		{ -1, 0, 2, 4, 2,     -20,    0,      10,    0 },   //  82
		{  2, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  83
		{  1, 1, 2,-2, 2,      10,    0,     -10,    0 },   //  84
		{  1, 0, 2, 2, 1,     -10,    0,      10,    0 },   //  85
		{ -2, 0, 2, 4, 2,     -10,    0,      10,    0 },   //  86
		{ -1, 0, 4, 0, 2,      10,    0,       0,    0 },   //  87
		{  1,-1, 0,-2, 0,      10,    0,       0,    0 },   //  88
		{  2, 0, 2,-2, 1,      10,    0,     -10,    0 },   //  89
		{  2, 0, 2, 2, 2,     -10,    0,       0,    0 },   //  90
		{  1, 0, 0, 2, 1,     -10,    0,       0,    0 },   //  91
		{  0, 0, 4,-2, 2,      10,    0,       0,    0 },   //  92
		{  3, 0, 2,-2, 2,      10,    0,       0,    0 },   //  93
		{  1, 0, 2,-2, 0,     -10,    0,       0,    0 },   //  94
		{  0, 1, 2, 0, 1,      10,    0,       0,    0 },   //  95
		{ -1,-1, 0, 2, 1,      10,    0,       0,    0 },   //  96
		{  0, 0,-2, 0, 1,     -10,    0,       0,    0 },   //  97
		{  0, 0, 2,-1, 2,     -10,    0,       0,    0 },   //  98
		{  0, 1, 0, 2, 0,     -10,    0,       0,    0 },   //  99
		{  1, 0,-2,-2, 0,     -10,    0,       0,    0 },   // 100
		{  0,-1, 2, 0, 1,     -10,    0,       0,    0 },   // 101
		{  1, 1, 0,-2, 1,     -10,    0,       0,    0 },   // 102
		{  1, 0,-2, 2, 0,     -10,    0,       0,    0 },   // 103
		{  2, 0, 0, 2, 0,      10,    0,       0,    0 },   // 104
		{  0, 0, 2, 4, 2,     -10,    0,       0,    0 },   // 105
		{  0, 1, 0, 1, 0,      10,    0,       0,    0 }    // 106
	};

	T  = ( Mjd_TT->Days + Mjd_TT->FracDay -MJD_J2000 ) / 36525.0;
	T2 = T  * T;
	T3 = T2 * T;
	rev = 360.0*3600.0;  

	/* Mean arguments of luni-solar motion

	l   mean anomaly of the Moon
	l'  mean anomaly of the Sun
	F   mean argument of latitude
	D   mean longitude elongation of the Moon from the Sun 
	Om  mean longitude of the ascending node
	*/

	l  = 485866.733 + (1325.0*rev +  715922.633)*T + 31.310*T2 + 0.064*T3;
	lp =1287099.804 + (  99.0*rev + 1292581.224)*T -  0.577*T2 - 0.012*T3;
	F  = 335778.877 + (1342.0*rev +  295263.137)*T - 13.257*T2 + 0.011*T3;
	D  = 1072261.307 + (1236.0*rev + 1105601.328)*T - 6.891*T2 + 0.019*T3;
	Om =  450160.280 - (   5.0*rev +  482890.539)*T + 7.455*T2 + 0.008*T3;

	/* Nutation in longitude and obliquity [rad] */

	*deps = *dpsi = 0.0;

	for ( i=0; i<106; i++ )
	{
		arg  =  ( C[i][0]*l+C[i][1]*lp+C[i][2]*F+C[i][3]*D+C[i][4]*Om ) / Arcs;
		*dpsi += ( C[i][5]+C[i][6]*T ) * sin(arg);
		*deps += ( C[i][7]+C[i][8]*T ) * cos(arg);
	}

	*dpsi = (1.0E-5 * *dpsi)/Arcs;
	*deps = (1.0E-5 * *deps)/Arcs;

}

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

void NutMatrix ( const MJDTIME* Mjd_TT, double Mat[])
{
	double dpsi, deps, eps;
	double M1[9], M2[9], M3[9];

	eps = MeanObliquity(Mjd_TT);  /* Mean obliquity of the ecliptic */

	NutAngles ( Mjd_TT, &dpsi, &deps ); /*Nutation in longitude and obliquity*/

	Rotation_x( -eps-deps, M1 );
	Rotation_z( -dpsi, M2 );
	MatrixMultiply( 3, 3, 3, 3, M1, M2, M3 );

	Rotation_x( eps, M1 );
	MatrixMultiply( 3, 3, 3, 3, M3, M1, Mat );    

}

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

void NutMatrixSimple ( const MJDTIME* Mjd_TT, double Mat[])
{
	double  T, ls, D, F, N;
	double  eps, dpsi, deps;
	double M1[9], M2[9], M3[9];

	T  = (Mjd_TT->Days +Mjd_TT->FracDay - MJD_J2000 ) / 36525.0;

	// Mean arguments of luni-solar motion

	ls = pi2 * fmod( 0.993133+  99.997306*T, 1.0 );   // mean anomaly Sun          
	D  = pi2 * fmod( 0.827362+1236.853087*T, 1.0 );   // diff. longitude Moon-Sun  
	F  = pi2 * fmod( 0.259089+1342.227826*T, 1.0 );   // mean argument of latitude 
	N  = pi2 * fmod( 0.347346-   5.372447*T, 1.0 );   // longit. ascending node    

	// Nutation angles

	dpsi = ( -17.200*sin(N)   - 1.319*sin(2*(F-D+N)) - 0.227*sin(2*(F+N))
		+ 0.206*sin(2*N) + 0.143*sin(ls) ) / Arcs;
	deps = ( + 9.203*cos(N)   + 0.574*cos(2*(F-D+N)) + 0.098*cos(2*(F+N))
		- 0.090*cos(2*N)                 ) / Arcs;

	// Mean obliquity of the ecliptic

	eps  = 0.4090928-2.2696E-4*T;   

	Rotation_x( -eps-deps, M1 );
	Rotation_z( -dpsi, M2 );
	MatrixMultiply( 3, 3, 3, 3, M1, M2, M3 );

	Rotation_x( eps, M1 );
	MatrixMultiply( 3, 3, 3, 3, M3, M1, Mat );    

}


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

double EqnEquinox ( const MJDTIME* Mjd_TT )
{
	double dpsi, deps;              /* Nutation angles */
	double T, T2, T3, rev;
	double Om, Angle;

	T  = ( Mjd_TT->Days + Mjd_TT->FracDay -MJD_J2000 ) / 36525.0;
	T2 = T*T;
	T3 = T2*T;
	rev = 360.0*3600.0;      /* arcsec/revolution */

	Om = ( 450160.280 - ( 5.0*rev +  482890.539)*T +  
		7.455*T2 + 0.008*T3 ) / Arcs;

	NutAngles (Mjd_TT, &dpsi, &deps );

	Angle=(0.00264 * sin(Om) + 0.000063 * sin(2*Om)) / Arcs ;

	return  dpsi * cos ( MeanObliquity(Mjd_TT) )+Angle ;

}


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

double GMST ( const MJDTIME* Mjd_UT1 )
{
	double T, T0, Res;
	T  = ( Mjd_UT1->Days - MJD_J2000 + Mjd_UT1->FracDay ) / 36525.0; 

	T0 = ( Mjd_UT1->Days - MJD_J2000 ) + Mjd_UT1->FracDay;

	Res = Mjd_UT1->FracDay + 0.2790572732640 + 0.00273781191135448 * T0;

	Res = fmod( Res*pi2, pi2 );
	Res = Res+ (0.014506+(4612.15739966+(1.39667721+ 
		( -0.00009344+0.00001882* T)*T)*T)*T) / Arcs;

	Res = fmod( Res, pi2 );

	return Res;   
}

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

double GAST ( const MJDTIME* Mjd_UT1)
{
	double Res;

	Res = GMST(Mjd_UT1) + EqnEquinox(Mjd_UT1);

	Res = fmod( Res, pi2 );

	return Res;
}


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

void GHAMatrix ( const MJDTIME* Mjd_UT1, double Mat[] )
{
	double Angle;

	Angle = GAST( Mjd_UT1 );

	Rotation_z( Angle, Mat );

}


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
//   x, y      Pole Motion para[s]
//   <return>  Pole matrix
//
***************************************************************************/

void PoleMatrix ( double x, double y, double Mat[])
{
	double M1[9], M2[9];

	Rotation_y( -x/Arcs, M1 );
	Rotation_x( -y/Arcs, M2 );

	MatrixMultiply( 3, 3, 3, 3, M1, M2, Mat );

}

/***************************************************************************
//
// ICRF_ITRF_MJD
//
// 目的: 天球系与地固系坐标间的相互转换,其中包含位置和速度分量, 输入时间
为MJD表示法的TT时.
//
// 参数:
//
//   Mjd1     Epoch given (Modified Julian Date TT)默认是MJD_2000
//   Mjd2     Epoch to precess to (Modified Julian Date TT)
//   Eop[2]   EOP参数
//   flag     1: 从天球系到地固系 
0: 从地固系到天球系
//   ICRF     天球系坐标，其中包含位置和速度
//   ITRF     地固系坐标
//
***************************************************************************/
void ICRF_ITRF_MJD( const double Mjd1, const MJDTIME* Mjd2,
	int flag, double ICRF[6], double ITRF[6] )

{
	int i;
	MJDTIME Mjd_UT1,Mjd_UTC;
	double  U[9], U_dot[9], UT[9], U_dotT[9];      
	double Prec[9], Nut[9], GH[9], Pole[9];
	double GH_dot[9];
	double tmp1[9], tmp2[9];
	double v1[3], v2[3];

	EOPPARA CurrEop;

	Mjd_UTC.Days = Mjd2->Days;
	Mjd_UTC.FracDay = Mjd2->FracDay - TT_UTC() / SECPERDAY;
	InterposeEOP( &Mjd_UTC, &CurrEop );

	Mjd_UT1.Days = Mjd_UTC.Days;
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	PrecMatrix( Mjd1, Mjd2, Prec );
	NutMatrixSimple ( Mjd2, Nut );
	GHAMatrix ( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );

	for( i=0; i<3; i++ )
	{
		GH_dot[i]   =  Omega_WGS * GH[3+i];
		GH_dot[3+i] = -Omega_WGS * GH[i];
		GH_dot[6+i] = 0.0;
	}

	MatrixMultiply( 3, 3, 3, 3, Nut, Prec, tmp1 );    
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, tmp2 );
	MatrixMultiply( 3, 3, 3, 3, tmp2, tmp1, U );   /*坐标转换矩阵*/

	MatrixMultiply( 3, 3, 3, 3, Pole, GH_dot, tmp2 );
	MatrixMultiply( 3, 3, 3, 3, tmp2, tmp1, U_dot );  /*速度转换矩阵 */

	if( flag == 1 )     /*天球系——>地固系 */
	{
		MatrixMultiply( 3, 3, 3, 1, U, ICRF, ITRF );

		MatrixMultiply( 3, 3, 3, 1, U, &ICRF[3], v1 );
		MatrixMultiply( 3, 3, 3, 1, U_dot, ICRF, v2 );

		MatrixAddition( 3, 1, v1, v2, &ITRF[3] );

	}
	else  /*地固系——>天球系 */
	{  
		MatrixTranspose( 3, 3, U, UT );
		MatrixTranspose( 3, 3, U_dot, U_dotT );

		MatrixMultiply( 3, 3, 3, 1, UT, ITRF, ICRF );

		MatrixMultiply( 3, 3, 3, 1, UT, &ITRF[3], v1 );
		MatrixMultiply( 3, 3, 3, 1, U_dotT, ITRF, v2 );

		MatrixAddition( 3, 1, v1, v2, &ICRF[3] );
	}
}

/***************************************************************************
//
// ICRF_ITRF_GPST
//
// 目的: 天球系与地固系坐标间的相互转换,其中包含位置和速度分量, 输入时间
为MJD表示法的GPS时.
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
	int flag, double ICRF[6], double ITRF[6] )
{
	MJDTIME Mjd_TT;

	GPSTimeToMJDTime( GT, &Mjd_TT );    /* 此处Mjd_TT为GPS时间的MJD表示 */
	Mjd_TT.FracDay = Mjd_TT.FracDay - ( GPST_TAI - TT_TAI )/SECPERDAY;

	ICRF_ITRF_MJD( Mjd1, &Mjd_TT, flag, ICRF, ITRF );

}

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
void ICRF_ITRF_Matrix( const double Mjd1, const GPSTIME* GT,
	int flag, double Mat[9] )
{
	int i;
	MJDTIME Mjd_UT1, Mjd_UTC, Mjd_TT;
	double Prec[9], Nut[9], GH[9], Pole[9];
	double tmp1[9], tmp2[9];

	EOPPARA CurrEop;

	GPSTimeToMJDTime( GT, &Mjd_TT );    /* 此处Mjd_TT为GPS时间的MJD表示 */
	Mjd_TT.FracDay = Mjd_TT.FracDay - ( GPST_TAI - TT_TAI )/SECPERDAY;

	Mjd_UTC.Days = Mjd_TT.Days;
	Mjd_UTC.FracDay = Mjd_TT.FracDay - TT_UTC() / SECPERDAY;
	InterposeEOP( &Mjd_UTC, &CurrEop );

	Mjd_UT1.Days = Mjd_UTC.Days;
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	PrecMatrix( Mjd1, &Mjd_TT, Prec );
	NutMatrixSimple ( &Mjd_TT, Nut );
	GHAMatrix ( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );

	MatrixMultiply( 3, 3, 3, 3, Nut, Prec, tmp1 );    
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, tmp2 );
	MatrixMultiply( 3, 3, 3, 3, tmp2, tmp1, Mat );   /*坐标转换矩阵*/

	if( flag == 0 )     /*天球系——>地固系 */
	{  
		MatrixTranspose( 3, 3, Mat, tmp1 );

		for( i=0; i<9; i++ )
		{
			Mat[i] = tmp1[i];
		}
	}


}
