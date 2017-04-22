/****************************************************************************
目的：    定义时间结构体及其相互转换函数,
空间直角坐标和大地坐标相互转换
编写时间：2008.11.22
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#include <math.h>
#include "GPSTime.h"
#include "RTOD_Const.h"
#include "CommonFuncs.h"
//2000年1月1日12点（J2000）的儒略日数为2451545.0
//sec是从J2000开始算的秒数
double SecTimeToJD(double sec)
{
	double jd=0.0,frac;
	frac = fmod(sec,1.0);
	int isec;
	if(frac < 0.0) 
		frac+=1.0;
	isec = sec - frac + 0.5;
	jd = isec/86400.0 + 2451545.0;
	return jd;
}
//儒略日到简化儒略日
double JDToMJD(double jd)
{
	return jd - 2400000.5;
}
//sec到通用时
void SecTimeToCT(double sec,COMMONTIME *CT)
{
	MJDTIME MJDT;
	double jd=SecTimeToJD(sec);
	double MJD=JDToMJD(jd);
	MJDT.FracDay=fmod(MJD,1.0);
	MJDT.Days=MJD-MJDT.FracDay;
	MJDTimeToCommonTime(&MJDT,CT);
}
/****************************************************************************
CommonTimeToMJDTime
目的：将通用时转换为简化儒略日

参数：
CT   通用时
MJDT 简化儒略日
****************************************************************************/

void CommonTimeToMJDTime( const COMMONTIME* CT, MJDTIME* MJDT)
{
	int y, m, temp;    /* 临时变量    */

	y = CT->Year;

	if( CT->Year < 100 )     /* 使用两位年记法 */
	{
		if( CT->Year < 80 )   /* 只支持1980~2079  */
		{
			y = y + 2000;  
		}
		else   
		{
			y = y + 1900;
		}
	}

	if( CT->Month <= 2 )   
	{
		y = y - 1;
		m = CT->Month + 12;
	}
	else
	{
		m = CT->Month;
	}

	temp = ( int )( 365.25 * y );
	temp += ( int )( 30.6001 * ( m+1 ) );
	temp += CT->Day;
	temp += -679019;

	MJDT->Days = temp;

	MJDT->FracDay = CT->Hour + CT->Minute / 60.0 + CT->Second / SECPERHOUR;
	MJDT->FracDay = MJDT->FracDay / 24.0;

}

/****************************************************************************
MJDTimeToCommonTime
目的：将简化儒略日转换为通用时

参数：
MJDT 简化儒略日
CT   通用时
****************************************************************************/

void MJDTimeToCommonTime( const MJDTIME* MJDT, COMMONTIME* CT )
{
	int a,b,c,d,e;

	a = ( int )( MJDT->Days + MJDT->FracDay + 2400000.5 + 0.5 );
	b = a + 1537;
	c = ( int )( ( b-122.1)/365.25 );
	d = ( int )( 365.25 * c );
	e = ( int )( ( b - d ) / 30.6001 );

	CT->Day = b - d - ( int )( 30.6001 * e );
	CT->Month = e - 1 - 12 * ( int )( e / 14 );
	CT->Year = c - 4715 - ( int )(( 7 + CT->Month ) / 10 );

	CT->Hour = ( int )( MJDT->FracDay * 24 );
	CT->Minute = ( int )( ( MJDT->FracDay * 24 - CT->Hour ) * 60 );
	CT->Second = ( ( MJDT->FracDay * 24 - CT->Hour ) * 60 - CT->Minute ) * 60;
}

/****************************************************************************
GPSTimeToMJDTime
目的：将GPS时间转换为简化儒略日

参数：
GT   GPS时间
MJDT 简化儒略日
****************************************************************************/

void GPSTimeToMJDTime( const GPSTIME* GT, MJDTIME* MJDT )
{
	int day;

	day = ( int )( GT->SecOfWeek / SECPERDAY );
	MJDT->FracDay = GT->SecOfWeek / SECPERDAY - day;

	MJDT->Days = JAN61980 + GT->Week * 7 + day;

}

/****************************************************************************
MJDTimeToGPSTime
目的：将简化儒略日转换为GPS时间

参数：
GT   GPS时间
MJDT 简化儒略日
****************************************************************************/
void MJDTimeToGPSTime ( const MJDTIME* MJDT, GPSTIME* GT )
{
	int RemainDay;

	GT->Week = (int)(( MJDT->Days - JAN61980) / 7 ); 

	RemainDay = MJDT->Days - GT->Week * 7 - JAN61980;

	GT->SecOfWeek = ( RemainDay + MJDT->FracDay ) * SECPERDAY;

}

/****************************************************************************
CommonTimeToGPSTime
目的：将通用时转换为GPS时间

参数：
CT   通用时
GT   GPS时间
****************************************************************************/
void CommonTimeToGPSTime ( const COMMONTIME* CT, GPSTIME* GT )
{
	MJDTIME mjd;

	CommonTimeToMJDTime( CT, &mjd );
	MJDTimeToGPSTime( &mjd, GT );

}

/****************************************************************************
GPSTimeToCommonTime
目的：将GPS时间转换为通用时

参数：
GT   GPS时间
CT   通用时
****************************************************************************/
void GPSTimeToCommonTime ( const GPSTIME* GT, COMMONTIME* CT )
{
	MJDTIME mjd;

	GPSTimeToMJDTime( GT, &mjd );
	MJDTimeToCommonTime( &mjd, CT );

}
/****************************************************************************
GPSTimeToep
目的：将GPS时间转换为通用时ep

参数：
GT   GPS时间
ep   rtklib下的通用时
****************************************************************************/
void GPSTimeToep ( const GPSTIME* GT, double ep[6])
{
	COMMONTIME CT;
	MJDTIME mjd;

	GPSTimeToMJDTime( GT, &mjd );
	MJDTimeToCommonTime( &mjd, &CT );
	ep[0] = (double)CT.Year;
	ep[1] = (double)CT.Month;
	ep[2] = (double)CT.Day;
	ep[3] = (double)CT.Hour;
	ep[4] = (double)CT.Minute;
	ep[5] = (double)CT.Second;

}

/****************************************************************************
CheckGPSTime
目的：检查GPS时间格式是否正确, 如周秒是否大于604800或小于0 

参数：
GT   GPS时间
****************************************************************************/
void CheckGPSTime( GPSTIME* GT )
{
	if( GT->SecOfWeek<0.0 )
	{
		GT->Week = GT->Week - 1;
		GT->SecOfWeek = GT->SecOfWeek + SECPERWEEK;
	}
	else if( GT->SecOfWeek>SECPERWEEK )
	{
		GT->Week = GT->Week + 1;
		GT->SecOfWeek = GT->SecOfWeek - SECPERWEEK;
	}

}

/****************************************************************************
XYZToBLH
目的：将空间直角坐标转换为大地坐标

参数：
XYZ   空间直角坐标[m]
BLH   大地坐标[Rad, m]
R     参考椭球的长半径[m]
F     扁率(1/f)
****************************************************************************/
void XYZToBLH( const double XYZ[], double BLH[], double R, double F )
{
	short   Iterator;
	double  e2, dZ, rho2, dZ_new, SinPhi;
	double  ZdZ, Nh, N;

	e2      = F * ( 2.0 - F );
	rho2 = XYZ[0] * XYZ[0] + XYZ[1] * XYZ[1]; 
	dZ = e2 * XYZ[2];

	Iterator = 0;
	for(;;)
	{
		ZdZ    =  XYZ[2] + dZ;
		Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 

		if( Nh < 1.0 )  /* 如果XYZ坐标为0.0的情况 */
		{
			BLH[0] = 0.0;
			BLH[1] = 0.0;
			BLH[2] = 0.0;
			return;
		}

		SinPhi =  ZdZ / Nh;                    
		N      =  R / sqrt(1.0 - e2*SinPhi*SinPhi); 
		dZ_new =  N * e2 * SinPhi;

		Iterator = Iterator + 1;
		if( (fabs(dZ-dZ_new) < 1E-10) || (Iterator >=10) )
		{
			break;
		}

		dZ = dZ_new;
	}

	BLH[1] = atan2 ( XYZ[1], XYZ[0] );
	BLH[0] = atan2 ( ZdZ, sqrt(rho2) );
	BLH[2] = Nh - N;
}

/****************************************************************************
BLHToXYZ
目的：将大地坐标转换为空间直角坐标

参数：
BLH   大地坐标[Rad, m]
XYZ   空间直角坐标[m]
R     参考椭球的长半径[m]
F     扁率(1/f)
****************************************************************************/
void BLHToXYZ( const double BLH[], double XYZ[], double R, double F )
{
	double  e2, CosLat, SinLat;
	double  N;

	e2      = F * ( 2.0 - F );
	CosLat = cos( BLH[0] );        
	SinLat = sin( BLH[0] );

	N = R_Earth / sqrt(1.0-e2*SinLat*SinLat);

	XYZ[0] =  (         N+BLH[2])*CosLat*cos(BLH[1]);
	XYZ[1] =  (         N+BLH[2])*CosLat*sin(BLH[1]);
	XYZ[2] =  ((1.0-e2)*N+BLH[2])*SinLat;
}

/****************************************************************************
CoorTranFromPZ90ToWGS84

目的：PZ90 与WGS84坐标转换
方法与参数来自陈俊勇<GPS和GLONASS定位成果的坐标转换>
测绘通报,2002年第七期

X(WGS84) = X(PZ90) + dX0 + R*X(PZ90)
|  D   -R3    R2 |
R = |  R3   D    -R1 |
| -R2   R1    D  |
dX0 = [ 7cm, 0, -77cm ]
D = -3E-9,  R1 = -19 mas, R2 = -4mas,  R3 = 353mas

MCC转换参数:
dX0 = [-0.47m, -0.51m, -1.56m]
D = 22E-9,   R3=1.728E-6, R2=-0.017E-6, R1 = -0.076E-6

MIT提供的转换参数:
dX0 = [0.0m, 2.5m, 0.0m]
D = 0,   R3=1.9E-6, R2=0.0, R1 = 0.0

2007.9.20后, 转换参数为
dX0 = [-0.36m, 0.08m, 0.180m]


参数：
BLH   大地坐标[Rad, m]
XYZ   空间直角坐标[m]
R     参考椭球的长半径[m]
F     扁率(1/f)
****************************************************************************/
/*  */

void CoorTranFromPZ90ToWGS84( double PZ90[3], double WGS84[3] )
{
	/*  double X[3], dX0[3] = { 0.07, 0.0, -0.77 };
	double R[9] = { 1.0-3E-9, -353/Arcs/1000.0, -4/Arcs/1000.0,
	353/Arcs/1000.0,  1.0-3E-9, 19/Arcs/1000.0,
	4/Arcs/1000.0,  -19/Arcs/1000.0, 1.0-3E-9 };
	*/ 
	/*    double X[3], dX0[3] = { -0.47, -0.51, -1.560 };
	double R[9] = { 1.0-3E-9,    -1.728E-6, -0.017E-6,
	1.728E-6,  1.0-3E-9,        0.076E-6,
	0.017E-6,  -0.076E-6,    1.0-3E-9 };
	*/   
	double X[3], dX0[3] = { 0.0, 2.5, 0.0 };
	double R[9] = { 1.0, -1.9E-6, 0.0,
		1.9E-6,  1.0,    0.0,
		0.0,  0.0,    1.0 };

	MatrixMultiply( 3, 3, 3, 1, R, PZ90, X );

	MatrixAddition( 3, 1, X, dX0, WGS84 );
}


/****************************************************************************
VeloTranFromPZ90ToWGS84

目的：PZ90 与WGS84速度转换
方法与参数来自陈俊勇<GPS和GLONASS定位成果的坐标转换>
测绘通报,2002年第七期

X(WGS84) = X(PZ90) + dX0 + R*X(PZ90)
|  D   -R3    R2 |
R = |  R3   D    -R1 |
| -R2   R1    D  |
dX0 = [ 7cm, 0, -77cm ]
D = -3E-9,  R1 = -19 mas, R2 = -4mas,  R3 = 353mas

MCC转换参数:
dX0 = [-0.47m, -0.51m, -1.56m]
D = 22E-9,   R3=1.728E-6, R2=-0.017E-6, R1 = -0.076E-6

MIT提供的转换参数:
dX0 = [0.0m, 2.5m, 0.0m]
D = 0,   R3=1.9E-6, R2=0.0, R1 = 0.0


参数：
BLH   大地坐标[Rad, m]
XYZ   空间直角坐标[m]
R     参考椭球的长半径[m]
F     扁率(1/f)
****************************************************************************/
/*  */

void VeloTranFromPZ90ToWGS84( double PZ90[3], double WGS84[3] )
{
	/*   double R[9] = { 1.0-3E-9,    -1.728E-6, -0.017E-6,
	1.728E-6,  1.0-3E-9,     0.076E-6,
	0.017E-6,  -0.076E-6,    1.0-3E-9 };
	double R[9] = { 1.0-3E-9, -353/Arcs/1000.0, -4/Arcs/1000.0,
	353/Arcs/1000.0,  1.0-3E-9, 19/Arcs/1000.0,
	4/Arcs/1000.0,  -19/Arcs/1000.0, 1.0-3E-9 };
	*/ 
	double R[9] = { 1.0, -1.9E-6, 0.0,
		1.9E-6,  1.0,    0.0,
		0.0,  0.0,    1.0 };

	MatrixMultiply( 3, 3, 3, 1, R, PZ90, WGS84 );

}


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

H        转换到NEU测站地平坐标的旋转矩阵

***************************************************************************/

void BLHToNEUMatrix( double BLH[3], double H[9])
{
	double sinB, cosB;
	double sinL, cosL;

	sinB = sin( BLH[0] );
	cosB = cos( BLH[0] );
	sinL = sin( BLH[1] );
	cosL = cos( BLH[1] );

	H[0]= -sinB * cosL;
	H[1]= -sinB * sinL;
	H[2]=  cosB;

	H[3]= -sinL;
	H[4]=  cosL;
	H[5]=  0.0;

	H[6]=  cosB * cosL;
	H[7]=  cosB * sinL;
	H[8]=  sinB;
}


/***************************************************************************
//
// PhaseCentToMassCent
//
// 目的: 将GNSS接收机天线相位中心坐标转换到质量中心

说明: 星固系定义(以CHAMP卫星为例)

原点: 卫星的质量中心
X轴:  指向卫星飞行方向
Z轴:  指向地球质量中心
Y轴:  与Z和X组成右手坐标系

//
// 输入参数:
//
Flag      1为相位中心转换到质量中心, 0为质量中心转换相位中心
Bias      星固系中的GNSS接收机天线相位中心与质量中心的偏差参数[m]

输出参数

State[6]  输入为星载接收机的位置和速度, 输出为经过改正后的位置

***************************************************************************/

void PhaseCentToMassCent( bool Flag, const double Bias[3], double State[6] )
{
	int    i;
	double ModR;             /* 卫星的地心矢径 */
	double ModVel;           /* 卫星速度的模   */
	double X[3] = { 0.0 }, Y[3] = { 0.0 }, Z[3] = { 0.0 };

	ModR   = sqrt( VectDot( 3, 3, State, State ) );
	ModVel = sqrt( VectDot( 3, 3, &State[3], &State[3] ) );

	if( (ModR <= (R_WGS84-1000.))
		&& (ModVel < 1000.0) )    /* 卫星位置不正确, 不改正 */
	{
		return; 
	}

	for( i=0; i<3; i++ )
	{
		Z[i] = State[i] / ( -1.0*ModR );
		X[i] = State[i+3] / ModVel;
	}

	CrossDot( 3, 3, Z, X, Y );

	if( Flag )   /* 相位中心到质量中心 */
	{
		for( i=0; i<3; i++ )
		{
			State[i] = State[i] - ( X[i]*Bias[0] + Y[i]*Bias[1] + Z[i]*Bias[2] );
		}
	}
	else
	{
		for( i=0; i<3; i++ )
		{
			State[i] = State[i] + ( X[i]*Bias[0] + Y[i]*Bias[1] + Z[i]*Bias[2] );
		}

	}
}

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
dXYZ      空间直角坐标系下的坐标分量
dRTN      RTN坐标系下的坐标分量

***************************************************************************/

void XYZToRTN( double State[6], double dXYZ[3], double dRTN[3] )   
{
	int i;
	double X[3],Y[3],Z[3], dis;
	for(i=0;i<3;i++)
		dRTN[i]=0.0;
	dis = sqrt( VectDot( 3, 3, State, State ) );
	X[0] = State[0]/dis;
	X[1] = State[1]/dis;
	X[2] = State[2]/dis;

	CrossDot( 3, 3, State, &State[3], Z );

	dis = sqrt( VectDot( 3, 3, Z, Z ) );
	Z[0] = Z[0]/dis;
	Z[1] = Z[1]/dis;
	Z[2] = Z[2]/dis;

	CrossDot( 3, 3, Z, X, Y );

	for( i=0; i<3; i++ )
	{
		dRTN[i] = X[i]*dXYZ[0]+Y[i]*dXYZ[1]+Z[i]*dXYZ[2];//原程序有误
		/*
		dRTN[0] += X[i]*dXYZ[i];//rj20160721
		dRTN[1] += Y[i]*dXYZ[i];
		dRTN[2] += Z[i]*dXYZ[i];
		*/
	}


}
//生成xyz转rtn的转换矩阵rj20160721
//flat=1  XYZtoRTN
//flat=0  RTNtoXYZ
void MatrixXYZ2RTN(double State[6],double transMat[9],int flag)
{
	int i;
	double X[3],Y[3],Z[3], dis;

	dis = sqrt( VectDot( 3, 3, State, State ) );
	X[0] = State[0]/dis;
	X[1] = State[1]/dis;
	X[2] = State[2]/dis;

	CrossDot( 3, 3, State, &State[3], Z );

	dis = sqrt( VectDot( 3, 3, Z, Z ) );
	Z[0] = Z[0]/dis;
	Z[1] = Z[1]/dis;
	Z[2] = Z[2]/dis;
	CrossDot( 3, 3, Z, X, Y );
	if(flag==1)//1：XYZ2RTN
		for( i=0; i<3; i++ )
		{
			transMat[i]=X[i];
			transMat[i+3]=Y[i];
			transMat[i+6]=Z[i];
		}
	else//0：RTN2XYZ
		for( i=0; i<3; i++ )
		{
			transMat[i*3]=X[i];
			transMat[i*3+1]=Y[i];
			transMat[i*3+2]=Z[i];
		}
}