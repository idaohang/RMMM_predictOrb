/****************************************************************************
目的：    单点实时定轨, 结合动力学模型滤波实时定轨


版本:     V1.1
版权：    武汉大学
****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include "RTOD_Const.h"
#include "CommonFuncs.h"
#include "RefSys.h"
#include "RTOD.h"

ONEPOCHOFSP3     *PRE_EPH;
pcvs_t pcvs;

void initKFState(EKFSTATE  *KFState)
{

	KFState->Step = 120.0;/* 积分步长 */
	//对于grace-A卫星
	KFState->Para[0].m_a = 50;//自主为10x10
	KFState->Para[0].n_a = 50;
	KFState->Para[0].GM  = GM_Earth;
	KFState->Para[0].R_ref = R_Earth;

	KFState->CentBias[0] =  0.0;
	KFState->CentBias[1] =  0.0;
	KFState->CentBias[2] = -0.444;

	KFState->Para[0].Mass = 487.2;
	KFState->Para[0].Cd   = 0.5;
	KFState->Para[0].Cr   = 1.1;
	KFState->Para[0].Area_H = 0.955;
	KFState->Para[0].Area_R = 5;

	//对于grace-B卫星
	
	KFState->Para[1].m_a = 50;//自主为10x10
	KFState->Para[1].n_a = 50;
	KFState->Para[1].GM  = GM_Earth;
	KFState->Para[1].R_ref = R_Earth;

	KFState->CentBias[3] =  0.0;
	KFState->CentBias[4] =  0.0;
	KFState->CentBias[5] = -0.444;

	KFState->Para[1].Mass = 487.2;
	KFState->Para[1].Cd   = 0.5;
	KFState->Para[1].Cr   = 1.1;
	KFState->Para[1].Area_H = 0.955;
	KFState->Para[1].Area_R = 5;

	KFState->Tao[0] = 1;KFState->Tao[1] = 1;       /* 补偿加速度的相关时间[s]  */
	KFState->Sigma[0] = 5E-6;KFState->Sigma[1] = 5E-6;  /*A、 B卫星补偿加速度的功率谱密度[m/s^2] */
	//KFState->Sigma[1] = 5E-7;
	//其余的在EmptyEKFSTATEStruct与EmptyAPRIORBITStruct中初始化
	KFState->eleMask=0;

	EmptyAPRIORBITStruct( &KFState->AprioriState[0]);
	EmptyAPRIORBITStruct( &KFState->AprioriState[1]);
	EmptyEKFSTATEStruct( KFState );
	EmptyCurComObsStruct(&KFState->CurComObs);
	EmptyddObsStruct(&KFState->ddOBS);

	KFState->wholeComObs=(Common11Obs *)malloc(MAXOBSNUMDAY*sizeof(Common11Obs));
}
void freeKFState(EKFSTATE KFState)
{
	free(KFState.wholeComObs);
}
/***************************************************************************
//
// EmptyAPRIORBITStruct
//
// 目的: 初始化APRIORBIT结构体

输入参数

ApriOrbit      待初始化的变量

***************************************************************************/

void EmptyAPRIORBITStruct( APRIORBIT* ApriOrbit )
{
	int i;
	ApriOrbit->LeoClk[0] = 0.0;
	ApriOrbit->LeoClk[1] = 0.0;

	ApriOrbit->OrbitAccuracy = 9999.999;
	ApriOrbit->Validate  = 0;

	for( i=0; i<3; i++ )
	{
		ApriOrbit->LeoPos[i] = 0.0;
	}

}

/***************************************************************************
//
// EmptySATLISTStruct
//
// 目的: 初始化SATLIST结构体

输入参数

SatList      待初始化的变量

***************************************************************************/

void EmptySATLISTStruct( SATLIST* SatList )
{
	int i;

	for( i=0; i<MAXCHANNUM; i++ )
	{
		SatList[i].Prn = 0;
		SatList[i].Status = 0;
		SatList[i].System = UNKS;
	}
}

/***************************************************************************
//
// EmptyPPRESULTStruct
//
// 目的: 初始化PPRESULT结构体

输入参数

PPResult      待初始化的变量

***************************************************************************/

void EmptyPPRESULTStruct( PPRESULT* PPResult )
{
	int i;

	for( i=0; i<3; i++ )
	{
		PPResult->Position[i] = 0.0;
		PPResult->PosBLH[i]   = 0.0;
		PPResult->Velocity[i] = 0.0;
		PPResult->VelNEU[i]   = 0.0;
		PPResult->MatrixFromBLHtoNEU[3*i+0] = 0.0;
		PPResult->MatrixFromBLHtoNEU[3*i+1] = 0.0;
		PPResult->MatrixFromBLHtoNEU[3*i+2] = 0.0;
	}

	PPResult->RcvClkOft[0] = 0.0;
	PPResult->RcvClkOft[1] = 0.0;
	PPResult->RcvClkSft = 0.0;
	PPResult->PDOP      = 999.0;
	PPResult->GDOP      = 999.0;
	PPResult->HDOP      = 999.0;
	PPResult->VDOP      = 999.0;
	PPResult->TDOP      = 999.0;
	PPResult->SigmaPos  = 999.0;
	PPResult->SigmaVel  = 999.0;
	PPResult->Coverage  = 999.0;
	PPResult->Iterator  = 20;
	PPResult->IsSuccess = false;

	for( i=0; i<MAXCHANNUM; i++ )
	{
		PPResult->Residual[i] = 0.0;

		PPResult->SatList[i].Prn = 0;
		PPResult->SatList[i].Status = 0;
		PPResult->SatList[i].System = UNKS;
	}
}                                  

/***************************************************************************
//
// EmptySATMIDRESULTStruct
//
// 目的: 初始化SATMIDRESULT结构体

输入参数

Num             SATMIDRESULT数组的维数
SatMidInfo      待初始化的变量

***************************************************************************/

void EmptySATMIDRESULTStruct( int Num, SATMIDRESULT* SatMidInfo )
{
	int i, j;

	for( j=0; j<Num; j++ )
	{
		for( i=0; i<3; i++ )
		{
			SatMidInfo[j].SatPos[i] = 0.0;
			SatMidInfo[j].SatVel[i] = 0.0;
		}

		SatMidInfo[j].Elevation = pi/2.0;
		SatMidInfo[j].IonoCorr  = 0.0;
		SatMidInfo[j].SatClkOft = 0.0;
		SatMidInfo[j].SatClkSft = 0.0;
		SatMidInfo[j].Relativity = 0.0;
		SatMidInfo[j].TropCorr   = 0.0;

	}

}
/*通过精密星历计算卫星位置和速度*/
void computeGPSfromSP3(const GPSTIME T_Tm,const int Prn,const int order,
	double * Pos0,double * Vel0,double * ClkOft,double * ClkSft)
{
	int k;
	ONESATEHP OneSatEph;//精密星历
	//精密钟差，为了算精确模糊度用

	//用于相位中心改正
	gtime_t time;//用于相位中心改正
	double ep[6];
	double rs[6];
	double dant[3];
	pcv_t *pcvt;
	pcvt=(pcv_t*)malloc(32*sizeof(pcv_t));

	/*精密星历*/
	Lagrange_SP3(T_Tm,Prn,order,PRE_EPH,&OneSatEph);
	//flag=Lagrange_CLK(T,(*OBS).obs[i].PRN,2,PRECLK,&OneSatClk);
	if (!OneSatEph.ISExist)
	{
		return;//如果精密星历内插失败，则返回,不对广播星历计算出的pos和vel进行赋值，即用广播星历计算得到的结果
	}
	Pos0[0]  =OneSatEph.satposx;
	Pos0[1]  =OneSatEph.satposy;
	Pos0[2]  =OneSatEph.satposz;
	Vel0[0]	 =OneSatEph.satvelx;
	Vel0[1]  =OneSatEph.satvely;
	Vel0[2]  =OneSatEph.satvelz;
	*ClkOft  =OneSatEph.satclk;
	//精密星历钟差的相对论改正
	*ClkOft   = *ClkOft - 2.0*VectDot(3,3,Pos0,Vel0)/C_Light/C_Light;
	*ClkSft   = 0.0;//避免后面钟速对广播星历钟差的补偿
	/*精密星历相位中心改正*/

	GPSTimeToep(&T_Tm,ep);
	time=epoch2time(ep);
	setpcv(time, pcvt, &pcvs);
	for(k=0;k<3;k++)
	{
		rs[k] = Pos0[k];
		rs[k+3] = Vel0[k];
	}
	/*rs就是地固系下卫星速度和坐标*/
	satantoff(time, rs,&pcvt[Prn-1],dant);
	/*精密星历下相位中心改正*/
	for(k=0;k<3;k++)
	{
		Pos0[k] += dant[k];
	}
	free(pcvt);
}
/***************************************************************************
//
// ComputeGPSSatOrbitAtSignalTrans
//
// 目的: 计算卫星发射时刻的卫星轨道,卫星钟差,相对论效应,卫星高度角
考虑地球自转改正(顾及卫星健康标记).

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
	SATMIDRESULT* SatMidInfo )
{
	int i,k;
	int Iterator;        /* 计算卫星轨道的迭代次数, 小于10次 */
	GPSTIME T_Tm;        /* 信号发射时刻的系统时间 */
	double  dt0, dt1;          /* 信号传播时间 [s]  */
	double  RotAng;      /* 信号传播时间内地球自转角度[Rad] */
	double  ClkOft, ClkSft;  /* 卫星钟差和钟速 */
	double  Mat[9], Pos0[3], Pos1[3], Vel0[3];
	double  dPos[3];

	

	dt1 = 0.075;     /*  初始概略传播时间 */
	Iterator = 0;
	T_Tm.Week = Time->Week;
	do {
		dt0 = dt1; 
		T_Tm.SecOfWeek = Time->SecOfWeek - dt0 - PreLeoClk/C_Light;

		if( ComputeGPSOrbitAndClockFullInfo( Prn, &T_Tm, GPSEph, IonPara, 
			Pos0, Vel0, &ClkOft, &ClkSft ) == false )
		{
			return false;
		}

		/*为了防止广播星历计算影响精密星历这里重新初始化*/
		/*
		T_Tm.SecOfWeek = Time->SecOfWeek - dt0 - PreLeoClk/C_Light;
		computeGPSfromSP3(T_Tm,Prn,5,Pos0,Vel0,&ClkOft,&ClkSft);
		computeGPSfromSP3(T_Tm,Prn,3,Pos0,Vel0,&ClkOft,&ClkSft);
		computeGPSfromSP3(T_Tm,Prn,4,Pos0,Vel0,&ClkOft,&ClkSft);*/
		/*  地球自转改正  */

		RotAng = Omega_WGS * dt0;
		Rotation_z( RotAng, Mat );
		MatrixMultiply( 3, 3, 3, 1, Mat, Pos0, Pos1 );

		for( i=0; i<3; i++ )
		{
			dPos[i] = Pos1[i] - PreLeoOrb[i];
		}

		dt1 = sqrt( VectDot( 3, 3, dPos, dPos ) ) / C_Light;

		Iterator++;  

	} while( (fabs( dt1-dt0 )>1E-8) && (Iterator < 10 ));

	SatMidInfo->SatClkOft = ClkOft;
	SatMidInfo->SatClkSft = ClkSft;

	SatMidInfo->Relativity = 0.0; /* 已在卫星钟差中改正, 此处为0 */

	/* 卫星速度改正 */

	MatrixMultiply( 3, 3, 3, 1, Mat, Vel0, SatMidInfo->SatVel );

	CopyArray( 3, SatMidInfo->SatPos, Pos1 );

	/*  卫星高度角计算 */

	SatMidInfo->Elevation = SatElev( SatMidInfo->SatPos, PreLeoOrb );

	/*  单频接收机的Tgd改正, 相对论效应已包含在钟差 

	SatMidInfo->SatClkOft = SatMidInfo->SatClkOft - GPSEph->TGD;  */

	/*   电离层延迟计算 */

	if( IonPara->IsValid )
	{
		SatMidInfo->IonoCorr = KlobucharForLEO( Time, SatMidInfo->SatPos, PreLeoOrb, IonPara );
	}
	else
	{
		SatMidInfo->IonoCorr = 0.0;
	}

	/*    对流层延迟改正, SPP使用概略改正模型  */

	if( Height <= 30000.0 )
	{
		SatMidInfo->TropCorr = 2.47/(sin(SatMidInfo->Elevation)+0.0121);
	}
	else
	{
		SatMidInfo->TropCorr = 0.0;
	}

	return true;
}


/***************************************************************************
//
// Klobuchar
//
// 目的: 使用Klobuchar模型, 计算单频接收机的电离层延迟改正量, 适用于地面测量

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
	double RcvPos[3], IONOPARA* IonoPara)
{
	int i;
	double correction;
	double RcvBLH[3];       /*  [Rad]  */
	double Mat[9], dPos[3], NEU[3];   /* To NEU matrix */
	double svE, svA;       /* satllite elevation and azimuth [ semi-circle] */
	double phi_u, lambda_u, psi, phi_i, lambda_i, phi_m;
	double iAMP, iPER, t, x, iF, t_iono;

	XYZToBLH( RcvPos, RcvBLH, R_WGS84, F_WGS84 );
	BLHToNEUMatrix( RcvBLH, Mat );

	for( i=0; i<3; i++ )
	{
		dPos[i] = SatPos[i] - RcvPos[i];
	}

	MatrixMultiply( 3, 3, 3, 1, Mat, dPos, NEU );

	svE = atan( NEU[2]/sqrt(NEU[0]*NEU[0]+NEU[1]*NEU[1]) ) / pi;
	svA = atan2( NEU[1], NEU[0] ) / pi;

	phi_u    = RcvBLH[0] / pi;
	lambda_u = RcvBLH[1] / pi;
	psi = 0.0137 / (svE + 0.11) - 0.022;

	phi_i = phi_u + psi * cos(svA*pi);
	if (phi_i > 0.416)
	{
		phi_i = 0.416;
	}
	if (phi_i < -0.416)
	{
		phi_i = -0.416;
	}

	lambda_i = lambda_u + psi * sin(svA*pi) / cos(phi_i*pi);

	phi_m = phi_i + 0.064 * cos((lambda_i - 1.617)*pi);

	iAMP = 0.0;
	iPER = 0.0;
	for (int n = 0; n < 4; n++)
	{
		iAMP += IonoPara->alpha[n] * pow(phi_m, n);
		iPER += IonoPara->beta[n] * pow(phi_m, n);
	}
	if (iAMP < 0.0)
	{
		iAMP = 0.0;
	}
	if (iPER < 72000.0)
	{
		iPER = 72000.0;
	}

	t = 43200.0 * lambda_i + fmod( Time->SecOfWeek, SECPERDAY );
	if (t >= 86400.0)
	{
		t = t - 86400.0;
	}
	if (t < 0)
	{
		t = t + 86400.0;
	}

	x = pi * (t - 50400.0) / iPER; // x is in radians

	/* 投影函数  */
	iF = 1.0 + 16.0 * pow(0.53 - svE, 3);    /* 地面测量使用 */

	t_iono = 0.0;
	if (fabs(x) < 1.57)
	{
		t_iono = iF * (5.0e-9 + iAMP * (1 - pow(x, 2)/2 + pow(x, 4)/24));
	}
	else
	{
		t_iono = iF * 5.0e-9;
	}

	correction = t_iono * C_Light;

	return correction;

}

/***************************************************************************
//
// KlobucharForLEO
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
	double RcvPos[3], IONOPARA* IonoPara)
{
	int i;
	double correction, alpha;
	double dPos[3], HIP, NormPos, Tmp; 
	double Mat[9], NEU[3], RcvBLH[3];      /* To NEU matrix */
	double phi_u, lambda_u, psi, phi_i, lambda_i, phi_m;
	double H0, dH, svE, svA;    /* satllite elevation and azimuth [ semi-circle] */
	double iAMP, iPER, t, x, iF, t_iono;

	H0 = 420000.0;                  /* 电离层电子密度最大的高度 */
	dH = 100000.0;

	NormPos = sqrt( VectDot( 3, 3, RcvPos, RcvPos ) );
	if( NormPos < 5000.0 || NormPos > ( R_WGS84 + 1E6 ) )     /* 接收机位置错误 */
	{
		return 0.0;
	}

	XYZToBLH( RcvPos, RcvBLH, R_WGS84, F_WGS84 );
	BLHToNEUMatrix( RcvBLH, Mat );

	for( i=0; i<3; i++ )
	{
		dPos[i] = SatPos[i] - RcvPos[i];
	}

	MatrixMultiply( 3, 3, 3, 1, Mat, dPos, NEU );

	svE = atan( NEU[2]/sqrt(NEU[0]*NEU[0]+NEU[1]*NEU[1]) )/pi;
	svA = atan2( NEU[1], NEU[0] )/pi;

	phi_u    = RcvBLH[0] / pi;
	lambda_u = RcvBLH[1] / pi;

	Tmp = exp( 1.0 - exp( ( -RcvBLH[2]+H0 )/dH ) );
	HIP = H0 - dH * log( 1.0 - log( ( 2.71828183 + Tmp )/2.0 ) ); /* e=2.71828183 */
	alpha = ( 1.359140915 - Tmp/2.0 ) / ( 2.71828183 - exp( 1.0 - exp( H0/dH ) ) );

	/* 计算地心夹角 zdh06411@163.com, caohaisheng@yahoo.com.cn */

	Tmp = asin( NormPos * cos(svE*pi) / ( R_WGS84 + HIP ) );
	psi = 0.5 - Tmp/pi - svE;
	iF = 1.0 / cos( Tmp );

	phi_i = phi_u + psi * cos(svA*pi);

	if (phi_i > 0.416)
	{
		phi_i = 0.416;
	}
	if (phi_i < -0.416)
	{
		phi_i = -0.416;
	}

	lambda_i = lambda_u + psi * sin(svA*pi) / cos(phi_i*pi);

	phi_m = phi_i + 0.064 * cos((lambda_i - 1.617)*pi);

	iAMP = 0.0;
	iPER = 0.0;
	for (int n = 0; n < 4; n++)
	{
		iAMP += IonoPara->alpha[n] * pow(phi_m, n);
		iPER += IonoPara->beta[n]  * pow(phi_m, n);
	}
	if (iAMP < 0.0)
	{
		iAMP = 0.0;
	}
	if (iPER < 72000.0)
	{
		iPER = 72000.0;
	}

	t = 43200.0 * lambda_i + fmod( Time->SecOfWeek, SECPERDAY );
	if (t >= 86400.0)
	{
		t = t - 86400.0;
	}
	if (t < 0)
	{
		t = t + 86400.0;
	}

	x = pi2 * (t - 50400.0) / iPER; // x is in radians

	t_iono = 0.0;
	if (fabs(x) < 1.57)
	{
		t_iono = iF * (5.0e-9 + iAMP * (1 - pow(x, 2)/2 + pow(x, 4)/24) );
	}
	else
	{
		t_iono = iF * 5.0e-9;
	}

	correction = t_iono * C_Light;

	return correction;

}

/***************************************************************************
//
// PointPositionRTOD
//
// 目的: 用GPS和GLONASS卫星的伪距, 进行单点定位, 如果有双频数据, 组合消除电离层延迟.
如果是单频数据,用改进的Klobuchar模型.         

//
// 输入参数:
//
GPSEph          GPS卫星星历[32]
IonPara         电离层参数
GLOEph          GLONASS卫星星历[32]
GloTmCorr       GLONASS时间与UTC(SU)系统间的改正值
EpochObs        历元的观测数据, 以此顺序计算卫星轨道
AprioriPos      接收机的先验坐标, 轨道积分预报轨道
AprioriClk      接收机的先验钟差[GPS,GLONASS]

输出参数

SatMidInfo      卫星轨道, 速度等中间计算过程信息, 用于滤波定轨
Result          单点定位结果(地固系)

返回值

参与定位的卫星数

算法流程:

先组织观测数据, 根据先验轨道信息( 如果是首次定位或重置后的定位, 先验信息为0.0),
计算导航卫星轨道和钟差信息, 然后对观测数据进行粗差探测, 如果是粗差数据, 该观测值
定义为不可用.


***************************************************************************/
int PointPositionRTOD( int graceType, GPSEPHREC* GPSEph, IONOPARA* IonPara,
	EpochObsData* Epoch, APRIORBIT* PreOrb, 
	SATMIDRESULT* SatMidInfo, PPRESULT* Result )
{
	int i, j, k, Iterator;                     /*  Iterator为单点定位迭代次数  */
	int SatNumUsed;      /*  单点定位使用的卫星数, GPS卫星数, GLONASS卫星数 */
	short Prn;
	short PosKind;            /* GPS:0, GLO:1, GPS&GLO:2  */

	int    PRValid;                            /* 观测值的伪距可用性和伪距观测值 */
	double PRange;
	double AprioriPos[3], AprioriClk[2], BLH[3];
	double Range, Ion, Height;                                 /*  接收机与导航卫星之间的距离  */
	double BPos[MAXCHANNUM];                      /* 伪距或多普勒的观测值-计算值 */
	double MeasA[MAXCHANNUM*4], MeasAT[MAXCHANNUM*4];          /* 观测矩阵 [MAXCHANNUM][4] */
	double Weight[MAXCHANNUM];                               /* 权矩阵, 单位阵只取对角线元素 */
	double Qvv[MAXCHANNUM];                               /* 观测值改正数的斜因数阵对角线元素 */
	double ATA[16], InvATA[16];
	double ATB[4];
	double dPos[4];               /* 定位结果  */
	double Residual[MAXCHANNUM] = {0.0};  /* 定位残差  */



	if( PreOrb->Validate )
	{
		CopyArray( 3, AprioriPos, PreOrb->LeoPos );
		CopyArray( 2, AprioriClk, PreOrb->LeoClk );
	}
	else
	{
		AprioriPos[0] = 0.0;
		AprioriPos[1] = 0.0;
		AprioriPos[2] = 0.0;
		AprioriClk[0] = 0.0;
		AprioriClk[1] = 0.0;
	}

	Iterator = 0;

	do{
		SatNumUsed = 0;
		XYZToBLH( AprioriPos, BLH, R_WGS84, F_WGS84 );
		Height = BLH[2];

		for( i=0; i<Epoch->SatNum; i++ )
		{
			Result->SatList[i].Prn = Epoch->SatObs[i].Prn;
			if(Epoch->SatObs[i].System!=GPS)
				continue;
			Result->SatList[i].System = Epoch->SatObs[i].System;

			PRValid = GetOneSatPseudoRange( Epoch->SatObs[i].System,
				&Epoch->SatObs[i], &PRange, &Ion );

			if(PRValid > 0 )
			{
				Prn = Epoch->SatObs[i].Prn;

				if( ComputeGPSSatOrbitAtSignalTrans( Prn, &Epoch->Time, AprioriPos, AprioriClk[0], 
					Height, &GPSEph[Prn-1], IonPara, &SatMidInfo[i])
					&& ( SatMidInfo[i].Elevation >= 10.0*Rad ))   /*  高度角大于10度 */
				{
					for(k=0;k<3;k++)
					{
						Epoch->SatObs[i].satPos[k] = SatMidInfo[i].SatPos[k];
					}

					Range = 0.0;
					for( k=0; k<3; k++ )
					{
						Range = Range + pow( Epoch->SatObs[i].satPos[k]-AprioriPos[k], 2.0 );
					}

					Range = sqrt( Range ); 
					MeasA[SatNumUsed*4+0] = ( AprioriPos[0] - Epoch->SatObs[i].satPos[0] )/Range;
					MeasA[SatNumUsed*4+1] = ( AprioriPos[1] - Epoch->SatObs[i].satPos[1] )/Range;
					MeasA[SatNumUsed*4+2] = ( AprioriPos[2] - Epoch->SatObs[i].satPos[2] )/Range;
					MeasA[SatNumUsed*4+3] = 1.0;             /* 接收机GPS系统钟差系数 */
					Weight[SatNumUsed]    = 1.0;  // * pow( sin(SatMidInfo[i].Elevation), 2.0 );

					BPos[SatNumUsed] = PRange - Range - AprioriClk[0] 
					+ SatMidInfo[i].SatClkOft * C_Light - SatMidInfo[i].TropCorr; 

					if( PRValid == 1 )
					{		
						BPos[SatNumUsed] = BPos[SatNumUsed] - SatMidInfo[i].IonoCorr;//如果是单频
					}

					Result->SatList[i].Status = 1;        /* 参与单点定位计算 */
					
					SatNumUsed++;
				}
				else    
				{
					Epoch->SatObs[i].Used = 0;     /* 无星历或高度角低 */                      

				}
				
			}
			Epoch->SatObs[i].Ion=Ion;			  /*保存由双频得到的电离层延迟，如果不是双频，则为0*/
		}

		/*  组成残差数据集合, 进行粗差检查,   */
		//这里要有先验位置和精度才能进行检验，通过观测值残差进行接收机钟差检验方法

		if( PreOrb->Validate )
		{
			if( DetectPseudoRangeOutliers(graceType, PreOrb->OrbitAccuracy, BPos, Epoch , Result->SatList ) )
			{
				DeleteOutlierFromObsList( &SatNumUsed, BPos, MeasA, Weight, Result->SatList);
			}

			PreOrb->Validate = 0;   /*  每个历元只探测一次  */
		}

		/*  组成误差方程  */

		for( k=SatNumUsed; k<MAXCHANNUM; k++)  /* 多余项清零 */
		{
			BPos[k] = 0.0;
			Weight[k] = 0.0;
			for( j=0; j<4; j++ )
			{
				MeasA[k*4+j] = 0.0;
			}
		}

		if( SatNumUsed>=4 )  /* 不考虑GLONASS卫星 */
		{

			MatrixTranspose( SatNumUsed, 4, MeasA, MeasAT );
			MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 4, MeasAT, Weight, MeasA, ATA );
			MatrixInv( 4, ATA, InvATA );

			MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 1, MeasAT, Weight, BPos, ATB );
			MatrixMultiply( 4, 4, 4, 1, InvATA, ATB, dPos );

			AprioriClk[0] = AprioriClk[0] + dPos[3];

		}
		else
		{
			//            printf("\nCannot position because of only %2d satellites at %10.3f epoch.\n",
			//                SatNumUsed, Epoch->Time.SecOfWeek );

			break;
		}

		for( k=0; k<3; k++ )
		{
			AprioriPos[k] = AprioriPos[k] + dPos[k];
		}

		Iterator++;
		Result->Iterator = Iterator;
		Result->Coverage = VectDot( 3, 3, dPos, dPos );

	}while ( (Result->Coverage>1E-6) && (Iterator < 10) );
	//循环结果计算高度角
	for(i=0;i<Epoch->SatNum;i++)
		Epoch->SatObs[i].elevation=SatMidInfo[i].Elevation;
	Result->SatNumUsed = SatNumUsed;

	CopyArray( 3, Result->Position, AprioriPos );
	CopyArray( 2, Result->RcvClkOft, AprioriClk );

	XYZToBLH( Result->Position, Result->PosBLH, R_WGS84, F_WGS84 );
	BLHToNEUMatrix( Result->PosBLH, Result->MatrixFromBLHtoNEU );

	/* 计算观测值残差 */

	MatrixMultiply( SatNumUsed, 4, 4, 1, MeasA, dPos, Residual );//
	for( i=0; i<SatNumUsed; i++ )
	{
		Residual[i] = Residual[i] - BPos[i];
	}

	/* 计算观测值改正数的协因数阵 */

	ComputeQvv( SatNumUsed, 4, MeasA, InvATA, MeasAT, Weight, Qvv );

	Result->PDOP = sqrt( InvATA[0] + InvATA[5] + InvATA[10] );
	Result->GDOP = sqrt( InvATA[0] + InvATA[5] + InvATA[10] + InvATA[15] );
	Result->HDOP = sqrt( InvATA[0] + InvATA[5] );
	Result->VDOP = sqrt( InvATA[10] );
	Result->TDOP = sqrt( InvATA[15] );

	if( SatNumUsed > 5 ) /* 多于5颗卫星才计算定位中误差, 增强可靠性 */
	{
		Result->SigmaPos = 0.0;

		for( j=0, i=0; i<Epoch->SatNum; i++ )
		{
			if( Result->SatList[i].Status == 1 )
			{
				Result->Residual[i] = Residual[j];
				Result->SigmaPos = Result->SigmaPos + Weight[i]*pow( Residual[i], 2.0 );
				j++;
			}
		}
		Result->SigmaPos = sqrt( Result->SigmaPos / (SatNumUsed-4) );
	}
	else
	{
		Result->SigmaPos = 999.99;
	}

	/* 验后粗差探测 */

	if( SatNumUsed > 5 )
	{
		CheckPostResidual_W( Epoch->SatNum, 3.0, Residual, Qvv , Result );
		//这里暂时不进行粗差检验，移到外面
	}
	/*RJ-2016-9-20
	**该模块是为了把单点定位中验出粗差的卫星作为不可用
	*/
	for(i=0;i<Epoch->SatNum;i++)
	{
		if(Result->SatList[i].Status!=1&&(Epoch->SatObs[i].Used==1))
			Epoch->SatObs[i].Used = 0;
	}
	if( Result->Coverage < 1E-6 && Result->Iterator < 10 &&
		Result->PDOP < 8.0 && Result->SigmaPos < 5.0 && SatNumUsed > 5 )
	{
		Result->IsSuccess = true;
	}
	else
	{
		Result->IsSuccess = false;
	}

	return SatNumUsed;
}

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
	double Pll[], double Qvv[] )
{
	int i;
	double AN[MAXCHANNUM*5], ANAT[MAXCHANNUM*MAXCHANNUM];

	MatrixMultiply( m, n, n, n, A, InvN, AN );
	MatrixMultiply( m, n, n, m, AN, AT, ANAT );

	for( i=0; i<m; i++ )
	{
		Qvv[i] = 1.0/Pll[i] - *(ANAT+i*m+i);

	}
}
/****************************************************************************
ComputeQLL

目的: 计算平差后观测值协因数阵
参数: int m
参数: int n
参数: const double A[]
参数: const double InvN[]
参数: const double AT[]
参数: const double Pll[]
参数: double QLL[]
****************************************************************************/
void ComputeQLL( int m, int n, const double A[],const double InvN[],const double AT[], 
	const double Pll[], double QLL[] )
{
	int i;
	double Q[MAXCHANNUM];
	for (i=0;i<m;i++)
	{
		Q[i]=1.0/Pll[i];
	}
	double AN[MAXCHANNUM*5], ANAT[MAXCHANNUM*MAXCHANNUM];
	double QANAT[MAXCHANNUM*MAXCHANNUM],QANATQ[MAXCHANNUM*MAXCHANNUM];
	MatrixMultiply( m, n, n, n, A, InvN, AN );
	MatrixMultiply( m, n, n, m, AN, AT, ANAT );
	MatrixMultiply2(1,m,m,m,ANAT,Q,QANAT);
	MatrixMultiply2(1,m,m,m,QANAT,Q,QANATQ);

	for( i=0; i<m; i++ )
	{
		QLL[i] = Q[i] - *(QANATQ+i*m+i);

	}
}
/***************************************************************************
//
// PointPositionVelocityDetermination
//
// 目的: 根据卫星星历和一个历元的观测数据的D1, 进行单点测速(包括仅GPS, 
仅GLONASS, GPS/GLONASS组合定位), 为滤波提供初始速度信息. 调试情况下,
CHAMP卫星等数据没有多普勒观测数据, 故先调用CreateDopplerObs函数, 生成
多普勒观测数据, 否则不能用该子程序.

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
	SATMIDRESULT* SatMidInfo, PPRESULT* Result )
{
	int i, j, k;    
	int Prn;
	int SatNumUsed;                               /*  单点定位使用的卫星数 */

	double Range;                                 /*  接收机与导航卫星之间的距离  */
	double BVel[MAXCHANNUM];                /* 多普勒的观测值-计算值 */
	double MeasA[MAXCHANNUM*4], MeasAT[MAXCHANNUM*4];          /* 观测矩阵 [MAXCHANNUM][5] */
	double Weight[MAXCHANNUM];                               /* 权矩阵, 单位阵只取对角线元素 */
	double ATA[16], InvATA[16];
	double ATB[4];
	double dPos[4];                                /* 定位结果  */
	double Residual[MAXCHANNUM];                   /* 定位残差  */
	double GLOL1WaveLen;                           /* 每颗GLONASS卫星的L1波长 */

	SatNumUsed = 0;

	/* 准备观测数据  */

	for( i=0; i<Epoch->SatNum; i++ )
	{
		if( Epoch->SatObs[i].Used != 1 )
		{
			continue;
		}

		if( Epoch->SatObs[i].System == GPS )   
		{
			for( j=0; j<MAXOBSTYPENUM; j++ )
			{
				if( (Epoch->SatObs[i].Data[j].Availability == true) 
					&& ( Epoch->SatObs[i].Data[j].Type == D1 ) )
				{
					for( k=0; k<3; k++ )
					{
						dPos[k] = SatMidInfo[i].SatPos[k] - Result->Position[k];
					}

					Range = sqrt( VectDot( 3, 3, dPos, dPos ) );

					MeasA[SatNumUsed*4+0] = -dPos[0] / Range;
					MeasA[SatNumUsed*4+1] = -dPos[1] / Range;
					MeasA[SatNumUsed*4+2] = -dPos[2] / Range;
					MeasA[SatNumUsed*4+3] = 1.0;             /* 接收机GPS系统频差系数 */

					Weight[SatNumUsed]    = 2 * pow( sin(SatMidInfo[i].Elevation), 2.0 );

					BVel[SatNumUsed] = -Epoch->SatObs[i].Data[j].Obs*C_Light/FG1_Freq  /* 注意多普勒观测值的符号 */
						- VectDot( 3, 3, dPos, SatMidInfo[i].SatVel )/Range
						+ SatMidInfo[i].SatClkSft*C_Light;

					SatNumUsed++;
				}
			}
		}

		if( Epoch->SatObs[i].System == GLONASS )  
			continue;
	}

	/*  组成误差方程  */

	for( k=SatNumUsed; k<MAXCHANNUM; k++)
	{
		BVel[k] = 0.0;
		for( j=0; j<4; j++ )
		{
			MeasA[k*4+j] = 0.0;
		}
	}

	for( k=0; k<4; k++ )
	{
		dPos[k] = 0.0;
	}

	if( SatNumUsed >= 4 )
	{
		MatrixTranspose( SatNumUsed, 4, MeasA, MeasAT );
		MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 4, MeasAT, Weight, MeasA, ATA );
		MatrixInv( 4, ATA, InvATA );

		MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 1, MeasAT, Weight, BVel, ATB );
		MatrixMultiply( 4, 4, 4, 1, InvATA, ATB, dPos );  /* dPos包含速度信息和频差信息 */

		/* 计算观测值残差和精度评定 */

		MatrixMultiply( SatNumUsed, 4, 4, 1, MeasA, dPos, Residual );

		Result->SigmaVel = 0.0;
		for( i=0; i<SatNumUsed; i++ )
		{
			Residual[i] = Residual[i] - BVel[i];
			Result->SigmaVel = Result->SigmaVel + Weight[i]*pow( Residual[i], 2.0 );
		}

	}

	/* 计算观测残差标准差 */

	if( SatNumUsed > 4 )
	{
		Result->SigmaVel = sqrt( Result->SigmaVel / (SatNumUsed-4) );
	}       
	else
	{
		Result->SigmaVel = 999.99;
	}

	/*  保存测速结果 */

	CopyArray( 3, Result->Velocity, dPos );
	Result->RcvClkSft = dPos[3];//解出接收机钟速

	MatrixMultiply( 3, 3, 3, 1, Result->MatrixFromBLHtoNEU, Result->Velocity, Result->VelNEU );

	return SatNumUsed;

}

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
	EpochObsData* Epoch, SATLIST* SatList )
{
	int  i, j;
	int  GPSSatNum, GLOSatNum;

	double GVal[MAXCHANNUM]={0.0};  /*  临时存储GPS, GLONASS的o_c */
	double GMean ;                   /* 均值 */
	double GStd,  GLim;        /* 标准差, 限差 */
	bool   HasOutlier = false;
	FILE   *ocFile;
	char   graceName[20]="";
	char   outName[30]="";
	//graceType for the out file name of o-c
	if(graceType==0)
		sprintf(graceName,"GRACE-A");
	else
		sprintf(graceName,"GRACE-B");
	
	GPSSatNum = 0;

	GLim   = Accuracy;

	/* 将GPS和GLONASS卫星的O_C数据分开 */

	for( j=0,i=0; i<MAXCHANNUM; i++ )
	{
		if( (SatList[i].Status==1) && (SatList[i].System==GPS) )
		{
			char cPrn[3];
			sprintf(outName,"OutFile\\OC\\%s_%d.txt",graceName,SatList[i].Prn);
			ocFile = fopen(outName,"a+");
			fprintf(ocFile,"%10d%14.3f%14.3f\n",Epoch->Time.Week,Epoch->Time.SecOfWeek,O_C[j]);
			fclose(ocFile);
			GVal[GPSSatNum] = O_C[j];
			GPSSatNum++;
			j++;
		}
	}

	mbbub( GPSSatNum, GVal );     /* 排序, 剔除两端点计算标准差 */
	GMean = GVal[GPSSatNum/2];    /* 取排序后的中间结果作为平均值 */

	GStd = 0.0;
	for( j=0, i=1; i<GPSSatNum-1; i++ )   /* 去掉端点计算标准差 */
	{
		if( fabs(GVal[i]-GMean) < 20*GLim )
		{
			GStd = GStd + (GVal[i]-GMean) * (GVal[i]-GMean);
			j++;
		}
	}
	if( j>0 )  GLim = sqrt( GStd/j );
	
	GLim = __max( GLim, Accuracy );

	if( GLim < 3.0 )    /* 防止标准差计算不准确,造成误删除 */
	{
		GLim = 3.0;
	}

	for( j=0,i=0; i<MAXCHANNUM; i++ )
	{
		if( (SatList[i].Status==1) && (SatList[i].System==GPS) )
		{
			if( fabs(O_C[j] - GMean ) >= 3*GLim)
			{
				// 				printf( "Outlier: %10.2f G%02d %12.3f %12.3f %12.3f\n", 
				// 					Epoch->Time.SecOfWeek, SatList[i].Prn, O_C[j], GMean, O_C[j]-GMean );
				Epoch->SatObs[i].Used = -1;
				SatList[i].Status = -1;
				HasOutlier = true;
			}
			j++;
		}
	}
	return HasOutlier;
}
/***************************************************************************
//
// DetectPseudoRangeOutliers
//
// 目的: 根据先验坐标或最小二乘估计的近似坐标, 计算出所有卫星的O-C
(观测值减去计算值), 比较每颗GPS卫星计算出来的接收机相对钟差是否一致.

//
// 输入参数:
//
Accuracy        先验轨道数据的精度指标, 用来确定粗差剔除的阈值.
O_C             LC观测值减去观测值的计算值, 即每颗卫星的接收机钟差.
Epoch           历元观测数据, 探测到粗差后将该卫星设置为粗差观测值,不可用
SatList         与O_C相对应, 以确知粗差观测值来自于那颗卫星

返回参数

当探测到粗差时, 返回1, 否则返回0.

***************************************************************************/
bool DetectOutliers_rel( double Accuracy, double O_C[], 
	Common11Obs * CurComObs, SATLIST* SatList, double sigma )
{
	int  i, j;
	int  GPSSatNum;

	double GVal[MAXCHANNUM]={0.0};						/*  临时存储GPSo_c */
	double GMean;								/* 均值 */
	double GStd, GLim;						/* 标准差, 限差 */
	bool   HasOutlier = false;

	GPSSatNum = 0;

	GLim = Accuracy;

	/* 将GPS和GLONASS卫星的O_C数据分开 */

	for( j=0,i=0; i<MAXCHANNUM; i++ )
	{
		if( SatList[i].Status==1 )
		{
			GVal[GPSSatNum] = O_C[j];
			GPSSatNum++;
			j++;
		}
	}

	mbbub( GPSSatNum, GVal );     /* 排序, 剔除两端点计算标准差 */

	GMean = GVal[GPSSatNum/2];    /* 取排序后的中间结果作为平均值 */

	GStd = 0.0;
	for( j=0, i=1; i<GPSSatNum-1; i++ )   /* 去掉端点计算标准差 */
	{
		if( fabs(GVal[i]-GMean) < 20*GLim )
		{
			GStd = GStd + (GVal[i]-GMean) * (GVal[i]-GMean);
			j++;
		}
	}
	if( j>0 )  GLim = sqrt( GStd/j );
	GLim = __max( GLim, Accuracy );
	
	if( GLim < sigma )    /* 防止标准差计算不准确,造成误删除 */
	{
		GLim = sigma;
	}

	for( j=0,i=0; i<MAXCHANNUM; i++ )
	{
		if( SatList[i].Status==1 )
		{
			if( fabs(O_C[j] - GMean ) >= 3*GLim)
			{
				// 				printf( "Outlier: %10.2f G%02d %12.3f %12.3f %12.3f\n", 
				// 					Epoch->Time.SecOfWeek, SatList[i].Prn, O_C[j], GMean, O_C[j]-GMean );
				SatList[i].Status = -1;
				if(sigma>1)
					CurComObs->comobs[SatList[i].index].used= -1;
				else
					CurComObs->comobs[SatList[i].index].flag= -1;
				HasOutlier = true;
			}
			j++;
		}
	}
	return HasOutlier;
}
/***************************************************************************
//
// DeleteOutlierFromLC
//
// 目的: 当探测出有粗差观测值后, 此函数将删除粗差卫星的观测矩阵/加权系数/残差数据等,
避免粗差数据参与最小二乘估计

//
// 输入参数:
//
GPSSatNum        GPS卫星数.
Resid            观测值-计算值
A                观测系数矩阵
P                观测值权矩阵

***************************************************************************/
void DeleteOutlierFrom_rel( int* GPSSatNum, double Resid[],
	double A[], double P[], SATLIST* SatList )
{
	int i, j, k;
	int satNum=*GPSSatNum;
	for( j=0,i=0; i<satNum; i++ )
	{
		if( SatList[i].Status == 1 )  /* 可用 */
		{
			j++;
		}
		else                /* 有粗差 */
		{
			*GPSSatNum = *GPSSatNum - 1;
			satNum=satNum-1;
			for( k=j; k<satNum-1; k++ )
			{
				Resid[k] = Resid[k+1];
				P[k]     = P[k+1];
				A[k*5+0] = A[(k+1)*5+0];
				A[k*5+1] = A[(k+1)*5+1];
				A[k*5+2] = A[(k+1)*5+2];
				A[k*5+3] = A[(k+1)*5+3];
				A[k*5+4] = A[(k+1)*5+4];
			}

		}
	}

}
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
	double A[], double P[], SATLIST* SatList )
{
	int i, j, k;

	for( j=0,i=0; i<MAXCHANNUM; i++ )
	{
		if( SatList[i].Status == 0 ) /* 没有星历, 跳过 */
		{
			continue;            
		}
		else if( SatList[i].Status == 1 )  /* 可用 */
		{
			j++;
		}
		else                /* 有粗差 */
		{
			if( SatList[i].System == GPS )
			{
				*GPSSatNum = *GPSSatNum - 1;
			}
			for( k=j; k<MAXCHANNUM-1; k++ )
			{
				Resid[k] = Resid[k+1];
				P[k]     = P[k+1];
				A[k*4+0] = A[(k+1)*4+0];
				A[k*4+1] = A[(k+1)*4+1];
				A[k*4+2] = A[(k+1)*4+2];
				A[k*4+3] = A[(k+1)*4+3];
			}

		}
	}

}

/***************************************************************************
//
// CheckPostResidual_W
//
// 目的: 在单点定位之后, 对验后残差进行假设检验_W检验. 检查结果的正确性和每颗卫星的残差
是否合理. 如果残差被确认为粗差, 则标识.

//
// 输入参数:
//
n                历元卫星总数
Residual         单点定位观测值残差
Qvv              观测值残差改正数
Result           单点定位的结果, 探测到粗差数据, 设置Result中SatList参数
***************************************************************************/
void CheckPostResidual_W( int n, double sigma0, double Residual[], double Qvv[], PPRESULT* Result )
{
	int i, j;
	double w;                    /* w检验的统计量 */
	//double sigma0 = 3.0;         /* 观测值的统计精度 */

	for( j=0, i=0; i<n; i++ )
	{
		if( Result->SatList[i].Status == 1 )
		{
			w = Residual[j] / sigma0 / sqrt( Qvv[j] );

			if( fabs(w) > 5.0 )
			{
				Result->SatList[i].Status = -2;  /* 该观测值在滤波中将不使用 */

				printf( "W-Test %5.2f Failed.\n", w );
			}

			j++;

		}
	}

}

/****************************************************************************
CheckPostResidual_t

目的: t分布检验-与w检验同样适用于独立观测值
参数: int n
参数: double sigma0
参数: double Residual[]
参数: double Qvv[]
参数: PPRESULT * Result
****************************************************************************/
void CheckPostResidual_t( int n, double postSigma, double Residual[], double QLL[], PPRESULT* Result )
{
	int i, j;
	double t;                    /* t检验的统计量 */

	for( j=0, i=0; i<n; i++ )
	{
		if( Result->SatList[i].Status == 1 )
		{
			t = Residual[j] / postSigma / sqrt( QLL[j] );

			if( fabs(t) > 7.0 )					 /*当自由度大于2时，7保证alpha/2<0.01，5保证alpha/2<0.02*/
			{
				Result->SatList[i].Status = -2;  /* 该观测值在滤波中将不使用 */

				printf( "t-Test %5.2f Failed.\n", t );
			}

			j++;

		}
	}

}

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
PPResult        单点定位结果(地固系)

输出参数

KFState         滤波定轨结果

返回值

参与定位的卫星数

***************************************************************************/
int EKFilterRTODA( double Mat[9], SATMIDRESULT* SatMidInfo_A, EpochObsData* Epoch_A,
	PPRESULT* PPResult, EKFSTATE* KFState )
{
	int i, k;
	int SatUsedNum[2];

	int  PRValid;
	double PRange, Ion,Phase;

	double Range;         /* 接收机与导航卫星之间的距离计算值 */
	double dPos[3];       /* 接收机坐标与导航卫星位置之差 */

	double ApriState[ADIMENSION], AECF[3], ARange; /*验前状态参数，用于状态估计信息使用*/
	double H[ADIMENSION];  /* 观测方程线性化系数向量 */
	double O_C;              /* 观测值减去计算值       */

	double O_C0, O_C1;  /* 滤波更新前/后的观测值-计算值 */

	PRValid = 0;
	SatUsedNum[0] = 0;
	SatUsedNum[1] = 0;
	if( PPResult[0].IsSuccess == true )
	{
		KFState->StateA[6] = PPResult[0].RcvClkOft[0];
	}
	if( PPResult[0].SatNumUsed < 4 )  /* 小于4颗卫星不更新 */
	{
		return 0 ;
	}

	CopyArray( ADIMENSION, ApriState, KFState->StateA );

	//A星单点测量更新
	for( i=0; i<Epoch_A->SatNum; i++ )
	{
		for( k=0; k<ADIMENSION; k++ )
		{
			H[k] = 0.0;
		}

		if( PPResult[0].SatList[i].Status == 1 )
		{
			MatrixMultiply( 3, 3, 3, 1, Mat, &KFState->StateA[0], &KFState->StateInECEFA[0] );

			PRValid = GetOneSatPseudoRange( Epoch_A->SatObs[i].System, 
				&Epoch_A->SatObs[i], &PRange, &Ion );
			Range = 0.0;
			for( k=0; k<3; k++ )
			{
				dPos[k] = KFState->StateInECEFA[k] - SatMidInfo_A[i].SatPos[k];
				Range = Range + dPos[k]*dPos[k];
			}
			Range = sqrt( Range );

			H[0] = ( dPos[0]*Mat[0] + dPos[1]*Mat[3] + dPos[2]*Mat[6] )/Range;//前面接收机位置已经转到地固系了这里为什么还要乘以转换矩阵。
			H[1] = ( dPos[0]*Mat[1] + dPos[1]*Mat[4] + dPos[2]*Mat[7] )/Range;//这里是把地固系转回惯性系，因为滤波状态是在惯性系下
			H[2] = ( dPos[0]*Mat[2] + dPos[1]*Mat[5] + dPos[2]*Mat[8] )/Range;//但是新息是标量，在地固系和惯性系下值不变

			if(PPResult[0].SatList[i].System == GPS )
			{
				H[6] = 1.0;
				O_C = PRange - Range - KFState->StateA[6] + SatMidInfo_A[i].SatClkOft*C_Light;

				if( PRValid == 1 )
				{
					O_C = O_C - SatMidInfo_A[i].IonoCorr;   // 
				}
			}

			//A星绝对测量更新权
			if( EKFMeasureUpdateA(O_C, 5, H, KFState ) )//这里5是观测噪声R，为啥设为5
			{
				SatUsedNum[0]++;				
			}
			else
			{
				printf( "A星 Kalman Filter Failed:%10.3f Prn%3d %10.3f\n", 
					Epoch_A->Time.SecOfWeek, Epoch_A->SatObs[i].Prn, O_C );
			}

		}
	}
	
	KFState->SatNumUsed[0] = SatUsedNum[0];
	KFState->ApriSigma[0] = 0.0;
	KFState->PostSigma[0] = 0.0;

	MatrixMultiply( 3, 3, 3, 1, Mat, &ApriState[0], AECF );   /*验前状态参数 */
	MatrixMultiply( 3, 3, 3, 1, Mat, &KFState->StateA[0], KFState->StateInECEFA );

	for( i=0; i<Epoch_A->SatNum; i++ )
	{
		if( PPResult[0].SatList[i].Status != 1 )
		{
			continue;
		}

		PRValid = GetOneSatPseudoRange( Epoch_A->SatObs[i].System, 
			&Epoch_A->SatObs[i], &PRange, &Ion );

		Range = 0.0;
		ARange = 0.0;
		for( k=0; k<3; k++ )
		{	
			ARange = ARange + pow( AECF[k]-SatMidInfo_A[i].SatPos[k], 2.0 );
			Range  = Range  + pow( KFState->StateInECEFA[k]-SatMidInfo_A[i].SatPos[k], 2.0 );
		}
		Range = sqrt( Range );
		ARange = sqrt( ARange );

		if ( Epoch_A->SatObs[i].System == GPS )
		{
			O_C0  = PRange - ARange - KFState->StateA[6]   /*验前残差,使用验后钟差 */
			+ SatMidInfo_A[i].SatClkOft*C_Light + SatMidInfo_A[i].Relativity;
			O_C1 = PRange - Range - KFState->StateA[6] /* 验后残差 */
			+ SatMidInfo_A[i].SatClkOft*C_Light + SatMidInfo_A[i].Relativity;

			if( PRValid == 1 )
			{
				O_C0  = O_C0  - SatMidInfo_A[i].IonoCorr; 
				O_C1  = O_C1  - SatMidInfo_A[i].IonoCorr; 
			}
		}

		KFState->ApriSigma[0] = KFState->ApriSigma[0] + pow( O_C0, 2.0);
		KFState->PostSigma[0] = KFState->PostSigma[0] + pow( O_C1, 2.0 );

		/* 时间 Prn SPP残差 预报残差 更新后残差 双频电离层改正 单频电离层改正 卫星高度角 */

		if( PPResult[0].SatList[i].Status == 1 )
		{
			/*fprintf( FRESIDUAL, " %10.1f %3d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %2d \n", 
				Epoch->Time.SecOfWeek, Epoch->SatObs[i].Prn, 
				PPResult->Residual[i], O_C0, O_C1, Ion,
				SatMidInfo[i].IonoCorr, SatMidInfo[i].Elevation*Deg,
				PPResult->PDOP,PPResult->SigmaPos, PPResult->SatList[i].Status );*/
		}
	}

	if( KFState->SatNumUsed[0] > 0 )
	{
		KFState->ApriSigma[0] = sqrt( KFState->ApriSigma[0] / SatUsedNum[0] );
		KFState->PostSigma[0] = sqrt( KFState->PostSigma[0] / SatUsedNum[0] );
	}
	else
	{
		KFState->ApriSigma[0] = 999.0;
		KFState->PostSigma[0] = 999.0;
	}

	// 	fprintf( FRESIDUAL, " %10.1f %8.3f %8.3f %8.3f %2d %2d %2d, 111111 \n", 
	// 				Epoch->Time.SecOfWeek, PPResult->SigmaPos, KFState->ApriSigma,
	// 				KFState->PostSigma, Epoch->SatNum, PPResult->GPSSatNum+PPResult->GLOSatNum,
	// 				KFState->SatNum );


	EKFConvergencyCheck(0,KFState);

	//B星单点测量更新

	return SatUsedNum[0];//A星B星共同卫星数
}
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
void EKFConvergencyCheck(int graceType, EKFSTATE* KFState )
{
	double limit;//判断发散sigma限值
	double Cov;
	if(graceType>=1){//相对测量更新
		graceType=1;
		limit=10.0;
		Cov=KFState->CovRel[0]+KFState->CovRel[ADIMENSION+1]+KFState->CovRel[2*RELDIMENSION+2];
	}
	else
	{//A星测量更新
		limit=20.0;
		Cov=KFState->CovA[0]+KFState->CovA[ADIMENSION+1]+KFState->CovA[2*RELDIMENSION+2];
		
	}
	switch( KFState->KFConvergency[graceType] )
	{
	case 0:
		if( sqrt(Cov)< 5.0 )  /* 表示已收敛，位置误差 */
		{
			KFState->KFConvergency[graceType] = 1;
		}
		break;

	case 1:
		if( KFState->SatNumUsed[graceType] > 0 )  /* 经过滤波更新以后, 残差仍然很大 */
		{
			if( KFState->PostSigma[graceType] > limit )
			{
				KFState->KFConvergency[graceType] = 2;
				KFState->IsInitial[graceType] = 0;      /* 需要重新初始化 */
			}

		}
		break;

	default:
		break;
	}

	if( KFState->KFConvergency[graceType] == 2 )
	{
		printf( "Epoch %10.1f Grace%d Kalman filtering diverged!.\n", KFState->Time.SecOfWeek,graceType+1 );
	}
}
