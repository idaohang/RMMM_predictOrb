/****************************************************************************
目的：    动力学轨道积分,使用RK4, RKF4,以及轨道差值等函数

编写时间：2008.11.26
版本:     V1.1
版权：    武汉大学

2009.3.11  修改了Hermite5函数中的一个错误
****************************************************************************/

#include <math.h>
#include "RTOD_Const.h"
#include "CommonFuncs.h"
#include "RefSys.h"
#include "OrbitInteg.h"


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
/*
void RK4Step( MJDTIME* Mjd_GPS, double step, double Y0[6], DYNMODELPARA* Para )
{
	int i;
	double Pos[3];
	double Vel[4][3], Acc[4][3];
	double Prec[9], Nut[9], GH[9], Pole[9], Mat[9], E[9], T[9];
	MJDTIME Mjd_TT, Mjd_UTC, Mjd_UT1;
	EOPPARA CurrEop;
	double h = step / SECPERDAY;

	Mjd_TT.Days = Mjd_GPS->Days;
	Mjd_TT.FracDay = Mjd_GPS->FracDay + (TT_TAI-GPST_TAI)/SECPERDAY;
	Mjd_UTC.Days = Mjd_TT.Days;
	Mjd_UTC.FracDay = Mjd_TT.FracDay - TT_UTC() / SECPERDAY;
	InterposeEOP( &Mjd_UTC, &CurrEop );

	Mjd_UT1.Days = Mjd_UTC.Days;
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	PrecMatrix( MJD_J2000, &Mjd_TT, Prec );
	NutMatrixSimple( &Mjd_TT, Nut );
	MatrixMultiply( 3, 3, 3, 3, Nut, Prec, T );

	GHAMatrix( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, Mat );
	MatrixMultiply( 3, 3, 3, 3, Mat, T, E );

	for( i=0; i<3; i++ )
	{
		Pos[i] = Y0[i];
		Vel[0][i] = Y0[3+i];
	}

	AccelMain(TODO, &Mjd_TT, Pos, Vel[0], E, T, Para , Acc[0] );

	//   compute right function at the second time 

	Mjd_TT.FracDay += h/2.0;  
	Mjd_UTC.FracDay += h/2.0;
	InterposeEOP( &Mjd_UTC, &CurrEop );
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	GHAMatrix( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, Mat );
	MatrixMultiply( 3, 3, 3, 3, Mat, T, E );      //短时间内没有考虑岁差和章动变化

	for( i=0; i<3; i++ )
	{
		Pos[i]    = Y0[i]   + Vel[0][i] * step / 2.0;
		Vel[1][i] = Y0[3+i] + Acc[0][i] * step / 2.0;
	}
	AccelMain(TODO, &Mjd_TT, Pos, Vel[1], E, T, Para , Acc[1] );

	//   compute right function at the third time 

	for( i=0; i<3; i++ )
	{
		Pos[i]    = Y0[i]   + Vel[1][i] * step / 2.0;
		Vel[2][i] = Y0[3+i] + Acc[1][i] * step / 2.0;
	}
	AccelMain(TODO, &Mjd_TT, Pos, Vel[2], E, T, Para , Acc[2] );

	//   compute right function at the fourth time 

	Mjd_TT.FracDay  += h/2.0;  
	Mjd_UTC.FracDay += h/2.0;
	InterposeEOP( &Mjd_UTC, &CurrEop );
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	GHAMatrix( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, Mat );
	MatrixMultiply( 3, 3, 3, 3, Mat, T, E );      //短时间内没有考虑岁差和章动变化

	for( i=0; i<3; i++ )
	{
		Pos[i]    = Y0[i]   + Vel[2][i] * step;
		Vel[3][i] = Y0[3+i] + Acc[2][i] * step;
	}
	AccelMain(TODO, &Mjd_TT, Pos, Vel[3], E, T, Para , Acc[3] ); 

	Mjd_GPS->FracDay += h;
	if( Mjd_GPS->FracDay >= 1.0 )    // 跨天检查 
	{
		Mjd_GPS->FracDay -= 1.0;
		Mjd_GPS->Days    += 1;
	}

	for( i=0; i<3; i++ )
	{
		Y0[i]   += ( Vel[0][i] + 2.0*Vel[1][i] + 2.0*Vel[2][i] + Vel[3][i] ) * step / 6.0;
		Y0[i+3] += ( Acc[0][i] + 2.0*Acc[1][i] + 2.0*Acc[2][i] + Acc[3][i] ) * step / 6.0;
	}

}
*/

/***************************************************************************
//
// RKF4Step
//
// Purpose:
//
//   龙格-库塔-Felberg4阶单步积分器函数, 只积分卫星位置和速度
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
/*
void RKF4Step( MJDTIME* Mjd_GPS, double step, double Y0[6], DYNMODELPARA* Para )
{
	int i;
	double Pos[3];
	double Vel[5][3], Acc[5][3];
	double Prec[9], Nut[9], GH[9], Pole[9], Mat[9], E[9], T[9];
	MJDTIME Mjd_TT, Mjd_UTC, Mjd_UT1;
	EOPPARA CurrEop;
	double h = step / SECPERDAY;

	Mjd_TT.Days = Mjd_GPS->Days;
	Mjd_TT.FracDay = Mjd_GPS->FracDay + (TT_TAI-GPST_TAI)/SECPERDAY;
	Mjd_UTC.Days = Mjd_TT.Days;
	Mjd_UTC.FracDay = Mjd_TT.FracDay - TT_UTC() / SECPERDAY;
	InterposeEOP( &Mjd_UTC, &CurrEop );

	Mjd_UT1.Days = Mjd_UTC.Days;
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	PrecMatrix( MJD_J2000, &Mjd_TT, Prec );
	NutMatrixSimple( &Mjd_TT, Nut );
	MatrixMultiply( 3, 3, 3, 3, Nut, Prec, T );

	GHAMatrix( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, Mat );
	MatrixMultiply( 3, 3, 3, 3, Mat, T, E );

	for( i=0; i<3; i++ )
	{
		Pos[i] = Y0[i];
		Vel[0][i] = Y0[3+i];
	}

	AccelMain(TODO, &Mjd_TT, Pos, Vel[0], E, T, Para , Acc[0] );

	//   compute right function at the second time 

	Mjd_TT.FracDay += h/4.0;  
	Mjd_UTC.FracDay += h/4.0;
	InterposeEOP( &Mjd_UTC, &CurrEop );
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	GHAMatrix( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, Mat );
	MatrixMultiply( 3, 3, 3, 3, Mat, T, E );      //短时间内没有考虑岁差和章动变化

	for( i=0; i<3; i++ )
	{
		Pos[i]    = Y0[i]   + Vel[0][i] * step / 4.0;
		Vel[1][i] = Y0[3+i] + Acc[0][i] * step / 4.0;
	}
	AccelMain(TODO, &Mjd_TT, Pos, Vel[1], E, T, Para , Acc[1] );

	//  compute right function at the third time 

	Mjd_TT.FracDay += h/8.0;  
	Mjd_UTC.FracDay += h/8.0;
	InterposeEOP( &Mjd_UTC, &CurrEop );
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	GHAMatrix( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, Mat );
	MatrixMultiply( 3, 3, 3, 3, Mat, T, E );      // 短时间内没有考虑岁差和章动变化

	for( i=0; i<3; i++ )
	{
		Pos[i]    = Y0[i]   + (Vel[0][i]*3/32 + Vel[1][i]*9/32.0) * step;
		Vel[2][i] = Y0[3+i] + (Acc[0][i]*3/32 + Acc[1][i]*9/32.0) * step;
	}
	AccelMain(TODO, &Mjd_TT, Pos, Vel[2], E, T, Para , Acc[2] );

	//   compute right function at the fourth time 

	Mjd_TT.FracDay  += h*57/104.0;  
	Mjd_UTC.FracDay += h*57/104.0;
	InterposeEOP( &Mjd_UTC, &CurrEop );
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	GHAMatrix( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, Mat );
	MatrixMultiply( 3, 3, 3, 3, Mat, T, E );      // 短时间内没有考虑岁差和章动变化

	for( i=0; i<3; i++ )
	{
		Pos[i]    = Y0[i]   + (Vel[0][i]*1932/2197.0 - Vel[1][i]*7200/2197.0
			+ Vel[2][i]*7296/2197.0 )* step;
		Vel[3][i] = Y0[3+i] + (Acc[0][i]*1932/2197.0 - Acc[1][i]*7200/2197.0
			+ Acc[2][i]*7296/2197.0 )* step;
	}
	AccelMain(TODO, &Mjd_TT, Pos, Vel[3], E, T, Para , Acc[3] ); 

	//   compute right function at the fourth time 

	Mjd_TT.FracDay  += h/13.0;  
	Mjd_UTC.FracDay += h/13.0;
	InterposeEOP( &Mjd_UTC, &CurrEop );
	Mjd_UT1.FracDay = Mjd_UTC.FracDay + CurrEop.dUT1 / SECPERDAY;

	GHAMatrix( &Mjd_UT1, GH );
	PoleMatrix( CurrEop.x, CurrEop.y, Pole );
	MatrixMultiply( 3, 3, 3, 3, Pole, GH, Mat );
	MatrixMultiply( 3, 3, 3, 3, Mat, T, E );      // 短时间内没有考虑岁差和章动变化

	for( i=0; i<3; i++ )
	{
		Pos[i]    = Y0[i]   + (Vel[0][i]*439/216.0 - Vel[1][i]*8.0
			+ Vel[2][i]*3680/513.0 - Vel[3][i]*845/4104.0 )* step;
		Vel[4][i] = Y0[3+i] + (Acc[0][i]*439/216.0 - Acc[1][i]*8.0
			+ Acc[2][i]*3680/513.0 - Acc[3][i]*845/4104.0 )* step;
	}
	AccelMain(TODO, &Mjd_TT, Pos, Vel[4], E, T, Para , Acc[4] ); 

	//  积分结果  

	Mjd_GPS->FracDay += h;
	if( Mjd_GPS->FracDay >= 1.0 )    // 跨天检查 
	{
		Mjd_GPS->FracDay -= 1.0;
		Mjd_GPS->Days    += 1;
	}

	for( i=0; i<3; i++ )
	{
		Y0[i]   += ( Vel[0][i]*25/216 + Vel[2][i]*1408/2565
			+ Vel[3][i]*2197/4104 - Vel[4][i]/5 ) * step;
		Y0[i+3] += ( Acc[0][i]*25/216 + Acc[2][i]*1408/2565
			+ Acc[3][i]*2197/4104 - Acc[4][i]/5 ) * step;
	}
}

*/
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
/*
void OrbitIntegToGivenTime( MJDTIME* Mjd_GPS, MJDTIME* Mjd_GivenTime,
	double Y0[6], DYNMODELPARA* Para )
{
	short   sign;
	double  Step;              // 积分步长[s] 
	double  dT;                //下一时刻与当前时间的时间间隔【s】 

	dT = ( Mjd_GivenTime->Days - Mjd_GPS->Days 
		+ Mjd_GivenTime->FracDay - Mjd_GPS->FracDay ) * SECPERDAY;

	if( dT >= 0.0 )
	{
		sign = 1;
	}
	else
	{
		sign = -1;
	}

	do{
		if( fabs(dT) > 30.0 )  //积分步长为30s 
		{
			Step = 30.0 * sign;    //  设置积分步长 
		}
		else
		{
			Step = dT;
		}

		RKF4Step( Mjd_GPS, Step, Y0, Para );

		dT = dT - Step;

	}while( fabs(dT) > 1E-8 );

}
*/
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
//                 积分后为下一时刻的卫星位置和速度
//   Y0         转移矩阵一定为对角单位阵，但是前面还有6个不为0的估计参数，
//   IntState   用Hermit5方法内插所需要的状态参数，在此处初始化[位置/速度/加速度]
//   Para       动力学模型参数
//
***************************************************************************/
void RKF4OrbitSTM( int graceType, MJDTIME* Mjd_GPS, double step, double Y0[54], 
	SCState Stat[2], DYNMODELPARA* Para )
{
	int i;
	double Y[54], dY[5][54];
	double h = step / SECPERDAY;

	for( i=0; i<54; i++ )
	{
		Y[i] = Y0[i];
	}

	VarEquation(graceType, Mjd_GPS, Y, dY[0] , Para );     /* first time */

	/* 初始化起点的内插状态参数 Stat[0] */

	for( i=0; i<3; i++ )
	{
		Stat[0].Pos[i]   = Y0[i];
		Stat[0].Pos[i+3] = *(Y0+3+i);
		Stat[0].Acc[i]   = *(dY[0]+3+i);
	}
	/**********************************/

	Mjd_GPS->FracDay += h/4.0;
	for( i=0; i<54; i++ )
	{
		Y[i] = Y0[i] + dY[0][i] * step / 4.0;
	}
	VarEquation(graceType, Mjd_GPS, Y, dY[1] , Para );     /* 2nd time */


	Mjd_GPS->FracDay += h/8.0;
	for( i=0; i<54; i++ )
	{
		Y[i] = Y0[i] + (dY[0][i]*3/32 + dY[1][i]*9/32.0) * step;
	}
	VarEquation(graceType, Mjd_GPS, Y, dY[2] , Para );     /* 3rd time */


	Mjd_GPS->FracDay += h*57.0/104.0;
	for( i=0; i<54; i++ )
	{
		Y[i] = Y0[i] + (dY[0][i]*1932/2197.0 - dY[1][i]*7200/2197.0
			+dY[2][i]*7296/2197.0) * step;
	}
	VarEquation(graceType, Mjd_GPS, Y, dY[3] , Para );     /* 4th time */


	Mjd_GPS->FracDay += h/13.0;
	for( i=0; i<54; i++ )
	{
		Y[i] = Y0[i] + (dY[0][i]*439/216.0 - dY[1][i]*8.0
			+ dY[2][i]*3680/513.0 - dY[3][i]*845/4104.0) * step;
	}
	VarEquation(graceType, Mjd_GPS, Y, dY[4] , Para );     /* 5th time */
	


	/*  积分结果  */

	if( Mjd_GPS->FracDay >= 1.0 )    /*  跨天检查 */
	{
		Mjd_GPS->FracDay -= 1.0;
		Mjd_GPS->Days    += 1;
	}

	for( i=0; i<54; i++ )
	{
		Y0[i]   += ( dY[0][i]*25/216 + dY[2][i]*1408/2565
			+ dY[3][i]*2197/4104 - dY[4][i]/5 ) * step;//
	}

	/* 初始化终点的内插状态参数 Stat[1] */

	for( i=0; i<3; i++ )
	{
		Stat[1].Pos[i]   = Y0[i];
		Stat[1].Pos[i+3] = *(Y0+3+i);
		Stat[1].Acc[i]   = dY[0][i+3]*25/216 + dY[2][i+3]*1408/2565
			+ dY[3][i+3]*2197/4104 - dY[4][i+3]/5 ;//这里还未进行内插，只是计算内插t2时刻的状态，这里没有
	}
	/**********************************/

}


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

void Hermite5( const SCState* S0, const SCState* S1, SCState* CurrState )
{
	int i;
	double delta, step;   /* delta is coefficent, step is time space */
	double dr[6], dv[6];  /* dr dv is coefficent for pos and vel interpolation */

	step = (S1->Mjd_GPS.Week-S0->Mjd_GPS.Week)*SECPERWEEK +
		S1->Mjd_GPS.SecOfWeek - S0->Mjd_GPS.SecOfWeek;

	/*      0<delta<1   */
	delta = ((CurrState->Mjd_GPS.Week-S0->Mjd_GPS.Week)*SECPERWEEK +
		CurrState->Mjd_GPS.SecOfWeek-S0->Mjd_GPS.SecOfWeek) / step;

	dr[0] = 1.0 + ( -10 + 15*delta - 6*delta*delta ) * pow(delta, 3.0 );
	dr[1] = delta + ( -6 + 8*delta - 3*delta*delta ) * pow(delta, 3.0 );
	dr[2] = (1 - 3*delta + 3*delta*delta - pow(delta,3.0)) * pow(delta, 2.0)/2.0;
	dr[3] = 1.0 - dr[0];
	dr[4] = ( -4.0 + 7*delta -3*delta*delta ) * pow(delta, 3.0);
	dr[5] = ( 1.0 -2*delta + delta*delta ) *pow(delta,3.0)/2.0;

	dv[0] = 30.0 * ( -1 + 2*delta - delta*delta ) * pow(delta, 2.0 );
	dv[1] = 1.0 + ( -18 + 32*delta - 15*delta*delta ) * pow(delta, 2.0 );
	dv[2] = delta + (-9 + 12*delta -5*delta*delta ) * pow(delta, 2.0) / 2.0;
	dv[3] = - dv[0];
	dv[4] = ( -12.0 + 28*delta -15*delta*delta ) * pow(delta, 2.0);
	dv[5] = ( 3.0 -8*delta + 5*delta*delta ) *pow(delta, 2.0) / 2.0;

	for( i=0; i<3; i++ )
	{
		CurrState->Pos[i] = dr[0]*S0->Pos[i] + dr[1]*step*S0->Pos[i+3] + dr[2]*step*step*S0->Acc[i]
		+ dr[3]*S1->Pos[i] + dr[4]*step*S1->Pos[i+3] + dr[5]*step*step*S1->Acc[i];
		CurrState->Pos[i+3] = dv[0]*S0->Pos[i]/step + dv[1]*S0->Pos[i+3] + dv[2]*step*S0->Acc[i]
		+ dv[3]*S1->Pos[i]/step + dv[4]*S1->Pos[i+3] + dv[5]*step*S1->Acc[i];
	}

}


/***************************************************************************
//
// InitStateTranMatrix
//
// Purpose:
//
//   初始化轨道积分向量中的状态转移矩阵, 初值为单位阵I
//
// Input/Output:
//
//   row    状态转移矩阵的行数
col    状态转移矩阵的列数
STM    状态转移矩阵向量

****************************************************************************/

void InitStateTranMatrix( int row, int col, double STM[] )
{
	int j, k;

	for( j=0; j<row; j++ )
	{
		for( k=0; k<col; k++ )
		{
			if(j==k)  
			{
				*(STM+j*8+k) = 1.0;
			}
			else
			{
				*(STM+j*8+k) = 0.0;
			}
		}
	}

}

