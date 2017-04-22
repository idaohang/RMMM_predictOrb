/****************************************************************************
目的：    考虑动力学补偿算法, 用动力学轨道积分预报轨道作为先验轨道, 用伪距观测
数据, 推广卡尔曼滤波更新卫星的运行轨道参数.

说明:     EKFSTATE中的卫星位置均以GNSS接收机天线相位中心为参考, 并进行推广卡尔曼
滤波. 轨道积分以卫星的质量中心为参考.

编写时间：2008.12.19
版本:     V1.1
版权：    武汉大学
****************************************************************************/

#include <math.h>
#include<string.h>
#include <stdlib.h>
#include "GPSTime.h"
#include "EKFilter.h"
#include "RefSys.h"
#include "RTOD_Const.h"
#include "CommonFuncs.h"
#include "ReadObs.h"
#include "RTOD.h"
#include "OrbitInteg.h"
#include "ReadLevel1B.h"


double diffTime(GPSTIME  Time_A,GPSTIME  Time_B)/*求两个历元的时间差diff=Time_A-TimeB*/
{
	double diff=0.0;
	diff=(Time_A.Week-Time_B.Week)*SECPERWEEK+
				Time_A.SecOfWeek-Time_B.SecOfWeek+1E-5;
	return diff;
}
void alignEpoch( FILE * Fobs_A,OBSTYPELIST * ObsTypeList_A,EpochObsData * Epoch_A,FILE * Fobs_B,OBSTYPELIST * ObsTypeList_B,EpochObsData * Epoch_B)
{
	while(fabs(diffTime(Epoch_A->Time,Epoch_B->Time))>1.0)
	{
		if(diffTime(Epoch_A->Time,Epoch_B->Time)>0.0)
			ReadEpochObsQC( Fobs_B, ObsTypeList_B, Epoch_B );
		else
			ReadEpochObsQC( Fobs_A, ObsTypeList_A, Epoch_A );
	}
			
}
//同时初始化主辅卫星
int init2Sat(FILE * Fobs_A,OBSTYPELIST * ObsTypeList_A,EpochObsData * Epoch_A,GPSEPHREC* GPSEph_A,SATMIDRESULT * SatMidInfo_A,
	FILE * Fobs_B,OBSTYPELIST * ObsTypeList_B,EpochObsData * Epoch_B,GPSEPHREC* GPSEph_B,SATMIDRESULT * SatMidInfo_B,
	double OiState[108],EKFSTATE * EkfState)
{
	double  BLH_A[3],BLH_B[3];             /* 获得接收机的高程，若小于300Km，不进行滤波处理 */
	int i,j,SatNumUsed_A,SatNumUsed_B;
	PPRESULT  Result[2];
	EpochObsData Data_A[3],Data_B[3];             /*用于生成Doppler观测值 */
	EmptyPPRESULTStruct( &Result[0]);
	EmptyPPRESULTStruct( &Result[1]);
	i = 0;j=0;
	do{
		if( ReadEpochObs( Fobs_A, ObsTypeList_A, Epoch_A )&&ReadEpochObs( Fobs_B, ObsTypeList_B, Epoch_B ))
		{
			//主辅卫星时间对准
			alignEpoch(Fobs_A, ObsTypeList_A, Epoch_A,Fobs_B, ObsTypeList_B, Epoch_B);
			if( i > 2 )
			{
				memcpy( &Data_A[0], &Data_A[1], sizeof(EpochObsData) );
				memcpy( &Data_A[1], &Data_A[2], sizeof(EpochObsData) );
				memcpy( &Data_A[2], Epoch_A, sizeof(EpochObsData) );

				memcpy( &Data_B[0], &Data_B[1], sizeof(EpochObsData) );
				memcpy( &Data_B[1], &Data_B[2], sizeof(EpochObsData) );
				memcpy( &Data_B[2], Epoch_B, sizeof(EpochObsData) );

				if( CreateDopplerObs( Data_A )&&CreateDopplerObs( Data_B ) )
				{
					SatNumUsed_A = PointPositionRTOD(0, GPSEph_A,&EkfState->IonPara, //计算GPS和GLONASS的钟差，未计算钟速
							&Data_A[1], &EkfState->AprioriState[0], SatMidInfo_A , &Result[0] );//Record为GPS星历结构体,这里还未读取，在计算卫星位置钟差的子函数里才读取，好处是储存空间少，只需要保存 最新的。
					SatNumUsed_B = PointPositionRTOD(1, GPSEph_B,
							&EkfState->IonPara, &Data_B[1], &EkfState->AprioriState[1], SatMidInfo_B , &Result[1] );
							
							
					PointPositionVelocityDetermination( &Data_A[1], SatMidInfo_A, &Result[0] );
					PointPositionVelocityDetermination( &Data_B[1], SatMidInfo_B, &Result[1] );
					XYZToBLH( Result[0].Position, BLH_A, R_WGS84, F_WGS84 );//
					XYZToBLH( Result[1].Position, BLH_B, R_WGS84, F_WGS84 );
					printf("SatA:%12.3f ", Data_A[1].Time.SecOfWeek );
					MatrixPrintf( 1, 3, Result[0].Position );//结果和单个卫星的一样

					printf("SatB:%12.3f ", Data_B[1].Time.SecOfWeek );
					MatrixPrintf( 1, 3, Result[1].Position );
	
					if(fabs(Data_A[1].Time.SecOfWeek-Data_B[1].Time.SecOfWeek)>1.0)
					{
						printf("SatA:%12.3f\nSatB:%12.3f ", Data_A[1].Time.SecOfWeek, Data_B[1].Time.SecOfWeek  );
						continue;
					}
					//A卫星的初始化
					if( fmod( Data_A[1].Time.SecOfWeek+0.1, 30.0)<0.2 && Result[0].SatNumUsed>=6 &&
						Result[0].Coverage<1E-6 && Result[0].SigmaPos<10.0 && Result[0].SigmaVel<0.1 
						&& BLH_A[2] > 300000.0 )
					{
						InitAStateAndCov( &Data_A[1].Time,EkfState->Para[0].Cd, EkfState->Para[0].Cr, EkfState->Tao[0], &Result[0], EkfState);
						printf("%12.3f ", Data_A[1].Time.SecOfWeek );
						MatrixPrintf( 1, ADIMENSION, EkfState->StateA );

						CopyArray( 6, &OiState[0], EkfState->StateA );
						PhaseCentToMassCent( true, &EkfState->CentBias[0], &OiState[0] );

						EkfState->IntSate[0].Mjd_GPS = EkfState->Time;
						EkfState->IntSate[1].Mjd_GPS = EkfState->Time;
						EkfState->IntSate[1].Mjd_GPS.SecOfWeek = EkfState->IntSate[1].Mjd_GPS.SecOfWeek +EkfState->Step;
						CheckGPSTime( &EkfState->IntSate[1].Mjd_GPS );
						//这里为什么要进行时间更新
						EkfState->IsInitial[0]= 1;//初始化成功由未初始化0转为初始化1
								
					}
					//B
					if( fmod( Data_B[1].Time.SecOfWeek+0.1, 30.0)<0.2 && Result[1].SatNumUsed>=6 &&
						Result[1].Coverage<1E-6 && Result[1].SigmaPos<10.0 && Result[1].SigmaVel<0.1 
						&& BLH_B[2] > 300000.0 )
					{
						//初始化相对位置
						InitRelStateAndCov( &Data_B[1].Time,Result[1], EkfState);
						printf("%12.3f ", Data_B[1].Time.SecOfWeek );
						MatrixPrintf( 1, RELDIMENSION, &EkfState->StateRel[RELDIMENSION]);
						/*A星+BA相对位置与速度*/
						getBpos(EkfState->StateA,EkfState->StateRel,&OiState[54]);
						PhaseCentToMassCent( true, &EkfState->CentBias[3], &OiState[54] );
					

						EkfState->IntSate[2].Mjd_GPS = EkfState->Time;
						EkfState->IntSate[3].Mjd_GPS = EkfState->Time;
						EkfState->IntSate[3].Mjd_GPS.SecOfWeek = EkfState->IntSate[3].Mjd_GPS.SecOfWeek +EkfState->Step;
						CheckGPSTime( &EkfState->IntSate[3].Mjd_GPS );

						EkfState->IsInitial[1]= 1;//初始化成功由未初始化0转为初始化1
								
					}
					if(EkfState->IsInitial[0]== 1 && EkfState->IsInitial[1]== 1)
					{
						ABTimeUpdate(OiState,EkfState);
						//initStateAndCov中模糊度协方差已清零，模糊度已清零
							EmptyCurComObsStruct(&EkfState->CurComObs);
						//当前历元共同观测值要清零
						return 1;
					}
					else{
						EkfState->IsInitial[0]=0;
						EkfState->IsInitial[1]=0;
					}
				}
			}
			else
			{
						
				memcpy( &Data_A[i], Epoch_A, sizeof( EpochObsData) );
				memcpy( &Data_B[i], Epoch_B, sizeof(EpochObsData) );
			}

			i++;
		}

	}while(1) ;
}

int intUpdate(FILE * Fout_A,EpochObsData * Epoch_A,
	FILE * Fout_B,
	double   OiState[108],EKFSTATE * EkfState)
{
	double h;
	FINALORBIT FinalOrb[2];            /* 输出的轨道平根数与瞬根值等参数 */
	MJDTIME  Mjd_GPS_A,Mjd_GPS_B;
	do{
		//A星
				CopyArray( 6, FinalOrb[0].ECI_Orb, EkfState->StateA );
				PhaseCentToMassCent( true,&EkfState->CentBias[0], FinalOrb[0].ECI_Orb );
				GPSTimeToMJDTime( &EkfState->Time, &Mjd_GPS_A );
				ICRF_ITRF_GPST( MJD_J2000, &EkfState->Time, true, FinalOrb[0].ECI_Orb, FinalOrb[0].ECF_Orb );//转到地固系
				fprintf( Fout_A, " %12.3f %12.3f %12.3f %12.3f %10.3f %10.3f %10.3f %12.3f %12.3f %12.3f\n",
					EkfState->Time.SecOfWeek, FinalOrb[0].ECF_Orb[0], FinalOrb[0].ECF_Orb[1], FinalOrb[0].ECF_Orb[2],		
					FinalOrb[0].ECF_Orb[3], FinalOrb[0].ECF_Orb[4], FinalOrb[0].ECF_Orb[5], 
					FinalOrb[0].ECF_Orb[0], FinalOrb[0].ECF_Orb[1], FinalOrb[0].ECF_Orb[2]);
				EkfState->Para[0].Cd = EkfState->StateA[8];
				EkfState->Para[0].Cr = EkfState->StateA[9];
				EkfState->Tao[0]     = EkfState->StateA[10];
				CopyArray( 6, OiState, FinalOrb[0].ECI_Orb );   
				InitStateTranMatrix( 6, 8, &OiState[6] );
				RKF4OrbitSTM(0, &Mjd_GPS_A, EkfState->Step, &OiState[0], &EkfState->IntSate[0] , &EkfState->Para[0] );//OiState包括前6个状态矢量，以及后面48转移矩阵
				EkfState->IntSate[0].Mjd_GPS = EkfState->IntSate[1].Mjd_GPS;
				EkfState->IntSate[1].Mjd_GPS.SecOfWeek = EkfState->IntSate[1].Mjd_GPS.SecOfWeek +EkfState->Step;
				CheckGPSTime( &EkfState->IntSate[1].Mjd_GPS );
				PhaseCentToMassCent( false, &EkfState->CentBias[0], &OiState[0] );

				//B星相位中心位置=A星相位中心位置-AB相位中心相对位置
				getBpos(EkfState->StateA,EkfState->StateRel,FinalOrb[1].ECI_Orb);
				PhaseCentToMassCent( true,&EkfState->CentBias[3], FinalOrb[1].ECI_Orb );
				GPSTimeToMJDTime( &EkfState->Time, &Mjd_GPS_B );
				ICRF_ITRF_GPST( MJD_J2000, &EkfState->Time, true, FinalOrb[1].ECI_Orb, FinalOrb[1].ECF_Orb );//转到地固系
				fprintf( Fout_B, " %12.3f %12.3f %12.3f %12.3f %10.3f %10.3f %10.3f %12.3f %12.3f %12.3f\n",
					EkfState->Time.SecOfWeek, FinalOrb[1].ECF_Orb[0], FinalOrb[1].ECF_Orb[1], FinalOrb[1].ECF_Orb[2],		
					FinalOrb[1].ECF_Orb[3], FinalOrb[1].ECF_Orb[4], FinalOrb[1].ECF_Orb[5], 
					FinalOrb[1].ECF_Orb[0], FinalOrb[1].ECF_Orb[1], FinalOrb[1].ECF_Orb[2]);
				EkfState->Para[1].Cd = EkfState->StateRel[7];
				EkfState->Para[1].Cr = EkfState->StateRel[8];
				CopyArray( 6, &OiState[54], FinalOrb[1].ECI_Orb );   
				InitStateTranMatrix( 6, 8, &OiState[60] );
				RKF4OrbitSTM(1, &Mjd_GPS_B, EkfState->Step, &OiState[54], &EkfState->IntSate[2] , &EkfState->Para[1] );//OiState包括前6个状态矢量，以及后面48转移矩阵
				EkfState->IntSate[2].Mjd_GPS = EkfState->IntSate[3].Mjd_GPS;
				EkfState->IntSate[3].Mjd_GPS.SecOfWeek = EkfState->IntSate[3].Mjd_GPS.SecOfWeek +EkfState->Step;
				CheckGPSTime( &EkfState->IntSate[3].Mjd_GPS );
				PhaseCentToMassCent( false, &EkfState->CentBias[3], &OiState[54] );

				ABTimeUpdate(OiState,EkfState);
				EkfState->Time.SecOfWeek = EkfState->Time.SecOfWeek + EkfState->Step;//时间只加一次Step
				CheckGPSTime( &EkfState->Time );

				h = (Epoch_A->Time.Week-EkfState->Time.Week)*SECPERWEEK+Epoch_A->Time.SecOfWeek-EkfState->Time.SecOfWeek;
			}while( h+1E-5>=0.0 );
	return 1;
			// 			time( &endt );    			
			// 			fprintf( ORB, "%5d%12.1f %10.3lf 1\n", Epoch.Time.Week, Epoch.Time.SecOfWeek, difftime( endt, start ) );
}
void EmptyOneSat11Obs(OneSat11Obs * oneSat11Obs)
{
	
		oneSat11Obs->CA1=0.0;
		oneSat11Obs->CA2=0.0;
		oneSat11Obs->CAApplyEpoch=0;
		oneSat11Obs->elevation[0]=0.0;
		oneSat11Obs->elevation[0]=0.0;//1
		//oneSat11Obs->map[0]=0.0;
		//oneSat11Obs->map[1]=0.0;
		oneSat11Obs->L11=0.0;
		oneSat11Obs->L12=0.0;
		oneSat11Obs->L21=0.0;
		oneSat11Obs->L22=0.0;
		oneSat11Obs->L1ApplyEpoch=0;
		oneSat11Obs->L2ApplyEpoch=0;
		oneSat11Obs->P11=0.0;
		oneSat11Obs->P12=0.0;
		oneSat11Obs->P1ApplyEpoch=0;
		oneSat11Obs->P21=0.0;
		oneSat11Obs->P22=0.0;
		oneSat11Obs->P2ApplyEpoch=0;
		oneSat11Obs->PRN=999;
		
		oneSat11Obs->satPos[0]=0.0;
		oneSat11Obs->satPos[1]=0.0;
		oneSat11Obs->satPos[2]=0.0;
		/*
		oneSat11Obs->satpos.spx2=0.0;
		oneSat11Obs->satpos.spy2=0.0;
		oneSat11Obs->satpos.spz2=0.0;*/
		oneSat11Obs->S11=0.0;
		oneSat11Obs->S12=0.0;
		oneSat11Obs->S21=0.0;
		oneSat11Obs->S22=0.0;
		//单差初始化
		oneSat11Obs->dCA=0.0;
		oneSat11Obs->dL1=0.0;
		oneSat11Obs->dL2=0.0;
		oneSat11Obs->dP1=0.0;
		oneSat11Obs->dP2=0.0;
		oneSat11Obs->flag=0;
		oneSat11Obs->used=0;
		//周跳相关初始化
		oneSat11Obs->Nw_SD	=0.0;
		oneSat11Obs->sigNw_SD=0.5;
		oneSat11Obs->Nw_num_SD=1;
		oneSat11Obs->index_satlist=-999;
		//
}
void EmptyCurComObsStruct(Common11Obs* ComObs)
{
	int j;
	ComObs->ComSatNum=0;
	ComObs->Time.SecOfWeek=0.0;
	ComObs->Time.Week=0;
	for(j=0;j<MAXCHANNUM;j++){
		EmptyOneSat11Obs(&ComObs->comobs[j]);
	}
	
}
void EmptyddObsStruct(DDPSEUDOOBS *ddObs)
{
	int i;
	ddObs->DDObsNum=0;
	for(i=0;i<MAXCHANNUM;i++)
	{
		ddObs->satddobs[i].nonRefPrn=999;
		ddObs->satddobs[i].ddF1=0;
		ddObs->satddobs[i].ddF2=0;
		/*
		ddObs->satddobs[i].satpos.spx1=0.0;
		ddObs->satddobs[i].satpos.spx2=0.0;
		ddObs->satddobs[i].satpos.spy1=0.0;
		ddObs->satddobs[i].satpos.spy2=0.0;
		ddObs->satddobs[i].satpos.spz1=0.0;
		ddObs->satddobs[i].satpos.spz2=0.0;
		*/
		
	}
}
void EmptyComObssStruct(Common11Obs * ComObss)
{
	int i;
	for(i=0;i<MaxValiEpochNum;i++)
	{
		ComObss[i].ComSatNum=0;
		ComObss[i].Time.SecOfWeek=0.0;
		ComObss[i].Time.Week=0;
		EmptyCurComObsStruct(&ComObss[i]);
	}
}
/***************************************************************************
//
// EmptyEKFSTATEStruct
//
// 目的: 清空EKFSTATE结构体
//
EkfState 推广卡尔曼滤波的状态参数, 协方差信息

***************************************************************************/

void EmptyEKFSTATEStruct( EKFSTATE* EkfState )
{
	int  i, j;

	EkfState->Time.Week      = 0;
	EkfState->Time.SecOfWeek = 0.0;
	/************************************************************************/
	/* 参考星的状态和协方差初始化                                           */
	/************************************************************************/
	for (i=0;i<ADIMENSION;i++)
	{
		if( i<6 )
		{
			EkfState->StateInECEFA[i] = 0.0;//A星
			EkfState->StateInECEFB[i] = 0.0;//B星
		}
		EkfState->StateA[i]=0.0;
		for( j=0; j<ADIMENSION; j++ )
		{
			EkfState->CovA[i*ADIMENSION+j] = 0.0;
		}
	}
	/************************************************************************/
	/* 相对状态和协方差初始化                                               */
	/************************************************************************/
	for (i=0;i<RELDIMENSION;i++)
	{
		EkfState->StateRel[i]=0.0;
		for( j=0; j<RELDIMENSION; j++ )
		{
			EkfState->CovRel[i*RELDIMENSION+j] = 0.0;
		}
	}
	for(int m=0;m<2*MAXCHANNUM;m++)
		EkfState->PRN[m]=0;
	EkfState->ApriSigma[0] = 999.999;
	EkfState->ApriSigma[1] = 999.999;
	EkfState->PostSigma[0] = 999.999;
	EkfState->PostSigma[1] = 999.999;
	EkfState->SatNumUsed[0]    = 0;
	EkfState->SatNumUsed[1]    = 0;
	EkfState->IsInitial[0] = 0;
	EkfState->IsInitial[1] = 0;
	EkfState->KFConvergency[0] = 0;
	EkfState->KFConvergency[1] = 0;
	EkfState->sigPosLC=0.0;
	EkfState->sigPosPC=0.0;
	//载波数据
	
}
/***************************************************************************
//
// InitEKFStateAndCov
//
// 目的: 用单点定位和测速的结果, 初始化EKF的状态和协方差参数, 在滤波重置时调用.
//
// 输入参数:
//
//  Time     卫星轨道参数对应的时间, 用GPS时间表示
Cd       大气阻力系数
Cr       光压系数
Tau      一阶高斯马尔可夫模型的相关时间参数
//  Result   单点定位结果, 包括位置/速度/接收机钟差及其变化率等

输出参数:

//  EkfState 推广卡尔曼滤波的状态参数, 协方差信息
//graceType=0表示A卫星，=1表示B卫星
***************************************************************************/

void InitAStateAndCov(const GPSTIME* Time, double Cd, double Cr, double Tau, 
	PPRESULT* Result, EKFSTATE* EkfState )
{
	int i, j;
	double ECEF[6];
	EkfState->Time.Week      = Time->Week;
	EkfState->Time.SecOfWeek = Time->SecOfWeek;
	//矩阵是一行一行存储
	for( i=0; i<ADIMENSION; i++ )
		{
			EkfState->StateA[i]=0.0;
			for( j=0; j<ADIMENSION; j++ )
			{
				EkfState->CovA[i*ADIMENSION+ADIMENSION+j]=0.0;
			}
	}
	EkfState->KFConvergency[0] = 0;

	for( i=0; i<3; i++ )
	{
		ECEF[i] = Result->Position[i];
		ECEF[i+3] = Result->Velocity[i];
	}

	ICRF_ITRF_GPST( MJD_J2000, Time, 0, EkfState->StateA, ECEF );

	for( i=0; i<3; i++ )
	{
		EkfState->StateA[ADIMENSION-3+i] = 0.0;         /* 补偿加速度  */

		EkfState->CovA[i*ADIMENSION+i] = 1E4;      /* 初始位置精度设置为100米 */
		EkfState->CovA[(i+3)*ADIMENSION+(i+3)] = 100.0;  /*初始速度精度设置为10m/s*/
		EkfState->CovA[(i+11)*ADIMENSION+(i+11)] = 0.01;   /*初始补偿加速度为0.1m/s^2 */
	}
	
	EkfState->IsInitial[0] = 1;                       /* 滤波器初始化 */
	EkfState->StateA[6] = Result->RcvClkOft[0];     /* GPS 钟差 */
	EkfState->StateA[7] = Result->RcvClkSft;        /* 钟差变化率  GPS、GLONASS变化率认为是一样的*/
	EkfState->StateA[8] = Cd;
	EkfState->StateA[9] = Cr;
	EkfState->StateA[10] = Tau;

	EkfState->CovA[6*ADIMENSION+6] = 1E4;  /*初始GPS接收机钟差精度100m */
	EkfState->CovA[7*ADIMENSION+7] = 100.0; /*初始GPS接收机钟速精度10m */
	EkfState->CovA[8*ADIMENSION+8] = 10;   /*初始Cd精度10 */
	EkfState->CovA[9*ADIMENSION+9] = 10;   /*初始Cr精度10 */
	EkfState->CovA[10*ADIMENSION+10]= 10;   /*初始Tau精度10 */
	//对于模糊度的初始化还未添加
}
/***************************************************************************
//
// InitEKFStateAndCov
//
// 目的: 用单点定位和测速的结果, 初始化EKF的状态和协方差参数, 在滤波重置时调用.
//
// 输入参数:
//
//  Time     卫星轨道参数对应的时间, 用GPS时间表示
Cd       大气阻力系数
Cr       光压系数
Tau      一阶高斯马尔可夫模型的相关时间参数
//  Result   单点定位结果, 包括位置/速度/接收机钟差及其变化率等

输出参数:

//  EkfState 推广卡尔曼滤波的状态参数, 协方差信息
//graceType=0表示A卫星，=1表示B卫星
***************************************************************************/

void InitRelStateAndCov(const GPSTIME* Time,
	PPRESULT Result, EKFSTATE* EkfState )
{
	int i, j;
	double ECEF[6];
	double ECI[6];
	
	/* 首先清空状态协方差矩阵 已根据不同卫星设置起始行列数*/
	for( i=0; i<RELDIMENSION; i++ )
		{
			EkfState->StateRel[i]=0.0;
			for( j=0; j<RELDIMENSION; j++ )
			{
				EkfState->CovRel[i*RELDIMENSION+j]=0.0;
			}
	}

	EkfState->KFConvergency[1] = 0;
		
	for( i=0; i<3; i++ )
	{
		ECEF[i]		= Result.Position[i] ;
		ECEF[i+3]	= Result.Velocity[i] ;
	}
	ICRF_ITRF_GPST( MJD_J2000, Time, 0, ECI,ECEF);//因为滤波是ECI系下
	getRelPos(EkfState->StateA,ECI,EkfState->StateRel);
	
	//变量顺序为x/y/z/Vx/Vy/Vz/Clk/Cd/Cr/Wx/Wy/Wz/L1  /WL
	//变量编号为0---2/ 3----5/ 6 / 7 /8 /9----11/12-23/24-35/
	for( i=0; i<3; i++ )
	{
		EkfState->StateRel[9+i] = 0.0;         /* 补偿加速度  */

		EkfState->CovRel[i*RELDIMENSION+i] = 1E4;      /* 初始位置精度设置为100米 */
		EkfState->CovRel[(i+3)*RELDIMENSION+(i+3)] = 100.0;  /*初始速度精度设置为10m/s*/
		EkfState->CovRel[(i+9)*RELDIMENSION+(i+9)] = 0.01;   /*初始补偿加速度为0.1m/s^2 */
	}
	
	EkfState->IsInitial[1] = 1;                       /* 滤波器初始化 */
	EkfState->StateRel[6]	= Result.RcvClkOft[0]-EkfState->StateA[6];     /* GPS 钟差 */
	EkfState->StateRel[7]	= EkfState->Para[1].Cd;
	EkfState->StateRel[8]	= EkfState->Para[1].Cr;

	EkfState->CovRel[6*RELDIMENSION+6]		= 1E4;  /*初始GPS接收机钟差精度100m */
	EkfState->CovRel[7*RELDIMENSION+7]		= 10;   /*初始Cd精度10 */
	EkfState->CovRel[8*RELDIMENSION+8]		= 10;   /*初始Cr精度10 */
	//对于模糊度的初始化还未添加
}


//补偿加速度在rtn方向的转移矩阵计算
void FormStateTransMatrixFromDynRTN4A( const double Tao, const double Step,
	double State[ADIMENSION], double STMDyn[48], double STM[] )
{
	int i, j;
	double ebdt, temp;
	double w[3],MatRTN2XYZ[9],drdwRTN[3],dvdwRTN[3],drdwXYZ[9],dvdwXYZ[9];
	for(i=0;i<3;i++)
		w[i]=State[11+i];//即补偿加速度

	/* 设置dr/dr0, dv/dv0, dr/dCd, dv/dCd  */
	for( i=0; i<6; i++ )
	{
		for( j=0; j<6; j++ )
		{
			STM[i*ADIMENSION+j] = STMDyn[i*8+j];
		}
		STM[i*ADIMENSION+8]  = STMDyn[i*8+6];           /* dr/dCd, dv/dCd  */
		STM[i*ADIMENSION+9] = STMDyn[i*8+7];           /* dr/dCr, dv/dCr  */
	}
	
	STM[8*ADIMENSION+8] = 1.0;          /* dCd/dCd  */
	STM[9*ADIMENSION+9] = 1.0;        /* dCr/dCr  */
	STM[10*ADIMENSION+10] = 1.0;        /* dTau/dTau */

	/* 设置dr/dw, dv/dw */

	ebdt =exp(-Step/Tao);//注意这里是三个方向相关时间一致，否则要用一个数组表示，但目前文献均是三个方向相关时间一致，调整的是功率谱
	temp = Step*Tao-(1.0-ebdt)*Tao*Tao;
	for(i=0;i<3;i++)
	{
		drdwRTN[i]=temp;
		dvdwRTN[i]=(1.0-ebdt)*Tao;
	}
	MatrixXYZ2RTN(State,MatRTN2XYZ,0);
	MatrixMultiply2(0,3,3,3,MatRTN2XYZ, drdwRTN, drdwXYZ );
	MatrixMultiply2(0,3,3,3,MatRTN2XYZ, dvdwRTN, dvdwXYZ );

	for( i=0; i<3; i++ )
		for(j=0;j<3;j++)
		{
			STM[i*ADIMENSION+11+j] = drdwXYZ[i*3+j];                 /* dr/dw */
			STM[(i+3)*ADIMENSION+11+j] = dvdwXYZ[i*3+j];				  /* dv/dw */
		}

	temp=Step*ebdt/Tao/Tao;                                  
	for( i=11; i<14; i++ )
	{
		STM[i*ADIMENSION+i] = ebdt;             /* dwdw */
		STM[i*ADIMENSION+10] = temp*w[i-11];    /* dw/dTau */
	}

	temp = (2.0*Tao*(ebdt-1)+Step*(ebdt+1));   
	for( i=0; i<3; i++ )
	{
		STM[i*ADIMENSION+10] = temp * w[i];    /* dr/dTao */
		STM[(i+3)*ADIMENSION+10] = ( (1.0-ebdt) - Step*ebdt/Tao ) * w[i]; /* dv/dTao */
		//因为w[i]一直等于0，所以其他不能影响到tao，tao值不变
	}
	
	STM[6*ADIMENSION+6] = 1.0;								// dclk/dclk 
	STM[6*ADIMENSION+7] = Step;								// dclk/dclksft 
	STM[7*ADIMENSION+7] = 1.0;								// dclksft/dclksft 
}

void FormStateTransMatrixFromDynRTN4rel( const double Tao, const double Step,
	double State[RELDIMENSION], double STMDyn[48], double STM[] )
{
	int i, j;
	double ebdt, temp;
	double w[3],MatRTN2XYZ[9],drdwRTN[3],dvdwRTN[3],drdwXYZ[9],dvdwXYZ[9];
	for(i=0;i<3;i++)
		w[i]=State[9+i];//即补偿加速度

	/* 设置dr/dr0, dv/dv0, dr/dCd, dv/dCd  */
	for( i=0; i<6; i++ )
	{
		for( j=0; j<6; j++ )
		{
			STM[i*RELDIMENSION+j] = STMDyn[i*8+j];
		}
		//因为前6个速度位置的积分没有带入该函数
		STM[i*RELDIMENSION+7] = STMDyn[i*8+6];           /* dr/dCd, dv/dCd  */
		STM[i*RELDIMENSION+8] = STMDyn[i*8+7];           /* dr/dCr, dv/dCr  */
	}
	STM[6*RELDIMENSION+6] = 1.0;				// dclk/dclk
	STM[7*RELDIMENSION+7] = 1.0;          /* dCd/dCd  */
	STM[8*RELDIMENSION+8] = 1.0;        /* dCr/dCr  */

	/* 设置dr/dw, dv/dw */

	ebdt =exp(-Step/Tao);//注意这里是三个方向相关时间一致，否则要用一个数组表示，但目前文献均是三个方向相关时间一致，调整的是功率谱
	temp = Step*Tao-(1.0-ebdt)*Tao*Tao;
	for(i=0;i<3;i++)
	{
		drdwRTN[i]=temp;
		dvdwRTN[i]=(1.0-ebdt)*Tao;
	}
	MatrixXYZ2RTN(State,MatRTN2XYZ,0);
	MatrixMultiply2(0,3,3,3,MatRTN2XYZ, drdwRTN, drdwXYZ );
	MatrixMultiply2(0,3,3,3,MatRTN2XYZ, dvdwRTN, dvdwXYZ );

	for( i=0; i<3; i++ )
		for(j=0;j<3;j++)
		{
			STM[i*RELDIMENSION+9+j] = drdwXYZ[i*3+j];                 /* dr/dw */
			STM[(i+3)*RELDIMENSION+9+j] = dvdwXYZ[i*3+j];				  /* dv/dw */
		}

	temp=Step*ebdt/Tao/Tao;                                  
	for( i=9; i<12; i++ )
	{
		STM[i*RELDIMENSION+i] = ebdt;             /* dwdw */
	}
	for( i=12;i<36;i++)
	{
		STM[i*RELDIMENSION+i] = 1;			/*模糊度*/
	}
}
//RTN方向的补偿加速度
void DMCrtn4A( const double Tao, const double Step,
	double State[])
{
	int i;
	double ebdt,temp;
	double MatXYZ2RTN[9]; 
	double deltRrtn[3],deltVrtn[3];
	double deltRxyz[3],deltVxyz[3];

	ebdt =exp(-Step/Tao);
	temp = Step*Tao-(1.0-ebdt)*Tao*Tao;
	for(i=0;i<3;i++)
	{
		deltRrtn[i]=temp*State[i+11];
		deltVrtn[i]=Tao*(1.0-ebdt)*State[i+11];
		deltRxyz[i]=0.0;
		deltVxyz[i]=0.0;
	}
	MatrixXYZ2RTN(State,MatXYZ2RTN,0);//生成的就是RTN2XYZ的转换矩阵
	for(i=0;i<3;i++)
	{
		//位置
		deltRxyz[0]+=MatXYZ2RTN[i]*deltRrtn[i];
		deltRxyz[1]+=MatXYZ2RTN[i+3]*deltRrtn[i];
		deltRxyz[2]+=MatXYZ2RTN[i+6]*deltRrtn[i];
		//速度
		deltVxyz[0]+=MatXYZ2RTN[i]*deltVrtn[i];
		deltVxyz[1]+=MatXYZ2RTN[i+3]*deltVrtn[i];
		deltVxyz[2]+=MatXYZ2RTN[i+6]*deltVrtn[i];
	}

	for( i=0; i<3; i++ )
	{
		State[i] = State[i] + deltRxyz[i];      /* 补偿加速度对位置的影响 */
		State[i+3] = State[i+3] + deltVxyz[i];  /* 补偿加速度对速度的影响 */
		State[i+11] = ebdt * State[i+11];         /* 补偿加速度自身的更新 */
	}

}
//RTN方向的补偿加速度
void DMCrtn4Rel( const double Tao, const double Step,
	double State[])
{
	int i;
	double ebdt,temp;
	double MatXYZ2RTN[9]; 
	double deltRrtn[3],deltVrtn[3];
	double deltRxyz[3],deltVxyz[3];

	ebdt =exp(-Step/Tao);
	temp = Step*Tao-(1.0-ebdt)*Tao*Tao;
	for(i=0;i<3;i++)
	{
		deltRrtn[i]=temp*State[i+9];
		deltVrtn[i]=Tao*(1.0-ebdt)*State[i+9];
		deltRxyz[i]=0.0;
		deltVxyz[i]=0.0;
	}
	MatrixXYZ2RTN(State,MatXYZ2RTN,0);//生成的就是RTN2XYZ的转换矩阵
	for(i=0;i<3;i++)
	{
		//位置
		deltRxyz[0]+=MatXYZ2RTN[i]*deltRrtn[i];
		deltRxyz[1]+=MatXYZ2RTN[i+3]*deltRrtn[i];
		deltRxyz[2]+=MatXYZ2RTN[i+6]*deltRrtn[i];
		//速度
		deltVxyz[0]+=MatXYZ2RTN[i]*deltVrtn[i];
		deltVxyz[1]+=MatXYZ2RTN[i+3]*deltVrtn[i];
		deltVxyz[2]+=MatXYZ2RTN[i+6]*deltVrtn[i];
	}

	for( i=0; i<3; i++ )
	{
		State[i] = State[i] + deltRxyz[i];      /* 补偿加速度对位置的影响 */
		State[i+3] = State[i+3] + deltVxyz[i];  /* 补偿加速度对速度的影响 */
		State[i+9] = ebdt * State[i+9];         /* 补偿加速度自身的更新 */
	}

}
//rtn方向补偿加速度时过程噪声的计算，上面是补偿加速度在xyz方向时的计算
void UdProcessNoiseCovRTN4A( double Tao, double Step, double Sigma,
	double State[6], double Qcov[] )
{
	int i, j;
	double Mrtn2xyz[9],Mxyz2rtn[9];//E为对角阵
	double qrr[3],qrv[3],qrw[3],qvv[3],qvw[3],qww[3];//rtn坐标系下对角噪声
	double Qrr[9],Qrv[9],Qrw[9],Qwr[9],Qvv[9],Qvw[9],Qwv[9];//转到xyz坐标系下则不是对角噪声了
	double ebdt, ebdt2;
	double sig2[3];//rtn三个方向可以不同，现在假定一样

	sig2[0]=0.01*Sigma*Sigma*Tao*Tao;
	sig2[1]=Sigma*Sigma*Tao*Tao;
	sig2[2]=0.09*Sigma*Sigma*Tao*Tao;


	ebdt  = exp(-Step/Tao);
	ebdt2 = ebdt * ebdt;
	//下面是力学模型补偿可见staticstical orbit determination P507
	
	//仅对角有数值
	for(i=0;i<3;i++){
		qrr[i]=(pow(Step,3.0)/3.0 - Step*Step*Tao + (1-ebdt2)*pow(Tao,3.0)/2.0 
		+ Step*(1-2*ebdt)*pow(Tao,2.0)) * sig2[i]; /*  Qrr  */
		qrv[i]=(Step*Step/2.0 + (1+ebdt2)*pow(Tao,2.0)/2 - Step*(1-ebdt)*Tao 
		- ebdt*pow(Tao,2.0)) * sig2[i];                      /*  Qrv Qvr */
		qrw[i]=(1 - ebdt2 -2*Step*ebdt/Tao) * sig2[i] * Tao / 2;   /* Qrw Qwr */  
		qvv[i]=(2*Step/Tao -3.0 + 4*ebdt -ebdt2)*sig2[i]*Tao/2;    /* Qvv  */  
		qvw[i]=(1-2*ebdt+ebdt2) * sig2[i] / 2;                     /* Qvw Qwv */ 
		qww[i]=(1-ebdt2)*sig2[i]/Tao/2.0;                          /* Qww  */  
	}
	MatrixXYZ2RTN(State,Mrtn2xyz,0);							//生成rtn转xyz矩阵
	MatrixXYZ2RTN(State,Mxyz2rtn,1);							//生成xyz转rtn矩阵
	//RTN转换到XYZ后,下列均是3x3矩阵，而Qww仍然表示为对角阵
	MatrixMultiply3( 3, 3, 3, 3, Mrtn2xyz, qrr,Mxyz2rtn, Qrr );
	MatrixMultiply3( 3, 3, 3, 3, Mrtn2xyz, qrv,Mxyz2rtn, Qrv );
	MatrixMultiply3( 3, 3, 3, 3, Mrtn2xyz, qvv,Mxyz2rtn, Qvv );
	MatrixMultiply2( 0, 3, 3, 3, Mrtn2xyz, qrw, Qrw );
	MatrixMultiply2( 0, 3, 3, 3, Mrtn2xyz, qvw, Qvw );
	MatrixMultiply2( 1, 3, 3, 3, Mxyz2rtn,qrw, Qwr );
	MatrixMultiply2( 1, 3, 3, 3, Mxyz2rtn,qvw, Qwv );

	for( i=0; i<ADIMENSION; i++ )
	{
		for( j=0; j<ADIMENSION; j++ )
		{
			Qcov[i*ADIMENSION+j] = 0.0;
		}
	}


	for( i=0; i<3; i++ )
	{	
		for(j=0;j<3;j++){
			//要考虑RTN不同方差的情况，以及转换矩阵后协方差不为0
		//变量顺序为x/y/z/Vx/Vy/Vz/Clk/sft/Cd/Cr/Tao/Wx/Wy/Wz/Vtec
		//变量编号为0/1/2/ 3/ 4/ 5/ 6 / 7 /8 /9 /10 /11/12/13/14
		Qcov[i*ADIMENSION+j] = Qrr[i*3+j];       /* Qrr */
		Qcov[i*ADIMENSION+j+3] = Qrv[i*3+j];     /* Qrv */
		Qcov[i*ADIMENSION+j+11] = Qrw[i*3+j];    /* Qrw */

		Qcov[(i+3)*ADIMENSION+j] = Qrv[i*3+j];    /* 因为转移矩阵相同Qvr=Qvr */
		Qcov[(i+3)*ADIMENSION+j+3] = Qvv[i*3+j];   /* Qvv */
		Qcov[(i+3)*ADIMENSION+j+11] = Qvw[i*3+j];   /* Qvw */

		Qcov[(i+11)*ADIMENSION+j] = Qwr[i*3+j];     /* Qwr */
		Qcov[(i+11)*ADIMENSION+j+3] = Qwv[i*3+j];   /* Qwv */
		}
		
		Qcov[(i+11)*ADIMENSION+i+11] = qww[i];   /* Qww */

	}
	// X Y Z Vx Vy Vz clk sft cd cr tao w1 w2 w3
	
	Qcov[6*ADIMENSION+6] = Sf*Step + Sg*pow(Step,3.0)/3.0;    //Qclk,clk for GPS 
	Qcov[6*ADIMENSION+7] = Sg*pow(Step,2.0)/2.0;              // Qclk,sft for GPS 

	Qcov[7*ADIMENSION+6] = Qcov[6*ADIMENSION+7];            // Qsft,clk for GPS 
	Qcov[7*ADIMENSION+7] = Sg*Step;                           // Qsft,sft  
	
	Qcov[8*ADIMENSION+8] = 1E-8*Step;                         /* QCd  1E-8*Step*/
	Qcov[9*ADIMENSION+9] = 1E-8*Step;                         /* QCr 1E-8*Step*/
	Qcov[10*ADIMENSION+10] = 1E-8*Step;                        /* QTao 1E-8*Step, Tao值随着时间逐渐增加或减小，可能引起滤波发散,设为0后，因为转移矩阵为1，不会改变 */
}
//rtn方向补偿加速度时过程噪声的计算，上面是补偿加速度在xyz方向时的计算
void UdProcessNoiseCovRTN4rel( double Tao, double Step, double Sigma,
	double State[6], double Qcov[] )
{
	int i, j;
	double Mrtn2xyz[9],Mxyz2rtn[9];//E为对角阵
	double qrr[3],qrv[3],qrw[3],qvv[3],qvw[3],qww[3];//rtn坐标系下对角噪声
	double Qrr[9],Qrv[9],Qrw[9],Qwr[9],Qvv[9],Qvw[9],Qwv[9];//转到xyz坐标系下则不是对角噪声了
	double ebdt, ebdt2;
	double sig2[3];//rtn三个方向可以不同，现在假定一样

	sig2[0]=0.01*Sigma*Sigma*Tao*Tao;
	sig2[1]=Sigma*Sigma*Tao*Tao;
	sig2[2]=0.09*Sigma*Sigma*Tao*Tao;


	ebdt  = exp(-Step/Tao);
	ebdt2 = ebdt * ebdt;
	//下面是力学模型补偿可见staticstical orbit determination P507
	
	//仅对角有数值
	for(i=0;i<3;i++){
		qrr[i]=(pow(Step,3.0)/3.0 - Step*Step*Tao + (1-ebdt2)*pow(Tao,3.0)/2.0 
		+ Step*(1-2*ebdt)*pow(Tao,2.0)) * sig2[i]; /*  Qrr  */
		qrv[i]=(Step*Step/2.0 + (1+ebdt2)*pow(Tao,2.0)/2 - Step*(1-ebdt)*Tao 
		- ebdt*pow(Tao,2.0)) * sig2[i];                      /*  Qrv Qvr */
		qrw[i]=(1 - ebdt2 -2*Step*ebdt/Tao) * sig2[i] * Tao / 2;   /* Qrw Qwr */  
		qvv[i]=(2*Step/Tao -3.0 + 4*ebdt -ebdt2)*sig2[i]*Tao/2;    /* Qvv  */  
		qvw[i]=(1-2*ebdt+ebdt2) * sig2[i] / 2;                     /* Qvw Qwv */ 
		qww[i]=(1-ebdt2)*sig2[i]/Tao/2.0;                          /* Qww  */  
	}
	MatrixXYZ2RTN(State,Mrtn2xyz,0);							//生成rtn转xyz矩阵
	MatrixXYZ2RTN(State,Mxyz2rtn,1);							//生成xyz转rtn矩阵
	//RTN转换到XYZ后,下列均是3x3矩阵，而Qww仍然表示为对角阵
	MatrixMultiply3( 3, 3, 3, 3, Mrtn2xyz, qrr,Mxyz2rtn, Qrr );
	MatrixMultiply3( 3, 3, 3, 3, Mrtn2xyz, qrv,Mxyz2rtn, Qrv );
	MatrixMultiply3( 3, 3, 3, 3, Mrtn2xyz, qvv,Mxyz2rtn, Qvv );
	MatrixMultiply2( 0, 3, 3, 3, Mrtn2xyz, qrw, Qrw );
	MatrixMultiply2( 0, 3, 3, 3, Mrtn2xyz, qvw, Qvw );
	MatrixMultiply2( 1, 3, 3, 3, Mxyz2rtn,qrw, Qwr );
	MatrixMultiply2( 1, 3, 3, 3, Mxyz2rtn,qvw, Qwv );

	for( i=0; i<RELDIMENSION; i++ )
	{
		for( j=0; j<RELDIMENSION; j++ )
		{
			Qcov[i*RELDIMENSION+j] = 0.0;
		}
	}


	for( i=0; i<3; i++ )
	{	
		for(j=0;j<3;j++){
			//要考虑RTN不同方差的情况，以及转换矩阵后协方差不为0
		//变量顺序为x/y/z/Vx/Vy/Vz/Clk/Cd/Cr/Wx/Wy/Wz/L1  /WL
		//变量编号为0---2/ 3----5/ 6 / 7 /8 /9----11/12-23/24-35/
		Qcov[i*RELDIMENSION+j] = Qrr[i*3+j];       /* Qrr */
		Qcov[i*RELDIMENSION+j+3] = Qrv[i*3+j];     /* Qrv */
		Qcov[i*RELDIMENSION+j+9] = Qrw[i*3+j];    /* Qrw */

		Qcov[(i+3)*RELDIMENSION+j] = Qrv[i*3+j];    /* 因为转移矩阵相同Qvr/Qvr应该互为转置？ */
		Qcov[(i+3)*RELDIMENSION+j+3] = Qvv[i*3+j];   /* Qvv */
		Qcov[(i+3)*RELDIMENSION+j+9] = Qvw[i*3+j];   /* Qvw */

		Qcov[(i+9)*RELDIMENSION+j] = Qwr[i*3+j];     /* Qwr */
		Qcov[(i+9)*RELDIMENSION+j+3] = Qwv[i*3+j];   /* Qwv */
		}
		
		Qcov[(i+9)*RELDIMENSION+i+9] = qww[i];   /* Qww */

	}
	//变量顺序为x/y/z/Vx/Vy/Vz/Clk/Cd/Cr/Wx/Wy/Wz/L1  /WL
	//变量编号为0---2/ 3----5/ 6 / 7 /8 /9----11/12-23/24-35/
	Qcov[6*RELDIMENSION+6] = t_nsd*Step;
	Qcov[7*RELDIMENSION+7] = 1E-8*Step;                         /* QCd  1E-8*Step*/
	Qcov[8*RELDIMENSION+8] = 1E-8*Step;                         /* QCr 1E-8*Step*/
	for (i=12;i<24;i++)
	{
		Qcov[i*RELDIMENSION+i]=1E-3;							/*QLc*/
	}
	for (i=24;i<36;i++)
	{
		Qcov[i*RELDIMENSION+i]=1E-4;							/*QMw*/
	}
}

//************************************
// 函数名称:  EKFTimeUpdateA
// 函数说明:  A星参考星时间更新
// 作者:      RJ
// 时间:   	  2016/12/01
// 返回值:    void
// 参数:      double OIState[54]
// 参数:      EKFSTATE * KFState
//************************************
void EKFTimeUpdateA(double OIState[54], EKFSTATE * KFState )
{
	int i,j;
	double STM[ADIMENSION*ADIMENSION] = {0.0};
	double Q[ADIMENSION*ADIMENSION] = {0.0};
	double TmpMat[ADIMENSION*ADIMENSION], STMT[ADIMENSION*ADIMENSION];
	double KFTmp[ADIMENSION*ADIMENSION];
	double Tao,Step,Sigma;
	Tao=KFState->Tao[0];
	Step=KFState->Step;
	Sigma=KFState->Sigma[0];
	for(i=0;i<ADIMENSION;i++)
		for(j=0;j<ADIMENSION;j++)
			KFTmp[i*ADIMENSION+j]=KFState->CovA[i*ADIMENSION+j];
	//状态更新为积分器中
	CopyArray( 6, KFState->StateA, OIState );
	/* P(-) = H * P(+) * H(T) + Q  */

	//rj20160721改为rtn方向动力学补偿噪声
	UdProcessNoiseCovRTN4A( Tao, Step, Sigma, KFState->StateA, Q);//A、B的过程噪声不相关，所以这里分开求
	FormStateTransMatrixFromDynRTN4A( Tao, Step, KFState->StateA, &OIState[6], STM );//A星检查一样，B星错误，OIState不对
	
	MatrixMultiply( ADIMENSION, ADIMENSION, ADIMENSION, ADIMENSION, 
		STM, KFTmp, TmpMat );
	
	MatrixTranspose( ADIMENSION, ADIMENSION, STM, STMT );
	MatrixMultiply( ADIMENSION, ADIMENSION, ADIMENSION, ADIMENSION,
		TmpMat, STMT, KFTmp );
	MatrixAddition2( ADIMENSION, ADIMENSION, Q, KFTmp );  
	//把协方差更新结果赋给滤波
	for(i=0;i<ADIMENSION;i++){
		for(j=0;j<ADIMENSION;j++){
			KFState->CovA[i*ADIMENSION+j]=KFTmp[i*ADIMENSION+j];
		}
	}
	/* 顾及估计的补偿加速度对卫星位置和速度的影响 */
	DMCrtn4A( Tao, Step, KFState->StateA);

}
/***************************************************************************
//
// EKFMeasureUpdate
//
// 目的: 使用UD滤波的原理, 将每个历元的多颗卫星的观测值组成序列, 依次对接收机的
状态进行测量更新.

说明： 测量更新后的卫星轨道是在惯性坐标系，以GPS接收机天线相位中心。

// 输入参数:

O_C      观测值减计算值[m]
R        观测值的精度指标[m]
H        观测值线性化的系数向量

输出参数

KFState  输入为经过时间更新后的卡尔曼滤波状态与协方差矩阵信息, 经过一颗GNSS
卫星的观测数据更新后, 输出更新后的状态和协方差矩阵信息

***************************************************************************/
bool EKFMeasureUpdateA( double O_C, double R, double H[], EKFSTATE* KFState )
{
	int i, j;
	double Error;
	double K[ADIMENSION]={0.0};		/* 卡尔曼滤波增益 */
	double TmpMat[ADIMENSION]={0.0};    
	double P_[ADIMENSION*ADIMENSION]={0.0};		/* 更新前的状态协方差矩阵 */
	double Mat[ADIMENSION*ADIMENSION]={0.0};
	for (i=0;i<ADIMENSION;i++)
	
	for( i=0; i<ADIMENSION; i++ )
    {
        for( j=0; j<ADIMENSION; j++ )
        {
            *(P_+i*ADIMENSION+j) = *(KFState->CovA+i*ADIMENSION+j);
        }
    }

	MatrixMultiply( ADIMENSION, ADIMENSION,ADIMENSION, 1,
		KFState->CovA, H, TmpMat );
	Error = R*R + VectDot( ADIMENSION, ADIMENSION, H, TmpMat );   /* HPH(T) + R  */
	
	if( fabs(O_C) > 3.0*Error)
	{
		return false;
	}

	/* 卡尔曼滤波增益计算 */

	for( i=0; i<ADIMENSION; i++ )
	{
		K[i] = TmpMat[i] / Error;//第一次更新A星一样
	}

	/* 卫星状态的测量更新 */

	for( i=0; i<ADIMENSION; i++ )
	{
		KFState->StateA[i] = KFState->StateA[i] + K[i] * O_C;//第一次更新A星一样，同一历元第二次更新方差全错
	}

	/* 状态协方差更新 */

	Dyadic( ADIMENSION, ADIMENSION, K, H, Mat );//实际上就是mx1矩阵乘以1xn矩阵，最后输出mxn矩阵

	for( i=0; i<ADIMENSION; i++ )
	{
		for( j=0; j<ADIMENSION; j++ )
		{
			*(Mat+i*ADIMENSION+j) = -1.0 * *(Mat+i*ADIMENSION+j);
			if( i==j )
			{
				*(Mat+i*ADIMENSION+j) = *(Mat+i*ADIMENSION+j) + 1.0;//因为1-KH实际上相当于E-KH
			}
		}
	}

	MatrixMultiply( ADIMENSION, ADIMENSION, ADIMENSION,ADIMENSION,
		Mat, P_, KFState->CovA );//这部出现错误

	if( KFState->StateA[10]<0.9 || KFState->StateA[10]>1.1 )   /* Tao  */
	{
		KFState->StateA[10] = 1.0;
	}
	return true;
}
bool EKFMeasureUpdateRel( double O_C, double R, double H[], EKFSTATE* KFState )
{
	int i, j;
	double Error;
	double K[RELDIMENSION]={0.0};
	double TmpMat[RELDIMENSION]={0.0};
	double P_[RELDIMENSION*RELDIMENSION]={0.0};
	double Mat[RELDIMENSION*RELDIMENSION]={0.0};
	
	for( i=0; i<RELDIMENSION; i++ )
    {
        for( j=0; j<RELDIMENSION; j++ )
        {
            *(P_+i*RELDIMENSION+j) = *(KFState->CovRel+i*RELDIMENSION+j);
        }
    }

	MatrixMultiply( RELDIMENSION, RELDIMENSION,RELDIMENSION, 1,
		KFState->CovRel, H, TmpMat );
	Error = R*R + VectDot( RELDIMENSION, RELDIMENSION, H, TmpMat );   /* HPH(T) + R  */
	
	if( fabs(O_C) > 3.0*Error)
	{
		/*
		//暂时不剔除粗差
		//return false;
		*/
	}

	/* 卡尔曼滤波增益计算 */

	for( i=0; i<RELDIMENSION; i++ )
	{
		K[i] = TmpMat[i] / Error;//第一次更新A星一样
	}

	/* 卫星状态的测量更新 */

	for( i=0; i<RELDIMENSION; i++ )
	{
		KFState->StateRel[i] = KFState->StateRel[i] + K[i] * O_C;//第一次更新A星一样，同一历元第二次更新方差全错
	}

	/* 状态协方差更新 */

	Dyadic( RELDIMENSION, RELDIMENSION, K, H, Mat );//实际上就是mx1矩阵乘以1xn矩阵，最后输出mxn矩阵

	for( i=0; i<RELDIMENSION; i++ )
	{
		for( j=0; j<RELDIMENSION; j++ )
		{
			*(Mat+i*RELDIMENSION+j) = -1.0 * *(Mat+i*RELDIMENSION+j);
			if( i==j )
			{
				*(Mat+i*RELDIMENSION+j) = *(Mat+i*RELDIMENSION+j) + 1.0;//因为1-KH实际上相当于E-KH
			}
		}
	}

	MatrixMultiply( RELDIMENSION, RELDIMENSION, RELDIMENSION,RELDIMENSION,
		Mat, P_, KFState->CovRel );//这部出现错误

	return true;
}
//总体测量更新
int WholeMeasureUpdate(FILE * Fout_A,EpochObsData * Epoch_A,GPSEPHREC* GPSEph_A,SATMIDRESULT * SatMidInfo_A,
	FILE * Fout_B,EpochObsData * Epoch_B,GPSEPHREC* GPSEph_B,SATMIDRESULT * SatMidInfo_B,
	double  OiState[108],EKFSTATE * KFState)
{
	int i,PCSatNum=0;
	double			TranMat[9];              /* 坐标转换矩阵       */
	PPRESULT		Result[2];
	PPRESULT		Result_LC;
	PPRESULT		Result_PC;   
	EOPPARA			CurrEop;
	FINALORBIT		FinalOrb[2];
	MJDTIME			Mjd_GPS_A;
	MJDTIME			Mjd_UT1;             /* 为计算格林尼治恒星时用 */
	Common11Obs		tempComObs;			 /*	当前历元共同卫星观测值*/
	int				flag;				 /* 是否所有模糊度初始化*/
	double			TranMatT[9];		/*地固转惯性系*/
	FILE * fout=fopen( "OutFile\\state.txt", "a+" );
	FILE * fcov=fopen( "OutFile\\cov.txt", "a+" );
	//--------------------------------------------------------------------------------
	// 转换到地固系
	ICRF_ITRF_Matrix( MJD_J2000, &Epoch_A->Time, true, TranMat );
	MatrixTranspose( 3, 3, TranMat, TranMatT);
	MatrixMultiply( 3, 3, 3, 1, TranMat, &OiState[0], KFState->AprioriState[0].LeoPos );
	MatrixMultiply( 3, 3, 3, 1, TranMat, &OiState[54], KFState->AprioriState[1].LeoPos );
	//--------------------------------------------------------------------------------

	//--------------------------------------------------------------------------------
	//初始状态精度
	if( KFState->KFConvergency[0] != 2 )
	{
		KFState->AprioriState[0].Validate = true;
		KFState->AprioriState[0].OrbitAccuracy = sqrt( KFState->CovA[0] + KFState->CovA[ADIMENSION+1] 
													 + KFState->CovA[2*ADIMENSION+2] );
	}
	else
	{
		KFState->AprioriState[0].Validate = 0;
		KFState->IsInitial[0] = 0;
	}
	if( KFState->KFConvergency[1] != 2 )
	{
		KFState->AprioriState[1].Validate = true;
		KFState->AprioriState[1].OrbitAccuracy = sqrt( KFState->CovRel[0] 
													 + KFState->CovRel[RELDIMENSION+1]
													 + KFState->CovRel[2*RELDIMENSION+2] );
	}
	else
	{
		KFState->AprioriState[1].Validate = 0;
		KFState->IsInitial[1] = 0;
	}
	//--------------------------------------------------------------------------------

	//--------------------------------------------------------------------------------
	//初始化几何定位结果
	EmptyPPRESULTStruct( &Result[0] );
	EmptyPPRESULTStruct( &Result[1] );
	EmptyPPRESULTStruct( &Result_LC );
	EmptyPPRESULTStruct( &Result_PC );
	//--------------------------------------------------------------------------------

	//------------------------------------------------------ --------------------------
	//单点定位
	PointPositionRTOD(0, GPSEph_A,&KFState->IonPara, Epoch_A, &KFState->AprioriState[0], SatMidInfo_A , &Result[0] );
	PointPositionRTOD(1, GPSEph_B,&KFState->IonPara, Epoch_B, &KFState->AprioriState[1], SatMidInfo_B , &Result[1] );
	/*
	if (!valResult(Result[0].SigmaPos,5,Result[0].SatNumUsed,4))
	{
		
		//Raim_spp();
	}
	if (!valResult(Result[1].SigmaPos,5,Result[1].SatNumUsed,4))
	{
		//Raim_spp();
	}
	*/
	double	Result1_ECI[3],Result2_ECI[3]; 
	MatrixMultiply( 3, 3, 3, 1, TranMatT, Result[0].Position, Result1_ECI );
	MatrixMultiply( 3, 3, 3, 1, TranMatT, Result[1].Position, Result2_ECI );
	EKFilterRTODA( TranMat, SatMidInfo_A, Epoch_A, &Result[0], KFState);						  /*A星伪距绝对测量更新*/
	//单点定位End------------------------------------------------------------------------	
	//--------------------------------------------------------------------------------
	//共同观测数据预处理
	EmptyCurComObsStruct( &tempComObs );						//先初始化
	GenCommomObs( Epoch_A, Epoch_B, &tempComObs );				//生成共同观测值				
	flag = CySlipDetection_LLI(TranMat, &KFState->CurComObs, &tempComObs , KFState );			/*探测周跳后初始化模糊度需要A、B星地固系初值*/
	//--------------------------------------------------------------------------------	
	//伪距滤波
	/*
	double RelPos_temp[3];
	PCSatNum=PPRTOD_PC( &KFState->CurComObs, KFState, &Result_PC, TranMat );
	if (PCSatNum>4)
	{
		MatrixMultiply(3,3,3,1,TranMatT,Result_PC.Position,RelPos_temp);
		for (i=0;i<3;i++)
		{
			KFState->StateRel[i]=RelPos_temp[i];
		}
		KFState->StateRel[6]=Result_PC.RcvClkOft[0];
	}*/
	relEKFilter_PC( TranMat, KFState );	
	//KFState->comSatNumUsed=PCSatNum;
	//载波滤波	
	if( flag )
	{
		//PPRTOD_LC(&KFState->CurComObs,KFState, &Result_LC, TranMat);
		//relEKFilter_LC( TranMat, KFState );
		
		//KFState->comSatNumUsed=KFState->SatNumUsed[1]+PCSatNum;
	}
	KFState->comSatNumUsed=KFState->CurComObs.ComSatNum;

	//将A星的测量更新结果代入积分器OiState[0]，使积分器的位置速度精度高一些，而不是靠初始化的位置速度一直积分
	CopyArray( 6, OiState, KFState->StateA );
	PhaseCentToMassCent( true, &KFState->CentBias[0], OiState );
	//将B星的测量更新结果代入积分器OiState[54]
	getBpos(KFState->StateA,KFState->StateRel,&OiState[54]);
	PhaseCentToMassCent( true, &KFState->CentBias[3], &OiState[54] );

	GPSTimeToMJDTime( &KFState->Time, &Mjd_GPS_A );

	/* 输出轨道 */
	//共同部分
	FinalOrb[0].Mjd_UTC = Mjd_GPS_A.Days + Mjd_GPS_A.FracDay - GetGPST_UTC()/SECPERDAY;
	Mjd_UT1 = Mjd_GPS_A;
	Mjd_UT1.FracDay = Mjd_UT1.FracDay - GetGPST_UTC()/SECPERDAY;
	MJDTimeToCommonTime( &Mjd_UT1, &FinalOrb[0].GT );
	MJDTimeToCommonTime( &Mjd_UT1, &FinalOrb[1].GT );
	FinalOrb[1].Mjd_UTC=FinalOrb[0].Mjd_UTC;
	//共同时间
	InterposeEOP( &Mjd_UT1, &CurrEop );
	Mjd_UT1.FracDay = Mjd_UT1.FracDay + CurrEop.dUT1/SECPERDAY;
	//A星
	CopyArray( 6, FinalOrb[0].ECI_Orb, OiState );
	ICRF_ITRF_GPST( MJD_J2000, &KFState->Time, true, OiState, FinalOrb[0].ECF_Orb );
	//B星
	CopyArray( 6, FinalOrb[1].ECI_Orb, &OiState[54] );
	ICRF_ITRF_GPST( MJD_J2000, &KFState->Time, true, &OiState[54], FinalOrb[1].ECF_Orb );//转到地固系

	fprintf(fout,"%12.3f",Epoch_A->Time.SecOfWeek);
	for(i=0;i<RELDIMENSION;i++)
		fprintf(fout,"%12.3f",KFState->StateRel[i]);
	fprintf(fout,"\n");
	fclose(fout);
	fprintf(fcov,"%12.3f",Epoch_A->Time.SecOfWeek);
	for(i=0;i<RELDIMENSION;i++)
		fprintf(fcov,"%12.3f",KFState->CovRel[i*RELDIMENSION+i]);
	fprintf(fcov,"\n");
	fclose(fcov);
	
	//A
	KFState->Para[0].Cd = KFState->StateA[8];
	KFState->Para[0].Cr = KFState->StateA[9];
	KFState->Tao[0]= KFState->StateA[10];
	//B
	KFState->Para[1].Cd = KFState->StateRel[7];
	KFState->Para[1].Cr = KFState->StateRel[8];
	return 1;
}
//************************************
// Method:    ABTimeUpdate
// FullName:  ABTimeUpdate
// Access:    public 
// Returns:   void
// Qualifier: 注意调用该函数之前，必须把测量更新的结果由相位中心转到质心
// Parameter: double OiState[108]
// Parameter: EKFSTATE * KFState
//************************************
void ABTimeUpdate(double  OiState[108],EKFSTATE* KFState)
{
	MJDTIME			Mjd_GPS_A;
	MJDTIME			Mjd_GPS_B;
	GPSTimeToMJDTime( &KFState->Time, &Mjd_GPS_A );
	Mjd_GPS_B.Days=Mjd_GPS_A.Days;
	Mjd_GPS_B.FracDay=Mjd_GPS_A.FracDay;
	//A
	InitStateTranMatrix( 6, 8, &OiState[6] );
	RKF4OrbitSTM(0, &Mjd_GPS_A, KFState->Step, &OiState[0], &KFState->IntSate[0], &KFState->Para[0]);
	KFState->IntSate[0].Mjd_GPS = KFState->IntSate[1].Mjd_GPS;
	KFState->IntSate[1].Mjd_GPS.SecOfWeek = KFState->IntSate[1].Mjd_GPS.SecOfWeek +KFState->Step;
	PhaseCentToMassCent( false, &KFState->CentBias[0], &OiState[0] );    /* 转换为相位中心 */
	//B

	InitStateTranMatrix( 6, 8, &OiState[60] );
	RKF4OrbitSTM(1, &Mjd_GPS_B, KFState->Step, &OiState[54], &KFState->IntSate[2] , &KFState->Para[1] );
	KFState->IntSate[2].Mjd_GPS = KFState->IntSate[3].Mjd_GPS;
	KFState->IntSate[3].Mjd_GPS.SecOfWeek = KFState->IntSate[3].Mjd_GPS.SecOfWeek +KFState->Step;
	PhaseCentToMassCent( false, &KFState->CentBias[3], &OiState[54] );    /* 转换为相位中心 */

	//要先时间更新卫星B，因为K时刻B星的位置要通过K时刻A星+K时刻相对位置，如果先更新A到K+1时刻，则不对
	EKFTimeUpdateRel(&OiState[54],KFState);
	EKFTimeUpdateA(&OiState[0],KFState);
	//最后更新后的相对状态又要通过更新后的A星与B星求差得到
	getRelPos(KFState->StateA,KFState->StateB,KFState->StateRel);//对相对状态速度和位置进更新
	KFState->Time.SecOfWeek = KFState->Time.SecOfWeek + KFState->Step;//时间只加一次Step
	CheckGPSTime( &KFState->Time );
	// 			time( &endt );    	
	printf("%5d%12.1f\n", KFState->Time.Week,KFState->Time.SecOfWeek);
}