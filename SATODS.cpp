/****************************************************************************
目的：    LEOSAOD软件的主函数

编写时间：2008.12.10
版本:     V1.1
版权：    武汉大学

修改记录：
2009.3.11    增加太阳光压摄动力，估计参数中增加Cr系数,
修改了长时间中断情况下的轨道输出。
2009.5.10    将重力场模型改为EGM2008，可改进重力场模型精度。

2009.8.3    如果使用单频数据进行实时定轨，需要修改以下两处代码
1. GetOneSatPseudoRange函数中，将P2改为其他类型中没有标记的符号
2. 将ComputeGPSSatOrbitAtSignalTrans中的Tgd加上，可以提高定轨精度。

****************************************************************************/
#pragma once
#include "Satods.h"
#include "RTOD.h"
#include "CommonFuncs.h"
FILE* FGPSEph;
FILE* FGLOEph;


void initFace()
{
	printf("\n        ********************************************************************\n");
	printf("        *                                                                  *\n");
	printf("        *                  LEO Spacecraft Spaceborne-GNSS                  *\n");
	printf("        *      Autonomous Orbit flying formation Determination Software    *\n");
	printf("        *                          LEOSAOD(DF)v1.1                         *\n");
	printf("        *                                                                  *\n");
	printf("        *                         Developed by WHU                         *\n");
	printf("        *                                                                  *\n");
	printf("        *                                   All Copyright reserved!        *\n");
	printf("        *                                                                  *\n");
	printf("        ********************************************************************\n\n\n");
}
void CheckDOY( int &year, int &doy )
{
	int days = 365;

	if( year%4 == 0 )
	{
		days++;
	}

	if( doy > days )
	{
		doy = doy - days;
		year++;
	}
}
void brushFile()
{
	//上一次计算刷新
	FILE* brush=fopen("OutFile\\graceA.txt","w+");
	fclose(brush);
	brush=fopen("OutFile\\graceB.txt","w+");
	fclose(brush);
	brush=fopen("OutFile\\line.txt","w+");
	fclose(brush);
	brush=fopen("OutFile\\state.txt","w+");
	fclose(brush);
	brush=fopen("OutFile\\cov.txt","w+");
	fclose(brush);
	brush=fopen("OutFile\\Out_obs.txt","w+");
	fclose(brush);
	brush=fopen("OutFile\\cycleSlip.txt","w+");
	fclose(brush);
	//上一次计算文件刷新结束
}
void InitPreOrbit(const int StartDOY,const int StartYear)
{
	int n_kbr;
	Lc_orbit=(LCPC_ORBIT*)malloc(17330*sizeof(struct LCPC_ORBIT));
	orbitB=(ORBIT*)malloc(17330*sizeof(struct ORBIT));
	orbitA=(ORBIT*)malloc(17330*sizeof(struct ORBIT));
	d_orbit=(D_ORBIT*)malloc(17330*sizeof(struct D_ORBIT));
	kbr=(KBR*)malloc(17330*sizeof(struct KBR));
	dOrbitA=(D_ORBIT*)malloc(17330*sizeof(struct D_ORBIT));
	dOrbitAB=(D_ORBIT*)malloc(17330*sizeof(struct D_ORBIT));
	dOrbitB=(D_ORBIT*)malloc(17330*sizeof(struct D_ORBIT));
	n_kbr=DataCheck(StartDOY,StartYear);//读取计算的参考值
	//基线坐标初始化
	for (int index=0;index<17330;index++)
		InitLc_orbit(Lc_orbit+index);
	//轨道坐标初始化
	for (int index=0;index<17330;index++)
		InitD_Orbit(dOrbitA+index);
	for (int index=0;index<17330;index++)
		InitD_Orbit(dOrbitB+index);
	for (int index=0;index<17330;index++)
		InitD_Orbit(dOrbitAB+index);
}
void freeAll()
{
	free(Lc_orbit);
	free(orbitB);
	free(orbitA);
	free(d_orbit);
	free(dOrbitA);
	free(dOrbitB);
	free(kbr);
	free(dOrbitAB);
}
void initFile(char * Satname_A, char * Satname_B, char * EphName, const int StartDOY, const int StartYear,
	char * ObsFileName_A, char * ObsFileName_B, char * GPSEphFileName, char * ResFileName_A, char * ResFileName_B,
	char * RsdFileName_A, char * RsdFileName_B, EKFSTATE KFState)
{
	//sprintf( ObsFileName_A, "InputFile\\%s%04d.%02do", Satname_A, StartDOY*10, StartYear );
	sprintf( ObsFileName_A, "InputFile\\%s%04d.%02do", Satname_A, StartDOY*10, StartYear );
	sprintf( ObsFileName_B, "InputFile\\%s%04d.%02do", Satname_B, StartDOY*10, StartYear );
	sprintf( GPSEphFileName, "InputFile\\%4s%04d.%02dn", EphName, StartDOY*10, StartYear );
	sprintf( ResFileName_A, "InputFile\\%s%04d.%02ds", Satname_A, StartDOY*10, StartYear );
	sprintf( ResFileName_B, "InputFile\\%s%04d.%02ds", Satname_B, StartDOY*10, StartYear );
	sprintf( RsdFileName_A, "OutFile\\graceA.txt_%d_%d",(int)KFState.Step,KFState.Para[0].m_a);
	sprintf( RsdFileName_B, "OutFile\\graceB.txt_%d_%d",(int)KFState.Step,KFState.Para[0].m_a);
}
void ReadAntFile(const char antexFile[],pcvs_t * pcvs)
{
	pcvs->n=0;
    pcvs->nmax=0;
    readantex(antexFile, pcvs);
}
int main( int args,char * argv[] )
{
	int  i,flag;
	time_t start, endt;
	double   transXYZV[6];//为了转换为RTN
	double   dXYZ[3];
	double   dRTN[3];
	/* 长时间实时定轨测试需要的文件目录等变量信息 */
	char ObsFileName_A[100],ObsFileName_B[100];
	char GPSEphFileName[100];
	char ResFileName_A[100],ResFileName_B[100],RsdFileName_A[100],RsdFileName_B[100];
	//char OrbFileName_A[100],OrbFileName_B[100];

	/* 卫星星历数据结构 */
	GPSEPHREC   Record_A[32],Record_B[32];

	/* 观测数据结构 */
	OBSTYPELIST ObsTypeList_A,ObsTypeList_B;
	EpochObsData Epoch_A,Epoch_B;
	PSPRange     PSPR_B[32],PSPR_A[32];    /* 用于存储相位平滑伪距的中间结果,适用于GPS */

	/* 单点定位数据结构 */
	SATMIDRESULT SatMidInfo_A[24],SatMidInfo_B[24];
	

	/* EKF定轨结果 */
	EKFSTATE KFState;
	
	/* 动力学轨道积分参数 */
	double   h;             
	SCState  CurrStat[2];              /* 内插轨道需要的参数 */
	double	 OiState[108]={0.0};

	/* 最终轨道输出 */
	char Satname_A[5], Satname_B[5],EphName[5];
	int  StartYear, StartDOY, EndYear, EndDOY;
	FINALORBIT FinalOrb[2];            /* 输出的轨道平根数与瞬根值等参数 */

	/*评价部分结构*/
	int numOfEpoch=0;
	int numOfA=0,numOfB=0;//当前A、B的GFZ真值索引

	double xA_pre[6]={0.0},xB_pre[6]={0.0};//前一个真值坐标的为了Kalman滤波坐标更新用
	initFace();

	printf("\nPlease input EPH abbr name (four alphabet: brdc)...\n" );
	//scanf("%s", EphName );
	sprintf(EphName,"brdc");
	printf("\nPlease input the year and DOY for first obs data...\n" );
	StartYear=10;StartDOY=152;//StartDOY=152;
	EndYear=10;EndDOY=152;//EndDOY=161;
	
	//刷新结果文件,初始化结构体
	brushFile();
	initKFState(&KFState);
	EmptySATMIDRESULTStruct( 24, SatMidInfo_A);
	EmptySATMIDRESULTStruct( 24, SatMidInfo_B);
	EmptyEpochObsDataStruct( &Epoch_A);
	EmptyEpochObsDataStruct( &Epoch_B);
	if(args!=3)
	{
		return 0;
	}
	else
	{
		KFState.Step = str2num(argv[1],0,4);
		KFState.Para[0].m_a = str2num(argv[2],0,4);
		KFState.Para[0].n_a = str2num(argv[2],0,4);
		KFState.Para[1].m_a = str2num(argv[2],0,4);
		KFState.Para[1].n_a = str2num(argv[2],0,4);
		//KFState.wuchaMod = str2num(argv[2],0,4);
	}

	time( &start );
	
	sprintf(Satname_A,"GPSA");
	sprintf(Satname_B,"GPSB");

	initFile(Satname_A, Satname_B, EphName, StartDOY, StartYear, ObsFileName_A, ObsFileName_B,GPSEphFileName,
		ResFileName_A,ResFileName_B,RsdFileName_A,RsdFileName_B, KFState);
	
	/************************************************************************/
	/* 读取文件*/
	/************************************************************************/
	FILE* FRESIDUAL_A= fopen( RsdFileName_A, "wt" );
	FILE* FRESIDUAL_B= fopen( RsdFileName_B, "wt" );
	FILE* Fout_A  = fopen( ResFileName_A, "a+" );
	FILE* Fout_B  = fopen( ResFileName_B, "a+" );
	/*
	FILE* Fobs_A= fopen( ObsFileName_A, "rt" );  
	FILE* Fobs_B= fopen( ObsFileName_B, "rt" );  
	*/
	
	/*读取精密星历文件,GPS天线相位中心改正文件*/
	/*
	PRE_EPH=(ONEPOCHOFSP3*)malloc(96*3*sizeof(ONEPOCHOFSP3));
	ReadSp3File("InputFile\\igs15861.sp3",PRE_EPH);
	ReadSp3File("InputFile\\igs15862.sp3",PRE_EPH);
	ReadSp3File("InputFile\\igs15863.sp3",PRE_EPH);
	//读取天线相位中心文件
	ReadAntFile("InputFile\\igs05_1627.atx",&pcvs);		//仅有GPS的
	//ReadAntFile("InputFile\\igs08_1758.atx",&pcvs);     //含有beidou等卫星的，读取需要额外判断
	
	if( (FGPSEph = fopen( GPSEphFileName, "rt" )) == NULL )
	{
		printf( "Cannot open %s GPS ephem file. \n", GPSEphFileName );
		return 0;
	}
	
	if( ReadObsHead( Fobs_A, &ObsTypeList_A ) ==0 )
	{
		printf( "Some problems in Reading the header of Obs file.\n" );
		return 0;
	}
	if( ReadObsHead( Fobs_B, &ObsTypeList_B ) ==0 )
	{
		printf( "Some problems in Reading the header of Obs file.\n" );
		return 0;
	}
	*/
	//精密低轨卫星轨道初始化
	int initFlag=0;
	double xA_pre_ECI[6],xB_pre_ECI[6];
	numOfEpoch=0;
	numOfA=0,numOfB=0;//当前A、B的GFZ真值索引
	InitPreOrbit(StartDOY,StartYear);
	if(orbitA[0].gpsTime.Week==orbitB[0].gpsTime.Week&&
		fabs(orbitA[0].gpsTime.SecOfWeek-orbitB[0].gpsTime.SecOfWeek)<1E-5)
	{
	
		KFState.Time.Week=orbitA[12].gpsTime.Week;
		KFState.Time.SecOfWeek=orbitA[12].gpsTime.SecOfWeek;
		xA_pre[0]=orbitA[12].x;
		xA_pre[1]=orbitA[12].y;
		xA_pre[2]=orbitA[12].z;
		xA_pre[3]=orbitA[12].dx;
		xA_pre[4]=orbitA[12].dy;
		xA_pre[5]=orbitA[12].dz;

		xB_pre[0]=orbitB[12].x;
		xB_pre[1]=orbitB[12].y;
		xB_pre[2]=orbitB[12].z;
		xB_pre[3]=orbitB[12].dx;
		xB_pre[4]=orbitB[12].dy;
		xB_pre[5]=orbitB[12].dz;
		ICRF_ITRF_GPST( MJD_J2000,  &KFState.Time, false, xA_pre_ECI, xA_pre);//转到惯性系
		ICRF_ITRF_GPST( MJD_J2000,  &KFState.Time, false, xB_pre_ECI, xB_pre);//转到惯性系
	}
	else
		return 0;
	do{
		double temptest=fmod(KFState.Time.SecOfWeek,86400.0);
		//if(fabs(KFState.Time.SecOfWeek-604770.0)<1E-5)
			//system("pause");
		//if( ReadEpochObs( Fobs_A, &ObsTypeList_A, &Epoch_A ) == 0||ReadEpochObs( Fobs_B, &ObsTypeList_B, &Epoch_B ) == 0 )//如果该天打开观测文件失败，则读取下一天数据
		if(fabs(fmod(KFState.Time.SecOfWeek,86400.0))<1E-5&&initFlag==1)//如果轨道位置比较到一天的最后一个历元，则读取下一个文件
		{
			
			//评价输出基线
			char baselineFile[60]="",baseVSFile[60]="",ASpeedFilePath[60]="",BSpeedFilePath[60]="";
			sprintf(baselineFile,"OutFile\\line.txt_%d_%d",(int)KFState.Step,KFState.Para[0].m_a);
			sprintf(baseVSFile,"OutFile\\speed.txt_%d_%d",(int)KFState.Step,KFState.Para[0].m_a);
			sprintf(ASpeedFilePath,"OutFile\\graceA_speed.txt_%d_%d",(int)KFState.Step,KFState.Para[0].m_a);
			sprintf(BSpeedFilePath,"OutFile\\graceB_speed.txt_%d_%d",(int)KFState.Step,KFState.Para[0].m_a);
			FILE * lineFile=fopen(baselineFile,"a+");
			FILE * speedFile=fopen(baseVSFile,"a+");
			FILE * ASpeedFile=fopen(ASpeedFilePath,"a+");
			FILE * BSpeedFile=fopen(BSpeedFilePath,"a+");
			double outDx,outDy,outDz,outDR;
			double outDDx,outDDy,outDDz,outDDR;
			for(i=0;i<numOfEpoch;i++){
				outDx=Lc_orbit[i].dx-dOrbitAB[i].dx;
				outDy=Lc_orbit[i].dy-dOrbitAB[i].dy;
				outDz=Lc_orbit[i].dz-dOrbitAB[i].dz; 
				outDR=sqrt(outDx*outDx+outDy*outDy+outDz*outDz);
				//输出相对位置误差文件
				fprintf( lineFile, " %10d %10d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
					Lc_orbit[i].week,(int)(Lc_orbit[i].sec+0.5), outDx,outDy,outDz,
					dOrbitAB[i].r,dOrbitAB[i].t,dOrbitAB[i].n ,outDR);
				outDDx=Lc_orbit[i].vx-dOrbitAB[i].ddx;
				outDDy=Lc_orbit[i].vy-dOrbitAB[i].ddy;
				outDDz=Lc_orbit[i].vz-dOrbitAB[i].ddz; 
				outDDR=sqrt(outDDx*outDDx+outDDy*outDDy+outDDz*outDDz);
				//输出相对速度误差文件
				fprintf(speedFile, "%10d%10d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
					Lc_orbit[i].week,(int)(Lc_orbit[i].sec+0.5), outDDx*1000,outDDy*1000,outDDz*1000,
					dOrbitAB[i].dr*1000,dOrbitAB[i].dt*1000,dOrbitAB[i].dn*1000 ,outDDR*1000);
			}
			fclose(lineFile);
			fclose(speedFile);
			//评价输出单星
			for(i=0;i<numOfEpoch;i++){
				//输出A星位置误差文件
				outDR=sqrt(dOrbitA[i].dx*dOrbitA[i].dx+dOrbitA[i].dy*dOrbitA[i].dy
					+dOrbitA[i].dz*dOrbitA[i].dz);
				fprintf( FRESIDUAL_A,"%10d%10d%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n",
					dOrbitA[i].week,(int)(dOrbitA[i].sec+0.5), dOrbitA[i].dx,dOrbitA[i].dy,dOrbitA[i].dz,
					dOrbitA[i].r,dOrbitA[i].t,dOrbitA[i].n,outDR);
				//输出A星速度误差文件
				outDDR=sqrt(dOrbitA[i].ddx*dOrbitA[i].ddx+dOrbitA[i].ddy*dOrbitA[i].ddy
					+dOrbitA[i].ddz*dOrbitA[i].ddz);
				fprintf( ASpeedFile,"%10d%10d%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n",
					dOrbitA[i].week,(int)(dOrbitA[i].sec+0.5), dOrbitA[i].ddx*1000,dOrbitA[i].ddy*1000,
					dOrbitA[i].ddz*1000,dOrbitA[i].dr*1000,dOrbitA[i].dt*1000,dOrbitA[i].dn*1000,outDDR*1E3);
				//输出B星位置误差文件
				outDR=sqrt(dOrbitB[i].dx*dOrbitB[i].dx+dOrbitB[i].dy*dOrbitB[i].dy
					+dOrbitB[i].dz*dOrbitB[i].dz);
				fprintf( FRESIDUAL_B,"%10d%10d%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n",
					dOrbitB[i].week,(int)(dOrbitB[i].sec+0.5),
					dOrbitB[i].dx,dOrbitB[i].dy,dOrbitB[i].dz,
					dOrbitB[i].r,dOrbitB[i].t,dOrbitB[i].n,outDR);
				//输出B星速度误差文件
				outDDR=sqrt(dOrbitB[i].ddx*dOrbitB[i].ddx+dOrbitB[i].ddy*dOrbitB[i].ddy
					+dOrbitB[i].ddz*dOrbitB[i].ddz);
				fprintf( BSpeedFile,"%10d%10d%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n",
					dOrbitB[i].week,(int)(dOrbitB[i].sec+0.5), dOrbitB[i].ddx*1E3,dOrbitB[i].ddy*1E3,dOrbitB[i].ddz*1E3,
					dOrbitB[i].dr*1E3,dOrbitB[i].dt*1E3,dOrbitB[i].dn*1E3,outDDR*1E3);
			}
			
			//Out_Residual(Lc_orbit,numOfEpoch,n_kbr);
			freeAll();
			/*
			fclose( Fobs_A);
			fclose( Fobs_B);
			fclose( FGPSEph );
			*/
			fclose( Fout_A );
			fclose( Fout_B );
			fclose( FRESIDUAL_A );
			fclose( FRESIDUAL_B );
			fclose( ASpeedFile);
			fclose( BSpeedFile);
			//评价结束
			time( &endt );    			
			printf("\nColapsed time: %lf\n", difftime( endt, start ) );
			//第二天
			StartDOY++;
			CheckDOY( StartYear, StartDOY );
			if( StartYear>=EndYear && StartDOY>EndDOY ) break;

			//重新初始化，前面已经free了
			//精密低轨卫星轨道初始化
			numOfEpoch=0;
			numOfA=0,numOfB=0;//当前A、B的GFZ真值索引
			InitPreOrbit(StartDOY,StartYear);
		
			initFile(Satname_A, Satname_B, EphName, StartDOY, StartYear, ObsFileName_A, ObsFileName_B,GPSEphFileName,
				ResFileName_A,ResFileName_B,RsdFileName_A,RsdFileName_B, KFState);

			printf( "Now processing %s......\n", ObsFileName_A );
			printf( "Now processing %s......\n", ObsFileName_B );
			/*
			Fout_A  = fopen( ResFileName_A, "wt" );
			Fout_B  = fopen( ResFileName_B, "wt" );

			if( (Fobs_A = fopen( ObsFileName_A, "rt" )) == NULL )
			{
				printf( "Cannot open GPSA Obs file. \n " );
				break;
			}
			if( (Fobs_B = fopen( ObsFileName_B, "rt" )) == NULL )
			{
				printf( "Cannot open GPSB Obs file. \n " );
				break;
			}
			if( (FGPSEph = fopen( GPSEphFileName, "rt" )) == NULL )
			{
				printf( "Cannot open %s GPS ephem file. \n", GPSEphFileName );
				break;
			}
			*/
			FRESIDUAL_A  = fopen( RsdFileName_A, "a+" );
			FRESIDUAL_B  = fopen( RsdFileName_B, "a+" );

			/* 观测文件读取 */
			/*
			if( ReadObsHead( Fobs_A, &ObsTypeList_A ) ==0 )
			{
				printf( "Some problems in Reading the header of ObsA file.\n" );
				break;
			}
			if( ReadObsHead( Fobs_B, &ObsTypeList_B ) ==0 )
			{
				printf( "Some problems in Reading the header of ObsB file.\n" );
				break;
			}
		
			ReadEpochObs( Fobs_A, &ObsTypeList_A, &Epoch_A );
			ReadEpochObs( Fobs_B, &ObsTypeList_B, &Epoch_B );
			*/
		}

		/*
		
		alignEpoch(Fobs_A, &ObsTypeList_A, &Epoch_A,Fobs_B, &ObsTypeList_B, &Epoch_B);//初始化成功后，对准后，可以直接用给后面测量更新，因为经过该步对准	
		
		/************************************************************************/
		//因为前面观测值已经对准，滤波也对准，故这里用A观测值的时间判断即可

		//注意KFState是在惯性系下的相心坐标，检验去掉转换直接质心结果是否相同

		//这里利用前一个历元assess的真值来作为初值
		initFlag=1;
		
		KFStateUpdate(xA_pre_ECI,xB_pre_ECI,&KFState);	
		for(int k=0;k<6;k++)
		{
			OiState[k] = KFState.StateA[k];
			OiState[k+54] = KFState.StateB[k];
		}
		ABTimeUpdate(OiState,&KFState);
		KFState.SatNumUsed[0]=0;
		KFState.SatNumUsed[1]=0;
		KFState.comSatNumUsed=0;
		assess(&numOfEpoch,&numOfA,&numOfB,&KFState,FinalOrb);
		for(int k=0;k<6;k++)
		{
			xA_pre_ECI[k] = KFState.StateA[k]; 
			xB_pre_ECI[k] = KFState.StateB[k]; 
		}
		/*
		//A星保持高精度，仅B星预报
		double delta[6];

		xA_pre[0] = orbitA[numOfA].x;
		xA_pre[1] = orbitA[numOfA].y;
		xA_pre[2] =	orbitA[numOfA].z;
		xA_pre[3] = orbitA[numOfA].dx;
		xA_pre[4] = orbitA[numOfA].dy;
		xA_pre[5] = orbitA[numOfA].dz;
		ICRF_ITRF_GPST( MJD_J2000,  &KFState.Time, false, xA_pre_ECI, xA_pre);//转到惯性系

		for(int k=0;k<6;k++)
		{
			delta[k] = KFState.StateB[k] - KFState.StateA[k];
			xB_pre_ECI[k] = xA_pre_ECI[k] + delta[k]; 
		}*/
		//这是利用JPL精密星历对滤波进行更新,去掉则是查看长期精度变化
	}while(1);
	//只能处理一天的，否则会因为内存问题错误
	freeKFState(KFState);
	time( &endt );    

	printf("\nColapsed time: %12.3f\n", difftime( endt, start ) );
	//free(PRE_EPH);
	printf( "Enter a key to exit...\n" );
//	getchar();
	//system("pause");
	return 0;

}
void KFStateUpdate(double * xA_pre,double * xB_pre,EKFSTATE *  KFState)
{

	PhaseCentToMassCent( false, &KFState->CentBias[0], xA_pre);    /* 转换为相位中心 */
	PhaseCentToMassCent( false, &KFState->CentBias[3], xB_pre);    /* 转换为相位中心 */
	for(int i=0;i<6;i++)
	{
		KFState->StateA[i] = xA_pre[i];
		KFState->StateB[i] = xB_pre[i];
	}
	getRelPos(KFState->StateA,KFState->StateB,KFState->StateRel);
}