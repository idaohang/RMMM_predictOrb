#include "qualityCheck.h"

const double chisqr[100]={      /* chi-sqr(n) (alpha=0.001) */
	10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
	31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
	46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
	61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
	74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
	88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
	101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
	113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
	126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
	138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};
/*绝对定轨粗差剔除*/
double LS_PC(const int satNum, int * PRN,int * banPRN,const double * satPos,double * ApriPos,const double * obs
	,double * AprioriClk, PPRESULT * Result)
{
	int i, j, k;                     /*  Iterator为单点定位迭代次数  */
	int SatNumUsed,Iterator;      /*  单点定位使用的卫星数 */
	double ATA[16], InvATA[16];
	double ATB[4];
	double dPos[4];               /* 定位结果  */
	double BPos[MAXCHANNUM];                      /* 伪距或多普勒的观测值-计算值 */
	double MeasA[MAXCHANNUM*4], MeasAT[MAXCHANNUM*4];          /* 观测矩阵 [MAXCHANNUM][5] */
	double Weight[MAXCHANNUM];                               /* 权矩阵, 单位阵只取对角线元素 */
	double Qvv[MAXCHANNUM];                               /* 观测值改正数的斜因数阵对角线元素 */
	double Range;                                 /*  接收机与导航卫星之间的距离  */
	double Coverage,SigmaPos=0.0;
	double Residual[MAXCHANNUM] = {0.0};  /* 定位残差  */
	int valid =0;
	Iterator = 0;
	do{
		SatNumUsed = 0;

		for( i=0; i<MAXCHANNUM; i++ )
		{
			valid=1;
			for (j=0;j<12;j++)
			{
				if(PRN[i]==999|PRN[i]==banPRN[j]){
					valid=0;
					break;
				}
			}
			if (valid==1)
			{
				Range = 0.0;
				for( k=0; k<3; k++ )
				{
					Range = Range + pow( satPos[i*3+k]-ApriPos[k], 2.0 );
				}

				Range = sqrt( Range );
				MeasA[SatNumUsed*4+0] = 1.0*( ApriPos[0] - satPos[i*3+0] )/Range;
				MeasA[SatNumUsed*4+1] = 1.0*( ApriPos[1] - satPos[i*3+1] )/Range;
				MeasA[SatNumUsed*4+2] = 1.0*( ApriPos[2] - satPos[i*3+2] )/Range;
				MeasA[SatNumUsed*4+3] = 1.0;             /* 接收机GPS系统钟差系数 */
				Weight[SatNumUsed]    = 1.0;  // * pow( sin(SatMidInfo[i].Elevation), 2.0 );

				BPos[SatNumUsed] =obs[i]-Range-*AprioriClk;
				SatNumUsed++;
			}


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

		if( SatNumUsed>=4 )  /* 肯定大于 */
		{
			MatrixTranspose( SatNumUsed, 4, MeasA, MeasAT );
			MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 4, MeasAT, Weight, MeasA, ATA );
			MatrixInv( 4, ATA, InvATA );

			MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 1, MeasAT, Weight, BPos, ATB );
			MatrixMultiply( 4, 4, 4, 1, InvATA, ATB, dPos );

			*AprioriClk = *AprioriClk + dPos[3];

		}
		for( k=0; k<3; k++ )
		{
			ApriPos[k] = ApriPos[k] + dPos[k];
		}
		Iterator++;
		Coverage = VectDot( 3, 3, dPos, dPos );


	}while ( (Coverage>1E-6) && (Iterator < 10) );
	//循环结果计算高度角

	/* 单独定位结果整理输出  */
	/* 计算观测值残差 */
	MatrixMultiply( SatNumUsed, 4, 4, 1, MeasA, dPos, Residual );//
	for( i=0; i<SatNumUsed; i++ )
	{
		Result->Residual[i] = Residual[i] - BPos[i];
	}

	/* 计算观测值改正数的协因数阵 */

	ComputeQvv( SatNumUsed, 4, MeasA, InvATA, MeasAT, Weight, Qvv );

	Result->PDOP = sqrt( InvATA[0] + InvATA[5] + InvATA[10] );
	Result->GDOP = sqrt( InvATA[0] + InvATA[5] + InvATA[10] + InvATA[15] );
	Result->HDOP = sqrt( InvATA[0] + InvATA[5] );
	Result->VDOP = sqrt( InvATA[10] );
	Result->TDOP = sqrt( InvATA[15] );
	Result->Iterator = Iterator;
	Result->Coverage = Coverage;
	if(Iterator<10)
	{
		for( i=0; i<SatNumUsed; i++ )
		{
			SigmaPos = SigmaPos + Weight[i]*pow( Result->Residual[i], 2.0 );
		}
		SigmaPos = sqrt( SigmaPos / (SatNumUsed-4) );
	}
	else 
		SigmaPos = 999.0;
	return SigmaPos;
}

double LSRel_PC(const int satNum, int * PRN,int * banPRN,const double * satPos,const double * StateInECEFA,const double * obs
	,double * AprioriPosRel,double * AprioriClk, PPRESULT * Result)
{
	int i, j, k;                     /*  Iterator为单点定位迭代次数  */
	int SatNumUsed,Iterator;      /*  单点定位使用的卫星数 */
	double ATA[16], InvATA[16];
	double ATB[4];
	double dPos[4];               /* 定位结果  */
	double BPos[MAXCHANNUM];                      /* 伪距或多普勒的观测值-计算值 */
	double MeasA[MAXCHANNUM*4], MeasAT[MAXCHANNUM*4];          /* 观测矩阵 [MAXCHANNUM][5] */
	double Weight[MAXCHANNUM];                               /* 权矩阵, 单位阵只取对角线元素 */
	double Qvv[MAXCHANNUM];                               /* 观测值改正数的斜因数阵对角线元素 */
	double RangeA,RangeB;                                 /*  接收机与导航卫星之间的距离  */
	double Coverage,SigmaPos=0.0;
	double Residual[MAXCHANNUM] = {0.0};  /* 定位残差  */
	int valid =0;
	Iterator = 0;
	do{
		SatNumUsed = 0;

		for( i=0; i<MAXCHANNUM; i++ )
		{
			valid=1;
			for (j=0;j<12;j++)
			{
				if(PRN[i]==999|PRN[i]==banPRN[j]){
					valid=0;
					break;
				}
			}
			if (valid==1)
			{
				RangeA = 0.0;
				RangeB = 0.0;
				for( k=0; k<3; k++ )
				{
					RangeA = RangeA + pow( satPos[i*3+k]-StateInECEFA[k], 2.0 );
					RangeB = RangeB + pow( satPos[i*3+k]-StateInECEFA[k]-AprioriPosRel[k], 2.0 );
				}

				RangeA = sqrt( RangeA );
				RangeB = sqrt( RangeB );
				MeasA[SatNumUsed*4+0] = 1.0*( AprioriPosRel[0] + StateInECEFA[0] - satPos[i*3+0] )/RangeB;
				MeasA[SatNumUsed*4+1] = 1.0*( AprioriPosRel[1] + StateInECEFA[1] - satPos[i*3+1] )/RangeB;
				MeasA[SatNumUsed*4+2] = 1.0*( AprioriPosRel[2] + StateInECEFA[2] - satPos[i*3+2] )/RangeB;
				MeasA[SatNumUsed*4+3] = 1.0;             /* 接收机GPS系统钟差系数 */
				Weight[SatNumUsed]    = 1.0;  // * pow( sin(SatMidInfo[i].Elevation), 2.0 );

				BPos[SatNumUsed] =obs[i]-(RangeB-RangeA)-*AprioriClk;
				SatNumUsed++;
			}


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

		if( SatNumUsed>=4 )  /* 肯定大于 */
		{
			MatrixTranspose( SatNumUsed, 4, MeasA, MeasAT );
			MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 4, MeasAT, Weight, MeasA, ATA );
			MatrixInv( 4, ATA, InvATA );

			MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 1, MeasAT, Weight, BPos, ATB );
			MatrixMultiply( 4, 4, 4, 1, InvATA, ATB, dPos );

			*AprioriClk = *AprioriClk + dPos[3];

		}
		for( k=0; k<3; k++ )
		{
			AprioriPosRel[k] = AprioriPosRel[k] + dPos[k];
		}
		Iterator++;
		Coverage = VectDot( 3, 3, dPos, dPos );


	}while ( (Coverage>1E-6) && (Iterator < 10) );
	//循环结果计算高度角

	/* 单独定位结果整理输出  */
	/* 计算观测值残差 */
	MatrixMultiply( SatNumUsed, 4, 4, 1, MeasA, dPos, Residual );//
	for( i=0; i<SatNumUsed; i++ )
	{
		Result->Residual[i] = Residual[i] - BPos[i];
	}

	/* 计算观测值改正数的协因数阵 */

	ComputeQvv( SatNumUsed, 4, MeasA, InvATA, MeasAT, Weight, Qvv );

	Result->PDOP = sqrt( InvATA[0] + InvATA[5] + InvATA[10] );
	Result->GDOP = sqrt( InvATA[0] + InvATA[5] + InvATA[10] + InvATA[15] );
	Result->HDOP = sqrt( InvATA[0] + InvATA[5] );
	Result->VDOP = sqrt( InvATA[10] );
	Result->TDOP = sqrt( InvATA[15] );
	Result->Iterator = Iterator;
	Result->Coverage = Coverage;
	if(Iterator<10)
	{
		for( i=0; i<SatNumUsed; i++ )
		{
			SigmaPos = SigmaPos + Weight[i]*pow( Result->Residual[i], 2.0 );
		}
		SigmaPos = sqrt( SigmaPos / (SatNumUsed-4) );
	}
	else 
		SigmaPos = 999.0;
	return SigmaPos;
}

double LEGE(GPSTIME Time, const int satNum, int * PRN,int * banPRN,const double * satPos,const double * StateInECEFA,const double * obs
	,double * AprioriPosRel,double * AprioriClk, PPRESULT * Result)
{
	int i,j;
	int SatNumUsed=satNum;
	int bestBan=999;
	double bestSigma0,bestSigma1,tempSigma=999.0;
	bestSigma0=LSRel_PC(SatNumUsed,PRN,banPRN,satPos,StateInECEFA,obs,AprioriPosRel,AprioriClk, Result);
	for (i=0;i<SatNumUsed-5;i++)
	{
		bestSigma1=999.0;
		for (j=0;j<SatNumUsed;j++)
		{
			banPRN[i]=PRN[j];
			tempSigma=LSRel_PC(SatNumUsed-i-1,PRN,banPRN,satPos,StateInECEFA,obs,AprioriPosRel,AprioriClk, Result);
			banPRN[i]=999;
			if (tempSigma<bestSigma1)
			{
				bestSigma1=tempSigma;
				bestBan=PRN[j];
			}
		}
		//如果又剔除一颗卫星后改善，则剔除最好的
		if (bestSigma0>bestSigma1)
		{
			banPRN[i]=bestBan;
		}
		//如果剔除后反而变差，则停止剔除。
		else
		{
			banPRN[i]=999;
			break;
		}
		//如果剔除后，改善到一定程度，则也停止剔除
		if (bestSigma0>3*bestSigma1)
		{
			FILE * lege=fopen("OutFile\\lege.txt","a+");
			fprintf(lege,"%14.3f,%14.3f\n",Time.SecOfWeek,bestSigma1);
			fclose(lege);
			SatNumUsed=SatNumUsed-i-1;//剔除后的卫星数
			Result->SigmaPos=LSRel_PC(SatNumUsed,PRN,banPRN,satPos,StateInECEFA,obs,AprioriPosRel,AprioriClk, Result);
			Result->SatNumUsed = SatNumUsed;
			break;
		}
		//没有达到理想程度则继续剔除。
		else{
			bestSigma0=bestSigma1;
		}
	}
	//剔除到只剩下5颗卫星时，则按原本所有卫星计算
	if(i==SatNumUsed-5)
	{
		for (j=0;j<MAXCHANNUM;j++)
		{
			banPRN[j]=999;
		}
		Result->SigmaPos=LSRel_PC(SatNumUsed,PRN,banPRN,satPos,StateInECEFA,obs,AprioriPosRel,AprioriClk, Result);
		Result->SatNumUsed = SatNumUsed;
	}
	return Result->SigmaPos;
}
/****************************************************************************
Raim_fde1

目的:
参数: GPSTIME Time
参数: const int satNum
参数: int * PRN
参数: int * banPRN
参数: const double * satPos
参数: const double * StateInECEFA
参数: const double * obs
参数: double * AprioriPosRel
参数: double * AprioriClk
参数: PPRESULT * Result
****************************************************************************/
double Raim_fde1(int * satNum, int * PRN,int * banPRN,const double * satPos,const double * StateInECEFA,const double * obs
	,double * AprioriPosRel,double * AprioriClk, PPRESULT * Result)
{
	int i,j;
	int SatNumUsed=*satNum;
	int bestBan=999;
	double bestSigma1,tempSigma=999.0;
	if(SatNumUsed>=6)
	{
		for (j=0;j<SatNumUsed;j++)
		{
			banPRN[0]=PRN[j];
			tempSigma=LSRel_PC(SatNumUsed-1,PRN,banPRN,satPos,StateInECEFA,obs,AprioriPosRel,AprioriClk, Result);
			banPRN[0]=999;
			if (j==0)
			{
				bestSigma1=tempSigma;
				continue;
			}
			if (tempSigma<bestSigma1)
			{
				bestSigma1=tempSigma;
				bestBan=PRN[j];
			}
		}
		banPRN[0]=bestBan;
		SatNumUsed=SatNumUsed-1;//剔除后的卫星数
		Result->SigmaPos=LSRel_PC(SatNumUsed,PRN,banPRN,satPos,StateInECEFA,obs,AprioriPosRel,AprioriClk, Result);
		*satNum=SatNumUsed;
	}
	else{
		Result->SigmaPos=999.0;
		Result->SatNumUsed=0;
	}
	return Result->SigmaPos;
}
//Raim_spp
void Raim_spp(EpochObsData * Epoch)
{
	int i,j;
	int SatNumUsed=0;
	int bestBan=999;
	double bestSigma1,tempSigma=999.0;
}
/****************************************************************************
valResult

目的:
参数: sigmaPos			验后中误差
参数: simga0			验前中误差  //与观测噪声有关，因为是前面设置的是单位权
参数: int nv			观测值数
参数: int nx			估计参数个数
****************************************************************************/
int valResult(const double sigmaPos,const double sigma0,int nv,int nx)
{
	double VTPV=0.0,GDop=0.0;
	if (nv>nx)
	{
		VTPV=sigmaPos * sigmaPos*(nv-nx)/(sigma0*sigma0);
		if(VTPV>chisqr[nv-nx-1]){
			printf("chi-square error nv=%d vv=%.1f cs=%.1f",nv,VTPV,chisqr[nv-nx-1]);
			return 0;
		}
	}
	return 1;
}