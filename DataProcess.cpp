#include "DataProcess.h"
int GenCommomObs(EpochObsData * Epoch_A,EpochObsData * Epoch_B,Common11Obs * tempComObs)
{
	int i,j,k;
	int commonObsNum=0;
	tempComObs->Time.Week		=Epoch_A->Time.Week;
	tempComObs->Time.SecOfWeek	=Epoch_A->Time.SecOfWeek;
	for(i=0;i<Epoch_A->SatNum;i++)
		for(j=0;j<Epoch_B->SatNum;j++)
		{
			if(Epoch_A->SatObs[i].System!=GPS||Epoch_A->SatObs[i].Used==0)//used==0表示高度角小于10度，used==-1表示伪距粗差检定
				break;
			if(Epoch_B->SatObs[j].System!=GPS||Epoch_B->SatObs[j].Used==0)//
				continue;
			if(Epoch_A->SatObs[i].Prn==Epoch_B->SatObs[j].Prn)
			{
				tempComObs->comobs[commonObsNum].PRN=Epoch_A->SatObs[i].Prn;
				for(k=0;k<3;k++)
					tempComObs->comobs[commonObsNum].satPos[k]=Epoch_A->SatObs[i].satPos[k];
				/*由单点定位计算的高度角传过来的，精度较低*/
				tempComObs->comobs[commonObsNum].elevation[0]=Epoch_A->SatObs[i].elevation;
				tempComObs->comobs[commonObsNum].elevation[1]=Epoch_B->SatObs[j].elevation;
				for(k=0;k<MAXOBSTYPENUM;k++){
					//A星
					if(Epoch_A->SatObs[i].Data[k].Availability==1&&Epoch_A->SatObs[i].Data[k].Type==C1
						&&Epoch_A->SatObs[i].Data[k].LLI!=9)
						tempComObs->comobs[commonObsNum].CA1=Epoch_A->SatObs[i].Data[k].Obs;
					if(Epoch_A->SatObs[i].Data[k].Availability==1&&Epoch_A->SatObs[i].Data[k].Type==P1
						&&Epoch_A->SatObs[i].Data[k].LLI!=9)
						tempComObs->comobs[commonObsNum].P11=Epoch_A->SatObs[i].Data[k].Obs;
					if(Epoch_A->SatObs[i].Data[k].Availability==1&&Epoch_A->SatObs[i].Data[k].Type==P2
						&&Epoch_A->SatObs[i].Data[k].LLI!=9)
						tempComObs->comobs[commonObsNum].P12=Epoch_A->SatObs[i].Data[k].Obs;
					if(Epoch_A->SatObs[i].Data[k].Availability==1&&Epoch_A->SatObs[i].Data[k].Type==L1){
						tempComObs->comobs[commonObsNum].L11=Epoch_A->SatObs[i].Data[k].Obs;
						tempComObs->comobs[commonObsNum].LLI_L1A=Epoch_A->SatObs[i].Data[k].LLI;
					}
					if(Epoch_A->SatObs[i].Data[k].Availability==1&&Epoch_A->SatObs[i].Data[k].Type==L2){
						tempComObs->comobs[commonObsNum].L12=Epoch_A->SatObs[i].Data[k].Obs;
						tempComObs->comobs[commonObsNum].LLI_L2A=Epoch_A->SatObs[i].Data[k].LLI;
					}
					if(Epoch_A->SatObs[i].Data[k].Availability==1&&Epoch_A->SatObs[i].Data[k].Type==S1)
						tempComObs->comobs[commonObsNum].S11=Epoch_A->SatObs[i].Data[k].Obs;
					if(Epoch_A->SatObs[i].Data[k].Availability==1&&Epoch_A->SatObs[i].Data[k].Type==S2)
						tempComObs->comobs[commonObsNum].S12=Epoch_A->SatObs[i].Data[k].Obs;
					//B星
					if(Epoch_B->SatObs[j].Data[k].Availability==1&&Epoch_B->SatObs[j].Data[k].Type==C1
						&&Epoch_B->SatObs[j].Data[k].LLI!=9)
						tempComObs->comobs[commonObsNum].CA2=Epoch_B->SatObs[j].Data[k].Obs;
					if(Epoch_B->SatObs[j].Data[k].Availability==1&&Epoch_B->SatObs[j].Data[k].Type==P1
						&&Epoch_B->SatObs[j].Data[k].LLI!=9)
						tempComObs->comobs[commonObsNum].P21=Epoch_B->SatObs[j].Data[k].Obs;
					if(Epoch_B->SatObs[j].Data[k].Availability==1&&Epoch_B->SatObs[j].Data[k].Type==P2
						&&Epoch_B->SatObs[j].Data[k].LLI!=9)
						tempComObs->comobs[commonObsNum].P22=Epoch_B->SatObs[j].Data[k].Obs;
					if(Epoch_B->SatObs[j].Data[k].Availability==1&&Epoch_B->SatObs[j].Data[k].Type==L1){
						tempComObs->comobs[commonObsNum].L21=Epoch_B->SatObs[j].Data[k].Obs;
						tempComObs->comobs[commonObsNum].LLI_L1B=Epoch_B->SatObs[j].Data[k].LLI;
					}
					if(Epoch_B->SatObs[j].Data[k].Availability==1&&Epoch_B->SatObs[j].Data[k].Type==L2){
						tempComObs->comobs[commonObsNum].L22=Epoch_B->SatObs[j].Data[k].Obs;
						tempComObs->comobs[commonObsNum].LLI_L2B=Epoch_B->SatObs[j].Data[k].LLI;
					}
					if(Epoch_B->SatObs[j].Data[k].Availability==1&&Epoch_B->SatObs[j].Data[k].Type==S1)
						tempComObs->comobs[commonObsNum].S21=Epoch_B->SatObs[j].Data[k].Obs;
					if(Epoch_B->SatObs[j].Data[k].Availability==1&&Epoch_B->SatObs[j].Data[k].Type==S2)
						tempComObs->comobs[commonObsNum].S22=Epoch_B->SatObs[j].Data[k].Obs;
				}
				if(fabs(tempComObs->comobs[commonObsNum].CA1)>0.1&&fabs(tempComObs->comobs[commonObsNum].CA2)>0.1){
						tempComObs->comobs[commonObsNum].CAApplyEpoch=1;
						tempComObs->comobs[commonObsNum].dCA=tempComObs->comobs[commonObsNum].CA2-tempComObs->comobs[commonObsNum].CA1;
				}
				if(fabs(tempComObs->comobs[commonObsNum].L11)>0.1&&fabs(tempComObs->comobs[commonObsNum].L21)>0.1){
						tempComObs->comobs[commonObsNum].L1ApplyEpoch=1;
						tempComObs->comobs[commonObsNum].dL1=tempComObs->comobs[commonObsNum].L21-tempComObs->comobs[commonObsNum].L11;
				}
				if(fabs(tempComObs->comobs[commonObsNum].L12)>0.1&&fabs(tempComObs->comobs[commonObsNum].L22)>0.1){
						tempComObs->comobs[commonObsNum].L2ApplyEpoch=1;
						tempComObs->comobs[commonObsNum].dL2=tempComObs->comobs[commonObsNum].L22-tempComObs->comobs[commonObsNum].L12;
				}
				if(fabs(tempComObs->comobs[commonObsNum].P11)>0.1&&fabs(tempComObs->comobs[commonObsNum].P21)>0.1){
						tempComObs->comobs[commonObsNum].P1ApplyEpoch=1;
						tempComObs->comobs[commonObsNum].dP1=tempComObs->comobs[commonObsNum].P21-tempComObs->comobs[commonObsNum].P11;
				}
				if(fabs(tempComObs->comobs[commonObsNum].P12)>0.1&&fabs(tempComObs->comobs[commonObsNum].P22)>0.1){
						tempComObs->comobs[commonObsNum].P2ApplyEpoch=1;
						tempComObs->comobs[commonObsNum].dP2=tempComObs->comobs[commonObsNum].P22-tempComObs->comobs[commonObsNum].P12;
				}
				tempComObs->comobs[commonObsNum].flag=0;//初始则认为都不是周跳
				tempComObs->comobs[commonObsNum].used=1;//初始认为都不是粗差
				//码观测值差分检验粗差
				if ( fabs(tempComObs->comobs[commonObsNum].P11-tempComObs->comobs[commonObsNum].P12) > 30.0||
					fabs(tempComObs->comobs[commonObsNum].P21-tempComObs->comobs[commonObsNum].P22) > 30.0
					)
				{
					tempComObs->comobs[commonObsNum].P2ApplyEpoch= 0;							//P2有粗差
				}
				if ( fabs(tempComObs->comobs[commonObsNum].P11-tempComObs->comobs[commonObsNum].CA1)>10.0||
					fabs(tempComObs->comobs[commonObsNum].P21-tempComObs->comobs[commonObsNum].CA2)>10.0)
				{
					tempComObs->comobs[commonObsNum].P1ApplyEpoch= 0;							//P1有粗差
				}
				if (Epoch_A->SatObs[i].Used==-1||Epoch_B->SatObs[j].Used==-1)			//如果前面伪距探测为粗差
				{
					tempComObs->comobs[commonObsNum].used=-1;
				}
				
				if(tempComObs->comobs[commonObsNum].L1ApplyEpoch==1&&tempComObs->comobs[commonObsNum].L2ApplyEpoch==1
					&&tempComObs->comobs[commonObsNum].P1ApplyEpoch==1&&tempComObs->comobs[commonObsNum].P2ApplyEpoch==1)
					commonObsNum++;//目前只考虑P1、P2、L1、L2观测值均存在其一情况下作为共同观测卫星
				else
					EmptyOneSat11Obs(&tempComObs->comobs[commonObsNum]);
				break;
			}
		}
		tempComObs->ComSatNum=commonObsNum;
		return 1;
}
void copyOneSat11Obs(OneSat11Obs * target,OneSat11Obs * src)
{
	int i=0;
	target->used		=	src->used;
	for(i=0;i<2;i++)
		target->elevation[i]=src->elevation[i];
	for(i=0;i<3;i++)
		target->satPos[i]=	src->satPos[i];
	target->map[0]		=	src->map[0];
	target->map[1]		=	src->map[1];
	target->CA1			=	src->CA1;
	target->CA2			=	src->CA2;
	target->CAApplyEpoch=	src->CAApplyEpoch;
	target->dCA			=	src->dCA;
	target->dL1			=	src->dL1;
	target->dL2			=	src->dL2;
	target->dP1			=	src->dP1;
	target->dP2			=	src->dP2;
	target->elevation[0]=	src->elevation[0];
	target->elevation[1]=	src->elevation[1];
	target->flag		=	src->flag;
	target->L11			=	src->L11;
	target->L12			=	src->L12;
	target->L1ApplyEpoch=	src->L1ApplyEpoch;
	target->L21			=	src->L21;
	target->L22			=	src->L22;
	target->L2ApplyEpoch=	src->L2ApplyEpoch;
	target->P11			=	src->P11;
	target->P12			=	src->P12;
	target->P1ApplyEpoch=	src->P1ApplyEpoch;
	target->P21			=	src->P21;
	target->P22			=	src->P22;
	target->P2ApplyEpoch=	src->P2ApplyEpoch;
	target->PRN			=	src->PRN;
	target->S11			=	src->S11;
	target->S12			=	src->S12;
	target->S21			=	src->S21;
	target->S22			=	src->S22;
	//MW周跳探测中间变量
	target->Nw_SD		=	src->Nw_SD;
	target->sigNw_SD	=	src->sigNw_SD;
	target->Nw_num_SD	=	src->Nw_num_SD;
	//MW单差探测中间变量
	target->Nw_A		=	src->Nw_A;
	target->sigNw_A		=	src->sigNw_A;
	target->Nw_num_A	=	src->Nw_num_A;
	//
	target->Nw_B		=	src->Nw_B;
	target->sigNw_B		=	src->sigNw_B;
	target->Nw_num_B	=	src->Nw_num_B;
}
void copyCommon11Obs(Common11Obs * targetComObs,Common11Obs * srcComObs)
{
	int i;
	targetComObs->ComSatNum		=	srcComObs->ComSatNum;
	targetComObs->Time.SecOfWeek=	srcComObs->Time.SecOfWeek;
	targetComObs->Time.Week		=	srcComObs->Time.Week;
	for(i=0;i<srcComObs->ComSatNum;i++){
		copyOneSat11Obs(&targetComObs->comobs[i],&srcComObs->comobs[i]);
	}
	for (i=srcComObs->ComSatNum;i<MAXCHANNUM;i++)
	{
		EmptyOneSat11Obs(&targetComObs->comobs[i]);
	}
}

//当失锁后所有卫星重新观测，初始化模糊度与方差模块
//这里初始化模糊度利用滤波相心位置，能加快收敛，但是也容易受到影响，如果位置精确，还可以加上相位中心偏差
void initialAllAmb(double * Mat, EKFSTATE * EkfState)
{
	int i,j;
	double rhoA,rhoB,dx,dy,dz,dtr;
	double B[48]={0.0},Dl[12]={0.0};
	double stateInECEFRel[6]={0.0};
	/************************************************************************/
	/* 初始化所需初值                                                       */
	/************************************************************************/
	dtr=EkfState->StateRel[6];
	MatrixMultiply( 3, 3, 3, 1, Mat, EkfState->StateA, EkfState->StateInECEFA );//每个观测值更新都要转到地固系中
	MatrixMultiply( 3, 3, 3, 1, Mat, EkfState->StateRel, stateInECEFRel );
	getBpos(EkfState->StateInECEFA,stateInECEFRel,EkfState->StateInECEFB);
	/************************************************************************/
	/* 初始化归零                                                           */
	/************************************************************************/
	for(i=0;i<12;i++){
		EkfState->StateRel[12+i]=0.0;
		EkfState->StateRel[24+i]=0.0;
		for(j=0;j<RELDIMENSION;j++)
		{
			EkfState->CovRel[(12+i)+RELDIMENSION*j]=0.0;		//第i个L1模糊度列全置为0.0;
			EkfState->CovRel[(12+i)*RELDIMENSION+j]=0.0;		//第i个L1模糊度行全置为0.0;
			EkfState->CovRel[(24+i)+RELDIMENSION*j]=0.0;		//第i个宽巷模糊度列全置为0.0;
			EkfState->CovRel[(24+i)*RELDIMENSION+j]=0.0;		//第15+i行全置为0.0;
		}
	}
	for(i=0;i<EkfState->CurComObs.ComSatNum;i++)
	{
		dx	=	EkfState->StateInECEFA[0]-EkfState->CurComObs.comobs[i].satPos[0];
		dy	=	EkfState->StateInECEFA[1]-EkfState->CurComObs.comobs[i].satPos[1];
		dz	=	EkfState->StateInECEFA[2]-EkfState->CurComObs.comobs[i].satPos[2];
		rhoA=	sqrt(dx*dx+dy*dy+dz*dz);
		dx=		EkfState->StateInECEFB[0]-EkfState->CurComObs.comobs[i].satPos[0];
		dy=		EkfState->StateInECEFB[1]-EkfState->CurComObs.comobs[i].satPos[1];
		dz=		EkfState->StateInECEFB[2]-EkfState->CurComObs.comobs[i].satPos[2];
		rhoB=	sqrt(dx*dx+dy*dy+dz*dz);
		EkfState->PRN[i]		=	EkfState->CurComObs.comobs[i].PRN;
		//模糊度初始化
		EkfState->StateRel[24+i] = ( EkfState->CurComObs.comobs[i].dL1 - EkfState->CurComObs.comobs[i].dL2 )
										   - LEMTA_NL/LEMTA_WL * ( EkfState->CurComObs.comobs[i].dP1/LEMTA_1 + EkfState->CurComObs.comobs[i].dP2/LEMTA_2 );
		EkfState->StateRel[12+i]				=	1.0/(1.0-RATIO2)*EkfState->CurComObs.comobs[i].dL1*LEMTA_1-RATIO2/(1.0-RATIO2)*EkfState->CurComObs.comobs[i].dL2*LEMTA_2
											   -(rhoB-rhoA)-dtr;//这里N=(l-Bx)/lemta，因此方差初始化通过误差传播
		//利用位置初始化LC									   
		/*EkfState->State[15+i]				=	1.0/(1.0-RATIO)*EkfState->CurComObs.comobs[i].dL1-RATIO/(1.0-RATIO)*EkfState->CurComObs.comobs[i].dL2
											   -(rhoB-rhoA)/LEMTA_NL-dtr/LEMTA_NL-EkfState->State[EKFDIMENSION+15+i]*RATIO/(1.0-RATIO);
		利用位置初始化MW与L1*/
		//模糊度方差初始化
		EkfState->CovRel[(12+i)*RELDIMENSION+(12+i)]=SIGL1;//模糊度初始化为3周
		EkfState->CovRel[(24+i)*RELDIMENSION+(24+i)]=SIGMW;//模糊度初始化为3周
		//方差初始化所需观测矩阵
		B[i*4+0] = dx/rhoB;
		B[i*4+1] = dy/rhoB;
		B[i*4+2] = dz/rhoB;
		B[i*4+3] = 1.0;
		Dl[i] = Me_dLC;
		//MW探测周跳中间变量初始化
		EkfState->CurComObs.comobs[i].Nw_SD = ( EkfState->CurComObs.comobs[i].dL1 - EkfState->CurComObs.comobs[i].dL2 )
									     - LEMTA_NL/LEMTA_WL*( EkfState->CurComObs.comobs[i].dP1/LEMTA_1 + EkfState->CurComObs.comobs[i].dP2/LEMTA_2 );
		EkfState->CurComObs.comobs[i].Nw_num_SD = 1;
		EkfState->CurComObs.comobs[i].sigNw_SD	 = MW_SIG;
		//MW非差模糊度初始化
		EkfState->CurComObs.comobs[i].Nw_A = ( EkfState->CurComObs.comobs[i].L11 - EkfState->CurComObs.comobs[i].L12 )
			- LEMTA_NL/LEMTA_WL*( EkfState->CurComObs.comobs[i].P11/LEMTA_1 + EkfState->CurComObs.comobs[i].P12/LEMTA_2 );
		EkfState->CurComObs.comobs[i].Nw_num_A = 1;
		EkfState->CurComObs.comobs[i].sigNw_A	 = MW_SIG;
		EkfState->CurComObs.comobs[i].Nw_B = ( EkfState->CurComObs.comobs[i].L21 - EkfState->CurComObs.comobs[i].L22 )
			- LEMTA_NL/LEMTA_WL*( EkfState->CurComObs.comobs[i].P21/LEMTA_1 + EkfState->CurComObs.comobs[i].P22/LEMTA_2 );
		EkfState->CurComObs.comobs[i].Nw_num_B = 1;
		EkfState->CurComObs.comobs[i].sigNw_B	 = MW_SIG;
	}
	//initialAmbsCov(B,Dl,EkfState,EkfState->Sigma[1]*EkfState->Sigma[1]);
	
}
void initialAmbsCov(const double B[48],const double Dl[12],EKFSTATE * EkfState,const double sigma)
{
	int i,j;
	double Qxx[16]={0.0},BT[48],BQ[48]={0.0},TempQ[144]={0.0};
	/************************************************************************/
	/* 方差初始化通过误差传播定律QNN=Dl+B*QXX*BT                            */
	/************************************************************************/
	//获得Qxx
	for(i=0;i<4;i++){
		for(j=0;j<4;j++)
		{
			Qxx[i*4+j]=EkfState->CovRel[i*RELDIMENSION+j];
		}
	}
	//获得TempQ=B*QXX*BT
	MatrixMultiply(EkfState->CurComObs.ComSatNum,4,4,4,B,Qxx,BQ);
	MatrixTranspose(EkfState->CurComObs.ComSatNum,4,B,BT);
	MatrixMultiply(EkfState->CurComObs.ComSatNum,4,4,EkfState->CurComObs.ComSatNum,BQ,BT,TempQ);
	MatrixMultiplyk(EkfState->CurComObs.ComSatNum,EkfState->CurComObs.ComSatNum,TempQ,sigma);
	MatrixAddition3(EkfState->CurComObs.ComSatNum,Dl,TempQ);
	//返回QNN
	for (i=0;i<EkfState->CurComObs.ComSatNum;i++)
	{
		for (j=0;j<EkfState->CurComObs.ComSatNum;j++)
		{
			EkfState->CovRel[(i+12)*RELDIMENSION+j+12]=TempQ[i*EkfState->CurComObs.ComSatNum+j];
		}
	}
}
/****************************************************************************
函数名:deletePRNCov
目的:删除索引为index的卫星模糊度相关协方差
思路;通过从index开始，先将，最终将最后一个清0
输入:
		index		需要删除的卫星在卫星序列中的索引
		ComSatNum	卫星序列中的总卫星数
		Cov			滤波的协方差
输出:
		Cov			删除目标卫星信息后的协方差
****************************************************************************/
void deletePRNCov(int index,int ComSatNum,double * Cov)
{
	int i,j;
	//先行覆盖
	for(i=index;i<ComSatNum-1;i++)
	{
		for(j=0;j<RELDIMENSION;j++)
		{
			Cov[(i+12)*RELDIMENSION+j]	=	Cov[(i+13)*RELDIMENSION+j];//L1模糊度协方差行覆盖
			Cov[(i+24)*RELDIMENSION+j]	=	Cov[(i+25)*RELDIMENSION+j];//宽巷模糊度协方差行覆盖
		}
	}
	//再列覆盖
	for(i=index;i<ComSatNum-1;i++)
	{
		for(j=0;j<RELDIMENSION;j++)
		{		
			Cov[j*RELDIMENSION+(i+12)]=Cov[j*RELDIMENSION+(i+13)];		//L1模糊度协方差列覆盖
			Cov[j*RELDIMENSION+(i+24)]=Cov[j*RELDIMENSION+(i+25)];		//宽巷模糊度协方差列覆盖
		}
	}
	//然后把最后一个卫星的模糊度协方差行列清0
	j=ComSatNum-1;
	for(i=0;i<RELDIMENSION;i++){
		Cov[(12+j)*RELDIMENSION+i]=0.0;		//L1模糊度协方差行清0
		Cov[(12+j)+RELDIMENSION*i]=0.0;		//L1模糊度协方差列清0
		Cov[(24+j)*RELDIMENSION+i]=0.0;		//宽巷模糊度协方差行清0		
		Cov[(24+j)+RELDIMENSION*i]=0.0;		//宽巷模糊度协方差列清0
	}
}
/****************************************************************************
函数名:deletePRN
目的:删除索引为index的卫星L1与宽巷模糊度、模糊度相关协方差
思路;通过从index开始，将后面一个覆盖前一个数据，最终将最后一个清0
输入:
		index		需要删除的卫星在卫星序列中的索引
		preComObs	前一个历元（当前滤波状态）的卫星序列及相关信息
		EkfState	滤波总结构体，为了调用其中的模糊度序列与协方差
输出:
		preComObs	删除目标卫星信息后的preComObs
		EkfState	删除目标卫星信息后的EkfState
****************************************************************************/
void deletePRN(int index,Common11Obs * preComObs,EKFSTATE * EkfState)
{
	int i,j;
	for(i=index;i<preComObs->ComSatNum-1;i++){
		copyOneSat11Obs(&preComObs->comobs[i],&preComObs->comobs[i+1]);
		EkfState->StateRel[12+i] = EkfState->StateRel[13+i];		//L1模糊度
		EkfState->StateRel[24+i] = EkfState->StateRel[25+i];		//宽巷模糊度
	}
	EkfState->StateRel[12+preComObs->ComSatNum-1] = 0.0;
	EkfState->StateRel[24+preComObs->ComSatNum-1] = 0.0;
	//清0最后一个模糊度的相关信息
	j=preComObs->ComSatNum-1;
	EmptyOneSat11Obs(&preComObs->comobs[j]);
	//删除相对应协方差Cov
	deletePRNCov(index,preComObs->ComSatNum,EkfState->CovRel);
	
	preComObs->ComSatNum=preComObs->ComSatNum-1;
}
//对同一历元的卫星序列进行调整，该函数是交换两卫星的序列
//注意要与模糊度和协方差同时调整
void SwapOneSat11Obs(OneSat11Obs * target,OneSat11Obs * src)
{
	OneSat11Obs tempSatObs;				//中间变量
	copyOneSat11Obs(&tempSatObs,target);
	copyOneSat11Obs(target,src);
	copyOneSat11Obs(src,&tempSatObs);
}
void getBpos(double * StateA,double * relState,double * StateB)/*前6个相加*/
{
	for (int i=0;i<6;i++)
	{
		StateB[i] = StateA[i] + relState[i];
	}
}
void getRelPos(double * StateA,double * StateB,double * relState)
{
	for (int i=0;i<6;i++)
	{
		relState[i] = StateB[i] - StateA[i];
	}
}
//对于发生周跳或者新出现的卫星初始化模糊度与协方差
void initialOneAmb(int index, double * Mat,EKFSTATE * EkfState,OneSat11Obs * obs)
{
	int i;
	double dx,dy,dz,rhoA,rhoB,dtr;
	double stateInECEFRel[6]={0.0};
	/************************************************************************/
	/* 初始化所需初值                                                       */
	/************************************************************************/
	dtr=EkfState->StateRel[6];
	MatrixMultiply( 3, 3, 3, 1, Mat, EkfState->StateA, EkfState->StateInECEFA );//每个观测值更新都要转到地固系中
	MatrixMultiply( 3, 3, 3, 1, Mat, EkfState->StateRel, stateInECEFRel );
	getBpos(EkfState->StateInECEFA,stateInECEFRel,EkfState->StateInECEFB);
	/************************************************************************/
	/* 初始化归零                                                           */
	/************************************************************************/
	for(i=0;i<RELDIMENSION;i++)
	{
		//L1模糊度
		EkfState->CovRel[(12+index)*RELDIMENSION+i]=0.0;
		EkfState->CovRel[(12+index)+RELDIMENSION*i]=0.0;
		//宽巷模糊度
		EkfState->CovRel[(24+index)*RELDIMENSION+i]=0.0;
		EkfState->CovRel[(24+index)+RELDIMENSION*i]=0.0;
	}
	//模糊度初始化
	dx	=	EkfState->StateInECEFA[0]-obs->satPos[0];
	dy	=	EkfState->StateInECEFA[1]-obs->satPos[1];
	dz	=	EkfState->StateInECEFA[2]-obs->satPos[2];
	rhoA=	sqrt(dx*dx+dy*dy+dz*dz);
	dx	=	EkfState->StateInECEFB[0]-obs->satPos[0];
	dy	=	EkfState->StateInECEFB[1]-obs->satPos[1];
	dz	=	EkfState->StateInECEFB[2]-obs->satPos[2];
	rhoB=	sqrt(dx*dx+dy*dy+dz*dz);
	EkfState->StateRel[24+index] = (obs->dL1-obs->dL2)-LEMTA_NL/LEMTA_WL*(obs->dP1/LEMTA_1+obs->dP2/LEMTA_2);

											
	EkfState->StateRel[12+index] = 1.0/(1.0-RATIO2)*obs->dL1*LEMTA_1-RATIO2/(1.0-RATIO2)*obs->dL2*LEMTA_2
								-(rhoB-rhoA)-dtr;
	
	/*EkfState->StateRel[12+index]				=1.0/(1.0-RATIO)*obs->dL1-RATIO/(1.0-RATIO)*obs->dL2
											   -(rhoB-rhoA)/LEMTA_NL-dtr/LEMTA_NL-EkfState->State[EKFDIMENSION+15+index]*RATIO/(1.0-RATIO);
	*/
	EkfState->CovRel[(12+index)+RELDIMENSION*(12+index)] = SIGL1;//初始化模糊度方差为周
	EkfState->CovRel[(24+index)+RELDIMENSION*(24+index)] = SIGMW;//初始化模糊度方差为周
	//周跳探测模糊度初始化
	obs->Nw_SD	=(obs->dL1-obs->dL2)
			-LEMTA_NL/LEMTA_WL*(obs->dP1/LEMTA_1+obs->dP2/LEMTA_2);
	obs->Nw_num_SD=1;
	obs->sigNw_SD=MW_SIG;
	//非差周跳探测模糊度初始化
	obs->Nw_A	=(obs->L11-obs->L12)
		-LEMTA_NL/LEMTA_WL*(obs->P11/LEMTA_1+obs->P12/LEMTA_2);
	obs->Nw_num_A=1;
	obs->sigNw_A=MW_SIG;

	obs->Nw_B	=(obs->L21-obs->L22)
		-LEMTA_NL/LEMTA_WL*(obs->P21/LEMTA_1+obs->P22/LEMTA_2);
	obs->Nw_num_B=1;
	obs->sigNw_B=MW_SIG;
}
//单差周跳探测主程序
int CySlipDetection_SD(double * Mat, Common11Obs * preComObs,Common11Obs * curComObs,EKFSTATE * EkfState)
{
	int i,j,isSearch,flag;
	int search[12]={0};
	char slipFilePath_GF[60]="OutFile\\cycleSlip.GF-SD";
	char slipFilePath_MW[60]="OutFile\\cycleSlip.MW-SD";
	char gfFilePath_SD[60]="OutFile\\gf-SD.txt";
	char mwFilePath_SD[60]="OutFile\\mw-SD.txt";
	FILE * cycleSlipFile_GF=fopen(slipFilePath_GF,"a+");
	FILE * cycleSlipFile_MW=fopen(slipFilePath_MW,"a+");
	FILE * gfFile_SD=fopen(gfFilePath_SD,"a+");
	FILE * mwFile_SD=fopen(mwFilePath_SD,"a+");
	int isHaveSlip_GF=0,isHaveSlip_MW=0;
	double tempLC=0.0,tempMW=0.0;
	/*周跳时间头*/
	double deltT=curComObs->Time.SecOfWeek-preComObs->Time.SecOfWeek+(curComObs->Time.Week-preComObs->Time.Week)*86400.0;
	if(preComObs->ComSatNum==0||fabs(deltT)>INTERVAL){
		copyCommon11Obs(&EkfState->CurComObs,curComObs);
		initialAllAmb(Mat,EkfState);
		fclose(cycleSlipFile_GF);
		fclose(cycleSlipFile_MW);
		fclose(gfFile_SD);
		fclose(mwFile_SD);
		return 0;//第一个历元，不用探测周跳
	}
	//TRACE("test");
	for(i=0;i<preComObs->ComSatNum;i++)
	{
		isSearch=0;
		for(j=0;j<curComObs->ComSatNum;j++)
		{
			if(preComObs->comobs[i].PRN==curComObs->comobs[j].PRN)
			{
				flag=0;
				isSearch=1;
				SwapOneSat11Obs(&curComObs->comobs[i],&curComObs->comobs[j]);
				search[i]=1;
				flag=GFDetection_SD(&preComObs->comobs[i],&curComObs->comobs[i], &tempLC);
				fprintf(gfFile_SD,"%12.3f %4d %12.3f\n",curComObs->Time.SecOfWeek,preComObs->comobs[i].PRN,tempLC);
				if(flag==1){
					fprintf(cycleSlipFile_GF,"%12.3f %4d\n",curComObs->Time.SecOfWeek,preComObs->comobs[i].PRN,tempLC);
					isHaveSlip_GF=1;
				}
				flag=MWDetection_SD(&preComObs->comobs[i],&curComObs->comobs[i], &tempMW);
				fprintf(mwFile_SD,"%12.3f %4d %12.3f\n",curComObs->Time.SecOfWeek,preComObs->comobs[i].PRN,tempMW);
				if(flag==1){
					fprintf(cycleSlipFile_MW,"%12.3f %4d %14.3f\n",curComObs->Time.SecOfWeek,preComObs->comobs[i].PRN,tempMW);
					isHaveSlip_MW=1;
				}
				break;
			}
		}
		if(isSearch==0)
		{
			deletePRN(i,preComObs,EkfState);//从preComObs->comobs[i]中删除该卫星，以及删除Ekf中该卫星的相关信息
			i--;
		}
	}
	for(i=0;i<curComObs->ComSatNum;i++)
	{
		if(curComObs->comobs[i].flag==1&&(curComObs->comobs[i].PRN==preComObs->comobs[i].PRN))
			initialOneAmb(i, Mat,EkfState,&curComObs->comobs[i]);
		if(search[i]==0){
			curComObs->comobs[i].flag=2;
			initialOneAmb(i, Mat,EkfState,&curComObs->comobs[i]);
		}
	}
	/*更新滤波中当前观测值*/
	copyCommon11Obs(preComObs,curComObs);

	fclose(cycleSlipFile_GF);
	fclose(cycleSlipFile_MW);
	fclose(gfFile_SD);
	fclose(mwFile_SD);
	return 1;
}
int CySlipDetection_LLI(double *Mat,Common11Obs * preComObs,Common11Obs * curComObs,EKFSTATE * EkfState){

	int i,j,isSearch,flag;
	int search[12]={0};
	int isHaveSlip_GF=0,isHaveSlip_MW=0;
    double tempMW,tempLC; 
	double deltT=curComObs->Time.SecOfWeek-preComObs->Time.SecOfWeek+(curComObs->Time.Week-preComObs->Time.Week)*86400.0;
	if(preComObs->ComSatNum==0||fabs(deltT)>INTERVAL){
		copyCommon11Obs(&EkfState->CurComObs,curComObs);
		initialAllAmb(Mat,EkfState);
		return 0;//第一个历元，不用探测周跳
	}
	//TRACE("test");
	for(i=0;i<preComObs->ComSatNum;i++)
	{
		isSearch=0;
		for(j=0;j<curComObs->ComSatNum;j++)
		{
			if(preComObs->comobs[i].PRN==curComObs->comobs[j].PRN)
			{
				flag=0;
				isSearch=1;
				SwapOneSat11Obs(&curComObs->comobs[i],&curComObs->comobs[j]);
				search[i]=1;
				if(curComObs->comobs[i].LLI_L1A==1||curComObs->comobs[i].LLI_L2A==1
					||curComObs->comobs[i].LLI_L1B==1||curComObs->comobs[i].LLI_L2B==1)
					curComObs->comobs[i].flag = -1;
				break;
			}
		}
		if(isSearch==0)
		{
			deletePRN(i,preComObs,EkfState);//从preComObs->comobs[i]中删除该卫星，以及删除Ekf中该卫星的相关信息
			i--;
		}
	}
	for(i=0;i<curComObs->ComSatNum;i++)
	{
		if(curComObs->comobs[i].flag==1&&(curComObs->comobs[i].PRN==preComObs->comobs[i].PRN))
			initialOneAmb(i, Mat,EkfState,&curComObs->comobs[i]);
		if(search[i]==0){
			curComObs->comobs[i].flag=2;
			initialOneAmb(i, Mat,EkfState,&curComObs->comobs[i]);
		}
	}
	/*更新滤波中当前观测值*/
	copyCommon11Obs(preComObs,curComObs);
	return 1;
}
int CySlipDetection_ND(double * Mat, Common11Obs * preComObs,Common11Obs * curComObs,EKFSTATE * EkfState)
{
	int i,j,isSearch,flag;
	int search[12]={0};
	char slipFilePath_GF[60]="OutFile\\cycleSlip.GF";
	char slipFilePath_MW[60]="OutFile\\cycleSlip.MW";
	char gfFilePath_ND[60]="OutFile\\gf-ND.txt";
	char mwFilePath_ND[60]="OutFile\\mw-ND.txt";
	FILE * cycleSlipFile_GF=fopen(slipFilePath_GF,"a+");
	FILE * cycleSlipFile_MW=fopen(slipFilePath_MW,"a+");
	FILE * gfFile_ND=fopen(gfFilePath_ND,"a+");
	FILE * mwFile_ND=fopen(mwFilePath_ND,"a+");
	int isHaveSlip_GF=0,isHaveSlip_MW=0;
	double tempLC[2]={0.0},tempMW[2]={0.0};
	/*周跳时间头*/
	double deltT=curComObs->Time.SecOfWeek-preComObs->Time.SecOfWeek+(curComObs->Time.Week-preComObs->Time.Week)*86400.0;
	if(preComObs->ComSatNum==0||fabs(deltT)>INTERVAL){
		copyCommon11Obs(&EkfState->CurComObs,curComObs);
		initialAllAmb(Mat,EkfState);
		fclose(cycleSlipFile_GF);
		fclose(cycleSlipFile_MW);
		fclose(gfFile_ND);
		fclose(mwFile_ND);
		return 0;//第一个历元，不用探测周跳
	}
	//TRACE("test");
	for(i=0;i<preComObs->ComSatNum;i++)
	{
		isSearch=0;
		for(j=0;j<curComObs->ComSatNum;j++)
		{
			if(preComObs->comobs[i].PRN==curComObs->comobs[j].PRN)
			{
				flag=0;
				isSearch=1;
				SwapOneSat11Obs(&curComObs->comobs[i],&curComObs->comobs[j]);
				search[i]=1;
				flag=GFDetection_ND(&preComObs->comobs[i],&curComObs->comobs[i], tempLC);
				fprintf(gfFile_ND,"%12.3f %4d %12.3f %12.3f\n",curComObs->Time.SecOfWeek,preComObs->comobs[i].PRN,tempLC[0],tempLC[1]);
				if(flag==1){
					fprintf(cycleSlipFile_GF,"%12.3f %4d %12.3f %12.3f\n",curComObs->Time.SecOfWeek,preComObs->comobs[i].PRN,tempLC[0],tempLC[1]);
				}
				flag=MWDetection_ND(&preComObs->comobs[i],&curComObs->comobs[i],tempMW);
				fprintf(mwFile_ND,"%12.3f %4d %12.3f %12.3f\n",curComObs->Time.SecOfWeek,preComObs->comobs[i].PRN,tempMW[0],tempMW[1]);
				if(flag==1){
					fprintf(cycleSlipFile_MW,"%12.3f %4d %12.3f %12.3f\n",curComObs->Time.SecOfWeek,preComObs->comobs[i].PRN,tempMW[0],tempMW[1]);
				}
				break;
			}
		}
		if(isSearch==0)
		{
			deletePRN(i,preComObs,EkfState);//从preComObs->comobs[i]中删除该卫星，以及删除Ekf中该卫星的相关信息
			i--;
		}
	}
	for(i=0;i<curComObs->ComSatNum;i++)
	{
		if(curComObs->comobs[i].flag==1&&(curComObs->comobs[i].PRN==preComObs->comobs[i].PRN))
			initialOneAmb(i, Mat,EkfState,&curComObs->comobs[i]);
		if(search[i]==0){
			curComObs->comobs[i].flag=2;
			initialOneAmb(i, Mat,EkfState,&curComObs->comobs[i]);
		}
	}
	/*更新滤波中当前观测值*/
	copyCommon11Obs(preComObs,curComObs);
	fclose(cycleSlipFile_GF);
	fclose(cycleSlipFile_MW);
	fclose(gfFile_ND);
	fclose(mwFile_ND);
	return 1;
}
/******************************************
COMBINED方法探测周跳
参考:秦显平博士论文，星载GPS低轨卫星定轨理论及方法研究
输入:
输出:
******************************************/
int combDetection()
{
	return 1;
}
/******************************************
电离层残差探测周跳法-单差观测值
输入:
		卫星i前一个历元dL1观测值，dL2观测值
		卫星i的当前历元dL1观测值，dL2观测值
输出：
		周跳标识flag
		-1:		数据不齐全，探测周跳失败
		 0:		未发现周跳
		 1:		发现周跳
******************************************/
int GFDetection_SD(OneSat11Obs * PreObs,OneSat11Obs * CurObs, double * dGF)
{
	double delta_GF;
	if(fabs(PreObs->dL1)>0.1&&fabs(PreObs->dL2)>0.1&&
		fabs(CurObs->dL1)>0.1&&fabs(CurObs->dL2)>0.1)
	{
		delta_GF=LEMTA_1*(CurObs->dL1-PreObs->dL1)+LEMTA_2*(PreObs->dL2-CurObs->dL2);//载波观测值单位为周
		*dGF=delta_GF;
		if(fabs(delta_GF)>GF_THRES)				//大于gf_thres周认为是周跳
		{
			CurObs->flag=1;
		}
		return CurObs->flag;
	}
	return -1;
}
/******************************************
电离层残差探测周跳法-非差观测值
输入:
		卫星i前一个历元dL1观测值，dL2观测值
		卫星i的当前历元dL1观测值，dL2观测值
输出：
		周跳标识flag
		-1:		数据不齐全，探测周跳失败
		 0:		未发现周跳
		 1:		发现周跳
******************************************/
int GFDetection_ND(OneSat11Obs * PreObs,OneSat11Obs * CurObs, double * dGF)
{
	double delta_GF;
	double preL1,preL2,curL1,curL2;
	double gflimitA=GF_THRES,gflimitB=GF_THRES;
	/*A星*/
	preL1=PreObs->L11;preL2=PreObs->L12;
	curL1=CurObs->L11;curL2=CurObs->L12;
	delta_GF=LEMTA_1*(curL1-preL1)+LEMTA_2*(preL2-curL2);//载波观测值单位为周
	dGF[0]=fabs(delta_GF);
	if(fabs(delta_GF)>gflimitA)				//大于gf_thres米认为是周跳
	{
		CurObs->flag=1;
	}
	
	/*B星*/
	preL1=PreObs->L21;preL2=PreObs->L22;
	curL1=CurObs->L21;curL2=CurObs->L22;
	delta_GF=LEMTA_1*(curL1-preL1)+LEMTA_2*(preL2-curL2);//载波观测值单位为周
	dGF[1]=fabs(delta_GF);
	if(fabs(delta_GF)>gflimitB)				//大于gf_thres周认为是周跳
	{
		CurObs->flag=1;
	}
	return CurObs->flag;
}
/***************************************
//MW法探测周跳-单差观测值
输入:
		卫星i的前一历元的dL1,dL2,dP1,dP2观测值
		卫星i的当前历元的dL1,dL2,dP1,dP2观测值
输出：
		周跳标识flag
		-1:		数据不齐全，探测周跳失败,生成共同卫星时对数据进行过要求，必须数据齐全
		 0:		未发现周跳
		 1:		发现周跳
****************************************/
int MWDetection_SD(OneSat11Obs * PreObs,OneSat11Obs * CurObs, double * tempMW)
{
	double temp,preSig2;
	temp=(CurObs->dL1-CurObs->dL2)-LEMTA_NL/LEMTA_WL*(CurObs->dP1/LEMTA_1+CurObs->dP2/LEMTA_2);
	*tempMW=fabs(temp-PreObs->Nw_SD);
	preSig2=PreObs->sigNw_SD*PreObs->sigNw_SD;
	if(fabs(temp-PreObs->Nw_SD)<(4*PreObs->sigNw_SD))
	{
		CurObs->Nw_num_SD=PreObs->Nw_num_SD+1;
		CurObs->Nw_SD=PreObs->Nw_SD+1.0/CurObs->Nw_num_SD*(temp-PreObs->Nw_SD);
		CurObs->sigNw_SD=sqrt(preSig2+1.0/CurObs->Nw_num_SD*((temp-PreObs->Nw_SD)*(temp-PreObs->Nw_SD)-preSig2));
		CurObs->flag	=0;
	}
	else
	{
		CurObs->flag	=1;
		CurObs->Nw_SD	=temp;
		CurObs->sigNw_SD=MW_SIG;
	}
	return CurObs->flag;
}
/***************************************
//MW法探测周跳-非差观测值
输入:
		卫星i的前一历元的dL1,dL2,dP1,dP2观测值
		卫星i的当前历元的dL1,dL2,dP1,dP2观测值
输出：
		周跳标识flag
		-1:		数据不齐全，探测周跳失败
		 0:		未发现周跳
		 1:		发现周跳
****************************************/
int MWDetection_ND(OneSat11Obs *  PreObs,OneSat11Obs *  CurObs,double * tempMW)
{
	double temp,preSig2;
	/*A星周跳探测*/
	temp=(CurObs->L11-CurObs->L12)-LEMTA_NL/LEMTA_WL*(CurObs->P11/LEMTA_1+CurObs->P12/LEMTA_2);
	preSig2=PreObs->sigNw_A*PreObs->sigNw_A;
	tempMW[0]=fabs(temp-PreObs->Nw_A);
	if(fabs(temp-PreObs->Nw_A)<(4*PreObs->sigNw_A))
	{
		CurObs->Nw_num_A=PreObs->Nw_num_A+1;
		CurObs->Nw_A=PreObs->Nw_A+1.0/CurObs->Nw_num_A*(temp-PreObs->Nw_A);
		CurObs->sigNw_A=sqrt(preSig2+1.0/CurObs->Nw_num_A*((temp-PreObs->Nw_A)*(temp-PreObs->Nw_A)-preSig2));
	}
	else
	{
		CurObs->flag	=1;
		CurObs->used	=-1;
		CurObs->Nw_A	=temp;
		CurObs->sigNw_A	=MW_SIG;
	}
	
	/*B星周跳探测*/
	
	temp=(CurObs->L21-CurObs->L22)-LEMTA_NL/LEMTA_WL*(CurObs->P21/LEMTA_1+CurObs->P22/LEMTA_2);
	preSig2=PreObs->sigNw_B*PreObs->sigNw_B;
	tempMW[1]=fabs(temp-PreObs->Nw_B);
	if(fabs(temp-PreObs->Nw_B)<(4*PreObs->sigNw_B))
	{
		CurObs->Nw_num_B=PreObs->Nw_num_B+1;
		CurObs->Nw_B=PreObs->Nw_B+1.0/CurObs->Nw_num_B*(temp-PreObs->Nw_B);
		CurObs->sigNw_B=sqrt(preSig2+1.0/CurObs->Nw_num_B*((temp-PreObs->Nw_B)*(temp-PreObs->Nw_B)-preSig2));
	}
	else
	{
		CurObs->flag	=1;
		CurObs->used	=-1;
		CurObs->Nw_B	=temp;
		CurObs->sigNw_B	=MW_SIG;
	}
	
	return CurObs->flag;
}
//利用高度角计算电离层投影函数
double getMap(double radEle)
{
	double map=0.0;
	//Lear模型
	double sinE=sin(radEle);
	map=2.037/(sqrt(sinE*sinE+0.076)+sinE);
	return map;
}

//************************************
// 函数名称:  relEKFilter_LP
// 函数说明:  伪距相对测量更新
// 作者:      RJ
// 时间:   	  2016/11/25
// 返回值:    void
// 参数:      double Mat[9]
// 参数:      EKFSTATE * KFState
//************************************
int relEKFilter_PC(double Mat[9],EKFSTATE * KFState)
{
	int    i, j, k;
	double stateInECEFRel[6]={0.0};
	double H[RELDIMENSION];	/* 观测方程线性化系数向量 */
	double RangeA,RangeB;				/* 接收机与导航卫星之间的距离计算值 */
	double dRange,dtr;				/*单差近似值*/
	double dPos_A[3],dPos_B[3];				/* 接收机坐标与导航卫星位置之差 */
	double Adtr,AECF[3], ARangeB; /*验前状态参数，用于状态估计信息使用*/
	double OC_LP;
	double OC_LP0[12],OC_LP1[12];
	int    satNumWhole=0;
	FILE * fPC;
	fPC=fopen( "OutFile\\LPOC.txt", "a+" );
	fprintf( fPC, "%12.3f\n", KFState->Time.SecOfWeek );
	//--------------------------------------------------------------------------------
	//先验残差
	//CopyArray( 2*EKFDIMENSION, ApriState, KFState->State );
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateA, KFState->StateInECEFA );//每个观测值更新都要转到地固系中
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateRel, stateInECEFRel );
	getBpos( KFState->StateInECEFA, stateInECEFRel, KFState->StateInECEFB);
	//验前状态，验前残差用
	CopyArray(3,AECF,KFState->StateInECEFB);
	Adtr=KFState->StateRel[6];
	//共同卫星数
	for(i=0;i<KFState->CurComObs.ComSatNum;i++)
	{
		if(KFState->CurComObs.comobs[i].used!=1)
			continue;
		fprintf( fPC, "%12d",KFState->CurComObs.comobs[i].PRN );
		satNumWhole++;
	}
	fprintf(fPC,"\n");
	//--------------------------------------------------------------------------------
	if ( satNumWhole<4 )
	{
		return 0;
	}
	
	for(i=0;i<KFState->CurComObs.ComSatNum;i++)
	{
		MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateRel, stateInECEFRel );
		getBpos( KFState->StateInECEFA, stateInECEFRel, KFState->StateInECEFB);
		for( j=0; j<RELDIMENSION; j++ )
		{
			H[j]	= 0.0;
		}
		//假设相对测量更新前，有相对定位或者卡尔曼新息粗差探测手段
		if(KFState->CurComObs.comobs[i].used == 1 )
		{
			RangeA=0.0;RangeB=0.0;
			for(k=0;k<3;k++){
				dPos_A[k]=(KFState->StateInECEFA[k]-KFState->CurComObs.comobs[i].satPos[k]);
				RangeA+=dPos_A[k]*dPos_A[k];
				dPos_B[k]=(KFState->StateInECEFB[k]-KFState->CurComObs.comobs[i].satPos[k]);
				RangeB+=dPos_B[k]*dPos_B[k];
			}

			RangeA	=	sqrt(RangeA);
			RangeB	=	sqrt(RangeB);
			dRange	=	RangeB-RangeA;
			dtr		=	KFState->StateRel[6];
			OC_LP	=	1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP1 - RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP2
							-dRange-dtr;
			//位置观测矩阵，观测值地固系下的
			H[0]	= 1.0 * ( dPos_B[0] * Mat[0] + dPos_B[1] * Mat[3] + dPos_B[2] * Mat[6] )/RangeB;
			H[1]	= 1.0 * ( dPos_B[0] * Mat[1] + dPos_B[1] * Mat[4] + dPos_B[2] * Mat[7] )/RangeB;
			H[2]	= 1.0 * ( dPos_B[0] * Mat[2] + dPos_B[1] * Mat[5] + dPos_B[2] * Mat[8] )/RangeB;
				//钟差观测矩阵
			H[6]	= 1.0;
			EKFMeasureUpdateRel(OC_LP, 1, H, KFState );	//Lp测量更新
		}
	}
	KFState->ApriSigma[1] = 0.0;
	KFState->PostSigma[1] = 0.0;
	//--------------------------------------------------------------------------------
	//残差
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateRel, stateInECEFRel );
	getBpos( KFState->StateInECEFA, stateInECEFRel, KFState->StateInECEFB);
	satNumWhole=0;
	for(i=0;i<KFState->CurComObs.ComSatNum;i++)
	{
		if(KFState->CurComObs.comobs[i].used!=1)
			continue;

		RangeA=0.0;RangeB=0.0;ARangeB=0.0;

		for(k=0;k<3;k++){
			RangeA = RangeA+pow(KFState->StateInECEFA[k]-KFState->CurComObs.comobs[i].satPos[k],2.0);
			RangeB = RangeB+pow(KFState->StateInECEFB[k]-KFState->CurComObs.comobs[i].satPos[k],2.0);
			ARangeB= ARangeB+pow(AECF[k]-KFState->CurComObs.comobs[i].satPos[k],2.0);
		}

		RangeA=sqrt(RangeA);
		RangeB=sqrt(RangeB);
		ARangeB=sqrt(ARangeB);
		dtr		=KFState->StateRel[6];
		OC_LP0[satNumWhole]	=1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP1-
			RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP2
			-(ARangeB-RangeA)-Adtr;
		OC_LP1[satNumWhole]	=1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP1-
			RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP2
			-(RangeB-RangeA)-dtr;
		KFState->ApriSigma[1] = KFState->ApriSigma[1] + pow( OC_LP0[satNumWhole], 2.0);
		KFState->PostSigma[1] = KFState->PostSigma[1] + pow( OC_LP1[satNumWhole], 2.0 );

		satNumWhole++;
	}
	if( satNumWhole > 0 )
	{
		KFState->ApriSigma[1] = sqrt( KFState->ApriSigma[1] / satNumWhole );
		KFState->PostSigma[1] = sqrt( KFState->PostSigma[1] / satNumWhole );
	}
	else
	{
		KFState->ApriSigma[1] = 999.0;
		KFState->PostSigma[1] = 999.0;
	}
	outRes(satNumWhole,fPC,OC_LP0,OC_LP1,KFState);
	outResEKFPC(satNumWhole,OC_LP0,OC_LP1,KFState);
	KFState->SatNumUsed[1] = satNumWhole;
	KFState->sigPosPC=KFState->PostSigma[1];
	EKFConvergencyCheck(1,KFState);
	return satNumWhole;
	//--------------------------------------------------------------------------------
}

//************************************
// 函数名称:  relEKFilter_LC
// 函数说明:  载波相对测量更新
// 作者:      RJ
// 时间:   	  2016/11/25
// 返回值:    void
// 参数:      double Mat[9]
// 参数:      EKFSTATE * KFState
//************************************
int relEKFilter_LC(double Mat[9],EKFSTATE * KFState)
{
	int    i, j, k;
	double stateInECEFRel[6]={0.0};
	double AL1[12]={0.0};
	double AW[12]={0.0};	
	double H[RELDIMENSION];	/* 观测方程线性化系数向量 */
	double RangeA,RangeB;				/* 接收机与导航卫星之间的距离计算值 */
	double dRange,dtr;				/*单差近似值*/
	double dPos_A[3],dPos_B[3];				/* 接收机坐标与导航卫星位置之差 */
	double AECF[3], ARangeB,Adtr; /*验前状态参数，用于状态估计信息使用*/
	double OC_Lc;
	double OC_LC0[12],OC_LC1[12];
	double LC_PC0=0.0,LC_PC1=0.0;//滤波前后无电离层载波减去无电离层伪距，均值为0说明周跳和模糊度处理正确
	FILE   *fLC=fopen( "OutFile\\LCOC.txt", "a+" );
	char   lc_str[60];
	int    satNumWhole=0;
	
	fprintf( fLC, "%12.3f\n", KFState->Time.SecOfWeek );

	//--------------------------------------------------------------------------------
	//先验残差
	//CopyArray( 2*EKFDIMENSION, ApriState, KFState->State );
	//traceMat(RELDIMENSION,RELDIMENSION,KFState->CovRel);
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateA, KFState->StateInECEFA );//每个观测值更新都要转到地固系中
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateRel, stateInECEFRel );
	getBpos( KFState->StateInECEFA, stateInECEFRel, KFState->StateInECEFB);
	CopyArray(3,AECF,KFState->StateInECEFB);
	Adtr=KFState->StateRel[6];
	for (i=0;i<12;i++)
	{
		AL1[i] = KFState->StateRel[12+i];
		AW[i] = KFState->StateRel[24+i];
	}

	for(i=0;i<KFState->CurComObs.ComSatNum;i++)
	{
		if(KFState->CurComObs.comobs[i].flag!=0||KFState->CurComObs.comobs[i].used!=1)
			continue;
		fprintf(fLC,"%12d",KFState->CurComObs.comobs[i].PRN);
		satNumWhole++;
	}
	//--------------------------------------------------------------------------------
	fprintf( fLC, "\n" );
	if ( satNumWhole<4 )
	{
		fclose(fLC);
		return 0;
	}
	
	for(i=0;i<KFState->CurComObs.ComSatNum;i++)
	{
		for( j=0; j<RELDIMENSION; j++ )
		{
			H[j]	= 0.0;
		}
		//假设相对测量更新前，有相对定位或者卡尔曼新息粗差探测手段
		if(KFState->CurComObs.comobs[i].flag==0&&KFState->CurComObs.comobs[i].used==1)
		{
			MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateRel, stateInECEFRel );
			getBpos( KFState->StateInECEFA, stateInECEFRel, KFState->StateInECEFB);
			RangeA=0.0;RangeB=0.0;
			for(k=0;k<3;k++){
				dPos_A[k]=(KFState->StateInECEFA[k]-KFState->CurComObs.comobs[i].satPos[k]);
				RangeA+=dPos_A[k]*dPos_A[k];
				dPos_B[k]=(KFState->StateInECEFB[k]-KFState->CurComObs.comobs[i].satPos[k]);
				RangeB+=dPos_B[k]*dPos_B[k];
			}

			RangeA	=	sqrt(RangeA);
			RangeB	=	sqrt(RangeB);
			dRange	=	RangeB-RangeA;
			dtr		=	KFState->StateRel[6];
			OC_Lc	=	1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL1*LEMTA_1 - RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL2*LEMTA_2
						-dRange-dtr-KFState->StateRel[12+i];//电离层模糊度
			//LEMTA_NL*KFState->State[15+i]-Hwl*KFState->State[15+i+EKFDIMENSION];
			/*
			Obs_MW	=LEMTA_WL*(KFState->CurComObs.comobs[i].dL1			-	KFState->CurComObs.comobs[i].dL2)-
					 LEMTA_NL*(KFState->CurComObs.comobs[i].dP1/LEMTA_1	+	KFState->CurComObs.comobs[i].dP2/LEMTA_2);
			OC_MW	=Obs_MW-LEMTA_WL*KFState->StateRel[24+i];
			OC_G1	=1.0/2*(LEMTA_1*KFState->CurComObs.comobs[i].dL1+KFState->CurComObs.comobs[i].dP1)
						-dRange-dtr-1.0/2*LEMTA_1*KFState->StateRel[12+i];
			Obs_MW_OC=Obs_MW/LEMTA_WL-KFState->StateRel[24+i];
			*/
			//位置观测矩阵，观测值地固系下的
			H[0]	= 1.0 * ( dPos_B[0] * Mat[0] + dPos_B[1] * Mat[3] + dPos_B[2] * Mat[6] )/RangeB;
			H[1]	= 1.0 * ( dPos_B[0] * Mat[1] + dPos_B[1] * Mat[4] + dPos_B[2] * Mat[7] )/RangeB;
			H[2]	= 1.0 * ( dPos_B[0] * Mat[2] + dPos_B[1] * Mat[5] + dPos_B[2] * Mat[8] )/RangeB;
			//Lc模糊度观测矩阵
				//钟差观测矩阵
			H[6]	= 1.0;
			H[12+i]	= 1.0;
			//H[12+i] = LEMTA_NL;
			//H[24+i] = Hwl;	//宽巷模糊度观测系数
			EKFMeasureUpdateRel(OC_Lc, Me_dLC, H, KFState );//Lc测量更新
			/*
			H[12+i]				= LEMTA_1/2.0;
			H[24+i]				= 0.0;
			EKFMeasureUpdate( OC_G1, 0.71, H, KFState );//G1测量更新
			H[24+i]=1.0*LEMTA_WL;
			EKFMeasureUpdate( OC_MW, 1.01, H, KFState );//MW测量更新,注意这里单位是m，而不是周
			*/
		}
	}
	KFState->ApriSigma[1] = 0.0;
	KFState->PostSigma[1] = 0.0;
	
	//--------------------------------------------------------------------------------
	//残差
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateRel, stateInECEFRel );
	getBpos( KFState->StateInECEFA, stateInECEFRel, KFState->StateInECEFB);
	satNumWhole=0;
	for(i=0;i<KFState->CurComObs.ComSatNum;i++)
	{
		if(KFState->CurComObs.comobs[i].flag!=0||KFState->CurComObs.comobs[i].used!=1)
			continue;

		RangeA=0.0;RangeB=0.0;ARangeB=0.0;

		for(k=0;k<3;k++){
			RangeA = RangeA+pow(KFState->StateInECEFA[k]-KFState->CurComObs.comobs[i].satPos[k],2.0);
			RangeB = RangeB+pow(KFState->StateInECEFB[k]-KFState->CurComObs.comobs[i].satPos[k],2.0);
			ARangeB= ARangeB+pow(AECF[k]-KFState->CurComObs.comobs[i].satPos[k],2.0);
		}

		RangeA = sqrt(RangeA);
		RangeB = sqrt(RangeB);
		ARangeB= sqrt(ARangeB);
		dtr	   = KFState->StateRel[6];
		OC_LC0[satNumWhole] = 1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL1*LEMTA_1 - RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL2*LEMTA_2
			-(ARangeB-RangeA)-Adtr-AL1[i];//电离层模糊度
		OC_LC1[satNumWhole] = 1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL1*LEMTA_1 - RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL2*LEMTA_2
							- (RangeB-RangeA)-dtr-KFState->StateRel[12+i];//电离层模糊度
		LC_PC0+= 1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL1*LEMTA_1 - RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL2*LEMTA_2
			-AL1[i]-(1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP1-RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP2);
		LC_PC1+= 1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL1*LEMTA_1 - RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dL2*LEMTA_2
			-KFState->StateRel[12+i]-(1.0/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP1-RATIO2/(1.0-RATIO2)*KFState->CurComObs.comobs[i].dP2);
		//LEMTA_NL*KFState->State[12+i]-Hwl*KFState->State[24+i];
		KFState->ApriSigma[1] = KFState->ApriSigma[1] + pow( OC_LC0[satNumWhole], 2.0);
		KFState->PostSigma[1] = KFState->PostSigma[1] + pow( OC_LC1[satNumWhole], 2.0 );
		satNumWhole++;
	}
	//--------------------------------------------------------------------------------
	if( satNumWhole > 0 )
	{
		KFState->ApriSigma[1] = sqrt( KFState->ApriSigma[1] / satNumWhole );
		KFState->PostSigma[1] = sqrt( KFState->PostSigma[1] / satNumWhole );
	}
	else
	{
		KFState->ApriSigma[1] = 999.0;
		KFState->PostSigma[1] = 999.0;
	}
	//输出评价
	outResAmbSD(KFState->Time.SecOfWeek,LC_PC0,LC_PC1);
	outRes(satNumWhole,fLC,OC_LC0,OC_LC1,KFState);
	outResEKFLC(satNumWhole,OC_LC0,OC_LC1,KFState);
	KFState->SatNumUsed[1] = satNumWhole;
	KFState->sigPosLC=KFState->PostSigma[1];
	EKFConvergencyCheck(1,KFState);
	return satNumWhole;
	
}

int PPRTOD_LC(Common11Obs * CurObs,EKFSTATE * KFState, PPRESULT * Result, double * Mat)
{
	int i, j, k, Iterator;                     /*  Iterator为单点定位迭代次数  */
	int SatNumUsed,GPSSatNum;      /*  单点定位使用的卫星数, GPS卫星数 */
	short Prn;
	short PosKind;            /* GPS:0, GLO:1, GPS&GLO:2  */		
	int    PRValid;                            /* 观测值的伪距可用性和伪距观测值 */
	double PRange;
	double AprioriPosRel[3], AprioriClk[2];
	double RangeA,RangeB, Ion, Height;                                 /*  接收机与导航卫星之间的距离  */
	double BPos[MAXCHANNUM];                      /* 伪距或多普勒的观测值-计算值 */
	double MeasA[MAXCHANNUM*5], MeasAT[MAXCHANNUM*5];          /* 观测矩阵 [MAXCHANNUM][5] */
	double Weight[MAXCHANNUM];                               /* 权矩阵, 单位阵只取对角线元素 */
	double Qvv[MAXCHANNUM];                               /* 观测值改正数的斜因数阵对角线元素 */
	double QLL[MAXCHANNUM*MAXCHANNUM]={0.0};			/*平差后观测值的协因数阵*/
	double ATA[25], InvATA[25];
	double ATB[5];
	double dPos[5];               /* 定位结果  */
	double Residual[MAXCHANNUM] = {0.0};  /* 定位残差  */
	int outlier=1;							/*第一次探测粗差*/
	double	OrbitAccuracy=0.0;
	PosKind=0;
	
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateA, KFState->StateInECEFA );//每个观测值更新都要转到地固系中
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateRel, AprioriPosRel );
	CopyArray(2,AprioriClk,&KFState->StateRel[6]);
	for(i=0;i<3;i++)
	{
		OrbitAccuracy = OrbitAccuracy + pow( KFState->CovRel[i*RELDIMENSION+i], 2.0 );
	}
	OrbitAccuracy=sqrt(OrbitAccuracy);
	Iterator = 0;
	
		do{
				SatNumUsed = 0;
				GPSSatNum  = 0;
	
				for( i=0; i<CurObs->ComSatNum; i++ )
				{
					if(CurObs->comobs[i].flag==0)
					{
					
						RangeA = 0.0;
						RangeB = 0.0;
						for( k=0; k<3; k++ )
						{
							RangeA = RangeA + pow( KFState->StateInECEFA[k]-CurObs->comobs[i].satPos[k], 2.0 );
							RangeB = RangeB + pow( KFState->StateInECEFA[k]+AprioriPosRel[k]-CurObs->comobs[i].satPos[k], 2.0 );
						}
	
						RangeA = sqrt( RangeA );
						RangeB = sqrt( RangeB );
						MeasA[SatNumUsed*5+0] = 1.0*( AprioriPosRel[0] + KFState->StateInECEFA[0] - CurObs->comobs[i].satPos[0] )/RangeB;
						MeasA[SatNumUsed*5+1] = 1.0*( AprioriPosRel[1] + KFState->StateInECEFA[1] - CurObs->comobs[i].satPos[1] )/RangeB;
						MeasA[SatNumUsed*5+2] = 1.0*( AprioriPosRel[2] + KFState->StateInECEFA[2]- CurObs->comobs[i].satPos[2] )/RangeB;
						MeasA[SatNumUsed*5+3] = 1.0;             /* 接收机GPS系统钟差系数 */
						MeasA[SatNumUsed*5+4] = 0.0;             /* 接收机GLONASS系统钟差系数 */
						Weight[SatNumUsed]    = 1.0;  // * pow( sin(SatMidInfo[i].Elevation), 2.0 );
	
						BPos[SatNumUsed] =1.0/(1.0-RATIO2)*CurObs->comobs[i].dL1*LEMTA_1-RATIO2/(1.0-RATIO2)*CurObs->comobs[i].dL2*LEMTA_2
							-(RangeB-RangeA)-AprioriClk[0]-LEMTA_NL*KFState->StateRel[12+i]-Hwl*KFState->StateRel[24+i];
						Result->SatList[SatNumUsed].Status=1;
						Result->SatList[SatNumUsed].index =i;
						CurObs->comobs[i].index_satlist=SatNumUsed;
						SatNumUsed++;
					
					}
				}
				/*  组成误差方程  */
	
				for( k=SatNumUsed; k<MAXCHANNUM; k++)  /* 多余项清零 */
				{
					BPos[k] = 0.0;
					Weight[k] = 0.0;
					for( j=0; j<5; j++ )
					{
						MeasA[k*5+j] = 0.0;
					}
				}
	
				if( SatNumUsed>=4 )  /* 不考虑GLONASS卫星 */
				{
					PosKind = 0;
					/*  将观测矩阵变换为MAXCHANNUM*4, 消除GLONASS卫星钟差参数 */
	
					for( i=0; i<SatNumUsed; i++ )
					{
						MeasA[i*4+0] = MeasA[i*5+0];
						MeasA[i*4+1] = MeasA[i*5+1];
						MeasA[i*4+2] = MeasA[i*5+2];
						MeasA[i*4+3] = MeasA[i*5+3];
					}
	
					MatrixTranspose( SatNumUsed, 4, MeasA, MeasAT );
					MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 4, MeasAT, Weight, MeasA, ATA );
					MatrixInv( 4, ATA, InvATA );
	
					MatrixMultiply3( 4, SatNumUsed, SatNumUsed, 1, MeasAT, Weight, BPos, ATB );
					MatrixMultiply( 4, 4, 4, 1, InvATA, ATB, dPos );
	
					AprioriClk[0] = AprioriClk[0] + dPos[3];
	
				}
				else
				{
					trace("\nCannot position because of only %2d satellites at %10.3f epoch.\n",SatNumUsed, CurObs->Time.SecOfWeek );
					break;
				}
	
				for( k=0; k<3; k++ )
				{
					AprioriPosRel[k] = AprioriPosRel[k] + dPos[k];
				}
				Iterator++;
				Result->Iterator = Iterator;
				Result->Coverage = VectDot( 3, 3, dPos, dPos );
	
			}while ( (Result->Coverage>1E-6) && (Iterator < 10) );
			//循环结果计算高度角
			Result->SatNumUsed = SatNumUsed;
			CopyArray(3,Result->Position,AprioriPosRel);
			CopyArray( 2, Result->RcvClkOft, AprioriClk );
	
	
			/* 单独定位结果整理输出  */
			/* 计算观测值残差 */
	
			MatrixMultiply( SatNumUsed, 4, 4, 1, MeasA, dPos, Residual );//
			for( i=0; i<SatNumUsed; i++ )
			{
				Residual[i] = Residual[i] - BPos[i];
			}
	
			/* 计算观测值改正数的协因数阵 */
			ComputeQvv( SatNumUsed, 4, MeasA, InvATA, MeasAT, Weight, Qvv );
			/* 计算平差后观测值的协因数阵*/
			ComputeQLL( SatNumUsed, 4, MeasA, InvATA, MeasAT, Weight, QLL);
			Result->PDOP = sqrt( InvATA[0] + InvATA[5] + InvATA[10] );
			Result->GDOP = sqrt( InvATA[0] + InvATA[5] + InvATA[10] + InvATA[15] );
			Result->HDOP = sqrt( InvATA[0] + InvATA[5] );
			Result->VDOP = sqrt( InvATA[10] );
			Result->TDOP = sqrt( InvATA[15] );
	
			if( SatNumUsed > 5 ) /* 多于5颗卫星才计算定位中误差, 增强可靠性 */
			{
				Result->SigmaPos = 0.0;
	
				for( j=0, i=0; i<CurObs->ComSatNum; i++ )
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
	
			if( SatNumUsed > 5 )			//多余观测值大于等于2时才能进行粗差检测，而W检验只有当观测值独立时才能进行
			{
				//CheckPostResidual_W( CurObs->ComSatNum, Me_dLC, Residual, Qvv , Result );//其中数值分别为w检验限值与观测值噪声
			}
			/*RJ-2016-9-20
			**该模块是为了把单点定位中验出粗差的卫星作为不可用
			*/
			for(i=0;i<CurObs->ComSatNum;i++)
			{
				if(CurObs->comobs[i].flag==0&&Result->SatList[CurObs->comobs[i].index_satlist].Status<0){
					CurObs->comobs[i].flag=1;
					//重新初始化模糊度
				}
			}

	if( Result->Coverage < 1E-6 && Result->Iterator < 5 &&
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
int PPRTOD_PC(Common11Obs * CurObs,EKFSTATE * KFState, PPRESULT * Result, double * Mat)
{
	int i, j, k,SatNumUsed=0;                     /*  Iterator为单点定位迭代次数  */
	double StateInECEFA[3];
	double AprioriPosRel[3], AprioriClk;
	double apriRes=0.0,posRes=0.0;
	double RangeA,RangeB;
	int    bestBan=999;
	int    PRN[MAXCHANNUM]={999};//只有第一个为999，后面的都是0
	int    banPRN[MAXCHANNUM]={999};
	double satPos[3*MAXCHANNUM]={0.0};
	double obs[MAXCHANNUM]={0.0};
	int    valFlag=0;
	FILE * apriRes_file,*posRes_file;
	for(i=0;i<MAXCHANNUM;i++)
	{
		PRN[i] = 999;
		banPRN[i] = 999;
	}
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateA, KFState->StateInECEFA );//每个观测值更新都要转到地固系中
	MatrixMultiply( 3, 3, 3, 1, Mat, KFState->StateRel, AprioriPosRel );
	AprioriClk=KFState->StateRel[6];
	for (i=0;i<CurObs->ComSatNum;i++)
	{
		RangeA=0.0;RangeB=0.0;
		if (CurObs->comobs[i].used==1)
		{
			for (j=0;j<3;j++)
			{
				satPos[SatNumUsed*3+j]=CurObs->comobs[i].satPos[j];
				RangeA+=pow(KFState->StateInECEFA[j]-CurObs->comobs[i].satPos[j],2.0);
				RangeB+=pow(KFState->StateInECEFA[j]+AprioriPosRel[j]-CurObs->comobs[i].satPos[j],2.0);
			}
			RangeA=sqrt(RangeA);
			RangeB=sqrt(RangeB);
			obs[SatNumUsed]=1.0/(1.0-RATIO2)*CurObs->comobs[i].dP1-RATIO2/(1.0-RATIO2)*CurObs->comobs[i].dP2;
			apriRes=apriRes+obs[SatNumUsed]-(RangeB-RangeA)-AprioriClk;
			PRN[SatNumUsed]=CurObs->comobs[i].PRN;
			SatNumUsed++;
		}
	}
	apriRes=apriRes/SatNumUsed;
	apriRes_file=fopen("OutFile\\apriRes.txt","a+");
	fprintf(apriRes_file,"%14.3f%14.3f\n",KFState->Time.SecOfWeek,apriRes);
	fclose(apriRes_file);
	if (SatNumUsed<5)
	{
		Result->IsSuccess = false;
		Result->SigmaPos=999.0;
		Result->SatNumUsed=0;
		return 0;
	}
	/*LEGE探测多维粗差*/
	Result->SigmaPos=LSRel_PC(SatNumUsed,PRN,banPRN,satPos,KFState->StateInECEFA,obs,AprioriPosRel,&AprioriClk, Result);
	valFlag=valResult(Result->SigmaPos,5,SatNumUsed,4);
	if (valFlag==1)
	{
		Result->SatNumUsed=SatNumUsed;
	}
	else
	{
		Raim_fde1(&SatNumUsed,PRN,banPRN,satPos,KFState->StateInECEFA,obs,AprioriPosRel,&AprioriClk,Result);		
	}
	valFlag=valResult(Result->SigmaPos,5,SatNumUsed,4);
	CopyArray(3,Result->Position,AprioriPosRel);
	Result->RcvClkOft[0]=AprioriClk;
	/*RJ-2016-9-20
	**该模块是为了把单点定位中验出粗差的卫星作为不可用
	*/
	for(i=0;i<CurObs->ComSatNum;i++)
	{
		for (j=0;j<MAXCHANNUM;j++)
		{
			if(CurObs->comobs[i].PRN==banPRN[j])
				CurObs->comobs[i].used=-2;
		}			
	}
	/*计算验后方差*/
	SatNumUsed=0;
	for (i=0;i<CurObs->ComSatNum;i++)
	{
		RangeA=0.0;RangeB=0.0;
		if (CurObs->comobs[i].used==1)
		{
			for (j=0;j<3;j++)
			{
				satPos[SatNumUsed*3+j]=CurObs->comobs[i].satPos[j];
				RangeA+=pow(KFState->StateInECEFA[j]-CurObs->comobs[i].satPos[j],2.0);
				RangeB+=pow(KFState->StateInECEFA[j]+AprioriPosRel[j]-CurObs->comobs[i].satPos[j],2.0);
			}
			RangeA=sqrt(RangeA);
			RangeB=sqrt(RangeB);
			obs[SatNumUsed]=1.0/(1.0-RATIO2)*CurObs->comobs[i].dP1-RATIO2/(1.0-RATIO2)*CurObs->comobs[i].dP2;
			posRes=posRes+obs[SatNumUsed]-(RangeB-RangeA)-AprioriClk;
			PRN[SatNumUsed]=CurObs->comobs[i].PRN;
			SatNumUsed++;
		}
	}
	posRes=posRes/SatNumUsed;
	posRes_file=fopen("OutFile\\posRes.txt","a+");
	fprintf(posRes_file,"%14.3f%14.3f\n",KFState->Time.SecOfWeek,posRes);
	fclose(posRes_file);
	if( Result->Coverage < 1E-6 && Result->Iterator < 5 &&
		Result->PDOP < 8.0 && Result->SigmaPos < 5.0 && SatNumUsed > 4 )//将至少有6颗卫星改为5颗
	{
		Result->IsSuccess = true;
	}
	else
	{
		Result->IsSuccess = false;
	}

	return SatNumUsed;
}
//************************************
// Method:    EKFTimeUpdateAB
// FullName:  EKFTimeUpdateAB
// Access:    public 
// Returns:   void
// Qualifier: ABTimeUpdate函数下面的子函数,用于滤波器时间更新
// Parameter: double OIState[108]
// Parameter: EKFSTATE * KFState
//************************************
void EKFTimeUpdateRel( double OIState[54], EKFSTATE* KFState )
{
	//对于cov矩阵  |A C|
	//             |D B|
	
	int i,j;
	double StateB[RELDIMENSION]={0.0};
	double STM[RELDIMENSION*RELDIMENSION] = {0.0};
	double Q[RELDIMENSION*RELDIMENSION] = {0.0};
	double TmpMat[RELDIMENSION*RELDIMENSION], STMT[RELDIMENSION*RELDIMENSION];
	double KFTmp[RELDIMENSION*RELDIMENSION];
	//下面两行是否有问题
	CopyArray(RELDIMENSION,StateB,KFState->StateRel);
	//除速度和位置外赋值给StateB
	CopyArray(6,StateB,OIState);						//直接利用积分器里面的
	/************************************************************************/
	/* 时间更新前的协方差                                                   */
	/************************************************************************/
	for(i=0;i<RELDIMENSION;i++)
	{
		for(j=0;j<RELDIMENSION;j++)
		{
			KFTmp[i*RELDIMENSION+j]=KFState->CovRel[i*RELDIMENSION+j];
		}
	}
	/************************************************************************/
	/* B星+宽巷模糊度的过程噪声，注意StateB是为了计算里面的坐标转移矩阵     */
	/************************************************************************/
	UdProcessNoiseCovRTN4rel( KFState->Tao[1], KFState->Step, KFState->Sigma[1], StateB, Q);
	
	FormStateTransMatrixFromDynRTN4rel( KFState->Tao[1], KFState->Step,StateB, &OIState[6], STM );
	
	MatrixMultiply( RELDIMENSION, RELDIMENSION, RELDIMENSION,RELDIMENSION, 
		STM, KFTmp, TmpMat );
	MatrixTranspose( RELDIMENSION, RELDIMENSION, STM, STMT );
	MatrixMultiply( RELDIMENSION, RELDIMENSION, RELDIMENSION, RELDIMENSION,
		TmpMat, STMT, KFTmp );
	MatrixAddition2( RELDIMENSION, RELDIMENSION, Q, KFTmp ); 
	for(i=0;i<RELDIMENSION;i++)
	{
		for(j=0;j<RELDIMENSION;j++)
		{
			KFState->CovRel[i*RELDIMENSION+j]=KFTmp[i*RELDIMENSION+j];
		}
	}
	/************************************************************************/
	/* 协方差更新完毕，接下来是补偿加速度对状态的补偿                       */
	/************************************************************************/
	DMCrtn4Rel( KFState->Tao[1], KFState->Step, StateB );			//B星的补偿
	//因为上述状态转移只是来求方差转移，位置速度转移是单独积分的，因此这里需要补偿
	CopyArray(RELDIMENSION,KFState->StateRel,StateB);//因为上述对相对状态的除速度位置外均进行了更新
	CopyArray(RELDIMENSION,KFState->StateB,StateB);//因为上述对相对状态的除速度位置外均进行了更新
}
void outResAmbSD(const double sow,const double LC_PC0,const double LC_PC1)
{
	FILE * fp;
	fp=fopen("OutFile\\LC-PC","a+");
	fprintf(fp,"%14.3f%14.3f%14.3f\n",sow,LC_PC0,LC_PC1);
	fclose(fp);
}
void outRes(int satNumWhole,FILE *fLC,double *OC_LC0,double *OC_LC1,EKFSTATE * KFState)
{
	int i;
	for (i=0;i<satNumWhole;i++)
	{
		fprintf(fLC,"%12.3f",OC_LC0[i]);
	}
	fprintf(fLC,"%12.3f\n",KFState->ApriSigma[1]);
	for (i=0;i<satNumWhole;i++)
	{
		fprintf(fLC,"%12.3f",OC_LC1[i]);
	}
	fprintf(fLC,"%12.3f\n",KFState->PostSigma[1]);
	fclose(fLC);
}
void outResEKFPC(int satNumWhole,double *OC_PC0,double *OC_PC1,EKFSTATE * KFState)
{
	int i;
	double apriRes=0.0,posRes=0.0;
	FILE *apriRes_file,*posRes_file;
	apriRes_file=fopen("OutFile\\apriRes_EKF_PC.txt","a+");
	posRes_file=fopen("OutFile\\posRes_EKF_PC.txt","a+");
	for (i=0;i<satNumWhole;i++)
	{
		apriRes+=OC_PC0[i];
		posRes+=OC_PC1[i];
	}
	fprintf(apriRes_file,"%12.3f%12.3f\n",KFState->Time.SecOfWeek,apriRes/satNumWhole);
	fprintf(posRes_file,"%12.3f%12.3f\n",KFState->Time.SecOfWeek,posRes/satNumWhole);
	fclose(apriRes_file);
	fclose(posRes_file);
}
void outResEKFLC(int satNumWhole,double *OC_LC0,double *OC_LC1,EKFSTATE * KFState)
{
	int i;
	double apriRes=0.0,posRes=0.0;
	FILE *apriRes_file,*posRes_file;
	apriRes_file=fopen("OutFile\\apriRes_EKF_LC.txt","a+");
	posRes_file=fopen("OutFile\\posRes_EKF_LC.txt","a+");
	for (i=0;i<satNumWhole;i++)
	{
		apriRes+=OC_LC0[i];
		posRes+=OC_LC1[i];
	}
	fprintf(apriRes_file,"%12.3f%12.3f\n",KFState->Time.SecOfWeek,apriRes/satNumWhole);
	fprintf(posRes_file,"%12.3f%12.3f\n",KFState->Time.SecOfWeek,posRes/satNumWhole);
	fclose(apriRes_file);
	fclose(posRes_file);
}