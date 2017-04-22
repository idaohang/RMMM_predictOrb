#include "ReadLevel1B.h"
#define ARCNUM  50
struct ORBIT* orbitA;
struct ORBIT* orbitB;
struct D_ORBIT* d_orbit;
struct KBR* kbr;
struct LCPC_ORBIT* Lc_orbit;
struct D_ORBIT* dOrbitA;
struct D_ORBIT* dOrbitB;
struct D_ORBIT* dOrbitAB;
void Readorbit(char* filename,ORBIT* orbit,int* num)
{
	FILE *obt;
	char line[512];
	COMMONTIME CurrCT={0};
	int m=0;
	if( ( obt = fopen(filename,"r") ) == NULL )
	{
		printf( "观测文件 %s 不能打开!\n", filename );
		return;
	}
	do
	{
		memset(line   ,0x00,sizeof(line)   );
		fgets(line,sizeof(line),obt);
	}while(!strstr(line,"END OF HEADER"));

	while( !feof( obt ) )
	{
		if( fgets( line, 512, obt ) != NULL )
		{
			sscanf( line, "%lf %*s %*s %lf %lf %lf %*lf %*lf %*lf %lf %lf %lf", 
				&orbit[m].sec,&orbit[m].x, &orbit[m].y, &orbit[m].z,&orbit[m].dx,&orbit[m].dy,&orbit[m].dz );
			SecTimeToCT(orbit[m].sec,&CurrCT);
			CommonTimeToGPSTime ( &CurrCT, &orbit[m].gpsTime );
			m++;
		}
	}
	(*num)=m;
	fclose(obt);
}
void Readkbr(char* filename,KBR* kbr,int* num)
{
	FILE *obt;
	char line[512];
	COMMONTIME CurrCT={0};
	int m=0;
	if( ( obt = fopen(filename,"r") ) == NULL )
	{
		printf( "观测文件 %s 不能打开!\n", filename );
		return;
	}
	do
	{
		memset(line   ,0x00,sizeof(line)   );
		fgets(line,sizeof(line),obt);
	}while(!strstr(line,"END OF HEADER"));
	char QF[10];

	while( !feof( obt ) )
	{
		if( fgets( line, 512, obt ) != NULL )
		{
			sscanf( &line[strlen(line)-10], "%s", QF);
			if (!strstr(QF,"1"))
			{
				sscanf( line, "%lf %lf %*lf %*lf %*lf %lf %*lf %*lf %lf", 
					&kbr[m].sec,&kbr[m].br, &kbr[m].ltc_br, &kbr[m].gc_br);//其中%*lf是跳过这个数据%*表示跳过
				kbr[m].cor_br=kbr[m].br+kbr[m].ltc_br+kbr[m].gc_br;
				SecTimeToCT(kbr[m].sec,&CurrCT);
				CommonTimeToGPSTime ( &CurrCT, &kbr[m].gpsTime );
				m++;
			}

		}
	}
	(*num)=m;
	fclose(obt);

}
void InitOrbit(ORBIT* orbit)
{
	orbit->gpsTime.SecOfWeek=0.0;
	orbit->gpsTime.Week=0;
	orbit->sec=0;
	orbit->x=0.;
	orbit->y=0.;
	orbit->z=0.;
	orbit->dx=0.;
	orbit->dy=0.;
	orbit->dz=0.;
}
void InitD_Orbit(D_ORBIT* orbit)
{
	orbit->week=0;
	orbit->sec=0;
	orbit->dx=0.;
	orbit->dy=0.;
	orbit->dz=0.;
	orbit->ddx=0.;
	orbit->ddy=0.;
	orbit->ddz=0.;
	orbit->R=0.;

}
void InitLc_orbit(LCPC_ORBIT* orbit)
{
	orbit->week=0;
	orbit->sec=0;
	orbit->dx=0.;
	orbit->dy=0.;
	orbit->dz=0.;
	orbit->R=0.;
	orbit->R0=0.;
	orbit->satnum=0;
	orbit->DOP=999;
}
void InitKBR(KBR* kbr)
{
	kbr->gpsTime.SecOfWeek=0.0;
	kbr->gpsTime.Week=0;
	kbr->sec=0;
	kbr->br=0.;
	kbr->ltc_br=0.;
	kbr->gc_br=0.;
	kbr->flag=0;
	kbr->cor_br=0.;
}
int DataCheck(int DOY,int year)
{
	if( year < 100 )     /* 使用两位年记法 */
	{
		if( year < 80 )   /* 只支持1980~2079  */
		{
			year = year + 2000;  
		}
		else   
		{
			year = year + 1900;
		}
	}
	int month,day;
	int monthDay[12]={31,28,31,30,31,30,31,31,30,31,30,29};//只需要前十一个月的即可，最后一个表示闰年的2月
	if(DOY>366)
		printf("日输入错误!");
	day=DOY;
	for(month=0;month<12;month++)
	{
		if(day<monthDay[month])
		{
			month=month+1;
			break;
		}
		day=day-monthDay[month];
	}
	//初始化
	for (int index=0;index<17330;index++)
		InitOrbit(orbitB+index);
	for (int index=0;index<17330;index++)
		InitOrbit(orbitA+index);
	for (int index=0;index<17330;index++)
		InitD_Orbit(d_orbit+index);
	for (int index=0;index<17330;index++)
		InitKBR(kbr+index);
	char GRACE_B[100],GRACE_A[100],kbrfile[100],ResFileName[100],ResFile[100];
	sprintf( GRACE_B, "InputFile\\GNV1B_%04d-%02d-%02d_B_02.asc",year, month,day);
	sprintf( GRACE_A, "InputFile\\GNV1B_%04d-%02d-%02d_A_02.asc",year, month,day);
	sprintf( kbrfile, "InputFile\\KBR1B_%04d-%02d-%02d_X_02.asc",year, month,day);
	
	int n_A,n_B,n_kbr;
	Readorbit(GRACE_B,orbitB,&n_B);
	Readorbit(GRACE_A,orbitA,&n_A);
	Readkbr(kbrfile,kbr,&n_kbr);

	double time,time1;
	int temp_num=0;
	int i,j;
	double dx,dy,dz,ddx,ddy,ddz,r;
	sprintf( ResFileName, "InputFile\\A_B.txt");
	FILE* Fout  = fopen( ResFileName, "a+" );
	for (i=0;i<n_B;i++)
	{
		time=orbitB[i].sec;
		for (j=0;j<n_A;j++)
		{
			time1=orbitA[j].sec;
			if (time==time1)
			{
				dx =orbitA[j].x -orbitB[i].x;
				dy =orbitA[j].y -orbitB[i].y;
				dz =orbitA[j].z -orbitB[i].z;
				ddx=orbitA[j].dx-orbitB[i].dx;
				ddy=orbitA[j].dy-orbitB[i].dy;
				ddz=orbitA[j].dz-orbitB[i].dz;
				r=sqrt(dx*dx+dy*dy+dz*dz);
				d_orbit[temp_num].sec=time;
				d_orbit[temp_num].dx=dx;
				d_orbit[temp_num].dy=dy;
				d_orbit[temp_num].dz=dz;
				d_orbit[temp_num].ddx=ddx;
				d_orbit[temp_num].ddy=ddy;
				d_orbit[temp_num].ddz=ddz;
				d_orbit[temp_num].R=r;

				fprintf( Fout, " %12.3f %12.3f %12.3f %12.3f %12.3f %10.3f %10.3f %10.3f\r\n",
					time,dx,dy,dz,r,ddx,ddy,ddz);
				temp_num++;
				break;
			}
		}
	}
	fclose(Fout);
	
	sprintf( ResFile, "InputFile\\kbr_R-%04d-%02d-%02d.txt",year, month,day);
	FILE* Ft  = fopen( ResFile, "w+" );
	double bias0,bias1;
	int FLAG=0;
	KBR_R* kbr_r;
	kbr_r=(KBR_R*)malloc(ARCNUM*sizeof(struct KBR_R));
	for (i=0;i<ARCNUM;i++)
	{
		kbr_r[i].kbr=(KBR*)malloc(17330*sizeof(struct KBR));
	}
	//当前弧段数
	int arc=-1;
	//当前弧段内观测值个数
	int num=0;
	double square[ARCNUM]={0.};//对应弧段的平方和
	double std[ARCNUM]={0.};//对应弧段的标准差
	for (i=0;i<temp_num;i++)
	{
		time=d_orbit[i].sec;
		for (j=0;j<n_kbr;j++)
		{
			time1=kbr[j].sec;
			if (time==time1)
			{
				if (FLAG==0)
				{
					bias0=d_orbit[i].R-kbr[j].cor_br;
					FLAG=1;
					arc++;
					square[arc]+=bias0*bias0;
					memcpy(&(kbr_r[arc].kbr[num]),&kbr[j],sizeof(struct KBR));
					num++;
				}
				else
				{
					bias1=d_orbit[i].R-kbr[j].cor_br;
					if ((bias1<bias0+2)&&(bias1>bias0-2))//说明为同一常值偏差
					{
						square[arc]+=bias1*bias1;
						memcpy(&(kbr_r[arc].kbr[num]),&kbr[j],sizeof(struct KBR));
						num++;
					}
					else
					{
						if (bias0>0)
							std[arc]=sqrt(square[arc]/num);
						else
							std[arc]=(-1)*sqrt(square[arc]/num);
						kbr_r[arc].bias=std[arc];
						kbr_r[arc].num=num;
						arc++;
						num=0;
						square[arc]+=bias1*bias1;
						memcpy(&(kbr_r[arc].kbr[num]),&kbr[j],sizeof(struct KBR));
						num++;
						bias0=bias1;
					}
				}
				break;
			}
		}
	}
	if (bias0>0)
		std[arc]=sqrt(square[arc]/num);
	else
		std[arc]=(-1)*sqrt(square[arc]/num);
	kbr_r[arc].bias=std[arc];
	kbr_r[arc].num=num;
	double res;
	int k;
	//数据输出
	/*char ResFile[100]="GNVresfile.txt";
	char Restime[100]="GNVrestime.txt";
	char Resbias[100]="GNVresbias.txt";

	FILE* Fres  = fopen( ResFile, "w+" );
	FILE* Ftime  = fopen( Restime, "w+" );
	FILE* Fbias  = fopen( Resbias, "w+" );*/

	for (i=0;i<temp_num;i++)
	{
		time=d_orbit[i].sec;
		for (k=0;k<=arc;k++)
		{
			for (j=0;j<kbr_r[k].num;j++)
			{
				time1=kbr_r[k].kbr[j].sec;
				if (time==time1)
				{
					res=d_orbit[i].R-kbr_r[k].kbr[j].cor_br-kbr_r[k].bias;
					/*fprintf( Ftime, " %d\r\n",time);
					fprintf( Fres, " %8.4f\r\n",res);*/
					fprintf( Ft, " %12.3f %8.4f %12.3f\r\n",time,res,kbr_r[k].bias);
				}
			}
		}
	}
	/*for (k=0;k<=arc;k++)
	{
		fprintf( Fbias, " %12.3f %d\r\n",kbr_r[k].bias,kbr_r[k].num);
	}*/

	/*fclose(Fres);
	fclose(Ftime);
	fclose(Fbias);*/
	fclose(Ft);

	free(kbr_r);
	return n_kbr;

}
void Out_Residual(LCPC_ORBIT* orbit,int numofepoch,int n_kbr)
{
	double time,time1;
	int i,j;
	
	double bias0,bias1;
	int FLAG=0;
	KBR_R* kbr_r2;
	kbr_r2=(KBR_R*)malloc(ARCNUM*sizeof(struct KBR_R));
	for (i=0;i<ARCNUM;i++)
	{
		kbr_r2[i].kbr=(KBR*)malloc(17330*sizeof(struct KBR));
	}
	//当前弧段数
	int arc=-1;
	//当前弧段内观测值个数
	int num=0;
	double square[ARCNUM]={0.};//对应弧段的平方和
	double std[ARCNUM]={0.};//对应弧段的标准差
	for (i=0;i<numofepoch;i++)
	{
		time=orbit[i].sec;
		for (j=0;j<n_kbr;j++)
		{
			time1=kbr[j].sec;
			if (time==time1)
			{
				if (FLAG==0)
				{
					bias0=orbit[i].R-kbr[j].cor_br;
					FLAG=1;
					arc++;
					square[arc]+=bias0*bias0;
					memcpy(&(kbr_r2[arc].kbr[num]),&kbr[j],sizeof(struct KBR));
					num++;
				}
				else
				{
					bias1=orbit[i].R-kbr[j].cor_br;
					if ((bias1<bias0+15)&&(bias1>bias0-15))//说明为同一常值偏差
					{
						square[arc]+=bias1*bias1;
						memcpy(&(kbr_r2[arc].kbr[num]),&kbr[j],sizeof(struct KBR));
						num++;
					}
					else
					{
						if (bias0>0)
							std[arc]=sqrt(square[arc]/num);
						else
							std[arc]=(-1)*sqrt(square[arc]/num);
						kbr_r2[arc].bias=std[arc];
						kbr_r2[arc].num=num;
						arc++;
						num=0;
						square[arc]+=bias1*bias1;
						memcpy(&(kbr_r2[arc].kbr[num]),&kbr[j],sizeof(struct KBR));
						num++;
						bias0=bias1;
					}
				}
				break;
			}
		}
	}
	if (bias0>0)
		std[arc]=sqrt(square[arc]/num);
	else
		std[arc]=(-1)*sqrt(square[arc]/num);
	kbr_r2[arc].bias=std[arc];
	kbr_r2[arc].num=num;
	double res;
	double rms,rms1,rms05,rms02,rms01;
	int numuse=0;
	int numuse1=0;
	int numuse05=0;
	int numuse02=0;
	int numuse01=0;
	int numusedop=0;
	rms=0.;rms1=0.;rms05=0.;rms02=0.;rms01=0.;
	int k;
	//残差数据输出
	char ResFile[100]="resfile.txt";
	char Restime[100]="restime.txt";
	char Resbias[100]="resbias.txt";
	char satnum[100]="satnum.txt";
	char sstd[100]="std0602.txt";
	FILE* Fres  = fopen( ResFile, "w+" );
	FILE* Ftime  = fopen( Restime, "w+" );
	FILE* Fbias  = fopen( Resbias, "w+" );
	FILE* Fsatnum  = fopen( satnum, "w+" );
	FILE* Fstd  = fopen( sstd, "w+" );

	char std_1[100]="std_1.txt";
	char std_05[100]="std_05.txt";
	char std_02[100]="std_02.txt";
	char std_01[100]="std_01.txt";
	FILE* Fstd_1=fopen(std_1,"a+");
	FILE* Fstd_05=fopen(std_05,"a+");
	FILE* Fstd_02=fopen(std_02,"a+");
	FILE* Fstd_01=fopen(std_01,"a+");

	char DOP[100]="DOP.txt";
	FILE* FDOP=fopen(DOP,"w+");
	char Eerror[100]="error.txt";
	FILE* Feror=fopen(Eerror,"w+");
	for (i=0;i<numofepoch;i++)
	{
		time=orbit[i].sec;
		for (k=0;k<=arc;k++)
		{
			for (j=0;j<kbr_r2[k].num;j++)
			{
				time1=kbr_r2[k].kbr[j].sec;
				if (time==time1)
				{
					res=orbit[i].R-kbr_r2[k].kbr[j].cor_br-kbr_r2[k].bias;
					numuse++;
					if (orbit[i].DOP>1)
					{
						fprintf(Feror," %8.4f %8.4f\r\n",res,orbit[i].DOP);
						numusedop++;
					}
					else
					{
						
						fprintf( Ftime, " %d\r\n",time);
						fprintf( Fres, " %8.4f\r\n",res);

						rms+=res*res;
						
						if (fabs(res)<=1)
						{
							rms1+=res*res;
							numuse1++;
						}
						if (fabs(res)<=0.5)
						{
							rms05+=res*res;
							numuse05++;
						}
						if (fabs(res)<=0.2)
						{
							rms02+=res*res;
							numuse02++;
						}
						if (fabs(res)<=0.1)
						{
							rms01+=res*res;
							numuse01++;
						}

					}	
					break;
				}
			}
		}
	}
	for (k=0;k<=arc;k++)
	{
		fprintf( Fbias, " %12.3f %d\r\n",kbr_r2[k].bias,kbr_r2[k].num);
	}
	rms=sqrt(rms/numuse)*100;
	rms1=sqrt(rms1/numuse1)*100;
	rms05=sqrt(rms05/numuse05)*100;
	rms02=sqrt(rms02/numuse02)*100;
	rms01=sqrt(rms01/numuse01)*100;
	fprintf( Fstd, " %12.3f %10d %12.3f\r\n",rms,numuse,1.0);
	fprintf( Fstd, " %12.3f %10d %12.3f\r\n",rms1,numuse1,numuse1/(double)numuse);
	fprintf( Fstd, " %12.3f %10d %12.3f\r\n",rms05,numuse05,numuse05/(double)numuse);
	fprintf( Fstd, " %12.3f %10d %12.3f\r\n",rms02,numuse02,numuse02/(double)numuse);
	fprintf( Fstd, " %12.3f %10d %12.3f\r\n",rms01,numuse01,numuse01/(double)numuse);
	fprintf( Fstd, " 粗差率： %12.3f\r\n",numusedop/(double)numuse);

	fprintf(Fstd_1," %12.3f\r\n",numuse1/(double)numuse);
	fprintf(Fstd_05," %12.3f\r\n",numuse05/(double)numuse);
	fprintf(Fstd_02," %12.3f\r\n",numuse02/(double)numuse);
	fprintf(Fstd_01," %12.3f\r\n",numuse01/(double)numuse);

	fclose(Fres);
	fclose(FDOP);
	fclose(Ftime);
	fclose(Fbias);
	fclose(Fsatnum);
	fclose(Fstd);

	fclose(Feror);

	fclose(Fstd_1);
	fclose(Fstd_05);
	fclose(Fstd_02);
	fclose(Fstd_01);
	free(kbr_r2);
}
/********************
统计结果模块
********************/
void assess(int * const numOfEpoch,int * const numOfA,int * const numOfB,EKFSTATE * KFState,FINALORBIT FinalOrb[2])
{
	int i;
	double stateB[6]={0.0};
	double transXYZV[6]	={0.0};
	double dXYZ[3]		={0.0};
	double dRTN[3]		={0.0};
	/*********************************************************
	********************将滤波值转换到ITRF质心****************
	**********************************************************/
	Lc_orbit[*numOfEpoch].week=KFState->Time.Week;
	Lc_orbit[*numOfEpoch].sec=(int)KFState->Time.SecOfWeek;
	CopyArray( 6, FinalOrb[0].ECI_Orb, KFState->StateA );
	getBpos(KFState->StateA,KFState->StateRel,stateB);
	CopyArray( 6, FinalOrb[1].ECI_Orb, stateB );
	PhaseCentToMassCent( true,&KFState->CentBias[0], FinalOrb[0].ECI_Orb );	
	PhaseCentToMassCent( true,&KFState->CentBias[3], FinalOrb[1].ECI_Orb );
	ICRF_ITRF_GPST( MJD_J2000,  &KFState->Time, true, FinalOrb[0].ECI_Orb, FinalOrb[0].ECF_Orb );//转到地固系
	ICRF_ITRF_GPST( MJD_J2000,  &KFState->Time, true, FinalOrb[1].ECI_Orb, FinalOrb[1].ECF_Orb );//转到地固系
	Lc_orbit[*numOfEpoch].satnum=KFState->comSatNumUsed;
	Lc_orbit[*numOfEpoch].dx=FinalOrb[0].ECF_Orb[0]-FinalOrb[1].ECF_Orb[0];//程序估计的位置差
	Lc_orbit[*numOfEpoch].dy=FinalOrb[0].ECF_Orb[1]-FinalOrb[1].ECF_Orb[1];
	Lc_orbit[*numOfEpoch].dz=FinalOrb[0].ECF_Orb[2]-FinalOrb[1].ECF_Orb[2];
	Lc_orbit[*numOfEpoch].vx=FinalOrb[0].ECF_Orb[3]-FinalOrb[1].ECF_Orb[3];//程序估计的速度差
	Lc_orbit[*numOfEpoch].vy=FinalOrb[0].ECF_Orb[4]-FinalOrb[1].ECF_Orb[4];
	Lc_orbit[*numOfEpoch].vz=FinalOrb[0].ECF_Orb[5]-FinalOrb[1].ECF_Orb[5];
	Lc_orbit[*numOfEpoch].R=sqrt(Lc_orbit[*numOfEpoch].dx*Lc_orbit[*numOfEpoch].dx+Lc_orbit[*numOfEpoch].dy*Lc_orbit[*numOfEpoch].dy+
				Lc_orbit[*numOfEpoch].dz*Lc_orbit[*numOfEpoch].dz);
	Lc_orbit[*numOfEpoch].sigPosPc=KFState->sigPosPC;	
	Lc_orbit[*numOfEpoch].sigPosLc=KFState->sigPosLC;


	/*********************************************************
	************************A星位置评估**********************
	**********************************************************/
	for(i=*numOfA;i<17330;i++)
		if(fabs((KFState->Time.Week-orbitA[i].gpsTime.Week)*604800.0+KFState->Time.SecOfWeek-orbitA[i].gpsTime.SecOfWeek)<1E-3)
			{
				*numOfA=i;
				break;
			}
	dOrbitA[*numOfEpoch].week=KFState->Time.Week;
	dOrbitA[*numOfEpoch].sec=(int)KFState->Time.SecOfWeek;

	dOrbitA[*numOfEpoch].dx= orbitA[*numOfA].x- FinalOrb[0].ECF_Orb[0];//A星估计的与真值的差异
	dOrbitA[*numOfEpoch].dy= orbitA[*numOfA].y- FinalOrb[0].ECF_Orb[1];
	dOrbitA[*numOfEpoch].dz= orbitA[*numOfA].z- FinalOrb[0].ECF_Orb[2];
	dOrbitA[*numOfEpoch].ddx=orbitA[*numOfA].dx-FinalOrb[0].ECF_Orb[3];
	dOrbitA[*numOfEpoch].ddy=orbitA[*numOfA].dy-FinalOrb[0].ECF_Orb[4];
	dOrbitA[*numOfEpoch].ddz=orbitA[*numOfA].dz-FinalOrb[0].ECF_Orb[5];
	dOrbitA[*numOfEpoch].R=sqrt(dOrbitA[*numOfEpoch].dx*dOrbitA[*numOfEpoch].dx+dOrbitA[*numOfEpoch].dy*dOrbitA[*numOfEpoch].dy+
		dOrbitA[*numOfEpoch].dz*dOrbitA[*numOfEpoch].dz);
				//A星转RTN
	transXYZV[0]=orbitA[*numOfA].x;
	transXYZV[1]=orbitA[*numOfA].y;
	transXYZV[2]=orbitA[*numOfA].z;
	transXYZV[3]=orbitA[*numOfA].dx;
	transXYZV[4]=orbitA[*numOfA].dy;
	transXYZV[5]=orbitA[*numOfA].dz;
	dXYZ[0]		=dOrbitA[*numOfEpoch].dx;
	dXYZ[1]		=dOrbitA[*numOfEpoch].dy;
	dXYZ[2]		=dOrbitA[*numOfEpoch].dz;
	XYZToRTN( transXYZV, dXYZ, dRTN );
	dOrbitA[*numOfEpoch].r=dRTN[0];
	dOrbitA[*numOfEpoch].t=dRTN[1];
	dOrbitA[*numOfEpoch].n=dRTN[2];
	
	dXYZ[0]		=dOrbitA[*numOfEpoch].ddx;
	dXYZ[1]		=dOrbitA[*numOfEpoch].ddy;
	dXYZ[2]		=dOrbitA[*numOfEpoch].ddz;
	XYZToRTN( transXYZV, dXYZ, dRTN );
	dOrbitA[*numOfEpoch].dr=dRTN[0];
	dOrbitA[*numOfEpoch].dt=dRTN[1];
	dOrbitA[*numOfEpoch].dn=dRTN[2];
	/*********************************************************
	************************B星位置评估**********************
	**********************************************************/
	for(i=*numOfB;i<17330;i++)
		if(fabs((KFState->Time.Week-orbitB[i].gpsTime.Week)*604800.0+KFState->Time.SecOfWeek-orbitB[i].gpsTime.SecOfWeek)<1E-3)
		{
			*numOfB=i;
			break;
		}
	dOrbitB[*numOfEpoch].week=KFState->Time.Week;
	dOrbitB[*numOfEpoch].sec=(int)KFState->Time.SecOfWeek;
	dOrbitB[*numOfEpoch].dx=orbitB[*numOfB].x- FinalOrb[1].ECF_Orb[0];//B星估计的与真值的差异
	dOrbitB[*numOfEpoch].dy=orbitB[*numOfB].y- FinalOrb[1].ECF_Orb[1];
	dOrbitB[*numOfEpoch].dz=orbitB[*numOfB].z- FinalOrb[1].ECF_Orb[2];
	dOrbitB[*numOfEpoch].ddx=orbitB[*numOfB].dx-FinalOrb[1].ECF_Orb[3];
	dOrbitB[*numOfEpoch].ddy=orbitB[*numOfB].dy-FinalOrb[1].ECF_Orb[4];
	dOrbitB[*numOfEpoch].ddz=orbitB[*numOfB].dz-FinalOrb[1].ECF_Orb[5];
	dOrbitB[*numOfEpoch].R=sqrt(dOrbitB[*numOfEpoch].dx*dOrbitB[*numOfEpoch].dx+dOrbitB[*numOfEpoch].dy*dOrbitB[*numOfEpoch].dy+
		dOrbitB[*numOfEpoch].dz*dOrbitB[*numOfEpoch].dz);
				//B星转RTN
	transXYZV[0]=orbitB[*numOfB].x;
	transXYZV[1]=orbitB[*numOfB].y;
	transXYZV[2]=orbitB[*numOfB].z;
	transXYZV[3]=orbitB[*numOfB].dx;
	transXYZV[4]=orbitB[*numOfB].dy;
	transXYZV[5]=orbitB[*numOfB].dz;

	dXYZ[0]		=dOrbitB[*numOfEpoch].dx;
	dXYZ[1]		=dOrbitB[*numOfEpoch].dy;
	dXYZ[2]		=dOrbitB[*numOfEpoch].dz;
	XYZToRTN( transXYZV, dXYZ, dRTN );
	dOrbitB[*numOfEpoch].r=dRTN[0];
	dOrbitB[*numOfEpoch].t=dRTN[1];
	dOrbitB[*numOfEpoch].n=dRTN[2];
	//速度
	dXYZ[0]		=dOrbitB[*numOfEpoch].ddx;
	dXYZ[1]		=dOrbitB[*numOfEpoch].ddy;
	dXYZ[2]		=dOrbitB[*numOfEpoch].ddz;
	XYZToRTN( transXYZV, dXYZ, dRTN );
	dOrbitB[*numOfEpoch].dr=dRTN[0];
	dOrbitB[*numOfEpoch].dt=dRTN[1];
	dOrbitB[*numOfEpoch].dn=dRTN[2];
	/*********************************************************
	************************相对位置评估**********************
	**********************************************************/
	dOrbitAB[*numOfEpoch].week=KFState->Time.Week;
	dOrbitAB[*numOfEpoch].sec=(int)KFState->Time.SecOfWeek;
	dOrbitAB[*numOfEpoch].dx=orbitA[*numOfA].x-orbitB[*numOfB].x;//真值的AB坐标差
	dOrbitAB[*numOfEpoch].dy=orbitA[*numOfA].y-orbitB[*numOfB].y;
	dOrbitAB[*numOfEpoch].dz=orbitA[*numOfA].z-orbitB[*numOfB].z;
	dOrbitAB[*numOfEpoch].ddx=orbitA[*numOfA].dx-orbitB[*numOfB].dx;
	dOrbitAB[*numOfEpoch].ddy=orbitA[*numOfA].dy-orbitB[*numOfB].dy;
	dOrbitAB[*numOfEpoch].ddz=orbitA[*numOfA].dz-orbitB[*numOfB].dz;
	dOrbitAB[*numOfEpoch].R=sqrt(dOrbitAB[*numOfEpoch].dx*dOrbitAB[*numOfEpoch].dx+
		dOrbitAB[*numOfEpoch].dy*dOrbitAB[*numOfEpoch].dy+
		dOrbitAB[*numOfEpoch].dz*dOrbitAB[*numOfEpoch].dz);
	//AB坐标差转为RTN，因为前面B星的参考为B真值，故这里不用再赋值参考点
	//若把RTN转到相对RTN上，则需要重新赋值参考点
	transXYZV[0]=dOrbitAB[*numOfEpoch].dx;
	transXYZV[1]=dOrbitAB[*numOfEpoch].dy;
	transXYZV[2]=dOrbitAB[*numOfEpoch].dz;
	transXYZV[3]=dOrbitAB[*numOfEpoch].ddx;
	transXYZV[4]=dOrbitAB[*numOfEpoch].ddy;
	transXYZV[5]=dOrbitAB[*numOfEpoch].ddz;
	dXYZ[0]		=Lc_orbit[*numOfEpoch].dx-dOrbitAB[*numOfEpoch].dx;
	dXYZ[1]		=Lc_orbit[*numOfEpoch].dy-dOrbitAB[*numOfEpoch].dy;
	dXYZ[2]		=Lc_orbit[*numOfEpoch].dz-dOrbitAB[*numOfEpoch].dz; 
	XYZToRTN( transXYZV, dXYZ, dRTN );
	dOrbitAB[*numOfEpoch].r=dRTN[0];
	dOrbitAB[*numOfEpoch].t=dRTN[1];
	dOrbitAB[*numOfEpoch].n=dRTN[2];

	dXYZ[0]		=Lc_orbit[*numOfEpoch].vx-dOrbitAB[*numOfEpoch].ddx;
	dXYZ[1]		=Lc_orbit[*numOfEpoch].vy-dOrbitAB[*numOfEpoch].ddy;
	dXYZ[2]		=Lc_orbit[*numOfEpoch].vz-dOrbitAB[*numOfEpoch].ddz; 
	XYZToRTN( transXYZV, dXYZ, dRTN );
	dOrbitAB[*numOfEpoch].dr=dRTN[0];
	dOrbitAB[*numOfEpoch].dt=dRTN[1];
	dOrbitAB[*numOfEpoch].dn=dRTN[2];
	*numOfEpoch=(*numOfEpoch)++;
}