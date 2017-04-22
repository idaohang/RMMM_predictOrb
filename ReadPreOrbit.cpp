#pragma once
#include "ReadPreOrbit.h"

void ReadSp3File(char* filename,ONEPOCHOFSP3 *OBSP)
{
	FILE *fpSp3;	
	if( ( fpSp3 = fopen(filename,"r") ) == NULL )
	{
		printf( "观测文件 %s 不能打开!\n", filename );
		return;
	}

	int					satnum;
	char				LineOfStr[256];
	COMMONTIME			sTtemp;
	GPSTIME				 pweek;
	double				 m_doy;
    //获取精密星历中有效卫星的数目
	for (int i=0;i<3;i++)
	{
		memset(LineOfStr   ,0x00,sizeof(LineOfStr)   );
		fgets(LineOfStr,sizeof(LineOfStr),fpSp3);
	}
    char tempstr[256] = {0x00};
	memset( tempstr, 0x00, sizeof(tempstr) );
	satnum = atoi( strncpy( tempstr, LineOfStr+4, 2 ) );
    //读完前22行
	for (int i=0;i<19;i++)
	{
		memset(LineOfStr   ,0x00,sizeof(LineOfStr)   );
		fgets(LineOfStr,sizeof(LineOfStr),fpSp3);
	}
	//开始读取有效数据信息
	static int sp3_num=0;
	memset(LineOfStr ,0x00,sizeof(LineOfStr));
	fgets(LineOfStr,sizeof(LineOfStr),fpSp3);
	while(!strstr(&LineOfStr[0],"EOF"))
	{
		//时间读取
		OBSP[sp3_num].Satnum      =satnum;
		memset( tempstr, 0x00, sizeof(tempstr) );
		OBSP[sp3_num].pT.Year   = atoi( strncpy( tempstr, LineOfStr+3, 4 ) );
		sTtemp.Year   = atoi( strncpy( tempstr, LineOfStr+3, 4 ) );
		memset( tempstr, 0x00, sizeof(tempstr) );
		OBSP[sp3_num].pT.Month  = atoi( strncpy( tempstr, LineOfStr+8 , 2 ) );
		sTtemp.Month  = atoi( strncpy( tempstr, LineOfStr+8 , 2 ) );
		memset( tempstr, 0x00, sizeof(tempstr) );
		OBSP[sp3_num].pT.Day    = atoi( strncpy( tempstr, LineOfStr+11, 2 ) );
		sTtemp.Day    = atoi( strncpy( tempstr, LineOfStr+11, 2 ) );
		memset( tempstr, 0x00, sizeof(tempstr) );
		OBSP[sp3_num].pT.Hour   = atoi( strncpy( tempstr, LineOfStr+14, 2 ) );
		sTtemp.Hour   = atoi( strncpy( tempstr, LineOfStr+14, 2 ) );
		memset( tempstr, 0x00, sizeof(tempstr) );
		OBSP[sp3_num].pT.Minute    = atoi( strncpy( tempstr, LineOfStr+17, 2 ) );
		sTtemp.Minute    = atoi( strncpy( tempstr, LineOfStr+17, 2 ) );
		memset( tempstr, 0x00, sizeof(tempstr) );
		OBSP[sp3_num].pT.Second = atof( strncpy( tempstr, LineOfStr+20, 2 ) );
		sTtemp.Second = atof( strncpy( tempstr, LineOfStr+20, 10 ) );

		CommonTimeToGPSTime(&sTtemp,&pweek);
		OBSP[sp3_num].dWeeksec    =pweek.SecOfWeek;
		OBSP[sp3_num].WeekNumber  =pweek.Week;
		OBSP[sp3_num].totalWeeksec=pweek.Week*604800.0+pweek.SecOfWeek;

		//星历数据读取
		for (int j=0;j<satnum;j++)
		{
			memset(LineOfStr ,0x00,sizeof(LineOfStr));
			fgets(LineOfStr,sizeof(LineOfStr),fpSp3);

			memset( tempstr, 0x00, sizeof(tempstr) );
			int m_PRN=atoi( strncpy( tempstr, LineOfStr+2, 2 ) );
			OBSP[sp3_num].obs[m_PRN-1].PRN      =m_PRN;
			
			memset( tempstr, 0x00, sizeof(tempstr) );
			OBSP[sp3_num].obs[m_PRN-1].satposx  =atof( strncpy( tempstr, LineOfStr+4, 14 ) )*1000.0;//单位为m
			
			memset( tempstr, 0x00, sizeof(tempstr) );
			OBSP[sp3_num].obs[m_PRN-1].satposy  =atof( strncpy( tempstr, LineOfStr+18, 14 ) )*1000.0;//单位为m
			memset( tempstr, 0x00, sizeof(tempstr) );
			OBSP[sp3_num].obs[m_PRN-1].satposz  =atof( strncpy( tempstr, LineOfStr+32, 14 ) )*1000.0;//单位为m
			memset( tempstr, 0x00, sizeof(tempstr) );
			OBSP[sp3_num].obs[m_PRN-1].satclk   =atof( strncpy( tempstr, LineOfStr+46, 14 ) )*1E-6;
			OBSP[sp3_num].obs[m_PRN-1].ISExist=true;
			if(fabs(atof( strncpy( tempstr, LineOfStr+4, 14 ) ))>999999||fabs(atof( strncpy( tempstr, LineOfStr+18, 14 ) ))>999999
				||fabs(atof( strncpy( tempstr, LineOfStr+32, 14 ) ))>999999||fabs(atof( strncpy( tempstr, LineOfStr+46, 14 ) ))>999999)
				OBSP[sp3_num].obs[m_PRN-1].ISExist=false;

		}
		sp3_num++;
		memset(LineOfStr ,0x00,sizeof(LineOfStr));
		fgets(LineOfStr,sizeof(LineOfStr),fpSp3);

	}

	fclose(fpSp3);

}

int FindFitPoint_SP3(const int PRN, GPSTIME OneEpoch,int n,struct ONEPOCHOFSP3 *OBSP, double t[], 
					 double  satPos[3][30], double clk[])//2n为拟合阶数，因为最大阶数即节点数减1
{
	double        totalweeksec,temp,res;
	int m=0;
	
	totalweeksec=OneEpoch.Week*604800.0+OneEpoch.SecOfWeek;
	
	do 
	{
		temp=OBSP[m].totalWeeksec;
		res=temp-totalweeksec;
		m++;
	} while(fabs(res)>450);
	m=m-1;
	for (int i=m-n;i<=m+n;i++)
	{
		if (!OBSP[i].obs[PRN-1].ISExist)
		{
			return 0;		//拟合序列中存在bad或absent的星历
		}
		t[i-m+n]  = OBSP[i].totalWeeksec;
		satPos[0][i-m+n] = OBSP[i].obs[PRN-1].satposx;
		satPos[1][i-m+n] = OBSP[i].obs[PRN-1].satposy;
		satPos[2][i-m+n] = OBSP[i].obs[PRN-1].satposz;
		clk[i-m+n] = OBSP[i].obs[PRN-1].satclk;
	}
	return m;
}
/* polynomial interpolation by Neville's algorithm ---------------------------*/

double interpPos(const double * xi,const double * yi,int n,const double x)
{
	double add = 0.0,multi;
	for (int k=0;k<2*n+1;k++)
	{
		multi=1.0;
		for (int i=0;i<2*n+1;i++)
		{
			if (i!=k)
			{
				multi*=(x-xi[i])/(xi[k]-xi[i]);
			}

		}
		add+=multi*yi[k];	
	}
	return add;
}
void Lagrange_SP3(GPSTIME OneEpoch,int PRN,int n,
                  struct ONEPOCHOFSP3 *OBSP,ONESATEHP *PREORBIT_SP3)
{
	int i;
	double t[30];//时间
	double satPos[3][30];		//X、Y、Z坐标
	double clk[30];				//钟差
	double pos0[3];				//TOTALSEC处的位置
	double pos1[3];				//微小时间差的位置
	double totalsec=OneEpoch.Week*604800.0+OneEpoch.SecOfWeek;
	double tt=1E-3;

	int mid = FindFitPoint_SP3(PRN,OneEpoch,n,OBSP, t,satPos, clk);
	if (mid == 0)
	{
		PREORBIT_SP3->ISExist=false;
		return;
	}
	//地球自转改正
	//Rotation(2*n+1,PRN,dt,FITPOINTSP3);
	//内插坐标
	for(i=0;i<3;i++)
		pos0[i]=interpPos(t,satPos[i],n,totalsec);
	//内插坐标
	for(i=0;i<3;i++)
		pos1[i]=interpPos(t,satPos[i],n,totalsec+tt);	
	//钟差
	PREORBIT_SP3->satclk=interpPos(t,clk,n,totalsec); 
	PREORBIT_SP3->satposx=pos0[0];
	PREORBIT_SP3->satposy=pos0[1];
	PREORBIT_SP3->satposz=pos0[2];
	PREORBIT_SP3->satvelx=(pos1[0]-pos0[0])/tt;
	PREORBIT_SP3->satvely=(pos1[1]-pos0[1])/tt;
	PREORBIT_SP3->satvelz=(pos1[2]-pos0[2])/tt;

	PREORBIT_SP3->ISExist=true;
	PREORBIT_SP3->PRN=PRN;
	
}
