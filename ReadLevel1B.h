#pragma once
#include "stdlib.h"
#include <cstring>
#include <math.h>
#include <stdio.h>
#include "satods.h"
//GNV轨道数据
struct ORBIT
{
	GPSTIME gpsTime;
	double sec;
	double x;
	double y;
	double z;
	double dx;
	double dy;
	double dz;
};
//GNV相对轨道数据
struct D_ORBIT
{
	int week;
	double sec;
	int comSatNum;
	double r;//坐标rtn
	double t;
	double n;

	double dr;//速度rtn
	double dt;
	double dn;

	double dx;//基线向量
	double dy;
	double dz;
	double ddx;//基线速度
	double ddy;
	double ddz;
	double R;//基线长
};
//相对定位数据
struct LCPC_ORBIT
{
	int week;
	double sec;//从J2000开始的秒数
	double x;//参考位置和速度，为了转换为rtn
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double dx;//基线向量残差，与GNV相对轨道进行比较得到
	double dy;
	double dz;
	int	satnum;//相对定位时的单差卫星数
	double sigPosPc;
	double sigPosLc;
	double R;//根据相对定位的数据算得的基线长
	double R0;//根据单点定轨的数据算得的基线长
	double DOP;
};
struct KBR
{
	GPSTIME gpsTime;
	double sec;
	double cor_br;
	double br;
	double ltc_br;
	double gc_br;
	int flag;//初始设置为0，若需要重新设置常值偏差则对应的时刻设置为1
};
//残差
struct KBR_R
{
	double bias;//常值偏差
	int num;
	KBR* kbr;

};
extern struct ORBIT* orbitA;
extern struct ORBIT* orbitB;
extern struct D_ORBIT* d_orbit;
extern struct KBR* kbr;
extern struct LCPC_ORBIT* Lc_orbit;
extern struct D_ORBIT* dOrbitA;
extern struct D_ORBIT* dOrbitB;
extern struct D_ORBIT* dOrbitAB;
void Readorbit(char* filename,ORBIT* orbit,int* num);
void Readkbr(char* filename,KBR* kbr,int* num);
void InitOrbit(ORBIT* orbit);
void InitD_Orbit(D_ORBIT* orbit);
void InitKBR(KBR* kbr);
void InitLc_orbit(LCPC_ORBIT* orbit);
void Out_Residual(LCPC_ORBIT* orbit,int numofepoch,int n_kbr);
int DataCheck(int doy,int year); //返回KBR文件观测值数目
void assess(int * const numOfEpoch,int * const numOfA,int * const numOfB,EKFSTATE * KFState,FINALORBIT FinalOrb[2]);