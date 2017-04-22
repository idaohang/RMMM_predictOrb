/****************************************************************************
目的：    定义LEOSAOD软件需要的所有结构体
编写时间：2008.21.42
版本:     V1.1
版权：    武汉大学
****************************************************************************/
#pragma once 
#ifndef _RTOD_STRUCTS_H_
#define _RTOD_STRUCTS_H_

#define ADIMENSION   15 /*position,velocity,GPS clk,GLONASS clk,clksft,Cd,Cr,tao,W*/
#define RELDIMENSION 36 /*position,velocity,GPS rel clk,Cd, Cr,W,L1模糊度12个,宽巷模糊度12个*/
#define MAXCHANNUM 12
#define MAXOBSTYPENUM 11
#define TT_TAI   32.184
#define GPST_TAI -19.0
#define MaxValiEpochNum 20
#define MAXOBSNUMDAY 8640
//   LEMTA_1  the wavelength of the primary carry frequency of GPS.
#define	LEMTA_1		0.190293672798         //单位为米
//   LEMTA_2  the wavelength of the secondary carry frequency.
#define LEMTA_2		0.2442102134           //单位为米
/* 导航卫星系统定义 */
enum GNSSSys { UNKS, GPS, GLONASS, GALILEO, GEO, COMPASS };

/* 观测数据类型定义 */
enum OBSTYPE  { UNKOBS, C1, P1, L1, D1, C2, P2, L2, D2, C5, L5, D5, S1, S2,Nw,Sw };  


typedef struct _COMMON_TIME_   /* 通用时间定义 */
{
	unsigned short Year;
	unsigned short Month;
	unsigned short Day;
	unsigned short Hour;
	unsigned short Minute;
	double Second;
} COMMONTIME;

typedef struct _GPS_TIME_              /* GPS时间定义 */
{
	unsigned short Week;          
	double         SecOfWeek;     
} GPSTIME;

typedef struct _MJD_TIME_             /* 简化儒略日 */
{
	int Days;             
	double FracDay;      
} MJDTIME;

/* 地球自转参数结构体定义 */
typedef struct _EOP_PARA_
{
	double Mjd;          /* 简化儒略日 */
	double x;            /* 极移参数x[角秒] */
	double y;            /* 极移参数y[角秒] */
	double dUT1;         /* UT1-UTC [s]     */
	double LeapSec;      /* 当前时刻的跳秒数[s] */
	char   Status;       /* 初始化标志 */
} EOPPARA;


/***************************************************************************
DYNMODELPARA

目的:  定义动力学模型使用的参数

***************************************************************************/
typedef struct _DYNAMICMODEL_
{
	int     n_a, m_a;   /* Degree and order of gravity field */
	double  Cd;         /* Coefficient of atmospheric drag */
	double  Cr;         /* Coefficient of solar radiation pressure */
	double  Mass;       /* Mass of spacecraft [kg] */
	double  Area_H;     /* surface area of spacecraft Cross-section[m^2]*/
	double  Area_R;     /* surface area of spacecraft solar radiation pressure[m^2]*/
	double  GM;         /* Gravitational coefficient [m^3/s^2] */
	double  R_ref;      /* Reference radius [m]   */
} DYNMODELPARA;

typedef struct _IONO_CORRECT_PARA_
{
	double alpha[4];        /* 电离层改正模型参数 alpha */
	double beta[4];         /* 电离层改正模型参数 beta  */
	bool   IsValid;         /* 电离层改正参数有效为true */
} IONOPARA;

typedef struct _GLONASS_TIMESYSTEM_CORRECTION_
{
	double  TauC;
	GPSTIME RefTime;
} GLOTIMECORR;


typedef struct _GPS_EPH_RECORD_
{
	short       PRN;
	GPSTIME  	TOC;
	double		ClkBias;
	double		ClkDrift;
	double		ClkDriftRate;
	double		IODE;
	double		Crs;
	double		DetlaN;
	double		M0;
	double		Cuc;
	double		e;
	double		Cus;
	double		SqrtA;
	GPSTIME	    TOE;
	double		Cic;
	double		OMEGA;
	double		Cis;
	double		i0;
	double		Crc;
	double		omega;
	double		OMEGADot;
	double		iDot;
	double		CodesOnL2Channel;
	double		L2PDataFlag;
	double		SVAccuracy;
	double		SVHealth;
	double		TGD;
	double		IODC;
	double		TransTimeOfMsg;
	double		FitInterval;
	double		Spare1;
	double		Spare2;
} GPSEPHREC;

typedef struct _GLONASS_EPH_RECORD
{
	short   SlotNum;          /*  Slot number in sat. constellation */
	GPSTIME RefTime;         /*  Epoch of ephemerides (UTC) */
	double  ClockBias;       /*  SV clock bias (sec) (-TauN)  */
	double  FreqBias;      /* SV relative frequency bias (+GammaN) */
	double  Tk;     /* message frame time (0 <= Tk < 86400 sec of day UTC)*/
	double  Pos[3];       /*  Satellite position [m]  */
	double  Vel[3];       /*  Satellite  velocity[m/s] */
	double  Acc[3];       /*  Satellite  velocity[m/s^2] */
	int     Health;       /*  health (0=OK) */
	double  FreqNum;      /*  frequency number (1-24) */
	double  InforAge;     /*  Age of oper. information (days) */
} GLONASSEPHREC;


/*  定义观测值类型的列表顺序  */
typedef struct _OBS_TYPE_LIST_
{
	int ObsTypeNum;
	OBSTYPE ObsType[MAXOBSTYPENUM];
} OBSTYPELIST;

/* 观测数据及其类型 */
typedef struct _OBSDATA_
{
	bool     Availability; /*  是否可用 */
	OBSTYPE  Type;     /*  观测数据类型  */
	double   Obs;      /*  观测值        */
	int		 LLI;		/*数据标识符，1为周跳，9为粗差*/
} ObsData;

/*  每颗卫星的观测数据定义  */
typedef struct _SATOBSDATA_
{
	short           Prn;
	GNSSSys         System;
	ObsData         Data[MAXOBSTYPENUM];     
	short           Used;            /* 没有星历为0, 正常为1, 粗差为-1, 初始化为1 */ 
	double			satPos[3];			 /*地固系下卫星坐标，这里保存的是B星算出来的，因为B星把A星的覆盖*/
	double			elevation;			/*卫星高度角*/
	double			Ion;				/*倾斜电离层延迟*/
	double			wValue;				/*w检验量*/
} SatObsData;

/*  每个历元的观测数据定义  */
typedef struct _EPOCHOBSDATA_
{
	GPSTIME           Time;
	short             EpochFlag;
	short             SatNum;
	SatObsData        SatObs[MAXCHANNUM];
} EpochObsData;

typedef struct ONESATTA11OBS{
	int		index_satlist;	//satlist的索引
	int		used;			//该卫星是否可用
	int		flag;			//单差观测值是否发生周跳，-1为发生周跳，0为不发生周跳，1为新卫星
	int     PRN;			//卫星的Prn号
	double  CA1;			//基准站的CA码
	double	CA2;			//流动站的CA码
	double  P11;			//流动站的P1码
	double  P21;			//流动站的P1码
	double  P12;			//流动站的P2码
	double  P22;			 //流动站的P2码
	//信噪比
	double  S11;			//1号卫星
	double  S21;			//2号卫星
	double  S12;			//1号卫星
	double  S22;			//2号卫星
	//L1载波
	double  L11;			//1号卫星
	double  L21;			//2号卫星

	//L2载波
	double  L12;			//1号卫星
	double  L22;			//2号卫星
	int	    L1ApplyEpoch;   //L1---观测数据的使用标志
	int	    L2ApplyEpoch;   //L2---观测数据的使用标志
	int     CAApplyEpoch;   //CA---观测数据的使用标志
	int	    P1ApplyEpoch;   //P1---观测数据的使用标志
	int	    P2ApplyEpoch;   //P2---观测数据的使用标志
	//单差观测值
	double	dL1;			//单差L1观测值,单位是周
	double	dL2;			//单差L2观测值,单位是周
	double  dCA;			//单差CA观测值
	double	dP1;			//单差P1观测值
	double	dP2;			//单差P2观测值
	double  satPos[3] ;     //卫星位置坐标
	double  elevation[2];   //elevation1---基准站的高度角标志  elevation2---流动站的高度角标志,单位是弧度rad
	double	map[2];			//基准站与流动站的投影函数
	//单差周跳探测中间变量
	double	Nw_SD;				//宽巷模糊度
	double	sigNw_SD;			//宽巷模糊度的方差
	int		Nw_num_SD;			//宽巷模糊度历元数
	//非差周跳探测中间变量
	double  Nw_A;				//A星宽巷模糊度
	double	sigNw_A;			//A星宽巷模糊度的方差
	int		Nw_num_A;			//A星宽巷模糊度历元数
	double  Nw_B;				//A星宽巷模糊度
	double	sigNw_B;			//A星宽巷模糊度的方差
	int		Nw_num_B;			//A星宽巷模糊度历元数
	//LLI周跳探测量
	int	   LLI_L1A;				//A星L1频段LLI
	int	   LLI_L2A;				//A星L2频段LLI
	int	   LLI_L1B;				//B星L1频段LLI
	int	   LLI_L2B;				//B星L2频段LLI
}OneSat11Obs;
/*  每个历元的共同观测值  */
typedef struct _COMMONOBS_
{
	GPSTIME           Time;//主卫星的时间
	short             ComSatNum;
	int				  PRN[MAXCHANNUM];
	OneSat11Obs       comobs[MAXCHANNUM];
}Common11Obs;
struct ONESATDDBIASE							//已固定的双差模糊度
{
	
	int nonRefPrn;								//非参考卫星的PRN值

	int     ddF1;							//双差载波L1已固定的模糊度
	int     ddF2;							//双差载波L2已固定的模糊度
	
	//SATPOSXYZ  satpos;							//当前的该星的位置

};
struct DDPSEUDOOBS                      //双差观测数据结构体
{
	int				refPrn;			   //参考卫星的PRN值
	int				DDObsNum;         //双差后的模糊度的数目，不一定是ComSatNum-1
    ONESATDDBIASE satddobs[MAXCHANNUM]; //双差观测数据
};
/*  载波相位平滑伪距的中间结构体定义  */
typedef struct _PHASESMOOTHPSEUDORANGE_
{
	short       CurrNum;   /* 当前平滑历元总数 */
	GPSTIME     CurrTime;  /* 当前历元的GPS时间 */

	double      PSC1;      /* 当前历元的CA码平滑结果  */
	double      PSP1;      /* 当前历元的P1码平滑结果  */
	double      PSP2;      /* 当前历元的P2码平滑结果  */
	double      CPL1;      /* 当前历元的L1  */
	double      CPL2;      /* 当前历元的L2  */
	bool        GoodC1;    /* 当前历元CA码平滑状态，true or false */
	bool        GoodP1;    /* 当前历元P1码平滑状态，true or false */
	bool        GoodP2;    /* 当前历元P2码平滑状态，true or false */
	bool        GoodL1;    /* 当前历元L1数据的质量，true or false */
	bool        GoodL2;    /* 当前历元L2数据的质量，true or false */
} PSPRange;


typedef struct _SPACECRAFT_STATE_
{
	GPSTIME   Mjd_GPS;     
	double    Pos[6];   /* position and velocity  */
	double    Acc[3];
} SCState;


/* 
单点定位的卫星列表
*/
typedef struct _SATLIST_PP_
{
	GNSSSys System;
	short   Prn;
	int     Status;        /*  所有卫星的状态: 1为可用, 0为无星历/无数据/高度角低,
						   -1为验前检验含有粗差, -2为验后检验有粗差*/
	int		index;			/*相对于观测值的索引*/
} SATLIST;

/* 每个历元单点定位和测速的结果及其精度指标 */
typedef struct _EPOCH_PP_RESULT_
{
	double Position[3];
	double PosBLH[3];
	double Velocity[3];
	double VelNEU[3];
	double RcvClkOft[2];         /* 0 为GPS钟差, 1为GLONASS钟差  */
	double MatrixFromBLHtoNEU[9];
	double RcvClkSft;
	double PDOP;
	double GDOP;
	double HDOP;
	double VDOP;
	double TDOP;
	double SigmaPos;
	double SigmaVel;
	double Residual[MAXCHANNUM];     /* 伪距观测值残差 */
	SATLIST SatList[MAXCHANNUM];     /* 单点定位卫星列表, 和观测值顺序保持一致 */
	int    Iterator;                 /* 单点定位迭代次数 */
	double Coverage;                 /* 单点定位收敛数据 */
	bool   IsSuccess;                /* 单点定位是否成功, 1为成功, 0为失败 */
	int    SatNumUsed;				/*用的GPS卫星数*/
} PPRESULT;

/* 定位计算(包括单点定位和滤波定位)需要用到的每颗卫星的中间计算结果 */
typedef struct _EPOCH_MIDLLE_RESULT_FOR_ONE_SV_
{
	double SatPos[3];
	double SatVel[3];
	double Elevation;
	double SatClkOft;
	double SatClkSft;
	double Relativity;
	double IonoCorr;
	double TropCorr;
} SATMIDRESULT;

typedef struct _APRIORI_LEO_ORBIT_CLK_
{
	double LeoPos[3];
	double LeoClk[2];

	double OrbitAccuracy;      /* PDOP*Sigma */
	int    Validate;           /* SatNum>=5, Coverage, Accuracy etc. */
} APRIORBIT; 


typedef struct _FINAL_RESULT_OUTPUT_STRUCT_
{
	COMMONTIME GT;              /* 轨道对应的GPS时间, 以通用时格式输出 */
	double     Mjd_UTC;         /* 轨道对应的UTC时间, 以简化儒略日格式输出 */
	double     ECF_Orb[6];      /* 滤波计算的卫星位置和速度,在WGS84地固系 */
	double     ECI_Orb[6];      /* 滤波计算的卫星位置和速度,以J2000.0惯性系 */
	double     MeanElements[6]; /* 平根轨道 a,e,i,OMEGA,omega,M */
	double     OscuElements[6]; /* 瞬时轨道 */
	double     FirstMeanElements[6]; /* 第一类无奇点平根 */
	double     FirstOscuElements[6]; /* 第一类无奇点瞬根 */
	double     da;                   /* 长半轴变化率 */
	double     c1, c2, c3;
	double     Sg0;                  /* 格林尼治恒星时 */
	double     Cov[6];               /* J2000惯性系滤波轨道的精度指标 */

	bool       Valid;               /* 轨道转换的有效性, 平根轨道是否正确 */
} FINALORBIT;
typedef struct _EKF_STATE_INFOMATION_
{
	double			Sigma[2],Step,wuchaMod;
	double			CentBias[6];								/*	相位偏差  A3B3*/
	double			Tao[2];										/*	A、B的Tao*/
	IONOPARA		IonPara;
	APRIORBIT		AprioriState[2];							/*地固系下的先验坐标*/
	DYNMODELPARA	Para[2];
	SCState			IntSate[4];									/*A2B2*/
	GPSTIME			Time;										/* Time       */
	double			StateA[ADIMENSION];							/*参考星A的状态*/
	double			StateB[RELDIMENSION];						/*B星状态，时间更新时用*/
	double			CovA[ADIMENSION*ADIMENSION];				/*参考星A的方差*/
	double			StateRel[RELDIMENSION];						/*相对状态*/
	double			CovRel[RELDIMENSION*RELDIMENSION];			/*相对状态的方差*/
	double			ApriSigma[2];								/* 先验信息的标准差 A1B1*/
	double			PostSigma[2];								/* 测量更新后的标准差 A1B1*/
	double			sigPosPC;
	double			sigPosLC;
	double			StateInECEFA[6];							/* 地固系下的滤波位置状态 */
	double			StateInECEFB[6];							/* 地固系下的滤波位置状态 */
	int				SatNumUsed[2];
	int				comSatNumUsed;
	short			KFConvergency[2];							/* 0表示没有收敛,1为收敛,2为发散  */
	short			IsInitial[2];								/* 0表示滤波没有初始化，1为已经初始化，轨道运行时使用 A1B1*/
	int				PRN[MAXCHANNUM];						    /*两个频段模糊度通道的卫星顺序应该一样未标志的为-999*/
	Common11Obs		CurComObs;									/*当前历元的共同观测数据，包括伪距和载波原始观测值*/ 
	Common11Obs     *wholeComObs;								/*一天总的观测值*/
	//int				ComEpochNum	;								/*上述的历元数*/
	int				eleMask;                                    /*高度角限值*/
	DDPSEUDOOBS     ddOBS;										/*虚拟模糊度观测值   rj 20160727*/
} EKFSTATE;

#endif
