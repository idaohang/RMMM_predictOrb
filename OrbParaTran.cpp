#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "OrbParaTran.h"
#include "RTOD_Const.h"

double GsiDSi, GsiDSi2, GsiDSi45, GsiDSi25;
double GsiDCi, GsiDCi2, DgfEe, DgfDe;
double DgfAa22, PubUb, DgfDee, PubAir, PubAir3;
double PubDSu, PubDCu, PdsDSu2, PdcDCu2, PdsDSu3, PdcDCu3, PdsDSu4, PdcDCu4, PdsDSu5, PdcDCu5;
double PubU, DgfP, DgfCf2, DgfCfb2;
double TsHocoSox[10][10], TsHocoCox[10][10];
double TsHocoClx[20], TsHocoSlx[20];
double ConTeAbh;   //some problem in defining the variable
double SatMoDfxh, SatMoGsbh;
double SatMoB1, SatMoB2, GdShoBt2, GdShoBt4;
int    GdShoN;
double GdShoAk[100], GdShoEk[100], GdShoWk[100], GdShoKuxk[100];
double SatMoBi[7];
double SatMoBi00, SatMoBi20;

const double SggDj[5][5] = { 
	{0.0,0.0,0.0,0.0,0.0},
	{0.0,0.0,0.0,0.0,0.0},
	{1.082626890e-003,0.0,2.807480e-006,0.0,0.0},
	{-2.5356350e-006,2.0441380e-006,1.0888760e-006,1.5765620e-006,0.0},
	{-1.6233590e-006,7.1179670e-007,7.516045000000001e-007,1.0088860e-006,3.5699490e-007}};

const double SggCjc[65]=
	{        0.1000000000000000e+01,  0.1000000000000000e+01,
     		 0.2000000000000000e+01,  0.6000000000000000e+01,
     	     0.2400000000000000e+02,  0.1200000000000000e+03,
     		 0.7200000000000000e+03,  0.5040000000000000e+04,
     	     0.4032000000000000e+05,  0.3628800000000000e+06,
     		 0.3628800000000000e+07,  0.3991680000000000e+08,
     	     0.4790016000000000e+09,  0.6227020800000000e+10,
     		 0.8717829120000000e+11,  0.1307674368000000e+13,
     	     0.2092278988800000e+14,  0.3556874280960000e+15,
     		 0.6402373705728000e+16,  0.1216451004088320e+18,
     	     0.2432902008176640e+19,  0.5109094217170944e+20,
     		 0.1124000727777608e+22,  0.2585201673888498e+23,
     	     0.6204484017332394e+24,  0.1551121004333099e+26,
     		 0.4032914611266056e+27,  0.1088886945041835e+29,
     	     0.3048883446117139e+30,  0.8841761993739702e+31,
     		 0.2652528598121911e+33,  0.8222838654177923e+34,
     	     0.2631308369336935e+36,  0.8683317618811886e+37,
     		 0.2952327990396041e+39,  0.1033314796638614e+41,
     	     0.3719933267899012e+42,  0.1376375309122635e+44,
     		 0.5230226174666011e+45,  0.2039788208119744e+47,
     	     0.8159152832478977e+48,  0.3345252661316381e+50,
     		 0.1405006117752880e+52,  0.6041526306337384e+53,
     	     0.2658271574788449e+55,  0.1196222208654802e+57,
     		 0.5502622159812089e+58,  0.2586232415111682e+60,
     	     0.1241391559253607e+62,  0.6082818640342676e+63,
     		 0.3041409320171338e+65,  0.1551118753287382e+67,
     	     0.8065817517094388e+68,  0.4274883284060026e+70,
     		 0.2308436973392414e+72,  0.1269640335365828e+74,
     	     0.7109985878048635e+75,  0.4052691950487722e+77,
     		 0.2350561331282879e+79,  0.1386831185456898e+81,
     	     0.8320987112741390e+82,  0.5075802138772248e+84,
     		 0.3146997326038794e+86,  0.1982608315404440e+88,
     	     0.1268869321858842e+90
};

const double	SggDbt[71]=
	{   0.1000000e+01,  0.2000000e+01,  0.4000000e+01,
     	0.8000000000e+01,  0.1600000000e+02,  0.3200000000e+02,	
     	0.6400000000e+02,  0.1280000000e+03,  0.2560000000e+03,
     	0.5120000000e+03,  0.1024000000e+04,  0.2048000000e+04,
     	0.4096000000e+04,  0.8192000000e+04,  0.1638400000e+05,
     	0.3276800000e+05,  0.6553600000e+05,  0.1310720000e+06,
     	0.2621440000e+06,  0.5242880000e+06,  0.1048576000e+07,
     	0.2097152000e+07,  0.4194304000e+07,  0.8388608000e+07,
     	0.1677721600e+08,  0.3355443200e+08,  0.6710886400e+08,
     	0.1342177280e+09,  0.2684354560e+09,  0.5368709120e+09,
     	0.1073741824e+10,  0.2147483648e+10,  0.4294967296e+10,
     	0.8589934592000000000e+10,  0.1717986918400000000e+11,
     	0.3435973836800000000e+11,  0.6871947673600000000e+11,
     	0.1374389534720000000e+12,  0.2748779069440000000e+12,
     	0.5497558138880000000e+12,  0.1099511627776000000e+13,
     	0.2199023255552000000e+13,  0.4398046511104000000e+13,
     	0.8796093022208000000e+13,  0.1759218604441600000e+14,
     	0.3518437208883200000e+14,  0.7036874417766400000e+14,
     	0.1407374883553280000e+15,  0.2814749767106560000e+15,
     	0.5629499534213120000e+15,  0.1125899906842624000e+16,
     	0.2251799813685248000e+16,  0.4503599627370496000e+16,
     	0.9007199254740992000e+16,  0.1801439850948198400e+17,
     	0.3602879701896396800e+17,  0.7205759403792793600e+17,
     	0.1441151880758558720e+18,  0.2882303761517117440e+18,
     	0.5764607523034234880e+18,  0.1152921504606846976e+19,
     	0.2305843009213693952e+19,  0.4611686018427387904e+19,
     	0.9223372036854775808e+19,  1.8446744073709551616e+19,
     	3.6893488147419103232e+19,  7.3786976294838206464e+19,
     	1.47573952589676412928e+20, 2.95147905179352825856e+20,
     	5.90295810358705651712e+20, 1.180591620717411303424e+21
};

const double SggDl[5][5]= {
		{0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,2.88075320e0,0.0,0.0},
		{0.0,0.12357890e0,2.83677290e0,0.37017860e0,0.0},
		{0.0,3.86142840e0,0.54163730e0,2.02727150e0,0.53742140e0}
	};

struct PODPara
{
	  double julianday;       // UTC儒略日，定轨时刻
    double t0;              // 北京时儒略日，定轨时刻
    double Sg0;             // 定轨时刻恒星时  单位 度 
    double x;               // J2000.0系X轴方向坐标，单位： 米	
    double y;	            // J2000.0系y轴方向坐标，单位： 米
    double z;	            // J2000.0系z轴方向坐标，单位： 米
    double vx;	            // J2000.0系X轴方向速度，单位： 米/秒
    double vy;	            // J2000.0系y轴方向速度，单位： 米/秒
    double vz;	            // J2000.0系z轴方向速度，单位： 米/秒
   
    ////////////  开普勒瞬根数定义///////////
    double a ;          // 半长轴，单位 米
    double e;           // 偏心率： 
    double i;           // 倾角，单位 度  
    double omiga;       // 近地点角距，单位 度
    double OMIGA;       // 升交点赤经，单位 度
    double M;           // 平近点角，单位 度
   
    ///////////// 开普勒平根数定义 ////////////
    double am ;         // 半长轴，单位 米
    double em;          // 偏心率： 
    double im;          // 倾角，单位 度  
    double omigam;      // 近地点角距，单位度
    double OMIGAm;      // 升交点赤经，单位 度
    double Mm;          // 平近点角，单位 度
   
    //  第一类无奇点瞬根数定义
    double a1;          // 半长轴，单位 米
    double i1;          // 倾角，单位 度
    double OMIGA1;      // 升交点赤经，单位 度 
    double kc1;         // 瞬根ξ
    double ita1;        // 瞬根η
    double lamda1;      // 瞬根λ 单位 度
   
    //  第一类无奇点平根数定义
    double a2 ;         // 半长轴，单位 米
    double i2;          //  倾角，单位 度 
    double OMIGA2;      //  升交点赤经，单位 度
    double kc2;         // 平根ξ
    double ita2;        // 平根η
    double lamda2;      // 平根λ 单位 度
  
    double da ;         // 半长轴变化率  米/天
    double c1;          //
    double c2;          //
    double c3;          //
};

void MMul31(double (*a)[3], double b[3], double c[3])
{
	c[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2];
	c[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2];
	c[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2];
}

//ZTSRMST
void MST(int Day, double SecOfDay, double& gt)
{
	int day0 = Day - 18262;

	if( fabs( SecOfDay ) > 86400. )
	{
		day0 = day0 + int( SecOfDay / 86400. );
		SecOfDay = SecOfDay - int( SecOfDay / 86400. ) * 86400.;
	}
	double ut1 = SecOfDay;

	double D = int( day0 - 1. + ut1 / 86400. - 5. / 6. );

	double x = day0 - 1. + ut1 / 86400. - 5. / 6. - D;

	gt = 280.4606184 + 0.9856122863 * ( D + x ) + 360.0 * x;
	gt = gt * Rad;
	RadianNorm(gt);
}

void RotationMatrix(int k, double ang, double mat[3][3])
{
	double c = cos(ang), s = sin(ang);

	if( ( k < 0 ) || ( k > 2 ) )
	{
		exit(1);
	}

	for(int i=0; i<3; i++)
	{
		mat[k][i] = 0.;
		mat[i][k] = 0.;
	}
	mat[k][k] = 1.;

	if( k == 0 )
	{
		mat[1][1] = c;
		mat[2][2] = c;
		mat[1][2] = s;
		mat[2][1] = -s;
		return;
	}
	else if( k == 1 )
	{
		mat[0][0] = c;
		mat[2][2] = c;
		mat[0][2] = -s;
		mat[2][0] = s;
		return;
	}
	else if( k == 2 )
	{
		mat[0][0] = c;
		mat[1][1] = c;
		mat[0][1] = s;
		mat[1][0] = -s;
		return;
	}
}

//ok
void PV_Keplerian(int flag, double* par1, double* par2, int& status)
{
	double ci, co, cw, si, so, sw, cohf, sihf;
	double tge, tgo, tgw;
	double px, py, pz, qx, qy, qz, rx, ry, rz;
	double a1, e1, e2, r1, r2, r3, r4;
	double r, vv, pp, rae, sn, ps;
	double da, de, ese, ece, ee, f;

	status = 0;
	if( flag == 1 )
	{
		r = sqrt(par1[0]*par1[0] + par1[1]*par1[1] + 
			par1[2]*par1[2]);
		vv = par1[3] * par1[3] + par1[4] * par1[4] + 
			par1[5] * par1[5];
		par2[0] = 1. / ( 2. / r - vv);
		if( par2[0] <= 0. )
		{
			par2[0] = - par2[0];
			pp =   pow( ( par1[0] * par1[4] - par1[1] * par1[3]), 2. )
     	         + pow( ( par1[1] * par1[5] - par1[2] * par1[4]), 2. )
    	         + pow( ( par1[2] * par1[3] - par1[0] * par1[5]), 2. );
			par2[1] = sqrt( pp / par2[0] + 1. );
			cohf = ( 1. + r / par2[0] ) / par2[1];
			sihf = sqrt( cohf * cohf - 1. );
			par2[5] = cohf + sihf;
			par2[5] = log( par2[5] );
			RadianNorm( par2[5] );
			de = sqrt( par2[1] * par2[1] - 1. );
			da = sqrt( par2[0] );
			pz = cohf / r * par1[2] - da * sihf * par1[5];
			qz = ( sihf / r * par1[2] + da * 
				( par2[1] - cohf) * par1[5] ) / de;
			rx = ( par1[1] * par1[5] - par1[2] * par1[4] ) / de / da;
			ry = ( par1[2] * par1[3] - par1[0] * par1[5] ) / de / da;
			rz = ( par1[0] * par1[4] - par1[1] * par1[3] ) / de / da;
		}
		else
		{
			a1 = sqrt(par2[0]);
			ece = 1. - r / par2[0];
			ese = ( par1[0] * par1[3] + par1[1] * par1[4] + par1[2] * par1[5] ) / a1;
			par2[1] = sqrt( ece * ece + ese * ese );
			tge = atan2( ese, ece );
			if( tge < 0. )
				tge = tge + pi2;
			par2[5] = tge - ese;
			RadianNorm( par2[5] );
			pz = ece / par2[1] / r * par1[2] - a1 * ( ese / par2[1] ) * par1[5];
			qz = 1. / sqrt( 1. - par2[1] * par2[1] ) * ( ese / par2[1] / r * par1[2] 
				+ a1 * ( ece / par2[1] - par2[1] ) * par1[5] );
			rae = sqrt( par2[0] * ( 1. - par2[1] * par2[1] ) );
			rx = ( par1[1] * par1[5] - par1[2] * par1[4] ) / rae;
			ry = ( par1[2] * par1[3] - par1[0] * par1[5] ) / rae;
			rz = ( par1[0] * par1[4] - par1[1] * par1[3] ) / rae;
		}
		par2[2] = acos( rz );
		sn = -ry;
		tgo = atan2(rx, sn );
		if( tgo < 0. )
			tgo = tgo + pi2;
		par2[3] = tgo;
		tgw = atan2( pz, qz );
		if( tgw < 0. )
			tgw = tgw + pi2;
		par2[4] = tgw;
	}
	else
	{
		ps = par1[0] * ( 1. - par1[1] * par1[1] );
		if( ps <= 0. )
		{
			status = -1;
			return;
		}
		ps = sqrt( ps );
		ci = cos( par1[2] );
		si = sin( par1[2] );
		co = cos( par1[3] );
		so = sin( par1[3] );
		cw = cos( par1[4] );
		sw = sin( par1[4] );
		px = co * cw - sw * so * ci;
		py = cw * so + sw * co * ci;
		pz = sw * si;
		qx = - sw * co - cw * so * ci;
		qy = - sw * so + cw * co * ci;
		qz = cw * si;
		if( par1[0] >= 0. )
		{
			e1 = par1[1];
			e2 = par1[5];
			FuncKepl(e1, e2 , ee, f, status);
			if( status != 0 )
				status = -3;
			r1 = par1[0] * ( cos(ee) - par1[1] );
			r2 = par1[0] * sqrt( 1. - par1[1] * par1[1] ) * sin(ee);
			r3 = - sin(f) / ps;
			r4 = ( cos(f) + par1[1] ) / ps;
		}
		else
		{
			da = sqrt( - par1[0] );
			cohf = cosh( par1[5] );
			sihf = sinh( par1[5] );
			r1 = - par1[0] * ( par1[1] - cohf );
			r2 = - par1[0] * sqrt( par1[1] * par1[1] - 1. ) * sihf;
			r3 = - sihf / da / ( par1[1] * cohf - 1. );
			r4 = sqrt( par1[1] * par1[1] - 1. ) * cohf / da / ( par1[1] * cohf - 1. );
		}
		par2[0] = r1 * px + r2 * qx;
		par2[1] = r1 * py + r2 * qy;
		par2[2] = r1 * pz + r2 * qz;
		par2[3] = r3 * px + r4 * qx;
		par2[4] = r3 * py + r4 * qy;
		par2[5] = r3 * pz + r4 * qz;
	}
}
//3乘3矩阵相乘
void MMul33(double (*a)[3], double (*b)[3], double (*c)[3])
{
	c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + 
		      a[0][2] * b[2][0];
	c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + 
		      a[0][2] * b[2][1];
	c[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + 
		      a[0][2] * b[2][2];

	c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + 
		      a[1][2] * b[2][0];
	c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + 
		      a[1][2] * b[2][1];
	c[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + 
		      a[1][2] * b[2][2];

	c[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + 
		      a[2][2] * b[2][0];
	c[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + 
		      a[2][2] * b[2][1];
	c[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + 
		      a[2][2] * b[2][2];
}

//ZTSRYMD
void YmdTo90(int year, int month, int day, int& Val)
{
	const static double limit = pow(2., 30.);
	if( month <= 2 )
	{
		month = month + 12;
		year  = year  - 1;
	}

	int Day_Y = int( (year - 1900) * 365.25 );
	int Day_M = int( (month * 153 - 1.5) / 5 );

	int Tol = Day_Y + Day_M + day - 18294;
	if( abs(Tol) <= limit )
		Val = int( Tol );
	else
		Val = 0;
}

//本程序为进行参数转换的主程序入口，此文件中其他函数为
//进行转化所必须的相关子函数
//ymdhm: 年月日时分，存储为整数类型
//sec: 当天中的秒数，需要指出的是这里输入的时间位于何
//时间系统之下还未明确，目前使用的是UTC系统进行测试，发现
//结果与给我们的参考结果基本可以符合，以后可以查询该问题
//osc: 该组参数为开普勒瞬根，长度以米为单位，e以1为单位，角度以角中度为单位
//mean：解算得到的开普勒平根，单位与osc的表达方式一样
//x1_osc：计算得到的第一类无奇点瞬根，单位与前相同
//x1_mean: 计算得到的第一类无奇点平根，单位与前相同
void OsculToMean(int ymdhm[5], double sec, double osc[6], double mean[6], 
				 double x1_osc[6], double x1_mean[6])
{
	int Day;
	double SecOfDay, Norm_A[8], Norm_B[8], Norm_B1[8], Sp[8], Sb[8], So[8];
	int is = 27;

	SecOfDay = sec + ymdhm[3] * 3600. + ymdhm[4] * 60.;
	YmdTo90(ymdhm[0], ymdhm[1], ymdhm[2], Day);

	SecOfDay = SecOfDay + 28800.0;
	if( SecOfDay > 86400. )
	{
		SecOfDay = SecOfDay - 86400.;
		Day = Day + 1;
	}

	///////////////////////////////////////////////////////////////
	//对输入值的单位进行调整，以进行后续计算
	Norm_B[0] = osc[0] / AE;
	Norm_B[1] = osc[1];
	for(int i=2; i<6; i++)
		Norm_B[i] = osc[i] * Rad;
	Norm_B[6] = 0.;
	Norm_B[7] = 0.;
	///////////////////////////////////////////////////////////////

	FuncMixIne(0, Day, SecOfDay, Norm_B, Sp);
	FuncTrKfe(1, Sp, Norm_B1);
	x1_osc[0] = Norm_B1[0] * AE;
	x1_osc[1] = Norm_B1[1] / Rad;
	x1_osc[2] = Norm_B1[2] / Rad;
	x1_osc[3] = Norm_B1[3];
	x1_osc[4] = Norm_B1[4];
	x1_osc[5] = Norm_B1[5] / Rad;

	FuncMast1(is, Day, SecOfDay, Norm_B1, Sb);
	x1_mean[0] = Sb[0] * AE;
	x1_mean[1] = Sb[1] / Rad;
	x1_mean[2] = Sb[2] / Rad;
	x1_mean[3] = Sb[3];
	x1_mean[4] = Sb[4];
	x1_mean[5] = Sb[5] / Rad;
	FuncTrKfe(2, Sb, So);
	FuncMixIne( 1, Day, SecOfDay, So, Norm_A );
	mean[0] = Norm_A[0] * AE;
	mean[1] = Norm_A[1];
	for(int i=2; i<6; i++)
		mean[i] = Norm_A[i] / Rad;
}

// 将角度转化到0到360度范围内
void AngleNorm(double &ang)
{
	ang = ang - int( ang / 360. ) * 360.;
	if( ang < 0. )
		ang = ang + 360.;
}

// 将弧度转化到0到2*PI范围内
void RadianNorm(double& ang)
{
	ang = ang - pi2 * int( ang / pi2 );
	if( ang < 0. )
		ang = ang + pi2;
}

// 将弧度转化到0到2*PI范围内
void RadianNorm2(double& ang)
{
	if( ang < 0. )
		ang = ang + pi2 * ( abs( int( ang / pi2 ) ) + 1 );
	ang = ang - int( ang / pi2 ) * pi2;
	if( ang >= pi )
		ang = ang - pi2;
}

// 具体函数意义未明白，计算结果目前正确
// 此函数在整个流程中调用了两次，包括flag为0和1各一次
// Day：整数天
// SecOfDay: 天内的秒数
// osc: 输入参数
// ggs: 输出结果
void FuncMixIne(int flag, int Day, double SecOfDay, double* osc, double* ggs)
{
	int status;
	double osc0[8], rv[6], r[3], v[3], r1[3], v1[3];

	for(int i=0; i<8; i++)
		osc0[i] = osc[i];

	double eps = Day + SecOfDay / 86400. - 4 / 3.;
	double Time = 2433282.5 + eps;
	double Yea  = ( Time - 2451545. ) / 36525.;

	double lm = 218.316656 + ( 481267.8813425 - 0.001329722E0 * Yea ) * Yea;
	AngleNorm(lm);
	lm = lm * Rad;

	double ls = 280.466447 + ( 36000.76981 + 0.000303611 * Yea ) * Yea;
	AngleNorm(ls);
	ls = ls * Rad;

	double dm = 125.044556 - ( 1934.136185 - 0.002076667 * Yea ) * Yea;
	AngleNorm(dm);
	dm = dm * Rad;

	if( flag != 0 )
		osc0[3] = osc0[3] - 0.6405 * Rad;

	PV_Keplerian(2, osc0, rv, status);

	for(int i=0; i<3; i++)
	{
		r[i] = rv[i];
		v[i] = rv[i+3];
	}

	double Pre[3];
	Pre[0] = ( 2306.2181 + 0.30188 * Yea + 0.017998 * Yea * Yea ) * Yea / 3600. * Rad;
	Pre[1] = ( 2306.2181 + 1.09468 * Yea + 0.018203 * Yea * Yea ) * Yea / 3600. * Rad;
	Pre[2] = ( 2004.3109 - 0.42665 * Yea - 0.041833 * Yea * Yea ) * Yea / 3600. * Rad;

	double sdm = sin( dm );
	double s2lm = sin( 2 * lm );
	double s2ls = sin( 2 * ls );

	double rad_miu = ( -15.813 * sdm + 0.191 * sin( 2 * dm ) - 1.166 * s2ls - 
		             0.187 * s2lm ) / 3600. * Rad;
	double rad_cta = ( -6.86 * sdm - 0.506 * s2ls + 0.083 * sin( 2 * dm ) - 
		             0.081 * s2lm ) / 3600. * Rad;
	double rad_eps = ( 9.21 * cos( dm ) + 0.551 * cos( 2 * ls ) - 0.09 * cos( 2 * dm ) + 
		             0.088 * cos( 2 * lm ) ) / 3600. * Rad;

	eps = eps / 365.2422 / 100.;

	double Par_A = 4609.896 * eps + 1.395 * eps * eps + 0.0362 * eps * eps * eps;
	Par_A = Par_A / 3600. * Rad + rad_miu;

	double rot[3][3], rot0[3][3], rot1[3][3], rot2[3][3];
	double angle;

	if( flag == 0 )
	{
		angle = -Pre[0];
		RotationMatrix(2, angle, rot0);
		angle = Pre[2];
		RotationMatrix(1, angle, rot2);
		MMul33(rot2, rot0, rot);

		angle = -Pre[1];
		RotationMatrix(2, angle, rot2);
		MMul33(rot2, rot, rot0);

		angle = -rad_miu;
		RotationMatrix(2, angle, rot2);
		RotationMatrix(1, rad_cta, rot);
		MMul33(rot, rot2, rot1);

		angle = -rad_eps;
		RotationMatrix(0, angle, rot2);
		MMul33(rot2, rot1, rot);
		MMul33(rot, rot0, rot1);

		RotationMatrix(2, Par_A, rot2);
		MMul33(rot2, rot1, rot0);
		MMul31(rot0, r, r1);
		MMul31(rot0, v, v1);

		for(int i=0; i<3; i++)
		{
			rv[i] = r1[i];
			rv[i+3] = v1[i];
		}
		PV_Keplerian(1, rv, ggs, status);
		ggs[3] = ggs[3] + 0.6405 * Rad;
	}
	else
	{
		angle = -Par_A;
		RotationMatrix(2, angle, rot2);
		RotationMatrix(0, rad_eps, rot);
		MMul33(rot, rot2, rot0);

		angle = -rad_cta;
		RotationMatrix(1, angle, rot1);
		MMul33(rot1, rot0, rot2);
		RotationMatrix(2, rad_miu, rot);
		MMul33(rot, rot2, rot0);

		angle = Pre[1];
		RotationMatrix(2, angle, rot1);
		MMul33(rot1, rot0, rot2);

		angle = -Pre[2];
		RotationMatrix(1, angle, rot);
		MMul33(rot, rot2, rot1);

		angle = Pre[0];
		RotationMatrix(2, angle, rot2);
		MMul33(rot2, rot1, rot0);
		MMul31(rot0, r, r1);
		MMul31(rot0, v, v1);
		for(int i=0; i<3; i++)
		{
			rv[i] = r1[i];
			rv[i+3] = v1[i];
		}
		PV_Keplerian(1, rv, ggs, status);
	}

	ggs[6] = osc[6];
	ggs[7] = osc[7];
}

//ZTSRTRKFE 
// 对函数已经进行了测试，计算结果正确
// 
void FuncTrKfe(int flag, double* Norm_A, double* Norm_B)
{
	double temp, Ang, Ang1, Ang2;

	const static double TM = 1.344686056927252e1;
	if( flag == 1 )
	{
		Norm_B[0] = 1. / pow( Norm_A[0], 1.5 );
		Norm_B[1] = Norm_A[2];
		Norm_B[2] = Norm_A[3];
		Norm_B[3] = Norm_A[1] * cos( Norm_A[4] );
		Norm_B[4] = - Norm_A[1] * sin( Norm_A[4] );
		Norm_B[5] = Norm_A[4] + Norm_A[5];
		RadianNorm( Norm_B[5] );
		Norm_B[6] = - Norm_B[0] * Norm_B[0] / pi2 * 
			Norm_A[6] / 86400.;
		temp = pi2 / Norm_B[0] * TM;
		Norm_B[7] = pi2 / temp / temp * ( 2. * Norm_A[6] * Norm_A[6] / temp 
			- 60. * Norm_A[7] ) * pow( TM, 3. ) / 86400. / 86400.;
	}
	else
	{
		Norm_B[0] = 1. / pow( Norm_A[0], 2./3. );
		Norm_B[1] = sqrt( Norm_A[3] * Norm_A[3] + Norm_A[4] * Norm_A[4] );
		Norm_B[2] = Norm_A[1];
		Norm_B[3] = Norm_A[2];
		Ang2 = Norm_A[3];
		Ang1 = -Norm_A[4];
		Ang = atan2(Ang1, Ang2);
		if( Ang < 0. )
			Ang = Ang + pi2;
		Norm_B[4] = Ang;
		Norm_B[5] = Norm_A[5] - Norm_B[4];
		RadianNorm( Norm_B[5] );

		Norm_B[6] = -86400. * pi2 * Norm_A[6] / Norm_A[0] / Norm_A[0];
		Norm_B[7] = 86400. * 1440. / TM * ( -Norm_A[7] + 2. / Norm_A[0] * Norm_A[6] * Norm_A[6])
			* pi2 / Norm_A[0] / Norm_A[0];
	}
}

//ZTSRMAST1
// 此转换过程中最重要的函数
// flag：具体意义不清楚，在整个流程中只调用了一次，输入值为27
// Day：整数天
// SecOfDay: 天内的秒数
// b1: 
// sb: 
void FuncMast1(int flag, int Day, double SecOfDay, double* b1, double* sb)
{
	double gt;
	double bmp[7], bsp[7];
	double s[8], sa[8], asa[32];
	double ss1[6], ss2[6], ss3[6], ss4[6];
	double six[17], cix[17];
	double nmp[5][5][5], nmpb[5][5][5];
	int N;

	MST(Day, SecOfDay, gt);
	BMoon(Day, SecOfDay, bmp);
	BSun(Day, SecOfDay, bsp);
	FuncTrKfe(2, b1, sa);

	for(int i=0; i<8; i++)
	{
		s[i] = b1[i];
		sb[i] = s[i];
	}

	int jmv = 0;
	asa[0] = s[5];
	double kgs = sa[0];

	do
	{
		FuncGkepl(s, sa);
		FuncGecu(s, sa);

		if( flag == 1 )
		{
			FuncSHO1(s, sa, ss1);
			FuncSHO2(s, sa, ss1[0], ss2);
		}
		else
		{
			FuncSHO1(s, sa, ss1);
			FuncSHO2(s, sa, ss1[0], ss2);
			N = 4;
			six[ 2*N ] = 1.;
			cix[ 2*N ] = 1.;
			for(int i=1; i<=2*N; i++)
			{
				six[ 2*N + i ] = six[ i - 1 + 2*N ] * GsiDSi;
				six[ 2*N - i ] = six[ 1 - i + 2*N ] / GsiDSi;
				cix[ 2*N + i ] = cix[ 2*N - 1 + i ] * ( 1. - GsiDCi );
				cix[ 2*N - i ] = cix[ 2*N - i + 1 ] / ( 1. - GsiDCi );
			}

			FuncQfNmp(N, six, cix, nmp, nmpb);
			FuncGtSco(0., gt, N, s);
			FuncGtSho(0., gt, N, nmp, nmpb, s, sa, ss3);
			FuncGdrto(s, sa);
			//分析发现由于变量s[6]未使用过，其值为0，使得下述变量值为0
			DgfAa22 = -2. * sa[0] / 3. / s[0] * s[6];
			FuncGdrsh(0., bsp[6], s, sa, ss4);
		}

		s[0] = kgs - ( ss1[0] + ss2[0] ) - ss3[0] - ss4[0];
		s[0] = 1. / sqrt( s[0] * s[0] * s[0] );
		for(int i=1; i<=5; i++)
			s[i] = sb[i] - ss1[i] - ss2[i] - ss3[i] - ss4[i];
		FuncTrKfe(2, s, sa);
		jmv = jmv + 1;
		asa[jmv] = s[5];
	}while( ( jmv != 31 ) && ( fabs(asa[jmv] - asa[jmv-1]) >= 1.e-10 ) );

	for(int i=0; i<8; i++)
		sb[i] = s[i];
}

//ZTSRJULIA ok
// 世纪数
double JulianCen(int Day, double SecOfDay)
{
	return ( Day - 18263. + SecOfDay / 86400. - 5./6. ) / 36525.;
}

//ZTSRTROCH ok
double TrochCen(int Day, double SecOfDay)
{
	return ( Day + SecOfDay/86400. - 5./6. - 0.4667 ) / 36524.22;
}

//ZTSRBMOON ok
//计算过程中月球相关参数
void BMoon(int Day, double SecOfDay, double* bm)
{
	double tu = JulianCen(Day, SecOfDay);
	double rad_eps = 23.44 * Rad;
	double rad_smj = 5.140 * Rad;
	double th = TrochCen(Day, SecOfDay);

	double dat50 = 50.24 * th / 3600. * Rad;
	double orm = 4.52360151198 - 33.7571462502 * tu - dat50;
	orm = orm - pi2 * int( orm / pi2 );

	double om = 4.71996656697 + 8399.70914498 * tu;
	om = om - pi2 * int( om / pi2 );

	bm[2] = cos( rad_eps ) * cos( rad_smj ) - sin( rad_eps ) * sin( rad_smj ) * cos( orm );
	bm[1] = sqrt( 1. - bm[2] * bm[2] );
	bm[5] = ( cos(rad_smj) - cos(rad_eps)*bm[2] ) / sin(rad_eps) / bm[1];
	bm[4] = sin(rad_smj) * sin(orm) / bm[1];
	double dsct = sin(rad_eps) * bm[4] / sin(rad_smj);
	double dcct = cos(orm) * bm[5] + sin(orm) * bm[4] * cos(rad_eps);
	double ct   = atan2(dsct, dcct);
	if( ct < 0. )
		ct = ct + pi2;

	bm[0] = atan2(bm[1], bm[2]);
	if( bm[0] < 0. )
		bm[0] = bm[0] + pi2;
	bm[3] = atan2(bm[4], bm[5]);
	if( bm[3] < 0. )
		bm[3] = bm[3] + pi2;

	double rlm = om - orm + ct;
	bm[6] = rlm - pi2 * int( rlm / pi2 );
}

//ZTSRBSUN ok
//计算过程中太阳相关参数
void BSun(int Day, double SecOfDay, double* bs)
{
	double cen = JulianCen(Day, SecOfDay);
	bs[0] = 23.44 * Rad;
	bs[1] = sin( bs[0] );
	bs[2] = cos( bs[0] );
	double mj = 5.14 * Rad;
	double th = ( Day + SecOfDay / 86400. - 5. / 6. - 0.42336 ) / 36524.22;
	double dat50 = 50.24 * th / 3600. * Rad;
	double sun = 4.88162793501 + 628.331950801 * cen - dat50;
	bs[6] = sun - pi2 * int( sun / pi2 );
	bs[3] = 0.;
	bs[4] = sin( bs[3] );
	bs[5] = cos( bs[3] );
}

//ZTSRGKEPL ok
void FuncGkepl(double* s, double* sb)
{
	double PubUb, PubDSub, PubDCub;
	double PubU2, PubU3, PubU4, PubU5;
	PubUb = s[5];
	do
	{
		PubU = PubUb;
		PubUb = s[5] + s[3] * sin(PubU) + s[4] * cos(PubU);
	}while( fabs(PubUb - PubU) > 1.e-10 );

	RadianNorm(PubUb);
	PubDSub = sin(PubUb);
	PubDCub = cos(PubUb);
	DgfDee = 1. / ( 1. + sqrt( 1. - sb[1] * sb[1] ) );
	PubAir = 1. / ( 1. - s[3] * PubDCub + s[4] * PubDSub);
	PubAir3 = pow(PubAir, 3);
	PubDSu = PubAir * ( s[3] * DgfDee * ( s[5] - PubUb ) + PubDSub + s[4] );
	PubDCu = PubAir * ( s[4] * DgfDee * ( s[5] - PubUb ) + PubDCub - s[3] );
	PubU   = atan2(PubDSu, PubDCu);
	RadianNorm(PubU);
	PubDSu = sin( PubU );
	PubDCu = cos( PubU );
	PubU2 = 2. * PubU;
	RadianNorm(PubU2);
	PdsDSu2 = sin(PubU2);
	PdcDCu2 = cos(PubU2);
	PubU3 = 3. * PubU;
	RadianNorm(PubU3);
	PdsDSu3 = sin(PubU3);
	PdcDCu3 = cos(PubU3);
	PubU4 = 4. * PubU;
	RadianNorm(PubU4);
	PdsDSu4 = sin(PubU4);
	PdcDCu4 = cos(PubU4);

	PubU5 = 5. * PubU;
	RadianNorm(PubU5);
	PdsDSu5 = sin(PubU5);
	PdcDCu5 = cos(PubU5);
}

//ZTSRGECU ok
void FuncGecu(double* s, double* sb)
{
	GsiDSi = sin(s[1]);
	GsiDSi2 = GsiDSi * GsiDSi;
	GsiDSi45 = 4. - 5. * GsiDSi2;
	GsiDSi25 = GsiDSi45 / 2.;
	GsiDCi = cos( s[1] );
	GsiDCi2 = GsiDCi * GsiDCi;
	DgfEe = sb[1] * sb[1];
	DgfDe = sqrt( 1. - DgfEe );
	DgfDee = ( 1. - DgfEe ) * DgfDe;
	DgfP   = sb[0] * ( 1. - DgfEe );
	double temp = SggDj[2][0] / DgfP / DgfP * s[0];
	//这些变量在整个程序中只出现了一次，目前并没有利用到 wkliu
//	double DgfDod = -1.5 * temp * GsiDCi;
//	double DgfWd  = 3. / 4. * temp * GsiDsi45;
//	double DgfDlmd = 3. / 4. * temp * ( GsiDsi45 + ( 2. - 3. * GsiDsi2 ) * DgfDe );
	DgfCf2  = ( 3. - 4. * DgfEe ) * DgfEe / ( DgfDee + 1. - 1.5 * DgfEe ) / 2.;
	DgfCfb2 = ( 1. + 9./8.*DgfEe ) / ( DgfDee + 1. - 1.5 * DgfEe + 3. / 8. * DgfEe * DgfEe ) / 4.;
}

//ZTSRGSHO1 ok
void FuncSHO1(double* s, double* sb, double* ss)
{
	double m33 = s[3] * s[3];
	double m44 = s[4] * s[4];
	double m34 = s[3] * s[4];
	double D3344 = m33 - m44;
	double ulm = PubU - s[5];
	RadianNorm2(ulm);
	double jp = 1.5 * SggDj[2][0] / DgfP / DgfP;
	ss[0] = SggDj[2][0] / sb[0] * ( ( 1. - 1.5 * GsiDSi2 ) * ( PubAir3 - 1. / DgfDee )
      	+ 1.5 * PubAir3 * GsiDSi2 * PdcDCu2 );

	ss[1] = jp / 12. * sin( 2. * s[1] ) * ( 3. * PdcDCu2
      	+ 3. * ( s[3] * PubDCu + s[4] * PubDSu ) + s[3] * PdcDCu3 - s[4] * PdsDSu3
        + D3344 * DgfCf2 / DgfEe );

	ss[2] = jp * GsiDCi * ( -ulm + 0.5 * PdsDSu2 - 0.5 * ( s[3] * PubDSu
      	+ 3. * s[4] * PubDCu ) + ( s[3] * PdsDSu3 + s[4] * PdcDCu3 ) / 6. 
      	- m34 / 3. * DgfCf2 / DgfEe );

	ss[3] = jp * (GsiDSi25 * ulm * s[4] + ( 1. + m33 / 4. + 2.25 * m44 ) * PubDCu
      	+ m34 * PubDSu + 0.5 * s[3] * PdcDCu2 - s[4] * PdsDSu2 + ( m33 - 3. * m44 ) / 12. * PdcDCu3
      	- m34 * PdsDSu3 / 3. + ( 1. - DgfCf2 / 12. ) * s[3]
      	- DgfCf2 / DgfEe / 12. * ( m33 - 3. * m44 ) * s[3] + GsiDSi2 * ( ( -5. / 4. 
      	+ 3. * m33 / 8. - 25. * m44 / 8. ) * PubDCu - m34 / 4. * PubDSu 
      	+ 0.5 * s[3] * PdcDCu2 + 2. * s[4] * PdsDSu2 + ( 7. / 12. + 11. / 48. * m33 
      	+ 25. / 48. * m44 ) * PdcDCu3 + 7. / 24. * m34 * PdsDSu3 + 3. / 8. * s[3] * PdcDCu4
      	- 3. / 8. * s[4] * PdsDSu4 + D3344 / 16. * PdcDCu5 - m34 / 8. * PdsDSu5 - ( 5. / 4. 
      	- DgfCf2 / 6. ) * s[3] + ( DgfCf2 / 4. / DgfEe - DgfCfb2 / 6. ) * ( m33 - 3. * m44 ) * s[3] ) );

	ss[4] = jp * ( -GsiDSi25 * ulm * s[3] - ( 1. + 1.25 * m33 + 0.25 * m44 ) * PubDSu 
      	- 2. * m34 * PubDCu - 0.5 * s[4] * PdcDCu2 + DgfEe / 12. * PdsDSu3 
      	+ ( 1. - DgfCf2 / 12. ) * s[4] - DgfCf2 / DgfEe / 12. * ( m44 - 3. * m33 ) * s[4]
      	+ GsiDSi2 * ( ( 7. / 4. + 9. / 8. * m33 + 9. / 8. * m44 ) * PubDSu 
      	+ 13. / 4. * m34 * PubDCu - 0.5 * s[3] * PdsDSu2 + 2. * s[4] * PdcDCu2 
      	- ( 7. / 12. + 13. / 48. * m33 + 23. / 48. * m44 ) * PdsDSu3 
      	+ 5. / 24. * m34 * PdcDCu3 - 3. / 8. * s[3] * PdsDSu4 - 3. / 8. * s[4] * PdcDCu4 
      	- D3344 / 16. * PdsDSu5 - m34 / 8. * PdcDCu5 - ( 5. / 4. - DgfCf2 / 6. ) * s[4] 
      	+ ( DgfCf2 / 4. / DgfEe - DgfCfb2 / 6. ) * ( 3. * m33 - m44 ) * s[4] ) );

	ss[5] = jp * ( GsiDSi25 * ulm - 0.5 * PdsDSu2 + 1.5 * s[3] * PubDSu 
      	+ 2.5 * s[4] * PubDCu - ( s[3] * PdsDSu3 + s[4] * PdcDCu3 ) / 6. 
      	+ ( ( 1 - DgfEe / 4. ) * ( s[3] * PubDSu + s[4] * PubDCu ) + D3344 * PdsDSu2 / 2. 
      	+ m34 * PdcDCu2 + ( m33 - 3. * m44 ) / 12. * s[3] * PdsDSu3 + ( 3. * m33 
      	- m44 ) / 12. * s[4] * PdcDCu3 ) / ( 1. + DgfDe ) + m34 / 3. * DgfCf2 / DgfEe );

	ss[5] = ss[5] + jp * GsiDSi2 * (1.25 * PdsDSu2 - 2.5 * ( s[3] * PubDSu 
      	+ s[4] * PubDCu ) - 5. / 6. * m34 * DgfCf2 / DgfEe + ( ( -0.5 + 1.25 * DgfDe + D3344 / 8. )
      	* s[3] * PubDSu + ( -2.5 - 1.25 * DgfDe + ( 7. * m33 + 5. * m44 ) / 8. )
      	* s[4] * PubDCu - 3. / 4. * D3344 * PdsDSu2 - 1.5 * m34 * PdcDCu2 + ( 1. + 5. / 12. * DgfDe 
      	- 7. / 48. * m33 + 17. / 48. * m44 ) * s[3] * PdsDSu3 + ( 1. + 5. / 12. * DgfDe 
      	- 19. / 48. * m33 + 5. / 48. * m44 ) * s[4] * PdcDCu3 
      	+ 3. / 8. * D3344 * PdsDSu4 + 3. / 4. * m34 * PdcDCu4 + ( m33 - 3. * m44 ) / 16. 
      	* s[3] * PdsDSu5 + ( 3. * m33 - m44 ) / 16. * s[4] * PdcDCu5
      	+ ( 0.25 + DgfCf2 / 3. / DgfEe + DgfCf2 / 6. ) * m34 ) / ( 1. + DgfDe ) );
}

//ZTSRGSHO2 ok
void FuncSHO2(double* s,double* sb,double as1,double* ss)
{
	double a = -0.5 + 1.5 * GsiDCi2;
	double b = 1.5 * GsiDSi2;
	double ab1 = 40.5 * a * a + 31. * a * b + 10. * b * b / 3.;
	double ab2 = -4.5 * a * a + 65. * a * b / 12. - 23. * b * b / 9.;
	double ab3 = -7. * a * b / 4. - b * b / 3.;
	double f1 = 1. - 8. * GsiDCi2 + 7. * GsiDCi2 * GsiDCi2;
	double f2 = 3. - 30. * GsiDCi2 + 35. * GsiDCi2 * GsiDCi2;
	double f3 = GsiDSi2 * GsiDSi2;
	double as2 = 9. / 32. * SggDj[2][0] * SggDj[2][0] / pow( sb[0], 3. ) * f3;
	double p2 = SggDj[2][0] * SggDj[2][0] / 8. / pow( DgfP, 4. );
	double p3 = SggDj[3][0] / 8. / DgfP / DgfP / DgfP;
	double p4 = SggDj[4][0] / pow( DgfP, 4. );

	ss[0] = p2 * DgfP * ( 4. * b * ( 5. - 6. * b ) * PdcDCu2 + 5. / 3. * b * b * PdcDCu4)
      	+ as1 * as1 / 4. / sb[0] - as2 + SggDj[3][0] / 4. / sb[0] / sb[0] * GsiDSi 
      	* ( 3. * ( 5. * GsiDCi2 - 1. ) * PubDSu + 5. * GsiDSi2 * PdsDSu3 )
      	- 5. / 8. * SggDj[4][0] / sb[0] / sb[0] / sb[0] 
      	* ( 3. / 20. * f2 * ( PubAir3 * PubAir * PubAir - ( 1. + 1.5 * DgfEe ) / ( 1. - DgfEe ) 
     	/ ( 1. - DgfEe ) / DgfDee ) - f1 * PdcDCu2 + 7. / 4. * f3 * PdcDCu4);

	ss[1] = p2 * b * b / tan( s[1] ) * ( -2. * PdcDCu2 - ( 4. / 3. 
      	+ 0.5 / tan( s[1] ) / tan( s[1] ) ) * PdcDCu4 ) + p3 * GsiDCi * ( ( 5. 
      	* GsiDCi2 - 1. ) * PubDSu + 5. * GsiDSi2 * PdsDSu3 ) + 5. / 16. * p4 / tan( s[1] ) 
      	* ( f1 * PdcDCu2 - 7. / 4. * f3 * PdcDCu4 );

	ss[2] = 3. * p2 * GsiDCi * ( 3. * GsiDSi2 * PdsDSu2 - ( 3. + GsiDSi2 ) / 4. * PdsDSu4 )
      	- p3 / tan( s[1] ) * ( 3. * ( 15. * GsiDCi2 - 11. ) * PubDCu + 5. * GsiDSi2 * PdcDCu3 ) 
      	- 5. / 8. * p4 * GsiDCi * ( ( 7. * GsiDCi2 - 4. ) * PdsDSu2 
      	+ 7. / 8. * GsiDSi2 * PdsDSu4 );

	ss[3] = p2 * ( ab1 * PubDCu + ab2 * PdcDCu3 + ab3 * PdcDCu5 ) 
      	+ SggDj[3][0] / 4. / sb[0] / sb[0] / sb[0] * GsiDSi * ( ( 10. * GsiDCi2 - 1. ) 
      	* PdsDSu2 + 25. / 8. * GsiDSi2 * PdsDSu4 ) + SggDj[4][0] / 32. / sb[0] / sb[0] 
      	/ sb[0] / sb[0] * ( -5. * ( f1 + 1.5 * f2 ) * PubDCu + 5. * ( 3. * f1 
      	- 7. / 4. * f3 ) * PdcDCu3 - 91. / 4. * f3 * PdcDCu5 );

	ss[4] = - p2 * ( ab1 * PubDSu + ab2 * PdsDSu3 + ab3 * PdsDSu5 ) 
      	+ SggDj[3][0] / 8. / sb[0] / sb[0] / sb[0] * GsiDSi * ( ( 25. * GsiDCi2 - 7. ) 
      	* PdcDCu2 + 25. / 4. * GsiDSi2 * PdcDCu4 ) + SggDj[4][0] / 32. / sb[0] / sb[0]
      	/ sb[0] / sb[0] * ( -5. * ( f1 - 1.5 * f2 ) * PubDSu - 5. * ( 3. * f1 
      	+ 7. / 4. * f3 ) * PdsDSu3 + 91. / 4. * f3 * PdsDSu5 );

	ss[5] = 3. * p2 / 4. * ( ( -7. - 26. * GsiDCi2 + 33. * GsiDCi2 * GsiDCi2 ) * PdsDSu2 
      	+ 0.5 * ( -8. + 27. * GsiDCi2 - 13. * GsiDCi2 * GsiDCi2 ) * PdsDSu4 ) - p3 / GsiDSi 
      	* ( -3. * ( 5. - 41. * GsiDCi2 + 40. * GsiDCi2 * GsiDCi2 ) * PubDCu + 5. / 3. 
      	* ( 5. - 13. * GsiDCi2 + 8. * GsiDCi2 * GsiDCi2 ) * PdcDCu3 ) + 5. / 32. * p4
      	* ( ( 7. - 72. * GsiDCi2 + 77. * GsiDCi2 * GsiDCi2 ) * PdsDSu2 - 7. / 8. * ( 7. - 18. 
      	* GsiDCi2 + 11. * GsiDCi2 * GsiDCi2 ) * PdsDSu4 );
}

//ZTSRQFNMP ok
void FuncQfNmp(int n,double six[17],double cix[17],double nmp[5][5][5],double nmpb[5][5][5])
{
	int num, kmi, kma;
	int j, k, lk;
	double ff, ffcc, fnm, fnmb, bz, cgm;
	for( int i=2; i<=n; i++ )
	{
		num = 1;
		if( i == 2 )
			num = 2;
		for( j=num; j<=i; j++)
		{
			for( k=0; k<=i; k++ )
			{
				ff = sqrt( 2. * SggCjc[i+j] * SggCjc[i-j] * ( 2 * i + 1 ) );
				ffcc = ff / SggDbt[i] / SggDbt[i] / SggCjc[k] / SggCjc[i-k] * SggCjc[2*(i-k)]
					* SggCjc[2*k];
				if( i-j-2*k <= 0 )
					kmi = 0;
				else
					kmi = i - j - 2 * k;
				if( ( 2 * i - 2 * k ) >= ( i - j ) )
					kma = i - j;
				else
					kma = 2 * i - 2 * k;

				fnm = 0.;
				fnmb = 0.;
				for(lk = kmi; lk <= kma; lk++)
				{
					if( lk&1 )  //?
						bz = -1;
					else
						bz = 1;
					cgm = bz / SggCjc[2*i-2*k-lk] / SggCjc[i-j-lk] / SggCjc[lk-i+j+2*k]
						/ SggCjc[lk] * six[3*i-j-2*k-2*lk+2*n] * cix[-2*i+j+2*k+2*lk+2*n];
					fnm = fnm + cgm;
					fnmb = fnmb + cgm * ( -i+j+2*k+2*lk-i*cix[2*n+1] ) / six[2*n+1];
				}
				nmp[i][j][k] = ffcc * fnm;
				nmpb[i][j][k] = ffcc * fnmb;
			}
		}
	}
}

//ZTSRGTSCO ok
void FuncGtSco(double dt, double s0, int n, double* s)
{
	double enw = 0.05883363336;
	double dob, dobj;

	for(int i=2; i<=n; i++)
	{
		int j1 = 1;
		if( i == 2 )
			j1 = 2;
		for(int j=j1; j<=i; j++)
		{
			dob = s[2] - SggDl[i][j] - enw * dt - s0;
			dobj = dob * j;
			RadianNorm(dobj);
			TsHocoSox[i][j] = sin(dobj);
			TsHocoCox[i][j] = cos(dobj);
		}
	}

	int n1 = - n - 1;
	double dlm;
	for(int  k=n1; k<=n+1; k++ )
	{
		dlm = k * s[5];
		RadianNorm(dlm);
		TsHocoClx[k-n1] = cos(dlm);
		TsHocoSlx[k-n1] = sin(dlm);
	}
}

//ZTSRGTSHO ok
void FuncGtSho(double dt, double s0, int n, double nmp[5][5][5], 
			   double nmpb[5][5][5], double* s, double* sb, double* ss)
{
	double enw = 0.05883363336;
	double ar  = enw / s[0];
	int j1, idt, idtcor, j, ijbz, ijdbz, l;
	double Dobj1, Dobj2, Dobj3, Dobj4, gcc;
	double DcGij1, DsGij1, DcGij2, DsGij2, DcGij3, DsGij3;
	double ab21, ab22, ab31, ab32, ab61, ab62, ab41, ab42, ab43, ab51, ab52, ab53;
	double abc1, abc2, abc3, abc4, abc5, abc6, abc7, abc8;

	for(int i=0; i<6; i++)
		ss[i] = 0.;

	for(int i=2; i<=n; i++)
	{
		j1 = 1;
		if( i ==  2 )
			j1 = 2;
		for( j=j1; j<=i; j++ )
		{
			Dobj1 = s[3] * TsHocoCox[i][j] - s[4] * TsHocoSox[i][j];
			Dobj2 = s[4] * TsHocoCox[i][j] + s[3] * TsHocoSox[i][j];
			Dobj3 = s[4] * TsHocoCox[i][j] - s[3] * TsHocoSox[i][j];
			Dobj4 = s[3] * TsHocoCox[i][j] + s[4] * TsHocoSox[i][j];
			ijbz = i - j;
			if( ijbz&1 )  //?
				idt = 1;
			else
				idt = 0;
			ijdbz = (i-j+idt) / 2;
			if(ijdbz&1)
				idtcor = -1;
			else
				idtcor = 1;
			gcc = SggDj[i][j] / pow( sb[0], i ) * idtcor;
			if( i&1 )
			{
				ab21 = ( GsiDCi - j ) * nmp[i][j][(i-1)/2] * ( (1-idt)*Dobj4 - idt*Dobj3 );
				ab22 = ( GsiDCi + j ) * nmp[i][j][(i+1)/2] * ( (1-idt)*Dobj1 + idt*Dobj2 );
				ab31 = nmpb[i][j][(i-1)/2] * ( (1-idt)*Dobj3 + idt*Dobj1 );
				ab32 = nmpb[i][j][(i+1)/2] * ( (1-idt)*Dobj2 - idt*Dobj1 );
				ab61 = nmp[i][j][(i-1)/2] * ( (1-idt)*Dobj3 + idt*Dobj1 );
				ab62 = nmp[i][j][(i+1)/2] * ( (1-idt)*Dobj2 - idt*Dobj1 );
				ab41 = 0.;
				ab42 = 0.;
				ab43 = 0.;
				ab51 = 0.;
				ab52 = 0.;
				ab53 = 0.;
			}
			else
			{
				ab21 = 0.;
				ab22 = 0.;
				ab31 = 0.;
				ab32 = 0.;
				ab61 = 0.;
				ab62 = 0.;
				ab41 = nmp[i][j][i/2+1] * ( (1-idt)*Dobj1 + idt*Dobj2 );
				ab42 = nmp[i][j][i/2-1] * ( -(1-idt)*Dobj4 + idt*Dobj3 );
				ab43 = nmp[i][j][i/2] * ( (1-idt)*TsHocoSox[i][j] - idt*TsHocoCox[i][j] ) * i * (i+1)
					/ j / ar;
				ab51 = nmp[i][j][i/2+1] * ( (1-idt)*Dobj2 - idt*Dobj1 );
				ab52 = nmp[i][j][i/2-1] * ( (1-idt)*Dobj3 + idt*Dobj4 );
			}
			abc1 = 0.;
			abc2 = 0.;
			abc3 = 0.;
			abc4 = 0.;
			abc5 = 0.;
			abc6 = 0.;
			abc7 = 0.;
			abc8 = 0.;

			for(int l=0; l<=i; l++)
			{
				DcGij1 = TsHocoClx[i-2*l+n+1] * TsHocoCox[i][j] - TsHocoSlx[i-2*l+n+1] * TsHocoSox[i][j];
				DsGij1 = TsHocoSlx[i-2*l+n+1] * TsHocoCox[i][j] + TsHocoClx[i-2*l+n+1] * TsHocoSox[i][j];
				DcGij2 = TsHocoClx[i-2*l+n+2] * TsHocoCox[i][j] - TsHocoSlx[i-2*l+n+2] * TsHocoSox[i][j];
				DsGij2 = TsHocoSlx[i-2*l+n+2] * TsHocoCox[i][j] + TsHocoClx[i-2*l+n+2] * TsHocoSox[i][j];
				DcGij3 = TsHocoClx[i-2*l+n] * TsHocoCox[i][j] - TsHocoSlx[i-2*l+n] * TsHocoSox[i][j];
				DsGij3 = TsHocoSlx[i-2*l+n] * TsHocoCox[i][j] + TsHocoClx[i-2*l+n] * TsHocoSox[i][j];

				abc1 = abc1 + 2. * (i-2*l) / (i-2*l-j*ar) * nmp[i][j][l] * ( (1-idt)*DcGij1 + idt*DsGij1);
				abc2 = abc2 + ( (i-2.*l)*GsiDCi - j ) / ( i - 2. * l - j * ar ) * nmp[i][j][l] * 
					( (1-idt)*DcGij1 + idt*DsGij1 );
				abc3 = abc3 + 1. / ( i - 2. * l - j * ar ) * nmpb[i][j][l] * 
					( (1-idt)*DsGij1 - idt * DcGij1 );
				abc4 = abc4 + (3.*i - 4.*l + 1) / (i - 2.*l + 1 - j*ar) * nmp[i][j][l] * 
					( (1-idt)*DcGij2 + idt*DsGij2 );
				abc5 = abc5 + (4.*l - i + 1) / (i - 2.*l - 1 - j*ar) * nmp[i][j][l] * 
					( (1-idt)*DcGij3 + idt*DsGij3 );
				abc6 = abc6 + (3.*i - 4.*l+1) / (i - 2.*l + 1 - j*ar) * nmp[i][j][l] * 
					( -(1-idt)*DsGij2 + idt*DcGij2 );
				abc7 = abc7 + (4.*l - i + 1) / (i - 2.*l - 1 - j*ar) * nmp[i][j][l] * 
					( (1-idt)*DsGij3 - idt*DcGij3 );
				abc8 = abc8 + ( 2. * (i + 1) / ( i - 2. * l - j * ar ) - 3. * ( i - 2 * l )
      			  / ( i - 2. * l - j * ar ) / ( i - 2. * l - j * ar ) ) * nmp[i][j][l]
      			  * ( ( 1 - idt ) * DsGij1 - idt * DcGij1 );
			}
			ss[0] = ss[0] + gcc * sb[0] * abc1;
			ss[1] = ss[1] + gcc / GsiDSi * ( abc2 - (i-1.) * (ab21-ab22) / 2. / j / ar);
			ss[2] = ss[2] + gcc / GsiDSi * ( abc3 - (i-1.) * (ab32-ab31) / 2. / j / ar);
			ss[3] = ss[3] + gcc / 2. * ( abc4 - abc5 - (i-1.)*(i-2.)*(ab41+ab42)/2./j/ar - s[4] * ab43);
			ss[4] = ss[4] + gcc / 2. * ( abc6 - abc7 + (i-1.)*(i-2.)*(ab51-ab52)/2./j/ar + s[3] * ab43);
			ss[5] = ss[5] + gcc * ( abc8 - 1. * (i+1) / j / ar * ( ab62 - ab61 ) );
		}
	}
	ss[3] = ss[3] - s[4] * GsiDCi * ss[2];
	ss[4] = ss[4] + s[3] * GsiDCi * ss[2];
	ss[5] = ss[5] - GsiDCi * ss[2];
}

//ZTSRGDRT0
void FuncGdrto(double* s,double* sb)
{
	double SatMoBh = ConTeAbh;   //maybe some problem
	SatMoDfxh = 0.;
	double ern = 0.5883363904e-1;
	double rp0 = sb[0] * ( 1. - sb[1] );
	double vp0 = sqrt( (1. + sb[1]) / rp0 );
	double df = 1. - rp0 * ern * cos( sb[2] ) / vp0;
	double bhb, fxh;
	double xbi[1004], dak[1000], dek[1000], dwk[1000], duk[1000];
	double kwk[1000], kuk[1000], uk[1000], bi[7];
	double SatMoGsbh2, m33, m44, D3344;
	double z, bt, bt3, obt, obtd, dbt, dbtd, dbtd2;
	int    k;
	
	if( SatMoBh == 0. )
	{
		bhb = ( rp0 - ( 1. - sin( sb[2] ) * sin( sb[2] )
      	    * sin( sb[4] ) * sin( sb[4] ) / 298.257 ) ) * 6378.140;

		if( bhb<260. ) 
			SatMoBh = 43.49 + 0.13679590 * ( bhb - 260. )
      	              - 0.86321009e-3 * ( bhb - 260. ) * ( bhb - 260. );

	    if( ( bhb >= 260. ) && ( bhb < 480. ) ) 
			SatMoBh = 43.49 + 0.13400589 * ( bhb - 260. ) 
			          - 0.24863434e-3 * ( bhb - 260. ) * ( bhb - 260. );

	    if( ( bhb >= 480. ) && ( bhb < 900. ) ) 
			SatMoBh = 61.85 + 5.9483895 * ( exp( ( bhb - 500. ) * log( 2. ) / 99.5 ) - 1.0 );

	    if( ( bhb >= 900. ) && ( bhb < 1000. ) ) 
			SatMoBh = 152.41 + 0.6072909 * ( bhb - 900. );
	    if( ( bhb >= 1000. ) && ( bhb < 1300. ) ) 
			SatMoBh = 213.14 - 121.04047 * ( exp( ( bhb - 1000. ) * 
			          log( 0.609275 ) / 100. ) - 1. );
	    if( bhb >= 1300. ) 
			SatMoBh = 306.79 + 0.12319071 * ( bhb - 1300. );
	    if( bhb < 512.8672 ) 
			fxh = 2.32 + 0.4 * ( bhb - 350. ) / 50. + 0.03 * pow( ( ( bhb - 350. ) / 50. ), 2. ) 
			     + 0.0017 * pow( ( ( bhb - 350. ) / 50. ), 3. );
	    if( bhb >= 512.8672 ) 
			fxh = 4.;
		SatMoDfxh = ( fxh - 1. ) / ( fxh + 1. );
	}

	SatMoBh = SatMoBh / AE * 1000.;
	SatMoGsbh = sb[0] / SatMoBh;
	SatMoGsbh2 = SatMoGsbh * SatMoGsbh;
	z = SatMoGsbh * sb[1];
	SatMoB1 = df * df;
	SatMoB2 = ern * df;
	bt = sb[1] / ( 1. + sqrt( 1. - sb[1] * sb[1] ) );
	GdShoBt2 = bt * bt;
	bt3 = GdShoBt2 * bt;
	GdShoBt4 = GdShoBt2 * GdShoBt2;
	obt = ( 1. - GdShoBt2 ) * ( 1. - GdShoBt2 );
	obtd = ( 1. - bt ) * ( 1. - bt );
	dbt = 2. - bt;
	dbtd = ( 1. + bt ) * ( 1. + bt );
	dbtd2 = dbtd * dbtd;

	//在这里开始，就体现出由于变量ConTeAbh未初始化，而使得
	//计算的GdShoN出现问题
	GdShoN = anint( sqrt( 2. * sb[0] * sb[1] / SatMoBh ) );
	int N = GdShoN;

	// 由于计算的N可能存在问题，xbi存在问题，通过women的检验
	//xbi变量并没有并初始化足够的长度，数据不足
	FuncGbese( N, z, xbi );

	if( z > 100. )
	{
		for( k=0; k<=GdShoN; k++ )
		{
			dak[k] = -8. * bt3 * dbt * ( 1. - k * k / z ) * SatMoBh / rp0 * xbi[k+3];
			dek[k] = -2. * GdShoBt2 * dbt * ( 1. - k * k / z ) * SatMoBh / rp0 * xbi[k+3];
			GdShoAk[k] = 2. * dbtd2 * xbi[k+3] + dak[k];
			GdShoEk[k] = 2. * dbtd * xbi[k+3] + dek[k];
		}
		for( k=1; k<=GdShoN; k++ )
		{
			dwk[k] = -2. / z * GdShoBt2 * dbt * ( 3. - k * k / z ) * SatMoBh / rp0 * xbi[k+3];
			duk[k] = -15. / z * GdShoBt2 * obtd / dbtd * ( 1. - 2. *k * k / 3. / z + pow(k, 4.) / 15. / z /z )
				* SatMoBh * SatMoBh / rp0 / rp0 * xbi[k+3];
			GdShoWk[k] = 2. / z * dbtd * xbi[k+3] + dwk[k];
			kwk[k] = k * GdShoWk[k];
			uk[k] = 2. / z / obt / xbi[k+3] + duk[k];
			kuk[k] = k * uk[k];
		}
	}
	else
	{
		for( k=0; k<=GdShoN; k++ )
		{
			//在这里由于xbi的问题，计算得到的GdShoAk，GdShoEk存在问题
			GdShoAk[k] = 2. * xbi[k+3] + 4. * bt * ( xbi[k+2] + xbi[k+4] )
      			+ 3. * GdShoBt2 * ( xbi[k+1] + 2. * xbi[k+3] + xbi[k+5] )
      			+ 2. * bt3 * ( xbi[k] + xbi[k+2] + xbi[k+4] + xbi[k+6])
      			+ 2. * GdShoBt4 * xbi[k+3];
			GdShoEk[k] = xbi[k+2] + xbi[k+4] + bt * ( xbi[k+1] + 2. * xbi[k+3] + xbi[k+5] )
      			+ GdShoBt2 / 2. * ( 2. * xbi[k] + xbi[k+2] + xbi[k+4] );
		}
		for( k=1; k<=GdShoN; k++ )
		{
			//计算得到GdShoWk存在问题
			kwk[k] = xbi[k+2] - xbi[k+4] + bt * ( xbi[k+1] - xbi[k+5] )
				+ GdShoBt2 / 2. * ( - xbi[k+2] + xbi[k+4] );
			GdShoWk[k] = kwk[k] / k;
			kuk[k] = xbi[k+2] - xbi[k+4] + GdShoBt4 * ( xbi[k+2] - xbi[k+4] )
				- GdShoBt2 / 2. * ( xbi[k] + xbi[k+2] - xbi[k+4] - xbi[k+6] );
		}
	}

	for( k=1; k<=GdShoN; k++ )
		GdShoKuxk[k] = sb[0] * ( kuk[k] - 0.5 * ( 1. - GdShoBt4 ) * kwk[k] );
	FuncGbesq(z, bi);
	for(k=0; k<=6; k++)
		SatMoBi[k] = bi[k];
	m33 = s[3] * s[3];
	m44 = s[4] * s[4];
	D3344 = m33 - m44;
	SatMoBi00 = SatMoBi[0] - SatMoBi[2];
	SatMoBi20 = 2. * ( 1. + SatMoGsbh ) * SatMoBi[0] + 2. * ( 3. - SatMoGsbh ) * SatMoBi[2];

	//variable nt used temporarily
//	SatMoBi20b = 2. * ( 1. + SatMoGsbh ) * SatMoBi[0] + 2. * ( 1. - SatMoGsbh ) * SatMoBi[2];
	double SatMoBi40 = SatMoBi[0] - 4. * SatMoBi[2] / 3. + SatMoBi[4] / 3.;
//	double SatMoBi42 = SatMoGsbh * D3344 * ( SatMobi[2] - SatMoBi[4] ) / 3. + 
//		        SatMoGsbh2 * D3344 * SatMoBi40 / 4.;
//	double SatMoBi44 = 2. * SatMoGsbh / 3. * ( SatMoBi[2] - SatMoBi[4] ) + SatMoGsbh2 / 2. * SatMoBi40;
//	double SatMoBi61 = SatMoBi[0] - 3. * SatMoBi[2] / 2. + 3. * SatMoBi[4] / 5. - SatMoBi[6] / 10.;
//	double SatMoBi60 = SatMoGsbh2 * ( m33 - 3. * m44 ) / 24. * SatMoBi61;
//	double SatMoBi60b = SatMoGsbh2 * ( m44 - 3. * m33 ) / 24. * SatMoBi61;
}

//ZTSRGBESE
void FuncGbese(int n, double z, double* bi)
{
	double Pi[100];
	double bik, bikb, bi1;
	int i, k, M;

	if( z == 0. )
	{
		for(i=0; i<8; i++)
			bi[i] = 0.;
		bi[3] = 1.;
	}
	else
	{
		if( z < 2. )
		{
			if( n >= 1 )
			{
				bi[3] = bi[4] = 0.;
				for( k=0; k<=64; k++ )
				{
					bik = pow( z/2., 2*k ) / SggCjc[k] / SggCjc[k];
					bikb = pow( z/2., 2*k+1 ) / SggCjc[k] / SggCjc[k] / (k+1);
					bi[3] = bi[3] + bik;
					bi[4] = bi[4] + bikb;
					if( ( fabs(bik) < 1.e-5 ) && ( fabs(bikb) < 1.e-5 ) )
						break;
				}
				bi[3] = bi[3] * exp(-z);
				bi[4] = bi[4] * exp(-z);
				for(i=2; i<=n; i++)
					bi[i+3] = bi[i+1] - bi[i+2] / z * 2. * (i-1);
				bi[2] = bi[4];
				bi[1] = 0.;
				bi[0] = 0.;
			}
			else
			{
				bi[3] = 0.;
				bi1 = 0.;
				for( k=0; k<=64; k++ )
				{
					bik = pow( z/2., 2*k ) / SggCjc[k] / SggCjc[k];
					bikb = pow( z/2., 2*k+1 ) / SggCjc[k] / SggCjc[k] / (k+1);
					bi[3] = bi[3] + bik;
					bi1 = bi1 + bikb;
					if( ( fabs(bik) < 1.e-5 ) && ( fabs(bikb) < 1.e-5 ) )
						break;
				}
				bi[3] = bi[3] * exp(-z);
				bi1 = bi1 * exp(-z);
				bi[2] = bi1;
				bi[1] = 0.;
				bi[0] = 0.;
			}
		}
		else
		{
			if( ( z >= 2. ) & ( z <= 100. ) )
			{
				M = anint( dmax1( z, n + 4. ) );
				Pi[M-1] = ( sqrt(z*z + M*M) - M ) / z;
				for( i=M-1; i>=1; i--)
					Pi[i-1] = 1. / ( 2. * i / z + Pi[i] );
				bi[3] = ( 1. + 1. / 8. / z + 9. / 128. / z / z + 25. / 1024. / z / z / z )
					/ sqrt( z * pi2 );
				for(i=1; i<=n; i++)
					bi[i+3] = bi[i+2] * Pi[i-1];
				bi[2] = bi[4];
				if( n<2 )
					bi[1] = 0.;
				else
					bi[1] = bi[5];
				if( n<3 )
					bi[0] = 0.;
				else
					bi[0] = bi[6];
			}
			else
			{
				if( z > 100. )
				{
					for(i=0; i<=n; i++)
						bi[i+3] = 1. * exp(-i*i/2./z) / sqrt(pi2*z);
					bi[2] = bi[4];
					bi[1] = bi[5];
					bi[0] = bi[4];
				}
			}
		}
	}
}

//ZTSRGBESQ else中的if部分分析结果正确，剩余部分通过比较代码无问题
void FuncGbesq(double z, double* bi)
{
	double bik, bikb, Pi[100];
	int i, k, n;

	if( z == 0. )
	{
		bi[0] = 1.;
		for( i=1; i<=6; i++ )
			bi[i] = 0.;
	}
	else
	{
		if( z < 2. )
		{
			bi[0] = 0.;
			bi[1] = 0.;
			for(k=0; k<=64; k++)
			{
				bik = pow( z/2., 2*k ) / SggCjc[k] / SggCjc[k];
				bikb = pow( z/2., 2*k+1 ) / SggCjc[k] / SggCjc[k] / (k+1);
				bi[0] = bi[0] + bik;
				bi[1] = bi[1] + bikb;
				if( ( fabs(bik) < 1.e-5 ) && ( fabs(bikb) < 1.e-5 ) )
					break;
			}
			bi[0] = bi[0] * exp(-z);
			bi[1] = bi[1] * exp(-z);
			for(i=2; i<=6; i++)
				bi[i] = bi[i-2] - bi[i-1] / z * 2. * (i-1);
		}
		else
		{
			if( ( z >= 2. ) & ( z <= 100. ) ) //根据分析这是&与&&其实结果相同 wkliu
			{
				n = anint( dmax1( z, 10. ) );
				Pi[n] = ( sqrt(z*z + n*n) - n ) / z;
				for( i=n-1; i>=1; i--)
					Pi[i-1] = 1. / ( 2. * i / z + Pi[i] );
				bi[0] = ( 1. + 0.125 / z + 9. / 128. / z / z + 25. / 1024. / z / z / z )
					/ sqrt( pi2 * z );
				for(i=1; i<=6; i++)
					bi[i] = bi[i-1] * Pi[i-1];
			}
			else
			{
				if( z > 100. )
				{
					for( i=0; i<=6; i++ )
						bi[i] = 1. * exp( - i * i / 2. / z ) / sqrt( pi2 * z );
				}
			}
		}
	}
}

//在测试中由于输入的DgfAa22的值可能存在问题，需要进一步分析 wkliu
void FuncGdrsh(double dt, double bl0, double* s, double* sb, double* ss)
{
	double sink[1000], cosk[1000];
	int status, k;
	double rad_eps = 23.44 * Rad;
	double bl = bl0 + 0.1606413502e-3 * dt / 2.;
	double axh0 = cos( s[2] - pi / 6. ) * cos( bl ) + sin( s[2] - pi / 6. ) * sin( bl ) * cos(rad_eps);
	double bxh0 = - ( sin( s[2] - pi / 6. ) * cos(bl) - cos( s[2] - pi / 6. ) * sin(bl) * cos( rad_eps ) )
		* GsiDCi + sin(bl) * sin( rad_eps ) * GsiDSi;
	double E, F, dker;
	FuncKepl(sb[1], sb[5], E, F, status);
	for(k=1; k<=GdShoN; k++)
	{
		dker = k * E;
		RadianNorm(dker);
		sink[k] = sin(dker);
		cosk[k] = cos(dker);
	}
	double axh2 = -SatMoB1 * sb[0] * sb[0] * s[0] * ( SatMoBi[0] + SatMoGsbh * DgfEe * SatMoBi00
		          + ( s[3] * axh0 - s[4] * bxh0 ) * SatMoDfxh * SatMoBi20 / 4. );
	double aa = DgfAa22 / axh2 * sb[0];
	double das = 0.;
	double des = 0.;
	double dws = 0.;
	double dls = 0.;
	for( k=1; k<=GdShoN; k++)
	{
		//这里将存在问题，由于相应GdShoAk、GdShoEk、GdShoWk、GdShoKuxk存在问题
		das += ( GdShoAk[k] / k * sink[k] );
		des += ( GdShoEk[k] / k * sink[k] );
		dws += ( GdShoWk[k] * cosk[k] );
		dls += ( GdShoKuxk[k] / k * cosk[k] );
	}

	//则计算的ss数组存在误差
	ss[0] = -aa * sb[0] / ( 1. - GdShoBt4 )
      	* ( 0.5 * GdShoAk[0] * sb[1] * sink[1] + das);
	ss[1] = 0.;
	ss[2] = 0.;
	ss[3] = -aa / pow( (1. + GdShoBt2), 2. ) *
      	( s[3] / sb[1] * ( 1. - GdShoBt2 ) * ( 0.5 * GdShoEk[0] * sb[1] * sink[1] + des )
      	- s[4] / sb[1] * ( 1. + GdShoBt2 ) * ( 0.5 * GdShoWk[1] * sb[1] + dws ) );
	ss[4] = -aa / pow( ( 1. + GdShoBt2 ), 2. ) * 
      	( s[4] / sb[1] * ( 1. - GdShoBt2 ) * ( 0.5 * GdShoEk[0] * sb[1] * sink[1] + des )
      	+ s[3] / sb[1] * ( 1. + GdShoBt2 ) * ( 0.5 * GdShoWk[1] * sb[1] + dws ) );
	ss[5] = -aa / ( 1. - GdShoBt4 ) * ( 0.5 * GdShoKuxk[1] * sb[1] + dls );
}

//ZTSRTREF ok
void FuncTref(int k, double e, double x, double& y)
{
	double e1 = e;
	if( k == -1 )
		e1 = -e;
	y = ( cos(x) - e1 ) / ( 1. - e1 * cos(x) );
	y = acos(y);
	double x1 = sin(x);
	if( x1 < 0. )
		y = pi2 - y;
}

//ZTSRKEPL ok
void FuncKepl(double e, double m, double& x, double& y, int& status)
{
	int i = 0;
	double x1, f, f1;

	x = m;

	do
	{
		i ++;
		f = x - e * sin(x) - m;
		f1 = 1. - e * cos(x);
		x1 = f / f1;
		x = x - x1;
		if( i > 32 )
		{
			FuncTref(1, e, x, y);
			status = -1;
			return;
		}
		if( fabs(x1) <= 1.e-12 )
		{
			FuncTref(1, e, x, y);
			status = 0;
			return;
		}
	}while(1);
}

// ok
int anint(double val)
{	
	if( val>0. )
	    return int( val + 0.5 );
	else
	    return int( val - 0.5 );
}

// 输出较大的值
// ok 
double dmax1(double x,double y)
{
    return ( x > y ) ? x : y ;
}

// *** wkliu 由于PV_Keplerian函数的输入输出与我们平时使用的数据的单位有区别，为了方便，这里提供一个二次接口
// 在本函数中flag的取值与PV_Keplerian相同，但flag为1时为由位置速度转开普勒参数
// 其它情况下由开普勒参数转卫星位置速度
// par1为输入的卫星参数，par2为输出结果
// 输入数据和输出结果的单位为：卫星位置速度为m和m/s，卫星开普勒参数中长度单位为m，角度单位为角度
void PV_Keplerian1(int flag, double* par1, double* par2, int& status)
{
	double TempP1[6], TempP2[6];
	int i;

	// 输入参数为位置速度
	if( flag == 1 )
	{
		for(i=0; i<3; i++)
			TempP1[i] = par1[i] / AE;  //位置单位转化 AE
		for(i=3; i<6; i++)
			TempP1[i] = par1[i] * CKV; // 速度单位转化
	}
	else  // 输入参数为开普勒参数
	{
		for(i=2; i<6; i++)
			TempP1[i] = par1[i] * Rad;
		TempP1[0] = par1[0] / AE;
		TempP1[1] = par1[1];
	}

	PV_Keplerian(flag, TempP1, TempP2, status); //调用函数进行计算

	//输出为开普勒参数
	if( flag == 1 )
	{
		par2[0] = TempP2[0] * AE;
		par2[1] = TempP2[1];
		for(i=2; i<6; i++)
			par2[i] = TempP2[i] / Rad;
	}
	else  //输出为位置速度
	{
		for(i=0; i<3; i++)
			par2[i] = TempP2[i] * AE;
		for(i=3; i<6; i++)
			par2[i] = TempP2[i] / CKV;
	}

}

// 本函数只是一个接口函数，其功能为输入相关的必须参数，包括时间(time, sec)以及位置速度，
// 解算以后将其结果返回，结果包括提供的参考结果文件中所有的相关参数，以后可以对接口根据需要进行修改
// utc: utc约化儒略日
// xin：解算得到的开普勒瞬根
// xin：返回值，开普勒瞬根
// xout：开普勒平根
// xit：第一类无奇点瞬根
// xot：第一类无奇点平根
// da：长半轴变化率
// c1、c2、c3
// 由于程序中存在判断，我们增加一个返回值，当最后通过判断，得到正确的输出结果时，返回值为
// 1，否则返回0
int OrbParaCon(double utc, double* xyz, double* xin, double* xout, double* xit, 
				double* xot, double& da, double& c1, double& c2, double& c3)
{
	int status;

	PV_Keplerian1(1, xyz, xin, status);  //*** wkliu  由于对接口进行了修改，这里首先该卫星的开普勒瞬根，程序其它部分不变

/*
	// *** wkliu  测试，检查逆向转化得到的结果是否相同
	xyz[0] = 0;
	PV_Keplerian1(2, xin, xyz, status);
*/

	int time[5]; 
	double sec;

	// 根据该时刻的utc儒略日（utc），计算得到对应的格林高历，
	// time中为年月日时分，sec为日内秒
	CalDaySec(utc, time, sec);

	double temp1[6], temp2[6];

	// 下述函数为进行计算的实际函数，在本函数中我们计算出所必须的全部相关参数信息
	// 开普勒瞬根xout在本程序中计算得到，另外两组参数结果在下面的代码中计算得到
	OsculToMean(time, sec, xin, xout, temp1, temp2);

	//***  这里是仿造提供的代码进行修改得到，其中的utcf在提供的代码中未说明该值
	//***  是如何得到的，必须进行求证
	//***  目前我们进行测试的时候未进行该判断，结果已经满足要求，但为了保证可信度
	//***  必须进行分析，这里将这部分代码仍然仿造参考代码，进行此判断
//	if( fabs( utc - utcf ) < 1.e-9 )   //由于暂时没有变量utcf的值，所以暂时屏蔽到这一行， liuwk 20090112
	if( 1 )
	{
		////////////////////////////////////////////////////////////
		//进行结果的保存，可以看到第一类无奇点平根和瞬根的结果
		//只存在第一个长半轴的结果保存的对应开普勒参数中的长半轴的值
		//其他的结果就是计算得到的temp1和temp2中的结果
		xit[0] = xin[0];
		xit[1] = temp1[1];
		xit[2] = temp1[2];
		xit[3] = temp1[3];
		xit[4] = temp1[4];
		xit[5] = temp1[5];
		xot[0] = xout[0];
		xot[1] = temp2[1];
		xot[2] = temp2[2];
		xot[3] = temp2[3];
		xot[4] = temp2[4];
		xot[5] = temp2[5];

		//结果的保存中我们未参考原始代码将保存所有参数，而是只截取了
		//相关的三组参数
		/////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////
		// 下面是计算da、c1、c2、c3的部分
		double hight = ( xout[0] * ( 1. - xout[1] ) - 6378140. ) / 1000.;
		da = -1500. / hight / 86400.;

		c1 = 1.5 * 0.00108263 * pow( 6378137., 2. ) / pow( xout[0], 2. );
		c2 = pow( sin( xout[2] * Rad ), 2. );
		c3 = sqrt( 398600.436e9 / pow( xout[0], 3. ) );
		////////////////////////////////////////////////////////////////////////

		return 1;
	}

	return 0;
}

void CalDaySec(double utc, int* ITime, double& stime)
{
	double utcp = utc + 1.e-16 + 2400000.5;
	int jd = int( utcp );
	double t = fmod( ( utcp - jd ) * 24. + 12., 24. );

	ITime[3] = int ( t );
	t = ( t - ITime[3] ) * 60.;
	ITime[4] = int( t );
	stime = ( t - ITime[4] ) * 60.;

	if( ITime[3] < 12 )
		jd ++;
	int lx = jd + 68569;
	int nx = 4 * lx / 146097;
	lx = lx - ( 146097 * nx + 3 ) / 4;

	ITime[0] = 4000 * ( lx + 1 ) / 1461001;
	lx = lx - 1461 * ITime[0] / 4 + 31;
	ITime[1] = 80 * lx / 2447;
	ITime[2] = lx - 2447 * ITime[1] / 80;
	lx = ITime[1] / 11;

	ITime[1] = ITime[1] + 2 - 12 * lx;
	ITime[0] = 100 * ( nx - 49 ) + ITime[0] + lx;
}

/*
// 本函数与函数OrbParaCon的区别在于本函数是为了将一个时刻卫星的所有相关参数
// 全部保存到一个PODPara类型的结构体中，与提供的参考文件中的每条数据中所包含的
// 全部数据信息相同
void OrbParaConAll(double utc, double* xin)
{
	PODPara para;

	double xout[6], xit[6], xot[6];
	double da, c1, c2, c3;

	if( OrbParaCon(utc, xin, xout, xit, xot, da, c1, c2, c3) == 1 )
	{
	}
}
*/

// 这是计算da的程序入口，其方法为输入一组数据及对应时刻的a值，然后通过数学拟合计算得到
// num：已知的元素对的数目
// tim：数组中各点的时间数组
// xtp：数组中个点的长半轴数组
// da：计算得到的da的结果
// 这部分代码基本与提供的参考代码无区别，如果以后需要可进行测试
void NomialWsh(int num, double* tim, double* xtp, double& da)
{
	double dt;
	double atpl[2], x[2], px[4], an[4];

	int n = 2 * num;
	int n2 = n * n;

	double* at  = (double*)malloc( n * sizeof(double) );
	double* atp = (double*)malloc( n * sizeof(double));
	double* pw  = (double*)malloc( n2 * sizeof(double));
	double* ea  = (double*)malloc( n * sizeof(double));
	double* el  = (double*)malloc( num * sizeof(double));

	int isw = 1, ii;
	for(int i=0; i<num; i++ )
	{
		dt = ( tim[i] - tim[0] ) * 86400.;
		ii = i * 2 - 1;
		ea[ ii + 1 ] = 1.;
		ea[ ii + 2 ] = dt;
		el[i] = xtp[i];
	}

	int j = num;
	int j0 = 0, m;

	for(int k=0; k<j; k++)
	{
		for(int l=0; l<j; l++)
		{
			m = k * num + l;
			pw[m] = 0.;
			if( k == l )
				pw[m] = 1.;
		}
	}

	VTVMul(2, j, ea, el, pw, x, px, atp, an, atpl, at, isw);

	da = x[1] * 86400.;

	free(at);
	free(atp);
	free(pw);
	free(ea);
	free(el);
}

int VTVMul(int n, int m, double* a, double* el, double* p, double* x, double* px, 
			double* atp, double* an, double* atpl, double* at, int& isw)
{
	int ii, jj;
	for (ii=0; ii<n; ii++)
	{
		for (jj=0; jj<m; jj++)
		{
			at[ii * m + jj] = a[jj * n + ii];
		}
	}

	FuncAbMul(a, p, atp, n, m, m);
	FuncAbMul(atp, at, an, n, m, n);
	FuncAbMul(atp, el, atpl, n, m, 1);
	FuncLdlt(n, 1, an, px, atpl, x, isw);

	if ( isw == 0 )
	{
		printf("Singular Equation");
		return 0;
	}
	return 1;
}

void FuncAbMul(double* Arra, double* Arrb, double* Arrr, int l, int m, int n)
{
	int i, j, k;
	for(i=0; i<l; i++)
	{
		for(j=0; j<n; j++)
		{
			Arrr[j * l + i] = 0.;

			for(k=0;k<m;k++)
				Arrr[j * l + i] = Arrr[j * l + i] + Arra[k * l + i] * Arrb[j * m + k];
		}
	}
}

void FuncLdlt(int n, int m, double* Arra, double* Arrp, 
		  double* Arrb, double* Arrx, int& isw)
{
	int i, j, jm1, k, l;
	double w, y, z;

	for(i=0; i<n; i++)
	{
		for(j=0; j<i+1; j++)
		{
			w = Arra[i * n + j];

			if( i == j )
			{
				if( i != 0 )
				{
					for ( k=0; k<i; k++ )
					{
						y = Arra[ k * n + i ];
						z = y * Arrp[k];
						Arra[ k * n + i ] = z;
						w = w - y * z;
					}
				}

				if ( w != 0. ) 
				{
					Arrp[i] = 1. / w;
					continue;
				}
				isw = 0;
				return;
			}

			if ( j == 0 ) 
			{
				Arra[ j * n + i ] = w;
				continue;
			}

			jm1 = j - 1;

			for ( k=0; k<jm1; k++ )
			{
				w = w - Arra[ k * n + i ] * Arra[ k * n + j ];
			}
		}
	}

	for( j=0; j<m; j++ )
	{
		for( i=0; i<n; i++ )
		{
			y = Arrb[ j * n + i ];

			if( i == 0 )
			{
				Arrx[ j * n + i ] = y;
				continue;
			}

			for ( k=0; k<i; k++ )
			{
				y = y - Arra[ k * n + i ] * Arrx[ j * n + k ];
			}
			Arrx[ j * n + i ] = y;
		}
		for ( l=0; l<n; l++ )
		{
			i = n - l - 1;
			y = Arrx[ j * n + i ] * Arrp[i];

			if ( i == n-1 ) 
			{
				Arrx[ j * n + i ] = y;
				continue;
			}

			for( k=i; k<n; k++ )
			{
				y = y - Arra[i * n + k] * Arrx[j * n + k];
			}
			Arrx[ j * n + i ] = y;
		}
	}

	isw = 1;
}