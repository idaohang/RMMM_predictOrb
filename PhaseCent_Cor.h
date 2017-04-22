#pragma once
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/* coordinate rotation matrix ------------------------------------------------*/
#define Rx(t,X) do { \
	(X)[0]=1.0; (X)[1]=(X)[2]=(X)[3]=(X)[6]=0.0; \
	(X)[4]=(X)[8]=cos(t); (X)[7]=sin(t); (X)[5]=-(X)[7]; \
} while (0)

#define Ry(t,X) do { \
	(X)[4]=1.0; (X)[1]=(X)[3]=(X)[5]=(X)[7]=0.0; \
	(X)[0]=(X)[8]=cos(t); (X)[2]=sin(t); (X)[6]=-(X)[2]; \
} while (0)

#define Rz(t,X) do { \
	(X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0; \
	(X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)
#define Au          149597870691.0  /* 1 AU (m) */
#define AS2R        (D2R/3600.0)    /* arc sec to radian */
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define Pi          3.1415926535897932  /* pi */
#define D2R         (Pi/180.0)          /* deg to rad */
#define R2D         (180.0/Pi)          /* rad to deg */
#define MAXANT      64                  /* max length of station name/antenna type */
#define NFREQ       3
#define SYS_GPS     0x01                /* navigation system: gps */
#define SYS_SBS     0x02                /* navigation system: sbas */
#define SYS_GLO     0x04                /* navigation system: glonass */
#define SYS_NONE    0x00                /* navigation system: none */
#define MAXSAT      32
#define MAXPRNGPS   32                  /* max satellite prn number of GPS */
#define MAXPRNGLO   24                  /* max satellite slot number of GLONASS */


const static double leaps[][7]={ /* leap seconds {y,m,d,h,m,s,utc-gpst,...} */
	{2012,7,1,0,0,0,-16},
	{2009,1,1,0,0,0,-15},
	{2006,1,1,0,0,0,-14},
	{1999,1,1,0,0,0,-13},
	{1997,7,1,0,0,0,-12},
	{1996,1,1,0,0,0,-11},
	{1994,7,1,0,0,0,-10},
	{1993,7,1,0,0,0, -9},
	{1992,7,1,0,0,0, -8},
	{1991,1,1,0,0,0, -7},
	{1990,1,1,0,0,0, -6},
	{1988,1,1,0,0,0, -5},
	{1985,7,1,0,0,0, -4},
	{1983,7,1,0,0,0, -3},
	{1982,7,1,0,0,0, -2},
	{1981,7,1,0,0,0, -1}
};


typedef struct gtime_t{     /* time struct */
	time_t time;			/* time (s) expressed by standard time_t */
	double sec;				/* fraction of second under 1 s */
}gtime_t;
typedef struct pcv_t{			/* antenna parameter type */
	int sat;					/* satellite number (0:receiver) */
	char type[MAXANT];			/* antenna type */
	char code[MAXANT];			/* serial number or satellite code */
	gtime_t ts;
	gtime_t te;				/* valid time start and end */
	double off[NFREQ][ 3];		/* phase center offset e/n/u or x/y/z (m) */

	//Panda 2011-04-02
	double dazi;				//Increment of the azimuth£º0 to 360 with increment 

	//'DAZI'(in degrees).
	double zen1,zen2,dzen;      //Receiver antenna:Definition of the grid in zenith angle.
	//Satellite antenna:Definition of the grid in nadir angle.

	//double var[NFREQ][80*20];			/* phase center variation (m) */
} pcv_t;

typedef struct pcvs_t{        /* antenna parameters type */
	int n,nmax;         /* number of data/allocated */
	pcv_t *pcv;         /* antenna parameters data */
}pcvs_t;
typedef struct erp_t{   /* earth rotation parameter type */
	double xp,yp;       /* pole offset (rad) */
	double ut1_utc;     /* ut1-utc (s) */
	double ddeps,ddpsi; /* corrections to iau 1980 nutation */
}erp_t;



extern void    satantoff(gtime_t time, const double *rs, const pcv_t *pcv,
					  double *dant);
extern int     readantex(const char *file, pcvs_t *pcvs);
extern void    setpcv(gtime_t time, pcv_t *pcvst, const pcvs_t *pcvs);
extern int     normv3(const double *a, double *b);
extern double  norm(const double *a, int n);
extern double  dot(const double *a, const double *b, int n);
extern void    cross3(const double *a, const double *b, double *c);
extern void    sunmoonpos(gtime_t tutc, const erp_t *erp, double *rsun,
					   double *rmoon, double *gmst);
static void    sunmoonpos_eci(gtime_t tut, double *rsun, double *rmoon);
static void    ast_args(double t, double *f);
extern gtime_t epoch2time(const double *ep);
extern void    eci2ecef(gtime_t tutc, const erp_t *erp, double *U, double *gmst);
extern void    matmul(const char *tr, int n, int k, int m, double alpha,
				   const double *A, const double *B, double beta, double *C);
static void    addpcv(const pcv_t *pcv, pcvs_t *pcvs);
static int     decodef(char *p, int n, double *v);
gtime_t        timeadd(gtime_t t, double sec);
double         timediff(gtime_t t1, gtime_t t2);
extern int     round(double dNum);
static double  time2sec(gtime_t time, gtime_t *day);
gtime_t        gpst2utc(gtime_t t);
gtime_t        utc2gpst(gtime_t t);
double         str2num(const char *s, int i, int n);
extern int     satno(int sys, int prn);
int            str2time(const char *s, int i, int n, gtime_t *t);
void           time2str(gtime_t t, char *s, int n);
extern void    time2epoch(gtime_t t, double *ep);
static void    nut_iau1980(double t, const double *f, double *dpsi, double *deps);
extern pcv_t  *searchpcv(int sat, const char *type, gtime_t time,
						const pcvs_t *pcvs);
extern int     satsys(int sat, int *prn);
extern void    satno2id(int sat, char *id);