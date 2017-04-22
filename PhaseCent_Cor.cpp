#include "PhaseCent_Cor.h"


/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          pcv_t  *pcv        I   satellite antenna parameter
*          double *dant       IO  satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m)
* return : none
*-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, const double *rs, const pcv_t *pcv,
                      double *dant)
{
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst;
    int i;

    /* sun position in ecef */
    sunmoonpos(gpst2utc(time),NULL,rsun,NULL,&gmst);

    /* unit vectors of satellite fixed coordinates */
    for (i=0; i<3; i++) r[i]=-rs[i];
    if (!normv3(r,ez)) return;
    for (i=0; i<3; i++) r[i]=rsun[i]-rs[i];
    if (!normv3(r,es)) return;
    cross3(ez,es,r);
    if (!normv3(r,ey)) return;
    cross3(ey,ez,ex);

    for (i=0; i<3; i++) /* use L1 value */
    {
        dant[i]=pcv->off[0][0]*ex[i]+pcv->off[0][1]*ey[i]+pcv->off[0][2]*ez[i];

    }
}
/* read antex file ----------------------------------------------------------*/
extern int readantex(const char *file, pcvs_t *pcvs)
{
    FILE *fp;
    static const pcv_t pcv0= {0};
    pcv_t pcv;
    double neu[3],dd;
    int f,prn,freq=0;
    char buff[256];

    if (!(fp=fopen(file,"r")))
    {
        return 0;
    }
    
    while (fgets(buff,sizeof(buff),fp))
    {
        if (strlen(buff)<60||strstr(buff+60,"COMMENT")) continue;

        if (strstr(buff+60,"START OF ANTENNA"))
        {
            pcv=pcv0;
        }
        if (strstr(buff+60,"END OF ANTENNA"))
        {
            addpcv(&pcv,pcvs);
            continue;
        }

        if (strstr(buff+60,"TYPE / SERIAL NO"))
        {
            strncpy(pcv.type,buff   ,20);
            pcv.type[20]='\0';
            strncpy(pcv.code,buff+20,20);
            pcv.code[20]='\0';
            if (!(prn=(int)str2num(pcv.code,1,2))) continue;
            switch (pcv.code[0])
            {
            case 'G':
                pcv.sat=satno(SYS_GPS,prn);
                break;
            case 'R':
                pcv.sat=satno(SYS_GLO,prn);
                break;
            }
        }
        else if (strstr(buff+60,"VALID FROM"))
        {
            if (!str2time(buff,0,43,&pcv.ts)) continue;
        }
        else if (strstr(buff+60,"VALID UNTIL"))
        {
            if (!str2time(buff,0,43,&pcv.te)) continue;
        }
        else if (strstr(buff+60,"DAZI"))
        {
            pcv.dazi=str2num(buff,2,6);
            continue;
        }
        else if (strstr(buff+60,"ZEN1 / ZEN2 / DZEN"))
        {
            pcv.zen1=str2num(buff,2,6);
            pcv.zen2=str2num(buff,8,6);
            pcv.dzen=str2num(buff,14,6);
            continue;
        }
        else if (strstr(buff+60,"START OF FREQUENCY"))
        {
            if (sscanf(buff+4,"%d",&f)<1) continue;

            freq=f==1?1:(f==2?2:(f==5?3:0)); /* L1/L2/L5 */
        }
        else if (strstr(buff+60,"END OF FREQUENCY"))
        {
            freq=0;
        }
        else if (strstr(buff+60,"NORTH / EAST / UP"))
        {
            if (freq<1||NFREQ<freq) continue;
            if (decodef(buff,3,neu)<3) continue;
            pcv.off[freq-1][0]=neu[pcv.sat?0:1]; /* x or e */
            pcv.off[freq-1][1]=neu[pcv.sat?1:0]; /* y or n */
            pcv.off[freq-1][2]=neu[2];           /* z or u */
        }

        //panda??    2010-12-01
        //接收机PCV考虑随方位角的变化
        else if (strstr(buff,"NOAZI"))
        {

            if (freq<1||NFREQ<freq) continue;

            dd=(pcv.zen2-pcv.zen1)/pcv.dzen+1;

            if (dd!=round(dd)||dd<=1)
            {
                printf("zen in atx file error (d!=round(d)||d<1)!");
                continue;
            }

            //if (pcv.dazi==0.0)
            //{
            //    //pcv.var[freq-1]=new double[int(dd)];
            //    i=decodef(buff+8,(int)dd,pcv.var[freq-1]);

            //    if (i<=0)
            //    {
            //        printf("error in reading atx (i<=0)!");
            //        continue;
            //    }
            //    else if (i!=(int)dd)
            //    {
            //        printf("error in reading atx (i!=(int)dd)!");
            //        continue;
            //    }
            //}
            //else
            //{
            //    id=(int)((360-0)/pcv.dazi)+1;
            //    //pcv.var[freq-1]=new double[int(dd)*id];

            //    for (i=0; i<id; i++)
            //    {
            //        fgets(buff,sizeof(buff),fp);

            //        j=decodef(buff+8,(int)dd,&pcv.var[freq-1][i*(int)dd]);

            //        if (j<=0)
            //        {
            //            printf("error in reading atx (j<=0)!");
            //            continue;
            //        }
            //        else if (j!=(int)dd)
            //        {
            //            printf("error in reading atx (j!=(int)dd)!");
            //            continue;
            //        }
            //    }
            //}
        }
    }
    fclose(fp);

    return 1;
}

/* set antenna parameters ----------------------------------------------------*/
extern void setpcv(gtime_t time, pcv_t *pcvst, const pcvs_t *pcvs)
{
    pcv_t *pcv;
    int i;
    char id[64];
    //popt->navsys==SYS_GPS;

    /* set satellite antenna parameters */
    for (i=0; i<MAXSAT; i++)
    {
        if (!(satsys(i+1,NULL)&SYS_GPS)) continue;
        if (!(pcv=searchpcv(i+1,"",time,pcvs)))
        {
            satno2id(i+1,id);
            continue;
        }
        pcvst[i]=*pcv;
    }
}
/* normalize 3d vector ---------------------------------------------------------
* normalize 3d vector
* args   : double *a        I   vector a (3 x 1)
*          double *b        O   normlized vector (3 x 1) || b || = 1
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int normv3(const double *a, double *b)
{
    double r;
    if ((r=norm(a,3))<=0.0) return 0;
    b[0]=a[0]/r;
    b[1]=a[1]/r;
    b[2]=a[2]/r;
    return 1;
}
/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern double norm(const double *a, int n)
{
    return sqrt(dot(a,a,n));
}
/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n)
{
    double c=0.0;

    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* outer product of 3d vectors -------------------------------------------------
* outer product of 3d vectors
* args   : double *a,*b     I   vector a,b (3 x 1)
*          double *c        O   outer product (a x b) (3 x 1)
* return : none
*-----------------------------------------------------------------------------*/
extern void cross3(const double *a, const double *b, double *c)
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}
/* sun and moon position -------------------------------------------------------
* get sun and moon position in ecef
* args   : gtime_t tut      I   time in ut1
*          erp_t  *erp      I   earth rotation parameters (NULL: not used)
*          double *rsun     O   sun position in ecef (m)
*          double *rmoon    O   moon position in ecef (m)
*          double *gmst     O   gmst (rad)
* return : none
*-----------------------------------------------------------------------------*/
extern void sunmoonpos(gtime_t tutc, const erp_t *erp, double *rsun,
                       double *rmoon, double *gmst)
{
    gtime_t tut;
    double rs[3],rm[3],U[9];

    tut=erp?timeadd(tutc,erp->ut1_utc):tutc; /* utc -> ut1 */

    /* sun and moon position in eci */
    sunmoonpos_eci(tut,rsun?rs:NULL,rmoon?rm:NULL);

    /* eci to ecef transformation matrix */
    eci2ecef(tutc,erp,U,gmst);

    /* sun and moon postion in ecef */
    if (rsun ) matmul("NN",3,1,3,1.0,U,rs,0.0,rsun );
    if (rmoon) matmul("NN",3,1,3,1.0,U,rm,0.0,rmoon);
}
/* sun and moon position in eci (ref [4] 5.1.1, 5.2.1) -----------------------*/
static void sunmoonpos_eci(gtime_t tut, double *rsun, double *rmoon)
{
    const double ep2000[]= {2000,1,1,12,0,0};
    double t,f[5],eps,Ms,ls,rs,lm,pm,rm,sine,cose,sinp,cosp,sinl,cosl;

    t=timediff(tut,epoch2time(ep2000))/86400.0/36525.0;

    /* astronomical arguments */
    ast_args(t,f);

    /* obliquity of the ecliptic */
    eps=23.439291-0.0130042*t;
    sine=sin(eps*D2R);
    cose=cos(eps*D2R);

    /* sun position in eci */
    if (rsun)
    {
        Ms=357.5277233+35999.05034*t;
        ls=280.460+36000.770*t+1.914666471*sin(Ms*D2R)+0.019994643*sin(2.0*Ms*D2R);
        rs=Au*(1.000140612-0.016708617*cos(Ms*D2R)-0.000139589*cos(2.0*Ms*D2R));
        sinl=sin(ls*D2R);
        cosl=cos(ls*D2R);
        rsun[0]=rs*cosl;
        rsun[1]=rs*cose*sinl;
        rsun[2]=rs*sine*sinl;
    }
    /* moon position in eci */
    if (rmoon)
    {
        lm=218.32+481267.883*t+6.29*sin(f[0])-1.27*sin(f[0]-2.0*f[3])+
           0.66*sin(2.0*f[3])+0.21*sin(2.0*f[0])-0.19*sin(f[1])-0.11*sin(2.0*f[2]);
        pm=5.13*sin(f[2])+0.28*sin(f[0]+f[2])-0.28*sin(f[2]-f[0])-
           0.17*sin(f[2]-2.0*f[3]);
        rm=RE_WGS84/sin((0.9508+0.0518*cos(f[0])+0.0095*cos(f[0]-2.0*f[3])+
                         0.0078*cos(2.0*f[3])+0.0028*cos(2.0*f[0]))*D2R);
        sinl=sin(lm*D2R);
        cosl=cos(lm*D2R);
        sinp=sin(pm*D2R);
        cosp=cos(pm*D2R);
        rmoon[0]=rm*cosp*cosl;
        rmoon[1]=rm*(cose*cosp*sinl-sine*sinp);
        rmoon[2]=rm*(sine*cosp*sinl+cose*sinp);
    }
}
/* astronomical arguments: f={l,l',F,D,OMG} (rad) ----------------------------*/
static void ast_args(double t, double *f)
{
    static const double fc[][5]=  /* coefficients for iau 1980 nutation */
    {
        { 134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470},
        { 357.52910918,  129596581.0481,  -0.5532,  0.000136, -0.00001149},
        {  93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417},
        { 297.85019547, 1602961601.2090,  -6.3706,  0.006593, -0.00003169},
        { 125.04455501,   -6962890.2665,   7.4722,  0.007702  -0.00005939}
    };
    double tt[4];
    int i,j;

    for (tt[0]=t,i=1; i<4; i++) tt[i]=tt[i-1]*t;
    for (i=0; i<5; i++)
    {
        f[i]=fc[i][0]*3600.0;
        for (j=0; j<4; j++) f[i]+=fc[i][j+1]*tt[j];
        f[i]=fmod(f[i]*AS2R,2.0*Pi);
    }
}
/* convert calendar day/time to time -------------------------------------------
* convert calendar day/time to gtime_t struct
* args   : double *ep       I   day/time {year,month,day,hour,min,sec}
* return : gtime_t struct
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern gtime_t epoch2time(const double *ep)
{
    const int doy[]= {1,32,60,91,121,152,182,213,244,274,305,335};
    gtime_t time= {0};
    int days,sec,year=(int)ep[0],mon=(int)ep[1],day=(int)ep[2];

    if (year<1970||2099<year||mon<1||12<mon) return time;

    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3?1:0);
    sec=(int)floor(ep[5]);
    time.time=(time_t)days*86400+(int)ep[3]*3600+(int)ep[4]*60+sec;
    time.sec=ep[5]-sec;
    return time;
}

/* eci to ecef transformation matrix -------------------------------------------
* compute eci to ecef transformation matrix
* args   : gtime_t tutc     I   time in utc
*          erp_t *erp       I   earth rotation parameters (NULL: not used)
*          double *U        O   eci to ecef transformation matrix (3 x 3)
*          double *gmst     IO  greenwich mean sidereal time (rad)
*                               (NULL: no output)
* return : none
* note   : see ref [3] chap 5
*-----------------------------------------------------------------------------*/
extern void eci2ecef(gtime_t tutc, const erp_t *erp, double *U, double *gmst)
{
    const double ep2000[]= {2000,1,1,12,0,0};
    static gtime_t tutc_;
    static double U_[9],gmst_;
    gtime_t tgps,tut,tut0;
    erp_t erp_= {0};
    double t,t2,t3,eps,ze,th,z,dpsi,deps,ut,gmst0,gast,f[5];
    double R1[9],R2[9],R3[9],R[9],W[9],N[9],P[9],NP[9];
    int i;

    if (fabs(timediff(tutc,tutc_))<0.01)   /* read cache */
    {
        for (i=0; i<9; i++) U[i]=U_[i];
        if (gmst) *gmst=gmst_;
        return;
    }
    tutc_=tutc;
    if (erp) erp_=*erp;

    /* terrestrial time */
    tgps=utc2gpst(tutc_);
    t=(timediff(tgps,epoch2time(ep2000))+19.0+32.184)/86400.0/36525.0;
    t2=t*t;
    t3=t2*t;

    /* astronomical arguments */
    ast_args(t,f);

    /* iau 1976 precession */
    ze=(2306.2181*t+0.30188*t2+0.017998*t3)*AS2R;
    th=(2004.3109*t-0.42665*t2-0.041833*t3)*AS2R;
    z =(2306.2181*t+1.09468*t2+0.018203*t3)*AS2R;
    eps=(84381.448-46.8150*t-0.00059*t2+0.001813*t3)*AS2R;
    Rz(-z,R1);
    Ry(th,R2);
    Rz(-ze,R3);
    matmul("NN",3,3,3,1.0,R1,R2,0.0,R);
    matmul("NN",3,3,3,1.0,R, R3,0.0,P); /* P=Rz(-z)*Ry(th)*Rz(-ze) */

    /* iau 1980 nutation */
    nut_iau1980(t,f,&dpsi,&deps);
    Rx(-eps-deps-erp_.ddeps,R1);
    Rz(-dpsi-erp_.ddpsi,R2);
    Rx(eps,R3);
    matmul("NN",3,3,3,1.0,R1,R2,0.0,R);
    matmul("NN",3,3,3,1.0,R ,R3,0.0,N); /* N=Rx(-eps)*Rz(-dspi)*Rx(eps) */

    /* greenwich mean/aparent sidereal time (rad) */
    tut=timeadd(tutc_,erp_.ut1_utc);
    ut=time2sec(tut,&tut0);
    t=timediff(tut0,epoch2time(ep2000))/86400.0/36525.0;
    t2=t*t;
    t3=t2*t;
    gmst0=24110.54841+8640184.812866*t+0.093104*t2-6.2E-6*t3;
    gmst_=gmst0+1.002737909350795*ut;
    gmst_=fmod(gmst_,86400.0)*Pi/43200.0;
    gast=gmst_+dpsi*cos(eps);
    gast+=(0.00264*sin(f[4])+0.000063*sin(2.0*f[4]))*AS2R;

    /* eci to ecef transformation matrix */
    Ry(-erp_.xp,R1);
    Rx(-erp_.yp,R2);
    Rz(gast,R3);
    matmul("NN",3,3,3,1.0,R1,R2,0.0,W );
    matmul("NN",3,3,3,1.0,W ,R3,0.0,R ); /* W=Ry(-xp)*Rx(-yp) */
    matmul("NN",3,3,3,1.0,N ,P ,0.0,NP);
    matmul("NN",3,3,3,1.0,R ,NP,0.0,U_); /* U=W*Rz(gast)*N*P */

    for (i=0; i<9; i++) U[i]=U_[i];
    if (gmst) *gmst=gmst_;

}


extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    double d;
    int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

    for (i=0; i<n; i++) for (j=0; j<k; j++)
        {
            d=0.0;
            switch (f)
            {
            case 1:
                for (x=0; x<m; x++) d+=A[i+x*n]*B[x+j*m];
                break;
            case 2:
                for (x=0; x<m; x++) d+=A[i+x*n]*B[j+x*k];
                break;
            case 3:
                for (x=0; x<m; x++) d+=A[x+i*m]*B[x+j*m];
                break;
            case 4:
                for (x=0; x<m; x++) d+=A[x+i*m]*B[j+x*k];
                break;
            }
            if (beta==0.0) C[i+j*n]=alpha*d;
            else C[i+j*n]=alpha*d+beta*C[i+j*n];
        }
}

/* add antenna parameter -----------------------------------------------------*/
static void addpcv(const pcv_t *pcv, pcvs_t *pcvs)
{
    if (pcvs->nmax<=pcvs->n)
    {
        pcvs->nmax+=256;
        if (!(pcvs->pcv=(pcv_t *)realloc(pcvs->pcv,sizeof(pcv_t)*pcvs->nmax)))
        {
            pcvs->n=pcvs->nmax=0;
            return;
        }
    }
    pcvs->pcv[pcvs->n++]=*pcv;
}
/* decode antenna parameter field --------------------------------------------*/
static int decodef(char *p, int n, double *v)
{
    int i;

    for (i=0; i<n; i++) v[i]=0.0;
    for (i=0,p=strtok(p," "); p&&i<n; p=strtok(NULL," "))
    {
        v[i++]=atof(p)*1E-3;
    }
    return i;
}
/* add time --------------------------------------------------------------------
* add time to gtime_t struct
* args   : gtime_t t        I   gtime_t struct
*          double sec       I   time to add (s)
* return : gtime_t struct (t+sec)
*-----------------------------------------------------------------------------*/
gtime_t timeadd(gtime_t t, double sec)
{
    double tt;

    t.sec+=sec;
    tt=floor(t.sec);
    t.time+=(int)tt;
    t.sec-=tt;
    return t;
}


/* time difference -------------------------------------------------------------
* difference between gtime_t structs
* args   : gtime_t t1,t2    I   gtime_t structs
* return : time difference (t1-t2) (s)
*-----------------------------------------------------------------------------*/
double timediff(gtime_t t1, gtime_t t2)
{
    return difftime(t1.time,t2.time)+t1.sec-t2.sec;
}

extern int round(double dNum)
{
    int iNum;
    if(dNum >= 0)
        iNum = (int)( dNum + 0.5 );
    else
        iNum = (int)( dNum - 0.5 );

    return iNum;
}
/* time to day and sec -------------------------------------------------------*/
static double time2sec(gtime_t time, gtime_t *day)
{
    double ep[6],sec;
    time2epoch(time,ep);
    sec=ep[3]*3600.0+ep[4]*60.0+ep[5];
    ep[3]=ep[4]=ep[5]=0.0;
    *day=epoch2time(ep);
    return sec;
}
/* gpstime to utc --------------------------------------------------------------
* convert gpstime to utc considering leap seconds
* args   : gtime_t t        I   time expressed in gpstime
* return : time expressed in utc
* notes  : ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
gtime_t gpst2utc(gtime_t t)
{
    gtime_t tu;
    int i;

    for (i=0; i<sizeof(leaps)/sizeof(*leaps); i++)
    {
        tu=timeadd(t,leaps[i][6]);
        if (timediff(tu,epoch2time(leaps[i]))>=0.0) return tu;
    }
    return t;
}


/* utc to gpstime --------------------------------------------------------------
* convert utc to gpstime considering leap seconds
* args   : gtime_t t        I   time expressed in utc
* return : time expressed in gpstime
* notes  : ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
gtime_t utc2gpst(gtime_t t)
{
    int i;

    for (i=0; i<sizeof(leaps)/sizeof(*leaps); i++)
    {
        if (timediff(t,epoch2time(leaps[i]))>=0.0) return timeadd(t,-leaps[i][6]);
    }
    return t;
}

/* string to number ------------------------------------------------------------
* convert substring in string to number
* args   : char   *s        I   string ("... nnn.nnn ...")
*          int    i,n       I   substring position and width
* return : converted number (0.0:error)
*-----------------------------------------------------------------------------*/
double str2num(const char *s, int i, int n)
{
    double value;
    char str[1024],*p=str;

    if (i<0||(int)strlen(s)<i||sizeof(str)-1<i) return 0.0;
    for (s+=i; *s&&--n>=0; s++) *p++=*s=='d'||*s=='D'?'E':*s;
    *p='\0';
    return sscanf(str,"%lf",&value)==1?value:0.0;
}


/* satellite system+prn/slot number to satellite number ------------------------
* convert satellite system+prn/slot number to satellite number
* args   : int    sys       I   satellite system (SYS_GPS,SYS_GLO,...)
*          int    prn       I   satellite prn/slot number
* return : satellite number (0:error)
*-----------------------------------------------------------------------------*/
extern int satno(int sys, int prn)
{
    if (prn<=0) return 0;
    switch (sys)
    {
    case SYS_GPS:
        return prn>32?0:prn;
    }
    return 0;
}
/* string to time --------------------------------------------------------------
* convert substring in string to gtime_t struct
* args   : char   *s        I   string ("... yyyy mm dd hh mm ss ...")
*          int    i,n       I   substring position and width
*          gtime_t *t       O   gtime_t struct
* return : status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
int str2time(const char *s, int i, int n, gtime_t *t)
{
    double ep[6];
    char str[256],*p=str;

    if (i<0||(int)strlen(s)<i||sizeof(str)-1<i) return -1;
    for (s+=i; *s&&--n>=0;) *p++=*s++;
    *p='\0';
    if (sscanf(str,"%lf %lf %lf %lf %lf %lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<6)
        return -1;
    if (ep[0]<100.0) ep[0]+=ep[0]<80.0?2000.0:1900.0;
    *t=epoch2time(ep);
    return 0;
}
void time2str(gtime_t t, char *s, int n)
{
    double ep[6];

    if (n<0) n=0;
    else if (n>12) n=12;
    if (1.0-t.sec<0.5/pow(10.0,n))
    {
        t.time++;
        t.sec=0.0;
    };
    time2epoch(t,ep);
    sprintf(s,"%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%0*.*f",ep[0],ep[1],ep[2],
            ep[3],ep[4],n<=0?2:n+3,n<=0?0:n,ep[5]);
}
/* time to calendar day/time ---------------------------------------------------
* convert gtime_t struct to calendar day/time
* args   : gtime_t t        I   gtime_t struct
*          double *ep       O   day/time {year,month,day,hour,min,sec}
* return : none
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern void time2epoch(gtime_t t, double *ep)
{
    const int mday[]=  /* # of days in a month */
    {
        31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
        31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    int days,sec,mon,day;

    /* leap year if year%4==0 in 1901-2099 */
    days=(int)(t.time/86400);
    sec=(int)(t.time-(time_t)days*86400);
    for (day=days%1461,mon=0; mon<48; mon++)
    {
        if (day>=mday[mon]) day-=mday[mon];
        else break;
    }
    ep[0]=1970+days/1461*4+mon/12;
    ep[1]=mon%12+1;
    ep[2]=day+1;
    ep[3]=sec/3600;
    ep[4]=sec%3600/60;
    ep[5]=sec%60+t.sec;
}

/* iau 1980 nutation ---------------------------------------------------------*/
static void nut_iau1980(double t, const double *f, double *dpsi, double *deps)
{
    static const double nut[106][10]=
    {
        {   0,   0,   0,   0,   1, -6798.4, -171996, -174.2, 92025,   8.9},
        {   0,   0,   2,  -2,   2,   182.6,  -13187,   -1.6,  5736,  -3.1},
        {   0,   0,   2,   0,   2,    13.7,   -2274,   -0.2,   977,  -0.5},
        {   0,   0,   0,   0,   2, -3399.2,    2062,    0.2,  -895,   0.5},
        {   0,  -1,   0,   0,   0,  -365.3,   -1426,    3.4,    54,  -0.1},
        {   1,   0,   0,   0,   0,    27.6,     712,    0.1,    -7,   0.0},
        {   0,   1,   2,  -2,   2,   121.7,    -517,    1.2,   224,  -0.6},
        {   0,   0,   2,   0,   1,    13.6,    -386,   -0.4,   200,   0.0},
        {   1,   0,   2,   0,   2,     9.1,    -301,    0.0,   129,  -0.1},
        {   0,  -1,   2,  -2,   2,   365.2,     217,   -0.5,   -95,   0.3},
        {  -1,   0,   0,   2,   0,    31.8,     158,    0.0,    -1,   0.0},
        {   0,   0,   2,  -2,   1,   177.8,     129,    0.1,   -70,   0.0},
        {  -1,   0,   2,   0,   2,    27.1,     123,    0.0,   -53,   0.0},
        {   1,   0,   0,   0,   1,    27.7,      63,    0.1,   -33,   0.0},
        {   0,   0,   0,   2,   0,    14.8,      63,    0.0,    -2,   0.0},
        {  -1,   0,   2,   2,   2,     9.6,     -59,    0.0,    26,   0.0},
        {  -1,   0,   0,   0,   1,   -27.4,     -58,   -0.1,    32,   0.0},
        {   1,   0,   2,   0,   1,     9.1,     -51,    0.0,    27,   0.0},
        {  -2,   0,   0,   2,   0,  -205.9,     -48,    0.0,     1,   0.0},
        {  -2,   0,   2,   0,   1,  1305.5,      46,    0.0,   -24,   0.0},
        {   0,   0,   2,   2,   2,     7.1,     -38,    0.0,    16,   0.0},
        {   2,   0,   2,   0,   2,     6.9,     -31,    0.0,    13,   0.0},
        {   2,   0,   0,   0,   0,    13.8,      29,    0.0,    -1,   0.0},
        {   1,   0,   2,  -2,   2,    23.9,      29,    0.0,   -12,   0.0},
        {   0,   0,   2,   0,   0,    13.6,      26,    0.0,    -1,   0.0},
        {   0,   0,   2,  -2,   0,   173.3,     -22,    0.0,     0,   0.0},
        {  -1,   0,   2,   0,   1,    27.0,      21,    0.0,   -10,   0.0},
        {   0,   2,   0,   0,   0,   182.6,      17,   -0.1,     0,   0.0},
        {   0,   2,   2,  -2,   2,    91.3,     -16,    0.1,     7,   0.0},
        {  -1,   0,   0,   2,   1,    32.0,      16,    0.0,    -8,   0.0},
        {   0,   1,   0,   0,   1,   386.0,     -15,    0.0,     9,   0.0},
        {   1,   0,   0,  -2,   1,   -31.7,     -13,    0.0,     7,   0.0},
        {   0,  -1,   0,   0,   1,  -346.6,     -12,    0.0,     6,   0.0},
        {   2,   0,  -2,   0,   0, -1095.2,      11,    0.0,     0,   0.0},
        {  -1,   0,   2,   2,   1,     9.5,     -10,    0.0,     5,   0.0},
        {   1,   0,   2,   2,   2,     5.6,      -8,    0.0,     3,   0.0},
        {   0,  -1,   2,   0,   2,    14.2,      -7,    0.0,     3,   0.0},
        {   0,   0,   2,   2,   1,     7.1,      -7,    0.0,     3,   0.0},
        {   1,   1,   0,  -2,   0,   -34.8,      -7,    0.0,     0,   0.0},
        {   0,   1,   2,   0,   2,    13.2,       7,    0.0,    -3,   0.0},
        {  -2,   0,   0,   2,   1,  -199.8,      -6,    0.0,     3,   0.0},
        {   0,   0,   0,   2,   1,    14.8,      -6,    0.0,     3,   0.0},
        {   2,   0,   2,  -2,   2,    12.8,       6,    0.0,    -3,   0.0},
        {   1,   0,   0,   2,   0,     9.6,       6,    0.0,     0,   0.0},
        {   1,   0,   2,  -2,   1,    23.9,       6,    0.0,    -3,   0.0},
        {   0,   0,   0,  -2,   1,   -14.7,      -5,    0.0,     3,   0.0},
        {   0,  -1,   2,  -2,   1,   346.6,      -5,    0.0,     3,   0.0},
        {   2,   0,   2,   0,   1,     6.9,      -5,    0.0,     3,   0.0},
        {   1,  -1,   0,   0,   0,    29.8,       5,    0.0,     0,   0.0},
        {   1,   0,   0,  -1,   0,   411.8,      -4,    0.0,     0,   0.0},
        {   0,   0,   0,   1,   0,    29.5,      -4,    0.0,     0,   0.0},
        {   0,   1,   0,  -2,   0,   -15.4,      -4,    0.0,     0,   0.0},
        {   1,   0,  -2,   0,   0,   -26.9,       4,    0.0,     0,   0.0},
        {   2,   0,   0,  -2,   1,   212.3,       4,    0.0,    -2,   0.0},
        {   0,   1,   2,  -2,   1,   119.6,       4,    0.0,    -2,   0.0},
        {   1,   1,   0,   0,   0,    25.6,      -3,    0.0,     0,   0.0},
        {   1,  -1,   0,  -1,   0, -3232.9,      -3,    0.0,     0,   0.0},
        {  -1,  -1,   2,   2,   2,     9.8,      -3,    0.0,     1,   0.0},
        {   0,  -1,   2,   2,   2,     7.2,      -3,    0.0,     1,   0.0},
        {   1,  -1,   2,   0,   2,     9.4,      -3,    0.0,     1,   0.0},
        {   3,   0,   2,   0,   2,     5.5,      -3,    0.0,     1,   0.0},
        {  -2,   0,   2,   0,   2,  1615.7,      -3,    0.0,     1,   0.0},
        {   1,   0,   2,   0,   0,     9.1,       3,    0.0,     0,   0.0},
        {  -1,   0,   2,   4,   2,     5.8,      -2,    0.0,     1,   0.0},
        {   1,   0,   0,   0,   2,    27.8,      -2,    0.0,     1,   0.0},
        {  -1,   0,   2,  -2,   1,   -32.6,      -2,    0.0,     1,   0.0},
        {   0,  -2,   2,  -2,   1,  6786.3,      -2,    0.0,     1,   0.0},
        {  -2,   0,   0,   0,   1,   -13.7,      -2,    0.0,     1,   0.0},
        {   2,   0,   0,   0,   1,    13.8,       2,    0.0,    -1,   0.0},
        {   3,   0,   0,   0,   0,     9.2,       2,    0.0,     0,   0.0},
        {   1,   1,   2,   0,   2,     8.9,       2,    0.0,    -1,   0.0},
        {   0,   0,   2,   1,   2,     9.3,       2,    0.0,    -1,   0.0},
        {   1,   0,   0,   2,   1,     9.6,      -1,    0.0,     0,   0.0},
        {   1,   0,   2,   2,   1,     5.6,      -1,    0.0,     1,   0.0},
        {   1,   1,   0,  -2,   1,   -34.7,      -1,    0.0,     0,   0.0},
        {   0,   1,   0,   2,   0,    14.2,      -1,    0.0,     0,   0.0},
        {   0,   1,   2,  -2,   0,   117.5,      -1,    0.0,     0,   0.0},
        {   0,   1,  -2,   2,   0,  -329.8,      -1,    0.0,     0,   0.0},
        {   1,   0,  -2,   2,   0,    23.8,      -1,    0.0,     0,   0.0},
        {   1,   0,  -2,  -2,   0,    -9.5,      -1,    0.0,     0,   0.0},
        {   1,   0,   2,  -2,   0,    32.8,      -1,    0.0,     0,   0.0},
        {   1,   0,   0,  -4,   0,   -10.1,      -1,    0.0,     0,   0.0},
        {   2,   0,   0,  -4,   0,   -15.9,      -1,    0.0,     0,   0.0},
        {   0,   0,   2,   4,   2,     4.8,      -1,    0.0,     0,   0.0},
        {   0,   0,   2,  -1,   2,    25.4,      -1,    0.0,     0,   0.0},
        {  -2,   0,   2,   4,   2,     7.3,      -1,    0.0,     1,   0.0},
        {   2,   0,   2,   2,   2,     4.7,      -1,    0.0,     0,   0.0},
        {   0,  -1,   2,   0,   1,    14.2,      -1,    0.0,     0,   0.0},
        {   0,   0,  -2,   0,   1,   -13.6,      -1,    0.0,     0,   0.0},
        {   0,   0,   4,  -2,   2,    12.7,       1,    0.0,     0,   0.0},
        {   0,   1,   0,   0,   2,   409.2,       1,    0.0,     0,   0.0},
        {   1,   1,   2,  -2,   2,    22.5,       1,    0.0,    -1,   0.0},
        {   3,   0,   2,  -2,   2,     8.7,       1,    0.0,     0,   0.0},
        {  -2,   0,   2,   2,   2,    14.6,       1,    0.0,    -1,   0.0},
        {  -1,   0,   0,   0,   2,   -27.3,       1,    0.0,    -1,   0.0},
        {   0,   0,  -2,   2,   1,  -169.0,       1,    0.0,     0,   0.0},
        {   0,   1,   2,   0,   1,    13.1,       1,    0.0,     0,   0.0},
        {  -1,   0,   4,   0,   2,     9.1,       1,    0.0,     0,   0.0},
        {   2,   1,   0,  -2,   0,   131.7,       1,    0.0,     0,   0.0},
        {   2,   0,   0,   2,   0,     7.1,       1,    0.0,     0,   0.0},
        {   2,   0,   2,  -2,   1,    12.8,       1,    0.0,    -1,   0.0},
        {   2,   0,  -2,   0,   1,  -943.2,       1,    0.0,     0,   0.0},
        {   1,  -1,   0,  -2,   0,   -29.3,       1,    0.0,     0,   0.0},
        {  -1,   0,   0,   1,   1,  -388.3,       1,    0.0,     0,   0.0},
        {  -1,  -1,   0,   2,   1,    35.0,       1,    0.0,     0,   0.0},
        {   0,   1,   0,   1,   0,    27.3,       1,    0.0,     0,   0.0}
    };
    double ang;
    int i,j;

    *dpsi=*deps=0.0;

    for (i=0; i<106; i++)
    {
        ang=0.0;
        for (j=0; j<5; j++) ang+=nut[i][j]*f[j];
        *dpsi+=(nut[i][6]+nut[i][7]*t)*sin(ang);
        *deps+=(nut[i][8]+nut[i][9]*t)*cos(ang);
    }
    *dpsi*=1E-4*AS2R; /* 0.1 mas -> rad */
    *deps*=1E-4*AS2R;
}

/* search antenna parameter ----------------------------------------------------
* read satellite antenna phase center position
* args   : int    sat         I   satellite number (0: receiver antenna)
*          char   *type       I   antenna type for receiver antenna
*          gtime_t time       I   time to search parameters
*          pcvs_t *pcvs       IO  antenna parameters
* return : antenna parameter (NULL: no antenna)
*-----------------------------------------------------------------------------*/
extern pcv_t *searchpcv(int sat, const char *type, gtime_t time,
                        const pcvs_t *pcvs)
{
    pcv_t *pcv;
    char buff[MAXANT],*types[2],*p;
    int i,j,n=0;

    if (sat)
    {
        /* search satellite antenna */
        for (i=0; i<pcvs->n; i++)
        {
            pcv=pcvs->pcv+i;
            if (pcv->sat!=sat) continue;
            if (pcv->ts.time!=0&&timediff(pcv->ts,time)>0.0) continue;
            if (pcv->te.time!=0&&timediff(pcv->te,time)<0.0) continue;
            return pcv;
        }
    }
    else
    {
        strcpy(buff,type);
        for (p=strtok(buff," "); p&&n<2; p=strtok(NULL," ")) types[n++]=p;
        if (n<=0) return NULL;

        /* search receiver antenna with radome at first */
        for (i=0; i<pcvs->n; i++)
        {
            pcv=pcvs->pcv+i;
            for (j=0; j<n; j++) if (!strstr(pcv->type,types[j])) break;
            if (j>=n) return pcv;
        }
        /* search receiver antenna without radome */
        for (i=0; i<pcvs->n; i++)
        {
            pcv=pcvs->pcv+i;
            if (strstr(pcv->type,types[0])!=pcv->type) continue;

            return pcv;
        }
    }
    return NULL;
}

/* satellite number to satellite system ----------------------------------------
* convert satellite number to satellite system
* args   : int    sat       I   satellite number (1-MAXSAT)
*          int    *prn      IO  satellite prn/slot number (NULL: no output)
* return : satellite system (SYS_GPS,SYS_GLO,...)
*-----------------------------------------------------------------------------*/
extern int satsys(int sat, int *prn)
{
    int sys=SYS_NONE;
    if (sat<=0||MAXSAT<sat) sat=0;
    else if ( sat            <=MAXPRNGPS) sys=SYS_GPS;
    else if ((sat-=MAXPRNGPS)<=MAXPRNGLO) sys=SYS_GLO;

    else sat=0;
    if (prn) *prn=sat;
    return sys;
}
/* satellite number to satellite id --------------------------------------------
* convert satellite number to satellite id
* args   : int    sat       I   satellite number
*          char   *id       O   satellite id ( nn,Gnn,Rnn,Enn,Jnn,Cnn or nnn)
* return : none
*-----------------------------------------------------------------------------*/
extern void satno2id(int sat, char *id)
{
    int prn;
    switch (satsys(sat,&prn))
    {
    case SYS_GPS:
        sprintf(id,"%d" ,prn);
        return;
    case SYS_GLO:
        sprintf(id,"R%d",prn);
        return;
    }
    strcpy(id,"");
}



