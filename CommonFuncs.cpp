/****************************************************************************
目的：    定义数组和矩阵运算的通用函数库

编写时间：2008.11.22
版本:     V1.1
版权：    武汉大学
****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include "commonfuncs.h"
#include "RTOD_Const.h"
double getModFirst(int length,double x[])
{
	double temp=0.0;
	for(int i=0;i<length;i++)
	{
		temp+=x[i]*x[i];
	}
	return sqrt(temp);
}
void randSed(){
	unsigned int seedVal;
	struct _timeb timeBuf ;
	_ftime (&timeBuf) ;
	seedVal = ( ( ( ( (unsigned int)timeBuf.time & 0xFFFF) +
		(unsigned int)timeBuf.millitm)^
		(unsigned int)timeBuf.millitm) );
	srand ((unsigned int)seedVal) ;
}
double randNum(double limitMin,double limitMax)
{
	double i ,j ,range ;
	range=limitMax-limitMin;    /* r is the range allowed; */
	i  =rand();    /* use the above example in this slot */
	/* Normalize the rand() output (scale to 0 to 1) */
	/* RAND_MAX is defined in stdlib, h */
	j  = ((double)i/(double)RAND_MAX) ;
	/* Scale the output to 1 to 44 */
	i  = (j * (double)range) ;
	i +=limitMin;
	return i;

}
/****************************************************************************
Rotation_x

目的：  计算绕X轴的旋转矩阵

参数:
Angle   旋转角[rad]
Mat     3*3阶旋转矩阵 
****************************************************************************/

void Rotation_x( double Angle, double Mat[] )
{
	double C, S;

	C = cos(Angle);
	S = sin(Angle);

	*(Mat+0*3+0) = 1.0; 
	*(Mat+0*3+1) = 0.0;  
	*(Mat+0*3+2) = 0.0;

	*(Mat+1*3+0) = 0.0; 
	*(Mat+1*3+1) =  +C;  
	*(Mat+1*3+2) =  +S;

	*(Mat+2*3+0) = 0.0; 
	*(Mat+2*3+1) =  -S;  
	*(Mat+2*3+2) =  +C;

}

/****************************************************************************
Rotation_y

目的：  计算绕X轴的旋转矩阵

参数:
Angle   旋转角[rad]
Mat     3*3阶旋转矩阵 
****************************************************************************/

void Rotation_y( double Angle, double Mat[] )
{
	double C, S;

	C = cos(Angle);
	S = sin(Angle);

	*(Mat+0*3+0) =  +C; 
	*(Mat+0*3+1) = 0.0;  
	*(Mat+0*3+2) =  -S;

	*(Mat+1*3+0) = 0.0; 
	*(Mat+1*3+1) = 1.0; 
	*(Mat+1*3+2) = 0.0;

	*(Mat+2*3+0) =  +S; 
	*(Mat+2*3+1) = 0.0;  
	*(Mat+2*3+2) =  +C;

}

/****************************************************************************
Rotation_z

目的：  计算绕X轴的旋转矩阵

参数:
Angle   旋转角[rad]
Mat     3*3阶旋转矩阵 
****************************************************************************/

void Rotation_z( double Angle, double Mat[] )
{
	double C, S;

	C = cos(Angle);
	S = sin(Angle);

	*(Mat+0*3+0) =  +C; 
	*(Mat+0*3+1) =  +S; 
	*(Mat+0*3+2) = 0.0;

	*(Mat+1*3+0) =  -S; 
	*(Mat+1*3+1) =  +C;  
	*(Mat+1*3+2) = 0.0;

	*(Mat+2*3+0) = 0.0;  
	*(Mat+2*3+1) = 0.0; 
	*(Mat+2*3+2) = 1.0;

}

/****************************************************************************
MatrixMultiply

目的：  矩阵相乘 M3 = M1*M2

参数:
m1      M1的行数
n1      M1的列数
m2      M2的行数
n2      M2的列数

****************************************************************************/

void MatrixMultiply( int m1, int n1, int m2, int n2, 
	const double M1[], const double M2[], double M3[] )
{    
	int i, j, k;
	double Sum;

	if( n1 != m2 )
	{
		printf("Matrix multiply fail!\n");
		exit(1);
	}

	for ( i=0; i<m1; i++) 
	{
		for ( j=0; j<n2; j++)
		{
			Sum = 0.0;

			for ( k=0; k<n1; k++)   
			{
				Sum += *(M1+i*n1+k) * *(M2+k*n2+j);
			}

			*(M3+i*n2+j) = Sum;
		}
	}
}

/****************************************************************************
MatrixMultiply2
功能：	矩阵和对角阵乘法
参数:
flag	表示相乘的顺序
flag=0:   M2 = M1*P
flag=1:    M2=P*M1
P       对角阵,用一维向量表示[m1]
m1      M1的行数
n1      M1的列数
length      M2矩阵的行列数[length=n1或length=m1]
M1      输入矩阵[m1*n1]

输出参数 

M2    输出矩阵[m1*n1]

****************************************************************************/
void MatrixMultiply2( int flag,int m1, int n1, int length, 
	const double M1[], const double P[],double M2[] )
{    
	int i, j;
	double Sum;
	if(0==flag){
		if(n1 != length)
			{
				printf("Matrix multiply fail!\n");
				exit(1);
			}
		for(i=0;i<m1;i++)
			for(j=0;j<n1;j++)
				M2[i*n1+j]=M1[i*n1+j]*P[j];

	}
	if(1==flag){
		if(m1 != length)
			{
				printf("Matrix multiply fail!\n");
				exit(1);
			}
		for(i=0;i<m1;i++)
			for(j=0;j<n1;j++)
				M2[i*n1+j]=M1[i*n1+j]*P[i];
	}
}
/****************************************************************************
MatrixMultiply3

目的：  矩阵相乘 M3 = M1*P*M2

参数:
m1      M1的行数
n1      M1的列数
m2      M2矩阵的行数[m2=n1]
n2      M2矩阵的列数
M1      输入矩阵[m1*n1]
P       对角阵,用一维向量表示[m1]
M2      输入矩阵[m2*n2]

输出参数 

M2    输出矩阵[n*n]

****************************************************************************/

void MatrixMultiply3( int m1, int n1, int m2, int n2, 
	const double M1[], const double P[], const double M2[], double M3[] )//没问题，是按行排列的矩阵的三个相乘，中间为对角矩阵用向量表示rj20160721
{    
	int i, j, k;
	double Sum;

	if( n1 != m2 )
	{
		printf("Matrix multiply fail!\n");
		exit(1);
	}

	for ( i=0; i<m1; i++) 
	{
		for ( j=0; j<n2; j++)
		{
			Sum = 0.0;

			for ( k=0; k<n1; k++)   
			{
				Sum += *(M1+i*n1+k) * *(P+k) * *(M2+k*n2+j);
			}

			*(M3+i*n2+j) = Sum;
		}
	}
}
/****************************************************************************
MatrixMultiplyk

目的: 常数乘上矩阵，通常为了单位权方差*权阵
参数: int m
参数: int n
参数: double M[]
参数: double k
****************************************************************************/
void MatrixMultiplyk(int m,int n,double M[],double k)
{
	int i,j;
	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			*(M+i*n+j) = *(M+i*n+j) * k; 
		}
	}
}
/****************************************************************************
ArrayAddition

目的：   向量相加 A3 = A1+A2

参数:
m      向量的的维数
A1      
A2
A3

****************************************************************************/
void ArrayAddition(int m,
	const double A1[], const double A2[], double A3[] )
{
	int i;
	for(i=0;i<m;i++)
	{
		A3[i]=A1[i]+A2[i];
	}
}
/****************************************************************************
MatrixAddition

目的：  矩阵相加 M3 = M1+M2

参数:
m      M1的行数
n      M1的列数

****************************************************************************/

void MatrixAddition( int m, int n,
	const double M1[], const double M2[], double M3[] )
{
	int i, j;

	for( i=0; i<m; i++ )
	{
		for( j=0; j<n; j++ )
		{
			*(M3+i*n+j) = *(M1+i*n+j) + *(M2+i*n+j);
		}
	}
}

/****************************************************************************
MatrixAddition2

目的：  矩阵相加 M2 = M1+M2

参数:
m      M1的行数
n      M1的列数

****************************************************************************/

void MatrixAddition2( int m, int n,
	const double M1[], double M2[] )
{
	int i, j;

	for( i=0; i<m; i++ )
	{
		for( j=0; j<n; j++ )
		{
			*(M2+i*n+j) = *(M1+i*n+j) + *(M2+i*n+j);
		}
	}
}

/****************************************************************************
MatrixAddition3

目的: 对角矩阵+方阵
参数: int n 方阵行列数
参数: const double M1[] 对角矩阵
参数: double M2[] 方阵
****************************************************************************/
void MatrixAddition3(int n,
	const double M1[], double M2[])
{
	int i;
	for (i=0;i<n;i++)
	{
		*(M2+i*n+i) = *(M1+i) + *(M2+i*n+i);
	}
}
/****************************************************************************
MatrixSubstract

目的：  矩阵相减 M3 = M1-M2

参数:
m      M1的行数
n      M1的列数

****************************************************************************/

void MatrixSubstract( int m, int n, 
	const double M1[], const double M2[], double M3[] )
{
	int i, j;

	for( i=0; i<m; i++ )
	{
		for( j=0; j<n; j++ )
		{
			*(M3+i*n+j) = *(M1+i*n+j) - *(M2+i*n+j);
		}
	}

}

/****************************************************************************
MatrixTranspose

目的：  矩阵转置 

参数:
m      M1的行数
n      M1的列数
M1     输入矩阵
MT     输出矩阵  MT = = M1(T)

待验证, 2008.12.8
****************************************************************************/

void MatrixTranspose( int m, int n, const double M1[], double MT[] )
{
	int i, j;
	for( i=0; i<m; i++ )
	{
		for( j=0; j<n; j++ )
		{
			*(MT+j*m+i) = *(M1+i*n+j);
		}
	}
}

/****************************************************************************
VectDot    a = A . B

目的：  计算两个向量的点积

参数:
m      A向量的元素个数
n      B向量的元素个数, 要求m=n

返回值  点积
****************************************************************************/

double VectDot( int m, int n, const double A[], const double B[] )
{
	int i;
	double Sum;

	if ( m != n )
	{
		printf("Vector dot fail!\n");
		exit(1);
	}

	Sum = 0.0;

	for( i=0; i<m; i++ )
	{
		Sum = Sum + *(A+i) * *(B+i);
	}

	return ( Sum );
}

/****************************************************************************
CrossDot    C = A X B

目的：  计算两个向量的叉积，只支持3向量叉乘3向量

参数:
m      A向量的元素个数
n      B向量的元素个数, 要求m=n


****************************************************************************/

void CrossDot( int m, int n, const double A[], const double B[], double C[] )
{ 
	if ( (n!=3) || (m!=3) )
	{
		printf("Cross dot fail!\n");
		exit(1);
	}

	C[0] = A[1]*B[2] - A[2]*B[1];
	C[1] = A[2]*B[0] - A[0]*B[2];
	C[2] = A[0]*B[1] - A[1]*B[0];

} 

/****************************************************************************
Dyadic    a = A . B

目的：  计算两个向量的Dyadic积

参数:
m      A向量的元素个数
n      B向量的元素个数

输出参数:
Mat    Dyadic积矩阵[m*n]

返回值  点积
****************************************************************************/

void Dyadic( int m, int n, const double A[], const double B[], double Mat[] )
{
	int i, j;

	for( i=0; i<m; i++ )
	{
		for( j=0; j<n; j++ )
		{
			*(Mat+i*n+j) = *(A+i) * *(B+j);
		}
	}
}

/****************************************************************************
SatElev

目的：  计算卫星高度角, 如果接收机位置没有初始化, 返回Pai/2.0

参数:

SatPos[3]  卫星位置[m]
RCVPos[3]  接收机位置[m]

返回值

卫星高度角[Rad]

****************************************************************************/

double SatElev( const double SatPos[], const double RCVPos[] )
{
	int    i;
	double Elev;
	double RcvR, SatRcvR;
	double dPos[3];       /* 卫星位置与接收机的位置差值 */


	for( i=0; i<3; i++ )
	{
		dPos[i] = SatPos[i] - RCVPos[i];
	}

	RcvR    = sqrt( VectDot(3, 3, RCVPos, RCVPos ) );
	SatRcvR = sqrt( VectDot(3, 3, dPos,   dPos ) );

	if( fabs( RcvR * SatRcvR ) <= 1.0 ) /* 检验卫星位置或接收机位置是否为0 */
	{
		Elev = pi / 2.0;
	}
	else
	{
		Elev = VectDot( 3, 3, RCVPos, dPos ) / ( RcvR * SatRcvR );

		Elev = pi/2.0 - acos( Elev );
	}

	return Elev;


}

/****************************************************************************
MatrixInv

目的：  矩阵求逆  

参数:
n      M1的行数和列数
a      输入矩阵
b      输出矩阵   b=inv(a)

****************************************************************************/

void MatrixInv(int n,double *a,double *b)
{
	double max0;
	int i,j,k,m;
	int r,s,ib[30],jb[30],ex;

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			b[i*n+j]=a[i*n+j];
		}
	}
	for(i=0;i<n;i++)
	{
		*(ib+i)=-1;
		*(jb+i)=-1;
	}
	for(k=0;k<n;k++)
	{
		max0=0.0;
		for(i=0;i<n;i++)
		{
			for(m=0;m<n;m++)
			{
				if (i==*(jb+m))
					goto loop1;
			}
			for(j=0;j<n;j++)
			{
				for (m=0;m<n;m++)
				{
					if(j==*(ib+m))   goto loop2;
				}
				if(max0<fabs(*(b+n*i+j)))
				{
					max0=fabs(*(b+n*i+j));
					r=i; s=j;
				}
loop2: ;
			}
loop1: ;
		}
		*(ib+r)=s; *(jb+s)=r;
		for(i=0;i<n;i++)
		{
			if(i!=r)
			{
				for(j=0;j<n;j++)
				{
					if(j!=s)
						*(b+n*i+j)=*(b+n*i+j)-*(b+i*n+s)*(*(b+n*r+j))/(*(b+n*r+s));
				}
				*(b+n*i+s)=*(b+i*n+s)/(*(b+n*r+s));
			}
		}
		for(j=0;j<n;j++)
		{
			if(j!=s)  *(b+n*r+j)=-(*(b+n*r+j))/(*(b+n*r+s));
		}
		*(b+n*r+s)=1/(*(b+n*r+s));
	}
	for(j=0;j<n-1;j++)
	{
		if(*(ib+j)!=j)
		{
			for (k=j+1;k<n;k++)
			{
				if(*(ib+k)==j)
				{
					for(i=0;i<n;i++)
					{
						max0=*(b+n*k+i);
						*(b+n*k+i)=*(b+n*j+i);
						*(b+n*j+i)=max0;
					}
					ex=*(ib+j);
					*(ib+j)=*(ib+k);
					*(ib+k)=ex;
				}
			}
		}
	}
	for(j=0;j<n-1;j++)
	{
		if(*(jb+j)!=j)
		{
			for(k=j+1;k<n;k++)
			{
				if(*(jb+k)==j)
				{
					for(i=0;i<n;i++)
					{
						max0=*(b+n*i+k);
						*(b+n*i+k)=*(b+n*i+j);
						*(b+n*i+j)=max0;
					}
					ex=*(jb+j);
					*(jb+j)=*(jb+k);
					*(jb+k)=ex;
				}
			}
		}
	}
}

/****************************************************************************
mbbub

目的：  实数冒泡排序  

参数:
n      待排序序列的长度
p      实数数组

****************************************************************************/
void mbbub( int n, double p[] )
{
	int m, k, i, j;
	double d;
	k=0;
	m=n-1;

	while (k<m)
	{
		j=m-1;
		m=0;
		for(i=k; i<=j; i++ )
		{
			if(p[i]>p[i+1])
			{
				d=p[i];
				p[i]=p[i+1];
				p[i+1]=d;
				m=i;
			}
		}

		j=k+1;
		k=0;
		for(i=m;i>=j;i--)
		{
			if(p[i-1]>p[i])
			{
				d=p[i];
				p[i]=p[i-1];
				p[i-1]=d;
				k=i;
			}
		}
	}  
}

/****************************************************************************
MatrixPrintf

目的：  将矩阵显示到屏幕上  

参数:
m, n   矩阵的行数和列数
Mat    待显示的矩阵

****************************************************************************/
void MatrixPrintf( int m, int n, double Mat[] )
{
	int i, j;

	for( i=0; i<m; i++ )
	{
		for( j=0; j<n; j++ )
		{
			printf( "  %lf ", *(Mat+i*n+j) );
		}
		printf("\n");
	}
}

/****************************************************************************
MatrixFprintf

目的：  将矩阵显示到文件中  

参数:
m, n   矩阵的行数和列数
Mat    待显示的矩阵

****************************************************************************/
void MatrixFprintf( int m, int n, FILE* Fout, double Mat[] )
{
	int i, j;

	for( i=0; i<m; i++ )
	{
		for( j=0; j<n; j++ )
		{
			fprintf( Fout, "  %lf ", *(Mat+i*n+j) );
		}
		fprintf( Fout, "\n" );
	}
}

/****************************************************************************
CopyArray

目的：  将一个数组拷贝到另一个数组中( 多维数组均以一维表示 )  

参数:
n      拷贝的数组元素个数
Dist   目标数组
Sour   源数组

****************************************************************************/
void CopyArray( int n, double Dist[], const double Sour[] )
{
	int i;

	for( i=0; i<n; i++ )
	{
		Dist[i] = Sour[i];
	}
}
/****************************************************************************
CopyMatrix

目的: 从Sour矩阵扣出子矩阵复制到Dist矩阵中
参数: int mLengh1   Dist矩阵行数
参数: int nLengh1   Dist矩阵列数
参数: int mstart1   Dist子矩阵起始行数
参数: int nstart1   Dist子矩阵起始列数
参数: double Dist[] 目的矩阵
参数: int mLengh2   Sour矩阵行数
参数: int nLengh2   Sour矩阵列数
参数: int mstart2   Sour子矩阵起始行数
参数: int nstart2   Sour子矩阵起始列数
参数: const double Sour[]
****************************************************************************/
void CopyMatrix(int mLengh1,int nLengh1,double Dist[],
	int mLengh2,int nLengh2,int mstart2,int nstart2,const double Sour[])
{
	int i,j,mstart,nstart;
	mstart=mstart2-1;
	nstart=nstart2-1;
	if (mstart2+mLengh1>mLengh2 || nstart2+nLengh1>nLengh2)
	{
		return;
	}
	for(i=0;i<mLengh1;i++)
	{
		for (j=0;j<nLengh1;j++)
		{
			Dist[i*nLengh1+j]=Sour[(mstart2+i)*nLengh2+nstart2+j];
		}
	}
}
void ArraySub(int m, const double A1[], const double A2[], double A3[])
{
	int i;
	for (i=0;i<m;i++)
	{
		A3[i]=A1[i]-A2[i];
	}
}
