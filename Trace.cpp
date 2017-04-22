#include <Windows.h>
#include "Trace.h"


void trace( const TCHAR * format,... )
{
	static const int BufferLen=1024;
	static int count;
	va_list pArg;
	TCHAR szMessageBuffer[BufferLen]={0};
	va_start(pArg,format);
	_vsntprintf(szMessageBuffer,BufferLen-1,format,pArg);
	va_end(pArg);
	OutputDebugString(szMessageBuffer);
}

void traceMat( int m,int n,double Mat[] )
{
	int i,j;
	trace("/**********************************************************************************************\n");
	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			trace("%12.3f",Mat[i*n+j]);	
		}
		trace("\n");
	}
	trace("**********************************************************************************************\\\n");
	return ;
}