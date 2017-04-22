#pragma once
#ifndef _TRACE_H_
#define _TRACE_H_
#include <stdio.h>
#include <tchar.h>

void trace( const TCHAR * format,... );
void traceMat( int m,int n,double Mat[] );
#endif // _DEBUG