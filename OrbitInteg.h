/****************************************************************************
Ŀ�ģ�    ����ѧ�������,ʹ��RK4, RKF4,�Լ������ֵ�Ⱥ���

��дʱ�䣺2008.11.26
�汾:     V1.1
��Ȩ��    �人��ѧ
****************************************************************************/
#pragma once 
#ifndef _ORBIT_INTEGRATION_H_
#define _ORBIT_INTEGRATION_H_

#include "DynaModel.h"
#include "GPSTime.h"



/***************************************************************************
//
// RK4Step
//
// Purpose:
//
//   ����-����4�׵�������������
//
// Input/Output:
//
//   Mjd_GPS    ��ʼʱ��MJDʱ��(GPS Time),���ֺ�Ϊmjd_gps+step
//   Step       ����ʱ����[s]
//   Y0         ���ֳ�ֵ,Mjd_GPSʱ�̵�����λ�ú��ٶ�[m, m/s],
//                 ���ֺ�Ϊ��һʱ�̵�����λ�ú��ٶ�
//   Para       ����ѧģ�Ͳ���
//
***************************************************************************/
//void RK4Step( MJDTIME* Mjd_GPS, double step, double Y0[6], DYNMODELPARA* Para );

/***************************************************************************
//
// RKF4Step
//
// Purpose:
//
//   ����-����-Felberg4�׵�������������
//
// Input/Output:
//
//   Mjd_GPS    ��ʼʱ��MJDʱ��(GPS Time),���ֺ�Ϊmjd_gps+step
//   Step       ����ʱ����[s]
//   Y0         ���ֳ�ֵ,Mjd_GPSʱ�̵�����λ�ú��ٶ�[m, m/s],
//                 ���ֺ�Ϊ��һʱ�̵�����λ�ú��ٶ�
//   Para       ����ѧģ�Ͳ���
//
***************************************************************************/
//void RKF4Step( MJDTIME* Mjd_GPS, double step, double Y0[6], DYNMODELPARA* Para );

/***************************************************************************
//
// RKF4OrbitSTM
//
// Purpose:
//
//   ����-����-Felberg4�׵�������������, ���й����״̬ת�ƾ������
//
// Input/Output:
//
//   Mjd_GPS    ��ʼʱ��MJDʱ��(GPS Time),���ֺ�Ϊmjd_gps+step
//   Step       ����ʱ����[s]
//   Y0         �����ת�ƾ����ֵ,Mjd_GPSʱ�̵�����λ�ú��ٶ�[m, m/s],
//                 ���ֺ�Ϊ��һʱ�̵�����λ�ú��ٶ�[48ά]
//   IntState   ��Hermit5�����ڲ�����Ҫ��״̬�������ڴ˴���ʼ��[λ��/�ٶ�/���ٶ�]
//   Para       ����ѧģ�Ͳ���
//
***************************************************************************/
void RKF4OrbitSTM( int graceType, MJDTIME* Mjd_GPS, double step, double Y0[54], SCState Stat[2], DYNMODELPARA* Para );

/***************************************************************************
//
// Hermite5
//
// Purpose:
//
//   Hermite 5�׶���ʽ�ڲ����ǵĹ��
//
// Input/Output:
//
//   S0         ��ʼʱ�����ǵ�״̬����[λ��/�ٶ�/���ٶ�]
//   S1         ��ֹʱ�����ǵ�״̬����
//   CurrState  �ڲ�ʱ�̵�����״̬
****************************************************************************/
void Hermite5( const SCState* S0, const SCState* S1, SCState* CurrState );

/***************************************************************************
//
// OrbitIntegToGivenTime
//
// Purpose:
//
//   ʹ��RKF4�������������������й�����֣�Ԥ������ʱ�����ǵĹ��
//
// Input/Output:
//
//   Mjd_GPS             ��ʼʱ������MJDʱ���ʾ��GPSʱ
//   Mjd_GivenTime       ��ҪԤ������Ĺ۲�ʱ�̣���Mjd_GPS��������ͬ
//   Y0                  ���ֳ�ֵ,Mjd_GPSʱ�̵�����λ�ú��ٶ�[m, m/s],
//                              ���ֺ�Ϊ��һʱ�̵�����λ�ú��ٶ�
//   Para                ����ѧģ�Ͳ���
//
***************************************************************************/
//void OrbitIntegToGivenTime( MJDTIME* Mjd_GPS, MJDTIME* Mjd_GivenTime,
//						   double Y0[6], DYNMODELPARA* Para );


/***************************************************************************
//
// InitStateTranMatrix
//
// Purpose:
//
//   ��ʼ��������������е�״̬ת�ƾ���
//
// Input/Output:
//
//   row    ״̬ת�ƾ��������
     col    ״̬ת�ƾ��������
	 STM    ״̬ת�ƾ�������

****************************************************************************/

void InitStateTranMatrix( int row, int col, double STM[] );


#endif