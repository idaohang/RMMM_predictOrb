#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ReadObs.h"
#include "RTOD_Const.h"

/***************************************************************************
//
// EmptyEpochObsDataStruct
//
// 目的: 初始化EpochObsData结构体, 在每次读数据之前调用
//
// 参数:
//
//  Epoch    待清空的观测数据

***************************************************************************/

void EmptyEpochObsDataStruct( EpochObsData* Epoch )
{
    int i, j;

    Epoch->EpochFlag = 0;
    Epoch->SatNum    = 0;
    Epoch->Time.Week = 0;
    Epoch->Time.SecOfWeek = 0.0;

    for( i=0; i<MAXCHANNUM; i++ )
    {
        Epoch->SatObs[i].Prn = 0;
        Epoch->SatObs[i].System = UNKS;
        Epoch->SatObs[i].Used   = 0;
        
        for( j=0; j<MAXOBSTYPENUM; j++ )
        {
            Epoch->SatObs[i].Data[j].Availability = 0;
            Epoch->SatObs[i].Data[j].Obs = 0.0;
            Epoch->SatObs[i].Data[j].Type = UNKOBS;
        }

    }

}

/***************************************************************************
//
// EmptyPSPRangeStruct
//
// 目的: 初始化PSPRange结构体, 在每次读数据之前调用
//
// 参数:
//
//  PSPR    待PSPRange结构体

***************************************************************************/

void EmptyPSPRangeStruct( PSPRange PSPR[] )
{
	int  i;

	for( i=0; i<32; i++ )
	{
		PSPR[i].CurrNum = 0;
		PSPR[i].CurrTime.Week = 0;
		PSPR[i].CurrTime.SecOfWeek = 0.0;

		PSPR[i].PSC1 = PSPR[i].PSP1 = PSPR[i].PSP2 = 0.0;
		PSPR[i].CPL1 = PSPR[i].CPL2 = 0.0;	

		PSPR[i].GoodC1 = PSPR[i].GoodP1 = PSPR[i].GoodP2 = false;
		PSPR[i].GoodL1 = PSPR[i].GoodL2 = false;

	}

}

/***************************************************************************
//
// ReadObsHead
//
// 目的: 读取Rinex文件头信息, 对于星载接收机的观测数据, 不使用文件头中的
         任何信息.
//
// 参数:
//
//  fin   观测数据文件指针
//
// 返回值:
// 
// 如果文件读取有问题, 返回0, 否则返回1
***************************************************************************/

int ReadObsHead( FILE* fin, OBSTYPELIST* ObsTypeList )
{
    int  i;
    char line[100], str[10], twoline[200];
    
    do{
        
        if( fgets( line, 100, fin )== NULL ) 
        {
            return 0;
        }
        
        if( strncmp( &line[60], "# / TYPES OF OBSERV", 18 ) == 0 )
        {
            strncpy( str, line, 6 );
            ObsTypeList->ObsTypeNum = atoi( str );

            strncpy( twoline, line, 60 );

            if( ObsTypeList->ObsTypeNum > 9 )
            {
                fgets( line, 100, fin );
                strncpy( &twoline[60], &line[6], 54 );
            }
            for( i=0; i<ObsTypeList->ObsTypeNum; i++ )
            { 
                strncpy( str, &twoline[10+6*i], 2 );
                
                if( strncmp( str, "C1", 2 ) == 0 )       ObsTypeList->ObsType[i] = C1;
                else if(strncmp( str, "P1", 2 ) == 0)    ObsTypeList->ObsType[i] = P1;
                else if(strncmp( str, "L1", 2 ) == 0)    ObsTypeList->ObsType[i] = L1;
                else if(strncmp( str, "D1", 2 ) == 0)    ObsTypeList->ObsType[i] = D1;
                else if(strncmp( str, "L2", 2 ) == 0)    ObsTypeList->ObsType[i] = L2;
                else if(strncmp( str, "P2", 2 ) == 0)    ObsTypeList->ObsType[i] = P2;
                else if(strncmp( str, "D2", 2 ) == 0)    ObsTypeList->ObsType[i] = D2;
                else if(strncmp( str, "S1", 2 ) == 0)    ObsTypeList->ObsType[i] = S1;
                else if(strncmp( str, "S2", 2 ) == 0)    ObsTypeList->ObsType[i] = S2;
                else if(strncmp( str, "C2", 2 ) == 0)    ObsTypeList->ObsType[i] = C2;
                else   ObsTypeList->ObsType[i] = UNKOBS;
            }
        } 
    }while( strncmp( &line[60],"END OF HEADER", 12 ) != 0 );
	
    return 1;
}

/***************************************************************************
//
// ReadEpochObs
//
// 目的: 读取一个历元的观测数据
//
// 参数:
//
//  fin   观测数据文件指针
//  Epoch 一个历元的观测数据
//
// 返回值:
// 
// 如果文件读取有问题, 返回0, 否则返回1
***************************************************************************/

int ReadEpochObs( FILE* fin, OBSTYPELIST* ObsTypeList, EpochObsData* Epoch )
{
    int i, j;
    char substr[60], line[100], twoline[200];
    char prnline[120];
    
    COMMONTIME CurrCT;
 
    do{
        if( fgets( line, 100, fin ) == NULL )
        {
            return 0;
        }
        sscanf( &line[26], "%d %d", &Epoch->EpochFlag, &Epoch->SatNum );

        if( Epoch->EpochFlag == 0 )
        {
            break;
        }
        else
        {
            for( i=0; i<Epoch->SatNum; i++ )
            {
                fgets( line, 100, fin );
            }
        }

    }while(1);

    
    sscanf( line, "%d %d %d %d %d %lf", &CurrCT.Year, &CurrCT.Month,
        &CurrCT.Day, &CurrCT.Hour, &CurrCT.Minute, &CurrCT.Second );

    sscanf( &line[26], "%d %d", &Epoch->EpochFlag, &Epoch->SatNum );
    
    CommonTimeToGPSTime ( &CurrCT, &Epoch->Time );
    
    strncpy( prnline, &line[32], 36 );
    
    if( Epoch->SatNum > 12 )
    {
        if( fgets( line, 100, fin ) == NULL )
        {
            return 0;
        }
        strncpy( &prnline[36], &line[32], 36 );
    }
    
    for( i=0; i<Epoch->SatNum; i++ )
    {
        if( fgets( line, 100, fin ) == NULL )
        {
            return 0;
        }
        
        strncpy( twoline, line, 80 );
        
        if( ObsTypeList->ObsTypeNum > 5 )
        {
            if( fgets( line, 100, fin ) == NULL )
            {
                return 0;
            }
            strncpy(&twoline[80], line, 80 );
        }
        
        if( prnline[3*i] == 'G'||prnline[3*i] == 'g'||prnline[3*i] == ' ')
        {
            Epoch->SatObs[i].System = GPS;
        }
        else if( prnline[3*i] == 'R'||prnline[3*i] == 'r' )
        {
            Epoch->SatObs[i].System = GLONASS;
        }

        strncpy( substr, &prnline[3*i+1], 2 );
        substr[2] = '\0';
        Epoch->SatObs[i].Prn    = atoi( substr );
        Epoch->SatObs[i].Used   = 1;             /* 开始均设置为可用 */
        
        for( j=0; j<ObsTypeList->ObsTypeNum; j++)
        {
            Epoch->SatObs[i].Data[j].Type = ObsTypeList->ObsType[j];
            
            strncpy( substr, &twoline[16*j], 14 );
            substr[14] = '\0';
            Epoch->SatObs[i].Data[j].Obs = atof( substr );

            if( fabs( Epoch->SatObs[i].Data[j].Obs ) < 0.01 )  /* 观测值为0或空 */
            {
                Epoch->SatObs[i].Data[j].Availability =  false;
            }
            else
            {
                Epoch->SatObs[i].Data[j].Availability =  true;
            }
        }
    }
    
    return 1;    
}
int ReadObsHeadQC( FILE* fin, OBSTYPELIST* ObsTypeList )
{
	int  i;
	char line[100], str[10], twoline[200];

	do{

		if( fgets( line, 100, fin )== NULL ) 
		{
			return 0;
		}

		if( strncmp( &line[60], "# / TYPES OF OBSERV", 18 ) == 0 )
		{
			strncpy( str, line, 6 );
			ObsTypeList->ObsTypeNum = atoi( str );

			strncpy( twoline, line, 60 );

			if( ObsTypeList->ObsTypeNum > 9 )
			{
				fgets( line, 100, fin );
				strncpy( &twoline[60], &line[6], 54 );
			}
			for( i=0; i<ObsTypeList->ObsTypeNum; i++ )
			{ 
				strncpy( str, &twoline[10+6*i], 2 );

				if( strncmp( str, "C1", 2 ) == 0 )       ObsTypeList->ObsType[i] = C1;
				else if(strncmp( str, "P1", 2 ) == 0)    ObsTypeList->ObsType[i] = P1;
				else if(strncmp( str, "L1", 2 ) == 0)    ObsTypeList->ObsType[i] = L1;
				else if(strncmp( str, "D1", 2 ) == 0)    ObsTypeList->ObsType[i] = D1;
				else if(strncmp( str, "L2", 2 ) == 0)    ObsTypeList->ObsType[i] = L2;
				else if(strncmp( str, "P2", 2 ) == 0)    ObsTypeList->ObsType[i] = P2;
				else if(strncmp( str, "D2", 2 ) == 0)    ObsTypeList->ObsType[i] = D2;
				else if(strncmp( str, "S1", 2 ) == 0)    ObsTypeList->ObsType[i] = S1;
				else if(strncmp( str, "S2", 2 ) == 0)    ObsTypeList->ObsType[i] = S2;
				else if(strncmp( str, "C2", 2 ) == 0)    ObsTypeList->ObsType[i] = C2;
				else if(strncmp( str, "Nw", 2 ) == 0)    ObsTypeList->ObsType[i] = Nw;
				else if(strncmp( str, "Sw", 2 ) == 0)    ObsTypeList->ObsType[i] = Sw;
				else   ObsTypeList->ObsType[i] = UNKOBS;
			}
		} 
	}while( strncmp( &line[60],"END OF HEADER", 12 ) != 0 );

	return 1;
}
int ReadEpochObsQC( FILE* fin, OBSTYPELIST* ObsTypeList, EpochObsData* Epoch )
{
	int i, j;
	char substr[60], line[200];
	char prnline[120];
	double tempTime;
	COMMONTIME CurrCT={0};

	do{
		if( fgets( line, 200, fin ) == NULL )
		{
			return 0;
		}
		sscanf( line, "%lf  %d", &tempTime, &Epoch->SatNum );
		Epoch->EpochFlag = 0;	
		if( Epoch->EpochFlag == 0 )
		{
			break;
		}
		else
		{
			for( i=0; i<Epoch->SatNum; i++ )
			{
				fgets( line, 200, fin );
			}
		}

	}while(1);

	SecTimeToCT(tempTime,&CurrCT);
	CommonTimeToGPSTime ( &CurrCT, &Epoch->Time );

	for( i=0; i<Epoch->SatNum; i++ )
	{
		if( fgets( line, 200, fin ) == NULL )
		{
			return 0;
		}
		if( line[1] == 'G'||line[1] == 'g'||line[1] == ' ')
		{
			Epoch->SatObs[i].System = GPS;
		}
		else if( prnline[1] == 'R'||prnline[1] == 'r' )
		{
			Epoch->SatObs[i].System = GLONASS;
		}

		strncpy( substr, &line[2], 2 );
		substr[2] = '\0';
		Epoch->SatObs[i].Prn    = atoi( substr );
		Epoch->SatObs[i].Used   = 1;             /* 开始均设置为可用 */

		for( j=0; j<ObsTypeList->ObsTypeNum; j++)
		{
			Epoch->SatObs[i].Data[j].Type = ObsTypeList->ObsType[j];

			strncpy( substr, &line[16*j+4], 14 );
			substr[14] = '\0';
			Epoch->SatObs[i].Data[j].Obs = atof( substr );
			strncpy(substr,&line[16*j+4+14],1);
			substr[1] = '\0';
			Epoch->SatObs[i].Data[j].LLI = atof( substr );
			if( fabs( Epoch->SatObs[i].Data[j].Obs ) < 0.01 )  /* 观测值为0或空 */
			{
				Epoch->SatObs[i].Data[j].Availability =  false;
				Epoch->SatObs[i].Data[j].LLI=9;
			}
			else
			{
				Epoch->SatObs[i].Data[j].Availability =  true;
			}
		}
	}

	return 1;    
}
/***************************************************************************
//
// SmoothPRangeWithPhase
//
// 目的: 用载波相位观测值平滑伪距观测值
//
// 参数:
//
//  PSPR  相位平滑伪距的中间结果
//  Epoch 一个历元的观测数据
//
// 注释:
// 
//   即使有双频观测数据，也采用单频伪距的平滑方法。考虑到电离层折射的影响，
//   平滑采用较短的时间段。超过时间限制，重新启动平滑器。

//    平滑的观测值覆盖历元观测数据的所有伪距，更新PSPR。
//    使用简单的周跳探测方法，检验相位数据和伪距数据的一致性。
***************************************************************************/

void SmoothPRangeWithPhase( PSPRange PSPR[], EpochObsData* Epoch )
{
	int i,j;
	int Prn;

	double l1, l2;
	double dt;

	bool bl1, bl2;  /* 各观测数据的可用性 */

	double dP, dL;
	bool   ResetFlagC1, ResetFlagP1, ResetFlagP2;

	for( i=0; i<Epoch->SatNum; i++ )
	{
		/* 先确定当前观测数据的卫星信息 */

		Prn = Epoch->SatObs[i].Prn;
		if( Epoch->SatObs[i].System != GPS )
		{
			continue;   /* GLONASS不平滑  */
		}
		
		/* 获取当前历元本卫星的L1、L2数据 */

		for( j=0; j<MAXOBSTYPENUM; j++ )
		{
			if( Epoch->SatObs[i].Data[j].Type == L1 )
			{
				l1 = Epoch->SatObs[i].Data[j].Obs;
				bl1 = false;
				if( Epoch->SatObs[i].Data[j].Availability == 1 )
				{
					bl1 = true;
				}
			}
			else if( Epoch->SatObs[i].Data[j].Type == L2 )
			{
				l2 = Epoch->SatObs[i].Data[j].Obs;
				bl2 = false;
				if( Epoch->SatObs[i].Data[j].Availability == 1 )
				{
					bl2 = true;
				}
			}
			else
			{
				continue;
			}
		}

		/* 对该卫星进行相位平滑，对每颗卫星进行分别平滑 */

		dt = ( Epoch->Time.Week - PSPR[Prn-1].CurrTime.Week ) * SECPERWEEK
			+ ( Epoch->Time.SecOfWeek - PSPR[Prn-1].CurrTime.SecOfWeek );

		if( dt > 3600.0 )   /* 超过3600秒，重新起步,保存数据  */
		{
			PSPR[Prn-1].CurrTime = Epoch->Time;

			PSPR[Prn-1].CPL1 = l1;
			PSPR[Prn-1].GoodL1 = bl1;

			PSPR[Prn-1].CPL2 = l2;
			PSPR[Prn-1].GoodL2 = bl2;

			ResetPhaseSmoothor( &Epoch->SatObs[i], &PSPR[Prn-1] );
		}
		else   /* 进行平滑  */
		{
			for( j=0; j<MAXOBSTYPENUM; j++ )
			{
				if( Epoch->SatObs[i].Data[j].Type == C1 
					&& Epoch->SatObs[i].Data[j].Availability == 1
					&& PSPR[Prn-1].GoodC1 == true )
				{
					ResetFlagC1 = false;
					dP = Epoch->SatObs[i].Data[j].Obs - PSPR[Prn-1].PSC1;
					dL = WAVELENGTHL1 * ( l1 - PSPR[Prn-1].CPL1 );
					if( fabs( dP - dL ) < 10.0 )
					{
						PSPR[Prn-1].PSC1 = Epoch->SatObs[i].Data[j].Obs / (PSPR[Prn-1].CurrNum+1) +
							PSPR[Prn-1].CurrNum*( PSPR[Prn-1].PSC1 + dL )/(PSPR[Prn-1].CurrNum+1);

						Epoch->SatObs[i].Data[j].Obs = PSPR[Prn-1].PSC1;
					}
					else  /* 重置平滑器 */
					{
						ResetFlagC1 = true;
					}
				}
				else if( Epoch->SatObs[i].Data[j].Type == P1 
					&& Epoch->SatObs[i].Data[j].Availability == 1
					&& PSPR[Prn-1].GoodP1 == true )
				{
					ResetFlagP1 = false;
			
					dP = Epoch->SatObs[i].Data[j].Obs - PSPR[Prn-1].PSP1;
					dL = WAVELENGTHL1 * ( l1 - PSPR[Prn-1].CPL1 );
					if( fabs( dP - dL ) < 10.0 )
					{
						PSPR[Prn-1].PSP1 = Epoch->SatObs[i].Data[j].Obs / (PSPR[Prn-1].CurrNum+1) +
							PSPR[Prn-1].CurrNum*( PSPR[Prn-1].PSP1 + dL )/(PSPR[Prn-1].CurrNum+1);
						
						Epoch->SatObs[i].Data[j].Obs = PSPR[Prn-1].PSP1;
					}
					else  /* 重置平滑器 */
					{
						ResetFlagP1 = true;
					}
				}
				else if( Epoch->SatObs[i].Data[j].Type == P2 
					&& Epoch->SatObs[i].Data[j].Availability == 1
					&& PSPR[Prn-1].GoodP2 == true )
				{
					ResetFlagP2 = false;
					dP = Epoch->SatObs[i].Data[j].Obs - PSPR[Prn-1].PSP2;
					dL = WAVELENGTHL2 * ( l2 - PSPR[Prn-1].CPL2 );
					if( fabs( dP - dL ) < 10.0 )
					{
						PSPR[Prn-1].PSP2 = Epoch->SatObs[i].Data[j].Obs / (PSPR[Prn-1].CurrNum+1) +
							PSPR[Prn-1].CurrNum*( PSPR[Prn-1].PSP2 + dL )/(PSPR[Prn-1].CurrNum+1);
						
						Epoch->SatObs[i].Data[j].Obs = PSPR[Prn-1].PSP2;
					}
					else  /* 重置平滑器 */
					{
						ResetFlagP2 = true;
					}
				}
				else
				{
					continue;
				}
			}

			if( ResetFlagC1 || ResetFlagP1 || ResetFlagP2 )
			{
				PSPR[Prn-1].CurrTime = Epoch->Time;
				
				PSPR[Prn-1].CPL1 = l1;
				PSPR[Prn-1].GoodL1 = bl1;
				
				PSPR[Prn-1].CPL2 = l2;
				PSPR[Prn-1].GoodL2 = bl2;
				
				ResetPhaseSmoothor( &Epoch->SatObs[i], &PSPR[Prn-1] );

			}
			else
			{
				PSPR[Prn-1].CurrNum = PSPR[Prn-1].CurrNum + 1;
				PSPR[Prn-1].CPL1 = l1;
				PSPR[Prn-1].GoodL1 = bl1;
				
				PSPR[Prn-1].CPL2 = l2;
				PSPR[Prn-1].GoodL2 = bl2;

			}
		}

	}

}

void ResetPhaseSmoothor( SatObsData* data, PSPRange* PSPR )
{
	int j;
	double c1, p1, p2;
	
	bool bc1, bp1, bp2;  /* 各观测数据的可用性 */
	PSPR->CurrNum = 1;
	
	for( j=0; j<MAXOBSTYPENUM; j++ )
	{
		if( data->Data[j].Type == C1 )
		{
			c1 = data->Data[j].Obs;
			bc1 = false;
			if( data->Data[j].Availability == 1 )
			{
				bc1 = true;
			}
			PSPR->PSC1 = c1;
			PSPR->GoodC1 = bc1;
		}
		else if( data->Data[j].Type == P1 )
		{
			p1 = data->Data[j].Obs;
			bp1 = false;
			if( data->Data[j].Availability == 1 )
			{
				bp1 = true;
			}
			PSPR->PSP1 = p1;
			PSPR->GoodP1 = bp1;
		}
		else if( data->Data[j].Type == P2 )
		{
			p2 = data->Data[j].Obs;
			bp2 = false;
			if( data->Data[j].Availability == 1 )
			{
				bp2 = true;
			}
			PSPR->PSP2 = p2;
			PSPR->GoodP2 = bp2;
		}
		else
		{
			continue;
		}
	}
}
/***************************************************************************
//
// CreateDopplerObs
//
// 目的: 用载波相位中心差分近似法计算Doppler观测值, 只适用用CHAMP卫星数据,其中没有
多普勒观测值
// 参数:
//
//
//  Epoch[3]         三个历元的观测数据
//
// 返回值:
// 
// 如果文件读取有问题, 返回0, 否则返回1
***************************************************************************/

int CreateDopplerObs( EpochObsData Epoch[3] )
{
	int i, j;
	int prn;
	double dt1, dt2;
	double l2, l1;

	dt1 = (Epoch[1].Time.Week-Epoch[0].Time.Week)*SECPERWEEK 
		+ Epoch[1].Time.SecOfWeek - Epoch[0].Time.SecOfWeek;
	dt2 = (Epoch[2].Time.Week-Epoch[1].Time.Week)*SECPERWEEK 
		+ Epoch[2].Time.SecOfWeek - Epoch[1].Time.SecOfWeek;

	if( fabs( dt2 - dt1 ) > 0.1 )
	{
		return 0;
	}

	for( i=0; i<Epoch[1].SatNum; i++ )
	{
		if( Epoch[1].SatObs[i].System == GPS )
		{
			prn = Epoch[1].SatObs[i].Prn;
			for( j=0; j<MAXOBSTYPENUM; j++ )
			{
				if( Epoch[1].SatObs[i].Data[j].Type == UNKOBS )
				{
					if( GetEpochSatL1( GPS, prn, &Epoch[0], &l1 ) && 
						GetEpochSatL1( GPS, prn, &Epoch[2], &l2) )
					{
						Epoch[1].SatObs[i].Data[j].Type = D1;
						Epoch[1].SatObs[i].Data[j].Obs = (l1-l2)/(dt2+dt1);
						Epoch[1].SatObs[i].Data[j].Availability = true;
						break;
					}
					
				}
			}
		}
				
	}

	return 1;
}

/***************************************************************************
//
// GetEpochSatL1
//
// 目的: 获取卫星的L1载波相位观测值, 临时函数, 仅用于用载波相位中心差分近似法
         计算Doppler观测值

***************************************************************************/

bool GetEpochSatL1( GNSSSys sys, short prn, EpochObsData* data, double* obs )
{
	int i, j;

	for( i=0; i<data->SatNum; i++ )
	{
		if( data->SatObs[i].System == sys && data->SatObs[i].Prn == prn )
		{
			for( j=0; j<MAXOBSTYPENUM; j++ )
			{
				if( data->SatObs[i].Data[j].Type == L1 )
				{
					*obs = data->SatObs[i].Data[j].Obs;
					return data->SatObs[i].Data[j].Availability;
				}
			}
		}		
	}

	return false;
}

 
/***************************************************************************
//
// GetOneSatPseudoRange
//
// 目的: 获取某颗卫星的伪距观测值，若是双频，使用组合观测值（C1+P2），单频返回C1
//
// 参数:
 
	data       卫星的观测数据
	pr         伪距观测值
	Ion        双频伪距组合计算的电离层改正值, 调试时用于检验单频电离层模型的
	           改正效率, 将来可以取消该参数，得出的是C1频段的电离层延迟
	
//
// 返回值:
// 
// 如果返回1, 单频C1码；返回2，双频组合伪距，返回0，表示有粗差
***************************************************************************/

int GetOneSatPseudoRange( GNSSSys Sys, SatObsData* data, double* pr, double *Ion )
{
	int i;
	double c1, p2;
	int ValidC1 = 0, ValidP2 = 0;
	int Val = 0;

	if( data->Used != 1 )
	{
		return 0;
	}
	for( i=0; i<MAXOBSTYPENUM; i++ )
	{
		if( (data->Data[i].Type == P1)
			&& ( data->Data[i].Availability == true ) )
		{
			ValidC1 = 1;
			c1 = data->Data[i].Obs;
		}
		if( ( data->Data[i].Type == P2 )
			&& ( data->Data[i].Availability == true ) )
		{
			ValidP2 = 1;
			p2 = data->Data[i].Obs;
		}
	}

	if( (ValidC1==1) && (ValidP2==1) )
	{
		if( fabs( p2 - c1 ) <= 30.0 )  /* P2有粗差，因为电离层延迟不会超过30m */
		{
			if( Sys == GPS )
			{
				*pr = ( FG12Ratio*FG12Ratio * c1 - p2 ) / ( FG12Ratio*FG12Ratio - 1.0 ); //双频无电离层组合
				*Ion = c1 - *pr;
			}
			else if( Sys == GLONASS )
			{
				printf("关于glonass清除错误!");
				system("pause");
			}
			Val = 2;
		}
		else
		{
			Val = 0;
		}
	}
	else if( ValidC1 == 1 )
	{
		*pr = c1;
		*Ion = 0.0;
		Val = 1;
	}
	else
	{
		Val = 0;
	}

	return Val;
}


