#pragma once 
#ifndef _SATODS_H_
#define _SATODS_H_

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "gpstime.h"
#include "CommonFuncs.h"
#include "rtod_const.h"
#include "RefSys.h"
#include "DynaModel.h"
#include "OrbitInteg.h"
#include "ReadObs.h"
#include "EphProc.h"
#include "RTOD.h"
#include "OrbParaTran.h"
#include "ReadLevel1B.h"
#include "PhaseCent_Cor.h"

void KFStateUpdate(double * xA_pre,double * xB_pre,EKFSTATE *  KFState);
#endif