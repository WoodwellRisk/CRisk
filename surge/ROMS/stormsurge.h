/*
** svn $Id: windbasin.h 1154 2023-02-17 20:52:30Z arango $
*******************************************************************************
** Copyright (c) 2002-2023 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** 
**
** Application flag:   STORMSURGE
** Input script:       roms.in

TIDAL SWITCHES:
#define SSH_TIDES
#define RAMP_TIDES
#

NO TIDES:

*/

#define UV_COR
#define UV_QDRAG
#define UV_ADV
#define DJ_GRADPS

#define ANA_FSOBC
#define ANA_M2OBC
#define ATM_PRESS
#define BODYFORCE
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define ANA_INITIAL
#define ANA_BTFLUX
#define ANA_STFLUX
#define ANA_SMFLUX
#define MASKING
#define CURVGRID
#define WET_DRY
#define LIMIT_BSTRESS