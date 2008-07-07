C=======================================================================
C @(#) $Id: ReadLong.F,v 1.9 2004/11/30 16:50:58 jcmorale Exp $
C-----------------------------------------------------------------------
      SUBROUTINE ReadLong(DATTYP,NERR,WHAT,WHERE,WHEN,NSIZE,IOPT,
     &     RARR,IARR,CARR,LARR,DARR,SARR)
CC
CD Read long record array
CC
CC --------------------
CC        Input
CC --------------------
CC
CC WHAT      : physcial info
CC WHERE     : space info (Gem location/mesh info)
CC             is WHICH for material and property
CC WHEN      : time info (time step #)
CC IOPT      : integer flag
CC                 = __stop_if_failed__
CC                 = __not_stop_if_failed__
CC                 failed means: dataset not found or IO error. 
CC
CC --------------------
CC      Modified
CC --------------------
CC
CC NSIZE     : Total number of enties in the current array
CC
CC --------------------
CC        Output
CC --------------------
CC
CC pVAR      :  pointer to the data area
CC
CC --------------------
CC       Details
CC --------------------
CC
CC --------------------
CC Modification History
CC --------------------
CC
CC --------------------
CC Variable dictionary
CC --------------------
CC
CC --------------------
CC Local variables & parameters
CC --------------------
CC
CC LCL_STRCT :  Number of structures to read (eg: number of vertices
CC              or number of block members for a parallel ru )
CC LENGTH    :  Size of each structure (eg: 3 for grid coordinates)
CC
CC======================================================================
C
C ------------------------------
C     Preprocessor includes
C ------------------------------
C
#include "MMS.h"
#include "stack_point.h"
#include "locale_params.h"
#include "cfd_io.h"
#include "parallel_partitioning.h"
C
C -----------------------------
C       Subroutine Arguments
C -----------------------------
C
C
      INTEGER             NSIZE,WHEN,IOPT,NERR
C
      CHARACTER*(*)       WHAT,WHERE
C
      INTEGER             IARR(*)
      REAL                RARR(*)
      CHARACTER*(*)       CARR
      __io_string_type__  SARR(*)
      LOGICAL             LARR(*)
      DOUBLE PRECISION    DARR(*)
C
C ------------------------------
C       External routines
C ------------------------------
C
      EXTERNAL ERRMSG
C
C ---------------------------
C       Local Variables
C ---------------------------
C
      INTEGER             LDATTYP(IO_MAXENT),STRLN0,NSTRCT,
     &                    DATNUM(IO_MAXENT),NATTRB,IBLOCK,
     &                    DUMMY_INTR,I,NUMBLOCK,INTR_ARR(1),
     &                    LCL_NUM,LCL_STRCT,DATTYP,LENGTH,
     &                    LCL_BEG,LCL_END,NERRCM,LCL_POS,MSGTAG,NP,
     &                    DIM1,DIM2,IVAL
      __io_offset_type__  OFFSET(2)
C
      REAL                REAL_ARR(1)
C
      DOUBLE PRECISION    DBLE_ARR(1)
C
      CHARACTER*(1)       LOGL_ARR(1)
      CHARACTER*(4)       CRESLT, CDTYPE_BLK
      CHARACTER*(10)      CNP
      __io_string_type__  DUMMY_CHAR,CLENGTH,WHERE_tmp,
     &                    ATTRBS(IO_MAXATTRBS),CSTRCT
C
      LOGICAL             LOFND, LNewIO
C
C ------------------------------
C        Stack pointers.
C ------------------------------
C
      __stack_point__     pTMP,pTMPD,pBLK,pBLKD
C
C ---------------------------
C      Local Parameters
C ---------------------------
C
      CHARACTER*(*) ROUTIN
      PARAMETER    (ROUTIN = 'ReadLong')
C
C
C ------------------------------
C Directory on entry:- any
C ------------------------------
C
C ---------------------------
C  Executable Statements
C ---------------------------
C
C-----------------------------------------------------------------------
C     Use new ReadLong (optimized parallel IO performance)
C-----------------------------------------------------------------------
C
C---- Initialization.
C
      NERR   = 0
      NATTRB = 0
      NSTRCT = 0
      IBLOCK = 0
C
      CALL SETARI(LDATTYP,0,IO_MAXENT)
      CALL SETARI(DATNUM,0,IO_MAXENT) 

C
C-----------------------------------------------------------------------
C     Sequential run or parallel run without communication
C-----------------------------------------------------------------------
C
      IF (.TRUE. ) THEN
C
C---- Call REDHDR to read the header
C
        CALL REDHDR(LDATTYP,DATNUM,NSTRCT,NATTRB,ATTRBS,IBLOCK,
     &              NERR,WHAT,WHERE,WHEN,IOPT,__no_communication__)
C
        IF (NERR .NE. __io_ok__) THEN
          IF(IOPT .EQ. __stop_if_failed__) THEN
            CALL ERRMSG('*B*')
            CALL ERRMSG(ROUTIN//': REDHDR failed: what=')
            CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
            CALL ERRMSG(' where=')
            CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
            CALL ERRMSG('*E*')
            CALL CFXSTP(ROUTIN)
          ENDIF
          GOTO 999
        ENDIF
C
C---- Now only support the single data type, check LDATTYP(2), if
C     it is not -1, report error
C     
        IF(LDATTYP(2) .NE. -1)THEN
          CALL ERRMSG('*B*')
          CALL ERRMSG(ROUTIN//': more than one data type found what=')
          CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
          CALL ERRMSG(' where=')
          CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
          CALL ERRMSG('*E*')
          CALL CFXSTP(ROUTIN)
        ENDIF
C
C---- Set the length
C
        LENGTH = DATNUM(1)
C 
C---- check length of the dataset
C
C
C---- Read the datatag
C
        DUMMY_INTR = 0
        CALL IOCNTF_WRAP(__io_read_datatag__, NERR, ERRMSG,
     &                   WHAT, WHERE, WHEN, NATTRB,DUMMY_CHAR,
     &                   NSTRCT,CSTRCT,IBLOCK,
     &                   OFFSET,DUMMY_INTR,
     &                   DUMMY_CHAR,DUMMY_INTR,
     &                   DUMMY_CHAR,DUMMY_INTR,
     &                   INTR_ARR,DUMMY_INTR,
     &                   REAL_ARR,DUMMY_INTR,
     &                   DBLE_ARR,DUMMY_INTR,
     &                   LOGL_ARR,DUMMY_INTR)
C
        IF( NERR .NE. __io_ok__) THEN
          CALL ERRMSG('*B*')
          CALL ERRMSG(ROUTIN//': read data  failed: what=')
          CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
          CALL ERRMSG(' where=')
          CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
          CALL ERRMSG('*E*')
          CALL CFXSTP(ROUTIN)
        ENDIF
C
C---- Set block size: LCL_STRCT
C
        IF (IBLOCK.EQ.1 .OR. NSTRCT.LT.IBLOCK) THEN
          LCL_STRCT = NSTRCT
        ELSE 
          LCL_STRCT = IBLOCK
        ENDIF
C
C---- Calculate number of blocks: NUMBLOCK
C
        IF (MOD(NSTRCT,LCL_STRCT).EQ.0) THEN
          NUMBLOCK = NSTRCT/LCL_STRCT
        ELSE
          NUMBLOCK = NSTRCT/LCL_STRCT + 1
        ENDIF
C
C---- Misc
C
        DUMMY_CHAR = ' '
        DUMMY_INTR = 0
        INTR_ARR(1) = 0
        CALL STRI2S(CLENGTH,LENGTH)
C
C---- Loop over all blocks to be read
C
        LCL_POS = 1
C
        DO 120 I=1,NUMBLOCK
          LCL_BEG = (I-1)*LCL_STRCT + 1
          LCL_END = MIN(I*LCL_STRCT,NSTRCT)
          LCL_NUM = LCL_END - LCL_BEG + 1
C     
C---- Read one block of data
C
          IF (DATTYP .EQ. __real_data_type__) THEN
C  
            CSTRCT = 'STRC{REAL:'//CLENGTH(1:STRLN0(CLENGTH))//'}'
C
            CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                       WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                       LCL_NUM,CSTRCT,IBLOCK,
     &                       OFFSET,DUMMY_INTR,
     &                       DUMMY_CHAR,DUMMY_INTR,
     &                       DUMMY_CHAR,DUMMY_INTR,
     &                       INTR_ARR,DUMMY_INTR,
     &                       RARR((I-1)*LCL_STRCT*LENGTH+1),
     &                       LCL_NUM*LENGTH,
     &                       DBLE_ARR,DUMMY_INTR,
     &                       LOGL_ARR,DUMMY_INTR)
C
          ELSE IF (DATTYP .EQ. __int_data_type__) THEN
C
            CSTRCT = 'STRC{INTR:'//CLENGTH(1:STRLN0(CLENGTH))//'}'
C
            CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                       WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                       LCL_NUM,CSTRCT,IBLOCK,
     &                       OFFSET,DUMMY_INTR,
     &                       DUMMY_CHAR,DUMMY_INTR,
     &                       DUMMY_CHAR,DUMMY_INTR,
     &                       IARR((I-1)*LCL_STRCT*LENGTH+1),
     &                       LCL_NUM*LENGTH,
     &                       REAL_ARR,DUMMY_INTR,
     &                       DBLE_ARR,DUMMY_INTR,
     &                       LOGL_ARR,DUMMY_INTR)
C
          ELSE IF (DATTYP .EQ. __double_data_type__) THEN
C
            CSTRCT = 'STRC{DBLE:'//CLENGTH(1:STRLN0(CLENGTH))//'}'
C
C---- Double precision solver runs: Read DBLE data directly to RZ
C     Single precision solver runs: Read DBLE data to DZ and convert 
C
            CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                       WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                       LCL_NUM,CSTRCT,IBLOCK,
     &                       OFFSET,DUMMY_INTR,
     &                       DUMMY_CHAR,DUMMY_INTR,
     &                       DUMMY_CHAR,DUMMY_INTR,
     &                       INTR_ARR,DUMMY_INTR,
     &                       REAL_ARR,DUMMY_INTR,
     &                       DARR((I-1)*LCL_STRCT*LENGTH+1),
     &                       LCL_NUM*LENGTH,
     &                       LOGL_ARR,DUMMY_INTR)
C
          ELSE
C
C---- Why error: Invalid data type for long record
C
            CALL ERRMSG('*B*')
            CALL ERRMSG('ReadLong: invalid data type')
            CALL ERRMSG('*E*')
            NERR = 2
            GOTO 999
          ENDIF
C
C---- Check the error flag from IOCNTF
C
          IF( NERR .NE. __io_ok__) THEN
            CALL ERRMSG('*B*')
            CALL ERRMSG(ROUTIN//': read data failed: what=')
            CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
            CALL ERRMSG(' where=')
            CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
            CALL ERRMSG('*E*')
            CALL CFXSTP(ROUTIN)
          ENDIF
C
 120    CONTINUE
C
      ENDIF
C
      NSIZE = NSTRCT
C
 999  CONTINUE
      END
