C -*- mode: Fortran; coding: utf-8; indent-tabs-mode: nil; -*-
C kate: space-indent on; indent-width 2; encoding utf-8;
C kate: auto-brackets on; mixedindent off; indent-mode cstyle;

C=======================================================================
C @(#) $Id: RedSht.F,v 1.13 2002/05/28 09:16:22 cfdasm Exp $
C-----------------------------------------------------------------------
      SUBROUTINE RedSht(DATTYP,LENGTH,NERR,
     &                  WHAT,WHERE,WHEN,
     &                  NSIZE,
     &                  IOPT,IOPTPAR,
     &                  RARR,IARR,CARR,LARR,DARR,SARR)

C
CSCDT Read short record array
CSCDB
C     This routine read short record array  from a 
C       dataset specified by WHAT,WHERE,WHEN
C
C     Note: logical type is not supported!!!!
C 
C      Inputs:
C        WHAT,WHERE,WHEN: The signature of the dataset
C                         WHAT: physcial info
C                         WHERE: space info (Gem location/mesh info)
C                                is WHICH for material and property
C                         WHEN: time info (time step #)
C        IOPT   : integer flag
C                 = 1, __stop_if_failed__
C                 = 2, __not_stop_if_failed__
C                 failed means: dataset not found or IO error. 
C        IOPTPAR: integer flag for a parallel run
C                 = 1, __cm_to_all_partitions__
C                 = 2, __cm_to_one_partition__
C
C      Modified:
C        NSIZE  : Total number of enties in the current array
C        LENGTH : vector length of one entry
C        RARR, IARR, CARR, DARR, LARR, SARR are data arrays
C
C      Outputs:
C        DATTYP : data type =1, real; =2, integer; =3 char; =4, logical
C                 data type =-1, struct
C        NERR   : integer,
C                      0 - no error.

C=======================================================================
C
C---------------------------
C Preprocessor includes
C---------------------------
#include "cfd_io.h"
C #include "cfd_variable_keys.h"
C #include "cfd_internal_lookup.h"
#include "parallel_partitioning.h"
C
C---------------------------
C Subroutine Arguments
C---------------------------
C
      INTEGER             DATTYP,NSIZE,LENGTH,IOPTPAR,
     &                    WHEN,IOPT,
     &                    NERR,IARR(*), NERRCM
C
      REAL                RARR(*)
C
      CHARACTER*(*)       CARR,WHAT,WHERE
      __io_string_type__  SARR(*)
C
      LOGICAL             LARR(*)
C
      DOUBLE PRECISION    DARR(*)
C
C---------------------------
C Local Parameters
C---------------------------
C
      CHARACTER*(*) ROUTIN
      PARAMETER    (ROUTIN = 'RedSht')
C
C---------------------------
C Local Variables
C---------------------------
C
      INTEGER             LDATTYP(IO_MAXENT),STRLN0,NSTRCT,INTR_ARR(1),
     &                    DATNUM(IO_MAXENT),NATTRB,IBLOCK, LN,
     &                    DUMMY_INTR,I,MSGTAG, PTR, NPA, N,
     &                    N_int,N_real,N_char,N_log,N_double,N_str
      __io_offset_type__  OFFSET(2)
C
      REAL                REAL_ARR(1)
C
      CHARACTER*(1)       LOGL_ARR(1),LARR2(1)
      CHARACTER*(256)     CSTRCT
      __io_string_type__  DUMMY_CHAR,CLENGTH,ATTRBS(IO_MAXATTRBS),
     &                    STRI_ARR(1)
C
      DOUBLE PRECISION    DBLE_ARR(1)
C
      EXTERNAL ERRMSG

C
C---------------------------
C Executable Statements
C---------------------------
C
C---- Initialization.
C
      NERR   = 0
      NATTRB = 0
      NSTRCT = 0
      IBLOCK = 0
      CSTRCT = ' '
C
      CALL SETARI(LDATTYP,-1,IO_MAXENT)
      CALL SETARI(DATNUM,-1,IO_MAXENT)
C
      MSGTAG = __redsht_message__
C
C---- Skip for slave in parallel run
C
      IF (RUNMOD.EQ.__parallel__ .AND.
     &    PARMOD.EQ.__slave__) GOTO 720
C
C---- Call REDHDR to read the header
C
      CALL REDHDR(LDATTYP,DATNUM,NSTRCT,NATTRB,
     &            ATTRBS,IBLOCK,NERR,
     &            WHAT,WHERE,WHEN,IOPT,__no_communication__)
C
      IF (NERR .NE. __io_ok__) THEN
        IF (IOPT .EQ. __stop_if_failed__) THEN
          CALL ERRMSG('*B*')
          CALL ERRMSG('RedSht: REDHDR failed: what=')
          CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
          CALL ERRMSG(' where=')
          CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
          CALL ERRMSG('*E*')
          GOTO 998
        ENDIF
        GOTO 999
      ENDIF        
C
C---- Check if the data array larger enough or not
C
      DUMMY_INTR = DATNUM(1)
      IF (LDATTYP(1) .EQ. __string_data_type__ .AND. 
     &    DATNUM(1) .EQ. 80)THEN
        DUMMY_INTR = 1
      ENDIF
C
      IF (NSTRCT*DUMMY_INTR .GT. NSIZE*LENGTH)THEN
        NERR = __io_memory_err__
        IF (IOPT .EQ. __stop_if_failed__) THEN
          CALL ERRMSG('*B*')
          CALL ERRMSG('RedSht: size of data array is too small: what=')
          CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
          CALL ERRMSG(' where=')
          CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
          CALL ERRMSG('*E*')
          GOTO 998
        END IF
        GOTO 999
      ENDIF        
C
C---- Set the data type 
C     --> If LDATTYP(2) is not -1,means the data type is a struct
C 
      DATTYP = LDATTYP(1)
      IF (LDATTYP(2) .NE. -1) DATTYP = -1
C
C---- Set the length
C
      LENGTH = DUMMY_INTR         
C 
C---- Get length of the dataset
C     
      NSIZE = NSTRCT
C     
C     read data body
C---- Get dataset type DATTYP
C     
      DUMMY_CHAR = ' '
      DUMMY_INTR = 0
      INTR_ARR(1)= 0
      CALL STRI2S(CLENGTH,DATNUM(1)) 
C
C---- Read datatag
C
      IF (NSTRCT .LE. 0) GOTO 999
      CALL IOCNTF_WRAP(__io_read_datatag__,NERR,ERRMSG,
     &                 WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                 NSTRCT,CSTRCT,IBLOCK,
     &                 OFFSET,DUMMY_INTR,
     &                 STRI_ARR,DUMMY_INTR,
     &                 DUMMY_CHAR,DUMMY_INTR,
     &                 INTR_ARR,DUMMY_INTR,
     &                 REAL_ARR,DUMMY_INTR,
     &                 DBLE_ARR,DUMMY_INTR,
     &                 LOGL_ARR,DUMMY_INTR)
C
      IF (NERR .NE. __io_ok__) THEN
        IF (IOPT .EQ. __stop_if_failed__) THEN
          CALL ERRMSG('*B*')
          CALL ERRMSG('RedSht: read datatag failed: what=')
          CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
          CALL ERRMSG(' where=')
          CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
          CALL ERRMSG('*E*')
          GOTO 998
        END IF
        GOTO 999
      ENDIF
C
      IF (DATTYP .EQ. __real_data_type__) THEN
        CSTRCT = 'STRC{REAL:'//CLENGTH(1:STRLN0(CLENGTH))//'}'
        CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                   WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                   NSTRCT,CSTRCT,IBLOCK,
     &                   OFFSET,DUMMY_INTR,
     &                   STRI_ARR,DUMMY_INTR,
     &                   DUMMY_CHAR,DUMMY_INTR,
     &                   INTR_ARR,DUMMY_INTR,
     &                   RARR,LENGTH*NSTRCT,
     &                   DBLE_ARR,DUMMY_INTR,
     &                   LOGL_ARR,DUMMY_INTR)
C
        IF (NERR .NE. __io_ok__) THEN
          IF (IOPT .EQ. __stop_if_failed__) THEN
            CALL ERRMSG('*B*')
            CALL ERRMSG('RedSht: read data failed: what=')
            CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
            CALL ERRMSG(' where=')
            CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
            CALL ERRMSG('*E*')
            GOTO 998
          END IF
          GOTO 999
        ENDIF
C
      ELSE IF (DATTYP .EQ. __int_data_type__) THEN
        CSTRCT = 'STRC{INTR:'//CLENGTH(1:STRLN0(CLENGTH))//'}'
        CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                   WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                   NSTRCT,CSTRCT,IBLOCK,
     &                   OFFSET,DUMMY_INTR,
     &                   STRI_ARR,DUMMY_INTR,
     &                   DUMMY_CHAR,DUMMY_INTR,
     &                   IARR,LENGTH*NSTRCT,
     &                   REAL_ARR,DUMMY_INTR,
     &                   DBLE_ARR,DUMMY_INTR,
     &                   LOGL_ARR,DUMMY_INTR)      
        IF (NERR .NE. __io_ok__) THEN
          IF (IOPT .EQ. __stop_if_failed__) THEN
            CALL ERRMSG('*B*')
            CALL ERRMSG('RedSht: read data failed: what=')
            CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
            CALL ERRMSG(' where=')
            CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
            CALL ERRMSG('*E*')
            GOTO 998
          END IF
          GOTO 999
        ENDIF
C
      ELSE IF (DATTYP .EQ. __char_data_type__) THEN
        CSTRCT = 'STRC{CHAR:'//CLENGTH(1:STRLN0(CLENGTH))//'}'
        CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                   WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                   NSTRCT,CSTRCT,IBLOCK,
     &                   OFFSET,DUMMY_INTR,
     &                   STRI_ARR,DUMMY_INTR,
     &                   CARR,LENGTH*NSTRCT,
     &                   INTR_ARR,DUMMY_INTR,
     &                   REAL_ARR,DUMMY_INTR,
     &                   DBLE_ARR,DUMMY_INTR,
     &                   LOGL_ARR,DUMMY_INTR)           
        IF (NERR .NE. __io_ok__) THEN
          IF (IOPT .EQ. __stop_if_failed__) THEN
            CALL ERRMSG('*B*')
            CALL ERRMSG('RedSht: read data failed: what=')
            CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
            CALL ERRMSG(' where=')
            CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
            CALL ERRMSG('*E*')
            GOTO 998
          END IF
          GOTO 999
        ENDIF       
      ELSE IF (DATTYP .EQ. __double_data_type__) THEN
c        write(*,*) 'double reader'
c        nstrct = 100
c        write(*,*) nstrct
        CSTRCT = 'STRC{DBLE:'//CLENGTH(1:STRLN0(CLENGTH))//'}'
        CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                   WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                   NSTRCT,CSTRCT,IBLOCK,
     &                   OFFSET,DUMMY_INTR,
     &                   STRI_ARR,DUMMY_INTR,
     &                   DUMMY_CHAR,DUMMY_INTR,
     &                   INTR_ARR,DUMMY_INTR,
     &                   REAL_ARR,DUMMY_INTR,
     &                   DARR,LENGTH*NSTRCT,
     &                   LOGL_ARR,DUMMY_INTR)           
        IF (NERR .NE. __io_ok__) THEN
          IF (IOPT .EQ. __stop_if_failed__) THEN
            CALL ERRMSG('*B*')
            CALL ERRMSG('RedSht: read data failed: what=')
            CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
            CALL ERRMSG(' where=')
            CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
            CALL ERRMSG('*E*')
            GOTO 998
          END IF
          GOTO 999
        ENDIF       
C
      ELSE IF (DATTYP .EQ. __logical_data_type__) THEN
        CSTRCT = 'STRC{LOGL:'//CLENGTH(1:STRLN0(CLENGTH))//'}'
        CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                   WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                   NSTRCT,CSTRCT,IBLOCK,
     &                   OFFSET,DUMMY_INTR,
     &                   STRI_ARR,DUMMY_INTR,
     &                   DUMMY_CHAR,DUMMY_INTR,
     &                   INTR_ARR,DUMMY_INTR,
     &                   REAL_ARR,DUMMY_INTR,
     &                   DBLE_ARR,DUMMY_INTR,
     &                   LARR2,LENGTH*NSTRCT)
        IF (NERR .NE. __io_ok__) THEN
          IF (IOPT .EQ. __stop_if_failed__) THEN
            CALL ERRMSG('*B*')
            CALL ERRMSG('RedSht: read data failed: what=')
            CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
            CALL ERRMSG(' where=')
            CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
            CALL ERRMSG('*E*')
            GOTO 998
          END IF
          GOTO 999
        ENDIF       
C
      ELSE IF (DATTYP .EQ. __string_data_type__) THEN
        CSTRCT = 'STRC{STRI:'//CLENGTH(1:STRLN0(CLENGTH))//'}'
        CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                   WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                   NSTRCT,CSTRCT,IBLOCK,
     &                   OFFSET,DUMMY_INTR,
     &                   SARR,NSTRCT,
     &                   DUMMY_CHAR,DUMMY_INTR,
     &                   INTR_ARR,DUMMY_INTR,
     &                   REAL_ARR,DUMMY_INTR,
     &                   DBLE_ARR,DUMMY_INTR,
     &                   LOGL_ARR,DUMMY_INTR)           
        IF (NERR .NE. __io_ok__) THEN
          IF (IOPT .EQ. __stop_if_failed__) THEN
            CALL ERRMSG('*B*')
            CALL ERRMSG('RedSht: read data failed: what=')
            CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
            CALL ERRMSG(' where=')
            CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
            CALL ERRMSG('*E*')
            GOTO 998
          END IF
          GOTO 999
        ENDIF       
C                     for struct
      ELSE IF (DATTYP .EQ. -1) THEN
        CSTRCT = 'STRC{'
        N_int = 0
        N_real = 0
        N_char = 0
        N_log = 0
        N_double = 0
        N_str = 0
        DO 200 I=1,IO_MAXENT
          CALL STRI2S(CLENGTH,DATNUM(I)) 
          IF (LDATTYP(I) .EQ. -1) GOTO 200
          IF (LDATTYP(I) .EQ. __real_data_type__) THEN
            IF (STRLN0(CSTRCT) .EQ. 5) THEN
              CALL STRCTN(CSTRCT,'REAL:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ELSE
              CALL STRCTN(CSTRCT,' REAL:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ENDIF
            N_real = N_real + DATNUM(I)
          ELSE IF (LDATTYP(I) .EQ. __int_data_type__) THEN
            IF (STRLN0(CSTRCT) .EQ. 5) THEN
              CALL STRCTN(CSTRCT,'INTR:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ELSE
              CALL STRCTN(CSTRCT,' INTR:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ENDIF
            N_int = N_int + DATNUM(I)
          ELSE IF (LDATTYP(I) .EQ. __char_data_type__) THEN
            IF (STRLN0(CSTRCT) .EQ. 5) THEN
              CALL STRCTN(CSTRCT,'CHAR:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ELSE
              CALL STRCTN(CSTRCT,' CHAR:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ENDIF
            N_char = N_char + DATNUM(I)
          ELSE IF (LDATTYP(I) .EQ. __double_data_type__) THEN
            IF (STRLN0(CSTRCT) .EQ. 5) THEN
              CALL STRCTN(CSTRCT,'DBLE:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ELSE
              CALL STRCTN(CSTRCT,' DBLE:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ENDIF
            N_double = N_double + DATNUM(I)
          ELSE IF (LDATTYP(I) .EQ. __logical_data_type__) THEN
            IF (STRLN0(CSTRCT) .EQ. 5) THEN
              CALL STRCTN(CSTRCT,'LOGL:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ELSE
              CALL STRCTN(CSTRCT,' LOGL:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ENDIF
            N_log = N_log + DATNUM(I)
          ELSE IF (LDATTYP(I) .EQ. __string_data_type__) THEN
            IF (STRLN0(CSTRCT) .EQ. 5) THEN
              CALL STRCTN(CSTRCT,'STRI:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ELSE
              CALL STRCTN(CSTRCT,' STRI:'//CLENGTH(1:STRLN0(CLENGTH)) )
            ENDIF
            N_str = N_str + DATNUM(I)
          ENDIF
C            
 200    CONTINUE
        CALL STRCTN(CSTRCT,'}')
C
        CALL IOCNTF_WRAP(__io_read_data__,NERR,ERRMSG,
     &                   WHAT,WHERE,WHEN,NATTRB,DUMMY_CHAR,
     &                   NSTRCT,CSTRCT,IBLOCK,
     &                   OFFSET,DUMMY_INTR,
     &                   SARR,N_str*NSTRCT,
     &                   CARR,N_char*NSTRCT,
     &                   IARR,N_int*NSTRCT,
     &                   RARR,N_real*NSTRCT,
     &                   DARR,N_double*NSTRCT,
     &                   LARR2,N_log*NSTRCT)
        IF (NERR .NE. __io_ok__) THEN
          IF (IOPT .EQ. __stop_if_failed__) THEN
            CALL ERRMSG('*B*')
            CALL ERRMSG('RedSht: read data failed: what=')
            CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
            CALL ERRMSG(' where=')
            CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
            CALL ERRMSG('*E*')
            GOTO 998
          END IF
          GOTO 999
        ENDIF       
C
      ELSE        
C     
C---> Why error: Invalid data type for short record
C
        PRINT*,'Invalid data type for short record'
        NERR = 2
        GOTO 999        
      ENDIF
C
 720  CONTINUE
C
C---- Exit.
C     
 999  CONTINUE
C
      RETURN
C
 998  CONTINUE
C       
      CALL CFXSTP(ROUTIN)
C
      END
