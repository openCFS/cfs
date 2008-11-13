C -*- mode: Fortran; coding: utf-8; indent-tabs-mode: nil; -*-
C kate: space-indent on; indent-width 2; encoding utf-8;
C kate: auto-brackets on; mixedindent off; indent-mode cstyle;

C=======================================================================
C @(#) $Id: RedHdr.F,v 1.10 2002/05/22 15:09:15 cfdasm Exp $
C-----------------------------------------------------------------------
      SUBROUTINE REDHDR(DATTYP,DATNUM,NSTRCT,NATTRB,
     &                  ATTRBS,IBLOCK,NERR,
     &                  WHAT,WHERE,WHEN,IOPT,IOPTPAR)
C
CSCDT Read the dataset header
CSCDB
C     This routine read dataset header from a 
C       dataset specified by WHAT,WHERE,WHEN
C           
C 
C      Inputs:
C        WHAT,WHERE,WHEN: The signature of the dataset
C                         WHAT: physcial info
C                         WHERE: space info (Gem location/mesh info)
C                                is WHICH for material and property
C                         WHEN: time info (time step #)
C        IOPTPAR: parallel option:
C                 = 0. __no_communication__
C                 = 2. __ds_to_partitions__
C      Outputs:
C        DATTYP : integer array,it contains the enumerate datatype
C                 for each entry, eg,
C                 1 if for real data, which is defined in cfd_io.h,
C                 __real_data_type__
C        DATNUM : integer array, number of data in each entry 
C                 of the structre,eg.
C                 if CSTRCT is STRC{INTR:3 REAL:2}, the DATNUN is 3 for the first
C                 entry, is 2 for the second entry
C        NATTRB : integer, is the number of the attrbutes
C        ATTRBS : character*80 array, contains all attributes strings
C        IBLOCK : integer, block factor
C        NERR   : integer, 
C                      0 - no error. 
C                      > 0 - an error code which defined in cfc_io.h
C
CSCDE
C=======================================================================
C
C---------------------------
C Preprocessor includes
C---------------------------
C
#include "cfd_io.h"
C #include "cfd_internal_lookup.h"
#include "parallel_partitioning.h"
C
C---------------------------
C Subroutine Arguments
C---------------------------
C
      INTEGER             DATTYP(IO_MAXENT), NSTRCT,IOPT,IOPTPAR,NERR,
     &                    DATNUM(IO_MAXENT), WHEN, NATTRB,IBLOCK
C
      CHARACTER*(*)       ATTRBS(IO_MAXATTRBS),WHAT,WHERE
C
C---------------------------
C Local Parameters
C---------------------------
C
      CHARACTER*(*) ROUTIN
      PARAMETER    (ROUTIN = 'REDHDR')
C
C---------------------------
C Local Variables
C---------------------------
C
      INTEGER             I,NUMINTR,NUMOFFSET,NUMSTRING,NUMCHAR,
     &                    INTR_ARR(1),DUMMY_INTR,
     &                    INDEX,START,END,LOCATION,STRLN0,
     &                    LOCAL_WHEN
C
      __io_offset_type__  OFFSET(2)
C
      REAL                REAL_ARR(1)
C
      CHARACTER*(1)       LOGL_ARR(1),CHAR_ARR(1)
      __io_string_type__  STRI_ARR(2),DUMMY2,LOCAL_WHAT,LOCAL_WHERE
      CHARACTER*(256)     CSTRCT,DUMMY,TMPSTR
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
      NATTRB = 0
      NERR   = 0
      CSTRCT = ' '
C
      CALL SETARI(DATTYP,-1,IO_MAXENT)
      CALL SETARI(DATNUM,-1,IO_MAXENT)
C
C---- Skip in case of parallel slave
C
      IF (PARMOD .EQ. __slave__) GOTO 999
C
C---- Call IOCNTF to locate dataset
C
      OFFSET(1) = 0
      OFFSET(2) = 0
      NUMOFFSET = 2
      NUMSTRING = 2
      DUMMY     = ' '
      DUMMY_INTR = 0
      NUMINTR   = 1
C
      CALL IOCNTF_WRAP(__io_locate_dataset__,NERR,ERRMSG,
     &                 WHAT,WHERE,WHEN,NATTRB,ATTRBS,
     &                 NSTRCT,CSTRCT,IBLOCK,
     &                 OFFSET,NUMOFFSET,
     &                 STRI_ARR,NUMSTRING,
     &                 CHAR_ARR,DUMMY_INTR,
     &                 INTR_ARR,NUMINTR,
     &                 REAL_ARR,DUMMY_INTR,
     &                 DBLE_ARR,DUMMY_INTR,
     &                 LOGL_ARR,DUMMY_INTR)
C
      IF( NERR .NE. __io_ok__ ) THEN
        IF(IOPT .EQ. __stop_if_failed__) THEN
          CALL ERRMSG('*B*')
          CALL ERRMSG('REDHDR: locating dataset failed: what=')
          CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
          CALL ERRMSG(' where=')
          CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
          CALL ERRMSG('*E*')
          GOTO 998
        ENDIF
        GOTO 999
      ENDIF
C
C---- Call IOCNTF to read the located  dataset header
C         
      NUMOFFSET = 1  
      CALL IOCNTF_WRAP(__io_read_header__,NERR,ERRMSG,
     &                 LOCAL_WHAT,LOCAL_WHERE,LOCAL_WHEN,NATTRB,ATTRBS,
     &                 NSTRCT,CSTRCT,IBLOCK,
     &                 OFFSET,NUMOFFSET,
     &                 STRI_ARR,NUMSTRING,
     &                 CHAR_ARR,NUMCHAR,
     &                 INTR_ARR,NUMINTR,
     &                 REAL_ARR,DUMMY_INTR,
     &                 DBLE_ARR,DUMMY_INTR,
     &                 LOGL_ARR,DUMMY_INTR)           
C
C---- Error checking
C
      IF( NERR .NE. __io_ok__ ) THEN
        CALL ERRMSG('*B*')
        CALL ERRMSG('REDHDR: reading header failed: what=')
        CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
        CALL ERRMSG(' where=')
        CALL ERRMSG(WHERE(1:STRLN0(WHERE)))
        CALL ERRMSG('*E*')
        GOTO 998
      ENDIF     
C
      IF( LOCAL_WHAT(1:STRLN0(LOCAL_WHAT)) .NE.
     &   WHAT(1:STRLN0(LOCAL_WHAT))) THEN
        CALL ERRMSG('*B*')
        CALL ERRMSG('REDHDR: what in the file is ')
        CALL ERRMSG(LOCAL_WHAT(1:STRLN0(LOCAL_WHAT)))
        CALL ERRMSG('what passed in is ')
        CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
        CALL ERRMSG('*E*')
        NERR = __io_read_err__
        GOTO 998
      ENDIF     
      IF(LOCAL_WHERE(1:STRLN0(LOCAL_WHERE)) .NE.
     &   WHERE(1:STRLN0(LOCAL_WHERE))) THEN
        CALL ERRMSG('*B*')
        CALL ERRMSG('REDHDR: where in the file is ')
        CALL ERRMSG(LOCAL_WHERE(1:STRLN0(LOCAL_WHERE)))
        CALL ERRMSG('what passed in is ')
        CALL ERRMSG(WHAT(1:STRLN0(WHAT)))
        CALL ERRMSG('*E*')
        NERR = __io_read_err__
        GOTO 998
      ENDIF     
      IF(INTR_ARR(1) .NE. LOCAL_WHEN)THEN
        CALL STRI2S(DUMMY,INTR_ARR(1))
        CALL STRI2S(DUMMY2,LOCAL_WHEN)
        CALL ERRMSG('REDHDR: when in the file is '
     &              //DUMMY2(1:STRLN0(DUMMY2))//
     &              'when in the index is '//DUMMY(1:STRLN0(DUMMY)))
        NERR = __io_read_err__
        GOTO 998
      ENDIF     
C
C---- Parse the datatype string CSTRCT to DATTYP AND DATNUM 
C
C      this code is only dealing with one datatype structure,since
C      in 5.2,only single datatype is allowed
C
      IF (NSTRCT .LE. 0) GOTO 999
      IF(CSTRCT(1:4) .EQ. 'REAL') THEN
        DATTYP(1) = __real_data_type__ 
        DATNUM(1) = 1
      ELSE IF(CSTRCT(1:4) .EQ. 'INTR') THEN
        DATTYP(1) = __int_data_type__ 
        DATNUM(1) = 1
      ELSE IF(CSTRCT(1:4) .EQ. 'CHAR') THEN
        DATTYP(1) = __char_data_type__
        DATNUM(1) = 1
      ELSE IF(CSTRCT(1:4) .EQ. 'LOGL') THEN
        DATTYP(1) = __logical_data_type__ 
        DATNUM(1) = 1
      ELSE IF(CSTRCT(1:4) .EQ. 'DBLE') THEN
        DATTYP(1) = __double_data_type__ 
        DATNUM(1) = 1     
      ELSE IF(CSTRCT(1:4) .EQ. 'STRI') THEN
        DATTYP(1) = __string_data_type__ 
        DATNUM(1) = 80
      ELSE IF(CSTRCT(1:4) .EQ. 'STRC') THEN
        DUMMY= CSTRCT(6:STRLN0(CSTRCT))
        INDEX = 1
        DO 400 I=1,IO_MAXENT
          IF(INDEX .GE. STRLN0(DUMMY)) THEN
            GOTO 700
          ENDIF
          START = INDEX
          END = INDEX+3
          IF(DUMMY(START:END) .EQ. 'REAL') THEN
            DATTYP(I) = __real_data_type__
          ELSE IF(DUMMY(START:END) .EQ. 'INTR') THEN
            DATTYP(I) = __int_data_type__
          ELSE IF(DUMMY(START:END) .EQ. 'CHAR') THEN
            DATTYP(I) = __char_data_type__
          ELSE IF(DUMMY(START:END) .EQ. 'LOGL') THEN
            DATTYP(I) = __logical_data_type__
          ELSE IF(DUMMY(START:END).EQ. 'DBLE') THEN
            DATTYP(I) = __double_data_type__
          ELSE IF(DUMMY(START:END) .EQ. 'STRI') THEN
            DATTYP(I) = __string_data_type__
          ENDIF
C
          TMPSTR = DUMMY(END+2:STRLN0(DUMMY))
          CALL STRTOI(DATNUM(I),DUMMY_INTR,TMPSTR)
C
          DUMMY = TMPSTR
          CALL STRBRK(DUMMY,LOCATION,' ')
          IF(INDEX .EQ. 0) THEN
            GOTO 700
          ENDIF
          INDEX = LOCATION+1
 400    CONTINUE
 700    CONTINUE
      ENDIF   
C
 999  CONTINUE
C
C---- Exit
C
      RETURN
 998  CONTINUE
      CALL CFXSTP(ROUTIN)
      END

