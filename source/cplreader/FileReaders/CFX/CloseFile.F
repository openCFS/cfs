C=======================================================================
C @(#) $Id: CloseFile.F,v 1.6 1999/06/11 14:00:30 cfdirh Exp $
C-----------------------------------------------------------------------
C
C
C
C=======================================================================
C
      SUBROUTINE CloseFile(NERR, WHATFILE)
C
CSCDT close a file 
CSCDB
C      Inputs:
C        WHATFILE : integer, __io_close_primaryfile__ or 
C                            __io_close_secondaryfile__
C      Oupputs:
C        NERR   : integer,
C                      0 - no error.
C
CSCDE
C=======================================================================
C---------------------------
C Preprocessor includes
C---------------------------
C
#include "cfd_io.h"
C
C
C---------------------------
C Subroutine Arguments
C---------------------------
C
      INTEGER           NERR, WHATFILE
C
C
C---------------------------
C Local Variables
C---------------------------
C
      INTEGER           DUMMY_NUM, DUMMY_VAL
C
      REAL              RARR(1)
C
      CHARACTER*1       LARR(1)
      CHARACTER*(80)     DUMMY_CHAR
C
      DOUBLE PRECISION  DARR(1)
C
      EXTERNAL ERRMSG
C
C---------------------------
C Executable Statements
C---------------------------
C
C------ Initialization.
C
      NERR   = 0
      DUMMY_CHAR = ' '
      DUMMY_NUM = 0
      DUMMY_VAL = 0 
      CALL IOCNTF_WRAP(WHATFILE,NERR, ERRMSG,
     &             DUMMY_CHAR,DUMMY_CHAR,DUMMY_NUM,
     &             DUMMY_NUM, DUMMY_CHAR,
     &             DUMMY_NUM,DUMMY_CHAR,DUMMY_NUM, 
     &             DUMMY_VAL,DUMMY_NUM,
     &             DUMMY_CHAR, DUMMY_NUM,
     &             DUMMY_CHAR, DUMMY_NUM,
     &             DUMMY_VAL,DUMMY_NUM,
     &             RARR, DUMMY_NUM,
     &             DARR, DUMMY_NUM,
     &             LARR, DUMMY_NUM)      
      IF( NERR .NE. __io_ok__) THEN
        IF(WHATFILE .EQ. __io_close_primaryfile__) THEN
          CALL ERRMSG('CloseFile: close primary file failed')
        ELSE
          CALL ERRMSG('CloseFile: close secondary file failed')
        ENDIF
      ENDIF
C
C------ Exit.
C
 999  CONTINUE
      RETURN
      END





