C=======================================================================
C @(#) $Id: OpenFile.F,v 1.6 1999/06/11 14:00:32 cfdirh Exp $
C-----------------------------------------------------------------------
C
C
C
C=======================================================================
C
      SUBROUTINE OpenFile(NERR, FILENAME, WHATFILE)
C
CSCDT Open a file 
CSCDB
C      Inputs:
C        FILENAME : charater string, file name
C        WHATFILE : integer,
C                    __io_close_primaryfile__
C                    __io_close_secondaryfile__
C      Oupputs:
C        NERR   : integer,
C                      0 - no error.
C
CSCDE
C=======================================================================
C
C---------------------------
C Preprocessor includes
C---------------------------
C
#include "cfd_io.h"
C
C---------------------------
C Subroutine Arguments
C---------------------------
C
      INTEGER           NERR, WHATFILE
C
      CHARACTER*(*)     FILENAME 
C
C
C---------------------------
C Local Variables
C---------------------------
C
      INTEGER           DUMMY_NUM, DUMMY_VAL,STRLN0,
     &                  LEN
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
      LEN = STRLN0(FILENAME) 

      CALL IOCNTF_WRAP(WHATFILE,NERR, ERRMSG,
     &            DUMMY_CHAR,DUMMY_CHAR,DUMMY_NUM,
     &            DUMMY_NUM, DUMMY_CHAR,
     &            DUMMY_NUM,DUMMY_CHAR,DUMMY_NUM, 
     &            DUMMY_VAL,DUMMY_NUM,
     &            DUMMY_CHAR, DUMMY_NUM,
     &            FILENAME, LEN,
     &            DUMMY_VAL,DUMMY_NUM,
     &            RARR, DUMMY_NUM,
     &            DARR, DUMMY_NUM,
     &            LARR, DUMMY_NUM)         
      IF( NERR .NE. __io_ok__) THEN
        IF(WHATFILE .EQ. __io_open_primaryfile__) THEN
          CALL ERRMSG('OpenFile: open primary file failed')
        ELSE
          CALL ERRMSG('OpenFile: open secondary file failed')
        ENDIF
      ENDIF
C
C------ Exit.
C
 999  CONTINUE
      RETURN
      END





