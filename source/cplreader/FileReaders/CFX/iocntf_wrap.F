      SUBROUTINE IOCNTF_WRAP(REQ,NERR,ERRMSG,
     &                  WHAT,WHERE,WHEN, NATTRIB,ATTRIB,
     &                  NSTRCT,STRCT,IBLOCK,
     &                  OFFSET,NOFFSET,
     &                  STRING_ARR,NSTRING,
     &                  STRING,NCHAR,
     &                  INTEGER_ARR,NINTEGER,
     &                  REAL_ARR,NREAL,
     &                  DBLE_ARR,NDBLE,
     &                  LOGL_ARR,NLOGL)
C
C     Converts *8 real arguments to *4 for calling IOCNTF when flow solver is
C     in double precision on DEC.  
C
C=======================================================================
C
C---------------------------
C Subroutine Arguments
C---------------------------
C
      INTEGER           REQ,NERR,WHEN,NATTRIB,NSTRCT,IBLOCK,
     &                  OFFSET(*),NOFFSET,
     &                  NSTRING,NCHAR,NINTEGER,NREAL,NDBLE,NLOGL,
     &                  INTEGER_ARR(*)
C
      REAL              REAL_ARR(*)
C
      CHARACTER*(*)     WHAT,WHERE,ATTRIB(*),STRCT,STRING_ARR(*),STRING
C
      CHARACTER*1       LOGL_ARR(*)
C
      DOUBLE PRECISION  DBLE_ARR(*)
C
C---------------------------
C Local Variables
C---------------------------
C
      EXTERNAL ERRMSG
C
C---------------------------
C Executable Statements
C---------------------------
C
C---  Convert word length to *4 from *8.
C        Set to zero on read as data array passed in may be unset
#ifdef DOUBLE_PRECISION
      IF (NREAL.GT.0) THEN
         IF (REQ.EQ.14.OR.REQ.EQ.6) THEN
            CALL R8_TO_R4(REAL_ARR,REAL_ARR,NREAL)        
         ELSE
            CALL ZERO_R4(REAL_ARR,NREAL)
         ENDIF
      ENDIF
#endif
C

      CALL IOCNTF(REQ,NERR, ERRMSG,
     &            WHAT, WHERE, WHEN, NATTRIB, ATTRIB,
     &            NSTRCT,STRCT,IBLOCK,
     &            OFFSET,NOFFSET,
     &            STRING_ARR, NSTRING,
     &            STRING, NCHAR,
     &            INTEGER_ARR,NINTEGER,
     &            REAL_ARR, NREAL,
     &            DBLE_ARR, NDBLE,
     &            LOGL_ARR, NLOGL)
C
C---  Convert word length to *8 from *4.
C
#ifdef DOUBLE_PRECISION
      IF (NREAL.GT.0) CALL R4_TO_R8(REAL_ARR,REAL_ARR,NREAL)
#endif
C
      END
