      SUBROUTINE CFXSTP(TEXT)
      CHARACTER*(*) TEXT
      PRINT*,TEXT
      END

      SUBROUTINE ERRMSG(TEXT)
      CHARACTER*(*) TEXT
      PRINT*,TEXT
      END

C=======================================================================
C @(#) $Id: stri2s.f,v 1.4 1999/02/13 11:06:21 cfdpsh Exp $
C-----------------------------------------------------------------------
C $Revision: 1.4 $
C
C
C
C=======================================================================
C
      SUBROUTINE STRI2S(STRING, N)
C
CSCDT Converts integer N to STRING.
C
CSCDB
C
CSCDE
C
C=======================================================================
C  
C----------------------
C Change Log
C----------------------
C  $Author: cfdpsh $
C  $Date: 1999/02/13 11:06:21 $
C  $Log: stri2s.f,v $
C  Revision 1.4  1999/02/13 11:06:21  cfdpsh
C  Correctly add the CVS Id tag to the top of the fortran files.
C
C  Revision 1.3  1999/02/11 18:55:52  cfdpsh
C  Imported from restructured solver project
C
C  Revision 1.1.1.1  1998/07/09 10:16:40  cfdlib
C  Imported source from CFX-5.2 (GE)
C
C  Revision 1.1  1998/04/24 13:01:16  rcs
C  Initial revision
C
C  Revision 1.1  1997/07/28 17:58:08  rcs
C  Initial revision
C
C  Revision 1.1  1997/03/26 15:34:04  wilson
C  Initial revision
C
C  
C     Rev 1.0   06/06/95 17:21:24   wilson
C  Initial revision.
C  
C     Rev 1.0   06/06/95 16:49:02   wilson
C  Initial revision.
C  
C------------------------
C Subroutine Arguments
C------------------------
C
      INTEGER       N
C
      CHARACTER*(*) STRING
C
C--------------------------
C Local Variables
C--------------------------
C
      INTEGER IBASE, M, I, J
C
      CHARACTER*15 NUMSTR
C
      LOGICAL NEGFLG
C
C---------------------------
C Executable Statements
C---------------------------
C
      IBASE = ICHAR('1')-1
C      L = LEN(STRING)
      STRING = ' '
C
      M = N
      I = 15
      IF (M.LT.0) THEN
        M = -M
        NEGFLG = .TRUE.
C
C  Test for case of the most negative representable number -- it may
C  not be handled properly by some machines' arithmetic in that M=-M !!
C  Bodge it by reducing its magnitude so we can take its negative.
C  NB Chances are the machine has simply rounded an even larger number
C  towards zero, so we aren't losing anything by doing this.
        IF (M.LT.0) M = -(M+1)
      ELSE
        NEGFLG = .FALSE.
      ENDIF
C
 600  J = M - (M/10)*10
      M = M/10
      IF(J.EQ.0) THEN
        NUMSTR(I:I) = '0'
      ELSE
        NUMSTR(I:I) = CHAR(IBASE+J)
      ENDIF
      I = I-1
      IF(M.NE.0) GOTO 600
      IF(NEGFLG) THEN
        NUMSTR(I:I) = '-'
        I = I-1
      ENDIF
      STRING = NUMSTR(I+1:15)
C
C------ Exit.
C
      RETURN
      END
C=======================================================================
C @(#) $Id: strln0.f,v 1.4 1999/02/13 11:06:24 cfdpsh Exp $
C-----------------------------------------------------------------------
C $Revision: 1.4 $
C
C
C
C=======================================================================
C
      INTEGER FUNCTION STRLN0(STRING)
C
CSCDT Returns length of string less trailing blanks.
C
CSCDB 
C      Inputs: 
C        STRING : character, string to return length of
C
C      Outputs:
C        STRLN0 : integer, length of STRING, minus trailing blanks
C
C      Local:
C        L      : integer, declared length of STRING
C        I      : integer, loop variable
CSCDE
C
C=======================================================================
C  
C---------------------------
C Function Arguments
C---------------------------
C
      CHARACTER*(*) STRING
C
C---------------------------
C Local Variables
C---------------------------
C
      INTEGER       L, I
C
C---------------------------
C Executable Statements
C---------------------------
C
      L = LEN(STRING)
C
      DO 100 I = L,1,-1
        IF( STRING(I:I) .NE. ' ' ) THEN
           STRLN0 = I
           GOTO 999
        ENDIF
  100 CONTINUE
C
      STRLN0 = 0
C
  999 CONTINUE
      RETURN
      END               
C=======================================================================
C @(#) $Id: strtoi.f,v 1.4 1999/02/13 11:06:30 cfdpsh Exp $
C-----------------------------------------------------------------------
C $Revision: 1.4 $
C
C
C                
C=======================================================================
C
      SUBROUTINE STRTOI(N, NI, STRING)
C                             
CSCDT Converts a string to an integer N.  
CSCDB 
C     NI is index to next character after number.
CSCDE
C
C=======================================================================
C  
      INTEGER         N, NI
      CHARACTER*(*)   STRING 
C
      INTEGER    I 
      LOGICAL    NEGFLG
C
      NEGFLG = .FALSE.                          
      CALL FLUSHB(I, STRING)
C
      IF (STRING(I:I).EQ.'-') THEN
        NEGFLG = .TRUE.  
        I=I+1
      ENDIF
C
      CALL STRTOW(N, NI, STRING(I:))
      NI = NI+I-1
      IF (NEGFLG) N = N*(-1)
C
      RETURN
      END
C=======================================================================
C @(#) $Id: strbrk.f,v 1.4 1999/02/13 11:06:18 cfdpsh Exp $
C-----------------------------------------------------------------------
C $Revision: 1.4 $
C
C
C
C=======================================================================
C
      SUBROUTINE STRBRK(STRING,INDEX,VALID)
C
CSCDT Return index of 1st valid character in STRING.
C
C=======================================================================
C  
      CHARACTER*(*) STRING,VALID
C
      INTEGER       NVALID, LSTR, I, J, INDEX
C
      NVALID=LEN(VALID)
      LSTR=LEN(STRING)
C
      DO 100 I=1,LSTR
        DO 100 J=1,NVALID
          IF (STRING(I:I).EQ.VALID(J:J)) THEN 
            INDEX = I
            GOTO 999
          ENDIF
 100  CONTINUE
C
      INDEX=0
C
  999 CONTINUE
C
      RETURN
      END
C=======================================================================
C @(#) $Id: strctn.f,v 1.4 1999/02/13 11:06:19 cfdpsh Exp $
C-----------------------------------------------------------------------
C $Revision: 1.4 $
C
C
C
C=======================================================================
C
      SUBROUTINE STRCTN( STR1,STR2 )
C
CSCDT Concatenate strings
CSCDB
C   Concatenates STR2 after the last non-blank character in STR1.
CSCDE
C=======================================================================
C  
C---------------------------
C Subroutine Arguments
C---------------------------
C
      CHARACTER*(*) STR1,STR2
C
C---------------------------
C Local Variables
C---------------------------
C
      INTEGER L, STRLN0
C
C---------------------------
C Executable Statements
C---------------------------
C
      L = STRLN0(STR1)
      IF (L .LT. LEN(STR1)) STR1(L+1:) = STR2
C
      RETURN
      END
      SUBROUTINE SETARI(IARR,IVALUE,NLEN)
C
C*****************************
C
C @(#) setari.f    VERSION 5.1    10/8/96 (mm/dd/yy)
C
C*****************************
C***********************************************************************
C
C  cfx maths module.
C
C  entry point subroutine to set integer array iarr(nlen) to a
C  constant ivalue.
C
C***********************************************************************
C
C   alteration record
C   -----------------
C
C   created
C     16/02/94 adb
C   modified
C
C***********************************************************************
C
C  declarations:
C  -------------
C
      INTEGER IARR(NLEN),IVALUE,NLEN,I
C
C***********************************************************************
C
      DO 100 I = 1, NLEN
         IARR(I) = IVALUE
  100 CONTINUE
C
      RETURN
      END
C=======================================================================
C @(#) $Id: flushb.f,v 1.4 1999/02/13 11:04:43 cfdpsh Exp $
C-----------------------------------------------------------------------
C $Revision: 1.4 $
C
C
C
C=======================================================================
C
      SUBROUTINE FLUSHB(J,STRING)
C
CSCDT Return index of first non-blank, non-TAB character in STRING.
CSCDE
C      Inputs:
C        STRING : input string
C
C      Outputs:
C        J      : index of first non-white space character.
C                 If none found, point to last character 
C                 in string.
C
C=======================================================================
C  
      CHARACTER*(*) STRING
C
      INTEGER       I, J, L
C
      L = LEN(STRING)
      DO 100 I = 1,L
        IF(STRING(I:I).NE.' ') THEN
C                                              Check for TAB (ASCII 011)
          IF(ICHAR(STRING(I:I)).NE.9) THEN
            J = I
            GOTO 999
          ENDIF
        ENDIF
  100 CONTINUE
C
      J = L
C
  999 CONTINUE
C
      RETURN
      END
C=======================================================================
C @(#) $Id: strtow.f,v 1.4 1999/02/13 11:06:32 cfdpsh Exp $
C-----------------------------------------------------------------------
C $Revision: 1.4 $
C
C
C
C=======================================================================
C
      SUBROUTINE STRTOW(N, NI, STRING)
C
CSCDT Convert a string to a whole number N.
CSCDB                             
C   Converts a string to a whole number N.  NI is index to next
C   character C after number. Allows leading '+' sign.
CSCDE
C
C=======================================================================
C  
      INTEGER         N, NI
      CHARACTER*(*)   STRING 
C
      INTEGER         LENSTR, I, J, L
      LOGICAL         ISNUM
C
      LENSTR = LEN(STRING)
      CALL FLUSHB(I,STRING)
C
C --- Identify Digit string.
C
      IF (STRING(I:I).EQ.'+') I=I+1
      J = I
 100  CONTINUE
      IF (ISNUM(STRING(J:J))) THEN
        J = J+1
        IF (J.LE.LENSTR) GOTO 100
      ENDIF
C
C --- Calculate value of number.
C
      N = 0
      DO 200 L = I,J-1
        N = N*10+ICHAR(STRING(L:L))-ICHAR('0')
 200  CONTINUE
C
 900  CONTINUE
      NI = J
C
      RETURN
      END
C=======================================================================
C @(#) $Id: isnum.f,v 1.4 1999/02/13 11:05:22 cfdpsh Exp $
C-----------------------------------------------------------------------
C $Revision: 1.4 $
C
C
C
C=======================================================================
C
      LOGICAL FUNCTION ISNUM(CHAR)
C
CSCDT Returns TRUE if char is numeric. '0'-'9'
C
C=======================================================================
C  
      CHARACTER*1  CHAR
C
      ISNUM = (CHAR.GE.'1'.AND.CHAR.LE.'9').OR.
     &        (CHAR.EQ.'0')
C
      RETURN
      END
      INTEGER FUNCTION LENACT(CSTR)
C
C Prototype.
C
C      Returns active length of a character string CSTR.
C
C
C***********************************************************************
C
C @(#) lenact.f    VERSION 5.1    10/8/96 (mm/dd/yy)
C
C***********************************************************************
C
      IMPLICIT REAL (A-H , O-Z)
C
C***********************************************************************
C
C Argument list.
C **************
      CHARACTER CSTR*(*)
C******************************************************************************
C Local variables.
C ****************
      INTEGER MAXLEN,I,ISPA
C******************************************************************************
      ISPA = ICHAR(' ')
      MAXLEN = LEN(CSTR)
C
      LENACT = 0
      DO 1 I=MAXLEN,1,-1
         IF ( ICHAR(CSTR(I:I)) .NE. ISPA ) THEN
            LENACT = I
            RETURN
         ENDIF
1     CONTINUE
C
      RETURN
      END
