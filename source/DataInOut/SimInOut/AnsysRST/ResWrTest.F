*deck,ResWrTest                  User
      function ResWrTest (FName, ncFName)
c  **********  Function writes a mini .rst file (just header
c              and extracts some infos from it  *****
#include "impcom.inc"

      external  ResWrBegin    , ResWrEnd
      integer   ResWrTest, ResWrBegin

      character*(*)  Fname
      integer ncFName

      integer  Nunit, Lunit, Units, UserCode, MaxResultSet
      integer  kerr 
      character*8   JobName

      data  Nunit,Lunit,Units,UserCode,MaxResultSet
     x     /  12,   6,    4,    10,        10      /
c      data  Fname /'ResWr.rst'/

      integer     MAXNODE,    MAXELEM,   NUMNODE,   NUMELEM, NUMDOF
      parameter  (MAXNODE=13, MAXELEM=3, NUMNODE=6, NUMELEM=2, NUMDOF=2)
      integer    KAN
      integer    DOF(NUMDOF)

      character*80  Title(5)

      data  Title / ' 2 element cantilever beam using plane42',
     x              'Test case for writing binary files',
     x              3*'       '/

      JobName  = 'ResWTest'

      kerr = ResWrBegin (Nunit,Lunit,Fname,ncFName,Title(1),JobName,
     x                   Units,NUMDOF,DOF(1),UserCode,
     x                   MAXNODE,NUMNODE,MAXELEM,NUMELEM,
     x                   MaxResultSet, KAN)

      call ResWrEnd

      if (kerr .ne. 0) then
         ResWrTest = 1
         return
      endif

      end
