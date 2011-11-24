*deck,bintst                                                            sam
      program bintst
c *** primary function: dumps a binary file with the name, file.rst,
c     to the screen. it then takes that file and copies it to a new
c     file, file2.rst. lastly, it dumps the file, file2.rst, to the
c     screen for comparison purposes
c
c *** Notice - This file contains ANSYS Confidential information ***
c
c *** ansys(r) copyright(c) 2000
c *** ansys, inc.
c
c  common variables:
c     iout      (int,sc,comm)     - output unit number
c     intpdp    (int,sc,comm)     - number of integers per double precision word
c     lenfnm    (int,sc,comm)     - number of characters in the filename
c     reclng    (int,sc,comm)     - system record length
c
c                           NOTE: bintst is not part of binlib.a.  it is 
c                                 included only as an aid to users.
c
#include "impcom.inc"
#include "unixpm.inc"

      external binini,bintrd,bintwr,biniqr,bintfo,modinfo,binsplit
      integer  biniqr
c
      integer        units, code, splitpoint

      character*80   title(2)
      character*260  pname,nname
      character*250  jobnam
c
      integer         iout,intpdp,lenfnm,reclng
      common /bintcm/ iout,intpdp,lenfnm,reclng
c
c          **********  define application exit mode   **********        system
      call modinfo(1,1)
c          **********  define the output unit number  **********        system
      iout   = 6
c          **********  define the number of integers per double precision
      intpdp = biniqr (0,2)
c          **********  define the number of characters in the file name
      lenfnm = biniqr (0,3)
c          **********  define the i/0 buffer page length (integer*4 words)
      reclng = biniqr (0,1)
c
c          **********  define the file names  **********                system
      pname = 'file.rst'
      nname = 'file2.rst'
c
c          **********  initialize the bin routines  **********
      call binini (iout)
c
c          **********  initialize the header common  **********
      title(1) = 'New title as given by BINTST'
      title(2) = 'New subtitle as given by BINTST'
      jobnam   = 'file    '
      units    = 0
c ---  the value for code below should be changed to appropriate 3rd party code
      code     = 200
      call bintfo (title,jobnam,units,code)

c          **********  set the file split point  **********
      splitpoint = verybig      ! split at 32768 MW
      call binsplit (splitpoint)

c          **********  read and list the file file.rst  **********
      call bintrd (pname)
c
c          **********  copy file file.rst to file file2.rst  **********
c      call bintwr (pname,nname)
c
c          **********  read and list the file file2.rst  **********
c      call bintrd (nname)
c
      stop
      end
