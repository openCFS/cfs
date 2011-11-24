*deck,bintrd                                                            sam
      subroutine bintrd (pname)
c *** primary function: bin file dump utility
c
c *** Notice - This file contains ANSYS Confidential information ***
c
c *** ansys(r) copyright(c) 2000
c *** ansys, inc.
c
c     typ=int,dp,log,chr,dcp   siz=sc,ar(n)   intent=in,out,inout
c
c  input arguments:
c     variable  (typ,siz,intent)    description
c     pname     (chr,sc,in)       - name of binary file which is to
c                                   be dumped to the screen
c
c  output arguments:
c     none.
c
c  common variables:
c     iout      (int,sc,comm)     - output unit number
c     intpdp    (int,sc,comm)     - number of integers per double precision word
c     lenfnm    (int,sc,comm)     - number of characters in the filename
c     reclng    (int,sc,comm)     - system record length
c
c                           NOTE: bintrd is not part of binlib.a.  it is 
c                                 included only as an aid to users.
c
#include "impcom.inc"
      external binset,binrd8,binclo
      integer  binset
c
c          *****  system record length - same as reclng  *****
      integer  ioleng
      parameter (ioleng=16384)
c
      integer  nunit,nstart,nstop,i,ir,n,k,
     x kext,lbuf,nbuf,j,kbf,npage,keyrw
      LONGINT locsav,iloc
c
c          *****  arrays for reading file data  *****
c            2 is the integers per double precision (intpdp)
c          100000 is the maximum expected ansys logical record size
c          this number may need to be increased for very large models
#if defined (FULL_64)
      integer  iarray(10000000)
#elif defined (LINUX_SYS)
      integer  iarray(100000*2)
#else
      integer  iarray(10000000*2)
#endif
#if defined (LINUX_SYS)
      double precision  darray(100000)
#else
      double precision  darray(10000000)
#endif
c
c          *****  buffer must be dimensioned at least
c             npage*ioleng*no.of files open
c     integer  buffer(1*ioleng*1)
      integer  buffer(ioleng)
c
c          *****  the file name  *****
      character*260  pname
c
c          *****  common for data from bintst  *****
      integer  iout,intpdp,lenfnm,reclng
      common /bintcm/ iout,intpdp,lenfnm,reclng
c
c          *****  make equivalent for listing integer or dp  *****
      equivalence (darray(1),iarray(1))
c
c
      write (iout,2001)
 2001 format (/'   ***** ANSYS BINARY FILE DUMP *****'//)
c
      write (iout,2003)
 2003 format (/' WRITE OUT ALL DATA ITEMS IN EACH RECORD')
c
c          **********  open the file  **********
c
c          *****  set the i/o unit number
      nunit = 12
c          *****  define buffer number *****
      nbuf = 1
c          *****  start address in the buffer  *****
      lbuf = 1
c          *****  number of pages for nbuf  *****
      npage = 1
c          *****  read key  *****
      keyrw = 1
c          *****  length of each page  *****
      j = reclng
c          *****  read external file format  *****
      kext = 1
c
      i = binset (nbuf,nunit,keyrw,lbuf,j,npage,pname,lenfnm,kext,
     x buffer(1))
c
c          *****  check for error  *****
      if (i .ne. 0) then
         write (iout,2004) pname
 2004    format (/' *** ERROR ***'/' DATA FILE ',a,' DOES NOT EXIST ')
         go to 999
      endif
c
c          ********** define the range of records to list  *********
      nstart = 2
      nstop = 50000000
      write (iout,2005)  pname,nstart,nstop
 2005 format (/' DUMP OF FILE= ',a/'   RECORDS',i8,' TO',i8/)
c
c          **********  list the file  **********
c
c          *****  file starting position  *****
      iloc = 0
c          *****  record starting number  *****
      ir = 0
c          *****  next record  *****
 100  ir = ir + 1
      if (ir .le. nstop) then
         locsav = iloc
c              *****  read a record  *****
         n = 10*intpdp
         call binrd8 (nbuf,iloc,n,iarray(1),kbf,buffer(1))
         if (n .le. 0) then
            write (iout,2006) pname
 2006       format (/' END-OF-FILE ON FILE= ',a//)
            go to 200
         endif
c             *****  skip if below first requested record  *****
         if (ir .lt. nstart) go to 100
c
         if (kbf .eq. 0) n = n/intpdp
         write (iout,2007) ir,n,locsav,kbf
 2007    format (/' RECORD=',i8,'   LENGTH=',i8,'   LOC=',i12,
     x     '   KBF=',i8)
            k = n
            if (kbf .eq. 0) then
               write (iout,2008) (darray(j),j=1,k)
 2008          FORMAT (1X,5G14.6)
            else
               write (iout,2009) (iarray(j),j=1,k)
 2009          FORMAT (1X,10I7)
            endif
         go to 100
      endif
c
c          **********  close the file  **********
  200 call binclo (nbuf,'KEEP',buffer(1))
c
  999 continue
c
      return
      end
