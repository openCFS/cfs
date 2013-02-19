*deck,bintwr                                                            sam
      subroutine bintwr (pname,nname)
c *** primary function: bin file copy utility
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
c     pname     (chr,sc,in)       - name of binary file which is to be copied
c
c  output arguments:
c     variable  (typ,siz,intent)    description
c     nname     (chr,sc,out)      - name of new binary file which is a copy
c                                   of pname
c  common variables:
c     iout      (int,sc,comm)     - output unit number
c     intpdp    (int,sc,comm)     - number of integers per double precision word
c     lenfnm    (int,sc,comm)     - number of characters in the filename
c     reclng    (int,sc,comm)     - system record length
c
c                           NOTE: bintwr is not part of binlib.a.  it is 
c                                 included only as an aid to users.
c
#include "impcom.inc"
      external binset,binrd8,binwrt8,binclo,binhed8
      integer  binset
c
c          *****  system record length - same as reclng  *****
      integer  ioleng
      parameter (ioleng=16384)
c
      integer  nunit,nstop,i,ir,n,nunit2,lbuf2,l,npage,
     x kext,lbuf,nbuf,j,kbf,nbuf2
      LONGINT iloc,loc2
c
c          *****  array for reading file data  *****
c            2 is the integers per double precision (integer*4)  *****
c          100000 is the maximum expected ansys logical record size
c          this number may need to be increased for very large models
#if defined (FULL_64)
      integer  iarray(10000000)
#elif defined (LINUX_SYS)
      integer  iarray(100000*2)
#else
      integer  iarray(10000000*2)
#endif
c
c          *****  buffer must be dimensioned at least
c             npage*ioleng*no.of files open
c     integer  buffer(1*ioleng*2)
      integer  buffer(2*ioleng)
c
c          *****  the file name  *****
      character*260  pname,nname
c
c          ******  common for data from bintst  *****
      integer  iout,intpdp,lenfnm,reclng
      common /bintcm/ iout,intpdp,lenfnm,reclng
c
c          *****  define file units *****
      nunit  = 13
      nunit2 = 12
c          *****  define buffer numbers  *****
      nbuf  = 1
      nbuf2 = 2
c          *****  start address in buffer  *****
      lbuf  = 1
      lbuf2 = ioleng + 1
c          *****  number of pages per file buffer  *****
      npage = 1
c          *****  length of each page  *****
      j = reclng
c          *****  external format  *****
      kext = 1
c
c          *****  open read file  *****
      i = binset (nbuf,nunit,1,lbuf,j,npage,pname,lenfnm,kext,
     x     buffer(1))
c
c          *****  check for error  *****
      if (i .ne. 0) then
         write (iout,2001) pname
 2001    format (/' *** ERROR ***'/' DATA FILE ',a,' DOES NOT EXIST')
         go to 999
      endif
c
c          *****  open write file  *****
      l = binset (nbuf2,nunit2,2,lbuf2,j,npage,nname,lenfnm,kext,
     x     buffer(1))
c
c          *****  check for error  *****
      if (l .ne. 0) then
         write (iout,2002) nname
 2002    format (/' *** ERROR ***'/' DATA FILE ',a,' DOES NOT EXIST')
         go to 999
      endif
c
c          ******  define range of records to copy  *****
      nstop = 50000000
c
c          **********  copy the file  **********
c
c          *****  file starting position  *****
      iloc = 0
      loc2 = 0
c          *****  record starting number  *****
      ir = 0
c          *****  next record  *****
 100  ir = ir + 1
      if (ir .le. nstop) then
c
c             *****  read record  *****
         n = 10*intpdp
         call binrd8 (nbuf,iloc,n,iarray(1),kbf,buffer(1))
c
c             *****  check for end-of-file  *****
         if (n .le. 0) go to 200
c
c         *****  rewrite header record  ****
         if (ir .eq. 1)  then
            call binhed8 (nbuf2,nunit2,loc2,buffer(1))
c
c             *****  write record  *****
         else
            call binwrt8 (nbuf2,loc2,n,iarray(1),kbf,buffer(1))
         endif
c
         go to 100
      endif
c
c          ***********  close the files  **********
  200 call binclo (nbuf,'KEEP',buffer(1))
      call binclo (nbuf2,'KEEP',buffer(1))
c
  999 continue
c
      return
      end
