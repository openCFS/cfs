*deck,ResWr                        User
c   **********   ANSYS RESULT FILE WRITING OBJECT   **********
c   Referenced libraries:  binlib.a
c   Include File:  ResWr.inc

c    This pair of files is to assist persons who would like to write a
c    result file which can be processed by the ANSYS program.


*deck,ResWrBegin
      function ResWrBegin (Nunit, Lunit, Fname, ncFname, Title, JobName,
     x                     Units, NumDOF, DOF, UserCode,
     x                     MaxNode, NumNode, MaxElem, NumElem,
     x                     MaxResultSet)
c primary function:    Open result file

c object/library: ResWr

c  input arguments:
c     Nunit    (int,sc,in)       - Fortran Unit number for file (ANSYS uses 12)
c     Lunit    (int,sc,in)       - Current print output unit (usually 6)
c     Fname    (ch*(ncFname),sc,in) - The name (with extension) for the file
c     ncFname  (int,sc,in)       - Number of characters in Fname
c     Title    (ch*80,ar(2),in)  - Title and First subtitle
c     JobName  (ch*8,sc,in)      - Jobname for file
c     Units    (int,sc,in)       - 0, unknown  1, SI  2, CSG 
c                                  3, U.S. Customary - foot
c                                  4, U.S. Customary - inch
c                                  6, MPA
c     NumDOF   (int,sc,in)       - Number of DOF per node
c     DOF      (int,ar(*),in)    - The DOFs per node
c     UserCode (int,sc,in)       - Code for this application (user defined for
c                                  now)
c     MaxNode  (int,sc,in)       - Maximum node number used
c     NumNode  (int,sc,in)       - Number of nodes used
c     MaxElem  (int,sc,in)       - Maximum element number used
c     NumElem  (int,sc,in)       - Number of elements used
c     MaxResultSet (int,sc,in)   - Maximum number of result sets (usually 1000)

c  output arguments:
c     ResWrBegin (int,sc,out)    - 0, successful  other, error in file open

#include "impcom.inc"
#include "ResWr.inc"

      external    fMallocI, fFree, fMallocD
      PTRFTN      fMallocI, fMallocD
      external    biniqr8, binini, binset, bintfo, binhed8, binwrt8,
     x            largeIntPut, binBasePut
      integer     binset
      LONGINT     biniqr8

      LONGINT     LongZeroL
      integer     ResWrBegin, Nunit, Lunit, ncFname, Units,
     x            NumDOF, DOF(NumDOF), UserCode,
     x            MaxNode, NumNode, MaxElem, NumElem, MaxResultSet, i
      character*8 JobName
      character*(*)  Fname
      character*80  Title(2)

c          *****  store arguments  *****
      iout   = Lunit
      nNode  = NumNode
      mxNode = MaxNode
      nElem  = NumElem
      mxElem = MaxElem
      nDOF   = NumDOF
      do  i = 1,nDOF
         curdof(i) = DOF(i)
      enddo
      MaxResult = MaxResultSet

c          ***** initialize binary file I/O system  *****
      call binini (iout)
      intpdp = biniqr8(0,2)
      intplong = biniqr8(0,5)
#if defined (FUJITSU_SYS)
      piIobufL = fMallocI(RECLENG*1_8)
#else
      piIobufL = fMallocI(RECLENG)
#endif


c          *****  open result file  *****
      ResWrBegin = binset(NBLK,Nunit,2,1,RECLENG,1,Fname,ncFname,
     x             1,iobuf(1))
      if (ResWrBegin .ne. 0) then
         write (iout,2000) Fname(1:ncFname)
 2000    format (/' *** ERROR ***'/' Could not open file= ',a,
     x    ' for writing.'/'   Check if the file already exists.',
     x    ' If so, delete it and re-run.')
         go to 999
      endif

      LongZeroL = 0
      call binBasePut (NBLK,LongZeroL)

c          *****  set up standard file header  *****
      call bintfo (Title(1),JobName,Units,UserCode)
      call binhed8 (NBLK,Nunit,iLocL,iobuf(1))

c          *****  set up results file header  *****
      do  i = 1,40
         Header(i) = 0
      enddo
      Header(1) = 12
      Header(2) = MaxNode
      Header(3) = nNode
      Header(4) = MaxResult
      Header(5) = nDOF
      Header(6) = MaxElem
      Header(7) = nElem
      Header(20)= Units
      Header(21)= 1
      LocHeaderL = iLocL
      call binwrt8 (NBLK,iLocL,40,Header(1),40,iobuf(1))

c          *****  write DOF list  *****
      call binwrt8 (NBLK,iLocL,nDOF,DOF(1),nDOF,iobuf(1))

c          *****  write Dummy Node List  *****
#if defined (FUJITSU_SYS)
      piEntitiesL = fMallocI(nNode*1_8)
#else
      piEntitiesL = fMallocI(nNode)
#endif
      LocNodeListL = iLocL
      Header(15)  = iLocL
      do  i = 1,nNode
         Entities(i) = i
      enddo
      call binwrt8 (NBLK,iLocL,nNode,Entities(1),nNode,iobuf(1))
      call fFree (piEntitiesL)


c          *****  write Dummy Element List  *****
#if defined (FUJITSU_SYS)
      piEntitiesL = fMallocI(nElem*1_8)
#else
      piEntitiesL = fMallocI(nElem)
#endif
      LocElemListL = iLocL
      Header(14)  = iLocL
      do  i = 1,nElem
         Entities(i) = i
      enddo
      call binwrt8 (NBLK,iLocL,nElem,Entities(1),nElem,iobuf(1))
      call fFree (piEntitiesL)

c          *****  initialize result indexes  *****
#if defined (FUJITSU_SYS)
      piResIndexL = fMallocI(2*MaxResult*1_8)
      piResLstepL = fMallocI(3*MaxResult*1_8)
      pdResTimeL  = fMallocD(MaxResult*1_8)
#else
      piResIndexL = fMallocI(2*MaxResult)
      piResLstepL = fMallocI(3*MaxResult)
      pdResTimeL  = fMallocD(MaxResult)
#endif

      nSets = 0
      LongZeroL = 0
      do  i = 1,MaxResult
         call largeIntPut (LongZeroL,iResIndex(i),
     x                              iResIndex(i+MaxResult))
         iResLstep(1,i) = 0
         iResLstep(2,i) = 0
         iResLstep(3,i) = 0
         dResTime(i)    = 0.0d0
      enddo

c          *****  write result indexes  *****
      LocResIndexL = iLocL
      Header(11) = iLocL
      i = 2*MaxResult
      call binwrt8 (NBLK,iLocL,i,iResIndex(1),i,iobuf(1))

      Header(12) = iLocL
      call binwrt8 (NBLK,iLocL,MaxResult*intpdp,dResTime(1),0,iobuf(1))

      Header(13) = iLocL
      i = 3*MaxResult
      call binwrt8 (NBLK,iLocL,i,iResLstep(1,1),i,iobuf(1))

 999  return
      end



*deck,ResWrGeomBegin
      subroutine ResWrGeomBegin (MaxType, NumType, MaxReal, NumReal,
     x                           MaxCsys, NumCsys)
c primary function:    Write Geometry Header Record

c object/library: ResWr

c  input arguments:
c     MaxType  (int,sc,in)       - Maximum element type
c     NumType  (int,sc,in)       - Number of element types
c     MaxReal  (int,sc,in)       - Maximum real constant set number
c     NumReal  (int,sc,in)       - Number of real constant sets
c     MaxCsys  (int,sc,in)       - Maximum coordinate system number
c     NumCsys  (int,sc,in)       - Number of user coordinate systems

c  output arguments:  none

#include "impcom.inc"
#include "echprm.inc"
#include "ResWr.inc"

      integer   MaxType, NumType, MaxReal, NumReal, MaxCsys, NumCsys, i

      mxType = MaxType
      nType  = NumType
      mxReal = MaxReal
      nReal  = NumReal
      mxCsys = MaxCsys
      nCsys  = NumCsys

c          **********  write out Geometry Header  *********
      do  i = 1,40
         GeomHeader(i) = 0
      enddo

      GeomHeader(2) = mxType
      GeomHeader(3) = mxReal
      GeomHeader(4) = nNode
      GeomHeader(5) = nElem
      GeomHeader(6) = mxCsys
      GeomHeader(17) = 24
      GeomHeader(19) = IELCSZ

      LocGeomHeaderL = iLocL
      Header(16)    = iLocL
      call binwrt8 (NBLK,iLocL,40,GeomHeader(1),40,iobuf(1))

      return
      end


*deck,ResWrTypeBegin
      subroutine ResWrTypeBegin
c primary function:    Start element type records

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8, largeIntPut, fMallocI
      PTRFTN    fMallocI

      integer   i


      if (mxType .le. 0) go to 999

#if defined (FUJITSU_SYS)
      piIndexToDataL = fMallocI(mxType*1_8)
#else
      piIndexToDataL = fMallocI(mxType)
#endif


      do  i = 1,mxType
         IndexToData(i) = 0
      enddo

c          *****  write out dummy index  *****
      LocIndexL = iLocL
      GeomHeader(7) = iLocL
      call largeIntPut (iLocL,GeomHeader(21),GeomHeader(22))
      call binwrt8 (NBLK,iLocL,mxType,IndexToData(1),mxType,iobuf(1))

 999  return
      end


*deck,ResWrType
      subroutine ResWrType (itype,ielc)
c primary function:    Write an element type record

c object/library: ResWr

c  input arguments:
c     itype    (int,sc,in)         - Element type number
c     ielc     (int,ar(IELCSZ),in) - Element characteristics

c  output arguments:  none


#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8

#include "echprm.inc"

      integer   itype, ielc(IELCSZ)

      IndexToData(itype) = iLocL - LocIndexL
      call binwrt8 (NBLK,iLocL,IELCSZ,ielc(1),IELCSZ,iobuf(1))

      return
      end


*deck,ResWrTypeEnd
      subroutine ResWrTypeEnd
c primary function:    Write out element type pointers

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external   binwrt8, fFree

      if (mxType .gt. 0) then
         call binwrt8 (NBLK,LocIndexL,mxType,IndexToData(1),mxType,
     x    iobuf(1))
         call fFree (piIndexToDataL)
      endif

      return
      end


*deck,ResWrRealBegin
      subroutine ResWrRealBegin
c primary function:    Start real constant records

c object/library: ResWr


#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8, fMallocI
      PTRFTN    fMallocI

      integer   i

      maxPerReal = 0

      if (mxReal .le. 0) go to 999

#if defined (FUJITSU_SYS)
      piIndexToDataL = fMallocI(mxReal*1_8)
#else
      piIndexToDataL = fMallocI(mxReal)
#endif


      do  i = 1,mxReal
         IndexToData(i) = 0
      enddo

c          *****  write out dummy IndexToData  *****
      LocIndexL = iLocL
      GeomHeader(8) = iLocL
      call largeIntPut (iLocL,GeomHeader(23),GeomHeader(24))
      call binwrt8 (NBLK,iLocL,mxReal,IndexToData(1),mxReal,iobuf(1))

 999  return
      end


*deck,ResWrReal
      subroutine ResWrReal (iReal,n,Rcon)
c primary function:    Write an real constant record

c object/library: ResWr

c  input arguments:
c     iReal    (int,sc,in)       - Real set number
c     n        (int,sc,in)       - Number of real constants in set
c     Rcon     (dp,ar(n),in)     - Real Constants

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8

      integer   iReal, n
      double precision Rcon(n)


      if (n .gt. maxPerReal) maxPerReal = n

      IndexToData(iReal) = iLocL - LocIndexL
      call binwrt8 (NBLK,iLocL,n*intpdp,Rcon(1),0,iobuf(1))

      return
      end


*deck,ResWrRealEnd
      subroutine ResWrRealEnd
c primary function:    Write out real constant pointers

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external   binwrt8, fFree

      if (mxReal .gt. 0) then
         GeomHeader(20) = maxPerReal
         call binwrt8 (NBLK,LocIndexL,mxReal,IndexToData(1),mxReal,
     x    iobuf(1))
         call fFree (piIndexToDataL)
      endif

      return
      end


*deck,ResWrCsysBegin
      subroutine ResWrCsysBegin
c primary function:    Start coordinate system records

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8, largeIntPut, fMallocI
      PTRFTN    fMallocI

      integer   i


      if (mxCsys .le. 0) go to 999

#if defined (FUJITSU_SYS)
      piIndexToDataL = fMallocI(mxCsys*1_8)
#else
      piIndexToDataL = fMallocI(mxCsys)
#endif


      do  i = 1,mxCsys
         IndexToData(i) = 0
      enddo

c          *****  write out dummy IndexToData  *****
      LocIndexL = iLocL
      GeomHeader(10) = iLocL
      call largeIntPut (iLocL,GeomHeader(25),GeomHeader(26))
      call binwrt8 (NBLK,iLocL,mxCsys,IndexToData(1),mxCsys,iobuf(1))

 999  return
      end


*deck,ResWrCsys
      subroutine ResWrCsys (iCsys,Csys)
c primary function:    Write an coordinate system record

c object/library: ResWr

c  input arguments:
c     iCsys    (int,sc,in)       - Coordinate system number
c     Csys     (dp,ar(24),in)    - Coordinate system description

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8

      integer   iCsys
      double precision  Csys(24)

      IndexToData(iCsys) = iLocL - LocIndexL
      call binwrt8 (NBLK,iLocL,24*intpdp,Csys(1),0,iobuf(1))

      return
      end


*deck,ResWrCsysEnd
      subroutine ResWrCsysEnd
c primary function:    Write out coordinate system pointers

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external   binwrt8, fFree

      if (mxCsys .gt. 0) then
         call binwrt8 (NBLK,LocIndexL,mxCsys,IndexToData(1),mxCsys,
     x    iobuf(1))
         call fFree (piIndexToDataL)
      endif

      return
      end


*deck,ResWrNodeBegin
      subroutine ResWrNodeBegin
c primary function:    Start node records

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  largeIntPut, fMallocI
      PTRFTN    fMallocI

      integer  i

      GeomHeader(9) = iLocL
      call largeIntPut (iLocL,GeomHeader(27),GeomHeader(28))

#if defined (FUJITSU_SYS)
         piIndexToDataL = fMallocI (nNode*1_8)
#else
         piIndexToDataL = fMallocI (nNode)
#endif


      nCount = 0
      do  i = 1,nNode
         IndexToData(i) = 0
      enddo

      return
      end


*deck,ResWrNode
      subroutine ResWrNode (iNode,xyzang)
c primary function:    Store a node

c object/library: ResWr

c  input arguments:
c     iNode    (int,sc,in)       - The node number
c     xyzang   (dp,ar(6),in)     - x,y,z,thxy,thyz,thzx for node

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      integer  iNode, i
      double precision  xyzang(6), dTemp(7)


      nCount = nCount + 1
      dTemp(1) = iNode
      do  i = 1,6
         dTemp(i+1) = xyzang(i)
      enddo
      IndexToData(nCount) = iNode
      call binwrt8 (NBLK,iLocL,7*intpdp,dTemp(1),0,iobuf(1))

      return
      end


*deck,ResWrNodeEnd
      subroutine ResWrNodeEnd
c primary function:    Finish Node Storage

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8, fMallocI, fFree
      PTRFTN    fMallocI

      integer   i, j


      call binwrt8 (NBLK,LocNodeListL,nNode,IndexToData(1),nNode,
     x iobuf(1))

#if defined (FUJITSU_SYS)
         piNinternalL = fMallocI(mxNode*1_8)
#else
         piNinternalL = fMallocI(mxNode)
#endif


      do  i = 1,mxNode
         iNinternal(i) = 0
      enddo

      do  i = 1,nNode
         j = IndexToData(i)
         iNinternal(j) = i
      enddo

      call fFree (piIndexToDataL)

      return
      end


*deck,ResWrElemBegin
      subroutine ResWrElemBegin
c primary function:   Start element output

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  fMallocL, fMallocI, binwrt8, largeIntPut
      PTRFTN    fMallocL, fMallocI

      integer  i


      maxPerElem = 0
#if defined (FUJITSU_SYS)
      plIndexL = fMallocL(nElem*1_8)
      piIndexToDataL  = fMallocI(nElem*1_8)
#else
      plIndexL = fMallocL(nElem)
      piIndexToDataL  = fMallocI(nElem)
#endif

      nCount = 0
      do  i = 1,nElem
         IndexToData(i)  = 0
         iIndexL(i) = 0
      enddo

      LocIndexL = iLocL
      GeomHeader(11) = iLocL
      call largeIntPut (iLocL,GeomHeader(29),GeomHeader(30))

      i = nElem*intplong
      call binwrt8 (NBLK,iLocL,i,IndexToData(1),i,iobuf(1))

      return
      end


*deck,ResWrElem
      subroutine ResWrElem (iElem, n, nodes, ElemData)
c primary function:    Write out an element

c object/library: ResWr

c  input arguments:
c     iElem    (int,sc,in)       - The element number
c     n        (int,sc,in)       - Number of nodes
c     nodes    (int,ar(n),in)    - Element nodes
c     ElemData (int,ar(10),in)   - Element information
c                                    mat    - material reference number
c                                    type   - element type number
c                                    real   - real constant reference number
c                                    elnum  - element number
c                                    esys   - element coordinate system
c                                    death  - death flag
c                                             = 0 - alive
c                                             = 1 - dead
c                                    solidm - solid model reference
c                                    shape  - coded shape key
c                                    pexcl  - P-Method exclude key
c                                    0      - Not used

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8

      integer   iElem, n, nodes(n), ElemData(10), i

      nCount = nCount + 1

      if (n .gt. maxPerElem) maxPerElem = n

      do  i = 1,10
         ElRecord(i) = ElemData(i)
      enddo

      ElRecord(9) = iElem
      do  i = 1,n
         ElRecord(i+10) = nodes(i)
      enddo

      iIndexL(nCount) = iLocL - LocIndexL
      call binwrt8 (NBLK,iLocL,n+10,ElRecord(1),n+10,iobuf(1))

      IndexToData(nCount) = iElem

      return
      end


*deck,ResWrElemEnd
      subroutine ResWrElemEnd
c primary function:   Finish element writing

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8, fMallocI, fFree, largeIntPutV
      PTRFTN    fMallocI

      integer   i, j

      GeomHeader(18) = maxPerElem

      call largeIntPutV (nElem,iIndexL(1))
      i = nElem*intplong
      call binwrt8 (NBLK,LocIndexL,i,iIndexL(1),i,iobuf(1))
      call binwrt8 (NBLK,LocElemListL,nElem,IndexToData(1),nElem,
     x iobuf(1))

      call fFree (plIndexL)

#if defined (FUJITSU_SYS)
      piEinternalL = fMallocI (mxElem*1_8)
#else
      piEinternalL = fMallocI (mxElem)
#endif

      do  i = 1,mxElem
         iEinternal(i) = 0
      enddo

      do  i = 1,nElem
         j = IndexToData(i)
         iEinternal(j) = i
      enddo

      call fFree (piIndexToDataL)

      return
      end



*deck,ResWrGeomEnd
      subroutine ResWrGeomEnd
c primary function:    Finish Geometry writing

c object/library: ResWr


#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8


      call binwrt8 (NBLK,LocGeomHeaderL,40,GeomHeader(1),40,iobuf(1))

      return
      end


*deck,ResWrSolBegin
      subroutine ResWrSolBegin (lstep,substep,ncumit,time,Title,DofLab)
c primary function:    Write the solution header records

c object/library: ResWr

c  input arguments:
c     lstep    (int,sc,in)       - Load step number
c     substep  (int,sc,in)       - Substep of this load step
c     ncumit   (int,sc,in)       - Cumulative iteration number
c     time     (int,sc,in)       - Current solution time
c     Title    (ch*80,ar(5),in)  - Title and 4 subtitles
c     DofLab   (ch*4,ar(nDOF),in)- Labels for DOFs

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      external  inexc4, binwrt8, binBasePut, binBaseGet, largeIntPut
      LONGINT   binBaseGet

      LONGINT  LongLocL, LongZeroL
      integer  lstep, substep, ncumit, i, ivect(200)
      double precision  time, dvect(100)
      character*80  Title(5)
      character*4   DofLab(*)

      equivalence (ivect(1),dvect(1))


      nSets = nSets + 1
      Header(9) = nSets
      dResTime(nSets) = time
      iResLstep(1,nSets) = lstep
      iResLstep(2,nSets) = substep
      iResLstep(3,nSets) = ncumit

      Mask = 0

      do  i = 1,200
         SolHeader(i) = 0
      enddo

      SolHeader(2) = nElem
      SolHeader(3) = nNode
      SolHeader(5) = lstep
      SolHeader(6) = substep
      SolHeader(7) = ncumit
      SolHeader(16)= 1
      SolHeader(18)= 1
      SolHeader(20)= nDOF
      do  i = 1,nDOF
         SolHeader(i+20) = curdof(i)
      enddo
      call inexc4 (Title(1),SolHeader(51),20)
      call inexc4 (Title(2),SolHeader(71),20)

      LongLocL = binBaseGet(NBLK) + iLocL
      call binBasePut (NBLK,LongLocL)
      iLocL = 0

      LocSolHeaderL = iLocL
      call largeIntPut (LongLocL,iResIndex(nSets),
     x                          iResIndex(nSets+MaxResult))
      call binwrt8 (NBLK,iLocL,200,SolHeader(1),200,iobuf(1))

      do  i = 1,100
         dvect(i) = 0.0d0
      enddo

      dvect(1) = time

      call binwrt8 (NBLK,iLocL,100*intpdp,dvect(1),0,iobuf(1))

      do  i = 1,200
         ivect(i) = 0
      enddo
      do  i = 1,nDOF
         ivect(i) = curdof(i)
      enddo
      call inexc4 (DofLab(1),ivect(33),nDOF)
      call inexc4 (Title(3),ivect(65),20)
      call inexc4 (Title(4),ivect(85),20)
      call inexc4 (Title(5),ivect(105),20)

      SolHeader(100)  = iLocL
      call largeIntPut (iLocL,SolHeader(103),SolHeader(104))

      call binwrt8 (NBLK,iLocL,200,ivect(1),200,iobuf(1))

c          *****  update the index  *****
      LongLocL = binBaseGet(NBLK)
      LongZeroL = 0
      call binBasePut (NBLK,LongZeroL)
      call binwrt8 (NBLK,LocResIndexL,
     x               2*MaxResult,iResIndex(1),2*MaxResult,iobuf(1))
      call binwrt8 (NBLK,LocResIndexL,
     x               MaxResult*intpdp,dResTime(1),0,iobuf(1))
      call binwrt8 (NBLK,LocResIndexL,
     x               3*MaxResult,iResLstep(1,1),3*MaxResult,iobuf(1))
      call binBasePut (NBLK,LongLocL)

      return
      end


*deck,ResWrDispBegin
      subroutine ResWrDispBegin
c primary function:     Begin output of displacement vector

c object/library: ResWr

c  input arguments:
c     NumDOF   (int,sc,in)       - Number of DOF per node

c  output arguments: none

#include "impcom.inc"
#include "ResWr.inc"

      external  fMallocD
      PTRFTN    fMallocD

      integer  i, j, k

#if defined (FUJITSU_SYS)
      pdDataL = fMallocD(nNode*nDOF*1_8)
#else
      pdDataL = fMallocD(nNode*nDOF)
#endif

      k = 0
      do  i = 1,nNode
         do  j = 1,nDOF
            k = k + 1
            dData(k) = 0.0d0
         enddo
      enddo

      return
      end

*deck,ResWrDisp
      subroutine ResWrDisp (node,Disp)
c primary function:    Store nodal displacements

c object/library: ResWr

c  input arguments:
c     node     (int,sc,in)       - Node number
c     Disp     (dp,ar(nDOF),in)  - Displacements
c     iNinternal(int,ar(mxNode),in) - External to internal node list

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      integer  node, i, j
      double precision  Disp(*)

      j = (iNinternal(node)-1)*nDOF
      do  i = 1,nDOF
         dData(j+i) = Disp(i)
      enddo

      return
      end

*deck,ResWrDispEnd
      subroutine ResWrDispEnd
c primary function:    Put displacement vector on result file

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  largeIntPut, fFree

      integer  n

      Mask = ibset(Mask,10)

      n = nNode*nDOF
      SolHeader(11) = iLocL
      call largeIntPut (iLocL,SolHeader(105),SolHeader(106))
      call binwrt8 (NBLK,iLocL,n*intpdp,dData(1),0,iobuf(1))

      call fFree (pdDataL)

      return
      end


*deck,ResWrRforBegin
      subroutine ResWrRforBegin (nRforce)
c primary function:   Start reaction force output

c object/library: ResWr

c  input arguments:
c     nRforce  (int,sc,in)       - Number of reaction forces

#include "impcom.inc"
#include "ResWr.inc"

      external  fMallocI, fMallocD
      PTRFTN    fMallocI, fMallocD

      integer  nRforce, i

#if defined (FUJITSU_SYS)
      piIndexToDataL = fMallocI(nRforce*1_8)
      pdDataL = fMallocD(nRforce*1_8)
#else
      piIndexToDataL = fMallocI(nRforce)
      pdDataL = fMallocD(nRforce)
#endif

      nCount = 0
      nRF    = nRforce

      do  i = 1,nRF
         IndexToData(i) = 0
         dData(i) = 0
      enddo


      return
      end

*deck,ResWrRfor
      subroutine ResWrRfor (node,idof,value)
c primary function:    Store a reaction force

c object/library: ResWr

c  input arguments:
c     node     (int,sc,in)       - External node number
c     idof     (int,sc,in)       - Internal dof number
c     value    (dp,sc,in)        - Value of reaction force

c  output arguments:  None

#include "impcom.inc"
#include "ResWr.inc"

      integer  node, idof
      double precision  value


      nCount = nCount + 1
      IndexToData(nCount) = (iNinternal(node)-1)*nDOF + idof
      dData(nCount) = value


      return
      end

*deck,ResWrRforEnd
      subroutine ResWrRforEnd
c primary function:   Write out the reaction forces

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  largeIntPut, fFree


      Mask = ibset(Mask,11)

      SolHeader(8) = nRF
      SolHeader(13) = iLocL
      call largeIntPut (iLocL,SolHeader(107),SolHeader(108))

      call binwrt8 (NBLK,iLocL,nRF,IndexToData(1),nRF,iobuf(1))
      call binwrt8 (NBLK,iLocL,nRF*intpdp,dData(1),0,iobuf(1))

      call fFree (pdDataL)
      call fFree (piIndexToDataL)


      return
      end


*deck,ResWrBCBegin
      subroutine ResWrBCBegin (nFixed,nForce)
c primary function:    Header for Boundary Conditions

c object/library: ResWr

c  input arguments:
c     nFixed   (int,sc,in)       - Number of Constraints
c     nForce   (int,sc,in)       - Number of Forces

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8, largeIntPut

      integer  nFixed, nForce, ivect(40), i

      do  i = 1,40
         ivect(i) = 0
      enddo

      nFix = nFixed
      nFor = nForce

      ivect(1) = nFix
      ivect(4) = nFor

      i = 40 + 3
      if (nFix .gt. 0) then
         ivect(2) = i
         i = i + nFix + 3
         ivect(3) = i
         i = i + 4*nFix*intpdp + 3
      endif
      if (nFor .gt. 0) then
         ivect(5) = i
         i = i + nFor + 3
         ivect(6) = i
      endif

      SolHeader(15) = iLocL
      call largeIntPut (iLocL,SolHeader(111),SolHeader(112))

      call binwrt8 (NBLK,iLocL,40,ivect(1),40,iobuf(1))

      Mask = ibset(Mask,23)

      return
      end


*deck,ResWrFixBegin
      subroutine ResWrFixBegin (nFixed)
c primary function:    Start constraint output

c object/library: ResWr

c  input arguments:
c     nFixed   (int,sc,in)          - Number of constraints

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      external  fMallocI, fMallocD
      PTRFTN    fMallocI, fMallocD

      integer   nFixed, i


#if defined (FUJITSU_SYS)
      piIndexToDataL = fMallocI(nFixed*1_8)
      pdDataL = fMallocD(4*nFixed*1_8)
#else
      piIndexToDataL = fMallocI(nFixed)
      pdDataL = fMallocD(4*nFixed)
#endif

      nCount = 0

      do  i = 1,nFixed
         IndexToData(i) = 0
      enddo
      do  i = 1,4*nFixed
         dData(i) = 0.0d0
      enddo

      return
      end

*deck,ResWrFix
      subroutine ResWrFix (node,idof,value)
c primary function:    Store a constraint value

c object/library: ResWr

c  input arguments:
c     node     (int,sc,in)       - External node number
c     idof     (int,sc,in)       - Internal dof number
c     value    (dp,ar(4),in)     - Real,Imag, RealOld,ImagOld

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      integer  node, idof, i, j
      double precision  value(4)


      j = 4*nCount
      nCount = nCount + 1

      IndexToData(nCount) = node*32 + idof
      do  i = 1,4
         j = j + 1
         dData(j) = value(i)
      enddo

      return
      end

*deck,ResWrFixEnd
      subroutine ResWrFixEnd
c primary function:    Write out the constraints

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8, fFree

      call binwrt8 (NBLK,iLocL,nFix,IndexToData(1),nFix,iobuf(1))
      call binwrt8 (NBLK,iLocL,4*nFix*intpdp,dData(1),0,iobuf(1))

      call fFree (pdDataL)
      call fFree (piIndexToDataL)

      return
      end


*deck,ResWrForcBegin
      subroutine ResWrForcBegin (nForce)
c primary function:    Start applied force output

c object/library: ResWr

c  input arguments:
c     nFixed   (int,sc,in)          - Number of constraints

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      external  fMallocI, fMallocD
      PTRFTN    fMallocI, fMallocD

      integer   nForce, i

#if defined (FUJITSU_SYS)
      piIndexToDataL = fMallocI(nForce*1_8)
      pdDataL = fMallocD(4*nForce*1_8)
#else
      piIndexToDataL = fMallocI(nForce)
      pdDataL = fMallocD(4*nForce)
#endif

      nCount = 0

      do  i = 1,nForce
         IndexToData(i) = 0
      enddo
      do  i = 1,4*nForce
         dData(i) = 0.0d0
      enddo


      return
      end

*deck,ResWrForc
      subroutine ResWrForc (node,idof,value)
c primary function:

c object/library: ResWr

c  input arguments:
c     node     (int,sc,in)       - External node number
c     idof     (int,sc,in)       - Internal dof number
c     value    (dp,ar(4),in)     - Real,Imag, RealOld,ImagOld

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      integer  node, idof, i, j
      double precision  value(4)


      j = 4*nCount
      nCount = nCount + 1

      IndexToData(nCount) = node*32 + idof
      do  i = 1,4
         j = j + 1
         dData(j) = value(i)
      enddo

      return
      end

*deck,ResWrForcEnd
      subroutine ResWrForcEnd
c primary function:    End of applied forces

c object/library: ResWr

c  input arguments:
c     IndexToData    (int,ar(nFor),in) - Dof of applied forces
c     dData    (dp,ar(4*nFor,in) - Values of applied forces

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8, fFree

      call binwrt8 (NBLK,iLocL,nFor,IndexToData(1),nFor,iobuf(1))
      call binwrt8 (NBLK,iLocL,4*nFor*intpdp,dData(1),0,iobuf(1))

      call fFree (pdDataL)
      call fFree (piIndexToDataL)

      return
      end


*deck,ResWrBCEnd
      subroutine ResWrBCEnd
c primary function:   End of BC output

c object/library: ResWr

c  input arguments:  none

c  output arguments:  none

#include "impcom.inc"


      return
      end


*deck,ResWrEresBegin
      subroutine ResWrEresBegin
c primary function:    Set up index array for element results

c object/library: ResWr

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8, largeIntPut, fMallocL
      PTRFTN    fMallocL

      integer  i

#if defined (FUJITSU_SYS)
      plIndexL = fMallocL(nElem*1_8)
#else
      plIndexL = fMallocL(nElem)
#endif

      do  i = 1,nElem
         iIndexL(i) = 0
      enddo

      SolHeader(12) = iLocL
      call largeIntPut (iLocL,SolHeader(119),SolHeader(120))
      LocIndexL = iLocL

      i = nElem*intplong
      call binwrt8 (NBLK,iLocL,i,iIndexL(1),i,iobuf(1))

      return
      end


*deck,ResWrEstrBegin
      subroutine ResWrEstrBegin  (iElem)
c primary function:    Start of element results

c object/library: ResWr


c  input arguments:
c     iElem    (int,sc,in)        - Element Number

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"


      integer  iElem, i, j


      Mask = ibset(Mask,12)

      do  i = 1,25
         ElemIndex(i) = 0
      enddo

      j = iEinternal(iElem)
      iIndexL(j) = iLocL - LocIndexL
      LocElemIndexL = iLocL

      call binwrt8 (NBLK,iLocL,25,ElemIndex(1),25,iobuf(1))

      return
      end

*deck,ResWrEstr
      subroutine ResWrEstr (iStr,nStr,Str)
c primary function:   Store element stress results

c object/library: ResWr


c  input arguments:
c     iStr    (int,sc,in)       - Stress record number (1-25)
c     nStr    (int,sc,in)       - Number of data words
c     Str     (dp,ar(nStr),in)  - Stress values

c  output arguments:  none

#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8

      integer  iStr, nStr, i
      double precision  Str(nStr), con

      con = 0.0d0
      do  i = 1,nStr
         con = con + abs(Str(i))
      enddo

      if (con .gt. 0.0d0) then
         ElemIndex(iStr) = iLocL - LocElemIndexL
         call binwrt8 (NBLK,iLocL,nStr*intpdp,Str(1),0,iobuf(1))
      else
         ElemIndex(iStr) = -nStr
      endif


      return
      end

*deck,ResWrEstrEnd
      subroutine ResWrEstrEnd
c primary function:     End of output for an element

c object/library: ResWr


#include "impcom.inc"
#include "ResWr.inc"


      call binwrt8 (NBLK,LocElemIndexL,25,ElemIndex(1),25,iobuf(1))

      return
      end


*deck,ResWrEresEnd
      subroutine ResWrEresEnd
c primary function:    Write out element index

c object/library: ResWr


#include "impcom.inc"
#include "ResWr.inc"

      external  fFree, largeIntPutV

      integer  i

      call largeIntPutV (nElem,iIndexL(1))
      i = nElem*intplong
      call binwrt8 (NBLK,LocIndexL,i,iIndexL(1),i,iobuf(1))

      call fFree (plIndexL)

      return
      end


*deck,ResWrSolEnd
      subroutine ResWrSolEnd
c primary function:    End of a Solution Result Set

c object/library: ResWr


#include "impcom.inc"
#include "ResWr.inc"

      external  binwrt8

      SolHeader(4) = Mask

      call binwrt8 (NBLK,LocSolHeaderL,200,SolHeader(1),200,iobuf(1))

      return
      end


*deck,ResWrEnd
      subroutine ResWrEnd
c primary function:    Close result file

c object/library: ResWr


#include "impcom.inc"
#include "ResWr.inc"

      external  fFree
      external  binwrt8, binclo, binBaseGet, binBasePut, largeIntPut
      LONGINT   binBaseGet

      LONGINT  LongLocL
      integer  i

      LongLocL = binBaseGet(NBLK) + iLocL
      call largeIntPut (LongLocL,Header(23),Header(24))
      Header(10) = LongLocL
      LongLocL = 0
      call binBasePut (NBLK,LongLocL)
      call binwrt8 (NBLK,LocHeaderL,40,Header(1),40,iobuf(1))

      i = NBLK
      call binclo (i,'KEEP',iobuf(1))

      call fFree (piEinternalL)
      call fFree (piNinternalL)


      call fFree (pdResTimeL)
      call fFree (piResLstepL)
      call fFree (piResIndexL)

      call fFree (piIobufL)

      return
      end


*deck,fMallocI
      function fMallocI (iLeng)
#include "impcom.inc"
#if defined(SOLARIS_SYS)
      external  malloc
#endif
      PTRFTN    malloc

#if defined (FUJITSU_SYS)
      integer*8  iLeng
#else
      integer  iLeng
#endif

      PTRFTN   fMallocI

#if defined(RS6000_SYS) || defined(RS64K_SYS)
      fMallocI = malloc(%val(iLeng*4))
#else 
      fMallocI = malloc(iLeng*4)
#endif

      return
      end

*deck,fMallocD
      function fMallocD (dLeng)
#include "impcom.inc"
#if defined(SOLARIS_SYS)
      external  malloc
#endif
      PTRFTN    malloc
#if defined (FUJITSU_SYS)
      integer*8  dLeng
#else
      integer  dLeng
#endif

      PTRFTN   fMallocD

#if defined(RS6000_SYS) || defined(RS64K_SYS)
      fMallocD = malloc(%val(dLeng*8))
#else 
      fMallocD = malloc(dLeng*8)
#endif

      return
      end

*deck,fMallocL
      function fMallocL (lLeng)
#include "impcom.inc"
#if defined(SOLARIS_SYS)
      external  malloc
#endif
      PTRFTN    malloc

#if defined (FUJITSU_SYS)
      integer*8  lLeng
#else
      integer  lLeng
#endif

      PTRFTN   fMallocL

#if defined(RS6000_SYS) || defined(RS64K_SYS)
      fMallocL = malloc(%val(lLeng*8))
#else 
      fMallocL = malloc(lLeng*8)
#endif

      return
      end

*deck,fFree
      subroutine fFree (piPointerL)
#include "impcom.inc"
#if defined(SOLARIS_SYS)
      external free
#endif
      PTRFTN  piPointerL

      call free (piPointerL)
      piPointerL = PTRFTNNULL

      return
      end
