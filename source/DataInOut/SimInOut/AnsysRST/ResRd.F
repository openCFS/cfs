*deck,ResRd        USERDISTRIB  User                                    ansys
c   **********   ANSYS RESULT FILE READING OBJECT   **********
c   Referenced libraries:  binlib.a
c   Include File:  ResRd.inc

c    This pair of files is to assist persons who would like to read a
c    result file which was written by the ANSYS program.


*deck,ResRdBegin
      function ResRdBegin (Nunit, Lunit, Fname, ncFname, Title, JobName,
     x                     Units, NumDOF, DOF, UserCode,
     x                     MaxNode, NumNode, MaxElem, NumElem,
     x                     MaxResultSet,NumResultSet)
c primary function:    Open result file and return global information

c object/library: ResRd

c  input arguments:
c     Nunit    (int,sc,in)       - Fortran Unit number for file (ANSYS uses 12)
c     Lunit    (int,sc,in)       - Current print output unit (usually 6 <STDOUT>)
c     Fname    (ch*(ncFname),sc,in) - The name (with extension) for the file
c     ncFname  (int,sc,in)       - Number of characters in Fname

c  output arguments:
c     Title    (ch*80,ar(2),out) - Title and First subtitle
c     JobName  (ch*32,sc,out)    - Jobname from file
c     Units    (int,sc,out)      - 0, unknown  1, SI  2, CSG
c                                  3, U.S. Customary - foot
c                                  4, U.S. Customary - inch
c                                  6, MPA
c     NumDOF   (int,sc,out)      - Number of DOF per node
c     DOF      (int,ar(*),out)   - The DOFs per node
c     UserCode (int,sc,out)      - Code for this application
c     MaxNode  (int,sc,out)      - Maximum node number used
c     NumNode  (int,sc,out)      - Number of nodes attached to elements
c     MaxElem  (int,sc,out)      - Maximum element number used
c     NumElem  (int,sc,out)      - Number of elements used
c     MaxResultSet (int,sc,out)  - Maximum number of result sets (usually 1000)
c     NumResultSet (int,sc,out)  - Number of result sets on file
c     ResRdBegin (int,sc,out)    - 0, successful  other, error in file open

#include "impcom.inc"
#include "ResRd.inc"

      external    fMallocI, fMallocD
      PTRFTN      fMallocI, fMallocD
      external    fFree
      external    biniqr8, binini, binset, binrd8, izero,
     x            binBasePut, exinc4, largeIntGet
      integer     binset
      LONGINT     biniqr8, largeIntGet

      LONGINT     LongZeroL
      integer     ResRdBegin, Nunit, Lunit, ncFname, Units, iTemp(100),
     x            NumDOF, DOF(*), UserCode, kbf, NumResultSet,
     x            MaxNode, NumNode, MaxElem, NumElem, MaxResultSet, i,j
      character*32 JobName
      character*(*)  Fname
      character*80  Title(2)

c          *****  set to .true. if want debug output  *****
      lDebug = .false.

c          *****  store arguments  *****
      iout   = Lunit

c          ***** initialize binary file I/O system  *****
      call binini (iout)
      intpdp = biniqr8(0,2)
      intplong = biniqr8(0,5)
      piIobufL = fMallocI(RECLENG)

c          *****  open result file  *****
      ResRdBegin = binset(NBLK,Nunit,1,1,RECLENG,1,Fname,ncFname,
     x             1,iobuf(1))
      if (ResRdBegin .ne. 0) then
         write (iout,2000) Fname(1:ncFname)
 2000    format (/' *** ERROR ***'/' Could not open file= ',a)
         go to 999
      endif

      LongZeroL = 0
      call binBasePut (NBLK,LongZeroL)

c          *****  read standard file header  *****
      iLocL = 0
      i = 100
      call binrd8 (NBLK,iLocL,i,iTemp(1),kbf,iobuf(1))
      JobName = ' '
      call exinc4 (iTemp(31),JobName,8)
      call exinc4 (iTemp(41),Title(1),20)
      call exinc4 (iTemp(61),Title(2),20)
      UserCode = iTemp(100)

c          *****  read results file header  *****
      LocHeaderL = iLocL
      i = 40
      call binrd8 (NBLK,iLocL,i,Header(1),kbf,iobuf(1))
      if (lDebug) write (*,3001) (Header(i),i=1,40)
 3001 format (' >Header=',10i6)
      mxNode    = Header(2)
      nNode     = Header(3)
      MaxResult = Header(4)
      nDOF      = Header(5)
      mxElem    = Header(6)
      nElem     = Header(7)
      NumResult = Header(9)
      Units     = Header(20)

      MaxNode = mxNode
      NumNode = nNode
      MaxElem = mxElem
      NumElem = nElem
      NumDOF  = nDOF
      MaxResultSet = MaxResult
      NumResultSet = NumResult

      LocElemListL = Header(14)
      LocNodeListL = Header(15)
      if (lDebug) write (*,3002) LocElemListL, LocNodeListL
 3002 format (' >LocElemListL=',i6,'  LocNodeListL=',i6)

c          *****  Read DOF list  *****
      call binrd8 (NBLK,iLocL,nDOF,curdof(1),kbf,iobuf(1))
      do  i = 1,nDOF
         DOF(i) = curdof(i)
      enddo
      if (lDebug) write (*,3003) (DOF(i),i=1,nDOF)
 3003 format (' >DOF=',10i6)

c          *****  Read nodal equivalence table  *****
      piIndexToDataL = fMallocI (nNode)
      iLocL = LocNodeListL
      i = nNode
      call binrd8 (NBLK,iLocL,i,IndexToData(1),kbf,iobuf(1))

      piNinternalL = fMallocI(mxNode)

c     *** create external-to-internal map
      do  i = 1,mxNode
         iNinternal(i) = 0
      enddo
      do  i = 1,nNode
         j = IndexToData(i)
         iNinternal(j) = i
      enddo

      if (lDebug) write (*,3007) (IndexToData(i),i=1,nNode)
 3007 format (' >IndexToData=',10i6)
      if (lDebug) write (*,3008) (iNinternal(i),i=1,mxNode)
 3008 format (' >iNinternal=',10i6)

      call fFree (piIndexToDataL)

c          *****  Read element equivalence table  *****
      piIndexToDataL = fMallocI (nElem)

      iLocL = LocElemListL
      i = nElem
      call binrd8 (NBLK,iLocL,i,IndexToData(1),kbf,iobuf(1))

      piEinternalL = fMallocI(mxElem)

c     *** create external-to-internal map
      do  i = 1,mxElem
         iEinternal(i) = 0
      enddo
      do  i = 1,nElem
         j = IndexToData(i)
         iEinternal(j) = i
      enddo

      if (lDebug) write (*,3009) (IndexToData(i),i=1,nElem)
 3009 format (' >IndexToData=',10i6)
      if (lDebug) write (*,3010) (iEinternal(i),i=1,mxElem)
 3010 format (' >iEinternal=',10i6)

      call fFree (piIndexToDataL)

c          *****  initialize result indexes  *****
      pdResTimeL  = fMallocD(MaxResult*1_8)
      piResIndexL = fMallocI(2*MaxResult)
      piResLstepL = fMallocI(3*MaxResult)
      pdResTimeL  = fMallocD(MaxResult)

c          *****  read result indexes  *****
      nSets = Header(9)

      iLocL = Header(11)
      LocResIndexL = iLocL
      i = 2*MaxResult
      call binrd8 (NBLK,iLocL,i,iResIndex(1),kbf,iobuf(1))
      if (i.lt.2*MaxResult) call izero (iResIndex(MaxResult+1),
     x  MaxResult)

      iLocL = Header(12)
      i = MaxResult*intpdp
      call binrd8 (NBLK,iLocL,i,dResTime(1),kbf,iobuf(1))

      iLocL = Header(13)
      i = 3*MaxResult
      call binrd8 (NBLK,iLocL,i,iResLstep(1,1),kbf,iobuf(1))

      if (lDebug) then
         write (*,'(a)') ' > Set    Ptr    LS   SS  NC     Time'
         do  i = 1,nSets
            write (*,3004) i,
     x         largeIntGet(iResIndex(i),iResIndex(i+MaxResult)),
     x         (iResLstep(j,i),j=1,3),dResTime(i)
 3004       format (' >',i4,i8,3i5,g12.5)
         enddo
      endif

c          *****  read geometry header to determine if long file  *****
      lLongRes = .true.
      iLocL  = Header(16)
      i = 40
      call binrd8 (NBLK,iLocL,i,iTemp(1),kbf,iobuf(1))
      if (i.lt.40) lLongRes = .false.
      if (lLongRes) then
         if (lDebug) write (*,3005)
 3005    format (' >>> Long result file format')
      else
         if (lDebug) write (*,3006)
 3006    format (' >>> Short result file format')
      endif

 999  return
      end



*deck,ResRdGeomBegin
      subroutine ResRdGeomBegin (MaxType, MaxReal, MaxCsys, nXYZ)
c primary function:    Read Geometry Header Record

c object/library: ResRd

c  input arguments:  none

c  output arguments:
c     MaxType  (int,sc,out)       - Maximum element type
c     MaxReal  (int,sc,out)       - Maximum real constant set number
c     MaxCsys  (int,sc,out)       - Maximum coordinate system number
c     nXYZ     (int,sc,out)       - number of nodes with coordinates

#include "impcom.inc"
#include "ResRd.inc"

      integer   MaxType, MaxReal, MaxCsys, nXYZ, i, kbf

c          *****  read geometry header  *****
      iLocL  = Header(16)
      LocGeomHeaderL = iLocL
      i = 40
      call binrd8 (NBLK,iLocL,i,GeomHeader(1),kbf,iobuf(1))

      if (lDebug) write (*,3001) (GeomHeader(i),i=1,40)
 3001 format (' >GeomHeader=',10i6)

      mxType = GeomHeader(2)
      mxReal = GeomHeader(3)
      nXYZ = GeomHeader(4)
      nXYZnode = nXYZ
      nElem  = GeomHeader(5)
      mxCsys = GeomHeader(6)

      MaxType = mxType
      MaxReal = mxReal
      MaxCsys = mxCsys

      return
      end


*deck,ResRdTypeBegin
      subroutine ResRdTypeBegin (NumType)
c primary function:    Start element type records

c object/library: ResRd

c  input arguments: none

c  output arguments:
c    NumType  (int,sc,out)   - Number of defined element types

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8, largeIntGet, fMallocI
      LONGINT   largeIntGet
      PTRFTN    fMallocI

      integer   NumType, i, kbf


      NumType = 0
      if (mxType .le. 0) go to 999

      piIndexToDataL = fMallocI(mxType)

c          *****  Read index  *****
      if (lLongRes) then
         iLocL = largeIntGet(GeomHeader(21),GeomHeader(22))
         LocIndexL = iLocL
      else
         iLocL = GeomHeader(7)
         LocIndexL = 0
      endif
      if (iLocL .gt. 0) then
         call binrd8 (NBLK,iLocL,mxType,IndexToData(1),kbf,iobuf(1))
         do  i = 1,mxType
            if (IndexToData(i) .ne. 0) NumType = NumType + 1
         enddo
      endif

      if (lDebug) write (*,3001) (IndexToData(i),i=1,mxType)
 3001 format (' >IndexToData=',10i6)

 999  return
      end


*deck,ResRdType
      function ResRdType (itype,ielc)
c primary function:    Read an element type record

c object/library: ResRd

c  input arguments:
c     itype    (int,sc,on)         - Element type number

c  output arguments:  none
c     ielc     (int,ar(IELCSZ),out)- Element characteristics
c     ResRdType (int,sc,out        - number of words read


#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8

#include "echprm.inc"

      integer   ResRdType, itype, ielc(IELCSZ), kbf, i

      ResRdType = 0

      if (IndexToData(itype) .ne. 0) then
         iLocL = abs(IndexToData(itype)) + LocIndexL
         i = IELCSZ
         call binrd8 (NBLK,iLocL,i,ielc(1),kbf,iobuf(1))
         ResRdType = i
      endif

      return
      end


*deck,ResRdTypeEnd
      subroutine ResRdTypeEnd
c primary function:    End element types

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external   fFree

      if (mxType .gt. 0) then
         call fFree (piIndexToDataL)
      endif

      return
      end


*deck,ResRdRealBegin
      subroutine ResRdRealBegin (NumReal,NumPerReal)
c primary function:    Start real constant records

c object/library: ResRd

c  input arguments: none

c  output arguments:
c    NumReal     (int,sc,out)   - Number of real constant sets defined
c    NumPerReal  (int,sc,out)   - Maximum number of real constant in a set

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8, largeIntGet,fMallocI
      LONGINT   largeIntGet
      PTRFTN    fMallocI

      integer   NumReal, NumPerReal, i, j, kbf
      double precision  dTemp(1000)

      NumPerReal = 0
      NumReal    = 0

      if (mxReal .le. 0) go to 999

      piIndexToDataL = fMallocI(mxReal)

c          *****  Read index  *****
      if (lLongRes) then
         iLocL = largeIntGet(GeomHeader(23),GeomHeader(24))
         LocIndexL = iLocL
      else
         iLocL = GeomHeader(8)
         LocIndexL = 0
      endif
      if (iLocL .gt. 0) then
         call binrd8 (NBLK,iLocL,mxReal,IndexToData(1),kbf,iobuf(1))
         do  i = 1,mxReal
            if (IndexToData(i) .ne. 0) then
               iLocL = abs(IndexToData(i)) + LocIndexL
               call binrd8 (NBLK,iLocL,j,dTemp(1),kbf,iobuf(1))
               j = j/intpdp
               if (j .gt. NumPerReal) NumPerReal = j
               NumReal = NumReal + 1
            endif
         enddo
      endif

      if (lDebug) write (*,3001) (IndexToData(i),i=1,mxReal)
 3001 format (' >IndexToData=',10i6)

 999  return
      end


*deck,ResRdReal
      function ResRdReal (iReal,Rcon)
c primary function:    Read real constant record

c object/library: ResRd

c  input arguments:
c     iReal    (int,sc,in)             - Real set number

c  output arguments:  none
c     Rcon      (dp,ar(ResRdReal),out) - Real Constants
c     ResRdReal (int,sc,out)           - Number of real constants in set

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8

      integer  ResRdReal, iReal, i, kbf
      double precision  Rcon(*)

      ResRdReal = 0

      if (IndexToData(iReal) .ne. 0) then
         iLocL = abs(IndexToData(iReal)) + LocIndexL
         i = 10
         call binrd8 (NBLK,iLocL,i,Rcon(1),kbf,iobuf(1))
         ResRdReal = i/intpdp
      endif

      return
      end


*deck,ResRdRealEnd
      subroutine ResRdRealEnd
c primary function:    End real constants

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external   fFree

      if (mxReal .gt. 0) then
         call fFree (piIndexToDataL)
      endif

      return
      end


*deck,ResRdCsysBegin
      subroutine ResRdCsysBegin (NumCsys)
c primary function:    Start coordinate system records

c object/library: ResRd

c  input arguments: none

c  output arguments:
c    NumCsys  (int,sc,out)   - Number of defined coordinate systems

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8, largeIntGet, fMallocI
      LONGINT   largeIntGet
      PTRFTN    fMallocI

      integer   NumCsys, i, kbf


      NumCsys = 0
      if (mxCsys .le. 0) go to 999

      piIndexToDataL = fMallocI(mxCsys)

c          *****  read index *****
      if (lLongRes) then
         iLocL = largeIntGet(GeomHeader(25),GeomHeader(26))
         LocIndexL = iLocL
      else
         iLocL = GeomHeader(10)
         LocIndexL = 0
      endif
      if (iLocL .gt. 0) then
         call binrd8 (NBLK,iLocL,mxCsys,IndexToData(1),kbf,iobuf(1))
         do  i = 1,mxCsys
            if (IndexToData(i) .ne. 0) NumCsys = NumCsys + 1
         enddo
      endif

      if (lDebug) write (*,3001) (IndexToData(i),i=1,mxCsys)
 3001 format (' >IndexToData=',10i6)

 999  return
      end


*deck,ResRdCsys
      function ResRdCsys (iCsys,Csys)
c primary function:    Read a coordinate system record

c object/library: ResRd

c  input arguments:
c     iCsys     (int,sc,in)           - Coordinate system number

c  output arguments:
c     Csys      (dp,ar(ResRdCsys),out)- Coordinate system description
c     ResRdCsys (int,sc,out)          - Number of values

c  output arguments:

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8

      integer   ResRdCsys, iCsys, kbf, i
      double precision  Csys(24)

      ResRdCsys = 0

      if (IndexToData(iCsys) .ne. 0) then
         iLocL = abs(IndexToData(iCsys)) + LocIndexL
          i = 24*intpdp
         call binrd8 (NBLK,iLocL,i,Csys(1),kbf,iobuf(1))
         ResRdCsys = i/intpdp
      endif

      return
      end


*deck,ResRdCsysEnd
      subroutine ResRdCsysEnd
c primary function:    End coordinate systems

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external   fFree

      if (mxCsys .gt. 0) then
         call fFree (piIndexToDataL)
      endif

      return
      end


*deck,ResRdNodeBegin
      subroutine ResRdNodeBegin
c primary function:    Start node records

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  largeIntGet
      LONGINT   largeIntGet


      if (lLongRes) then
         LocIndexL = largeIntGet(GeomHeader(27),GeomHeader(28))
      else
         LocIndexL = GeomHeader(9)
      endif

      pdNodesL = PTRFTNNULL

      return
      end


*deck,ResRdNode
      function ResRdNode (iNode,xyzang)
c primary function:    Get a node

c object/library: ResRd

c  input arguments:
c     iNode     (int,sc,in)     - node sequence number
c                                  (1 to nXYZnode)

c  output arguments:
c     xyzang    (dp,ar(6),out)  - x,y,z,thxy,thyz,thzx for node
c     ResRdNode (int,sc,out)    - Node number

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8, fMallocD
      PTRFTN    fMallocD

      integer  ResRdNode, iNode, i, kbf, loc
      double precision  xyzang(6), dTemp(7)


      ResRdNode = 0

      if (lLongRes) then
         iLocL = LocIndexL + (7*intpdp+3)*(iNode-1)
         i = 7*intpdp
         call binrd8 (NBLK,iLocL,i,dTemp(1),kbf,iobuf(1))

         ResRdNode = dTemp(1)

         do  i = 1,6
            xyzang(i) = dTemp(i+1)
         enddo

      else
         if (pdNodesL.eq.PTRFTNNULL) then
            pdNodesL = fMallocD (7*nXYZnode)
            i = 7*nXYZnode*intpdp
            call binrd8 (NBLK,LocIndexL,i,dNodes(1),kbf,iobuf(1))
         endif

         loc = 7*(iNode-1) + 1
         ResRdNode = dNodes(loc)

         do  i = 1,6
            xyzang(i) = dNodes(loc+i)
         enddo

      endif

      return
      end


*deck,ResRdNodeEnd
      subroutine ResRdNodeEnd
c primary function:    End Nodes

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fFree

      call fFree (pdNodesL)

      return
      end


*deck,ResRdElemBegin
      subroutine ResRdElemBegin
c primary function:   Start element input

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fMallocL, fMallocI, binrd8, largeIntGet, largeIntGetV
      external  fFree
      PTRFTN    fMallocL, fMallocI
      LONGINT   largeIntGet

      integer  i, kbf

      pointer (piWorkL,iWork)
      integer iWork(*)

      plIndexL = PTRFTNNULL

      if (lLongRes) then
         LocIndexL = largeIntGet(GeomHeader(29),GeomHeader(30))
      else
         LocIndexL = GeomHeader(11)
      endif

      if (LocIndexL .gt. 0) then

c          *****  read element description index  *****
         plIndexL = fMallocL(nElem)
         iLocL = LocIndexL
         if (lDebug) write (*,3000) LocIndexL
 3000    format (' >LocIndexL=',i6)

         if (lLongRes) then
            i = nElem*intplong
            call binrd8 (NBLK,iLocL,i,iIndexL(1),kbf,iobuf(1))
            call largeIntGetV (nElem,iIndexL(1))
         else
            piWorkL = fMallocI(nElem)
            i = nElem
            call binrd8 (NBLK,iLocL,i,iWork(1),kbf,iobuf(1))
            do i = 1,nElem
               iIndexL(i) = iWork(i)
            enddo
            call fFree (piWorkL)
            LocIndexL = 0
         endif
         if (lDebug) write (*,3003) (iIndexL(i),i=1,nElem)
 3003    format (' >iIndexL=',10i6)

      endif

      return
      end


*deck,ResRdElem
      function ResRdElem (iElem, nodes, ElemData)
c primary function:    Read an element

c object/library: ResRd

c  input arguments:
c     iElem    (int,sc,in)       - The element number

c  output arguments:
c     ResRdElem(int,sc,out)      - Number of nodes
c     nodes    (int,ar(n),out)   - Element nodes
c     ElemData (int,ar(10),out)  - Element information
c                                    mat    - material reference number
c                                    type   - element type number
c                                    real   - real constant reference number
c                                    secnum - section number
c                                    esys   - element coordinate system
c                                    death  - death flag
c                                             = 0 - alive
c                                             = 1 - dead
c                                    solidm - solid model reference
c                                    shape  - coded shape key
c                                    elnum  - element number
c                                    pexcl  - P-Method exclude key

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8

      integer   ResRdElem, iElem, n, nodes(*), ElemData(10), i, kbf,
     x          iEl

      ResRdElem = 0

      iEl = iEinternal(iElem)
      if (iEl .gt. 0) then
         if (lDebug) write (*,3001) iElem,iEl,iIndexL(iEl)
 3001    format (' >iElem=',i4,'  iEl=',i4,'  iIndexL(iEl)=',i6)
         if (iIndexL(iEl) .ne. 0) then
            iLocL = iIndexL(iEl) + LocIndexL
            n = 10
            call binrd8 (NBLK,iLocL,n,ElRecord(1),kbf,iobuf(1))

            do  i = 1,10
               ElemData(i) = ElRecord(i)
            enddo

            do  i = 11,n
               nodes(i-10) = ElRecord(i)
            enddo

            ResRdElem = n - 10
         endif
      endif

      return
      end


*deck,ResRdElemEnd
      subroutine ResRdElemEnd
c primary function:   Finish element writing

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fFree

      call fFree (plIndexL)

      return
      end



*deck,ResRdGeomEnd           parallel
      subroutine ResRdGeomEnd
c primary function:    Finish Geometry Input

c object/library: ResRd


#include "impcom.inc"

      return
      end


*deck,ResRdSolBegin
      function ResRdSolBegin (key,lstep,substep,ncumit,kcmplx,time,
     x                        Title,DofLab)
c primary function:    Read the solution header records

c object/library: ResRd

c  input arguments:
c     key      (int,sc,in)       - 1, find by lstep/substep
c                                  2, find by ncumit
c                                  3, find by time
c     lstep    (int,sc,in/out)   - Load step number
c     substep  (int,sc,in/out)   - Substep of this load step
c     ncumit   (int,sc,in/out)   - Cumulative iteration number
c     kcmplx   (int,sc,in)       - 0, Real solution   1, Imaginary solution
c     time     (int,sc,in/out)   - Current solution time

c  output arguments:
c     Title    (ch*80,ar(5),out)  - Title and 4 subtitles
c     DofLab   (ch*4,ar(nDOF),out)- Labels for DOFs
c     ResRdSolBegin (int,sc,out)  - 0, requested solution set found
c                                   1, not found

#include "impcom.inc"
#include "ResRd.inc"

      external  exinc4, binrd8, binBasePut, largeIntGet, izero
      LONGINT   largeIntGet

      integer  ResRdSolBegin, key, lstep, substep, ncumit, kcmplx,
     x         i, ivect(200), kbf
      double precision  time, dvect(100)
      character*80  Title(5)
      character*4   DofLab(*)

      ResRdSolBegin = 1

c *** find data set no. by load step and substep

      if (key .eq. 1)  then

         do  nSets = 1,NumResult
            if (iResLstep(1,nSets).ne.lstep .or.
     x          iResLstep(2,nSets).ne.substep) cycle
            if (kcmplx .eq. 1)  then
               if (iResLstep(1,nSets+1).ne.lstep .or.
     x             iResLstep(2,nSets+1).ne.substep)  cycle
            endif
            go to 10
         enddo
         go to 999

c ***  find data set by ncumit

      elseif (key .eq. 2) then

         do  nSets = 1,NumResult
            if (iResLstep(3,nSets) .ne. ncumit)  cycle
            if (kcmplx .eq. 1) then
               if (iResLstep(3,nSets+1) .ne. ncumit)  cycle
            endif
            go to 10
         enddo
         go to 999

c *** find data set. no. by time/freq

      elseif (key .eq. 3)  then

c        *** check for 0 time
         if (time .eq. 0.0d0)  go to 999
c        *** get time/freq index
         do  nSets = 1, NumResult
            if (dResTime(nSets) .lt. time)  cycle
            if (kcmplx .eq. 1)  then
               if (dResTime(nSets+1) .lt. time) cycle
            endif
            go to 10
         enddo
         go to 999

      endif

c          *****  result set found - read load set headers  *****
  10  ResRdSolBegin = 0

      LocSolHeaderL = largeIntGet(iResIndex(nSets),
     x                            iResIndex(nSets+MaxResult))
      call binBasePut(NBLK,LocSolHeaderL)

      iLocL = 0
      i = 200
      call binrd8 (NBLK,iLocL,i,SolHeader(1),kbf,iobuf(1))
      if (i.lt.200) then
         call izero (SolHeader(101),100)
         lLongRes = .false.
      endif
      if (lDebug) write (*,2001) SolHeader
 2001 format (' >SolHeader=',10i6)

      lstep   = SolHeader(5)
      substep = SolHeader(6)
      ncumit  = SolHeader(7)
      call exinc4 (SolHeader(51),Title(1),20)
      call exinc4 (SolHeader(71),Title(2),20)

      i = 100*intpdp
      call binrd8 (NBLK,iLocL,i,dvect(1),kbf,iobuf(1))

      time = dvect(1)

      i = 200
      call binrd8 (NBLK,iLocL,i,ivect(1),kbf,iobuf(1))

      call exinc4 (ivect(33), DofLab(1),nDOF)
      call exinc4 (ivect(65), Title(3),20)
      call exinc4 (ivect(85), Title(4),20)
      call exinc4 (ivect(105),Title(5),20)

 999  return
      end


*deck,ResRdDispBegin
      subroutine ResRdDispBegin
c primary function:     Read displacement vector

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fMallocD, largeIntGet
      PTRFTN    fMallocD
      LONGINT   largeIntGet

      integer  i, kbf

      i = nNode*nDOF
      pdDataL = fMallocD(i)

      if (lLongRes) then
         iLocL = largeIntGet(SolHeader(105),SolHeader(106))
      else
         iLocL = SolHeader(11)
      endif
      if (iLocL .gt. 0) then
         i = i*intpdp
         call binrd8 (NBLK,iLocL,i,Data(1),kbf,iobuf(1))
      endif

      return
      end

*deck,ResRdDisp
      function ResRdDisp (node,Disp)
c primary function:    Retrieve a nodal displacement

c object/library: ResRd

c  input arguments:
c     node     (int,sc,in)       - Node number

c  output arguments:  none
c     Disp     (dp,ar(nDOF),out) - Displacements
c     ResRdDisp(int,sc,out)      - Number of displacements

#include "impcom.inc"
#include "ResRd.inc"

      integer  ResRdDisp, node, i, j
      double precision  Disp(*)

      ResRdDisp = 0
      if (iNinternal(node) .gt. 0) then
         j = (iNinternal(node)-1)*nDOF
         do  i = 1,nDOF
            Disp(i) = Data(j+i)
         enddo
         ResRdDisp = nDOF
      endif

      return
      end

*deck,ResRdDispEnd
      subroutine ResRdDispEnd
c primary function:    Finish displacements

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fFree

      call fFree (pdDataL)

      return
      end


*deck,ResRdRforBegin
      subroutine ResRdRforBegin (nRForce)
c primary function:   Read reaction forces

c object/library: ResRd

c  input arguments: none

c  output arguments:
c     nRForce    (int,sc,out)     - Number of reactions

#include "impcom.inc"
#include "ResRd.inc"

      external  fMallocL, fMallocD, largeIntGet, binrd8
      PTRFTN    fMallocL, fMallocD
      LONGINT   largeIntGet

      integer  nRForce, kbf, i

      nRF = SolHeader(8)
      nRForce = nRF

      if (nRF .gt. 0) then
         plIndexToDataL = fMallocL(nRF)
         piIndexToDataL = plIndexToDataL
         pdDataL = fMallocD(nRF)

         if (lLongRes) then
            iLocL = largeIntGet(SolHeader(107),SolHeader(108))
         else
            iLocL = SolHeader(13)
         endif

         i = nRF
         call binrd8 (NBLK,iLocL,i,IndexToData(1),kbf,iobuf(1))
         if (i .eq. nRF) then
            do  i = nRF,1,-1
               IndexToDataL(i) = IndexToData(i)
            enddo
         endif
         i = nRF*intpdp
         call binrd8 (NBLK,iLocL,i,Data(1),kbf,iobuf(1))

      endif

      return
      end

*deck,ResRdRfor
      function ResRdRfor (node,idof,value)
c primary function:    Retrieve a reaction force

c object/library: ResRd

c  input arguments:
c     node     (int,sc,in)       - External node number
c     idof     (int,sc,in)       - Internal dof number

c  output arguments:
c     value     (dp,sc,in)       - Value of reaction force
c     ResRdRfor (int,sc,out)     - Number of returned values (0 or 1)

#include "impcom.inc"
#include "ResRd.inc"

      LONGINT  kL
      integer  ResRdRfor, node, idof, iNode, i
      double precision  value

      ResRdRfor = 0

      if (nRF .gt. 0) then

         iNode = iNinternal(node)
         if (iNode .gt. 0) then
            kL = iNode - 1
            kL = kL*nDOF + idof
            do  i = 1,nRF
               if (kL .eq. IndexToDataL(i)) then
                  ResRdRfor = 1
                  value = Data(i)
                  exit
               endif
            enddo
         endif

      endif


      return
      end

*deck,ResRdRforEnd
      subroutine ResRdRforEnd
c primary function:   Finish reaction forces

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fFree


      call fFree (pdDataL)
      call fFree (plIndexToDataL)


      return
      end


*deck,ResRdBCBegin
      subroutine ResRdBCBegin (BCHeader)
c primary function:    Read Header for Boundary Conditions

c object/library: ResRd

c  input arguments: none

c  output arguments:
c     BCHeader   (int,ar(40),out)    - BC Header

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8, izero, largeIntGet
      LONGINT   largeIntGet

      integer  BCHeader(40), i, kbf

      if (lLongRes) then
         iLocL = largeIntGet(SolHeader(111),SolHeader(112))
         LocIndexL = iLocL
      else
         iLocL = SolHeader(15)
         LocIndexL = 0
      endif
      if (iLocL .gt. 0) then
         i = 40
         call binrd8 (NBLK,iLocL,i,BCHeader(1),kbf,iobuf(1))
      else
         call izero (BCHeader(1),40)
      endif

      return
      end


*deck,ResRdFixBegin
      subroutine ResRdFixBegin (BCHeader,nFixed)
c primary function:    Read constraint information

c object/library: ResRd

c  input arguments:
c     BCHeader   (int,ar(40),in)     - BC Header

c  output arguments:
c     nFixed     (int,sc,out)        - number of constraints

#include "impcom.inc"
#include "ResRd.inc"

      external  fMallocI, fMallocD
      PTRFTN    fMallocI, fMallocD

      integer   BCHeader(40), nFixed, i, kbf, format

      piIndexToDataL = PTRFTNNULL
      pibcNodesL = PTRFTNNULL
      pibcDofsL = PTRFTNNULL
      pdDataL = PTRFTNNULL

      nFixed = BCHeader(1)
      nFix = nFixed
      format = BCHeader(7)

      if (nFixed .gt. 0) then
         if (format .eq. 0) then
            piIndexToDataL = fMallocI(nFixed)
         else
            pibcNodesL = fMallocI(nFixed)
            pibcDofsL  = fMallocI(nFixed)
         endif
         pdDataL = fMallocD(4*nFixed)

         iLocL = BCHeader(2) + LocIndexL
         if (format .eq. 0) then
            i = nFixed
            call binrd8 (NBLK,iLocL,i,IndexToData(1),kbf,iobuf(1))
         else
            i = nFixed
            call binrd8 (NBLK,iLocL,i,bcNodes(1),kbf,iobuf(1))
            i = nFixed
            call binrd8 (NBLK,iLocL,i,bcDofs(1),kbf,iobuf(1))
         endif

         iLocL = BCHeader(3) + LocIndexL
         i = 4*nFixed*intpdp
         call binrd8 (NBLK,iLocL,i,Data(1),kbf,iobuf(1))

      endif

      return
      end

*deck,ResRdFixOld
      function ResRdFixOld (node,idof,value)
c primary function:    Retrieve a constraint value (in old format)

c object/library: ResRd

c  input arguments:
c     node        (int,sc,in)       - External node number
c     idof        (int,sc,in)       - Internal dof number

c  output arguments:
c     value       (dp,ar(4),in)     - Real,Imag, RealOld,ImagOld
c     ResRdFixOld (int,sc,out)      - Number of returned values (0 or 4)

#include "impcom.inc"
#include "ResRd.inc"

      integer  node, idof, i, k, iNode, ResRdFixOld
      double precision  value(4)

      ResRdFixOld = 0

      if (nFix .gt. 0) then

         iNode = iNinternal(node)
         if (iNode .gt. 0) then
            k = node*32 + idof
            do  i = 1,nFix
               if (k .eq. IndexToData(i)) then
                  ResRdFixOld = 4
                  value(1) = Data(4*(i-1)+1)
                  value(2) = Data(4*(i-1)+2)
                  value(3) = Data(4*(i-1)+3)
                  value(4) = Data(4*(i-1)+4)
                  exit
               endif
            enddo
         endif

      endif

      return
      end

*deck,ResRdFix
      function ResRdFix (node,idof,value)
c primary function:    Retrieve a constraint value

c object/library: ResRd

c  input arguments:
c     node     (int,sc,in)       - External node number
c     idof     (int,sc,in)       - Internal dof number

c  output arguments:
c     value    (dp,ar(4),in)     - Real,Imag, RealOld,ImagOld
c     ResRdFix (int,sc,out)      - Number of returned values (0 or 4)

#include "impcom.inc"
#include "ResRd.inc"

      integer  node, idof, i, iNode, ResRdFix, bcNode, bcDof
      double precision  value(4)

      ResRdFix = 0

      if (nFix .gt. 0) then

         iNode = iNinternal(node)
         if (iNode .gt. 0) then
            do  i = 1,nFix
               bcNode = bcNodes(i)
               bcDof = bcDofs(i)
               if (iNode.eq.bcNode .and. idof.eq.bcDof) then
                  ResRdFix = 4
                  value(1) = Data(4*(i-1)+1)
                  value(2) = Data(4*(i-1)+2)
                  value(3) = Data(4*(i-1)+3)
                  value(4) = Data(4*(i-1)+4)
                  exit
               endif
            enddo
         endif

      endif

      return
      end

*deck,ResRdFixEnd
      subroutine ResRdFixEnd
c primary function:    Finish the constraints

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fFree


      call fFree (pdDataL)
      call fFree (piIndexToDataL)
      call fFree (pibcNodesL)
      call fFree (pibcDofsL)

      return
      end


*deck,ResRdForcBegin
      subroutine ResRdForcBegin (BCHeader,nForces)
c primary function:    Read applied force information

c object/library: ResRd

c  input arguments:
c     BCHeader   (int,ar(40),in)     - BC Header

c  output arguments:
c     nForces    (int,sc,out)        - number of forces

#include "impcom.inc"
#include "ResRd.inc"

      external  fMallocI, fMallocD
      PTRFTN    fMallocI, fMallocD

      integer   BCHeader(40), nForces, i, kbf, format

      piIndexToDataL = PTRFTNNULL
      pibcNodesL = PTRFTNNULL
      pibcDofsL = PTRFTNNULL
      pdDataL = PTRFTNNULL

      nForces = BCHeader(4)
      nFor = nForces
      format = BCHeader(7)

      if (nForces .gt. 0) then
         if (format .eq. 0) then
            piIndexToDataL = fMallocI(nForces)
         else
            pibcNodesL = fMallocI(nForces)
            pibcDofsL  = fMallocI(nForces)
         endif
         pdDataL = fMallocD(4*nForces)

         iLocL = BCHeader(5) + LocIndexL
         if (format .eq. 0) then
            i = nForces
            call binrd8 (NBLK,iLocL,i,IndexToData(1),kbf,iobuf(1))
         else
            i = nForces
            call binrd8 (NBLK,iLocL,i,bcNodes(1),kbf,iobuf(1))
            i = nForces
            call binrd8 (NBLK,iLocL,i,bcDofs(1),kbf,iobuf(1))
         endif

         iLocL = BCHeader(6) + LocIndexL
         i = 4*nForces*intpdp
         call binrd8 (NBLK,iLocL,i,Data(1),kbf,iobuf(1))

      endif

      return
      end

*deck,ResRdForcOld
      function ResRdForcOld (node,idof,value)
c primary function:    Retrieve an applied force value (in old format)

c object/library: ResRd

c  input arguments:
c     node        (int,sc,in)       - External node number
c     idof        (int,sc,in)       - Internal dof number

c  output arguments:
c     value        (dp,ar(4),in)    - Real,Imag, RealOld,ImagOld
c     ResRdForcOld (int,sc,out)     - Number of returned values (0 or 4)

#include "impcom.inc"
#include "ResRd.inc"

      integer  node, idof, i, k, iNode, ResRdForcOld
      double precision  value(4)

      ResRdForcOld = 0

      if (nFor .gt. 0) then

         iNode = iNinternal(node)
         if (iNode .gt. 0) then
            k = node*32 + idof
            do  i = 1,nFor
               if (k .eq. IndexToData(i)) then
                  ResRdForcOld = 4
                  value(1) = Data(4*(i-1)+1)
                  value(2) = Data(4*(i-1)+2)
                  value(3) = Data(4*(i-1)+3)
                  value(4) = Data(4*(i-1)+4)
                  exit
               endif
            enddo
         endif

      endif

      return
      end

*deck,ResRdForc
      function ResRdForc (node,idof,value)
c primary function:    Retrieve an applied force value

c object/library: ResRd

c  input arguments:
c     node      (int,sc,in)       - External node number
c     idof      (int,sc,in)       - Internal dof number

c  output arguments:
c     value     (dp,ar(4),in)     - Real,Imag, RealOld,ImagOld
c     ResRdForc (int,sc,out)      - Number of returned values (0 or 4)

#include "impcom.inc"
#include "ResRd.inc"

      integer  node, idof, i, iNode, ResRdForc, bcNode, bcDof
      double precision  value(4)

      ResRdForc = 0

      if (nFor .gt. 0) then

         iNode = iNinternal(node)
         if (iNode .gt. 0) then
            do  i = 1,nFor
               bcNode = bcNodes(i)
               bcDof = bcDofs(i)
               if (iNode.eq.bcNode .and. idof.eq.bcDof) then
                  ResRdForc = 4
                  value(1) = Data(4*(i-1)+1)
                  value(2) = Data(4*(i-1)+2)
                  value(3) = Data(4*(i-1)+3)
                  value(4) = Data(4*(i-1)+4)
                  exit
               endif
            enddo
         endif

      endif

      return
      end

*deck,ResRdForcEnd
      subroutine ResRdForcEnd
c primary function:    Finish the applied forces

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fFree


      call fFree (pdDataL)
      call fFree (piIndexToDataL)
      call fFree (pibcNodesL)
      call fFree (pibcDofsL)

      return
      end


*deck,ResRdBCEnd           parallel
      subroutine ResRdBCEnd
c primary function:   End of BC output

c object/library: ResRd

c  input arguments:  none

c  output arguments:  none

#include "impcom.inc"


      return
      end


*deck,ResRdEresBegin
      subroutine ResRdEresBegin
c primary function:    Retrieve index array for element results

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8, largeIntGet, fMallocL, largeIntGetV, fMallocI
      external  fFree
      PTRFTN    fMallocL, fMallocI
      LONGINT   largeIntGet

      integer  i, kbf

      pointer (piWorkL,iWork)
      integer iWork(*)


      plIndexL = fMallocL(nElem)

      if (lLongRes) then
         iLocL = largeIntGet(SolHeader(119),SolHeader(120))
         LocIndexL = iLocL
      else
         iLocL = SolHeader(12)
         LocIndexL = 0
      endif
      if (iLocL .gt. 0) then
         if (lLongRes) then
            i = nElem*intplong
            call binrd8 (NBLK,iLocL,i,iIndexL(1),kbf,iobuf(1))
            call largeIntGetV (nElem,iIndexL(1))
         else
            piWorkL = fMallocI(nElem)
            i = nElem
            call binrd8 (NBLK,iLocL,i,iWork(1),kbf,iobuf(1))
            do i = 1,nElem
               iIndexL(i) = iWork(i)
            enddo
            call fFree (piWorkL)
         endif
      else
         do  i = 1,nElem
            iIndexL(i) = 0
         enddo
      endif

      return
      end


*deck,ResRdEstrBegin
      function ResRdEstrBegin  (iElem)
c primary function:    Retrieve results index for an element

c object/library: ResRd


c  input arguments:
c     iElem    (int,sc,in)        - Element Number

c  output arguments:  none

#include "impcom.inc"
#include "ResRd.inc"


      integer  ResRdEstrBegin, iElem, i, j, kbf

      ResRdEstrBegin = 1

      LocElemIndexL = 0
      j = iEinternal(iElem)
      if (j .gt. 0) then
         if (iIndexL(j) .gt. 0) then
            iLocL = iIndexL(j) + LocIndexL
            if (lLongRes) then
               LocElemIndexL = iLocL
            else
               LocElemIndexL = 0
            endif
            i = 25
            call binrd8 (NBLK,iLocL,i,ElemIndex(1),kbf,iobuf(1))
            ResRdEstrBegin = 0
         endif
      else
         do i = 1,25
            ElemIndex(i) = 0
         enddo
      endif

      return
      end

*deck,ResRdEstr
      function ResRdEstr (iStr,Str)
c primary function:   Get an element's results

c object/library: ResRd


c  input arguments:
c     iStr    (int,sc,in)       - element record number (1-25)

c  output arguments:
c     ResRdEstr (int,sc,out)       - Number of element values
c     Str       (dp,ar(nStr),out)  - element values

#include "impcom.inc"
#include "ResRd.inc"

      external  binrd8

      integer  ResRdEstr, iStr, n, i, kbf
      double precision  Str(*)

      ResRdEstr = 0

      n = ElemIndex(iStr)
      if (n .ne. 0) then
         if (n .lt. 0) then
            n = - n
            do  i = 1,n
               Str(i) = 0.0d0
            enddo
         else
            iLocL = ElemIndex(iStr) + LocElemIndexL
            i = 10
            call binrd8 (NBLK,iLocL,i,Str(1),kbf,iobuf(1))
            n = i/intpdp
         endif
         ResRdEstr = n
      endif


      return
      end

*deck,ResRdEstrEnd     parallel
      subroutine ResRdEstrEnd
c primary function:     End of output for an element

c object/library: ResRd


#include "impcom.inc"

      return
      end


*deck,ResRdEresEnd
      subroutine ResRdEresEnd
c primary function:    End of element stresses

c object/library: ResRd


#include "impcom.inc"
#include "ResRd.inc"

      external  fFree

      call fFree (plIndexL)

      return
      end



*deck,ResRdNstrBegin
      subroutine ResRdNstrBegin (kNodStr)
c primary function:     Read nodal stress vector(s)

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fMallocD, largeIntGet
      PTRFTN    fMallocD
      LONGINT   largeIntGet

      integer   kNodStr, i, j, k, kbf


      kNodStr = 0
      if (lLongRes) then
         iLocL = largeIntGet(SolHeader(145),SolHeader(146))
         kNodStr = SolHeader(143)
         if (kNodStr.gt.0 .and. iLocL.gt.0) then
            i = nNode*6*kNodStr
            pdDataL = fMallocD(i)
            j = 1
            do  k = 1,kNodStr
               i = 6*nNode*intpdp
               call binrd8 (NBLK,iLocL,i,Data(j),kbf,iobuf(1))
               j = j + 6*nNode
            enddo
         endif
      endif

      return
      end


*deck,ResRdNstr
      function ResRdNstr (node,Nstr)
c primary function:    Retrieve nodal stresses

c object/library: ResRd

c  input arguments:
c     node     (int,sc,in)            - Node number

c  output arguments:  none
c     Nstr     (dp,ar(6,kNodstr),out) - Nodal stresses
c     ResRdNstr(int,sc,out)           - kNodStr

#include "impcom.inc"
#include "ResRd.inc"

      integer  ResRdNstr, node, kNodStr, i, j, n
      double precision  Nstr(6,*)

      ResRdNstr = 0
      
      kNodStr = SolHeader(143)
      if (kNodStr .gt. 0) then
         if (iNinternal(node) .gt. 0) then
            j = (iNinternal(node)-1)*6
            do  n = 1,kNodStr
               do  i = 1,6
                  Nstr(i,n) = Data(j+i)
               enddo
               j = j + nNode*6
            enddo
            ResRdNstr = kNodStr
         endif
      endif

      return
      end


*deck,ResRdNstrEnd
      subroutine ResRdNstrEnd
c primary function:    Finish displacements

c object/library: ResRd

#include "impcom.inc"
#include "ResRd.inc"

      external  fFree

      call fFree (pdDataL)

      return
      end


*deck,ResRdSolEnd      parallel
      subroutine ResRdSolEnd
c primary function:    End of a Solution Result Set

c object/library: ResRd


#include "impcom.inc"

      return
      end


*deck,ResRdEnd
      subroutine ResRdEnd
c primary function:    Close result file

c object/library: ResRd


#include "impcom.inc"
#include "ResRd.inc"

      external  fFree
      external  binclo

      integer  i

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

      integer  iLeng
      PTRFTN   fMallocI
#if defined(FUJITSU_SYS)
      integer*8  iLeng8
#endif

#if defined(RS6000_SYS) || defined(RS64K_SYS)
      fMallocI = malloc(%val(iLeng*4))
#elif defined(FUJITSU_SYS)
      iLeng8 = iLeng*4*1_8
      fMallocI = malloc(iLeng8)
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

      integer  dLeng
      PTRFTN   fMallocD
#if defined(FUJITSU_SYS)
      integer*8  dLeng8
#endif

#if defined(RS6000_SYS) || defined(RS64K_SYS)
      fMallocD = malloc(%val(dLeng*8))
#elif defined(FUJITSU_SYS)
      dLeng8 = dLeng*8*1_8
      fMallocD = malloc(dLeng8)
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

      integer  lLeng
      PTRFTN   fMallocL
#if defined(FUJITSU_SYS)
      integer*8  lLeng8
#endif

#if defined(RS6000_SYS) || defined(RS64K_SYS)
      fMallocL = malloc(%val(lLeng*8))
#elif defined(FUJITSU_SYS)
      lLeng8 = lLeng*8*1_8
      fMallocL = malloc(lLeng8)
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

      if (piPointerL.ne.PTRFTNNULL) call free (piPointerL)
      piPointerL = PTRFTNNULL

      return
      end
