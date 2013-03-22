*deck,ResRdDemo    USERDISTRIB  User                                    ansys
      program ResRdDemo
c  **********  test program for ResRd.F routines  *****
#include "impcom.inc"

      external  ResRdBegin, ResRdGeomBegin, ResRdGeomEnd,
     x          ResRdTypeBegin, ResRdType, ResRdTypeEnd,
     x          ResRdRealBegin, ResRdReal, ResRdRealEnd,
     x          ResRdCsysBegin, ResRdCsys, ResRdCsysEnd,
     x          ResRdNodeBegin, ResRdNode, ResRdNodeEnd,
     x          ResRdElemBegin, ResRdElem, ResRdElemEnd,
     x          ResRdSolBegin , ResRdSolEnd, ResRdEnd,
     x          ResRdDispBegin, ResRdDisp, ResRdDispEnd,
     x          ResRdRforBegin, ResRdRfor, ResRdRforEnd,
     x          ResRdEresBegin, ResRdEresEnd,
     x          ResRdEstrBegin, ResRdEstr, ResRdEstrEnd,
     x          ResRdBCBegin, ResRdBCEnd,
     x          ResRdFixBegin, ResRdFix, ResRdFixOld, ResRdFixEnd,
     x          ResRdForcBegin, ResRdForc, ResRdForcOld, ResRdForcEnd,
     x          ResRdNstrBegin, ResRdNstr, ResRdNstrEnd
      integer   ResRdBegin, ResRdType, ResRdReal, ResRdCsys,
     x          ResRdNode,  ResRdElem, ResRdSolBegin,
     x          ResRdDisp,  ResRdRfor, ResRdEstr, ResRdEstrBegin,
     x          ResRdFix,   ResRdFixOld, ResRdForc, ResRdForcOld,
     x          ResRdNstr

c *** define max array sizes
c ***     !!! NOTE: may need to be increased for other models !!!
      integer  MAXDOFSPERNODE, MAXNODESPERELEM, MAXREALS, MAXELERECORD
      parameter (MAXDOFSPERNODE=32, MAXNODESPERELEM=24, MAXREALS=100,
     x           MAXELERECORD=MAXNODESPERELEM*11)

      integer  Nunit,Lunit,Units,NumDof,DOF(MAXDOFSPERNODE),UserCode,
     x         MaxNode, NumNode, MaxElem, NumElem,MaxResultSet,
     x         MaxType, MaxReal, MaxCsys, NumType, NumReal,
     x         NumCsys, NumPerReal, NumResultSet, nRForce, nFixed,
     x         nForces, kNodStr
      integer  ielc(200), i, j, k, n, kerr, nodes(MAXNODESPERELEM),
     x         ElemData(10), iSet, ncumit, lstep, substep, kcmplx, ii,
     x         BCHeader(40), nXYZnode
      character*4   DofLab(MAXDOFSPERNODE)
      character*32  JobName
      character*12  Fname
      character*80  Title(5)
      double precision  Rcon(MAXREALS), xyz(6), time,
     x                  Disp(MAXDOFSPERNODE), Str(MAXELERECORD),
     x                  Value, values(4), NodStr(6,3)
      character*3   EstrLab(25)
      data EstrLab /'EMS','ENF','ENS','ENG','EGR',
     x              'EEL','EPL','ECR','ETH','EUL',
     x              'EFX','ELF','EMN','ECD','ENL',
     x              'EHC','EPT','ESF','UNK','ETB',
     x              'ECT','EXY','EBA','ESV','UNK'/

      Nunit = 12
      Lunit = 6
      Fname = 'test.rst'

      write (*,3000) Nunit,Lunit,Fname
 3000 format (/'  *****  ResRdDemo  *****'/
     x         ' Nunit=',i3,' Lunit=',i3,'  Fname=',a)

      kerr = ResRdBegin (Nunit,Lunit,Fname,8,Title(1),JobName,
     x                 Units,NumDof,DOF(1),UserCode,
     x                 MaxNode, NumNode, MaxElem, NumElem,
     x                 MaxResultSet,NumResultSet)
      if (kerr .ne. 0) go to 999
      if (NumDof.gt.MAXDOFSPERNODE) then
         write (*,9100) NumDof,MAXDOFSPERNODE
 9100    format (/' *** ERROR ***'/' Number of DOFs per node on the ',
     x     'file exceeds the currently dimensioned size.'/
     x     '   DOFs per node on the file               =',i5/
     x     '   DOFs per node dimensioned in ResRdDemo.F=',i5)
         go to 999
      endif

      write (*,3001) Title(1),Title(2)
 3001 format (' *** Header'/' Title=',a/' Stitle1=',a)
      write (*,3002) JobName, Units, NumDof
 3002 format (' JobName=',a,'  Units=',i3,'  NumDof=',i4)
      write (*,3003) (DOF(i),i=1,NumDof)
 3003 format (' DOF=',10i4)
      write (*,3004) UserCode,MaxResultSet
 3004 format (' UserCode=',i8,'  MaxResultSet=',i6)
      write (*,3005) MaxNode,NumNode,MaxElem,NumElem
 3005 format (' MaxNode=',i6,'  NumNode=',i6,'  MaxElem=',i6,
     x        '  NumElem=',i6)


c          **********  Geometry  **********
      call ResRdGeomBegin (MaxType, MaxReal, MaxCsys, nXYZnode)

      write (*,3011) MaxType,MaxReal,MaxCsys,nXYZnode
 3011 format (/'  **********  Geometry  **********'/
     x  ' MaxType=',i4,'  MaxReal=',i4,'  MaxCsys=',i4,'  nXYZnode=',i4)


c          *****  Element Types  *****
      call ResRdTypeBegin (NumType)
      write (*,3012) NumType
 3012 format (/' *** Element Types'/' NumType=',i4)

      do  i = 1,MaxType
         n = ResRdType(i,ielc(1))
         if (n .gt. 0) then
            write (*,3013) i,n,ielc(1),ielc(2)
 3013       format ('  iType=',i4,'  n=',i4,'  ielc(1:2)=',2i4)
         endif
      enddo

      call ResRdTypeEnd


c          *****  Real Constants  *****
      call ResRdRealBegin (NumReal,NumPerReal)
      write (*,3022) NumReal, NumPerReal
 3022 format (/' *** Real Constants'/' NumReal=',i4,'  NumPerReal=',i4)

      do  i = 1,MaxReal
         n = ResRdReal(i,Rcon)
         if (n .gt. 0) then
            if (n.gt.MAXREALS) then
               write (*,9200) n,MAXREALS
 9200          format (/' *** ERROR ***'/' Number of Real Constants ',
     x           'on the file exceeds the currently dimensioned size.'/
     x           '   Number of Reals on the file               =',i6/
     x           '   Number of Reals dimensioned in ResRdDemo.F=',i6)
               go to 999
            endif
            write (*,3023) i,n,(Rcon(j),j=1,n)
 3023       format ('  iReal=',i4,'  n=',i4/' Rcon=',6g12.5)
         endif
      enddo

      call ResRdRealEnd


c          *****  Coordinate Systems  *****
      call ResRdCsysBegin (NumCsys)
      write (*,3032) NumCsys
 3032 format (/' *** Coordinate Systems'/' NumCsys=',i4)

      do  i = 1,MaxCsys
         n = ResRdCsys(i,Rcon)
         if (n .gt. 0) then
            write (*,3033) i,n,(Rcon(j),j=1,n)
 3033       format ('  iCsys=',i4,'  n=',i4/(' Data=',6g12.5))
         endif
      enddo

      call ResRdCsysEnd


c          *****  Nodal Coordinates  *****
      call ResRdNodeBegin

      write (*,3040)
 3040 format (/' *** Nodes')
      do  i = 1,nXYZnode
         n = ResRdNode(i,xyz(1))
         if (n .gt. 0) then
            write (*,3041) n,xyz
 3041       format (' Node=',i4,'  xyz=',3g12.5,3x,3f8.2)
         endif
      enddo

      call ResRdNodeEnd


c          **********  Elements  **********
      call ResRdElemBegin

      write (*,3050)
 3050 format (/' *** Elements')
      do  i = 1,MaxElem
         n = ResRdElem(i,nodes,ElemData)
         if (n .gt. 0) then
            if (n.gt.MAXNODESPERELEM) then
               write (*,9300) n,MAXNODESPERELEM
 9300          format (/' *** ERROR ***'/' Number of nodes per element',
     x           ' on the file exceeds the currently dimensioned size.'/
     x           '   Number of nodes on the file               =',i6/
     x           '   Number of nodes dimensioned in ResRdDemo.F=',i6)
               go to 999
            endif
            write (*,3051) i,n,ElemData
 3051       format (' iElem=',i4,'  nNode=',i4,'  ElemData=',10i4)
            write (*,3052) (nodes(j),j=1,n)
 3052       format ('   Nodes=',10i6)
         endif
      enddo

      call ResRdElemEnd


      call ResRdGeomEnd



c          **********  Load Sets  **********
      do  iSet = 1,NumResultSet

         ncumit = iSet
         kcmplx = 0
         kerr = ResRdSolBegin (2,lstep,substep,ncumit,kcmplx,time,
     x                         Title(1),DofLab(1))
         if (kerr.ne.0) then
            write (*,9250) iSet
 9250       format (/' *** WARNING ***'/' Load set=',i6,' not found',
     x        ' on the file.')
            cycle
         endif
         write (*,3061) iSet, lstep,substep,ncumit,kcmplx,time
 3061    format (/'  **********  Load Set =',i3,'  **********'/
     x            ' lstep=',i3,'  substep=',i3,'  ncumit=',i3,
     x            '  kcmplx=',i2,'  time=',g12.5)
         write (*,3062) Title(1),Title(2)
 3062    format (' Title=',a/' Stitle1=',a)
         write (*,3063) (DofLab(i),i=1,NumDof)
 3063    format (' DofLab=',6a6)


c          *****  Displacements  *****
         call ResRdDispBegin
         write (*,3071)
 3071    format (/'  *****  Displacements  *****')

         do  i = 1,MaxNode
            k = ResRdDisp(i,Disp(1))
            if (k .gt. 0) then
               write (*,3073) i,(Disp(j),j=1,k)
 3073          format (' Node=',i4,'  Disp=',6g12.5/(10x,6g12.5))
            endif
         enddo
         call ResRdDispEnd


c          ******  Reaction Forces  *****
         call ResRdRforBegin (nRForce)
         write (*,3075)
 3075    format (/'     *****  Reaction Forces  *****')

         if (nRForce .gt. 0) then

            do  i = 1,MaxNode
               do  j = 1,NumDof
                  k = ResRdRfor (i,j,Value)
                  if (k .gt. 0) then
                     write (*,3076) i,j,Value
 3076                format (' Node=',i6,'  DOF=',i3,'  Value=',g12.5)
                  endif
               enddo
            enddo

            call ResRdRforEnd
         endif


c          ******  Boundary Conditions  *****
         call ResRdBCBegin (BCHeader(1))

c          ******  Constraints  *****
         call ResRdFixBegin (BCHeader(1),nFixed)
         write (*,3077)
 3077    format (/'     *****  Constraints  *****')

         if (nFixed .gt. 0) then

            do  i = 1,MaxNode
               do  j = 1,NumDof
                  if (BCHeader(7) .eq. 0) then
                     k = ResRdFixOld (i,j,values(1))
                  else
                     k = ResRdFix (i,j,values(1))
                  endif
                  if (k .gt. 0) then
                     write (*,3078) i,j,values
 3078                format (' Node=',i6,'  DOF=',i3,'  Values=',4g12.5)
                  endif
               enddo
            enddo

            call ResRdFixEnd
         endif

c          ******  Applied Forces  *****
         call ResRdForcBegin (BCHeader(1),nForces)
         write (*,3079)
 3079    format (/'     *****  Applied Forces  *****')

         if (nForces .gt. 0) then

            do  i = 1,MaxNode
               do  j = 1,NumDof
                  if (BCHeader(7) .eq. 0) then
                     k = ResRdForcOld (i,j,values(1))
                  else
                     k = ResRdForc (i,j,values(1))
                  endif
                  if (k .gt. 0) then
                     write (*,3078) i,j,values
                  endif
               enddo
            enddo

            call ResRdForcEnd
         endif

         call ResRdBCEnd


c          *****  Element Stresses  *****
         call ResRdEresBegin

         write (*,3081)
 3081    format (/'     *****  Element Stresses  *****')
         do  i = 1,MaxElem
            kerr = ResRdEstrBegin(i)
            if (kerr .eq. 0) then
               write (*,3082) i
 3082          format (/' Element=',i4)
               do  j = 1,25
                  k = ResRdEstr(j,Str(1))
                  if (k .gt. 0) then
                     if (k.gt.MAXELERECORD) then
                        write (*,9400) EstrLab(j),k,MAXELERECORD
 9400                   format (/' *** ERROR ***'/' Number of element ',
     x                    'results on the file exceeds the currently ',
     x                    'dimensioned size'/' for record type ',a/
     x              '   Number of items on the file               =',i6/
     x              '   Number of items dimensioned in ResRdDemo.F=',i6)
                        go to 999
                     endif
                     write (*,3083) EstrLab(j),(Str(ii),ii=1,k)
 3083                format (2x,a,'=',6g12.5/(6x,6g12.5))
                  endif
               enddo
               call ResRdEstrEnd
            endif
         enddo

         call ResRdEresEnd


c          **********  Nodal Stresses  **********
         call ResRdNstrBegin (kNodStr)

         if (kNodStr .gt. 0) then
            write (*,3091) kNodStr
 3091       format (/'     *****  Nodal Stresses  *****  kNodStr=',i2)
            do  i = 1,MaxNode
               n = ResRdNstr(i,NodStr(1,1))
               if (n .gt. 0) then
                  do  j = 1,n
                     write (*,3092) i,(NodStr(k,j),k=1,6)
 3092                format (' Node=',i7,'  Str=',6f10.0)
                  enddo
               endif
            enddo
         endif

         call ResRdNstrEnd

         call ResRdSolEnd

      enddo

      call ResRdEnd

 999  continue
      end
