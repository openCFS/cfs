*deck,ResWrDemo                  User
      program ResWrDemo
c  **********  test program for ResWr.F routines  *****
#include "impcom.inc"

      external  ResWrBegin    , ResWrEnd    ,
     x          ResWrGeomBegin, ResWrGeomEnd,
     x          ResWrTypeBegin, ResWrType, ResWrTypeEnd,
     x          ResWrRealBegin, ResWrReal, ResWrRealEnd,
     x          ResWrCsysBegin, ResWrCsys, ResWrCsysEnd,
     x          ResWrNodeBegin, ResWrNode, ResWrNodeEnd,
     x          ResWrElemBegin, ResWrElem, ResWrElemEnd,
     x          ResWrSolBegin , ResWrSolEnd,
     x          ResWrRforBegin, ResWrRfor, ResWrRforEnd, 
     x          ResWrBCBegin  , ResWrBCEnd,
     x          ResWrFixBegin , ResWrFix , ResWrFixEnd, 
     x          ResWrForcBegin, ResWrForc, ResWrForcEnd, 
     x          ResWrDispBegin, ResWrDisp, ResWrDispEnd,
     x          ResWrEresBegin, ResWrEresEnd,
     x          ResWrEstrBegin, ResWrEstr, ResWrEstrEnd
      integer   ResWrBegin

      integer  Nunit, Lunit, ncFname, Units, UserCode, MaxResultSet
      integer  i, j, kerr, iElem, iSet, kcmplx, ii, iNode
      double precision value(4)
      character*8   JobName
      character*12  Fname

      data  Nunit,Lunit,ncFname,Units,UserCode,MaxResultSet
     x     /  12,   6,     9,     4,    10,        10      /
      data  Fname /'ResWr.rst   '/
      character*3   EstrLab(25)
      data EstrLab /'EMS','ENF','ENS','ENG','EGR',
     x              'EEL','EPL','ECR','ETH','EUL',
     x              'EFX','ELF','EMN','ECD','ENL',
     x              'EHC','EPT','ESF','UNK','ETB',
     x              'ECT','EXY','EBA','ESV','UNK'/

c========================  Example Data ======================================
c   /batch,list
c   /com     EXAMPLE PROBLEM TO BE USED - 2 ELEMENT CANTILEVER BEAM USING PLANE42
c   /com
c   /com             11        12        13
c   /com     FIXED  >*.........*.........*
c   /com             .         .         .
c   /com             .         .         .
c   /com             .    1    .    3    .
c   /com             .         .         .
c   /com     FIXED  >*.........*.........* 3
c   /com             1         2         ! 
c   /com                                 !  
c   /com                                 V  10 #
c   
c   /BATCH,list
c   /config,nres,10       ! to limit size
c   /title, 2 element cantilever beam using plane42
c   /stitle,1,Test case for writing binary files
c   /units,bin
c   /prep7
c   et,1,42
c   mp,ex,1,3e7
c   n,1,0,0
c   n,2,1,0
c   n,3,2,0
c   ngen,2,10,1,3,1,0,1
c   en,1,1,2,12,11
c   en,3,2,3,13,12
c   finish
c   /solution
c   d,1,all
c   d,11,all
c   f,3,fy,-10
c   outres,all,none
c   outres,nsol,1
c   outres,rsol,1
c   outres,strs,1
c   outres,nload,1
c   outpr,all,none
c   outpr,nsol,all
c   outpr,esol,all
c   outpr,nload,all
c   solve
c   finish
c   /post1
c   set,1,1
c   nlist,all
c   elist,all
c   etlist,all
c   rlist,all
c   flist,all
c   dlist,all
c   prnsol,dof
c   presol,s
c   presol,forc
c   prrsol,
c   finish
c   /aux02
c   file,file,rst
c   form,long
c   dump,all
c   finish
c   /exit,nosave

      integer    NSETS,    MAXNODE,    MAXELEM,   NUMNODE,   NUMELEM
      parameter (NSETS=1,  MAXNODE=13, MAXELEM=3, NUMNODE=6, NUMELEM=2)
      integer    NUMDOF,   NFIXED,     NFORCE
      parameter (NUMDOF=2, NFIXED=4,   NFORCE=1)
      integer    MAXTYPE,  MAXREAL,    MAXCSYS,   NRCON,     NPEREL
      parameter (MAXTYPE=1,MAXREAL=0,  MAXCSYS=0, NRCON=1,   NPEREL=4)

      integer  LSET(NSETS), ielc(150),
     x         DOF(NUMDOF), SUBST(NSETS), NCUM(NSETS),
     x         NumType, NumReal, NumCsys,
     x         nodeRF(NFIXED), dofRF(NFIXED),
     x         nodeF(NFORCE), dofF(NFORCE),
     x         Norder(NUMNODE), Eorder(NUMELEM), Enode(NPEREL,MAXELEM),
     x         ElemData(10), NumResultSet
      double precision  Nxyz(6,MAXNODE), Time(NSETS),
     x                  Disp(NUMDOF,NUMNODE),Estress(11,NPEREL,NUMELEM),
     x                  Eforce(NUMDOF,NPEREL,NUMELEM), Etemps(NPEREL),
     x                  Rforce(NFIXED), Dvalue(NFIXED), Fvalue(NFORCE),
     x                  Rcon(NRCON), Csys(24), Area(NUMELEM),
     x                  Energy(NUMELEM), Evalue(11)
      character*4   DofLab(NUMDOF)
      character*80  Title(5)

      data  Title / 'CFS++ GFortran standalone ANSYS RST test.',
     x              'Test case for writing binary files',
     x              3*'       '/

c --- Geometry
      data  NumType, NumReal, NumCsys
     x     /    1,       0,       0  /

c      see elccmt.inc and elcpar.inc for information on the following
c      This was optained with an /aux02 file dump of the test input data
      data  ielc / 
     x            1,42,0,0,0, 0,0,0,0,0,
     x            0,0,0,0,0, 0,0,0,0,0,
     x            0,4,0,0,0, 1,0,1,0,0,
     x            0,0,3,3,1, 0,0,0,0,0,
     x            0,0,0,0,0, 0,1,0,0,4,
     x            0,4,0,0,0, 2,0,0,0,0,
     x            4,4,4,2,4, 22,4,0,0,0,
     x            0,0,0,2,4, 0,0,0,0,0,
     x            0,0,4,4,0, 0,0,0,0,0,
     x            4,0,0,4,1, 0,0,0,0,0,
     x            0,0,0,0,0, 1,4,0,8,25,
     x            12,139,147,0,0, 0,1,0,0,0,
     x            0,122,0,0,0, 152,24,108,0,0,
     x            0,0,0,1347174734,1161048608, 0,0,0,1,0,
     x            55,0,0,0,1, 0,0,0,0,0/

      data  DOF    / 1, 2 /
      data  DofLab / 'UX  ', 'SP01' /
  
      data  Rcon   / 1.0d0 /

      data Csys    / 24*0.0d0 /

      data  Norder  /  1,  2, 12, 11, 3, 13/     ! Node equivalence table
      data  Eorder  /  1,  3 /                   ! Element equivalence table

      data  Nxyz   /  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  ! Node 1 
     x                1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  ! Node 2
     x                2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  ! Node 3
     x             42*0.0d0,
     x                0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  ! Node 11
     x                1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  ! Node 12
     x                2.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 / ! Node 13

      data  Enode  /  1,  2, 12, 11,             ! Element nodes
     x                0,  0,  0,  0,
     x                2,  3, 13, 12 /

      data ElemData / 1,        ! material number
     x                1,        ! type number
     x                1,        ! real set number
     x                1,        ! section number
     x                0,        ! element coordinate system
     x                0,        ! birth-death flag  0, alive  1, dead
     x                0,        ! solid model reference
     x                0,        ! coded shape key
     x                0,        ! p-method exclude key
     x                0 /       ! not used


c --- Results
      data NumResultSet / 1 /

      data  Time  / 1.0d0 /

      data  LSET  / 1 /
      data  SUBST / 1 /
      data  NCUM  / 1 /

c --- Reactions
      data  nodeRF / 1, 1, 11, 11 /
      data  dofRF  / 1, 2,  1, 2 /
c                      FX        FY
      data  Rforce / 20.0d0, 5.7142858d0,           ! Node 1
     x              -20.0d0, 4.2857142d0 /          ! Node 11

c --- Constraints
      data  Dvalue /  0.0d0, 0.0d0, 0.0d0, 0.0d0 /

c --- Applied forces
      data  nodeF  / 3 /
      data  dofF   / 2 /
      data  Fvalue / -10.0d0 /

c --- Displacements
c           UX            UY
      data  Disp /
     x   0.000000   ,  0.000000,           ! Node 1
     x  -0.297857D-5, -0.379524D-5,        ! Node 2
     x   0.302143D-5, -0.393810D-5,        ! Node 12
     x   0.000000   ,  0.000000,           ! Node 11
     x  -0.404286D-5, -0.120190D-4,        ! Node 3
     x   0.395714D-5, -0.114476D-4 /       ! Node 13

c --- Element temperatures
      data  Etemps / 4*0.0d0 /

c --- Element nodal forces
      data  Eforce /
     x -20.0000,  -5.71429, 10.0000,  3.57143,  ! Element 1
     x -10.0000,   6.42857, 20.0000, -4.28571,

     x -10.0000,  -3.57143,  0.0000, 10.00000,  ! Element 3
     x   0.0000,   0.00000, 10.0000, -6.42857 /

c --- Element nodal stresses
      data  Estress /
     x -90.0000,  0.0000,   0.0000, -10.0000,   0.0000, 0.0000,  ! Element 1
     x   1.0977,  0.0000, -91.1954,  92.1954,  91.6515,
     x -90.0000, -4.2857,   0.0000, -10.0000,   0.0000, 0.0000,
     x   0.0000, -3.13451,-91.1512,  91.1512,  89.6251,
     x  90.0000, -4.28571,  0.0000, -10.0000,   0.0000, 0.0000,
     x  91.0489,  0.0000,  -5.3346,  96.3836,  93.8301,
     x  90.0000,  0.0000,   0.0000, -10.0000,   0.0000, 0.0000,
     x  91.0977,  0.0000,  -1.0977,  92.1954,  91.6515,

     x -30.0000, -4.2857,   0.0000, -10.0000,   0.0000, 0.0000,  ! Element 3
     x   0.0000, -0.8546, -33.4310,  33.4311,  33.0121,
     x -30.0000, 17.1429,   0.0000, -10.0000,   0.0000, 0.0000,
     x  19.1764,  0.0000, -32.0335,  51.2099,  44.8125,
     x  30.0000, 17.1429,   0.0000, -10.0000,   0.0000, 0.0000,
     x  35.4595, 11.6833,   0.0000,  35.4595,  31.2984,
     x  30.0000, -4.2857,   0.0000, -10.0000,   0.0000, 0.0000,
     x  32.7035,  0.0000,  -6.9892,  39.6927,  36.7007 /

c --- Element energies
      data  Area /  1.0d0, 1.0d0 /
      data  Energy / 0.494354d-4, 0.106599d-4 /
      data  Evalue / 0.0d0, 0.0d0, TINY, TINY, 0.0d0, 0.0d0, 5*TINY/
c======================== End Example Data =====================================


      JobName  = 'ResWDemo'


      write (*,3000) Nunit,Lunit,Fname
 3000 format (/'  *****  ResWrDemo  *****'/
     x         ' Nunit=',i3,' Lunit=',i3,'  Fname=',a)
      write (*,3001) Title(1),Title(2)
 3001 format (' Title=',a/' Stitle1=',a)
      write (*,3002) JobName, Units, NUMDOF
 3002 format (' JobName=',a,'  Units=',i3,'  NumDof=',i4)
      write (*,3003) (DOF(i),i=1,NUMDOF)
 3003 format (' DOF=',10i4)
      write (*,3004) UserCode,MaxResultSet
 3004 format (' UserCode=',i8,'  MaxResultSet=',i6)
      write (*,3005) MAXNODE,NUMNODE,MAXELEM,NUMELEM
 3005 format (' MaxNode=',i6,'  NumNode=',i6,'  MaxElem=',i6,
     x        '  NumElem=',i6)


      kerr = ResWrBegin (Nunit,Lunit,Fname,ncFname,Title(1),JobName,
     x                   Units,NUMDOF,DOF(1),UserCode,
     x                   MAXNODE,NUMNODE,MAXELEM,NUMELEM,
     x                   MaxResultSet)

      if (kerr .ne. 0) go to 999


c          **********  Geometry  **********

      write (*,3011) MAXTYPE,MAXREAL,MAXCSYS
 3011 format (/'  *****  Geometry  *****'/
     x        ' MaxType=',i4,'  MaxReal=',i4,'  MaxCsys=',i4)
      write (*,3012) NumType,NumReal,NumCsys
 3012 format (' NumType=',i4,'  NumReal=',i4,'  NumCsys=',i4/)

      NumReal = 0
      NumCsys = 0

      call ResWrGeomBegin (MAXTYPE,NumType,MAXREAL,NumReal,
     x                     MAXCSYS,NumCsys)


c          *****  element types  *****
      call ResWrTypeBegin

      do  i = 1,NumType
         write (*,3013) i,ielc(1),ielc(2)
 3013    format (' iType=',i4,'  ielc(1:2)=',2i4)
         call ResWrType(i,ielc(1))
      enddo

      call ResWrTypeEnd


c          *****  real constants  *****
c      call ResWrRealBegin

c      do  i = 1,NumReal
c         write (*,3023) i,NRCON,(Rcon(j),j=1,NRCON)
c 3023    format (' iReal=',i4,'  n=',i4/' Rcon=',6g12.5)
c         call ResWrReal (i,NRCON,Rcon(1))
c      enddo

c      call ResWrRealEnd


c          *****  Coordinate systems  *****
c      call ResWrCsysBegin
c
c      do  i = 1,NumCsys
c         write (*,3033) i,(Csys(j),j=1,24)
c 3033    format (' iCsys=',i4/(' Data=',6g12.5))
c         call ResWrCsys (i,Csys)
c      enddo
c
c      call ResWrCsysEnd


c          *****  nodal coordinates  *****
      call ResWrNodeBegin

      write (*,3040)
 3040 format (/'  Node locations  ')
      do  i = 1,NUMNODE
         iNode = Norder(i)
         write (*,3041) iNode,(Nxyz(j,iNode),j=1,6)
 3041    format (' Node=',i4,'  xyz=',6f8.2)
         call ResWrNode(iNode,Nxyz(1,iNode))
      enddo

      call ResWrNodeEnd


c          **********  Elements  **********
      call ResWrElemBegin

      write (*,3050)
 3050 format (/'  Elements  ')
      do  i = 1,NUMELEM
         iElem = Eorder(i)
         write (*,3051) iElem,NPEREL,ElemData
 3051    format (' iElem=',i4,'  nNode=',i4,'  ElemData=',10i4)
         write (*,3052) (Enode(j,iElem),j=1,NPEREL)
 3052    format (' Nodes=',10i6)
         call ResWrElem(iElem,NPEREL,Enode(1,iElem),ElemData(1))
      enddo

      call ResWrElemEnd


      call ResWrGeomEnd



      call ResWrEnd

 999  continue
      end
