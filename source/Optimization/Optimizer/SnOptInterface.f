*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     snfilewrapper   snopenappend   snclose   snopenread
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snfilewrapper
     &   ( name, ispec, inform, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     snSpec
      character*(*)
     &     name
      integer
     &     inform, lencw, leniw, lenrw, iw(leniw), ispec
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8
*     ==================================================================
*     Read options for snopt from the file named name. inform .eq 0 if
*     successful.
*
*     09 Jan 2000: First version of snfilewrapper.
*     ==================================================================
      integer
     &     iostat

      integer ln 
      open( ispec, iostat=iostat, file=name, status='old' )
      if ( 0 .ne. iostat ) then
         inform = 2 + iostat
      else
         call snSpec( ispec, inform, cw, lencw, iw, leniw, rw, lenrw )
         close( ispec )
      end if

      end ! subroutine snfilewrapper

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snopenappend
     &   ( iunit, name, inform )

      implicit
     &     none
      integer
     &     iunit
      character*(*)
     &     name
      integer
     &     inform
*     ==================================================================
*     Open file named name to Fortran unit iunit. inform .eq. 0 if
*     sucessful.  Although opening for appending is not in the f 77
*     standard, it is understood by f2c.
*
*     09 Jan 2000: First version of snoptopenappend
*     ==================================================================
      open( iunit, iostat=inform, file=name, access='append' )

      end ! subroutine snoptopenappend

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snclose
     &   ( iunit )
      implicit
     &     none
      integer
     &     iunit
*     ==================================================================
*     Close unit iunit.
*
*     09 Jan 2000: First version of snoptclose
*     ==================================================================
      close( iunit )

      end ! subroutine snoptclose

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snopenread
     &   ( iunit, name, inform )

      implicit
     &     none
      integer
     &     iunit
      character*(*)
     &     name
      integer
     &     inform
*     ==================================================================
*     Open file called name to Fortran unit iunit. inform .eq. 0 if
*     sucessful.  Although opening for appending is not in the f77
*     standard, it is understood by f2c.
*
*     09 Jan 2000: First version of npopenread
*     ==================================================================
      open ( iunit, iostat=inform, file=name, form = 'FORMATTED',
     *     status = 'OLD' )

      end ! subroutine snopenread
