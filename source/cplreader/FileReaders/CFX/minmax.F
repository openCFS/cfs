      subroutine minmax(text,nvx,ncompt,var)
c
      character*(*) text
      integer nvx, ncompt, ivx, ic
      real    var(ncompt,nvx), varmin(3), varmax(3)
c
      ivx = 1
      do ic=1,ncompt
        varmin(ic) = var(ic,ivx)
        varmax(ic) = var(ic,ivx)
      enddo
c
      do ivx=2,nvx
        do ic=1,ncompt
          varmin(ic) = min( var(ic,ivx), varmin(ic) )
          varmax(ic) = max( var(ic,ivx), varmax(ic) )
        enddo
      enddo
c
      write (*,100) text
      do ic=1,ncompt
        write (*,110) varmin(ic), varmax(ic)
      enddo
c
  100 format(/,' Variable: ',A)
  110 format('   Minimum:',1pe15.3,10x,'Maximum:',1pe15.3)
c
      end
