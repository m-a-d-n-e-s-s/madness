      subroutine fci_print_vector(c, thresh, nwalko, nactive, statesym)
      implicit none
#include "fciparam.fh"
c      
c     Print large components of the CI vector out
c
      double precision c(*)   ! [input] The determinant CI vector
      double precision thresh ! [input] Print coeffs > thresh
      integer nwalko          ! [input] No. of orbital occupations
      integer nactive         ! [input] No. of active orbitals
      integer statesym        ! [input] Symmetry of the CI state (0...7)
c
      integer iwalk, mwalk, ns, ioff, numdet, i, testo, testm, sym
      integer iocc(nactmx), mocc(nsmax)
      character*80 fmt
      character*2 spin(2)
      data spin/'+ ','- '/
c
      write(6,1)
 1    format(/
     $     '    Index     Owalk Mwalk      Coeff       Orbital and',
     $     ' Spin occupation'/
     $     ' --------    ------ -----    ---------     -----------',
     $     '----------------')
      do iwalk = 1, nwalko
         call fci_owalk_info(iwalk, iocc, ns, ioff, numdet,sym)
         if (sym .eq. statesym) then
            write(fmt,12) nactive, max(ns,1)
 12         format('(1x,i8,4x,i6,2x,i4,4x,f9.6,4x,',i2,
     $           'i2,4x,',i2,'a2)')
            do mwalk = 1, numdet
               if (abs(c(ioff+mwalk)).gt.thresh) then
                  call fci_mfind_occ(mwalk, mocc, ns)
                  write(6,fmt) ioff+mwalk,iwalk,mwalk,
     $                 c(ioff+mwalk),
     $                 (iocc(i)-1,i=1,nactive),
     $                 (spin(mocc(i)),i=1,ns)
c     
c     This just to check that the det_to_walk routine is working
c     
                  call fci_det_to_walk(ioff+mwalk, testo, testm)
                  if (testo.ne.iwalk .or. testm.ne.mwalk)
     $                 call fci_err('fci_det_to_walk failed',0)
               endif
            enddo
         endif
      enddo
c
      write(6,*)
c
      end
