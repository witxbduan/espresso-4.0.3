  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                   
  
module tetrahedron
USE kinds, ONLY : DP
implicit none
type point
real(DP) :: w(10),c(10),d(10),e(10)
integer :: i
end type point

type tetra
type(point) :: p(4)
end type tetra

end module tetrahedron
                                           
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_elec_bo
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the electron self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the electron linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
use tetrahedron
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  USE phcom, ONLY : lgamma, nmodes
  USE epwcom, ONLY : nbndsub, lrepmatf, iunepmatf, &
      fsthick, eptemp, ngaussw, degaussw, iuetf,   &
      nbndskip, ecutse, parallel_k, &
      parallel_q, epf_mem, etf_mem, eig_read, eps_acustic, &
nqf1, nqf2, nqf3 !! extra
  USE pwcom, ONLY : nelec, ef, isk
  USE el_phon, ONLY : etf, ibndmin, ibndmax, nksf, etfq, &
      epf17, wkf, nksqf, nxqf, wf, wqf, xkf, nkstotf, xqf
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : me_pool, inter_pool_comm, my_pool_id
#endif
  implicit none
  !
  
  integer :: ncorn, ntet, ierr, itet, jtet
type(tetra), allocatable :: tet(:)
real(kind=DP), allocatable :: wkt(:), eprim(:,:,:), Fqnu(:,:,:,:), ekj(:,:)
  
  REAL(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq,  &
       weight, wgauss, dosef, dos_ef, sigmar(nbndsub, nksf), &
       sigmai(nbndsub, nksf), zi(nbndsub, nksf), eptemp0, &
sigmai_chk(nbndsub, nksf)  !! extra
  logical :: already_skipped
  REAL(kind=DP), external :: efermig
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf
  REAL(kind=DP), allocatable :: xkf_all(:,:) , etf_all(:,:), &
                                sigmar_all (:,:), sigmai_all (:,:), zi_all (:,:), sigmai_chk_all (:,:)
  !
  WRITE(6,'(/5x,a)') repeat('=',67)
  WRITE(6,'(5x,"Electron (Imaginary) Self-Energy in the Migdal Approximation")') 
  WRITE(6,'(5x,a/)') repeat('=',67)
  !
  ! loop over temperatures can be introduced
  !
  eptemp0 = eptemp(1)
  !
  IF ( fsthick .lt. 1.d3 ) &
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
      'Fermi Surface thickness = ', fsthick, ' Ry'
  WRITE(stdout, '(/5x,a,e18.9,a)' ) &
    'Golden Rule strictly enforced with T = ',eptemp0, ' Ry'
  !
  ! here we take into account that we may skip bands when we wannierize
  ! (spin-unpolarized)
  already_skipped = .false.
  IF ( nbndskip .gt. 0 ) THEN
     IF ( .not. already_skipped ) THEN 
        nelec = nelec - two * nbndskip
        already_skipped = .true.
        WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
        WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec
     ENDIF
  ENDIF
  !
  ! Fermi level and corresponding DOS
  !
  ef0 = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw,ngaussw,0,isk)
  !
  !   if 'fine' Fermi level differs by more than 250 meV, there is probably
  !   something wrong with the wannier functions
  IF ( abs(ef0 - ef) * ryd2eV .gt. 0.5 .and. (.not.eig_read) ) &
       CALL errore ('selfen_elec', 'Something wrong with Ef, check MLWFs', 1)
  !
  dosef = dos_ef (ngaussw, degaussw, ef0, etf, wkf, nksf, nbndsub)
  !   N(Ef) in the equation for lambda is the DOS per spin
  dosef = dosef / two
  !
  WRITE (6, 100) degaussw, ngaussw
  WRITE (6, 101) dosef, ef0 * ryd2ev
  WRITE (6, 101) dosef, ef  * ryd2ev
  WRITE (6,'(a)') ' '
  !
  sigmar = zero
  sigmai = zero
  sigmai_chk = zero
  zi = zero
  !
  ! loop over all k points of the fine mesh
  !
  fermicount = 0 
  
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Tetrahedron method for sigmai
ncorn = 4 !! four corners of a tetrahedron
if(nxqf .ne. nqf1*nqf2*nqf3) then
	write(*,*) "Error in q-mesh: nxqf!= nqf1*nqf2*nqf3"
	call exit(1)
endif
ntet=nqf1*nqf2*nqf3*6
allocate(tet(ntet), stat=ierr)
allocate(wkt(nxqf), eprim(2*nmodes, nxqf, ibndmax-ibndmin+1), Fqnu(nxqf, 2*nmodes, ibndmax-ibndmin+1, ibndmax-ibndmin+1))
allocate(ekj(2*nmodes,nxqf ))
if (ierr/=0) print*, "tet : Allocation failed"


CALL make_kp_reg(nqf1,nqf2,nqf3,wkt,ntet, tet) !! associate q points with corners of tetrahedron
! DO iq = 1, nxqf
! 	write(*,*) xqf(:, iq)
! ENDDO  !! from the print outs, the arrangement of q points does follow the same way for tetra generation


DO ik = 1, nksqf
	!
	IF (lgamma) THEN
		ikk = ik
		ikq = ik
	ELSE
		ikk = 2 * ik - 1
		ikq = ikk + 1
	ENDIF
	!
	DO iq = 1, nxqf
		!
		! we read the hamiltonian eigenvalues (those at k+q depend on q!)
		!
		IF (etf_mem) THEN
			etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk, iq)
			etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq, iq)
		ELSE
			nrec = (iq-1) * nksf + ikk
			CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
			nrec = (iq-1) * nksf + ikq
			CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
		ENDIF
		!
		! here we must have ef, not ef0, to be consistent with ephwann_shuffle
		!
 !!		IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
 !!		( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
 		IF (1) THEN
			!
			DO imode= 1, nmodes
				!
				! the phonon frequency and Bose occupation
				wq = wf (imode, iq)
				wgq = wgauss( -wq/eptemp0, -99)
				wgq = wgq / ( one - two * wgq )
				!
				!  we read the e-p matrix
				!
				IF (etf_mem) THEN
					epf(:,:) = epf17 ( ik, iq, :, :, imode)
				ELSE
					nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
					CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
				ENDIF
				!
				DO  jbnd = 1, ibndmax-ibndmin+1
					!  the fermi occupation for k+q
					ekq = etf (ibndmin-1+jbnd, ikq) - ef0
					wgkq = wgauss( -ekq/eptemp0, -99)
					
					eprim( imode, iq , jbnd ) = ekq - wq
					eprim( nmodes + imode, iq, jbnd ) = ekq + wq
					DO ibnd = 1, ibndmax-ibndmin+1
						IF (wq .gt. eps_acustic) THEN
							g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
						ELSE
							g2 = 0.d0
						ENDIF
						Fqnu( iq, imode, ibnd, jbnd ) = g2 * ( wgq + wgkq )
						Fqnu( iq, nmodes+imode, ibnd, jbnd ) = g2 * ( 1 + wgq - wgkq )
					ENDDO ! ibnd
				ENDDO ! jbnd
			ENDDO ! imode
		ENDIF
	ENDDO ! iq
	DO jbnd = 1, ibndmax-ibndmin+1
		ekj(:,:)=eprim(:,:,jbnd)
		CALL eigen_tet(ntet,ekj,tet,2*nmodes,nxqf)
		DO ibnd = 1, ibndmax-ibndmin+1
			!  the energy of the electron at k (relative to Ef)
			ekk = etf (ibndmin-1+ibnd, ikk) - ef0
			!
			CALL weight_tet(nxqf, ntet,2*nmodes,ekk,tet,wkt)
			DO itet= 1, ntet
				DO jtet = 1, ncorn
					DO imode = 1, 2*nmodes
						! Gamma(ikk, ibnd, imode) = Gamma(ikk, ibnd, imode) + Fqnu( iq, imode, ibnd, jbnd ) * tet_weight(imode) * 2 * pi
						sigmai_chk(ibndmin-1+ibnd, ikk ) = sigmai_chk(ibndmin-1+ibnd, ikk ) + Fqnu( tet(itet)%p(jtet)%i, imode, ibnd, jbnd ) *tet(itet)%p(jtet)%c(imode) * pi
					  WRITE(6,'(/5x," wt: ", f14.10)') tet(itet)%p(jtet)%c(imode)
          ENDDO ! imode
				ENDDO ! jtet
			ENDDO ! itet
		ENDDO ! ibnd
	ENDDO ! jbnd
	write(stdout,*) "The ",ik,"th K point done..."
ENDDO ! ik

write(stdout,*) "At the end of Tetra method"
deallocate(tet)
deallocate(wkt, eprim, Fqnu)
deallocate(ekj)
  

  
  
  
  
  
  
  DO ik = 1, nksqf
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
     ! loop over the q-points
     !
     DO iq = 1, nxqf
        !
        ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
        !
        IF (etf_mem) THEN
           etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk, iq)
           etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq, iq)
        ELSE
           nrec = (iq-1) * nksf + ikk
           CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
           nrec = (iq-1) * nksf + ikq
           CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        ENDIF
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        !
        IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
             ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
           !
           fermicount = fermicount + 1
           DO imode = 1, nmodes
              !
              ! the phonon frequency and Bose occupation
              wq = wf (imode, iq)
              wgq = wgauss( -wq/eptemp0, -99)
              wgq = wgq / ( one - two * wgq )
              !
              !  we read the e-p matrix
              !
              IF (etf_mem) THEN
                 epf(:,:) = epf17 ( ik, iq, :, :, imode)
              ELSE
                 nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
                 CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
              ENDIF
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                 !
                 !  the energy of the electron at k (relative to Ef)
                 ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                 !
                 DO jbnd = 1, ibndmax-ibndmin+1
                    !
                    !  the fermi occupation for k+q
                    ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                    wgkq = wgauss( -ekq/eptemp0, -99)  
                    !
                    ! here we take into account the zero-point sqrt(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    IF (wq .gt. eps_acustic) THEN
                       g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                    ELSE
                       g2 = 0.d0
                    ENDIF
                    !
                    ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
                    ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book 
                    ! (Many-Particle Physics, 3rd edition)
                    ! 
                    weight = wqf(iq) * real (                                        &
                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) )
!                   ecutse needs to be defined if it's used 
!@                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                    !
                    sigmar(ibndmin-1+ibnd,ikk) = sigmar(ibndmin-1+ibnd,ikk) + g2 * weight
                    !
                    weight = wqf(iq) * aimag (                                        &
                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) ) 
!@                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                    !
                    sigmai(ibndmin-1+ibnd,ikk) = sigmai(ibndmin-1+ibnd,ikk) + g2 * weight
                    !
                    ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                    !
                    weight = wqf(iq) * &
                         ( (       wgkq + wgq ) * ( (ekk - ( ekq - wq ))**two - degaussw**two ) /       &
                                                  ( (ekk - ( ekq - wq ))**two + degaussw**two )**two +  &
                           ( one - wgkq + wgq ) * ( (ekk - ( ekq + wq ))**two - degaussw**two ) /       &
                                                  ( (ekk - ( ekq + wq ))**two + degaussw**two )**two )  
!@                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                    !
                    zi(ibndmin-1+ibnd,ikk) = zi(ibndmin-1+ibnd,ikk) + g2 * weight
                    !
                 ENDDO !jbnd
                 !
              ENDDO !ibnd
              !
           ENDDO !imode
           !
        ENDIF ! endif  fsthick
        !
     ENDDO  ! iq's
     !
  ENDDO ! end loop on k
  !
  ! The k points are distributed among pools: here we collect them
  !
  nksqtotf = nkstotf/2 ! odd-even for k,k+q
  !
  ALLOCATE ( xkf_all    ( 3,       nkstotf ), &
             etf_all    ( nbndsub, nkstotf ), &
             sigmar_all ( nbndsub, nkstotf ), &  
             sigmai_all ( nbndsub, nkstotf ), &  
             sigmai_chk_all ( nbndsub, nkstotf ), &
             zi_all     ( nbndsub, nkstotf )  )
  !
#ifdef __PARA
  !
  ! note that poolgather2 works with the doubled grid (k and k+q)
  ! therefore we need to use the sigma array with both grids, even
  ! though one of them is useless. This should be fixed by modifying
  ! poolgather2 (it's a waste of memory).
  !
  CALL poolgather2 ( 3,       nkstotf, nksf, xkf,    xkf_all  )
  CALL poolgather2 ( nbndsub, nkstotf, nksf, etf,    etf_all  )
  CALL poolgather2 ( nbndsub, nkstotf, nksf, sigmar, sigmar_all)
  CALL poolgather2 ( nbndsub, nkstotf, nksf, sigmai, sigmai_all)
  CALL poolgather2 ( nbndsub, nkstotf, nksf, sigmai_chk, sigmai_chk_all)
  CALL poolgather2 ( nbndsub, nkstotf, nksf, zi,     zi_all)
  CALL mp_sum(fermicount, inter_pool_comm)
  !
  ! test output from each pool
  ! DO ik = 1, nksqf
  !    IF (lgamma) THEN
  !       ikk = ik
  !    ELSE
  !       ikk = 2 * ik - 1
  !    ENDIF
  !    WRITE(1000+my_pool_id,'(/5x,"ik = ",i5," coord.: ", 3f9.5)') ik, xkf(:,ikk)
  !    WRITE(1000+my_pool_id,'(6(2x,f12.6))') ( ryd2mev*sigmar(ibnd,ikk), ibnd=ibndmin,ibndmax )
  ! ENDDO
  ! CALL mp_barrier()
  !
#else
  !
  xkf_all = xkf
  etf_all = etf
  sigmar_all = sigmar
  sigmai_all = sigmai
  sigmai_chk_all = sigmai_chk
  zi_all     = zi
  !
#endif
  !
  WRITE(6,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")') 
  !
  WRITE( stdout, '(/5x,a,i5,a,i5/)' ) &
    'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nksqtotf*nxqf
  !
  DO ik = 1, nksqtotf
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
     WRITE(stdout,'(/5x,"ik = ",i5," coord.: ", 3f9.5)') ik, xkf_all(:,ikk)
     WRITE(stdout,'(5x,a)') repeat('-',67)
     !
     DO ibnd = ibndmin, ibndmax
        !
        ! note that ekk does not depend on q 
        ekk = etf_all (ibnd, ikk) - ef0
        !
        ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
        zi_all (ibnd,ikk) = one / ( one + zi_all (ibnd,ikk) )
        !
!        WRITE(stdout, 102) ibnd, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk), &
!            ryd2mev * sigmai_all (ibnd,ikk),  zi_all (ibnd,ikk)
        WRITE(stdout, 103) ik, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk), &
              ryd2mev * sigmai_all (ibnd,ikk), ryd2mev * sigmai_chk_all (ibnd,ikk) ! ,  zi_all (ibnd,ikk)
       !
    ENDDO
    WRITE(stdout,'(5x,a/)') repeat('-',67)
    !
  ENDDO
  !
100 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'E( ',i3,' )=',f9.3,' eV   Re[Sigma]=',f9.3,' meV   Im[Sigma]=',f9.3,' meV     Z=',f9.3)
103 format(5x,'k( ',i6,' )=',f10.4,' eV   Re[Sigma]=',f10.4,' meV   Im[Sigma]=',f10.4,' meV     Im[Sigma]1=',f10.4)
  !
  RETURN
  !
  END SUBROUTINE selfen_elec_bo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine make_kp_reg(nx,ny,nz,wkt,ntet, tet)
! generate a regular mesh from 0 to g_i, with eventual shift of origin only
! Tetrahedron are generated. How are Tetra near BZ boundaries treated correctly?
use tetrahedron
USE kinds, ONLY : DP
implicit none
integer :: i,j,k,l,nx,ny,nz,nk,nt,np,pt(8), ntet
real(DP) :: q(3),ep(3),kpt(3,nx*ny*nz),wkt(nx*ny*nz)
type(tetra) :: tet(ntet)

ep = 0d0  ! shift by a small amount so that only one boundary K remains in the FBZ after folding
! shft = (-0.5d0)*(g1+g2+g3)

do i=1,nx*ny*nz
	wkt(i)=1d0/(nx*ny*nz)
enddo

np = 0
nt = 0
do i = 1,nx
	do j = 1,ny
		do k = 1,nz			
			do l = 1,6
				pt(1) = (i-1)*ny*nz + (j-1)*nz + (k-1) +1
				pt(2) = pt(1)+1					
				pt(3) = pt(1)+nz
				pt(4) = pt(1)+nz+1
				pt(5) = pt(1)+nz*ny
				pt(6) = pt(1)+nz*ny+1
				pt(7) = pt(1)+nz*ny+nz
				pt(8) = pt(1)+nz*ny+nz+1
				if(k .eq. nz) then
					pt(2)=pt(2)-nz
					pt(4)=pt(4)-nz
					pt(6)=pt(6)-nz
					pt(8)=pt(8)-nz
				endif
				if(j .eq. ny) then
					pt(3)=pt(3)-ny*nz
					pt(4)=pt(4)-ny*nz
					pt(7)=pt(7)-ny*nz
					pt(8)=pt(8)-ny*nz
				endif
				if(i .eq. nx) then
					pt(5)=pt(5)-nx*ny*nz
					pt(6)=pt(6)-nx*ny*nz
					pt(7)=pt(7)-nx*ny*nz
					pt(8)=pt(8)-nx*ny*nz
				endif

				nt = nt+1
				np = np+1

				if (l.eq.1) then   			! 1,2,3,6
					tet(nt)%p(1)%i = pt(1)
					tet(nt)%p(2)%i = pt(2)
					tet(nt)%p(3)%i = pt(3)
					tet(nt)%p(4)%i = pt(6)
				elseif (l.eq.2) then			! 2,3,4,6
					tet(nt)%p(1)%i = pt(2)
					tet(nt)%p(2)%i = pt(3)
					tet(nt)%p(3)%i = pt(4)
					tet(nt)%p(4)%i = pt(6)
				elseif (l.eq.3) then			! 1,3,5,6
					tet(nt)%p(1)%i = pt(1)
					tet(nt)%p(2)%i = pt(3)
					tet(nt)%p(3)%i = pt(5)
					tet(nt)%p(4)%i = pt(6)
				elseif (l.eq.4) then			! 3,4,6,8
					tet(nt)%p(1)%i = pt(3)
					tet(nt)%p(2)%i = pt(4)
					tet(nt)%p(3)%i = pt(6)
					tet(nt)%p(4)%i = pt(8)
				elseif (l.eq.5) then			! 3,5,6,7
					tet(nt)%p(1)%i = pt(3)
					tet(nt)%p(2)%i = pt(5)
					tet(nt)%p(3)%i = pt(6)
					tet(nt)%p(4)%i = pt(7)
				else           			! 3,6,7,8
					tet(nt)%p(1)%i = pt(3)
					tet(nt)%p(2)%i = pt(6)
					tet(nt)%p(3)%i = pt(7)
					tet(nt)%p(4)%i = pt(8)
				endif
			enddo	
		enddo
	enddo
enddo

! write(ulog,*)'KP_REG: Number of regular kpoints generated is=',nk
! 
! if (mod(nx,2).eq.0 .and. mod(ny,2).eq.0 .and. mod(nz,2).eq.0 ) then
! 	do i=1,nk
! 		kpt(:,i) = kpt(:,i)+shft(:)
! 	enddo
! endif
! 
! 2  format(i7,2x,3(1x,f12.5),5x,f9.5)
! 3  format(3(i3),2x,i6,2x,3(1x,f12.5),5x,f9.5)
! close(126)

end subroutine make_kp_reg


subroutine eigen_tet(ntet,eival,tet,n,nkp)
use tetrahedron
USE kinds, ONLY : DP
implicit none
integer j,k,l,ntet, n,nkp,w,cnd
real(DP), intent(in) :: eival(n,nkp)
type(tetra), intent(inout) :: tet(ntet)

do j=1,ntet !! this many tetra
	do k=1,4   !! corners
		do l=1,n   !! number of bands
			w = tet(j)%p(k)%i  !! i is the index of k points in the original kmesh
			tet(j)%p(k)%w(l)=eival(l,w) !! eigenvalues (band energies of band l, kpoint "w"
		enddo
	enddo
enddo

end subroutine eigen_tet

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine weight_tet(nktot, ntet,nb,e,tet,wkt)
!! input:  number of bands, the energy e that the weights are going to be computed, tet: tetrahedron
!! the first 5 parameters are easy. Just need to get tet right.
use tetrahedron
USE kinds, ONLY : DP
implicit none
integer i,j,k,l,m,nb,nktot, ntet,uns(4),zero,kkk,ii
real(DP) :: e1,e2,e3,e4,e21,e31,e41,e32,e42,e43,en(4),en1(4),a1,a2,a3,a4,b11,b12,b21,b22,b31,b32,b41,b42
real(DP) :: a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44
real(DP) :: d11,d12,d13,d21,d22,d23,d31,d32,d33,d41,d42,d43
real(DP) :: c11,c12,c21,c22,c31,c32,c41,c42
real(DP) :: s1,s2,s3,s4,s5,s6
real(DP) :: e,wkt(nktot)
type(tetra), intent(inout) :: tet(ntet)

zero=0

do i=1,ntet    ! allocate weights for imaginary terms
	do j=1,nb
		do k=1,4
			en(k)=tet(i)%p(k)%w(j)
		enddo
		en1=en
		call ssort(en,uns,4)  ! 
		if (en(1).ne.en(2) .and. en(2).ne.en(3) .and. en(3).ne.en(4)) then
			do k=1,4
				if (en(k).eq.en1(1)) then
					uns(k) = 1
				elseif (en(k).eq.en1(2)) then
					uns(k) = 2
				elseif (en(k).eq.en1(3)) then
					uns(k) = 3
				elseif (en(k).eq.en1(4)) then
					uns(k) = 4
				endif
			enddo
		else
			do k=1,4
				if (en(k).eq.en1(1)) then
					uns(k) = 1
				elseif (en(k).eq.en1(2)) then
					uns(k) = 2
				elseif (en(k).eq.en1(3)) then
					uns(k) = 3
				elseif (en(k).eq.en1(4)) then
					uns(k) = 4
				endif
				if (k.ne.1) then
					do l=1,k-1
						if (uns(l).eq.uns(k) .and. uns(k).eq.1) then
							if (en(k).eq.en1(2)) then
								uns(k) = 2
							elseif (en(k).eq.en1(3)) then
								uns(k) = 3
							elseif (en(k).eq.en1(4)) then
								uns(k) = 4
							endif
						elseif (uns(l).eq.uns(k) .and. uns(k).eq.2) then
							if (en(k).eq.en1(3)) then
								uns(k) = 3
							elseif (en(k).eq.en1(4)) then
								uns(k) = 4
							endif
						elseif (uns(l).eq.uns(k) .and. uns(k).eq.3) then
							uns(k) = 4
						endif
					enddo
				endif
			enddo  !! k = 1,4
		endif
		e1 = en(4)
		e2 = en(3)
		e3 = en(2)
		e4 = en(1)
		e21 = e2-e1
		e31 = e3-e1
		e41 = e4-e1
		e32 = e3-e2
		e42 = e4-e2
		e43 = e4-e3
! 		if(e21 .eq. 0. .or. e31 .eq. 0. .or. e41 .eq. 0. .or. e32 .eq. 0. .or. e42 .eq. 0. .or. e43 .eq. 0.) then 
! 			write(*,*) "Equal e tetra corners: ", i
! 		endif
		if (e.lt.e1 .or. e.gt.e4) then
			tet(i)%p(uns(4))%c(j) = 0
			tet(i)%p(uns(3))%c(j) = 0
			tet(i)%p(uns(2))%c(j) = 0
			tet(i)%p(uns(1))%c(j) = 0
		elseif (e1.le.e .and. e.le.e2) then
			tet(i)%p(uns(4))%c(j) = ((e2-e)/e21+(e3-e)/e31+(e4-e)/e41)*(e-e1)**2/(e41*e31*e21)
			tet(i)%p(uns(3))%c(j) = (e-e1)**3/(e21**2*e31*e41)
			tet(i)%p(uns(2))%c(j) = (e-e1)**3/(e21*e31**2*e41)
			tet(i)%p(uns(1))%c(j) = (e-e1)**3/(e21*e31*e41**2)
		elseif (e2.le.e .and. e.le.e3) then
			c11 = (e3-e)/e31**2
			c12 = (e4-e)/e41**2
			c21 = (e3-e)/e32**2
			c22 = (e4-e)/e42**2
			c31 = (e-e2)/e32**2
			c32 = (e-e1)/e31**2
			c41 = (e-e2)/e42**2
			c42 = (e-e1)/e41**2
			b11 = (e3-e)*(e-e2)/(e42*e32)+(e4-e)*(e-e1)/(e41*e42)+(e3-e)*(e-e1)/(e32*e41)
			b12 = (e4-e)*(e-e1)/(e42*e31)+(e4-e)*(e-e2)/(e42*e32)+(e3-e)*(e-e1)/(e31*e32)
			b21 = (e3-e)*(e-e2)/(e42*e31)+(e4-e)*(e-e2)/(e42*e41)+(e3-e)*(e-e1)/(e31*e41)
			b22 = (e3-e)*(e-e2)/(e32*e31)+(e4-e)*(e-e1)/(e41*e31)+(e4-e)*(e-e2)/(e32*e41)
			b31 = (e3-e)*(e-e2)/(e42*e31)+(e4-e)*(e-e2)/(e42*e41)+(e3-e)*(e-e1)/(e31*e41)
			b32 = (e3-e)*(e-e2)/(e42*e32)+(e4-e)*(e-e1)/(e41*e42)+(e3-e)*(e-e1)/(e32*e41)
			b41 = (e3-e)*(e-e2)/(e32*e31)+(e4-e)*(e-e1)/(e41*e31)+(e4-e)*(e-e2)/(e32*e41)
			b42 = (e4-e)*(e-e1)/(e42*e31)+(e4-e)*(e-e2)/(e42*e32)+(e3-e)*(e-e1)/(e31*e32)
			tet(i)%p(uns(4))%c(j) = .5*(c11*b11+c12*b12)
			tet(i)%p(uns(3))%c(j) = .5*(c21*b21+c22*b22)
			tet(i)%p(uns(2))%c(j) = .5*(c31*b31+c32*b32)
			tet(i)%p(uns(1))%c(j) = .5*(c41*b41+c42*b42)
		elseif (e3.le.e .and. e.le.e4) then
			tet(i)%p(uns(4))%c(j) = (e4-e)**3/(e41**2*e42*e43)
			tet(i)%p(uns(3))%c(j) = (e4-e)**3/(e41*e42**2*e43)
			tet(i)%p(uns(2))%c(j) = (e4-e)**3/(e41*e42*e43**2)
			tet(i)%p(uns(1))%c(j) = ((e-e3)/e43+(e-e2)/e42+(e-e1)/e41)*(e4-e)**2/(e41*e42*e43)
		endif
		tet(i)%p(uns(1))%c(j)=tet(i)%p(uns(1))%c(j)*wkt(tet(i)%p(uns(1))%i)*1.0d0/6.0d0 !! the weight is V_t/V_G which is the inverse of # tetra
		tet(i)%p(uns(2))%c(j)=tet(i)%p(uns(2))%c(j)*wkt(tet(i)%p(uns(2))%i)*1.0d0/6.0d0 !! spin=2 should not be included here, it's not relevant
		tet(i)%p(uns(3))%c(j)=tet(i)%p(uns(3))%c(j)*wkt(tet(i)%p(uns(3))%i)*1.0d0/6.0d0
		tet(i)%p(uns(4))%c(j)=tet(i)%p(uns(4))%c(j)*wkt(tet(i)%p(uns(4))%i)*1.0d0/6.0d0
		do ii=1,4
			if(isnan(tet(i)%p(uns(ii))%c(j) ) ) then
		 		tet(i)%p(uns(ii))%c(j)=0.0d0
		 	endif
		enddo
	enddo
enddo

end subroutine weight_tet

SUBROUTINE SSORT (X, IY, N)
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER N
REAL(DP) :: X(1:N)
INTEGER IY(N)
REAL(DP) :: TEMP
INTEGER I, ISWAP(1), ITEMP, ISWAP1
INTRINSIC MAXLOC
DO 200 I=1,N-1
	ISWAP=MAXLOC(X(I:N))
	ISWAP1=ISWAP(1)+I-1
	IF(ISWAP1.NE.I) THEN
		TEMP=X(I)
		X(I)=X(ISWAP1)
		X(ISWAP1)=TEMP
		ITEMP=IY(I)
		IY(I)=IY(ISWAP1)
		IY(ISWAP1)=ITEMP
	ENDIF
200 CONTINUE
RETURN
END