  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !   
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_elec
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
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  USE phcom, ONLY : lgamma, nmodes
  USE epwcom, ONLY : nbndsub, lrepmatf, iunepmatf, &
      fsthick, eptemp, ngaussw, degaussw, iuetf,   &
      nbndskip, ecutse, parallel_k, &
      parallel_q, epf_mem, etf_mem, eig_read, eps_acustic
  USE pwcom, ONLY : nelec, ef, isk
  USE el_phon, ONLY : etf, ibndmin, ibndmax, nksf, etfq, &
      epf17, wkf, nksqf, nxqf, wf, wqf, xkf, nkstotf, &
      sigmar_all, sigmai_all, zi_all
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : me_pool, inter_pool_comm, my_pool_id
#endif
  implicit none
  !
  REAL(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq,  &
       weight, wgauss, dosef, dos_ef, sigmar(nbndsub, nksf), &
       sigmai(nbndsub, nksf), zi(nbndsub, nksf), eptemp0
  logical :: already_skipped
  REAL(kind=DP), external :: efermig
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf
  REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:), etf_all(:,:)
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
  zi = zero
  !
  ! loop over all k points of the fine mesh
  !
  fermicount = 0 
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
             zi_all     ( nbndsub, nkstotf )  )
  xkf_all(:,:) = zero
  etf_all(:,:) = zero
  sigmar_all(:,:) = zero
  sigmai_all(:,:) = zero
  zi_all(:,:) = zero
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
  CALL poolgather2 ( nbndsub, nkstotf, nksf, zi,     zi_all)
  CALL mp_sum(fermicount, inter_pool_comm)
  CALL mp_barrier()
  !

#else
  !
  xkf_all = xkf
  etf_all = etf
  sigmar_all = sigmar
  sigmai_all = sigmai
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
              ryd2mev * sigmai_all (ibnd,ikk),  zi_all (ibnd,ikk)
       !
    ENDDO
    WRITE(stdout,'(5x,a/)') repeat('-',67)
    !
  ENDDO
  !
  IF ( ALLOCATED(xkf_all) )    DEALLOCATE( xkf_all )
  IF ( ALLOCATED(etf_all) )    DEALLOCATE( etf_all )
  IF ( ALLOCATED(sigmar_all) ) DEALLOCATE( sigmar_all )
  IF ( ALLOCATED(sigmai_all) ) DEALLOCATE( sigmai_all )
  IF ( ALLOCATED(zi_all) )     DEALLOCATE( zi_all )
  !
100 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'E( ',i3,' )=',f9.3,' eV   Re[Sigma]=',f9.3,' meV   Im[Sigma]=',f9.3,' meV     Z=',f9.3)
103 format(5x,'k( ',i6,' )=',f10.4,' eV   Re[Sigma]=',f10.4,' meV   Im[Sigma]=',f10.4,' meV     Z=',f9.3)
  !
  RETURN
  !
  END SUBROUTINE selfen_elec
  !   
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_elec_fly ( iq )
  !-----------------------------------------------------------------------
  !
  !  On-the-fly calculation of electron self energy
  !
  !  compute the imaginary part of the electron self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the electron linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !  Modified by Bo Qiu based on the original subroutine selfen_elec
  !  Email: 200210qb@gmail.com
  !  March 11th, 2013
  !
  !  02/06/2013 Modified by Roxana Margine 
  !
  !  This subroutine computes the contribution from phonon iq to all k-points
  !  The outer loop in ephwann_shuffle.f90 will loop over all iq points
  !  The contribution from each iq is summed at the end of this subroutine for iq=nxqf 
  !  to recover the per-ik electron self energy
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  USE phcom, ONLY : lgamma, nmodes
  USE epwcom, ONLY : nbndsub, lrepmatf, iunepmatf, &
      fsthick, eptemp, ngaussw, degaussw, iuetf,   &
      nbndskip, ecutse, parallel_k, &
      parallel_q, epf_mem, etf_mem, eig_read, eps_acustic
  USE pwcom, ONLY : nelec, ef, isk
  USE el_phon, ONLY : etf, ibndmin, ibndmax, nksf, etfq, &
      epf17, wkf, nksqf, nxqf, wf, wqf, xkf, nkstotf, &
      sigmar_all, sigmai_all, zi_all
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : me_pool, inter_pool_comm, my_pool_id
#endif
  implicit none
  !
  REAL(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq,  &
       weight, wgauss, dosef, dos_ef, sigmar(nbndsub, nksf), &
       sigmai(nbndsub, nksf), zi(nbndsub, nksf), eptemp0
  logical :: already_skipped
  REAL(kind=DP), external :: efermig
  integer, parameter :: out_unit=20
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf
  REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:), etf_all(:,:), &
                                sigmar_all_q(:,:), sigmai_all_q(:,:), zi_all_q(:,:)          
  !
  !
  ! loop over temperatures can be introduced
  !
  !WRITE(6,'(5x,"made nrec_end3= ")')
  !open (unit=out_unit,file="test.txt",action="write",status="replace")
  eptemp0 = eptemp(1)

  !
  IF (iq.eq.1) THEN
     !
     !WRITE(stdout,'(/5x,a)') repeat('=',67)
     !WRITE(6,'(5x,"Electron (Imaginary) Self-Energy in the Migdal Approximation (on the fly)")')
     !WRITE(6,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick .lt. 1.d3 ) &
        !WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick, ' Ry'
     !WRITE(stdout, '(/5x,a,e18.9,a)' ) &
     !     'Golden Rule strictly enforced with T = ',eptemp0, ' Ry'
     !
     ! here we take into account that we may skip bands when we wannierize
     ! (spin-unpolarized)
     !
     already_skipped = .false.
     IF ( nbndskip .gt. 0 ) THEN
        IF ( .not. already_skipped ) THEN
           nelec = nelec - two * nbndskip
           already_skipped = .true.
           !WRITE(stdout,'(/5x,"Skipping the first ",i4," bands:")') nbndskip
           !WRITE(stdout,'(/5x,"The Fermi level will be determined with ",f9.5," electrons")') nelec
        ENDIF
     ENDIF
     !
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
  IF ( iq .eq. 1 ) THEN 
     !WRITE (6, 100) degaussw, ngaussw
     !WRITE (6, 101) dosef, ef0 * ryd2ev
     !WRITE (6, 101) dosef, ef  * ryd2ev
     !WRITE (6,'(a)') ' '
100 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
  ENDIF
  sigmar = zero
  sigmai = zero
  zi = zero
  !
  ! loop over all k points of the fine mesh
  !
  fermicount = 0 
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
     ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
     ! when we see references to iq, it is always = 1 for on the fly calculations
     !
     IF (etf_mem) THEN
        etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk, 1)
        etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq, 1)
     ELSE
        nrec = (1-1) * nksf + ikk
        !WRITE(6,'(5x,"made nrec= ",i5, " ikk =", i5," ik =", i5)') nrec, ikk,ik
        CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        nrec = (1-1) * nksf + ikq
        !WRITE(6,'(5x,"made ikq= ",i5)') nrec
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
              epf(:,:) = epf17 ( ik, 1, :, :, imode)
           ELSE
              nrec = (1-1) * nmodes * nksqf + (imode-1) * nksqf + ik
              !WRITE(6,'(5x,"made it= ",i5)') nrec
              CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
              !CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, 76, nrec, -1)
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
!@                  if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                 !
                 sigmar(ibndmin-1+ibnd,ikk) = sigmar(ibndmin-1+ibnd,ikk) + g2 * weight
                 !
                 weight = wqf(iq) * aimag (                                        &
                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) ) 
!@                 if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
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
!@                 if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
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
  ENDDO ! end loop on k
  !
  ! The k points are distributed among pools: here we collect them
  !
  nksqtotf = nkstotf/2 ! odd-even for k,k+q
  !
  ALLOCATE ( xkf_all      ( 3,       nkstotf ), &
             etf_all      ( nbndsub, nkstotf ), &
             sigmar_all_q ( nbndsub, nkstotf ), &
             sigmai_all_q ( nbndsub, nkstotf ), &
             zi_all_q     ( nbndsub, nkstotf )  )
  xkf_all(:,:) = zero
  etf_all(:,:) = zero
  sigmar_all_q(:,:) = zero
  sigmai_all_q(:,:) = zero
  zi_all_q(:,:) = zero
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
  CALL poolgather2 ( nbndsub, nkstotf, nksf, sigmar, sigmar_all_q)
  CALL poolgather2 ( nbndsub, nkstotf, nksf, sigmai, sigmai_all_q)
  CALL poolgather2 ( nbndsub, nkstotf, nksf, zi,     zi_all_q)
  CALL mp_sum(fermicount, inter_pool_comm)
  CALL mp_barrier()
  !
#else
  !
  xkf_all = xkf
  etf_all = etf
  sigmar_all_q = sigmar
  sigmai_all_q = sigmai
  zi_all_q     = zi
  !
#endif
  !
  IF ( iq .eq. 1 ) THEN 
     IF (.not. ALLOCATED (sigmar_all)) then
        ALLOCATE(sigmar_all(nbndsub, nkstotf)) !! for the sum of selfen_elec contribution from each iq
        sigmar_all = zero
     ENDIF
     IF (.not. ALLOCATED (sigmai_all)) then
        ALLOCATE(sigmai_all(nbndsub, nkstotf)) !! for the sum of selfen_elec contribution from each iq
        sigmai_all = zero
     ENDIF
     IF (.not. ALLOCATED (zi_all)) then
        ALLOCATE(zi_all(nbndsub, nkstotf)) !! for the sum of selfen_elec contribution from each iq
        zi_all = zero
     ENDIF
  ENDIF
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
     DO ibnd = ibndmin, ibndmax
        !
        ! note that ekk does not depend on q 
        ekk = etf_all (ibnd, ikk) - ef0
        !
        ! sum contribution from each iq to global sigmar, sigmai, zi
        sigmar_all (ibnd,ikk) = sigmar_all (ibnd,ikk) + sigmar_all_q (ibnd,ikk) 
        sigmai_all (ibnd,ikk) = sigmai_all (ibnd,ikk) + sigmai_all_q (ibnd,ikk) 
        zi_all (ibnd,ikk) = zi_all (ibnd,ikk) + zi_all_q (ibnd,ikk) 
       !
    ENDDO
!     WRITE(1000+my_pool_id,'(/5x,"iq = ",i5," ik = ",i5)') iq, ik
!     WRITE(1000+my_pool_id,'(6(2x,f12.6))') (sigmar_all(ibnd,ikk),ibnd=ibndmin,ibndmax )
    !
  ENDDO
  !
  ! Output electron SE here after looping over all q-points (with their contributions 
  ! summed in sigmar_all, etc.)
  IF ( iq .eq. nxqf ) THEN
     !
     WRITE(6,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
     !
     open (unit=out_unit,file="test.txt",action="write",status="replace")
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
!           WRITE(stdout, 102) ibnd, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk), &
!                               ryd2mev * sigmai_all (ibnd,ikk),  zi_all (ibnd,ikk)
           WRITE(stdout, 103) ik, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk), &
                              ryd2mev * sigmai_all (ibnd,ikk), zi_all (ibnd,ikk)
           write (out_unit,*)  ryd2ev * ekk, ryd2mev * sigmai_all (ibnd,ikk)
        ENDDO
        WRITE(stdout,'(5x,a/)') repeat('-',67)
        !
     ENDDO
     !
     IF ( ALLOCATED(sigmar_all) ) DEALLOCATE( sigmar_all )
     IF ( ALLOCATED(sigmai_all) ) DEALLOCATE( sigmai_all )
     IF ( ALLOCATED(zi_all) )     DEALLOCATE( zi_all )
     !
     close (out_unit)
     102 FORMAT(5x,'E( ',i3,' )=',f9.3,' eV   Re[Sigma]=',f9.3,' meV   Im[Sigma]=',f9.3,' meV     Z=',f9.3)
     103 format(5x,'k( ',i6,' )=',f10.4,' eV   Re[Sigma]=',f10.4,' meV   Im[Sigma]=',f10.4,' meV     Z=',f9.3)
     !
  ENDIF
  !
  IF ( ALLOCATED(xkf_all) )      DEALLOCATE( xkf_all )
  IF ( ALLOCATED(etf_all) )      DEALLOCATE( etf_all )
  IF ( ALLOCATED(sigmar_all_q) ) DEALLOCATE( sigmar_all_q )
  IF ( ALLOCATED(sigmai_all_q) ) DEALLOCATE( sigmai_all_q )
  IF ( ALLOCATED(zi_all_q) )     DEALLOCATE( zi_all_q )
  !
  RETURN
  !
  END SUBROUTINE selfen_elec_fly
