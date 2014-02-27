 !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine selfen_phon_tet
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  !USE phcom,     ONLY : nmodes
  USE epwcom,    ONLY : nbndsub, lrepmatf, iunepmatf, fsthick, &
                        eptemp, ngaussw, degaussw, iuetf,     &
                        wmin, wmax, nw, nbndskip, a2f, epf_mem, etf_mem, &
                        nsmear, delta_smear, eig_read, eps_acustic,  &
                        nqf1, nqf2, nqf3, nkf1, nkf2, nkf3 !! extra from schuberm 
  USE pwcom,     ONLY : nelec, ef, isk, nbnd, wg, et, nspin, s, t_rev, irt, ftau, nsym, invsym, d1,d2,d3, &
                                 time_reversal
 
  USE el_phon,   ONLY : epf17, ibndmax, ibndmin, etf, &
                        etfq, wkf, xqf, wqf, nksf, nxqf,   &
                        nksqf, wf, nkstotf, xkf, xqf, dmef, &
                        lambda_all, lambda_v_all, nrottet

  !Modules added for tetrahedra-related subroutines
  !USE wvfct,              ONLY : nbnd, wg, et
  !USE lsda_mod,           ONLY : nspin, isk
  !USE klist,              ONLY : xk, wk, nks
  use phcom
  !USE ktetra,             ONLY : k1, k2, k3
  !USE symme,              ONLY : s, t_rev, irt, ftau, nsym, invsym, d1,d2,d3, &
  !                               time_reversal
#ifdef __PARA
  USE mp,        ONLY : mp_barrier,mp_sum
  USE mp_global, ONLY : me_pool,inter_pool_comm,my_pool_id
#endif
  !
  implicit none
  !
  real(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0), czero=(0.d0,0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount, ismear
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgkk, wgkq, lambda, lambda_v, &
                   weight, wgauss, dosef, dos_ef, w0g1, w0g2, w0gauss,    &
                   gamma(nmodes),gamma_v(nmodes), degaussw0, eptemp0, lambda_tot
  !
  real(kind=DP), external :: efermig
  real(kind=DP), external :: efermit
  real(kind=DP) :: gamma_all  ( nmodes, nxqf, 10 )
  real(kind=DP) :: gamma_all_v  ( nmodes, nxqf, 10 )
  real(kind=DP) :: coskkq(nbndsub, nbndsub)
  real(kind=DP) :: DDOT,  vkk(3,nbndsub), vkq(3,nbndsub)
  !
  !Variables added for tetrahedra-related subroutines
  integer :: is
  real(kind=DP):: nks_
  real(kind=DP):: wgf (nbnd, nksf),wgff (nbnd, nksf), eff, gamma_chk(nmodes),Fqnu(nksqf, nmodes, ibndmax-ibndmin+1)
  integer :: ntetraf, itet, jtet,kp1, kp2, kp3, kp4, itetra (4)
  integer, allocatable :: tetraf(:,:)
  real(kind=DP):: xkf_(3,nkf1*nkf2*nkf3), wkf_(nkf1*nkf2*nkf3)
  !real(kind=DP), ALLOCATABLE :: xkf_(:,:), wkf_(:)
  integer, parameter :: out_unit=20
  !
  WRITE(6,'(/5x,a)') repeat('=',67)
  WRITE(6,'(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation")') 
  WRITE(6,'(5x,a/)') repeat('=',67)

  open (unit=out_unit,file="tet.txt",action="write",status="replace")
  !
  IF ( fsthick .lt. 1.d3 ) &
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
     'Fermi Surface thickness = ', fsthick, ' Ry'
  WRITE(stdout, '(/5x,a,e18.9,a)' ) &
     'Golden Rule strictly enforced with T = ',eptemp(1), ' Ry'
  ! here we take into account that we may skip bands when we wannierize
  ! (spin-unpolarized)
  ! 
  IF (nbndskip.gt.0) then
     nelec = nelec - two * nbndskip
     WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
     WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec     
  ENDIF
  !
  ! here we loop on smearing values JN - this may be an option later on 
  !
  !

  !DO ismear = 1, nsmear
     !
     degaussw0 = (ismear-1) * delta_smear + degaussw
     eptemp0 = (ismear-1) * delta_smear + eptemp(1)
     !
     IF (.not.ALLOCATED (lambda_all) .and. .not.ALLOCATED(lambda_v_all)) then
        ALLOCATE(lambda_all(nmodes, nxqf, nsmear))
        ALLOCATE(lambda_v_all(nmodes, nxqf, nsmear))
     ENDIF
     gamma_all =  0.d0
     gamma_all_v = 0.d0
     lambda_all = 0.d0
     lambda_v_all = 0.d0
     ! 
     !
     ! Fermi level and corresponding DOS
     !
     !   Note that the weights of k+q points must be set to zero here
     !   no spin-polarized calculation here
    ef0 = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw0,ngaussw,0,isk)
    !Number of tetrahedra for fine mesh
    ntetraf=6*nkf1*nkf2*nkf3
    ALLOCATE ( tetraf  ( 4,  ntetraf))
    !ALLOCATE ( xkf_  ( 3,  nkf1*nkf2*nkf3))
    !ALLOCATE ( wkf_  ( nkf1*nkf2*nkf3))
    
    !defined kpoint_grid ( nrot, time_reversal, s, t_rev, bg, npk, &
    !                       k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
    call kpoint_grid ( nrottet, time_reversal, s, t_rev, bg, nkf1*nkf2*nkf3, &
                          0,0,0, nkf1,nkf2,nkf3, nksf, xkf_, wkf_)

    !defined tetrahedra ( nsym, s, minus_q, at, bg, npk, k1,k2,k3, &
    !   nk1,nk2,nk3, nks, xk, wk, ntetra, tetra )
    call tetrahedra ( nrottet, s, minus_q, at, bg, nkf1*nkf2*nkf3, 0,0,0, &
       nkf1,nkf2,nkf3, nksf, xkf_, wkf_, ntetraf, tetraf )

     !Calculate the Fermi energy ef
     !defined efermit (et, nbnd, nks, nelec, nspin, ntetra, tetra, is, isk)
    eff = efermit (etf, nbndsub, nksf, nelec, nspin, ntetraf, tetraf, 0, isk)

    !defined tweights (nks, nspin, nbnd, nelec, ntetra, tetra, et, &
    !   ef, wg, is, isk )
    call tweights (nksf, nspin, nbnd, nelec, ntetraf, tetraf, etf, &
       eff, wgf, 0, isk )

    WRITE(6,'(/5x,"wgf = ",i5," nks = ",i5," nksf = ",i5," nksqf = ",i5, " ntetra =", i5)') size(wgf,2),nks_,nksf,nksqf, ntetraf
    WRITE(6,'(/5x,"eff = ",f14.10," ef0 = ",f14.10," ef = ",f14.10)') eff* ryd2eV, ef0* ryd2eV, ef* ryd2eV
    wgff=wgf
    !ALLOCATE (Fqnu(nksqf, nmodes, ibndmax-ibndmin+1, ibndmax-ibndmin+1))
    !DO ik = 1, nksqf
     ! ikk = 2 * ik - 1
      !DO ibnd = 1, nbndsub
       !     WRITE(6,'(/5x,"ikk = ",i5," ibnd: ", i5, " wt: ", f14.10)') ikk, ibnd, wgf(ibnd,ikk)
      !ENDDO
    !ENDDO
     !   if 'fine' Fermi level differs by more than 250 meV, there is probably
     !   something wrong with the wannier functions
     IF (abs(ef0 - ef) * ryd2eV .gt. 0.5 .and. (.not.eig_read) ) &
        CALL errore ('selfen_phon', 'Something wrong with Ef, check MLWFs', 1)
     !
     dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
     !   N(Ef) in the equation for lambda is the DOS per spin
     dosef = dosef / two
     !
     WRITE (6, 100) degaussw0, ngaussw
     WRITE (6, 101) dosef, ef0 * ryd2ev
     
     !
     !
     ! loop over all q points of the fine mesh (this is in the k-para case 
     ! it should always be k-para for selfen_phon)
     DO iq = 1, nxqf
        !
        CALL start_clock('PH SELF-ENERGY')
        !
        fermicount = 0
        gamma = zero
        gamma_v = zero
        gamma_chk = zero
        !
        DO ik = 1, nksqf
           !
           ikk = 2 * ik - 1
           ikq = ikk + 1
           !
           coskkq = 0.d0
           DO ibnd = 1, nbndsub
              DO jbnd = 1, nbndsub
                 ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
                 ! 1/m  = 2 in Rydberg atomic units
                 vkk(:, ibnd ) = 2.0 * REAL (dmef (:, ibnd, ibnd, ikk ) )
                 vkq(:, jbnd ) = 2.0 * REAL (dmef (:, jbnd, jbnd, ikq ) )
                 IF ( abs ( DDOT(3,vkk(:,ibnd), 1, vkk(:,ibnd), 1) ) .gt. 1.d-4 ) &
                    coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1) / &
                                          DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
              ENDDO
           ENDDO
           !
           ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
           !
           IF (etf_mem) then
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
           IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
                ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) then
              !
              fermicount = fermicount + 1
              !
              DO imode = 1, nmodes
                 !
                 ! the phonon frequency
                 wq = wf (imode, iq)
                 !
                 !  we read the e-p matrix from disk / memory
                 !
                 IF (etf_mem) then
                    epf(:,:) = epf17 ( ik, iq, :, :, imode)
                 ELSE
                    nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
                    CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
                 ENDIF
                 !
                 DO ibnd = 1, ibndmax-ibndmin+1
                    !
                    !  the fermi occupation for k
                    ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                    wgkk = wgauss( -ekk/eptemp0, -99)
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
                       ! NON FUNZIONA SCAMBIANDO i,j
                       ! the coupling from Gamma acoustic phonons is negligible
                       IF ( wq .gt. eps_acustic ) THEN
                          g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                       ELSE
                          g2 = 0.d0
                       ENDIF
                       !
                       ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
                       ! This is the imaginary part of the phonon self-energy, sans the matrix elements
                       !
                       !weight = wkf (ikk) * (wgkk - wgkq) * &
                       !   aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) ) 
                       !
                       ! the below expression is positive-definite, but also an approximation
                       ! which neglects some fine features
                       !
                       w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                       w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                       weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                       !weight = pi * wq*wkf (ikk) *abs((wgkk-wgkq)/(ekq-ekk-wq+degaussw0))
                       !
                       !if (weight.gt. 1.0d-8) WRITE(6,'(/5x,"ikk = ",i5," ibnd: ", i5, " wt: ", f9.5)') ikk, ibnd, weight
                       !gamma(imode) =   gamma   (imode) + weight * g2 
                       !gamma_v(imode) = gamma_v (imode) + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
                       gamma(imode) =   gamma   (imode) + g2 * wgff(jbnd,ikk)!*wkf(ikk)*pi
                       gamma_v(imode) = gamma_v (imode) + weight * g2 * (1-coskkq(ibnd, jbnd) )
                       !WRITE(6,'(/5x,"wgf = ", f14.10," weight = ", f14.10)') wgff(jbnd,ikk)-wgff(ibnd,ikk),weight
                       !
                       !Fqnu( ik, imode, ibnd, jbnd ) = g2 * ( wgkk -wgkq)
                       !Fqnu( ik, imode, ibnd) = Fqnu( ik, imode, ibnd)+g2 *wgff(jbnd,ik)*( wgkq - wgkk )
                    ENDDO ! jbnd
                    !gamma(imode) =   gamma   (imode)* wgf(ibnd,ikk)
                    !gamma(imode) =   gamma   (imode) + g2 *wgff(ibnd,ik)*pi 
                 ENDDO   ! ibnd
                 !
              ENDDO ! loop on q-modes
              !
              !
           ENDIF ! endif fsthick
           !
           CALL stop_clock('PH SELF-ENERGY')
           !
        ENDDO ! loop on k 

        !DO jbnd = 1, ibndmax-ibndmin+1
        !  DO ibnd = 1, ibndmax-ibndmin+1
        !    DO itet= 1, ntetraf
         !      DO imode = 1, nmodes
        !        ! Gamma(ikk, ibnd, imode) = Gamma(ikk, ibnd, imode) + Fqnu( iq, imode, ibnd, jbnd ) * tet_weight(imode) * 2 * pi
                  !gamma_chk(imode) = gamma_chk(imode) + Fqnu( itet, imode, ibnd)*pi 
         !      ENDDO ! imode
         !    ENDDO ! itet
         ! ENDDO ! ibnd
        !ENDDO ! jbnd       !
#ifdef __PARA
        !
        ! collect contributions from all pools (sum over k-points)
        ! this finishes the integral over the BZ  (k)
        !
        CALL mp_sum(gamma,inter_pool_comm) 
        CALL mp_sum(gamma_v,inter_pool_comm) 
        CALL mp_sum(fermicount, inter_pool_comm)
        !
#endif
        !
        WRITE(6,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:,iq) , wqf(iq)
        WRITE(6,'(5x,a)') repeat('-',67)
        lambda_tot = zero
        DO imode = 1, nmodes
           ! 
           wq = wf (imode, iq)
           lambda = zero
           lambda_v = zero
           IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda   = gamma(imode) / pi / wq**two / dosef
           IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda_v = gamma_v(imode) / pi / wq**two / dosef
           lambda_tot = lambda_tot + lambda
           !
           gamma_all( imode, iq, ismear ) = gamma(imode)
           gamma_all_v( imode, iq, ismear ) = gamma_v(imode) / pi / wq**two / dosef
           lambda_all( imode, iq, ismear ) = lambda
           lambda_v_all( imode, iq, ismear ) = lambda_v
           !
           WRITE(6, 102) imode, lambda, ryd2mev * gamma(imode), ryd2mev * wq
          write (out_unit,*)  ryd2mev * wq, ryd2mev * gamma(imode)
  call flush(6)
        ENDDO
        !
        WRITE(6,103) lambda_tot
        WRITE(6,'(5x,a/)') repeat('-',67)
        !
        ! test ONLY
#ifdef __PARA
!        if (me.eq.1) & 
        IF (me_pool == 0) &
#endif
        !
        !     
        WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
             'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkstotf/2
        !
     ENDDO ! loop on q
     close (out_unit)
     open (unit=out_unit,file="smear.txt",action="write",status="replace")
     DO iq = 1, nxqf
        !
        CALL start_clock('PH SELF-ENERGY')
        !
        fermicount = 0
        gamma = zero
        gamma_v = zero
        !
        DO ik = 1, nksqf
           !
           ikk = 2 * ik - 1
           ikq = ikk + 1
           !
           coskkq = 0.d0
           DO ibnd = 1, nbndsub
              DO jbnd = 1, nbndsub
                 ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
                 ! 1/m  = 2 in Rydberg atomic units
                 vkk(:, ibnd ) = 2.0 * REAL (dmef (:, ibnd, ibnd, ikk ) )
                 vkq(:, jbnd ) = 2.0 * REAL (dmef (:, jbnd, jbnd, ikq ) )
                 IF ( abs ( DDOT(3,vkk(:,ibnd), 1, vkk(:,ibnd), 1) ) .gt. 1.d-4 ) &
                    coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1) / &
                                          DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
              ENDDO
           ENDDO
           !
           ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
           !
           IF (etf_mem) then
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
           IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
                ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) then
              !
              fermicount = fermicount + 1
              !
              DO imode = 1, nmodes
                 !
                 ! the phonon frequency
                 wq = wf (imode, iq)
                 !
                 !  we read the e-p matrix from disk / memory
                 !
                 IF (etf_mem) then
                    epf(:,:) = epf17 ( ik, iq, :, :, imode)
                 ELSE
                    nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
                    CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
                 ENDIF
                 !
                 DO ibnd = 1, ibndmax-ibndmin+1
                    !
                    !  the fermi occupation for k
                    ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                    wgkk = wgauss( -ekk/eptemp0, -99)
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
                       ! NON FUNZIONA SCAMBIANDO i,j
                       ! the coupling from Gamma acoustic phonons is negligible
                       IF ( wq .gt. eps_acustic ) THEN
                          g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                       ELSE
                          g2 = 0.d0
                       ENDIF
                       !
                       ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
                       ! This is the imaginary part of the phonon self-energy, sans the matrix elements
                       !
                       !weight = wkf (ikk) * (wgkk - wgkq) * &
                       !   aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) ) 
                       !
                       ! the below expression is positive-definite, but also an approximation
                       ! which neglects some fine features
                       !
                       w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                       w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                       weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                       !
                       gamma(imode) =   gamma   (imode) + weight * g2 
                       gamma_v(imode) = gamma_v (imode) + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
                       !
                    ENDDO ! jbnd
                 ENDDO   ! ibnd
                 !
              ENDDO ! loop on q-modes
              !
              !
           ENDIF ! endif fsthick
           !
           CALL stop_clock('PH SELF-ENERGY')
           !
        ENDDO ! loop on k
#ifdef __PARA
        !
        ! collect contributions from all pools (sum over k-points)
        ! this finishes the integral over the BZ  (k)
        !
        CALL mp_sum(gamma,inter_pool_comm) 
        CALL mp_sum(gamma_v,inter_pool_comm) 
        CALL mp_sum(fermicount, inter_pool_comm)
        !
#endif
        !
        WRITE(6,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:,iq) , wqf(iq)
        WRITE(6,'(5x,a)') repeat('-',67)
        lambda_tot = zero
        DO imode = 1, nmodes
           ! 
           wq = wf (imode, iq)
           lambda = zero
           lambda_v = zero
           IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda   = gamma(imode) / pi / wq**two / dosef
           IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda_v = gamma_v(imode) / pi / wq**two / dosef
           lambda_tot = lambda_tot + lambda
           !
           gamma_all( imode, iq, ismear ) = gamma(imode)
           gamma_all_v( imode, iq, ismear ) = gamma_v(imode) / pi / wq**two / dosef
           lambda_all( imode, iq, ismear ) = lambda
           lambda_v_all( imode, iq, ismear ) = lambda_v
           !
           WRITE(6, 102) imode, lambda, ryd2mev * gamma(imode), ryd2mev * wq
           write (out_unit,*)  ryd2mev * wq, ryd2mev * gamma(imode)
  call flush(6)
        ENDDO
        !
        WRITE(6,103) lambda_tot
        WRITE(6,'(5x,a/)') repeat('-',67)
        !
        ! test ONLY
#ifdef __PARA
!        if (me.eq.1) & 
        IF (me_pool == 0) &
#endif
        !
        !     
        WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
             'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkstotf/2
        !
     ENDDO ! loop on q
     close (out_unit)

  !ENDDO !smears
  ! generate the Eliashberg spectral function
  !
  IF (a2f) call eliashberg_a2f( lambda_all(:,:,1), lambda_v_all(:,:,1))
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 format(5x,'lambda( ',i3,' )=',f9.3,'   gamma=',f9.3,' meV','   omega=',f9.3,' meV')
103 format(5x,'lambda( tot )=',f9.3)
  !
  end subroutine selfen_phon_tet
  !
  !-----------------------------------------------------------------------
  subroutine selfen_phon_fly_tet (iq )
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation.  This routine is similar to the one above
  !  but it is ONLY called from within ephwann_shuffle and calculates 
  !  the selfenergy for one phonon at a time.  Much smaller footprint on
  !  the disk
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  use phcom,     ONLY : nmodes
  use epwcom,    ONLY : nbndsub, lrepmatf, iunepmatf, fsthick, &
                        eptemp, ngaussw, degaussw, iuetf,     &
                        wmin, wmax, nw, nbndskip, epf_mem, etf_mem, &
                        nsmear, delta_smear, eig_read, eps_acustic
  use pwcom,     ONLY : nelec, ef, isk, nbnd
  use el_phon,   ONLY : epf17, ibndmax, ibndmin, etf, &
                        etfq, wkf, xqf, wqf, nksf, nxqf,   &
                        nksqf, wf, nkstotf, xkf, xqf, &
                        lambda_all, lambda_v_all, nrr_k, &
                        dmef, ndegen_k, irvec
#ifdef __PARA
  use mp,        ONLY : mp_barrier,mp_sum
  use mp_global, ONLY : me_pool,inter_pool_comm,my_pool_id
#endif
  !
  implicit none
  !
  real(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0), czero=(0.d0,0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount, ismear
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgkk, wgkq, lambda, lambda_v, &
     weight, wgauss, dosef, dos_ef, w0g1, w0g2, w0gauss,   &
     gamma(nmodes),gamma_v(nmodes), degaussw0, eptemp0, lambda_tot
  !
  real(kind=DP), external :: efermig
  real(kind=DP) :: coskkq(nbndsub, nbndsub)
  real(kind=DP) :: DDOT,  vkk(3,nbndsub), vkq(3,nbndsub)
  !
  !
  IF (iq.eq.1) then 
     WRITE(6,'(/5x,a)') repeat('=',67)
     WRITE(6,'(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation (on the fly)")') 
     WRITE(6,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick.lt.1.d3 ) &
          WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Fermi Surface thickness = ', fsthick, ' Ry'
     WRITE(stdout, '(/5x,a,e18.9,a)' ) &
          'Golden Rule strictly enforced with T = ',eptemp(1), ' Ry'
     !
     ! here we take into account that we may skip bands when we wannierize
     ! (spin-unpolarized)
     ! 
     ! could be done for iq = 1
     IF (.not.ALLOCATED (lambda_all) .and. .not.ALLOCATED(lambda_v_all)) then
        ALLOCATE(lambda_all(nmodes, nxqf, nsmear))
        ALLOCATE(lambda_v_all(nmodes, nxqf, nsmear))
        lambda_all   = zero
        lambda_v_all = zero
     ENDIF
     !
     IF (nbndskip.gt.0) then
        nelec = nelec - two * nbndskip
        WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
        WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec     
     ENDIF
     !
  ENDIF
  !
  DO ismear = 1, nsmear
     !
     degaussw0 = (ismear-1) * delta_smear + degaussw
     eptemp0   = (ismear-1) * delta_smear + eptemp(1)
     !
     ! Fermi level and corresponding DOS
     !
     !   Note that the weights of k+q points must be set to zero here
     !   no spin-polarized calculation here
     ef0 = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw0,ngaussw,0,isk)
     dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
     !
     !   N(Ef) in the equation for lambda is the DOS per spin
     dosef = dosef / two
     !
     IF ( iq .eq. 1 ) THEN 
        WRITE (6, 100) degaussw0, ngaussw
        WRITE (6, 101) dosef, ef0 * ryd2ev
        WRITE (6, 101) dosef, ef  * ryd2ev
     ENDIF
     !   if 'fine' Fermi level differs by more than 250 meV, there is probably
     !   something wrong with the wannier functions
     IF (abs(ef0 - ef) * ryd2eV .gt. 0.5 .and. (.not.eig_read) ) &
        CALL errore ('selfen_phon', 'Something wrong with Ef, check MLWFs', 1)
     !
     CALL start_clock('PH SELF-ENERGY')
     !
     fermicount = 0
     gamma = zero
     gamma_v = zero
     !
     DO ik = 1, nksqf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        ! 
        coskkq = 0.d0
        DO ibnd = 1, nbndsub
           DO jbnd = 1, nbndsub
              ! coskkq = (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
              ! this is different from :   coskkq = (vk dot vkq) / |vk||vkq|
              ! In principle the only coskkq contributing to lambda_tr are both near the
              ! Fermi surface and the magnitudes will not differ greatly between vk and vkq
              ! we may implement the approximation to the angle between k and k+q vectors also 
              ! listed in Grimvall
              !
              ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
              vkk(:, ibnd ) = 2.0 * REAL (dmef (:, ibnd, ibnd, ikk ) )
              vkq(:, jbnd ) = 2.0 * REAL (dmef (:, jbnd, jbnd, ikq ) )
              if ( abs ( vkk(1,ibnd)**2 + vkk(2,ibnd)**2 + vkk(3,ibnd)**2) .gt. 1.d-4) &
                   coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1)  / &
                   DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
           ENDDO
        ENDDO
        !
        ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
        !
        ! when we see references to iq for file readinq, it is always = 1 for on the fly calculations
        IF (etf_mem) then
           etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk,  1)
           etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq,  1)
        ELSE
           nrec = (iq-1) * nksf + ikk
           nrec = ikk
           CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
           nrec = (iq-1) * nksf + ikq
           nrec = ikq
           CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        ENDIF
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
             ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) then
           !
           fermicount = fermicount + 1
           !
           DO imode = 1, nmodes
              !
              ! the phonon frequency
              wq = wf (imode, iq)
              !
              !  we read the e-p matrix from disk / memory
              !
              IF (etf_mem) then
                 epf(:,:) = epf17 ( ik,  1, :, :, imode)
              ELSE
                 nrec = (imode-1) * nksqf + ik
                 CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
              ENDIF
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                 !
                 !  the fermi occupation for k
                 ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                 wgkk = wgauss( -ekk/eptemp0, -99)
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
                    ! NON FUNZIONA SCAMBIANDO i,j
                    ! the coupling from Gamma acoustic phonons is negligible
                    IF ( wq .gt. eps_acustic ) THEN
                       g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                    ELSE
                       g2 = 0.d0
                    ENDIF
                    !
                    ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
                    ! This is the imaginary part of the phonon self-energy, sans the matrix elements
                    !
                    !weight = wkf (ikk) * (wgkk - wgkq) * &
                    !     aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) ) 
                    !
                    ! the below expression is positive-definite, but also an approximation
                    ! which neglects some fine features
                    !
                    w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                    w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                    weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                    !
                    gamma   (imode) =   gamma   (imode) + weight * g2 
                    gamma_v (imode) =   gamma_v (imode) + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
                    !
                    !
                 ENDDO ! jbnd
              ENDDO   ! ibnd
              !
           ENDDO ! loop on q-modes
           !
           !
        ENDIF ! endif fsthick
        !
        CALL stop_clock('PH SELF-ENERGY')
        !
     ENDDO ! loop on k
#ifdef __PARA
     !
     ! collect contributions from all pools (sum over k-points)
     ! this finishes the integral over the BZ  (k)
     !
     CALL mp_sum(gamma,inter_pool_comm) 
     CALL mp_sum(gamma_v,inter_pool_comm) 
     CALL mp_sum(fermicount, inter_pool_comm)
     !
#endif
     !
     WRITE(6,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5, " ismear: ",i3)') iq, xqf(:,iq) , wqf(iq), ismear
     WRITE(6,'(5x,a)') repeat('-',67)
     lambda_tot = 0
     DO imode = 1, nmodes
        ! 
        wq = wf (imode, iq)
        lambda   = zero
        lambda_v = zero
        IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda   = gamma(imode)   / pi / wq**two / dosef
        IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda_v = gamma_v(imode) / pi / wq**two / dosef
        lambda_tot = lambda_tot + lambda
        !
        lambda_all( imode, iq, ismear ) = lambda
        lambda_v_all( imode, iq, ismear ) = lambda_v
        !
        !
        WRITE(6, 102) imode, lambda, ryd2mev * gamma(imode), ryd2mev * wq
     ENDDO
     WRITE(6, 103) lambda_tot
     WRITE(6,'(5x,a/)') repeat('-',67)
     !

     IF (ismear .eq. 1) write( stdout, '(/5x,a,i8,a,i8/)' ) &
          'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkstotf/2
     !
     !
  ENDDO !smears
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 format(5x,'lambda( ',i3,' )=',f9.3,'   gamma=',f9.3,' meV','   omega=',f9.3,' meV')
103 format(5x,'lambda( sum )=',f9.3)
  !
  return
end subroutine selfen_phon_fly_tet
!

