!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine q2qstar_ph (dyn, at, bg, nat, nsym, s, invs, irt, rtau, &
     nq, sxq, isq, imq, iudyn)
  !-----------------------------------------------------------------------
  ! Generates the dynamical matrices for the star of q and writes them on
  ! disk for later use.
  ! If there is a symmetry operation such that q -> -q +G then imposes on
  ! dynamical matrix those conditions related to time reversal symmetry.
  !
#include "f_defs.h"
  USE kinds, only : DP
  implicit none
  ! input variables
  integer :: nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
       nq, isq (48), imq, iudyn
  ! number of atoms in the unit cell
  ! number of symmetry operations
  ! the symmetry operations
  ! index of the inverse operations
  ! index of the rotated atom
  ! degeneracy of the star of q
  ! symmetry op. giving the rotated q
  ! index of -q in the star (0 if non present)
  ! unit number
  complex(DP) :: dyn (3 * nat, 3 * nat)
  ! the input dynamical matrix. if imq.ne.0 the
  ! output matrix is symmetrized w.r.t. time-reversal

  real(DP) :: at (3, 3), bg (3, 3), rtau (3, 48, nat), sxq (3, 48)
  ! direct lattice vectors
  ! reciprocal lattice vectors
  ! for each atom and rotation gives the R vector involved
  ! list of q in the star
  !
  !  local variables
  integer :: na, nb, iq, nsq, isym, icar, jcar, i, j
  ! counters
  ! nsq: number of sym.op. giving each q in the list

  complex(DP) :: phi (3, 3, nat, nat), phi2 (3, 3, nat, nat)
  ! work space
  !
  ! Sets number of symmetry operations giving each q in the list
  !
  nsq = nsym / nq
  if (nsq * nq /= nsym) call errore ('q2star_ph', 'wrong degeneracy', 1)
  !
  ! Writes dyn.mat. dyn(3*nat,3*nat) on the 4-index array phi(3,3,nat,nta)
  !
  do i = 1, 3 * nat
     na = (i - 1) / 3 + 1
     icar = i - 3 * (na - 1)
     do j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcar = j - 3 * (nb - 1)
        phi (icar, jcar, na, nb) = dyn (i, j)
     enddo
  enddo
  !
  ! Go to crystal coordinates
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi (1, 1, na, nb), at, bg, - 1)
     enddo
  enddo
  !
  ! If -q is in the list impose first of all the conditions coming from
  ! time reversal symmetry
  !
  if (imq /= 0) then
     phi2 (:,:,:,:) = (0.d0, 0.d0)
     isym = 1
     do while (isq (isym) /= imq)
        isym = isym + 1
     enddo
     call rotate_and_add_dyn (phi, phi2, nat, isym, s, invs, irt, &
          rtau, sxq (1, imq) )
     do na = 1, nat
        do nb = 1, nat
           do i = 1, 3
              do j = 1, 3
                 phi (i, j, na, nb) = 0.5d0 * (phi (i, j, na, nb) + &
                                        CONJG(phi2(i, j, na, nb) ) )
              enddo
           enddo
        enddo
     enddo
     phi2 (:,:,:,:) = phi (:,:,:,:)
     !
     ! Back to cartesian coordinates
     !
     do na = 1, nat
        do nb = 1, nat
           call trntnsc (phi2 (1, 1, na, nb), at, bg, + 1)
        enddo
     enddo
     !
     ! Saves 4-index array phi(3,3,nat,nta) on the dyn.mat. dyn(3*nat,3*nat)
     !
     do i = 1, 3 * nat
        na = (i - 1) / 3 + 1
        icar = i - 3 * (na - 1)
        do j = 1, 3 * nat
           nb = (j - 1) / 3 + 1
           jcar = j - 3 * (nb - 1)
           dyn (i, j) = phi2 (icar, jcar, na, nb)
        enddo
     enddo
  endif
  !
  ! For each q of the star rotates phi with the appropriate sym.op. -> phi
  !
  do iq = 1, nq
     phi2 (:,:,:,:) = (0.d0, 0.d0)
     do isym = 1, nsym
        if (isq (isym) == iq) then
           call rotate_and_add_dyn (phi, phi2, nat, isym, s, invs, irt, &
                rtau, sxq (1, iq) )
        endif
     enddo
     phi2 (:,:,:,:) = phi2 (:,:,:,:) / DBLE (nsq)
     !
     ! Back to cartesian coordinates
     !
     do na = 1, nat
        do nb = 1, nat
           call trntnsc (phi2 (1, 1, na, nb), at, bg, + 1)
        enddo
     enddo
     !
     ! Writes the dynamical matrix in cartesian coordinates on file
     !
     call write_dyn_on_file (sxq (1, iq), phi2, nat, iudyn)
     if (imq == 0) then
        !
        ! if -q is not in the star recovers its matrix by time reversal
        !
        do na = 1, nat
           do nb = 1, nat
              do i = 1, 3
                 do j = 1, 3
                    phi2 (i, j, na, nb) = CONJG(phi2 (i, j, na, nb) )
                 enddo
              enddo
           enddo
        enddo
        !
        ! and writes it (changing temporarily sign to q)
        !
        sxq (:, iq) = - sxq (:, iq)
        call write_dyn_on_file (sxq (1, iq), phi2, nat, iudyn)
        sxq (:, iq) = - sxq (:, iq)
     endif
  enddo
  !
  return
end subroutine q2qstar_ph
