!
! Copyright (C) 2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
PROGRAM read_upf_tofile 
  !---------------------------------------------------------------------
  !
  !   This small program reads the pseudopotential in the Unified 
  !   Pseudopotential Format and writes three files
  !   in a format which can be plotted. The files are:
  !
  !   filewfc    with the pseudo-wavefunctions
  !   filebeta   with the beta functions
  !   filepot    with the local potential, the valence and core charge.
  !
  !
  ! PWSCF modules
  !
  !
  USE constants, only : fpi
  USE pseudo_types
  USE upf_module
  USE radial_grids, only : radial_grid_type, nullify_radial_grid
  !
  IMPLICIT NONE
  !
  INTEGER :: iunps, ierr
  !
  CHARACTER(30) :: file_pseudo
  !
  !     Local variables
  !
  INTEGER :: ios, n, j
  TYPE (pseudo_upf) :: upf
  TYPE (radial_grid_type) :: grid
  !
  WRITE(6,'("Name of the upf file > ", $)') 
  READ(5,'(a)') file_pseudo 

  !  nullify objects as soon as they are instantiated
 
  CALL nullify_pseudo_upf( upf )
  CALL nullify_radial_grid( grid )

  iunps=2
  OPEN(UNIT=iunps,FILE=file_pseudo,STATUS='old',FORM='formatted', &
       ERR=100, IOSTAT=ios)
100   CALL errore('read_upf_tofile','open error on file '//file_pseudo,ios)

  CALL read_upf(upf, grid, ierr, unit=iunps)
  !
  IF (ierr .NE. 0) &
     CALL errore('read_upf_tofile','reading pseudo upf', ABS(ierr))
  !
  CLOSE(iunps)
  !
  OPEN(UNIT=iunps,FILE='filewfc',STATUS='unknown',FORM='formatted', &
       ERR=200, IOSTAT=ios)
200   CALL errore('read_upf_tofile','open error on file filewfc',ABS(ios))

  DO n=1,upf%mesh
     WRITE(iunps,'(30f12.6)') upf%r(n), (upf%chi(n,j), j=1,upf%nwfc)
  END DO

  CLOSE(iunps)

  OPEN(UNIT=iunps,FILE='filebeta',STATUS='unknown',FORM='formatted', &
       ERR=300, IOSTAT=ios)
300   CALL errore('read_upf_tofile','open error on file filebeta',ABS(ios))

  DO n=1,upf%mesh
     WRITE(iunps,'(30f12.6)') upf%r(n), (upf%beta(n,j), j=1,upf%nbeta)
  END DO

  CLOSE(iunps)

  OPEN(UNIT=iunps,FILE='filepot',STATUS='unknown',FORM='formatted', &
       ERR=400, IOSTAT=ios)
400   CALL errore('read_upf_tofile','open error on file filepot',ABS(ios))

  DO n=1,upf%mesh
     WRITE(iunps,'(4f12.6)') upf%r(n), upf%vloc(n), &
               upf%rho_at(n), upf%rho_atc(n)*fpi*upf%r(n)**2 
  END DO

  CLOSE(iunps)

  CALL deallocate_pseudo_upf( upf )

END PROGRAM read_upf_tofile
