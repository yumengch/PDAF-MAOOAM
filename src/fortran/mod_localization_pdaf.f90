! This file is part of pyPDAF

! Copyright (C) 2022 University of Reading and
! National Centre for Earth Observation

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module mod_localization_pdaf
use mod_kind_pdaf, only: wp
implicit none

integer :: loc_weight
real(wp) :: local_range
real(wp) :: srange
integer, allocatable :: id_lstate_in_pstate(:)
real(wp) :: coords_l(2)
namelist /local_nml/ loc_weight, local_range, srange
contains
   SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)
      USE mod_model_pdaf, &             ! Model variables
         ONLY: dim_state_p, nx, ny, pi, maooam_model
      USE mod_parallel_pdaf, &     ! assimilation parallelization variables
         ONLY: mype_filter

      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in)  :: step     !< Current time step
      INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
      INTEGER, INTENT(out) :: dim_l    !< Local state dimension

      ! *** local variables ***
      INTEGER :: i                       ! Counters
      INTEGER :: off_p                   ! Process-local offset in global state vector
      real(wp) :: n
      real(wp) :: dx, dy

      ! ****************************************
      ! *** Initialize local state dimension ***
      ! ****************************************

      dim_l = 1

      ! **********************************************
      ! *** Initialize coordinates of local domain ***
      ! **********************************************

      ! Global coordinates of local analysis domain
      ! We use grid point indices as coordinates, but could e.g. use meters
      off_p = 0
      DO i = 1, mype_filter
         off_p = off_p + dim_state_p
      END DO
      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)
      coords_l(2) = int((mod(domain_p - 1, nx*ny)) / nx)
      coords_l(1) = (int(mod(domain_p - 1, nx*ny) -  coords_l(2)*nx))*dx
      coords_l(2) = coords_l(2)*dy
      ! ******************************************************
      ! *** Initialize array of indices of the local state ***
      ! ***  vector elements in the global state vector.   ***
      ! ******************************************************

      ! Allocate array
      IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
      ALLOCATE(id_lstate_in_pstate(dim_l))

      ! Here the local domain is a single grid point and variable given by DOMAIN_P
      id_lstate_in_pstate(1) = domain_p
   END SUBROUTINE init_dim_l_pdaf

   SUBROUTINE init_n_domains_pdaf(step, n_domains_p)
      USE mod_model_pdaf, &      ! Assimilation variables
         ONLY: dim_state_p

      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in)  :: step        !< Current time step
      INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains

      ! ************************************
      ! *** Initialize number of domains ***
      ! ************************************

      ! Here simply the state dimension
      n_domains_p = dim_state_p
   END SUBROUTINE init_n_domains_pdaf

   SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)
      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in) :: step           !< Current time step
      INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
      INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
      INTEGER, INTENT(in) :: dim_l          !< Local state dimension
      REAL(wp), INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
      REAL(wp), INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain

      ! *** local variables ***
      INTEGER :: i                          ! Counter

      ! *************************************
      ! *** Initialize local state vector ***
      ! *************************************

      ! Generic initialization using ID_LSTATE_IN_PSTATE set in INIT_DIM_L_PDAF
      DO i = 1, dim_l
         state_l(i) = state_p(id_lstate_in_pstate(i))
      END DO
   END SUBROUTINE g2l_state_pdaf

   SUBROUTINE l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p)
      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in) :: step           !< Current time step
      INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
      INTEGER, INTENT(in) :: dim_l          !< Local state dimension
      INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
      REAL(wp), INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
      REAL(wp), INTENT(inout) :: state_p(dim_p) !< PE-local full state vector 

      ! *** local variables ***
      INTEGER :: i                          ! Counter


      ! **************************************************
      ! *** Initialize elements of global state vector ***
      ! **************************************************

      ! Generic initialization using ID_LSTATE_IN_PSTATE set in INIT_DIM_L_PDAF
      DO i = 1, dim_l
         state_p(id_lstate_in_pstate(i)) = state_l(i)
      END DO
   END SUBROUTINE l2g_state_pdaf

   subroutine finalize
      if (allocated(id_lstate_in_pstate)) deallocate(id_lstate_in_pstate)
   end subroutine finalize
end module mod_localization_pdaf
