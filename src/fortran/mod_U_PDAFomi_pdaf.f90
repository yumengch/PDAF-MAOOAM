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

module mod_U_PDAFomi_pdaf
use mod_kind_pdaf, only: wp
use mod_observations_pdaf, only: n_obs
implicit none

contains
   subroutine get_obs_f(step, dim_obs_f, observation_f)
      use mod_model_pdaf, only: dim_state
      use mod_obswriter_pdaf, only: writeObs
      INTEGER, INTENT(in)  :: step                ! Current time step
      INTEGER, INTENT(in)  :: dim_obs_f           ! Dimension of obs. vector
      REAL(wp), INTENT(in) :: observation_f(dim_obs_f) ! Observation vector

      integer :: i_obs, istart, iend

      istart = 1
      do i_obs = 1, n_obs
         iend = istart + dim_state - 1
         call writeObs(i_obs, step, observation_f(istart:iend))
         istart = iend + 1
      end do
   end subroutine get_obs_f


   subroutine init_dim_obs_gen_pdafomi(step, dim_obs)
      use mod_observations_pdaf, only: init_dim_obs_gen
      INTEGER, INTENT(in)  :: step     !< Current time step
      INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

      integer :: this_dim_obs
      integer :: i_obs

      dim_obs = 0
      
      do i_obs = 1, n_obs
         call init_dim_obs_gen(i_obs, this_dim_obs)
         dim_obs = dim_obs + this_dim_obs
      end do
   end subroutine init_dim_obs_gen_pdafomi


   subroutine init_dim_obs_pdafomi(step, dim_obs)
      use mod_observations_pdaf, only: init_dim_obs
      INTEGER, INTENT(in)  :: step     !< Current time step
      INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

      integer :: this_dim_obs
      integer :: i_obs

      dim_obs = 0
      
      do i_obs = 1, n_obs
         call init_dim_obs(i_obs, this_dim_obs)
         dim_obs = dim_obs + this_dim_obs
      end do
   end subroutine init_dim_obs_pdafomi


   SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)
      ! Include functions for different observations
      USE mod_observations_pdaf, ONLY: obs_op_gridpoint

      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in)     :: step                 !< Current time step
      INTEGER, INTENT(in)     :: dim_p                !< PE-local state dimension
      INTEGER, INTENT(in)     :: dim_obs              !< Dimension of full observed state
      REAL(wp), INTENT(in)    :: state_p(dim_p)       !< PE-local model state
      REAL(wp), INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state

      integer :: i_obs
      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************

      ! The order of these calls is not relevant as the setup
      ! of the overall observation vector is defined by the
      ! order of the calls in init_dim_obs_pdafomi
      do i_obs = 1, n_obs
         CALL obs_op_gridpoint(i_obs, dim_p, dim_obs, state_p, ostate)
      end do
   END SUBROUTINE obs_op_pdafomi

   SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)
      ! Include functions for different observations
      USE mod_observations_pdaf, ONLY: init_dim_obs_l

      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
      INTEGER, INTENT(in)  :: step       !< Current time step
      INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
      INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

      integer :: i_obs
      ! **********************************************
      ! *** Initialize local observation dimension ***
      ! **********************************************

      ! Call init_dim_obs_l specific for each observation
      do i_obs = 1, n_obs
         CALL init_dim_obs_l(i_obs, dim_obs_l)
      enddo

   END SUBROUTINE init_dim_obs_l_pdafomi

   SUBROUTINE localize_covar_pdafomi(dim_p, dim_obs, HP_p, HPH)

      ! Include functions for different observations
      USE mod_observations_pdaf, ONLY: localize_covar
      USE mod_model_pdaf, &          ! Include information on model grid
          ONLY: dim_state_p
      USE mod_parallel_pdaf, &  ! Include rank of filter process
          ONLY: mype_filter

      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in)     :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in)     :: dim_obs               !< number of observations
      REAL(wp), INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
      REAL(wp), INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH

      ! *** local variables ***
      INTEGER :: i_obs, j, cnt_p         ! Counters
      INTEGER :: off_nx                  ! Offset in x-direction for parallelization
      REAL(wp), ALLOCATABLE :: coords_p(:,:) ! Coordinates of PE-local state vector entries


      ! **********************
      ! *** INITIALIZATION ***
      ! **********************

      ! *** Initialize coordinate array ***

      ALLOCATE(coords_p(1, dim_p))

      ! Get offset of local domain in global domain in x-direction
      off_nx = 0
      DO j = 1, mype_filter
         off_nx = off_nx + dim_state_p
      END DO

      cnt_p = 0
      DO j = 1 + off_nx, dim_state_p + off_nx
         cnt_p = cnt_p + 1
         coords_p(1, cnt_p) = REAL(j)
      END DO


      ! *************************************
      ! *** Apply covariance localization ***
      ! *************************************

      ! Call localize_covar specific for each observation
      do i_obs = 1, n_obs
         CALL localize_covar(i_obs, dim_p, dim_obs, HP_p, HPH, coords_p)
      enddo


      ! ****************
      ! *** Clean up ***
      ! ****************

      DEALLOCATE(coords_p)

   END SUBROUTINE localize_covar_pdafomi
end module mod_U_PDAFomi_pdaf