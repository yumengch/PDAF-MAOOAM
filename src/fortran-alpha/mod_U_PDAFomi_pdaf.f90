! this file is part of pypdaf

! copyright (c) 2022 university of reading and
! national centre for earth observation

! this program is free software: you can redistribute it and/or modify
! it under the terms of the gnu general public license as published by
! the free software foundation, either version 3 of the license, or
! (at your option) any later version.

! this program is distributed in the hope that it will be useful,
! but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose.  see the
! gnu general public license for more details.

! you should have received a copy of the gnu general public license
! along with this program.  if not, see <http://www.gnu.org/licenses/>.

module mod_u_pdafomi_pdaf
use mod_kind_pdaf, only: wp
use mod_observations_pdaf, only: n_obs, obs
implicit none

integer :: timer_getobs_start, timer_getobs_end, t_rate
integer :: timer_dimomi_start, timer_dimomi_end
integer :: timer_op_start, timer_op_end
real(wp) :: getobs_dur, dimomi_dur, op_dur

integer :: op_cnt

contains
   !! set up the observation information and obtain the
   !! size of the observation vector for assimilation
   subroutine init_dim_obs_pdafomi(step, dim_obs)
      use mod_observations_pdaf, only: init_dim_obs

      integer, intent(in)  :: step     !< current time step
      integer, intent(out) :: dim_obs  !< dimension of full observation vector

      integer :: this_dim_obs
      integer :: i_obs

      call system_clock(timer_dimomi_start)

      dim_obs = 0
      do i_obs = 1, n_obs
         if (mod(step, obs(i_obs)%delt_obs) == 0) then
            call init_dim_obs(i_obs, step, this_dim_obs)
            dim_obs = dim_obs + this_dim_obs
         end if
      end do

      call system_clock(timer_dimomi_end, t_rate)
      dimomi_dur = dimomi_dur + &
          (real(timer_dimomi_end, wp) - real(timer_dimomi_start, wp))/real(t_rate, wp)
   end subroutine init_dim_obs_pdafomi

   ! observation operator for obs. on gridpoints
   subroutine obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate, ismean)
      ! include functions for different observations
      use mod_observations_pdaf, only: obs_op_A, obs_op_Aonly, obs_op_O
      use mod_statevector_pdaf, only: update_ocean, update_both
      use mod_U_pdaf, only: ens_p_noupdate
      implicit none
      ! *** arguments ***
      integer,  intent(in)    :: step                 !< current time step
      integer,  intent(in)    :: dim_p                !< pe-local state dimension
      integer,  intent(in)    :: dim_obs              !< dimension of full observed state
      real(wp), intent(in)    :: state_p(dim_p)       !< pe-local model state
      real(wp), intent(out) :: ostate(dim_obs)      !< pe-local full observed state
      logical, intent(in)  :: ismean

      ! ******************************************************
      ! *** apply observation operator h on a state vector ***
      ! ******************************************************
      call system_clock(timer_op_start)
      if (update_both) then
         op_cnt = op_cnt + 1
         if (update_ocean) then
            call obs_op_O(dim_p, dim_obs, state_p, ostate, ens_p_noupdate, op_cnt)
         else
            call obs_op_A(dim_p, dim_obs, state_p, ostate, ens_p_noupdate, op_cnt)
         end if
      else
         call obs_op_Aonly(dim_p, dim_obs, state_p, ostate)
      end if
      call system_clock(timer_op_end, t_rate)
      op_dur = op_dur + &
          (real(timer_op_end, wp) - real(timer_op_start, wp))/real(t_rate, wp)
   end subroutine obs_op_pdafomi

   SUBROUTINE prodRinvA_pdafomi(step, dim_obs_p, dim_ens, obs_p, A_p, C_p)
      use mod_observations_pdaf, only: prodRinvA_A, prodRinvA_Aonly, prodRinvA_O
      use mod_statevector_pdaf, only: update_ocean, update_both
      INTEGER,  INTENT(in)  :: step                    ! Current time step
      INTEGER,  INTENT(in)  :: dim_obs_p               ! PE-local dimension of obs. vector
      INTEGER,  INTENT(in)  :: dim_ens                 ! Ensemble size
      REAL(wp), INTENT(in)  :: obs_p(dim_obs_p)        ! PE-local vector of observations
      REAL(wp), INTENT(in)  :: A_p(dim_obs_p, dim_ens) ! Input matrix from analysis routine
      REAL(wp), INTENT(out) :: C_p(dim_obs_p, dim_ens) ! Output matrix
      if (update_both) then
         if (update_ocean) then
            call prodRinvA_O(step, dim_obs_p, dim_ens, obs_p, A_p, C_p)
         else
            call prodRinvA_A(step, dim_obs_p, dim_ens, obs_p, A_p, C_p)
         end if
      else
         call prodRinvA_Aonly(step, dim_obs_p, dim_ens, obs_p, A_p, C_p)
      end if
   END SUBROUTINE prodRinvA_pdafomi

   SUBROUTINE init_obs_pdafomi(step, dim_obs_p, observation_p)
      use mod_observations_pdaf, only: init_obs, n_obs
      INTEGER,  INTENT(in)   :: step             ! Current time step
      INTEGER,  INTENT(in)   :: dim_obs_p        ! PE-local dimension of obs. vector
      REAL(wp), INTENT(out)  :: observation_p(dim_obs_p) ! PE-local observation vector

      integer :: i_obs
      do i_obs = 1, n_obs
         if (mod(step, obs(i_obs)%delt_obs) == 0) &
            call init_obs(i_obs, step, dim_obs_p, observation_p)
      end do
   end subroutine init_obs_pdafomi

   SUBROUTINE init_obsvar_pdafomi(step, dim_obs_p, obs_p, meanvar)
      INTEGER, INTENT(in) :: step          ! Current time step
      INTEGER, INTENT(in) :: dim_obs_p     ! PE-local dimension of observation vector
      REAL, INTENT(in) :: obs_p(dim_obs_p) ! PE-local observation vector
      REAL, INTENT(out)   :: meanvar       ! Mean observation error variance
   END SUBROUTINE init_obsvar_pdafomi
end module mod_u_pdafomi_pdaf

