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
use mod_config_pdaf, only: is_strong
use mod_statevector_pdaf, only: component
implicit none

integer :: timer_getobs_start, timer_getobs_end, t_rate
integer :: timer_dimomi_start, timer_dimomi_end
integer :: timer_op_start, timer_op_end
real(wp) :: getobs_dur, dimomi_dur, op_dur

contains
   ! write the synthetic observations into netcdf file
   subroutine get_obs_f(step, dim_obs_f, observation_f)
      use mod_obswriter_pdaf, only: writeobs
      ! arguments
      integer,  intent(in)  :: step                     ! current time step
      integer,  intent(in)  :: dim_obs_f                ! dimension of obs. vector
      real(wp), intent(in)  :: observation_f(dim_obs_f) ! observation vector
      ! local variables
      integer :: i_obs, istart, iend

      call system_clock(timer_getobs_start)

      istart = 1
      do i_obs = 1, n_obs
         iend = istart + obs(i_obs)%dim_obs - 1
         call writeobs(i_obs, step, observation_f(istart:iend))
         istart = iend + 1
      end do
      call system_clock(timer_getobs_end, t_rate)
      getobs_dur = getobs_dur + &
          (real(timer_getobs_end, wp) - real(timer_getobs_start, wp))/real(t_rate, wp)
   end subroutine get_obs_f

   !! set up the observation information and obtain the 
   !! size of the observation vector for synthetic obs.
   subroutine init_dim_obs_gen_pdafomi(step, dim_obs)
      use mod_observations_pdaf, only: init_dim_obs_gen
      ! arguments
      integer, intent(in)  :: step     !< current time step
      integer, intent(out) :: dim_obs  !< dimension of full observation vector
      ! local variables
      integer :: this_dim_obs
      integer :: i_obs

      call system_clock(timer_dimomi_start)

      dim_obs = 0
      do i_obs = 1, n_obs
         call init_dim_obs_gen(i_obs, this_dim_obs)
         dim_obs = dim_obs + this_dim_obs
      end do

      call system_clock(timer_dimomi_end, t_rate)
      dimomi_dur = dimomi_dur + &
          (real(timer_dimomi_end, wp) - real(timer_dimomi_start, wp))/real(t_rate, wp)
   end subroutine init_dim_obs_gen_pdafomi

   !! set up the observation information and obtain the 
   !! size of the observation vector for assimilation
   subroutine init_dim_obs_pdafomi(step, dim_obs)
      use mod_observations_pdaf, only: init_dim_obs
      integer, intent(in)  :: step     !< current time step
      integer, intent(out) :: dim_obs  !< dimension of full observation vector

      integer :: this_dim_obs
      integer :: i_obs
      integer :: true_step, shift, factor
      logical :: condition

      call system_clock(timer_dimomi_start)

      dim_obs = 0
      condition = (.not. is_strong) .and. (component == 'b')
      factor = 1
      shift = 0
      if (condition) then
         factor = 2
         shift = -1
      end if

      true_step = (step - shift)/factor
      do i_obs = 1, n_obs
         call init_dim_obs(i_obs, step, this_dim_obs)
         dim_obs = dim_obs + this_dim_obs
      end do

      call system_clock(timer_dimomi_end, t_rate)
      dimomi_dur = dimomi_dur + &
          (real(timer_dimomi_end, wp) - real(timer_dimomi_start, wp))/real(t_rate, wp)
   end subroutine init_dim_obs_pdafomi

   ! observation operator for obs. on gridpoints
   subroutine obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)
      ! include functions for different observations
      use mod_observations_pdaf, only: obs_op_gridpoint
      implicit none
      ! *** arguments ***
      integer,  intent(in)    :: step                 !< current time step
      integer,  intent(in)    :: dim_p                !< pe-local state dimension
      integer,  intent(in)    :: dim_obs              !< dimension of full observed state
      real(wp), intent(in)    :: state_p(dim_p)       !< pe-local model state
      real(wp), intent(inout) :: ostate(dim_obs)      !< pe-local full observed state
      ! local variables
      integer :: i_obs

      ! ******************************************************
      ! *** apply observation operator h on a state vector ***
      ! ******************************************************
      call system_clock(timer_op_start)

      do i_obs = 1, n_obs
         if (obs(i_obs)%doassim == 1) &
            call obs_op_gridpoint(i_obs, dim_p, dim_obs, state_p, ostate)
      end do

      call system_clock(timer_op_end, t_rate)
      op_dur = op_dur + &
          (real(timer_op_end, wp) - real(timer_op_start, wp))/real(t_rate, wp)
   end subroutine obs_op_pdafomi

   ! set up the localisation information
   subroutine init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)
      ! include functions for different observations
      use mod_observations_pdaf, only: init_dim_obs_l

      implicit none

      ! *** arguments ***
      integer, intent(in)  :: domain_p   !< index of current local analysis domain
      integer, intent(in)  :: step       !< current time step
      integer, intent(in)  :: dim_obs    !< full dimension of observation vector
      integer, intent(out) :: dim_obs_l  !< local dimension of observation vector

      integer :: i_obs
      ! **********************************************
      ! *** initialize local observation dimension ***
      ! **********************************************

      ! call init_dim_obs_l specific for each observation
      do i_obs = 1, n_obs
         if (obs(i_obs)%doassim == 1) &
            call init_dim_obs_l(i_obs, dim_obs_l)
      enddo

   end subroutine init_dim_obs_l_pdafomi

   ! apply covariance localisation
   subroutine localize_covar_pdafomi(dim_p, dim_obs, hp_p, hph)
      use mod_model_pdaf, only: nx, ny, pi, maooam_model
      ! include functions for different observations
      use mod_observations_pdaf, only: localize_covar
      use mod_statevector_pdaf, only: nVar

      implicit none

      ! *** arguments ***
      integer, intent(in)     :: dim_p                 !< pe-local state dimension
      integer, intent(in)     :: dim_obs               !< number of observations
      real(wp), intent(inout) :: hp_p(dim_obs, dim_p)  !< pe local part of matrix hp
      real(wp), intent(inout) :: hph(dim_obs, dim_obs) !< matrix hph

      ! *** local variables ***
      integer   :: i_obs, i, j, k, cnt_p         ! counters
      real(wp)  :: dx, dy
      real(wp)  :: n
      real(wp), allocatable :: coords_p(:,:) ! coordinates of pe-local state vector entries

      ! **********************
      ! *** initialization ***
      ! **********************

      ! *** initialize coordinate array ***

      allocate(coords_p(2, dim_p))

      ! get offset of local domain in global domain in x-direction
      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)
      cnt_p = 0
      do k = 1, nVar
         do j = 1, ny
            do i = 1, nx
               cnt_p = cnt_p + 1
               coords_p(1, cnt_p) = (i-1)*dx
               coords_p(2, cnt_p) = (j-1)*dy
            end do
         end do
      end do

      ! *************************************
      ! *** apply covariance localization ***
      ! *************************************

      ! call localize_covar specific for each observation
      do i_obs = 1, n_obs
         if (obs(i_obs)%doassim == 1) &
            call localize_covar(i_obs, dim_p, dim_obs, hp_p, hph, coords_p)
      enddo

      ! ****************
      ! *** clean up ***
      ! ****************

      deallocate(coords_p)

   end subroutine localize_covar_pdafomi
end module mod_u_pdafomi_pdaf

