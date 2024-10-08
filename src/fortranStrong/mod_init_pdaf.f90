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
module mod_init_pdaf
use mod_kind_pdaf, only: wp
use mod_parallel_pdaf, &
   only: comm_model, comm_couple, comm_filter, &
         task_id, filterpe, dim_ens => n_modeltasks, mype_world, &
         abort_parallel
use mod_filteroptions_pdaf, only: filtertype, subtype, &
                                  type_sqrt, type_trans, &
                                  incremental
use mod_inflation_pdaf, only: type_forget, forget
implicit none

contains
   subroutine setETKFOptions(filter_param_i, filter_param_r)
      use mod_statevector_pdaf, only: dim_state_p
      integer,   intent(inout) :: filter_param_i(:)
      real(wp),  intent(inout) :: filter_param_r(:)

      filter_param_i(1) = dim_state_p
      filter_param_i(2) = dim_ens
      filter_param_i(3) = 0
      filter_param_i(4) = incremental
      filter_param_i(5) = type_forget
      filter_param_i(6) = type_trans
      filter_param_i(7) = type_sqrt

      filter_param_r(1) = forget
   end subroutine setETKFOptions

   subroutine init_pdaf(screen)
      use pdaf_interfaces_module, only: PDAF_init, PDAF_get_state
      use mod_localisation_pdaf, only: set_lim_coords
      use mod_U_pdaf, only: init_ens_pdaf, &
                            distribute_state_pdaf, &
                            next_observation_pdaf,prepoststep_ens_pdaf
      integer, intent(in) :: screen

      integer           :: filter_param_i(7) ! Integer parameter array for filter
      real(wp)          :: filter_param_r(2) ! Real parameter array for filter
      integer           :: n_param_i
      integer           :: status_pdaf
      integer           :: steps = 0, doexit
      real(wp)          :: timenow = 0.

      ! All other filters
      n_param_i = 7
      call setETKFOptions(filter_param_i, filter_param_r)

      call PDAF_init(filtertype, subtype, 0, filter_param_i, n_param_i, &
                     filter_param_r, 2, comm_model, comm_filter, &
                     comm_couple, task_id, dim_ens, filterpe, &
                     init_ens_pdaf, screen, status_pdaf)

      if (status_pdaf /= 0) then
         write (*,'(/1x,a6,i3,a43,i4,a1/)') &
             'error ', status_pdaf, &
             ' in initialization of pdaf - stopping! (pe ', mype_world,')'
         call abort_parallel()
      end if

      call PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, distribute_state_pdaf, &
                          prepoststep_ens_pdaf, status_pdaf)

      ! Set the domain limits for the localisation in OMI
      call set_lim_coords()
   end subroutine init_pdaf

   subroutine finalize_pdaf()
      use mod_parallel_pdaf, &
         only: mype_world, mpierr
      use pdaf_interfaces_module, only: PDAF_print_info, PDAF_deallocate
      use mod_StateWriter_pdaf, only: finalize_state_writer
      use mod_obswriter_pdaf, only: finalizeObs
      use mod_observations_pdaf, only: finalize_obs
      use mpi

      if (filterpe) call finalize_state_writer()
      ! *** Show allocated memory for PDAF ***
      call PDAF_print_info(11)

      ! *** Print PDAF timings onto screen ***
      if (mype_world==0) call PDAF_print_info(3)

      if (filtertype == 100)  call finalizeObs()

      call finalize_obs()

      ! Deallocate PDAF arrays
      call PDAF_deallocate()
      call mpi_barrier(mpi_comm_world, mpierr)
   end subroutine finalize_pdaf
end module mod_init_pdaf
