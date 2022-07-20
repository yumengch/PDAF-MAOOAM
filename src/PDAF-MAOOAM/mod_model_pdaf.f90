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

module mod_model_pdaf
use mod_kind_pdaf, only: wp
use model_def, only: model
use rk4_integrator, only: RK4Integrator
use mod_ModelWriter_pdaf, only: init_model_writer, finalize_model_writer
implicit none

! Dimension of state vector and ensemble size
integer :: dim_ens
integer :: dim_state
integer :: dim_state_p
integer :: noc
integer :: natm

integer :: total_steps
real(wp) :: total_time
real(wp), allocatable :: field(:)
real(wp), allocatable :: field_new(:)

type(Model), TARGET :: maooam_model
type(RK4Integrator) :: integr
contains
   subroutine initialize_model()
      print *, 'Model MAOOAM v1.4'
      print *, 'Loading information...'
      ! initialise model configurations 
      CALL maooam_model%init

      natm = maooam_model%model_configuration%modes%natm
      noc = maooam_model%model_configuration%modes%noc
      dim_state_p = maooam_model%model_configuration%modes%ndim

      total_time = maooam_model%model_configuration%integration%t_run
      dim_state = dim_state_p
      total_steps = int(total_time/maooam_model%model_configuration%integration%dt)

      CALL integr%init(maooam_model)

      ! initialise the model writer
      call init_model_writer('maooam.nc', natm, noc, dim_ens)

      ! initialise initial condition
      ALLOCATE(field(0:dim_state_p),field_new(0:dim_state_p))
      field = maooam_model%load_IC()
   end subroutine initialize_model

   subroutine finalize_model
      DEALLOCATE(field, field_new)
      CALL maooam_model%clean
      CALL integr%clean
      call finalize_model_writer()
   end subroutine finalize_model
end module mod_model_pdaf