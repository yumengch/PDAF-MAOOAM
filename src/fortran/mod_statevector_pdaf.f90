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

module mod_statevector_pdaf
use mod_parallel_pdaf, only: n_modeltasks
implicit none

integer :: dim_state_p
integer :: dim_state
integer :: dim_ens
integer :: nVar 
logical :: sv_atm = .false.
logical :: sv_ocean = .false.
CHARACTER :: component
namelist /state_vector_nml/ component
contains
   subroutine initSV(is_strong)
      use mod_model_pdaf, only: nx, ny
      logical, intent(in) :: is_strong
      if (component == 'a') then
         sv_atm = .true.
      else if (component == 'o') then
         sv_ocean = .true.
      else if (component == 'b') then
         sv_atm = .true.
         sv_ocean = .true.
      end if
      nVar = 2
      if (is_strong) nVar = 4
      dim_state_p = nx*ny*nVar
      dim_state = nx*ny*nVar
      dim_ens = n_modeltasks
   end subroutine initSV

   subroutine setField(vartype)
      CHARACTER, intent(in) :: vartype
      if (vartype == 'a') then
         sv_atm = .true.
         sv_ocean = .false.
      else if (component == 'o') then
         sv_atm = .false.
         sv_ocean = .true.
      else if (component == 'b') then
         sv_atm = .true.
         sv_ocean = .true.
      end if
   end subroutine setField

end module mod_statevector_pdaf