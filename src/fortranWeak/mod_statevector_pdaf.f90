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

contains
   subroutine initSV()
      use mod_model_pdaf, only: nx, ny
      nVar = 2
      dim_state_p = nx*ny*nVar
      dim_state = nx*ny*nVar
      dim_ens = n_modeltasks
   end subroutine initSV

   subroutine setField(vartype)
      CHARACTER, intent(in) :: vartype
      if (vartype == 'a') then
         sv_atm = .true.
         sv_ocean = .false.
      else if (vartype == 'o') then
         sv_atm = .false.
         sv_ocean = .true.
      end if
   end subroutine setField

end module mod_statevector_pdaf
