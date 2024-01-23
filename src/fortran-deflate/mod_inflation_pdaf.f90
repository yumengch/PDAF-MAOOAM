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



module mod_inflation_pdaf
use mod_kind_pdaf
implicit none
! Inflation options

! Attributes
! ----------
! forget : float
!     forgetting factor
! type_forget : int
!     Type of forget factor
!     type of forgetting factor
!     - (0) fixed
!     - (1) global adaptive
!     - (2) local adaptive for LSEIK/LETKF/LESTKF
integer :: type_forget = 0
! forgeting factor
real(wp) :: forget = 1
real(wp) :: alpha = 0.
namelist /infl_nml/ type_forget, forget, alpha
end module mod_inflation_pdaf