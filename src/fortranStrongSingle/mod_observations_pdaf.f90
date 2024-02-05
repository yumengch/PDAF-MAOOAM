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

module mod_observations_pdaf
use mod_kind_pdaf, only: wp
use PDAFomi, only: obs_f
use mod_parallel_pdaf, only: abort_parallel
implicit none

integer   :: n_obs

type :: obs_t
   character(len=50) :: obsname
   character(len=50) :: filename
   character(len=50) :: filename_var
   CHARACTER         :: obsvar

   integer  :: doassim
   integer  :: delt_obs
   integer  :: obs_den = 8
   integer  :: nrows
   integer  :: dim_obs

   real(wp) :: rms_obs
   real(wp) :: missing_value

   real(wp),    allocatable :: obs_field_p(:, :, :)
   real(wp),    allocatable :: var_obs(:, :, :)
   type(obs_f)              :: thisobs

   integer :: file_timecount = 0
   integer :: file_timestep
end type obs_t

type(obs_t),   ALLOCATABLE, target :: obs(:)

contains
   subroutine init()
      use mod_model_pdaf, only: nx, ny
      use mod_filteroptions_pdaf, only: filtertype
      integer                        :: i_obs
      character(len=50), allocatable :: namelist_names(:)
      namelist /n_obs_nml/ n_obs
      namelist /obs_nml/ namelist_names

      ! read number of observation types
      open (20, file='PDAF_config.nml')
      read(20, nml=n_obs_nml)
      rewind(20)

      ! allocate observation type-specific options
      allocate(namelist_names(n_obs))
      allocate(obs(n_obs))

      ! read filename of observation-related namelist
      read(20, nml=obs_nml)
      rewind(20)
      close(20)

      ! initialize observations
      do i_obs = 1, n_obs
         call init_single_obs(i_obs, trim(namelist_names(i_obs)))
      end do
      deallocate(namelist_names)
   end subroutine init

   subroutine init_single_obs(i_obs, nmlname)
      use mod_model_pdaf   , only: nx, ny, pi, maooam_model
      use mod_parallel_pdaf, only: mype_world

      integer     , intent(in) :: i_obs
      character(*), intent(in) :: nmlname

      character(len=50) :: obsname
      character(len=50) :: filename
      character(len=50) :: filename_var

      integer  :: doassim
      integer  :: delt_obs
      integer  :: file_timestep
      integer  :: file_timecount

      real(wp) :: rms_obs
      integer  :: disttype
      integer  :: ncoord
      integer  :: nrows
      integer  :: obs_err_type
      integer  :: use_global_obs
      integer  :: obs_den

      character :: obsvar

      integer  :: nVar, nxo, nyo

      namelist /setup_nml/ obsname, filename, filename_var, &
                           doassim, delt_obs, file_timestep, &
                           rms_obs, disttype, ncoord, &
                           nrows, obs_err_type, &
                           use_global_obs, obs_den, obsvar, &
                           file_timecount

      ! read options for the observation
      open (20, file=nmlname)
      read(20, nml=setup_nml)
      close(20)

      ! assign observation options
      obs(i_obs)%obsvar = obsvar
      obs(i_obs)%obsname = obsname
      obs(i_obs)%filename = filename
      obs(i_obs)%filename_var = filename_var

      obs(i_obs)%delt_obs = delt_obs
      obs(i_obs)%rms_obs = rms_obs

      ! observation types
      obs(i_obs)%nrows = nrows
      obs(i_obs)%doassim = doassim
      obs(i_obs)%thisobs%doassim = doassim
      obs(i_obs)%thisobs%disttype = disttype
      obs(i_obs)%thisobs%ncoord = ncoord
      obs(i_obs)%thisobs%obs_err_type = obs_err_type
      obs(i_obs)%thisobs%use_global_obs = use_global_obs

      obs(i_obs)%obs_den = obs_den
      obs(i_obs)%file_timecount = file_timecount
      obs(i_obs)%file_timestep = file_timestep

      if (mype_world == 0) &
          print *, 'Assimilate observations:', obsname

      ! Size of domain for periodicity for disttype=1
      ! (<0 for no periodicity)
      allocate(obs(i_obs)%thisobs%domainsize(ncoord))
      obs(i_obs)%thisobs%domainsize(1) = 2*pi/maooam_model%model_configuration%physics%n
      obs(i_obs)%thisobs%domainsize(2) = pi
      obs(i_obs)%missing_value = -99999

      ! allocate array for observation and its variance
      nVar = 2
      if (obs(i_obs)%obsvar == 'b') nVar = 4
      nxo = nx/obs(i_obs)%obs_den + 1
      nyo = ny/obs(i_obs)%obs_den + 1
      allocate(obs(i_obs)%obs_field_p(nxo, nyo, nVar))
      allocate(obs(i_obs)%var_obs(nxo, nyo, nVar))
      obs(i_obs)%obs_field_p = 0.
   end subroutine init_single_obs

   subroutine init_dim_obs(i_obs, dim_obs)
      use mod_model_pdaf, only: nx, ny, pi, maooam_model
      use PDAFomi, only: PDAFomi_gather_obs

      integer, intent(in)  :: i_obs
      integer, intent(out) :: dim_obs

      integer  :: i, j, k
      integer  :: cnt
      integer  :: offset
      integer  :: dim_obs_p
      integer  :: nvar, nxo, nyo

      real(wp) :: n
      real(wp) :: dx, dy

      real(wp), allocatable :: obs_p(:)
      real(wp), allocatable :: ivar_obs_p(:)
      real(wp), allocatable :: ocoord_p(:, :)

      nvar = 2

      call get_obs_field(i_obs)
      ! read observation variance
      call get_var_obs(i_obs)

      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)
      nxo = nx/obs(i_obs)%obs_den + 1
      nyo = ny/obs(i_obs)%obs_den + 1

      ! count valid observations
      dim_obs_p = nvar*nxo*nyo
      obs(i_obs)%dim_obs = dim_obs_p

      ! Initialize vector of observations on the process sub-domain
      ! Initialize coordinate array of observations
      ! on the process sub-domain
      print *, 'assimilate ', obs(i_obs)%obsvar, ' component'
      offset = 0
      if (obs(i_obs)%obsvar == 'o') offset = 2*nx*ny
      if (dim_obs_p > 0) then
         allocate(obs_p(dim_obs_p))
         allocate(obs(i_obs)%thisobs%id_obs_p(obs(i_obs)%nrows, dim_obs_p))
         allocate(ocoord_p(obs(i_obs)%thisobs%ncoord, dim_obs_p))
         allocate(ivar_obs_p(dim_obs_p))

         cnt = 1
         ! loop over 4 variables
         do k = 1, nvar
            ! loop over spatial domain
            do j = 1, nyo
               do i = 1, nxo
                  obs_p(cnt) = obs(i_obs)%obs_field_p(i, j, k)
                  obs(i_obs)%thisobs%id_obs_p(1, cnt) = offset + &
                     (k-1)*nx*ny + nx*obs(i_obs)%obs_den*(j-1) + &
                     (i-1)*obs(i_obs)%obs_den + 1
                  ocoord_p(1, cnt) = (i-1)*obs(i_obs)%obs_den*dx
                  ocoord_p(2, cnt) = (j-1)*obs(i_obs)%obs_den*dy
                  ivar_obs_p(cnt) = 1._wp/obs(i_obs)%rms_obs/obs(i_obs)%var_obs(i, j, k)
                  if (obs(i_obs)%var_obs(i, j, k) < 1e-12) ivar_obs_p(cnt) = 1e14
                  cnt = cnt + 1
               end do
            end do
         end do
      else
         allocate(obs_p(1))
         allocate(obs(i_obs)%thisobs%id_obs_p(obs(i_obs)%nrows, 1))
         allocate(ocoord_p(obs(i_obs)%thisobs%ncoord, 1))
         allocate(ivar_obs_p(1))
      end if

      CALL PDAFomi_gather_obs(obs(i_obs)%thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              obs(i_obs)%thisobs%ncoord, 0._wp, dim_obs)

      deallocate(obs_p, ocoord_p, ivar_obs_p)
   end subroutine init_dim_obs

   subroutine get_obs_field(i_obs)
      use mod_nfcheck_pdaf, only: check
      use mod_model_pdaf, only: nx, ny
      use netcdf
      integer,      intent(in)    :: i_obs

      integer :: ncid, varid
      integer :: i, cnt
      integer :: i_start, i_end
      integer :: nxo, nyo
      character(len=5) :: varname(4)

      obs(i_obs)%file_timecount = obs(i_obs)%file_timecount + obs(i_obs)%file_timestep
      varname = [character(len=5) :: 'psi_a', 'T_a', 'psi_o', 'T_o']

      call check( nf90_open(trim(obs(i_obs)%filename), nf90_nowrite, ncid) )
      ! choose observed variables
      i_start = 1
      i_end = 4
      if (obs(i_obs)%obsvar == 'a') then
         i_end = 2
      else if (obs(i_obs)%obsvar == 'o') then
         i_start = 3
      else if (obs(i_obs)%obsvar == 'b') then
      else
         print *, '...obsvar is not correct...', i_obs, obs(i_obs)%obsvar
         call abort_parallel()
      endif
      nxo = nx/obs(i_obs)%obs_den + 1
      nyo = ny/obs(i_obs)%obs_den + 1

      ! read observations
      cnt = 1
      do i = i_start, i_end
         call check( nf90_inq_varid(ncid, trim(varname(i)), varid) )
         call check( nf90_get_var(ncid, varid, &
                                  obs(i_obs)%obs_field_p(:, :, cnt), &
                                  start=[1, 1, obs(i_obs)%file_timecount], &
                                  count=[nxo, nyo, 1] &
                                  ) &
                    )
         cnt = cnt + 1
      end do
      call check( nf90_close(ncid) )
   end subroutine get_obs_field

   subroutine get_var_obs(i_obs)
      use mod_model_pdaf, only: nx, ny
      use mod_nfcheck_pdaf, only: check
      use netcdf
      integer,      intent(in)    :: i_obs

      integer :: ncid, varid
      integer :: i, cnt
      integer :: i_start, i_end
      integer :: nxo, nyo
      character(len=10) :: varname(4)

      varname = [character(len=10) :: 'psi_a_var', 'T_a_var', 'psi_o_var', 'T_o_var']

      ! open NC file
      call check( nf90_open(trim(obs(i_obs)%filename_var), nf90_nowrite, ncid) )

      ! choose observed variables
      i_start = 1
      i_end = 4
      if (obs(i_obs)%obsvar == 'a') then
         i_end = 2
      else if (obs(i_obs)%obsvar == 'o') then
         i_start = 3
      else if (obs(i_obs)%obsvar == 'b') then
      else
         print *, '...obsvar is not correct...', obs(i_obs)%obsvar, i_obs
         call abort_parallel()
      endif
      nxo = nx/obs(i_obs)%obs_den + 1
      nyo = ny/obs(i_obs)%obs_den + 1

      ! read observations
      cnt = 1
      do i = i_start, i_end
         call check( nf90_inq_varid(ncid, trim(varname(i)), varid) )
         call check( nf90_get_var(ncid, varid, &
                             obs(i_obs)%var_obs(:, :, cnt), &
                             start=[1, 1], count=[nxo, nyo]) )
         cnt = cnt + 1
      end do
      call check( nf90_close(ncid) )
   end subroutine get_var_obs

   SUBROUTINE obs_op_gridpoint(i_obs, dim_p, dim_obs, state_p, ostate)
      USE PDAFomi, ONLY: PDAFomi_obs_op_gridpoint
      IMPLICIT NONE

      ! *** Arguments ***
      integer,  intent(in)    :: i_obs
      INTEGER,  INTENT(in)    :: dim_p                 !< PE-local state dimension
      INTEGER,  INTENT(in)    :: dim_obs               !< Dimension of full observed state (all observed fields)
      REAL(wp), INTENT(in)    :: state_p(dim_p)        !< PE-local model state
      REAL(wp), INTENT(inout) :: ostate(dim_obs)       !< Full observed state

      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************
      ! observation operator for observed grid point values
      CALL PDAFomi_obs_op_gridpoint(obs(i_obs)%thisobs, state_p, ostate)
   END SUBROUTINE obs_op_gridpoint

   subroutine finalize_obs()
      integer :: i_obs

      do i_obs = 1, n_obs
         deallocate(obs(i_obs)%obs_field_p)
         deallocate(obs(i_obs)%var_obs)
      end do
      deallocate(obs)
   end subroutine finalize_obs
end module mod_observations_pdaf
