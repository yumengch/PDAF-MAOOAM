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
use PDAFomi, only: obs_f, obs_l
use mod_parallel_pdaf, only: abort_parallel
implicit none

integer   :: n_obs
integer   :: nxo, nyo
integer   :: obs_den = 8

CHARACTER, dimension(:), ALLOCATABLE :: obsvar_all(:)

integer ,            allocatable :: delt_obs_all(:)
integer ,            allocatable :: nrows_all(:)
real(wp),            allocatable :: missing_value(:)
real(wp),            allocatable :: rms_obs_all(:)
real(wp),            allocatable :: obs_field_p_all(:, :, :, :)
real(wp),            allocatable :: var_obs(:, :, :, :)
type(obs_f), target, allocatable :: thisobs(:)
type(obs_l), target, allocatable :: thisobs_l(:)

integer,             ALLOCATABLE :: loc_weight_all(:)
real(wp),            ALLOCATABLE :: local_range_all(:)
real(wp),            ALLOCATABLE :: srange_all(:)


integer,             allocatable :: time_count(:)
integer,             allocatable :: time_interval(:)
character(len=50),   allocatable :: obsnames(:)
character(len=50),   allocatable :: filenames(:)

contains
   subroutine init()
      use mod_model_pdaf, only: nx, ny
      use mod_filteroptions_pdaf, only: filtertype
      integer                        :: i_obs
      integer                        :: nvar
      character(len=50), allocatable :: namelist_names(:)
      namelist /n_obs_nml/ n_obs
      namelist /obs_nml/ namelist_names

      ! read number of observation types
      open (20, file='PDAF_config.nml')
      read(20, nml=n_obs_nml)
      rewind(20)

      ! allocate observation type-specific options
      allocate(namelist_names(n_obs))
      allocate(delt_obs_all(n_obs))
      allocate(nrows_all(n_obs), missing_value(n_obs))
      allocate(rms_obs_all(n_obs), thisobs(n_obs), thisobs_l(n_obs))
      allocate(filenames(n_obs), obsnames(n_obs))
      allocate(loc_weight_all(n_obs), local_range_all(n_obs), srange_all(n_obs))
      allocate(obsvar_all(n_obs))
      allocate(time_count(n_obs))
      allocate(time_interval(n_obs))
      time_count = 0
      time_interval = 1

      ! read filename of observation-related namelist
      read(20, nml=obs_nml)
      rewind(20)
      close(20)

      ! get the size of the observation vector
      nxo = nx/obs_den + 1
      nyo = ny/obs_den + 1

      nvar = 2
      if (filtertype == 100) nvar = 4

      ! observation vector and error variance
      allocate(obs_field_p_all(nxo, nyo, nvar, n_obs))
      allocate(var_obs(nxo, nyo, nvar, n_obs))
      ! initialize observations
      do i_obs = 1, n_obs
         call init_single_obs(i_obs, trim(namelist_names(i_obs)))
      end do
      deallocate(namelist_names)
   end subroutine init

   subroutine init_single_obs(i_obs, nmlname)
      use mod_model_pdaf   , only: pi, maooam_model
      use mod_parallel_pdaf, only: mype_filter

      integer     , intent(in) :: i_obs
      character(*), intent(in) :: nmlname

      character(len=50) :: obsname
      character(len=50) :: filename
      character(len=50) :: filename_var

      integer  :: doassim
      integer  :: delt_obs
      integer  :: disttype
      integer  :: ncoord
      integer  :: nrows
      integer  :: obs_err_type
      integer  :: use_global_obs
      real(wp) :: rms_obs
      
      integer  :: loc_weight
      real(wp) :: local_range
      real(wp) :: srange
      character :: obsvar
      integer   :: interval

      namelist /setup_nml/ obsname, filename, filename_var, &
                           doassim, delt_obs, &
                           rms_obs, disttype, ncoord, &
                           nrows, obs_err_type, &
                           use_global_obs, &
                           loc_weight, local_range, srange, &
                           obsvar, interval 

      ! read options for the observation
      open (20, file=nmlname)
      read(20, nml=setup_nml)
      rewind(20)
      close(20)

      ! assign observation options
      thisobs(i_obs)%doassim = doassim
      delt_obs_all(i_obs) = delt_obs
      nrows_all(i_obs) = nrows
      rms_obs_all(i_obs) = rms_obs
      filenames(i_obs) = filename
      obsnames(i_obs) = obsname

      loc_weight_all(i_obs) = loc_weight
      local_range_all(i_obs) = local_range
      srange_all(i_obs) = srange

      obsvar_all(i_obs) = obsvar

      thisobs(i_obs)%disttype = disttype
      thisobs(i_obs)%ncoord = ncoord
      thisobs(i_obs)%obs_err_type = obs_err_type
      thisobs(i_obs)%use_global_obs = use_global_obs

      time_interval(i_obs) = interval
      if (mype_filter == 0) &
          print *, 'Assimilate observations:', obsnames(i_obs)

      ! read observation variance
      call get_var_obs(i_obs, trim(filename_var), var_obs(:, :, :, i_obs))

      ! Size of domain for periodicity for disttype=1
      ! (<0 for no periodicity)
      allocate(thisobs(i_obs)%domainsize(ncoord))
      thisobs(i_obs)%domainsize(1) = 2*pi/maooam_model%model_configuration%physics%n
      thisobs(i_obs)%domainsize(2) = pi
      missing_value(i_obs) = -99999
   end subroutine init_single_obs

   subroutine init_dim_obs(i_obs, step, dim_obs)
      use mod_model_pdaf, only: nx, ny, pi, maooam_model
      use PDAFomi, only: PDAFomi_gather_obs

      integer, intent(in)  :: i_obs
      integer, intent(in)  :: step
      integer, intent(out) :: dim_obs

      integer  :: i, j, k
      integer  :: cnt
      integer  :: offset
      integer  :: dim_obs_p
      integer  :: nvar

      real(wp) :: n
      real(wp) :: dx, dy

      real(wp), allocatable :: obs_p(:)
      real(wp), allocatable :: ivar_obs_p(:)
      real(wp), allocatable :: ocoord_p(:, :)

      nvar = 2
      if (mod ( step, delt_obs_all(i_obs) ) /= 0) then
         print *, 'no assim', i_obs, obsvar_all(i_obs), step
         dim_obs = 0
         dim_obs_p = 0

         allocate(obs_p(1))
         allocate(ivar_obs_p(1))
         allocate(ocoord_p(thisobs(i_obs)%ncoord, 1))
         allocate(thisobs(i_obs)%id_obs_p(1, dim_obs_p))

         call PDAFomi_gather_obs(thisobs(i_obs), dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
            thisobs(i_obs)%ncoord, local_range_all(i_obs), dim_obs)

         deallocate(obs_p, ocoord_p, ivar_obs_p)

         return
      end if
      call get_obs_field(i_obs, trim(filenames(i_obs)))

      offset = 0
      if (obsvar_all(i_obs) == 'o') offset = 2*nx*ny

      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)

      ! count valid observations
      cnt = 0
      do k = 1, nvar
         ! loop over spatial domain
         do j = 1, nyo
            do i = 1, nxo
               if (obs_field_p_all(i, j, k, i_obs) > missing_value(i_obs)) cnt = cnt + 1
            enddo
         enddo
      end do
      dim_obs_p = cnt

      ! Initialize vector of observations on the process sub-domain
      ! Initialize coordinate array of observations
      ! on the process sub-domain
      if (dim_obs_p > 0) then
         allocate(obs_p(dim_obs_p))
         allocate(thisobs(i_obs)%id_obs_p(nrows_all(i_obs), dim_obs_p))
         allocate(ocoord_p(thisobs(i_obs)%ncoord, dim_obs_p))
         allocate(ivar_obs_p(dim_obs_p))

         cnt = 1
         ! loop over 4 variables
         do k = 1, nvar
            ! loop over spatial domain
            do j = 1, nyo
               do i = 1, nxo
                  if (obs_field_p_all(i, j, k, i_obs) > missing_value(i_obs)) then
                     obs_p(cnt) = obs_field_p_all(i, j, k, i_obs)
                     thisobs(i_obs)%id_obs_p(1, cnt) = offset + (k-1)*nx*ny + nx*obs_den*(j-1) + (i-1)*obs_den + 1
                     ocoord_p(1, cnt) = (i-1)*obs_den*dx
                     ocoord_p(2, cnt) = (j-1)*obs_den*dy
                     ivar_obs_p(cnt) = 1._wp/rms_obs_all(i_obs)/var_obs(i, j, k, i_obs)
                     if (var_obs(i, j, k, i_obs) < 1e-12) ivar_obs_p(cnt) = 1e14
                     cnt = cnt + 1
                  end if
               end do
            end do
         end do
      else
         allocate(obs_p(1))
         allocate(thisobs(i_obs)%id_obs_p(nrows_all(i_obs), 1))
         allocate(ocoord_p(thisobs(i_obs)%ncoord, 1))
         allocate(ivar_obs_p(1))
         thisobs(i_obs)%id_obs_p = 0._wp
      end if

      CALL PDAFomi_gather_obs(thisobs(i_obs), dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs(i_obs)%ncoord, local_range_all(i_obs), dim_obs)

      deallocate(obs_p, ocoord_p, ivar_obs_p)
   end subroutine init_dim_obs

   subroutine init_dim_obs_gen(i_obs, dim_obs)
      use mod_model_pdaf, only: nx, ny, pi, maooam_model
      use PDAFomi, only: PDAFomi_gather_obs
      integer, intent(in)  :: i_obs
      integer, intent(out) :: dim_obs

      integer  :: i, j, k
      integer  :: cnt
      integer  :: offset
      integer  :: dim_obs_p
      integer  :: nvar

      real(wp) :: n
      real(wp) :: dx, dy

      real(wp), allocatable :: obs_p(:)
      real(wp), allocatable :: ivar_obs_p(:)
      real(wp), allocatable :: ocoord_p(:, :)


      nvar = 4

      offset = 0
      if (obsvar_all(i_obs) == 'o') offset = 2*nx*ny

      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)

      ! count valid observations
      dim_obs_p = nvar*nxo*nyo

      ! Initialize vector of observations on the process sub-domain
      ! Initialize coordinate array of observations
      ! on the process sub-domain
      if (dim_obs_p > 0) then
         allocate(obs_p(dim_obs_p))
         allocate(thisobs(i_obs)%id_obs_p(nrows_all(i_obs), dim_obs_p))
         allocate(ocoord_p(thisobs(i_obs)%ncoord, dim_obs_p))
         allocate(ivar_obs_p(dim_obs_p))
         cnt = 1
         ! loop over 4 variables
         do k = 1, nvar
            ! loop over spatial domain
            do j = 1, nyo
               do i = 1, nxo
                  obs_p(cnt) = obs_field_p_all(i, j, k, i_obs)
                  thisobs(i_obs)%id_obs_p(1, cnt) = offset + (k-1)*nx*ny + nx*obs_den*(j-1) + (i-1)*obs_den + 1
                  ocoord_p(1, cnt) = (i-1)*obs_den*dx
                  ocoord_p(2, cnt) = (j-1)*obs_den*dy
                  ivar_obs_p(cnt) = 1._wp/rms_obs_all(i_obs)/var_obs(i, j, k, i_obs)
                  if (var_obs(i, j, k, i_obs) < 1e-12) ivar_obs_p(i) = 1e14
                  cnt = cnt + 1
               end do
            end do
         end do
      else
         allocate(obs_p(1))
         allocate(thisobs(i_obs)%id_obs_p(nrows_all(i_obs), 1))
         allocate(ocoord_p(thisobs(i_obs)%ncoord, 1))
         allocate(ivar_obs_p(1))
         thisobs(i_obs)%id_obs_p = 0._wp
      end if

      CALL PDAFomi_gather_obs(thisobs(i_obs), dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs(i_obs)%ncoord, local_range_all(i_obs), dim_obs)

      deallocate(obs_p, ocoord_p, ivar_obs_p)
   end subroutine init_dim_obs_gen

   subroutine get_obs_field(i_obs, filename)
      use mod_nfcheck_pdaf, only: check
      use netcdf
      integer,      intent(in)    :: i_obs
      character(*), intent(in)    :: filename

      integer :: ncid, varid
      integer :: i, cnt
      integer :: i_start, i_end
      character(len=5) :: varname(4)

      time_count(i_obs) = time_count(i_obs) + time_interval(i_obs)
      varname = [character(len=5) :: 'psi_a', 'T_a', 'psi_o', 'T_o']

      call check( nf90_open(filename, nf90_nowrite, ncid) )
      ! choose observed variables
      i_start = 1
      i_end = 4
      if (obsvar_all(i_obs) == 'a') then
         i_end = 2
      else if (obsvar_all(i_obs) == 'o') then
         i_start = 3
      else if (obsvar_all(i_obs) == 'b') then
      else
         print *, '...obsvar is not correct...', i_obs, obsvar_all(i_obs)
         call abort_parallel()
      endif

      ! read observations
      cnt = 1
      do i = i_start, i_end
         call check( nf90_inq_varid(ncid, trim(varname(i)), varid) )
         call check( nf90_get_var(ncid, varid, &
                             obs_field_p_all(:, :, cnt, i_obs), &
                             start=[1, 1, time_count(i_obs)], count=[nxo, nyo, 1]) )
         cnt = cnt + 1
      end do
      call check( nf90_close(ncid) )
   end subroutine get_obs_field

   subroutine get_var_obs(i_obs, filename, var_obs_i)
      use mod_nfcheck_pdaf, only: check
      use netcdf
      integer,      intent(in)    :: i_obs
      character(*), intent(in)    :: filename
      real(wp),     intent(inout) :: var_obs_i(:, :, :)

      integer :: ncid, varid
      integer :: i, cnt
      integer :: i_start, i_end
      character(len=10) :: varname(4)

      varname = [character(len=10) :: 'psi_a_var', 'T_a_var', 'psi_o_var', 'T_o_var']

      ! open NC file
      call check( nf90_open(filename, nf90_nowrite, ncid) )

      ! choose observed variables
      i_start = 1
      i_end = 4
      if (obsvar_all(i_obs) == 'a') then
         i_end = 2
      else if (obsvar_all(i_obs) == 'o') then
         i_start = 3
      else if (obsvar_all(i_obs) == 'b') then
      else
         print *, '...obsvar is not correct...', obsvar_all(i_obs), i_obs
         call abort_parallel()
      endif

      ! read observations
      cnt = 1
      do i = i_start, i_end
         call check( nf90_inq_varid(ncid, trim(varname(i)), varid) )
         call check( nf90_get_var(ncid, varid, &
                             var_obs_i(:, :, cnt), &
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

      IF (thisobs(i_obs)%doassim==1) THEN
         ! observation operator for observed grid point values
         CALL PDAFomi_obs_op_gridpoint(thisobs(i_obs), state_p, ostate)
      END IF
   END SUBROUTINE obs_op_gridpoint

   SUBROUTINE init_dim_obs_l(i_obs, dim_obs_l)
      ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l
      ! Include localization radius and local coordinates
      USE mod_localization_pdaf, &   
            ONLY: coords_l

      IMPLICIT NONE

      ! *** Arguments ***
      integer, intent(in)    :: i_obs
      INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector

      ! **********************************************
      ! *** Initialize local observation dimension ***
      ! **********************************************

      CALL PDAFomi_init_dim_obs_l(thisobs_l(i_obs), thisobs(i_obs), coords_l, &
                                  loc_weight_all(i_obs), local_range_all(i_obs), &
                                  srange_all(i_obs), dim_obs_l)
   END SUBROUTINE init_dim_obs_l

   SUBROUTINE localize_covar(i_obs, dim_p, dim_obs, HP_p, HPH, coords_p)
      ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_localize_covar

      IMPLICIT NONE

      ! *** Arguments ***
      integer, intent(in)  :: i_obs
      INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
      REAL(wp), INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
      REAL(wp), INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
      REAL(wp), INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


      ! *************************************
      ! *** Apply covariance localization ***
      ! *************************************

      CALL PDAFomi_localize_covar(thisobs(i_obs), dim_p, loc_weight_all(i_obs), &
                                  local_range_all(i_obs), srange_all(i_obs), &
                                  coords_p, HP_p, HPH)
   END SUBROUTINE localize_covar

   subroutine finalize_obs()
      deallocate(delt_obs_all)
      deallocate(nrows_all, missing_value)
      deallocate(rms_obs_all, thisobs, thisobs_l)
      deallocate(filenames, obsnames)

      deallocate(obs_field_p_all)
      deallocate(var_obs)
      deallocate(time_count)
   end subroutine finalize_obs
end module mod_observations_pdaf
