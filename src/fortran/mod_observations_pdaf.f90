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
integer :: n_obs
integer :: obssize
integer, allocatable :: delt_obs_all(:)
integer, allocatable :: nrows_all(:)
real(wp), allocatable :: missing_value(:)
real(wp), allocatable :: rms_obs_all(:)
real(wp), allocatable :: obs_field_p_all(:, :)
real(wp), allocatable :: var_obs(:, :)
type(obs_f), target, allocatable :: thisobs(:)
type(obs_l), target, allocatable :: thisobs_l(:)

character :: obsvar
character(len=50), allocatable :: obsnames(:)
character(len=50), allocatable :: filenames(:)

contains
   subroutine init()
      use mod_model_pdaf, only: dim_state_p, dim_state, &
                                 nx, ny
      integer :: i_obs
      character(len=50), allocatable :: namelist_names(:)
      namelist /n_obs_nml/ n_obs
      namelist /obs_nml/ namelist_names, obsvar

      open (20, file='PDAF_config.nml')
      read(20, nml=n_obs_nml)
      rewind(20)

      allocate(namelist_names(n_obs))
      allocate(delt_obs_all(n_obs))
      allocate(nrows_all(n_obs), missing_value(n_obs))
      allocate(rms_obs_all(n_obs), thisobs(n_obs), thisobs_l(n_obs))
      allocate(filenames(n_obs), obsnames(n_obs))

      read(20, nml=obs_nml)
      rewind(20)
      close(20)

      if (obsvar == 'b') then
         obssize = 4*nx*ny
      else if ((obsvar == 'a') .or. (obsvar == 'o')) then
         obssize = 2*nx*ny
      else
         print *, 'Not implemented for', obsvar
         call abort_parallel()
      end if
      allocate(obs_field_p_all(obssize, n_obs))
      allocate(var_obs(obssize, n_obs))

      do i_obs = 1, n_obs
         call init_single_obs(i_obs, trim(namelist_names(i_obs)))
      end do
      deallocate(namelist_names)
   end subroutine init

   subroutine init_single_obs(i_obs, nmlname)
      use mod_model_pdaf, only: pi, maooam_model
      use mod_parallel_pdaf, only: mype_filter
      integer, intent(in) :: i_obs
      character(*), intent(in) :: nmlname
      character(len=50) :: obsname
      character(len=50) :: filename
      integer :: doassim
      integer :: delt_obs
      real(wp) :: rms_obs
      integer :: disttype, ncoord, nrows
      integer :: obs_err_type, use_global_obs
      
      namelist /setup_nml/ obsname, filename, &
                           doassim, delt_obs, &
                           rms_obs, disttype, ncoord, &
                           nrows, obs_err_type, &
                           use_global_obs

      open (20, file=nmlname)
      read(20, nml=setup_nml)
      rewind(20)
      close(20)
      thisobs(i_obs)%doassim = doassim
      delt_obs_all(i_obs) = delt_obs
      nrows_all(i_obs) = nrows
      rms_obs_all(i_obs) = rms_obs
      filenames(i_obs) = filename
      obsnames(i_obs) = obsname
      thisobs(i_obs)%disttype = disttype
      thisobs(i_obs)%ncoord = ncoord
      thisobs(i_obs)%obs_err_type = obs_err_type
      thisobs(i_obs)%use_global_obs = use_global_obs
      if (mype_filter == 0) &
          print *, 'Assimilate observations:', obsnames(i_obs)

      call get_var_obs(i_obs, trim(filenames(i_obs)), var_obs(:, i_obs))
      ! Size of domain for periodicity for disttype=1
      ! (<0 for no periodicity)
      allocate(thisobs(i_obs)%domainsize(ncoord))
      thisobs(i_obs)%domainsize(1) = 2*pi/maooam_model%model_configuration%physics%n
      thisobs(i_obs)%domainsize(2) = pi
      missing_value(i_obs) = -99999
   end subroutine init_single_obs

   subroutine init_dim_obs(i_obs, dim_obs)
      use mod_localization_pdaf, only: local_range
      use mod_parallel_pdaf, only: mype_filter
      use mod_model_pdaf, only: dim_state_p, dim_state
      use PDAFomi, only: PDAFomi_gather_obs
      integer, intent(in) :: i_obs
      integer, intent(out) :: dim_obs

      integer  :: pe_start, pe_end
      integer  :: i
      integer  :: cnt_p
      integer  :: dim_obs_p
      
      ! real(wp) :: obs_field(dim_state)
      real(wp), allocatable :: obs_p(:)
      real(wp), allocatable :: ivar_obs_p(:)
      real(wp), allocatable :: var_obs_p(:)
      real(wp), allocatable :: ocoord_p(:, :)

      ! call get_obs_field(i_obs, trim(filenames(i_obs)), obs_field)

      ! ! Count valid observations that
      ! ! lie within the process sub-domain
      ! pe_start = dim_state_p*mype_filter + 1
      ! pe_end = dim_state_p*(mype_filter+1)
      ! obs_field_p = obs_field(pe_start:pe_end)

      ! count valid observations
      cnt_p = 0
      do i = 1, dim_state_p
         if (obs_field_p_all(i, i_obs) > missing_value(i_obs)) cnt_p = cnt_p + 1
      end do
      dim_obs_p = cnt_p

      ! Initialize vector of observations on the process sub-domain
      ! Initialize coordinate array of observations
      ! on the process sub-domain
      if (dim_obs_p > 0) then
         allocate(obs_p(dim_obs_p), var_obs_p(dim_obs_p))
         call set_obs_p(i_obs, obs_field_p_all(:, i_obs), var_obs(:, i_obs), obs_p, var_obs_p)
         allocate(thisobs(i_obs)%id_obs_p(nrows_all(i_obs), dim_obs_p))
         call set_id_obs_p(i_obs, obs_field_p_all(:, i_obs), thisobs(i_obs)%id_obs_p)
         allocate(ocoord_p(thisobs(i_obs)%ncoord, dim_obs_p))
         call set_ocoord_p(i_obs, pe_start, obs_field_p_all(:, i_obs), ocoord_p)
         allocate(ivar_obs_p(dim_obs_p))
         ivar_obs_p = 1._wp/rms_obs_all(i_obs)/var_obs_p
      else
         allocate(obs_p(1))
         allocate(thisobs(i_obs)%id_obs_p(nrows_all(i_obs), 1))
         allocate(ocoord_p(thisobs(i_obs)%ncoord, 1))
         allocate(ivar_obs_p(1))
         thisobs(i_obs)%id_obs_p = 0._wp
      end if

      CALL PDAFomi_gather_obs(thisobs(i_obs), dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs(i_obs)%ncoord, local_range, dim_obs)

      deallocate(obs_p, ocoord_p, ivar_obs_p)
   end subroutine init_dim_obs

   subroutine init_dim_obs_gen(i_obs, dim_obs)
      use mod_localization_pdaf, only: local_range
      use mod_parallel_pdaf, only: mype_filter
      use mod_model_pdaf, only: dim_state_p, dim_state
      use PDAFomi, only: PDAFomi_gather_obs
      integer, intent(in)  :: i_obs
      integer, intent(out) :: dim_obs

      integer  :: pe_start, pe_end
      integer  :: i
      integer  :: cnt_p
      integer  :: dim_obs_p
      ! real(wp) :: obs_field_p(dim_state_p)
      ! real(wp) :: obs_field(dim_state)
      real(wp), allocatable :: obs_p(:)
      real(wp), allocatable :: ivar_obs_p(:)
      real(wp), allocatable :: var_obs_p(:)
      real(wp), allocatable :: ocoord_p(:, :)

      ! obs_field = 0._wp

      ! Count valid observations that
      ! lie within the process sub-domain
      ! pe_start = dim_state_p*mype_filter + 1
      ! pe_end = dim_state_p*(mype_filter+1)
      ! obs_field_p = obs_field(pe_start:pe_end)

      ! count valid observations
      cnt_p = 0
      do i = 1, dim_state_p
         if (obs_field_p_all(i, i_obs) > missing_value(i_obs)) cnt_p = cnt_p + 1
      end do
      dim_obs_p = cnt_p

      ! Initialize vector of observations on the process sub-domain
      ! Initialize coordinate array of observations
      ! on the process sub-domain
      if (dim_obs_p > 0) then
         allocate(obs_p(dim_obs_p), var_obs_p(dim_obs_p))
         call set_obs_p(i_obs, obs_field_p_all(:, i_obs), var_obs(:, i_obs), obs_p, var_obs_p)
         allocate(thisobs(i_obs)%id_obs_p(nrows_all(i_obs), dim_obs_p))
         call set_id_obs_p(i_obs, obs_field_p_all(:, i_obs), thisobs(i_obs)%id_obs_p)
         allocate(ocoord_p(thisobs(i_obs)%ncoord, dim_obs_p))
         call set_ocoord_p(i_obs, pe_start, obs_field_p_all(:, i_obs), ocoord_p)
         allocate(ivar_obs_p(dim_obs_p))
         ivar_obs_p = 1._wp/rms_obs_all(i_obs)/rms_obs_all(i_obs)
      else
         allocate(obs_p(1))
         allocate(thisobs(i_obs)%id_obs_p(nrows_all(i_obs), 1))
         allocate(ocoord_p(thisobs(i_obs)%ncoord, 1))
         allocate(ivar_obs_p(1))
         thisobs(i_obs)%id_obs_p = 0._wp
      end if

      CALL PDAFomi_gather_obs(thisobs(i_obs), dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs(i_obs)%ncoord, local_range, dim_obs)

      deallocate(obs_p, ocoord_p, ivar_obs_p)
   end subroutine init_dim_obs_gen

   subroutine set_obs_p(i_obs, obs_field_p, var_obs, obs_p, var_obs_p)
      integer,  intent(in)    :: i_obs
      real(wp), intent(in)    :: obs_field_p(:)
      real(wp), intent(in)    :: var_obs(:)
      real(wp), intent(inout) :: obs_p(:)
      real(wp), intent(inout) :: var_obs_p(:)

      integer :: cnt
      integer :: i
      integer :: dim_p

      ! set up PE-local observation vector
      cnt = 1
      dim_p = size(obs_field_p)
      do i = 1, dim_p
         if (obs_field_p(i) > missing_value(i_obs)) then
            obs_p(cnt) = obs_field_p(i)
            var_obs_p(cnt) = var_obs(i)
            cnt = cnt + 1
         end if
      end do
   end subroutine set_obs_p

   subroutine set_id_obs_p(i_obs, obs_field_p, id_obs_p)
      use mod_model_pdaf, only: nx, ny
      integer,  intent(in)    :: i_obs
      real(wp), intent(in)    :: obs_field_p(:)
      integer,  intent(inout) :: id_obs_p(:, :)

      integer :: cnt
      integer :: dim_p
      integer :: i
      integer :: offset

      cnt = 1
      dim_p = size(obs_field_p)
      if ((obsvar == 'a') .or. (obsvar == 'b')) then
         offset = 0
      else
         offset = 2*nx*ny
      end if
      do i = 1, dim_p
         if (obs_field_p(i) > missing_value(i_obs)) then
            id_obs_p(1, cnt) = i + offset
            cnt = cnt + 1
         end if
      end do
   end subroutine set_id_obs_p

   subroutine set_ocoord_p(i_obs, offset, obs_field_p, ocoord_p)
      use mod_model_pdaf, only: nx, ny, pi, maooam_model
      integer,  intent(in)    :: i_obs
      integer,  intent(in)    :: offset
      real(wp), intent(in)    :: obs_field_p(:)
      real(wp),  intent(inout) :: ocoord_p(:, :)

      integer :: cnt
      integer :: i
      integer :: dim_p
      real(wp) :: n
      real(wp) :: dx, dy

      cnt = 1
      dim_p = size(obs_field_p)
      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)
      do i = 1, dim_p
         if (obs_field_p(i) > missing_value(i_obs)) then
            ocoord_p(2, cnt) = int((mod(i-1, nx*ny)) / nx)
            ocoord_p(1, cnt) = (int(mod(i-1, nx*ny) -  ocoord_p(2, cnt)*nx))*dx
            ocoord_p(2, cnt) = ocoord_p(2, cnt)*dy
            cnt = cnt + 1
         end if
      end do
   end subroutine set_ocoord_p

   subroutine get_obs_field(i_obs, filename, obs_field)
      use mod_model_pdaf, only: nx, ny
      use netcdf
      integer,      intent(in)    :: i_obs
      character(*), intent(in)    :: filename
      real(wp),     intent(inout) :: obs_field(:)

      integer, save, allocatable :: time_count(:)
      integer :: ncid, varid
      integer :: ierr
      integer :: i
      character(len=5) :: varname(4)

      if (.not. allocated(time_count)) allocate(time_count(n_obs))

      time_count(i_obs) = time_count(i_obs) + 1
      ierr = nf90_open(filename, nf90_nowrite, ncid)

      varname = [character(len=5) :: 'psi_a', 'T_a', 'psi_o', 'T_o']
      do i = 1, 4
         ierr = nf90_inq_varid(ncid, trim(varname(i)), varid)
         ierr = nf90_get_var(ncid, varid, &
                             obs_field((i-1)*nx*ny + 1: i*nx*ny), &
                             start=[1, 1, time_count(i_obs)], count=[nx, ny, 1])
      end do
      ierr = nf90_close(ncid)
   end subroutine get_obs_field

   subroutine get_var_obs(i_obs, filename, var_obs)
      use mod_model_pdaf, only: nx, ny
      use netcdf
      integer,      intent(in)    :: i_obs
      character(*), intent(in)    :: filename
      real(wp),     intent(inout) :: var_obs(:)

      integer :: ncid, varid
      integer :: ierr
      integer :: i
      integer :: i_start, i_end
      character(len=5) :: varname(4)

      ierr = nf90_open(filename, nf90_nowrite, ncid)

      varname = [character(len=5) :: 'psi_a', 'T_a_v', 'psi_o', 'T_o_v']
      i_start = 1
      i_end = 4
      if (obsvar == 'a') then
         i_end = 2
      else if (obsvar == 'o') then
         i_start = 3
      endif
 
      do i = i_start, i_end
         ierr = nf90_inq_varid(ncid, trim(varname(i)), varid)
         ierr = nf90_get_var(ncid, varid, &
                             var_obs((i-1)*nx*ny + 1:i*nx*ny), &
                             start=[1, 1], count=[nx, ny])
      end do
      ierr = nf90_close(ncid)
   end subroutine get_var_obs

   SUBROUTINE obs_op_gridpoint(i_obs, dim_p, dim_obs, state_p, ostate)
      USE PDAFomi, ONLY: PDAFomi_obs_op_gridpoint
      IMPLICIT NONE

      ! *** Arguments ***
      integer, intent(in) :: i_obs
      INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
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
            ONLY: coords_l, local_range, loc_weight, srange

      IMPLICIT NONE

      ! *** Arguments ***
      integer, intent(in)  :: i_obs
      INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector

      ! **********************************************
      ! *** Initialize local observation dimension ***
      ! **********************************************

      CALL PDAFomi_init_dim_obs_l(thisobs_l(i_obs), thisobs(i_obs), coords_l, &
                                  loc_weight, local_range, srange, dim_obs_l)
   END SUBROUTINE init_dim_obs_l

   SUBROUTINE localize_covar(i_obs, dim_p, dim_obs, HP_p, HPH, coords_p)
      ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_localize_covar

      ! Include localization radius and local coordinates
      USE mod_localization_pdaf, &   
            ONLY: local_range, loc_weight, srange

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

      CALL PDAFomi_localize_covar(thisobs(i_obs), dim_p, loc_weight, local_range, srange, &
                                  coords_p, HP_p, HPH)
   END SUBROUTINE localize_covar

   subroutine finalize_obs()
      deallocate(delt_obs_all)
      deallocate(nrows_all, missing_value)
      deallocate(rms_obs_all, thisobs, thisobs_l)
      deallocate(filenames, obsnames)
      deallocate(obs_field_p_all)
      deallocate(var_obs)
   end subroutine finalize_obs
end module mod_observations_pdaf

