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
use mod_parallel_pdaf, only: abort_parallel
implicit none

integer   :: n_obs

type :: obs_t
   character ::  obsvar
   character(len=50) :: filename
   character(len=50) :: filename_var

   integer  :: doassim
   integer  :: delt_obs
   integer  :: obs_den = 8
   integer  :: dim_obs

   real(wp) :: rms_obs
   real(wp) :: missing_value

   real(wp),    allocatable :: obs_field_p(:, :, :)
   real(wp),    allocatable :: var_obs(:, :, :)
   real(wp),    ALLOCATABLE :: A(:, :)
   real(wp),    ALLOCATABLE :: Ainv(:, :)

   integer :: file_timecount = 0
   integer :: file_timestep
end type obs_t


type(obs_t),   ALLOCATABLE, target :: obs(:)
REAL(wp), ALLOCATABLE :: svals(:)      ! Singular values of Ainv
REAL(wp), ALLOCATABLE :: work(:)       ! Work array for SYEVTYPE
integer :: ocean_interval

contains
   subroutine init()
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
      ocean_interval = maxval(obs(:)%delt_obs)
   end subroutine init

   subroutine init_single_obs(i_obs, nmlname)
      use mod_model_pdaf   , only: nx, ny
      use mod_parallel_pdaf, only: dim_ens => n_modeltasks

      integer     , intent(in) :: i_obs
      character(*), intent(in) :: nmlname

      character         :: obsvar
      character(len=50) :: filename
      character(len=50) :: filename_var

      integer  :: doassim
      integer  :: delt_obs
      integer  :: obs_den

      real(wp) :: rms_obs

      integer  :: file_timestep
      integer  :: file_timecount

      integer  :: nVar, nxo, nyo

      namelist /setup_nml/ obsvar, filename, filename_var, &
                           doassim, delt_obs, file_timestep, &
                           rms_obs, obs_den, file_timecount

      ! read options for the observation
      open (20, file=nmlname)
      read(20, nml=setup_nml)
      rewind(20)
      close(20)

      ! assign observation options
      obs(i_obs)%obsvar = obsvar  
      obs(i_obs)%filename = filename
      obs(i_obs)%filename_var = filename_var

      obs(i_obs)%delt_obs = delt_obs
      obs(i_obs)%rms_obs = rms_obs

      ! observation types
      obs(i_obs)%doassim = doassim

      obs(i_obs)%obs_den = obs_den
      obs(i_obs)%file_timecount = file_timecount
      obs(i_obs)%file_timestep = file_timestep

      ! Size of domain for periodicity for disttype=1
      ! (<0 for no periodicity)
      obs(i_obs)%missing_value = -99999

      ! allocate array for observation and its variance
      nVar = 2
      nxo = nx/obs(i_obs)%obs_den + 1
      nyo = ny/obs(i_obs)%obs_den + 1
      allocate(obs(i_obs)%obs_field_p(nxo, nyo, nVar))
      allocate(obs(i_obs)%var_obs(nxo, nyo, nVar))
      allocate(obs(i_obs)%A(2*nxo*nyo, dim_ens))
      allocate(obs(i_obs)%Ainv(2*nxo*nyo, 2*nxo*nyo))

      if (.not. allocated(svals)) allocate(svals(2*nxo*nyo), work(6*nxo*nyo))
      obs(i_obs)%obs_field_p = 0.
   end subroutine init_single_obs

   subroutine init_dim_obs(i_obs, step, dim_obs)
      use mod_model_pdaf, only: nx, ny

      integer, intent(in)  :: i_obs
      integer, intent(in)  :: step
      integer, intent(out) :: dim_obs

      integer  :: nxo, nyo

      print *, 'assimilate ', obs(i_obs)%obsvar, ' component'

      nxo = nx/obs(i_obs)%obs_den + 1
      nyo = ny/obs(i_obs)%obs_den + 1

      call get_obs_field(step, i_obs)
      ! read observation variance
      call get_var_obs(i_obs)
      ! count valid observations
      dim_obs = 2*nxo*nyo
      obs(i_obs)%dim_obs = dim_obs
   end subroutine init_dim_obs

   subroutine get_obs_field(step, i_obs)
      use mod_nfcheck_pdaf, only: check
      use mod_model_pdaf, only: nx, ny
      use netcdf
      integer, intent(in)  :: step
      integer,      intent(in)    :: i_obs

      integer :: ncid, varid
      integer :: i, cnt
      integer :: i_start, i_end
      integer :: nxo, nyo
      integer :: timecount
      character(len=5) :: varname(4)

      timecount = obs(i_obs)%file_timecount + step/obs(i_obs)%file_timestep

      varname = [character(len=5) :: 'psi_a', 'T_a', 'psi_o', 'T_o']

      call check( nf90_open(trim(obs(i_obs)%filename), nf90_nowrite, ncid) )
      ! choose observed variables
      i_start = 1
      i_end = 4
      if (obs(i_obs)%obsvar == 'a') then
         i_end = 2
      else if (obs(i_obs)%obsvar == 'o') then
         i_start = 3
      else
         print *, '...obsvar is not correct...', i_obs, obs(i_obs)%obsvar
         call abort_parallel()
      endif

      print *, 'reading', i_obs, obs(i_obs)%obsvar, timecount
      nxo = nx/obs(i_obs)%obs_den + 1
      nyo = ny/obs(i_obs)%obs_den + 1

      ! read observations
      cnt = 1
      do i = i_start, i_end
         call check( nf90_inq_varid(ncid, trim(varname(i)), varid) )
         call check( nf90_get_var(ncid, varid, &
                                  obs(i_obs)%obs_field_p(:, :, cnt), &
                                  start=[1, 1, timecount], &
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

   SUBROUTINE obs_op_Aonly(dim_p, dim_obs, state_p, ostate)
      use mod_model_pdaf, only: nx, ny
      IMPLICIT NONE
      ! *** Arguments ***
      INTEGER,  INTENT(in)    :: dim_p                 !< PE-local state dimension
      INTEGER,  INTENT(in)    :: dim_obs               !< Dimension of full observed state (all observed fields)
      REAL(wp), INTENT(in)    :: state_p(dim_p)        !< PE-local model state
      REAL(wp), INTENT(inout) :: ostate(dim_obs)       !< Full observed state

      integer :: nxo, nyo, nvar
      integer :: i, j, k, l, cnt
      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************
      nxo = nx/obs(1)%obs_den + 1
      nyo = ny/obs(1)%obs_den + 1
      nvar = 2
      cnt = 1
      do k = 1, nvar
         ! loop over spatial domain
         do j = 1, nyo
            do i = 1, nxo
               l = (k-1)*nx*ny + nx*obs(1)%obs_den*(j-1) + &
                        (i-1)*obs(1)%obs_den + 1
               ostate(cnt) = state_p(l)
               cnt = cnt + 1
            end do
         end do
      end do
   END SUBROUTINE obs_op_Aonly

   SUBROUTINE obs_op_A(dim_p, dim_obs, state_p, ostate, is_mean)
      use mod_model_pdaf, only: nx, ny
      use mod_model_pdaf, only: psi_o_f, T_o_f
      use mod_inflation_pdaf, only: alpha_O
      IMPLICIT NONE
      ! *** Arguments ***
      INTEGER,  INTENT(in)    :: dim_p                 !< PE-local state dimension
      INTEGER,  INTENT(in)    :: dim_obs               !< Dimension of full observed state (all observed fields)
      REAL(wp), INTENT(in)    :: state_p(dim_p)        !< PE-local model state
      REAL(wp), INTENT(inout) :: ostate(dim_obs)       !< Full observed state
      logical, intent(in)  :: is_mean

      integer :: nxo, nyo, nvar
      integer :: i, j, k, l, cnt
      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************
      nxo = nx/obs(1)%obs_den + 1
      nyo = ny/obs(1)%obs_den + 1
      nvar = 2
      cnt = 1
      do k = 1, nvar
         ! loop over spatial domain
         do j = 1, nyo
            do i = 1, nxo
               l = (k-1)*nx*ny + nx*obs(1)%obs_den*(j-1) + &
                        (i-1)*obs(1)%obs_den + 1
               ostate(cnt) = state_p(l)
               cnt = cnt + 1
            end do
         end do
      end do

      cnt = 2*nxo*nyo + 1
      ! loop over spatial domain
      do j = 1, nyo
         do i = 1, nxo
            ostate(cnt) = psi_o_f((i-1)*obs(2)%obs_den + 1, (j-1)*obs(2)%obs_den + 1)
            ostate(nxo*nyo + cnt) = T_o_f((i-1)*obs(2)%obs_den + 1, (j-1)*obs(2)%obs_den + 1)
            cnt = cnt + 1
         end do
      end do

      if (.not. is_mean) ostate(2*nxo*nyo+1:) = alpha_O*ostate(2*nxo*nyo+1:)
   END SUBROUTINE obs_op_A

   SUBROUTINE obs_op_O(dim_p, dim_obs, state_p, ostate, is_mean)
      use mod_model_pdaf, only: nx, ny
      use mod_model_pdaf, only: psi_a_f, T_a_f
      use mod_inflation_pdaf, only: alpha_A
      IMPLICIT NONE
      ! *** Arguments ***
      INTEGER,  INTENT(in)    :: dim_p                 !< PE-local state dimension
      INTEGER,  INTENT(in)    :: dim_obs               !< Dimension of full observed state (all observed fields)
      REAL(wp), INTENT(in)    :: state_p(dim_p)        !< PE-local model state
      REAL(wp), INTENT(inout) :: ostate(dim_obs)       !< Full observed state
      logical, intent(in)  :: is_mean

      integer :: nxo, nyo, nvar
      integer :: i, j, k, l, cnt
      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************
      nxo = nx/obs(1)%obs_den + 1
      nyo = ny/obs(1)%obs_den + 1

      cnt = 1
      ! loop over spatial domain
      do j = 1, nyo
         do i = 1, nxo
            ostate(cnt) = psi_a_f((i-1)*obs(1)%obs_den + 1, (j-1)*obs(1)%obs_den + 1)
            ostate(cnt + nxo*nyo) = T_a_f((i-1)*obs(1)%obs_den + 1, (j-1)*obs(1)%obs_den + 1)
            cnt = cnt + 1
         end do
      end do

      if (.not. is_mean) ostate(:2*nxo*nyo) = alpha_A*ostate(:2*nxo*nyo)

      nvar = 2
      cnt = 2*nxo*nyo + 1
      do k = 1, nvar
         ! loop over spatial domain
         do j = 1, nyo
            do i = 1, nxo
               l = (k-1)*nx*ny + nx*obs(2)%obs_den*(j-1) + &
                        (i-1)*obs(2)%obs_den + 1
               ostate(cnt) = state_p(l)
               cnt = cnt + 1
            end do
         end do
      end do
   END SUBROUTINE obs_op_O

   SUBROUTINE prodRinvA_Aonly(step, dim_obs_p, dim_ens, obs_p, A_p, C_p)
      use mod_model_pdaf, only: nx, ny
      INTEGER,  INTENT(in)  :: step                    ! Current time step
      INTEGER,  INTENT(in)  :: dim_obs_p               ! PE-local dimension of obs. vector
      INTEGER,  INTENT(in)  :: dim_ens                 ! Ensemble size
      REAL(wp), INTENT(in)  :: obs_p(dim_obs_p)        ! PE-local vector of observations
      REAL(wp), INTENT(in)  :: A_p(dim_obs_p, dim_ens) ! Input matrix from analysis routine
      REAL(wp), INTENT(out) :: C_p(dim_obs_p, dim_ens) ! Output matrix

      integer :: nxo, nyo
      integer :: cnt, i, j, k

      nxo = nx/obs(1)%obs_den + 1
      nyo = ny/obs(1)%obs_den + 1
      ! first 2*nxo*nyo rows
      cnt = 1
      do k = 1, 2
         ! loop over spatial domain
         do j = 1, nyo
            do i = 1, nxo
               if (obs(1)%var_obs(i, j, k) < 1e-12) then
                  C_p(cnt, :) = A_p(cnt, :)*1e14
               else
                  C_p(cnt, :) = A_p(cnt, :)/obs(1)%var_obs(i, j, k)/obs(1)%rms_obs
               end if
               cnt = cnt + 1
            end do
         end do
      end do
   END SUBROUTINE prodRinvA_Aonly

   SUBROUTINE prodRinvA_A(step, dim_obs_p, dim_ens, obs_p, A_p, C_p)
      use mod_model_pdaf, only: nx, ny
      use mod_inflation_pdaf, only: forget, alpha_O
      INTEGER,  INTENT(in)  :: step                    ! Current time step
      INTEGER,  INTENT(in)  :: dim_obs_p               ! PE-local dimension of obs. vector
      INTEGER,  INTENT(in)  :: dim_ens                 ! Ensemble size
      REAL(wp), INTENT(in)  :: obs_p(dim_obs_p)        ! PE-local vector of observations
      REAL(wp), INTENT(in)  :: A_p(dim_obs_p, dim_ens) ! Input matrix from analysis routine
      REAL(wp), INTENT(out) :: C_p(dim_obs_p, dim_ens) ! Output matrix

      integer :: nxo, nyo
      integer :: cnt, i, j, k
      real(wp) :: coeff

      INTEGER :: lib_info                ! Status flag for LAPACK calls
      INTEGER :: ldwork                  ! Size of work array for syevTYPE

      ! A_p is basically Y^s matrix
      nxo = nx/obs(1)%obs_den + 1
      nyo = ny/obs(1)%obs_den + 1
      ! first 2*nxo*nyo rows
      cnt = 1
      do k = 1, 2
         ! loop over spatial domain
         do j = 1, nyo
            do i = 1, nxo
               if (obs(1)%var_obs(i, j, k) < 1e-12) then
                  C_p(cnt, :) = A_p(cnt, :)*1e14
               else
                  C_p(cnt, :) = A_p(cnt, :)/obs(1)%var_obs(i, j, k)/obs(1)%rms_obs
               end if
               cnt = cnt + 1
            end do
         end do
      end do

      ! second 2*nxo*nyo rows
      if (abs(alpha_O - 1.) > 1e-6 .and. abs(alpha_O) > 1e-6) then
         coeff = (1. - alpha_O*alpha_O)/forget/(dim_ens - 1)/alpha_O/alpha_O
         obs(2)%A = coeff*matmul(A_p(2*nxo*nyo + 1:, :), transpose(A_p(2*nxo*nyo + 1:, :)))
         cnt = 1
         do k = 1, 2
            ! loop over spatial domain
            do j = 1, nyo
               do i = 1, nxo
                  if (obs(2)%var_obs(i, j, k) < 1e-12) obs(2)%var_obs(i, j, k) = 1e-14
                  obs(2)%A(cnt, cnt) = obs(2)%A(cnt, cnt) + obs(2)%var_obs(i, j, k)*obs(2)%rms_obs
                  cnt = cnt + 1
               end do
            end do
         end do
         ldwork = 6*nxo*nyo
         CALL dsyev('V', 'L', 2*nxo*nyo, obs(2)%A, 2*nxo*nyo, svals, work, ldwork, lib_info)
         do i = 1, 2*nxo*nyo
            obs(2)%Ainv(:, i) = obs(2)%A(:, i)/svals(i)
         end do
         C_p(2*nxo*nyo + 1:, :) = matmul(matmul(obs(2)%Ainv, &
                                                transpose(obs(2)%A)), &
                                                A_p(2*nxo*nyo + 1:, :))
      else if (abs(alpha_O) < 1e-6) then
         C_p(2*nxo*nyo + 1:, :) = 0.
      else
         cnt = 2*nxo*nyo + 1
         do k = 1, 2
            ! loop over spatial domain
            do j = 1, nyo
               do i = 1, nxo
                  if (obs(2)%var_obs(i, j, k) < 1e-12) then
                     C_p(cnt, :) = A_p(cnt, :)*1e14
                  else
                     C_p(cnt, :) = A_p(cnt, :)/obs(2)%var_obs(i, j, k)/obs(2)%rms_obs
                  end if
                  cnt = cnt + 1
               end do
            end do
         end do
      end if
   END SUBROUTINE prodRinvA_A

   SUBROUTINE prodRinvA_O(step, dim_obs_p, dim_ens, obs_p, A_p, C_p)
      use mod_model_pdaf, only: nx, ny
      use mod_inflation_pdaf, only: forget, alpha_A
      INTEGER, INTENT(in) :: step                ! Current time step
      INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
      INTEGER, INTENT(in) :: dim_ens             ! Ensemble size
      REAL(wp), INTENT(in)    :: obs_p(dim_obs_p)    ! PE-local vector of observations
      REAL(wp), INTENT(in)    :: A_p(dim_obs_p, dim_ens) ! Input matrix from analysis routine
      REAL(wp), INTENT(out)   :: C_p(dim_obs_p, dim_ens) ! Output matrix

      integer :: nxo, nyo
      integer :: cnt, i, j, k
      real(wp) :: coeff

      INTEGER :: lib_info                ! Status flag for LAPACK calls
      INTEGER :: ldwork                  ! Size of work array for syevTYPE
      ! A_p is basically Y^s matrix

      nxo = nx/obs(1)%obs_den + 1
      nyo = ny/obs(1)%obs_den + 1
      ! first 2*nxo*nyo rows
      cnt = 2*nxo*nyo + 1
      do k = 1, 2
         ! loop over spatial domain
         do j = 1, nyo
            do i = 1, nxo
               if (obs(2)%var_obs(i, j, k) < 1e-12) then
                  C_p(cnt, :) = A_p(cnt, :)*1e14
               else
                  C_p(cnt, :) = A_p(cnt, :)/obs(2)%var_obs(i, j, k)/obs(2)%rms_obs
               end if
               cnt = cnt + 1
            end do
         end do
      end do
      ! second 2*nxo*nyo rows
      if (abs(alpha_A - 1.) > 1e-6 .and. abs(alpha_A) > 1e-6) then
         coeff = (1. - alpha_A*alpha_A)/forget/(dim_ens - 1)/alpha_A/alpha_A
         obs(1)%A = coeff*matmul(A_p(:2*nxo*nyo, :), transpose(A_p(:2*nxo*nyo, :)))
         cnt = 1
         do k = 1, 2
            ! loop over spatial domain
            do j = 1, nyo
               do i = 1, nxo
                  if (obs(1)%var_obs(i, j, k) < 1e-12) obs(1)%var_obs(i, j, k) = 1e-14
                  obs(1)%A(cnt, cnt) = obs(1)%A(cnt, cnt) + obs(1)%var_obs(i, j, k)*obs(1)%rms_obs
                  cnt = cnt + 1
               end do
            end do
         end do
         ldwork = 6*nxo*nyo
         CALL dsyev('V', 'L', 2*nxo*nyo, obs(1)%A, 2*nxo*nyo, svals, work, ldwork, lib_info)
         do i = 1, 2*nxo*nyo
            obs(1)%Ainv(:, i) = obs(1)%A(:, i)/svals(i)
         end do
         C_p(:2*nxo*nyo, :) = matmul(matmul(obs(1)%Ainv, &
                                            transpose(obs(1)%A)), &
                                     A_p(:2*nxo*nyo, :))
      else if (abs(alpha_A) < 1e-6 ) then
         C_p(:2*nxo*nyo, :) = 0.
      else
         cnt = 1
         do k = 1, 2
            ! loop over spatial domain
            do j = 1, nyo
               do i = 1, nxo
                  if (obs(1)%var_obs(i, j, k) < 1e-12) then
                     C_p(cnt, :) = A_p(cnt, :)*1e14
                  else
                     C_p(cnt, :) = A_p(cnt, :)/obs(1)%var_obs(i, j, k)/obs(1)%rms_obs
                  end if
                  cnt = cnt + 1
               end do
            end do
         end do
      end if
   END SUBROUTINE prodRinvA_O

   SUBROUTINE init_obs(i_obs, step, dim_obs_p, observation_p)
      use mod_model_pdaf, only: nx, ny
      integer,  intent(in)   :: i_obs
      INTEGER,  INTENT(in)   :: step             ! Current time step
      INTEGER,  INTENT(in)   :: dim_obs_p        ! PE-local dimension of obs. vector
      REAL(wp), INTENT(out)  :: observation_p(dim_obs_p) ! PE-local observation vector

      integer :: nxo, nyo, nvar
      integer :: i, j, k, cnt

      nxo = nx/obs(i_obs)%obs_den + 1
      nyo = ny/obs(i_obs)%obs_den + 1
      nvar = 2
      cnt = 1 + (i_obs - 1)*nvar*nxo*nyo
      do k = 1, nvar
         ! loop over spatial domain
         do j = 1, nyo
            do i = 1, nxo
               observation_p(cnt) = obs(i_obs)%obs_field_p(i, j, k)
               cnt = cnt + 1
            end do
         end do
      end do
   end subroutine init_obs

   subroutine finalize_obs()
      integer :: i_obs

      do i_obs = 1, n_obs
         deallocate(obs(i_obs)%obs_field_p)
         deallocate(obs(i_obs)%var_obs)
         deallocate(obs(i_obs)%A)
         deallocate(obs(i_obs)%Ainv)
      end do
      deallocate(obs)
      deallocate(svals, work)
   end subroutine finalize_obs
end module mod_observations_pdaf
