!$Id: next_observation_pdaf.F90 872 2021-11-22 16:45:59Z lnerger $
!>  Initialize information on next observation
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters
!!
!! The subroutine is called before each forecast phase
!! by PDAF_get_state. It has to initialize the number
!! of time steps until the next available observation
!! (nsteps) and the current model time (time). In
!! addition the exit flag (exit) has to be initialized.
!! It indicates if the data assimilation process is
!! completed such that the ensemble loop in the model
!! routine can be exited.
!!
!! The routine is called from PDAF_get_state by all processes
!!
!! Version for the 2D tutorial model.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
module mod_U_pdaf
use mod_kind_pdaf, only: wp
USE mod_parallel_pdaf, ONLY: mype_world
implicit none

contains
   subroutine init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, uinv, ens_p, status_pdaf)
      use netcdf
      use mod_nfcheck_pdaf, only: check
      USE mod_model_pdaf, &             ! Model variables
         ONLY: nx, ny, &
               psi_a, T_a, psi_o, T_o, &
               toPhysical_A, toPhysical_O, &
               ensscale
      use pdaf_interfaces_module, only: PDAF_sampleens
      implicit none
      ! type of filter to initialize
      integer, intent(in) :: filtertype
      ! pe-local state dimension
      integer, intent(in) :: dim_p
      ! size of ensemble
      integer, intent(in) :: dim_ens
      ! pe-local model state
      real(wp), intent(inout) :: state_p(dim_p)
      ! array not referenced for ensemble filters
      real(wp), intent(inout) :: uinv(dim_ens - 1,dim_ens - 1)
      ! pe-local state ensemble
      real(wp), intent(inout) :: ens_p(dim_p, dim_ens)
      ! pdaf status flag
      integer, intent(inout) :: status_pdaf
      
      ! local variables
      integer :: ncid, varid
      integer :: i, j
      real(wp), allocatable :: eofV(:, :)
      real(wp), allocatable :: svals(:)
      character(len=5) :: varname(4)

      allocate(svals(dim_ens - 1))
      allocate(eofV(dim_p, dim_ens - 1))

      varname = [character(len=5) :: 'psi_a', 'T_a', 'psi_o', 'T_o']
      call check( nf90_open('covariance.nc', nf90_nowrite, ncid) )

      ! read singular values
      call check( nf90_inq_varid(ncid, 'sigma', varid) )
      call check( nf90_get_var(ncid, varid, svals, start=[1], count=[dim_ens - 1]) )

      ! read singular vectors
      do j = 1, dim_ens - 1
         do i = 1, 4
            call check( nf90_inq_varid(ncid, trim(varname(i))//'_svd', varid) )
            call check( nf90_get_var(ncid, varid, eofV((i-1)*nx*ny+1:i*nx*ny, j), start=[1, 1, j], count=[nx, ny, 1]) )
         end do
      end do

      ! convert ensemble mean from spectral space to grid point space
      call toPhysical_A()
      call toPhysical_O()
      state_p(:nx*ny) = reshape(psi_a, [nx*ny])
      state_p(nx*ny+1:2*nx*ny) = reshape(T_a, [nx*ny])
      state_p(2*nx*ny+1:3*nx*ny) = reshape(psi_o, [nx*ny])
      state_p(3*nx*ny+1:4*nx*ny) = reshape(T_o, [nx*ny])

      ! sample ensemble
      ens_p = 0.
      call PDAF_sampleens(dim_p, dim_ens, eofV, svals, state_p, ens_p, verbose=1, flag=status_pdaf)

      ! get ensemble mean
      state_p = 0.
      do j = 1, dim_ens
         state_p = state_p + ens_p(:, j)/dim_ens
      end do

      ! get inflate ensemble spread
      do i = 0, 3
         print *, 'The ensscale of the ', i+1, '-th variable', ensscale(i+1)
         do j = 1, dim_ens
            ens_p(1 + i*nx*ny:(i+1)*nx*ny, j) = ensscale(i+1)*(ens_p(1 + i*nx*ny:(i+1)*nx*ny, j) - &
                                                               state_p(1 + i*nx*ny:(i+1)*nx*ny)) &
                                                + state_p(1 + i*nx*ny:(i+1)*nx*ny)
         end do
      enddo

      call check( nf90_close(ncid) )
      deallocate(eofV, svals)
   end subroutine init_ens_pdaf

   SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)
      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in)   :: stepnow  !< Number of the current time step
      INTEGER, INTENT(out)  :: nsteps   !< Number of time steps until next obs
      INTEGER, INTENT(out)  :: doexit   !< Whether to exit forecasting (1 for exit)
      REAL(wp), INTENT(out) :: time     !< Current model (physical) time

      ! *******************************************************
      ! *** Set number of time steps until next observation ***
      ! *******************************************************
      ! *** End of assimilation process ***
      nsteps = 1         ! No more steps
      doexit = 0         ! Exit assimilation

   END SUBROUTINE next_observation_pdaf

   !!
   SUBROUTINE distribute_state_pdaf(dim_p, state_p)
      USE mod_model_pdaf, &             ! Model variables
         only: nx, ny, psi_a, T_a, psi_o, T_o, &
               toFourier_A, toFourier_O
      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
      REAL(wp), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

      ! *************************************************
      ! *** Initialize model fields from state vector ***
      ! *** for process-local model domain            ***
      !**************************************************
      if (mype_world == 0) print *, 'distribute to atmosphere component'
      psi_a = reshape(state_p(:nx*ny)           , [nx, ny])
      T_a   = reshape(state_p(nx*ny+1:2*nx*ny)  , [nx, ny])
      call toFourier_A(nx, ny)

      if (mype_world == 0) print *, 'distribute to ocean component'
      psi_o = reshape(state_p(2*nx*ny + 1:3*nx*ny), [nx, ny])
      T_o   = reshape(state_p(3*nx*ny+1:4*nx*ny), [nx, ny])
      call toFourier_O(nx, ny)
   END SUBROUTINE distribute_state_pdaf


   SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, &
                                   dim_obs_p, state_p, uinv, ens_p, flag)
      use mod_parallel_pdaf, only: mype_filter
      include 'mpif.h'

      INTEGER, INTENT(in) :: step        !< Current time step (negative for call after forecast)
      INTEGER, INTENT(in) :: dim_p       !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
      INTEGER, INTENT(in) :: dim_ens_p   !< PE-local size of ensemble
      INTEGER, INTENT(in) :: dim_obs_p   !< PE-local dimension of observation vector
      REAL(wp), INTENT(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
      !< (The array 'state_p' is not generally not initialized in the case of SEIK.
      !< It can be used freely here.)
      REAL(wp), INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
      REAL(wp), INTENT(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
      INTEGER, INTENT(in) :: flag        !< PDAF status flag

      real(wp) :: variance_p(dim_p)
      real(wp) :: variance(dim_p)
      real(wp) :: inv_dim_ens, inv_dim_ens1, rmserror_est
      integer :: i

      ! pre- and post-processing of ensemble
      print *, 'Analyze initial state ensemble'

      ! ensemble mean
      inv_dim_ens = 1._wp/dim_ens
      if (dim_ens > 1) then
         inv_dim_ens1 = 1._wp/(dim_ens - 1)
      else
         inv_dim_ens1 = 0._wp
      end if

      state_p = 0._wp
      do i = 1, dim_ens
         state_p = state_p + ens_p(:, i)
      end do
      state_p = state_p*inv_dim_ens

      ! ensemble variance
      variance_p = 0._wp
      do i = 1, dim_ens
         variance_p = variance_p + (ens_p(:, i) - state_p)*(ens_p(:, i) - state_p)
      end do
      variance_p = variance_p*inv_dim_ens1

      variance(:dim_p) = variance_p

      rmserror_est = sqrt(sum(variance)/dim_p)

      if (mype_filter == 0) print*, 'RMS error: ', rmserror_est
   END SUBROUTINE prepoststep_ens_pdaf
end module mod_U_pdaf
