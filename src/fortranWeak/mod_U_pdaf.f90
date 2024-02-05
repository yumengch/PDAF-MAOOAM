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

logical :: firsttime = .true.
logical :: firsttime_distribute = .true.
integer :: timer_collect_start, timer_collect_end, t_rate
integer :: timer_distr_start, timer_distr_end
integer :: timer_next_start, timer_next_end
integer :: timer_prepost_start, timer_prepost_end
real(wp) :: collect_dur, distr_dur, next_dur, prepost_dur

contains
   subroutine init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, uinv, ens_p, status_pdaf)
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
      
      ens_p = 0.
   end subroutine init_ens_pdaf

   SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)
      USE mod_model_pdaf, ONLY: total_steps
      USE mod_observations_pdaf, ONLY: obs
      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in)   :: stepnow  !< Number of the current time step
      INTEGER, INTENT(out)  :: nsteps   !< Number of time steps until next obs
      INTEGER, INTENT(out)  :: doexit   !< Whether to exit forecasting (1 for exit)
      REAL(wp), INTENT(out) :: time     !< Current model (physical) time

      call SYSTEM_CLOCK(timer_next_start)
      ! *******************************************************
      ! *** Set number of time steps until next observation ***
      ! *******************************************************
      nsteps = minval(obs(:)%delt_obs)

      IF (stepnow + nsteps <= total_steps) THEN
         ! *** During the assimilation process ***
         doexit = 0          ! Not used in this impl

         IF (mype_world == 0) WRITE (*, '(i7, 3x, a, i7)') &
            stepnow, 'Next observation at time step', stepnow + nsteps
      ELSE
         ! *** End of assimilation process ***
         doexit = 1          ! Exit assimilation

         IF (mype_world == 0) WRITE (*, '(i7, 3x, a)') &
            stepnow, 'No more observations - end assimilation'
      END IF

      call SYSTEM_CLOCK(timer_next_end, t_rate)
      next_dur = next_dur + &
        (real(timer_next_end, wp) - real(timer_next_start, wp))/real(t_rate, wp)
   END SUBROUTINE next_observation_pdaf

   !!
   SUBROUTINE distribute_state_pdaf(dim_p, state_p)
      USE mod_model_pdaf, &             ! Model variables
         only: nx, ny, psi_a, T_a, psi_o, T_o, &
               toFourier_A, toFourier_O
      use mod_statevector_pdaf, only: sv_atm, sv_ocean
      USE mod_filteroptions_pdaf, only: filtertype
      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
      REAL(wp), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

      call SYSTEM_CLOCK(timer_distr_start)

      ! *************************************************
      ! *** Initialize model fields from state vector ***
      ! *** for process-local model domain            ***
      !**************************************************
      if (firsttime_distribute) then
         if (mype_world == 0) &
             print *, 'distribute_state_pdaf: starting from restart files'
         firsttime_distribute = .false.
         return
      end if

      if (filtertype == 100) return

      if (mype_world == 0) print *, 'sv_ocean', sv_ocean, 'sv_atm', sv_atm
      if (sv_atm) then
         if (mype_world == 0) print *, 'distribute to atmosphere component'
         psi_a = reshape(state_p(:nx*ny)           , [nx, ny])
         T_a   = reshape(state_p(nx*ny+1:2*nx*ny)  , [nx, ny])
         call toFourier_A(nx, ny)
      end if

      if (sv_ocean) then
         if (mype_world == 0) print *, 'distribute to ocean component'
         psi_o = reshape(state_p(1:nx*ny), [nx, ny])
         T_o   = reshape(state_p(nx*ny+1:2*nx*ny), [nx, ny])
         call toFourier_O(nx, ny)
      end if

      call SYSTEM_CLOCK(timer_distr_end, t_rate)
      distr_dur = distr_dur + &
        (real(timer_distr_end, wp) - real(timer_distr_start, wp))/real(t_rate, wp)
   END SUBROUTINE distribute_state_pdaf


   SUBROUTINE collect_state_pdaf(dim_p, state_p)

      USE mod_model_pdaf, &             ! Model variables
         ONLY: psi_a, T_a, psi_o, T_o, toPhysical_A, toPhysical_O, nx, ny
      use mod_statevector_pdaf, only: sv_atm, sv_ocean
      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
      REAL(wp), INTENT(inout) :: state_p(dim_p)  !< local state vector

      call SYSTEM_CLOCK(timer_collect_start)
      ! *************************************************
      ! *** Initialize state vector from model fields ***
      ! *** for process-local model domain            ***
      ! *************************************************
      state_p = 0.
      if (sv_atm) then
         if (mype_world == 0) print *, 'collect atmosphere'
         call toPhysical_A()
         state_p(:nx*ny) = reshape(psi_a, [nx*ny])
         state_p(nx*ny+1:2*nx*ny) = reshape(T_a, [nx*ny])
      endif

      if (sv_ocean) then
         call toPhysical_O()
         if (mype_world == 0) print *, 'collect ocean'
         state_p(1:nx*ny) = reshape(psi_o, [nx*ny])
         state_p(nx*ny+1:2*nx*ny) = reshape(T_o, [nx*ny])
      endif

      call SYSTEM_CLOCK(timer_collect_end, t_rate)
      collect_dur = collect_dur + &
        (real(timer_collect_end, wp) - real(timer_collect_start, wp))/real(t_rate, wp)
   END SUBROUTINE collect_state_pdaf


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

      call SYSTEM_CLOCK(timer_prepost_start)

      ! pre- and post-processing of ensemble
      if (firsttime) then
         print *, 'Analyze initial state ensemble'
      else
         if (step < 0) then
            print *, 'Analyze forecasted state ensemble'
         else
            print *, 'Analyze assimilated state ensemble'
         endif
      endif

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

      firsttime = .false.
      call SYSTEM_CLOCK(timer_prepost_end, t_rate)

      prepost_dur = prepost_dur + &
          (real(timer_prepost_end, wp) - real(timer_prepost_start, wp))/real(t_rate, wp)
   END SUBROUTINE prepoststep_ens_pdaf
end module mod_U_pdaf
