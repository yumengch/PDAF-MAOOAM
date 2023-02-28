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
implicit none

logical :: firsttime = .true.
integer :: timer_collect_start, timer_collect_end, t_rate
integer :: timer_distr_start, timer_distr_end
integer :: timer_next_start, timer_next_end
integer :: timer_prepost_start, timer_prepost_end
real(wp) :: collect_dur, distr_dur, next_dur, prepost_dur
contains
   subroutine init_ens_pdaf_freerun(filtertype, dim_p, dim_ens, state_p, uinv, ens_p, status_pdaf)
      use netcdf
      USE mod_model_pdaf, &             ! Model variables
         ONLY: nx, ny, psi_a, T_a, psi_o, T_o, toPhysical
      use mod_nfcheck_pdaf, only: check
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
      integer :: ncid, dimid, varid
      integer :: i, j
      integer :: rank
      real(wp), allocatable :: eofV(:, :)
      real(wp), allocatable :: svals(:)
      character(len=5) :: varname(4)

      varname = [character(len=5) :: 'psi_a', 'T_a', 'psi_o', 'T_o']
      call check( nf90_open('covariance.nc', nf90_nowrite, ncid) )
      
      call check( nf90_inq_dimid(ncid, 'rank', dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=rank) )

      allocate(svals(rank))
      call check( nf90_inq_varid(ncid, 'sigma', varid) )
      call check( nf90_get_var(ncid, varid, svals, start=[1], count=[rank]) )

      allocate(eofV(dim_ens - 1, dim_p))

      call toPhysical()
      state_p(:nx*ny) = reshape(psi_a, [nx*ny])
      state_p(nx*ny+1:2*nx*ny) = reshape(T_a, [nx*ny])
      state_p(2*nx*ny+1:3*nx*ny) = reshape(psi_o, [nx*ny])
      state_p(3*nx*ny+1:4*nx*ny) = reshape(T_o, [nx*ny])

      if (dim_ens > 1) then
         do j = 1, dim_ens - 1
            do i = 1, 4
               call check( nf90_inq_varid(ncid, trim(varname(i))//'_svd', varid) )
               call check( nf90_get_var(ncid, varid, eofV(j, (i-1)*nx*ny+1:i*nx*ny), start=[1, 1, 1], count=[nx, ny, 1]) )
            end do
         end do

         call PDAF_sampleens(dim_p, dim_ens, eofV, svals, state_p, ens_p, verbose=1, flag=status_pdaf)
      else
         ens_p(:, 1) = state_p
      end if
      call check( nf90_close(ncid) )
      deallocate(eofV, svals)
   end subroutine init_ens_pdaf_freerun

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
      USE mod_parallel_pdaf, &    ! Parallelization variables
         ONLY: mype_world
      USE mod_model_pdaf, &            ! Model variables
         ONLY: total_steps, is_freerun
      USE mod_observations_pdaf, &
         ONLY: delt_obs_all
      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in)  :: stepnow  !< Number of the current time step
      INTEGER, INTENT(out) :: nsteps   !< Number of time steps until next obs
      INTEGER, INTENT(out) :: doexit   !< Whether to exit forecasting (1 for exit)
      REAL(wp), INTENT(out)    :: time     !< Current model (physical) time

      integer :: delt_obs

      call SYSTEM_CLOCK(timer_next_start)      
      ! *******************************************************
      ! *** Set number of time steps until next observation ***
      ! *******************************************************
      time = 0.0          ! Not used in this implementation
      delt_obs = 90
      if (.not. is_freerun) delt_obs = minval(delt_obs_all)
      IF (stepnow + delt_obs <= total_steps) THEN
         ! *** During the assimilation process ***
         nsteps = delt_obs   ! This assumes a constant time step interval
         doexit = 0          ! Not used in this impl

         IF (mype_world == 0) WRITE (*, '(i7, 3x, a, i7)') &
            stepnow, 'Next observation at time step', stepnow + nsteps
      ELSE
         ! *** End of assimilation process ***
         nsteps = 0          ! No more steps
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
         ONLY: nx, ny, psi_a, T_a, psi_o, T_o, toFourier, toPhysical

      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
      REAL(wp), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

      call SYSTEM_CLOCK(timer_distr_start)

      ! *************************************************
      ! *** Initialize model fields from state vector ***
      ! *** for process-local model domain            ***
      !**************************************************
      psi_a = reshape(state_p(:nx*ny)           , [nx, ny])
      T_a   = reshape(state_p(nx*ny+1:2*nx*ny)  , [nx, ny])
      psi_o = reshape(state_p(2*nx*ny+1:3*nx*ny), [nx, ny])
      T_o   = reshape(state_p(3*nx*ny+1:4*nx*ny), [nx, ny])
      call toFourier(nx, ny)

      call SYSTEM_CLOCK(timer_distr_end, t_rate)
      distr_dur = distr_dur + &
        (real(timer_distr_end, wp) - real(timer_distr_start, wp))/real(t_rate, wp)
   END SUBROUTINE distribute_state_pdaf

   SUBROUTINE distribute_state_pdaf_init(dim_p, state_p)
      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
      REAL(wp), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

   END SUBROUTINE distribute_state_pdaf_init

   SUBROUTINE collect_state_pdaf(dim_p, state_p)

      USE mod_model_pdaf, &             ! Model variables
         ONLY: psi_a, T_a, psi_o, T_o, toPhysical, nx, ny

      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
      REAL(wp), INTENT(inout) :: state_p(dim_p)  !< local state vector

      call SYSTEM_CLOCK(timer_collect_start)
      ! *************************************************
      ! *** Initialize state vector from model fields ***
      ! *** for process-local model domain            ***
      ! *************************************************
      call toPhysical()
      state_p(:nx*ny) = reshape(psi_a, [nx*ny])
      state_p(nx*ny+1:2*nx*ny) = reshape(T_a, [nx*ny])
      state_p(2*nx*ny+1:3*nx*ny) = reshape(psi_o, [nx*ny])
      state_p(3*nx*ny+1:4*nx*ny) = reshape(T_o, [nx*ny])

      call SYSTEM_CLOCK(timer_collect_end, t_rate)
      collect_dur = collect_dur + &
        (real(timer_collect_end, wp) - real(timer_collect_start, wp))/real(t_rate, wp)
   END SUBROUTINE collect_state_pdaf

   SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, &
                                   dim_obs_p, state_p, uinv, ens_p, flag)
      use mod_parallel_pdaf, only: mype_filter, comm_filter, &
                                   npes_filter, MPIerr, MPIstatus
      use mod_model_pdaf, only: nx, ny, dim_state_p, integr
      use mod_StateWriter_pdaf, only: write_state

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
      integer :: off_p
      integer :: i

      call SYSTEM_CLOCK(timer_prepost_start)
      ! pre- and post-processing of ensemble
      if (firsttime) then
         print *, 'Analyze initial state ensemble'
      else
         if (step < 0) then
            print *, 'Analyze and write forecasted state ensemble'
         else
            print *, 'Analyze and write assimilated state ensemble'
         endif
      endif

      ! ensemble mean    
      state_p = 0._wp
      inv_dim_ens = 1._wp/dim_ens
      if (dim_ens > 1) then
         inv_dim_ens1 = 1._wp/(dim_ens - 1)
      else
         inv_dim_ens1 = 0._wp
      end if
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
    

      if (mype_filter /= 0) then
         call mpi_send(variance_p, dim_p, MPI_DOUBLE_PRECISION, 0, &
                       mype_filter, COMM_filter, MPIerr)
      else
         variance(:dim_p) = variance_p

         off_p = 0
         do i = 2, npes_filter
            off_p = off_p + dim_p
            CALL MPI_recv(variance(1 + off_p: 1 + off_p + dim_p), dim_p,  MPI_DOUBLE_PRECISION, i - 1, &
                          i - 1, COMM_filter, MPIstatus, MPIerr)
         end do
      endif

      rmserror_est = sqrt(sum(variance)/dim_state_p)

      if (mype_filter == 0) print*, 'RMS error: ', rmserror_est

      ! if (step <= 0) then
      !    print *, '---writing ensemble forecast---'
      !    call write_state(-step*integr%dt, 'f', ens_p, nx, ny, dim_ens)
      ! else
      !    print *, '---writing ensemble analysis---'
      !    call write_state(step*integr%dt, 'a', ens_p, nx, ny, dim_ens)
      ! endif

      firsttime = .false.
      call SYSTEM_CLOCK(timer_prepost_end, t_rate)
      prepost_dur = prepost_dur + &
          (real(timer_prepost_end, wp) - real(timer_prepost_start, wp))/real(t_rate, wp)
   END SUBROUTINE prepoststep_ens_pdaf
end module mod_U_pdaf
