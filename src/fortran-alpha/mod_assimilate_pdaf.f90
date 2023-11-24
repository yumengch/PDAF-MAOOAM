!$Id: assimilate_pdaf.F90 870 2021-11-22 14:02:55Z lnerger $
!>  Routine to call PDAF for analysis step
!!
!! This routine is called during the model integrations at each time
!! step. It calls the filter-specific assimilation routine of PDAF
!! (PDAF_assimilate_X), which checks whether the forecast phase is
!! completed. If so, the analysis step is computed inside PDAF
!!
!! __Revision history:__
!! * 2013-08 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!

module mod_assimilate_pdaf
implicit none

contains
   SUBROUTINE assimilate_pdaf(it)
      USE PDAFomi, ONLY: PDAFomi_dealloc
      USE PDAF_mod_filter,ONLY: cnt_steps, step_obs, nsteps
      USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
           ONLY: PDAF_assimilate_etkf
      USE mod_parallel_pdaf, &       ! Parallelization
           ONLY: mype_world, abort_parallel
      USE mod_U_pdaf, only: collect_oceanstate_pdaf, &    ! Collect a state vector from model fields
                            distribute_oceanstate_pdaf, &  ! Distribute a state vector to model fields
                            collect_atmospherestate_pdaf, &  ! Distribute a state vector to model fields
                            distribute_atmospherestate_pdaf, &  ! Distribute a state vector to model fields
                            next_observation_pdaf, &  ! Provide time step of next observation
                            prepoststep_ens_pdaf
      USE mod_U_PDAFomi_pdaf, only: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
                             obs_op_pdafomi, &         ! Obs. operator for full obs. vector for PE-local domain
                             init_obs_pdafomi, & ! Get dimension of obs. vector for local analysis domain
                             prodRinvA_pdafomi, &    ! Apply localization to covariance matrix in LEnKF
                             init_obsvar_pdafomi
      use mod_statevector_pdaf, only: update_ocean, update_both
      use mod_model_pdaf, only: toPhysical_A, toPhysical_O

      IMPLICIT NONE
      integer, intent(in) :: it
      ! *** Local variables ***
      INTEGER :: status_pdaf          ! PDAF status flag
      integer :: steps, time, doexit
      ! *********************************
      ! *** Call assimilation routine ***
      ! *********************************

      ! Call assimilate routine for global or local filter
      call toPhysical_A()
      call toPhysical_O()
      if (update_both .and. cnt_steps + 1 == nsteps) then
         update_ocean = .true.
         CALL PDAF_put_state_etkf(collect_oceanstate_pdaf, &
                                  init_dim_obs_pdafomi, &
                                  obs_op_pdafomi, &
                                  init_obs_pdafomi, &
                                  prepoststep_ens_pdaf, prodRinvA_pdafomi, &
                                  init_obsvar_pdafomi, status_pdaf)

         ! *** Prepare start of next ensemble forecast ***
         IF (status_pdaf==0) then
            CALL PDAF_get_state(steps, time, doexit, &
                                next_observation_pdaf, &
                                distribute_oceanstate_pdaf, &
                                prepoststep_ens_pdaf, status_pdaf)
            step_obs = step_obs - nsteps
            nsteps = 1
         else
            WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
                  'ERROR ', status_pdaf, &
                  ' in PDAF_put_state_etkf - stopping! (PE ', mype_world,')'
            CALL  abort_parallel()
         end if
         CALL PDAFomi_dealloc()
         update_ocean = .false.
      end if
      
      call PDAF_assimilate_etkf(collect_atmospherestate_pdaf, &
                                distribute_atmospherestate_pdaf, &
                                init_dim_obs_pdafomi, &
                                obs_op_pdafomi, init_obs_pdafomi, &
                                prepoststep_ens_pdaf, prodRinvA_pdafomi, &
                                init_obsvar_pdafomi, next_observation_pdaf, status_pdaf)
      ! Check for errors during execution of PDAF

      IF (status_pdaf /= 0) THEN
         WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
              'ERROR ', status_pdaf, &
              ' in PDAF_assimilate_etkf - stopping! (PE ', mype_world,')'
         CALL  abort_parallel()
      END IF

   END SUBROUTINE assimilate_pdaf

end module mod_assimilate_pdaf
