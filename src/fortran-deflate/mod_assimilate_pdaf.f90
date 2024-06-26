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
           ONLY: PDAFomi_assimilate_local, PDAFomi_assimilate_global, &
           PDAFomi_assimilate_lenkf, PDAF_get_localfilter
      USE mod_parallel_pdaf, &       ! Parallelization
           ONLY: mype_world, abort_parallel
      USE mod_filteroptions_pdaf, &         ! Variables for assimilation
           ONLY: filtertype
      USE mod_U_pdaf, only: collect_state_pdaf, &    ! Collect a state vector from model fields
                            distribute_state_pdaf, &  ! Distribute a state vector to model fields
                            ! distribute_atmospherestate_pdaf, &  ! Distribute a state vector to model fields
                            next_observation_pdaf, &  ! Provide time step of next observation
                            prepoststep_ens_pdaf
      USE mod_localization_pdaf, only: init_n_domains_pdaf, &  ! Provide number of local analysis domains
                                  init_dim_l_pdaf, & ! Initialize state dimension for local analysis domain
                                  g2l_state_pdaf, &  ! Get state on local analysis domain from global state
                                  l2g_state_pdaf     ! Update global state from state on local analysis domain
      USE mod_U_PDAFomi_pdaf, only: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
                             obs_op_pdafomi, &         ! Obs. operator for full obs. vector for PE-local domain
                             init_dim_obs_l_pdafomi, & ! Get dimension of obs. vector for local analysis domain
                             localize_covar_pdafomi    ! Apply localization to covariance matrix in LEnKF
      use mod_statevector_pdaf, only: update_ocean

      IMPLICIT NONE
      integer, intent(in) :: it
      ! *** Local variables ***
      INTEGER :: status_pdaf          ! PDAF status flag
      integer :: steps, time, doexit
      INTEGER :: localfilter          ! Flag for domain-localized filter (1=true)
      EXTERNAL :: PDAFomi_init_obs_f_cb, & ! Initialize observation vector
      PDAFomi_init_obsvar_cb, &       ! Initialize mean observation error variance
      PDAFomi_prodRinvA_cb        ! Provide product R^-1 A
      ! *********************************
      ! *** Call assimilation routine ***
      ! *********************************

      ! Check  whether the filter is domain-localized
      CALL PDAF_get_localfilter(localfilter)
      CALL PDAFomi_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
           init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_ens_pdaf, &
           next_observation_pdaf, status_pdaf)
      ! Check for errors during execution of PDAF

      IF (status_pdaf /= 0) THEN
         WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
              'ERROR ', status_pdaf, &
              ' in PDAFomi_assimilate - stopping! (PE ', mype_world,')'
         CALL  abort_parallel()
      END IF

   END SUBROUTINE assimilate_pdaf

end module mod_assimilate_pdaf
