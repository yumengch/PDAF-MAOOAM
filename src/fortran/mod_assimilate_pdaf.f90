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
   SUBROUTINE assimilate_pdaf()
      USE PDAFomi, ONLY: PDAFomi_dealloc
      USE PDAF_mod_filter,ONLY: cnt_steps, step_obs, nsteps
      USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
           ONLY: PDAFomi_assimilate_local, PDAFomi_assimilate_global, &
           PDAFomi_assimilate_lenkf, PDAF_get_localfilter
      USE mod_parallel_pdaf, &       ! Parallelization
           ONLY: mype_world, abort_parallel
      USE mod_filteroptions_pdaf, &         ! Variables for assimilation
           ONLY: filtertype
      use mod_config_pdaf, only: is_strong
      use mod_statevector_pdaf, only: setField, component
      USE mod_U_pdaf, only: collect_state_pdaf, &    ! Collect a state vector from model fields
                            distribute_state_pdaf, &  ! Distribute a state vector to model fields
                            next_observation_pdaf, &  ! Provide time step of next observation
                            prepoststep_ens_pdaf
      USE mod_localization_pdaf, only: init_n_domains_pdaf, &  ! Provide number of local analysis domains
                                  init_dim_l_pdaf, & ! Initialize state dimension for local analysis domain
                                  g2l_state_pdaf, &  ! Get state on local analysis domain from global state
                                  l2g_state_pdaf     ! Update global state from state on local analysis domain
      USE mod_U_PDAFomi_pdaf, only: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
                             obs_op_pdafomi, &         ! Obs. operator for full obs. vector for PE-local domain
                             init_dim_obs_l_pdafomi, & ! Get dimension of obs. vector for local analysis domain
                             localize_covar_pdafomi, & ! Apply localization to covariance matrix in LEnKF
                             init_dim_obs_gen_pdafomi, &
                             get_obs_f
      use mod_observations_pdaf, only: obs
      IMPLICIT NONE
      ! *** Local variables ***
      INTEGER :: status_pdaf          ! PDAF status flag
      integer :: steps, time, doexit
      INTEGER :: localfilter          ! Flag for domain-localized filter (1=true)
      integer :: iobs(1)
      EXTERNAL :: PDAFomi_init_obs_f_cb, & ! Initialize observation vector
      PDAFomi_init_obsvar_cb, &       ! Initialize mean observation error variance
      PDAFomi_prodRinvA_cb        ! Provide product R^-1 A

      ! *********************************
      ! *** Call assimilation routine ***
      ! *********************************

      ! Check  whether the filter is domain-localized
      CALL PDAF_get_localfilter(localfilter)

      ! Call assimilate routine for global or local filter
      IF (localfilter==1) THEN
         CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
              init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_ens_pdaf, init_n_domains_pdaf, &
              init_dim_l_pdaf, init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, &
              next_observation_pdaf, status_pdaf)
      ELSE
         IF (filtertype==8) THEN
            ! localized EnKF has its own OMI interface routine
            CALL PDAFomi_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
                 init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_ens_pdaf, &
                 localize_covar_pdafomi, next_observation_pdaf, status_pdaf)
         else if (filtertype == 100) then
            ! generate observations
            call PDAFomi_generate_obs(collect_state_pdaf, &
                                         distribute_state_pdaf, &
                                         init_dim_obs_gen_pdafomi, &
                                         obs_op_pdafomi, &
                                         get_obs_f, &
                                         prepoststep_ens_pdaf, &
                                         next_observation_pdaf, status_pdaf)
         ELSE
            ! All global filters, except LEnKF
            if ((.not. is_strong) .and. (cnt_steps + 1 == nsteps) .and. (component == 'b')) then
               call setField('o')
               iobs = findloc(obs(:)%obsvar, 'o')
               if (mod(step_obs, obs(iobs(1))%delt_obs) == 0) then
                  CALL PDAF_put_state_estkf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
                                          PDAFomi_init_obs_f_cb, prepoststep_ens_pdaf, &
                                          PDAFomi_prodRinvA_cb, &
                                          PDAFomi_init_obsvar_cb, status_pdaf)

                  ! *** Prepare start of next ensemble forecast ***
                  IF (status_pdaf==0) then
                     CALL PDAF_get_state(steps, time, doexit, next_observation_pdaf, distribute_state_pdaf, &
                     prepoststep_ens_pdaf, status_pdaf)
                     step_obs = step_obs - nsteps
                     nsteps = 1
                  end if
                  CALL PDAFomi_dealloc()
               end if
               call setField('a')
            end if

            CALL PDAFomi_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
                 init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_ens_pdaf, &
                 next_observation_pdaf, status_pdaf)

         END IF
      END IF

      ! Check for errors during execution of PDAF

      IF (status_pdaf /= 0) THEN
         WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
              'ERROR ', status_pdaf, &
              ' in PDAFomi_assimilate - stopping! (PE ', mype_world,')'
         CALL  abort_parallel()
      END IF

   END SUBROUTINE assimilate_pdaf

end module mod_assimilate_pdaf
