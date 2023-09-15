module mod_config_pdaf
   use mod_kind_pdaf, only: wp
implicit none

integer :: screen
logical :: is_freerun, is_strong, do_hybrid_ens
real(wp) :: hybrid_coeff

contains
   subroutine read_namelist
      use mod_filteroptions_pdaf, only: filter_nml
      use mod_inflation_pdaf, only: infl_nml
      use mod_model_pdaf, only: model_nml
      use mod_statevector_pdaf, only: state_vector_nml

      namelist /PDAF_nml/ screen, is_strong, is_freerun, do_hybrid_ens, hybrid_coeff
      open (20, file='PDAF_config.nml')
      read(20, nml=filter_nml)
      rewind(20)
      read(20, nml=PDAF_nml)
      rewind(20)
      read(20, nml=infl_nml)
      rewind(20)
      read(20, nml=model_nml)
      rewind(20)
      read(20, nml=state_vector_nml)
      rewind(20)
      close(20)
   end subroutine read_namelist
end module mod_config_pdaf
