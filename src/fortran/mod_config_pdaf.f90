module mod_config_pdaf
implicit none

contains
   subroutine read_namelist
      use mod_filteroptions_pdaf, only: filter_nml
      use mod_inflation_pdaf, only: infl_nml
      use mod_localization_pdaf, only: local_nml
      use mod_model_pdaf, only: model_nml

      open (20, file='PDAF_config.nml')
      read(20, nml=filter_nml)
      rewind(20)
      read(20, nml=infl_nml)
      rewind(20)
      read(20, nml=local_nml)
      rewind(20)
      read(20, nml=model_nml)
      close(20)
   end subroutine read_namelist
end module mod_config_pdaf
