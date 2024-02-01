module mod_config_pdaf
use mod_kind_pdaf, only: wp
implicit none

integer :: screen

contains
   subroutine read_namelist
      use mod_model_pdaf, only: model_nml

      namelist /PDAF_nml/ screen
      open (20, file='PDAF_config.nml')
      read(20, nml=PDAF_nml)
      rewind(20)
      read(20, nml=model_nml)
      close(20)
   end subroutine read_namelist
end module mod_config_pdaf
