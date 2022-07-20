PROGRAM maooam_pdaf
  use mod_kind_pdaf, only: wp
  use mod_model_pdaf, &
      only: dim_ens, &
            initialize_model, &
            finalize_model, total_steps, &
            integr, field, field_new
  use mod_parallel_pdaf, only: init_parallel_pdaf, finalize_parallel_pdaf
  use mod_init_pdaf, only: init_pdaf, finalize_pdaf
  use mod_assimilate_pdaf, only: assimilate_pdaf
  implicit none
  REAL(wp) :: t=0.D0  !< Time variable
  integer :: it
  integer :: screen = 0 !< verbosity of screen output

  ! Initialise parallization
  call init_parallel_pdaf(dim_ens, screen)
  ! initialise model
  call initialize_model()
  ! initialise PDAF
  call init_pdaf(screen)

  t=0.D0
  print *, 'Starting the time evolution...'
  DO it = 1, total_steps
     CALL integr%step(field,t, field_new)
     field = field_new 
     call assimilate_pdaf()
  END DO
  print *, 'Evolution finished.'

  call finalize_model()
  CALL finalize_pdaf()
  call finalize_parallel_pdaf()
END PROGRAM maooam_pdaf
