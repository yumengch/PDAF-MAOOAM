PROGRAM maooam_pdaf
   use mod_kind_pdaf, only: wp
   use mod_model_pdaf, &
              only: dim_ens, natm, noc, &
                    initialize_model, &
                    finalize_model, total_steps, &
                    writeout, tw, &
                    integr, field, field_new
   use mod_config_pdaf, only: read_namelist
   use mod_parallel_pdaf, only: init_parallel_pdaf, finalize_parallel_pdaf
   use mod_init_pdaf, only: init_pdaf, finalize_pdaf
   use mod_ModelWriter_pdaf, only: write_model
   use mod_assimilate_pdaf, only: assimilate_pdaf
   implicit none
   REAL(wp) :: t=0.D0  !< Time variable
   integer  :: it
   integer  :: screen = 0 !< verbosity of screen output

   ! Initialise parallization
   call init_parallel_pdaf(dim_ens, screen)
   call read_namelist()
   ! initialise model
   call initialize_model()
   ! initialise PDAF
   call init_pdaf(screen)

   t=0.D0
   if (writeout) &
      call write_model(0._wp, 'f', field(1:), natm, noc)
   print *, 'Starting the time evolution...'
   DO it = 1, total_steps
      IF (writeout .AND. mod(t,tw) < integr%dt) THEN
         call write_model(t, 'a', field(1:), natm, noc)
      endif
      CALL integr%step(field, t, field_new)
      field = field_new
      IF (writeout .AND. mod(t,tw) < integr%dt) THEN
         call write_model(t, 'f', field(1:), natm, noc)
      end if
      call assimilate_pdaf()
   END DO
   print *, 'Evolution finished.'

   call finalize_model()
   CALL finalize_pdaf()
   call finalize_parallel_pdaf()
END PROGRAM maooam_pdaf
