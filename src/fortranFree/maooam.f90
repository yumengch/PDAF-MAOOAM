PROGRAM maooam_pdaf
   use mod_kind_pdaf, only: wp
   use mod_config_pdaf, only: screen, read_namelist
   use mod_parallel_pdaf, only: initialize_parallel_pdaf, mype_world, &
                                init_parallel_pdaf, finalize_parallel_pdaf
   use mod_model_pdaf, only: initialize_model, &
                             natm, noc, &
                             finalize_model, total_steps, &
                             writeout, tw, restart_it,&
                             integr, field, field_new, current_time
   use mod_statevector_pdaf, only: initSV
   use mod_init_pdaf, only: init_pdaf, finalize_pdaf
   use mod_ModelWriter_pdaf, only: write_model

   implicit none
   REAL(wp) :: t=0.D0  !< Time variable
   integer  :: it

   ! Initialise parallization
   call initialize_parallel_pdaf()
   call read_namelist()
   call init_parallel_pdaf(screen)
   ! initialise model
   call initialize_model()
   ! initialise state vector
   call initSV()
   ! initialise PDAF
   call init_pdaf(screen)

   if (mype_world == 0) print *, 'Starting the time evolution...'

   t= current_time(1)
   DO it = 1, total_steps
      ! write analysis data into netCDF
      IF ((writeout) .AND. (mod(it-1, tw) < integr%dt)) THEN
         if (mype_world == 0 ) print *, 'a', it
         call write_model(t, 'a', field(1:), natm, noc)
      endif
      ! model integration
      CALL integr%step(field, t, field_new)
      field = field_new
      ! write forcast data into netCDF
      IF ((writeout) .AND. (mod(it, tw) < integr%dt)) THEN
         if (mype_world == 0) print *, 'f', it
         call write_model(t, 'f', field(1:), natm, noc)
      end if

   END DO

   if (mype_world == 0)  then
      print *, 'Evolution finished.'
   end if

   call finalize_model()
   CALL finalize_pdaf()
   call finalize_parallel_pdaf()
END PROGRAM maooam_pdaf
