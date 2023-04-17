PROGRAM maooam_pdaf
   use mod_kind_pdaf, only: wp
   use mod_config_pdaf, only: screen, is_freerun, is_strong, read_namelist
   use mod_parallel_pdaf, only: initialize_parallel_pdaf, mype_world, &
                                init_parallel_pdaf, finalize_parallel_pdaf
   use mod_statevector_pdaf, only: component, setField
   use mod_model_pdaf, only: initialize_model, &
                             natm, noc, &
                             finalize_model, total_steps, &
                             writeout, tw, restart_it,&
                             integr, field, field_new, &
                             current_time
   use mod_observations_pdaf, only: initObs => init, obs
   use mod_FilterOptions_pdaf, only: filtertype
   use mod_statevector_pdaf, only: initSV
   use mod_obswriter_pdaf, only: init_obs_writer
   use mod_init_pdaf, only: init_pdaf, finalize_pdaf
   use mod_ModelWriter_pdaf, only: write_model
   use mod_assimilate_pdaf, only: assimilate_pdaf
   use mod_U_pdaf, only: collect_dur, distr_dur, next_dur, prepost_dur
   use mod_U_pdafomi_pdaf, only: getobs_dur, dimomi_dur, op_dur

   implicit none
   REAL(wp) :: t=0.D0  !< Time variable
   integer  :: timer_model_start, timer_model_end, t_rate
   integer  :: timer_PDAF_start, timer_PDAF_end
   real(wp) :: timer_model_dur
   real(wp) :: timer_PDAF_dur

   integer  :: it

   ! Initialise parallization
   call initialize_parallel_pdaf()
   call read_namelist()
   call init_parallel_pdaf(screen)
   ! initialise model
   call initialize_model()
   ! initialise state vector
   call initSV(is_strong)
   ! initialise observations
   call initObs()
   if (filtertype == 100) call init_obs_writer()
   ! initialise PDAF
   call init_pdaf(screen)

   if (mype_world == 0) print *, 'Starting the time evolution...'
   timer_model_dur = 0.
   timer_PDAF_dur = 0.   
   collect_dur = 0._wp
   distr_dur = 0._wp
   next_dur = 0._wp
   prepost_dur = 0._wp
   getobs_dur = 0._wp
   dimomi_dur = 0._wp
   op_dur = 0._wp

   t= current_time(1)
   DO it = 1, total_steps
      IF ((writeout) .AND. (mod(it-1, tw) < integr%dt)) THEN
         if (mype_world == 0 ) print *, 'a', it
         call write_model(t, 'a', field(1:), natm, noc)
      endif

      call SYSTEM_CLOCK(timer_model_start)
      CALL integr%step(field, t, field_new)
      field = field_new
      call SYSTEM_CLOCK(timer_model_end, t_rate)
      timer_model_dur = timer_model_dur + &
         (real(timer_model_end, wp) - real(timer_model_start, wp))/real(t_rate, wp)

      IF ((writeout) .AND. (mod(it, tw) < integr%dt)) THEN
         if (mype_world == 0) print *, 'f', it
         call write_model(t, 'f', field(1:), natm, noc)
      end if

      call SYSTEM_CLOCK(timer_PDAF_start)
      if ((.not. is_Strong) .and. (component == 'b')) then
          call setField('a')
          obs(1)%doassim = 1
          obs(2)%doassim = 0
      end if
      if (.not. is_freerun) call assimilate_pdaf()

      if ((.not. is_strong) .and. (component == 'b')) then
          call setField('o')
          obs(1)%doassim = 0
          obs(2)%doassim = 1
          call assimilate_pdaf()
      end if
      call SYSTEM_CLOCK(timer_PDAF_end, t_rate)
      timer_pdaf_dur = timer_pdaf_dur + &
         (real(timer_pdaf_end, wp) - real(timer_pdaf_start, wp))/real(t_rate, wp)



   END DO

   if (mype_world == 0)  then
      print *, 'Evolution finished.'

      print *, 'model evolution:', timer_model_dur
      print *, 'total assimilation time', timer_pdaf_dur
      print *, 'user-defined state collection:'       , collect_dur
      print *, 'user-defined state distribution:'     , distr_dur
      print *, 'user-defined next observation:'       , next_dur
      print *, 'user-defined prepost processing:'     , prepost_dur
      print *, 'user-defined observation output:'     , getobs_dur
      print *, 'user-defined initialise observation:' , dimomi_dur
      print *, 'user-defined observation operator:'   , op_dur
   end if

   call finalize_model()
   CALL finalize_pdaf()
   call finalize_parallel_pdaf()
END PROGRAM maooam_pdaf
