PROGRAM maooam_pdaf
   use mod_kind_pdaf, only: wp
   use mod_model_pdaf, &
              only: dim_ens, natm, noc, &
                    initialize_model, &
                    finalize_model, total_steps, &
                    writeout, tw, restart_it,&
                    integr, field, field_new
   use mod_config_pdaf, only: verbose, read_namelist
   use mod_parallel_pdaf, only: init_parallel_pdaf, finalize_parallel_pdaf
   use mod_init_pdaf, only: init_pdaf, finalize_pdaf
   use mod_ModelWriter_pdaf, only: write_model
   use mod_assimilate_pdaf, only: assimilate_pdaf
   use mod_U_pdaf, only: collect_dur, distr_dur, next_dur, prepost_dur
   use mod_U_pdafomi_pdaf, only: getobs_dur, dimomi_dur, op_dur

   implicit none
   REAL(wp) :: t=0.D0  !< Time variable
   integer :: timer_model_start, timer_model_end, t_rate
   integer :: timer_PDAF_start, timer_PDAF_end
   real(wp) :: timer_model_dur
   real(wp) :: timer_PDAF_dur

   integer  :: it

   ! Initialise parallization
   call read_namelist()
   call init_parallel_pdaf(dim_ens, verbose)
   ! initialise model
   call initialize_model()
   ! initialise PDAF
   call init_pdaf(verbose)

   t= (restart_it-1)*tw
   if (writeout) &
      call write_model(t, 'f', field(1:), natm, noc)
   print *, 'Starting the time evolution...'
   timer_model_dur = 0.
   timer_PDAF_dur = 0.   
   collect_dur = 0._wp
   distr_dur = 0._wp
   next_dur = 0._wp
   prepost_dur = 0._wp
   getobs_dur = 0._wp
   dimomi_dur = 0._wp
   op_dur = 0._wp
   DO it = 1, total_steps
      IF (writeout .AND. mod(t,tw) < integr%dt) THEN
         call write_model(t, 'a', field(1:), natm, noc)
      endif
      call SYSTEM_CLOCK(timer_model_start)
      CALL integr%step(field, t, field_new)
      field = field_new
      call SYSTEM_CLOCK(timer_model_end, t_rate)
      timer_model_dur = timer_model_dur + &
         (real(timer_model_end, wp) - real(timer_model_start, wp))/real(t_rate, wp)
      IF (writeout .AND. mod(t,tw) < integr%dt) THEN
         call write_model(t, 'f', field(1:), natm, noc)
      end if
      call SYSTEM_CLOCK(timer_PDAF_start)
      call assimilate_pdaf()
      call SYSTEM_CLOCK(timer_PDAF_end, t_rate)
      timer_pdaf_dur = timer_pdaf_dur + &
         (real(timer_pdaf_end, wp) - real(timer_pdaf_start, wp))/real(t_rate, wp)
 
   END DO
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
   call finalize_model()
   CALL finalize_pdaf()
   call finalize_parallel_pdaf()
END PROGRAM maooam_pdaf
