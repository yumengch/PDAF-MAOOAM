module mod_nfcheck_pdaf
implicit none
contains
  subroutine check(status)
    use mod_parallel_pdaf, only: abort_parallel
    use netcdf

    ! *** Aruments ***
    integer, intent ( in) :: status   ! Reading status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       call abort_parallel()
    end if

  end subroutine check
end module mod_nfcheck_pdaf
