module mod_StateWriter_pdaf
use mod_kind_pdaf, only: wp
use netcdf
use mod_nfcheck_pdaf, only: check
implicit none

integer :: time_count
integer :: ncid
integer :: dimid(4)
integer :: varid_time
integer :: varid(8)

contains
   subroutine init_state_writer(filename, nx, ny, dim_ens)
      character(*), intent(in) :: filename
      integer, intent(in)      :: nx, ny, dim_ens

      integer :: i, j
      integer :: dimids(4, 4)

      character         :: vartype(2) = ['f', 'a']
      character(len=8)  :: typename(2) = ['forecast', 'analysis']
      character(len=5)  :: varname(4)
      character(len=40) :: standard_name(4)
      character(len=50) :: long_name(4)

      call check( nf90_create(filename, nf90_netcdf4, ncid) )
      ! set global attributes
      call setAttrs
      ! initialise dimensions
      call check( nf90_def_dim( ncid, 'time', nf90_unlimited, dimid(1) ) )
      call check( nf90_def_var( ncid, 'time', nf90_double, dimid(1), varid_time ) )
      call check( nf90_put_att( ncid, varid_time, 'long_name', 'time' ) )
      call check( nf90_put_att( ncid, varid_time, 'units', 'seconds since 1900-1-1 0:0:0' ) )
      call check( nf90_def_dim( ncid, 'ens', dim_ens, dimid(2) ) )
      call check( nf90_def_dim( ncid, 'nx', nx, dimid(3) ) )
      call check( nf90_def_dim( ncid, 'ny', ny, dimid(4) ) )
      ! initialise output variables
      call getVarAttrs(varname, standard_name, long_name, dimids)
      do i = 1, 2
         do j = 1, 4
            call check( nf90_def_var(ncid, trim(varname(j))//'_'//vartype(i), nf90_double, dimids(j, :), varid((i-1)*4+j)) )
            call check( nf90_put_att(ncid, varid((i-1)*4+j), 'standard_name', trim(standard_name(j))//'_'//trim(typename(i))) )
            call check( nf90_put_att(ncid, varid((i-1)*4+j), 'long_name', trim(typename(i))//' of '//trim(long_name(j))) )
         end do
      end do
      call check( NF90_ENDDEF(ncid) )
      time_count = 0
   end subroutine init_state_writer

   subroutine setAttrs()
      integer :: ierr
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'Conventions', 'CF-1.8')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'title', &
      'NetCDF output from MAOOAM-pyPDAF')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'institution', &
      'NCEO-AWI-UoR')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'source', 'MAOOAM-pyPDAF')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'history', iso8601()//': Data created')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'reference', 'https://github.com/yumengch/pyPDAF')
   end subroutine setAttrs

   subroutine getVarAttrs(fieldnames, standard_name, long_name, dims)
      integer, intent(out) :: dims(4, 4)
      character(len=5), intent(out) :: fieldnames(4)
      character(len=40), intent(out) :: standard_name(4)
      character(len=50), intent(out) :: long_name(4)

      fieldnames = [character(len=5) :: 'psi_a', 'T_a', 'psi_o', 'T_o']
      standard_name = [character(len=40):: 'atmosphere_streamfunction', &
                       'atmosphere_temperature', &
                       'ocean_streamfunction', &
                       'ocean_temperature' &
                       ]
      long_name = [character(len=50) :: 'streamfunction in the atmosphere', &
                   'temperature in the atmosphere', &
                   'streamfunction in the ocean', &
                   'temperature in the ocean' &
                   ]

      dims = reshape([dimid(2), dimid(3), dimid(4), dimid(1), &
              dimid(2), dimid(3), dimid(4), dimid(1), &
              dimid(2), dimid(3), dimid(4), dimid(1), &
              dimid(2), dimid(3), dimid(4), dimid(1)], shape(dims), order=[2, 1])
   end subroutine getVarAttrs

   function iso8601() result(datetime)
      !! Returns current date and time in ISO 8601 format.
      !! from https://cyber.dabamos.de/programming/modernfortran/date-and-time.html
      character(len=*), parameter :: ISO8601_FMT = '(i4, 2("-", i2.2), "T", 2(i0.2, ":"), i0.2, ".", i0.3, a, ":", a)'
      character(len=29) :: datetime
      character(len=5)  :: zone
      integer           :: dt(8)

      call date_and_time(values=dt, zone=zone)

      write (datetime, ISO8601_FMT) dt(1), dt(2), dt(3), dt(5), dt(6), &
          dt(7), dt(8), zone(1:3), zone(4:5)
   end function iso8601

   subroutine write_state(step, vartype, inputData, nx, ny, dim_ens)
      real(wp), intent(in) :: step
      character, intent(in) :: vartype
      real(wp), intent(in) :: inputData(:, :)
      integer, intent(in)      :: nx, ny, dim_ens

      integer :: i, ierr, i_ens
      integer :: i0,tc

      print *, step, vartype, nx, ny, dim_ens
      if (vartype == 'f') then
         i0 = 0
         time_count = time_count + 1
         tc = time_count
         ierr = nf90_put_var(ncid, varid_time, [abs(step)], start=[tc], count=[1])
      else
         i0 = 4
         tc = time_count
      endif

      do i = 1, 4
         do i_ens = 1, dim_ens
            ierr = nf90_put_var(ncid, varid(i0 + i), &
                                inputData((i-1)*nx*ny + 1: i*nx*ny, i_ens), &
                                        start=[i_ens, 1, 1, tc], &
                                count=[1, nx, ny, 1] &
                                )
         end do
      end do
   end subroutine write_state

   subroutine finalize_state_writer()
      integer :: ierr
      ierr = nf90_close(ncid)
   end subroutine finalize_state_writer
end module mod_StateWriter_pdaf
