module mod_ModelWriter_pdaf
use mod_kind_pdaf, only: wp
use netcdf
implicit none

integer :: time_count = 0
integer :: ncid
integer :: dimid(4)
integer :: varid_time
integer :: varid(8)

contains
   subroutine init_model_writer(filename, natm, noc, dim_ens)
      character(*), intent(in) :: filename
      integer, intent(in)      :: natm, noc, dim_ens

      integer :: ierr
      integer :: i, j
      integer :: dimids(4, 3)

      character         :: vartype(2) = ['f', 'a']
      character(len=8)  :: typename(2) = ['forecast', 'analysis']
      character(len=5)  :: varname(4)
      character(len=40) :: standard_name(4)
      character(len=50) :: long_name(4)

      ierr = nf90_create(filename, nf90_netcdf4, ncid)
      ! set global attributes
      call setAttrs
      ! initialise dimensions
      ierr = nf90_def_dim( ncid, 'time', nf90_unlimited, dimid(1) )
      ierr = nf90_def_var( ncid, 'time', nf90_double, dimid(1), varid_time )
      ierr = nf90_put_att( ncid, varid_time, 'long_name', 'time' )
      ierr = nf90_put_att( ncid, varid_time, 'units', 'days since 1900-1-1 0:0:0' )
      ierr = nf90_def_dim( ncid, 'ens', dim_ens, dimid(2) )
      ierr = nf90_def_dim( ncid, 'natm', natm, dimid(3) )
      ierr = nf90_def_dim( ncid, 'noc', noc, dimid(4) )
      ! initialise output variables
      call getVarAttrs(varname, standard_name, long_name, dimids)
      do i = 1, 2
         do j = 1, 4
            ierr = nf90_def_var(ncid, trim(varname(j))//'_'//vartype(i), nf90_double, dimids(j, :), varid((i-1)*4+j))
            ierr = nf90_put_att(ncid, varid((i-1)*4+j), 'standard_name', trim(standard_name(j))//'_'//trim(typename(i)))
            ierr = nf90_put_att(ncid, varid((i-1)*4+j), 'long_name', trim(typename(i))//' of '//trim(long_name(j)))
         end do
      end do
      ierr = NF90_ENDDEF(ncid)
   end subroutine init_model_writer

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
      integer, intent(out) :: dims(4, 3)
      character(len=5), intent(out) :: fieldnames(4)
      character(len=40), intent(out) :: standard_name(4)
      character(len=50), intent(out) :: long_name(4)

      fieldnames = [character(len=5) :: 'psi_a', 'T_a', 'psi_o', 'T_o']
      standard_name = [character(len=40):: 'atmosphere_streamfunction_coefficient', &
                       'atmosphere_temperature_coefficient', &
                       'ocean_streamfunction_coefficient', &
                       'ocean_temperature_coefficient' &
                       ]
      long_name = [character(len=50) :: 'coefficient of streamfunction in the atmosphere', &
                   'coefficient of temperature in the atmosphere', &
                   'coefficient of streamfunction in the ocean', &
                   'coefficient of temperature in the ocean' &
                   ]

      dims = reshape([dimid(2), dimid(3),dimid(1), &
              dimid(2), dimid(3),dimid(1), &
              dimid(2), dimid(4),dimid(1), &
              dimid(2), dimid(4),dimid(1)], shape(dims), order=[2, 1])
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

   subroutine write_model(step, vartype, inputData, natm, noc, dim_ens)
      integer, intent(in) :: step
      character, intent(in) :: vartype
      real(wp), intent(in) :: inputData(:, :)
      integer, intent(in)      :: natm, noc, dim_ens

      integer :: offsets(5)
      integer :: dims(4)
      integer :: i, ierr
      integer :: i0

      time_count = time_count + 1
      ierr = nf90_put_var(ncid, varid_time, [real(step, wp)], start=[time_count], count=[1])

      if (vartype == 'f') then
         i0 = 0
      else
         i0 = 4
      endif
      offsets = [0, natm, 2*natm, 2*natm + noc, 2*(natm + noc)]
      dims = [natm, natm, noc, noc]

      do i = 1, 4
         ierr = nf90_put_var(ncid, varid(i0 + i), &
                             inputData(offsets(i)+1:offsets(i+1), :), &
                             start=[1, 1, time_count], &
                             count=[dim_ens, dims(i), 1] &
                             )
      end do
   end subroutine write_model

   subroutine finalize_model_writer()
      integer :: ierr
      ierr = nf90_close(ncid)
   end subroutine finalize_model_writer
end module mod_ModelWriter_pdaf