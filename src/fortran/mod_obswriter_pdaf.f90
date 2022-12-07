module mod_ObsWriter_pdaf
use mod_kind_pdaf, only: wp
use mod_model_pdaf, only: nx, ny, dim_ens
use mod_observations_pdaf, only: n_obs, nxo, nyo
use netcdf
implicit none

integer, allocatable :: time_count(:)
integer, allocatable :: ncid(:)
integer, allocatable :: dimid(:, :)
integer, allocatable :: varid_time(:)
integer, allocatable :: varid(:, :)

contains
   subroutine init_obs_writer(filenames)
      character(*), intent(in) :: filenames(:)

      integer :: ierr
      integer :: j
      integer :: i_obs
      integer :: dimids(4, 3)

      character(len=5)  :: varname(4)
      character(len=40) :: standard_name(4)
      character(len=50) :: long_name(4)

      allocate(time_count(n_obs))
      allocate(ncid(n_obs), dimid(n_obs, 3), varid_time(n_obs))
      allocate(varid(n_obs, 4))

      time_count = 0
      do i_obs = 1, n_obs
         ierr = nf90_create(trim(filenames(i_obs)), nf90_netcdf4, ncid(i_obs))
         ! set global attributes
         call setAttrs(i_obs)
         ! initialise dimensions
         ierr = nf90_def_dim( ncid(i_obs), 'time', nf90_unlimited, dimid(i_obs, 1) )
         ierr = nf90_def_var( ncid(i_obs), 'time', nf90_double, dimid(i_obs, 1), varid_time(i_obs) )
         ierr = nf90_put_att( ncid(i_obs), varid_time(i_obs), 'long_name', 'time' )
         ierr = nf90_put_att( ncid(i_obs), varid_time(i_obs), 'units', 'days since 1900-1-1 0:0:0' )
         ierr = nf90_def_dim( ncid(i_obs), 'nx', nxo, dimid(i_obs, 2) )
         ierr = nf90_def_dim( ncid(i_obs), 'ny', nyo, dimid(i_obs, 3) )
         ! initialise output variables
         call getVarAttrs(i_obs, varname, standard_name, long_name, dimids)

         do j = 1, 4
            ierr = nf90_def_var(ncid(i_obs), trim(varname(j)), nf90_double, dimids(j, :), varid(i_obs, j))
            ierr = nf90_put_att(ncid(i_obs), varid(i_obs, j), 'standard_name', trim(standard_name(j)))
            ierr = nf90_put_att(ncid(i_obs), varid(i_obs, j), 'long_name', trim(long_name(j)))
         end do
         ierr = NF90_ENDDEF(ncid(i_obs))
      end do
   end subroutine init_obs_writer

   subroutine setAttrs(i_obs)
      integer, intent(in) :: i_obs
      integer :: ierr
      
      ierr = nf90_put_att(ncid(i_obs),  NF90_GLOBAL, 'Conventions', 'CF-1.8')
      ierr = nf90_put_att(ncid(i_obs),  NF90_GLOBAL, 'title', &
      'NetCDF output of synthetic observations from MAOOAM-pyPDAF')
      ierr = nf90_put_att(ncid(i_obs),  NF90_GLOBAL, 'institution', &
      'NCEO-AWI-UoR')
      ierr = nf90_put_att(ncid(i_obs),  NF90_GLOBAL, 'source', 'MAOOAM-pyPDAF')
      ierr = nf90_put_att(ncid(i_obs),  NF90_GLOBAL, 'history', iso8601()//': Data created')
      ierr = nf90_put_att(ncid(i_obs),  NF90_GLOBAL, 'reference', 'https://github.com/yumengch/pyPDAF')
   end subroutine setAttrs

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

   subroutine getVarAttrs(i_obs, fieldnames, standard_name, long_name, dims)
      integer, intent(in)  :: i_obs
      integer, intent(out) :: dims(4, 3)
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

      dims = reshape([dimid(i_obs, 2),dimid(i_obs, 3), dimid(i_obs, 1), &
              dimid(i_obs, 2),dimid(i_obs, 3), dimid(i_obs, 1), &
              dimid(i_obs, 2),dimid(i_obs, 3), dimid(i_obs, 1), &
              dimid(i_obs, 2),dimid(i_obs, 3), dimid(i_obs, 1)], shape(dims), order=[2, 1])
   end subroutine getVarAttrs

   subroutine writeObs(i_obs, step, inputData)
      integer, intent(in) :: i_obs
      integer, intent(in) :: step
      real(wp), intent(in) :: inputData(:)

      integer :: i, ierr

      time_count(i_obs) = time_count(i_obs) + 1
      ierr = nf90_put_var(ncid(i_obs), varid_time(i_obs), &
                          [real(step, wp)], &
                          start=[time_count(i_obs)], count=[1])

      do i = 1, 4
         ierr = nf90_put_var(ncid(i_obs), varid(i_obs, i), &
                             inputData((i-1)*nxo*nyo + 1: i*nxo*nyo), &
                             start=[1, 1, time_count(i_obs)], &
                             count=[nxo, nyo, 1] &
                             )
      end do
   end subroutine writeObs

   subroutine finalizeObs()
      integer :: ierr
      integer :: i_obs
      do i_obs = 1, n_obs
         ierr = nf90_close(ncid(i_obs))
      end do
   end subroutine finalizeObs
end module mod_ObsWriter_pdaf
