program genCovar
use mod_kind_pdaf, only: wp
use netcdf
use model_def, only: model
implicit none
type(model) :: maooam_model
integer :: ncid
integer :: varid(9)

integer :: nx, ny, nt
integer :: natm, noc, nmod
integer :: dim_state
integer :: dim_field
integer :: status
integer :: i
real(wp), allocatable :: coeffs(:, :)
real(wp), allocatable :: psi_a(:, :)
real(wp), allocatable :: T_a(:, :)
real(wp), allocatable :: psi_o(:, :)
real(wp), allocatable :: T_o(:, :)

real(wp), allocatable :: states(:, :)
real(wp), allocatable :: meanstate(:)
real(wp), allocatable :: svdU(:, :)
real(wp), allocatable :: svals(:)
real(wp), allocatable :: stddev(:)

! timers
integer :: t_rate
integer :: t_reader_start, t_reader_end
integer :: t_transform_start, t_transform_end
integer :: t_gen_start, t_gen_end
integer :: t_write_start, t_write_end
real(wp) :: t_reader_dur, t_transform_dur, t_gen_dur, t_write_dur

namelist /genCov_nml/ nx, ny

! read the namelist
open (20, file='genCov_pdaf.nml')
read(20, nml=genCov_nml)
close(20)
! initialise the model parameters for Fourier/physical state transformation
CALL maooam_model%init
natm = maooam_model%model_configuration%modes%natm
noc = maooam_model%model_configuration%modes%noc
nmod = maooam_model%model_configuration%modes%ndim
call SYSTEM_CLOCK(t_reader_start)
! read model trajectory
call trajectory_reader(natm, noc, nmod, nt, coeffs)
call system_clock(t_reader_end, t_rate)
t_reader_dur = (real(t_reader_end, wp) - real(t_reader_start,wp))/real(t_rate, wp)

dim_state = 4*nx*ny
if (.not. allocated(psi_a)) allocate(psi_a(nx, ny))
if (.not. allocated(T_a)) allocate(T_a(nx, ny))
if (.not. allocated(psi_o)) allocate(psi_o(nx, ny))
if (.not. allocated(T_o)) allocate(T_o(nx, ny))
if (.not. allocated(states))  ALLOCATE(states(dim_state, nt))

call system_clock(t_transform_start)
do i = 1, nt
   ! transform from Fourier to physical space
   call toPhysical(maooam_model, nx, ny, coeffs(:, i), &
                   psi_a, T_a, psi_o, T_o)
   ! transform model field to state vector
   call FieldtoState(nx, ny, psi_a, T_a, psi_o, T_o, states(:, i))
end do
call system_clock(t_transform_end, t_rate)
t_transform_dur = (real(t_transform_end, wp) - real(t_transform_start,wp))/real(t_rate,wp)

! compute the temporal mean of the state vector
if (.not. allocated(meanstate))  ALLOCATE(meanstate(dim_state))
meanstate = sum(states, dim=2)/nt

! Allocate arrays for singular values and vectors
dim_field = nx*ny
ALLOCATE(svals(nt))
ALLOCATE(svdU(dim_state, nt))
allocate(stddev(4))
call system_clock(t_gen_start)
! Call routine generating matrix decomposition
CALL PDAF_eofcovar(dim_state, nt, 4, &
                   [dim_field, dim_field, dim_field, dim_field], &
                   [0, dim_field, 2*dim_field, 3*dim_field], &
                    1, 1, states, stddev, svals, svdU, meanstate, 3, status)
call system_clock(t_gen_end, t_rate)
t_gen_dur = (real(t_gen_end, wp) - real(t_gen_start ,wp))/real(t_rate, wp)

call system_clock(t_write_start)
! save the EOF results
call init_covar_writer('covariance.nc', nt, nx, ny, ncid, varid)
call write_covar(ncid, varid, nt, nx, ny, svals, svdU, meanstate)
call finalize_covar_writer()
call system_clock(t_write_end, t_rate)
t_write_dur = (real(t_write_end, wp) - real(t_write_start,wp))/real(t_rate,wp)

print *, 'reader time', t_reader_dur
print *, 'transform time', t_transform_dur
print *, 'gencov time', t_gen_dur
print *, 'write time', t_write_dur

deallocate(coeffs, psi_a,  T_a,  psi_o,  T_o)
deallocate(states, meanstate, svdU, svals, stddev)

contains
   subroutine trajectory_reader(natm, noc, nmod, nt, coeffs)
      integer, intent(in)                 :: natm, noc, nmod
      real(wp),intent(inout), allocatable :: coeffs(:, :)
      integer, intent(inout)              :: nt

      ! local variables
      integer          :: ncid
      integer          :: varid
      integer          :: ierr
      integer          :: dimid
      integer          :: i
      integer          :: dim(4)
      integer          :: offset(5)
      character(len=7) :: varnames(4)

      varnames = [character(len=7) :: 'psi_a_f', 'T_a_f', 'psi_o_f', 'T_o_f']
      dim = [natm, natm, noc, noc]
      offset = [0, natm, 2*natm, 2*natm+noc, 2*natm+2*noc]

      ierr = nf90_open('trajectory.nc', nf90_nowrite, ncid)
      ierr = nf90_inq_dimid(ncid, 'time', dimid)
      ierr = nf90_inquire_dimension(ncid, dimid, len=nt)

      allocate(coeffs(nmod, nt))
      do i = 1, 4
         ierr = nf90_inq_varid(ncid, trim(varnames(i)), varid)
         ierr = nf90_get_var(ncid, varid, coeffs(offset(i)+1:offset(i+1), :), &
                             start=[1, 1], count=[dim(i), nt]) 
      end do
      ierr = nf90_close(ncid)
   end subroutine trajectory_reader

   real(wp) function Fa(M, H, P, n, nx, ny, x, y)
      integer,  intent(in) :: nx, ny
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: n
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y

      Fa = sqrt(2.)*cos(P*y)
   end function Fa

   real(wp) function Fk(M, H, P, n, nx, ny, x, y)
      integer,  intent(in) :: nx, ny
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: n
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y

      Fk = 2*cos(M*x*n)*sin(P*y)
   end function Fk

   real(wp) function Fl(M, H, P, n, nx, ny, x, y)
      integer,  intent(in) :: nx, ny
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: n
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y

      Fl = 2*sin(H*x*n)*sin(P*y)
   end function Fl

   real(wp) function phi(H, P, n, nx, ny, x, y)
      integer,  intent(in) :: nx, ny
      integer,  intent(in) :: H, P
      real(wp), intent(in) :: n
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y

      phi = 2*sin(0.5*H*x*n)*sin(P*y)
   end function phi

   subroutine toPhysical(maooam_model, nx, ny, coeffs, psi_a, T_a, psi_o, T_o)
      type(model), intent(in) :: maooam_model
      integer,  intent(in) :: nx, ny
      real(wp), intent(in) :: coeffs(:)
      real(wp), intent(inout) :: psi_a(:, :)
      real(wp), intent(inout) :: T_a(:, :)
      real(wp), intent(inout) :: psi_o(:, :)
      real(wp), intent(inout) :: T_o(:, :)

      real(wp), parameter :: pi = dacos(-1.D0)
      integer   :: H, M, P
      real(wp)  :: n
      integer   :: i, j, k
      real(wp)  :: dx, dy
      CHARACTER :: typ=" "

      interface basis
         function basis(M, H, P, n, nx, ny, x, y) result(res)
            import :: wp
            integer,  intent(in) :: nx, ny
            integer,  intent(in) :: M, H, P
            real(wp), intent(in) :: n
            real(wp), intent(in) :: x
            real(wp), intent(in) :: y

            real(wp) :: res
         end function basis
      end interface basis
      procedure(basis), pointer :: f => null()

      ! set up grid
      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)

      ! get atmospheric components
      psi_a = 0.
      T_a = 0.
      do k = 1, natm
         H = maooam_model%inner_products%awavenum(k)%H
         M = maooam_model%inner_products%awavenum(k)%M
         P = maooam_model%inner_products%awavenum(k)%P
         typ = maooam_model%inner_products%awavenum(k)%typ
         if (typ == "A") then
            f => Fa
         else if (typ == "K") then
            f => Fk
         else if (typ == "L") then
            f => Fl
         else
            print *, "error in function type"
            stop
         end if
         do j = 1, ny
            do i = 1, nx
               psi_a(i, j) = psi_a(i, j) + coeffs(k)*f(M, H, P, n, nx, ny, (i-1)*dx, (j-1)*dy)
               T_a(i, j) = T_a(i, j) + coeffs(k + natm)*f(M, H, P, n, nx, ny, (i-1)*dx, (j-1)*dy)
            enddo
         enddo
      end do

      ! get ocean components
      psi_o = 0.
      T_o = 0.
      do k = 1, noc
         H = maooam_model%inner_products%owavenum(k)%H
         P = maooam_model%inner_products%owavenum(k)%P
         do j = 1, ny
            do i = 1, nx
               psi_o(i, j) = psi_o(i, j) + coeffs(k + 2*natm)*phi(H, P, n, nx, ny, (i-1)*dx, (j-1)*dy)
               T_o(i, j) = T_o(i, j) + coeffs(k + 2*natm + noc)*phi(H, P, n, nx, ny, (i-1)*dx, (j-1)*dy)
            enddo
         enddo
      end do
   end subroutine toPhysical

   subroutine init_covar_writer(filename, nrank, nx, ny, ncid, varid)
      character(*), intent(in) :: filename
      integer, intent(in)      :: nrank, nx, ny
      integer, intent(inout)   :: ncid
      integer, intent(inout)   :: varid(9)

      integer :: ierr
      integer :: i, j
      integer :: dimid(3)
      integer :: dimids(4, 3)

      character(len=5)  :: vartype(2)
      character(len=5)  :: fieldnames(4)
      character(len=15) :: std_typename(2)
      character(len=15) :: long_typename(2)
      character(len=40) :: standard_name(4)
      character(len=50) :: long_name(4)

      vartype = [character(len=5)  :: '_svd', '_mean']
      std_typename = [character(len=15) :: 'singular_vector', 'temporal_mean']
      long_typename = [character(len=15) :: 'singular vector', 'temporal mean']

      ierr = nf90_create(filename, nf90_netcdf4, ncid)

      ! set global attributes
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'Conventions', 'CF-1.8')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'title', &
      'NetCDF output for covariance matrix based on model trajectory using PDAF_eofcovar')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'institution', &
      'NCEO-AWI-UoR')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'source', 'MAOOAM-pyPDAF')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'history', iso8601()//': Data created')
      ierr = nf90_put_att(ncid,  NF90_GLOBAL, 'reference', 'https://github.com/yumengch/pyPDAF')

      ! initialise dimensions
      ierr = nf90_def_dim( ncid, 'nx', nx, dimid(1) )
      ierr = nf90_def_dim( ncid, 'ny', ny, dimid(2) )
      ierr = nf90_def_dim( ncid, 'rank', nrank, dimid(3) )
      ! initialise output variables
      ierr = nf90_def_var(ncid, 'sigma', nf90_double, dimid(3), varid(1))
      ierr = nf90_put_att(ncid, varid(1), 'standard_name', 'singular_value' )
      ierr = nf90_put_att(ncid, varid(1), 'long_name', 'singular value of the state vector')

      ! set variable attributes
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

      dimids = reshape([dimid(1), dimid(2),dimid(3), &
              dimid(1), dimid(2),dimid(3), &
              dimid(1), dimid(2),dimid(3), &
              dimid(1), dimid(2),dimid(3)], shape(dimids), order=[2, 1])

      do i = 1, 2
         do j = 1, 4
            ierr = nf90_def_var(ncid, trim(fieldnames(j))//trim(vartype(i)), nf90_double, dimids(j, :4-i), varid((i-1)*4+j+1))
            ierr = nf90_put_att(ncid, varid((i-1)*4+j+1), 'standard_name', trim( trim(std_typename(i))//'_'//standard_name(j)) )
            ierr = nf90_put_att(ncid, varid((i-1)*4+j+1), 'long_name', trim(long_typename(i))//' of '//trim(long_name(j)))
         end do
      end do
      ierr = NF90_ENDDEF(ncid)
   end subroutine init_covar_writer

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

   subroutine write_covar(ncid, varid, nrank, nx, ny, svals, svdU, meanstate)
      integer,   intent(in) :: ncid
      integer,   intent(in) :: varid(9)
      integer,   intent(in) :: nrank, nx, ny
      real(wp),  intent(in) :: svals(nrank)
      real(wp),  intent(in) :: meanstate(ny*nx*4)
      real(wp),  intent(in) :: svdU(ny*nx*4, nrank)

      integer :: ierr, i, j

      ierr = nf90_put_var( ncid, varid(1), svals, start=[1], count=[nrank] )
      do j = 1, 4
         do i = 1, nrank
            ierr = nf90_put_var(ncid, varid(j+1), &
                                reshape(svdU((j-1)*nx*ny + 1: j*nx*ny, i), [nx, ny]), &
                                start=[1, 1, i], &
                                count=[nx, ny, 1] &
                                )
         end do
      end do

      do j = 1, 4
         ierr = nf90_put_var(ncid, varid(j+5), &
                             reshape(meanstate((j-1)*nx*ny + 1: j*nx*ny), [nx, ny]), &
                             start=[1, 1], &
                             count=[nx, ny] &
                             )
      end do
   end subroutine write_covar

   subroutine finalize_covar_writer()
      integer :: ierr
      ierr = nf90_close(ncid)
   end subroutine finalize_covar_writer

   subroutine FieldtoState(nx, ny, psi_a, T_a, psi_o, T_o, states)
      integer, intent(in)  :: nx, ny
      real(wp), intent(in) :: psi_a(nx, ny)
      real(wp), intent(in) :: T_a(nx, ny)
      real(wp), intent(in) :: psi_o(nx, ny)
      real(wp), intent(in) :: T_o(nx, ny)

      real(wp), intent(inout) :: states(4*nx*ny)

      states(:nx*ny) = reshape(psi_a, [nx*ny])
      states(nx*ny+1:2*nx*ny) = reshape(T_a, [nx*ny])
      states(2*nx*ny+1:3*nx*ny) = reshape(psi_o, [nx*ny])
      states(3*nx*ny+1:4*nx*ny) = reshape(T_o, [nx*ny])
   end subroutine FieldtoState
end program genCovar
