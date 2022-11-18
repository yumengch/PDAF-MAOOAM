program genCovar
use mod_kind_pdaf, only: wp
use netcdf
use model_def, only: model
use mod_romb_pdaf,only: romb
implicit none
type(model)      :: maooam_model
integer          :: ncid
integer          :: varid
integer          :: ierr
integer          :: dimid
integer          :: dims(4)
integer          :: offset(5)
character(len=7) :: varnames(4)

integer :: ncid_o
integer :: varid_o(4)

integer :: nx, ny, nt
integer :: natm, noc, nmod
integer :: i, it

real(wp) :: dt
real(wp), allocatable :: coeffs(:)
real(wp), allocatable :: psi_a(:, :)
real(wp), allocatable :: T_a(:, :)
real(wp), allocatable :: psi_o(:, :)
real(wp), allocatable :: T_o(:, :)

real(wp), allocatable :: psi_a_std(:, :)
real(wp), allocatable :: T_a_std(:, :)
real(wp), allocatable :: psi_o_std(:, :)
real(wp), allocatable :: T_o_std(:, :)

real(wp), allocatable :: psi_a_mean(:, :)
real(wp), allocatable :: T_a_mean(:, :)
real(wp), allocatable :: psi_o_mean(:, :)
real(wp), allocatable :: T_o_mean(:, :)


! initialise the model parameters for Fourier/physical state
! transformation
CALL maooam_model%init
natm = maooam_model%model_configuration%modes%natm
noc = maooam_model%model_configuration%modes%noc
nmod = maooam_model%model_configuration%modes%ndim

! read model trajectory
if (.not. allocated(coeffs)) allocate(coeffs(nmod))

ierr = nf90_open('trajectory.nc', nf90_nowrite, ncid)
varnames = [character(len=7) :: 'psi_a_f', 'T_a_f', 'psi_o_f', 'T_o_f']
dims = [natm, natm, noc, noc]
offset = [0, natm, 2*natm, 2*natm+noc, 2*natm+2*noc]

ierr = nf90_inq_dimid(ncid, 'time', dimid)
ierr = nf90_inquire_dimension(ncid, dimid, len=nt)
print *, nt
nx = 257
ny = nx
if (.not. allocated(psi_a)) allocate(psi_a(nx, ny))
if (.not. allocated(T_a)  ) allocate(T_a(nx, ny))
if (.not. allocated(psi_o)) allocate(psi_o(nx, ny))
if (.not. allocated(T_o)  ) allocate(T_o(nx, ny))

if (.not. allocated(psi_a_std)) allocate(psi_a_std(nx, ny))
if (.not. allocated(T_a_std)  ) allocate(T_a_std(nx, ny))
if (.not. allocated(psi_o_std)) allocate(psi_o_std(nx, ny))
if (.not. allocated(T_o_std)  ) allocate(T_o_std(nx, ny))

if (.not. allocated(psi_a_mean)) allocate(psi_a_mean(nx, ny))
if (.not. allocated(T_a_mean)  ) allocate(T_a_mean(nx, ny))
if (.not. allocated(psi_o_mean)) allocate(psi_o_mean(nx, ny))
if (.not. allocated(T_o_mean)  ) allocate(T_o_mean(nx, ny))

psi_a_mean = 0._wp
T_a_mean = 0._wp
psi_o_mean = 0._wp
T_o_mean = 0._wp
dt = 1._wp/nt
do it = 1, nt 
   do i = 1, 4
      ierr = nf90_inq_varid(ncid, trim(varnames(i)), varid)
      ierr = nf90_get_var(ncid, varid, coeffs(offset(i)+1:offset(i+1)),&
                          start=[1, it], count=[dims(i), 1]) 
   end do
   ! transform from Fourier to physical space
   call toPhysical(maooam_model, nx, ny, coeffs, &
                   psi_a, T_a, psi_o, T_o)
   psi_a_mean = psi_a_mean + psi_a*dt
   T_a_mean = T_a_mean + T_a*dt
   psi_o_mean = psi_o_mean + psi_o*dt
   T_o_mean = T_o_mean + T_o*dt
end do

print *, minval(psi_a_mean), maxval(psi_a_mean)

psi_a_std = 0.
T_a_std = 0.
psi_o_std = 0.
T_o_std = 0.
do it = 1, nt
   do i = 1, 4
      ierr = nf90_inq_varid(ncid, trim(varnames(i)), varid)
      ierr = nf90_get_var(ncid, varid, coeffs(offset(i)+1:offset(i+1)),&
                          start=[1, it], count=[dims(i), 1]) 
   end do
   ! transform from Fourier to physical space
   call toPhysical(maooam_model, nx, ny, coeffs, &
                   psi_a, T_a, psi_o, T_o)
   psi_a_std = psi_a_std + dt*(psi_a - psi_a_mean)*(psi_a - psi_a_mean) 
   T_a_std = T_a_std + dt*(T_a - T_a_mean)*(T_a - T_a_mean) 
   psi_o_std = psi_o_std + dt*(psi_o - psi_o_mean)*(psi_o - psi_o_mean) 
   T_o_std = T_o_std + dt*(T_o - T_o_mean)*(T_o - T_o_mean) 
end do

print *, minval(psi_a_std), maxval(psi_a_std)
deallocate(psi_a, T_a, psi_o, T_o)
deallocate(psi_a_mean, T_a_mean, psi_o_mean, T_o_mean)

ierr = nf90_close(ncid)

call init_writer('traj_var.nc', nx, ny, ncid_o, varid_o)
call writestd(ncid_o, varid_o, nx, ny, psi_a_std, T_a_std, psi_o_std, T_o_std)
call finalize_writer(ncid_o)

deallocate(coeffs)

contains
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

   subroutine toFourier(maooam_model, nx, ny, coeffs, psi_a, T_a, psi_o, T_o)
      type(model), intent(in)    :: maooam_model
      integer,     intent(in)    :: nx, ny
      real(wp),    intent(inout) :: coeffs(:)
      real(wp),    intent(in)    :: psi_a(:, :)
      real(wp),    intent(in)    :: T_a(:, :)
      real(wp),    intent(in)    :: psi_o(:, :)
      real(wp),    intent(in)    :: T_o(:, :)

      real(wp), parameter :: pi = dacos(-1.D0)
      integer :: H, M, P
      integer :: i, j, k
      integer :: nk(2)
      CHARACTER :: typ=" "
      real(wp) :: n
      real(wp) :: dx, dy
      real(wp) :: integrand(nx, ny)
      real(wp) :: integral(nx)
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

      ! set up numerical integration
      nk(1) = int(log(real(nx, wp))/log(2.) + 1)
      nk(2) = int(log(real(ny, wp))/log(2.) + 1)

      ! get atmospheric components
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
               integrand(i, j) = psi_a(i, j)*f(M, H, P, n, nx, ny, (i-1)*dx, (j-1)*dy)
            enddo
         enddo
         
         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dy)
         end do
         coeffs(k) = romb(nx, nk(1), integral, dx)

         do j = 1, ny
            do i = 1, nx
               integrand(i, j) = T_a(i, j)*f(M, H, P, n, nx, ny, (i-1)*dx, (j-1)*dy)
            enddo
         enddo

         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dy)
         end do
         coeffs(k+natm) = romb(nx, nk(1), integral, dx)
      end do

      ! get ocean components
      do k = 1, noc
         H = maooam_model%inner_products%owavenum(k)%H
         P = maooam_model%inner_products%owavenum(k)%P

         do j = 1, ny
            do i = 1, nx
               integrand(i, j) = psi_o(i, j)*phi(H, P, n, nx, ny, (i-1)*dx, (j-1)*dy)
            enddo
         enddo
         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dy)
         end do
         coeffs(k + 2*natm) = romb(nx, nk(1), integral, dx)

         do j = 1, ny
            do i = 1, nx
               integrand(i, j) = T_o(i, j)*phi(H, P, n, nx, ny, (i-1)*dx, (j-1)*dy)
            enddo
         enddo

         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dy)
         end do
         coeffs(k + 2*natm + noc) = romb(nx, nk(1), integral, dx)
      end do

      coeffs(1:) =  coeffs(1:)*n/2/pi/pi
   end subroutine toFourier

   function iso8601() result(datetime)
      !! Returns current date and time in ISO 8601 format.
      !! from
      !https://cyber.dabamos.de/programming/modernfortran/date-and-time.html
      character(len=*), parameter :: ISO8601_FMT = '(i4, 2("-", i2.2), "T", 2(i0.2, ":"), i0.2, ".", i0.3, a, ":", a)'
      character(len=29) :: datetime
      character(len=5)  :: zone
      integer           :: dt(8)

      call date_and_time(values=dt, zone=zone)

      write (datetime, ISO8601_FMT) dt(1), dt(2), dt(3), dt(5), dt(6), &
          dt(7), dt(8), zone(1:3), zone(4:5)
   end function iso8601

   subroutine init_writer(filename, nx, ny, ncid_o, varid_o)
      character(*), intent(in) :: filename
      integer, intent(in)      :: nx, ny
      integer, intent(inout)   :: ncid_o
      integer, intent(inout)   :: varid_o(4)

      integer :: ierr
      integer :: i, j
      integer :: dimid(2)
      integer :: dimids(4, 2)

      character(len=5)  :: fieldnames(4)
      character(len=40) :: standard_name(4)
      character(len=50) :: long_name(4)


      ierr = nf90_create(filename, nf90_netcdf4, ncid_o)

      ! set global attributes
      ierr = nf90_put_att(ncid_o,  NF90_GLOBAL, 'Conventions', 'CF-1.8')
      ierr = nf90_put_att(ncid_o,  NF90_GLOBAL, 'title', &
      'NetCDF output for covariance matrix based on model trajectory using PDAF_eofcovar')
      ierr = nf90_put_att(ncid_o,  NF90_GLOBAL, 'institution', &
      'NCEO-AWI-UoR')
      ierr = nf90_put_att(ncid_o,  NF90_GLOBAL, 'source', 'MAOOAM-pyPDAF')
      ierr = nf90_put_att(ncid_o,  NF90_GLOBAL, 'history', iso8601()//': Data created')
      ierr = nf90_put_att(ncid_o,  NF90_GLOBAL, 'reference', 'https://github.com/yumengch/pyPDAF')

      ! initialise dimensions
      ierr = nf90_def_dim( ncid_o, 'nx', nx, dimid(1) )
      ierr = nf90_def_dim( ncid_o, 'ny', ny, dimid(2) )

      ! initialise output variables
      ! set variable attributes
      fieldnames = [character(len=5) :: 'psi_a_var', 'T_a_var', 'psi_o_var', 'T_o_var']
      standard_name = [character(len=40):: 'variance_atmosphere_streamfunction', &
                       'variance_atmosphere_temperature', &
                       'variance_ocean_streamfunction', &
                       'variance_ocean_temperature' &
                       ]
      long_name = [character(len=50) :: 'variance of streamfunction in the atmosphere', &
                   'variance of temperature in the atmosphere', &
                   'variance of streamfunction in the ocean', &
                   'variance of temperature in the ocean' &
                   ]

      dimids = reshape([dimid(1), dimid(2), &
              dimid(1), dimid(2), &
              dimid(1), dimid(2), &
              dimid(1), dimid(2)], shape(dimids), order=[2, 1])

      do j = 1, 4
         ierr = nf90_def_var(ncid_o, trim(fieldnames(j)), nf90_double, dimids(j, :), varid_o(j))
         ierr = nf90_put_att(ncid_o, varid_o(j), 'standard_name', trim( standard_name(j) ) )
         ierr = nf90_put_att(ncid_o, varid_o(j), 'long_name', trim(long_name(j)))
      end do
      ierr = NF90_ENDDEF(ncid_o)
   end subroutine init_writer

   subroutine writestd(ncid_o, varid_o, nx, ny, psi_a_std, T_a_std, psi_o_std, T_o_std)
      integer,   intent(in) :: ncid_o
      integer,   intent(in) :: varid_o(4)
      integer,   intent(in) :: nx, ny
      real(wp),  intent(in) :: psi_a_std(:, :)
      real(wp),  intent(in) :: T_a_std(:, :)
      real(wp),  intent(in) :: psi_o_std(:, :)
      real(wp),  intent(in) :: T_o_std(:, :)


      integer :: ierr

      ierr = nf90_put_var(ncid_o, varid_o(1), &
                          psi_a_std, &
                          start=[1, 1], &
                          count=[nx, ny] &
                          )
      ierr = nf90_put_var(ncid_o, varid_o(2), &
                          T_a_std, &
                          start=[1, 1], &
                          count=[nx, ny] &
                          )
      ierr = nf90_put_var(ncid_o, varid_o(3), &
                          psi_o_std, &
                          start=[1, 1], &
                          count=[nx, ny] &
                          )
      ierr = nf90_put_var(ncid_o, varid_o(4), &
                          T_o_std, &
                          start=[1, 1], &
                          count=[nx, ny] &
                          )
   end subroutine writestd

   subroutine finalize_writer(ncid_o)
      integer, intent(in) :: ncid_o
      integer :: ierr
      ierr = nf90_close(ncid_o)
   end subroutine finalize_writer
end program genCovar

