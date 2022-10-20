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

integer :: nx, ny, nt
integer :: natm, noc, nmod
integer :: i, it, i_exp, i_sample
integer :: n_exp, n_sample
real(wp), allocatable :: time(:)
real(wp), allocatable :: error(:)
real(wp), allocatable :: coeffs(:)
real(wp), allocatable :: coeffs_trans(:)
real(wp), allocatable :: psi_a(:, :)
real(wp), allocatable :: T_a(:, :)
real(wp), allocatable :: psi_o(:, :)
real(wp), allocatable :: T_o(:, :)

namelist /grid_nml/ n_exp, nt, n_sample

! read the namelist
open (20, file='compare.nml')
read(20, nml=grid_nml)
close(20)

print *, 'n_exp', n_exp, 'nt', nt, 'n_sample', n_sample

! initialise the model parameters for Fourier/physical state transformation
CALL maooam_model%init
natm = maooam_model%model_configuration%modes%natm
noc = maooam_model%model_configuration%modes%noc
nmod = maooam_model%model_configuration%modes%ndim

! read model trajectory
if (.not. allocated(coeffs)) allocate(coeffs(nmod))
if (.not. allocated(coeffs_trans)) allocate(coeffs_trans(nmod))

ierr = nf90_open('trajectory.nc', nf90_nowrite, ncid)
varnames = [character(len=7) :: 'psi_a_f', 'T_a_f', 'psi_o_f', 'T_o_f']
dims = [natm, natm, noc, noc]
offset = [0, natm, 2*natm, 2*natm+noc, 2*natm+2*noc]

if (.not. allocated(time)) allocate(time(nt))
if (.not. allocated(error)) allocate(error(n_exp))

error = 0._wp
do i_sample = 1, n_sample
   call random_seed(size=nt)
   call random_seed(put=[(i + 20000*i_sample,  i = 1, nt)])
   call random_number(time)
   do i_exp = 1, n_exp 
      nx = 2**i_exp + 1
      ny = nx
      if (.not. allocated(psi_a)) allocate(psi_a(nx, ny))
      if (.not. allocated(T_a)) allocate(T_a(nx, ny))
      if (.not. allocated(psi_o)) allocate(psi_o(nx, ny))
      if (.not. allocated(T_o)) allocate(T_o(nx, ny))
      
      do it = 1, nt
         do i = 1, 4
            ierr = nf90_inq_varid(ncid, trim(varnames(i)), varid)
            ierr = nf90_get_var(ncid, varid, coeffs(offset(i)+1:offset(i+1)), &
                                start=[1, 1+floor(20224719*time(it))], count=[dims(i), 1]) 
         end do
         ! transform from Fourier to physical space
         call toPhysical(maooam_model, nx, ny, coeffs, &
                         psi_a, T_a, psi_o, T_o)
         ! transform model field to state vector
         call toFourier(maooam_model, nx, ny, coeffs_trans, &
            psi_a, T_a, psi_o, T_o)
         error(i_exp) = error(i_exp) + maxval(abs(coeffs - coeffs_trans))
      end do
      deallocate(psi_a, T_a, psi_o, T_o)
   end do
end do
error = error/nt/n_sample
print *, error
deallocate(time, error)
ierr = nf90_close(ncid)

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
