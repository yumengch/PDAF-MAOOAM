! This file is part of pyPDAF

! Copyright (C) 2022 University of Reading and
! National Centre for Earth Observation

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module mod_model
use mod_maooam, only: maooam_model, integr
use mod_romb_pdaf, only: romb
implicit none

real(kind=8), parameter   :: pi =  3.14159265358979323846
real(kind=8), ALLOCATABLE :: field_g(:), field_new_g(:)
integer :: natm, noc, ndim

contains
   !> init model
   subroutine initialize_model()
      ! initialise model configurations
      CALL maooam_model%init()
      CALL integr%init(maooam_model)

      ALLOCATE(field_g(0:integr%ndim))
      ALLOCATE(field_new_g(0:integr%ndim))
      field_g = 1.
      field_new_g = 1.

      natm = maooam_model%model_configuration%modes%natm
      noc = maooam_model%model_configuration%modes%noc
      ndim = maooam_model%model_configuration%modes%ndim
   end subroutine initialize_model

   subroutine step(field, t, field_out)
      real(kind=8), intent(in) :: field(:)
      REAL(kind=8), INTENT(inout) :: t
      REAL(kind=8), INTENT(inout) :: field_out(:)

      field_g(1:) = field(:)
      CALL integr%step(field_g, t, field_new_g)
      field_out(:) = field_new_g(1:)
   end subroutine step

   function Fa(P, y) result(c)
      integer,  intent(in) :: P
      real(kind=8), intent(in) :: y

      real(kind=8) :: c

      c = sqrt(2.)*cos(P*y)
   end function Fa

   function Fk(M, P, x, y) result(c)
      integer,  intent(in) :: M, P
      real(kind=8), intent(in) :: x, y

      real(kind=8) :: c
      real(kind=8) :: n

      n = maooam_model%model_configuration%physics%n
      c = 2*cos(M*x*n)*sin(P*y)
   end function Fk

   function Fl(H, P, x, y) result(c)
      integer,  intent(in) :: H, P
      real(kind=8), intent(in) :: x, y

      real(kind=8) :: c
      real(kind=8) :: n

      n = maooam_model%model_configuration%physics%n
      c = 2*sin(H*x*n)*sin(P*y)
   end function Fl

   function phi(H, P, x, y) result(c)
      integer,  intent(in) :: H, P
      real(kind=8), intent(in) :: x, y

      real(kind=8) :: c
      real(kind=8) :: n

      n = maooam_model%model_configuration%physics%n
      c = 2*sin(0.5*H*x*n)*sin(P*y)
   end function phi

   subroutine toPhysical_A(field, psi_a, T_a)
      real(kind=8), intent(in) :: field(:)
      real(kind=8), intent(inout) :: psi_a(:, :), T_a(:, :)

      integer :: H, M, P
      integer :: i, j, k
      integer :: nx, ny
      real(kind=8) :: dx, dy, n
      CHARACTER :: typ=" "

      ! set up grid
      n = maooam_model%model_configuration%physics%n
      nx = size(psi_a(:, 1))
      ny = size(psi_a(1, :))
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
            do j = 1, ny
               do i = 1, nx
                  psi_a(i, j) = psi_a(i, j) + field(k)       *Fa(P, (j-1)*dy)
                  T_a(i, j)   = T_a(i, j)   + field(k + natm)*Fa(P, (j-1)*dy)
               enddo
            enddo
         else if (typ == "K") then
            do j = 1, ny
               do i = 1, nx
                  psi_a(i, j) = psi_a(i, j) + field(k)       *Fk(M, P, (i-1)*dx, (j-1)*dy)
                  T_a(i, j)   = T_a(i, j)   + field(k + natm)*Fk(M, P, (i-1)*dx, (j-1)*dy)
               enddo
            enddo
         else if (typ == "L") then
            do j = 1, ny
               do i = 1, nx
                  psi_a(i, j) = psi_a(i, j) + field(k)       *Fl(H, P, (i-1)*dx, (j-1)*dy)
                  T_a(i, j)   = T_a(i, j)   + field(k + natm)*Fl(H, P, (i-1)*dx, (j-1)*dy)
               enddo
            enddo
         else
            print *, "error in function type"
            stop
         end if
      end do
   end subroutine toPhysical_A


   subroutine toPhysical_O(field, psi_o, T_o)
      real(kind=8), intent(in) :: field(:)
      real(kind=8), intent(inout) :: psi_o(:, :), T_o(:, :)
      integer :: H, P
      integer :: i, j, k
      integer :: nx, ny
      real(kind=8) :: dx, dy, n

      nx = size(psi_o(:, 1))
      ny = size(psi_o(1, :))
      ! set up grid
      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)

      ! get ocean components
      psi_o = 0.
      T_o = 0.
      do k = 1, noc
         H = maooam_model%inner_products%owavenum(k)%H
         P = maooam_model%inner_products%owavenum(k)%P
         do j = 1, ny
            do i = 1, nx
               psi_o(i, j) = psi_o(i, j) + field(k + 2*natm)*phi(H, P, (i-1)*dx, (j-1)*dy)
               T_o(i, j) = T_o(i, j) + field(k + 2*natm + noc)*phi(H, P, (i-1)*dx, (j-1)*dy)
            enddo
         enddo
      end do
   end subroutine toPhysical_O


   subroutine toPhysical_A_multistep(field, psi_a, T_a)
      real(kind=8), intent(in) :: field(:, :)
      real(kind=8), intent(inout) :: psi_a(:, :, :), T_a(:, :, :)

      integer :: H, M, P
      integer :: i, j, k, it
      integer :: nx, ny, nt
      real(kind=8) :: dx, dy, n
      CHARACTER :: typ=" "

      ! set up grid
      n = maooam_model%model_configuration%physics%n
      nt = size(psi_a(:, 1, 1))
      nx = size(psi_a(1, :, 1))
      ny = size(psi_a(1, 1, :))

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
            do j = 1, ny
               do i = 1, nx
                  psi_a(:, i, j) = psi_a(:, i, j) + field(:, k)       *Fa(P, (j-1)*dy)
                  T_a(:, i, j)   = T_a(:, i, j)   + field(:, k + natm)*Fa(P, (j-1)*dy)
               enddo
            enddo
         else if (typ == "K") then
            do j = 1, ny
               do i = 1, nx
                  psi_a(:, i, j) = psi_a(:, i, j) + field(:, k)       *Fk(M, P, (i-1)*dx, (j-1)*dy)
                  T_a(:, i, j)   = T_a(:, i, j)   + field(:, k + natm)*Fk(M, P, (i-1)*dx, (j-1)*dy)
               enddo
            enddo
         else if (typ == "L") then
            do j = 1, ny
               do i = 1, nx
                  psi_a(:, i, j) = psi_a(:, i, j) + field(:, k)       *Fl(H, P, (i-1)*dx, (j-1)*dy)
                  T_a(:, i, j)   = T_a(:, i, j)   + field(:, k + natm)*Fl(H, P, (i-1)*dx, (j-1)*dy)
               enddo
            enddo
         else
            print *, "error in function type"
            stop
         end if
      end do

   end subroutine toPhysical_A_multistep


   subroutine toPhysical_O_multistep(field, psi_o, T_o)
      real(kind=8), intent(in) :: field(:, :)
      real(kind=8), intent(inout) :: psi_o(:, :, :), T_o(:, :, :)
      integer :: H, P
      integer :: i, j, k, it
      integer :: nx, ny, nt
      real(kind=8) :: dx, dy, n

      nt = size(psi_o(:, 1, 1))
      nx = size(psi_o(1, :, 1))
      ny = size(psi_o(1, 1, :))
      ! set up grid
      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)

      ! get ocean components
      psi_o = 0.
      T_o = 0.

      do k = 1, noc
         H = maooam_model%inner_products%owavenum(k)%H
         P = maooam_model%inner_products%owavenum(k)%P
         do j = 1, ny
            do i = 1, nx
               psi_o(:, i, j) = psi_o(:, i, j) + field(:, k + 2*natm)*phi(H, P, (i-1)*dx, (j-1)*dy)
               T_o(:, i, j) = T_o(:, i, j) + field(:, k + 2*natm + noc)*phi(H, P, (i-1)*dx, (j-1)*dy)
            enddo
         enddo
      end do
   end subroutine toPhysical_O_multistep


   subroutine toFourier_A(nx, ny, psi_a, T_a, field)
      integer, intent(in) :: nx, ny
      real(kind=8), intent(in) :: psi_a(:, :), T_a(:, :)
      real(kind=8), intent(inout) :: field(:)
      integer :: H, M, P
      integer :: i, j, k
      integer :: nk(2)
      CHARACTER :: typ=" "
      real(kind=8) :: n
      real(kind=8) :: dx, dy
      real(kind=8) :: integrand(nx, ny)
      real(kind=8) :: integral(nx)

      ! set up grid
      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)

      ! set up numerical integration
      nk(1) = int(log(real(nx, kind=8))/log(2.) + 1)
      nk(2) = int(log(real(ny, kind=8))/log(2.) + 1)

      ! get atmospheric components
      do k = 1, natm
         H = maooam_model%inner_products%awavenum(k)%H
         M = maooam_model%inner_products%awavenum(k)%M
         P = maooam_model%inner_products%awavenum(k)%P
         typ = maooam_model%inner_products%awavenum(k)%typ
         if (typ == "A") then
            do j = 1, ny
               do i = 1, nx
                  integrand(i, j) = psi_a(i, j)*Fa(P, (j-1)*dy)
               enddo
            enddo
            
            do j = 1, nx
               integral(j) = romb(ny, nk(2), integrand(j, :), dy)
            end do
            field(k) = romb(nx, nk(1), integral, dx)

            do j = 1, ny
               do i = 1, nx
                  integrand(i, j) = T_a(i, j)*Fa(P, (j-1)*dy)
               enddo
            enddo
         else if (typ == "K") then
            do j = 1, ny
               do i = 1, nx
                  integrand(i, j) = psi_a(i, j)*Fk(M, P, (i-1)*dx, (j-1)*dy)
               enddo
            enddo
            
            do j = 1, nx
               integral(j) = romb(ny, nk(2), integrand(j, :), dy)
            end do
            field(k) = romb(nx, nk(1), integral, dx)

            do j = 1, ny
               do i = 1, nx
                  integrand(i, j) = T_a(i, j)*Fk(M, P, (i-1)*dx, (j-1)*dy)
               enddo
            enddo
         else if (typ == "L") then
             do j = 1, ny
               do i = 1, nx
                  integrand(i, j) = psi_a(i, j)*Fl(H, P, (i-1)*dx, (j-1)*dy)
               enddo
            enddo
            
            do j = 1, nx
               integral(j) = romb(ny, nk(2), integrand(j, :), dy)
            end do
            field(k) = romb(nx, nk(1), integral, dx)

            do j = 1, ny
               do i = 1, nx
                  integrand(i, j) = T_a(i, j)*Fl(H, P, (i-1)*dx, (j-1)*dy)
               enddo
            enddo
         else
            print *, "error in function type"
            stop
         end if

         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dy)
         end do
         field(k+natm) = romb(nx, nk(1), integral, dx)
      end do

      field(1:2*natm) =  field(1:2*natm)*n/2/pi/pi
   end subroutine toFourier_A


   subroutine toFourier_O(nx, ny, psi_o, T_o, field)
      integer, intent(in) :: nx, ny
      real(kind=8), intent(in) :: psi_o(:, :), T_o(:, :)
      real(kind=8), intent(inout) :: field(:)
      integer :: H, P
      integer :: i, j, k
      integer :: nk(2)
      real(kind=8) :: n
      real(kind=8) :: dx, dy
      real(kind=8) :: integrand(nx, ny)
      real(kind=8) :: integral(nx)

      ! set up grid
      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)

      ! set up numerical integration
      nk(1) = int(log(real(nx, kind=8))/log(2.) + 1)
      nk(2) = int(log(real(ny, kind=8))/log(2.) + 1)

      ! get ocean components
      do k = 1, noc
         H = maooam_model%inner_products%owavenum(k)%H
         P = maooam_model%inner_products%owavenum(k)%P

         do j = 1, ny
            do i = 1, nx
               integrand(i, j) = psi_o(i, j)*phi(H, P, (i-1)*dx, (j-1)*dy)
            enddo
         enddo
         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dy)
         end do
         field(k + 2*natm) = romb(nx, nk(1), integral, dx)

         do j = 1, ny
            do i = 1, nx
               integrand(i, j) = T_o(i, j)*phi(H, P, (i-1)*dx, (j-1)*dy)
            enddo
         enddo

         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dy)
         end do
         field(k + 2*natm + noc) = romb(nx, nk(1), integral, dx)
      end do
      field(2*natm+1:) =  field(2*natm+1:)*n/2/pi/pi
   end subroutine toFourier_O
   
   subroutine getMoving(nt, window, coeffs, &
                      psi_a, T_a, psi_o, T_o, &
                      psi_a_out, T_a_out, psi_o_out, T_o_out)
      integer, intent(in) :: nt, window
      real(kind=8), intent(in) :: coeffs(:, :)
      real(kind=8), intent(inout) :: psi_a(:, :), T_a(:, :)
      real(kind=8), intent(inout) :: psi_o(:, :), T_o(:, :)
      real(kind=8), intent(inout) :: psi_a_out(:, :, :), T_a_out(:, :, :)
      real(kind=8), intent(inout) :: psi_o_out(:, :, :), T_o_out(:, :, :)
      integer :: i, j, nmax, nmin, nstart

      psi_a_out = 0.
      T_a_out = 0.
      psi_o_out = 0.
      T_o_out = 0.
      do i = 1, window
         call toPhysical_A(coeffs(:, i), psi_a, T_a)
         call toPhysical_O(coeffs(:, i), psi_o, T_o)
         nmax = min(i, nt)
         do j = 1, nmax
            psi_a_out(:, :, j) = psi_a_out(:, :, j) + psi_a/window
            T_a_out  (:, :, j) = T_a_out  (:, :, j) + T_a/window
            psi_o_out(:, :, j) = psi_o_out(:, :, j) + psi_o/window
            T_o_out  (:, :, j) = T_o_out  (:, :, j) + T_o/window
         end do
      end do
      do i = window + 1, nt
         call toPhysical_A(coeffs(:, i), psi_a, T_a)
         call toPhysical_O(coeffs(:, i), psi_o, T_o)
         do j = i - window + 1, i
            psi_a_out(:, :, j) = psi_a_out(:, :, j) + psi_a/window
            T_a_out  (:, :, j) = T_a_out  (:, :, j) + T_a/window
            psi_o_out(:, :, j) = psi_o_out(:, :, j) + psi_o/window
            T_o_out  (:, :, j) = T_o_out  (:, :, j) + T_o/window
         end do
      end do

      nstart = max(nt + 1, window + 1)
      do i = nstart, nt + window - 1
         call toPhysical_A(coeffs(:, i), psi_a, T_a)
         call toPhysical_O(coeffs(:, i), psi_o, T_o)
         nmin = max(i - window + 1, 1)
         do j = nmin, nt
            psi_a_out(:, :, j) = psi_a_out(:, :, j) + psi_a/window
            T_a_out  (:, :, j) = T_a_out  (:, :, j) + T_a/window
            psi_o_out(:, :, j) = psi_o_out(:, :, j) + psi_o/window
            T_o_out  (:, :, j) = T_o_out  (:, :, j) + T_o/window
         end do
      end do
   end subroutine getMoving

   subroutine finalize_model
      DEALLOCATE(field_g, field_new_g)
      CALL maooam_model%clean
      CALL integr%clean
   end subroutine finalize_model
end module mod_model
