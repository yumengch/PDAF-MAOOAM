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

module mod_model_pdaf
use mod_kind_pdaf, only: wp
use model_def, only: model
use rk4_integrator, only: RK4Integrator
use mod_romb_pdaf, only: romb
use mod_ModelWriter_pdaf, only: init_model_writer, finalize_model_writer
implicit none

! Dimension of state vector and ensemble size
integer  :: dim_ens
integer  :: dim_state
integer  :: dim_state_p
integer  :: noc
integer  :: natm

integer  :: total_steps
real(wp) :: total_time
real(wp), allocatable :: field(:, :)
real(wp), allocatable :: field_new(:, :)

real(wp), allocatable :: psi_a(:, :)
real(wp), allocatable :: T_a(:, :)

real(wp), allocatable :: psi_o(:, :)
real(wp), allocatable :: T_o(:, :)

logical :: writeout
real(wp) :: tw
type(Model), TARGET :: maooam_model
type(RK4Integrator) :: integr

interface basis
   function basis(M, H, P, nx, ny, x, y) result(res)
      import :: wp
      integer,  intent(in) :: nx, ny
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: x(nx, ny)
      real(wp), intent(in) :: y(nx, ny)

      real(wp) :: res(nx, ny)
   end function basis
end interface basis

contains
   subroutine initialize_model()

      print *, 'Model MAOOAM v1.4'
      print *, 'Loading information...'
      ! initialise model configurations 
      CALL maooam_model%init

      natm = maooam_model%model_configuration%modes%natm
      noc = maooam_model%model_configuration%modes%noc
      dim_state_p = maooam_model%model_configuration%modes%ndim

      total_time = maooam_model%model_configuration%integration%t_run
      dim_state = dim_state_p
      total_steps = int(total_time/maooam_model%model_configuration%integration%dt)

      writeout = maooam_model%model_configuration%integration%writeout
      tw = maooam_model%model_configuration%integration%tw

      CALL integr%init(maooam_model)

      ! initialise the model writer
      call init_model_writer('maooam.nc', natm, noc, dim_ens)

      ! initialise initial condition
      ALLOCATE(field(0:dim_state_p, 1),field_new(0:dim_state_p, 1))
      field(0:, 1) = maooam_model%load_IC()

   end subroutine initialize_model

   function Fa(M, H, P, nx, ny, x, y) result(res)
      integer,  intent(in) :: nx, ny
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: x(nx, ny)
      real(wp), intent(in) :: y(nx, ny)

      real(wp) :: res(nx, ny)

      res = sqrt(2.)*cos(P*y)
   end function Fa

   function Fk(M, H, P, nx, ny, x, y) result(res)
      integer,  intent(in) :: nx, ny
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: x(nx, ny)
      real(wp), intent(in) :: y(nx, ny)

      real(wp) :: res(nx, ny)

      integer :: n
      n = maooam_model%model_configuration%physics%n

      res = 2*cos(M*x*n)*sin(P*y)
   end function Fk

   function Fl(M, H, P, nx, ny, x, y) result(res)
      integer,  intent(in) :: nx, ny
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: x(nx, ny)
      real(wp), intent(in) :: y(nx, ny)

      real(wp) :: res(nx, ny)

      integer :: n
      n = maooam_model%model_configuration%physics%n

      res = 2*sin(H*x*n)*sin(P*y)
   end function Fl

   function phi(H, P, nx, ny, x, y) result(res)
      integer,  intent(in) :: nx, ny
      integer,  intent(in) :: H, P
      real(wp), intent(in) :: x(nx, ny)
      real(wp), intent(in) :: y(nx, ny)

      real(wp) :: res(nx, ny)

      integer :: n
      n = maooam_model%model_configuration%physics%n

      res = 2*sin(0.5*H*x*n)*sin(P*y)
   end function phi

   subroutine toPhysical(nx, ny, xc, yc)
      integer,  intent(in) :: nx, ny
      real(wp), intent(in) :: xc(nx, ny), yc(nx, ny)
      integer :: H, M, P
      integer :: i
      CHARACTER :: typ=" "
      procedure(basis), pointer :: f => null()

      ! get atmospheric components
      if (.not. allocated(psi_a)) allocate(psi_a(nx, ny))
      if (.not. allocated(T_a)) allocate(T_a(nx, ny))

      psi_a = 0.
      T_a = 0.
      do i = 1, natm
         H = maooam_model%inner_products%awavenum(i)%H
         M = maooam_model%inner_products%awavenum(i)%M
         P = maooam_model%inner_products%awavenum(i)%P
         typ = maooam_model%inner_products%awavenum(i)%typ
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
         psi_a = psi_a + field(i, 1)*f(M, H, P, nx, ny, xc, yc)
         T_a = T_a + field(i + natm, 1)*f(M, H, P, nx, ny, xc, yc)
      end do

      ! get ocean components
      if (.not. allocated(psi_o)) allocate(psi_o(nx, ny))
      if (.not. allocated(T_o)) allocate(T_o(nx, ny))
      psi_o = 0.
      T_o = 0.
      do i = 1, noc
         H = maooam_model%inner_products%owavenum(i)%H
         P = maooam_model%inner_products%owavenum(i)%P
         psi_o = psi_o + field(i + 2*natm, 1)*phi(H, P, nx, ny, xc, yc)
         T_o = T_o + field(i + 2*natm + noc, 1)*phi(H, P, nx, ny, xc, yc)
      end do
   end subroutine toPhysical


   subroutine toFourier(nx, ny, xc, yc)
      integer,  intent(in) :: nx, ny
      real(wp), intent(in) :: xc(nx, ny), yc(nx, ny)
      integer :: H, M, P, n
      integer :: i, j
      integer :: nk(2)
      real(wp) :: pi =  3.14159265358979323846
      real(wp) :: dx(2)
      CHARACTER :: typ=" "
      real(wp) :: integrand(nx, ny)
      real(wp) :: integral(nx)
      procedure(basis), pointer :: f => null()

      if (.not. allocated(field)) allocate(field(0:dim_state_p, 1))

      ! set up numerical integration
      nk(1) = log(real(nx, wp))/log(2.) + 1
      nk(2) = log(real(ny, wp))/log(2.) + 1
      n = maooam_model%model_configuration%physics%n
      dx(1) = 2*pi/n/nx
      dx(2) = pi/ny

      ! get atmospheric components
      do i = 1, natm
         H = maooam_model%inner_products%awavenum(i)%H
         M = maooam_model%inner_products%awavenum(i)%M
         P = maooam_model%inner_products%awavenum(i)%P
         typ = maooam_model%inner_products%awavenum(i)%typ
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

         integrand = psi_a*f(M, H, P, nx, ny, xc, yc)
         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dx(2))
         end do
         field(i, 1) = romb(nx, nk(1), integral, dx(1))

         integrand = T_a*f(M, H, P, nx, ny, xc, yc)
         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dx(2))
         end do
         field(i+natm, 1) = romb(nx, nk(1), integral, dx(1))
      end do

      ! get ocean components
      do i = 1, noc
         H = maooam_model%inner_products%owavenum(i)%H
         P = maooam_model%inner_products%owavenum(i)%P

         integrand = psi_o*phi(H, P, nx, ny, xc, yc)
         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dx(2))
         end do
         field(i + 2*natm, 1) = romb(nx, nk(1), integral, dx(1))

         integrand = T_o*phi(H, P, nx, ny, xc, yc)
         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dx(2))
         end do
         field(i + 2*natm + noc, 1) = romb(nx, nk(1), integral, dx(1))
      end do
   end subroutine toFourier


   subroutine finalize_model

      DEALLOCATE(field, field_new)
      CALL maooam_model%clean
      CALL integr%clean
      call finalize_model_writer()

   end subroutine finalize_model
end module mod_model_pdaf