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
use mod_parallel_pdaf, only: task_id, abort_parallel
use mod_ModelWriter_pdaf, only: init_model_writer, finalize_model_writer
implicit none

! Dimension of state vector and ensemble size
integer  :: dim_ens
integer  :: dim_state
integer  :: dim_state_p
integer  :: noc
integer  :: natm

integer  :: total_steps
real(wp) :: total_time, current_time(1)
integer :: nx, ny

real(wp), parameter   :: pi =  3.14159265358979323846
real(wp), allocatable :: field(:)
real(wp), allocatable :: field_new(:)

real(wp), allocatable :: psi_a(:, :)
real(wp), allocatable :: T_a(:, :)

real(wp), allocatable :: psi_o(:, :)
real(wp), allocatable :: T_o(:, :)

logical :: writeout
logical :: ln_restart
logical :: is_freerun
integer :: restart_it
integer :: tw
type(Model), TARGET :: maooam_model
type(RK4Integrator) :: integr

namelist /model_nml/ nx, ny, ln_restart, restart_it, is_freerun, tw

interface basis
   real(wp) function basis(M, H, P, x, y)
      import :: wp
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y
   end function basis
end interface basis

contains
   subroutine initialize_model()
      integer :: ndim
      character(len=3)  :: task_id_str

      print *, 'Model MAOOAM v1.4'
      print *, 'Loading information...'
      ! initialise model configurations 
      CALL maooam_model%init
      CALL integr%init(maooam_model)

      natm = maooam_model%model_configuration%modes%natm
      noc = maooam_model%model_configuration%modes%noc
      ndim = maooam_model%model_configuration%modes%ndim
      dim_state_p = nx*ny*4

      total_time = maooam_model%model_configuration%integration%t_run
      dim_state = dim_state_p

      writeout = maooam_model%model_configuration%integration%writeout
      ! tw = maooam_model%model_configuration%integration%tw

      ! initialise the model writer
      if (writeout) &
         call init_model_writer(natm, noc, dim_ens)
         
      ! initialise initial condition
      ALLOCATE(field(0:ndim),field_new(0:ndim))
      field = maooam_model%load_IC()
      
      if (ln_restart) then
         write(task_id_str, '(I3.3)') task_id
         call read_restart('restart/maooam_'//trim(task_id_str)//'.nc', natm, noc, field(1:), restart_it)
      else
         restart_it = 1
         current_time(1) = 0
      endif

      total_steps = int(total_time/maooam_model%model_configuration%integration%dt)
      print *, 'nt', total_steps

      ! get atmospheric components
      if (.not. allocated(psi_a)) allocate(psi_a(nx, ny))
      if (.not. allocated(T_a)) allocate(T_a(nx, ny))
      if (.not. allocated(psi_o)) allocate(psi_o(nx, ny))
      if (.not. allocated(T_o)) allocate(T_o(nx, ny))
   end subroutine initialize_model

   subroutine read_restart(filename, natm, noc, fields, it)
      use netcdf
      character(*), intent(in)    :: filename
      integer,      intent(in)    :: natm, noc
      real(wp),     intent(inout) :: fields(:)
      integer,      intent(inout) :: it

      ! local variables
      integer          :: ncid
      integer          :: varid
      integer          :: dimid
      integer          :: i, nt
      integer          :: ierr
      integer          :: dim(4)
      integer          :: offset(5)
      character(len=7) :: varnames(4)

      varnames = [character(len=7) :: 'psi_a_f', 'T_a_f', 'psi_o_f', 'T_o_f']
      dim = [natm, natm, noc, noc]
      offset = [0, natm, 2*natm, 2*natm+noc, 2*natm+2*noc]

      ierr = nf90_open(filename, nf90_nowrite, ncid)
      ierr = nf90_inq_dimid(ncid, 'time', dimid)
      ierr = nf90_inquire_dimension(ncid, dimid, len=nt)

      if (it > nt) then
         print *, 'The selected restart time step is larger than the time steps in the restart file.'
         call abort_parallel()
      end if

      ierr = nf90_inq_varid(ncid, 'time', varid)
      ierr = nf90_get_var(ncid, varid, current_time, &
                          start=[it], count=[1])

      do i = 1, 4
         ierr = nf90_inq_varid(ncid, trim(varnames(i)), varid)
         ierr = nf90_get_var(ncid, varid, fields(offset(i)+1:offset(i+1)), &
                             start=[1, it], count=[dim(i), 1]) 
      end do
      ierr = nf90_close(ncid)

   end subroutine read_restart

   real(wp) function Fa(M, H, P, x, y)
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y

      Fa = sqrt(2.)*cos(P*y)
   end function Fa

   real(wp) function Fk(M, H, P, x, y)
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y

      real(wp) :: n
      n = maooam_model%model_configuration%physics%n

      Fk = 2*cos(M*x*n)*sin(P*y)
   end function Fk

   real(wp) function Fl(M, H, P, x, y)
      integer,  intent(in) :: M, H, P
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y

      real(wp) :: n
      n = maooam_model%model_configuration%physics%n

      Fl = 2*sin(H*x*n)*sin(P*y)
   end function Fl

   real(wp) function phi(H, P, x, y)
      integer,  intent(in) :: H, P
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y

      real(wp) :: n
      n = maooam_model%model_configuration%physics%n

      phi = 2*sin(0.5*H*x*n)*sin(P*y)
   end function phi

   subroutine toPhysical()

      integer :: H, M, P
      integer :: i, j, k
      real(wp) :: dx, dy, n
      CHARACTER :: typ=" "
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
               psi_a(i, j) = psi_a(i, j) + field(k)*f(M, H, P, (i-1)*dx, (j-1)*dy)
               T_a(i, j) = T_a(i, j) + field(k + natm)*f(M, H, P, (i-1)*dx, (j-1)*dy)
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
               psi_o(i, j) = psi_o(i, j) + field(k + 2*natm)*phi(H, P, (i-1)*dx, (j-1)*dy)
               T_o(i, j) = T_o(i, j) + field(k + 2*natm + noc)*phi(H, P, (i-1)*dx, (j-1)*dy)
            enddo
         enddo
      end do
   end subroutine toPhysical


   subroutine toFourier(nx, ny)
      integer, intent(in) :: nx, ny
      integer :: H, M, P
      integer :: i, j, k
      integer :: nk(2)
      CHARACTER :: typ=" "
      real(wp) :: n
      real(wp) :: dx, dy
      real(wp) :: integrand(nx, ny)
      real(wp) :: integral(nx)
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
               integrand(i, j) = psi_a(i, j)*f(M, H, P, (i-1)*dx, (j-1)*dy)
            enddo
         enddo
         
         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dy)
         end do
         field(k) = romb(nx, nk(1), integral, dx)

         do j = 1, ny
            do i = 1, nx
               integrand(i, j) = T_a(i, j)*f(M, H, P, (i-1)*dx, (j-1)*dy)
            enddo
         enddo

         do j = 1, nx
            integral(j) = romb(ny, nk(2), integrand(j, :), dy)
         end do
         field(k+natm) = romb(nx, nk(1), integral, dx)
      end do

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

      field(1:) =  field(1:)*n/2/pi/pi
   end subroutine toFourier


   subroutine finalize_model
      if (allocated(psi_a)) deallocate(psi_a)
      if (allocated(T_a)) deallocate(T_a)
      if (allocated(psi_o)) deallocate(psi_o)
      if (allocated(T_o)) deallocate(T_o)
      DEALLOCATE(field, field_new)
      CALL maooam_model%clean
      CALL integr%clean
      if (writeout) &
         call finalize_model_writer()

   end subroutine finalize_model
end module mod_model_pdaf
