module mod_romb_pdaf
use mod_kind_pdaf, only: wp
implicit none

contains
   function romb(nx, nk, y, dx) result(integral)
      integer,  intent(in)  :: nx, nk
      real(wp), intent(in)  :: y(0:nx-1)
      real(wp), intent(in)  :: dx
      real(wp)              :: integral

      integer :: Nsamps, Ninterv
      integer :: start, step, limit
      integer :: n
      integer :: k
      integer :: i, j
      real(wp) :: h
      real(wp) :: prev
      real(wp) :: R(0:nk, 0:nk)


      Nsamps = nx
      Ninterv = Nsamps-1
      n = 1
      k = 0
      do while (n < Ninterv)
          n = n * 2
          k = k + 1
      end do

      if (n /= Ninterv) then
         print *, "Number of samples must be one plus a non-negative power of 2."
         stop
      end if

      h = Ninterv * dx
      R(0, 0) = (y(0) + y(nx-1))/2.0*h

      start = Ninterv
      limit = Ninterv
      step = Ninterv
      do i = 1, k
         start = start/2
         R(i, 0) = 0.5*(R(i-1, 0) + h*sum(y(start:limit:step)))
         step = step/2
         do j = 1, i
            prev = R(i, j-1)
            R(i, j) = prev + (prev-R(i-1, j-1)) / ((2**(2*j))-1)
         end do
         h = h / 2.0
      end do
      integral = R(k, k)
   end function romb
end module mod_romb_pdaf