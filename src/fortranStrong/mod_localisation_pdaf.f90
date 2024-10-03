module mod_localisation_pdaf
   use mod_kind_pdaf, only: wp
   use PDAFomi, only: PDAFomi_set_domain_limits
implicit none

real(wp) :: cradius
real(wp) :: sradius
integer :: locweight


namelist /local_nml/ cradius, sradius, locweight

contains
   subroutine set_lim_coords()
      use mod_model_pdaf, only: nx, ny, pi, maooam_model

      real(wp) :: lim_coords(2, 2)
      real(wp) :: n

      n = maooam_model%model_configuration%physics%n
      lim_coords(1, 1) = 0
      lim_coords(1, 2) = pi*2/n
      lim_coords(2, 1) = pi
      lim_coords(2, 2) = 0

      call PDAFomi_set_domain_limits(lim_coords)
   end subroutine set_lim_coords
end module mod_localisation_pdaf
