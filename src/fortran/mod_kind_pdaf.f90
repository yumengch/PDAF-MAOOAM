module mod_kind_pdaf
implicit none
INTEGER, PUBLIC, PARAMETER ::   sp = SELECTED_REAL_KIND( 6, 37)   !< single precision (real 4)
INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !< double precision (real 8)
INTEGER, PUBLIC, PARAMETER ::   wp = dp                           !< working precision
end module mod_kind_pdaf