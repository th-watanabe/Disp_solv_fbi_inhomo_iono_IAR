MODULE re_header

  implicit none

  integer, parameter :: DP = selected_real_kind(14)
  integer, parameter :: QP = selected_real_kind(28)

  integer :: olog = 6

  real(kind=DP) ,parameter :: epsilon = 0.1_DP
  real(kind=DP) ,parameter :: epsilon_i = 1._DP / epsilon

  real(kind=DP) ,parameter :: ell    = 30._DP   !磁気圏の長さ
  real(kind=DP)            :: height            !電離圏の長さ

  complex(kind=DP), parameter :: ui = ( 0._DP, 1._DP )
  real(kind=DP),    parameter :: pi = 3.141592653589793_DP


  public

END MODULE re_header
