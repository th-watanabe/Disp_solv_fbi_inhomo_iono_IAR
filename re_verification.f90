MODULE re_verify

  use re_header

  implicit none
  integer, parameter :: nz = 4000

  private

  public verify

contains

SUBROUTINE verify( zz, z, zbx, zey, valp, mp_z, nu_z, mp_z_hat, ff_z, alpha_l, sigma_l, recomb, isw )

  real(kind=DP), dimension(0:nz), intent(in)    :: zz, mp_z, nu_z, mp_z_hat, ff_z, valp
  real(kind=DP), intent(in)                     :: alpha_l, sigma_l, recomb
  complex(kind=DP), dimension(0:nz), intent(in) :: zbx, zey
  complex(kind=DP), intent(in) :: z
  integer         , intent(in)                  :: isw


  complex(kind=DP), dimension(0:nz) :: eq1_l, eq1_r, eq2_l, eq2_r, eq3_l, eq3_r
  integer :: i, j

!
!
! !--- parameters
!
!   complex(kind=DP), dimension(0:nz) :: zbx, zey
!   real(kind=DP),    dimension(0:nz) :: zz, valp
!
!   real(kind=DP),    dimension(0:nz) :: mp_z, nu_z, dn_z, ff_z, al_z, et_z, mp_z_hat, sigma_l, alpha_l, recomb
!
!   complex(kind=DP):: z
!   real(kind=DP)   :: x
  ! character(15)   :: dummy





    ! read(unit=30,fmt="(1p, 2e15.7)") z
    !   do i = 0, nz-1
    !     read(unit=30,fmt="(1p, 13e15.7)") &
    !     zz(i), zbx(i), zey(i), valp(i), mp_z(i), nu_z(i), mp_z_hat(i), ff_z(i), alpha_l, sigma_l, recomb
    !   end do
    !   ! zeta(i), zpsi(i), zphi(i), zdns(i), valp(i), mp_z(i), nu_z(i), dn_z(i), ff_z(i), al_z(i)

   ! open(30)
   write(unit=30+isw,fmt="(1p, 2e15.7)") z
  do i = 0, nz, 2
   write(unit=30+isw,fmt="(1p, 14e15.7)") &
     zz(i), zbx(i), zey(i), valp(i), mp_z(i), nu_z(i), mp_z_hat(i), ff_z(i), alpha_l, sigma_l, recomb, height
  end do
   write(30+isw,*)
   ! close(30)

    eq1_l(:) = ( 0._DP, 0._DP)
    eq1_r(:) = ( 0._DP, 0._DP)
    eq2_l(:) = ( 0._DP, 0._DP)
    eq2_r(:) = ( 0._DP, 0._DP)
    eq3_l(:) = ( 0._DP, 0._DP)
    eq3_r(:) = ( 0._DP, 0._DP)


      ! do i = 2, nz-2, 2
      !   eq1_l(i) = ( zbx(i+2) - zbx(i-2) ) / ( zz(i+2) - zz(i-2) )
      !   ! eq1_r(i) = (-ui * z * ff_z(i) / valp(i)**2 + epsilon_i**2 * alpha_l * mp_z_hat(i) ) * zey(i) &
      !   !           +( ui * z * ff_z(i) / valp(i)**2 * sigma_l * mp_z_hat(i) &
      !   !            - epsilon_i**2 * alpha_l * sigma_l * mp_z_hat(i)**2 ) &
      !   !            / ( z - sigma_l * mp_z_hat(i) + 2._DP * ui * recomb )
      !
      !   eq1_r(i) = ( - ui * z * ff_z(i) / valp(i)**2 &
      !               + epsilon_i**2 * alpha_l * mp_z_hat(i)) &
      !            / ( 1._DP - sigma_l * mp_z_hat(i) / ( z + 2._DP * ui * recomb ) ) * zey(i)
      !
      !
      !   eq2_l(i) = ( zey(i+2) - zey(i-2) ) / ( zz(i+2) - zz(i-2) )
      !   eq2_r(i) = - ui * z * zbx(i)
      !   ! eq3_l(i) = ( - ui * z + 2._DP * al_z(i) ) *  zdns(i)
      !   ! eq3_r(i) =-( zpsi(i+1) - zpsi(i-1) ) / ( zeta(i+1) - zeta(i-1) ) * k_perp**2 / dn_z(i)
      ! end do
      !
      !   open(40)
      !   write(unit=40,fmt="(1p, 2e15.7)") z
      !   do i = 2, nz-2, 2
      !     write(unit=40,fmt="(1p, 13e15.7)") &
      !       zz(i), eq1_l(i), eq1_r(i), eq2_l(i), eq2_r(i), eq3_l(i), eq3_r(i)
      !   enddo
      !   close(40)

    return

END SUBROUTINE verify

END MODULE re_verify
