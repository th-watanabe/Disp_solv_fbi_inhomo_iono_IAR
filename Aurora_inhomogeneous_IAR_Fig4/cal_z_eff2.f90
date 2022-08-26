PROGRAM cal_z_eff

  use re_header

  implicit none

  real(kind=DP), save :: k_perp, sigma_l, alpha_l, recomb, nu_in

  real(kind=DP) :: integral_mp, average_mp, integral_mp_z_tilde_sq

  integer, save :: isw

  complex(kind=DP) :: omega, domega
  integer :: i, j, n=10000

  ! initial values
  data nu_in   / 1._DP  /
  data alpha_l / 1.0_DP  /
  data sigma_l / 5.0_DP  /
  data recomb  / 0.01_DP /
  data isw       / 0     /

    height = 0.03_DP

    domega = ( 0.01_DP, 0._DP )

    nu_in = 0._DP
    do j = 1, 5
      call fbi_set_nu_in( real(j, kind=DP) )
      call fbi_init
      omega = ( 0._DP, 0._DP )
    do i = 1, n
      omega = omega + domega
      call fbi_z_eff_model( omega )
    end do
      write(6,*)
    end do


CONTAINS

  SUBROUTINE fbi_init

!    write(6,fmt="(a,1p,1e15.7)")
    write(6,fmt="(a,1p,1e15.7)") "# k_perp = ", k_perp
    write(6,fmt="(a,1p,1e15.7)") "# height = ", height
    write(6,fmt="(a,1p,1e15.7)") "# nu_in  = ", nu_in
    write(6,fmt="(a,1p,1e15.7)") "# alpha_l= ", alpha_l
    write(6,fmt="(a,1p,1e15.7)") "# sigma_l= ", sigma_l

  END SUBROUTINE


  SUBROUTINE fbi_set_nu_in( nu_in_new )
    real(kind=DP) :: nu_in_new
      nu_in = nu_in_new
  END SUBROUTINE

  SUBROUTINE fbi_set_alpha_l( alpha_l_new )
    real(kind=DP) :: alpha_l_new
      alpha_l = alpha_l_new
  END SUBROUTINE

  SUBROUTINE fbi_set_sigma_l( sigma_l_new )
    real(kind=DP) :: sigma_l_new
      sigma_l = sigma_l_new
  END SUBROUTINE

  SUBROUTINE fbi_set_isw( isw_new )
    integer :: isw_new
      isw = isw_new
  END SUBROUTINE fbi_set_isw


  SUBROUTINE fbi_z_eff_model( omega )

    integer, parameter :: nz = 40000

    complex(kind=DP), intent(in)  :: omega
!    complex(kind=DP), intent(out) :: z_eff
    complex(kind=DP) :: z_eff, z_eff2, beta

    real(kind=DP),    dimension(0:nz) :: zz

    real(kind=DP),    dimension(0:nz) :: mp_z, nu_z, dn_z, ff_z, mp_z_hat, mp_z_tilde, mp_z_tilde_sq

    real(kind=DP)    :: dz

    integer       :: i, j


    integral_mp = 0._DP
    dz          = height / dble( nz )
    zz(0)       = - height
! slab geometry
    do i = 0, nz
      zz(i)     = zz(0) + dz * dble(i)
    end do


    do i = 0, nz
      nu_z(i) = nu_in * ( 1._DP - cos( pi * zz(i) / height ) )
      mp_z(i) = nu_z(i) / ( 1 + nu_z(i)**2 )
      ! ff_z(i) = ( 1._DP - nu_z(i)**2 ) / ( 1._DP + nu_z(i)**2 )**2
      integral_mp = integral_mp + mp_z(i) * dz
!        write(unit=90,fmt="(1p,11e15.7)") zeta(i),mp_z(i), nu_z(i)
    end do

   average_mp = integral_mp / height

   do i = 0, nz
     mp_z_hat(i) = mp_z(i) / average_mp
!     mp_z(i)     = average_mp
!     mp_z_hat(i) = 1._DP

!     mp_z_tilde(i) = ( mp_z(i) - average_mp ) / average_mp
     mp_z_tilde(i) = mp_z_hat(i) - 1._DP

!     if(nu_in == 1.0)then
!        write(99,fmt="(4e15.7)") zz(i), mp_z_hat(i), mp_z_tilde(i), mp_z_tilde_sq(i)
!     endif
   enddo

   z_eff = ( 0._DP, 0._DP )
   z_eff2= ( 0._DP, 0._DP )

   beta  = ( omega + ui*2._DP*recomb ) / sigma_l

   do i = 0, nz
     z_eff = z_eff + mp_z_tilde(i) * ( 1._DP + mp_z_tilde(i) ) &
                   / ( ( 1._DP + mp_z_tilde(i) )*sigma_l - omega - ui*2._DP*recomb )
     z_eff2= z_eff2 + 1._DP / ( mp_z_tilde(i) + 1._DP - beta )
   enddo

     z_eff = z_eff * sigma_l * dz / height

     z_eff2= beta * ( 1._DP - z_eff2 * dz / height * ( 1._DP - beta ) )

     write(6,fmt="(1p,6e15.7)") omega, z_eff, z_eff2


  END SUBROUTINE fbi_z_eff_model


END PROGRAM cal_z_eff
