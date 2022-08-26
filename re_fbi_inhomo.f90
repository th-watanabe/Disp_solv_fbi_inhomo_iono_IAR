MODULE re_fbi

  use re_header
  use re_verify

  implicit none

  real(kind=DP), save :: k_perp, sigma_l, alpha_l, recomb, nu_in

  real(kind=DP) :: integral_mp, average_mp

  integer, save :: isw

  ! initial values
  data nu_in   / 5._DP  /
  data alpha_l / 5.0_DP  /
  data sigma_l / 5.0_DP  /
!  data recomb  / 0.0_DP  / !mi_coupleにも書いてある
  data recomb  / 1.d-2  / !mi_coupleにも書いてある
  data isw       / 0     /

  private

  public fbi_init, fbi_set_isw, fbi_set_alpha_l, fbi_set_nu_in
  public fbi_set_sigma_l, fbi_func_inhomo, fbi_func

CONTAINS

  SUBROUTINE fbi_init

    write(olog,fmt="(a,1p,1e15.7)")
    write(olog,fmt="(a,1p,1e15.7)") "# k_perp = ", k_perp
    write(olog,fmt="(a,1p,1e15.7)") "# height = ", height
    write(olog,fmt="(a,1p,1e15.7)") "# nu_in  = ", nu_in
    write(olog,fmt="(a,1p,1e15.7)") "# alpha_l= ", alpha_l
    write(olog,fmt="(a,1p,1e15.7)") "# sigma_l= ", sigma_l

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

  SUBROUTINE fbi_func( omega, bx_m, ey_m )
    complex(kind=DP), intent(out) :: omega
    complex(kind=DP), intent(in)  :: bx_m, ey_m

! THW 2021.09.03
!    omega = sigma_l / ( 1._DP + alpha_l * ey_m / bx_m ) - 2._DP * ui * recomb
    omega = sigma_l / ( 1._DP - alpha_l * ey_m / bx_m ) - 2._DP * ui * recomb
! THW 2021.09.03


  END SUBROUTINE fbi_func



  SUBROUTINE fbi_func_inhomo( z, bx_m, ey_m, bx_bottom )

    integer, parameter :: nz = 4000

    complex(kind=DP), intent(in)  :: z, bx_m, ey_m
    complex(kind=DP), intent(out) :: bx_bottom

    complex(kind=DP)                  :: ey, bx

    complex(kind=DP), dimension(0:nz) :: zbx, zey
    real(kind=DP),    dimension(0:nz) :: zz, valp

    real(kind=DP),    dimension(0:nz) :: mp_z, nu_z, dn_z, ff_z, mp_z_hat

    real(kind=DP)    :: dz

    complex(kind=DP), dimension(4) :: bx_k, bx_w, ey_k, ey_w
    real(kind=DP),    dimension(4) :: ck, ci

    real(kind=DP)    :: dz2
    complex(kind=DP) :: fc1, fc2

    integer       :: i, j


    integral_mp = 0
    dz          = height / dble( nz )
    zz(0)       = - height
! slab geometry
    do i = 0, nz
      zz(i)     = zz(0) + dz * dble(i)
      valp(i)     = sqrt( 1._DP / (epsilon**2 + exp(- 2._DP * zz(i) ) ) )
    end do


    do i = 0, nz
      nu_z(i) = nu_in * ( 1._DP - cos( pi * zz(i) / height ) )
      mp_z(i) = nu_z(i) / ( 1 + nu_z(i)**2 )
      ! ff_z(i) = ( 1._DP - nu_z(i)**2 ) / ( 1._DP + nu_z(i)**2 )**2
      ff_z(i) = 0._DP
      integral_mp = integral_mp + mp_z(i) * dz
!        write(unit=90,fmt="(1p,11e15.7)") zeta(i),mp_z(i), nu_z(i)
    end do
   average_mp = integral_mp / height

   do i = 0, nz
     ! mp_z_hat(i) = mp_z(i) / average_mp
     mp_z(i)     = average_mp
     mp_z_hat(i) = 1._DP
   enddo


! coefficients used in the RK method
         ck(1)   = 1._DP/6._DP
         ck(2)   = 1._DP/3._DP
         ck(3)   = 1._DP/3._DP
         ck(4)   = 1._DP/6._DP
         ci(1)   = 0._DP
         ci(2)   = 0.5_DP
         ci(3)   = 0.5_DP
         ci(4)   = 1._DP

       bx_k = ( 0._DP, 0._DP )
       bx_w = ( 0._DP, 0._DP )
       ey_k = ( 0._DP, 0._DP )
       ey_w = ( 0._DP, 0._DP )

! boundary codition given by the magnetosphere
!        bph   = ( 1._DP, 0._DP )
!        psi   = - bph / phi_jpara
       bx      = bx_m
       ey      = ey_m
       zbx(:)  = ( 0._DP, 0._DP )
       zey(:)  = ( 0._DP, 0._DP )


       zbx(nz) = bx
       zey(nz) = ey

! THW 2021.09.03
!       dz2 = dz * 2._DP
       dz2 =-dz * 2._DP
! THW 2021.09.03

       do i = nz, 2, -2
! Runge-Kutta starts
         do j = 1, 4
           if( j .eq. 1 ) then
             bx_w(j) = bx
             ey_w(j) = ey
           else
             bx_w(j) = bx + ci(j) * bx_k(j-1) * dz2
             ey_w(j) = ey + ci(j) * ey_k(j-1) * dz2
           end if

           ! if( j == 1 ) then
           !
           !   fc1 = ( - ui * z * ff_z(i  ) / valp(i  )**2 &
           !        + epsilon_i**2 * alpha_l * mp_z_hat(i  )) &
           !        / ( 1._DP - sigma_l * mp_z_hat(i  ) / ( z + 2._DP * ui * recomb ) )
           !
           ! else if ( j == 4 ) then
           !   fc1 = ( - ui * z * ff_z(i-2) / valp(i-2)**2 &
           !        + epsilon_i**2 * alpha_l * mp_z_hat(i-2)) &
           !        / ( 1._DP - sigma_l * mp_z_hat(i-2) / ( z + 2._DP * ui * recomb ) )
           !
           ! else
           !   fc1 = ( - ui * z * ff_z(i-1) / valp(i-1)**2 &
           !        + epsilon_i**2 * alpha_l * mp_z_hat(i-1)) &
           !        / ( 1._DP - sigma_l * mp_z_hat(i-1) / ( z + 2._DP * ui * recomb ) )
           !
           ! endif

           if( j == 1 ) then

             fc1 = ( - ui * z * ff_z(i  ) / valp(i  )**2 &
                  + alpha_l * mp_z_hat(i  ) / height ) &
                  / ( 1._DP - sigma_l * mp_z_hat(i  ) / ( z + 2._DP * ui * recomb ) )

           else if ( j == 4 ) then
             fc1 = ( - ui * z * ff_z(i-2) / valp(i-2)**2 &
                  + alpha_l * mp_z_hat(i-2) / height ) &
                  / ( 1._DP - sigma_l * mp_z_hat(i-2) / ( z + 2._DP * ui * recomb ) )

           else
             fc1 = ( - ui * z * ff_z(i-1) / valp(i-1)**2 &
                  + alpha_l * mp_z_hat(i-1) / height ) &
                  / ( 1._DP - sigma_l * mp_z_hat(i-1) / ( z + 2._DP * ui * recomb ) )

           endif


           fc2  = - ui * z

           bx_k(j) = ey_w(j) *   fc1
           ey_k(j) = bx_w(j) *   fc2
         end do

         do j = 1, 4
           bx   = bx + dz2 * ck(j) * bx_k(j)
           ey   = ey + dz2 * ck(j) * ey_k(j)
         end do

         zbx(i-2) = bx
         zey(i-2) = ey
!          zdns(i-2) = bph / b0(i-2) * fc1 * ui / z
          ! zdns(i-2) = bph / b0(i-2) * fc1 * ui / ( z + ui*2._DP*recomb )

       end do

       bx_bottom = zbx(0)

       if( isw .gt. 0 ) then
         call  verify( zz, z, zbx, zey, valp, mp_z, nu_z, mp_z_hat, ff_z, alpha_l, sigma_l, recomb, isw )
         isw = 0
       endif


       return

END SUBROUTINE fbi_func_inhomo


END MODULE re_fbi
