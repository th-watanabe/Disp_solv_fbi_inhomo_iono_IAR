MODULE re_mi_couple

  use re_header
  use re_slv
  use re_fbi
  use re_mag

  implicit none
  private
  public mi_couple_find, mi_couple_func

CONTAINS


SUBROUTINE mi_couple_find

  complex(kind=DP) :: zin, zout, omega
  real(kind=DP)    :: alpha_l, sigma_l, nu_in
  complex(kind=DP) :: bx_bottom, bx_m, ey_m

  real(kind=DP)    :: dalpha_l, dnu_in, dheight
  integer          :: i, j, k, m, aa, ss, nn, hh, clabel

  integer, parameter :: jloop=5
  integer, parameter :: nalpha_l = 1000
  integer, parameter :: nnu_in   = 50
!  integer, parameter :: nheight  = 30
  integer, parameter :: nheight  = 1

  complex(kind=DP)  :: zsave1, zsave2

  dalpha_l = 5._DP  / real( nalpha_l,kind=DP )
  dnu_in   = 5._DP / real( nnu_in,kind=DP )
  dheight  = 0.03_DP / real( nheight,kind=DP )
!  dheight  = 1.d-6 / real( nheight,kind=DP )

  sigma_l  = 5.0_DP
  ! nu_in    = 5.0_DP
  ! height   = 0.03_DP

  ! call fbi_set_nu_in( nu_in )
  call fbi_set_sigma_l( sigma_l )
  ! call fbi_init
  ! call mag_init
  ! call mi_couple_init

  call mi_couple_file_open

  do j = 0, jloop

    if(j .eq. 0) then
            ! zin=cmplx(0.15_DP, 0.001_DP,kind=DP )
            zin=cmplx(0.15_DP, 0._DP,kind=DP )
    elseif(j .eq. 1) then
            ! zin=cmplx(1.0_DP, 0.001_DP,kind=DP )
            zin=cmplx(1.0_DP, 0._DP,kind=DP )
            call mi_couple_file_close(1)
    elseif(j .eq. 2) then
            ! zin=cmplx(1.95_DP, 0.05_DP,kind=DP )
            zin=cmplx(1.95_DP, 0._DP,kind=DP )
            call mi_couple_file_close(2)
    elseif(j .eq. 3) then
            ! zin=cmplx(2.59_DP, 0.05_DP,kind=DP )
            zin=cmplx(2.59_DP, 0._DP,kind=DP )
            call mi_couple_file_close(3)
    elseif(j .eq. 4) then
            zin=cmplx(3.3_DP, 0._DP,kind=DP )
            call mi_couple_file_close(4)
    elseif(j .eq. 5) then
            zin=cmplx( 4.3_DP, 0._DP,kind=DP )
            call mi_couple_file_close(5)
    end if

    zsave1 = zin
    zsave2 = zin

    do hh = nheight, 1, -1
    ! do hh = nheight, nheight
      height = dheight * dble(hh)
      zin = zsave2

!      do nn = 0, nnu_in-10
       do nn = 0, 0
        nu_in = 1.0 + dnu_in * dble(nn)
        call fbi_set_nu_in( nu_in )
        zin = zsave1

        do aa = nalpha_l, 10, -1
         alpha_l = dalpha_l * real( aa,kind=DP )
         call fbi_set_alpha_l( alpha_l )

         call slv_newton( mi_couple_func, zin, zout )

         omega  = zout

         ! write(olog,fmt="(1p,5e15.7)") omega, alpha_l, sigma_l, nu_in
         select case(j)
         case(0)
           write(71, fmt="(6e15.7)") omega, alpha_l, sigma_l, nu_in, height
         case(1)
           write(72, fmt="(6e15.7)") omega, alpha_l, sigma_l, nu_in, height
         case(2)l
           write(73, fmt="(6e15.7)") omega, alpha_l, sigma_l, nu_in, height
         case(3)
           write(74, fmt="(6e15.7)") omega, alpha_l, sigma_l, nu_in, height
         case(4)
           write(75, fmt="(6e15.7)") omega, alpha_l, sigma_l, nu_in, height
         case(5)
           write(76, fmt="(6e15.7)") omega, alpha_l, sigma_l, nu_in, height
         end select

         zin = zout

         if( aa .eq. nalpha_l )then
           zsave1 = zout
           if( nn .eq. 1)then
             zsave2 = zout
           end if
         endif

         ! integer, parameter :: nalpha_l = 1000
         ! integer, parameter :: nnu_in   = 100
         ! integer, parameter :: nheight  = 100

         if( aa == nalpha_l .or. aa == nalpha_l*2/5 .or. aa == nalpha_l/5 .or. aa == nalpha_l/10 ) then
           if( nn == 0 .or. nn == 10 .or.nn == 20 .or. nn == 30 .or. nn == 40 )then
             if( hh == nheight .or. hh == 10 .or. hh == 5 .or. hh == 1 )then
               call fbi_set_isw( 1 )
               call mag_set_isw( 1 )
               call mag_response( omega, bx_m, ey_m )
               call fbi_func_inhomo( omega, bx_m, ey_m, bx_bottom )
               call fbi_set_isw( 0 )
               call mag_set_isw( 0 )
             endif
           endif
         endif

         enddo !aa
         ! write(olog,fmt="(1p,7e15.7)")
         select case(j)
         case(0)
           write(71, fmt="(1e15.7)")
         case(1)
           write(72, fmt="(1e15.7)")
         case(2)
           write(73, fmt="(1e15.7)")
         case(3)
           write(74, fmt="(1e15.7)")
         case(4)
           write(75, fmt="(1e15.7)")
         case(5)
           write(76, fmt="(1e15.7)")
         end select

      enddo !nn
    enddo !hh
  enddo ! j

END SUBROUTINE mi_couple_find

SUBROUTINE mi_couple_file_open
  open( 71, file="result_sigma5.0_1st_mode.dat")
  open( 72, file="result_sigma5.0_2nd_mode.dat")
  open( 73, file="result_sigma5.0_3rd_mode.dat")
  open( 74, file="result_sigma5.0_4th_mode.dat")
  open( 75, file="result_sigma5.0_5th_mode.dat")
  open( 76, file="result_sigma5.0_6th_mode.dat")
END SUBROUTINE mi_couple_file_open


SUBROUTINE mi_couple_file_close( clabel )
  integer, intent(in) :: clabel

  select case(clabel)
  case(1)
    close(71)
  case(2)
    close(72)
  case(3)
    close(73)
  case(4)
    close(74)
  case(5)
    close(75)
  case(6)
    close(76)
  end select

END SUBROUTINE mi_couple_file_close


  SUBROUTINE mi_couple_init

      ! write(olog,fmt="(a,1p,e15.7)")

  END SUBROUTINE mi_couple_init


  FUNCTION mi_couple_func( zin )

    complex(kind=DP) :: mi_couple_func
    complex(kind=DP), intent(in) :: zin

    complex(kind=DP) :: omega_m

    complex(kind=DP) :: bx_bottom, bx_m, ey_m

      omega_m = zin
      call mag_response( omega_m, bx_m, ey_m )

!! for inhomogeneous ionosphere
      call fbi_func_inhomo( omega_m, bx_m, ey_m, bx_bottom )

      mi_couple_func = bx_bottom

  END FUNCTION mi_couple_func


END MODULE re_mi_couple
