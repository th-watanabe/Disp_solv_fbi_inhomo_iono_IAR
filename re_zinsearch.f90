PROGRAM zinsearch

      use re_header
      use re_fbi
      use re_mag

      implicit none

      integer, save :: series
      data series / 0 /

      complex(kind=DP) :: fun_search
      complex(kind=DP) :: somega
      complex(kind=DP) :: ey_m, bx_m, bx_bottom
      integer :: re, ima, i
      ! real(kind=DP),   parameter :: remax=5._DP / 20._DP, imamax=0.5_DP / 20._DP
      ! real(kind=DP),   parameter :: remin=0._DP / 20._DP, imamin=-0.5_DP / 20._DP
      real(kind=DP),   parameter :: remax=5._DP, imamax=0.5_DP
      real(kind=DP),   parameter :: remin=0._DP, imamin=-0.5_DP

!      real(kind=DP),   parameter :: remax=0.14_DP , imamax=0.0005_DP
!      real(kind=DP),   parameter :: remin=0.12_DP , imamin=0._DP
      real(kind=DP)              :: dre, dima
      integer      ,   parameter :: nre = 200  , nima = 200  , nalpha_l = 10
      real(kind=DP),   parameter :: malpha = 500._DP
      real(kind=DP)              :: k_perp, alpha_l, dalpha_l
!      real(kind=DP),   parameter :: k_perp = 1._DP, alpha_l = 5._DP

      if (series .eq. 1 ) then

      dalpha_l = malpha / real( nalpha_l,kind=DP )

      dre  = ( remax  - remin  ) / dble( nre )
      dima = ( imamax - imamin ) / dble( nima )

      do i   = 0, nalpha_l
      alpha_l = dalpha_l * dble( i )
       if( i.eq.0)then
       alpha_l = 0.0001
       endif
      call fbi_set_alpha_l(alpha_l)
      do re  = 0, nre-1
      do ima = 0, nima-1

         somega = cmplx( remin+dble(re)*dre, imamin+dble(ima)*dima,kind=DP )

         call mag_response( somega, bx_m, ey_m )
!         call fbi_set_phi_jpara( phi_jpara )
         call fbi_func_inhomo( somega, bx_m, ey_m, bx_bottom )

!         fun_search = somega - omega

!         write (*,*) real( somega), aimag(somega), cabs( cmplx(fun_search,kind=4))
         write(*,*) alpha_l, real(somega), aimag(somega), cabs( cmplx(bx_bottom,kind=4) )

      enddo
      write(*,*)
      enddo
!      write(*,*)
      enddo

    elseif(series .eq. 0)then

      dre  = ( remax  - remin  ) / dble( nre )
      dima = ( imamax - imamin ) / dble( nima )

      alpha_l = malpha

      call fbi_set_alpha_l(alpha_l)

      do re  = 0, nre-1
      do ima = 0, nima-1

         somega = cmplx( remin+dble(re)*dre, imamin+dble(ima)*dima,kind=DP )

         call mag_response( somega, bx_m, ey_m )
!         call fbi_set_phi_jpara( phi_jpara )
         call fbi_func_inhomo( somega, bx_m, ey_m, bx_bottom )

!         fun_search = somega - omega

!         write (*,*) real( somega), aimag(somega), cabs( cmplx(fun_search,kind=4))
         write(6,fmt="(4e15.7)") alpha_l, real(somega), aimag(somega), cabs( cmplx(bx_bottom) )

      enddo
      write(*,*)
      enddo

    endif

END PROGRAM
