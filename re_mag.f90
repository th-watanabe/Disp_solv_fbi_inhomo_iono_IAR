MODULE re_mag

  use re_header

  implicit none

  integer,          parameter :: nz = 1000
  integer, save :: isw


data isw / 0 /

private

public mag_response, mag_init, mag_set_isw

CONTAINS

  SUBROUTINE mag_init

      write(olog,fmt="(a,1p,1e15.7)")
      write(olog,fmt="(a,1p,1e15.7)") "# zl     = ", ell

  END SUBROUTINE mag_init

  SUBROUTINE mag_set_isw( isw_new )
    integer :: isw_new
      isw = isw_new
  END SUBROUTINE mag_set_isw


SUBROUTINE mag_response( z, bx_m, ey_m )

  complex(kind=DP), intent(in)  :: z
  complex(kind=DP), intent(out) :: ey_m, bx_m

  complex(kind=DP)                  :: ey, bx

  complex(kind=DP), dimension(0:nz) :: zey, zbx
  real(kind=DP),    dimension(0:nz) :: zz, valp

  real(kind=DP)    :: dz

  complex(kind=DP), dimension(4) :: bx_k, bx_w, ey_k, ey_w
  real(kind=DP),    dimension(4) :: ck, ci

  real(kind=DP)    :: dz2, fc1, fc2, fc3

  integer       :: i, j

  dz   = ell / dble( nz )

  do i = 0, nz
    zz(i) = dz * dble(i)
    valp(i) = sqrt(1._DP / ( epsilon**2 + exp(- 2._DP * zz(i) ) ) )
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

          zbx(:)   = ( 0._DP, 0._DP )
          zey(:)   = ( 0._DP, 0._DP )


!!!! anti-symmetric for phi
          bx = ( 1._DP, 0._DP )
          ey = ( 0._DP, 0._DP )

          zbx( nz )   = bx
          zey( nz )   = ey

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

            if( j == 1 ) then
              fc1  = - 1._DP / valp(i  )**2
            else if ( j == 4 ) then
              fc1  = - 1._DP / valp(i-2)**2
            else
              fc1  = - 1._DP / valp(i-1)**2
            end if

            fc2 = 1._DP
            fc3 = 0._DP

            bx_k(j) = ui * ey_w(j) * fc1 * z
            ey_k(j) = - ui * bx_w(j) * fc2 * z
          end do

          do j = 1, 4
            bx = bx + dz2 * ck(j) * bx_k(j)
            ey = ey + dz2 * ck(j) * ey_k(j)
          end do

          zbx(i-2) = bx
          zey(i-2) = ey

        end do

      bx_m = bx
      ey_m = ey



      if( isw .gt. 0 ) then
          write(unit=20+isw,fmt="(a, 1p, 3e15.7)") "# z = ", z
        do i = 0, nz, 2
          write(unit=20+isw,fmt="(1p, 6e15.7)") zz(i), zbx(i), zey(i), valp(i)
        end do
          write(20+isw,*)
          isw = 0
      end if

   return

END SUBROUTINE

END MODULE re_mag
