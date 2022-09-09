function kron_delta(i, j) result(kd)

    implicit none
    integer, intent(in) :: i, j
    integer :: kd

    kd = 0
    if (i == j) kd = 1

end function kron_delta

function step_func(x) result(hev)

    implicit none
    integer, intent(in) :: x
    integer :: hev

    hev = 0
    if (x > 0) hev = 1

end function step_func

subroutine gaunt(l1, l2, l3, m1, m2, m3, gnt) ! One-Way Schulten-Gordon-Cruzan algorithm to calculate the integral of three spherical harmonics.

    implicit none
    integer, parameter :: dp=kind(0.d0)
    integer, intent(in) :: l1, l2, l3, m1, m2, m3
    integer :: l3t, l3t1, i
    integer(kind=16) :: fact(0:l1+l2+l3)
    real(dp) :: pi=4*atan(1.0_dp), A, A1, B, K, K1, gnt1, gnt2, gnt3, gnt3_temp
    real(dp), intent(out) :: gnt

    ! Apply the angular momentum selection rules to speedup the algorithm.
    gnt = 0.0_dp
    if ((l1 < abs(m1)) .or. (l2 < abs(m2)) .or. (l3 < abs(m3))) then
        gnt = 0.0_dp
    else if (mod(l1+l2+l3, 2) /= 0) then
        gnt = 0.0_dp
    else if (l2+l1 < l3) then
        gnt = 0.0_dp
    else if (l1-l2+l3 < 0) then
        gnt = 0.0_dp
    else if (-l1+l2+l3 < 0) then
        gnt = 0.0_dp
    else if (m1+m2+m3 /= 0) then
        gnt = 0.0_dp
    else

        ! Compute all factorials needed beforehand and store the results in a list.
        fact(0) = 1
        do i=1, l1+l2+l3
            fact(i) = fact(i-1)*i
        end do

        gnt1 = (-1)**(l2-l1)*sqrt((fact(l1+l2)*gamma(l1+0.5_dp)*gamma(l2+0.5_dp))/&
                (sqrt(pi)*fact(l1)*fact(l2)*(2*(l1+l2)+1)*gamma(l1+l2+0.5_dp)))
        gnt2 = (-1)**(l2-l1+m3)*sqrt((fact(l1)*fact(l2)*gamma(l1+0.5_dp)*gamma(l2+0.5_dp)*fact(l1+l2+m3)*fact(l1+l2-m3)&
                )/(sqrt(pi)*fact(l1+l2)*(2*(l1+l2)+1)*gamma(l1+l2+0.5_dp)*fact(l1+m1)*fact(l1-m1)*fact(l2+m2)*fact(l2-&
                m2)))
        gnt3 = 2*(l2*m1-l1*m2)*gnt2*sqrt((2*(l1+l2)+1.0_dp)/((2*l1)*(2*l2)*(l1+l2-m3)*(l1+l2+m3)))

        ! If l3 = l2 + l1 or l3 = |l2-l1|, calculate the Wigner 3-j symbols from a special value of the reccurence relations.
        if (l2+l1 == l3) then
            gnt = gnt1*gnt2*sqrt(((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*pi))
        else if (l2-l1 == l3) then
            gnt1 = (-1)**(2*l1+l2)*sqrt((fact(l2)*gamma(l1+0.5_dp)*gamma(l2-l1+0.5_dp))/&
                    (sqrt(pi)*fact(l2-l1)*fact(l1)*(2*l2+1)*gamma(l2+0.5_dp)))
            gnt2 = (-1)**(2*l1+m1+m3+l2)*sqrt((fact(l1)*fact(l2-l1)*gamma(l1+0.5_dp)*gamma(l2-l1+0.5_dp)*fact(l2+m2)*&
                    fact(l2-m2))/(sqrt(pi)*fact(l2)*(2*l2+1)*gamma(l2+0.5_dp)*fact(l2-l1+m3)*fact(l2-l1-m3)*&
                    fact(l1+m1)*fact(l1-m1)))
            gnt = gnt1*gnt2*sqrt(((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*pi))
        else if (l1-l2 == l3) then
            gnt1 = (-1)**(2*l2-l1)*sqrt((fact(l1)*gamma(l2+0.5_dp)*gamma(l1-l2+0.5_dp))/&
                    (sqrt(pi)*fact(l1-l2)*fact(l2)*(2*l1+1)*gamma(l1+0.5_dp)))
            gnt2 = (-1)**(2*l2+m2+m3-l1)*sqrt((fact(l2)*fact(l1-l2)*gamma(l2+0.5_dp)*gamma(l1-l2+0.5_dp)*fact(l1+m1)*&
                    fact(l1-m1))/(sqrt(pi)*fact(l1)*(2*l1+1)*gamma(l1+0.5_dp)*fact(l1-l2+m3)*fact(l1-l2-m3)*&
                    fact(l2+m2)*fact(l2-m2)))
            gnt = gnt1*gnt2*sqrt(((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*pi))
        else
        ! Recursively calculate the Wigner 3-j symbols starting from l3 = l2 + l1.
            l3t = l2+l1-1
            l3t1 = l2+l1
            !B = -(2*l3t+1)*(l1*(l1+1)*m3-l2*(l2+1)*m3-l3t*(l3t+1)*(m2-m1t))
            do while ((l3t > abs(l2-l1)) .and. (l3t /= l3))
                A = sqrt(real((l3t**2-m3**2)*(l3t**2-(l1-l2)**2)*((l1+l2+1)**2-l3t**2), dp))
                A1 = sqrt(real(((l3t+1)**2-m3**2)*((l3t+1)**2-(l1-l2)**2)*((l1+l2+1)**2-(l3t+1)**2), dp))
                B = (2*l3t+1)*(l1*(l1+1)*m3-l2*(l2+1)*m3+l3t*(l3t+1)*(m1-m2))
                gnt3_temp = gnt3
                gnt3 = (B*gnt3 - l3t*A1*gnt2)/((l3t+1)*A)
                gnt2 = gnt3_temp
                l3t = l3t-1
                !B = (2*l3t1-1)*(B/(2*l3t1+1)-2*l3t1*(m2-m1t))
            end do
            if (mod(l1+l2+l3t, 2) /= 0) then
                gnt1 = 0.0_dp
            else
                do while ((l3t1 > abs(l2-l1)) .and. (l3t1 /= l3))
                    K = real((l3t1**2-(l1-l2)**2)*((l1+l2+1)**2-l3t1**2), dp)
                    K1 = real(((l3t1-1)**2-(l1-l2)**2)*((l1+l2+1)**2-(l3t1-1)**2), dp)
                    gnt1 = -sqrt(K/K1)*gnt1
                    l3t1 = l3t1-2
                end do
            end if
            gnt = gnt1*gnt3*sqrt(((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*pi))
        end if
    end if

end subroutine gaunt

subroutine unit_real_gaunt(l1, l2, l3, m1, m2, m3, ugnt) ! Algorithm to calculate the integral of three real spherical harmonics, by treating the real
                                                         ! spherical harmonics as a unitary transformation of the standard spherical harmonics.
    implicit none
    integer, parameter :: dp=kind(0.d0)
    integer, intent(in) :: l1, l2, l3, m1, m2, m3
    integer :: i, j, k, kron_delta, step_func
    real(dp) :: gauntval, sumgnt
    complex(dp) :: U1, U2, U3
    real(dp), intent(out) :: ugnt

    ! Apply the angular momentum selection rules to speedup the algorithm.
    if (mod(l1+l2+l3, 2) /= 0) then
        ugnt = 0.0_dp
    else
        ! Compute unitary transformation matrices.
        do i=-l1, l1
            U1 = complex(kron_delta(abs(m1), abs(i))*kron_delta(i, 0)*kron_delta(m1, 0)+kron_delta(abs(m1), abs(i))*&
                    step_func(m1)*(kron_delta(i, m1)+(-1)**i*kron_delta(i, -m1))/sqrt(2.0_dp), kron_delta(abs(m1)&
                    , abs(i))*step_func(-m1)*(-kron_delta(i, m1)*(-1)**(i-m1)+kron_delta(i, -m1)*(-1)**(-m1))/&
                    sqrt(2.0_dp))
            do j=-l2, l2
                U2 = complex(kron_delta(abs(m2), abs(j))*kron_delta(j, 0)*kron_delta(m2, 0)+kron_delta(abs(m2), abs(j))&
                        *step_func(m2)*(kron_delta(j, m2)+(-1)**j*kron_delta(j, -m2))/sqrt(2.0_dp), kron_delta(&
                        abs(m2), abs(j))*step_func(-m2)*(-kron_delta(j, m2)*(-1)**(j-m2)+kron_delta(j, -m2)*(-1)&
                        **(-m2))/sqrt(2.0_dp))
                U3 = complex(kron_delta(abs(m3), abs(-i-j))*kron_delta(-i-j, 0)*kron_delta(m3, 0)+kron_delta(abs(m3),&
                        abs(-i-j))*step_func(m3)*(kron_delta(-i-j, m3)+(-1)**(-i-j)*kron_delta(-i-j, -m3))/sqrt(&
                        2.0_dp), kron_delta(abs(m3), abs(-i-j))*step_func(-m3)*(-kron_delta(-i-j, m3)*(-1)**(-i-j&
                        -m3)+kron_delta(-i-j, -m3)*(-1)**(-m3))/sqrt(2.0_dp))
                call gaunt(l1, l2, l3, i, j, -i-j, gauntval)
                sumgnt = sumgnt + real(U1*U2*U3)*gauntval
            end do
        end do
        ugnt = sumgnt
    end if

end subroutine unit_real_gaunt
