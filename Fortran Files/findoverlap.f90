function kron_delta(i, j) result(kd) ! Kronecker delta function.

    implicit none
    integer, intent(in) :: i, j
    integer :: kd

    kd = 0
    if (i == j) kd = 1

end function kron_delta

subroutine findoverlap(cvt, lvt, mvt, base_num, overlap) ! Calculates the overlap matrix.

    implicit none
    integer, parameter :: dp=kind(0.d0) ! Specifies double precision more robustly.
    double precision, intent(in) :: cvt(:)
    integer, intent(in) :: lvt(:), mvt(:)
    integer, intent(in) :: base_num
    integer :: i, j, kron_delta
    real(dp) :: gam
    double precision, intent(out) :: overlap(base_num, base_num)

    overlap = 0.0_dp
    do j=1, size(lvt)
        do i=1, size(lvt)
            gam = gamma(0.5_dp*(lvt(i)+lvt(j))+1.5_dp)/sqrt(gamma(lvt(i)+1.5_dp)*gamma(lvt(j)+1.5_dp))
            overlap(i, j) = (kron_delta(lvt(i), lvt(j))*kron_delta(mvt(i), mvt(j))*gam*(2.0_dp*cvt(j))**(0.5_dp*lvt(j)&
                             +0.75_dp)*(2.0_dp*cvt(i))**(0.5_dp*lvt(i)+0.75_dp))/(cvt(j)+cvt(i))**(0.5_dp*(lvt(i)+&
                             lvt(j))+1.5_dp)
        end do
    end do

end subroutine findoverlap
