program dead_simple
    use exact_banfor_module
    implicit none

    complex(8), dimension(2)    :: psi
    real(8), dimension(2, 2)    :: rho
    real(8), parameter          :: hbar = 6.582119569D-16
    real(8)                     :: Dm, vel, theta0, Vopt, Wsc, Wabs, lambda
    real(8)                     :: A, B, tStep, d, R, wStart, wEnd
    integer                     :: w, N, wDivs

    psi = (/ cmplx(1.D0, 0.D0), cmplx(0.D0, 0.D0)/)
    data rho / 1.D0, 0.D0, 0.D0, 0.D0 /
    Dm = 300.0D-9
    vel = 1000 * 100    ! m/s -> cm / s 
    theta0 = 0.1
    Vopt = 166D-9
    Wsc = 0.D0
    lambda = 0.D0
    A = 0.D0
    B = 0.D0
    d = 2.2
    tStep = d / vel
!    N = 200

    wStart = -15.D0 
    wEnd = -6.D0
    wDivs= 20

    N = int(-(wStart - wEnd)*float(wDivs))

    write(unit = 6, fmt = '(A1, A16, 3A17)') "#", "Wabs", "rho_11", "rho_22", "abs"
    do w = 0, N
        psi = (/ cmplx(1.D0, 0.D0), cmplx(0.D0, 0.D0)/)
        data rho / 1.D0, 0.D0, 0.D0, 0.D0 /
        Wabs = 10**(wStart + w / float(wDivs))
        !Wabs = 3.0D-9
        call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, lambda, A, B, tStep, psi, rho)

        R = (1.D0 - rho(1, 1) - rho(2, 2)) / &
                &(1.D0 - exp(-2.D0 * d * wabs / (hbar * vel)))
        write(unit = 6, fmt = '(4ES17.8E3)') Wabs, rho(1, 1), rho(2, 2), R
    end do
end program dead_simple