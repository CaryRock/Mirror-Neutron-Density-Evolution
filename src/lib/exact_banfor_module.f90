module exact_banfor_module
contains
    recursive subroutine exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, tStep, rho)
        implicit none
        ! Input argument declarations
        real    Dm, vel, theta0, Vopt, Wsc, Wabs, tStep
        complex, dimension(2, 2)  :: s
        real, dimension(2, 2)     :: rho
        
        ! Required variable declarations - not from invocation
        real        hbar, TOF, time, eps   
        real        V, W, U1, U2, WW1, WW2, DE
!        real        unitar
!        real        lambda, A, B, nmass, Qe, theta, omega
!        complex     s2ze, D, arg2
        complex     i, zeta, zeta2, arg, cze, cze2, sze, sze2, H1, H2
        complex     czestar, szestar, cze2star, sze2star, uroo, roo             
        complex     H1C, H2C
        complex     aee1, aee2, ee1, ee2, ee1C, ee2C
        
        !lambda  = 0.0
        !A       = 0.0
        !B       = 0.0
        !nmass   = 0.0
        !Qe      = 0.0
        !s2ze    = cmplx(0.0, 0.0)
        !theta   = 0.0
        !omega   = 0.0
        !D       = cmplx(0.0, 0.0)
        !arg2    = cmplx(0.0, 0.0)

        ! Maximum angle possible - 45 degrees
        if (theta0 .ge. 0.785398163) theta0 = 0.785398163
        hbar = 6.582119569E-16

        V       = Vopt - Dm
        i       = cmplx(0.0, 1.0)
        W       = Wabs + vel * Wsc
        eps     = 0.5 * abs(Dm) * tan(2. * theta0)
        arg     = 2. * eps / (V - i*W)
        zeta2   = atan(arg)
        zeta    = 0.5 * zeta2
        cze     = cos(zeta)
        sze     = sin(zeta)
        cze2    = cze * cze
        sze2    = sze * sze
        
        czestar = conjg(cze)
        szestar = conjg(sze)
        cze2star= conjg(cze2)
        sze2star= conjg(sze2)
        
        uroo    = (V - i*W)**2 + 4. * eps * eps
        roo     = sqrt(uroo)
        if (real(zeta) .lt. 0.0) roo = -roo
        H1      = 0.5 * (V - i*W + roo)
        H2      = 0.5 * (V - i*W - roo)
        H1C     = conjg(H1)
        H2C     = conjg(H2)
        
        U1      = real(H1)
        U2      = real(H2)
        DE      = U1 - U2
        WW1     = aimag(H1)
        WW2     = aimag(H2)

        time    = tStep
        TOF     = time / hbar

        aee1    = -i * H1 * TOF
        ee1     = cexp(aee1)
        aee2    = -i * H2 * TOF
        ee2     = cexp(aee2)
        ee1C    = conjg(ee1)
        ee2C    = conjg(ee2)

        s(1, 1) = cze2 * ee1 + sze2 * ee2
        s(2, 2) = sze2 * ee1 + cze2 * ee2
        s(1, 2) = cze * sze * (ee1 - ee2)
        s(2, 1) = cze * sze * (ee1 - ee2)
        
        !psi = matmul(s, psi)

        !rho(1, 1) = real(psi(1) * conjg(psi(1)))
        !rho(2, 2) = real(psi(2) * conjg(psi(2)))
        !rho(1, 2) = real(psi(1) * conjg(psi(2)))
        !rho(2, 1) = real(psi(2) * conjg(psi(1)))
        rho(1, 1) = real(s(1, 1) * conjg(s(1, 1)))
        rho(1, 2) = real(s(1, 2) * conjg(s(2, 1)))
        rho(2, 1) = real(s(2, 1) * conjg(s(1, 2)))
        rho(2, 2) = real(s(2, 2) * conjg(s(2, 2)))

        !print *, "rho = ", rho
        return 
    end subroutine exactBanfor
end module exact_banfor_module