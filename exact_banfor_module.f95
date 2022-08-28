module exact_banfor_module
contains
    subroutine exactBanfor(Dm, vel, theta0, Vopt, W1, W2, lambda, A, B, tStep, psi, rho)
        implicit none
        ! Input argument declarations
        real(8) Dm, vel, theta0, Vopt, W1, W2, lambda, A, B, tStep
        complex(8), dimension(2, 2) :: s
        real(8), dimension(2, 2)    :: rho
        complex(8), dimension(2)    :: psi, psi2
        
        ! Required variable declarations - not from invocation
        real(8)     hbar, nmass, Qe, s2ze, TOF, time, eps   ! originally real
        real(8)     theta, omega, V, W, U1, U2, WW1, WW2, DE! originally real
!        real(8)     unitar
        complex(8)  i, zeta, D, zeta2, arg, cze, cze2, sze, sze2, H1, H2, arg2  ! originally complex
        complex(8)  czestar, szestar, cze2star, sze2star, uroo, roo             ! originally complex
        complex(8)  H1C, H2C
        complex(8)  aee1, aee2, ee1, ee2, ee1C, ee2C
!        integer     nti, ntf, ntdiff
        
        ! Variables required for exact solution
        !complex(8) psiN, psiM, psiNC, psiMC
        complex(8)  psiN, psiM
        
        lambda  = 0.D0
        A       = 0.D0
        B       = 0.D0
        nmass   = 0.D0
        Qe      = 0.D0
        s2ze    = 0.D0
        theta   = 0.D0
        omega   = 0.D0
        D       = 0.D0
        arg2    = cmplx(0.D0, 0.D0, kind = 8)

        ! Maximum angle possible - 45 degrees
        if (theta0 .ge. 0.785) theta0 = 0.785398163
        hbar = 6.582119569E-16

        V       = Vopt - Dm - B
        i       = cmplx(0.0, 1.0)
        W       = W2 + vel * W1
        eps     = 0.5 * abs(Dm) * tan(2. * theta0)
        
        arg     = 2. * eps / (V - i*W)
        zeta2   = atan(arg)
        zeta    = 0.5 * zeta2
        
        psiN    = psi(1)
        psiM    = psi(2)

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

        !time   = lambda / vel
        time    = tStep
        TOF     = time / hbar

        aee1    = -i * H1 * TOF
        ee1     = cdexp(aee1)
        aee2    = -i * H2 * TOF
        ee2     = cdexp(aee2)
        ee1C    = conjg(ee1)
        ee2C    = conjg(ee2)

        s(1, 1) = cze2 * ee1 + sze2 * ee2
        s(2, 2) = sze2 * ee1 + cze2 * ee2
        s(1, 2) = cze * sze * (ee1 - ee2)
        s(2, 1) = cze * sze * (ee1 - ee2)
        
        psi2 = matmul(s, psi)

        rho(1, 1) = realpart(psi2(1) * conjg(psi2(1)))
        rho(2, 2) = realpart(psi2(2) * conjg(psi2(2)))
        rho(1, 2) = realpart(psi2(1) * conjg(psi2(2)))
        rho(2, 1) = realpart(psi2(2) * conjg(psi2(1)))

        !print *, "rho = ", rho
        return 
    end subroutine exactBanfor
end module exact_banfor_module