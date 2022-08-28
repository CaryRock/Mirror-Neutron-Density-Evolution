! This will be imported to the "main" function to do the computation. It takes
! in as parameters:
!   Dm, vel, theta0, Vopt, W1, W2, A, and B
!
! Everything else is calculated from those.
! It returns an array of weights, rho_ij, i & j in [1, 4], representing the 
! probability that the neutron is in one of those states.
module banfor_module
contains
    subroutine banfor(Dm, vel, theta0, Vopt, W1, W2, lambda, A, B, tStep, rho)
        implicit none
    
        real hbar, nmass, Qe, Dm, vel, s2ze, TOF, theta0, time, eps, tStep
        real theta, omega, Vopt, V, W, W1, W2, A, B, U1, U2, WW1, WW2, DE, lambda
        complex(8) pnn8, pmn8, pnm8, pmm8
        real(8) pnn, pmn, pnm, pmm, unitar
        complex i, zeta, D, zeta2, arg, cze, cze2, sze, sze2, H1, H2, arg2
        complex czestar, szestar, cze2star, sze2star, uroo, roo
        complex(8) snn, smn, snm, smm, H1C, H2C, snnC, smnC, snmC, smmC
        complex(8) aee1, aee2, ee1, ee2, ee1C, ee2C, ee1CC, ee2CC
        integer nti, ntf, ntdiff
        
        real(8), dimension(2,2) :: rho
        
        B = 0.D0

        if (theta0 .ge. 0.8) theta0 = 0.785398163
        hbar = 6.582119569E-16
        
        V       = Vopt - Dm - B                        ! delta_m in eV
        i       = cmplx(0.0, 1.0)                   ! sqrt(-1.0)
        W       = W2 + vel*W1
        eps     = 0.5 * abs(Dm) * tan(2. * theta0)  ! epsilon for given Dm and theta0
        
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
        
        !time    = lambda / vel
        time    = tStep
        TOF     = time / hbar
        
        aee1    = -i * H1 * TOF
        ee1     = cdexp(aee1)
        aee2    = -i * H2 * TOF
        ee2     = cdexp(aee2)
        ee1C    = conjg(ee1)
        ee2C    = conjg(ee2)
        
        snn     = cze2 * ee1 + sze2 * ee2
        snnC    = cze2star * ee1C + sze2star * ee2C
        smn     = cze * sze * (ee1 - ee2)
        smnC    = czestar * szestar * (ee1C - ee2C)
        snm     = cze * sze * (ee1 - ee2)!smn
        snmC    = czestar * szestar * (ee1C - ee2C)!snmC
        smm     = sze2 * ee1 + cze2 * ee2
        smmC    = sze2star * ee1C + cze2star * ee2C
        
        pnn     = realpart(snn * snnC)
        pnm     = realpart(snm * snmC)
        pmn     = realpart(smn * smnC)
        pmm     = realpart(smm * smmC)
        
!        unitar  = pnn + pmm
        
        rho(1, 1)   = pnn
        rho(1, 2)   = pnm
        rho(2, 1)   = pmn
        rho(2, 2)   = pmm
        
        return
    end subroutine banfor
end module banfor_module