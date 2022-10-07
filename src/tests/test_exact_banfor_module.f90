program test_exact_banfor
        !use material_list
        use exact_banfor_module
        implicit none

        ! Needs Dm, vel, theta0, Vopt, Wsc, Wabs, lambda, A, B, tStep, Psi, and rho
        
        ! Dm:           real(8) Mass difference
        ! vel:          real(8) velocity
        ! theta0:       real(8) mixing angle
        ! Vopt:         real(8) real part of the optical potential of the material
        ! Wsc:          real(8) Scattering potential of the material
        ! Wabs:         real(8) Absorption potential of the material
        ! lambda:       real(8) I don't remember what lambda is, but it's here for future use
        ! A:            real(8) I assume this is the scalar potential? Not used
        ! B:            real(8) Magnetic field strength. Not currently assigned to nonzero.
        ! tStep:        real(8) The time part of the ToF; time particle has been travelling
        ! psi:          complex((2), 8) spinor of current particle phase (e.g., neutron, mirror-neutron, or some mix)
        ! rho:          real((2, 2), 8) density matrix

        real(8) :: Dm, vel, theta0, Vopt, Wsc, Wabs, lambda, A, B, tStep
        complex(8), dimension(2)        :: psi
        real(8), dimension(2, 2)        :: rho
        integer :: jo, jn
        ! Copying the basic details of another file for Cd

        Dm = 100.0E-9
        vel = 1000.0
        tStep = 0.0035 / vel
        Vopt = 5.377E-8
        Wsc = 9.914E-15
        Wabs = 8.4558E-9

        lambda = 0.0
        A = 0.0
        B = 0.0

        psi = (/ cmplx(1.0, 0.0), cmplx(0.0, 0.0) /)
        rho(1, 1) = 1.0
        rho(1, 2) = 0.0
        rho(2, 1) = 0.0
        rho(2, 2) = 0.0


    do 500 jo = 1, 6
        do 501 jn = 1, 9
            theta0 = jn * 10. ** (-7 + jo)
            call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, lambda, A, B, &
                            &tStep, psi, rho)
            write(6, 57) theta0, rho(1, 1), rho(2, 2)
501     end do
500 end do

57  format(F10.8, 1P2E16.6)
end program test_exact_banfor
