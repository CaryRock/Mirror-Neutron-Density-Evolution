program test
  use exact_banfor_module
  implicit none
  real, dimension(2, 2) :: rho
  complex, dimension(2)    :: psi
  real                  :: Dm, vel, theta0, Vopt, Wsc, Wabs, tStep

  psi(1) = cmplx(1.0, 0.0)
  psi(2) = cmplx(0.0, 0.0)

  rho = 0.0
  rho(1, 1) = 1.0

  print *, "rho before: ", rho

  open(unit = 1, file = "banfor_vars.in", status = "old")
  read(1, *)  Dm
  read(1, *)  vel
  read(1, *)  theta0
  read(1, *)  Vopt
  read(1, *)  Wsc
  read(1, *)  Wabs
  read(1, *)  tStep
  close(1)

  call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, tStep, psi, rho)

  print *, "rho: ", rho
  print *, ""
end program test
