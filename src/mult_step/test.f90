program test
  use f90getopt
  use material_list
  use exact_banfor_module
  use get_parameters
  use distributions
  implicit none
  ! KISS!
  ! Just hard-code most of these things - it doesn't have to be *complete*, it
  ! just has to demonstrate correct working of the new implementation of the
  ! step sampling method

  type(materiallist), dimension(:), allocatable  :: inventory
  integer, parameter                :: prec = 18

  ! Arrays
  complex, dimension(2) :: psi

  real, dimension(2, 2) :: rho
  real, dimension(2, 2), parameter  :: rhoN = &
    &reshape((/1.0, 0.0, 0.0, 0.0/), shape(rho))
  real, dimension(:, :), allocatable:: meshNavg, meshMavg
  real, dimension(:), allocatable   :: vel_list, Masses, Angles
  real, dimension(4)                :: O, P

  ! Scalars
  real  :: PN, PM, dlina, Dm, theta0, vel, Vopt, Wabs, Wsc, tStep, x
  real  :: lambda, Navg, Mavg, PNvels, PMvels, vAvg, temp
  real  :: dx, matD, remThick, elLambda

  integer :: i, j, k, l, m, numSteps, nMaterials, nMasses, nAngles, nVels
  integer :: kmin, kmax, totSteps, accumulatedSteps

  ! Setup inventory array
  nMaterials = 1
  allocate(inventory(1))
  inventory(1)%matName = "b4c_test"
  inventory(1)%d = 0.008          ! 8 mm
  inventory(1)%V = 199.2E-9       ! eV
  inventory(1)%Wsc = 0.0          ! eV
  inventory(1)%Wabs= 6.102E-9     ! eV
  inventory(1)%elscatl = .01373   ! m
  inventory(1)%absscatl= 1.187E-4 ! m
  inventory(1)%steps = 2          !

  ! Setup deltaMs
  !nMasses = 261
  nMasses = 1
  allocate(Masses(nMasses))
  if (nMasses .gt. 1) then
    Masses(1) = 1e-9
    do i = 2, 20*13 + 1
      Masses(i) = 10**(-20. + float(i)*.05)
    end do
  else
    Masses(1) = 2e-7
  end if
  ! Setup thetas
!  nAngles = 159
  nAngles = 1
  allocate(Angles(nAngles))
  if (nAngles .gt. 1) then
    Angles(1) = 1e-8
    do i = 2, 159
      theta0 = 10**(-8 + float(i)*.05)
      if (theta0 .lt. 0.785) Angles(i) = theta0
    end do
  else
    Angles(1) = 0.001
  end if

  ! Set number of velocities to use
  nVels = 1
  temp = 340.
  allocate(vel_list(nVels))

massLoop: do i = 1, nMasses
    Dm = Masses(i)

anglLoop: do j = 1, nAngles
      theta0 = Angles(j)
      
      do k = 1, nVels
        call YK_MAXW(temp, vel_list(k))
      end do

      vAvg = vAvg + sum(vel_list)
      vAvg = vAvg / nVels

      PNvels = 0.0
      PMvels = 0.0

velsLoop: do k = 1, nVels
        PN = 1.0
        PM = 0.0

        vel = vel_list(k)
        rho = rhoN

! This whole mtrlLoop section needs to be rewritten to accomodate the fact that
! the step size basically changes every step of the way. Some things do remain
! constants, but most do end up changing via random number.
mtrlLoop: do l = 1, nMaterials
          matD = inventory(l)%d
          numSteps = inventory(l)%steps
          elLambda = inventory(l)%elscatl
          ! tStep = inventory(l)%elscatl / 100. / vel
          ! Set elscatl to m in material file
          Vopt = inventory(l)%V     ! Set these to
          Wabs = inventory(l)%Wabs  !              eV in material file

          remThick = matD

            m = 1
matSteps: do while (remThick .gt. 0.)
            print *, remThick
            rho = rhoN  ! Reset rho

            call scttStep(elLambda, remThick, dx) ! Get scatt. length
            tStep = dx / vel          ! Get dt for that length at given vel
              ! Calculate BANFOR
            call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, tStep, psi, rho)
            O(1) = rho(1, 1)
            O(2) = rho(1, 2)
            O(3) = rho(2, 1)
            O(4) = rho(2, 2)
            print *, dx, theta0, O
          
            if (l .eq. 1) then
              P(1) = PN * O(1)
              P(2) = PN * O(2)
              P(3) = 0.0
              P(4) = 0.0
            end if
            
            print *, remThick
            PN = P(1) + P(3)
            PM = P(2) + P(4)

            P(1) = PN * O(1)
            P(2) = PN * O(2)
            P(3) = PM * O(3)
            P(4) = PM * O(4)

!            if (remThick .eq. 0.) exit
            m = m + 1
          end do matSteps
        end do mtrlLoop

        PNvels = PNvels + PN
        PMvels = PMvels + PM
      end do velsLoop

      !meshNavg(i, j) = PNvels / float(nVels)
      !meshMavg(i, j) = PMvels / float(nVels)
    end do anglLoop
  end do massLoop

contains
  ! Takes in a scattering length and (remaining) material length, returns an
  ! appropriate step based on those lengths
  subroutine scttStep(lambda, remThick, variate)
    ! -l * ln(rand)
    implicit none
    real, intent(in)  :: lambda
    real              :: remThick, variate

    call random_number(variate)
    variate = -lambda * log(variate)
    if (variate .ge. remThick) variate = remThick
    remThick = remThick - variate
    ! return variate, altered remThick
  end subroutine scttStep
end program test
