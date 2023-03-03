module NNp_Loop
  use material_list
  use distributions
  use exact_banfor_module
  !$ use OMP_LIB
  implicit none

  type loop_variables
    type(materiallist), allocatable :: inventory(:)
    real  :: rhoReset(2, 2), rho(2, 2)
    real  :: O(4), P(4)
    real, allocatable :: vel_list(:), matAveN(:), matAVeM(:)

    real  :: vel, Dm, theta0, tStep, PN, PM, dlina, Vopt, PNvels, PMvels
    real  :: Wabs, x, dx, matD, remThick, elLambda, Wsc

    integer :: nMaterials, nVels, numSteps, i, j, k, l, m
  end type loop_variables

contains
  recursive subroutine scttStep(lambda, remThick, variate)
    implicit none
    real, intent(in)  :: lambda
    real              :: remThick, variate
    
    call random_number(variate)
    variate = -lambda * log(variate)
    if (variate .ge. remThick) variate = remThick
    remThick = remThick - variate
  end subroutine scttStep

  subroutine get_log_data(inventory, rhoReset, rho, O, P, &
    &vel_list, PNvels, PMvels, Dm, theta0, nMaterials, nVels, &
    &header, N_initial, matAveN, matAveM, no_scattering)
    implicit none
    type(materiallist), dimension(:), allocatable :: inventory
    real, dimension(2, 2), intent(in) :: rhoReset
    real, dimension(2, 2) :: rho
    real, dimension(4)    :: O, P
    real, dimension(:), allocatable :: vel_list
    real, dimension(:), allocatable :: matAveN, matAveM

    real  :: vel, Dm, theta0, tStep, PN, PM, dlina, Vopt, PNvels, PMvels
    real  :: Wabs, x, dx, matD, remThick, elLambda, Wsc

    integer, intent(in) :: nMaterials, nVels
    integer :: i, j, k, l, m, numSteps

    character(32), dimension(6) :: header
    logical, intent(in) :: N_initial, no_scattering
    ! Call the passed function loop - probably main_vel_loop below - and
    ! return the desired 

    call main_vel_loop(inventory, rhoReset, rho, O, P, &
    &vel_list, PNvels, PMvels, Dm, theta0, nMaterials, nVels, &
    &header, N_initial, matAveN, matAveM, no_scattering)
    
  end subroutine get_log_data

  recursive subroutine main_vel_loop(inventory, rhoReset, rho, O, P, &
    &vel_list, PNvels, PMvels, Dm, theta0, nMaterials, nVels, &
    &header, N_initial, matAveN, matAveM, no_scattering)
    implicit none
    type(materiallist), dimension(:), allocatable :: inventory
    real, dimension(2, 2), intent(in) :: rhoReset
    real, dimension(2, 2) :: rho
    real, dimension(4)    :: O, P
    real, dimension(:), allocatable :: vel_list
    real, dimension(:), allocatable :: matAveN, matAveM

    real  :: vel, Dm, theta0, tStep, PN, PM, dlina, Vopt, PNvels, PMvels
    real  :: Wabs, x, dx, matD, remThick, elLambda, Wsc

    integer, intent(in) :: nMaterials, nVels
    integer :: i, j, k, l, m, numSteps

    character(32), dimension(6) :: header
    logical, intent(in) :: N_initial, no_scattering

        ! # Mtrl MtrlThck G.Dist PN PM
!    !$OMP MASTER
!    print '(A2, A35, A28, A30, A30, A27)', "# ", header(2), header(3), &
!      &header(4), header(5), header(6)
!    print *, adjustl(header(1)), ",", 20.0, ",", 20.0,",", 1.0,",", 0.0
!    !$OMP END MASTER
    !!$OMP PARALLEL DO Private(i, j, k, l, m, vel, tStep, rho, O, P, &
    !!$OMP& PN, PM, dlina, Vopt, Wabs, numSteps, x, dx, matD, remThick,&
    !!$OMP& elLambda, temp)

vl: do k = 1, nVels
      if (N_initial) then
        PN = 1.0
        PM = 0.0
      else
        PN = 0.0
        PM = 1.0
      end if

      vel = vel_list(k)
      rho = rhoReset

nMat: do l = 1, 1!nMaterials
        matD = inventory(l)%d           ! m
        numSteps = inventory(l)%steps   !
        Vopt = inventory(l)%V           ! eV
        Wabs = inventory(l)%Wabs        ! eV
        if (.not. no_scattering) Wsc  = inventory(l)%Wsc         ! eV
        elLambda = inventory(l)%elscatl ! m

        remThick = matD
        m = 0

mStep:  do while(remThick .gt. 0)
          m = m + 1
          rho = rhoReset
          call scttStep(elLambda, remThick, dx)
          tStep = dx / vel
          call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, tStep, rho)
          O(1) = rho(1, 1)
          O(2) = rho(1, 2)
          O(3) = rho(2, 1)
          O(4) = rho(2, 2)

          if (l .eq. 1 .and. m .eq. 1) then
            if (N_initial) then
              P(1) = PN * O(1)
              P(2) = PN * O(2)
              P(3) = 0.0
              P(4) = 0.0
            else
              P(1) = 0.0
              P(2) = 0.0
              P(3) = PM * O(3)
              P(4) = PM * O(4)
            end if
          end if

          PN = P(1) + P(3)
          PM = P(2) + P(4)

          P(1) = PN * O(1)
          P(2) = PN * O(2)
          P(3) = PM * O(3)
          P(4) = PM * O(4)
        end do mStep  ! matSteps

        ! For later outputting the desired probabilities as a function of material
        ! for a given mass and angle.
        matAveN(l) = matAveN(l) + PN
        matAveM(l) = matAveM(l) + PM
      end do nMat     ! nMaterials

      PNvels = PNvels + PN
      PMvels = PMvels + PM
    end do vl   ! nVels

  end subroutine main_vel_loop
end module NNp_Loop