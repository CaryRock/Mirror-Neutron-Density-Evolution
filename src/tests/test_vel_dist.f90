program vel_dist_test
    use distributions
    implicit none
    real, dimension(:), allocatable :: velocities, PDF, CDF, sampled
    character(80), dimension(4)     :: file_names
    real                            :: temp, vmin, vmax
    integer, parameter              :: num_samples = int(1e6)
    integer                         :: i

    allocate(velocities(num_samples + 1))
    allocate(PDF(num_samples + 1))
    allocate(CDF(num_samples))
    allocate(sampled(num_samples))

    temp = 340  ! K
    vmin = 300  ! m/s
    vmax = 12000! m/s

    velocities  = 0.0
    PDF         = 0.0
    CDF         = 0.0
    sampled     = -1.0

    file_names(1) = "pdf.txt"
    file_names(2) = "cdf.txt"
    file_names(3) = "velocities.txt"
    file_names(4) = "sampled.txt"

    call setup_MBDist(num_samples, temp, velocities, PDF, CDF, vmin, vmax)

    do i = 1, num_samples
        CDF(i) = CDF(i) / CDF(num_samples)
    end do

    do i = 1, num_samples
        sampled(i) = sample_MBDist(velocities, CDF, num_samples)
    end do

    call write_to_file(file_names(1), PDF, num_samples)
    call write_to_file(file_names(2), CDF, num_samples + 1)
    call write_to_file(file_names(3), velocities, num_samples + 1)
   ! open(unit = 21, file = trim(adjustl(file_names(3))))
   ! do i = 1, num_samples + 1
   !     write(21, *) velocities(i)
   ! end do
   ! close(unit = 21)
    call write_to_file(file_names(4), sampled, num_samples)

contains
    subroutine write_to_file(file_name, data, limit)
        implicit none
        character(80), intent(in)   :: file_name
        real, dimension(:)          :: data
        integer, intent(in)         :: limit
        integer                     :: i

        open(unit = 21, file = trim(adjustl(file_name)))
        do i = 1, num_samples
            write(21, *) data(i)
        end do
        close(unit = 21)
    end subroutine write_to_file
end program vel_dist_test