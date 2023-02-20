module material_list
    implicit none
    !type materiallist
    !    character(32)   :: matName, matRange  ! Cols 1, 3
    !        ! cols          2, 4,          5, 6,   7,   8,    9,     
    !    real            ::  d, numDensity, V, Wth, Wsc, Wabs, sigabs
    !        ! cols          10,      11,      12,
    !    real            ::  sigel, elscatl, abslngth
    !        ! cols            13
    !    integer         ::  steps
    !end type materiallist

  type materiallist
    character(32) :: matName
    real          :: d, V, Wsc, Wabs, elscatl, absscatl
    integer       :: steps
  end type materiallist
contains

! === get_lines ================================================================
    ! Simply gets the number of lines in the given file file
    ! Returns the number of lines contained in the file
    !
    ! Input:
    !   character INFILE_2: Filename - the list of materials
    !
    ! Output:
    !   integer nlines: the number of lines in INFILE_2
    integer function get_lines(INFILE_2) result(nlines)
        implicit none
        character(64),  intent(in)              :: INFILE_2
        character(256)                          :: line
        integer                                 :: error

        nlines = 0
        error = 0
        
        open(unit = 1, file = trim(adjustl(INFILE_2)), status = "old")
        read(1, *)  ! Skip the filename
        read(1, *)  ! Skip the header line
        do while(error == 0)
            nlines = nlines + 1
            read(1, *, iostat = error) line
        end do
        nlines = nlines - 1 ! Comment at the beginning - don't want that
        
        write(6, '(A,I0)') "Total number of materials: ", nlines
        close(1)

        return
    end function get_lines
! === get_lines ================================================================

! === get_num_masses ===========================================================
    integer function get_num_masses(INFILE_4) result(nMasses)
        implicit none
        character(64), intent(in)   :: INFILE_4
        character(256)              :: line
        integer                     :: error

        nMasses = 0
        error = 0

        open(unit = 1, file = trim(adjustl(INFILE_4)), status = "old")
        do while(error == 0)
            nMasses = nMasses + 1
            read(1, *, iostat = error) line
        end do
        close(1)

        nMasses = nMasses - 1   ! Overcounting
        return
    end function get_num_masses
! === get_num_masses ===========================================================

! === get_num_angles ===========================================================
    integer function get_num_angles(INFILE_5) result(nAngles)
        implicit none
        character(64), intent(in)   :: INFILE_5
        character(256)              :: line
        integer                     :: error

        nAngles = 0
        error = 0

        open(unit = 1, file = trim(adjustl(INFILE_5)), status = "old")
        do while(error == 0)
            nAngles = nAngles + 1
            read(1, *, iostat = error) line
        end do
        close(1)

        nAngles = nAngles - 1   ! Overcounting
        return
    end function get_num_angles
! === get_num_angles ===========================================================

! === get_vel_lines ============================================================
    integer function get_vel_lines(INFILE_3) result(nlines)
        implicit none
        character(64), intent(in)               :: INFILE_3
        character(256)                          :: line
        integer                                 :: error

        nlines = 0
        error = 0

        open(unit = 1, file = trim(INFILE_3), status = "old")
        do while(error == 0)
            nlines = nlines + 1
            read(1, *, iostat = error) line
        end do
        close(1)

        nlines = nlines - 1 ! Don't need the extra line
        return
    end function get_vel_lines

    subroutine get_velocities(INFILE_3, num_vels, vel_list)
        implicit none
        character(64), intent(in):: INFILE_3
        real, dimension(:)       :: vel_list
        integer                  :: num_vels, l, error

        error = 0
        open(unit = 1, file = trim(INFILE_3), status = "old")
        do 101 l = 1, num_vels
            read(1, *, iostat = error) vel_list(l)
101     end do
        close(1)
        end subroutine get_velocities
! === get_vel_lines ============================================================

! === get_materials ============================================================
    ! getMaterials(...) is called by the main function to read in the data file
    ! "material.list".
    !
    ! The output is to create an array of some meaningful data type which will 
    ! then be iterated over in the main loop of the file/program.
    !
    ! Input: 
    !   type(...), dimension(:) :: inventory 
    !   integer, intent(in)     :: nlines
    !   character(64), intent(in)   :: INFILE_2

    ! Output:
    !   type(...), dimension(:) :: inventory
    !       Inventory will have been "filled in" appropriately
    subroutine get_materials(inventory, name, nlines, INFILE_2)
        implicit none
        character(128),      intent(in)     :: INFILE_2
        character(64)                       :: name
        integer,            intent(in)      :: nlines
        type(materiallist), dimension(:)    :: inventory
        integer                             :: n = 0, error = 0
        

        open(unit = 1, file = INFILE_2, status = "old")
        read(1, *, iostat = error) name
        read(1, *, iostat = error) ! Mulch the header

        do 100 n = 1, nlines
            read(1, *, iostat = error) inventory(n)%matName, &
                &inventory(n)%d, inventory(n)%V, &
                &inventory(n)%Wsc, inventory(n)%Wabs,&
                &inventory(n)%elscatl, inventory(n)%absscatl, inventory(n)%steps
100     end do
        
        close(1)
    end subroutine get_materials
! === get_materials ============================================================

end module material_list