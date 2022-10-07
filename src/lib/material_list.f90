module material_list
    implicit none
    type materiallist
        character(32)   :: matName, matRange  ! Cols 1, 3
            ! cols          2, 4,          5, 6,   7,   8,    9,      10,  
        real(8)        ::   d, numDensity, V, Wth, Wsc, Wabs, sigabs, sigel
            ! cols          11       12
        real(8)        ::   elscatl, abslngth
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
        
        !write(*, '(A,I0)') "Total number of materials: ", nlines
        close(1)

        return
    end function get_lines
! === get_lines ================================================================

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
        character(64), intent(in)   :: INFILE_3
        real(8), dimension(:)       :: vel_list
        integer                     :: num_vels, l, error

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
    subroutine get_materials(inventory, name, nlines, INFILE_2, &
                            &skip_Wsc, skip_Wabs)
        implicit none
        character(128),      intent(in)     :: INFILE_2
        character(64)                       :: name
        integer,            intent(in)      :: nlines
        type(materiallist), dimension(:)    :: inventory
        logical, intent(in)                 :: skip_Wsc, skip_Wabs
        integer                             :: n = 0, error = 0
        

        open(unit = 1, file = INFILE_2, status = "old")
        read(1, *, iostat = error) name
        read(1, *, iostat = error) ! Mulch the header

        do 100 n = 1, nlines
            read(1, *, iostat = error) inventory(n)%matName, &
                &inventory(n)%d, inventory(n)%matRange, &
                &inventory(n)%numDensity, inventory(n)%V, inventory(n)%Wth, &
                &inventory(n)%Wsc, inventory(n)%Wabs, inventory(n)%sigabs, &
                &inventory(n)%sigel, inventory(n)%elscatl, inventory(n)%abslngth
            ! Not all of the above are useful - some are only there for ease
            ! of file reading because kludging gets the job done-ing
                if (skip_Wsc .eqv. .true.) then
                    inventory(n)%Wsc = 0.D0
                end if

                if (skip_Wabs .eqv. .true.) then
                    inventory(n)%Wsc = 0.D0
                end if
100     end do
        
        close(1)
    end subroutine get_materials
! === get_materials ============================================================

end module material_list
