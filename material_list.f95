module material_list
    implicit none
    type materiallist
        character(32)   :: matName, matRange  ! Cols 1, 3
            ! cols     2, 4,          5, 6,   7,   8,    9,      10,    11
        real(8)        :: d, numDensity, V, Wth, Wsc, Wabs, sigabs, sigel, elscatl
    end type materiallist

contains
    ! Simply gets the number of lines in the given file file
    ! Returns the number of lines contained in the file
    !
    ! Input:
    !   character INFILE_2: Filename - assumed to be the list of materials
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
        
        open(1, file = trim(INFILE_2), status = "old")
        read(1, *)  ! Skip the filename
        read(1, *)  ! Skip the velocity
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
    subroutine get_materials(inventory, name, velocity, nlines, INFILE_2)
        implicit none
        character(64),      intent(in)                  :: INFILE_2
        character(64)                                   :: name
        real(8)                                         :: velocity
        integer,            intent(in)                  :: nlines
        type(materiallist), dimension(:), intent(inout) :: inventory
        character(256)                                  :: line
        integer                                         :: n = 0, error
        
        open(unit = 1, file = INFILE_2, status = "old")
        read(1, *, iostat = error) name
        read(1, *, iostat = error) velocity
        read(1, *, iostat = error) line  ! Grab header line
        
        do 100 n = 1, nlines
            read(1, *, iostat = error) inventory(n)%matName, &
                inventory(n)%d, inventory(n)%matRange, &
                inventory(n)%numDensity, inventory(n)%V, inventory(n)%Wth, &
                inventory(n)%Wsc, inventory(n)%Wabs, inventory(n)%sigabs, &
                inventory(n)%sigel, inventory(n)%elscatl
            ! Not all of the above are useful - some are only there for ease
            ! of file reading because kludging gets the job done-ing
100     end do
        
        close(1)
        return
    end subroutine get_materials
end module material_list