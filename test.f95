program arg_test
    implicit none
    character(100)  :: num1char, num2char
    real num1, num2, numsum

    if (command_argument_count() .ne. 2) then
        write(*, *) "Error, two command-line arguments required. Stopping."
        stop 1
    end if

    call get_command_argument(1, num1char)
    call get_command_argument(2, num2char)

    read(num1char, *) num1
    read(num2char, *) num2

    numsum = num1 + num2
    write(*, *) numsum

end program arg_test
