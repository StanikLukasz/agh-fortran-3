program main
    use sequiential_lib
    use coarray_lib
    implicit none
    
    ! variables
    integer (kind = 4) :: status, N
    real (kind = 8), dimension(:,:), allocatable :: A, B, X
    real (kind = 16) :: before, after ! for time measuring
    character (len=16) :: arg ! for getting from command line
    
    ! getting argument - size of matrices
    call GET_COMMAND_ARGUMENT(1, arg)
    read(arg, *) N
    
    ! preparing matrices
    allocate(A(N,N))
    allocate(B(N,N))
    allocate(X(N,N))
    
    A = 1.6d0
    B = 4.2d0
    
    ! opening log file
    if (THIS_IMAGE() .EQ. 1) then
        open (unit=19, file="results.txt", position="append", &
            form="formatted", action="write")
    end if
    
    !!! getting times
    
    ! SEQUENTIAL - multiplying
    if (THIS_IMAGE() .EQ. 1) then
        call CPU_TIME(before)
        call mult_seq(A,B,X,status)
        call CPU_TIME(after)
        
        write(19,*)"mult_seq:   ", NUM_IMAGES(), N, (after-before)
    end if
    
    syncall()
    
    ! COARRAY - multiplying
    call CPU_TIME(before)
    call mult_coarr(A,B,X,status)
    call CPU_TIME(after)
    
    if (THIS_IMAGE() .EQ. 1) then
        write(19,*)"mult_coarr: ", NUM_IMAGES(), N, (after-before)
    end if
    
    syncall()
    
    ! SEQUENTIAL - Gauss elimination
    if (THIS_IMAGE() .EQ. 1) then
        call CPU_TIME(before)
        call gauss_seq(A,B(1,:),N-1)
        call CPU_TIME(after)
        
        write(19,*)"gauss_seq:  ", NUM_IMAGES(), N, (after-before)
    end if
    
    syncall()
    
    
    
    

end program main
    
    
    
    