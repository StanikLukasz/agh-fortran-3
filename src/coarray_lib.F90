!------------------------------------------------------------------------------
! MODULE: coarray_lib
!
!> @author
!> Lukasz Stanik
!
! DESCRIPTION: 
!> The module provides implementations
!> of some mathematical functions using coarray.
!
! REVISION HISTORY:
! 23 06 2018 - Initial version
!------------------------------------------------------------------------------
    
module coarray_lib
    
    contains
  !------------------------------------------------------------------------------
  !> @author
  !> Lukasz Stanik
  !
  ! DESCRIPTION: 
  !> Multiplies two matrices and returns result in the `X` param.
  !> Works with coarray to boost computings.
  !
  !> @param[in] A First matrix
  !> @param[in] B Second matrix
  !> @param[out] X Result matrix
  !> @param[out] status Status, 0 = success
  !------------------------------------------------------------------------------
    subroutine mult_coarr(A, B, X, status)
     implicit none
        
        ! variables
        real (kind = 8), intent(in) :: A(:,:) 
        real (kind = 8), intent(in) :: B(:,:) 
        real (kind = 8), intent(out) :: X(:,:)
        
        real (kind = 8), codimension[:], dimension(:,:), allocatable :: buffer
        integer (kind = 4) :: width, image
        integer (kind = 4), codimension[:], allocatable :: first_row, last_row
            
        integer (kind = 4), intent(out) :: status
        integer (kind = 4) :: Ar, Ac, Br, Bc, Xr, Xc
        integer (kind = 4) :: i, j, k
        
        logical :: dim_condition
        
        allocate(first_row[*])
        allocate( last_row[*])

        ! checking dimensions condition
        Ar = size(A(:,1))
        Ac = size(A(1,:))
        Br = size(B(:,1))
        Bc = size(B(1,:))
        Xr = size(X(:,1))
        Xc = size(X(1,:))
        
        ! initializing coarray variables
          
        width = CEILING(real(Ar)/NUM_IMAGES())
        allocate(buffer(width, Xc)[*])
        first_row = MIN(Ar, (THIS_IMAGE() - 1) * width) + 1
        last_row = MIN(Ar, THIS_IMAGE() * width)
        
        dim_condition = (Ac .EQ. Br) .AND. (Ar .EQ. Xr) .AND. (Bc .EQ. Xc)
        
        ! returning with status = -1 if dim_condition is not true
        if (.NOT. dim_condition) then
            status = -1
            return
        end if
        
        ! multiplying itself     
        do i = first_row, last_row
            do j = 1, Xc         
                    buffer(i - first_row + 1,j) = 0.d0
                    do k = 1, Ac
                        buffer(i - first_row + 1, j) = buffer(i - first_row + 1,j) + A(i,k) * B(k,j)
                    end do
            end do
        end do

        syncall()
        
        ! getting results together
        if (THIS_IMAGE() .EQ. 1) then
            do image = 1, NUM_IMAGES()
                X(first_row[image]:last_row[image],:) = buffer(1:(last_row[image] - first_row[image] + 1),:)[image]
            end do
        end if
        
        deallocate(buffer)
        deallocate(first_row)
        deallocate(last_row)
        
        status = 0
        
    end subroutine mult_coarr
    
  !------------------------------------------------------------------------------
  !> @author
  !> Lukasz Stanik
  !
  ! DESCRIPTION: 
  !> Provides a gauss elimination mechanism.
  !> Works with coarray to boost computings.
  !
  !> @param[inout] A Matrix of coefficients
  !> @param[inout] X Matrix of values
  !> @param[in] N Max row number of the from-0-numbered matrices
  !------------------------------------------------------------------------------
    subroutine gauss_coarr (A, X, N)
        implicit none
        
        
        ! variables
        integer (kind = 8), intent(in) :: N
        real (kind = 8), intent(inout) :: A(N,N), X(N)
        real (kind = 8), codimension[:], allocatable :: cA(:,:), cX(:)
        
        integer (kind = 8) :: I, J
        real (kind = 8) :: C
        
        allocate(cA(0:N,0:N)[*])
        allocate(cX(0:N)[*])

        if (THIS_IMAGE() .EQ. 1) then
            cA(:,:)[1] = A(:,:)
            cX(:)[1] = X(:)
        end if
        
        
        ! algorithm
        do I = 1,N-1
            do J = THIS_IMAGE()-1, N-1, NUM_IMAGES()
                if (I .NE. J) then
                    C = cA(I, J+1) / cA(I, I+1)
                    cA(:, J+1) = cA(:, J+1) - C*cA(:,I+1)
                    cX(J+1) = cX(J+1) - C*cX(I+1)
                end if
            end do
        end do
        
        ! getting results together
        if (THIS_IMAGE() .EQ. 1) then
            A(:,:) = cA(:,:)[1]
            X(:) = cX(:)[1]
        end if
        
    end subroutine gauss_coarr
    
    
end module coarray_lib