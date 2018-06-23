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
        allocate(buffer(width, Bc)[*])
        
        width = CEILING(real(Ar)/NUM_IMAGES())
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
                    X(i,j) = 0.d0
                    do k = 1, Ac
                        X(i - first_row + 1, j) = X(i - first_row + 1,j) + A(i,k) * B(k,j)
                    end do
            end do
        end do

        sync all
        
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
    
end module coarray_lib