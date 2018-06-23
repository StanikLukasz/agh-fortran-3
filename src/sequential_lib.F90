!------------------------------------------------------------------------------
! MODULE: sequential_lib
!
!> @author
!> Lukasz Stanik
!
! DESCRIPTION: 
!> The module provides sequential implementations
!> of some mathematical functions.
!
! REVISION HISTORY:
! 23 06 2018 - Initial version
!------------------------------------------------------------------------------
    
module sequiential_lib
    contains
  !------------------------------------------------------------------------------
  !> @author
  !> Lukasz Stanik
  !
  ! DESCRIPTION: 
  !> Multiplies two matrices and returns result in the `X` param.
  !
  !> @param[in] A First matrix
  !> @param[in] B Second matrix
  !> @param[out] X Result matrix
  !> @param[out] status Status, 0 = success
  !------------------------------------------------------------------------------
     subroutine mult_seq (A, B, X, status)
        implicit none
        
        ! variables
        real (kind = 8), intent(in) :: A(:,:) 
        real (kind = 8), intent(in) :: B(:,:) 
        real (kind = 8), intent(out) :: X(:,:)
        
        integer (kind = 4), intent(out) :: status
        integer (kind = 4) :: Ar, Ac, Br, Bc, Xr, Xc
        integer (kind = 4) :: i, j, k
        
        
        ! checking dimensions condition
        logical :: dim_condition
        
        Ar = size(A(:,1))
        Ac = size(A(1,:))
        Br = size(B(:,1))
        Bc = size(B(1,:))
        Xr = size(X(:,1))
        Xc = size(X(1,:))
        
        dim_condition = (Ac .EQ. Br) .AND. (Ar .EQ. Xr) .AND. (Bc .EQ. Xc)
        
        ! returning with status = -1 if dim_condition is not true
        if (.NOT. dim_condition) then
            status = -1
            return
        end if
        
        ! multiplying itself     
        do i = 1, Xr
            do j = 1, Xc         
                    X(i,j) = 0.d0
                    do k = 1, Ac
                        X(i,j) = X(i,j) + A(i,k) * B(k,j)
                    end do
            end do
        end do

        status = 0
        
     end subroutine mult_seq
    
 !------------------------------------------------------------------------------
  !> @author
  !> Lukasz Stanik
  !
  ! DESCRIPTION: 
  !> `mult` version with explicit sizes of matrices.
  !> Needed for use in python code.
  !> Accepts only square matrices.
  !>
  !> @param[in] A First matrix
  !> @param[in] B Second matrix
  !> @param[out] X Result matrix
  !> @param[out] status Status, 0 = success
  !> @param[in] size Size of the square matrices
  !------------------------------------------------------------------------------
  subroutine mult_seq_explicit(A, B, X, status, size)
    implicit none
   
    ! variables
    integer ( kind = 4), intent(in) :: size
    real ( kind = 8), intent(in) :: A(size,size)
    real ( kind = 8), intent(in) :: B(size,size)
    real ( kind = 8), intent(out) :: X(size,size)
    integer ( kind = 4), intent(out) :: status
    
    !f2py intent(in) :: A, B, size
    !f2py intent(out) :: X, status
    
    ! calling a mult_seq subroutine
    call mult_seq(A, B, X, status)
 
  end subroutine mult_seq_explicit

 !------------------------------------------------------------------------------
  !> @author
  !> Lukasz Stanik
  !
  ! DESCRIPTION: 
  !> Provides a gauss elimination mechanism.
  !>
  !> @param[inout] A Matrix of coefficients
  !> @param[inout] X Matrix of values
  !> @param[in] N Size of the matrices
  !------------------------------------------------------------------------------
  
   subroutine gauss_seq (A, X, N)
        implicit none
        
        
        ! variables
        integer (kind = 8), intent(in) :: N
        real (kind = 8), intent(inout) :: A(N,N), X(N)
        
        integer (kind = 8) :: I, J
        real (kind = 8) :: C
        
        !f2py intent(in) :: A, X, N
        !f2py intent(out) :: A, X
        
        ! algorithm
        do I = 1,N-1
            do J = 0,N-1
                if (I .NE. J) then
                    C = A(I, J+1) / A(I, I+1)
                    A(:, J+1) = A(:, J+1) - C*A(:,I+1)
                    X(J+1) = X(J+1) - C*X(I+1)
                end if
            end do
        end do
        
    end subroutine gauss_seq
  
end module sequiential_lib
    
    
        