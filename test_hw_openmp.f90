program test_hw_openmp 
  ! Use iso_fortran_env module for portable 
  ! declaration of floating point double 
  ! precision kind size parameter DP and 
  ! standard error and output unit files
  ! stderr and stdout, resp.
  use, intrinsic :: iso_fortran_env
  !$ use omp_lib
  
  implicit none
  integer, parameter :: DP     = REAL64
  integer, parameter :: stdout = OUTPUT_UNIT
  integer, parameter :: stderr = ERROR_UNIT
  logical, parameter :: DEBUG  = .true.
  ! Size of A, x and y arrays
  integer, parameter   :: N    = 10000 

  ! Number of times matrix-vector product is repeated
  ! for reproducibility-related reasons
  integer, parameter   :: R    = 10

  ! Declare allocatable arrays storing the entries of A, x, y 
  real(DP), allocatable :: A(:,:)
  real(DP), allocatable :: x(:)
  real(DP), allocatable :: y(:) 

  ! Rest of miscellaneous program variables  
  double precision :: t_start, t_stop, t_current, t_min
  integer          :: i, j

  ! Allocate arrays
  allocate(A(N,N))
  allocate(x(N)) 
  allocate(y(N))
   
  do j=1,R
    call initialize(A,x,y)

    t_start = 0.0d0
    !$ t_start   = omp_get_wtime()
    call matvec(A,x,y)
    t_stop  = 0.0d0
    !$ t_stop    = omp_get_wtime()
    t_current = t_stop - t_start 
    if (DEBUG) write(stdout,'(A,I3,A,E32.25,A,E32.25)') 'Repetition #'          , j, & 
                                                        ' Elapsed Time (secs.) ', t_current, &
                                                        ' Relative Error'       , compute_relative_error(y) 
    if (j==1 .or. t_current < t_min) then
      t_min = t_current
    end if 
  end do  
  write(stdout,'(A,E32.25,A,E32.25)') ' Elapsed Time (secs.)', t_min, &
                                      ' Relative Error'      , compute_relative_error(y)

  ! Deallocate arrays
  deallocate(A)
  deallocate(x) 
  deallocate(y)
contains

  subroutine matvec(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   
   write (*,'(a)') '! *** Master on Numerical Methods in Engineering (MNME) ***'
   write (*,'(a)') '! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***'        
   write (*,'(a)') '! *** OpenMP HomeWork'
   write (*,'(a)') '! *** CODE DEVELOPMENT REQUIRED ***'        
   write (*,'(a)') '!     THE BODY OF THE FOLLOWING'        
   write (*,'(a)') '!     SUBROUTINE MUST BE DEVELOPED'        
   write (*,'(a)') '!     FOR THE OpenMP HW ASSIGNMENT'
   stop

  end subroutine matvec

  subroutine initialize(A,x,y)
    implicit none
    real(DP), intent(inout) :: A(:,:)
    real(DP), intent(inout) :: x(:)
    real(DP), intent(inout) :: y(:)
    integer :: i,j   
 
    ! Check size compatibility (if DEBUG parameter is .true.)
    if (DEBUG) call check_size_compatibility(A,x,y)

   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO schedule(runtime)
   do i=1,size(A,1)
     do j=1,size(A,2)
       A(i,j) = real(i+j,DP) 
     end do
     x(i) = 1.0_DP
     y(i) = 1.0_DP 
   end do
   !$OMP END DO 
   !$OMP END PARALLEL 
  end subroutine initialize

  subroutine check_size_compatibility(A,x,y)
    implicit none
    real(DP), intent(in) :: A(:,:)
    real(DP), intent(in) :: x(:)
    real(DP), intent(in) :: y(:)
    logical :: compatible_sizes
    compatible_sizes = ((size(A,1) == size(A,2)) .and. &
                        (size(x)   == size(A,2)) .and. &
                        (size(y)   == size(A,1))) 
    if (.not. compatible_sizes) then
      write(stderr,'(A)') 'check_size_compatibility::sizes of [A,x,y] triplet provided not compatible'
      stop 
    endif  
  end subroutine check_size_compatibility
 
  ! Computes the relative error ||y_exact-y_computed||_2/||y_exact||_2
  ! IMPORTANT NOTE: this function assumes that A,x,y have been
  !                 initialized as A_ij = i+j, x_i = 1, y_i =1
  !                 within the initialize subroutine above. 
  function compute_relative_error(y) 
    implicit none
    real(DP), intent(in) :: y(:)
    real(DP) :: compute_relative_error

    real(DP) :: nrm2_y_exact
    real(DP) :: nrm2_y_exact_sub_y_computed
    real(DP) :: aux
    integer :: i,j 

    nrm2_y_exact = 0.0_DP
    nrm2_y_exact_sub_y_computed = 0.0_DP

    do i=1,size(y)
      aux = 1.0_DP
      do j=1,size(y)
        aux = aux + real(i+j,DP)
      end do
      nrm2_y_exact = nrm2_y_exact + aux*aux
      nrm2_y_exact_sub_y_computed = nrm2_y_exact_sub_y_computed + & 
                                      (aux-y(i))*(aux-y(i))
   end do    
   nrm2_y_exact                = sqrt(nrm2_y_exact)
   nrm2_y_exact_sub_y_computed = sqrt(nrm2_y_exact_sub_y_computed) 
   compute_relative_error      = nrm2_y_exact_sub_y_computed/nrm2_y_exact 
  end function compute_relative_error  

end program test_hw_openmp
