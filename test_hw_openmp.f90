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
  logical, parameter :: DEBUG  = .false.
  integer, parameter :: out_unit=20
  character(len=1024) :: out_name
  CHARACTER(LEN=30) :: rowfmt,in_itefmt,im_itefmt
  integer, parameter :: N_div = 15
  integer :: num_threads, i_num_threads
  
  ! Size of A, x and y arrays
  integer, parameter   :: N    = 2**(N_div-1)

  ! Number of times matrix-vector product is repeated
  ! for reproducibility-related reasons
  integer, parameter   :: R    = 10

  ! Declare allocatable arrays storing the entries of A, x, y 
  real(DP), allocatable :: A(:,:)
  real(DP), allocatable :: x(:)
  real(DP), allocatable :: y(:)
  !real(DP), allocatable :: y_average(:) 

  ! Rest of miscellaneous program variables  
  double precision :: t_start, t_stop, t_current, t_average, e_current, e_average
  integer          :: i, j, n0, i_n, i_m, i2
  integer, dimension(1:N_div):: n0_vec

  !$OMP PARALLEL
  !$ num_threads = omp_get_num_threads()
  !$OMP END PARALLEL
  i_num_threads=int(log(real(num_threads)) / log(2.)+1.0)
  write (out_name, "(A,I2.2,A)") "results", num_threads , ".m"
  
  do n0=1,N_div
	n0_vec(n0)=2**(n0-1)
  end do

open (unit=out_unit,file=out_name,action="write",status="replace")
  
  WRITE(in_itefmt,'(A,I10,A)') '(A,',N_div,'(1X,I10),A)'
  WRITE(im_itefmt,'(A,I10,A)') '(A,',6,'(1X,I2),A)'
do i_n=1,N_div
  n0=n0_vec(i_n)
  write(*,*)'N = ',n0
  !Writting N
  write(out_unit,*)
  write(out_unit,*)'%#### N = ',n0,'#############'
  ! Allocate arrays
  allocate(A(n0,n0))
  allocate(x(n0)) 
  allocate(y(n0))
!  allocate(y_average(n0))
  
  

  WRITE(rowfmt,'(A,I10,A)') '(',n0,'(1X,D20.15),A)'
  do i_m=1,9
  write(*,*)'Method = ',i_m
  write(out_unit,*)
  write(out_unit,*)'%######## Method = ',i_m,'#########'
  write(out_unit,*)
  t_average = 0.0d0
  e_average = 0.0d0
!   !$OMP PARALLEL default(none) shared(y_average,n0,num_threads) private(i)
!   !$OMP DO schedule(runtime)
!  do i=1,n0
!  y_average(i) = 0.0d0
!  end do
!   !$OMP END DO
!   !$OMP END PARALLEL
  do j=1,R
    call initialize(A,x,y)
	!Write(*,*) "Calculating y=y+A*x"
	call matvec_type(A,x,y,i_m,t_start, t_stop)
    t_current = t_stop - t_start 
	!Write(*,*) "Calculating error norm"
	e_current=compute_relative_error(y) 
    if (DEBUG) write(stdout,'(A,I3,A,E32.25,A,E32.25)') 'Repetition #'          , j, & 
                                                        ' Elapsed Time (secs.) ', t_current, &
                                                        ' Relative Error'       , e_current
    !if (j==1 .or. t_current < t_min) then
    !  t_min = t_current
    !end if
	!WRITE(*,FMT=rowfmt)(y(i), i=1,n0)
!	!$OMP PARALLEL default(none) shared(y_average,n0,y) private(i)
!	!$OMP DO schedule(runtime)
!	do i=1,n0
!	y_average(i)=y_average(i)+y(i)/R
!	end do
!	!$OMP END DO
!	!$OMP END PARALLEL
	t_average=t_average+t_current/R
	e_average=e_average+e_current/R
  end do  
  write(stdout,'(A,E32.25,A,E32.25)') ' Elapsed Time (secs.)', t_average, &
                                      ' Relative Error'      , e_average
									  
  !Writting E_current & T_current
  WRITE(out_unit,"(A,I10,A,I2,A,I2,A)")'E{',i_n,',',i_m,',',i_num_threads,'}=['
  WRITE(out_unit,"(E32.25,A)")e_average,'];'
  
  WRITE(out_unit,"(A,I10,A,I2,A,I2,A)")'T{',i_n,',',i_m,',',i_num_threads,'}=['
  WRITE(out_unit,"(E32.25,A)")t_average,'];'
  
!  !Writting Y
!  WRITE(out_unit,*)'Y{',i_n,',',i_m,',',i_num_threads,'}=['
!  WRITE(out_unit,FMT=rowfmt)(y_average(j), j=1,n0),'];'
  end do
!  !Writting X
!  WRITE(out_unit,*)'X{',i_n,'}=['
!  WRITE(out_unit,FMT=rowfmt)(x(j), j=1,n0),'];'
!  
!  !Writting A
!  WRITE(out_unit,*)'A{',i_n,'}=['
!  do i=1,n0
!  if (i==n0) then 
!  WRITE(out_unit,FMT=rowfmt)(A(i,j), j=1,n0),'];' 
!  else 
!  WRITE(out_unit,FMT=rowfmt)(A(i,j), j=1,n0),';...' 
!  endif
!  end do
  
  ! Deallocate arrays
  deallocate(A)
  deallocate(x) 
  deallocate(y)
!  deallocate(y_average)
end do
!Writting Array Iterators:
  WRITE(out_unit,FMT=in_itefmt)'i_n=[',n0_vec,'];'
  WRITE(out_unit,FMT=im_itefmt)'i_m=[',1,2,3,4,5,6,'];'

close(out_unit)
  
contains

  subroutine matvec_type(A,x,y,method,t_start, t_stop)
	implicit none
	real(DP), intent(in)     		:: A(:,:)
	real(DP), intent(in)     		:: x(:)
	real(DP), intent(inout)  		:: y(:)
	integer, intent(in) 	 		:: method
	double precision , intent(out) 	:: t_start, t_stop
	
	if (method==1) then
	
		t_start = 0.0d0
		!$ t_start   = omp_get_wtime()
		call matvec1(A,x,y)
		t_stop  = 0.0d0
		!$ t_stop    = omp_get_wtime()
		
	elseif (method==2) then
  
  		t_start = 0.0d0
		!$ t_start   = omp_get_wtime()
		call matvec2(A,x,y)
		t_stop  = 0.0d0
		!$ t_stop    = omp_get_wtime()
		
	elseif (method==3) then
  
  		t_start = 0.0d0
		!$ t_start   = omp_get_wtime()
		call matvec3(A,x,y)
		t_stop  = 0.0d0
		!$ t_stop    = omp_get_wtime()
		
	elseif (method==4) then
  
  		t_start = 0.0d0
		!$ t_start   = omp_get_wtime()
		call matvec4(A,x,y)
		t_stop  = 0.0d0
		!$ t_stop    = omp_get_wtime()
		
	elseif (method==5) then
  
  		t_start = 0.0d0
		!$ t_start   = omp_get_wtime()
		call matvec5(A,x,y)
		t_stop  = 0.0d0
		!$ t_stop    = omp_get_wtime()
		
	elseif (method==6) then
  
  		t_start = 0.0d0
		!$ t_start   = omp_get_wtime()
		call matvec6(A,x,y)
		t_stop  = 0.0d0
		!$ t_stop    = omp_get_wtime()
		
	elseif (method==7) then
  
  		t_start = 0.0d0
		!$ t_start   = omp_get_wtime()
		call matvec6(A,x,y)
		t_stop  = 0.0d0
		!$ t_stop    = omp_get_wtime()
		
	elseif (method==8) then
  
  		t_start = 0.0d0
		!$ t_start   = omp_get_wtime()
		call matvec6(A,x,y)
		t_stop  = 0.0d0
		!$ t_stop    = omp_get_wtime()
		
	elseif (method==9) then
  
  		t_start = 0.0d0
		!$ t_start   = omp_get_wtime()
		call matvec6(A,x,y)
		t_stop  = 0.0d0
		!$ t_stop    = omp_get_wtime()
  
	endif

  end subroutine matvec_type

  subroutine matvec1(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   !write(stdout,*) 'Debug_variable: ',DEBUG
   
   !axpy
   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO schedule(runtime)
   do i=1,size(A,1) !COLUMN
     do j=1,size(A,2) !ROW
       y(i) = y(i) + A(i,j)*x(j)
     end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL
   
   !write (*,'(a)') '! *** Master on Numerical Methods in Engineering (MNME) ***'
   !write (*,'(a)') '! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***'        
   !write (*,'(a)') '! *** OpenMP HomeWork'
   !write (*,'(a)') '! *** CODE DEVELOPMENT REQUIRED ***'        
   !write (*,'(a)') '!     THE BODY OF THE FOLLOWING'        
   !write (*,'(a)') '!     SUBROUTINE MUST BE DEVELOPED'        
   !write (*,'(a)') '!     FOR THE OpenMP HW ASSIGNMENT'
   !stop

  end subroutine matvec1

  
  subroutine matvec2(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   !write(stdout,*) 'Debug_variable: ',DEBUG

   !axpy
   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO  schedule(runtime)
   do i=1,size(A,1) !COLUMN
     do j=1,size(A,2) !ROW
       y(j) = y(j) + A(j,i)*x(i)
     end do
	 
   end do
   !$OMP END DO
   !$OMP END PARALLEL

  end subroutine matvec2
  
  
   subroutine matvec3(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
   real(DP) :: s
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   !write(stdout,*) 'Debug_variable: ',DEBUG
   
   !axpy
   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO REDUCTION(+:y) schedule(runtime)
   do i=1,size(A,1) !COLUMN
     do j=1,size(A,2) !ROW
       y(j) = y(j) + A(j,i)*x(i)
     end do
	 
   end do
   !$OMP END DO
   !$OMP END PARALLEL

  end subroutine matvec3

  
  subroutine matvec4(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   !write(stdout,*) 'Debug_variable: ',DEBUG
   
   !axpy
   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO schedule(DYNAMIC)
   do i=1,size(A,1) !COLUMN
     do j=1,size(A,2) !ROW
       y(i) = y(i) + A(i,j)*x(j)
     end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL
   
  end subroutine matvec4

  
  subroutine matvec5(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   !write(stdout,*) 'Debug_variable: ',DEBUG

   !axpy
   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO  schedule(DYNAMIC)
   do i=1,size(A,1) !COLUMN
     do j=1,size(A,2) !ROW
       y(j) = y(j) + A(j,i)*x(i)
     end do
	 
   end do
   !$OMP END DO
   !$OMP END PARALLEL

  end subroutine matvec5
  
  
   subroutine matvec6(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
   real(DP) :: s
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   !write(stdout,*) 'Debug_variable: ',DEBUG
   
   !axpy
   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO REDUCTION(+:y) schedule(DYNAMIC)
   do i=1,size(A,1) !COLUMN
     do j=1,size(A,2) !ROW
       y(j) = y(j) + A(j,i)*x(i)
     end do
	 
   end do
   !$OMP END DO
   !$OMP END PARALLEL

  end subroutine matvec6
  
  
  subroutine matvec7(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   !write(stdout,*) 'Debug_variable: ',DEBUG
   
   !axpy
   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO schedule(GUIDED,8192)
   do i=1,size(A,1) !COLUMN
     do j=1,size(A,2) !ROW
       y(i) = y(i) + A(i,j)*x(j)
     end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL
   
  end subroutine matvec7

  
  subroutine matvec8(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   !write(stdout,*) 'Debug_variable: ',DEBUG

   !axpy
   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO  schedule(GUIDED,8192)
   do i=1,size(A,1) !COLUMN
     do j=1,size(A,2) !ROW
       y(j) = y(j) + A(j,i)*x(i)
     end do
	 
   end do
   !$OMP END DO
   !$OMP END PARALLEL

  end subroutine matvec8
  
  
   subroutine matvec9(A,x,y)
   implicit none
   real(DP), intent(in)     :: A(:,:)
   real(DP), intent(in)     :: x(:)
   real(DP), intent(inout)  :: y(:)
   integer :: i, j
   real(DP) :: s
  
   ! Check size compatibility (if DEBUG parameter is .true.)
   if (DEBUG) call check_size_compatibility(A,x,y)
   !write(stdout,*) 'Debug_variable: ',DEBUG
   
   !axpy
   !$OMP PARALLEL default(none) shared(A,x,y) private(i,j)
   !$OMP DO REDUCTION(+:y) schedule(GUIDED,8192)
   do i=1,size(A,1) !COLUMN
     do j=1,size(A,2) !ROW
       y(j) = y(j) + A(j,i)*x(i)
     end do
	 
   end do
   !$OMP END DO
   !$OMP END PARALLEL

  end subroutine matvec9
  
  
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
