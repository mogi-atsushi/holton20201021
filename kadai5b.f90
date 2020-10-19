program kadai5b
  implicit none
  
  ! kind parameters
  integer(4), parameter :: kindi=4
  integer(4), parameter :: kindr=8
  
  ! model settings
  integer(kindi), parameter :: jmax=10
  real(kindr), parameter :: dt=1.0e-2_kindr
  real(kindr), parameter :: nu=1.0_kindr, sigma=0.1_kindr
  real(kindr), parameter :: time_max=0.1_kindr
  real(kindr), parameter :: time_int=0.01_kindr ! output interval
  character(len=*), parameter :: outputname='kadai5b.csv'
  character(len=2) :: scheme
  
  ! work
  integer(kindi), parameter :: nmax=int(time_max/dt,kind=kindi)
  integer(kindi), parameter :: nint=int(time_int/dt,kind=kindi)
  real(kindr), parameter :: dx=1.0_kindr/real(jmax,kind=kindr)
  real(kindr), parameter :: mu=nu*dt*real(jmax,kind=kindr)**2
  real(kindr), dimension(-jmax:jmax,0:nmax) :: temperature
  real(kindr), dimension(-jmax:jmax) :: x
  real(kindr), dimension(0:nmax) :: time

  ! if (jmax==5) then
  !
  ! x -1                   0                  +1
  !    *---*---*---*---*---*---*---*---*---*---*
  ! j -5  -4  -3  -2  -1   0   1   2   3   4   5

  ! --- start --- !
  write(*,'("   nu=",e12.6)') nu
  write(*,'("sigma=",e12.6)') sigma
  write(*,'("   dt=",e12.6)') dt
  write(*,'("   dx=",e12.6)') dx
  write(*,'("   mu=",e12.6)') mu

  ! initialization
  call init(x, time, temperature)

  ! time integration
  write(*,'("Which scheme? (f, b, cn): ")',advance='no')
  read(*,*) scheme
  if (scheme=='f') then
     call forward(temperature)
  else if (scheme=='b') then
     call backward(temperature)
  else if (scheme=='cn') then
     call cn(temperature)
  else
     write(*,*) 'Please answer f, b, or cn.'
     stop
  end if

  ! output
  call output(x, time, temperature)
  
  stop
contains
  
  subroutine init(x, time, temperature)
    implicit none
    real(kindr), dimension(-jmax:jmax,0:nmax), intent(inout) :: temperature
    real(kindr), dimension(-jmax:jmax), intent(inout) :: x
    real(kindr), dimension(0:nmax), intent(inout) :: time
    ! work
    integer (kindi) :: j,n
    ! --- start --- !
    do j = -jmax, jmax
       x(j) = real(j,kind=kindr) / real(jmax,kind=kindr)
    end do
    do n = 0, nmax
       time(n) = dt*n
    end do
    ! initial condition
    do j = -jmax, jmax
       temperature(j,0) = exp(-1.0_kindr*x(j)**2 / (2.0_kindr*sigma**2))
    end do
  end subroutine init

  subroutine forward(temperature)
    implicit none
    real(kindr), dimension(-jmax:jmax,0:nmax), intent(inout) :: temperature
    ! work
    integer(kindi) :: j,n
    ! --- start --- !
    do n = 0, nmax-1
       do j = -jmax+1, jmax-1
          temperature(j,n+1) = temperature(j,n) &
               & + mu*(temperature(j+1,n)-2.0_kindr*temperature(j,n)+temperature(j-1,n))
       end do
       ! boundary condition
       temperature(-jmax,n+1) = temperature(-jmax+1,n+1)
       temperature(+jmax,n+1) = temperature(+jmax-1,n+1)
    end do
  end subroutine forward

  subroutine inverse(A, dimA)
    implicit none
    integer(kindi), intent(in) :: dimA
    real(kindr), dimension(dimA,dimA), intent(inout) :: A
    ! work
    integer(kindi) :: j,info
    integer(kindi), dimension(dimA) :: ipiv
    real(kindr), dimension(dimA) :: work
    ! --- start --- !
    ! LU decomposition
    call dgetrf(dimA, dimA, A, dimA, ipiv, info)
    if (info/=0) then
       write(*,*) 'Error (LU decomposition)', info
       stop
    end if
    ! inverse matrix
    call dgetri(dimA, A, dimA, ipiv, work, dimA, info)
    if (info/=0) then
       write(*,*) 'Error (inverse matrix)', info
       stop
    end if
  end subroutine inverse
  
  subroutine backward(temperature)
    implicit none
    real(kindr), dimension(-jmax:jmax,0:nmax), intent(inout) :: temperature
    ! work
    integer(kindi) :: j,n
    integer(kindi), parameter :: dimA=2*jmax-1
    real(kindr) :: diagonal
    real(kindr), dimension(dimA,dimA) :: A
    ! --- start --- !
    ! T(n) = A T(n+1)
    A = 0.0_kindr
    diagonal = 1.0_kindr + 2.0_kindr * mu
    do j = 2, dimA-1
       A(j,j-1) = -mu
       A(j,j) = diagonal
       A(j,j+1) = -mu
    end do
    A(1,1) = 1.0_kindr + mu
    A(1,2) = -mu
    A(dimA,dimA) = A(1,1)
    A(dimA,dimA-1) = A(1,2)
    ! T(n+1) = A T(n)
    call inverse(A ,dimA)
    ! time integration
    do n = 0, nmax-1
       temperature(-jmax+1:jmax-1,n+1) = matmul(A, temperature(-jmax+1:jmax-1,n))
    end do
    temperature(-jmax,:) = temperature(-jmax+1,:)
    temperature(+jmax,:) = temperature(+jmax-1,:)
  end subroutine backward
  
  subroutine cn(temperature) ! Crank-Nicolson scheme
    implicit none
    real(kindr), dimension(-jmax:jmax,0:nmax), intent(inout) :: temperature
    ! work
    integer(kindi) :: j,n
    integer(kindi), parameter :: dimA=2*jmax-1
    real(kindr) :: diagonal, halfmu=mu*0.5_kindr
    real(kindr), dimension(dimA,dimA) :: A, B
    ! --- start --- !
    ! A T(n+1) = B T(n)
    ! --- A --- !
    A = 0.0_kindr
    diagonal = 1.0_kindr + mu
    do j = 2, dimA-1
       A(j,j-1) = -halfmu
       A(j,j) = diagonal
       A(j,j+1) = -halfmu
    end do
    A(1,1) = 1.0_kindr + halfmu
    A(1,2) = -halfmu
    A(dimA,dimA) = A(1,1)
    A(dimA,dimA-1) = A(1,2)
    ! --- B --- !
    B = 0.0_kindr
    diagonal = 1.0_kindr - mu
    do j = 2, dimA-1
       B(j,j-1) = halfmu
       B(j,j) = diagonal
       B(j,j+1) = halfmu
    end do
    B(1,1) = 1.0_kindr - halfmu
    B(1,2) = halfmu
    B(dimA,dimA) = B(1,1)
    B(dimA,dimA-1) = B(1,2)
    ! T(n+1) = A B T(n)
    call inverse(A, dimA)
    ! T(n+1) = A T(n)
    A = matmul(A, B)
    ! time integration
    do n = 0, nmax-1
       temperature(-jmax+1:jmax-1,n+1) = matmul(A, temperature(-jmax+1:jmax-1,n))
    end do
    temperature(-jmax,:) = temperature(-jmax+1,:)
    temperature(+jmax,:) = temperature(+jmax-1,:)
  end subroutine cn

  subroutine output(x, time, temperature)
    implicit none
    real(kindr), dimension(-jmax:jmax,0:nmax), intent(in) :: temperature
    real(kindr), dimension(-jmax:jmax), intent(in) :: x
    real(kindr), dimension(0:nmax), intent(in) :: time
    ! work
    integer(kindi) :: j,n
    ! --- start --- !
    open(10, file=outputname, action='write', form='formatted', status='replace')
    do n = 0, nmax, nint
       write(10,'(",",f5.2)',advance='no') time(n)
    end do
    write(10,*)
    do j = -jmax, jmax
       write(10,'(f5.2)',advance='no') x(j)
       do n = 0, nmax, nint
          write(10,'(",",e12.6)',advance='no') temperature(j,n)
       end do
       write(10,*)
    end do
    close(10)
  end subroutine output
end program kadai5b
