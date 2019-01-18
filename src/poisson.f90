program poisson
use omp_lib
implicit none

  integer,parameter :: rk=4, N=100, t=1000
  real(kind=rk), parameter :: pi = 4.d0*atan(1.d0)
  integer :: i,j,tt,nt,ia
  real(kind=rk) :: x,y,t1,t2,wt0,wt1, gama
  real(kind=rk),dimension(N,N) :: f, fnew, freal, g
  character(len=10) :: arg

  call getarg(1,arg)
  read(arg,*) gama

  f=0
  fnew=0
  
  !creating g
  do i = 1,N
    do j = 1,N
      x = real(i)/N
      y = real(j)/N
      g(i,j) = pi*exp(x+y)*(cos(pi/2*exp(x+y))-pi/2*exp(x+y)*sin(pi/2*exp(x+y))) & 
              -pi*exp(y-x)*(sin(pi/2*exp(y-x))-pi/2*exp(y-x)*cos(pi/2*exp(y-x)))
    end do
  end do
  
  !print '(10(f6.2))', transpose(g)
  
  !computing real f   
  do i = 1,N
    do j = 1,N
      x = real(i)/N
      y = real(j)/N
      freal(i,j) = cos(pi/2*exp(y-x))+sin(pi/2*exp(y-x))
    end do
  end do
  
  !boundary conditions
  f(1,1:N) = freal(1,1:N)
  f(N,1:N) = freal(N,1:N)
  f(1:N,1) = freal(1:N,1)
  f(1:N,N) = freal(1:N,N)
  
  fnew = f
  
  
  !print '(10(f6.2))', transpose(freal)
  !print *, '-------------------------------------------------------'
  
  !computing f with JOR without parallelization and measuring time it takes
  call cpu_time(t1)
  do tt = 1,t
    do i = 2,N-1
      do j = 2,N-1
        x = real(i)/N
        y = real(j)/N
        fnew(i,j) = (1-gama)*f(i,j)+gama/4*(f(i+1,j)+f(i,j+1)+f(i-1,j)+f(i,j-1)) & 
                    -(gama/(4*N**2)*g(i,j))
      end do
    end do
    f = fnew
  end do
  
  call cpu_time(t2)
  print *,'JOR time without parallelization: ',t2-t1
  
  
  !computing f with SOR without parallelization and measuring time it takes
  call cpu_time(t1)
  do tt = 1,t
    
    !Black points
    do i = 2,N-1,2
      do j = 2,N-1,2
        fnew(i,j) = (1-gama)*f(i,j)+gama/4*(f(i+1,j)+f(i,j+1)+f(i-1,j)+f(i,j-1)) & 
                    -(gama/(4*N**2)*g(i,j))
      end do
    end do
    
    !Red points
    do i = 3,N-1,2
      do j = 3,N-1,2
        fnew(i,j) = (1-gama)*f(i,j)+gama/4*(fnew(i+1,j)+fnew(i,j+1)+fnew(i-1,j)+ &
                     fnew(i,j-1)) -(gama/(4*N**2)*g(i,j))
      end do
    end do   
    f = fnew
    
  end do    
  
  call cpu_time(t2)
  print *,'SOR time without parallelization: ',t2-t1
  print *,'difference: ',meandiff(freal,f,N)
  
  !computing f with parallelization using SOR
  !$call omp_set_num_threads(4)
  !$setenv OMP_STACKSIZE 10M
  wt0=omp_get_wtime()
  
  !$omp parallel
  
    !$omp master
    !$ print *, "Running with ", omp_get_num_threads(), " threads"
    !$omp end master
  
    do tt = 1,t
   
      !Black points
      !$omp do collapse(2)
      do i = 2,N-1,2
        do j = 2,N-1,2
          fnew(i,j) = (1-gama)*f(i,j)+gama/4*(f(i+1,j)+f(i,j+1)+f(i-1,j)+f(i,j-1)) & 
                      -(gama/(4*N**2)*g(i,j))
        end do
      end do
      !$omp end do
      
      !Red points
      !$omp do collapse(2)
      do i = 3,N-1,2
        do j = 3,N-1,2
          fnew(i,j) = (1-gama)*f(i,j)+gama/4*(fnew(i+1,j)+fnew(i,j+1)+fnew(i-1,j)+ &
                       fnew(i,j-1)) -(gama/(4*N**2)*g(i,j))
        end do
      end do
      !$omp end do
      f = fnew
    
    end do
    
  !$omp end parallel
  
   wt1=omp_get_wtime()
   print *,'wall clock time using SOR with omp: ',wt1-wt0
  
contains

  function meandiff(a,b,n) result (y)
    real(kind=4), dimension(n,n), intent(in) :: a,b
    integer, intent(in) :: n
    integer :: i,j
    real(kind=4) :: y
    real(kind=4) :: dif(n,n)
    
    do i = 1,n
      do j= 1,n
        dif(i,j) = a(i,j) -b(i,j)
      end do
    end do
    
    y = sum(dif)/(n*n)
  end function  

  
end program
  
  
  
  
  




