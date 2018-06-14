program pvals2chains 
  implicit none 
  integer, parameter :: ndim_pso=6, ndim_chains=38,n_max=80000 
  real, dimension(:), allocatable ::  y
  real, dimension(:,:), allocatable ::  x1, x2
  real, dimension(ndim_pso+1) :: tt 
  character (len=240)inputfile,str
  integer :: i,j,k, NT, N, ios    
  real :: h 
  CHARACTER(LEN=1) :: junk 
   
   if (iargc() .lt. 1) then
      write(*,*)"ERROR ! Plase give use as follows"
    !  write(*,*)"cat | awk '{printf(/"%s %s %s %s %s %s  %s\n",$3,$4,$5,$6,$7,$8,$9)}' | sort -n -r  > pvals1.dat"
      write(*,*)"./pvals2chain pvals1.dat > test_1.txt"
      stop
   end if
   

   CALL GETARG(1, STR)
   read (STR,'(A)')inputfile
   
   open(unit=1,file=trim(inputfile))
   
   DO 100 i = 1,n_max     ! We know there are fewer than 10000 lines  
     READ(1,*,END=200) (tt(j),j=1,ndim_pso+1)
100  CONTINUE
200  N  = i - 1     
     
     REWIND(1)
     
     allocate(Y(N))
     allocate(X1(N,NDIM_PSO))
     allocate(X2(N,NDIM_CHAINS))
     
     do i = 1, N 
      read(1,*)y(i),(x1(i,j),j=1,ndim_pso)
   end do
   close(1)
   
   do i = 1, N  
      do j = 1, ndim_chains 
         x2(i,j) = 0.0
      end do
 end do
 
 do i = 1, N 
   x2(i,1) = x1(i,1) ! ombh2
   x2(i,2) = x1(i,2) ! omch2
   x2(i,7) = x1(i,3) ! H0 
   x2(i,18) = x1(i,4)!  10^9 A_s 
   x2(i,6) = x1(i,5) ! ns 
   x2(i,4) = x1(i,6) ! tau   
   h =  x1(i,3) /100.0 
   x2(i,8) = 1.0- (x1(i,1)+x1(i,2))/(h*h)        ! omegal  
   x2(i,5) =log(10.0* x1(i,4))                   ! log[10^10As]
end do

do i = 1, N 
   write(*,'(41(E14.7))')real(i),y(i),(x2(i,j),j=1,ndim_chains)
end do

end program pvals2chains 
