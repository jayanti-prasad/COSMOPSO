!-----------------------------------------------------------------------
!-This is the main progrm. It calls camb module for computing Cls and 
! WMAP likelihood module for computing likelihood function.
!-For sampling it uses PSO module 
!-There is io interface which reads parameters for CAMB and PSO module.
!-The current module has automatic termination criterian in place i.e.,
! if gbest does not improve after a certain number steps, exploration is stopped.
!                                     Jayanti Prasad, July 01, 2011 
!------------------------------------------------------------------------

program cmbrpso_main
  use defs
  use pso
  use set_para
  use wmap_likelihood_9yr
  use wmap_options 
  use wmap_util
  use clik
  implicit none
  include "mpif.h"
  integer :: i, j, l,rank,ierr
  logical :: converge,nstop 
  character (len=MAX_LEN) :: str, paramfile, CAMspec,  commander,lowLike
  
  if (iargc() .lt. 1) then
     write(*,*)"ERROR ! Plase give use as follows"
     write(*,*)"./mpirun -np <np> ./cosmopso param.in"
     stop
 end if
 
 call getarg(1,str)
 
 read (str,'(A)')paramfile
 
 call MPI_INIT(ierr)
 
 call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
 
 call read_para(paramfile)

 
 if (USE_WMAP)  then 
    call wmap_likelihood_init
 end if
 
 if (USE_PLANCK) then 
    CAMspec=trim(planck_data)//'CAMspec_v6.2TN_2013_02_26_dist.clik/'
    commander=trim(planck_data)//'commander_v4.1_lm49.clik/'
    lowLike=trim(planck_data)//'lowlike_v222.clik/'   
 
    call clik_init(pself_commander,trim(commander))
    call clik_init(pself_camspec,trim(CAMspec))
    call clik_init(pself_lowl,trim(lowLike))
    
 end if
 
 call pso_init()
 
 write(*,*)"pso_init done "
 
 nstop=.false. 
 
 do i=1,ntimes_max      
    if (i .ge. 2) then 
       call pso_main(i)
    end if
    
    if (rank .eq. 0) then
       
       call write_output(i)
       
       if (converge(i)) then         
          write(*,*)"convergence found after",i,"steps"
          nstop=.true.
       end if
       
    end if
    
    call  MPI_Bcast(nstop,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    
    if(nstop) go to 10  
 end do
 
10 call MPI_FINALIZE(ierr)
 
 stop 
 
end program cmbrpso_main

subroutine write_output(i)
  use defs 
  implicit none 
  integer, intent(in) :: i 
  integer :: j, k 
  logical :: ex

  inquire(file="pvals.dat", exist=ex)
  if (ex) then
     open(10,file="pvals.dat",status='old',position='append')
  else
     open(10,file="pvals.dat",status='new',position='append')
  endif

  inquire(file="pbest.dat", exist=ex)
  if (ex) then
     open(11,file="pbest.dat",status='old',position='append')
  else
     open(11,file="pbest.dat",status='new',position='append')
  endif

  inquire(file="gbest.dat", exist=ex)
  if (ex) then
     open(12,file="gbest.dat",status='old',position='append')
  else
     open(12,file="gbest.dat",status='new',position='append')
  endif

  do j = 1, npart
      write(10,'(I8,I8,F18.8,6(F12.8))')i,j,- f(j),(x(j,k),k = 1, ndim)
     write(11,'(I8,I8,F18.8,6(F12.8))')i,j,- pbest(j),(xpbest(j,k),k = 1, ndim)
  end do

  write(*,*)"---------------------gbest-----------------------"
  write(12,'(I8,F18.8,6(F12.8))')i,-gbest(i),(xgbest(i,k),k= 1, ndim)
  write(*,'(I8,F18.8,6(F12.8))')i,-gbest(i),(xgbest(i,k),k= 1, ndim)
  write(*,*)"------------------------------------------------"
 
end subroutine write_output 



logical function converge(istep)
  use defs
  implicit none
  integer, intent(in) :: istep 
  logical :: cond, fcond 
  integer :: i   
  real :: diff 
  
  if (istep .lt. nstops) then 
     converge = .false.
  else
     fcond  = .true.
     do i = istep-nstops, istep
        diff = abs(gbest(i)-gbest(i-1))
        if (diff .lt. eps) then 
           cond  = .true.
        else
           cond =  .false.   
        end if
        fcond = fcond.and.cond 
     end do
     converge = fcond
  end if
  
end function converge
