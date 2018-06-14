!--------------------------------------------------------------------------------------
!                                COSMOPSO V2.0
!-------------------------------------------------------------------------------------
!  *  This is the Particle Swarm Optimization function written for the
!     COSMOPSO package by Jayanti Prasad.  
!  *   Since this follows the most basic algorithm so I am not giving references 
!     but it is worth to check James Kennedy and Russell Eberhart's 1995 paper:
!    << Kennedy, J. and Eberhart, R. C. Particle swarm optimization. 
!     Proceedings of IEEE International Conference on Neural Networks, Piscataway, 
!     NJ. pp. 1942-1948, 1995>>
!  * Details of our implimentaion are given in 
!   Jayanti Prasad, Tarun Souradeep (2012), Phys. Rev. D 85, 123008 (2012) [arXiv:1108.5600v2 ]
!   Cosmological parameter estimation using Particle Swarm Optimization (PSO) 
!                                                   -----Jayanti Prasad, Sept 20, 2014
!------------------------------------------------------------------------------------

module pso   
  use defs
contains 
  
  subroutine pso_main(itime) 
    use defs 
    implicit none 
    include 'mpif.h'
    integer, intent(in) :: itime   
    real :: xi_1,xi_2 
    integer ::rank,ierr,numtasks,itag,i,j  
    double precision :: r8_uniform, a_uran,b_uran 
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
    
    if (rank.eq.0) then 
       
       gbest(itime) = gbest(itime-1)
       xgbest(itime,:) = xgbest(itime-1,:)
   
        a_uran = 0.d0  
        b_uran = 1.d0 

        do i = 1, npart 
          do j =  1, ndim    
   
           xi_1  = r8_uniform (a_uran,b_uran,seed_uran )
           xi_2  = r8_uniform (a_uran,b_uran,seed_uran )
           
           v(i,j) = w * v(i,j) + c_pbest * xi_1 * (xpbest(i,j)-x(i,j)) + c_gbest * &
                xi_2 * (xgbest(itime,j)-x(i,j))  
      
          end do 

        end do
        
        call maxvel()
        
        x(:,:) = x(:,:) + v(:,:)  
        
       call  boundary()
     
 
    end if
    
    f = 0.d0 
 
    call pso_func_mpi()
    
    if (rank.eq.0) then 
        
       call compute_best(itime)
 
    end if
    
  end subroutine pso_main
  
  subroutine pso_init()  
    ! initilize the postion and the velocities of particles
    ! position and velocities are initilized randomly  
    implicit none 
    include 'mpif.h'
    integer :: i,j,itime,rank,numtasks,ierr
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
    
    if (rank.eq.0) then     
       call  pso_initpos_random()
    end if
    
    call pso_func_mpi() ! This will return f(i),i=1,npart 
    
    if (rank.eq.0) then 
       
       pbest(:) = f(:) 
     
       xpbest(:,:) = x(:,:) 
  
       itime  = 1 
 
       gbest(itime) = -1.0E38       

       do i =  1,  npart   
         if (pbest(i) .gt. gbest(1)) then 
           gbest(itime)  = pbest(i)
           xgbest(itime,:) = xpbest(i,:)        
         end if    
       end do 
      
 
       call compute_best(itime) 
   
    end if
    
  end subroutine pso_init
  
  subroutine  boundary() 
    implicit none 
    integer :: i,j  
    
    do i  =  1, npart  
       do j =  1, ndim 
          if (x(i,j) .gt. xmax(j)) then 
             x(i,j) =  xmax(j)
             v(i,j) = -v(i,j)    
          end if
          
          if (x(i,j) .lt. xmin(j)) then
             x(i,j) =  xmin(j)
             v(i,j) = -v(i,j)
          end if
       end do
    end do
    
  end subroutine  boundary
  
  subroutine maxvel()
    implicit none 
    integer :: i, j
    
    do i =  1, npart    
       
       do j =  1, ndim  
          if (v(i,j) .lt. vmin(j)) then 
             v(i,j) = vmin(j) 
          end if
          
          if (v(i,j) .gt. vmax(j)) then
             v(i,j) = vmax(j)
          end if
       end do
       
    end do
    
  end subroutine maxvel
  
  
  ! initially particles positions can be assigned randomly 
  
  subroutine pso_initpos_random()
    use defs 
    implicit none 
    integer :: i,j 
    double precision :: r8_uniform, a_uran,b_uran 
    
    ! set the maximum velicity here 
    
    vmax(:) = abs(xmax(:)-xmin(:)) * c_vmax  

    vmin(:) = -vmax(:)
 
    a_uran = 0.d0
    b_uran = 1.d0  

    do i  =  1, npart 
       do j =  1, ndim 
       x(i,j) =  xmin(j) + (xmax(j)-xmin(j)) * &
        r8_uniform (a_uran,b_uran, seed_uran)
       
       v(i,j) = c_vmax  * (r8_uniform (a_uran,b_uran, seed_uran) - 0.5d0) 
       end do   

    end do
  
    call boundary()
    
  end subroutine pso_initpos_random
  
  subroutine compute_best(itime)
    implicit none 
    integer, intent(in) :: itime 
    integer :: i
   
    if (itime .gt. 1) then
       gbest(itime) = gbest(itime-1)
       xgbest(itime,:) = xgbest(itime-1,:)
     end if 
 
    do i = 1,npart 
       
       if (f(i).gt.pbest(i)) then ! update pbest 
          
          pbest(i)  = f(i) 
          xpbest(i,:) = x(i,:) 
       end if
       
    end do 

    do i =  1, npart  
       if(pbest(i) .gt. gbest(itime))  then  ! update gbest  
          gbest(itime) = pbest(i)
          xgbest(itime,:) = xpbest(i,:) 
         end if
    end do
    
  end subroutine compute_best
  
end module pso

