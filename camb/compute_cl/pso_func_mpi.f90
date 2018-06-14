subroutine pso_func_mpi() 
  use defs
  use pso 
  implicit none
  include 'mpif.h'
  integer :: numtasks, rank,ierr,itag,error,I,J  
  double precision,dimension(npart,ndim) :: xtmp 
  double precision,dimension(ndim) :: p
  real(8) :: y
  interface 
     subroutine compute_likelihood_new(p,y)
       Use wmap_likelihood_9yr
       Use wmap_options
       Use wmap_util
       implicit none
       double precision, intent(in) :: p(:)
       double precision, intent(out) :: y
      end subroutine compute_likelihood_new
  end interface

  
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
  
  if (rank.eq.0) then 
     do i = 1, npart
        do j = 1,  ndim 
           xtmp(i,j)= x(i,j)
        end do
     end do
  end if
  
  ! this is a serious bug fix ! in older versions 'ierr' may be missing which 
  ! makes the code difficult to export - jayanti, july 14, 2011
  
  call   MPI_Bcast(xtmp,npart *ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
  do i = rank+1, npart,numtasks
     
     do j= 1, ndim
        p(j) =  xtmp(i,j)
     end do
     
     call compute_likelihood_new(p,y)  
     
     write(*,'(7(F14.6))')y,(p(j),j=1, ndim)  
     
     call MPI_GATHER (-y,1,MPI_DOUBLE_PRECISION,f(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     
  end do
  
end subroutine pso_func_mpi

