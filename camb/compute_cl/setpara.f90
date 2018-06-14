module set_para
  use defs
  use IniFile 
 
contains 
  
  subroutine read_para(paramfile)
    implicit none
    character (len=MAX_LEN),intent(in)  :: paramfile  
    character(len=MAX_LEN) :: prior
    integer :: i  
    logical :: bad 

    call Ini_Open(trim(paramfile), 1, bad, .false.)
    if (bad) stop 'Error opening parameter file'

     seed_uran= Ini_Read_Int('seed_uran',seed_uran)
     ndim     = Ini_Read_Int('ndim',ndim)
     npart    = Ini_Read_Int('npart',npart)
     nstops   = Ini_Read_Int('nstops',nstops)
     eps      = Ini_Read_Double('eps',eps)   
     c_gbest  = Ini_Read_Double('c_gbest',c_gbest);
     c_pbest  = Ini_Read_Double('c_pbest',c_pbest);
     c_cbest  = Ini_Read_Double('c_cbest',c_cbest);
     w        = Ini_Read_Double('w',w);
     c_vmax   = Ini_Read_Double('c_vmax',c_vmax);
     planck_data   = Ini_Read_String('planck_data')
     prior    = Ini_Read_String("prior_file")

     USE_WMAP   = Ini_Read_Logical('USE_WMAP')
     USE_PLANCK = Ini_Read_Logical('USE_PLANCK')
 
     write(*,'(4(A12,I8))')"seed=",seed_uran,"ndim=",ndim,"npart=",npart,"nstops=",nstops
     write(*,'(5(A12,D12.6))')"c1=",c_pbest,"c2=",c_gbest,"w=",w,"c_vmax=",c_vmax,"eps=",eps
     write(*,*)"Planck data=",trim(planck_data)
     write(*,*)"USE WMAP:",USE_WMAP,"USE_PLANCK=",USE_PLANCK 

    call Ini_Close

    open(11,file=trim(prior))
    do i = 1, ndim
       read(11,*) xdef(i),xmin(i),xmax(i)
       !write(*,*) "param",xdef(i),xmin(i),xmax(i)
    end do
    close(11)
   
  end subroutine read_para
  
end module set_para
