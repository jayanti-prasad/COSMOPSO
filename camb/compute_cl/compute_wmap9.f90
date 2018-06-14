   function test_wmap9(cl_tt, cl_te, cl_ee, cl_bb)
    Use wmap_likelihood_9yr
    Use wmap_options
    Use wmap_util
    implicit none 
    double precision, intent(in) :: cl_tt(:), cl_te(:), cl_ee(:), cl_bb(:) 
    double precision :: test_wmap9
    Real(8) :: like(num_WMAP), like_tot, expected_like_tot
    Integer :: lun, l, i, olun
    Integer :: tt_npix, teeebb_npix

    use_TT = .True.
    use_TE = .True.
    use_lowl_TT = .True.
    use_lowl_pol = .True.

    !---------------------------------------------------
    ! get likelihoods
    !---------------------------------------------------
    like = 0.d0
    !Call wmap_likelihood_init
    Call wmap_likelihood_dof(tt_npix, teeebb_npix)
    Call wmap_likelihood_compute(cl_tt, cl_te, cl_ee, cl_bb, like)
    Call wmap_likelihood_error_report

    like_tot = Sum(like(1:num_WMAP))

    test_wmap9=like_tot 

    !write(*,*)"wmap likelihood=",like_tot

end function test_wmap9 
