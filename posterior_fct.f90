subroutine posterior_fct_mu(p,sample_k,posterior_mu)
    use nrtype; use global_var
    real(DP),intent(out)::posterior_mu
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(adls,clusters),intent(in)::p
    real(DP),dimension(indv)::posterior_mu_i
    integer::i_l,e_l,h_l,c_l,g_l,ge_l,age,ge_d,health_d
    double precision::likelihood,prior

    posterior_mu_i=1.0d0
    !!$OMP PARALLEL NUM_THREADS(10) DEFAULT(SHARED) PRIVATE(i_l,g_l,likelihood)
    !!$OMP DO 
    do i_l=1,indv; do g_l=first_age(i_l),last_age(i_l)
        if (data_adls(i_l,1,g_l)/=-1) then
            call likelihood_i(i_l,sample_k(i_l,g_l),g_l,p,likelihood)
            posterior_mu_i(i_l)=posterior_mu_i(i_l)*likelihood
        end if
    end do; end do
    !!$OMP END DO
    !!$OMP END PARALLEL
    posterior_mu=sum(log(posterior_mu_i))

    
end subroutine