subroutine likelihood_all(p,likelihood)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(indv,clusters,generations),intent(out)::likelihood
    integer::g_l,i_l,v_l,c_l
    real(DP),dimension(adls,clusters),intent(in)::p
    
    likelihood(:,:,:)=1.0d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(g_l,i_l,v_l,c_l)
    !$OMP DO
    do i_l=1,indv;do g_l=first_age(i_l),last_age(i_l)    
        if (data_adls(i_l,1,g_l)/=-1) then
          do v_l=1,adls; do c_l=1,clusters
            if (data_adls(i_l,v_l,g_l)==1) then
                likelihood(i_l,c_l,g_l)=likelihood(i_l,c_l,g_l)*p(v_l,c_l)
            elseif (data_adls(i_l,v_l,g_l)==0) then
                likelihood(i_l,c_l,g_l)=likelihood(i_l,c_l,g_l)*(1-p(v_l,c_l))
            end if 
          end do; end do      
        end if
    end do; end do
    !$OMP END DO
    !$OMP END PARALLEL
end subroutine
    
subroutine likelihood_i(i_l,c_l,g_l,p,likelihood)
    use global_var; use nrtype
    implicit none
    real(DP),intent(out)::likelihood
    integer,intent(in)::i_l,c_l,g_l
    real(DP),dimension(adls,clusters),intent(in)::p
    integer::v_l
    
    likelihood=1.0d0
    do v_l=1,adls
        if (data_adls(i_l,v_l,g_l)==1) then
            likelihood=likelihood*p(v_l,c_l)
        elseif (data_adls(i_l,v_l,g_l)==0) then
            likelihood=likelihood*(1-p(v_l,c_l))
        end if
    end do
end subroutine