subroutine fraction_h_e_g(sample_k,share_h) 
    use global_var; use nrtype
    implicit none
    integer,dimension(indv,generations),intent(in)::sample_k
    real(DP),dimension(clusters,L_gender,L_educ),intent(out)::share_h
    integer::h_l,g_l,ge_l,age,ge_d,it,i_l,d_l,h_l2,t_l,e_l,y_l,c_l
    integer::ind

    share_h=0.0d0
    do i_l=1,indv;
        g_l=first_age(i_l)
        if (sample_k(i_l,g_l)>=1 .and. race(i_l)==1 .and. first_age(i_l)<5) then
            share_h(sample_k(i_l,g_l),gender(i_l),educ(i_l))=share_h(sample_k(i_l,g_l),gender(i_l),educ(i_l))+1.0d0
        end if
    end do
    
    do e_l=1,L_educ;do ge_l=1,L_gender
        share_h(:,ge_l,e_l)=share_h(:,ge_l,e_l)/sum(share_h(:,ge_l,e_l))
    end do;end do

end subroutine