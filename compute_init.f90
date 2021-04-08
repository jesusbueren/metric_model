subroutine compute_init(sample_k,y,init_cond)  
    use global_var; use nrtype
    implicit none
    integer,dimension(indv,generations),intent(in)::sample_k
    integer,dimension(indv,1),intent(in)::y
    real(DP),dimension(clusters,types,L_gender,L_educ)::init_cond 
    integer,dimension(types,L_gender,L_educ)::counter
    integer::c_l,i_l,g_l
    
    counter=0.0d0
    init_cond=0.0d0
    do i_l=1,indv;do g_l=1,3
        if (sample_k(i_l,g_l)>=1 .and. sample_k(i_l,g_l)<=2 ) then
            counter(y(i_l,1),gender(i_l),educ(i_l))=counter(y(i_l,1),gender(i_l),educ(i_l))+1
            init_cond(sample_k(i_l,g_l),y(i_l,1),gender(i_l),educ(i_l))=init_cond(sample_k(i_l,g_l),y(i_l,1),gender(i_l),educ(i_l))+1.0d0
        end if
    end do;end do
    
    do c_l=1,clusters
        init_cond(c_l,:,:,:)=init_cond(c_l,:,:,:)/dble(counter)
    end do
    
    end subroutine