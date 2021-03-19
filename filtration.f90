subroutine filtration(H,p,smoothed_states,L_i,y,filtered_states)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    real(DP),dimension(adls,clusters),intent(in)::p
    double precision,dimension(indv,clusters+1,generations),intent(in)::smoothed_states
    real(DP),dimension(indv,clusters,generations),intent(in)::L_i
    integer,dimension(indv,1),intent(in)::y
    real(DP),dimension(indv,clusters,generations),intent(out)::filtered_states
    real(DP),dimension(clusters,clusters)::joint_states
    real(DP)::A
    integer::g_l,i_l,c_l1,c_l2
    
    filtered_states=0.0d0
    !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(g_l,i_l,c_l1,c_l2,A,joint_states)
    !$OMP DO 
    do i_l=1,indv
        do g_l=1,last_age(i_l)
            do c_l1=1,clusters; do c_l2=1,clusters
                if (g_l==1) then
                    joint_states(c_l2,c_l1)=smoothed_states(i_l,c_l2,g_l)*L_i(i_l,c_l2,g_l)
                else
                    joint_states(c_l2,c_l1)=filtered_states(i_l,c_l1,g_l-1)*H(c_l1,c_l2,g_l-1,y(i_l,1),gender(i_l),educ(i_l))*L_i(i_l,c_l2,g_l)
                end if
            end do;end do
            A=sum(joint_states)
            do c_l2=1,clusters
                filtered_states(i_l,c_l2,g_l)=sum(joint_states(c_l2,:))/A
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

end subroutine