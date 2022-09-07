subroutine compute_weights(fraction_y,H,share_h,weights,joint_yh) 
    use global_var; use nrtype
    implicit none
    real(DP),dimension(clusters,L_gender,L_educ,types,cohorts),intent(in)::fraction_y
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    real(DP),dimension(clusters,L_gender,L_educ),intent(in)::share_h
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(out)::weights,joint_yh
    real(DP),dimension(clusters,types)::joint_yh_new
    real(DP),dimension(covariates_mixture,1)::x
    integer::e_l,ge_l,g_l,h_l,y_l,y_l2,i_l,h_l2,co_l
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    weights=0.0d0
    joint_yh=0.0d0
    do e_l=1,L_educ;do ge_l=1,L_gender;do co_l=1,cohorts
        g_l=1
        do h_l=1,clusters
            do y_l=1,types
                weights(g_l,h_l,ge_l,e_l,y_l,co_l)=fraction_y(h_l,ge_l,e_l,y_l,co_l)
                joint_yh(g_l,h_l,ge_l,e_l,y_l,co_l)=weights(g_l,h_l,ge_l,e_l,y_l,co_l)*share_h(h_l,ge_l,e_l)
            end do 
        end do
        do g_l=2,generations
            joint_yh_new=0.0d0
            do y_l=1,types;do h_l=1,clusters
                do h_l2=1,clusters
                    joint_yh_new(h_l,y_l)=joint_yh_new(h_l,y_l)+H(h_l2,h_l,g_l-1,y_l,ge_l,e_l)*joint_yh(g_l-1,h_l2,ge_l,e_l,y_l,co_l) 
                end do
            end do; end do
            joint_yh(g_l,:,ge_l,e_l,:,co_l)=joint_yh_new/sum(joint_yh_new)
            do h_l=1,clusters;do y_l=1,types
                weights(g_l,h_l,ge_l,e_l,y_l,co_l)=joint_yh(g_l,h_l,ge_l,e_l,y_l,co_l)/sum(joint_yh(g_l,h_l,ge_l,e_l,:,co_l)) 
                if (isnan(weights(g_l,h_l,ge_l,e_l,y_l,co_l))) then
                    print*,'problem computing weights'
                end if
            end do; end do
        end do
    end do;end do;end do
    


end subroutine
    
