subroutine compute_weights(delta,H,share_h,weights,joint_yh) 
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types),intent(in)::delta
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
    real(DP),dimension(clusters,L_gender,L_educ),intent(in)::share_h
    real(DP),dimension(generations,clusters,L_gender,L_educ,types),intent(out)::weights,joint_yh
    real(DP),dimension(types)::y_star
    real(DP),dimension(clusters,types)::joint_yh_new
    real(DP)::age,health_d
    real(DP),dimension(covariates_mixture,1)::x
    integer::e_l,ge_l,g_l,h_l,y_l,y_l2,i_l,h_l2
    integer,parameter::nodes=5
    real(DP),dimension(nodes):: xs=(/0.117581320211778,	1.0745620124369,	3.08593744371755,	6.41472973366203,	11.8071894899717/), &
                                node_weight=(/1.22172526747065,	0.480277222164629,	0.0677487889109621,	0.00268729149356246,	1.52808657104652E-05/),prod1,prod2
    
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    weights=0.0d0
    joint_yh=0.0d0
    do e_l=1,L_educ;do ge_l=1,L_gender
        g_l=1
        do h_l=1,clusters
            health_d=dble(h_l-1)
            x(1:covariates_mixture,1)=(/1.0_dp,health_d/)
            do y_l=1,clusters
                y_star(y_l)=sum(x(:,1)*delta(:,ge_l,e_l,y_l))
            end do
            do y_l=1,types
                prod1=1
                prod2=1
                do y_l2=1,types
                    if (y_l2/=y_l) then
                        prod1=prod1*0.5d0*(1.0d0+erf((-sqrt(2.0d0*xs)-(y_star(y_l2)-y_star(y_l)))/sqrt(2.0d0)))
                        prod2=prod2*0.5d0*(1.0d0+erf(( sqrt(2.0d0*xs)-(y_star(y_l2)-y_star(y_l)))/sqrt(2.0d0)))
                    end if
                end do
                weights(g_l,h_l,ge_l,e_l,y_l)=0.5d0/sqrt(pi)*sum(node_weight*(prod1+prod2))  !weights(:,1,1,1,1)
                joint_yh(g_l,h_l,ge_l,e_l,y_l)=joint_yh(g_l,h_l,ge_l,e_l,y_l)+weights(g_l,h_l,ge_l,e_l,y_l)*share_h(h_l,ge_l,e_l)
            end do 
        end do
        do g_l=2,generations
            joint_yh_new=0.0d0
            do h_l=1,clusters;do y_l=1,types
                do h_l2=1,clusters
                    joint_yh_new(h_l,y_l)=joint_yh_new(h_l,y_l)+H(h_l2,h_l,g_l-1,y_l,ge_l,e_l)*joint_yh(g_l-1,h_l2,ge_l,e_l,y_l)
                end do
            end do; end do
            joint_yh(g_l,:,ge_l,e_l,:)=joint_yh_new/sum(joint_yh_new)
            do h_l=1,clusters;do y_l=1,types
                weights(g_l,h_l,ge_l,e_l,y_l)=joint_yh(g_l,h_l,ge_l,e_l,y_l)/sum(joint_yh(g_l,h_l,ge_l,e_l,:)) 
                if (isnan(weights(g_l,h_l,ge_l,e_l,y_l))) then
                    print*,''
                end if
            end do; end do
        end do
    end do;end do
    


end subroutine
    
