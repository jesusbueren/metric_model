subroutine transitions(beta_h,beta_d,H,LE,joint_yh)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_h
    real(DP),dimension(covariates,clusters,L_gender,L_educ),intent(in)::beta_d
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(in)::joint_yh
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(out)::H
    real(DP),dimension(types,L_gender,L_educ,clusters+1),intent(out)::LE
    integer::e_l,c_l,c_l2,g_l,ge_l,age,it,max_loc,d_l,c_l3,t_l
    double precision,dimension(L_educ-1)::educ_d
    real(DP),dimension(covariates,1)::x
    real(DP),dimension(types-1)::dummy_type,dummy_type_x_age
    real(DP),dimension(clusters,clusters,generations,types,L_gender,L_educ)::H_new
    integer,parameter::sims=1000
    integer,dimension(clusters+1)::counter_h
    real(DP),dimension(clusters)::h_star
    double precision,dimension(clusters+1,generations)::p
    integer,parameter::nodes=5
    real(DP),dimension(nodes):: xs=(/0.117581320211778,	1.0745620124369,	3.08593744371755,	6.41472973366203,	11.8071894899717/), &
                                weight=(/1.22172526747065,	0.480277222164629,	0.0677487889109621,	0.00268729149356246,	1.52808657104652E-05/),prod1,prod2
    integer::ind
    real(DP)::gender_d
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    
    H=-9.0d0
    !!$OMP PARALLEL  DEFAULT(PRIVATE) SHARED(H,beta)
    !!$OMP  DO collapse(4)
    do t_l=1,types; do c_l=1,clusters; do g_l=generations,1,-1; do ge_l=1,L_gender;do e_l=1,L_educ
        age=initial_age+(g_l-1)*2-70
        dummy_type=0.0d0
        dummy_type_x_age=0.0d0
        if (t_l>1)then
            dummy_type(t_l-1)=1.0d0
            dummy_type_x_age(t_l-1)=dble(age)
        end if
        x(:,1)=[(/1.0_dp,dble(age)/),dummy_type,dummy_type_x_age]!,dble(age)**2.0d0
        if (clusters==2) then
            H(c_l,1,g_l,t_l,ge_l,e_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_h(:,c_l,ge_l,e_l))/sqrt(2.0_dp)))
            H(c_l,2,g_l,t_l,ge_l,e_l)=1.0d0-H(c_l,1,g_l,t_l,ge_l,e_l)
        else
            print*,'need to do smthg else'
            !!By numerical integration
            !do c_l2=1,clusters
            !    h_star(c_l2)=sum(x(:,1)*beta_h(:,t_l,c_l,c_l2,ge_l,e_l))
            !end do
            !do c_l2=1,clusters
            !    prod1=1
            !    prod2=1
            !    do c_l3=1,clusters
            !        if (c_l2/=c_l3) then
            !            prod1=prod1*0.5d0*(1.0d0+erf((-sqrt(2.0d0*xs)-(h_star(c_l3)-h_star(c_l2)))/sqrt(2.0d0)))
            !            prod2=prod2*0.5d0*(1.0d0+erf(( sqrt(2.0d0*xs)-(h_star(c_l3)-h_star(c_l2)))/sqrt(2.0d0)))
            !        end if
            !    end do
            !    H(c_l,c_l2,g_l,t_l,ge_l,e_l)=0.5d0/sqrt(pi)*sum(weight*(prod1+prod2))
            !end do   
        end if
            
        H(c_l,clusters+1,g_l,t_l,ge_l,e_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_d(:,c_l,ge_l,e_l))/sqrt(2.0_dp)))



        if (isnan(sum(H(c_l,:,g_l,t_l,ge_l,e_l)))) then
            print*,'error in transitions'
        end if
        H(c_l,1:clusters,g_l,t_l,ge_l,e_l)=H(c_l,1:2,g_l,t_l,ge_l,e_l)*(1.0d0-H(c_l,clusters+1,g_l,t_l,ge_l,e_l))
    end do; end do; end do; end do;end do
    !!$OMP END DO
    !!$OMP END PARALLEL !H(1,3,:,3,1,3)

    !No resurection
    H(clusters+1,1:clusters,:,:,:,:)=0.0_dp 
    H(clusters+1,clusters+1,:,:,:,:)=1.0_dp

    !Compute LE in each health status at age 50
    LE=0.0d0

    do ge_l=1,L_gender; do t_l=1,types;do e_l=1,L_educ
        p=-9.0d0
        if (cohorts==5) then
            p(1:clusters,1)=joint_yh(1,:,ge_l,e_l,t_l,3)/sum(joint_yh(1,:,ge_l,e_l,t_l,3))
        else
            p(1:clusters,1)=joint_yh(1,:,ge_l,e_l,t_l,5)/sum(joint_yh(1,:,ge_l,e_l,t_l,5))
        end if
        if (isnan(sum(p)))then
            print*,'error in transitions: initial cond. Don t worry if it=1'
        end if
        do g_l=2,generations-4
            if (g_l>12) then
                if (g_l==13) then
                    p(3,g_l-1)=0.0d0
                    p(:,g_l-1)=p(:,g_l-1)/sum(p(:,g_l-1)) 
                end if                    
                do c_l=1,clusters;do c_l2=1,clusters
                    if (c_l==c_l2) then
                        LE(t_l,ge_l,e_l,c_l)=LE(t_l,ge_l,e_l,c_l)+2.0d0*H(c_l,c_l,g_l-1,t_l,ge_l,e_l)*p(c_l,g_l-1)
                    else
                        LE(t_l,ge_l,e_l,c_l)=LE(t_l,ge_l,e_l,c_l)+1.0d0*H(c_l,c_l2,g_l-1,t_l,ge_l,e_l)*p(c_l,g_l-1)
                        LE(t_l,ge_l,e_l,c_l)=LE(t_l,ge_l,e_l,c_l)+1.0d0*H(c_l2,c_l,g_l-1,t_l,ge_l,e_l)*p(c_l2,g_l-1)
                    end if
                end do;end do
            end if
            p(:,g_l)=matmul(transpose(H(:,:,g_l-1,t_l,ge_l,e_l)),p(:,g_l-1))    
        end do
        LE(t_l,ge_l,e_l,clusters+1)=sum(LE(t_l,ge_l,e_l,1:clusters))
    end do;end do;end do

end subroutine