subroutine transitions(beta_h,beta_d,init_cond,H,LE)
    use global_var; use nrtype
    implicit none
    real(DP),dimension(covariates,types,clusters,clusters),intent(in)::beta_h
    real(DP),dimension(covariates,types,clusters),intent(in)::beta_d
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(out)::H
    real(DP),dimension(types,L_gender,L_educ,clusters+1),intent(out)::LE
    real(DP),dimension(clusters,types,L_gender,L_educ),intent(in)::init_cond
    integer::e_l,c_l,c_l2,g_l,ge_l,age,it,max_loc,d_l,c_l3,t_l
    double precision,dimension(L_educ-1)::educ_d
    real(DP),dimension(covariates,1)::x
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
    
    !!$OMP PARALLEL  DEFAULT(PRIVATE) SHARED(H,beta)
    !!$OMP  DO collapse(4)
    do t_l=1,types; do c_l=1,clusters; do g_l=1,generations; do ge_l=1,2;do e_l=1,L_educ
        age=initial_age+(g_l-1)*2-70
        gender_d=ge_l-1
        x(1:4,1)=(/1.0_dp,dble(age),gender_d,dble(age)*gender_d/)
        educ_d=0.0d0
        if (e_l>1) then
            educ_d(e_l-1)=1.0d0
        end if        
        x(5:6,1)=educ_d
        x(7:8,1)=educ_d*dble(age)

        if (clusters==2) then
            H(c_l,1,g_l,t_l,ge_l,e_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_h(:,t_l,c_l,1))/sqrt(2.0_dp)))
            H(c_l,2,g_l,t_l,ge_l,e_l)=0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_h(:,t_l,c_l,1))/sqrt(2.0_dp)))
        else
            !By numerical integration
            do c_l2=1,clusters
                h_star(c_l2)=sum(x(:,1)*beta_h(:,t_l,c_l,c_l2))
            end do
            do c_l2=1,clusters
                prod1=1
                prod2=1
                do c_l3=1,clusters
                    if (c_l2/=c_l3) then
                        prod1=prod1*0.5d0*(1.0d0+erf((-sqrt(2.0d0*xs)-(h_star(c_l3)-h_star(c_l2)))/sqrt(2.0d0)))
                        prod2=prod2*0.5d0*(1.0d0+erf(( sqrt(2.0d0*xs)-(h_star(c_l3)-h_star(c_l2)))/sqrt(2.0d0)))
                    end if
                end do
                H(c_l,c_l2,g_l,t_l,ge_l,e_l)=0.5d0/sqrt(pi)*sum(weight*(prod1+prod2))
            end do
            
        end if
        H(c_l,clusters+1,g_l,t_l,ge_l,e_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*beta_d(:,t_l,c_l))/sqrt(2.0_dp)))
        H(c_l,1:clusters,g_l,t_l,ge_l,e_l)=H(c_l,1:2,g_l,t_l,ge_l,e_l)/sum(H(c_l,1:clusters,g_l,t_l,ge_l,e_l))*(1.0d0-H(c_l,clusters+1,g_l,t_l,ge_l,e_l))
    end do; end do; end do; end do;end do
    !!$OMP END DO
    !!$OMP END PARALLEL

    !No resurection
    H(clusters+1,1:clusters,:,:,:,:)=0.0_dp 
    H(clusters+1,clusters+1,:,:,:,:)=1.0_dp

    !Compute LE in each health status
    LE=0.0_dp
    do ge_l=1,2; do t_l=1,types;do e_l=1,L_educ
        p(1:clusters,1)=init_cond(:,t_l,ge_l,e_l)
        do g_l=2,generations
            do c_l=1,clusters;do c_l2=1,clusters
                if (c_l==c_l2) then
                    LE(t_l,ge_l,e_l,c_l)=LE(t_l,ge_l,e_l,c_l)+2.0d0*H(c_l,c_l,g_l-1,t_l,ge_l,e_l)*p(c_l,g_l-1)
                else
                    LE(t_l,ge_l,e_l,c_l)=LE(t_l,ge_l,e_l,c_l)+1.0d0*H(c_l,c_l2,g_l-1,t_l,ge_l,e_l)*p(c_l,g_l-1)
                    LE(t_l,ge_l,e_l,c_l)=LE(t_l,ge_l,e_l,c_l)+1.0d0*H(c_l,c_l2,g_l-1,t_l,ge_l,e_l)*p(c_l,g_l-1)
                end if
                p(:,g_l)=matmul(transpose(H(:,:,g_l-1,t_l,ge_l,e_l)),p(:,g_l-1))    
            end do;end do
        end do
        LE(t_l,ge_l,e_l,clusters+1)=sum(LE(t_l,ge_l,e_l,1:clusters))
    end do;end do;end do

end subroutine