subroutine compute_pr_surv(type_pr,sample_k)
    use global_var;use nrtype
    implicit none
    double precision,dimension(indv,types),intent(in)::type_pr
    integer,dimension(indv,generations),intent(in)::sample_k
    integer,parameter::iterations=2728
    real(dp),dimension(covariates,clusters,L_gender,L_educ,iterations)::c_tr_all
    real(dp),dimension(covariates,clusters,L_gender,L_educ,iterations)::c_tr_d_all
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_h
    real(DP),dimension(covariates,clusters,L_gender,L_educ)::beta_d
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts)::joint_yh
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE 
    real(DP),dimension(clusters,generations,types,L_gender,L_educ,4)::LIV10_model
    real(DP),dimension(2,generations,indv_HRS)::sur_pr_data
    real(DP),dimension(generations,indv_HRS)::sur_pr_model
    integer,dimension(indv,1)::y
    
    open(unit=9,file=path_s//'c_tr.txt')
                read(9,'(F20.8)') c_tr_all
        close(9)
    open(unit=9,file=path_s//'c_tr_d.txt')
            read(9,'(F20.8)') c_tr_d_all
    close(9)
    open(unit=10,file=path//"Data\subj_exp.csv")
        read(10,*) sur_pr_data 
    close(10)
    beta_h=sum(c_tr_all,5)/dble(iterations)
    beta_d=sum(c_tr_d_all,5)/dble(iterations)
    
    joint_yh=1.0d0
    call transitions(beta_h,beta_d,H,LE,joint_yh)
    call simulate_tr(H,LIV10_model)
    
    
    !Sample health behavior type
    call sample_health_behavior(type_pr,y) 
    
    !Compute model predicted survival probabilities
    call model_survival(LIV10_model,y,sur_pr_data,sample_k,sur_pr_model)
    
    
    
end subroutine
    
subroutine simulate_tr(H,LIV10)
use global_var;use nrtype
implicit none
real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ),intent(in)::H
real(DP),dimension(clusters,generations,types,L_gender,L_educ,4),intent(out)::LIV10
integer,parameter::indv_sim=100000
integer::i_l,h_l,e_l,ge_l,y_l,g_l,h_l2,ind,h_ini,t_l
real(DP)::u_h

    LIV10=-9.0d0

    ge_l=1 !gender
    do e_l=1,L_educ !education level
    do y_l=1,2 !lifestyle
    do g_l=12,29 !initial age
    do h_ini=1,clusters !initial health
        LIV10(h_ini,g_l,y_l,ge_l,e_l,:)=0.0d0
        do i_l=1,indv_sim
            do t_l=g_l,generations-1
                if (t_l==g_l) then
                    h_l=1 
                end if
                h_l2=-9
                ind=1
                call RANDOM_NUMBER(u_h)
                do while (h_l2==-9)
                    if (u_h<sum(H(h_l,1:ind,t_l,y_l,ge_l,e_l))) then
                        h_l2=ind
                    else
                        ind=ind+1
                    end if        
                end do
                h_l=h_l2
                if (h_l2==clusters+1 .or. t_l==generations-1) then
                    if (t_l>=27) then
                        LIV10(h_ini,g_l,y_l,ge_l,e_l,1)=LIV10(h_ini,g_l,y_l,ge_l,e_l,1)+1.0d0
                    end if
                    if (t_l>29) then
                        LIV10(h_ini,g_l,y_l,ge_l,e_l,2)=LIV10(h_ini,g_l,y_l,ge_l,e_l,2)+1.0d0
                    end if
                    if (t_l>32) then
                        LIV10(h_ini,g_l,y_l,ge_l,e_l,3)=LIV10(h_ini,g_l,y_l,ge_l,e_l,3)+1.0d0
                    end if
                    if (t_l>34) then
                        LIV10(h_ini,g_l,y_l,ge_l,e_l,4)=LIV10(h_ini,g_l,y_l,ge_l,e_l,4)+1.0d0
                    end if 
                    exit
                end if
            end do    
        end do
        LIV10(h_ini,g_l,y_l,ge_l,e_l,:)=LIV10(h_ini,g_l,y_l,ge_l,e_l,:)/dble(indv_sim)
    end do
    end do
    end do
    end do


end subroutine

    
subroutine model_survival(LIV10_model,y,sur_pr_data,sample_k,sur_pr_model)
use global_var;use nrtype
implicit none
real(DP),dimension(clusters,generations,types,L_gender,L_educ,4),intent(in)::LIV10_model
integer,dimension(indv,1),intent(in)::y
real(DP),dimension(2,generations,indv_HRS),intent(in)::sur_pr_data
integer,dimension(indv,generations),intent(in)::sample_k
real(DP),dimension(generations,indv_HRS),intent(out)::sur_pr_model
integer::g_l,i_l

sur_pr_model=-9.0d0

open(unit=9,file=path_s//'sur_pr_model_vs_data.txt')

        
sur_pr_model=-9.0d0
do i_l=1,indv_HRS;do g_l=first_age(i_l),last_age(i_l)
    if (gender(i_l)==1 .and. sur_pr_data(1,g_l,i_l)/=-9.0d0 .and. sample_k(i_l,g_l)>=1) then
        if (sur_pr_data(2,g_l,i_l)==80.0d0 .and. LIV10_model(sample_k(i_l,g_l),g_l,y(i_l,1),gender(i_l),educ(i_l),1)/=-9.0d0) then
             sur_pr_model(g_l,i_l)=LIV10_model(sample_k(i_l,g_l),g_l,y(i_l,1),gender(i_l),educ(i_l),1)*100.0d0
        elseif (sur_pr_data(2,g_l,i_l)==85.0d0 .and. LIV10_model(sample_k(i_l,g_l),g_l,y(i_l,1),gender(i_l),educ(i_l),2)/=-9.0d0) then
            sur_pr_model(g_l,i_l)=LIV10_model(sample_k(i_l,g_l),g_l,y(i_l,1),gender(i_l),educ(i_l),2)*100.0d0
        elseif (sur_pr_data(2,g_l,i_l)==90.0d0 .and. LIV10_model(sample_k(i_l,g_l),g_l,y(i_l,1),gender(i_l),educ(i_l),3)/=-9.0d0) then
            sur_pr_model(g_l,i_l)=LIV10_model(sample_k(i_l,g_l),g_l,y(i_l,1),gender(i_l),educ(i_l),3)*100.0d0
        elseif (sur_pr_data(2,g_l,i_l)==95.0d0 .and. LIV10_model(sample_k(i_l,g_l),g_l,y(i_l,1),gender(i_l),educ(i_l),4)/=-9.0d0) then
            sur_pr_model(g_l,i_l)=LIV10_model(sample_k(i_l,g_l),g_l,y(i_l,1),gender(i_l),educ(i_l),4)*100.0d0
        end if
        
        if (sur_pr_model(g_l,i_l)/=-9.0d0 ) then
            write(9,'(<6>I5,<2>F8.3,I5)'), i_l,g_l,birth_cohort(i_l),y(i_l,1),educ(i_l),sample_k(i_l,g_l),sur_pr_model(g_l,i_l),sur_pr_data(1,g_l,i_l),int(sur_pr_data(2,g_l,i_l))
        end if
        
    end if
end do;end do

close(9)

end subroutine
    
