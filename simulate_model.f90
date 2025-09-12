subroutine simulate_model()
    use global_var; use nrtype
    implicit none
    integer,parameter:: iterations=2999
    real(dp),dimension(covariates,clusters,L_gender,L_educ,iterations)::beta_h_all,beta_d_all
    real(dp),dimension(covariates_habits,habits,types,iterations)::gamma_all
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types,iterations)::delta_all
    real(dp),dimension(covariates,clusters,L_gender,L_educ)::beta_h,beta_d
    real(dp),dimension(covariates_habits,habits,types)::gamma
    real(DP),dimension(covariates_mixture,L_gender,L_educ,types)::delta
    real(DP),dimension(clusters,L_gender,L_educ,types,cohorts)::fraction
    real(DP),dimension(clusters,L_gender,L_educ)::share_h
    real(DP),dimension(clusters+1,clusters+1,generations,types,L_gender,L_educ)::H
    real(DP),dimension(types,L_gender,L_educ,clusters+1)::LE !LE(:,1,3,3)
    real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts)::weights,joint_yh,joint_yh2 !p(y|h) p(y,h)
    integer,dimension(indv,generations)::sample_k
    integer::h_l,c_l,g_l,e_d,age,h2_l,i_l,e_l
    real(DP),dimension(habits,generations,types,clusters)::alphas
    real(DP),dimension(covariates_habits,1)::x
    real(dp)::health_d,u
    integer,dimension(indv,1)::y
    real(DP),dimension(indv,types)::type_pr
    integer,dimension(indv_HRS,habits,generations)::data_habits_sim
    integer,dimension(1)::seed=456
    real(DP),dimension(habits,L_educ)::pr_h_e_sim,pr_h_e_data
    real(DP),dimension(habits,habits,2)::cov_h1_h2_sim,cov_h1_h2_data,obs_h1_h2
    integer,dimension(habits,L_educ)::obs_h_e
    real(DP),dimension(habits,2)::pr_h_a_sim,pr_h_a_data,var_h_a_sim,var_h_a_data,cov_h_a_sim,cov_h_a_data
    integer,dimension(habits,2)::obs_h_a
    real(DP),dimension(habits)::mean_h_sim,mean_h_data
    character(len=4)::table_name
    integer,dimension(4)::burn=(/0,0,1000,2000/)
        call random_seed(PUT=seed)
    
        !Compute share of indv in good and bad health in the initial period across education and gender
        sample_k=data_shlt
        call fraction_h_e_g(sample_k,share_h)
    
        open(unit=9,file=path_s//'c_tr_'//types_s//'.txt')
            read(9,'(F20.8)') beta_h_all
        close(9)
        open(unit=10,file=path_s//'c_tr_d_'//types_s//'.txt')
            read(10,'(F20.8)') beta_d_all
        close(10)
        open(unit=11,file=path_s//'c_habits_'//types_s//'.txt')
            read(11,'(F20.8)') gamma_all
        close(11)
         open(unit=16,file=path_s//'delta_'//types_s//'.txt')
            read(16,'(F20.8)') delta_all
        close(16)
        
        beta_h=sum(beta_h_all(:,:,:,:,burn(types)+1:iterations),5)/dble(iterations-burn(types))
        beta_d=sum(beta_d_all(:,:,:,:,burn(types)+1:iterations),5)/dble(iterations-burn(types))
        ! Given health and survival parameters compute health transitions
        joint_yh=1.0d0/dble(clusters*types) 
        call transitions(beta_h,beta_d,H,LE,joint_yh) 
        
        !Compute the joint pr(y,h|a,e,s,c)
        delta=sum(delta_all(:,:,:,:,burn(types)+1:iterations),5)/dble(iterations-burn(types))
        call delta_2_fraction(delta,fraction)
        call compute_weights(fraction,H,share_h,weights,joint_yh) 
        compute_LE=1
        call transitions(beta_h,beta_d,H,LE,joint_yh) 
        call create_table_3(LE,joint_yh)
        
        !Compute pr of behaviors given type
        gamma=sum(gamma_all(:,:,:,burn(types)+1:iterations),4)/dble(iterations-burn(types))
        do e_l=1,types; do h_l=1,habits;do c_l=1,clusters; do g_l=1,generations
            age=initial_age+(g_l-1)*2-70
            health_d=dble(c_l-1)
            x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),health_d/)
            alphas(h_l,g_l,e_l,c_l)=1.0_dp-0.5_dp*(1.0_dp+erf(-sum(x(:,1)*gamma(:,h_l,e_l))/sqrt(2.0_dp)))
        end do;end do; end do;end do
        
        !Sample a type for each individual in the data
        y=1
        call sample_y(gamma,y,sample_k,H,weights,type_pr)
        
        !Compute artificial set of behaviors for individuals for all indv in the HRS
        data_habits_sim=-9
        do i_l=1,indv_HRS;
        do g_l=first_age(i_l),last_age(i_l);do h_l=1,habits
            if (sample_k(i_l,g_l)/=-1) then
                health_d=dble(sample_k(i_l,g_l)-1)
                age=initial_age+(g_l-1)*2-70
                x(:,1)=(/1.0_dp,dble(age),dble(age**2.0_dp-1.0_dp),health_d/)
                if ((data_habits(i_l,h_l,g_l)==1 .or. data_habits(i_l,h_l,g_l)==0) .and. race(i_l)==1) then
                    call RANDOM_NUMBER(u)
                    if (u<alphas(h_l,g_l,y(i_l,1),sample_k(i_l,g_l))) then
                        data_habits_sim(i_l,h_l,g_l)=1
                    else
                        data_habits_sim(i_l,h_l,g_l)=0
                    end if
                end if
            end if
        end do; end do
        end do
       
        !Mean behavior across education groups and young and old: model vs data
        pr_h_e_sim=0.0d0
        pr_h_e_data=0.0d0
        obs_h_e=0
        pr_h_a_sim=0.0d0
        pr_h_a_data=0.0d0
        obs_h_a=0
        do i_l=1,indv_HRS;
        do g_l=first_age(i_l),last_age(i_l);do h_l=1,habits
            if (sample_k(i_l,g_l)/=-1 .and. (data_habits(i_l,h_l,g_l)==1 .or. data_habits(i_l,h_l,g_l)==0) .and. race(i_l)==1) then
               pr_h_e_sim(h_l,educ(i_l))=pr_h_e_sim(h_l,educ(i_l))+data_habits_sim(i_l,h_l,g_l)
               pr_h_e_data(h_l,educ(i_l))=pr_h_e_data(h_l,educ(i_l))+data_habits(i_l,h_l,g_l)
               obs_h_e(h_l,educ(i_l))=obs_h_e(h_l,educ(i_l))+1
               age=initial_age+(g_l-1)*2
               if (age>=65 .and. age<=69) then
                   pr_h_a_sim(h_l,1)=pr_h_a_sim(h_l,1)+data_habits_sim(i_l,h_l,g_l)
                   pr_h_a_data(h_l,1)=pr_h_a_data(h_l,1)+data_habits(i_l,h_l,g_l)
                   obs_h_a(h_l,1)=obs_h_a(h_l,1)+1
               elseif (age>=75 .and. age<=79) then
                   pr_h_a_sim(h_l,2)=pr_h_a_sim(h_l,2)+data_habits_sim(i_l,h_l,g_l)
                   pr_h_a_data(h_l,2)=pr_h_a_data(h_l,2)+data_habits(i_l,h_l,g_l)
                   obs_h_a(h_l,2)=obs_h_a(h_l,2)+1
               end if
            end if
        end do; end do
        end do
        pr_h_e_sim=pr_h_e_sim/dble(obs_h_e)
        pr_h_e_data=pr_h_e_data/dble(obs_h_e)
        pr_h_a_sim=pr_h_a_sim/dble(obs_h_a)
        pr_h_a_data=pr_h_a_data/dble(obs_h_a)
        !Variance
        var_h_a_sim=pr_h_a_sim*(1.0d0-pr_h_a_sim)
        var_h_a_data=pr_h_a_data*(1.0d0-pr_h_a_data) 
        
        !covariance between t and t+4
        cov_h_a_sim=0.0d0
        cov_h_a_data=0.0d0
        obs_h_a=0
        do i_l=1,indv_HRS;
            if (last_age(i_l)>=first_age(i_l)+2) then
                do g_l=first_age(i_l),last_age(i_l)-2;do h_l=1,habits
                    if (sample_k(i_l,g_l)/=-1 .and. (data_habits(i_l,h_l,g_l)==1 .or. data_habits(i_l,h_l,g_l)==0) .and. race(i_l)==1 .and. &
                        sample_k(i_l,g_l+2)/=-1 .and. (data_habits(i_l,h_l,g_l+2)==1 .or. data_habits(i_l,h_l,g_l+2)==0)) then
                       age=initial_age+(g_l-1)*2
                       if (age>=65 .and. age<=69) then
                           cov_h_a_sim(h_l,1)=cov_h_a_sim(h_l,1)+(data_habits_sim(i_l,h_l,g_l)-pr_h_a_sim(h_l,1))*(data_habits_sim(i_l,h_l,g_l+2)-pr_h_a_sim(h_l,1))
                           cov_h_a_data(h_l,1)=cov_h_a_data(h_l,1)+(data_habits(i_l,h_l,g_l)-pr_h_a_data(h_l,1))*(data_habits(i_l,h_l,g_l+2)-pr_h_a_data(h_l,1))
                           obs_h_a(h_l,1)=obs_h_a(h_l,1)+1
                       elseif (age>=75 .and. age<=79) then
                           cov_h_a_sim(h_l,2)=cov_h_a_sim(h_l,2)+(data_habits_sim(i_l,h_l,g_l)-pr_h_a_sim(h_l,2))*(data_habits_sim(i_l,h_l,g_l+2)-pr_h_a_sim(h_l,2))
                           cov_h_a_data(h_l,2)=cov_h_a_data(h_l,2)+(data_habits(i_l,h_l,g_l)-pr_h_a_data(h_l,2))*(data_habits(i_l,h_l,g_l+2)-pr_h_a_data(h_l,2))
                           obs_h_a(h_l,2)=obs_h_a(h_l,2)+1
                       end if
                    end if
                end do; end do
            end if
        end do
        cov_h_a_sim=cov_h_a_sim/obs_h_a
        cov_h_a_data=cov_h_a_data/obs_h_a
        cov_h_a_sim=cov_h_a_sim/var_h_a_sim
        cov_h_a_data=cov_h_a_data/var_h_a_data
        
        !Cross-correlations
        cov_h1_h2_sim=0.0d0
        cov_h1_h2_data=0.0d0
        obs_h1_h2=0
        
        do i_l=1,indv_HRS;
            do g_l=first_age(i_l),last_age(i_l);do h_l=1,habits;do h2_l=1,habits
                if (sample_k(i_l,g_l)/=-1 .and. (data_habits(i_l,h_l,g_l)==1 .or. data_habits(i_l,h_l,g_l)==0) .and. race(i_l)==1 .and. &
                    (data_habits(i_l,h2_l,g_l)==1 .or. data_habits(i_l,h2_l,g_l)==0)) then
                    age=initial_age+(g_l-1)*2
                    if (age>=65 .and. age<=69) then
                        cov_h1_h2_sim(h_l,h2_l,1)=cov_h1_h2_sim(h_l,h2_l,1)+(data_habits_sim(i_l,h_l,g_l)-pr_h_a_sim(h_l,1))*(data_habits_sim(i_l,h2_l,g_l)-pr_h_a_sim(h2_l,1))
                        cov_h1_h2_data(h_l,h2_l,1)=cov_h1_h2_data(h_l,h2_l,1)+(data_habits(i_l,h_l,g_l)-pr_h_a_data(h_l,1))*(data_habits(i_l,h2_l,g_l)-pr_h_a_data(h2_l,1))
                        obs_h1_h2(h_l,h2_l,1)=obs_h1_h2(h_l,h2_l,1)+1
                    elseif (age>=75 .and. age<=79) then
                        cov_h1_h2_sim(h_l,h2_l,2)=cov_h1_h2_sim(h_l,h2_l,2)+(data_habits_sim(i_l,h_l,g_l)-pr_h_a_sim(h_l,2))*(data_habits_sim(i_l,h2_l,g_l)-pr_h_a_sim(h2_l,2))
                        cov_h1_h2_data(h_l,h2_l,2)=cov_h1_h2_data(h_l,h2_l,2)+(data_habits(i_l,h_l,g_l)-pr_h_a_data(h_l,2))*(data_habits(i_l,h2_l,g_l)-pr_h_a_data(h2_l,2))
                        obs_h1_h2(h_l,h2_l,2)=obs_h1_h2(h_l,h2_l,2)+1
                    end if
                end if
            end do; end do;end do
        end do
        cov_h1_h2_sim=cov_h1_h2_sim/obs_h1_h2
        cov_h1_h2_data=cov_h1_h2_data/obs_h1_h2
        
        do h2_l=1,habits;do h_l=1,habits
            cov_h1_h2_sim(h_l,h2_l,:)=cov_h1_h2_sim(h_l,h2_l,:)/sqrt(var_h_a_sim(h_l,:)*var_h_a_sim(h2_l,:))
            cov_h1_h2_data(h_l,h2_l,:)=cov_h1_h2_data(h_l,h2_l,:)/sqrt(var_h_a_data(h_l,:)*var_h_a_data(h2_l,:))
        end do;end do
        
        !Create Tables
        table_name='data'
        call create_table_1(pr_h_e_data,pr_h_a_data,cov_h_a_data,table_name)
        call create_table_2(cov_h1_h2_data,table_name)
        table_name='sim'//types_s
        call create_table_1(pr_h_e_sim,pr_h_a_sim,cov_h_a_sim,table_name)
        call create_table_2(cov_h1_h2_sim,table_name)
        
        !Simulate lives for computing sd
        call sim_lives(joint_yh,H)
        
        
end subroutine
    
subroutine create_table_1(pr_h_e,pr_h_a,cov_h_a,table_name) 
use global_var;use nrtype
implicit none
real(DP),dimension(habits,L_educ),intent(in)::pr_h_e
real(DP),dimension(habits,2),intent(in)::pr_h_a,cov_h_a
character(len=4),intent(in)::table_name
integer,dimension(habits)::order
character(len=11),dimension(habits)::str_habit
integer::h_l,h_l2

order=(/1,4,5,2,3,6/) !(/2,3,1,4,5,6/)
str_habit=(/'Cancer test','Drinking','Smoking','Cholesterol','Flu shot','Exercise'/)
open(unit=10,file=path_s2//'table_1_'//table_name//'.txt')
    write(10,'(A100)'),' & \multicolumn{5}{c}{Mean} & \multicolumn{2}{c}{AC} \\'
    write(10,'(A100)'),'\cmidrule(lr){2-6}\cmidrule(lr){7-8} '
    write(10,'(A100)'),'  & \textsc{hsd} & \textsc{hsg} & \textsc{cg} &  65-70 & 75-80 &  65-70 & 75-80 \\'
    write(10,'(A100)'),'\cmidrule(lr){2-4}\cmidrule(lr){5-6} \cmidrule(lr){7-8} '
    do h_l2=1,habits
        h_l=order(h_l2)
        write(10,'(A8,A15,A2,F5.2,A2,F5.2,A2,F5.2,A2,F5.2,A2,F5.2,A2,F5.2,A2,F5.2,A2)'), &
        '\text{',str_habit(h_l),'}&',pr_h_e(h_l,1),'&',pr_h_e(h_l,2),'&',pr_h_e(h_l,3),'&' &
            ,pr_h_a(h_l,1),'&',pr_h_a(h_l,2),'&',cov_h_a(h_l,1),'&',cov_h_a(h_l,2),'\\'
    end do
    write(10,'(A100)'),' \bottomrule'
close(10)


end subroutine

subroutine create_table_2(cov_h1_h2,table_name) 
use global_var;use nrtype
implicit none
real(DP),dimension(habits,habits,2),intent(in)::cov_h1_h2
character(len=4),intent(in)::table_name
integer,dimension(habits)::order
character(len=11),dimension(habits)::str_habit
integer::h_l2,h_l2o,h_l1,h_l1o

order=(/2,3,1,4,5,6/)
str_habit=(/'Cancer test','Drinking','Smoking','Cholesterol','Flu shot','Exercise'/)
open(unit=10,file=path_s2//'table_2_'//table_name//'.txt')
do h_l2=1,habits
    h_l2o=order(h_l2)
    if (h_l2<habits) then
        write(10,'(A10,A11,A1)'),'& \text{',str_habit(h_l2o),'} '  
    else
        write(10,'(A10,A11,A4)'),'& \text{',str_habit(h_l2o),'} \\'
    end if
end do
write(10,'(A20)'),'\midrule' 
do h_l1=1,habits
    h_l1o=order(h_l1)
    write(10,'(A8,A11,A1)'),'\text{',str_habit(h_l1o),'} '  
    do h_l2=1,habits
        h_l2o=order(h_l2) 
        if (h_l1>h_l2) then
            if (h_l2<habits) then
                write(10,'(A2,F5.2)'),'&', cov_h1_h2(h_l1o,h_l2o,1) 
            else
                write(10,'(A2,F5.2,A2)'),'&', cov_h1_h2(h_l1o,h_l2o,1),'\\'
            end if
        elseif (h_l1==h_l2) then
            if (h_l2<habits) then
                write(10,'(A22,F5.2,A1)'),'& \textcolor{gray!90}{', cov_h1_h2(h_l1o,h_l2o,2),'}'
            else
                write(10,'(A22,F5.2,A3)'),'& \textcolor{gray!90}{', cov_h1_h2(h_l1o,h_l2o,2),'}\\'
            end if
        else
            if (h_l2<habits) then
                write(10,'(A21,F5.2)'),'& \cellcolor{gray!20}', cov_h1_h2(h_l1o,h_l2o,2) 
            else
                write(10,'(A21,F5.2,A2)'),'& \cellcolor{gray!20}', cov_h1_h2(h_l1o,h_l2o,2),'\\'
            end if
        end if
    end do
end do
write(10,'(A20)'),'\bottomrule' 

close(10)


end subroutine    
    
    
subroutine create_table_3(LE,joint_yh)
use global_var;use nrtype
implicit none 
real(DP),dimension(types,L_gender,L_educ,clusters+1),intent(in)::LE
real(DP),dimension(generations,clusters,L_gender,L_educ,types,cohorts),intent(in)::joint_yh
character(len=30),dimension(types)::group_name
character(len=1)::y_l_s
character(len=200),dimension(types)::group_s
character(len=200)::all_s,delta_s
integer::y_l,e_l
real(DP),dimension(2,L_educ,types)::columns
real(DP),dimension(2,types)::columns_all
real(DP),dimension(L_educ)::pr_e
real(DP),dimension(L_educ,types)::Pr_e_y
character(len=150) :: fmt

print*,'Adjust pr of education as needed'
pr_e=(/0.11d0,0.52d0,0.37d0/)   
if (types==2) then
    group_name(1)='\ \ \textsc{pro}'
    group_name(2)='\ \ \textsc{det}'
else
    do y_l=1,types
        write(y_l_s, '(I1)') y_l
        group_name(y_l)='\ \ Group '//y_l_s
    end do
end if

fmt = '(A30,A3,F10.1,A3,F10.1,A3,F10.1,A3,F10.1,A3,F10.1,A3,F10.1,A3,F10.1,A3,F10.1,A3,F10.1,A3,F10.1,A3,F10.1,A3)'

open(unit=10,file=path_s2//'table_3_new'//types_s//'.txt')
write(10,'(A220)'),'& \multicolumn{2}{c}{All} & \multicolumn{2}{c}{\textsc{hsd}} &	\multicolumn{2}{c}{\textsc{hsg}} & \multicolumn{2}{c}{\textsc{cg}} & \multicolumn{3}{c}{$\Delta_{\text{e}}$LE (\textsc{cg}-\textsc{hsd})} \\'
write(10,'(A100)'),'\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){6-7} \cmidrule(lr){8-9}\cmidrule(lr){10-12} '
write(10,'(A200)'),'   &\% 	&	LE  &\% 	&	LE 	&	\% 	&	LE 	&	\% 	&	LE 	& Data	 	& (a)	& (b) 	\\'
write(10,'(A100)'),'\midrule'
write(10,'(A100)'),' \\[-2ex]'

do y_l=1,types
    do e_l=1,L_educ
        columns(1,e_l,y_l)=sum(joint_yh(1,:,1,e_l,y_l,4))
        columns(2,e_l,y_l)=LE(y_l,1,e_l,3) 
    end do
    columns_all(1,y_l)=sum(columns(1,:,y_l)*Pr_e)
    Pr_e_y(:,y_l)=columns(1,:,y_l)*Pr_e/sum(columns(1,:,y_l)*Pr_e)
    columns_all(2,y_l)=sum(Pr_e_y(:,y_l)*columns(2,:,y_l))
    write(group_s(y_l),fmt), group_name(y_l) ,'&', nint(columns_all(1,y_l)*1000.0)/10.0d0, ' & ', nint(columns_all(2,y_l)*10.0)/10.0, &
                                    ' & ', nint(columns(1,1,y_l)*1000.0)/10.0d0, ' & ', nint(columns(2,1,y_l)*10.0)/10.0, &
                                    ' & ', nint(columns(1,2,y_l)*1000.0)/10.0d0, ' & ', nint(columns(2,2,y_l)*10.0)/10.0, &
                                    ' & ', nint(columns(1,3,y_l)*1000.0)/10.0d0, ' & ', nint(columns(2,3,y_l)*10.0)/10.0,'\\'
end do
if (types==2) then
    write(delta_s,fmt), '\ \ $\Delta_y$ ' ,'&', nint((columns_all(1,1)-columns_all(1,2))*1000.0)/10.0d0, ' & ', nint((columns_all(2,1)-columns_all(2,2))*10.0)/10.0, &
                                    ' & ', nint((columns(1,1,1)-columns(1,1,2))*1000.0)/10.0d0, ' & ', nint((columns(2,1,1)-columns(2,1,2))*10.0)/10.0, &
                                    ' & ', nint((columns(1,2,1)-columns(1,2,2))*1000.0)/10.0d0, ' & ', nint((columns(2,2,1)-columns(2,2,2))*10.0)/10.0, &
                                    ' & ', nint((columns(1,3,1)-columns(1,3,2))*1000.0)/10.0d0, ' & ', nint((columns(2,3,1)-columns(2,3,2))*10.0)/10.0,'\\'
end if
    
    write(all_s,fmt), 'All' ,' & ', 100.0, ' & ', nint(sum(columns_all(2,:)*columns_all(1,:))*10.0)/10.0, &
                            ' & ', 100.0, ' & ', nint(sum(columns(2,1,:)*columns(1,1,:))*10.0)/10.0, &
                            ' & ', 100.0, ' & ', nint(sum(columns(2,2,:)*columns(1,2,:))*10.0)/10.0, &
                            ' & ', 100.0, ' & ', nint(sum(columns(2,3,:)*columns(1,3,:))*10.0)/10.0, &
                            ' & ',nint((sum(columns(2,3,:)*columns(1,3,:))-sum(columns(2,1,:)*columns(1,1,:)))*10.0)/10.0, &
                            ' & ',nint((sum(columns(2,3,:)*columns(1,3,:))-sum(columns(2,1,:)*columns(1,3,:)))*10.0)/10.0 , & 
                            ' & ',nint(((sum(columns(2,3,:)*columns(1,3,:))-sum(columns(2,1,:)*columns(1,1,:)))-(sum(columns(2,3,:)*columns(1,3,:))-sum(columns(2,1,:)*columns(1,3,:))))*10.0)/10.0 , '\\'
    write(10,'(A)'),all_s
    do y_l=1,types
        write(10,'(A)'),group_s(y_l)
    end do
    if (types==2) then
        write(10,'(A)'),delta_s
    end if
        
    write(10,'(A100)'),'\bottomrule'
close(10)


end subroutine
    
