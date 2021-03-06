subroutine charge_data()
    use global_var
    implicit none 
    integer,dimension(indv_HRS*habits*generations)::data_vec_habits
    integer,dimension(2,indv)::ages
    integer,dimension(generations,indv_HRS)::data_shlt_srh
    integer,dimension(generations,indv_psid)::data_shlt_psid
    
    integer::i_l,g_l,index
    
    
    open(unit=10,file=path//"data\ages_all.csv")
        read(10,*) ages(:,1:indv_HRS)
    close(10)
    open(unit=10,file=path//"data\ages_all_psid.csv")
        read(10,*) ages(:,indv_HRS+1:indv)
    close(10)
    open(unit=10,file=path//"data\gender_all.csv")
        read(10,*) gender(1:indv_HRS)
    close(10)
    open(unit=10,file=path//"data\gender_all_psid.csv")
        read(10,*) gender(indv_HRS+1:indv)
    close(10)
    open(unit=10,file=path//"Data\habits.csv")
        read(10,*) data_vec_habits
    close(10)
    open(unit=10,file=path//"data\educ_all.csv")
        read(10,*) educ(1:indv_HRS)
    close(10)
    open(unit=10,file=path//"data\educ_all_psid.csv")
        read(10,*) educ(indv_HRS+1:indv)
    close(10)
    open(unit=10,file=path//"Data\shlt.csv")
        read(10,*) data_shlt_srh
    close(10)
    open(unit=10,file=path//"Data\shlt_psid.csv")
        read(10,*) data_shlt_psid
    close(10)
    
    data_shlt(1:indv_HRS,:)=reshape(data_shlt_srh,(/indv_HRS,generations/),order=(/2,1/))
    data_shlt(indv_HRS+1:indv,:)=reshape(data_shlt_psid,(/indv_psid,generations/),order=(/2,1/))
    
    do i_l=1,indv
        do g_l=1,generations 
            if (data_shlt(i_l,g_l)>=1 .and. data_shlt(i_l,g_l)<=3) then
                data_shlt(i_l,g_l)=1
            elseif (data_shlt(i_l,g_l)>=4 .and. data_shlt(i_l,g_l)<=5) then
                data_shlt(i_l,g_l)=2
            elseif (data_shlt(i_l,g_l)>=6) then
                data_shlt(i_l,g_l)=3
            end if                
        end do    
    end do
    
    
    do i_l=1,indv
        if (i_l<=indv_HRS) then
            do g_l=1,generations 
                index=(i_l-1)*(generations*habits)+(g_l-1)*habits+1
                data_habits(i_l,:,g_l)=data_vec_habits(index:index+habits-1)
            end do
        else
            data_habits(i_l,:,:)=-9
        end if
    end do

    high_school=0
    college=0
    do i_l=1,indv
        first_age(i_l)=ages(1,i_l)
        last_age(i_l)=ages(2,i_l)
        if (educ(i_l)==2) then
            high_school(i_l)=1
        elseif (educ(i_l)==3) then
            college(i_l)=1
        end if
    end do
end subroutine
