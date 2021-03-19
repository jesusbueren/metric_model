subroutine charge_data()
    use global_var
    implicit none 
    integer,dimension(indv*adls*generations)::data_vec_adls
    integer,dimension(indv*habits*generations)::data_vec_habits
    integer,dimension(indv*2)::ages
    integer::i_l,g_l,index
    integer,dimension(generations,indv)::data_shlt_rsh
    
    open(unit=10,file=path//"data\ages_all.csv")
        read(10,*) ages
    close(10)
    open(unit=10,file=path//"data\gender_all.csv")
        read(10,*) gender
    close(10)
    open(unit=10,file=path//"Data\adls.csv")
        read(10,*) data_vec_adls
    close(10)
    open(unit=10,file=path//"Data\habits.csv")
        read(10,*) data_vec_habits
    close(10)
    open(unit=10,file=path//"data\educ_all.csv")
        read(10,*) educ
    close(10)
    open(unit=10,file=path//"Data\shlt.csv")
        read(10,*) data_shlt_rsh
    close(10)
    data_shlt=reshape(data_shlt_rsh,(/indv,generations/),order=(/2,1/))
    
    do i_l=1,indv
        do g_l=1,generations 
            if (data_shlt(i_l,g_l)>=1 .and. data_shlt(i_l,g_l)<=3) then
                data_shlt(i_l,g_l)=1
            elseif (data_shlt(i_l,g_l)>=4 .and. data_shlt(i_l,g_l)<=5) then
                data_shlt(i_l,g_l)=2
            end if                
        end do    
    end do
    
    do i_l=1,indv
        do g_l=1,generations 
            index=(i_l-1)*(generations*adls)+(g_l-1)*adls+1
            data_adls(i_l,:,g_l)=data_vec_adls(index:index+adls-1)
        end do
    end do
    
    do i_l=1,indv
        do g_l=1,generations 
            index=(i_l-1)*(generations*habits)+(g_l-1)*habits+1
            data_habits(i_l,:,g_l)=data_vec_habits(index:index+habits-1)
        end do
    end do
    
    high_school=0
    college=0
    do i_l=1,indv
        index=(i_l-1)*2+1
        first_age(i_l)=ages(index)
        last_age(i_l)=ages(index+1)
        if (educ(i_l)==2) then
            high_school(i_l)=1
        elseif (educ(i_l)==3) then
            college(i_l)=1
        end if
    end do
end subroutine
