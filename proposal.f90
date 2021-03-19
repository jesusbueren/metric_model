subroutine proposal(c,c_g,C_mat,dim)
    implicit none
    integer,intent(in)::dim
    double precision,intent(in),dimension(dim,1)::c
    double precision,intent(out),dimension(dim,1)::c_g
    double precision,intent(in),dimension(dim,dim)::C_mat
    double precision,dimension(dim,1)::normal_d
    integer::v_l
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    character::continue_program
    
    do v_l=1,dim
        normal_d(v_l,1)=c4_normal_01()
    end do
    
    c_g=c+matmul(C_mat,normal_d)  
    
end subroutine
    
double precision function c4_normal_01 (  )
!------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  double precision, parameter :: r4_pi=3.14159265358979323846264338327950288419716939937510
  double precision:: v1
  double precision:: v2
  double precision:: x_c
  double precision:: x_r 
  call random_number(v1)
  call random_number(v2)
  x_r =sqrt(-2.0d0*log(v1))*cos(2.0d0*r4_pi*v2)
  x_c =sqrt(-2.0d0*log(v1))*sin(2.0d0*r4_pi*v2)
  c4_normal_01=x_r
  return
end