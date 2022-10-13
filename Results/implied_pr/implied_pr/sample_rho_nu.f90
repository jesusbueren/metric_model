subroutine sample_rho_nu(u_draw,s2_nu,rho,s2_nu_new)
    use global_var;use nrtype
    implicit none
    double precision,dimension(indv,generations),intent(in)::u_draw
    double precision,dimension(L_educ),intent(in)::s2_nu
    double precision,dimension(L_educ,1),intent(out)::rho
    double precision,dimension(L_educ),intent(out)::s2_nu_new
    integer::ind,t_l,i_l,e_l
    double precision,dimension(indv*generations,1)::x_u,y_u,e
    double precision,dimension(1,1)::rho_hat
    double precision::v,s2,shape,scale
    interface
        double precision function c4_normal_01( )
            implicit none
        end function c4_normal_01
    end interface
    interface
        double precision function r8_gamma_01_sample ( shape )
            implicit none
            double precision,intent(in)::shape
        end function r8_gamma_01_sample
    end interface
    
    !Sample rho
    do e_l=1,L_educ
    ind=0
    do i_l=indv_HRS+1,indv
        do t_l=first_age(i_l),last_age(i_l)-1
            if (u_draw(i_l,t_l)/=-1.0d0 .and. u_draw(i_l,t_l+1)/=-1.0d0 .and. gender(i_l)==1 .and. initial_age+(t_l-1)*2<64 .and. educ(i_l)==e_l ) then 
                ind=ind+1
                x_u(ind,1)=u_draw(i_l,t_l)
                y_u(ind,1)=u_draw(i_l,t_l+1)
            end if
        end do
    end do
    
    rho_hat=matmul(1/matmul(transpose(x_u(1:ind,1:1)),x_u(1:ind,1:1)),matmul(transpose(x_u(1:ind,1:1)),y_u(1:ind,1:1)))
    rho(e_l:e_l,1:1)=rho_hat(1,1)+c4_normal_01( )*sqrt(1/matmul(transpose(x_u(1:ind,1:1)),x_u(1:ind,1:1))*s2_nu(e_l))
    
    !Sample s2_nu
    v=dble(ind-1)
    e(1:ind,1:1)=y_u(1:ind,1:1)-rho(e_l,1)*x_u(1:ind,1:1) 
    s2=sum(e(1:ind,1:1)**2)/v
    shape=v/2
    scale=1/(v*s2/2)
    s2_nu_new(e_l)=1/(r8_gamma_01_sample(shape)*scale)
    end do
    
    end subroutine

    
double precision function r8_gamma_01_sample ( a )

!*****************************************************************************80
!
!! R8_GAMMA_01_SAMPLE samples the standard Gamma distribution.
!
!  Discussion:
!
!    This procedure corresponds to algorithm GD in the reference.
!
!    pdf ( a; x ) = 1/gamma(a) * x^(a-1) * exp ( - x )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by Barry Brown, James Lovato.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Joachim Ahrens, Ulrich Dieter,
!    Generating Gamma Variates by a Modified Rejection Technique,
!    Communications of the ACM,
!    Volume 25, Number 1, January 1982, pages 47-54.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the shape parameter. 
!    0.0 < A.
!
!    Output, real ( kind = 8 ) R8_GAMMA_01_SAMPLE, a random deviate 
!    from the distribution.
!
  implicit none
  double precision,intent(in):: a
  double precision, parameter :: a1 =  0.3333333D+00
  double precision, parameter :: a2 = -0.2500030D+00
  double precision, parameter :: a3 =  0.2000062D+00
  double precision, parameter :: a4 = -0.1662921D+00
  double precision, parameter :: a5 =  0.1423657D+00
  double precision, parameter :: a6 = -0.1367177D+00
  double precision, parameter :: a7 =  0.1233795D+00
  double precision:: b
  double precision:: c
  double precision:: d
  double precision:: e
  double precision, parameter :: e1 = 1.0D+00
  double precision, parameter :: e2 = 0.4999897D+00
  double precision, parameter :: e3 = 0.1668290D+00
  double precision, parameter :: e4 = 0.0407753D+00
  double precision, parameter :: e5 = 0.0102930D+00
  double precision:: p
  double precision:: q
  double precision:: q0
  double precision, parameter :: q1 =  0.04166669D+00
  double precision, parameter :: q2 =  0.02083148D+00
  double precision, parameter :: q3 =  0.00801191D+00
  double precision, parameter :: q4 =  0.00144121D+00
  double precision, parameter :: q5 = -0.00007388D+00
  double precision, parameter :: q6 =  0.00024511D+00
  double precision, parameter :: q7 =  0.00024240D+00
  double precision:: r
  double precision:: r8_exponential_01_sample
  double precision:: r8_normal_01_sample
  double precision:: r8_uniform_01_sample
  double precision:: s
  double precision:: s2
  double precision:: si
  double precision, parameter :: sqrt32 = 5.656854D+00
  double precision:: t
  double precision:: u
  double precision:: v
  double precision:: w
  double precision:: x

  if ( 1.0D+00 <= a ) then

    s2 = a - 0.5D+00
    s = sqrt ( s2 )
    d = sqrt32 - 12.0D+00 * s
!
!  Immediate acceptance.
!
    t = r8_normal_01_sample ( )
    x = s + 0.5D+00 * t
    r8_gamma_01_sample = x * x

    if ( 0.0D+00 <= t ) then
      return
    end if
!
!  Squeeze acceptance.
!
    u = r8_uniform_01_sample ( )
    if ( d * u <= t * t * t ) then
      return
    end if

    r = 1.0D+00 / a
    q0 = (((((( q7 &
      * r + q6 ) &
      * r + q5 ) &
      * r + q4 ) &
      * r + q3 ) &
      * r + q2 ) &
      * r + q1 ) &
      * r
!
!  Approximation depending on size of parameter A.
!
    if ( 13.022D+00 < a ) then
      b = 1.77D+00
      si = 0.75D+00
      c = 0.1515D+00 / s
    else if ( 3.686D+00 < a ) then
      b = 1.654D+00 + 0.0076D+00 * s2
      si = 1.68D+00 / s + 0.275D+00
      c = 0.062D+00 / s + 0.024D+00
    else
      b = 0.463D+00 + s + 0.178D+00 * s2
      si = 1.235D+00
      c = 0.195D+00 / s - 0.079D+00 + 0.16D+00 * s
    end if
!
!  Quotient test.
!
    if ( 0.0D+00 < x ) then

      v = 0.5D+00 * t / s

      if ( 0.25D+00 < abs ( v ) ) then
        q = q0 - s * t + 0.25D+00 * t * t + 2.0D+00 * s2 * log ( 1.0D+00 + v )
      else
        q = q0 + 0.5D+00 * t * t * (((((( a7 &
          * v + a6 ) &
          * v + a5 ) &
          * v + a4 ) &
          * v + a3 ) &
          * v + a2 ) &
          * v + a1 ) &
          * v
      end if

      if ( log ( 1.0D+00 - u ) <= q ) then
        return
      end if

    end if

    do

      e = r8_exponential_01_sample ( )
      u = 2.0D+00 * r8_uniform_01_sample ( ) - 1.0D+00
 
      if ( 0.0D+00 <= u ) then
        t = b + abs ( si * e )
      else
        t = b - abs ( si * e )
      end if
!
!  Possible rejection.
!
      if ( t < -0.7187449D+00 ) then
        cycle
      end if
!
!  Calculate V and quotient Q.
!
      v = 0.5D+00 * t / s

      if ( 0.25D+00 < abs ( v ) ) then
        q = q0 - s * t + 0.25D+00 * t * t + 2.0D+00 * s2 * log ( 1.0D+00 + v )
      else
        q = q0 + 0.5D+00 * t * t * (((((( a7 &
          * v + a6 ) &
          * v + a5 ) &
          * v + a4 ) &
          * v + a3 ) &
          * v + a2 ) &
          * v + a1 ) &
          *  v
      end if
!
!  Hat acceptance.
!
      if ( q <= 0.0D+00 ) then
        cycle
      end if

      if ( 0.5D+00 < q ) then
        w = exp ( q ) - 1.0D+00
      else
        w = (((( e5 * q + e4 ) * q + e3 ) * q + e2 ) * q + e1 ) * q
      end if
!
!  May have to sample again.
!
      if ( c * abs ( u ) <= w * exp ( e - 0.5D+00 * t * t ) ) then
        exit
      end if

    end do

    x = s + 0.5D+00 * t
    r8_gamma_01_sample = x * x

    return
!
!  Method for A < 1.
!
  else

    b = 1.0D+00 + 0.3678794D+00 * a

    do

      p = b * r8_uniform_01_sample ( )

      if ( p < 1.0D+00 ) then

        r8_gamma_01_sample = exp ( log ( p ) / a )

        if ( r8_gamma_01_sample <= r8_exponential_01_sample ( ) ) then
          return
        end if

        cycle

      end if

      r8_gamma_01_sample = - log ( ( b - p ) / a )

      if ( ( 1.0D+00 - a ) * log ( r8_gamma_01_sample ) <= &
        r8_exponential_01_sample ( ) ) then
        exit
      end if

    end do

  end if

  return
    end
    
    function r8_uniform_01_sample ( )

!*****************************************************************************80
!
!! R8_UNIFORM_01_SAMPLE generates a uniform random deviate from [0,1].
!
!  Discussion:
!
!    This function should be the only way that the package accesses random
!    numbers.
!
!    Setting OPTION to 0 accesses the R8_UNI_01() function in RNGLIB,
!    for which there are versions in various languages, which should result
!    in the same values being returned.
!
!    Setting OPTION to 1 in the FORTRAN90 version calls the system
!    RNG "random_number()".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2013
!
!  Author:
!
!    Original FORTRAN77 version by Barry Brown, James Lovato.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01_SAMPLE, a random deviate 
!    from the distribution.
!
  implicit none

  integer ( kind = 4 ), parameter :: option = 0
  real ( kind = 8 ) r8_uniform_01_sample
  real ( kind = 8 ) value

 call random_number ( harvest = value )
 
  r8_uniform_01_sample = value

  return
end
    
    function r8_normal_01_sample ( )

!*****************************************************************************80
!
!! R8_NORMAL_01_SAMPLE returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    The Box-Muller method is used, which is efficient, but
!    generates two values at a time.
!
!    Typically, we would use one value and save the other for the next call.
!    However, the fact that this function has saved memory makes it difficult
!    to correctly handle cases where we want to re-initialize the code,
!    or to run in parallel.  Therefore, we will instead use the first value
!    and DISCARD the second.
!
!    EFFICIENCY must defer to SIMPLICITY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_NORMAL_01_SAMPLE, a sample of the standard
!    normal PDF.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01_sample
  real ( kind = 8 ) r8_uniform_01_sample
  real ( kind = 8 ) x

  r1 = r8_uniform_01_sample ( )
  r2 = r8_uniform_01_sample ( )

  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )

  r8_normal_01_sample = x

  return
    end
    
    function r8_exponential_01_sample ( )

!*****************************************************************************80
!
!! R8_EXPONENTIAL_01_SAMPLE samples the standard exponential PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EXPONENTIAL_01_SAMPLE, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) r8_exponential_01_sample
  real ( kind = 8 ) r8_uniform_01_sample

  r = r8_uniform_01_sample ( )

  r8_exponential_01_sample = - log ( r )

  return
end