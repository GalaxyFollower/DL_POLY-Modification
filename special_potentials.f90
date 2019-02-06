!**********************************************************************************************************************!
!                           Modification by Mian Muhammad Faheem and Mehdi Zare                                        !
!**********************************************************************************************************************!
module special_potentials

  use kinds_f90
  implicit none

  public

  !Variables for modified energy/gradient.
  real(kind=wp) :: feng, fgam

  interface spohr_h_numerical
    module procedure spohr_h_numerical_wpX3
  end interface spohr_h_numerical

  interface spohr_o_numerical
    module procedure spohr_o_numerical_wpX5
  end interface spohr_o_numerical

  contains

  !******************************************************************************************************************!
  !                                           INTERFACE SPOHR_H_NUMERICAL                                            !
  !******************************************************************************************************************!
  subroutine spohr_h_numerical_wpX3 (eng, gam, rr)
    use kinds_f90
    implicit none

    real(kind=wp), intent(in)  :: rr
    real(kind=wp), intent(out) :: eng, gam

    real(kind=wp) :: a4, b4

    !Set model parameters in DLPOLY internal units.
   !*************************************************************************!
   !                Modification by Mehdi Zare                               !
   !  DL_POLY internal unit for energy: 1.6605402x10^-23 jouls               !
   !  SH unit of energy: 10^-19 jouls                                        !
   !  SH to DL_POLY unit Example for a4: (1.7142x10^-19)/(1.6605402x10^-23)  !
   !*************************************************************************!                 
    a4  = 10323.1466_wp
    b4  = 1.2777_wp
    !Calculate energy at given point.
    eng = a4 * exp(-b4*rr)  
 
    !Calculate gam = -1/r(dU/dr) 
    gam = (a4*b4/rr)*exp(-b4*rr)
   !*************************************************************************!
   !                 End of  Modification by Mehdi Zare                      !
   !*************************************************************************!
  end subroutine spohr_h_numerical_wpX3
  !******************************************************************************************************************!

  !******************************************************************************************************************!
  !                                           INTERFACE SPOHR_O_NUMERICAL                                            !
  !******************************************************************************************************************!
  subroutine spohr_o_numerical_wpX5 (eng, gam, gamrho, rr, dx, dy)
    use kinds_f90
    implicit none
    
    real(kind=wp), intent(in)  :: rr, dx, dy
    real(kind=wp), intent(out) :: eng, gam
   !*************************************************************************!
   !                     Modification by Mehdi Zare                          !
   !*************************************************************************! 
    real(kind=wp), intent(out) :: gamrho
   !*************************************************************************!
   !                 End of  Modification by Mehdi Zare                      !
   !*************************************************************************!
    real(kind=wp) :: a1, a2, a3, b1, b2, b3, cc
    real(kind=wp) :: t1, xp

    !Set model parameters in DLPOLY internal units.
    a1  = 11407131.2456_wp
    a2  = 11359556.3661_wp
    a3  = 6022136651.6752_wp
    b1  = 1.1004_wp
    b2  = 1.0966_wp
    b3  = 5.3568_wp
    cc  = 0.5208_wp
   !*************************************************************************!
   !                     Modification by Mehdi Zare                          !
   !    For Pt-O the potential is a function of both r and rho, So for force !
   !              calculation we perform some algebra:                       !
   !           Fx = -dU/dx = -(dU/dr)*(dr/dx) - (dU/drho)*(drho/dx)          !
   !                  ... Fx = gam*xij + gamrho*xij                          !
   !                  ... Fy = gam*yij + gamrho*yij                          !
   !                  ... Fz = gam*zij + 0                                   !
   !*************************************************************************!
    ! t1 is rho square in the SH paper
    ! xp is f(rho) in the SH paper
    t1  = dx*dx + dy*dy
    xp  = exp(-cc*t1)
    
    ! Calculare energy at a given point
    eng = (a1*exp(-b1*rr)-a2*exp(-b2*rr))*xp + a3*exp(-b3*rr)*(1.0_wp-xp)
 
    !Calculate -1/r(dU/dr) : Modification by Mehdi Zare:for Ex.-dU/dx=-(dU/dr)(dr/dx) = -(x/r)(dU/dr)! 
    gam = (-a1*b1*exp(-b1*rr)+a2*b2*exp(-b2*rr))*xp + (-a3*b3*exp(-b3*rr))*(1.0_wp-xp)
    gam = -gam/rr
    
    !Calculate -1/rho(dU/drho)
    gamrho = 2.0_wp*cc*sqrt(t1)*xp * (a3*exp(-b3*rr) + a2*exp(-b2*rr) - a1*exp(-b1*rr))
    gamrho = -gamrho/sqrt(t1)
   !*************************************************************************!
   !                 End of  Modification by Mehdi Zare                      !
   !*************************************************************************!    
  end subroutine spohr_o_numerical_wpX5

end module special_potentials
!**********************************************************************************************************************!
!                                                 End of Modification                                                  !
!**********************************************************************************************************************!
