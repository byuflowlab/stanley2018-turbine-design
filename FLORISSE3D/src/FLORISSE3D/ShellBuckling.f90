subroutine shellBucklingEurocode(nPoints, d, t, sigma_z, sigma_t, tau_zt, L_reinforced, E, sigma_y, &
  &gamma_f, gamma_b, EU_utilization)
  ! Estimate shell buckling utilization along tower.

  ! Arguments:
  ! npt - number of locations at each node at which stress is evaluated.
  ! sigma_z - axial stress at npt*node locations.  must be in order
  !               [(node1_pts1-npt), (node2_pts1-npt), ...]
  ! sigma_t - azimuthal stress given at npt*node locations
  ! tau_zt - shear stress (z, theta) at npt*node locations
  ! E - modulus of elasticity
  ! sigma_y - yield stress
  ! L_reinforced - reinforcement length - structure is re-discretized with this spacing
  ! gamma_f - safety factor for stresses
  ! gamma_b - safety factor for buckling
  !
  ! Returns:
  ! z
  ! EU_utilization: - array of shell buckling utilizations evaluted at (z[0] at npt locations, \n
  !                   z[0]+L_reinforced at npt locations, ...). \n
  !                   Each utilization must be < 1 to avoid failure.

  implicit none

  ! define precision to be the standard for a double precision ! on local system
  integer, parameter :: dp = kind(0.d0)

  ! in
  integer, intent(in) :: nPoints
  real(dp), intent(in) :: gamma_f, gamma_b
  real(dp), dimension(nPoints), intent(in) :: d, t, sigma_z, sigma_t, L_reinforced, E, sigma_y
  real(dp), dimension(nPoints), intent(in) :: tau_zt

  !out
  real(dp), dimension(nPoints), intent(out) :: EU_utilization

  !local
  integer :: i
  real(dp) :: h, r1, r2, t1, t2, sigma_z_shell, sigma_t_shell, tau_zt_shell, utilization
  real(dp), dimension(nPoints) :: sigma_z_sh, sigma_t_sh, tau_zt_sh

  do i = 1, nPoints
    h = L_reinforced(i)

    r1 = d(i)/2.0_dp - t(i)/2.0_dp
    r2 = d(i)/2.0_dp - t(i)/2.0_dp
    t1 = t(i)
    t2 = t(i)

    sigma_z_shell = sigma_z(i)
    sigma_t_shell = sigma_t(i)
    tau_zt_shell = tau_zt(i)

    ! TO DO: the following is non-smooth, although in general its probably OK
    ! change to magnitudes and add safety factor
    sigma_z_shell = gamma_f*abs(sigma_z_shell)
    sigma_t_shell = gamma_f*abs(sigma_t_shell)
    tau_zt_shell = gamma_f*abs(tau_zt_shell)

    call shellBucklingOneSection(h, r1, r2, t1, t2, gamma_b, sigma_z_shell, sigma_t_shell, &
                                &tau_zt_shell, E(i), sigma_y(i), utilization)
    EU_utilization(i) = utilization

    ! make them into vectors
    ! TODO is this necessary?
    sigma_z_sh(i)=sigma_z_shell
    sigma_t_sh(i)=sigma_t_shell
    tau_zt_sh(i)=tau_zt_shell

  end do

end subroutine shellBucklingEurocode


subroutine shellBucklingOneSection(h, r1, r2, t1, t2, gamma_b, sigma_z, sigma_t, tau_zt, E, sigma_y, utilization)

  ! Estimate shell buckling for one tapered cylindrical shell section.
  !
  ! Arguments:
  ! h - height of conical section
  ! r1 - radius at bottom
  ! r2 - radius at top
  ! t1 - shell thickness at bottom
  ! t2 - shell thickness at top
  ! E - modulus of elasticity
  ! sigma_y - yield stress
  ! gamma_b - buckling reduction safety factor
  ! sigma_z - axial stress component
  ! sigma_t - azimuthal stress component
  ! tau_zt - shear stress component (z, theta)
  !
  ! Returns:
  ! EU_utilization, shell buckling utilization which must be < 1 to avoid failure

  ! NOTE: definition of r1, r2 switched from Eurocode document to be consistent with FEM.


  implicit none

  ! define precision to be the standard for a double precision ! on local system
  integer, parameter :: dp = kind(0.d0)

  !in
  real(dp), intent(in) :: h, r1, r2, t1, t2, gamma_b, sigma_z, sigma_t, tau_zt, E, sigma_y

  !out
  real(dp), intent(out) :: utilization

  !local
  real(dp) :: beta, L, t, le, re, omega, rovert, Cx, sigma_z_Rcr, lambda_z0, beta_z, &
              &eta_z, Q, lambda_z, delta_wk, alpha_z, chi_z, sigma_z_Rk, sigma_z_Rd, &
              &sigma_t_Rcr, alpha_t, lambda_t0, beta_t, eta_t, lambda_t, chi_theta, &
              &sigma_t_Rk, sigma_t_Rd, rho, C_tau, tau_zt_Rcr, alpha_tau, beta_tau, &
              &lambda_tau0, eta_tau, lambda_tau, chi_tau, tau_zt_Rk, tau_zt_Rd, k_z, &
              &k_theta, k_tau, k_i

  ! ----- geometric parameters --------
  beta = atan2(r1-r2, h)
  L = h/cos(beta)
  t = 0.5_dp*(t1+t2)

  ! ------------- axial stress -------------
  ! length parameter
  le = L
  re = 0.5_dp*(r1+r2)/cos(beta)
  omega = le/sqrt(re*t)
  rovert = re/t

  ! compute Cx
  call cxsmooth(omega, rovert, Cx)


  ! if omega <= 1.7:
  !     Cx = 1.36 - 1.83/omega + 2.07/omega/omega
  ! elif omega > 0.5*rovert:
  !     Cxb = 6.0  ! clamped-clamped
  !     Cx = max(0.6, 1 + 0.2/Cxb*(1-2.0*omega/rovert))
  ! else:
  !     Cx = 1.0

  ! critical axial buckling stress
  sigma_z_Rcr = 0.605_dp*E*Cx/rovert

  ! compute buckling reduction factors
  lambda_z0 = 0.2_dp
  beta_z = 0.6_dp
  eta_z = 1.0_dp
  Q = 25.0_dp  ! quality parameter - high
  lambda_z = sqrt(sigma_y/sigma_z_Rcr)
  delta_wk = 1.0_dp/Q*sqrt(rovert)*t
  alpha_z = 0.62_dp/(1.0_dp + 1.91_dp*(delta_wk/t)**1.44)

  call buckling_reduction_factor(alpha_z, beta_z, eta_z, lambda_z0, lambda_z, chi_z)

  ! design buckling stress
  sigma_z_Rk = chi_z*sigma_y
  sigma_z_Rd = sigma_z_Rk/gamma_b

  ! ---------------- hoop stress ------------------

  ! length parameter
  le = L
  re = 0.5_dp*(r1+r2)/(cos(beta))
  omega = le/sqrt(re*t)
  rovert = re/t

  ! Ctheta = 1.5  ! clamped-clamped
  ! CthetaS = 1.5 + 10.0/omega**2 - 5.0/omega**3

  ! ! critical hoop buckling stress
  ! if (omega/Ctheta < 20.0):
  !     sigma_t_Rcr = 0.92*E*CthetaS/omega/rovert
  ! elif (omega/Ctheta > 1.63*rovert):
  !     sigma_t_Rcr = E*(1.0/rovert)**2*(0.275 + 2.03*(Ctheta/omega*rovert)**4)
  ! else:
  !     sigma_t_Rcr = 0.92*E*Ctheta/omega/rovert

  call sigmasmooth(omega, E, rovert, sigma_t_Rcr)

  ! buckling reduction factor
  alpha_t = 0.65_dp  ! high fabrication quality
  lambda_t0 = 0.4_dp
  beta_t = 0.6_dp
  eta_t = 1.0_dp
  lambda_t = sqrt(sigma_y/sigma_t_Rcr)

  call buckling_reduction_factor(alpha_t, beta_t, eta_t, lambda_t0, lambda_t, chi_theta)

  sigma_t_Rk = chi_theta*sigma_y
  sigma_t_Rd = sigma_t_Rk/gamma_b

  ! ----------------- shear stress ----------------------

  ! length parameter
  le = h
  rho = sqrt((r1+r2)/(2.0_dp*r2))
  re = (1.0_dp + rho - 1.0_dp/rho)*r2*cos(beta)
  omega = le/sqrt(re*t)
  rovert = re/t

  ! if (omega < 10):
  !     C_tau = sqrt(1.0 + 42.0/omega**3)
  ! elif (omega > 8.7*rovert):
  !     C_tau = 1.0/3.0*sqrt(omega/rovert)
  ! else:
  !     C_tau = 1.0
  call tausmooth(omega, rovert, C_tau)

  tau_zt_Rcr = 0.75_dp*E*C_tau*sqrt(1.0_dp/omega)/rovert

  ! reduction factor
  alpha_tau = 0.65_dp  ! high fabrifaction quality
  beta_tau = 0.6_dp
  lambda_tau0 = 0.4_dp
  eta_tau = 1.0_dp
  lambda_tau = sqrt(sigma_y/sqrt(3.0_dp)/tau_zt_Rcr)

  call buckling_reduction_factor(alpha_tau, beta_tau, eta_tau, lambda_tau0, lambda_tau, chi_tau)

  tau_zt_Rk = chi_tau*sigma_y/sqrt(3.0_dp)
  tau_zt_Rd = tau_zt_Rk/gamma_b

  ! buckling interaction parameters

  k_z = 1.25_dp + 0.75_dp*chi_z
  k_theta = 1.25_dp + 0.75_dp*chi_theta
  k_tau = 1.75_dp + 0.25_dp*chi_tau
  k_i = (chi_z*chi_theta)**2

  ! shell buckling utilization

  utilization = (sigma_z/sigma_z_Rd)**k_z + (sigma_t/sigma_t_Rd)**k_theta - &
  &k_i*(sigma_z*sigma_t/sigma_z_Rd/sigma_t_Rd) + (tau_zt/tau_zt_Rd)**k_tau

end subroutine shellBucklingOneSection


subroutine cxsmooth(omega, rovert, Cx)

  implicit none

  ! define precision to be the standard for a double precision ! on local system
  integer, parameter :: dp = kind(0.d0)

  !in
  real(dp), intent(in) :: omega, rovert

  !out
  real(dp), intent(out) :: Cx

  !local
  real(dp) :: Cxb, constant, ptL1, ptR1, ptL2, ptR2, ptL3, ptR3, fL, fR, gL, gR

  !evaluate the function
  Cxb = 6.0_dp  ! clamped-clamped
  constant = 1.0_dp + 1.83_dp/1.7_dp - 2.07_dp/1.7_dp**2

  ptL1 = 1.7_dp-0.25_dp
  ptR1 = 1.7_dp+0.25_dp

  ptL2 = 0.5_dp*rovert - 1.0_dp
  ptR2 = 0.5_dp*rovert + 1.0_dp

  ptL3 = (0.5_dp+Cxb)*rovert - 1.0_dp
  ptR3 = (0.5_dp+Cxb)*rovert + 1.0_dp


  if (omega < ptL1) then
    Cx = constant - 1.83_dp/omega + 2.07_dp/omega**2

  else if (omega >= ptL1 .and. omega <= ptR1) then

    fL = constant - 1.83_dp/ptL1 + 2.07_dp/ptL1**2
    fR = 1.0_dp
    gL = 1.83_dp/ptL1**2 - 4.14_dp/ptL1**3
    gR = 0.0_dp

    call cubic_spline_eval(ptL1, ptR1, fL, fR, gL, gR, omega, Cx)

  else if (omega > ptR1 .and. omega < ptL2) then
    Cx = 1.0_dp

  else if (omega >= ptL2 .and. omega <= ptR2) then

    fL = 1.0_dp
    fR = 1.0_dp + 0.2_dp/Cxb*(1.0_dp-2.0_dp*ptR2/rovert)
    gL = 0.0_dp
    gR = -0.4_dp/Cxb/rovert

    call cubic_spline_eval(ptL2, ptR2, fL, fR, gL, gR, omega, Cx)

  else if (omega > ptR2 .and. omega < ptL3) then
    Cx = 1.0_dp + 0.2_dp/Cxb*(1.0_dp-2.0_dp*omega/rovert)

  else if (omega >= ptL3 .and. omega <= ptR3) then

    fL = 1.0_dp + 0.2_dp/Cxb*(1.0_dp-2.0_dp*ptL3/rovert)
    fR = 0.6_dp
    gL = -0.4_dp/Cxb/rovert
    gR = 0.0_dp

    call cubic_spline_eval(ptL3, ptR3, fL, fR, gL, gR, omega, Cx)

  else
    Cx = 0.6_dp

  end if

end subroutine cxsmooth



subroutine sigmasmooth(omega, E, rovert, sigma)

  implicit none

  ! define precision to be the standard for a double precision ! on local system
  integer, parameter :: dp = kind(0.d0)

  !in
  real(dp), intent(in) :: omega, E, rovert

  !out
  real(dp), intent(out) :: sigma

  !local
  real(dp) :: Ctheta, ptL, ptR, offset, Cthetas, fL, fR, gL, gR, alpha1


  Ctheta = 1.5_dp  ! clamped-clamped

  ptL = 1.63_dp*rovert*Ctheta - 1.0_dp
  ptR = 1.63_dp*rovert*Ctheta + 1.0_dp

  if (omega < 20.0_dp*Ctheta) then
    offset = (10.0_dp/(20.0_dp*Ctheta)**2 - 5.0_dp/(20.0_dp*Ctheta)**3)
    Cthetas = 1.5_dp + 10.0_dp/omega**2 - 5.0_dp/omega**3 - offset
    sigma = 0.92_dp*E*Cthetas/omega/rovert

  else if (omega >= 20.0_dp*Ctheta .and. omega < ptL) then

    sigma = 0.92_dp*E*Ctheta/omega/rovert

  else if (omega >= ptL .and. omega <= ptR) then

    alpha1 = 0.92_dp/1.63_dp - 2.03_dp/1.63_dp**4

    fL = 0.92_dp*E*Ctheta/ptL/rovert
    fR = E*(1.0_dp/rovert)**2*(alpha1 + 2.03_dp*(Ctheta/ptR*rovert)**4)
    gL = -0.92_dp*E*Ctheta/rovert/ptL**2
    gR = -E*(1.0_dp/rovert)*2.03_dp*4.0_dp*(Ctheta/ptR*rovert)**3*Ctheta/ptR**2

    call cubic_spline_eval(ptL, ptR, fL, fR, gL, gR, omega, sigma)

  else

    alpha1 = 0.92_dp/1.63_dp - 2.03_dp/1.63_dp**4
    sigma = E*(1.0_dp/rovert)**2*(alpha1 + 2.03_dp*(Ctheta/omega*rovert)**4)

  end if

end subroutine sigmasmooth


subroutine tausmooth(omega, rovert, C_tau)

  implicit none

  ! define precision to be the standard for a double precision ! on local system
  integer, parameter :: dp = kind(0.d0)

  !in
  real(dp), intent(in) :: omega, rovert

  !out
  real(dp), intent(out) :: C_tau

  !local
  real(dp) :: ptL1, ptR1, ptL2, ptR2, fL, fR, gL, gR


  ptL1 = 9.0_dp
  ptR1 = 11.0_dp

  ptL2 = 8.7_dp*rovert - 1.0_dp
  ptR2 = 8.7_dp*rovert + 1.0_dp

  if (omega < ptL1) then
    C_tau = sqrt(1.0_dp + 42.0_dp/omega**3 - 42.0_dp/10.0_dp**3)

  else if (omega >= ptL1 .and. omega <= ptR1) then
    fL = sqrt(1.0_dp + 42.0_dp/ptL1**3 - 42.0_dp/10.0_dp**3)
    fR = 1.0_dp
    gL = -63.0_dp/ptL1**4/fL
    gR = 0.0_dp

    call cubic_spline_eval(ptL1, ptR1, fL, fR, gL, gR, omega, C_tau)

  else if (omega > ptR1 .and. omega < ptL2) then
    C_tau = 1.0_dp

  else if (omega >= ptL2 .and. omega <= ptR2) then
    fL = 1.0_dp
    fR = 1.0_dp/3.0_dp*sqrt(ptR2/rovert) + 1.0_dp - sqrt(8.7_dp)/3.0_dp
    gL = 0.0_dp
    gR = 1.0_dp/6.0_dp/sqrt(ptR2*rovert)

    call cubic_spline_eval(ptL2, ptR2, fL, fR, gL, gR, omega, C_tau)

  else
    C_tau = 1.0_dp/3.0_dp*sqrt(omega/rovert) + 1.0_dp - sqrt(8.7_dp)/3.0_dp

  end if

end subroutine tausmooth


subroutine buckling_reduction_factor(alpha, beta, eta, lambda_0, lambda_bar, chi)
  ! Computes a buckling reduction factor used in Eurocode shell buckling formula.

  implicit none

  ! define precision to be the standard for a double precision ! on local system
  integer, parameter :: dp = kind(0.d0)

  !in
  real(dp), intent(in) :: alpha, beta, eta, lambda_0, lambda_bar

  !out
  real(dp), intent(out) :: chi

  !local
  real(dp) :: lambda_p, ptL, ptR, fracR, fL, fR, gL, gR

  lambda_p = sqrt(alpha/(1.0_dp-beta))

  ptL = 0.9_dp*lambda_0
  ptR = 1.1_dp*lambda_0

  if (lambda_bar < ptL) then
    chi = 1.0_dp

  else if (lambda_bar >= ptL .and. lambda_bar <= ptR) then  ! cubic spline section

    fracR = (ptR-lambda_0)/(lambda_p-lambda_0)
    fL = 1.0_dp
    fR = 1.0_dp-beta*fracR**eta
    gL = 0.0_dp
    gR = -beta*eta*fracR**(eta-1)/(lambda_p-lambda_0)

    call cubic_spline_eval(ptL, ptR, fL, fR, gL, gR, lambda_bar, chi)

  else if (lambda_bar > ptR .and. lambda_bar < lambda_p) then
    chi = 1.0_dp - beta*((lambda_bar-lambda_0)/(lambda_p-lambda_0))**eta

  else
    chi = alpha/lambda_bar**2

  end if

  ! if (lambda_bar <= lambda_0):
  !     chi = 1.0
  ! elif (lambda_bar >= lambda_p):
  !     chi = alpha/lambda_bar**2
  ! else:
  !     chi = 1.0 - beta*((lambda_bar-lambda_0)/(lambda_p-lambda_0))**eta

end subroutine buckling_reduction_factor




subroutine cubic_spline_eval(x1, x2, f1, f2, g1, g2, x, poly)

  implicit none

  ! define precision to be the standard for a double precision ! on local system
  integer, parameter :: dp = kind(0.d0)

  !in
  real(dp), intent(in) :: x1, x2, f1, f2, g1, g2, x

  !out
  real(dp), intent(out) :: poly

  !local
  real(dp) :: det, B11, B12, B13, B14, B21, B22, B23, B24,&
              &B31, B32, B33, B34, B41, B42, B43, B44
  real(dp), dimension(4) :: A1, A2, A3, A4, b, coeff

  A1(1) = x1**3
  A1(2) = x1**2
  A1(3) = x1
  A1(4) = 1.0_dp

  A2(1) = x2**3
  A2(2) = x2**2
  A2(3) = x2
  A2(4) = 1.0_dp

  A3(1) = 3.0_dp*x1**2
  A3(2) = 2.0_dp*x1
  A3(3) = 1.0_dp
  A3(4) = 0.0_dp

  A4(1) = 3.0_dp*x2**2
  A4(2) = 2.0_dp*x2
  A4(3) = 1.0_dp
  A4(4) = 0.0_dp

  b(1) = f1
  b(2) = f2
  b(3) = g1
  b(4) = g2

  det = A1(1)*(A2(2)*A3(3)*A4(4)+A2(3)*A3(4)*A4(2)+A2(4)*A3(2)*A4(3)-A2(4)*A3(3)*&
        &A4(2)-A2(3)*A3(2)*A4(4)-A2(2)*A3(4)*A4(3))&
        & - A1(2)*(A2(1)*A3(3)*A4(4)+A2(3)*A3(4)*A4(1)+A2(4)*A3(1)*A4(3)-A2(4)*A3(3)*&
        &A4(1)-A2(3)*A3(1)*A4(4)-A2(1)*A3(4)*A4(3))&
        & + A1(3)*(A2(1)*A3(2)*A4(4)+A2(2)*A3(4)*A4(1)+A2(4)*A3(1)*A4(2)-A2(4)*A3(2)*&
        &A4(1)-A2(2)*A3(1)*A4(4)-A2(1)*A3(4)*A4(2))&
        & - A1(4)*(A2(1)*A3(2)*A4(3)+A2(2)*A3(3)*A4(1)+A2(3)*A3(1)*A4(2)-A2(3)*A3(2)*&
        &A4(1)-A2(2)*A3(1)*A4(3)-A2(1)*A3(3)*A4(2))

  ! entries of cof(A)
  B11 = A2(2)*A3(3)*A4(4)+A2(3)*A3(4)*A4(2)+A2(4)*A3(2)*A4(3)-A2(2)*A3(4)*A4(3)-A2(3)*A3(2)*A4(4)-A2(4)*A3(3)*A4(2)
  B12 = A1(2)*A3(4)*A4(3)+A1(3)*A3(2)*A4(4)+A1(4)*A3(3)*A4(2)-A1(2)*A3(3)*A4(4)-A1(3)*A3(4)*A4(2)-A1(4)*A3(2)*A4(3)
  B13 = A1(2)*A2(3)*A4(4)+A1(3)*A2(4)*A4(2)+A1(4)*A2(2)*A4(3)-A1(2)*A2(4)*A4(3)-A1(3)*A2(2)*A4(4)-A1(4)*A2(3)*A4(2)
  B14 = A1(2)*A2(4)*A3(3)+A1(3)*A2(2)*A3(4)+A1(4)*A2(3)*A3(2)-A1(2)*A2(3)*A3(4)-A1(3)*A2(4)*A3(2)-A1(4)*A2(2)*A3(3)

  B21 = A2(1)*A3(4)*A4(3)+A2(3)*A3(1)*A4(4)+A2(4)*A3(3)*A4(1)-A2(1)*A3(3)*A4(4)-A2(3)*A3(4)*A4(1)-A2(4)*A3(1)*A4(3)
  B22 = A1(1)*A3(3)*A4(4)+A1(3)*A3(4)*A4(1)+A1(4)*A3(1)*A4(3)-A1(1)*A3(4)*A4(3)-A1(3)*A3(1)*A4(4)-A1(4)*A3(3)*A4(1)
  B23 = A1(1)*A2(4)*A4(3)+A1(3)*A2(1)*A4(4)+A1(4)*A2(3)*A4(1)-A1(1)*A2(3)*A4(4)-A1(3)*A2(4)*A4(1)-A1(4)*A2(1)*A4(3)
  B24 = A1(1)*A2(3)*A3(4)+A1(3)*A2(4)*A3(1)+A1(4)*A2(1)*A3(3)-A1(1)*A2(4)*A3(3)-A1(3)*A2(1)*A3(4)-A1(4)*A2(3)*A3(1)

  B31 = A2(1)*A3(2)*A4(4)+A2(2)*A3(4)*A4(1)+A2(4)*A3(1)*A4(2)-A2(1)*A3(4)*A4(2)-A2(2)*A3(1)*A4(4)-A2(4)*A3(2)*A4(1)
  B32 = A1(1)*A3(4)*A4(2)+A1(2)*A3(1)*A4(4)+A1(4)*A3(2)*A4(1)-A1(1)*A3(2)*A4(4)-A1(2)*A3(4)*A4(1)-A1(4)*A3(1)*A4(2)
  B33 = A1(1)*A2(2)*A4(4)+A1(2)*A2(4)*A4(1)+A1(4)*A2(1)*A4(2)-A1(1)*A2(4)*A4(2)-A1(2)*A2(1)*A4(4)-A1(4)*A2(2)*A4(1)
  B34 = A1(1)*A2(4)*A3(2)+A1(2)*A2(1)*A3(4)+A1(4)*A2(2)*A3(1)-A1(1)*A2(2)*A3(4)-A1(2)*A2(4)*A3(1)-A1(4)*A2(1)*A3(2)

  B41 = A2(1)*A3(3)*A4(2)+A2(2)*A3(1)*A4(3)+A2(3)*A3(2)*A4(1)-A2(1)*A3(2)*A4(3)-A2(2)*A3(3)*A4(1)-A2(3)*A3(1)*A4(2)
  B42 = A1(1)*A3(2)*A4(3)+A1(2)*A3(3)*A4(1)+A1(3)*A3(1)*A4(2)-A1(1)*A3(3)*A4(2)-A1(2)*A3(1)*A4(3)-A1(3)*A3(2)*A4(1)
  B43 = A1(1)*A2(3)*A4(2)+A1(2)*A2(1)*A4(3)+A1(3)*A2(2)*A4(1)-A1(1)*A2(2)*A4(3)-A1(2)*A2(3)*A4(1)-A1(3)*A2(1)*A4(2)
  B44 = A1(1)*A2(2)*A3(3)+A1(2)*A2(3)*A3(1)+A1(3)*A2(1)*A3(2)-A1(1)*A2(3)*A3(2)-A1(2)*A2(1)*A3(3)-A1(3)*A2(2)*A3(1)

  !solve the equation Ax = b; x = A^(-1)b
  coeff(1) = (B11/det)*b(1)+(B12/det)*b(2)+(B13/det)*b(3)+(B14/det)*b(4)
  coeff(2) = (B21/det)*b(1)+(B22/det)*b(2)+(B23/det)*b(3)+(B24/det)*b(4)
  coeff(3) = (B31/det)*b(1)+(B32/det)*b(2)+(B33/det)*b(3)+(B34/det)*b(4)
  coeff(4) = (B41/det)*b(1)+(B42/det)*b(2)+(B43/det)*b(3)+(B44/det)*b(4)

  poly = coeff(1)*x**3+coeff(2)*x**2+coeff(3)*x+coeff(4)

end subroutine cubic_spline_eval








!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of shellbucklingeurocode in forward (tangent) mode:
!   variations   of useful results: eu_utilization
!   with respect to varying inputs: d t sigma_t sigma_z tau_zt
!   RW status of diff variables: d:in t:in eu_utilization:out sigma_t:in
!                sigma_z:in tau_zt:in
SUBROUTINE SHELLBUCKLINGEUROCODE_DV(npoints, d, dd, t, td, sigma_z, &
& sigma_zd, sigma_t, sigma_td, tau_zt, tau_ztd, l_reinforced, e, sigma_y&
& , gamma_f, gamma_b, eu_utilization, eu_utilizationd, nbdirs)
  ! USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: npoints
  REAL(dp), INTENT(IN) :: gamma_f, gamma_b
  REAL(dp), DIMENSION(npoints), INTENT(IN) :: d, t, sigma_z, sigma_t, &
& l_reinforced, e, sigma_y
  REAL(dp), DIMENSION(nbdirs, npoints), INTENT(IN) :: dd, td, &
& sigma_zd, sigma_td
  REAL(dp), DIMENSION(npoints), INTENT(IN) :: tau_zt
  REAL(dp), DIMENSION(nbdirs, npoints), INTENT(IN) :: tau_ztd
!out
  REAL(dp), DIMENSION(npoints), INTENT(OUT) :: eu_utilization
  REAL(dp), DIMENSION(nbdirs, npoints), INTENT(OUT) :: &
& eu_utilizationd
!local
  INTEGER :: i
  REAL(dp) :: h, r1, r2, t1, t2, sigma_z_shell, sigma_t_shell, &
& tau_zt_shell, utilization
  REAL(dp), DIMENSION(nbdirs) :: r1d, r2d, t1d, t2d, sigma_z_shelld, &
& sigma_t_shelld, tau_zt_shelld, utilizationd
  REAL(dp), DIMENSION(npoints) :: sigma_z_sh, sigma_t_sh, tau_zt_sh
  INTRINSIC KIND
  INTRINSIC ABS
  REAL(dp) :: abs0
  REAL(dp), DIMENSION(nbdirs) :: abs0d
  REAL(dp) :: abs1
  REAL(dp), DIMENSION(nbdirs) :: abs1d
  REAL(dp) :: abs2
  REAL(dp), DIMENSION(nbdirs) :: abs2d
  INTEGER :: nd
  INTEGER :: nbdirs
  DO nd=1,nbdirs
    eu_utilizationd(nd, :) = 0.0
  END DO
  DO i=1,npoints
    h = l_reinforced(i)
    DO nd=1,nbdirs
      r1d(nd) = dd(nd, i)/2.0_dp - td(nd, i)/2.0_dp
      r2d(nd) = dd(nd, i)/2.0_dp - td(nd, i)/2.0_dp
      t1d(nd) = td(nd, i)
      t2d(nd) = td(nd, i)
      sigma_z_shelld(nd) = sigma_zd(nd, i)
      sigma_t_shelld(nd) = sigma_td(nd, i)
      tau_zt_shelld(nd) = tau_ztd(nd, i)
    END DO
    r1 = d(i)/2.0_dp - t(i)/2.0_dp
    r2 = d(i)/2.0_dp - t(i)/2.0_dp
    t1 = t(i)
    t2 = t(i)
    sigma_z_shell = sigma_z(i)
    sigma_t_shell = sigma_t(i)
    tau_zt_shell = tau_zt(i)
    IF (sigma_z_shell .GE. 0.) THEN
      DO nd=1,nbdirs
        abs0d(nd) = sigma_z_shelld(nd)
      END DO
      abs0 = sigma_z_shell
    ELSE
      DO nd=1,nbdirs
        abs0d(nd) = -sigma_z_shelld(nd)
      END DO
      abs0 = -sigma_z_shell
    END IF
    DO nd=1,nbdirs
! TO DO: the following is non-smooth, although in general its probably OK
! change to magnitudes and add safety factor
      sigma_z_shelld(nd) = gamma_f*abs0d(nd)
    END DO
    sigma_z_shell = gamma_f*abs0
    IF (sigma_t_shell .GE. 0.) THEN
      DO nd=1,nbdirs
        abs1d(nd) = sigma_t_shelld(nd)
      END DO
      abs1 = sigma_t_shell
    ELSE
      DO nd=1,nbdirs
        abs1d(nd) = -sigma_t_shelld(nd)
      END DO
      abs1 = -sigma_t_shell
    END IF
    DO nd=1,nbdirs
      sigma_t_shelld(nd) = gamma_f*abs1d(nd)
    END DO
    sigma_t_shell = gamma_f*abs1
    IF (tau_zt_shell .GE. 0.) THEN
      DO nd=1,nbdirs
        abs2d(nd) = tau_zt_shelld(nd)
      END DO
      abs2 = tau_zt_shell
    ELSE
      DO nd=1,nbdirs
        abs2d(nd) = -tau_zt_shelld(nd)
      END DO
      abs2 = -tau_zt_shell
    END IF
    DO nd=1,nbdirs
      tau_zt_shelld(nd) = gamma_f*abs2d(nd)
    END DO
    tau_zt_shell = gamma_f*abs2
    CALL SHELLBUCKLINGONESECTION_DV(h, r1, r1d, r2, r2d, t1, t1d, t2, &
&                             t2d, gamma_b, sigma_z_shell, &
&                             sigma_z_shelld, sigma_t_shell, &
&                             sigma_t_shelld, tau_zt_shell, &
&                             tau_zt_shelld, e(i), sigma_y(i), &
&                             utilization, utilizationd, nbdirs)
    DO nd=1,nbdirs
      eu_utilizationd(nd, i) = utilizationd(nd)
    END DO
    eu_utilization(i) = utilization
! make them into vectors
! TODO is this necessary?
    sigma_z_sh(i) = sigma_z_shell
    sigma_t_sh(i) = sigma_t_shell
    tau_zt_sh(i) = tau_zt_shell
  END DO
END SUBROUTINE SHELLBUCKLINGEUROCODE_DV








!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of shellbucklingonesection in forward (tangent) mode:
!   variations   of useful results: utilization
!   with respect to varying inputs: sigma_t sigma_z t1 t2 r1 r2
!                tau_zt
SUBROUTINE SHELLBUCKLINGONESECTION_DV(h, r1, r1d, r2, r2d, t1, t1d, t2, &
& t2d, gamma_b, sigma_z, sigma_zd, sigma_t, sigma_td, tau_zt, tau_ztd, e&
& , sigma_y, utilization, utilizationd, nbdirs)
  ! USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
!in
  REAL(dp), INTENT(IN) :: h, r1, r2, t1, t2, gamma_b, sigma_z, sigma_t, &
& tau_zt, e, sigma_y
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: r1d, r2d, t1d, t2d, &
& sigma_zd, sigma_td, tau_ztd
!out
  REAL(dp), INTENT(OUT) :: utilization
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: utilizationd
!local
  REAL(dp) :: beta, l, t, le, re, omega, rovert, cx, sigma_z_rcr, &
& lambda_z0, beta_z, eta_z, q, lambda_z, delta_wk, alpha_z, chi_z, &
& sigma_z_rk, sigma_z_rd, sigma_t_rcr, alpha_t, lambda_t0, beta_t, eta_t&
& , lambda_t, chi_theta, sigma_t_rk, sigma_t_rd, rho, c_tau, tau_zt_rcr&
& , alpha_tau, beta_tau, lambda_tau0, eta_tau, lambda_tau, chi_tau, &
& tau_zt_rk, tau_zt_rd, k_z, k_theta, k_tau, k_i
  REAL(dp), DIMENSION(nbdirs) :: betad, ld, td, led, red, omegad, &
& rovertd, cxd, sigma_z_rcrd, lambda_zd, delta_wkd, alpha_zd, chi_zd, &
& sigma_z_rkd, sigma_z_rdd, sigma_t_rcrd, alpha_td, lambda_td, &
& chi_thetad, sigma_t_rkd, sigma_t_rdd, rhod, c_taud, tau_zt_rcrd, &
& alpha_taud, lambda_taud, chi_taud, tau_zt_rkd, tau_zt_rdd, k_zd, &
& k_thetad, k_taud, k_id
  INTRINSIC KIND
  INTRINSIC ATAN2
  INTRINSIC COS
  INTRINSIC SQRT
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  REAL(dp) :: pwx1
  REAL(dp), DIMENSION(nbdirs) :: pwx1d
  REAL(dp) :: pwr1
  REAL(dp), DIMENSION(nbdirs) :: pwr1d
  REAL(dp) :: pwx2
  REAL(dp), DIMENSION(nbdirs) :: pwx2d
  REAL(dp) :: pwr2
  REAL(dp), DIMENSION(nbdirs) :: pwr2d
  REAL(dp) :: pwx3
  REAL(dp), DIMENSION(nbdirs) :: pwx3d
  REAL(dp) :: pwr3
  REAL(dp), DIMENSION(nbdirs) :: pwr3d
  INTEGER :: nd
  INTEGER :: nbdirs
  beta = ATAN2(r1 - r2, h)
  l = h/COS(beta)
  t = 0.5_dp*(t1+t2)
  le = l
  re = 0.5_dp*(r1+r2)/COS(beta)
  arg1 = re*t
  result1 = SQRT(arg1)
  rovert = re/t
  DO nd=1,nbdirs
! ----- geometric parameters --------
    betad(nd) = (r1d(nd)-r2d(nd))*h/((r1-r2)**2+h**2)
    ld(nd) = -((-(h*betad(nd)*SIN(beta)))/COS(beta)**2)
    td(nd) = 0.5_dp*(t1d(nd)+t2d(nd))
! ------------- axial stress -------------
! length parameter
    led(nd) = ld(nd)
    red(nd) = (0.5_dp*(r1d(nd)+r2d(nd))*COS(beta)+0.5_dp*(r1+r2)*betad(&
&     nd)*SIN(beta))/COS(beta)**2
    arg1d(nd) = red(nd)*t + re*td(nd)
    IF (arg1 .EQ. 0.0) THEN
      result1d(nd) = 0.0
    ELSE
      result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    omegad(nd) = (led(nd)*result1-le*result1d(nd))/result1**2
    rovertd(nd) = (red(nd)*t-re*td(nd))/t**2
    IF (rovert .EQ. 0.0) THEN
      result1d(nd) = 0.0
    ELSE
      result1d(nd) = rovertd(nd)/(2.0*SQRT(rovert))
    END IF
! ---------------- hoop stress ------------------
! length parameter
    led(nd) = ld(nd)
    red(nd) = (0.5_dp*(r1d(nd)+r2d(nd))*COS(beta)+0.5_dp*(r1+r2)*betad(&
&     nd)*SIN(beta))/COS(beta)**2
    alpha_td(nd) = 0.0
    alpha_taud(nd) = 0.0
  END DO
  omega = le/result1
! compute Cx
  CALL CXSMOOTH_DV(omega, omegad, rovert, rovertd, cx, cxd, nbdirs)
  sigma_z_rcr = 0.605_dp*e*cx/rovert
! quality parameter - high
  q = 25.0_dp
  arg1 = sigma_y/sigma_z_rcr
  result1 = SQRT(rovert)
  delta_wk = 1.0_dp/q*result1*t
  re = 0.5_dp*(r1+r2)/COS(beta)
  DO nd=1,nbdirs
! if omega <= 1.7:
!     Cx = 1.36 - 1.83/omega + 2.07/omega/omega
! elif omega > 0.5*rovert:
!     Cxb = 6.0  ! clamped-clamped
!     Cx = max(0.6, 1 + 0.2/Cxb*(1-2.0*omega/rovert))
! else:
!     Cx = 1.0
! critical axial buckling stress
    sigma_z_rcrd(nd) = (0.605_dp*e*cxd(nd)*rovert-0.605_dp*e*cx*rovertd(&
&     nd))/rovert**2
    arg1d(nd) = -(sigma_y*sigma_z_rcrd(nd)/sigma_z_rcr**2)
    IF (arg1 .EQ. 0.0) THEN
      lambda_zd(nd) = 0.0
    ELSE
      lambda_zd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    delta_wkd(nd) = (result1d(nd)*t+result1*td(nd))/q
    alpha_zd(nd) = -(0.62_dp*1.91_dp*1.44*(delta_wk/t)**0.44*(delta_wkd(&
&     nd)*t-delta_wk*td(nd))/t**2/(1.0_dp+1.91_dp*(delta_wk/t)**1.44)**2&
&     )
    arg1d(nd) = red(nd)*t + re*td(nd)
    rovertd(nd) = (red(nd)*t-re*td(nd))/t**2
  END DO
! compute buckling reduction factors
  lambda_z0 = 0.2_dp
  beta_z = 0.6_dp
  eta_z = 1.0_dp
  lambda_z = SQRT(arg1)
  alpha_z = 0.62_dp/(1.0_dp+1.91_dp*(delta_wk/t)**1.44)
  CALL BUCKLING_REDUCTION_FACTOR_DV(alpha_z, alpha_zd, beta_z, eta_z, &
&                             lambda_z0, lambda_z, lambda_zd, chi_z, &
&                             chi_zd, nbdirs)
  sigma_z_rk = chi_z*sigma_y
  sigma_z_rd = sigma_z_rk/gamma_b
  le = l
  arg1 = re*t
  result1 = SQRT(arg1)
  k_z = 1.25_dp + 0.75_dp*chi_z
  pwx1 = sigma_z/sigma_z_rd
  DO nd=1,nbdirs
! design buckling stress
    sigma_z_rkd(nd) = sigma_y*chi_zd(nd)
    sigma_z_rdd(nd) = sigma_z_rkd(nd)/gamma_b
    IF (arg1 .EQ. 0.0) THEN
      result1d(nd) = 0.0
    ELSE
      result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    omegad(nd) = (led(nd)*result1-le*result1d(nd))/result1**2
! buckling interaction parameters
    k_zd(nd) = 0.75_dp*chi_zd(nd)
! shell buckling utilization
    pwx1d(nd) = (sigma_zd(nd)*sigma_z_rd-sigma_z*sigma_z_rdd(nd))/&
&     sigma_z_rd**2
    IF (pwx1 .GT. 0.0) THEN
      pwr1d(nd) = pwx1**k_z*(LOG(pwx1)*k_zd(nd)+k_z*pwx1d(nd)/pwx1)
    ELSE IF (pwx1 .EQ. 0.0) THEN
      IF (k_z .EQ. 1.0) THEN
        pwr1d(nd) = pwx1d(nd)
      ELSE
        pwr1d(nd) = 0.0
      END IF
    ELSE IF (k_z .EQ. INT(k_z)) THEN
      pwr1d(nd) = k_z*pwx1**(k_z-1)*pwx1d(nd)
    ELSE
      pwr1d(nd) = 0.0
    END IF
  END DO
  omega = le/result1
  rovert = re/t
! Ctheta = 1.5  ! clamped-clamped
! CthetaS = 1.5 + 10.0/omega**2 - 5.0/omega**3
! ! critical hoop buckling stress
! if (omega/Ctheta < 20.0):
!     sigma_t_Rcr = 0.92*E*CthetaS/omega/rovert
! elif (omega/Ctheta > 1.63*rovert):
!     sigma_t_Rcr = E*(1.0/rovert)**2*(0.275 + 2.03*(Ctheta/omega*rovert)**4)
! else:
!     sigma_t_Rcr = 0.92*E*Ctheta/omega/rovert
  CALL SIGMASMOOTH_DV(omega, omegad, e, rovert, rovertd, sigma_t_rcr, &
&               sigma_t_rcrd, nbdirs)
! buckling reduction factor
! high fabrication quality
  alpha_t = 0.65_dp
  lambda_t0 = 0.4_dp
  beta_t = 0.6_dp
  eta_t = 1.0_dp
  arg1 = sigma_y/sigma_t_rcr
  DO nd=1,nbdirs
    arg1d(nd) = -(sigma_y*sigma_t_rcrd(nd)/sigma_t_rcr**2)
    IF (arg1 .EQ. 0.0) THEN
      lambda_td(nd) = 0.0
    ELSE
      lambda_td(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    arg1d(nd) = ((r1d(nd)+r2d(nd))*2.0_dp*r2-(r1+r2)*2.0_dp*r2d(nd))/(&
&     2.0_dp*r2)**2
  END DO
  lambda_t = SQRT(arg1)
  CALL BUCKLING_REDUCTION_FACTOR_DV(alpha_t, alpha_td, beta_t, eta_t, &
&                             lambda_t0, lambda_t, lambda_td, chi_theta&
&                             , chi_thetad, nbdirs)
  sigma_t_rk = chi_theta*sigma_y
  sigma_t_rd = sigma_t_rk/gamma_b
  arg1 = (r1+r2)/(2.0_dp*r2)
  rho = SQRT(arg1)
  re = (1.0_dp+rho-1.0_dp/rho)*r2*COS(beta)
  k_theta = 1.25_dp + 0.75_dp*chi_theta
  pwx2 = sigma_t/sigma_t_rd
  DO nd=1,nbdirs
    sigma_t_rkd(nd) = sigma_y*chi_thetad(nd)
    sigma_t_rdd(nd) = sigma_t_rkd(nd)/gamma_b
    IF (arg1 .EQ. 0.0) THEN
      rhod(nd) = 0.0
    ELSE
      rhod(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    red(nd) = ((rhod(nd)+rhod(nd)/rho**2)*r2+(1.0_dp+rho-1.0_dp/rho)*r2d&
&     (nd))*COS(beta) - (1.0_dp+rho-1.0_dp/rho)*r2*betad(nd)*SIN(beta)
    arg1d(nd) = red(nd)*t + re*td(nd)
    rovertd(nd) = (red(nd)*t-re*td(nd))/t**2
    k_thetad(nd) = 0.75_dp*chi_thetad(nd)
    k_id(nd) = 2*chi_z*chi_theta*(chi_zd(nd)*chi_theta+chi_z*chi_thetad(&
&     nd))
    pwx2d(nd) = (sigma_td(nd)*sigma_t_rd-sigma_t*sigma_t_rdd(nd))/&
&     sigma_t_rd**2
    IF (pwx2 .GT. 0.0) THEN
      pwr2d(nd) = pwx2**k_theta*(LOG(pwx2)*k_thetad(nd)+k_theta*pwx2d(nd&
&       )/pwx2)
    ELSE IF (pwx2 .EQ. 0.0) THEN
      IF (k_theta .EQ. 1.0) THEN
        pwr2d(nd) = pwx2d(nd)
      ELSE
        pwr2d(nd) = 0.0
      END IF
    ELSE IF (k_theta .EQ. INT(k_theta)) THEN
      pwr2d(nd) = k_theta*pwx2**(k_theta-1)*pwx2d(nd)
    ELSE
      pwr2d(nd) = 0.0
    END IF
  END DO
! ----------------- shear stress ----------------------
! length parameter
  le = h
  arg1 = re*t
  result1 = SQRT(arg1)
  omega = le/result1
  DO nd=1,nbdirs
    IF (arg1 .EQ. 0.0) THEN
      result1d(nd) = 0.0
    ELSE
      result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    omegad(nd) = -(le*result1d(nd)/result1**2)
    arg1d(nd) = -(omegad(nd)/omega**2)
  END DO
  rovert = re/t
! if (omega < 10):
!     C_tau = sqrt(1.0 + 42.0/omega**3)
! elif (omega > 8.7*rovert):
!     C_tau = 1.0/3.0*sqrt(omega/rovert)
! else:
!     C_tau = 1.0
  CALL TAUSMOOTH_DV(omega, omegad, rovert, rovertd, c_tau, c_taud, &
&             nbdirs)
  arg1 = 1.0_dp/omega
  result1 = SQRT(arg1)
  DO nd=1,nbdirs
    IF (arg1 .EQ. 0.0) THEN
      result1d(nd) = 0.0
    ELSE
      result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    tau_zt_rcrd(nd) = (0.75_dp*e*(c_taud(nd)*result1+c_tau*result1d(nd))&
&     *rovert-0.75_dp*e*c_tau*result1*rovertd(nd))/rovert**2
  END DO
  tau_zt_rcr = 0.75_dp*e*c_tau*result1/rovert
! reduction factor
! high fabrifaction quality
  alpha_tau = 0.65_dp
  beta_tau = 0.6_dp
  lambda_tau0 = 0.4_dp
  eta_tau = 1.0_dp
  result1 = SQRT(3.0_dp)
  arg1 = sigma_y/result1/tau_zt_rcr
  DO nd=1,nbdirs
    arg1d(nd) = -(sigma_y*tau_zt_rcrd(nd)/result1/tau_zt_rcr**2)
    IF (arg1 .EQ. 0.0) THEN
      lambda_taud(nd) = 0.0
    ELSE
      lambda_taud(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
  END DO
  lambda_tau = SQRT(arg1)
  CALL BUCKLING_REDUCTION_FACTOR_DV(alpha_tau, alpha_taud, beta_tau, &
&                             eta_tau, lambda_tau0, lambda_tau, &
&                             lambda_taud, chi_tau, chi_taud, nbdirs)
  result1 = SQRT(3.0_dp)
  tau_zt_rk = chi_tau*sigma_y/result1
  tau_zt_rd = tau_zt_rk/gamma_b
  k_tau = 1.75_dp + 0.25_dp*chi_tau
  k_i = (chi_z*chi_theta)**2
  pwx3 = tau_zt/tau_zt_rd
  DO nd=1,nbdirs
    tau_zt_rkd(nd) = sigma_y*chi_taud(nd)/result1
    tau_zt_rdd(nd) = tau_zt_rkd(nd)/gamma_b
    k_taud(nd) = 0.25_dp*chi_taud(nd)
    pwx3d(nd) = (tau_ztd(nd)*tau_zt_rd-tau_zt*tau_zt_rdd(nd))/tau_zt_rd&
&     **2
    IF (pwx3 .GT. 0.0) THEN
      pwr3d(nd) = pwx3**k_tau*(LOG(pwx3)*k_taud(nd)+k_tau*pwx3d(nd)/pwx3&
&       )
    ELSE IF (pwx3 .EQ. 0.0) THEN
      IF (k_tau .EQ. 1.0) THEN
        pwr3d(nd) = pwx3d(nd)
      ELSE
        pwr3d(nd) = 0.0
      END IF
    ELSE IF (k_tau .EQ. INT(k_tau)) THEN
      pwr3d(nd) = k_tau*pwx3**(k_tau-1)*pwx3d(nd)
    ELSE
      pwr3d(nd) = 0.0
    END IF
    utilizationd(nd) = pwr1d(nd) + pwr2d(nd) - k_id(nd)*sigma_z*sigma_t/&
&     (sigma_z_rd*sigma_t_rd) - k_i*(((sigma_zd(nd)*sigma_t+sigma_z*&
&     sigma_td(nd))*sigma_z_rd-sigma_z*sigma_t*sigma_z_rdd(nd))*&
&     sigma_t_rd/sigma_z_rd**2-sigma_z*sigma_t*sigma_t_rdd(nd)/&
&     sigma_z_rd)/sigma_t_rd**2 + pwr3d(nd)
  END DO
  pwr1 = pwx1**k_z
  pwr2 = pwx2**k_theta
  pwr3 = pwx3**k_tau
  utilization = pwr1 + pwr2 - k_i*(sigma_z*sigma_t/sigma_z_rd/sigma_t_rd&
&   ) + pwr3
END SUBROUTINE SHELLBUCKLINGONESECTION_DV









!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of tausmooth in forward (tangent) mode:
!   variations   of useful results: c_tau
!   with respect to varying inputs: omega rovert
SUBROUTINE TAUSMOOTH_DV(omega, omegad, rovert, rovertd, c_tau, c_taud, &
& nbdirs)
  ! USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
!in
  REAL(dp), INTENT(IN) :: omega, rovert
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: omegad, rovertd
!out
  REAL(dp), INTENT(OUT) :: c_tau
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: c_taud
!local
  REAL(dp) :: ptl1, ptr1, ptl2, ptr2, fl, fr, gl, gr
  REAL(dp), DIMENSION(nbdirs) :: ptl1d, ptr1d, ptl2d, ptr2d, fld, frd&
& , gld, grd
  INTRINSIC KIND
  INTRINSIC SQRT
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  REAL(dp) :: result2
  INTEGER :: nd
  INTEGER :: nbdirs
  ptl1 = 9.0_dp
  ptr1 = 11.0_dp
  DO nd=1,nbdirs
    ptl2d(nd) = 8.7_dp*rovertd(nd)
    ptr2d(nd) = 8.7_dp*rovertd(nd)
  END DO
  ptl2 = 8.7_dp*rovert - 1.0_dp
  ptr2 = 8.7_dp*rovert + 1.0_dp
  IF (omega .LT. ptl1) THEN
    arg1 = 1.0_dp + 42.0_dp/omega**3 - 42.0_dp/10.0_dp**3
    DO nd=1,nbdirs
      arg1d(nd) = -(42.0_dp*3*omega**2*omegad(nd)/(omega**3)**2)
      IF (arg1 .EQ. 0.0) THEN
        c_taud(nd) = 0.0
      ELSE
        c_taud(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
    END DO
    c_tau = SQRT(arg1)
  ELSE IF (omega .GE. ptl1 .AND. omega .LE. ptr1) THEN
    arg1 = 1.0_dp + 42.0_dp/ptl1**3 - 42.0_dp/10.0_dp**3
    fl = SQRT(arg1)
    fr = 1.0_dp
    gl = -(63.0_dp/ptl1**4/fl)
    gr = 0.0_dp
    DO nd=1,nbdirs
      grd(nd) = 0.0
      gld(nd) = 0.0
      frd(nd) = 0.0
      fld(nd) = 0.0
      ptr1d(nd) = 0.0
      ptl1d(nd) = 0.0
    END DO
    CALL CUBIC_SPLINE_EVAL_DV(ptl1, ptl1d, ptr1, ptr1d, fl, fld, fr, frd&
&                       , gl, gld, gr, grd, omega, omegad, c_tau, c_taud&
&                       , nbdirs)
  ELSE IF (omega .GT. ptr1 .AND. omega .LT. ptl2) THEN
    c_tau = 1.0_dp
    DO nd=1,nbdirs
      c_taud(nd) = 0.0
    END DO
  ELSE IF (omega .GE. ptl2 .AND. omega .LE. ptr2) THEN
    fl = 1.0_dp
    arg1 = ptr2/rovert
    DO nd=1,nbdirs
      arg1d(nd) = (ptr2d(nd)*rovert-ptr2*rovertd(nd))/rovert**2
      IF (arg1 .EQ. 0.0) THEN
        result1d(nd) = 0.0
      ELSE
        result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
      frd(nd) = result1d(nd)/3.0_dp
      arg1d(nd) = ptr2d(nd)*rovert + ptr2*rovertd(nd)
      gld(nd) = 0.0
      fld(nd) = 0.0
    END DO
    result1 = SQRT(arg1)
    result2 = SQRT(8.7_dp)
    fr = 1.0_dp/3.0_dp*result1 + 1.0_dp - result2/3.0_dp
    gl = 0.0_dp
    arg1 = ptr2*rovert
    result1 = SQRT(arg1)
    DO nd=1,nbdirs
      IF (arg1 .EQ. 0.0) THEN
        result1d(nd) = 0.0
      ELSE
        result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
      grd(nd) = -(result1d(nd)/6.0_dp/result1**2)
    END DO
    gr = 1.0_dp/6.0_dp/result1
    CALL CUBIC_SPLINE_EVAL_DV(ptl2, ptl2d, ptr2, ptr2d, fl, fld, fr, frd&
&                       , gl, gld, gr, grd, omega, omegad, c_tau, c_taud&
&                       , nbdirs)
  ELSE
    arg1 = omega/rovert
    DO nd=1,nbdirs
      arg1d(nd) = (omegad(nd)*rovert-omega*rovertd(nd))/rovert**2
      IF (arg1 .EQ. 0.0) THEN
        result1d(nd) = 0.0
      ELSE
        result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
      c_taud(nd) = result1d(nd)/3.0_dp
    END DO
    result1 = SQRT(arg1)
    result2 = SQRT(8.7_dp)
    c_tau = 1.0_dp/3.0_dp*result1 + 1.0_dp - result2/3.0_dp
  END IF
END SUBROUTINE TAUSMOOTH_DV









!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of sigmasmooth in forward (tangent) mode:
!   variations   of useful results: sigma
!   with respect to varying inputs: omega rovert
SUBROUTINE SIGMASMOOTH_DV(omega, omegad, e, rovert, rovertd, sigma, &
& sigmad, nbdirs)
  ! USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
!in
  REAL(dp), INTENT(IN) :: omega, e, rovert
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: omegad, rovertd
!out
  REAL(dp), INTENT(OUT) :: sigma
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: sigmad
!local
  REAL(dp) :: ctheta, ptl, ptr, offset, cthetas, fl, fr, gl, gr, alpha1
  REAL(dp), DIMENSION(nbdirs) :: ptld, ptrd, cthetasd, fld, frd, gld&
& , grd
  INTRINSIC KIND
  INTEGER :: nd
  INTEGER :: nbdirs
! clamped-clamped
  ctheta = 1.5_dp
  DO nd=1,nbdirs
    ptld(nd) = 1.63_dp*ctheta*rovertd(nd)
    ptrd(nd) = 1.63_dp*ctheta*rovertd(nd)
  END DO
  ptl = 1.63_dp*rovert*ctheta - 1.0_dp
  ptr = 1.63_dp*rovert*ctheta + 1.0_dp
  IF (omega .LT. 20.0_dp*ctheta) THEN
    offset = 10.0_dp/(20.0_dp*ctheta)**2 - 5.0_dp/(20.0_dp*ctheta)**3
    cthetas = 1.5_dp + 10.0_dp/omega**2 - 5.0_dp/omega**3 - offset
    DO nd=1,nbdirs
      cthetasd(nd) = 5.0_dp*3*omega**2*omegad(nd)/(omega**3)**2 - &
&       10.0_dp*2*omega*omegad(nd)/(omega**2)**2
      sigmad(nd) = ((0.92_dp*e*cthetasd(nd)*omega-0.92_dp*e*cthetas*&
&       omegad(nd))*rovert/omega**2-0.92_dp*e*cthetas*rovertd(nd)/omega)&
&       /rovert**2
    END DO
    sigma = 0.92_dp*e*cthetas/omega/rovert
  ELSE IF (omega .GE. 20.0_dp*ctheta .AND. omega .LT. ptl) THEN
    DO nd=1,nbdirs
      sigmad(nd) = (-(0.92_dp*e*ctheta*omegad(nd)*rovert/omega**2)-&
&       0.92_dp*e*ctheta*rovertd(nd)/omega)/rovert**2
    END DO
    sigma = 0.92_dp*e*ctheta/omega/rovert
  ELSE IF (omega .GE. ptl .AND. omega .LE. ptr) THEN
    alpha1 = 0.92_dp/1.63_dp - 2.03_dp/1.63_dp**4
    DO nd=1,nbdirs
      fld(nd) = (-(0.92_dp*e*ctheta*ptld(nd)*rovert/ptl**2)-0.92_dp*e*&
&       ctheta*rovertd(nd)/ptl)/rovert**2
      frd(nd) = e*(rovert*2.03_dp*4*ctheta**3*(ctheta*rovertd(nd)/ptr-&
&       ctheta*ptrd(nd)*rovert/ptr**2)/ptr**3-2*rovertd(nd)*(alpha1+&
&       2.03_dp*(ctheta/ptr*rovert)**4)/rovert**3)
      gld(nd) = -((-(0.92_dp*e*ctheta*rovertd(nd)*ptl**2/rovert**2)-&
&       0.92_dp*e*ctheta*2*ptl*ptld(nd)/rovert)/(ptl**2)**2)
      grd(nd) = -((e*2.03_dp*4.0_dp*ctheta*(rovert*3*ctheta**2*(ctheta*&
&       rovertd(nd)/ptr-ctheta*ptrd(nd)*rovert/ptr**2)/ptr**2-rovertd(nd&
&       )*rovert*ctheta**3/ptr**3)*ptr**2-e*rovert**2*2.03_dp*4.0_dp*&
&       ctheta**4*2*ptrd(nd)/ptr**2)/(ptr**2)**2)
    END DO
    fl = 0.92_dp*e*ctheta/ptl/rovert
    fr = e*(1.0_dp/rovert)**2*(alpha1+2.03_dp*(ctheta/ptr*rovert)**4)
    gl = -(0.92_dp*e*ctheta/rovert/ptl**2)
    gr = -(e*(1.0_dp/rovert)*2.03_dp*4.0_dp*(ctheta/ptr*rovert)**3*&
&     ctheta/ptr**2)
    CALL CUBIC_SPLINE_EVAL_DV(ptl, ptld, ptr, ptrd, fl, fld, fr, frd, gl&
&                       , gld, gr, grd, omega, omegad, sigma, sigmad, &
&                       nbdirs)
  ELSE
    alpha1 = 0.92_dp/1.63_dp - 2.03_dp/1.63_dp**4
    DO nd=1,nbdirs
      sigmad(nd) = e*(rovert*2.03_dp*4*ctheta**3*(ctheta*rovertd(nd)/&
&       omega-ctheta*omegad(nd)*rovert/omega**2)/omega**3-2*rovertd(nd)*&
&       (alpha1+2.03_dp*(ctheta/omega*rovert)**4)/rovert**3)
    END DO
    sigma = e*(1.0_dp/rovert)**2*(alpha1+2.03_dp*(ctheta/omega*rovert)**&
&     4)
  END IF
END SUBROUTINE SIGMASMOOTH_DV









!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of buckling_reduction_factor in forward (tangent) mode:
!   variations   of useful results: chi
!   with respect to varying inputs: alpha lambda_bar
SUBROUTINE BUCKLING_REDUCTION_FACTOR_DV(alpha, alphad, beta, eta, &
& lambda_0, lambda_bar, lambda_bard, chi, chid, nbdirs)
  ! USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! if (lambda_bar <= lambda_0):
!     chi = 1.0
! elif (lambda_bar >= lambda_p):
!     chi = alpha/lambda_bar**2
! else:
!     chi = 1.0 - beta*((lambda_bar-lambda_0)/(lambda_p-lambda_0))**eta
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
!in
  REAL(dp), INTENT(IN) :: alpha, beta, eta, lambda_0, lambda_bar
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: alphad, lambda_bard
!out
  REAL(dp), INTENT(OUT) :: chi
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: chid
!local
  REAL(dp) :: lambda_p, ptl, ptr, fracr, fl, fr, gl, gr
  REAL(dp), DIMENSION(nbdirs) :: lambda_pd, ptld, ptrd, fracrd, fld, &
& frd, gld, grd
  INTRINSIC KIND
  INTRINSIC SQRT
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: pwr1
  REAL(dp), DIMENSION(nbdirs) :: pwr1d
  REAL(dp) :: pwy1
  REAL(dp) :: pwx1
  REAL(dp), DIMENSION(nbdirs) :: pwx1d
  INTEGER :: nd
  INTEGER :: nbdirs
  arg1 = alpha/(1.0_dp-beta)
  DO nd=1,nbdirs
    arg1d(nd) = alphad(nd)/(1.0_dp-beta)
    IF (arg1 .EQ. 0.0) THEN
      lambda_pd(nd) = 0.0
    ELSE
      lambda_pd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
  END DO
  lambda_p = SQRT(arg1)
  ptl = 0.9_dp*lambda_0
  ptr = 1.1_dp*lambda_0
  IF (lambda_bar .LT. ptl) THEN
    chi = 1.0_dp
    DO nd=1,nbdirs
      chid(nd) = 0.0
    END DO
  ELSE IF (lambda_bar .GE. ptl .AND. lambda_bar .LE. ptr) THEN
    fracr = (ptr-lambda_0)/(lambda_p-lambda_0)
    pwr1 = fracr**eta
    fr = 1.0_dp - beta*pwr1
    pwy1 = eta - 1
    pwr1 = fracr**pwy1
    DO nd=1,nbdirs
! cubic spline section
      fracrd(nd) = -((ptr-lambda_0)*lambda_pd(nd)/(lambda_p-lambda_0)**2&
&       )
      IF (fracr .GT. 0.0 .OR. (fracr .LT. 0.0 .AND. eta .EQ. INT(eta))) &
&     THEN
        pwr1d(nd) = eta*fracr**(eta-1)*fracrd(nd)
      ELSE IF (fracr .EQ. 0.0 .AND. eta .EQ. 1.0) THEN
        pwr1d(nd) = fracrd(nd)
      ELSE
        pwr1d(nd) = 0.0
      END IF
      frd(nd) = -(beta*pwr1d(nd))
      IF (fracr .GT. 0.0 .OR. (fracr .LT. 0.0 .AND. pwy1 .EQ. INT(pwy1))&
&     ) THEN
        pwr1d(nd) = pwy1*fracr**(pwy1-1)*fracrd(nd)
      ELSE IF (fracr .EQ. 0.0 .AND. pwy1 .EQ. 1.0) THEN
        pwr1d(nd) = fracrd(nd)
      ELSE
        pwr1d(nd) = 0.0
      END IF
      grd(nd) = -((beta*eta*pwr1d(nd)*(lambda_p-lambda_0)-beta*eta*pwr1*&
&       lambda_pd(nd))/(lambda_p-lambda_0)**2)
      gld(nd) = 0.0
      fld(nd) = 0.0
      ptrd(nd) = 0.0
      ptld(nd) = 0.0
    END DO
    fl = 1.0_dp
    gl = 0.0_dp
    gr = -(beta*eta*pwr1/(lambda_p-lambda_0))
    CALL CUBIC_SPLINE_EVAL_DV(ptl, ptld, ptr, ptrd, fl, fld, fr, frd, gl&
&                       , gld, gr, grd, lambda_bar, lambda_bard, chi, &
&                       chid, nbdirs)
  ELSE IF (lambda_bar .GT. ptr .AND. lambda_bar .LT. lambda_p) THEN
    pwx1 = (lambda_bar-lambda_0)/(lambda_p-lambda_0)
    DO nd=1,nbdirs
      pwx1d(nd) = (lambda_bard(nd)*(lambda_p-lambda_0)-(lambda_bar-&
&       lambda_0)*lambda_pd(nd))/(lambda_p-lambda_0)**2
      IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. eta .EQ. INT(eta))) &
&     THEN
        pwr1d(nd) = eta*pwx1**(eta-1)*pwx1d(nd)
      ELSE IF (pwx1 .EQ. 0.0 .AND. eta .EQ. 1.0) THEN
        pwr1d(nd) = pwx1d(nd)
      ELSE
        pwr1d(nd) = 0.0
      END IF
      chid(nd) = -(beta*pwr1d(nd))
    END DO
    pwr1 = pwx1**eta
    chi = 1.0_dp - beta*pwr1
  ELSE
    DO nd=1,nbdirs
      chid(nd) = (alphad(nd)*lambda_bar**2-alpha*2*lambda_bar*&
&       lambda_bard(nd))/(lambda_bar**2)**2
    END DO
    chi = alpha/lambda_bar**2
  END IF
END SUBROUTINE BUCKLING_REDUCTION_FACTOR_DV









!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of cxsmooth in forward (tangent) mode:
!   variations   of useful results: cx
!   with respect to varying inputs: omega rovert
SUBROUTINE CXSMOOTH_DV(omega, omegad, rovert, rovertd, cx, cxd, nbdirs)
  ! USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
!in
  REAL(dp), INTENT(IN) :: omega, rovert
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: omegad, rovertd
!out
  REAL(dp), INTENT(OUT) :: cx
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: cxd
!local
  REAL(dp) :: cxb, constant, ptl1, ptr1, ptl2, ptr2, ptl3, ptr3, fl, fr&
& , gl, gr
  REAL(dp), DIMENSION(nbdirs) :: ptl1d, ptr1d, ptl2d, ptr2d, ptl3d, &
& ptr3d, fld, frd, gld, grd
  INTRINSIC KIND
  INTEGER :: nd
  INTEGER :: nbdirs
!evaluate the function
! clamped-clamped
  cxb = 6.0_dp
  constant = 1.0_dp + 1.83_dp/1.7_dp - 2.07_dp/1.7_dp**2
  ptl1 = 1.7_dp - 0.25_dp
  ptr1 = 1.7_dp + 0.25_dp
  DO nd=1,nbdirs
    ptl2d(nd) = 0.5_dp*rovertd(nd)
    ptr2d(nd) = 0.5_dp*rovertd(nd)
    ptl3d(nd) = (0.5_dp+cxb)*rovertd(nd)
    ptr3d(nd) = (0.5_dp+cxb)*rovertd(nd)
  END DO
  ptl2 = 0.5_dp*rovert - 1.0_dp
  ptr2 = 0.5_dp*rovert + 1.0_dp
  ptl3 = (0.5_dp+cxb)*rovert - 1.0_dp
  ptr3 = (0.5_dp+cxb)*rovert + 1.0_dp
  IF (omega .LT. ptl1) THEN
    DO nd=1,nbdirs
      cxd(nd) = 1.83_dp*omegad(nd)/omega**2 - 2.07_dp*2*omega*omegad(nd)&
&       /(omega**2)**2
    END DO
    cx = constant - 1.83_dp/omega + 2.07_dp/omega**2
  ELSE IF (omega .GE. ptl1 .AND. omega .LE. ptr1) THEN
    fl = constant - 1.83_dp/ptl1 + 2.07_dp/ptl1**2
    fr = 1.0_dp
    gl = 1.83_dp/ptl1**2 - 4.14_dp/ptl1**3
    gr = 0.0_dp
    DO nd=1,nbdirs
      grd(nd) = 0.0
      gld(nd) = 0.0
      frd(nd) = 0.0
      fld(nd) = 0.0
      ptr1d(nd) = 0.0
      ptl1d(nd) = 0.0
    END DO
    CALL CUBIC_SPLINE_EVAL_DV(ptl1, ptl1d, ptr1, ptr1d, fl, fld, fr, frd&
&                       , gl, gld, gr, grd, omega, omegad, cx, cxd, &
&                       nbdirs)
  ELSE IF (omega .GT. ptr1 .AND. omega .LT. ptl2) THEN
    cx = 1.0_dp
    DO nd=1,nbdirs
      cxd(nd) = 0.0
    END DO
  ELSE IF (omega .GE. ptl2 .AND. omega .LE. ptr2) THEN
    fl = 1.0_dp
    DO nd=1,nbdirs
      frd(nd) = -(0.2_dp*(2.0_dp*ptr2d(nd)*rovert-2.0_dp*ptr2*rovertd(nd&
&       ))/(cxb*rovert**2))
      grd(nd) = 0.4_dp*rovertd(nd)/cxb/rovert**2
      gld(nd) = 0.0
      fld(nd) = 0.0
    END DO
    fr = 1.0_dp + 0.2_dp/cxb*(1.0_dp-2.0_dp*ptr2/rovert)
    gl = 0.0_dp
    gr = -(0.4_dp/cxb/rovert)
    CALL CUBIC_SPLINE_EVAL_DV(ptl2, ptl2d, ptr2, ptr2d, fl, fld, fr, frd&
&                       , gl, gld, gr, grd, omega, omegad, cx, cxd, &
&                       nbdirs)
  ELSE IF (omega .GT. ptr2 .AND. omega .LT. ptl3) THEN
    DO nd=1,nbdirs
      cxd(nd) = -(0.2_dp*(2.0_dp*omegad(nd)*rovert-2.0_dp*omega*rovertd(&
&       nd))/(cxb*rovert**2))
    END DO
    cx = 1.0_dp + 0.2_dp/cxb*(1.0_dp-2.0_dp*omega/rovert)
  ELSE IF (omega .GE. ptl3 .AND. omega .LE. ptr3) THEN
    DO nd=1,nbdirs
      fld(nd) = -(0.2_dp*(2.0_dp*ptl3d(nd)*rovert-2.0_dp*ptl3*rovertd(nd&
&       ))/(cxb*rovert**2))
      gld(nd) = 0.4_dp*rovertd(nd)/cxb/rovert**2
      grd(nd) = 0.0
      frd(nd) = 0.0
    END DO
    fl = 1.0_dp + 0.2_dp/cxb*(1.0_dp-2.0_dp*ptl3/rovert)
    fr = 0.6_dp
    gl = -(0.4_dp/cxb/rovert)
    gr = 0.0_dp
    CALL CUBIC_SPLINE_EVAL_DV(ptl3, ptl3d, ptr3, ptr3d, fl, fld, fr, frd&
&                       , gl, gld, gr, grd, omega, omegad, cx, cxd, &
&                       nbdirs)
  ELSE
    cx = 0.6_dp
    DO nd=1,nbdirs
      cxd(nd) = 0.0
    END DO
  END IF
END SUBROUTINE CXSMOOTH_DV









!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of cubic_spline_eval in forward (tangent) mode:
!   variations   of useful results: poly
!   with respect to varying inputs: f1 f2 x g1 g2 x1 x2
SUBROUTINE CUBIC_SPLINE_EVAL_DV(x1, x1d, x2, x2d, f1, f1d, f2, f2d, g1, &
& g1d, g2, g2d, x, xd, poly, polyd, nbdirs)
  ! USE DIFFSIZES
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
!in
  REAL(dp), INTENT(IN) :: x1, x2, f1, f2, g1, g2, x
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: x1d, x2d, f1d, f2d, g1d&
& , g2d, xd
!out
  REAL(dp), INTENT(OUT) :: poly
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: polyd
!local
  REAL(dp) :: det, b11, b12, b13, b14, b21, b22, b23, b24, b31, b32, b33&
& , b34, b41, b42, b43, b44
  REAL(dp), DIMENSION(nbdirs) :: detd, b11d, b12d, b13d, b14d, b21d, &
& b22d, b23d, b24d, b31d, b32d, b33d, b34d, b41d, b42d, b43d, b44d
  REAL(dp), DIMENSION(4) :: a1, a2, a3, a4, b, coeff
  REAL(dp), DIMENSION(nbdirs, 4) :: a1d, a2d, a3d, a4d, bd, coeffd
  INTRINSIC KIND
  INTEGER :: nd
  INTEGER :: nbdirs
  a1(1) = x1**3
  a1(2) = x1**2
  a1(3) = x1
  a1(4) = 1.0_dp
  a2(1) = x2**3
  a2(2) = x2**2
  a2(3) = x2
  a2(4) = 1.0_dp
  a3(1) = 3.0_dp*x1**2
  a3(2) = 2.0_dp*x1
  a3(3) = 1.0_dp
  a3(4) = 0.0_dp
  a4(1) = 3.0_dp*x2**2
  a4(2) = 2.0_dp*x2
  a4(3) = 1.0_dp
  a4(4) = 0.0_dp
  b(1) = f1
  b(2) = f2
  b(3) = g1
  b(4) = g2
  det = a1(1)*(a2(2)*a3(3)*a4(4)+a2(3)*a3(4)*a4(2)+a2(4)*a3(2)*a4(3)-a2(&
&   4)*a3(3)*a4(2)-a2(3)*a3(2)*a4(4)-a2(2)*a3(4)*a4(3)) - a1(2)*(a2(1)*&
&   a3(3)*a4(4)+a2(3)*a3(4)*a4(1)+a2(4)*a3(1)*a4(3)-a2(4)*a3(3)*a4(1)-a2&
&   (3)*a3(1)*a4(4)-a2(1)*a3(4)*a4(3)) + a1(3)*(a2(1)*a3(2)*a4(4)+a2(2)*&
&   a3(4)*a4(1)+a2(4)*a3(1)*a4(2)-a2(4)*a3(2)*a4(1)-a2(2)*a3(1)*a4(4)-a2&
&   (1)*a3(4)*a4(2)) - a1(4)*(a2(1)*a3(2)*a4(3)+a2(2)*a3(3)*a4(1)+a2(3)*&
&   a3(1)*a4(2)-a2(3)*a3(2)*a4(1)-a2(2)*a3(1)*a4(3)-a2(1)*a3(3)*a4(2))
  b11 = a2(2)*a3(3)*a4(4) + a2(3)*a3(4)*a4(2) + a2(4)*a3(2)*a4(3) - a2(2&
&   )*a3(4)*a4(3) - a2(3)*a3(2)*a4(4) - a2(4)*a3(3)*a4(2)
  b12 = a1(2)*a3(4)*a4(3) + a1(3)*a3(2)*a4(4) + a1(4)*a3(3)*a4(2) - a1(2&
&   )*a3(3)*a4(4) - a1(3)*a3(4)*a4(2) - a1(4)*a3(2)*a4(3)
  b13 = a1(2)*a2(3)*a4(4) + a1(3)*a2(4)*a4(2) + a1(4)*a2(2)*a4(3) - a1(2&
&   )*a2(4)*a4(3) - a1(3)*a2(2)*a4(4) - a1(4)*a2(3)*a4(2)
  b14 = a1(2)*a2(4)*a3(3) + a1(3)*a2(2)*a3(4) + a1(4)*a2(3)*a3(2) - a1(2&
&   )*a2(3)*a3(4) - a1(3)*a2(4)*a3(2) - a1(4)*a2(2)*a3(3)
  b21 = a2(1)*a3(4)*a4(3) + a2(3)*a3(1)*a4(4) + a2(4)*a3(3)*a4(1) - a2(1&
&   )*a3(3)*a4(4) - a2(3)*a3(4)*a4(1) - a2(4)*a3(1)*a4(3)
  b22 = a1(1)*a3(3)*a4(4) + a1(3)*a3(4)*a4(1) + a1(4)*a3(1)*a4(3) - a1(1&
&   )*a3(4)*a4(3) - a1(3)*a3(1)*a4(4) - a1(4)*a3(3)*a4(1)
  b23 = a1(1)*a2(4)*a4(3) + a1(3)*a2(1)*a4(4) + a1(4)*a2(3)*a4(1) - a1(1&
&   )*a2(3)*a4(4) - a1(3)*a2(4)*a4(1) - a1(4)*a2(1)*a4(3)
  b24 = a1(1)*a2(3)*a3(4) + a1(3)*a2(4)*a3(1) + a1(4)*a2(1)*a3(3) - a1(1&
&   )*a2(4)*a3(3) - a1(3)*a2(1)*a3(4) - a1(4)*a2(3)*a3(1)
  b31 = a2(1)*a3(2)*a4(4) + a2(2)*a3(4)*a4(1) + a2(4)*a3(1)*a4(2) - a2(1&
&   )*a3(4)*a4(2) - a2(2)*a3(1)*a4(4) - a2(4)*a3(2)*a4(1)
  b32 = a1(1)*a3(4)*a4(2) + a1(2)*a3(1)*a4(4) + a1(4)*a3(2)*a4(1) - a1(1&
&   )*a3(2)*a4(4) - a1(2)*a3(4)*a4(1) - a1(4)*a3(1)*a4(2)
  b33 = a1(1)*a2(2)*a4(4) + a1(2)*a2(4)*a4(1) + a1(4)*a2(1)*a4(2) - a1(1&
&   )*a2(4)*a4(2) - a1(2)*a2(1)*a4(4) - a1(4)*a2(2)*a4(1)
  b34 = a1(1)*a2(4)*a3(2) + a1(2)*a2(1)*a3(4) + a1(4)*a2(2)*a3(1) - a1(1&
&   )*a2(2)*a3(4) - a1(2)*a2(4)*a3(1) - a1(4)*a2(1)*a3(2)
  b41 = a2(1)*a3(3)*a4(2) + a2(2)*a3(1)*a4(3) + a2(3)*a3(2)*a4(1) - a2(1&
&   )*a3(2)*a4(3) - a2(2)*a3(3)*a4(1) - a2(3)*a3(1)*a4(2)
  b42 = a1(1)*a3(2)*a4(3) + a1(2)*a3(3)*a4(1) + a1(3)*a3(1)*a4(2) - a1(1&
&   )*a3(3)*a4(2) - a1(2)*a3(1)*a4(3) - a1(3)*a3(2)*a4(1)
  b43 = a1(1)*a2(3)*a4(2) + a1(2)*a2(1)*a4(3) + a1(3)*a2(2)*a4(1) - a1(1&
&   )*a2(2)*a4(3) - a1(2)*a2(3)*a4(1) - a1(3)*a2(1)*a4(2)
  b44 = a1(1)*a2(2)*a3(3) + a1(2)*a2(3)*a3(1) + a1(3)*a2(1)*a3(2) - a1(1&
&   )*a2(3)*a3(2) - a1(2)*a2(1)*a3(3) - a1(3)*a2(2)*a3(1)
  coeff(1) = b11/det*b(1) + b12/det*b(2) + b13/det*b(3) + b14/det*b(4)
  coeff(2) = b21/det*b(1) + b22/det*b(2) + b23/det*b(3) + b24/det*b(4)
  coeff(3) = b31/det*b(1) + b32/det*b(2) + b33/det*b(3) + b34/det*b(4)
  coeff(4) = b41/det*b(1) + b42/det*b(2) + b43/det*b(3) + b44/det*b(4)
  DO nd=1,nbdirs
    a1d(nd, :) = 0.0
    a1d(nd, 1) = 3*x1**2*x1d(nd)
    a1d(nd, 2) = 2*x1*x1d(nd)
    a1d(nd, 3) = x1d(nd)
    a1d(nd, 4) = 0.0
    a2d(nd, :) = 0.0
    a2d(nd, 1) = 3*x2**2*x2d(nd)
    a2d(nd, 2) = 2*x2*x2d(nd)
    a2d(nd, 3) = x2d(nd)
    a2d(nd, 4) = 0.0
    a3d(nd, :) = 0.0
    a3d(nd, 1) = 3.0_dp*2*x1*x1d(nd)
    a3d(nd, 2) = 2.0_dp*x1d(nd)
    a3d(nd, 3) = 0.0
    a3d(nd, 4) = 0.0
    a4d(nd, :) = 0.0
    a4d(nd, 1) = 3.0_dp*2*x2*x2d(nd)
    a4d(nd, 2) = 2.0_dp*x2d(nd)
    a4d(nd, 3) = 0.0
    a4d(nd, 4) = 0.0
    bd(nd, :) = 0.0
    bd(nd, 1) = f1d(nd)
    bd(nd, 2) = f2d(nd)
    bd(nd, 3) = g1d(nd)
    bd(nd, 4) = g2d(nd)
    detd(nd) = a1d(nd, 1)*(a2(2)*a3(3)*a4(4)+a2(3)*a3(4)*a4(2)+a2(4)*a3(&
&     2)*a4(3)-a2(4)*a3(3)*a4(2)-a2(3)*a3(2)*a4(4)-a2(2)*a3(4)*a4(3)) + &
&     a1(1)*((a2d(nd, 2)*a3(3)+a2(2)*a3d(nd, 3))*a4(4)+a2(2)*a3(3)*a4d(&
&     nd, 4)+(a2d(nd, 3)*a3(4)+a2(3)*a3d(nd, 4))*a4(2)+a2(3)*a3(4)*a4d(&
&     nd, 2)+(a2d(nd, 4)*a3(2)+a2(4)*a3d(nd, 2))*a4(3)+a2(4)*a3(2)*a4d(&
&     nd, 3)-(a2d(nd, 4)*a3(3)+a2(4)*a3d(nd, 3))*a4(2)-a2(4)*a3(3)*a4d(&
&     nd, 2)-(a2d(nd, 3)*a3(2)+a2(3)*a3d(nd, 2))*a4(4)-a2(3)*a3(2)*a4d(&
&     nd, 4)-(a2d(nd, 2)*a3(4)+a2(2)*a3d(nd, 4))*a4(3)-a2(2)*a3(4)*a4d(&
&     nd, 3)) - a1d(nd, 2)*(a2(1)*a3(3)*a4(4)+a2(3)*a3(4)*a4(1)+a2(4)*a3&
&     (1)*a4(3)-a2(4)*a3(3)*a4(1)-a2(3)*a3(1)*a4(4)-a2(1)*a3(4)*a4(3)) -&
&     a1(2)*((a2d(nd, 1)*a3(3)+a2(1)*a3d(nd, 3))*a4(4)+a2(1)*a3(3)*a4d(&
&     nd, 4)+(a2d(nd, 3)*a3(4)+a2(3)*a3d(nd, 4))*a4(1)+a2(3)*a3(4)*a4d(&
&     nd, 1)+(a2d(nd, 4)*a3(1)+a2(4)*a3d(nd, 1))*a4(3)+a2(4)*a3(1)*a4d(&
&     nd, 3)-(a2d(nd, 4)*a3(3)+a2(4)*a3d(nd, 3))*a4(1)-a2(4)*a3(3)*a4d(&
&     nd, 1)-(a2d(nd, 3)*a3(1)+a2(3)*a3d(nd, 1))*a4(4)-a2(3)*a3(1)*a4d(&
&     nd, 4)-(a2d(nd, 1)*a3(4)+a2(1)*a3d(nd, 4))*a4(3)-a2(1)*a3(4)*a4d(&
&     nd, 3)) + a1d(nd, 3)*(a2(1)*a3(2)*a4(4)+a2(2)*a3(4)*a4(1)+a2(4)*a3&
&     (1)*a4(2)-a2(4)*a3(2)*a4(1)-a2(2)*a3(1)*a4(4)-a2(1)*a3(4)*a4(2)) +&
&     a1(3)*((a2d(nd, 1)*a3(2)+a2(1)*a3d(nd, 2))*a4(4)+a2(1)*a3(2)*a4d(&
&     nd, 4)+(a2d(nd, 2)*a3(4)+a2(2)*a3d(nd, 4))*a4(1)+a2(2)*a3(4)*a4d(&
&     nd, 1)+(a2d(nd, 4)*a3(1)+a2(4)*a3d(nd, 1))*a4(2)+a2(4)*a3(1)*a4d(&
&     nd, 2)-(a2d(nd, 4)*a3(2)+a2(4)*a3d(nd, 2))*a4(1)-a2(4)*a3(2)*a4d(&
&     nd, 1)-(a2d(nd, 2)*a3(1)+a2(2)*a3d(nd, 1))*a4(4)-a2(2)*a3(1)*a4d(&
&     nd, 4)-(a2d(nd, 1)*a3(4)+a2(1)*a3d(nd, 4))*a4(2)-a2(1)*a3(4)*a4d(&
&     nd, 2)) - a1d(nd, 4)*(a2(1)*a3(2)*a4(3)+a2(2)*a3(3)*a4(1)+a2(3)*a3&
&     (1)*a4(2)-a2(3)*a3(2)*a4(1)-a2(2)*a3(1)*a4(3)-a2(1)*a3(3)*a4(2)) -&
&     a1(4)*((a2d(nd, 1)*a3(2)+a2(1)*a3d(nd, 2))*a4(3)+a2(1)*a3(2)*a4d(&
&     nd, 3)+(a2d(nd, 2)*a3(3)+a2(2)*a3d(nd, 3))*a4(1)+a2(2)*a3(3)*a4d(&
&     nd, 1)+(a2d(nd, 3)*a3(1)+a2(3)*a3d(nd, 1))*a4(2)+a2(3)*a3(1)*a4d(&
&     nd, 2)-(a2d(nd, 3)*a3(2)+a2(3)*a3d(nd, 2))*a4(1)-a2(3)*a3(2)*a4d(&
&     nd, 1)-(a2d(nd, 2)*a3(1)+a2(2)*a3d(nd, 1))*a4(3)-a2(2)*a3(1)*a4d(&
&     nd, 3)-(a2d(nd, 1)*a3(3)+a2(1)*a3d(nd, 3))*a4(2)-a2(1)*a3(3)*a4d(&
&     nd, 2))
! entries of cof(A)
    b11d(nd) = (a2d(nd, 2)*a3(3)+a2(2)*a3d(nd, 3))*a4(4) + a2(2)*a3(3)*&
&     a4d(nd, 4) + (a2d(nd, 3)*a3(4)+a2(3)*a3d(nd, 4))*a4(2) + a2(3)*a3(&
&     4)*a4d(nd, 2) + (a2d(nd, 4)*a3(2)+a2(4)*a3d(nd, 2))*a4(3) + a2(4)*&
&     a3(2)*a4d(nd, 3) - (a2d(nd, 2)*a3(4)+a2(2)*a3d(nd, 4))*a4(3) - a2(&
&     2)*a3(4)*a4d(nd, 3) - (a2d(nd, 3)*a3(2)+a2(3)*a3d(nd, 2))*a4(4) - &
&     a2(3)*a3(2)*a4d(nd, 4) - (a2d(nd, 4)*a3(3)+a2(4)*a3d(nd, 3))*a4(2)&
&     - a2(4)*a3(3)*a4d(nd, 2)
    b12d(nd) = (a1d(nd, 2)*a3(4)+a1(2)*a3d(nd, 4))*a4(3) + a1(2)*a3(4)*&
&     a4d(nd, 3) + (a1d(nd, 3)*a3(2)+a1(3)*a3d(nd, 2))*a4(4) + a1(3)*a3(&
&     2)*a4d(nd, 4) + (a1d(nd, 4)*a3(3)+a1(4)*a3d(nd, 3))*a4(2) + a1(4)*&
&     a3(3)*a4d(nd, 2) - (a1d(nd, 2)*a3(3)+a1(2)*a3d(nd, 3))*a4(4) - a1(&
&     2)*a3(3)*a4d(nd, 4) - (a1d(nd, 3)*a3(4)+a1(3)*a3d(nd, 4))*a4(2) - &
&     a1(3)*a3(4)*a4d(nd, 2) - (a1d(nd, 4)*a3(2)+a1(4)*a3d(nd, 2))*a4(3)&
&     - a1(4)*a3(2)*a4d(nd, 3)
    b13d(nd) = (a1d(nd, 2)*a2(3)+a1(2)*a2d(nd, 3))*a4(4) + a1(2)*a2(3)*&
&     a4d(nd, 4) + (a1d(nd, 3)*a2(4)+a1(3)*a2d(nd, 4))*a4(2) + a1(3)*a2(&
&     4)*a4d(nd, 2) + (a1d(nd, 4)*a2(2)+a1(4)*a2d(nd, 2))*a4(3) + a1(4)*&
&     a2(2)*a4d(nd, 3) - (a1d(nd, 2)*a2(4)+a1(2)*a2d(nd, 4))*a4(3) - a1(&
&     2)*a2(4)*a4d(nd, 3) - (a1d(nd, 3)*a2(2)+a1(3)*a2d(nd, 2))*a4(4) - &
&     a1(3)*a2(2)*a4d(nd, 4) - (a1d(nd, 4)*a2(3)+a1(4)*a2d(nd, 3))*a4(2)&
&     - a1(4)*a2(3)*a4d(nd, 2)
    b14d(nd) = (a1d(nd, 2)*a2(4)+a1(2)*a2d(nd, 4))*a3(3) + a1(2)*a2(4)*&
&     a3d(nd, 3) + (a1d(nd, 3)*a2(2)+a1(3)*a2d(nd, 2))*a3(4) + a1(3)*a2(&
&     2)*a3d(nd, 4) + (a1d(nd, 4)*a2(3)+a1(4)*a2d(nd, 3))*a3(2) + a1(4)*&
&     a2(3)*a3d(nd, 2) - (a1d(nd, 2)*a2(3)+a1(2)*a2d(nd, 3))*a3(4) - a1(&
&     2)*a2(3)*a3d(nd, 4) - (a1d(nd, 3)*a2(4)+a1(3)*a2d(nd, 4))*a3(2) - &
&     a1(3)*a2(4)*a3d(nd, 2) - (a1d(nd, 4)*a2(2)+a1(4)*a2d(nd, 2))*a3(3)&
&     - a1(4)*a2(2)*a3d(nd, 3)
    b21d(nd) = (a2d(nd, 1)*a3(4)+a2(1)*a3d(nd, 4))*a4(3) + a2(1)*a3(4)*&
&     a4d(nd, 3) + (a2d(nd, 3)*a3(1)+a2(3)*a3d(nd, 1))*a4(4) + a2(3)*a3(&
&     1)*a4d(nd, 4) + (a2d(nd, 4)*a3(3)+a2(4)*a3d(nd, 3))*a4(1) + a2(4)*&
&     a3(3)*a4d(nd, 1) - (a2d(nd, 1)*a3(3)+a2(1)*a3d(nd, 3))*a4(4) - a2(&
&     1)*a3(3)*a4d(nd, 4) - (a2d(nd, 3)*a3(4)+a2(3)*a3d(nd, 4))*a4(1) - &
&     a2(3)*a3(4)*a4d(nd, 1) - (a2d(nd, 4)*a3(1)+a2(4)*a3d(nd, 1))*a4(3)&
&     - a2(4)*a3(1)*a4d(nd, 3)
    b22d(nd) = (a1d(nd, 1)*a3(3)+a1(1)*a3d(nd, 3))*a4(4) + a1(1)*a3(3)*&
&     a4d(nd, 4) + (a1d(nd, 3)*a3(4)+a1(3)*a3d(nd, 4))*a4(1) + a1(3)*a3(&
&     4)*a4d(nd, 1) + (a1d(nd, 4)*a3(1)+a1(4)*a3d(nd, 1))*a4(3) + a1(4)*&
&     a3(1)*a4d(nd, 3) - (a1d(nd, 1)*a3(4)+a1(1)*a3d(nd, 4))*a4(3) - a1(&
&     1)*a3(4)*a4d(nd, 3) - (a1d(nd, 3)*a3(1)+a1(3)*a3d(nd, 1))*a4(4) - &
&     a1(3)*a3(1)*a4d(nd, 4) - (a1d(nd, 4)*a3(3)+a1(4)*a3d(nd, 3))*a4(1)&
&     - a1(4)*a3(3)*a4d(nd, 1)
    b23d(nd) = (a1d(nd, 1)*a2(4)+a1(1)*a2d(nd, 4))*a4(3) + a1(1)*a2(4)*&
&     a4d(nd, 3) + (a1d(nd, 3)*a2(1)+a1(3)*a2d(nd, 1))*a4(4) + a1(3)*a2(&
&     1)*a4d(nd, 4) + (a1d(nd, 4)*a2(3)+a1(4)*a2d(nd, 3))*a4(1) + a1(4)*&
&     a2(3)*a4d(nd, 1) - (a1d(nd, 1)*a2(3)+a1(1)*a2d(nd, 3))*a4(4) - a1(&
&     1)*a2(3)*a4d(nd, 4) - (a1d(nd, 3)*a2(4)+a1(3)*a2d(nd, 4))*a4(1) - &
&     a1(3)*a2(4)*a4d(nd, 1) - (a1d(nd, 4)*a2(1)+a1(4)*a2d(nd, 1))*a4(3)&
&     - a1(4)*a2(1)*a4d(nd, 3)
    b24d(nd) = (a1d(nd, 1)*a2(3)+a1(1)*a2d(nd, 3))*a3(4) + a1(1)*a2(3)*&
&     a3d(nd, 4) + (a1d(nd, 3)*a2(4)+a1(3)*a2d(nd, 4))*a3(1) + a1(3)*a2(&
&     4)*a3d(nd, 1) + (a1d(nd, 4)*a2(1)+a1(4)*a2d(nd, 1))*a3(3) + a1(4)*&
&     a2(1)*a3d(nd, 3) - (a1d(nd, 1)*a2(4)+a1(1)*a2d(nd, 4))*a3(3) - a1(&
&     1)*a2(4)*a3d(nd, 3) - (a1d(nd, 3)*a2(1)+a1(3)*a2d(nd, 1))*a3(4) - &
&     a1(3)*a2(1)*a3d(nd, 4) - (a1d(nd, 4)*a2(3)+a1(4)*a2d(nd, 3))*a3(1)&
&     - a1(4)*a2(3)*a3d(nd, 1)
    b31d(nd) = (a2d(nd, 1)*a3(2)+a2(1)*a3d(nd, 2))*a4(4) + a2(1)*a3(2)*&
&     a4d(nd, 4) + (a2d(nd, 2)*a3(4)+a2(2)*a3d(nd, 4))*a4(1) + a2(2)*a3(&
&     4)*a4d(nd, 1) + (a2d(nd, 4)*a3(1)+a2(4)*a3d(nd, 1))*a4(2) + a2(4)*&
&     a3(1)*a4d(nd, 2) - (a2d(nd, 1)*a3(4)+a2(1)*a3d(nd, 4))*a4(2) - a2(&
&     1)*a3(4)*a4d(nd, 2) - (a2d(nd, 2)*a3(1)+a2(2)*a3d(nd, 1))*a4(4) - &
&     a2(2)*a3(1)*a4d(nd, 4) - (a2d(nd, 4)*a3(2)+a2(4)*a3d(nd, 2))*a4(1)&
&     - a2(4)*a3(2)*a4d(nd, 1)
    b32d(nd) = (a1d(nd, 1)*a3(4)+a1(1)*a3d(nd, 4))*a4(2) + a1(1)*a3(4)*&
&     a4d(nd, 2) + (a1d(nd, 2)*a3(1)+a1(2)*a3d(nd, 1))*a4(4) + a1(2)*a3(&
&     1)*a4d(nd, 4) + (a1d(nd, 4)*a3(2)+a1(4)*a3d(nd, 2))*a4(1) + a1(4)*&
&     a3(2)*a4d(nd, 1) - (a1d(nd, 1)*a3(2)+a1(1)*a3d(nd, 2))*a4(4) - a1(&
&     1)*a3(2)*a4d(nd, 4) - (a1d(nd, 2)*a3(4)+a1(2)*a3d(nd, 4))*a4(1) - &
&     a1(2)*a3(4)*a4d(nd, 1) - (a1d(nd, 4)*a3(1)+a1(4)*a3d(nd, 1))*a4(2)&
&     - a1(4)*a3(1)*a4d(nd, 2)
    b33d(nd) = (a1d(nd, 1)*a2(2)+a1(1)*a2d(nd, 2))*a4(4) + a1(1)*a2(2)*&
&     a4d(nd, 4) + (a1d(nd, 2)*a2(4)+a1(2)*a2d(nd, 4))*a4(1) + a1(2)*a2(&
&     4)*a4d(nd, 1) + (a1d(nd, 4)*a2(1)+a1(4)*a2d(nd, 1))*a4(2) + a1(4)*&
&     a2(1)*a4d(nd, 2) - (a1d(nd, 1)*a2(4)+a1(1)*a2d(nd, 4))*a4(2) - a1(&
&     1)*a2(4)*a4d(nd, 2) - (a1d(nd, 2)*a2(1)+a1(2)*a2d(nd, 1))*a4(4) - &
&     a1(2)*a2(1)*a4d(nd, 4) - (a1d(nd, 4)*a2(2)+a1(4)*a2d(nd, 2))*a4(1)&
&     - a1(4)*a2(2)*a4d(nd, 1)
    b34d(nd) = (a1d(nd, 1)*a2(4)+a1(1)*a2d(nd, 4))*a3(2) + a1(1)*a2(4)*&
&     a3d(nd, 2) + (a1d(nd, 2)*a2(1)+a1(2)*a2d(nd, 1))*a3(4) + a1(2)*a2(&
&     1)*a3d(nd, 4) + (a1d(nd, 4)*a2(2)+a1(4)*a2d(nd, 2))*a3(1) + a1(4)*&
&     a2(2)*a3d(nd, 1) - (a1d(nd, 1)*a2(2)+a1(1)*a2d(nd, 2))*a3(4) - a1(&
&     1)*a2(2)*a3d(nd, 4) - (a1d(nd, 2)*a2(4)+a1(2)*a2d(nd, 4))*a3(1) - &
&     a1(2)*a2(4)*a3d(nd, 1) - (a1d(nd, 4)*a2(1)+a1(4)*a2d(nd, 1))*a3(2)&
&     - a1(4)*a2(1)*a3d(nd, 2)
    b41d(nd) = (a2d(nd, 1)*a3(3)+a2(1)*a3d(nd, 3))*a4(2) + a2(1)*a3(3)*&
&     a4d(nd, 2) + (a2d(nd, 2)*a3(1)+a2(2)*a3d(nd, 1))*a4(3) + a2(2)*a3(&
&     1)*a4d(nd, 3) + (a2d(nd, 3)*a3(2)+a2(3)*a3d(nd, 2))*a4(1) + a2(3)*&
&     a3(2)*a4d(nd, 1) - (a2d(nd, 1)*a3(2)+a2(1)*a3d(nd, 2))*a4(3) - a2(&
&     1)*a3(2)*a4d(nd, 3) - (a2d(nd, 2)*a3(3)+a2(2)*a3d(nd, 3))*a4(1) - &
&     a2(2)*a3(3)*a4d(nd, 1) - (a2d(nd, 3)*a3(1)+a2(3)*a3d(nd, 1))*a4(2)&
&     - a2(3)*a3(1)*a4d(nd, 2)
    b42d(nd) = (a1d(nd, 1)*a3(2)+a1(1)*a3d(nd, 2))*a4(3) + a1(1)*a3(2)*&
&     a4d(nd, 3) + (a1d(nd, 2)*a3(3)+a1(2)*a3d(nd, 3))*a4(1) + a1(2)*a3(&
&     3)*a4d(nd, 1) + (a1d(nd, 3)*a3(1)+a1(3)*a3d(nd, 1))*a4(2) + a1(3)*&
&     a3(1)*a4d(nd, 2) - (a1d(nd, 1)*a3(3)+a1(1)*a3d(nd, 3))*a4(2) - a1(&
&     1)*a3(3)*a4d(nd, 2) - (a1d(nd, 2)*a3(1)+a1(2)*a3d(nd, 1))*a4(3) - &
&     a1(2)*a3(1)*a4d(nd, 3) - (a1d(nd, 3)*a3(2)+a1(3)*a3d(nd, 2))*a4(1)&
&     - a1(3)*a3(2)*a4d(nd, 1)
    b43d(nd) = (a1d(nd, 1)*a2(3)+a1(1)*a2d(nd, 3))*a4(2) + a1(1)*a2(3)*&
&     a4d(nd, 2) + (a1d(nd, 2)*a2(1)+a1(2)*a2d(nd, 1))*a4(3) + a1(2)*a2(&
&     1)*a4d(nd, 3) + (a1d(nd, 3)*a2(2)+a1(3)*a2d(nd, 2))*a4(1) + a1(3)*&
&     a2(2)*a4d(nd, 1) - (a1d(nd, 1)*a2(2)+a1(1)*a2d(nd, 2))*a4(3) - a1(&
&     1)*a2(2)*a4d(nd, 3) - (a1d(nd, 2)*a2(3)+a1(2)*a2d(nd, 3))*a4(1) - &
&     a1(2)*a2(3)*a4d(nd, 1) - (a1d(nd, 3)*a2(1)+a1(3)*a2d(nd, 1))*a4(2)&
&     - a1(3)*a2(1)*a4d(nd, 2)
    b44d(nd) = (a1d(nd, 1)*a2(2)+a1(1)*a2d(nd, 2))*a3(3) + a1(1)*a2(2)*&
&     a3d(nd, 3) + (a1d(nd, 2)*a2(3)+a1(2)*a2d(nd, 3))*a3(1) + a1(2)*a2(&
&     3)*a3d(nd, 1) + (a1d(nd, 3)*a2(1)+a1(3)*a2d(nd, 1))*a3(2) + a1(3)*&
&     a2(1)*a3d(nd, 2) - (a1d(nd, 1)*a2(3)+a1(1)*a2d(nd, 3))*a3(2) - a1(&
&     1)*a2(3)*a3d(nd, 2) - (a1d(nd, 2)*a2(1)+a1(2)*a2d(nd, 1))*a3(3) - &
&     a1(2)*a2(1)*a3d(nd, 3) - (a1d(nd, 3)*a2(2)+a1(3)*a2d(nd, 2))*a3(1)&
&     - a1(3)*a2(2)*a3d(nd, 1)
!solve the equation Ax = b; x = A^(-1)b
    coeffd(nd, :) = 0.0
    coeffd(nd, 1) = (b11d(nd)*det-b11*detd(nd))*b(1)/det**2 + b11*bd(nd&
&     , 1)/det + (b12d(nd)*det-b12*detd(nd))*b(2)/det**2 + b12*bd(nd, 2)&
&     /det + (b13d(nd)*det-b13*detd(nd))*b(3)/det**2 + b13*bd(nd, 3)/det&
&     + (b14d(nd)*det-b14*detd(nd))*b(4)/det**2 + b14*bd(nd, 4)/det
    coeffd(nd, 2) = (b21d(nd)*det-b21*detd(nd))*b(1)/det**2 + b21*bd(nd&
&     , 1)/det + (b22d(nd)*det-b22*detd(nd))*b(2)/det**2 + b22*bd(nd, 2)&
&     /det + (b23d(nd)*det-b23*detd(nd))*b(3)/det**2 + b23*bd(nd, 3)/det&
&     + (b24d(nd)*det-b24*detd(nd))*b(4)/det**2 + b24*bd(nd, 4)/det
    coeffd(nd, 3) = (b31d(nd)*det-b31*detd(nd))*b(1)/det**2 + b31*bd(nd&
&     , 1)/det + (b32d(nd)*det-b32*detd(nd))*b(2)/det**2 + b32*bd(nd, 2)&
&     /det + (b33d(nd)*det-b33*detd(nd))*b(3)/det**2 + b33*bd(nd, 3)/det&
&     + (b34d(nd)*det-b34*detd(nd))*b(4)/det**2 + b34*bd(nd, 4)/det
    coeffd(nd, 4) = (b41d(nd)*det-b41*detd(nd))*b(1)/det**2 + b41*bd(nd&
&     , 1)/det + (b42d(nd)*det-b42*detd(nd))*b(2)/det**2 + b42*bd(nd, 2)&
&     /det + (b43d(nd)*det-b43*detd(nd))*b(3)/det**2 + b43*bd(nd, 3)/det&
&     + (b44d(nd)*det-b44*detd(nd))*b(4)/det**2 + b44*bd(nd, 4)/det
    polyd(nd) = coeffd(nd, 1)*x**3 + coeff(1)*3*x**2*xd(nd) + coeffd(nd&
&     , 2)*x**2 + coeff(2)*2*x*xd(nd) + coeffd(nd, 3)*x + coeff(3)*xd(nd&
&     ) + coeffd(nd, 4)
  END DO
  poly = coeff(1)*x**3 + coeff(2)*x**2 + coeff(3)*x + coeff(4)
END SUBROUTINE CUBIC_SPLINE_EVAL_DV
