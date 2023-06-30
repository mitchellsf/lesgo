!!
!!  Copyright (C) 2016  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
module qeqwm
!*******************************************************************************
use types, only : rprec

private
public qeqwm_initialize, qeqwm_finalize, lagrangian_rewm, eqwm_compute, &
    qeqwm_write_checkpoint, qeqwm_read_checkpoint, &
    utau, utx, uty, ssx, ssy, Ts, Us, vtx, vty, utau_filt

!Lagrangian relaxation wall model variables
! friction velocities
real(rprec), dimension(:,:), allocatable :: utau, utx, uty
! wall stress angle components
real(rprec), dimension(:,:), allocatable :: ssx, ssy
! relaxation time scale
real(rprec), dimension(:,:), allocatable :: Ts
! velocity at first grid point computed from assumed profile and utau
real(rprec), dimension(:,:), allocatable :: Us
! advecting velocities in lagrangian rewm
real(rprec), dimension(:,:), allocatable :: vtx, vty
!! reynolds stress evaluated at Delta
!real(rprec), dimension(:,:), allocatable :: upwp, vpwp
real(rprec), dimension(:,:), allocatable :: utau_filt

contains

!*******************************************************************************
subroutine qeqwm_initialize
!*******************************************************************************
use param, only : nx, ny, nz, dz, coord, wmpt
use sim_param, only : u

implicit none

allocate(utau(nx,ny))
allocate(utx(nx,ny))
allocate(uty(nx,ny))
allocate(ssx(nx,ny))
allocate(ssy(nx,ny))
allocate(Ts(nx,ny))
allocate(Us(nx,ny))
allocate(vtx(nx,ny))
allocate(vty(nx,ny))
!allocate(upwp(nx,ny))
!allocate(vpwp(nx,ny))
allocate(utau_filt(nx,ny))

utau = 1._rprec
utx = 1._rprec
uty = 0._rprec
ssx = 1._rprec
ssy = 1._rprec
if (coord==0) then
    Us(1:nx,1:ny) = sum(u(1:nx,1:ny,wmpt))/nx/ny
else
    Us(1:nx,1:ny) = sum(u(1:nx,1:ny,nz-wmpt))/nx/ny
endif
Ts(1:nx,1:ny) = Us(1:nx,1:ny)*(wmpt-0.5_rprec)*dz
vtx = 0._rprec
vty = 0._rprec
!upwp = -utx*utau
!vpwp = -uty*utau
utau_filt = 1._rprec

end subroutine qeqwm_initialize
          
!*******************************************************************************
subroutine qeqwm_finalize
!*******************************************************************************
          
implicit none

deallocate(utau)
deallocate(utx)
deallocate(uty)
deallocate(ssx)
deallocate(ssy)
deallocate(Ts)
deallocate(Us)
deallocate(vtx)
deallocate(vty)
!deallocate(upwp)
!deallocate(vpwp)
deallocate(utau_filt)

end subroutine qeqwm_finalize

!*******************************************************************************
subroutine lagrangian_rewm(Deltay,Ud,theta_d,dpdxbar1,dpdybar1,&
    dpdxbar2,dpdybar2,twxbar,twybar)
!*******************************************************************************
use param, only : ld, nx, ny, nu_molec, dt, dt_f, jt, coord, jt_total, &
    total_time, vonk
use grid_m
use param, only : path
!use eqwm_dyn, only : dynamic_smag_model

implicit none
real(rprec), dimension(nx,ny), intent(in) :: Ud, theta_d
real(rprec), dimension(nx,ny), intent(in) :: dpdxbar1, dpdybar1,&
    dpdxbar2, dpdybar2
real(rprec), intent(in) :: Deltay
real(rprec), dimension(nx,ny), intent(out) :: twxbar, twybar
real(rprec), dimension(nx,ny) :: sxp, syp, Deltap, fu, reds, delta_star, retdu,&
        taud, redu, theta_delta, rhsx, rhsy, rhsx0, rhsy0, dpdxu, psi_p,&
        taudx, taudy
real(rprec), dimension(ld,ny) :: temp1, temp2
real(rprec), dimension(nx,ny) :: tmp
!real(rprec), dimension(nx,ny):: lp, dudz_tot, upwp_eq, vpwp_eq, Tt, eps
real(rprec), pointer, dimension(:) :: x, y
real(rprec) :: x0, y0
integer :: i, j, imax, jmax, fid
integer, dimension(2) :: temp
character*50 :: fname
real(rprec), dimension(nx,ny) :: fu_const
real(rprec), dimension(nx,ny) :: utxp, utyp

nullify(x,y)
x => grid % x
y => grid % y

! Lagrangian relaxation wall model
sxp = ssx
syp = ssy
ssx = utx/utau
ssy = uty/utau
Deltap = utau*Deltay/nu_molec ! Delta_y in plus units (from previous tstep)
call velocity_fit(Deltap,fu)
Ts = fu*Deltay/utau
!! slowly varying time scale ---------------------------------
!call velocity_fit(Deltay/nu_molec*utau_filt,fu_const)
!Ts = fu_const*Deltay/utau_filt
!! -----------------------------------------------------
Us = utau*fu
reds = Us*Deltay/nu_molec
call integral_scale_fits(Deltap,fu,delta_star,theta_delta)
vtx = (1-delta_star-theta_delta)*utx*fu
vty = (1-delta_star-theta_delta)*uty*fu
redu = Ud*Deltay/nu_molec
dpdxu = dpdxbar2*cos(theta_d) + dpdybar2*sin(theta_d)
psi_p = abs(dpdxu)*Deltay**3.0/nu_molec**2.0
call retd_pres_calc(redu,psi_p,retdu)
!call retd_fit_calc(redu,retdu)
taud = (retdu*nu_molec/Deltay)**2
taudx = taud*cos(theta_d)
taudy = taud*sin(theta_d)

!call rs_model(Deltay,taud,taudx,taudy,theta_d)

!call dynamic_smag_model(temp1,temp2)
!taudx = temp1(1:nx,1:ny)
!taudy = temp2(1:nx,1:ny)
!taud = sqrt(taudx**2.0+taudy**2.0)

rhsx = utx + dt*(utau*delta_star*(ssx-sxp)/dt_f &
    + (-Deltay/utau*(dpdxbar1-dpdxbar2) &
    + taudx/utau - utx)/Ts)
rhsy = uty + dt*(utau*delta_star*(ssy-syp)/dt_f &
    + (-Deltay/utau*(dpdybar1-dpdybar2) &
    + taudy/utau - uty)/Ts)
!rhsx = utx + dt*(utau*delta_star*(ssx-sxp)/dt_f &
!    + (-Deltay/utau*(dpdxbar1-dpdxu*cos(theta_d)) &
!    + taud*cos(theta_d)/utau - utx)/Ts)
!rhsy = uty + dt*(utau*delta_star*(ssy-syp)/dt_f &
!    + (-Deltay/utau*(dpdybar1-dpdxu*sin(theta_d)) &
!    + taud*sin(theta_d)/utau - uty)/Ts)
!rhsx = utx + dt*(utau*delta_star*(ssx-sxp)/dt_f &
!    + (taudx/utau - utx)/Ts)
!rhsy = uty + dt*(utau*delta_star*(ssy-syp)/dt_f &
!    + (taudy/utau - uty)/Ts)

do i = 1,nx
do j = 1,ny
    x0 = x(i) - vtx(i,j)*dt
    y0 = y(j) - vty(i,j)*dt
    call bilinear_interp(rhsx,x0,y0,rhsx0(i,j))
    call bilinear_interp(rhsy,x0,y0,rhsy0(i,j))
enddo
enddo

!store utx and uty at previous time step
!(not required, just `compute_larte_terms' function)
utxp = utx
utyp = uty

! Lagrangian version
utx = rhsx0
uty = rhsy0

!! Eulerian version
!utx = rhsx
!uty = rhsy

utau = sqrt(utx**2 + uty**2)

twxbar = utau*utx
twybar = utau*uty

!if (coord==0 .and. mod(jt_total,100)==0) then
!    call dynamic_smag_model2
!endif


!if (coord==0) then
!    call compute_larte_terms(utxp,utyp,taudx,taudy,delta_star)
!endif

end subroutine lagrangian_rewm

!*******************************************************************************
subroutine compute_larte_terms(utxp,utyp,taudx,taudy,delta_star)
!*******************************************************************************
use param, only : jt_total, total_time, path, dt, ld, nx, ny
use types, only : rprec
use derivatives, only : ddxy_plane
implicit none

real(rprec), dimension(nx,ny), intent(in) :: utxp, utyp, taudx, taudy, delta_star
real(rprec), dimension(nx,ny) :: utaup
real(rprec), dimension(nx,ny) :: dutxdt,dutydt,taudx_larte,taudy_larte,utx_larte,&
    uty_larte,dsxdt_larte,dsydt_larte,advx_larte,advy_larte
real(rprec), dimension(ld,ny) :: temp,dutxdx,dutxdy,dutydx,dutydy
integer :: fid
character*50 :: fname

dutxdt = (utx-utxp)/dt
dutydt = (uty-utyp)/dt

temp(1:nx,1:ny) = utx
call ddxy_plane(temp,dutxdx,dutxdy)
temp(1:nx,1:ny) = uty
call ddxy_plane(temp,dutydx,dutydy)
!dutxdx = 1
!dutxdy = 1
!dutydx = 1
!dutydy = 1
advx_larte = vtx*dutxdx(1:nx,1:ny) + vty*dutxdy(1:nx,1:ny)
advy_larte = vtx*dutydx(1:nx,1:ny) + vty*dutydy(1:nx,1:ny)

taudx_larte = taudx/utau/Ts
taudy_larte = taudy/utau/Ts

utx_larte = utx/Ts
uty_larte = uty/Ts

utaup = sqrt(utxp**2.0+utyp**2.0)
dsxdt_larte = utau*delta_star*(utx/utau-utxp/utaup)/dt
dsydt_larte = utau*delta_star*(uty/utau-utyp/utaup)/dt

fname = path // 'output/larte_terms.dat'
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,&
    total_time,&
    dt,&
    sum(dutxdt)/nx/ny,&
    sum(dutydt)/nx/ny,&
    sum(advx_larte)/nx/ny,&
    sum(advy_larte)/nx/ny,&
    sum(taudx_larte)/nx/ny,&
    sum(taudy_larte)/nx/ny,&
    sum(utx_larte)/nx/ny,&
    sum(uty_larte)/nx/ny,&
    sum(dsxdt_larte)/nx/ny,&
    sum(dsydt_larte)/nx/ny
close(fid)

end subroutine compute_larte_terms

!*******************************************************************************
subroutine eqwm_compute(Deltay,Ud,theta_d,dpdxbar1,dpdybar1,&
    dpdxbar2,dpdybar2,twxbar,twybar)
!*******************************************************************************
use param, only : nx, ny, nu_molec

implicit none
real(rprec), dimension(nx,ny), intent(in) :: Ud, theta_d
real(rprec), dimension(nx,ny), intent(in) :: dpdxbar1, dpdybar1,&
    dpdxbar2, dpdybar2
real(rprec), intent(in) :: Deltay
real(rprec), dimension(nx,ny), intent(out) :: twxbar, twybar
real(rprec), dimension(nx,ny) :: Deltap, fu, redu, retdu, taud, taudx, taudy, &
    dpdxu, psi_p

Deltap = utau*Deltay/nu_molec ! Delta_y in plus units (from previous tstep)
call velocity_fit(Deltap,fu)
Ts = fu*Deltay/utau
redu = Ud*Deltay/nu_molec
dpdxu = dpdxbar2*cos(theta_d) + dpdybar2*sin(theta_d)
psi_p = abs(dpdxu)*Deltay**4.0/nu_molec**2.0
call retd_pres_calc(redu,psi_p,retdu)
taud = (retdu*nu_molec/Deltay)**2
taudx = taud*cos(theta_d)
taudy = taud*sin(theta_d)

twxbar = -Deltay*(dpdxbar1-dpdxbar2) + taudx
twybar = -Deltay*(dpdybar1-dpdybar2) + taudy

utau = (twxbar**2+twybar**2)**0.25
utx = twxbar/utau
uty = twybar/utau

end subroutine eqwm_compute

!*******************************************************************************
subroutine retd_fit_calc(redelta_in,retd_out)
!*******************************************************************************
use param, only : nx, ny

implicit none
real(rprec), dimension(nx,ny), intent(in) :: redelta_in
real(rprec), dimension(nx,ny), intent(out) :: retd_out
real(rprec) :: kappa3
real(rprec), dimension(nx,ny) :: kappa4, beta1, beta2

kappa3 = 0.005_rprec
beta1 = (1._rprec + 0.155_rprec*redelta_in**-0.03_rprec)**-1
beta2 = 1.7 - (1._rprec + 36._rprec*redelta_in**-0.75_rprec)**-1
kappa4 = kappa3**(beta1 - 0.5_rprec)

retd_out = kappa4*redelta_in**beta1*                            &
    (1._rprec + (kappa3*redelta_in)**-beta2)**((beta1-0.5_rprec)/beta2)

end subroutine retd_fit_calc

!*******************************************************************************
subroutine retd_fit_calc2(redelta_in,retd_out)
!*******************************************************************************
use param, only : ld, ny

implicit none
real(rprec), dimension(ld,ny), intent(in) :: redelta_in
real(rprec), dimension(ld,ny), intent(out) :: retd_out
real(rprec) :: kappa3
real(rprec), dimension(ld,ny) :: kappa4, beta1, beta2

kappa3 = 0.005_rprec
beta1 = (1._rprec + 0.155_rprec*redelta_in**-0.03_rprec)**-1
beta2 = 1.7 - (1._rprec + 36._rprec*redelta_in**-0.75_rprec)**-1
kappa4 = kappa3**(beta1 - 0.5_rprec)

retd_out = kappa4*redelta_in**beta1*                            &
    (1._rprec + (kappa3*redelta_in)**-beta2)**((beta1-0.5_rprec)/beta2)

end subroutine retd_fit_calc2

!*******************************************************************************
subroutine retd_pres_calc(redelta_in,psi_p,retd_out)
!*******************************************************************************
use param, only : nx, ny, nu_molec

implicit none
real(rprec), dimension(nx,ny), intent(in) :: redelta_in, psi_p
real(rprec), dimension(nx,ny), intent(out) :: retd_out
real(rprec), dimension(nx,ny) :: retd_eq, retd_min, pp

call retd_fit_calc(redelta_in,retd_eq)
retd_min = 1.5*psi_p**0.39*(1.0+(1000.0/psi_p)**2.0)**-0.055
pp = 2.5 - 0.6*(1.0 + tanh(2.0*log10(psi_p)**-6.0))
retd_out = (retd_min**pp + retd_eq**pp)**(1.0/pp)

end subroutine retd_pres_calc

!*******************************************************************************
subroutine integral_scale_fits(Deltap,fu,delta_star,theta_delta)
!*******************************************************************************
use param, only : nx, ny

implicit none
real(rprec) :: kappa, c1, c2, c3, c4, c5, c6, c7, c8
real(rprec), dimension(nx,ny) :: gamma1, gamma2, delta_log, theta_log, red
real(rprec), dimension(nx,ny), intent(in) :: Deltap, fu
real(rprec), dimension(nx,ny), intent(out) :: delta_star, theta_delta

kappa = 0.4_rprec
c1 = 23.664_rprec
c2 = 0.0016_rprec
c3 = 1.516_rprec
c4 = 1.177_rprec
c5 = -103.5_rprec
c6 = 2586._rprec
c7 = 0.00154_rprec
c8 = 2.475_rprec
red = Deltap*fu
gamma1 = 1._rprec/(1._rprec+c2*red**c4)
gamma2 = 1._rprec/(1._rprec+c7*red)
delta_log = c1/red + Deltap/red/kappa
theta_log = (c5+Deltap/kappa)/red + Deltap/red**2.0*(c6-2.0*Deltap/kappa**2.0)
delta_star = 0.5_rprec*gamma1 + (1._rprec-gamma1)**c3*delta_log
theta_delta = 1._rprec/6._rprec*gamma2 + (1._rprec-gamma2)**c8*theta_log

end subroutine integral_scale_fits

!*******************************************************************************
subroutine velocity_fit(Deltap,fu)
!*******************************************************************************
use param, only : nx, ny

implicit none
real(rprec) :: kappa,B,k2,k1,beta
real(rprec), dimension(nx,ny), intent(in) :: Deltap
real(rprec), dimension(nx,ny), intent(out) :: fu

kappa = 0.4_rprec
B = 4.95_rprec
k2 = 9.753_rprec
k1 = 1._rprec/kappa*log(k2) + B
beta = 1.903_rprec
fu = (1._rprec/kappa*log(k2+Deltap)+B)* &
    (1._rprec + ((k1**-1._rprec)*Deltap)**-beta)**(-1._rprec/beta)

end subroutine velocity_fit

!*******************************************************************************
subroutine bilinear_interp(u,x0,y0,u0)
!*******************************************************************************
use param, only : nx, ny, L_x, dx, L_y, dy
use grid_m

implicit none
real(rprec), dimension(nx,ny), intent(in) :: u
real(rprec), intent(in) :: x0, y0
real(rprec), intent(out) :: u0
real(rprec) :: eps, xp, yp, x0p, y0p
integer :: i1, i2, j1, j2
real(rprec), pointer, dimension(:) :: x, y

nullify(x,y)
x => grid % x
y => grid % y

eps = 10._rprec**-5

xp = modulo(x0,L_x)
i1 = ceiling(abs(xp-eps)/dx)
i2 = ceiling(modulo(real(i1+1),nx+eps))

yp = modulo(y0,L_y)
j1 = ceiling(abs(yp-eps)/dy)
j2 = ceiling(modulo(real(j1+1),ny+eps))

x0p = (xp - x(i1))/dx
y0p = (yp - y(j1))/dy

u0 = u(i1,j1)*(1._rprec-x0p)*(1._rprec-y0p) + u(i2,j1)*x0p*(1._rprec-y0p)&
    + u(i1,j2)*(1._rprec-x0p)*y0p + u(i2,j2)*x0p*y0p

end subroutine bilinear_interp

!*******************************************************************************
subroutine dynamic_model
!*******************************************************************************
use param, only : ld, ny, dz, jt_total, nu_molec
use sim_param, only : u,v,w
use test_filtermodule
use types, only : rprec
use functions, only : dealias_mult

implicit none

real(rprec) :: Deltay
real(rprec), dimension(ld,ny) :: u_filt,v_filt,w_filt,uw_filt,vw_filt,L13,L23,&
    Ud,retd,taud,theta_d,F13,F23,F13_ufilt,F23_ufilt,F13_filt,F23_filt,&
    M13,M23,alpha
integer :: jz

do jz = 1, 3

    Deltay = (jz-0.5_rprec)*dz
    
    u_filt = u(:,:,jz)
    v_filt = v(:,:,jz)
    if (jz==1) then
        w_filt = w(:,:,jz+1)*0.25_rprec
    else
        w_filt = 0.5_rprec*(w(:,:,jz)+w(:,:,jz+1))
    endif
    uw_filt = dealias_mult(u_filt,w_filt)
    vw_filt = dealias_mult(v_filt,w_filt)
    
    call test_filter(u_filt)
    call test_filter(v_filt)
    call test_filter(w_filt)
    call test_filter(uw_filt)
    call test_filter(vw_filt)
    
    L13 = uw_filt - dealias_mult(u_filt,w_filt)
    L23 = vw_filt - dealias_mult(v_filt,w_filt)
    
    Ud = sqrt(u(:,:,jz)**2._rprec+v(:,:,jz)**2._rprec)
    call retd_fit_calc2(Ud*Deltay/nu_molec,retd)
    taud = (retd*nu_molec/Deltay)**2._rprec
    theta_d = atan2(v(:,:,jz),u(:,:,jz))
    F13 = taud*cos(theta_d)
    F23 = taud*sin(theta_d)
    
    Ud = sqrt(u_filt**2._rprec+v_filt**2._rprec)
    call retd_fit_calc2(Ud*Deltay/nu_molec,retd)
    taud = (retd*nu_molec/Deltay)**2._rprec
    theta_d = atan2(v_filt,u_filt)
    F13_ufilt = taud*cos(theta_d)
    F23_ufilt = taud*sin(theta_d)
    
    F13_filt = F13
    F23_filt = F23
    call test_filter(F13_filt)
    call test_filter(F23_filt)
    
    M13 = F13_filt - F13_ufilt
    M23 = F23_filt - F23_ufilt
    
    alpha = (dealias_mult(L13,M13)+dealias_mult(L23,M23))/&
        (dealias_mult(M13,M13)+dealias_mult(M23,M23))
    
    if (jt_total==1000) then
        call write_dynamic_model(jz,u_filt,v_filt,w_filt,uw_filt,vw_filt,L13,L23,&
            F13,F23,F13_filt,F23_filt,F13_ufilt,F23_ufilt,M13,M23,alpha)
    endif

enddo

end subroutine dynamic_model

!*******************************************************************************
subroutine dynamic_smag_model2
!*******************************************************************************
use param, only : ld,nx,ny,dx,dy,dz,vonk,jt_total,nu_molec
use sim_param, only : u,v,w,dudx,dudy,dvdx,dvdy,dwdx,dwdy
use test_filtermodule
use types, only : rprec
use functions, only : dealias_mult

implicit none

real(rprec) :: Deltay,Delta_tilde,Delta_hat,Cs,taudx_tot,taudy_tot
real(rprec), dimension(ld,ny) :: u_filt,v_filt,w_filt,&
    uu_filt,vv_filt,ww_filt,uv_filt,uw_filt,vw_filt,&
    L11,L22,L33,L12,L13,L23,&
    M11,M22,M33,M12,M13,M23,&
    u_w,v_w,Ud,Deltap,taud,lp,dudz_tot,nu_T,taudx_dyn,taudy_dyn,&
    ux,uy,uz,vx,vy,vz,wz,&
    S11,S22,S33,S12,S13,S23,S,&
    S11_filt,S22_filt,S33_filt,S12_filt,S13_filt,S23_filt,S_filt,&
    S_S11_filt,S_S22_filt,S_S33_filt,S_S12_filt,S_S13_filt,S_S23_filt
integer :: jz

! use first LES grid point on ***w grid***
Deltay = dz
Delta_tilde = sqrt(dx*dy)
Delta_hat = 2._rprec*Delta_tilde
u_filt = (u(:,:,1)+u(:,:,2))/2._rprec
v_filt = (v(:,:,1)+v(:,:,2))/2._rprec
w_filt = w(:,:,2)
uu_filt = u_filt*u_filt
vv_filt = v_filt*v_filt
ww_filt = w_filt*w_filt
uv_filt = u_filt*v_filt
uw_filt = u_filt*w_filt
vw_filt = v_filt*w_filt

call test_filter(u_filt)
call test_filter(v_filt)
call test_filter(w_filt)
call test_filter(uu_filt)
call test_filter(vv_filt)
call test_filter(ww_filt)
call test_filter(uv_filt)
call test_filter(uw_filt)
call test_filter(vw_filt)

L11 = uu_filt - u_filt*u_filt
L22 = vv_filt - v_filt*v_filt
L33 = ww_filt - w_filt*w_filt
L12 = uv_filt - u_filt*v_filt
L13 = uw_filt - u_filt*w_filt
L23 = vw_filt - v_filt*w_filt

! compute dudz,dvdz using EQWM
! normally this is computed using FD
u_w = 0.5_rprec*(u(:,:,1)+u(:,:,2))
v_w = 0.5_rprec*(v(:,:,1)+v(:,:,2))
Ud = sqrt(u_w**2._rprec + v_w**2._rprec)
call retd_fit_calc2(Ud*Deltay/nu_molec,Deltap)
taud = (Deltap*nu_molec/Deltay)**2._rprec
lp = vonk*Deltap*(1._rprec-exp(-Deltap/25._rprec))
dudz_tot = taud/(2._rprec*nu_molec*lp**2._rprec)*(-1._rprec &
    + sqrt(1._rprec+4._rprec*lp**2._rprec))
uz = dudz_tot*u_w/Ud
vz = dudz_tot*v_w/Ud

! Calculate Sij for jz=2 on w-node
! uv-grid: dudx,dudy,dvdx,dvdy,dwdz
! w-grid: dudz,dvdz,dwdx,dwdy
ux = 0.5_rprec*(dudx(:,:,1)+dudx(:,:,2))
uy = 0.5_rprec*(dudy(:,:,1)+dudy(:,:,2))
vx = 0.5_rprec*(dvdx(:,:,1)+dvdx(:,:,2))
vy = 0.5_rprec*(dvdy(:,:,1)+dvdy(:,:,2))
! compute dwdz using continuity
wz = -(ux+vy)
S11 = ux
S12 = 0.5_rprec*(uy+vx)
S13 = 0.5_rprec*(uz + dwdx(:,:,2))
S22 = vy
S23 = 0.5_rprec*(vz + dwdy(:,:,2))
S33 = wz
S = sqrt(2._rprec*((S11**2.0 + S22**2.0 + S33**2.0)&
    + 2._rprec*(S12**2.0 + S13**2.0 + S23**2.0)))

! compute filtered Sij
S11_filt = S11
S12_filt = S12
S13_filt = S13
S22_filt = S22
S23_filt = S23
S33_filt = S33
call test_filter(S11_filt)
call test_filter(S12_filt)
call test_filter(S13_filt)
call test_filter(S22_filt)
call test_filter(S23_filt)
call test_filter(S33_filt)
S_filt = sqrt(2._rprec*((S11_filt**2.0 + S22_filt**2.0 &
    + S33_filt**2.0) + 2._rprec*(S12_filt**2.0 &
    + S13_filt**2.0 + S23_filt**2.0)))

! filter S*Sij
S_S11_filt = S*S11
S_S12_filt = S*S12
S_S13_filt = S*S13
S_S22_filt = S*S22
S_S23_filt = S*S23
S_S33_filt = S*S33
call test_filter(S_S11_filt)
call test_filter(S_S12_filt)
call test_filter(S_S13_filt)
call test_filter(S_S22_filt)
call test_filter(S_S23_filt)
call test_filter(S_S33_filt)

! compute Mij
M11 = 2._rprec*(Delta_tilde**2.0*S_S11_filt - Delta_hat**2.0*S_filt*S11_filt)
M12 = 2._rprec*(Delta_tilde**2.0*S_S12_filt - Delta_hat**2.0*S_filt*S12_filt)
M13 = 2._rprec*(Delta_tilde**2.0*S_S13_filt - Delta_hat**2.0*S_filt*S13_filt)
M22 = 2._rprec*(Delta_tilde**2.0*S_S22_filt - Delta_hat**2.0*S_filt*S22_filt)
M23 = 2._rprec*(Delta_tilde**2.0*S_S23_filt - Delta_hat**2.0*S_filt*S23_filt)
M33 = 2._rprec*(Delta_tilde**2.0*S_S33_filt - Delta_hat**2.0*S_filt*S33_filt)

Cs = sum(L11*M11+L22*M22+L33*M33 + 2._rprec*(L12*M12+L13*M13+L23*M23))/ &
     sum(M11*M11+M22*M22+M33*M33 + 2._rprec*(M12*M12+M13*M13+M23*M23))

nu_T = Cs*Delta_tilde**2.0*S
taudx_dyn = 2._rprec*(nu_molec + nu_T)*S13
taudy_dyn = 2._rprec*(nu_molec + nu_T)*S23
taudx_tot = sum(taudx_dyn(1:nx,1:ny))/nx/ny &
          + sum(u_w(1:nx,1:ny))/nx/ny*sum(w(1:nx,1:ny,2))/nx/ny &
          - sum(u_w(1:nx,1:ny)*w(1:nx,1:ny,2))/nx/ny
taudy_tot = sum(taudy_dyn(1:nx,1:ny))/nx/ny &
          + sum(v_w(1:nx,1:ny))/nx/ny*sum(w(1:nx,1:ny,2))/nx/ny &
          - sum(v_w(1:nx,1:ny)*w(1:nx,1:ny,2))/nx/ny

if (mod(jt_total,100)==0) then
    write(*,*) jt_total,taudx_tot,taudy_tot,&
               Cs,sum(taudx_dyn(1:nx,1:ny))/nx/ny,sum(taudy_dyn(1:nx,1:ny))/nx/ny
endif
if (jt_total==1000) then
    call write_dynamic_smag_model(u_w,v_w,w(:,:,2),u_filt,v_filt,w_filt,uz,vz,&
    S11,S12,S13,S22,S23,S33,L11,L12,L13,L22,L23,L33,M11,M12,M13,M22,M23,M33)
endif

end subroutine dynamic_smag_model2

!!*******************************************************************************
!subroutine rs_model(Deltay,taud,taudx,taudy,theta_d)
!!*******************************************************************************
!use param, only : nx, ny, nu_molec, dt, vonk, coord
!
!implicit none
!real(rprec), intent(in) :: Deltay
!real(rprec), dimension(nx,ny), intent(in) :: taud, theta_d
!real(rprec), dimension(nx,ny), intent(out) :: taudx, taudy
!real(rprec), dimension(nx,ny):: Deltap, lp, dudz_tot, upwp_eq, vpwp_eq, Tt, eps
!
!Deltap = Deltay*sqrt(taud)/nu_molec
!lp = vonk*Deltap*(1._rprec - exp(-Deltap/25._rprec))
!dudz_tot = taud/nu_molec/2._rprec/lp**2*(-1 + sqrt(1 + 4._rprec*lp**2))
!upwp_eq = -(lp*nu_molec/sqrt(taud))**2*dudz_tot**2*cos(theta_d)
!vpwp_eq = -(lp*nu_molec/sqrt(taud))**2*dudz_tot**2*sin(theta_d)
!Tt = vonk*Deltay/utau
!eps = min(abs(dt/Tt),1._rprec)
!upwp = eps*upwp_eq + (1._rprec-eps)*upwp
!vpwp = eps*vpwp_eq + (1._rprec-eps)*vpwp
!
!taudx = nu_molec*dudz_tot*cos(theta_d) - upwp
!taudy = nu_molec*dudz_tot*sin(theta_d) - vpwp
!
!!if (coord==0) then
!!    call rs_monitor(Deltay,Deltap,lp,dudz_tot,upwp_eq,vpwp_eq,Tt,&
!!        taudx,taudy,taud,theta_d)
!!endif
!
!end subroutine rs_model

!*******************************************************************************
subroutine qeqwm_write_checkpoint
!*******************************************************************************

use param, only : nx, ny, coord, nproc, write_endian
use types, only : rprec

implicit none

integer :: fid1, fid2

if (coord==0) then
    open(newunit=fid1, file='qeqwm_checkPoint_bot.bin', convert=write_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    write(fid1,rec=1) utau(1:nx,1:ny)
    write(fid1,rec=2) utx(1:nx,1:ny)
    write(fid1,rec=3) uty(1:nx,1:ny)
    write(fid1,rec=4) ssx(1:nx,1:ny)
    write(fid1,rec=5) ssy(1:nx,1:ny)
    write(fid1,rec=6) Ts(1:nx,1:ny)
    write(fid1,rec=7) Us(1:nx,1:ny)
!    write(fid1,rec=8) upwp(1:nx,1:ny)
!    write(fid1,rec=9) vpwp(1:nx,1:ny)
    close(fid1)
else
    open(newunit=fid2, file='qeqwm_checkPoint_top.bin', convert=write_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    write(fid2,rec=1) utau(1:nx,1:ny)
    write(fid2,rec=2) utx(1:nx,1:ny)
    write(fid2,rec=3) uty(1:nx,1:ny)
    write(fid2,rec=4) ssx(1:nx,1:ny)
    write(fid2,rec=5) ssy(1:nx,1:ny)
    write(fid2,rec=6) Ts(1:nx,1:ny)
    write(fid2,rec=7) Us(1:nx,1:ny)
!    write(fid2,rec=8) upwp(1:nx,1:ny)
!    write(fid2,rec=9) vpwp(1:nx,1:ny)
    close(fid2)
end if

end subroutine qeqwm_write_checkpoint

!*******************************************************************************
subroutine qeqwm_read_checkpoint
!*******************************************************************************

use param, only : nx, ny, coord, nproc, read_endian
use types, only : rprec

implicit none

integer :: fid1, fid2

if (coord==0) then
    open(newunit=fid1, file='qeqwm_checkPoint_bot.bin', convert=read_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    read(fid1,rec=1) utau(1:nx,1:ny)
    read(fid1,rec=2) utx(1:nx,1:ny)
    read(fid1,rec=3) uty(1:nx,1:ny)
    read(fid1,rec=4) ssx(1:nx,1:ny)
    read(fid1,rec=5) ssy(1:nx,1:ny)
    read(fid1,rec=6) Ts(1:nx,1:ny)
    read(fid1,rec=7) Us(1:nx,1:ny)
!    read(fid1,rec=8) upwp(1:nx,1:ny)
!    read(fid1,rec=9) vpwp(1:nx,1:ny)
    close(fid1)
else if (coord==nproc-1) then
    open(newunit=fid2, file='qeqwm_checkPoint_top.bin', convert=read_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    read(fid2,rec=1) utau(1:nx,1:ny)
    read(fid2,rec=2) utx(1:nx,1:ny)
    read(fid2,rec=3) uty(1:nx,1:ny)
    read(fid2,rec=4) ssx(1:nx,1:ny)
    read(fid2,rec=5) ssy(1:nx,1:ny)
    read(fid2,rec=6) Ts(1:nx,1:ny)
    read(fid2,rec=7) Us(1:nx,1:ny)
!    read(fid2,rec=8) upwp(1:nx,1:ny)
!    read(fid2,rec=9) vpwp(1:nx,1:ny)
    close(fid2)
end if

end subroutine qeqwm_read_checkpoint

!*******************************************************************************
subroutine write_dynamic_model(jz,u_filt,v_filt,w_filt,uw_filt,vw_filt,L13,L23,&
    F13,F23,F13_filt,F23_filt,F13_ufilt,F23_ufilt,M13,M23,alpha)
!*******************************************************************************
use sim_param, only : u,v,w
use param, only : ld, nx, ny, coord, write_endian
use types, only : rprec
use string_util, only : string_splice

implicit none

real(rprec), dimension(ld,ny) :: u_filt,v_filt,w_filt,uw_filt,vw_filt,L13,L23,&
    F13,F23,F13_ufilt,F23_ufilt,F13_filt,F23_filt,M13,M23,alpha
integer :: fid, jz
character(64) :: fname

if (coord==0) then
    call string_splice(fname,'output/dynamic_model_wmpt', jz, '.bin')
    open(newunit=fid, file=fname, convert=write_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    write(fid,rec=1) u(1:nx,1:ny,jz)
    write(fid,rec=2) v(1:nx,1:ny,jz)
    if (jz==1) then
        write(fid,rec=3) w(1:nx,1:ny,jz+1)*0.25_rprec
    else
        write(fid,rec=3) (w(1:nx,1:ny,jz)+w(1:nx,1:ny,jz+1))*0.5_rprec
    endif
    write(fid,rec=4) u_filt(1:nx,1:ny)
    write(fid,rec=5) v_filt(1:nx,1:ny)
    write(fid,rec=6) w_filt(1:nx,1:ny)
    write(fid,rec=7) uw_filt(1:nx,1:ny)
    write(fid,rec=8) vw_filt(1:nx,1:ny)
    write(fid,rec=9) L13(1:nx,1:ny)
    write(fid,rec=10) L23(1:nx,1:ny)
    write(fid,rec=11) F13(1:nx,1:ny)
    write(fid,rec=12) F23(1:nx,1:ny)
    write(fid,rec=13) F13_filt(1:nx,1:ny)
    write(fid,rec=14) F23_filt(1:nx,1:ny)
    write(fid,rec=15) F13_ufilt(1:nx,1:ny)
    write(fid,rec=16) F23_ufilt(1:nx,1:ny)
    write(fid,rec=17) M13(1:nx,1:ny)
    write(fid,rec=18) M23(1:nx,1:ny)
    write(fid,rec=19) alpha(1:nx,1:ny)
    close(fid)
end if

end subroutine write_dynamic_model

!*******************************************************************************
subroutine write_dynamic_smag_model(u,v,w,u_filt,v_filt,w_filt,dudz,dvdz,&
    S11,S12,S13,S22,S23,S33,L11,L12,L13,L22,L23,L33,M11,M12,M13,M22,M23,M33)
!*******************************************************************************
use param, only : ld, nx, ny, coord, write_endian
use types, only : rprec
use string_util, only : string_splice

implicit none

real(rprec), dimension(ld,ny) :: u,v,w,u_filt,v_filt,w_filt,dudz,dvdz,&
    S11,S12,S13,S22,S23,S33,L11,L12,L13,L22,L23,L33,M11,M12,M13,M22,M23,M33
integer :: fid

if (coord==0) then
    open(newunit=fid, file='dynamic_smag_model.bin', convert=write_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    write(fid,rec=1) u(1:nx,1:ny)
    write(fid,rec=2) v(1:nx,1:ny)
    write(fid,rec=3) w(1:nx,1:ny)
    write(fid,rec=4) u_filt(1:nx,1:ny)
    write(fid,rec=5) v_filt(1:nx,1:ny)
    write(fid,rec=6) w_filt(1:nx,1:ny)
    write(fid,rec=7) dudz(1:nx,1:ny)
    write(fid,rec=8) dvdz(1:nx,1:ny)
    write(fid,rec=9) S11(1:nx,1:ny)
    write(fid,rec=10) S12(1:nx,1:ny)
    write(fid,rec=11) S13(1:nx,1:ny)
    write(fid,rec=12) S22(1:nx,1:ny)
    write(fid,rec=13) S23(1:nx,1:ny)
    write(fid,rec=14) S33(1:nx,1:ny)
    write(fid,rec=15) L11(1:nx,1:ny)
    write(fid,rec=16) L12(1:nx,1:ny)
    write(fid,rec=17) L13(1:nx,1:ny)
    write(fid,rec=18) L22(1:nx,1:ny)
    write(fid,rec=19) L23(1:nx,1:ny)
    write(fid,rec=20) L33(1:nx,1:ny)
    write(fid,rec=21) M11(1:nx,1:ny)
    write(fid,rec=22) M12(1:nx,1:ny)
    write(fid,rec=23) M13(1:nx,1:ny)
    write(fid,rec=24) M22(1:nx,1:ny)
    write(fid,rec=25) M23(1:nx,1:ny)
    write(fid,rec=26) M33(1:nx,1:ny)
    close(fid)
end if

end subroutine write_dynamic_smag_model

!!*******************************************************************************
!subroutine rs_monitor(Deltay,Deltap,lp,dudz_tot,upwp_eq,vpwp_eq,Tt,&
!    taudx,taudy,taud,theta_d)
!!*******************************************************************************
!!
!! This subroutine monitors the temporal evolution of neqwm quantities
!!
!use param, only : jt_total, total_time, path, nx, ny, dt
!
!implicit none
!
!integer :: fid, i, j
!character*50 :: fname
!
!real(rprec), intent(in) :: Deltay
!real(rprec), dimension(nx,ny), intent(in) :: taud, theta_d, taudx, taudy,&
!     Deltap, lp, dudz_tot, upwp_eq, vpwp_eq, Tt
!
!i = int(nx/2._rprec)
!j = int(ny/2._rprec)
!
!fname = path // 'output/rs_monitor.dat'
!open(newunit=fid, file=fname, status='unknown', position='append')
!write(fid,*) jt_total,&
!    total_time,&
!    dt,&
!    Deltay,&
!    Deltap(i,j),&
!    lp(i,j),&
!    dudz_tot(i,j),&
!    upwp_eq(i,j),&
!    vpwp_eq(i,j),&
!    Tt(i,j),&
!    upwp(i,j),&
!    vpwp(i,j),&
!    taudx(i,j),&
!    taudy(i,j),&
!    taud(i,j),&
!    theta_d(i,j),&
!    utx(i,j),&
!    uty(i,j)
!close(fid)
!
!end subroutine rs_monitor

end module qeqwm 
