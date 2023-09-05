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
use wm_param

private
public qeqwm_initialize, qeqwm_finalize, lagrangian_rewm, eqwm_compute, &
    qeqwm_write_checkpoint, qeqwm_read_checkpoint, &
!    utau, utx, uty, ssx, ssy, Ts, Us, vtx, vty, utau_filt, &
    velocity_fit, retd_fit

!!Lagrangian relaxation wall model variables
!! friction velocities
!real(rprec), dimension(:,:), allocatable :: utau, utx, uty
!! wall stress angle components
!real(rprec), dimension(:,:), allocatable :: ssx, ssy
!! relaxation time scale
!real(rprec), dimension(:,:), allocatable :: Ts
!! velocity at first grid point computed from assumed profile and utau
!real(rprec), dimension(:,:), allocatable :: Us
!! advecting velocities in lagrangian rewm
!real(rprec), dimension(:,:), allocatable :: vtx, vty
!!! reynolds stress evaluated at Delta
!!real(rprec), dimension(:,:), allocatable :: upwp, vpwp
!real(rprec), dimension(:,:), allocatable :: utau_filt

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
allocate(utau_filt(nx,ny))

utau = sqrt(twxbar)
utx = sqrt(twxbar)
uty = 0._rprec
ssx = 1._rprec
ssy = 0._rprec
Us = ubar
Ts = Us*Deltay/utau**2.0
vtx = 0._rprec
vty = 0._rprec
utau_filt = utau

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
!subroutine lagrangian_rewm(Deltay,Ud,theta_d,dpdxbar1,dpdybar1,&
!    dpdxbar2,dpdybar2,twxbar,twybar)
subroutine lagrangian_rewm()
!*******************************************************************************
use param, only : ld, nx, ny, nu_molec, dt, dt_f, jt, coord, jt_total, &
    total_time, vonk
use grid_m
use param, only : path
use param, only : inflow_type, fringe_region_end, fringe_region_len, L_x
!use eqwm_dyn, only : dynamic_smag_model

implicit none
!real(rprec), dimension(nx,ny), intent(in) :: Ud, theta_d
!real(rprec), dimension(nx,ny), intent(in) :: dpdxbar1, dpdybar1,&
!    dpdxbar2, dpdybar2
!real(rprec), intent(in) :: Deltay
!real(rprec), dimension(nx,ny), intent(out) :: twxbar, twybar
real(rprec), dimension(nx,ny) :: sxp, syp, Deltap, fu, reds, delta_star, retdu,&
        taud, theta_delta, rhsx, rhsy, rhsx0, rhsy0, dpdxu, &
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
real(rprec) :: xf1

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
vtx = (1.0-delta_star-theta_delta)*utx*fu
vty = (1.0-delta_star-theta_delta)*uty*fu
!redu = Ud*Deltay/nu_molec
!dpdxu = dpdxbar2*cos(theta_d) + dpdybar2*sin(theta_d)
!psi_p = dpdxu*Deltay**3.0/nu_molec**2.0
!call retd_pres_fit(redelta,psi_p,retdu)
call retd_fit(redelta,retdu)
!taud = (retdu*nu_molec/Deltay)**2
!taudx = taud*cos(theta_d)
!taudy = taud*sin(theta_d)
twx_eq = (retdu*nu_molec/Deltay)**2.0*cos(theta_d) + 10._rprec**-10._rprec
twy_eq = (retdu*nu_molec/Deltay)**2.0*sin(theta_d)

rhsx = utx + dt*(utau*delta_star*(ssx-sxp)/dt_f &
    + (twx_eq/utau - utx)/Ts)
rhsy = uty + dt*(utau*delta_star*(ssy-syp)/dt_f &
    + (twy_eq/utau - uty)/Ts)

xf1 = modulo(fringe_region_end-fringe_region_len,1._rprec)*L_x

!do i = 1,100
!do j = 1,ny
!!    rhsx0(i,j) = twx_eq(i,j)/(twx_eq(i,j)**2+twy_eq(i,j)**2)**0.25 
!!    rhsy0(i,j) = twy_eq(i,j)/(twx_eq(i,j)**2+twy_eq(i,j)**2)**0.25
!    rhsx0(i,j) = rhsx(i,j)
!    rhsy0(i,j) = rhsy(i,j)
!enddo
!enddo

do i = 1,nx
do j = 1,ny
    x0 = x(i) - vtx(i,j)*dt
    y0 = y(j) - vty(i,j)*dt
!    if (inflow_type .ne. 0 .and. &
!        modulo(x0,L_x) >= xf1 .and. &
!        modulo(x0,L_x) <= xf1 + fringe_region_len*L_x) then
!        rhsx0(i,j) = rhsx(i,j)
!        rhsy0(i,j) = rhsy(i,j)
!!        rhsx0(i,j) = twx_eq(i,j)/(twx_eq(i,j)**2+twy_eq(i,j)**2)**0.25 
!!        rhsy0(i,j) = twy_eq(i,j)/(twx_eq(i,j)**2+twy_eq(i,j)**2)**0.25
!!        rhsx0(i,j) = twx_eq(i,j)/utau(i,j) + 10._rprec**-10.0
!!        rhsy0(i,j) = twy_eq(i,j)/utau(i,j) + 10._rprec**-10.0
!!        write(*,*) jt_total,twx_eq(i,j),twy_eq(i,j),rhsx0(i,j)
!    else
        call bilinear_interp(rhsx,x0,y0,rhsx0(i,j))
        call bilinear_interp(rhsy,x0,y0,rhsy0(i,j))
!    endif
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

!write(*,*) 'larte: ',jt_total, sum(twxbar)/nx/ny, sum(twybar)/nx/ny

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
!subroutine eqwm_compute(Deltay,Ud,theta_d,dpdxbar1,dpdybar1,&
!    dpdxbar2,dpdybar2,twxbar,twybar)
subroutine eqwm_compute()
!*******************************************************************************
use param, only : nx, ny, nu_molec
use param, only : jt_total

implicit none
!real(rprec), dimension(nx,ny), intent(in) :: Ud, theta_d
!real(rprec), dimension(nx,ny), intent(in) :: dpdxbar1, dpdybar1,&
!    dpdxbar2, dpdybar2
!real(rprec), intent(in) :: Deltay
!real(rprec), dimension(nx,ny), intent(out) :: twxbar, twybar
real(rprec), dimension(nx,ny) :: Deltap, fu, retdu, taud, taudx, taudy, &
    dpdxu

Deltap = utau*Deltay/nu_molec ! Delta_y in plus units (from previous tstep)
call velocity_fit(Deltap,fu)
Ts = fu*Deltay/utau
!redu = Ud*Deltay/nu_molec
!dpdxu = dpdxbar2*cos(theta_d) + dpdybar2*sin(theta_d)
!psi_p = dpdxu*Deltay**3.0/nu_molec**2.0
!call retd_pres_fit(redelta,psi_p,retdu)
call retd_fit(redelta,retdu)
!taud = (retdu*nu_molec/Deltay)**2
!taudx = taud*cos(theta_d)
!taudy = taud*sin(theta_d)
twx_eq = (retdu*nu_molec/Deltay)**2.0*cos(theta_d)
twy_eq = (retdu*nu_molec/Deltay)**2.0*sin(theta_d)

twxbar = twx_eq
twybar = twy_eq

utau = (twxbar**2+twybar**2)**0.25
utx = twxbar/utau
uty = twybar/utau

!write(*,*) 'eqwm : ',jt_total, sum(twxbar)/nx/ny, sum(twybar)/nx/ny
!write(*,*) 'mts: ',jt_total, twxbar(1,1), redelta(1,1), psi_p(1,1)
!write(*,*) 'mts        : ',jt_total, dplesdx(1,1),dplesdy(1,1), psi_p(1,1)

end subroutine eqwm_compute

!*******************************************************************************
subroutine retd_pres_fit(redelta_in,psi_in,retd_out)
!*******************************************************************************
use param, only : nx, ny, nu_molec

implicit none
real(rprec), dimension(nx,ny), intent(in) :: redelta_in, psi_in
real(rprec), dimension(nx,ny), intent(out) :: retd_out
real(rprec), dimension(nx,ny) :: retd_eq
real(rprec) :: retd_min, pp, redelta_min
integer :: jx, jy

call retd_fit(redelta_in,retd_eq)
do jx = 1,nx
do jy = 1,ny
if (psi_in(jx,jy) < 0) then
    retd_min = 1.5*(-psi_in(jx,jy))**0.39*&
        (1.0+(1000.0/-psi_in(jx,jy))**2.0)**-0.055
    pp = 2.5 - 0.6*(1.0 + tanh(2.0*(log10(-psi_in(jx,jy))-6.0)))
    retd_out(jx,jy) = (retd_min**pp + retd_eq(jx,jy)**pp)**(1.0/pp)
else
    redelta_min = 2.5*psi_in(jx,jy)**0.54*(1.0+(30.0/psi_in(jx,jy))**0.5)**-0.88
    if (redelta_in(jx,jy) > redelta_min) then
        retd_out(jx,jy) = retd_eq(jx,jy)*&
            (1.0-(1.0+log(redelta_in(jx,jy)/redelta_min))**-1.9)
    else
        retd_out(jx,jy) = 0._rprec
    endif
endif
enddo
enddo

end subroutine retd_pres_fit

!*******************************************************************************
subroutine retd_fit(redelta_in,retd_out)
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

end subroutine retd_fit

!!*******************************************************************************
!subroutine retd_fit_calc(redelta_in,retd_out)
!!*******************************************************************************
!use param, only : nx, ny
!
!implicit none
!real(rprec), dimension(nx,ny), intent(in) :: redelta_in
!real(rprec), dimension(nx,ny), intent(out) :: retd_out
!real(rprec) :: kappa3
!real(rprec), dimension(nx,ny) :: kappa4, beta1, beta2
!
!kappa3 = 0.005_rprec
!beta1 = (1._rprec + 0.155_rprec*redelta_in**-0.03_rprec)**-1
!beta2 = 1.7 - (1._rprec + 36._rprec*redelta_in**-0.75_rprec)**-1
!kappa4 = kappa3**(beta1 - 0.5_rprec)
!
!retd_out = kappa4*redelta_in**beta1*                            &
!    (1._rprec + (kappa3*redelta_in)**-beta2)**((beta1-0.5_rprec)/beta2)
!
!end subroutine retd_fit_calc
!
!!*******************************************************************************
!subroutine retd_fit_calc2(redelta_in,retd_out)
!!*******************************************************************************
!use param, only : ld, ny
!
!implicit none
!real(rprec), dimension(ld,ny), intent(in) :: redelta_in
!real(rprec), dimension(ld,ny), intent(out) :: retd_out
!real(rprec) :: kappa3
!real(rprec), dimension(ld,ny) :: kappa4, beta1, beta2
!
!kappa3 = 0.005_rprec
!beta1 = (1._rprec + 0.155_rprec*redelta_in**-0.03_rprec)**-1
!beta2 = 1.7 - (1._rprec + 36._rprec*redelta_in**-0.75_rprec)**-1
!kappa4 = kappa3**(beta1 - 0.5_rprec)
!
!retd_out = kappa4*redelta_in**beta1*                            &
!    (1._rprec + (kappa3*redelta_in)**-beta2)**((beta1-0.5_rprec)/beta2)
!
!end subroutine retd_fit_calc2
!
!!*******************************************************************************
!subroutine retd_pres_calc(redelta_in,psi_p,retd_out)
!!*******************************************************************************
!use param, only : nx, ny, nu_molec
!
!implicit none
!real(rprec), dimension(nx,ny), intent(in) :: redelta_in, psi_p
!real(rprec), dimension(nx,ny), intent(out) :: retd_out
!real(rprec), dimension(nx,ny) :: retd_eq, retd_min, pp
!
!call retd_fit_calc(redelta_in,retd_eq)
!retd_min = 1.5*psi_p**0.39*(1.0+(1000.0/psi_p)**2.0)**-0.055
!pp = 2.5 - 0.6*(1.0 + tanh(2.0*log10(psi_p)**-6.0))
!retd_out = (retd_min**pp + retd_eq**pp)**(1.0/pp)
!
!end subroutine retd_pres_calc

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

end module qeqwm 
