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
module mts_wm
!*******************************************************************************
use types, only : rprec
use wm_param

implicit none

private
public mts_initialize, mts_finalize, mts_wallstress_calc,&
    mts_monitor, mts_monitor_plane, mts_write_checkpoint,&
    mts_read_checkpoint, get_wallmodelstress
!    twxbar, twybar, twxpp, twypp, twxpp_turb, twypp_turb,&
!    ubar, vbar, ls, Deltay

!! wall stress from quasi-equilibrium and non-equilibrium models
!real(rprec), dimension(:,:), allocatable :: twxbar, twybar, twxpp, twypp,&
!    twxpp_turb,twypp_turb
!! velocities at first grid point (non-filtered)
!real(rprec), dimension(:,:), allocatable :: uinst, vinst, Ud
!! velocity angle(measured from x-axis): d = Deltay
!real(rprec), dimension(:,:), allocatable :: theta_d
!! pressure gradients at first grid point
!real(rprec), dimension(:,:), allocatable :: dplesdx, dplesdy, dpdxbar1, dpdybar1,&
!    dpdxbar2, dpdybar2, dpdxpp, dpdypp, dpdxpp_m, dpdypp_m, dpdxpp_mm, dpdypp_mm
!! wall model height
!real(rprec) :: Deltay
!! temporary files for debug_wm
!real(rprec), dimension(:,:,:), allocatable :: utemp2,vtemp2,dpdxtemp,dpdytemp,&
!    utxtemp, utytemp
!! time scale for filtering the pressure
!real(rprec), dimension(:,:), allocatable :: Tnu
!! filtered LES velocity
!real(rprec), dimension(:,:), allocatable :: ubar, vbar
!real(rprec), dimension(:,:), allocatable :: uinf, vinf
!real(rprec), dimension(:,:), allocatable :: uinf_lam, vinf_lam
!real(rprec), dimension(:,:), allocatable :: Tdelta
!real(rprec), dimension(:,:), allocatable :: dpdxbar3,dpdybar3
!real(rprec), dimension(:,:), allocatable :: dpdxppdelta,dpdyppdelta
!real(rprec) :: ls=12._rprec ! thickness of Stokes layer in inner units
!
!real(rprec) :: PI=4*atan(1._rprec)

contains

!*******************************************************************************
subroutine mts_initialize
!*******************************************************************************
use param, only : nx, ny, nz, dz, coord, mean_p_force_x, mean_p_force_y, &
    nu_molec, wmpt
use sim_param, only : u, v
use qeqwm, only : qeqwm_initialize, retd_fit
use neqwm, only : neqwm_initialize

implicit none
real(rprec), dimension(nx,ny) :: retd

allocate(twxbar(nx,ny))
allocate(twybar(nx,ny))
allocate(twxpp(nx,ny))
allocate(twypp(nx,ny))
allocate(twxp(nx,ny))
allocate(twyp(nx,ny))
allocate(uinst(nx,ny))
allocate(vinst(nx,ny))
allocate(Ud(nx,ny))
allocate(dplesdx(nx,ny))
allocate(dplesdy(nx,ny))
allocate(dpdxbar1(nx,ny))
allocate(dpdybar1(nx,ny))
allocate(dpdxbar2(nx,ny))
allocate(dpdybar2(nx,ny))
allocate(dpdxpp(nx,ny))
allocate(dpdypp(nx,ny))
allocate(dpdxpp_m(nx,ny))
allocate(dpdypp_m(nx,ny))
allocate(dpdxpp_mm(nx,ny))
allocate(dpdypp_mm(nx,ny))
allocate(theta_d(nx,ny))
allocate(Tnu(nx,ny))
allocate(ubar(nx,ny))
allocate(vbar(nx,ny))
allocate(uinf(nx,ny))
allocate(vinf(nx,ny))
allocate(uinf_lam(nx,ny))
allocate(vinf_lam(nx,ny))
allocate(Tdelta(nx,ny))
allocate(dpdxbar3(nx,ny))
allocate(dpdybar3(nx,ny))
allocate(dpdxppdelta(nx,ny))
allocate(dpdyppdelta(nx,ny))
allocate(dpdx_fit(nx,ny))
allocate(dpdy_fit(nx,ny))
allocate(twx_eq(nx,ny))
allocate(twy_eq(nx,ny))
allocate(redelta(nx,ny))
allocate(psi_p(nx,ny))

allocate(utemp2(nx,ny,1000))
allocate(vtemp2(nx,ny,1000))
allocate(dpdxtemp(nx,ny,1000))
allocate(dpdytemp(nx,ny,1000))
allocate(utxtemp(nx,ny,1000))
allocate(utytemp(nx,ny,1000))

! provide initial values
if (coord==0) then
    uinst = u(1:nx,1:ny,wmpt)
    vinst = v(1:nx,1:ny,wmpt)
else
    uinst = u(1:nx,1:ny,nz-wmpt)
    vinst = v(1:nx,1:ny,nz-wmpt)
endif
ubar = uinst
vbar = vinst
Ud = sqrt(uinst**2+vinst**2)
redelta = Ud*Deltay/nu_molec
call retd_fit(redelta,retd)
twx_eq = (retd*nu_molec/Deltay)**2.0
twy_eq = 0._rprec
twxbar = twx_eq
twybar = twy_eq
twxpp = 0._rprec
twypp = 0._rprec
dplesdx = - mean_p_force_x
dplesdy = - mean_p_force_y
dpdxbar1 = dplesdx
dpdybar1 = dplesdy
dpdxbar2 = dplesdx
dpdybar2 = dplesdy
dpdx_fit = dplesdx
dpdy_fit = dplesdy
dpdxpp = 0._rprec
dpdypp = 0._rprec
dpdxpp_m = 0._rprec
dpdypp_m = 0._rprec
dpdxpp_mm = 0._rprec
dpdypp_mm = 0._rprec
!Deltay = dz/2._rprec
Deltay = (wmpt-0.5_rprec)*dz
Tnu = 100._rprec*nu_molec/twxbar
uinf = 0._rprec
vinf = 0._rprec
uinf_lam = 0._rprec
vinf_lam = 0._rprec
Tdelta = Tnu
dpdxbar3 = dplesdx
dpdybar3 = dplesdy
dpdxppdelta = 0._rprec
dpdyppdelta = 0._rprec
twxp = 0._rprec
twyp = 0._rprec
psi_p = dplesdx*Deltay**3.0/nu_molec**2.0

call qeqwm_initialize
call neqwm_initialize

end subroutine mts_initialize

!*******************************************************************************
subroutine mts_finalize
!*******************************************************************************

use qeqwm, only : qeqwm_finalize
use neqwm, only : neqwm_finalize

implicit none

deallocate(twxbar)
deallocate(twybar)
deallocate(twxpp)
deallocate(twypp)
deallocate(uinst)
deallocate(vinst)
deallocate(Ud)
deallocate(dplesdx)
deallocate(dplesdy)
deallocate(dpdxbar1)
deallocate(dpdybar1)
deallocate(dpdxbar2)
deallocate(dpdybar2)
deallocate(dpdxpp)
deallocate(dpdypp)
deallocate(dpdxpp_m)
deallocate(dpdypp_m)
deallocate(dpdxpp_mm)
deallocate(dpdypp_mm)
deallocate(theta_d)
deallocate(Tnu)
deallocate(ubar)
deallocate(vbar)
deallocate(uinf)
deallocate(vinf)
deallocate(uinf_lam)
deallocate(vinf_lam)
deallocate(Tdelta)
deallocate(dpdxbar3)
deallocate(dpdybar3)
deallocate(dpdxppdelta)
deallocate(dpdyppdelta)
deallocate(dpdx_fit)
deallocate(dpdy_fit)
deallocate(twx_eq)
deallocate(twy_eq)

deallocate(utemp2)
deallocate(vtemp2)
deallocate(dpdxtemp)
deallocate(dpdytemp)
deallocate(utxtemp)
deallocate(utytemp)

call qeqwm_finalize
call neqwm_finalize

end subroutine mts_finalize

!*******************************************************************************
subroutine mts_wallstress_calc
!*******************************************************************************

use param, only : nx, ny, nz, dz, ld, nz, coord, dt, vonk, nu_molec, jt, &
    mean_p_force_x, mean_p_force_y, nu_molec, wmpt, path
use param, only : jt_total, total_time
use param, only : qeq_case, lamNEQ_flag, turbNEQ_flag, velocity_correction_flag
use sim_param, only : u, v, w, dpdx, dpdy, txz, tyz, dudz, dvdz
use qeqwm, only : lagrangian_rewm, eqwm_compute, retd_fit
!, utau, utx, uty, Ts, utau_filt
use neqwm, only : neq_laminar_calc, neq_turb_calc
use test_filtermodule

!use param, only : dx, dy, dz, lbc_mom, ubc_mom, coord, nproc, jt_total, nu_molec, &
!    vonk, ld, nx, ny, nz, dt, mean_p_force_x, mean_p_force_y, jt, wmpt, &
!    path, total_time, write_endian, read_endian
!use sim_param, only : u, v, w, p, txz, tyz, dudz, dvdz
!use sim_param, only : dpdx, dpdy, dudx, dvdx, dwdx, dudy, dvdy, dwdy
!use test_filtermodule
!!use fft
!!use emul_complex, only : OPERATOR(.MULI.)
!use string_util
use derivatives

implicit none

real(rprec), dimension(nx,ny) :: eps1, eps2, eps3, eps4

integer :: i, j
real(rprec), dimension(nx,ny) :: lp, dudz_tot, theta_w, retd2, utau2
!real(rprec), dimension(nx,ny) :: uinst_m, vinst_m, Ules_m, Vles_m, advx, advy
!real(rprec), dimension(nx,ny) :: Usbar
!real(rprec), dimension(ld,ny) :: dummy1, dummy2, dummy3, dummy4
!real(rprec), dimension(ld,ny) :: ududx, vdvdx, wdwdx
!real(rprec), dimension(ld,ny) :: ududy, vdvdy, wdwdy
real(rprec), dimension(ld,ny) :: utemp, vtemp, wtemp
real(rprec), dimension(ld,ny) :: ke_dealiased, dkedx, dkedy
!real(rprec), dimension(ld,ny) :: ke_dealiased, ke1, ke2, ke3, dkedx, dkedy
!real(rprec), dimension(ld,ny) :: dkedx1, dkedx2, dkedx3, dkedy1, dkedy2, dkedy3
!real(rprec), dimension(ld,ny) :: dpdx_ave, dpdy_ave
integer :: fid
character*50 :: fname
!real(rprec), dimension(nx,ny) :: uinf_turb,vinf_turb,fu
real(rprec), dimension(nx,ny) :: upp, utaupp
real(rprec) :: Tlam
real(rprec), dimension(nx,ny) :: taudx_eq, taudy_eq, retd, utau_eq, fu_eq

!Usbar = sqrt(Ules**2+Vles**2)
!
!Ts = abs(Usbar*Deltay/2._rprec/tauws)
!Tn = abs(Us*Deltay/tauws*(1._rprec - delta_star))

!eps = dt/maxval(abs(Ts))
!
!! store velocity at previous time step to calculate time derivative
!uinst_m = uinst
!vinst_m = vinst
!Ules_m = Ules
!Vles_m = Vles


! collect LES variables
if (coord==0) then
    uinst(1:nx,1:ny) = u(1:nx,1:ny,wmpt)
    vinst(1:nx,1:ny) = v(1:nx,1:ny,wmpt)
    utemp = u(:,:,wmpt)
    vtemp = v(:,:,wmpt)
    if (wmpt==1) then
        wtemp = 0.25_rprec*w(:,:,2)
    else
        wtemp = 0.5_rprec*(w(:,:,wmpt)+w(:,:,wmpt+1))
    endif
    if (jt .ne. 1) then
        dplesdx = dpdx(1:nx,1:ny,wmpt)
        dplesdy = dpdy(1:nx,1:ny,wmpt)
    endif
else
    uinst(1:nx,1:ny) = u(1:nx,1:ny,nz-wmpt)
    vinst(1:nx,1:ny) = v(1:nx,1:ny,nz-wmpt)
    utemp = u(:,:,nz-wmpt)
    vtemp = v(:,:,nz-wmpt)
    if (wmpt==1) then
        wtemp = 0.25_rprec*w(:,:,nz-wmpt)
    else
        wtemp = 0.5_rprec*(w(:,:,nz-wmpt)+w(:,:,nz-wmpt-1))
    endif
    if (jt .ne. 1) then
        dplesdx = dpdx(1:nx,1:ny,nz-wmpt)
        dplesdy = dpdy(1:nx,1:ny,nz-wmpt)
    endif
end if
! p is set to zero after restarting
if (jt .ne. 1) then
    ! calculate real pressure gradient 
    ! (dealiased and gradient computed spectrally
    ! for accurate estimate of PG)
    call ke_dealiased_calc(ke_dealiased,utemp,vtemp,wtemp)
    call ddxy_plane(ke_dealiased, dkedx, dkedy)
    dplesdx(1:nx,1:ny) = dplesdx(1:nx,1:ny) - dkedx(1:nx,1:ny) &
        - mean_p_force_x
    dplesdy(1:nx,1:ny) = dplesdy(1:nx,1:ny) - dkedy(1:nx,1:ny) &
        - mean_p_force_y
endif

! filtered friction velocity for time scales
!utau_filt = eps2*sqrt(twxbar**2+twybar**2) + (1._rprec-eps2)*utau_filt
Tnu = ls**2._rprec*nu_molec/utau**2._rprec
Tdelta = Tnu + Deltay/0.4_rprec/utau
Tlam = 8._rprec*Deltay**2.0/4._rprec/nu_molec

eps1 = min(abs(dt/Tnu),1._rprec)
eps2 = min(abs(dt/Ts/3._rprec),1._rprec)
eps3 = min(abs(dt/Ts),1._rprec)
eps4 = min(abs(dt/Tdelta),1._rprec)

! filtered LES velocity
!ubar = eps3*uinst + (1._rprec-eps3)*ubar
!vbar = eps3*vinst + (1._rprec-eps3)*vbar
!ubar = sum(uinst)/nx/ny
!vbar = sum(vinst)/nx/ny
!call test_filter(utemp)
!call test_filter(vtemp)
ubar = utemp(1:nx,1:ny)
vbar = vtemp(1:nx,1:ny)

!! temporal filter applied to pressure gradient
! 1 should be faster (shorter time scale) than 2
dpdxbar1(1:nx,1:ny) = eps1(1:nx,1:ny)*dplesdx(1:nx,1:ny) &
    + (1._rprec-eps1(1:nx,1:ny))*dpdxbar1(1:nx,1:ny)
dpdybar1(1:nx,1:ny) = eps1(1:nx,1:ny)*dplesdy(1:nx,1:ny) &
    + (1._rprec-eps1(1:nx,1:ny))*dpdybar1(1:nx,1:ny)
dpdxbar2(1:nx,1:ny) = eps2(1:nx,1:ny)*dplesdx(1:nx,1:ny) &
    + (1._rprec-eps2(1:nx,1:ny))*dpdxbar2(1:nx,1:ny)
dpdybar2(1:nx,1:ny) = eps2(1:nx,1:ny)*dplesdy(1:nx,1:ny) &
    + (1._rprec-eps2(1:nx,1:ny))*dpdybar2(1:nx,1:ny)
dpdxbar3 = eps4*dplesdx + (1._rprec-eps4)*dpdxbar3
dpdybar3 = eps4*dplesdy + (1._rprec-eps4)*dpdybar3

dpdxppdelta = dplesdx - dpdxbar3
dpdyppdelta = dplesdy - dpdybar3

if (velocity_correction_flag) then
! subtracting Stokes layer non-equilibrium velocity
!uinf = -eps4*Tdelta*dpdxppdelta + (1._rprec-eps4)*uinf
!vinf = -eps4*Tdelta*dpdyppdelta + (1._rprec-eps4)*vinf
!uinf = -eps1*Tnu*(dplesdx-dpdxbar1) + (1._rprec-eps1)*uinf
!vinf = -eps1*Tnu*(dplesdy-dpdybar1) + (1._rprec-eps1)*vinf
uinf = -dt*(dplesdx-dpdxbar1) + (1._rprec-eps4)*uinf
vinf = -dt*(dplesdy-dpdybar1) + (1._rprec-eps4)*vinf
ubar = ubar - uinf
vbar = vbar - vinf
endif

! non-equilibrium PG
dpdxpp_mm = dpdxpp_m
dpdypp_mm = dpdypp_m
dpdxpp_m = dpdxpp
dpdypp_m = dpdypp
dpdxpp(1:nx,1:ny) = dplesdx(1:nx,1:ny) - dpdxbar1(1:nx,1:ny)
dpdypp(1:nx,1:ny) = dplesdy(1:nx,1:ny) - dpdybar1(1:nx,1:ny)
!dpdxpp = dpdxppdelta + uinf/Tdelta
!dpdypp = dpdyppdelta + vinf/Tdelta

!dpdxbar1 = 0._rprec
!dpdybar1 = 0._rprec
!dpdxbar2 = 0._rprec
!dpdybar2 = 0._rprec

!Ud = sqrt(uinst**2+vinst**2)
!theta_d = atan2(vinst,uinst)
Ud = sqrt(ubar**2 + vbar**2)
theta_d = atan2(vbar,ubar)

! inputs to eqwm fit
redelta = Ud*Deltay/nu_molec
if (velocity_correction_flag) then
dpdx_fit = dpdxbar1
dpdy_fit = dpdybar1
else
dpdx_fit = dplesdx
dpdy_fit = dplesdy
endif
psi_p = (dpdx_fit*cos(theta_d)+dpdy_fit*sin(theta_d))*Deltay**3.0/nu_molec**2.0

! ******** calling wall models **********

!if (coord==0) then
!    write(*,*) 'before: ', jt_total, utau(1,1), utx(1,1), uty(1,1)
!endif

! Quasi-equilibrium
select case (qeq_case)
    case (0)
        call eqwm_compute()
    case (1)
        call lagrangian_rewm()
        !call eqwm_compute()
end select

!if (coord==0) then
!    write(*,*) 'after : ', jt_total, utau(1,1), utx(1,1), uty(1,1)
!endif

! Laminar non-equilibrium
if (lamNEQ_flag) then
    call neq_laminar_calc()
else
    twxpp = 0._rprec
    twypp = 0._rprec
endif

! Turbulent non-equilibrium
if (turbNEQ_flag) then
    call neq_turb_calc()
else
    twxp = 0._rprec
    twyp = 0._rprec
endif

! velocity gradient at first grid point computed without PG in fit
! for simplicity, friction velocity computed using LES velocity at wmpt
call retd_fit(redelta,retd2)
utau2 = nu_molec*retd2/Deltay
lp = vonk*dz/2._rprec*utau2/nu_molec                          &
    *(1 - exp(-dz/2._rprec*utau2/nu_molec/25._rprec))
! total derivative at first grid point
dudz_tot = utau2**2/nu_molec/2._rprec/lp**2*(-1 + sqrt(1 + 4._rprec*lp**2))

theta_w = atan2(twybar,twxbar)

!if (coord==0 .and. mod(jt,1000) .eq. 0) then
!    call write_velocity_decomp(u1,u2,u3,u4,v1,v2,v3,v4)
!endif

!if (coord==0) then
!    call write_velocity_pt1_pt2
!endif

!i = int(nx/2._rprec)
!j = int(ny/2._rprec)
!if (coord==0 .and. mod(jt,1) .eq. 0) then
!    write(*,*) jt_total, uinf_turb(i,j), vinf_turb(i,j), twxpp_turb(i,j), twypp_turb(i,j)
!endif


if (coord==0) then
    txz(1:nx,1:ny,1) = - twxbar - twxpp - twxp
    tyz(1:nx,1:ny,1) = - twybar - twypp - twyp
    dudz(1:nx,1:ny,1) = dudz_tot*cos(theta_w)
    dvdz(1:nx,1:ny,1) = dudz_tot*sin(theta_w)
else
    txz(1:nx,1:ny,nz) = twxbar + twxpp + twxp
    tyz(1:nx,1:ny,nz) = twybar + twypp + twyp
    dudz(1:nx,1:ny,nz) = -dudz_tot*cos(theta_w)
    dvdz(1:nx,1:ny,nz) = -dudz_tot*sin(theta_w)
end if

!if (coord==0) then
!i = int(nx/2._rprec)
!j = int(ny/2._rprec)
!fname = path // 'output/mts_monitor_point_extra.dat'
!open(newunit=fid, file=fname, status='unknown', position='append')
!write(fid,*) jt_total,&
!    total_time,&
!    utemp(i,j),&
!    vtemp(i,j),&
!    uinf_lam(i,j),&
!    vinf_lam(i,j),&
!    uinf_turb(i,j),&
!    vinf_turb(i,j),&
!    txz(i,j,1),&
!    tyz(i,j,1),&
!    twxpp_turb(i,j),&
!    twypp_turb(i,j),&
!    utau_filt(i,j),&
!    Tnu(i,j)
!close(fid)
!endif

!if (coord==0) then
!    write(*,*) jt, 'after:',txz(nx/2,ny/2,1)
!endif
!if (coord==0) then
!    call mts_monitor()
!endif
!if (coord==0) then
!    call mts_ws_plane
!    call mts_test_point
!endif

end subroutine mts_wallstress_calc

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
subroutine dealias_mult(fg_out,f_in,g_in)
!*******************************************************************************
! multiplies f and g with dealiasing

use types, only : rprec
use param, only : ld, ld_big, nx, nx2, ny, ny2
use fft

implicit none

real(rprec), dimension(ld,ny), intent(in) :: f_in, g_in
real(rprec), dimension(ld,ny), intent(out) :: fg_out
real(rprec), dimension(ld,ny) :: dummy1, dummy2
real(rprec), dimension(ld_big,ny2) :: f_big, g_big, fg_big

dummy1 = f_in/(nx*ny)
dummy2 = g_in/(nx*ny)

call dfftw_execute_dft_r2c(forw, dummy1, dummy1)
call dfftw_execute_dft_r2c(forw, dummy2, dummy2)

call padd(f_big,dummy1)
call padd(g_big,dummy2)

call dfftw_execute_dft_c2r(back_big, f_big, f_big)
call dfftw_execute_dft_c2r(back_big, g_big, g_big)

fg_big = f_big*g_big/(nx2*ny2)
call dfftw_execute_dft_r2c(forw_big, fg_big, fg_big)
call unpadd(fg_out,fg_big)
call dfftw_execute_dft_c2r(back, fg_out, fg_out)

end subroutine dealias_mult

!*******************************************************************************
subroutine ke_dealiased_calc(ke_dealiased,u_in,v_in,w_in)
!*******************************************************************************

use types, only : rprec
use param, only : coord, ld, ld_big, nx, nx2, ny, ny2, nz
use sim_param, only : u, v, w
use fft

implicit none

real(rprec), dimension(ld,ny), intent(out) :: ke_dealiased
real(rprec), dimension(ld_big,ny2) :: u_big, v_big, w_big, ke_big
real(rprec), dimension(ld,ny), intent(in) :: u_in, v_in, w_in

! pad to big domain for dealiasing
! (u,v,w assumed to be on uv grid)
call to_big(u_in,u_big)
call to_big(v_in,v_big)
call to_big(w_in,w_big)

! compute kinetic energy after padding
ke_big = u_big**2 + v_big**2 + w_big**2
ke_big = 0.5_rprec*ke_big/(nx2*ny2)

! return to regular domain
call dfftw_execute_dft_r2c(forw_big, ke_big, ke_big)
call unpadd(ke_dealiased,ke_big)
call dfftw_execute_dft_c2r(back, ke_dealiased, ke_dealiased)

end subroutine ke_dealiased_calc


!*******************************************************************************
subroutine to_big(f,f_big)
!*******************************************************************************
! returns padded variable for dealiasing

use types, only : rprec
use param, only : ld, ld_big, nx, nx2, ny, ny2
use fft

implicit none

real(rprec), dimension(ld,ny), intent(in) :: f
real(rprec), dimension(ld_big,ny2), intent(out) :: f_big
real(rprec), dimension(ld,ny) :: dummy

dummy = f/(nx*ny)
call dfftw_execute_dft_r2c(forw, dummy, dummy)
call padd(f_big,dummy)
call dfftw_execute_dft_c2r(back_big, f_big, f_big)

end subroutine to_big

!*******************************************************************************
subroutine get_wallmodelstress(twx_rewm,twy_rewm,twx_fluc,twy_fluc)
!*******************************************************************************
use param, only : nx, ny, coord
use types, only : rprec

implicit none

real(rprec), dimension(nx,ny), intent(inout) :: twx_rewm,twy_rewm,&
    twx_fluc,twy_fluc

if (coord==0) then

    twx_rewm = -twxbar
    twy_rewm = -twybar
    twx_fluc = -twxpp
    twy_fluc = -twypp

else

    twx_rewm = twxbar
    twy_rewm = twybar
    twx_fluc = twxpp
    twy_fluc = twypp

endif

end subroutine get_wallmodelstress

!*******************************************************************************
subroutine mts_monitor
!*******************************************************************************
!
! This subroutine monitors the temporal evolution of neqwm quantities
!
use param, only : jt_total, total_time, path, nx, ny, dt
use sim_param, only : u, v
!use qeqwm, only : utau, utx, uty, ssx, ssy, Ts, Us, vtx, vty

implicit none

integer :: fid, i, j
character*50 :: fname

i = int(nx/2._rprec)
j = int(ny/2._rprec)

fname = path // 'output/mts_monitor_point.dat'
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,&
    total_time,&
    dt,&
    Deltay,&
    twxbar(i,j),&
    twybar(i,j),&
    twxpp(i,j),&
    twypp(i,j),&
    uinst(i,j),&
    vinst(i,j),&
    Ud(i,j),&
    dplesdx(i,j),&
    dplesdy(i,j),&
    dpdxbar1(i,j),&
    dpdybar1(i,j),&
    dpdxbar2(i,j),&
    dpdybar2(i,j),&
    dpdxpp(i,j),&
    dpdypp(i,j),&
    theta_d(i,j),&
    ! qeqwm variables
    utau(i,j),&
    utx(i,j),&
    uty(i,j),&
    ssx(i,j),&
    ssy(i,j),&
    Ts(i,j),&
    Us(i,j),&
    vtx(i,j),&
    vty(i,j),&
    uinf(i,j),&
    vinf(i,j),&
    Tdelta(i,j),&
    dpdxbar3(i,j),&
    dpdybar3(i,j),&
    dpdxppdelta(i,j),&
    dpdyppdelta(i,j)
close(fid)

end subroutine mts_monitor

!*******************************************************************************
subroutine mts_monitor_plane
!*******************************************************************************
!
! This subroutine monitors the temporal evolution of wm quantities (plane
! averaged)
!
use param, only : jt_total, total_time, path, nx, ny, dt
!use qeqwm, only : utau, utx, uty, ssx, ssy, Ts, Us, vtx, vty

implicit none

integer :: fid
character*50 :: fname

fname = path // 'output/mts_monitor_plane.dat'
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,&
    total_time,&
    dt,&
    sum(twxbar)/nx/ny,&
    sum(twybar)/nx/ny,&
    sum(twxpp)/nx/ny,&
    sum(twypp)/nx/ny,&
    sum(uinst)/nx/ny,&
    sum(vinst)/nx/ny,&
    sum(Ud)/nx/ny,&
    sum(dplesdx)/nx/ny,&
    sum(dplesdy)/nx/ny,&
    sum(dpdxbar1)/nx/ny,&
    sum(dpdybar1)/nx/ny,&
    sum(dpdxbar2)/nx/ny,&
    sum(dpdybar2)/nx/ny,&
    sum(dpdxpp)/nx/ny,&
    sum(dpdypp)/nx/ny,&
    sum(theta_d)/nx/ny,&
    ! qeqwm variables
    sum(utau)/nx/ny,&
    sum(utx)/nx/ny,&
    sum(uty)/nx/ny,&
    sum(ssx)/nx/ny,&
    sum(ssy)/nx/ny,&
    sum(Ts)/nx/ny,&
    sum(Us)/nx/ny,&
    sum(vtx)/nx/ny,&
    sum(vty)/nx/ny,&
    sum(uinf)/nx/ny,&
    sum(vinf)/nx/ny,&
    sum(Tdelta)/nx/ny,&
    sum(dpdxbar3)/nx/ny,&
    sum(dpdybar3)/nx/ny,&
    sum(dpdxppdelta)/nx/ny,&
    sum(dpdyppdelta)/nx/ny
close(fid)

end subroutine mts_monitor_plane

!*******************************************************************************
subroutine mts_write_checkpoint
!*******************************************************************************

use param, only : nx, ny, coord, nproc, write_endian
use types, only : rprec
use qeqwm, only : qeqwm_write_checkpoint
use neqwm, only : neqwm_write_I_checkpoint

implicit none

integer :: fid1, fid2

if (coord==0) then
    open(newunit=fid1, file='mts_checkPoint_bot.bin', convert=write_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    write(fid1,rec=1) twxbar(1:nx,1:ny)
    write(fid1,rec=2) twybar(1:nx,1:ny)
    write(fid1,rec=3) twxpp(1:nx,1:ny)
    write(fid1,rec=4) twypp(1:nx,1:ny)
    write(fid1,rec=5) uinst(1:nx,1:ny)
    write(fid1,rec=6) vinst(1:nx,1:ny)
    write(fid1,rec=7) Ud(1:nx,1:ny)
    write(fid1,rec=8) dplesdx(1:nx,1:ny)
    write(fid1,rec=9) dplesdy(1:nx,1:ny)
    write(fid1,rec=10) dpdxbar1(1:nx,1:ny)
    write(fid1,rec=11) dpdybar1(1:nx,1:ny)
    write(fid1,rec=12) dpdxbar2(1:nx,1:ny)
    write(fid1,rec=13) dpdybar2(1:nx,1:ny)
    write(fid1,rec=14) dpdxpp(1:nx,1:ny)
    write(fid1,rec=15) dpdypp(1:nx,1:ny)
    write(fid1,rec=16) theta_d(1:nx,1:ny)
    write(fid1,rec=17) ubar(1:nx,1:ny)
    write(fid1,rec=18) vbar(1:nx,1:ny)
    write(fid1,rec=19) uinf(1:nx,1:ny)
    write(fid1,rec=20) vinf(1:nx,1:ny)
    write(fid1,rec=21) Tdelta(1:nx,1:ny)
    write(fid1,rec=22) dpdxbar3(1:nx,1:ny)
    write(fid1,rec=23) dpdybar3(1:nx,1:ny)
    write(fid1,rec=24) dpdxppdelta(1:nx,1:ny)
    write(fid1,rec=25) dpdyppdelta(1:nx,1:ny)
    close(fid1)
else
    open(newunit=fid2, file='mts_checkPoint_top.bin', convert=write_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    write(fid2,rec=1) twxbar(1:nx,1:ny)
    write(fid2,rec=2) twybar(1:nx,1:ny)
    write(fid2,rec=3) twxpp(1:nx,1:ny)
    write(fid2,rec=4) twypp(1:nx,1:ny)
    write(fid2,rec=5) uinst(1:nx,1:ny)
    write(fid2,rec=6) vinst(1:nx,1:ny)
    write(fid2,rec=7) Ud(1:nx,1:ny)
    write(fid2,rec=8) dplesdx(1:nx,1:ny)
    write(fid2,rec=9) dplesdy(1:nx,1:ny)
    write(fid2,rec=10) dpdxbar1(1:nx,1:ny)
    write(fid2,rec=11) dpdybar1(1:nx,1:ny)
    write(fid2,rec=12) dpdxbar2(1:nx,1:ny)
    write(fid2,rec=13) dpdybar2(1:nx,1:ny)
    write(fid2,rec=14) dpdxpp(1:nx,1:ny)
    write(fid2,rec=15) dpdypp(1:nx,1:ny)
    write(fid2,rec=16) theta_d(1:nx,1:ny)
    write(fid2,rec=17) ubar(1:nx,1:ny)
    write(fid2,rec=18) vbar(1:nx,1:ny)
    write(fid2,rec=19) uinf(1:nx,1:ny)
    write(fid2,rec=20) vinf(1:nx,1:ny)
    write(fid2,rec=21) Tdelta(1:nx,1:ny)
    write(fid2,rec=22) dpdxbar3(1:nx,1:ny)
    write(fid2,rec=23) dpdybar3(1:nx,1:ny)
    write(fid2,rec=24) dpdxppdelta(1:nx,1:ny)
    write(fid2,rec=25) dpdyppdelta(1:nx,1:ny)
    close(fid2)
end if

call qeqwm_write_checkpoint
call neqwm_write_I_checkpoint

end subroutine mts_write_checkpoint

!*******************************************************************************
subroutine mts_read_checkpoint
!*******************************************************************************

use param, only : nx, ny, coord, nproc, read_endian, lbc_mom, ubc_mom
use types, only : rprec
use qeqwm, only : qeqwm_read_checkpoint
use neqwm, only : neqwm_read_I_checkpoint

implicit none

integer :: fid1, fid2
logical :: file_flag

if (coord==0) then
    open(newunit=fid1, file='mts_checkPoint_bot.bin', convert=read_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    read(fid1,rec=1) twxbar(1:nx,1:ny)
    read(fid1,rec=2) twybar(1:nx,1:ny)
    read(fid1,rec=3) twxpp(1:nx,1:ny)
    read(fid1,rec=4) twypp(1:nx,1:ny)
    read(fid1,rec=5) uinst(1:nx,1:ny)
    read(fid1,rec=6) vinst(1:nx,1:ny)
    read(fid1,rec=7) Ud(1:nx,1:ny)
    read(fid1,rec=8) dplesdx(1:nx,1:ny)
    read(fid1,rec=9) dplesdy(1:nx,1:ny)
    read(fid1,rec=10) dpdxbar1(1:nx,1:ny)
    read(fid1,rec=11) dpdybar1(1:nx,1:ny)
    read(fid1,rec=12) dpdxbar2(1:nx,1:ny)
    read(fid1,rec=13) dpdybar2(1:nx,1:ny)
    read(fid1,rec=14) dpdxpp(1:nx,1:ny)
    read(fid1,rec=15) dpdypp(1:nx,1:ny)
    read(fid1,rec=16) theta_d(1:nx,1:ny)
    read(fid1,rec=17) ubar(1:nx,1:ny)
    read(fid1,rec=18) vbar(1:nx,1:ny)
    read(fid1,rec=19) uinf(1:nx,1:ny)
    read(fid1,rec=20) vinf(1:nx,1:ny)
    read(fid1,rec=21) Tdelta(1:nx,1:ny)
    read(fid1,rec=22) dpdxbar3(1:nx,1:ny)
    read(fid1,rec=23) dpdybar3(1:nx,1:ny)
    read(fid1,rec=24) dpdxppdelta(1:nx,1:ny)
    read(fid1,rec=25) dpdyppdelta(1:nx,1:ny)
    close(fid1)
else if (coord==nproc-1) then
    open(newunit=fid2, file='mts_checkPoint_top.bin', convert=read_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec)
    read(fid2,rec=1) twxbar(1:nx,1:ny)
    read(fid2,rec=2) twybar(1:nx,1:ny)
    read(fid2,rec=3) twxpp(1:nx,1:ny)
    read(fid2,rec=4) twypp(1:nx,1:ny)
    read(fid2,rec=5) uinst(1:nx,1:ny)
    read(fid2,rec=6) vinst(1:nx,1:ny)
    read(fid2,rec=7) Ud(1:nx,1:ny)
    read(fid2,rec=8) dplesdx(1:nx,1:ny)
    read(fid2,rec=9) dplesdy(1:nx,1:ny)
    read(fid2,rec=10) dpdxbar1(1:nx,1:ny)
    read(fid2,rec=11) dpdybar1(1:nx,1:ny)
    read(fid2,rec=12) dpdxbar2(1:nx,1:ny)
    read(fid2,rec=13) dpdybar2(1:nx,1:ny)
    read(fid2,rec=14) dpdxpp(1:nx,1:ny)
    read(fid2,rec=15) dpdypp(1:nx,1:ny)
    read(fid2,rec=16) theta_d(1:nx,1:ny)
    read(fid2,rec=17) ubar(1:nx,1:ny)
    read(fid2,rec=18) vbar(1:nx,1:ny)
    read(fid2,rec=19) uinf(1:nx,1:ny)
    read(fid2,rec=20) vinf(1:nx,1:ny)
    read(fid2,rec=21) Tdelta(1:nx,1:ny)
    read(fid2,rec=22) dpdxbar3(1:nx,1:ny)
    read(fid2,rec=23) dpdybar3(1:nx,1:ny)
    read(fid2,rec=24) dpdxppdelta(1:nx,1:ny)
    read(fid2,rec=25) dpdyppdelta(1:nx,1:ny)
    close(fid2)
end if

if (coord==0) then
    inquire (file='qeqwm_checkPoint_bot.bin', exist=file_flag)
    if (file_flag .and. lbc_mom==4) then
        call qeqwm_read_checkPoint()
        write(*,*) 'reading qeqwm_checkpoint_bot'
    end if
    inquire (file='neqwm_I_checkPoint_bot.bin', exist=file_flag)
    if (file_flag .and. lbc_mom==4) then
        call neqwm_read_I_checkPoint()
        write(*,*) 'reading neqwm_I_checkpoint_bot'
    end if
else
    inquire (file='qeqwm_checkPoint_top.bin', exist=file_flag)
    if (file_flag .and. ubc_mom==4) then
        call qeqwm_read_checkPoint()
        write(*,*) 'reading qeqwm_checkpoint_top'
    end if
    inquire (file='neqwm_I_checkPoint_top.bin', exist=file_flag)
    if (file_flag .and. ubc_mom==4) then
        call neqwm_read_I_checkPoint()
        write(*,*) 'reading neqwm_I_checkpoint_top'
    end if
endif

end subroutine mts_read_checkpoint

!*******************************************************************************
subroutine debug_wm
!*******************************************************************************
use param, only : nx, ny, jt, coord, write_endian
!use qeqwm, only : utx, uty

implicit none

integer fid, nt
character*50 :: fname

nt = 1000

if (jt .le. nt) then
    utemp2(1:nx,1:ny,jt) = uinst(1:nx,1:ny)
    vtemp2(1:nx,1:ny,jt) = vinst(1:nx,1:ny)
    dpdxtemp(1:nx,1:ny,jt) = dplesdx(1:nx,1:ny)
    dpdytemp(1:nx,1:ny,jt) = dplesdy(1:nx,1:ny)
    utxtemp(1:nx,1:ny,jt) = utx(1:nx,1:ny)
    utytemp(1:nx,1:ny,jt) = uty(1:nx,1:ny)
endif

if (jt == nt) then
    fname = 'output/debug_wm.bin'
    open(newunit=fid, file=fname, convert=write_endian,&
        form='unformatted', access='direct', recl=nx*ny*rprec*nt)
    write(fid,rec=1) utemp2(1:nx,1:ny,1:nt)
    write(fid,rec=2) vtemp2(1:nx,1:ny,1:nt)
    write(fid,rec=3) dpdxtemp(1:nx,1:ny,1:nt)
    write(fid,rec=4) dpdytemp(1:nx,1:ny,1:nt)
    write(fid,rec=5) utxtemp(1:nx,1:ny,1:nt)
    write(fid,rec=6) utytemp(1:nx,1:ny,1:nt)
endif

end subroutine debug_wm

!*******************************************************************************
subroutine write_velocity_plane
!*******************************************************************************
use param, only : nx, ny, path, jt_total, write_endian
use types, only : rprec
use string_util, only : string_splice

character*50 :: fname

call string_splice(fname, path // 'output/velocity_', jt_total, '.bin')
open(unit=13,file=fname,form='unformatted',convert=write_endian,   &
                access='direct',recl=nx*ny*rprec)
write(13,rec=1) ubar(1:nx,1:ny)
write(13,rec=2) vbar(1:nx,1:ny)
close(13)

end subroutine write_velocity_plane

!*******************************************************************************
subroutine write_velocity_pt1_pt2
!*******************************************************************************
use param, only : nx, ny, path, jt_total, total_time
use sim_param, only : u, v
use types, only : rprec

integer :: fid, i, j
character*50 :: fname

i = int(nx/2._rprec)
j = int(ny/2._rprec)

fname = path // 'output/velocity_pt1_pt2.dat'
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,&
    total_time,&
    u(i,j,1),&
    v(i,j,1),&
    u(i,j,2),&
    v(i,j,2)
close(fid)

end subroutine write_velocity_pt1_pt2

!*******************************************************************************
subroutine mts_ws_plane
!*******************************************************************************
use param, only : nx, ny, path, jt_total
use sim_param, only : txz, tyz
use types, only : rprec

integer :: fid

character*50 :: fname

fname = path // 'output/mts_ws_plane.dat'
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,&
    sum(txz(1:nx,1:ny,1))/nx/ny,&
    sum(twxbar(1:nx,1:ny))/nx/ny,&
    sum(twxpp(1:nx,1:ny))/nx/ny,&
    sum(twxp(1:nx,1:ny))/nx/ny,&
    sum(tyz(1:nx,1:ny,1))/nx/ny,&
    sum(twybar(1:nx,1:ny))/nx/ny,&
    sum(twypp(1:nx,1:ny))/nx/ny,&
    sum(twyp(1:nx,1:ny))/nx/ny
close(fid)

end subroutine mts_ws_plane

!*******************************************************************************
subroutine mts_test_point
!*******************************************************************************
use param, only : nx, ny, path, jt_total, total_time, dt
use sim_param, only : u, v
use types, only : rprec

integer :: fid, i, j
character*50 :: fname

i = int(nx/2._rprec)
j = int(ny/2._rprec)

fname = path // 'output/mts_test_point.dat'
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,&
!    uinst(i,j),&
!    vinst(i,j),&
!    ubar(i,j),&
!    vbar(i,j),&
!    dplesdx(i,j),&
!    dplesdy(i,j),&
!    dpdx_fit(i,j),&
!    dpdy_fit(i,j),&
!    dpdxpp(i,j),&
!    dpdypp(i,j)
    dplesdx(i,j),&
    dpdxbar1(i,j),&
    dpdxpp(i,j),&
    Tnu(i,j),&
    Ts(i,j),&
    utau(i,j),&
    utau_filt(i,j),&
    dt
    
close(fid)

end subroutine mts_test_point

end module mts_wm

