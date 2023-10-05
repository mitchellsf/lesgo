!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
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
module functions_bl
!*******************************************************************************
use param, only : nx
use types, only : rprec
implicit none

save
private
public initialize_bl, bl_mean_velocity, wtop_bl, ws_inter

real(rprec), dimension(:), allocatable :: wtop_bl
real(rprec), dimension(:,:), allocatable :: gamma_r
real(rprec) :: total_time_ws_inter

contains

!*******************************************************************************
subroutine initialize_bl()
!*******************************************************************************
use param, only : nx, ny, nz, nz_tot, nu_molec, coord, nproc
use sim_param, only : u, v, w

implicit none
real(rprec), dimension(nx) :: utau_init
!real(rprec), dimension(nz_tot) :: wbar_inlt
real(rprec), dimension(nx,nz_tot) :: ubar_init, wbar_init
real(rprec), dimension(nx,ny,nz) :: ufluc, vfluc, wfluc
integer :: jx, jz, kstart

kstart = coord*(nz-1) + 1
u = 0._rprec
v = 0._rprec
w = 0._rprec

! mean velocity calculations
call bl_mean_velocity(ubar_init,utau_init,wbar_init)
!wbar_inlt = wbar_init(1,1:nz_tot)
!do jz = 1,nz
!    w(1,:,jz) = wbar_inlt(kstart+jz-1)
!enddo

! fluctuations (random noise)
call bl_fluctuations(utau_init,ufluc,vfluc,wfluc)
!ufluc = 0._rprec
!vfluc = 0._rprec
!wfluc = 0._rprec

do jx = 1,nx
do jz = 1,nz
    u(jx,:,jz) = ubar_init(jx,kstart+jz-1) + ufluc(jx,:,jz)
    v(jx,:,jz) = vfluc(jx,:,jz)
    w(jx,:,jz) = wbar_init(jx,kstart+jz-1) + wfluc(jx,:,jz)
enddo
enddo

! call afterwards to make sure top boundary has the correct velocity
if (coord==nproc-1) then
    do jx = 1,nx
        w(jx,:,nz) = wtop_bl(jx)
    enddo
endif

end subroutine initialize_bl

!*******************************************************************************
subroutine bl_mean_velocity(ubar,utau,wbar)
!*******************************************************************************
use param, only : nu_molec, dx, dz, nx, ny, nz, nz_tot, coord, nproc,&
    L_x, fringe_region_len, pi, suction_blowing, phi_top_sb

implicit none
real(rprec), dimension(nx,nz_tot) :: up, ubar, wbar
real(rprec), dimension(nz_tot) :: z_uv
real(rprec), dimension(nx) :: utau, Retheta, Redeltas, ddx_deltas, &
    xx, w1, w2
real(rprec) :: ddx_ubar,&
    qout,qtot,a1,a2,k1,k2,x1,x2,xp,wtop1,wtop2
integer :: jx, jy, jz, ix1, ix2, ixp, ixf, jxm, jxp

ixf = modulo(floor((1._rprec-fringe_region_len)*nx + 1._rprec) - 1, nx) + 1

do jz = 1, nz_tot
    z_uv(jz) = (jz-0.5_rprec)*dz
enddo

Retheta(1) = 1._rprec/nu_molec
up(1,:) = iterate_velocity(z_uv,1._rprec/nu_molec)
Redeltas(1) = compute_Redeltas_corrected(up(1,:),z_uv/nu_molec/up(1,nz_tot))
do jx = 2,ixf-1
    ! Based on Cf = 2*(d theta/dx) -> integrated momentum equation for ZPGTBL
    Retheta(jx) = Retheta(jx-1) + dx/nu_molec/(up(jx-1,nz_tot))**2.0
    ! Find velocity profile (based on Monkewitz fit) that matches this Retheta
    up(jx,:) = iterate_velocity(z_uv,Retheta(jx))
    Redeltas(jx) = compute_Redeltas_corrected(up(jx,:),z_uv/nu_molec/up(jx,nz_tot))
    ddx_deltas(jx) = (Redeltas(jx)-Redeltas(jx-1))/dx*nu_molec
enddo
! backwards integration in fringe region
Retheta(nx) = Retheta(1) - dx/nu_molec/(up(1,nz_tot))**2.0
up(nx,:) = iterate_velocity(z_uv,Retheta(nx))
Redeltas(nx) = compute_Redeltas_corrected(up(nx,:),z_uv/nu_molec/up(nx,nz_tot))
ddx_deltas(nx) = (Redeltas(1)-Redeltas(nx))/dx*nu_molec
do jx = nx-1,ixf,-1
    ! Based on Cf = 2*(d theta/dx) -> integrated momentum equation for ZPGTBL
    Retheta(jx) = Retheta(jx+1) - dx/nu_molec/(up(jx+1,nz_tot))**2.0
    ! Find velocity profile (based on Monkewitz fit) that matches this Retheta
    up(jx,:) = iterate_velocity(z_uv,Retheta(jx))
    Redeltas(jx) = compute_Redeltas_corrected(up(jx,:),z_uv/nu_molec/up(jx,nz_tot))
    ddx_deltas(jx) = (Redeltas(jx+1)-Redeltas(jx))/dx*nu_molec
enddo
ddx_deltas(1) = (Redeltas(1)-Redeltas(nx))/dx*nu_molec
utau = 1._rprec/up(:,nz_tot)
do jx = 1,nx
    ubar(jx,:) = up(jx,:)*utau(jx)
enddo

! vertical velocity
wbar = 0._rprec
do jx = 1,ixf-1
    jxm = modulo(jx-2,nx)+1
do jz = 2,nz_tot
    ddx_ubar = (ubar(jx,jz)-ubar(jxm,jz))/dx
    wbar(jx,jz) = wbar(jx,jz-1) - ddx_ubar*dz
enddo
enddo
! fringe region
do jx = ixf,nx
    jxp = modulo(jx,nx)+1
do jz = 2,nz_tot
    ddx_ubar = (ubar(jxp,jz)-ubar(jx,jz))/dx
    wbar(jx,jz) = wbar(jx,jz-1) - ddx_ubar*dz
enddo
enddo

! modify velocity on top surface to account for BL growth
!if (coord==nproc-1) then

ix1 = modulo(floor((1._rprec-2._rprec*fringe_region_len)*nx + 1._rprec) - 1, nx) + 1
ix2 = modulo(floor((1._rprec-fringe_region_len)*nx + 1._rprec) - 1, nx) + 1
!ix1 = modulo(floor((1._rprec-fringe_region_len)*nx + 1._rprec) - 1, nx) + 1
!ix2 = modulo(floor((1._rprec-0.75_rprec*fringe_region_len)*nx + 1._rprec) - 1, nx) + 1
if (modulo(ix2-ix1,2)/=0) then
    ix2 = ix2-1
endif
ixp = ix1+(ix2-ix1)/2
x1 = (ix1-1)*dx
x2 = (ix2-1)*dx
xp = (ixp-1)*dx
if (suction_blowing) then
    wtop_bl = ddx_deltas
    do jx = 1,nx
        wtop_bl(jx) = phi_top_sb
    enddo
    call suction_blowing_bc()
else
    wtop_bl = ddx_deltas
    wtop_bl(ix2:nx) = wtop_bl(1) ! constant in fringe region
endif
qout = wtop_bl(1)*(L_x-x2)
!qout = 0._rprec
!qtot = 0.5_rprec*(wtop_bl(1)+wtop_bl(nx))*dx
do jx = 2,ix1
    qout = qout + 0.5_rprec*(wtop_bl(jx)+wtop_bl(jx-1))*dx
enddo
!do jx = ix2+1,nx
!    qout = qout + 0.5_rprec*(wtop_bl(jx)+wtop_bl(jx-1))*dx
!enddo
wtop1 = wtop_bl(ix1)
wtop2 = wtop_bl(ix2)
a1 = qout/(x2-x1) + 0.75_rprec*wtop1 + 0.25_rprec*wtop2
a2 = qout/(x2-x1) + 0.25_rprec*wtop1 + 0.75_rprec*wtop2
k1 = 2._rprec*pi/(x2-x1)
k2 = 2._rprec*pi/(x2-x1)
do jx = 1,nx
    xx(jx) = (jx-1)*dx
enddo
w1 = a1*(cos(k1*(xx-x1))-1._rprec)+wtop1
w2 = a2*(cos(k2*(xx-x2))-1._rprec)+wtop2
wtop_bl(ix1:ixp) = w1(ix1:ixp)
wtop_bl(ixp:ix2) = w2(ixp:ix2)

!wtop_bl = 0._rprec

qtot = 0.5_rprec*(wtop_bl(1)+wtop_bl(nx))*dx
do jx = 2,nx
    qtot = qtot + 0.5_rprec*(wtop_bl(jx)+wtop_bl(jx-1))*dx
enddo
if (coord==0) then
    write(*,*) 'Qtot: ',qtot
endif




!endif

end subroutine bl_mean_velocity

!*******************************************************************************
subroutine suction_blowing_bc()
!*******************************************************************************
use param, only : nx, dx, fringe_region_len, wmax_sb, sigma_sb, xc_sb, phi_top_sb

implicit none
real(rprec) :: xx
integer :: jx, ix1

ix1 = modulo(floor((1._rprec-fringe_region_len)*nx + 1._rprec) - 1, nx) + 1
do jx = 2,ix1
    xx = (jx-1._rprec)*dx
    wtop_bl(jx) = -sqrt(2._rprec)*wmax_sb*(xx-xc_sb)/sigma_sb*&
        exp(0.5_rprec-((xx-xc_sb)/sigma_sb)**2._rprec) + phi_top_sb
enddo

end subroutine suction_blowing_bc

!*******************************************************************************
subroutine bl_fluctuations(utau,ufluc,vfluc,wfluc)
!*******************************************************************************
use param, only : nx, ny, nz, dz, coord, nu_molec

implicit none
real(rprec), dimension(nx), intent(in) :: utau
real(rprec), dimension(nx,ny,nz), intent(out) :: ufluc, vfluc, wfluc
real(rprec), dimension(nx) :: delta
real(rprec) :: z, rms, sigma_rv
integer :: jx, jz

rms = 3._rprec
sigma_rv = 0.289_rprec

! Fill u, v, and w with uniformly distributed random numbers between 0 and 1
call init_random_seed
call random_number(ufluc)
call random_number(vfluc)
call random_number(wfluc)

! Center random number about 0 and rescale
ufluc = rms / sigma_rv * (ufluc - 0.5_rprec)*utau(1)
vfluc = rms / sigma_rv * (vfluc - 0.5_rprec)*utau(1)
wfluc = rms / sigma_rv * (wfluc - 0.5_rprec)*utau(1)

! fit used to estimate delta
delta = (0.3851_rprec/nu_molec + 72.1823_rprec)*nu_molec/utau(:)

! Rescale noise depending on distance from wall and mean log profile
do jz = 1, nz
do jx = 1, nx
    z = (coord*(nz-1) + jz - 0.5_rprec) * dz
    if (z <= delta(jx)) then
        ufluc(jx,:,jz) = ufluc(jx,:,jz) * (1._rprec-z/delta(jx))
        vfluc(jx,:,jz) = wfluc(jx,:,jz) * (1._rprec-z/delta(jx))
        wfluc(jx,:,jz) = wfluc(jx,:,jz) * (1._rprec-z/delta(jx))
    else
        ufluc(jx,:,jz) = 0._rprec
        vfluc(jx,:,jz) = 0._rprec
        wfluc(jx,:,jz) = 0._rprec
    end if
enddo
enddo


end subroutine bl_fluctuations

!*******************************************************************************
function iterate_velocity(z_uv,Retheta_tar) result(up)
!*******************************************************************************
use param, only : nu_molec,wmpt, nz_tot

implicit none
real(rprec) :: deltap, up_inf, deltap_eps, up_inf_eps, Retheta, Retheta_tar,&
    Retheta_eps, eps, tol
real(rprec), dimension(nz_tot) :: z_uv, up, up_eps
integer :: iter

! initial guess for deltap based on linear fit
deltap = 0.3851_rprec/nu_molec + 72.1823_rprec
up_inf = monkewitz_fit(deltap,deltap)
up = generate_velocity(z_uv/nu_molec/up_inf,deltap)
Retheta = compute_Retheta_corrected(up,z_uv/nu_molec/up_inf)

! iterate until target Retheta is achieved
eps = 0.000001_rprec
tol = 0.000001_rprec
iter = 0
do while (abs(Retheta-Retheta_tar)/Retheta_tar>tol .and. iter .le. 100)
    ! incremental change in deltap
    deltap_eps = deltap + eps
    up_inf_eps = monkewitz_fit(deltap_eps,deltap_eps)
    up_eps = generate_velocity(z_uv/nu_molec/up_inf_eps,deltap_eps)
    Retheta_eps = compute_Retheta_corrected(up_eps,z_uv/nu_molec/up_inf_eps)

    ! Update with Newton-Raphson
    deltap = deltap - (Retheta-Retheta_tar)/(Retheta_eps-Retheta)*eps
    up_inf = monkewitz_fit(deltap,deltap)
    up = generate_velocity(z_uv/nu_molec/up_inf,deltap)
    Retheta = compute_Retheta_corrected(up,z_uv/nu_molec/up_inf)
    iter = iter+1
enddo

end function iterate_velocity

!*******************************************************************************
function generate_velocity(zp,deltap) result(up)
!*******************************************************************************
use param, only : nz_tot

implicit none
real(rprec), dimension(nz_tot) ::  zp, up
real(rprec) :: up_inf, deltap
integer :: jz

up_inf = monkewitz_fit(deltap,deltap)
do jz = 1,nz_tot
if (zp(jz) < deltap) then
    if (zp(jz) == 0) then
        up(jz) = 0
    else
        up(jz) = monkewitz_fit(zp(jz),deltap)
    endif
else
    up(jz) = up_inf
endif
enddo
end function generate_velocity

!*******************************************************************************
function compute_Retheta_corrected(up,zp) result(Retheta)
!*******************************************************************************
use param, only : nz_tot,wmpt

implicit none
real(rprec), dimension(nz_tot) :: up, zp, dummy
real(rprec) :: a, Retheta
integer :: jz

a = up(wmpt)/up(nz_tot)
Retheta = zp(wmpt)*up(nz_tot)*a*(&
    (1.0-a)*(1.0-deltas_fit(zp(wmpt))) + a*theta_fit(zp(wmpt)) )
dummy = up*(1._rprec-up/up(nz_tot))
do jz = wmpt+1,nz_tot
    Retheta = Retheta + 0.5_rprec*(dummy(jz)+dummy(jz-1))*(zp(jz)-zp(jz-1))
enddo

end function compute_Retheta_corrected

!*******************************************************************************
function compute_Redeltas_corrected(up,zp) result(Redeltas)
!*******************************************************************************
use param, only : nz_tot,wmpt

implicit none
real(rprec), dimension(nz_tot) :: up, zp, dummy
real(rprec) :: a, Redeltas
integer :: jz

a = up(wmpt)/up(nz_tot)
Redeltas = zp(wmpt)*up(nz_tot)*(1.0-a*(1.0-deltas_fit(zp(wmpt))))
dummy = up(nz_tot)*(1._rprec-up/up(nz_tot))
do jz = wmpt+1,nz_tot
    Redeltas = Redeltas + 0.5_rprec*(dummy(jz)+dummy(jz-1))*(zp(jz)-zp(jz-1))
enddo

end function compute_Redeltas_corrected

!*******************************************************************************
function monkewitz_fit(zp,deltap) result(up)
!*******************************************************************************

implicit none
real(rprec) :: zp,deltap,up
real(rprec) :: wake,kappa,a,alpha,beta,a1,a2,a3,up_inner,Z,Wexp

wake = 0.446_rprec
kappa = 0.384_rprec
a = 10.306_rprec
alpha = (a-1._rprec/kappa)/2._rprec
beta = sqrt(2._rprec*a*alpha-alpha**2._rprec)
a1 = 132.8410_rprec
a2 = -166.2041_rprec
a3 = 71.9114_rprec
up_inner = 1._rprec/kappa*log((zp+a)/a) + alpha/(a + 4._rprec*alpha)*(&
    (a-4._rprec*alpha)*log(a*((zp-alpha)**2.0+beta**2.0)/(2._rprec*alpha*(zp+a)**2.0))&
    + 2._rprec*alpha*(5._rprec*a-4._rprec*alpha)/beta*(atan((zp-alpha)/beta)&
    + atan(alpha/beta)) )
Z = zp/deltap
Wexp = (1._rprec-0.5_rprec/wake*log(Z))*&
    (1._rprec-exp(Z**4.0*(&
    a1*(Z-1.25_rprec)+a2*(Z**2.0-1.5_rprec)+a3*(Z**3.0-1.75_rprec) )))/&
    (1._rprec-exp(-0.25_rprec*(a1+2._rprec*a2+3._rprec*a3)))
up = up_inner + 2._rprec*wake/kappa*Wexp

end function monkewitz_fit

!*******************************************************************************
function deltas_fit(Deltap) result(deltas)
!*******************************************************************************

implicit none
real(rprec) :: Deltap
real(rprec) :: deltas
real(rprec) :: Redelta, k, c1, c2, c3, c4, gamma1, deltas_log

Redelta = Deltap*f_fit(Deltap)
k = 0.4_rprec
c1 = 23.664_rprec
c2 = 0.0016_rprec
c3 = 1.516_rprec
c4 = 1.177_rprec
gamma1 = 1._rprec/(1._rprec+c2*Redelta**c4)
deltas_log = c1/Redelta+Deltap/Redelta/k
deltas = gamma1/2._rprec + (1._rprec-gamma1)**c3*deltas_log

end function deltas_fit

!*******************************************************************************
function theta_fit(Deltap) result(theta)
!*******************************************************************************

implicit none
real(rprec) :: Deltap
real(rprec) :: theta
real(rprec) :: Redelta, k, c5, c6, c7, c8, gamma2, theta_log

Redelta = Deltap*f_fit(Deltap)
k = 0.4_rprec
c5 = -103.5_rprec
c6 = 2586._rprec
c7 = 0.00154_rprec
c8 = 2.475_rprec
gamma2 = 1._rprec/(1._rprec+c7*Redelta)
theta_log = 1._rprec/Redelta*(c5+Deltap/k) + &
    Deltap/Redelta**2._rprec*(c6-2._rprec*Deltap/k**2._rprec)
theta = gamma2/6._rprec + (1._rprec-gamma2)**c8*theta_log

end function theta_fit

!*******************************************************************************
function f_fit(Deltap) result(func)
!*******************************************************************************

implicit none
real(rprec) :: Deltap
real(rprec) :: func
real(rprec) :: B, k, k2, beta, k1

B = 4.95_rprec
k = 0.4_rprec
k2 = 9.753_rprec
beta = 1.903_rprec
k1 = 1._rprec/k*log(k2)+B
func = (1._rprec/k*log(k2+Deltap)+B)*&
    (1._rprec+(Deltap/k1)**-beta)**(-1._rprec/beta)

end function f_fit


!*******************************************************************************
subroutine ws_inter()
!*******************************************************************************

use param, only : jt_total, nx, ny, dt, write_endian
use param, only : ws_inter_nstart, ws_inter_nend
use sim_param, only : txz

implicit none
integer :: jx, jy
character(64) :: fname

if (jt_total == ws_inter_nstart) then

    allocate(gamma_r(nx,ny))
    gamma_r = 0._rprec
    total_time_ws_inter = 0._rprec

else
    do jx = 1,nx
    do jy = 1,ny
        if (-txz(jx,jy,1) < 0) then
            gamma_r(jx,jy) = gamma_r(jx,jy) + dt
        endif
    enddo
    enddo
    total_time_ws_inter = total_time_ws_inter + dt

endif

if (jt_total == ws_inter_nend) then

    gamma_r = gamma_r/total_time_ws_inter
    fname = 'ws_inter.bin'
    open(unit=13,file=fname,action='write',access='direct',&
        form='unformatted',convert=write_endian,recl=nx*ny*rprec)
    write(13,rec=1) gamma_r(1:nx,1:ny)
    close(13)

endif

end subroutine ws_inter

end module functions_bl
