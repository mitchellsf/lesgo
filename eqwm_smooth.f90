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
module eqwm_smooth
!*******************************************************************************
use types, only : rprec
use param, only : dz, ld, nx, ny, nz, vonk, nu_molec, jt_total, coord, &
    nproc, wmpt, mean_p_force_x, mean_p_force_y
use sim_param, only : u, v, w, dpdx, dpdy, txz, tyz, dudz, dvdz
use test_filtermodule

private
public eqwm_calc

contains

!*******************************************************************************
subroutine eqwm_calc
!*******************************************************************************

implicit none
integer :: i, j
real(rprec), dimension(nx,ny) :: Us, red, retd, retd2, utau2, theta_d, twx, twy
real(rprec), dimension(nx,ny) :: dpds, psi_p
real(rprec), dimension(nx,ny) :: lp, dudz_tot
real(rprec), dimension(ld,ny) :: Ules, Vles
real(rprec), dimension(nx,ny) :: dplesdx, dplesdy
real(rprec) :: Deltay, kappa3 
real(rprec), dimension(ld,ny) :: utemp, vtemp, wtemp
real(rprec), dimension(ld,ny) :: ke_dealiased, dkedx, dkedy

Deltay = dz*(wmpt - 0.5_rprec)

if (coord==0) then
    Ules = u(:,:,wmpt)
    Vles = v(:,:,wmpt)
    utemp = u(:,:,wmpt)
    vtemp = v(:,:,wmpt)
    if (wmpt==1) then
        wtemp = 0.25_rprec*w(:,:,2)
    else
        wtemp = 0.5_rprec*(w(:,:,wmpt)+w(:,:,wmpt+1))
    endif
    dplesdx = dpdx(1:nx,1:ny,wmpt)
    dplesdy = dpdy(1:nx,1:ny,wmpt)
else
    Ules = u(:,:,nz-wmpt)
    Vles = v(:,:,nz-wmpt)
    utemp = u(:,:,nz-wmpt)
    vtemp = v(:,:,nz-wmpt)
    if (wmpt==1) then
        wtemp = 0.25_rprec*w(:,:,nz-wmpt)
    else
        wtemp = 0.5_rprec*(w(:,:,nz-wmpt)+w(:,:,nz-wmpt-1))
    endif
    dplesdx = dpdx(1:nx,1:ny,nz-wmpt)
    dplesdy = dpdy(1:nx,1:ny,nz-wmpt)
end if
! calculate real pressure gradient 
! (dealiased and gradient computed spectrally
! for accurate estimate of PG)
call ke_dealiased_calc(ke_dealiased,utemp,vtemp,wtemp)
call ddxy_plane(ke_dealiased, dkedx, dkedy)
dplesdx(1:nx,1:ny) = dplesdx(1:nx,1:ny) - dkedx(1:nx,1:ny) &
    - mean_p_force_x
dplesdy(1:nx,1:ny) = dplesdy(1:nx,1:ny) - dkedy(1:nx,1:ny) &
    - mean_p_force_y

!call test_filter(Ules)
!call test_filter(Vles)
!Ules = sum(Ules(1:nx,1:ny))/nx/ny
!Vles = sum(Vles(1:nx,1:ny))/nx/ny
Us = sqrt(Ules(1:nx,1:ny)**2+Vles(1:nx,1:ny)**2)
theta_d = atan2(Vles(1:nx,1:ny),Ules(1:nx,1:ny))

! EQWM with PG effects
red = Us*Deltay/nu_molec
dpds = dplesdx*cos(theta_d) + dplesdy*sin(theta_d)
psi_p = dpds*Deltay**3.0/nu_molec**2.0
call retd_pres_fit(red,psi_p,retd)
twx = (retd*nu_molec/Deltay)**2.0*cos(theta_d)
twy = (retd*nu_molec/Deltay)**2.0*sin(theta_d)

! velocity gradient at first grid point computed without PG in fit
! for simplicity, friction velocity computed using LES velocity at wmpt
call retd_fit(red,retd2)
utau2 = nu_molec*retd2/Deltay
lp = vonk*dz/2._rprec*utau2/nu_molec                          &
    *(1 - exp(-dz/2._rprec*utau2/nu_molec/25._rprec))
! total derivative at first grid point
dudz_tot = utau2**2/nu_molec/2._rprec/lp**2*(-1 + sqrt(1 + 4._rprec*lp**2))

if (coord==0) then
    txz(1:nx,1:ny,1) = -twx
    tyz(1:nx,1:ny,1) = -twy
    dudz(1:nx,1:ny,1) = dudz_tot*cos(theta_d)
    dvdz(1:nx,1:ny,1) = dudz_tot*sin(theta_d)
else
    txz(1:nx,1:ny,nz) = twx
    tyz(1:nx,1:ny,nz) = twy
    dudz(1:nx,1:ny,nz) = -dudz_tot*cos(theta_d)
    dvdz(1:nx,1:ny,nz) = -dudz_tot*sin(theta_d)
end if

!if ( mod(jt_total,1)==0 ) then
!    call eqwm_monitor()
!end if

!write(*,*) 'eqwm2: ',jt_total, sum(twx)/nx/ny, sum(twy)/nx/ny
!write(*,*) 'eqwm_smooth: ',jt_total, dplesdx(1,1), dplesdy(1,1), psi_p(1,1)

end subroutine eqwm_calc

!*******************************************************************************
subroutine retd_pres_fit(redelta_in,psi_p,retd_out)
!*******************************************************************************
use param, only : nx, ny, nu_molec

implicit none
real(rprec), dimension(nx,ny), intent(in) :: redelta_in, psi_p
real(rprec), dimension(nx,ny), intent(out) :: retd_out
real(rprec), dimension(nx,ny) :: retd_eq
real(rprec) :: retd_min, pp, redelta_min
integer :: jx, jy

call retd_fit(redelta_in,retd_eq)
do jx = 1,nx
do jy = 1,ny
if (psi_p(jx,jy) < 0) then
    retd_min = 1.5*(-psi_p(jx,jy))**0.39*&
        (1.0+(1000.0/-psi_p(jx,jy))**2.0)**-0.055
    pp = 2.5 - 0.6*(1.0 + tanh(2.0*(log10(-psi_p(jx,jy))-6.0)))
    retd_out(jx,jy) = (retd_min**pp + retd_eq(jx,jy)**pp)**(1.0/pp)
else
    redelta_min = 2.5*psi_p(jx,jy)**0.54*(1.0+(30.0/psi_p(jx,jy))**0.5)**-0.88
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
subroutine ddxy_plane (f, dfdx, dfdy)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! x and y using spectral decomposition. Done for single z-plane.
!
use types, only : rprec
use param, only : ld, nx, ny
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

real(rprec), dimension(:,:), intent(in) :: f
real(rprec), dimension(:,:), intent(inout) :: dfdx, dfdy
real(rprec) :: const

const = 1._rprec / ( nx * ny )

! Use dfdy to hold f; since we are doing in place FFTs this is required
dfdx = const*f
call dfftw_execute_dft_r2c(forw, dfdx, dfdx)

! Zero padded region and Nyquist frequency
dfdx(ld-1:ld,:) = 0._rprec
dfdx(:,ny/2+1) = 0._rprec

! Derivatives: must to y's first here, because we're using dfdx as storage
! Use complex emulation of dfdy to perform complex multiplication
! Optimized version for real(eye*ky)=0
! only passing imaginary part of eye*ky
dfdy = dfdx .MULI. ky
dfdx = dfdx .MULI. kx

! Perform inverse transform to get pseudospectral derivative
call dfftw_execute_dft_c2r(back, dfdx, dfdx)
call dfftw_execute_dft_c2r(back, dfdy, dfdy)

end subroutine ddxy_plane

!*******************************************************************************
subroutine eqwm_monitor
!*******************************************************************************
!
! This subroutine monitors the temporal evolution of eqwm quantities
!
use param, only : jt_total, total_time, path

implicit none

real(rprec) :: txzbar, tyzbar
integer :: fid, i, j, k
character*50 :: fname

if (coord==0) then
!    txzbar = 0._rprec
!    tyzbar = 0._rprec
!    do i = 1, nx
!    do j = 1, ny
!        txzbar = txzbar + txz(i,j,1)
!        tyzbar = tyzbar + tyz(i,j,1)
!    enddo
!    enddo
!    txzbar = txzbar/nx/ny
!    tyzbar = tyzbar/nx/ny
    txzbar = 0._rprec
    tyzbar = 0._rprec
    do j = 1, ny
        txzbar = txzbar + txz(1,j,1)
        tyzbar = tyzbar + tyz(1,j,1)
    enddo
    txzbar = txzbar/ny
    tyzbar = tyzbar/ny
    fname = path // 'output/eqwm_track_bot.dat'
    open(newunit=fid, file=fname, status='unknown',   &
        position='append')
    if (jt_total==100) then
        write(fid,*) 'time step    total time    txz    tyz    '
    end if
    write(fid,*) jt_total, total_time, txzbar, tyzbar
    close(fid)
else
    txzbar = 0._rprec
    tyzbar = 0._rprec
    do i = 1, nx
    do j = 1, ny
        txzbar = txzbar + txz(i,j,nz)
        tyzbar = tyzbar + tyz(i,j,nz)
    enddo
    enddo
    txzbar = txzbar/nx/ny
    tyzbar = tyzbar/nx/ny
    fname = path // 'output/eqwm_track_top.dat'
    open(newunit=fid, file=fname, status='unknown',   &
        position='append')
    if (jt_total==100) then
        write(fid,*) 'time step    total time    txz    tyz    '
    end if
    write(fid,*) jt_total, total_time, txzbar, tyzbar
    close(fid)
end if


end subroutine eqwm_monitor

end module eqwm_smooth
