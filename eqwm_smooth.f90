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
    nproc, wmpt
use sim_param, only : u, v, txz, tyz, dudz, dvdz
use test_filtermodule

private
public eqwm_calc

contains

!*******************************************************************************
subroutine eqwm_calc
!*******************************************************************************

implicit none
integer :: i, j
real(rprec), dimension(nx,ny) :: Us, red, utaus, tauws, retaud_eq, theta_d
real(rprec), dimension(nx,ny) :: kappa4, beta1, beta2, lp, dudz_tot
real(rprec), dimension(ld,ny) :: Ules, Vles
real(rprec) :: Deltay, kappa3 

Deltay = dz*(wmpt - 0.5_rprec)

if (coord==0) then
    Ules = u(:,:,wmpt)
    Vles = v(:,:,wmpt)
else
    Ules = u(:,:,nz-wmpt)
    Vles = v(:,:,nz-wmpt)
end if

call test_filter(Ules)
call test_filter(Vles)
!Ules = sum(Ules(1:nx,1:ny))/nx/ny
!Vles = sum(Vles(1:nx,1:ny))/nx/ny
Us = sqrt(Ules(1:nx,1:ny)**2+Vles(1:nx,1:ny)**2)
theta_d = atan2(Vles(1:nx,1:ny),Ules(1:nx,1:ny))

! Re_\Deltay
red = Us*Deltay/nu_molec
kappa3 = 0.005_rprec
beta1 = (1._rprec + 0.155_rprec*red**-0.03_rprec)**-1
beta2 = 1.7 - (1._rprec + 36._rprec*red**-0.75_rprec)**-1
kappa4 = kappa3**(beta1 - 0.5_rprec)
retaud_eq = kappa4*red**beta1*                            &
    (1._rprec + (kappa3*red)**-beta2)**((beta1-0.5_rprec)/beta2)

utaus = nu_molec*retaud_eq/Deltay

! During equilibrium
tauws = utaus**2

! mixing length at first grid point in inner units
lp = vonk*dz/2._rprec*utaus/nu_molec                          &
    *(1 - exp(-dz/2._rprec*utaus/nu_molec/25._rprec))
! total derivative at first grid point
dudz_tot = utaus**2/nu_molec/2._rprec/lp**2*(-1 + sqrt(1 + 4._rprec*lp**2))

if (coord==0) then
    txz(1:nx,1:ny,1) = -tauws*cos(theta_d)
    tyz(1:nx,1:ny,1) = -tauws*sin(theta_d)
    dudz(1:nx,1:ny,1) = dudz_tot*cos(theta_d)
    dvdz(1:nx,1:ny,1) = dudz_tot*sin(theta_d)
else
    txz(1:nx,1:ny,nz) = tauws*cos(theta_d)
    tyz(1:nx,1:ny,nz) = tauws*sin(theta_d)
    dudz(1:nx,1:ny,nz) = -dudz_tot*cos(theta_d)
    dvdz(1:nx,1:ny,nz) = -dudz_tot*sin(theta_d)
end if

!if ( mod(jt_total,100)==0 ) then
!    call eqwm_monitor()
!end if

end subroutine eqwm_calc

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
    txzbar = 0._rprec
    tyzbar = 0._rprec
    do i = 1, nx
    do j = 1, ny
        txzbar = txzbar + txz(i,j,1)
        tyzbar = tyzbar + tyz(i,j,1)
    enddo
    enddo
    txzbar = txzbar/nx/ny
    tyzbar = tyzbar/nx/ny
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
