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
module neqwm
!*******************************************************************************
use types, only : rprec

private
public neqwm_initialize, neqwm_finalize, &
    neqwm_write_I_checkpoint, neqwm_read_I_checkpoint, neq_laminar_calc

! sum-of-exponentials approximation coefficients
real(rprec), dimension(:), allocatable :: ws, ss
! integrals for history part of non-equilibrium wall model
real(rprec), dimension(:,:,:), allocatable :: Ix, Iy
! number of exponential terms
integer :: nexp=0

contains

!*******************************************************************************
subroutine neqwm_initialize
!*******************************************************************************
use param, only : nx, ny, nz, dz, coord
use sim_param, only : u

implicit none
integer :: i, io

nexp=0
open(1,file="soe.dat")
do
    read(1,*,iostat=io)
    if (io/=0) exit
    nexp = nexp + 1
enddo
close(1)

if (coord==0) write(*,*) nexp, ' exponential terms in SOE approx.'

allocate(ws(nexp))
allocate(ss(nexp))
allocate(Ix(nx,ny,nexp))
allocate(Iy(nx,ny,nexp))

open(1,file="soe.dat")
do i = 1,nexp
    read(1,*) ws(i), ss(i)
end do
close(1)

end subroutine neqwm_initialize
          
!*******************************************************************************
subroutine neqwm_finalize
!*******************************************************************************
          
implicit none

deallocate(ws)
deallocate(ss)
deallocate(Ix)
deallocate(Iy)

end subroutine neqwm_finalize

!*******************************************************************************
subroutine neq_laminar_calc(dpdxpp,dpdxpp_m,dpdxpp_mm,dpdypp,dpdypp_m,&
    dpdypp_mm,twxpp,twypp)
!*******************************************************************************
use param, only : dt, nx, ny, nu_molec, dt_f

implicit none
real(rprec), dimension(nx,ny), intent(in) :: dpdxpp,dpdxpp_m,dpdxpp_mm,dpdypp,&
    dpdypp_m,dpdypp_mm
real(rprec), dimension(nx,ny), intent(out) :: twxpp,twypp
real(rprec), dimension(nx,ny) :: dpdx12, dpdx32, dpdy12, dpdy32, sumx, sumy
real(rprec), dimension(nx,ny) :: tauwxpp_l, tauwxpp_h, tauwypp_l, tauwypp_h
real(rprec) :: PI=4*atan(1._rprec)
integer :: i

dpdx12 = 0.5_rprec*(dpdxpp + dpdxpp_m)
dpdx32 = 0.5_rprec*(dpdxpp_m + dpdxpp_mm)

dpdy12 = 0.5_rprec*(dpdypp + dpdypp_m)
dpdy32 = 0.5_rprec*(dpdypp_m + dpdypp_mm)

tauwxpp_l = -2*sqrt(nu_molec*dt/PI)*dpdx12
tauwypp_l = -2*sqrt(nu_molec*dt/PI)*dpdy12

sumx = 0._rprec
sumy = 0._rprec

do i = 1,nexp

    Ix(:,:,i) = exp(-ss(i)*dt)*(Ix(:,:,i) - dpdx32/ss(i)*(1-exp(-ss(i)*dt_f)))
    sumx = sumx + ws(i)*Ix(:,:,i)

    Iy(:,:,i) = exp(-ss(i)*dt)*(Iy(:,:,i) - dpdy32/ss(i)*(1-exp(-ss(i)*dt_f)))
    sumy = sumy + ws(i)*Iy(:,:,i)

end do

tauwxpp_h = sqrt(nu_molec/PI)*sumx
tauwypp_h = sqrt(nu_molec/PI)*sumy

twxpp = tauwxpp_l + tauwxpp_h
twypp = tauwypp_l + tauwypp_h

end subroutine neq_laminar_calc

!*******************************************************************************
subroutine neqwm_write_I_checkpoint
!*******************************************************************************

use param, only : nx, ny, coord, nproc, write_endian
use types, only : rprec

implicit none

integer :: fid1, fid2

if (coord==0) then
    open(newunit=fid1, file='neqwm_I_checkPoint_bot.bin', convert=write_endian,&
        form='unformatted', access='direct', recl=nx*ny*nexp*rprec)
    write(fid1,rec=1) Ix(1:nx,1:ny,1:nexp)
    write(fid1,rec=2) Iy(1:nx,1:ny,1:nexp)
    close(fid1)
else
    open(newunit=fid2, file='neqwm_I_checkPoint_top.bin', convert=write_endian,&
        form='unformatted', access='direct', recl=nx*ny*nexp*rprec)
    write(fid2,rec=1) Ix(1:nx,1:ny,1:nexp)
    write(fid2,rec=2) Iy(1:nx,1:ny,1:nexp)
    close(fid2)
end if

end subroutine neqwm_write_I_checkpoint

!*******************************************************************************
subroutine neqwm_read_I_checkpoint
!*******************************************************************************

use param, only : nx, ny, coord, nproc, read_endian
use types, only : rprec

implicit none

integer :: fid1, fid2

if (coord==0) then
    open(newunit=fid1, file='neqwm_I_checkPoint_bot.bin', convert=read_endian,&
        form='unformatted', access='direct', recl=nx*ny*nexp*rprec)
    read(fid1,rec=1) Ix(1:nx,1:ny,1:nexp)
    read(fid1,rec=2) Iy(1:nx,1:ny,1:nexp)
    close(fid1)
else
    open(newunit=fid2, file='neqwm_I_checkPoint_top.bin', convert=read_endian,&
        form='unformatted', access='direct', recl=nx*ny*nexp*rprec)
    read(fid2,rec=1) Ix(1:nx,1:ny,1:nexp)
    read(fid2,rec=2) Iy(1:nx,1:ny,1:nexp)
    close(fid2)
end if

end subroutine neqwm_read_I_checkpoint

end module neqwm 
