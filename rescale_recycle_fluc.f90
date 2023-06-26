!!
!!  Copyright (C) 2020  Johns Hopkins University
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
module rescale_recycle_fluc
!*******************************************************************************
use types, only : rprec
use fringe

implicit none

private
public rescale_recycle_fluc_init, rescale_recycle_fluc_calc, apply_fringe

integer :: isample, iter

! velocity decomposition
type velocity_t
    real(rprec), dimension(:,:), allocatable :: tot, fluc
    real(rprec), dimension(:), allocatable :: avg
end type velocity_t

! all variables at inlet and sample location
! velocities, thicknesses, coordinates
type sample_inlet_t
    type(velocity_t) :: u, v, w
    real(rprec) :: delta, theta, utau, uinf, tauw
    character*50 :: loc_name
end type sample_inlet_t

real(rprec), dimension(:), allocatable :: z_uv, z_w
type(sample_inlet_t) :: inlt, samp
type(fringe_t) :: apply_fringe
real(rprec) :: tavg
real(rprec), dimension(:,:,:), allocatable :: u_fringe,v_fringe,w_fringe
real(rprec), dimension(:,:), allocatable :: ubar_init, wbar_init
real(rprec), dimension(:), allocatable :: utau_init
real(rprec), dimension(:), allocatable :: u_fluc_avg, v_fluc_avg, w_fluc_avg

contains
!*******************************************************************************
subroutine rescale_recycle_fluc_init
!*******************************************************************************
use param, only : fringe_region_end, fringe_region_len, sampling_region_end
use param, only : nx, ny, nz, dx
use param, only : nz_tot, dz, coord
use functions_bl, only : bl_mean_velocity

!real(rprec), dimension(nx,nz_tot) :: ubar_init, wbar_init
!real(rprec), dimension(nx) :: utau_init
integer :: jx,jz

apply_fringe = fringe_t(fringe_region_end, fringe_region_len)
call sample_inlet_init(inlt)
call sample_inlet_init(samp)
isample = modulo(floor(sampling_region_end*nx + 1._rprec) - 1, nx) + 1
samp%loc_name = 'sample'
inlt%loc_name = 'inlet'

allocate(z_uv(nz_tot))
allocate(z_w(nz_tot))
do jz = 1, nz_tot
    z_uv(jz) = (jz-0.5_rprec)*dz
    z_w(jz) = (jz-1._rprec)*dz
enddo

! initial mean velocities
allocate(ubar_init(nx,nz_tot))
allocate(wbar_init(nx,nz_tot))
allocate(utau_init(nx))
call bl_mean_velocity(ubar_init,utau_init,wbar_init)
samp%u%avg = ubar_init(isample,:)
samp%utau = utau_init(isample)
samp%tauw = samp%utau**2._rprec
samp%w%avg = wbar_init(isample,:)
samp%v%avg = 0._rprec
inlt%u%avg = ubar_init(1,:)
inlt%utau = utau_init(1)
inlt%tauw = inlt%utau**2._rprec
inlt%w%avg = wbar_init(1,:)
inlt%v%avg = 0._rprec
call compute_bl_thickness(inlt)
call compute_bl_thickness(samp)
call compute_momentum_thickness2(samp)
call compute_momentum_thickness2(inlt)

tavg = 0._rprec

allocate(u_fringe(apply_fringe%nx,ny,nz_tot))
allocate(v_fringe(apply_fringe%nx,ny,nz_tot))
allocate(w_fringe(apply_fringe%nx,ny,nz_tot))
!do jx = 1, apply_fringe%nx
!do jz = 1, nz_tot
!    u_fringe(jx,:,jz) = ubar_init(apply_fringe%iwrap(jx),jz)
!    w_fringe(jx,:,jz) = wbar_init(apply_fringe%iwrap(jx),jz)
!enddo
!enddo
do jz = 1,nz_tot
    u_fringe(:,:,jz) = inlt%u%avg(jz)
    w_fringe(:,:,jz) = inlt%w%avg(jz)
enddo
v_fringe = 0._rprec

allocate(u_fluc_avg(nz_tot)); u_fluc_avg = 0._rprec
allocate(v_fluc_avg(nz_tot)); v_fluc_avg = 0._rprec
allocate(w_fluc_avg(nz_tot)); w_fluc_avg = 0._rprec

! Time averaging
! reads time averaged velocity profiles if they exist, tavg & utau read
call read_sample_velocity_average()

end subroutine rescale_recycle_fluc_init

!*******************************************************************************
subroutine sample_inlet_init(this)
!*******************************************************************************

class(sample_inlet_t), intent(inout) :: this

call velocity_init(this%u)
call velocity_init(this%v)
call velocity_init(this%w)

this%delta = 10._rprec
this%theta = 1._rprec
this%utau = 1._rprec
this%uinf = 1._rprec
this%tauw = 1._rprec

end subroutine sample_inlet_init

!*******************************************************************************
subroutine velocity_init(this)
!*******************************************************************************
use param, only : ny, nz_tot

class(velocity_t), intent(inout) :: this

allocate(this%tot(ny,nz_tot))
allocate(this%avg(nz_tot))
allocate(this%fluc(ny,nz_tot))

end subroutine velocity_init

!*******************************************************************************
subroutine rescale_recycle_fluc_calc
!*******************************************************************************
use param
use sim_param, only : u, v, w, txz, fxa, fya, fza, fx, fy, fz

integer :: i, i_w, kstart, kend, jx, jy, jz
real(rprec), dimension(ny,nz_tot) :: dummy1, dummy2, dummy3, force_tscale

!! update fringe weights
!apply_fringe = fringe_t(fringe_region_end, fringe_region_len)

! averaging time scale
if (total_time < start_tavg_bl) then
    if (total_time < start_tavg2_bl) then
        tavg = tavg1_bl
    else
        tavg = tavg2_bl
    endif
else
    tavg = tavg + dt
endif

! Sampling plane
kstart = coord*(nz-1) + 1
kend = (coord+1)*(nz-1)
dummy1 = 0._rprec
dummy2 = 0._rprec
dummy3 = 0._rprec
dummy1(:,kstart:kend) = u(isample,1:ny,1:nz-1)
dummy2(:,kstart:kend) = v(isample,1:ny,1:nz-1)
dummy3(:,kstart:kend) = w(isample,1:ny,1:nz-1)
if (coord==nproc-1) then
    dummy1(:,nz_tot) = u(isample,1:ny,nz)
    dummy2(:,nz_tot) = v(isample,1:ny,nz)
    dummy3(:,nz_tot) = w(isample,1:ny,nz)
endif
call mpi_allreduce(dummy1,samp%u%tot,ny*nz_tot,MPI_RPREC,MPI_SUM, &
     MPI_COMM_WORLD,ierr)
call mpi_allreduce(dummy2,samp%v%tot,ny*nz_tot,MPI_RPREC,MPI_SUM, &
     MPI_COMM_WORLD,ierr)
call mpi_allreduce(dummy3,samp%w%tot,ny*nz_tot,MPI_RPREC,MPI_SUM, &
     MPI_COMM_WORLD,ierr)
! Compute time and spanwise average
samp%u%avg = dt/tavg*sum(samp%u%tot,1)/ny + (1._rprec-dt/tavg)*samp%u%avg
samp%v%avg = dt/tavg*sum(samp%v%tot,1)/ny + (1._rprec-dt/tavg)*samp%v%avg
samp%w%avg = dt/tavg*sum(samp%w%tot,1)/ny + (1._rprec-dt/tavg)*samp%w%avg
samp%uinf = samp%u%avg(nz_tot-1)
do jz = 1,nz_tot
    samp%u%fluc(:,jz) = samp%u%tot(:,jz) - samp%u%avg(jz)
    samp%v%fluc(:,jz) = samp%v%tot(:,jz) - samp%v%avg(jz)
    samp%w%fluc(:,jz) = samp%w%tot(:,jz) - samp%w%avg(jz)
enddo
samp%utau = utau_init(isample)
samp%tauw = dt/tavg*sum(abs(txz(isample,1:ny,1)))/ny &
    + (1._rprec-dt/tavg)*samp%tauw
!samp%utau = sqrt(samp%tauw)
!call compute_bl_thickness(samp)
call compute_momentum_thickness2(samp)

! Mean flow from initialized velocity profile
inlt%u%avg = ubar_init(1,:)
inlt%utau = utau_init(1)
inlt%tauw = utau_init(1)**2._rprec
inlt%w%avg = wbar_init(1,:)
inlt%v%avg = 0._rprec
call compute_bl_thickness(inlt)
call compute_momentum_thickness2(inlt)
! Only rescales the fluctuations (total velocity updated here)
call rescale_fluctuations()
!check fluctuations average to zero
u_fluc_avg = dt/tavg*sum(inlt%u%fluc,1)/ny + (1._rprec-dt/tavg)*u_fluc_avg
v_fluc_avg = dt/tavg*sum(inlt%v%fluc,1)/ny + (1._rprec-dt/tavg)*v_fluc_avg
w_fluc_avg = dt/tavg*sum(inlt%w%fluc,1)/ny + (1._rprec-dt/tavg)*w_fluc_avg

!! outputs everything at a given time step
!if (coord==0 .and. mod(jt_total,1000)==0) then
!    call check_rescale_velocity()
!endif

! simple mirroring technique to reduce error accumulation
! could also add spanwise shift
do jx = 1,apply_fringe%nx
    u_fringe(jx,1,:) = inlt%u%tot(1,:)
    v_fringe(jx,1,:) = -inlt%v%tot(1,:)
    w_fringe(jx,1,:) = inlt%w%tot(1,:)
do jy = 2,ny
    u_fringe(jx,jy,:) = inlt%u%tot(ny-jy+2,:)
    v_fringe(jx,jy,:) = -inlt%v%tot(ny-jy+2,:)
    w_fringe(jx,jy,:) = inlt%w%tot(ny-jy+2,:)
enddo
enddo


!do jx = 1,apply_fringe%nx
!do jz = 1,nz_tot
!!    u_fringe(jx,:,jz) = ubar_init(apply_fringe%iwrap(jx),jz) + inlt%u%fluc(:,jz)
!!    v_fringe(jx,:,jz) = inlt%v%fluc(:,jz)
!!    w_fringe(jx,:,jz) = wbar_init(apply_fringe%iwrap(jx),jz) + inlt%w%fluc(:,jz)
!    u_fringe(jx,:,jz) = inlt%u%avg(jz)
!    v_fringe(jx,:,jz) = inlt%v%avg(jz)
!    w_fringe(jx,:,jz) = inlt%w%avg(jz)
!enddo
!enddo
!u_fringe(apply_fringe%nx,:,:) = inlt%u%tot
!v_fringe(apply_fringe%nx,:,:) = inlt%v%tot
!w_fringe(apply_fringe%nx,:,:) = inlt%w%tot

! Apply inflow conditions
!do i = 1, apply_fringe%nx
!    i_w = apply_fringe%iwrap(i)
!    u(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * u(i_w,1:ny,1:nz)&
!        + apply_fringe%beta(i) * inlt%u%tot(1:ny,kstart:kend+1)
!    v(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * v(i_w,1:ny,1:nz)&
!        + apply_fringe%beta(i) * inlt%v%tot(1:ny,kstart:kend+1)
!    w(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * w(i_w,1:ny,1:nz)&
!        + apply_fringe%beta(i) * inlt%w%tot(1:ny,kstart:kend+1)
!end do
force_tscale = 5._rprec*dt
do i = 1, apply_fringe%nx
    i_w = apply_fringe%iwrap(i)
    fxa(i_w,1:ny,1:nz) = apply_fringe%beta(i)/force_tscale*&
        (u_fringe(i,1:ny,kstart:kend+1) - u(i_w,1:ny,1:nz))
    fya(i_w,1:ny,1:nz) = apply_fringe%beta(i)/force_tscale*&
        (v_fringe(i,1:ny,kstart:kend+1) - v(i_w,1:ny,1:nz))
    fza(i_w,1:ny,1:nz) = apply_fringe%beta(i)/force_tscale*&
        (w_fringe(i,1:ny,kstart:kend+1) - w(i_w,1:ny,1:nz))
end do
!do i = 1, apply_fringe%nx
!    i_w = apply_fringe%iwrap(i)
!    u(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * u(i_w,1:ny,1:nz)&
!        + apply_fringe%beta(i) * u_fringe(i,1:ny,kstart:kend+1)
!    v(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * v(i_w,1:ny,1:nz)&
!        + apply_fringe%beta(i) * v_fringe(i,1:ny,kstart:kend+1)
!    w(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * w(i_w,1:ny,1:nz)&
!        + apply_fringe%beta(i) * w_fringe(i,1:ny,kstart:kend+1)
!end do

!if (mod(jt_total,10)==0) then
!    call write_fringe_force()
!endif

if (coord==0 .and. mod(jt_total,1)==0) then
    call write_sample_velocity_average()
endif

if (coord==0 .and. mod(jt_total,1)==0) then
    call monitor_rescale_recycle()
endif


end subroutine rescale_recycle_fluc_calc

!*******************************************************************************
subroutine rescale_fluctuations()
!*******************************************************************************
use param

implicit none
real(rprec), dimension(ny,nz_tot) :: up_inner, up_outer, vp_inner, vp_outer, &
    wp_inner, wp_outer
real(rprec) :: gg, a, b, weight_uv, weight_w, fluc_damp
integer :: jz

fluc_damp = 2._rprec
a = 4._rprec
b = 0.2_rprec
gg = inlt%utau/samp%utau
do jz = 1,nz_tot
    up_inner(:,jz) = gg*lin_interp2(samp%u%fluc,z_uv*samp%utau/nu_molec,&
        z_uv(jz)*inlt%utau/nu_molec)
    vp_inner(:,jz) = gg*lin_interp2(samp%v%fluc,z_uv*samp%utau/nu_molec,&
        z_uv(jz)*inlt%utau/nu_molec)
    wp_inner(:,jz) = gg*lin_interp2(samp%w%fluc,z_w*samp%utau/nu_molec,&
        z_w(jz)*inlt%utau/nu_molec)
!    up_outer(:,jz) = gg*lin_interp2(samp%u%fluc,z_uv/samp%theta,z_uv(jz))
!    vp_outer(:,jz) = gg*lin_interp2(samp%v%fluc,z_uv/samp%theta,z_uv(jz))
!    wp_outer(:,jz) = gg*lin_interp2(samp%w%fluc,z_w/samp%theta,z_w(jz))
    up_outer(:,jz) = gg*lin_interp2(samp%u%fluc,z_uv/samp%delta,z_uv(jz)/inlt%delta)
    vp_outer(:,jz) = gg*lin_interp2(samp%v%fluc,z_uv/samp%delta,z_uv(jz)/inlt%delta)
    wp_outer(:,jz) = gg*lin_interp2(samp%w%fluc,z_w/samp%delta,z_w(jz)/inlt%delta)
    weight_uv = 0.5_rprec*(1.0+tanh( a*(z_uv(jz)/inlt%delta-b)/&
        ((1.0-2.0*b)*z_uv(jz)/inlt%delta+b) )/tanh(a))
    weight_w = 0.5_rprec*(1.0+tanh( a*(z_w(jz)/inlt%delta-b)/&
        ((1.0-2.0*b)*z_w(jz)/inlt%delta+b) )/tanh(a))
    if (z_uv(jz)>inlt%delta) then
        weight_uv = 1._rprec
        if (z_uv(jz)<fluc_damp*inlt%delta) then
            up_outer(:,jz) = up_outer(:,jz)*cos((z_uv(jz)/inlt%delta-1.0)*&
                pi/2._rprec/(fluc_damp-1._rprec))
            vp_outer(:,jz) = vp_outer(:,jz)*cos((z_uv(jz)/inlt%delta-1.0)*&
                pi/2._rprec/(fluc_damp-1._rprec))
        else
            up_outer(:,jz) = 0._rprec
            vp_outer(:,jz) = 0._rprec
        endif
    endif
    if (z_w(jz)>inlt%delta) then
        weight_w = 1._rprec
        if (z_w(jz)<fluc_damp*inlt%delta) then
            wp_outer(:,jz) = wp_outer(:,jz)*cos((z_w(jz)/inlt%delta-1.0)*&
                pi/2._rprec/(fluc_damp-1._rprec))
        else
            wp_outer(:,jz) = 0._rprec
        endif
    endif
    inlt%u%tot(:,jz) = inlt%u%avg(jz) + up_inner(:,jz)*(1._rprec-weight_uv) +&
        (up_outer(:,jz))*weight_uv
    inlt%v%tot(:,jz) = vp_inner(:,jz)*(1._rprec-weight_uv) +&
        (vp_outer(:,jz))*weight_uv
    inlt%w%tot(:,jz) = inlt%w%avg(jz) + wp_inner(:,jz)*(1._rprec-weight_w) +&
        (wp_outer(:,jz))*weight_w
    inlt%u%fluc(:,jz) = inlt%u%tot(:,jz) - inlt%u%avg(jz)
    inlt%v%fluc(:,jz) = inlt%v%tot(:,jz) - inlt%v%avg(jz)
    inlt%w%fluc(:,jz) = inlt%w%tot(:,jz) - inlt%w%avg(jz)
enddo

end subroutine rescale_fluctuations

!*******************************************************************************
function lin_interp1(phi,grid,pt) result(phi_i)
!*******************************************************************************
use param, only : nz_tot
use functions, only : binary_search

implicit none
real(rprec), dimension(:), intent(in) :: phi, grid
real(rprec), intent(in) :: pt
real(rprec) :: phi_i
integer :: low, nn

nn = size(grid)

! find lower index for pt
low = binary_search(grid,pt)

!perform linear interpolation (or extrapolation if outside bounds of array)
!note: nz_tot is out of domain bounds for uv grid
if (low>nn-2) then
    phi_i = phi(nn-1)
else
    if (low==0) then
        low = 1
    endif
    phi_i = phi(low) + (pt-grid(low))*&
        (phi(low+1)-phi(low))/(grid(low+1)-grid(low))    
endif

end function lin_interp1

!*******************************************************************************
function lin_interp2(phi,grid,pt) result(phi_i)
!*******************************************************************************
use param, only : ny, nz_tot
use functions, only : binary_search

implicit none
real(rprec), dimension(ny,nz_tot) :: phi
real(rprec), dimension(nz_tot) :: grid
real(rprec), dimension(ny) :: phi_i
real(rprec) :: pt
integer :: low

! find lower index for pt
low = binary_search(grid,pt)

!perform linear interpolation (or extrapolation if outside bounds of array)
!note: nz_tot is out of domain bounds for uv grid
if (low>nz_tot-2) then
    phi_i(:) = phi(:,nz_tot-1)
else
    if (low==0) then
        low = 1
    endif
    phi_i(:) = phi(:,low) + (pt-grid(low))*&
        (phi(:,low+1)-phi(:,low))/(grid(low+1)-grid(low))
endif
!if (low>nz_tot-2) then
!    phi_i(:) = phi(:,nz_tot-1)
!elseif (low==0) then
!    phi_i(:) = phi(:,1)
!else
!!    if (low==0) then
!!        low = 1
!!    endif
!    phi_i(:) = phi(:,low) + (pt-grid(low))*&
!        (phi(:,low+1)-phi(:,low))/(grid(low+1)-grid(low))
!endif

end function lin_interp2

!*******************************************************************************
subroutine compute_bl_thickness(this)
!*******************************************************************************
use param, only : nz_tot, dz

implicit none
class(sample_inlet_t), intent(inout) :: this
real(rprec), dimension(nz_tot) :: ubar
real(rprec) :: uinf
integer :: jz

ubar = this%u%avg
uinf = this%uinf
jz = 1
do while (ubar(jz)/uinf<0.99_rprec)
    jz = jz + 1
enddo
! linear interpolation to estimate delta99
this%delta = z_uv(jz-1) + dz*(0.99_rprec*uinf-ubar(jz-1))/(ubar(jz)-ubar(jz-1))

end subroutine compute_bl_thickness

!*******************************************************************************
subroutine compute_momentum_thickness2(this)
!*******************************************************************************
use param

implicit none
class(sample_inlet_t), intent(inout) :: this
real(rprec), dimension(nz_tot) :: ubar
real(rprec) :: utau
real(rprec), dimension(nz_tot) :: dummy
real(rprec) :: uinf, a, Deltap, hwm
integer :: jz,jz_delta

ubar = this%u%avg
uinf = this%uinf
dummy = ubar/uinf*(1._rprec-ubar/uinf)
a = ubar(wmpt)/uinf
hwm = (wmpt-0.5_rprec)*dz
Deltap = hwm*this%utau/nu_molec
this%theta = hwm*a*(&
    (1._rprec-a)*(1._rprec-deltas_fit(Deltap)) + a*theta_fit(Deltap) )
jz = wmpt+1
do while (ubar(jz)<uinf .and. jz < nz_tot)
    this%theta = this%theta + 0.5_rprec*(dummy(jz)+dummy(jz-1))*dz
    jz = jz+1
enddo

end subroutine compute_momentum_thickness2

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
subroutine monitor_rescale_recycle()
!*******************************************************************************
use param, only : path, jt_total, total_time

integer :: fid
character*50 :: fname

fname = path // 'output/sample.dat'
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,&
    total_time,&
    samp%delta,&
    samp%utau,&
    samp%theta,&
    samp%uinf,&
    iter
close(fid)

fname = path // 'output/inlet.dat'
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,&
    total_time,&
    inlt%delta,&
    inlt%utau,&
    inlt%theta,&
    inlt%uinf,&
    iter
close(fid)

end subroutine monitor_rescale_recycle

!*******************************************************************************
subroutine write_sample_velocity_average()
!*******************************************************************************
use param, only : path, ny, nz_tot, jt_total, total_time
use string_util, only : string_splice, string_concat

integer :: fid, jz
character*50 :: fname

fname = trim(adjustl(path)) // trim(adjustl('output/')) // &
    trim(adjustl(samp%loc_name)) // trim(adjustl('.vel-avg.dat'))
open(newunit=fid, file=fname, status='unknown', position='rewind')
write(fid,*) jt_total,&
    total_time,&
    Tavg, &
    samp%utau, &
    samp%tauw, &
    samp%delta
    do jz = 1,nz_tot
        write(fid,*) samp%u%avg(jz),&
        samp%v%avg(jz),&
        samp%w%avg(jz)
    enddo
close(fid)

fname = trim(adjustl(path)) // trim(adjustl('output/')) // &
    trim(adjustl(inlt%loc_name)) // trim(adjustl('.fluc-avg.dat'))
open(newunit=fid, file=fname, status='unknown', position='rewind')
do jz = 1,nz_tot
    write(fid,*) u_fluc_avg(jz),&
    v_fluc_avg(jz),&
    w_fluc_avg(jz)
enddo
close(fid)
if (mod(jt_total,100)==0) then
fname = trim(adjustl(path)) // trim(adjustl('output/')) // &
    trim(adjustl(inlt%loc_name)) // trim(adjustl('.fluc-avg-tot.dat'))
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,sum(u_fluc_avg)/nz_tot,sum(v_fluc_avg)/nz_tot,&
    sum(w_fluc_avg)/nz_tot
endif


end subroutine write_sample_velocity_average

!*******************************************************************************
subroutine read_sample_velocity_average()
!*******************************************************************************
use param, only : coord, path, ny, nz_tot, jt_total, total_time
use string_util, only : string_splice, string_concat

real(rprec) :: dummy
integer :: fid, jz
character*50 :: fname
logical :: file_flag

fname = trim(adjustl(path)) // trim(adjustl('output/')) // &
    trim(adjustl(samp%loc_name)) // trim(adjustl('.vel-avg.dat'))
inquire(file=fname,exist=file_flag)
if (file_flag) then
    if (coord==0) then
        write(*,*) 'reading ',fname
    endif
    open(newunit=fid, file=fname, position='rewind')
    read(fid,*) dummy,&
        dummy,&
        Tavg,&
        samp%utau,&
        samp%tauw,&
        samp%delta
    do jz = 1,nz_tot
        read(fid,*) samp%u%avg(jz),&
            samp%v%avg(jz),&
            samp%w%avg(jz)
    enddo
    close(fid)
endif

fname = trim(adjustl(path)) // trim(adjustl('output/')) // &
    trim(adjustl(inlt%loc_name)) // trim(adjustl('.fluc-avg.dat'))
inquire(file=fname,exist=file_flag)
if (file_flag) then
    if (coord==0) then
        write(*,*) 'reading ',fname
    endif
    open(newunit=fid, file=fname, position='rewind')
    do jz = 1,nz_tot
        read(fid,*) u_fluc_avg(jz),&
            v_fluc_avg(jz),&
            w_fluc_avg(jz)
    enddo
    close(fid)
endif

end subroutine read_sample_velocity_average

!*******************************************************************************
subroutine write_fringe_force()
!*******************************************************************************
use param, only : path, coord, jt_total, nx, ny, nz, write_endian
use sim_param, only : fx,fy,fz
use string_util

character (64) :: fname, bin_ext

call string_splice(fname, path //'output/force.', jt_total)
call string_splice(bin_ext, '.c', coord, '.bin')
call string_concat(fname, bin_ext)
open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) fx(:nx,:ny,1:nz)
write(13,rec=2) fy(:nx,:ny,1:nz)
write(13,rec=3) fz(:nx,:ny,1:nz)
close(13)

end subroutine write_fringe_force

!*******************************************************************************
subroutine check_rescale_velocity()
!*******************************************************************************
use param, only : path, coord, jt_total, nx, ny, nz_tot, write_endian
use string_util

integer :: fid, jz
character (64) :: fname

call string_splice(fname, path //'output/sample.vel-tot.', jt_total)
call string_concat(fname, '.bin')
open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
    access='direct', recl=ny*nz_tot*rprec)
write(13,rec=1) samp%u%tot(1:ny,1:nz_tot)
write(13,rec=2) samp%v%tot(1:ny,1:nz_tot)
write(13,rec=3) samp%w%tot(1:ny,1:nz_tot)
close(13)

call string_splice(fname, path //'output/inlet.vel-tot.', jt_total)
call string_concat(fname, '.bin')
open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
    access='direct', recl=ny*nz_tot*rprec)
write(13,rec=1) inlt%u%tot(1:ny,1:nz_tot)
write(13,rec=2) inlt%v%tot(1:ny,1:nz_tot)
write(13,rec=3) inlt%w%tot(1:ny,1:nz_tot)
close(13)

call string_splice(fname, path //'output/sample.vel-fluc.', jt_total)
call string_concat(fname, '.bin')
open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
    access='direct', recl=ny*nz_tot*rprec)
write(13,rec=1) samp%u%fluc(1:ny,1:nz_tot)
write(13,rec=2) samp%v%fluc(1:ny,1:nz_tot)
write(13,rec=3) samp%w%fluc(1:ny,1:nz_tot)
close(13)

call string_splice(fname, path //'output/inlet.vel-fluc.', jt_total)
call string_concat(fname, '.bin')
open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
    access='direct', recl=ny*nz_tot*rprec)
write(13,rec=1) inlt%u%fluc(1:ny,1:nz_tot)
write(13,rec=2) inlt%v%fluc(1:ny,1:nz_tot)
write(13,rec=3) inlt%w%fluc(1:ny,1:nz_tot)
close(13)

call string_splice(fname, path //'output/sample.vel-avg.', jt_total)
call string_concat(fname, '.dat')
open(newunit=fid, file=fname, status='unknown', position='rewind')
write(fid,*) samp%utau, &
    samp%delta
    do jz = 1,nz_tot
        write(fid,*) samp%u%avg(jz),&
        samp%v%avg(jz),&
        samp%w%avg(jz)
    enddo
close(fid)

call string_splice(fname, path //'output/inlet.vel-avg.', jt_total)
call string_concat(fname, '.dat')
open(newunit=fid, file=fname, status='unknown', position='rewind')
write(fid,*) inlt%utau, &
    inlt%delta
    do jz = 1,nz_tot
        write(fid,*) inlt%u%avg(jz),&
        inlt%v%avg(jz),&
        inlt%w%avg(jz)
    enddo
close(fid)

end subroutine check_rescale_velocity

end module rescale_recycle_fluc
