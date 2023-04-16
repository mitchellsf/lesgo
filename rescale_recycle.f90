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
module rescale_recycle
!*******************************************************************************
use types, only : rprec
use fringe

implicit none

private
public rescale_recycle_init, rescale_recycle_calc, compute_momentum_thickness

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
    real(rprec) :: delta, theta, utau, uinf
    character*50 :: loc_name
end type sample_inlet_t

real(rprec), dimension(:), allocatable :: z_uv, z_w
type(sample_inlet_t) :: inlt, samp
type(fringe_t) :: apply_fringe

contains
!*******************************************************************************
subroutine rescale_recycle_init
!*******************************************************************************
use param, only : fringe_region_end, fringe_region_len, sampling_region_end
use param, only : nx, ny, nz
use param, only : nz_tot, dz

integer :: jz

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

end subroutine rescale_recycle_init

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
subroutine rescale_recycle_calc
!*******************************************************************************
use param
use sim_param, only : u, v, w, txz

integer :: i, i_w, kstart, kend, jz
real(rprec), dimension(ny,nz_tot) :: dummy1, dummy2, dummy3
type(sample_inlet_t) :: inlt_eps
real(rprec) :: tol, eps

call sample_inlet_init(inlt_eps)
tol = 0.000001_rprec
eps = 0.000001_rprec

dummy1 = 0._rprec
dummy2 = 0._rprec
dummy3 = 0._rprec

! Sampling plane
kstart = coord*(nz-1) + 1
kend = (coord+1)*(nz-1)
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
samp%u%avg = sum(samp%u%tot,1)/ny
samp%v%avg = sum(samp%v%tot,1)/ny
samp%w%avg = sum(samp%w%tot,1)/ny
samp%uinf = samp%u%avg(nz_tot-1)
do jz = 1,nz_tot
    samp%u%fluc(:,jz) = samp%u%tot(:,jz) - samp%u%avg(jz)
    samp%v%fluc(:,jz) = samp%v%tot(:,jz) - samp%v%avg(jz)
    samp%w%fluc(:,jz) = samp%w%tot(:,jz) - samp%w%avg(jz)
enddo
samp%utau = sqrt(sum(sqrt(txz(isample,1:ny,1)**2.0+txz(isample,1:ny,1)**2.0))/ny)
call mpi_bcast(samp%utau,1,MPI_RPREC,0,MPI_COMM_WORLD,ierr)
call compute_momentum_thickness2(samp%u%avg,samp%theta,samp%utau)
call compute_bl_thickness(samp%u%avg,samp%delta)

! Inlet plane computation (sampling plane rescaled)
!call rescale_velocity(samp,inlt)
call iterate_utau(samp,inlt)
iter = 1
do while (abs(inlt%theta-1._rprec) .gt. tol .and. iter .lt. 100)
    iter = iter + 1
    inlt_eps%delta = inlt%delta + eps
!    call rescale_velocity(samp,inlt_eps)
    call iterate_utau(samp,inlt_eps)
    inlt%delta = inlt%delta - (inlt%theta-1._rprec)/(inlt_eps%theta-inlt%theta)*eps
!    call rescale_velocity(samp,inlt)
    if (inlt%delta < 0._rprec .or. inlt%delta > L_z) then
        inlt%delta = 5._rprec
    endif
    call iterate_utau(samp,inlt)
enddo

if (coord==0 .and. mod(jt_total,1)==0) then
    call monitor_velocity_rescale_recycle(samp)
    call monitor_velocity_rescale_recycle(inlt)
endif

if (coord==0 .and. mod(jt_total,1)==0) then
    call monitor_rescale_recycle()
endif

! Apply inflow conditions
do i = 1, apply_fringe%nx
    i_w = apply_fringe%iwrap(i)
    u(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * u(i_w,1:ny,1:nz)&
        + apply_fringe%beta(i) * inlt%u%tot(1:ny,kstart:kend+1)
    v(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * v(i_w,1:ny,1:nz)&
        + apply_fringe%beta(i) * inlt%v%tot(1:ny,kstart:kend+1)
    w(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * w(i_w,1:ny,1:nz)&
        + apply_fringe%beta(i) * inlt%w%tot(1:ny,kstart:kend+1)
end do

end subroutine rescale_recycle_calc

!*******************************************************************************
subroutine rescale_velocity(sam,inl)
!*******************************************************************************
use param

implicit none
class(sample_inlet_t), intent(in) :: sam
class(sample_inlet_t), intent(inout) :: inl
real(rprec), dimension(nz_tot) :: ubar_inner, ubar_outer
real(rprec), dimension(ny,nz_tot) :: up_inner, up_outer, vp_inner, vp_outer, &
    wp_inner, wp_outer
real(rprec) :: gg, a, b, weight_uv, weight_w
integer :: jz


a = 4._rprec
b = 0.2_rprec
gg = inl%utau/sam%utau
do jz = 1,nz_tot
!    ubar_inner(jz) = inl%utau*f_fit(z_uv(jz)*inl%utau/nu_molec)
    ubar_inner(jz) = gg*lin_interp1(sam%u%avg,z_uv*sam%utau/nu_molec,&
        z_uv(jz)*inl%utau/nu_molec)
    ubar_outer(jz) = gg*lin_interp1(sam%u%avg,z_uv/sam%delta,&
        z_uv(jz)/inl%delta) + (1._rprec-gg*sam%uinf)*(1._rprec)  
!    ubar_outer(jz) = gg*lin_interp1(sam%u%avg,z_uv/sam%delta,&
!        z_uv(jz)/inl%delta) + (1._rprec-gg)*(1._rprec)
    up_inner(:,jz) = gg*lin_interp2(sam%u%fluc,z_uv*sam%utau/nu_molec,&
        z_uv(jz)*inl%utau/nu_molec)
    vp_inner(:,jz) = gg*lin_interp2(sam%v%fluc,z_uv*sam%utau/nu_molec,&
        z_uv(jz)*inl%utau/nu_molec)
    wp_inner(:,jz) = gg*lin_interp2(sam%w%fluc,z_w*sam%utau/nu_molec,&
        z_w(jz)*inl%utau/nu_molec)
    up_outer(:,jz) = gg*lin_interp2(sam%u%fluc,z_uv/sam%delta,z_uv(jz)/inl%delta)
    vp_outer(:,jz) = gg*lin_interp2(sam%v%fluc,z_uv/sam%delta,z_uv(jz)/inl%delta)
    wp_outer(:,jz) = gg*lin_interp2(sam%w%fluc,z_w/sam%delta,z_w(jz)/inl%delta)
    weight_uv = 0.5_rprec*(1.0+tanh( a*(z_uv(jz)/inl%delta-b)/&
        ((1.0-2.0*b)*z_uv(jz)/inl%delta+b) )/tanh(a))
    weight_w = 0.5_rprec*(1.0+tanh( a*(z_w(jz)/inl%delta-b)/&
        ((1.0-2.0*b)*z_w(jz)/inl%delta+b) )/tanh(a))
    if (z_uv(jz)>inl%delta) then
        weight_uv = 1._rprec
        if (z_uv(jz)<2.0*inl%delta) then
            up_outer(:,jz) = up_outer(:,jz)*cos((z_uv(jz)/inl%delta-1.0)*pi/2.0)
            vp_outer(:,jz) = vp_outer(:,jz)*cos((z_uv(jz)/inl%delta-1.0)*pi/2.0)
        else
            up_outer(:,jz) = 0._rprec
            vp_outer(:,jz) = 0._rprec
        endif
    endif
    if (z_w(jz)>inl%delta) then
        weight_w = 1._rprec
        if (z_w(jz)<2.0*inl%delta) then
            wp_outer(:,jz) = wp_outer(:,jz)*cos((z_w(jz)/inl%delta-1.0)*pi/2.0)
        else
            wp_outer(:,jz) = 0._rprec
        endif
    endif

    inl%u%tot(:,jz) = (ubar_inner(jz)+up_inner(:,jz))*(1._rprec-weight_uv) +&
        (ubar_outer(jz)+up_outer(:,jz))*weight_uv
    inl%v%tot(:,jz) = (vp_inner(:,jz))*(1._rprec-weight_uv) +&
        (vp_outer(:,jz))*weight_uv
    inl%w%tot(:,jz) = (wp_inner(:,jz))*(1._rprec-weight_w) +&
        (wp_outer(:,jz))*weight_w
enddo

inl%u%avg = sum(inl%u%tot,1)/ny
inl%uinf = inl%u%avg(nz_tot-1)
call compute_momentum_thickness2(inl%u%avg,inl%theta,inl%utau)
!inl%utau = retd_fit(inl%u%avg(1)*z_uv(1)/nu_molec)*nu_molec/z_uv(1)

end subroutine rescale_velocity

!*******************************************************************************
subroutine iterate_utau(samp,inlt)
!*******************************************************************************
use param, only : coord, nu_molec

implicit none
real(rprec) :: tol, utau_before, utau_after, retd_inlt
type(sample_inlet_t), intent(inout) :: samp, inlt
integer :: iter

!tol = 0.01_rprec
!utau_before = 1._rprec
!utau_after = -1._rprec
!iter = 1
!do while (abs((utau_before-utau_after)/utau_after) .gt. tol .and. iter .lt. 1000)
!    utau_before = inlt%utau
!    call rescale_velocity(samp,inlt)
!    utau_after = inlt%utau
!    iter = iter + 1
!enddo
!if (coord==0) then
!    write(*,*) iter
!endif

inlt%utau = samp%utau*(samp%theta/1._rprec)**0.125_rprec
!inlt%utau = samp%utau*(samp%theta/inlt%theta)**0.125_rprec
!inlt%utau = 0.9_rprec*samp%utau
!retd_inlt = 0.2988_rprec/nu_molec + 70.2624_rprec
!inlt%utau = retd_inlt*nu_molec/inlt%delta
call rescale_velocity(samp,inlt)

end subroutine iterate_utau

!*******************************************************************************
function retd_fit(red) result(retd)
!*******************************************************************************

implicit none
real(rprec) :: red, retd
real(rprec) :: kappa3, beta1, beta2, kappa4

kappa3 = 0.005_rprec
beta1 = (1._rprec + 0.155_rprec*red**-0.03_rprec)**-1
beta2 = 1.7 - (1._rprec + 36._rprec*red**-0.75_rprec)**-1
kappa4 = kappa3**(beta1 - 0.5_rprec)
retd = kappa4*red**beta1*                            &
    (1._rprec + (kappa3*red)**-beta2)**((beta1-0.5_rprec)/beta2)

end function retd_fit

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

end function lin_interp2

!*******************************************************************************
subroutine compute_momentum_thickness(ubar,theta)
!*******************************************************************************
use param

implicit none
real(rprec) :: theta_local
real(rprec), dimension(nz), intent(in) :: ubar
real(rprec), intent(out) :: theta
real(rprec), dimension(nz) :: dummy
integer :: jz

dummy = ubar*(1._rprec-ubar)
theta_local = 0.25_rprec*dummy(1)*dz
do jz = 2, nz-1
    theta_local = theta_local + 0.5_rprec*(dummy(jz)+dummy(jz-1))*dz
enddo
call mpi_allreduce(theta_local,theta,1,MPI_RPREC,MPI_SUM, &
     MPI_COMM_WORLD,ierr)

end subroutine compute_momentum_thickness

!*******************************************************************************
subroutine compute_momentum_thickness2(ubar,theta,utau)
!*******************************************************************************
use param

implicit none
real(rprec), dimension(nz_tot), intent(in) :: ubar
real(rprec), intent(in) :: utau
real(rprec), intent(out) :: theta
real(rprec), dimension(nz_tot) :: dummy
real(rprec) :: uinf, a, Deltap, hwm
integer :: jz

uinf = ubar(nz_tot-1)
dummy = ubar/uinf*(1._rprec-ubar/uinf)
!theta = 0.25_rprec*dummy(1)*dz
a = ubar(wmpt)/uinf
hwm = (wmpt-0.5_rprec)*dz
Deltap = hwm*utau/nu_molec
theta = hwm*a*(&
    (1._rprec-a)*(1._rprec-deltas_fit(Deltap)) + a*theta_fit(Deltap) );
do jz = wmpt+1, nz_tot-1
    theta = theta + 0.5_rprec*(dummy(jz)+dummy(jz-1))*dz
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
subroutine compute_bl_thickness(ubar,delta)
!*******************************************************************************
use param

implicit none
real(rprec), dimension(nz_tot), intent(in) :: ubar
real(rprec), intent(out) :: delta
real(rprec) :: dz_f, z_f, ubar_f
integer :: jz

dz_f = dz/10._rprec
jz = 1
z_f = jz*dz_f
ubar_f = lin_interp1(ubar,z_uv,z_f)
do while (ubar_f/ubar(nz_tot-1)<0.99_rprec)
    jz = jz + 1
    z_f = jz*dz_f
    ubar_f = lin_interp1(ubar,z_uv,z_f)
enddo
delta = z_f

end subroutine compute_bl_thickness

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
subroutine monitor_velocity_rescale_recycle(this)
!*******************************************************************************
use param, only : path, ny, nz_tot, jt_total, total_time, write_endian
use string_util, only : string_splice, string_concat

class(sample_inlet_t), intent(in) :: this
integer :: fid
character*50 :: fname

!call string_splice(fname,'.vel-tot.', jt_total)
!fname = trim(adjustl(path)) // trim(adjustl('output/')) // &
!    trim(adjustl(this%loc_name)) // trim(adjustl(fname))
!open(unit=fid, file=fname, form='unformatted', convert=write_endian,        &
!    access='direct', recl=ny*nz_tot*rprec)
!write(fid,rec=1) this%u%tot
!write(fid,rec=2) this%v%tot
!write(fid,rec=3) this%w%tot
!close(fid)

!call string_splice(fname,'.vel-avg.', jt_total)
!fname = trim(adjustl(path)) // trim(adjustl('output/')) // &
!    trim(adjustl(this%loc_name)) // trim(adjustl(fname))
!open(unit=fid, file=fname, form='unformatted', convert=write_endian,        &
!    access='direct', recl=nz_tot*rprec)
!write(fid,rec=1) this%u%avg
!close(fid)

fname = trim(adjustl(path)) // trim(adjustl('output/')) // &
    trim(adjustl(this%loc_name)) // trim(adjustl('.vel-avg.dat'))
open(newunit=fid, file=fname, status='unknown', position='append')
write(fid,*) jt_total,&
    total_time,&
    this%u%avg
close(fid)

end subroutine monitor_velocity_rescale_recycle

end module rescale_recycle
