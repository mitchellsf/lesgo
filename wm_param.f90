!!
!!  Copyright (C) 2011-2017  Johns Hopkins University
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
module wm_param
!*******************************************************************************
use types, only : rprec

save
private rprec
public

! variables for all wall models
! wall stress from quasi-equilibrium and non-equilibrium models
real(rprec), dimension(:,:), allocatable :: twxbar, twybar, twxpp, twypp,&
    twxpp_turb,twypp_turb
real(rprec), dimension(:,:), allocatable :: twx_eq, twy_eq
! velocities at first grid point (non-filtered)
real(rprec), dimension(:,:), allocatable :: uinst, vinst, Ud
! velocity angle(measured from x-axis): d = Deltay
real(rprec), dimension(:,:), allocatable :: theta_d
! pressure gradients at first grid point
real(rprec), dimension(:,:), allocatable :: dplesdx, dplesdy, dpdxbar1, dpdybar1,&
    dpdxbar2, dpdybar2, dpdxpp, dpdypp, dpdxpp_m, dpdypp_m, dpdxpp_mm, dpdypp_mm
! wall model height
real(rprec) :: Deltay
! temporary files for debug_wm
real(rprec), dimension(:,:,:), allocatable :: utemp2,vtemp2,dpdxtemp,dpdytemp,&
    utxtemp, utytemp
! time scale for filtering the pressure
real(rprec), dimension(:,:), allocatable :: Tnu
! filtered LES velocity
real(rprec), dimension(:,:), allocatable :: ubar, vbar
real(rprec), dimension(:,:), allocatable :: uinf, vinf
real(rprec), dimension(:,:), allocatable :: uinf_lam, vinf_lam
real(rprec), dimension(:,:), allocatable :: Tdelta
real(rprec), dimension(:,:), allocatable :: dpdxbar3,dpdybar3
real(rprec), dimension(:,:), allocatable :: dpdxppdelta,dpdyppdelta
real(rprec) :: ls=12._rprec ! thickness of Stokes layer in inner units
real(rprec), dimension(:,:), allocatable :: dpdx_fit,dpdy_fit
real(rprec), dimension(:,:), allocatable :: redelta, psi_p
real(rprec) :: PI=4*atan(1._rprec)

!LaRTE wall model variables
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

! lamNEQ variables
! sum-of-exponentials approximation coefficients
real(rprec), dimension(:), allocatable :: ws, ss
! integrals for history part of non-equilibrium wall model
real(rprec), dimension(:,:,:), allocatable :: Ix, Iy
! number of exponential terms
integer :: nexp=0

end module wm_param
