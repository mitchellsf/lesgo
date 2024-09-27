!!
!!  Copyright (C) 2010-2016  Johns Hopkins University
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
module adjust_mean_dpdx
!*******************************************************************************

use types, only : rprec
use param
use grid_m
use messages
use string_util
#ifdef PPMPI
use mpi_defs, only : MPI_SYNC_DOWNUP, mpi_sync_real_array 
#endif
use sim_param, only : u 

implicit none

save
private

public :: adjust_dpdx

real(rprec) :: ubulk_old

! Commonly used indices
integer :: i,j,k

contains

!*******************************************************************************
subroutine adjust_dpdx()
!*******************************************************************************

implicit none

real(rprec) :: ubulk_target, ubulk
real(rprec) :: ubulk_global
real(rprec) :: Rem, Retau

if (lin_accel) then

    ! initialize if restarting
    if (jt == 1 .and. jt_total /= 1) then
       call read_restart
    endif

    if (jt_total==lin_accel_nstart) then
        lin_accel_tstart = total_time
    endif

    if (jt_total < lin_accel_nstart) then
        mean_p_force_x = 1._rprec
    elseif (total_time-lin_accel_tstart <= lin_accel_time) then
        Rem = (total_time-lin_accel_tstart)*lin_accel_dRemdt &
            + lin_accel_Rem_i
        ! Dean's correlation
        Retau = 0.0955248659_rprec*Rem**0.875_rprec
        mean_p_force_x = 0.5_rprec/lin_accel_Retau_i*lin_accel_dRemdt &
            + (Retau/lin_accel_Retau_i)**2
    else
        mean_p_force_x = (lin_accel_Retau_f/lin_accel_Retau_i)**2
    endif

    ubulk = 0.0_rprec
    do k=1, nz-1
       do j=1, ny
          do i= 1, nx
             ubulk = ubulk + u(i,j,k)*dz
          enddo
       enddo
    enddo
    ubulk = ubulk / (real(nx)*real(ny)*L_z)
    
    call mpi_allreduce(ubulk,ubulk_global,1,MPI_RPREC,MPI_SUM, &
         MPI_COMM_WORLD,ierr)

    ubulk_old = ubulk_global

    if (coord == 0 .and. mod(jt_total,100)==0) then
       call write_mean_dpdx_mean_velocity
    endif

    if (jt_total == nsteps .and. coord == 0) then
       call write_restart
    endif

elseif (pulse) then

   call pulse_flow 

else

    ubulk_target = 1.0_rprec
    
    ubulk = 0.0_rprec
    do k=1, nz-1
       do j=1, ny
          do i= 1, nx
             ubulk = ubulk + u(i,j,k)*dz
          enddo
       enddo
    enddo
    ubulk = ubulk / (real(nx)*real(ny)*L_z)
    
    call mpi_allreduce(ubulk,ubulk_global,1,MPI_RPREC,MPI_SUM, &
         MPI_COMM_WORLD,ierr)
    
    ! initialize if restarting
    if (jt == 1 .and. jt_total /= 1) then
       call read_restart
    endif
    
    if (jt_total == 1) then
       mean_p_force_x = 0._rprec
       ubulk_old = ubulk_global
    endif
    
    mean_p_force_x = mean_p_force_x - ( 1.5_rprec * ( ubulk_global - ubulk_target ) / dt &
         - 0.5_rprec * ( ubulk_old - ubulk_target ) / dt )
    
    ubulk_old = ubulk_global
    
    if (jt_total == nsteps .and. coord == 0) then
       call write_restart
    endif
    
endif

end subroutine adjust_dpdx

!*******************************************************************************
subroutine pulse_flow()
!*******************************************************************************

implicit none

real(rprec) :: Re0, beta

! initialize if restarting
if (jt == 1 .and. jt_total /= 1) then
    call read_restart
endif

if (jt_total==pulse_nstart) then
    pulse_tstart = total_time
endif

if (jt_total > pulse_nstart) then
    Re0 = 9078._rprec ! constant for now, should have Re dependence
    beta = pulse_amp*pulse_freq*Re0
    mean_p_force_x = 1._rprec + &
        beta*sin(pulse_freq*(total_time-pulse_tstart)/nu_molec)
endif

if (jt_total == nsteps .and. coord == 0) then
   call write_restart
endif

end subroutine pulse_flow

!*******************************************************************************
subroutine write_restart()
!*******************************************************************************

open(2,file='adjust_mean_dpdx_restart.txt')
if (lin_accel) then
    write(2,*) lin_accel_tstart
elseif(pulse) then
    write(2,*) pulse_tstart
else
    write(2,*) ubulk_old, mean_p_force_x
endif
close(2)

end subroutine write_restart

!*******************************************************************************
subroutine write_mean_dpdx_mean_velocity()
!*******************************************************************************

open(2,file='mean_dpdx_mean_velocity.dat',position="append")
write(2,*) jt_total, total_time, ubulk_old, mean_p_force_x
close(2)

end subroutine write_mean_dpdx_mean_velocity

!*******************************************************************************
subroutine read_restart()
!*******************************************************************************

open(2,file='adjust_mean_dpdx_restart.txt')
if (lin_accel) then
    read(2,*) lin_accel_tstart
elseif(pulse) then
    read(2,*) pulse_tstart
else
    read(2,*) ubulk_old, mean_p_force_x
endif
close(2)

end subroutine read_restart


end module adjust_mean_dpdx

