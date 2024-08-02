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
program main
!*******************************************************************************
!
! Main file for lesgo solver
! Contains main time loop
!

use types, only : rprec
use clock_m
use param
use sim_param
use grid_m
use io, only : energy, output_loop, output_final, jt_total
use io, only : write_tau_wall_bot, write_tau_wall_top
use fft
use derivatives, only : filt_da, ddz_uv, ddz_w
use test_filtermodule
use cfl_util
use sgs_stag_util, only : sgs_stag
use forcing
use functions, only: get_tau_wall_bot, get_tau_wall_top

#ifdef PPMPI
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif

#ifdef PPLVLSET
use level_set, only : level_set_global_CA, level_set_vel_err
use level_set_base, only : global_CA_calc
#endif

#ifdef PPTURBINES
use turbines, only : turbines_forcing, turbine_vel_init
#endif

#ifdef PPSCALARS
use scalars, only : buoyancy_force, scalars_transport, scalars_deriv
#endif

use sponge
use coriolis, only : coriolis_calc, coriolis_forcing, alpha, G, phi_actual
use messages

implicit none

character (*), parameter :: prog_name = 'main'
integer :: nca
character(:), allocatable :: ca

integer :: jt_step, nstart
real(rprec) :: rmsdivvel, ke, maxcfl, tt

type(clock_t) :: clock, clock_total, clock_forcing

! Measure total time in forcing function
real(rprec) :: clock_total_f = 0.0

#ifdef PPMPI
! Buffers used for MPI communication
real(rprec) :: rbuffer
real(rprec) :: maxdummy ! Used to calculate maximum with mpi_allreduce
real(rprec) :: tau_top   ! Used to write top wall stress at first proc
#endif

! Initialize MPI
#ifdef PPMPI
call mpi_init (ierr)
#endif

! Get path if needed
nca = COMMAND_ARGUMENT_COUNT()
if (nca == 1) then
    call GET_COMMAND_ARGUMENT(1, length=nca)
    allocate(character(nca) :: ca)
    call GET_COMMAND_ARGUMENT(1, value=ca)
    allocate(path, source = './' // ca // '/')
else
    allocate(path, source='./')
endif

! Start the clocks, both local and total
call clock%start

! Initialize time variable
tt = 0
jt = 0
jt_total = 0

! Initialize all data
call initialize()

if(coord == 0) then
    call clock%stop
#ifdef PPMPI
    write(*,'(1a,E15.7)') 'Initialization wall time: ', clock % time
#else
    write(*,'(1a,E15.7)') 'Initialization cpu time: ', clock % time
#endif
endif

call clock_total%start

! Initialize starting loop index
! If new simulation jt_total=0 by definition, if restarting jt_total
! provided by total_time.dat
nstart = jt_total+1

! BEGIN TIME LOOP
time_loop: do jt_step = nstart, nsteps

    ! Get the starting time for the iteration
    call clock%start

    if (use_cfl_dt) then

        dt_f = dt
        dt = get_cfl_dt()
        dt_dim = dt * z_i / u_star

        tadv1 = 1._rprec + 0.5_rprec * dt / dt_f
        tadv2 = 1._rprec - tadv1

    end if

   ! Advance time
   jt_total = jt_step
   jt = jt + 1
   total_time = total_time + dt
   total_time_dim = total_time_dim + dt_dim
   tt = tt+dt

    ! Save previous time's right-hand-sides for Adams-Bashforth Integration
    ! NOTE: RHS does not contain the pressure gradient
    RHSx_f = RHSx
    RHSy_f = RHSy
    RHSz_f = RHSz

    ! Calculate velocity derivatives
    ! Calculate dudx, dudy, dvdx, dvdy, dwdx, dwdy (in Fourier space)
    call filt_da(u, dudx, dudy, lbz)
    call filt_da(v, dvdx, dvdy, lbz)
    call filt_da(w, dwdx, dwdy, lbz)

    ! Calculate dudz, dvdz using finite differences (for 1:nz on uv-nodes)
    !  except bottom coord, only 2:nz
    call ddz_uv(u, dudz, lbz)
    call ddz_uv(v, dvdz, lbz)

    ! Calculate dwdz using finite differences (for 0:nz-1 on w-nodes)
    !  except bottom coord, only 1:nz-1
    call ddz_w(w, dwdz, lbz)

#ifdef PPSCALARS
    call scalars_deriv()
#endif

    ! Calculate wall stress and derivatives at the wall
    ! (txz, tyz, dudz, dvdz at jz=1)
    ! using the velocity log-law
    ! MPI: bottom and top processes only
    if (coord == 0 .or. coord == nproc-1) then
        call wallstress()
    end if

    ! Calculate turbulent (subgrid) stress for entire domain
    !   using the model specified in param.f90 (Smag, LASD, etc)
    !   MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz (stress-free lid)
    call sgs_stag()

    ! Exchange ghost node information (since coords overlap) for tau_zz
    !   send info up (from nz-1 below to 0 above)
#ifdef PPMPI
    call mpi_sendrecv (tzz(:,:,nz-1), ld*ny, MPI_RPREC, up, 6,                 &
                       tzz(:,:,0), ld*ny, MPI_RPREC, down, 6,                  &
                       comm, status, ierr)
#endif

    ! Compute divergence of SGS shear stresses
    ! the divt's and the diagonal elements of t are not equivalenced
    ! in this version. Provides divtz 1:nz-1, except 1:nz at top process
    call divstress_uv(divtx, divty, txx, txy, txz, tyy, tyz)
    call divstress_w(divtz, txz, tyz, tzz)

    ! Calculates u x (omega) term in physical space. Uses 3/2 rule for
    ! dealiasing. Stores this term in RHS (right hand side) variable
    call convec()

    ! Add div-tau term to RHS variable
    !   this will be used for pressure calculation
    RHSx(:,:,1:nz-1) = -RHSx(:,:,1:nz-1) - divtx(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = -RHSy(:,:,1:nz-1) - divty(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = -RHSz(:,:,1:nz-1) - divtz(:,:,1:nz-1)
    if (coord == nproc-1) RHSz(:,:,nz) = -RHSz(:,:,nz)-divtz(:,:,nz)

    call coriolis_calc()

#ifdef PPSCALARS
    call scalars_transport()
    call buoyancy_force()
#endif
    call sponge_force()

    !--calculate u^(*) (intermediate vel field)
    !  at this stage, p, dpdx_i are from previous time step
    !  (assumes old dpdx has NOT been added to RHSx_f, etc)
    !  we add force (mean press forcing) here so that u^(*) is as close
    !  to the final velocity as possible
    if (use_mean_p_force) then
        RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + mean_p_force_x
        RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + mean_p_force_y
    end if

    ! Optional random forcing, i.e. to help prevent relaminarization
    if (use_random_force .and. jt_total < stop_random_force) then
        call forcing_random()
    end if

    !//////////////////////////////////////////////////////
    !/// APPLIED FORCES                                 ///
    !//////////////////////////////////////////////////////
    !  In order to save memory the arrays fxa, fya, and fza are now only defined when needed.
    !  For Levelset RNS all three arrays are assigned.
    !  For turbines at the moment only fxa is assigned.
    !  Look in forcing_applied for calculation of forces.
    !  Look in sim_param.f90 for the assignment of the arrays.

    !  Applied forcing (forces are added to RHS{x,y,z})

    ! Calculate forcing time
    call clock_forcing%start

    ! Apply forcing. These forces will later go into RHS
    call forcing_applied()

    ! Calculate forcing time
    call clock_forcing%stop

    ! Calculate the total time of the forcing
    clock_total_f = clock_total_f + clock_forcing % time

    !  Update RHS with applied forcing
#if defined(PPTURBINES) || defined(PPATM)
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + fya(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) + fza(:,:,1:nz-1)
#endif

    !//////////////////////////////////////////////////////
    !/// EULER INTEGRATION CHECK                        ///
    !//////////////////////////////////////////////////////
    ! Set RHS*_f if necessary (first timestep)
    if ((jt_total == 1) .and. (.not. initu)) then
        ! if initu, then this is read from the initialization file
        ! else for the first step put RHS_f=RHS
        !--i.e. at first step, take an Euler step
        RHSx_f = RHSx
        RHSy_f = RHSy
        RHSz_f = RHSz
    end if

    !//////////////////////////////////////////////////////
    !/// INTERMEDIATE VELOCITY                          ///
    !//////////////////////////////////////////////////////
    ! Calculate intermediate velocity field
    !   only 1:nz-1 are valid
    u(:,:,1:nz-1) = u(:,:,1:nz-1) +                                            &
        dt * ( tadv1 * RHSx(:,:,1:nz-1) + tadv2 * RHSx_f(:,:,1:nz-1) )
    v(:,:,1:nz-1) = v(:,:,1:nz-1) +                                            &
        dt * ( tadv1 * RHSy(:,:,1:nz-1) + tadv2 * RHSy_f(:,:,1:nz-1) )
    w(:,:,1:nz-1) = w(:,:,1:nz-1) +                                            &
        dt * ( tadv1 * RHSz(:,:,1:nz-1) + tadv2 * RHSz_f(:,:,1:nz-1) )
    if (coord == nproc-1) then
        w(:,:,nz) = w(:,:,nz) +                                                &
            dt * ( tadv1 * RHSz(:,:,nz) + tadv2 * RHSz_f(:,:,nz) )
    end if

    ! Set unused values to BOGUS so unintended uses will be noticable
#ifdef PPSAFETYMODE
#ifdef PPMPI
    u(:,:,0) = BOGUS
    v(:,:,0) = BOGUS
    w(:,:,0) = BOGUS
#endif
    u(:,:,nz) = BOGUS
    v(:,:,nz) = BOGUS
    if(coord < nproc-1) w(:,:,nz) = BOGUS
#endif

    !//////////////////////////////////////////////////////
    !/// PRESSURE SOLUTION                              ///
    !//////////////////////////////////////////////////////
    ! Solve Poisson equation for pressure
    !   div of momentum eqn + continuity (div-vel=0) yields Poisson eqn
    !   do not need to store p --> only need gradient
    !   provides p, dpdx, dpdy at 0:nz-1 and dpdz at 1:nz-1
    call press_stag_array()

    ! Add pressure gradients to RHS variables (for next time step)
    !   could avoid storing pressure gradients - add directly to RHS
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) - dpdx(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) - dpdy(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) - dpdz(:,:,1:nz-1)
    if(coord == nproc-1) then
        RHSz(:,:,nz) = RHSz(:,:,nz) - dpdz(:,:,nz)
    end if

    !//////////////////////////////////////////////////////
    !/// INDUCED FORCES                                 ///
    !//////////////////////////////////////////////////////
    ! Calculate external forces induced forces. These are
    ! stored in fx,fy,fz arrays. We are calling induced
    ! forces before applied forces as some of the applied
    ! forces (RNS) depend on the induced forces and the
    ! two are assumed independent
    call forcing_induced()

    !//////////////////////////////////////////////////////
    !/// PROJECTION STEP                                ///
    !//////////////////////////////////////////////////////
    ! Projection method provides u,v,w for jz=1:nz
    !   uses fx,fy,fz calculated above
    !   for MPI: syncs 1 -> Nz and Nz-1 -> 0 nodes info for u,v,w
    call project ()

    ! Write ke to file
    if (modulo (jt_total, nenergy) == 0) call energy(ke)

#ifdef PPLVLSET
    if (global_CA_calc) call level_set_global_CA()
#endif

    ! Write output files
    call output_loop()

    ! Check the total time of the simulation up to this point on the master
    ! node and send this to all

    if (modulo (jt_total, wbase) == 0) then

        ! Get the ending time for the iteration
        call clock%stop
        call clock_total%stop

        ! Calculate rms divergence of velocity
        ! only written to screen, not used otherwise
        call rmsdiv(rmsdivvel)
        maxcfl = get_max_cfl()

        ! This takes care of the clock times, to obtain the quantities based
        ! on all the processors, not just processor 0
#ifdef PPMPI
        call mpi_allreduce(clock % time, maxdummy, 1, mpi_rprec,               &
            MPI_MAX, comm, ierr)
        clock % time = maxdummy
        call mpi_allreduce(clock_total % time, maxdummy, 1, mpi_rprec,         &
            MPI_MAX, comm, ierr)
        clock_total % time = maxdummy
        call mpi_allreduce(clock_forcing % time, maxdummy, 1, mpi_rprec,       &
            MPI_MAX, comm, ierr)
        clock_forcing % time = maxdummy
        call mpi_allreduce(clock_total_f , maxdummy, 1, mpi_rprec,             &
            MPI_MAX, comm, ierr)
        clock_total_f = maxdummy
#endif

        ! Send top wall stress to bottom process
#ifdef PPMPI
        if (coord == nproc-1) then
            tau_top = get_tau_wall_top()
        else
            tau_top = 0._rprec
        endif

        call mpi_allreduce(tau_top, maxdummy, 1, mpi_rprec,               &
            MPI_SUM, comm, ierr)
        tau_top = maxdummy
#endif

            if (coord == 0) then
            write(*,*)
            write(*,'(a)') '==================================================='
            write(*,'(a)') 'Time step information:'
            write(*,'(a,i9)') '  Iteration: ', jt_total
            write(*,'(a,E15.7)') '  Time step: ', dt
            write(*,'(a,E15.7)') '  Dimensional time: ', total_time_dim
            write(*,'(a,E15.7)') '  CFL: ', maxcfl
            write(*,'(a,2E15.7)') '  AB2 TADV1, TADV2: ', tadv1, tadv2
            write(*,*)
            write(*,'(a)') 'Flow field information:'
            write(*,'(a,E15.7)') '  Velocity divergence metric: ', rmsdivvel
            write(*,'(a,E15.7)') '  Kinetic energy: ', ke
            write(*,'(a,E15.7)') '  Bot wall stress: ', get_tau_wall_bot()
#ifdef PPMPI
            write(*,'(a,E15.7)') '  Top wall stress: ', tau_top
#else
            write(*,'(a,E15.7)') '  Top wall stress: ', get_tau_wall_top()
#endif
            write(*,*)
            write(*,'(1a)') 'Simulation wall times (s): '
            write(*,'(1a,E15.7)') '  Iteration: ', clock % time
            write(*,'(1a,E15.7)') '  Cumulative: ', clock_total % time
            write(*,'(1a,E15.7)') '  Forcing: ', clock_forcing % time
            write(*,'(1a,E15.7)') '  Cumulative Forcing: ', clock_total_f
            write(*,'(1a,E15.7)') '  Forcing %: ',                             &
                clock_total_f /clock_total % time
            if (coriolis_forcing > 0) then
                write(*,*)
                write(*,'(1a)') 'Coriolis parameters: '
                write(*,'(1a,E15.7)') '  G: ', G
                write(*,'(1a,E15.7)') '  alpha: ', alpha
                if (coriolis_forcing == 2) then
                    write(*,'(1a,E15.7)') '  wind direction: ', phi_actual
                end if
            end if
            write(*,'(a)') '==================================================='
            call write_tau_wall_bot()
        end if
        if(coord == nproc-1) then
            call write_tau_wall_top()
        end if

#ifdef PPMPI
        call mpi_barrier(comm, ierr)
#endif

        ! Check if we are to check the allowable runtime
        if (runtime > 0) then

#ifdef PPMPI
            ! Determine the processor that has used most time and communicate
            ! this. Needed to make sure that all processors stop at the same
            ! time and not just some of them
            call mpi_allreduce(clock_total % time, rbuffer, 1, MPI_RPREC,      &
                MPI_MAX, MPI_COMM_WORLD, ierr)
            clock_total % time = rbuffer
#endif

            ! If maximum time is surpassed go to the end of the program
            if ( clock_total % time >= real(runtime,rprec) ) then
                call mesg( prog_name,                                          &
                    'Specified runtime exceeded. Exiting simulation.')
                exit time_loop
            endif

       endif

    end if

end do time_loop
! END TIME LOOP

! Finalize
close(2)

! Write total_time.dat and tavg files
call output_final()

! Stop wall clock
call clock_total%stop
#ifdef PPMPI
if (coord == 0) write(*,"(a,e15.7)") 'Simulation wall time (s) : ',            &
    clock_total % time
#else
if (coord == 0) write(*,"(a,e15.7)") 'Simulation cpu time (s) : ',             &
    clock_total % time
#endif

call finalize()

if(coord == 0 ) write(*,'(a)') 'Simulation complete'

end program main
