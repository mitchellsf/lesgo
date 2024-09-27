cmake_minimum_required(VERSION 2.8.9)

##############################################################################
#                               User Inputs
##############################################################################

# Specify Compiler
# Options: Intel and GNU
set(CMAKE_Fortran_COMPILER_ID "GNU")

# Set the name of the computer you are compiling on.
# Supported computers: rockfish, panoramix, njord, marcc, darwin-brew
# Otherwise, default settings are used
set(hostname "personal")

# Read and write binary files as big_endian or little_endian
# Options: DEFAULT, LITTLE, or BIG
set(WRITE_ENDIAN "DEFAULT")
set(READ_ENDIAN "DEFAULT")

# Options to enable various features in the code
option(USE_MPI "Enable parallelization support." ON)
option(USE_CPS "Enable concurrent precursor simulations." OFF)
option(USE_HIT "Enable homogeneous isotropic turbulence input." OFF)
option(USE_STREAKS "Shift the domain in the y (spanwise) direction to eliminate streaks effects on averaging." OFF)
option(USE_LVLSET "Enable level set." OFF)
option(USE_DYN_TN "Dynamically update T, the timescale for Lagrangian averaging, using Taylor timescale." OFF)
option(USE_SAFETYMODE "Set extra safety measures." ON)
option(USE_CGNS "Enable CGNS data output." OFF)
option(USE_MAPPING "Enable grid mapping transformation." OFF)
option(USE_LVLSET_STRETCH "Enable simplified level set with mapping." OFF)
option(USE_CNDIFF "Enable Crank-Nicolson time marching for wall-normal diffusion." ON)

option(OUTPUT_SGS "Output time-averaging statistics on SGS model." OFF)
option(OUTPUT_BUDGET "Output time-averaging statistics needed for budgets." OFF)
option(OUTPUT_TURBSPEC "Output time-averaging spectra data." OFF)
option(OUTPUT_SPECBUDG "Output time-averaging spectral budget data." OFF)
option(OUTPUT_VELVORT "Output additional velocity/vorticity statistics." OFF)
option(OUTPUT_WMLES "Output time-averaging statistics on wall model." OFF)
option(OUTPUT_SSP "Output self-sustaining process analysis." OFF)
option(OUTPUT_CLOCK "Output time data for each simulation step." OFF)
option(USE_BL "Enable boundary layer simulations." ON)

# Enable verbose mode
set(CMAKE_VERBOSE_MAKEFILE ON)

##############################################################################
#                         Set up build (do not change)
##############################################################################
# Set compiler versions for Cmake
#  Intel
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(compiler_name "ifort")
#  GNU
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(compiler_name "gfortran")
endif ()
find_program(WHICH NAMES which)
exec_program(${WHICH} ARGS ${compiler_name} OUTPUT_VARIABLE compiler_path)
set(CMAKE_Fortran_COMPILER ${compiler_path})

# Project details
enable_language (Fortran)
project (lesgo)
set(exec_name "lesgo")

# Library links
if (${hostname} STREQUAL "marcc")
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} -mkl)
    else()
        link_libraries(fftw3)
    endif()
else()
    link_libraries(fftw3)
endif()

if (${hostname} STREQUAL "rockfish")
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        set(CMAKE_Fortran_COMPILER ${compiler_path})
	set(CMAKE_C_COMPILER gcc)
	set(CMAKE_CXX_COMPILER g++)
    else()
	link_libraries(fftw3)
    endif()
endif()

# Common include and link directories
include_directories(/usr/local/include)
link_directories(/usr/local/lib)
if ( NOT(${hostname} STREQUAL "darwin-brew") )
    # include_directories($ENV{HOME}/include /usr/include)
    include_directories(/usr/include)
    link_directories($ENV{HOME}/lib /usr/lib64)
endif ()

# Include host-specific paths
#  panoramix
if (${hostname} STREQUAL "panoramix")
    include_directories(/opt/cgns-3.2.1-hdf5-mpi/include
    /opt/hdf5-1.8.12-mpi/include /usr/local/hdf5/include /usr/local/cgns/include
    /act/fftw3/mpich/gcc-4.7.2/include /act/hdf5-1.8.15/mpich/gcc-4.7.2/include
    /act/cgns-3.3.0/mpich/gcc-4.7.2/include)
    link_directories(/opt/hdf5-1.8.12-mpi/lib  /opt/cgns-3.2.1-hdf5-mpi/lib
    /usr/local/hdf5/lib  /usr/local/cgns/lib /act/fftw3/mpich/gcc-4.7.2/lib
    /act/hdf5-1.8.15/mpich/gcc-4.7.2/lib /act/cgns-3.3.0/mpich/gcc-4.7.2/lib)
#  marcc
elseif (${hostname} STREQUAL "marcc")
    # Intel
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        include_directories(/software/apps/compilers/intel/mkl/include/fftw/)
        link_directories(/software/apps/compilers/intel/mkl/lib/intel64/)
    # GNU -- These have to be installed locally for now.
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        include_directories(/software/apps/fftw3/3.3.8/gcc/5.5.0/include/)
        link_directories(/software/apps/fftw3/3.3.8/gcc/5.5.0/lib/)
        if (USE_CGNS)
            include_directories($ENV{HOME}/cgns/include/)
            link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
        endif ()
    endif ()
#  rockfish
elseif (${hostname} STREQUAL "rockfish")
    # Intel
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        include_directories(/data/apps/linux-centos8-cascadelake/intel-19.1.2.254/fftw-3.3.8-wxxony3sfn4vctui2jomkwlzs6jfnvqr/include/)
        link_directories(/data/apps/linux-centos8-cascadelake/intel-19.1.2.254/fftw-3.3.8-wxxony3sfn4vctui2jomkwlzs6jfnvqr/lib/)
    # GNU -- These have to be installed locally for now.
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        include_directories(/data/apps/linux-centos8-cascadelake/gcc-9.3.0/fftw-3.3.8-uxkh4w5nrgeviv2ooaorhcuolbtbosyi/include/)
        link_directories(/data/apps/linux-centos8-cascadelake/gcc-9.3.0/fftw-3.3.8-uxkh4w5nrgeviv2ooaorhcuolbtbosyi/lib/)
        if (USE_CGNS)
            include_directories($ENV{HOME}/cgns/include/)
            link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
        endif ()
    endif ()
# Cheyenne
elseif (${hostname} STREQUAL "cheyenne")
    # Intel
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        include_directories(/glade/u/apps/ch/opt/fftw3-mpi/3.3.6-pl2/mpt/2.15f/intel/18.0.1/include/)
        link_directories(/glade/u/apps/ch/opt/fftw3-mpi/3.3.6-pl2/mpt/2.15f/intel/18.0.1/lib/)
    # GNU -- These have to be installed locally for now.
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        include_directories(/software/apps/fftw3/3.3.8/gcc/5.5.0/include/)
        link_directories(/software/apps/fftw3/3.3.8/gcc/5.5.0/lib/)
        if (USE_CGNS)
            include_directories($ENV{HOME}/cgns/include/)
            link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
        endif ()
    endif ()
# stampede
elseif (${hostname} STREQUAL "stampede")
    # Intel
        include_directories(/opt/apps/intel18/impi18_0/fftw3/3.3.8/include /opt/intel/mkl/include/intel64/lp64)
        link_directories(/opt/apps/intel18/impi18_0/fftw3/3.3.8/lib /opt/intel/mkl/lib/intel64_lin)
#  njord
elseif (${hostname} STREQUAL "njord")
    # nothing yet, but can add if needed
# Personal computer
elseif (${hostname} STREQUAL "personal")
    # GNU
        include_directories(/usr/include/)
        link_directories(/usr/lib/x86_64-linux-gnu/)
#  all other systems
elseif (${hostname} STREQUAL "default")
    # include_directories($ENV{HOME}/codes/fftw3/include)
    # link_directories($ENV{HOME}/codes/fftw3/lib)
    # nothing yet, but can add if needed
    if (USE_CGNS)
        include_directories($ENV{HOME}/cgns/include/)
        link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
    endif ()
endif ()

# List common sources
set (Sources cfl_util.f90 clocks.f90 convec.f90 coriolis.f90 derivatives.f90 divstress_uv.f90
     divstress_w.f90 emul_complex.f90 fft.f90 finalize.f90
     forcing.f90 fringe.f90 functions.f90 grid.f90 inflow.f90
     initial.f90 initialize.f90 input_util.f90 interpolag_Sdep.f90
     interpolag_Ssim.f90 io.f90 lagrange_Sdep.f90 lagrange_Ssim.f90 main.f90
     messages.f90 param.f90 param_output.f90 pid.f90
     press_stag_array.f90 rmsdiv.f90 scaledep_dynamic.f90 shifted_inflow.f90
     sgs_param.f90 sgs_stag_util.f90 sim_param.f90 sponge.f90 stat_defs.f90
     std_dynamic.f90 string_util.f90 test_filtermodule.f90 tridag_array.f90
     time_average.f90 types.f90 wallstress.f90 init_random_seed.f90
     iwmles.f90 mts_wm.f90 qeqwm.f90 neqwm.f90 wm_param.f90 eqwm_smooth.f90)
#set (Sources cfl_util.f90 clocks.f90 convec.f90 derivatives.f90 divstress.f90
#     emul_complex.f90 fft.f90 finalize.f90 forcing.f90 fringe_util.f90
#     functions.f90 grid.f90 initial.f90 initialize.f90 input_util.f90
#     io.f90 main.f90 messages.f90 param.f90 param_output.f90
#     press_solve.f90 rmsdiv.f90 sgs_param.f90 sgs_stag_util.f90 
#     sim_param.f90 stat_defs.f90 string_util.f90 test_filtermodule.f90
#     types.f90 wallstress.f90 init_random_seed.f90 iwmles.f90
#     streamwise_bc.f90 bl_calcs.f90)

if (USE_MPI)
    add_definitions(-DPPMPI)
    set (CMAKE_Fortran_COMPILER mpif90)
    set(Sources ${Sources} mpi_transpose_mod.f90 mpi_defs.f90)
    set(exec_name "${exec_name}-mpi")
endif (USE_MPI)

if (USE_CPS)
    add_definitions(-DPPCPS)
    set(Sources ${Sources} concurrent_precursor.f90)
    set(exec_name "${exec_name}-cps")
endif (USE_CPS)

if (USE_HIT)
    add_definitions(-DPPHIT)
    set(Sources ${Sources} hit_inflow.f90)
    set(exec_name "${exec_name}-HIT")
endif (USE_HIT)

if (USE_STREAKS)
    add_definitions(-DPPSTREAKS)
    set(exec_name "${exec_name}-streaks")
endif (USE_STREAKS)

if (USE_LVLSET)
    add_definitions(-DPPLVLSET)
    set(Sources ${Sources} level_set_base.f90 level_set.f90 linear_simple.f90
        trees_pre_ls.f90 trees_base_ls.f90 trees_setup_ls.f90 trees_io_ls.f90
        trees_global_fmask_ls.f90)
    set(exec_name "${exec_name}-ls")
endif (USE_LVLSET)

if (USE_DYN_TN)
    add_definitions(-DPPDYN_TN)
    set(exec_name "${exec_name}-dyntn")
endif (USE_DYN_TN)

if (USE_SAFETYMODE)
    add_definitions(-DPPSAFETYMODE)
else (USE_SAFETYMODE)
    set(exec_name "${exec_name}-safety_off")
endif (USE_SAFETYMODE)

if (USE_CGNS)
    add_definitions(-DPPCGNS)
    link_libraries(cgns hdf5)
    set(exec_name "${exec_name}-cgns")
endif (USE_CGNS)

if (USE_MAPPING)
    add_definitions(-DPPMAPPING)
    set(Sources ${Sources} load_jacobian.f90)
endif (USE_MAPPING)

if (USE_LVLSET_STRETCH)
    add_definitions(-DPPLVLSET_STRETCH)
    set(Sources ${Sources} load_topography.f90 level_set2.f90)
    set(exec_name "${exec_name}-ls2")
endif (USE_LVLSET_STRETCH)

if (USE_TLWM_LVLSET)
    add_definitions(-DPPTLWM_LVLSET)
    set(exec_name "${exec_name}-ls2")
endif (USE_TLWM_LVLSET)

if (OUTPUT_SGS)
    add_definitions(-DPPOUTPUT_SGS)
endif (OUTPUT_SGS)

if (OUTPUT_BUDGET)
    add_definitions(-DPPOUTPUT_BUDGET)
endif (OUTPUT_BUDGET)

if (OUTPUT_TURBSPEC)
    add_definitions(-DPPOUTPUT_TURBSPEC)
endif (OUTPUT_TURBSPEC)

if (OUTPUT_SPECBUDG)
    add_definitions(-DPPOUTPUT_SPECBUDG)
    add_definitions(-DPPOUTPUT_TURBSPEC)
    add_definitions(-DPPOUTPUT_BUDGET)
endif (OUTPUT_SPECBUDG)

if (OUTPUT_VELVORT)
    add_definitions(-DPPOUTPUT_VELVORT)
endif (OUTPUT_VELVORT)

if (OUTPUT_WMLES)
    add_definitions(-DPPOUTPUT_WMLES)
endif (OUTPUT_WMLES)

if (OUTPUT_SSP)
    add_definitions(-DPPOUTPUT_SSP)
endif (OUTPUT_SSP)

if (OUTPUT_CLOCK)
    add_definitions(-DPPOUTPUT_CLOCK)
endif (OUTPUT_CLOCK)

if (USE_BL)
    add_definitions(-DPPBL)
    set(Sources ${Sources} rescale_recycle.f90 rescale_recycle_fluc.f90	functions_bl.f90 )
    set(exec_name "${exec_name}-bl")
endif (USE_BL)

if (WRITE_ENDIAN STREQUAL "LITTLE")
    add_definitions(-DPPWRITE_LITTLE_ENDIAN)
endif (WRITE_ENDIAN STREQUAL "LITTLE")

if (WRITE_ENDIAN STREQUAL "BIG")
    add_definitions(-DPPWRITE_BIG_ENDIAN)
endif (WRITE_ENDIAN STREQUAL "BIG")

if (READ_ENDIAN STREQUAL "LITTLE")
    add_definitions(-DPPREAD_LITTLE_ENDIAN)
endif (READ_ENDIAN STREQUAL "LITTLE")

if (READ_ENDIAN STREQUAL "BIG")
    add_definitions(-DPPREAD_BIG_ENDIAN)
endif (READ_ENDIAN STREQUAL "BIG")

# Intel
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O2 -funroll-loops -no-wrap-margin -assume byterecl")
    #set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O0 -g -traceback -no-wrap-margin -assume byterecl")
# GNU
else (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    #set (CMAKE_Fortran_FLAGS "-cpp -O2 -ffast-math -funroll-loops -fall-intrinsics -Wall -fcheck=all")
    set (CMAKE_Fortran_FLAGS "-cpp -O0 -g -fbacktrace -fall-intrinsics -Wall -fcheck=all")
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

# main lesgo executable
add_executable(${exec_name} ${Sources})
