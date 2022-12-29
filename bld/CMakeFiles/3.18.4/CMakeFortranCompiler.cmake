set(CMAKE_Fortran_COMPILER "/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/bin/gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "9.3.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/lib/gcc/x86_64-pc-linux-gnu/9.3.0/finclude;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/ipp/include;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/mkl/include;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/pstl/include;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/pstl/stdlib;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/tbb/include;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/daal/include;/cm/shared/apps/slurm/current/include;/data/apps/linux-centos8-cascadelake/intel-19.1.2.254/openmpi-3.1.6-d5iqzgx77lspks5kuha5bexs6ikhjrl7/include;/cm/shared/apps/Intel/2020/include;/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/lib/gcc/x86_64-pc-linux-gnu/9.3.0/include;/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/lib/gcc/x86_64-pc-linux-gnu/9.3.0/include-fixed;/usr/local/include;/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/lib/gcc/x86_64-pc-linux-gnu/9.3.0;/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/lib/gcc;/cm/shared/apps/slurm/current/lib64;/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/lib64;/lib64;/usr/lib64;/data/apps/linux-centos8-cascadelake/intel-19.1.2.254/openmpi-3.1.6-d5iqzgx77lspks5kuha5bexs6ikhjrl7/lib;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/mpi/intel64/libfabric/lib;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/ipp/lib/intel64;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/compiler/lib/intel64_lin;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/mkl/lib/intel64_lin;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/tbb/lib/intel64/gcc4.8;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/daal/lib/intel64_lin;/cm/shared/apps/Intel/2020/compilers_and_libraries_2020.2.254/linux/tbb/lib/intel64_lin/gcc4.8;/cm/shared/apps/Intel/2020/lib;/cm/shared/apps/slurm/current/lib64/slurm;/data/apps/linux-centos8-cascadelake/gcc-9.2.0/gcc-9.3.0-bnvby67rgbqevwsd264rgz44xucnkhpm/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
