proj = yes
FC = mpif90 
LINKER = mpic++ 
#LINKER = nvcc -arch=sm_70
CXX = $(PREP) nvcc
ckernel=yes 
openmp=no
raja_cuda=yes

LAPACK_DIR = /gpfs/alpine/csc340/proj-shared/larsen/sw4_uberenv/uberenv_libs/spack/opt/spack/linux-rhel7-ppc64le/gcc-5.4.0/netlib-lapack-3.8.0-2zkztxddylrbg4yua4bgm7afwcsgyh34
BLAS_DIR = /gpfs/alpine/csc340/proj-shared/larsen/sw4_uberenv/uberenv_libs/spack/opt/spack/linux-rhel7-ppc64le/gcc-5.4.0/netlib-lapack-3.8.0-2zkztxddylrbg4yua4bgm7afwcsgyh34
RAJA_LOCATION = /gpfs/alpine/csc340/proj-shared/larsen/sw4_uberenv/uberenv_libs/spack/opt/spack/linux-rhel7-ppc64le/gcc-5.4.0/raja-0.7.0-b37dcelfsjdjsesyyuwrb3l73xgptsze/
PROJ_HOME = /gpfs/alpine/csc340/proj-shared/larsen/sw4_uberenv/uberenv_libs/spack/opt/spack/linux-rhel7-ppc64le/gcc-5.4.0/proj-5.0.1-bdtgbpiylu6n3gmad54rpeam6lpb2bia/
UMPIRE_HOME = /gpfs/alpine/csc340/proj-shared/larsen/sw4_uberenv/uberenv_libs/spack/opt/spack/linux-rhel7-ppc64le/gcc-5.4.0/umpire-0.3.3-cqgdzeiuw2cl66vdhggzy5pltflsxuy2
CUDA_HOME = $(OLCF_CUDA_ROOT)
EXTRA_FORT_FLAGS = 
MORE_FLAGS = -DENABLE_MPI_TIMING_BARRIER=1 -DSW4_USE_UMPIRE=1 -DSW4_RAJA_VERSION=7
#LINKFLAGS =  -arch=sm_70
EXTRA_CXX_FLAGS =  -O3  -use_fast_math -DDISABLE_PREFETCH  -ccbin mpic++ -std=c++11 -Xcompiler "-fopenmp" --expt-extended-lambda -restrict -arch=sm_70 -I $(CUDA_HOME)/include -I$(RAJA_LOCATION)/include  -I$(CALIPER_LOCATION)/include --x cu -DUSE_NVTX -DRAJA_USE_CUDA -DSW4_CROUTINES -DRAJA_USE_RESTRICT_PTR -DCUDA_CODE -DENABLE_CUDA -dc -I$(PROJ_HOME)/include -I$(UMPIRE_HOME)/include $(MORE_FLAGS)

EXTRA_CXX_FLAGS2 =  -O3  -use_fast_math -DDISABLE_PREFETCH -ccbin mpic++ -std=c++11 --expt-extended-lambda -Xcompiler "-fopenmp" -restrict -arch=sm_70 -I $(CUDA_HOME)/include -I$(RAJA_LOCATION)/include  --x cu  -DUSE_NVTX -DRAJA_USE_CUDA -DSW4_CROUTINES -DRAJA_USE_RESTRICT_PTR -DCUDA_CODE -DENABLE_CUDA   -I$(UMPIRE_HOME)/include $(MORE_FLAGS)

EXTRA_LINK_FLAGS = -O3  -dc -L $(OLCF_XLF_ROOT)/lib/ -L $(LAPACK_DIR)/lib64 -llapack -lblas  -L $(CUDA_HOME)/lib64 -lcudart -lnvToolsExt -lcuda -L $(RAJA_LOCATION)/lib -lRAJA  -lcudadevrt -L $(PROJ_HOME)/lib -lproj -L /usr/lib64/nvidia -lnvidia-ml -L$(UMPIRE_HOME)/lib -lcudadevrt -lumpire -lgfortran -lgomp -lcudadevrt -lcuda
