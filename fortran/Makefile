CC=gfortran -fopenmp -O3
FFLAGS=-c

OBJECTS_COMMON1 = misc.f90 photon_packet.f90 dust.f90 sources.f90

OBJECTS_COMMON2 = grid_general.f90 lucy_iteration.f90 bw_iteration.f90 thermal_mc.f90 main.f90

OBJECTS_COMMON2_OMP = grid_general.f90 lucy_iteration_omp.f90 bw_iteration_omp.f90 thermal_mc.f90 main.f90

OBJECTS_CAR = $(OBJECTS_COMMON1) grid_cartesian.f90 $(OBJECTS_COMMON2)

OBJECTS_CYL = $(OBJECTS_COMMON1) grid_cylindrical.f90 $(OBJECTS_COMMON2)

OBJECTS_SPH = $(OBJECTS_COMMON1) grid_spherical.f90 $(OBJECTS_COMMON2)

OBJECTS_CAR_OMP = $(OBJECTS_COMMON1) grid_cartesian.f90 $(OBJECTS_COMMON2_OMP)

OBJECTS_CYL_OMP = $(OBJECTS_COMMON1) grid_cylindrical.f90 $(OBJECTS_COMMON2_OMP)

OBJECTS_SPH_OMP = $(OBJECTS_COMMON1) grid_spherical.f90 $(OBJECTS_COMMON2_OMP)

all: car polarmc_car cyl polarmc_cyl sph polarmc_sph

# Single thread versions.

polarmc_car: $(OBJECTS_CAR)
	$(CC) $(OBJECTS_CAR) -o polarmc_car

car: $(OBJECTS_CAR)
	$(CC) $(FFLAGS) $(OBJECTS_CAR)

polarmc_cyl: $(OBJECTS_CYL)
	$(CC) $(OBJECTS_CYL) -o polarmc_cyl

cyl: $(OBJECTS_CYL)
	$(CC) $(FFLAGS) $(OBJECTS_CYL)

polarmc_sph: $(OBJECTS_SPH)
	$(CC) $(OBJECTS_SPH) -o polarmc_sph

sph: $(OBJECTS_SPH)
	$(CC) $(FFLAGS) $(OBJECTS_SPH)

# OpenMP versions.

polarmc_car_omp: $(OBJECTS_CAR_OMP)
	$(CC) $(OBJECTS_CAR_OMP) -o polarmc_car_omp

car_omp: $(OBJECTS_CAR_OMP)
	$(CC) $(FFLAGS) $(OBJECTS_CAR_OMP)

polarmc_cyl_omp: $(OBJECTS_CYL_OMP)
	$(CC) $(OBJECTS_CYL_OMP) -o polarmc_cyl_omp

cyl_omp: $(OBJECTS_CYL_OMP)
	$(CC) $(FFLAGS) $(OBJECTS_CYL_OMP)

polarmc_sph_omp: $(OBJECTS_SPH_OMP)
	$(CC) $(OBJECTS_SPH) -o polarmc_sph_omp

sph_omp: $(OBJECTS_SPH)
	$(CC) $(FFLAGS) $(OBJECTS_SPH_OMP)

clean:
	rm -rf *.o *.mod polarmc_*
