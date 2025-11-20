#========================================================================
#
#                   S P E C F E M 2 D  Version 7 . 0
#                   --------------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                              CNRS, France
#                       and Princeton University, USA
#                 (there are currently many more authors!)
#                           (c) October 2017
#
# This software is a computer program whose purpose is to solve
# the two-dimensional viscoelastic anisotropic or poroelastic wave equation
# using a spectral-element method (SEM).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# The full text of the license is available in file "LICENSE".
#
#========================================================================
#
# Makefile.  Generated from Makefile.in by configure.
#######################################

FC = gfortran
FCFLAGS = -g -O2
FC_DEFINE = -D
MPIFC = mpif90
MPILIBS = 

FLAGS_CHECK = -std=f2008 -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused -Werror -O3 -finline-functions

FCFLAGS_f90 = -J./obj -I./obj -I.  -I${SETUP}

FC_MODEXT = mod
FC_MODDIR = ./obj

FCCOMPILE_CHECK = ${FC} ${FCFLAGS} $(FLAGS_CHECK)

MPIFCCOMPILE_CHECK = ${MPIFC} ${FCFLAGS} $(FLAGS_CHECK)

CC = gcc
CFLAGS = -g -O2 $(CPPFLAGS) -I${SETUP}
CPPFLAGS =  $(COND_MPI_CPPFLAGS)

# all linker flags
MPILIBS +=  

#######################################
####
#### MPI
####
#######################################

## serial or parallel
MPI = yes
#MPI = no

FCLINK = $(MPIFCCOMPILE_CHECK)
#FCLINK = $(FCCOMPILE_CHECK)

COND_MPI_CPPFLAGS = $(FC_DEFINE)WITH_MPI
#COND_MPI_CPPFLAGS =

MPI_INCLUDES = 

# fortran compiler setting
F90 = $(MPIFCCOMPILE_CHECK) $(FC_DEFINE)USE_MPI
#F90 = $(FCCOMPILE_CHECK)

#######################################
####
#### SCOTCH
####
#######################################

SCOTCH=yes
#SCOTCH=no

USE_BUNDLED_SCOTCH = 1

SCOTCH_DIR = ./external_libs/scotch
SCOTCH_INCDIR = ./external_libs/scotch/include
SCOTCH_LIBDIR = ./external_libs/scotch/lib

SCOTCH_INC = -I${SCOTCH_INCDIR}
SCOTCH_LIBS = -L${SCOTCH_LIBDIR} -lscotch -lscotcherr

SCOTCH_FLAGS = -DUSE_SCOTCH $(SCOTCH_INC)

#######################################
####
#### CUDA
#### with configure: ./configure --with-cuda=cuda5 CUDA_FLAGS=.. CUDA_LIB=.. CUDA_INC=.. MPI_INC=.. ..
####
#######################################

#CUDA = yes
CUDA = no

#CUDA5 = yes
CUDA5 = no

#CUDA6 = yes
CUDA6 = no

#CUDA7 = yes
CUDA7 = no

#CUDA8 = yes
CUDA8 = no

#CUDA9 = yes
CUDA9 = no

# CUDA compilation with linking
#CUDA_PLUS = yes
CUDA_PLUS = no

# default cuda libraries
# runtime library -lcudart needed, others are optional -lcuda -lcublas

CUDA_FLAGS = 
CUDA_INC =  -I${SETUP}
CUDA_LINK =   -lstdc++

#NVCC = nvcc
NVCC = gcc

##
## GPU architecture
##
# CUDA architecture / code version
# Fermi: -gencode=arch=compute_10,code=sm_10 not supported
# Tesla (default): -gencode=arch=compute_20,code=sm_20
# Geforce GT 650m: -gencode=arch=compute_30,code=sm_30
# Kepler (cuda5) : -gencode=arch=compute_35,code=sm_35
# Kepler (cuda6.5,K80): -gencode=arch=compute_37,code=sm_37
# Maxwell (cuda6.5+/cuda7,Quadro K2200): -gencode=arch=compute_50,code=sm_50
# Pascal (cuda8,P100): -gencode=arch=compute_60,code=sm_60
# Volta (cuda9,V100): -gencode=arch=compute_70,code=sm_70
GENCODE_20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE_30 = -gencode=arch=compute_30,code=\"sm_30,compute_30\"
GENCODE_35 = -gencode=arch=compute_35,code=\"sm_35,compute_35\"
GENCODE_37 = -gencode=arch=compute_37,code=\"sm_37\"
GENCODE_50 = -gencode=arch=compute_50,code=\"sm_50\"
GENCODE_60 = -gencode=arch=compute_60,code=\"sm_60,compute_60\"
GENCODE_70 = -gencode=arch=compute_70,code=\"sm_70,compute_70\"

# cuda preprocessor flag
# CUDA version 9.0
##GENCODE = $(GENCODE_70) $(FC_DEFINE)GPU_DEVICE_Volta
# CUDA version 8.0
##GENCODE = $(GENCODE_60) $(FC_DEFINE)GPU_DEVICE_Pascal
# CUDA version 7.x
##GENCODE = $(GENCODE_50) $(FC_DEFINE)GPU_DEVICE_Maxwell
# CUDA version 6.5
##GENCODE = $(GENCODE_37) $(FC_DEFINE)GPU_DEVICE_K80
# CUDA version 5.x
##GENCODE = $(GENCODE_35) $(FC_DEFINE)GPU_DEVICE_K20
# CUDA version 4.x
#GENCODE = $(GENCODE_20)

# CUDA flags and linking
#NVCC_FLAGS_BASE = $(CUDA_FLAGS) $(CUDA_INC) $(MPI_INCLUDES) $(COND_MPI_CPPFLAGS)
##NVCC_FLAGS = $(NVCC_FLAGS_BASE) -dc -DCUDA $(GENCODE) -ftz=true
#NVCC_FLAGS = $(NVCC_FLAGS_BASE) -DCUDA -DUSE_OLDER_CUDA4_GPU $(GENCODE)

##NVCCLINK_BASE = $(NVCC) $(CUDA_INC) $(MPI_INCLUDES) $(COND_MPI_CPPFLAGS) -DCUDA
##NVCCLINK = $(NVCCLINK_BASE) -dlink $(GENCODE)
#NVCCLINK = $(NVCCLINK_BASE) -DUSE_OLDER_CUDA4_GPU $(GENCODE)

NVCC_FLAGS = $(MPI_INCLUDES) $(COND_MPI_CPPFLAGS)
NVCCLINK = $(NVCC) $(NVCC_FLAGS)


#######################################
####
#### directories
####
#######################################

## compilation directories
# B : build directory
B = .
# E : executables directory
E = $B/bin
# O : objects directory
O = $B/obj
# S_TOP : source file root directory
S_TOP = .
# L : libraries directory
L = $B/lib
# setup file directory
SETUP = $B/setup
# output file directory
OUTPUT = $B/OUTPUT_FILES


#######################################
####
#### targets
####
#######################################

# code subdirectories
SUBDIRS = \
	auxiliaries \
	cuda \
	meshfem2D \
	shared \
	specfem2D \
	tomography/postprocess_sensitivity_kernels \
	tomography \
	$(EMPTY_MACRO)

# default targets for the pure Fortran version
DEFAULT = \
	xmeshfem2D \
	xspecfem2D \
	xadj_seismogram \
	xcheck_quality_external_mesh \
	xconvolve_source_timefunction \
	$(EMPTY_MACRO)


default: $(DEFAULT)

all: default auxiliaries postprocess tomography


ifdef CLEAN
clean:
	@echo "cleaning by CLEAN"
	-rm -f $(foreach dir, $(CLEAN), $($(dir)_OBJECTS) $($(dir)_MODULES) $($(dir)_SHARED_OBJECTS) $($(dir)_TARGETS))
	-rm -f ${E}/*__genmod.*
	-rm -f ${O}/*__genmod.*
else
clean:
	@echo "cleaning all"
	-rm -f $(foreach dir, $(SUBDIRS), $($(dir)_OBJECTS) $($(dir)_MODULES) $($(dir)_TARGETS))
	-rm -f ${E}/*__genmod.*
	-rm -f ${O}/*__genmod.*
endif

realclean: clean
ifeq (${USE_BUNDLED_SCOTCH},1)
	@echo "cleaning bundled Scotch in directory: ${SCOTCH_DIR}/src"
	$(MAKE) -C ${SCOTCH_DIR}/src realclean
endif
	-rm -rf $E/* $O/*

# unit testing
# If the first argument is "test"...
ifeq (test,$(findstring test,firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  TEST_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # turn them into do-nothing targets
  $(eval $(TEST_ARGS):;@:)
endif

tests:
	@echo "testing in directory: ${S_TOP}/tests/"
	cd ${S_TOP}/tests; ./run_all_tests.sh $(TEST_ARGS)
	@echo ""


help:
	@echo "usage: make [executable]"
	@echo ""
	@echo "supported main executables:"
	@echo "    xspecfem2D"
	@echo "    xmeshfem2D"
	@echo ""
	@echo "additional executables:"
	@echo "- auxiliary executables: [make aux]"
	@echo "    xadj_seismogram"
	@echo "    xcheck_quality_external_mesh"
	@echo "    xconvolve_source_timefunction"
	@echo ""
	@echo "- sensitivity kernel postprocessing tools: [make postprocess]"
	@echo "    xcombine_sem"
	@echo "    xsmooth_sem"
	@echo ""
	@echo "- tomography tools: [make tomography]"
	@echo "    xsum_kernels"
	@echo ""
	@echo "for unit testing:"
	@echo "    tests"
	@echo ""


.PHONY: all default backup clean realclean help tests

#######################################

${SETUP}/version.fh: 
	@echo "GEN $@"
	@echo "! This file is generated by Make. Do not edit this file!" > $@
	@echo "character(len=*), parameter :: git_package_version = \"$$(cd ${S_TOP} && git describe --tags)\"" >> $@
	@echo "character(len=*), parameter :: git_commit_version = \"$$(cd ${S_TOP} && git log -n 1 | head -1 | /bin/grep commit)\"" >> $@
	@echo "character(len=*), parameter :: git_date_version = \"From $$(cd ${S_TOP} && git log -n 1 | head -5 | /bin/grep Date:)\"" >> $@

#######################################

# Get dependencies and rules for building stuff
include $(patsubst %, ${S_TOP}/src/%/rules.mk, $(SUBDIRS))

#######################################

##
## Shortcuts
##

# Shortcut for: <prog>/<xprog> -> bin/<xprog>
define target_shortcut
$(patsubst $E/%, %, $(1)): $(1)
.PHONY: $(patsubst $E/%, %, $(1))
$(patsubst $E/x%, %, $(1)): $(1)
.PHONY: $(patsubst $E/x%, %, $(1))
endef

# Shortcut for: dir -> src/dir/<targets in here>
define shortcut
$(1): $($(1)_TARGETS)
.PHONY: $(1)
$$(foreach target, $$(filter $E/%,$$($(1)_TARGETS)), $$(eval $$(call target_shortcut,$$(target))))
endef

$(foreach dir, $(SUBDIRS), $(eval $(call shortcut,$(dir))))

# shortcut testing
test : tests

# Other old shortcuts
mesh: $E/xmeshfem2D
spec: $E/xspecfem2D
.PHONY: mesh spec
