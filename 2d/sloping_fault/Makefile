
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
# Note Seismic version takes setfault.py into account
CLAWMAKE = $(CLAW)/seismic/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
SETFAULT_FILE = setfault.py					# File containing function to make fault
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot_surface.py   # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= -fopenmp

# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
SEISLIB = $(CLAW)/seismic/src/2d
include $(AMRLIB)/Makefile.amr_2d

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \
  $(SEISLIB)/fault_module.f90

SOURCES = \
  qinit.f \
  setprob.f90 \
  bc2amr.f \
  setaux.f90 \
  mapc2p.f90 \
  $(SEISLIB)/b4step2.f90 \
	$(SEISLIB)/tick.f \
  $(CLAW)/seismic/src/rp/rpn2_vc_elasticity_mapped_fault_abl.f90 \
  $(CLAW)/seismic/src/rp/rpt2_vc_elasticity_mapped_fault_abl.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

elasticity_vc_mapped_fault_abl.so: ../../src/rp/rpn2_vc_elasticity_mapped_fault_abl.f90
	f2py -c ../../src/rp/rpn2_vc_elasticity_mapped_fault_abl.f90 -m elasticity_vc_mapped_fault_abl

sloping_fault.so: setaux.f90 mapc2p.f90 ../../src/2d/fault_module_no_types.f90 ../../src/2d/b4step2.f90
	f2py -c setaux.f90 mapc2p.f90 ../../src/2d/fault_module_no_types.f90 ../../src/2d/b4step2.f90 ${CLAW}/amrclaw/src/2d/opendatafile.f -m sloping_fault


