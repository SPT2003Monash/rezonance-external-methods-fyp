
#------------------------------------------------------------------------------
# Project 'GFMBESS_PSCAD502_20250513' make using the 'Intel_ Fortran Compiler Classic 2021.12.0 (64-bit)' compiler.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# All project
#------------------------------------------------------------------------------

all: targets
	@echo !--Make: succeeded.



#------------------------------------------------------------------------------
# Directories, Platform, and Version
#------------------------------------------------------------------------------

Arch        = windows
EmtdcDir    = C:\Program Files (x86)\PSCAD50\emtdc\if18
EmtdcInc    = $(EmtdcDir)\inc
EmtdcBin    = $(EmtdcDir)\$(Arch)
EmtdcMain   = $(EmtdcBin)\main.obj
EmtdcLib    = $(EmtdcBin)\emtdc.lib
SolverLib    = $(EmtdcBin)\Solver.lib


#------------------------------------------------------------------------------
# Fortran Compiler
#------------------------------------------------------------------------------

FC_Name         = ifort.exe
FC_Suffix       = obj
FC_Args         = /nologo /c /free /real_size:64 /fpconstant /warn:declarations /iface:default /align:dcommons /fpe:0
FC_Debug        =  /O2
FC_Preprocess   = 
FC_Preproswitch = 
FC_Warn         = 
FC_Checks       = 
FC_Includes     = /include:"$(EmtdcInc)" /include:"$(EmtdcDir)" /include:"$(EmtdcBin)"
FC_Compile      = $(FC_Name) $(FC_Args) $(FC_Includes) $(FC_Debug) $(FC_Warn) $(FC_Checks)

#------------------------------------------------------------------------------
# C Compiler
#------------------------------------------------------------------------------

CC_Name     = cl.exe
CC_Suffix   = obj
CC_Args     = /nologo /MT /W3 /EHsc /c
CC_Debug    =  /O2
CC_Includes = 
CC_Compile  = $(CC_Name) $(CC_Args) $(CC_Includes) $(CC_Debug)

#------------------------------------------------------------------------------
# Linker
#------------------------------------------------------------------------------

Link_Name   = link.exe
Link_Debug  = 
Link_Args   = /out:$@ /nologo /nodefaultlib:libc.lib /nodefaultlib:libcmtd.lib /subsystem:console
Link        = $(Link_Name) $(Link_Args) $(Link_Debug)

#------------------------------------------------------------------------------
# Build rules for generated files
#------------------------------------------------------------------------------


.f.$(FC_Suffix):
	@echo !--Compile: $<
	$(FC_Compile) $<



.c.$(CC_Suffix):
	@echo !--Compile: $<
	$(CC_Compile) $<



#------------------------------------------------------------------------------
# Build rules for file references
#------------------------------------------------------------------------------


rezonance_if18_1.lib: 
	@echo !--Copy: "C:\Users\erict\Downloads\drive-download-20260413T095038Z-3-001 (1)\pscad_module\lib\if18\rezonance_if18.lib"
	copy "C:\Users\erict\Downloads\drive-download-20260413T095038Z-3-001 (1)\pscad_module\lib\if18\rezonance_if18.lib" "rezonance_if18_1.lib"

#------------------------------------------------------------------------------
# Dependencies
#------------------------------------------------------------------------------


FC_Objects = \
 Station.$(FC_Suffix) \
 Main.$(FC_Suffix) \
 BESS.$(FC_Suffix) \
 PPC.$(FC_Suffix) \
 Startup.$(FC_Suffix) \
 AVR.$(FC_Suffix) \
 PI_AntiWindUp.$(FC_Suffix) \
 Battery.$(FC_Suffix) \
 UC.$(FC_Suffix) \
 Parameters_Bases.$(FC_Suffix) \
 Measurement_Filter.$(FC_Suffix) \
 GSC.$(FC_Suffix) \
 Primary_Controller.$(FC_Suffix) \
 APC_Droop.$(FC_Suffix) \
 RPC_Droop.$(FC_Suffix) \
 RPC_PI.$(FC_Suffix) \
 APC_VSG.$(FC_Suffix) \
 APC_GCVSG.$(FC_Suffix) \
 RPC_RotorFluxModel.$(FC_Suffix) \
 Calculations.$(FC_Suffix) \
 ABC_to_DQ0.$(FC_Suffix) \
 GSC_Voltage_Controller.$(FC_Suffix) \
 Current_Limiter.$(FC_Suffix) \
 GSC_Current_Controller.$(FC_Suffix) \
 Uref_Generator.$(FC_Suffix) \
 DQ0_to_ABC.$(FC_Suffix) \
 Virtual_Impedance.$(FC_Suffix) \
 Parameters_Filters.$(FC_Suffix) \
 ConvBridge_IGBT.$(FC_Suffix) \
 PWM.$(FC_Suffix) \
 Harmonic_Filter.$(FC_Suffix)

FC_ObjectsLong = \
 "Station.$(FC_Suffix)" \
 "Main.$(FC_Suffix)" \
 "BESS.$(FC_Suffix)" \
 "PPC.$(FC_Suffix)" \
 "Startup.$(FC_Suffix)" \
 "AVR.$(FC_Suffix)" \
 "PI_AntiWindUp.$(FC_Suffix)" \
 "Battery.$(FC_Suffix)" \
 "UC.$(FC_Suffix)" \
 "Parameters_Bases.$(FC_Suffix)" \
 "Measurement_Filter.$(FC_Suffix)" \
 "GSC.$(FC_Suffix)" \
 "Primary_Controller.$(FC_Suffix)" \
 "APC_Droop.$(FC_Suffix)" \
 "RPC_Droop.$(FC_Suffix)" \
 "RPC_PI.$(FC_Suffix)" \
 "APC_VSG.$(FC_Suffix)" \
 "APC_GCVSG.$(FC_Suffix)" \
 "RPC_RotorFluxModel.$(FC_Suffix)" \
 "Calculations.$(FC_Suffix)" \
 "ABC_to_DQ0.$(FC_Suffix)" \
 "GSC_Voltage_Controller.$(FC_Suffix)" \
 "Current_Limiter.$(FC_Suffix)" \
 "GSC_Current_Controller.$(FC_Suffix)" \
 "Uref_Generator.$(FC_Suffix)" \
 "DQ0_to_ABC.$(FC_Suffix)" \
 "Virtual_Impedance.$(FC_Suffix)" \
 "Parameters_Filters.$(FC_Suffix)" \
 "ConvBridge_IGBT.$(FC_Suffix)" \
 "PWM.$(FC_Suffix)" \
 "Harmonic_Filter.$(FC_Suffix)"

CC_Objects =

CC_ObjectsLong =

LK_Objects = \
  rezonance_if18_1.lib

LK_ObjectsLong = \
  "rezonance_if18_1.lib"

SysLibs  = ws2_32.lib

Binary   = GFMBESS_PSCAD502_20250513.exe

$(Binary): $(FC_Objects) $(CC_Objects) $(LK_Objects) 
	@echo !--Link: $@
	$(Link) "$(EmtdcMain)" $(FC_ObjectsLong) $(CC_ObjectsLong) $(LK_ObjectsLong) "$(EmtdcLib)" "$(SolverLib)" $(SysLibs)

targets: $(Binary)


clean:
	-del EMTDC_V*
	-del *.obj
	-del *.o
	-del *.exe
	@echo !--Make clean: succeeded.



