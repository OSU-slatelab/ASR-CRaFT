################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include src/subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: CRFDecode

# Tool invocations
CRFDecode: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L/u/drspeech/src/quicknet-v3_20/H-i586-linux -L"/u/morrijer/workspaces/CRF_cpp_new/CRF/Debug" -L/u/drspeech/src/quicknet-v3_20 -L/u/drspeech/src/ATLAS-3.7.8/lib/Linux_P4SSE2_4 -o"CRFDecode" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(EXECUTABLES)$(CC_DEPS)$(C_DEPS)$(CPP_DEPS)$(CXX_DEPS)$(C_UPPER_DEPS) CRFDecode
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
