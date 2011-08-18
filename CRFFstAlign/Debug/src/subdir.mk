################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Main.cpp 

OBJS += \
./src/Main.o 

CPP_DEPS += \
./src/Main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	/u/drspeech/opt/gcc-4.4.0/x86_64/bin/g++ -I"/home3/hey/segmental-branch/CRF/src" -I/u/drspeech/opt/OpenFst-1.1/include -I/u/drspeech/src/quicknet-v3_20/H-x86_64 -I/u/drspeech/src/ATLAS-3.8.2/Linux_XE64SSE3/include -I/u/drspeech/src/ATLAS-3.8.2/include -I/u/drspeech/src/quicknet-v3_20 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


