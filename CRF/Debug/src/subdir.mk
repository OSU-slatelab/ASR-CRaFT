################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CRF_Model.cpp 

OBJS += \
./src/CRF_Model.o 

CPP_DEPS += \
./src/CRF_Model.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	/usr/bin/g++ -I/u/drspeech/src/quicknet-v3_20/H-x86_64 -I/u/drspeech/src/ATLAS-3.8.2/Linux_XE64SSE3/include -I/u/drspeech/src/ATLAS-3.8.2/include -I/u/drspeech/src/quicknet-v3_20 -I/u/drspeech/opt/OpenFst-1.1/include -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


