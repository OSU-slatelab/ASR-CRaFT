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
src/Main.o: ../src/Main.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DUSE_SSE -I"/u/fosler/research/workspace-cpp/ASR-CRaFT/trunk/CRF/src" -I/u/drspeech/opt/OpenFst-beta-20080317/ -I/u/drspeech/src/quicknet-v3_20/H-i586-linux -I/u/drspeech/src/ATLAS/include -I/u/drspeech/src/quicknet-v3_20 -O3 -g3 -Wall -c -fmessage-length=0 -msse2 -ffast-math -MMD -MP -MF"$(@:%.o=%.d)" -MT"src/Main.d" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


