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
	g++ -I"/u/morrijer/workspaces/CRF_cpp_new/CRFDecode/src" -I"/u/morrijer/workspaces/CRF_cpp_new/CRF/src" -I/u/drspeech/src/quicknet-v3_20/H-i586-linux -I/u/drspeech/src/quicknet-v3_20/ -I/u/drspeech/src/ATLAS/include -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

