################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CRF_FeatureMap.cpp \
../src/CRF_FeatureStream.cpp \
../src/CRF_FeatureStreamManager.cpp \
../src/CRF_GradBuilder.cpp \
../src/CRF_InFtrStream_RandPresent.cpp \
../src/CRF_InLabStream_RandPresent.cpp \
../src/CRF_LabelPath.cpp \
../src/CRF_LatticeBuilder.cpp \
../src/CRF_LocalPosteriorBuilder.cpp \
../src/CRF_LogMath.cpp \
../src/CRF_Model.cpp \
../src/CRF_NewGradBuilder.cpp \
../src/CRF_NewGradBuilderLog.cpp \
../src/CRF_NewGradBuilderSparseLog.cpp \
../src/CRF_NewLocalPosteriorBuilder.cpp \
../src/CRF_NewViterbi.cpp \
../src/CRF_SGTrainer.cpp \
../src/CRF_Seq.cpp \
../src/CRF_SparseFeatureMap.cpp \
../src/CRF_StateNode.cpp \
../src/CRF_StateVector.cpp \
../src/CRF_StdFeatureMap.cpp \
../src/CRF_StdGradBuilder.cpp \
../src/CRF_StdGradBuilderLog.cpp \
../src/CRF_StdNStateNode.cpp \
../src/CRF_StdNStateNodeLog.cpp \
../src/CRF_StdNStateVector.cpp \
../src/CRF_StdNStateVectorLog.cpp \
../src/CRF_StdSparseFeatureMap.cpp \
../src/CRF_StdStateNode.cpp \
../src/CRF_StdStateNodeLog.cpp \
../src/CRF_StdStateVector.cpp \
../src/CRF_StdStateVectorLog.cpp \
../src/CRF_StdTransFeatureMap.cpp \
../src/CRF_Trainer.cpp 

OBJS += \
./src/CRF_FeatureMap.o \
./src/CRF_FeatureStream.o \
./src/CRF_FeatureStreamManager.o \
./src/CRF_GradBuilder.o \
./src/CRF_InFtrStream_RandPresent.o \
./src/CRF_InLabStream_RandPresent.o \
./src/CRF_LabelPath.o \
./src/CRF_LatticeBuilder.o \
./src/CRF_LocalPosteriorBuilder.o \
./src/CRF_LogMath.o \
./src/CRF_Model.o \
./src/CRF_NewGradBuilder.o \
./src/CRF_NewGradBuilderLog.o \
./src/CRF_NewGradBuilderSparseLog.o \
./src/CRF_NewLocalPosteriorBuilder.o \
./src/CRF_NewViterbi.o \
./src/CRF_SGTrainer.o \
./src/CRF_Seq.o \
./src/CRF_SparseFeatureMap.o \
./src/CRF_StateNode.o \
./src/CRF_StateVector.o \
./src/CRF_StdFeatureMap.o \
./src/CRF_StdGradBuilder.o \
./src/CRF_StdGradBuilderLog.o \
./src/CRF_StdNStateNode.o \
./src/CRF_StdNStateNodeLog.o \
./src/CRF_StdNStateVector.o \
./src/CRF_StdNStateVectorLog.o \
./src/CRF_StdSparseFeatureMap.o \
./src/CRF_StdStateNode.o \
./src/CRF_StdStateNodeLog.o \
./src/CRF_StdStateVector.o \
./src/CRF_StdStateVectorLog.o \
./src/CRF_StdTransFeatureMap.o \
./src/CRF_Trainer.o 

CPP_DEPS += \
./src/CRF_FeatureMap.d \
./src/CRF_FeatureStream.d \
./src/CRF_FeatureStreamManager.d \
./src/CRF_GradBuilder.d \
./src/CRF_InFtrStream_RandPresent.d \
./src/CRF_InLabStream_RandPresent.d \
./src/CRF_LabelPath.d \
./src/CRF_LatticeBuilder.d \
./src/CRF_LocalPosteriorBuilder.d \
./src/CRF_LogMath.d \
./src/CRF_Model.d \
./src/CRF_NewGradBuilder.d \
./src/CRF_NewGradBuilderLog.d \
./src/CRF_NewGradBuilderSparseLog.d \
./src/CRF_NewLocalPosteriorBuilder.d \
./src/CRF_NewViterbi.d \
./src/CRF_SGTrainer.d \
./src/CRF_Seq.d \
./src/CRF_SparseFeatureMap.d \
./src/CRF_StateNode.d \
./src/CRF_StateVector.d \
./src/CRF_StdFeatureMap.d \
./src/CRF_StdGradBuilder.d \
./src/CRF_StdGradBuilderLog.d \
./src/CRF_StdNStateNode.d \
./src/CRF_StdNStateNodeLog.d \
./src/CRF_StdNStateVector.d \
./src/CRF_StdNStateVectorLog.d \
./src/CRF_StdSparseFeatureMap.d \
./src/CRF_StdStateNode.d \
./src/CRF_StdStateNodeLog.d \
./src/CRF_StdStateVector.d \
./src/CRF_StdStateVectorLog.d \
./src/CRF_StdTransFeatureMap.d \
./src/CRF_Trainer.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


