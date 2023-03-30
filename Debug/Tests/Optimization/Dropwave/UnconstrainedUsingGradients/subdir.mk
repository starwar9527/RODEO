################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Dropwave/UnconstrainedUsingGradients/Dropwave.cpp 

CPP_DEPS += \
./Tests/Optimization/Dropwave/UnconstrainedUsingGradients/Dropwave.d 

OBJS += \
./Tests/Optimization/Dropwave/UnconstrainedUsingGradients/Dropwave.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Dropwave/UnconstrainedUsingGradients/%.o: ../Tests/Optimization/Dropwave/UnconstrainedUsingGradients/%.cpp Tests/Optimization/Dropwave/UnconstrainedUsingGradients/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Dropwave-2f-UnconstrainedUsingGradients

clean-Tests-2f-Optimization-2f-Dropwave-2f-UnconstrainedUsingGradients:
	-$(RM) ./Tests/Optimization/Dropwave/UnconstrainedUsingGradients/Dropwave.d ./Tests/Optimization/Dropwave/UnconstrainedUsingGradients/Dropwave.o

.PHONY: clean-Tests-2f-Optimization-2f-Dropwave-2f-UnconstrainedUsingGradients

