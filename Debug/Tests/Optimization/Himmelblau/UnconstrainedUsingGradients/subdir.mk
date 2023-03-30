################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Himmelblau/UnconstrainedUsingGradients/himmelblau.cpp 

CPP_DEPS += \
./Tests/Optimization/Himmelblau/UnconstrainedUsingGradients/himmelblau.d 

OBJS += \
./Tests/Optimization/Himmelblau/UnconstrainedUsingGradients/himmelblau.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Himmelblau/UnconstrainedUsingGradients/%.o: ../Tests/Optimization/Himmelblau/UnconstrainedUsingGradients/%.cpp Tests/Optimization/Himmelblau/UnconstrainedUsingGradients/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Himmelblau-2f-UnconstrainedUsingGradients

clean-Tests-2f-Optimization-2f-Himmelblau-2f-UnconstrainedUsingGradients:
	-$(RM) ./Tests/Optimization/Himmelblau/UnconstrainedUsingGradients/himmelblau.d ./Tests/Optimization/Himmelblau/UnconstrainedUsingGradients/himmelblau.o

.PHONY: clean-Tests-2f-Optimization-2f-Himmelblau-2f-UnconstrainedUsingGradients

