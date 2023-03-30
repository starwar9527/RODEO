################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Rosenbrock/UnconstrainedUsingGradients/Rosenbrock.cpp 

CPP_DEPS += \
./Tests/Optimization/Rosenbrock/UnconstrainedUsingGradients/Rosenbrock.d 

OBJS += \
./Tests/Optimization/Rosenbrock/UnconstrainedUsingGradients/Rosenbrock.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Rosenbrock/UnconstrainedUsingGradients/%.o: ../Tests/Optimization/Rosenbrock/UnconstrainedUsingGradients/%.cpp Tests/Optimization/Rosenbrock/UnconstrainedUsingGradients/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-UnconstrainedUsingGradients

clean-Tests-2f-Optimization-2f-Rosenbrock-2f-UnconstrainedUsingGradients:
	-$(RM) ./Tests/Optimization/Rosenbrock/UnconstrainedUsingGradients/Rosenbrock.d ./Tests/Optimization/Rosenbrock/UnconstrainedUsingGradients/Rosenbrock.o

.PHONY: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-UnconstrainedUsingGradients

