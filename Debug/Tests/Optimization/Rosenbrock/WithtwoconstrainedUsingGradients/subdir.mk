################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock.cpp \
../Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock1.cpp 

CPP_DEPS += \
./Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock.d \
./Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock1.d 

OBJS += \
./Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock.o \
./Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock1.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/%.o: ../Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/%.cpp Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithtwoconstrainedUsingGradients

clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithtwoconstrainedUsingGradients:
	-$(RM) ./Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock.d ./Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock.o ./Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock1.d ./Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/Rosenbrock1.o

.PHONY: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithtwoconstrainedUsingGradients

