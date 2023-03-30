################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock.cpp \
../Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock1.cpp 

CPP_DEPS += \
./Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock.d \
./Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock1.d 

OBJS += \
./Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock.o \
./Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock1.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/%.o: ../Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/%.cpp Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithoneconstrainedUsingGradients

clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithoneconstrainedUsingGradients:
	-$(RM) ./Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock.d ./Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock.o ./Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock1.d ./Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/Rosenbrock1.o

.PHONY: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithoneconstrainedUsingGradients

