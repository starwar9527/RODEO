################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Rosenbrock/WithTwoConstraint/Rosenbrock.cpp 

CPP_DEPS += \
./Tests/Optimization/Rosenbrock/WithTwoConstraint/Rosenbrock.d 

OBJS += \
./Tests/Optimization/Rosenbrock/WithTwoConstraint/Rosenbrock.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Rosenbrock/WithTwoConstraint/%.o: ../Tests/Optimization/Rosenbrock/WithTwoConstraint/%.cpp Tests/Optimization/Rosenbrock/WithTwoConstraint/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithTwoConstraint

clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithTwoConstraint:
	-$(RM) ./Tests/Optimization/Rosenbrock/WithTwoConstraint/Rosenbrock.d ./Tests/Optimization/Rosenbrock/WithTwoConstraint/Rosenbrock.o

.PHONY: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithTwoConstraint

