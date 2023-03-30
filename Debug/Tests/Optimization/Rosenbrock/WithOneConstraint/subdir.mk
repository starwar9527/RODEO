################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Rosenbrock/WithOneConstraint/Rosenbrock.cpp 

CPP_DEPS += \
./Tests/Optimization/Rosenbrock/WithOneConstraint/Rosenbrock.d 

OBJS += \
./Tests/Optimization/Rosenbrock/WithOneConstraint/Rosenbrock.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Rosenbrock/WithOneConstraint/%.o: ../Tests/Optimization/Rosenbrock/WithOneConstraint/%.cpp Tests/Optimization/Rosenbrock/WithOneConstraint/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithOneConstraint

clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithOneConstraint:
	-$(RM) ./Tests/Optimization/Rosenbrock/WithOneConstraint/Rosenbrock.d ./Tests/Optimization/Rosenbrock/WithOneConstraint/Rosenbrock.o

.PHONY: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-WithOneConstraint

