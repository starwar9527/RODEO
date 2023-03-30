################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Rosenbrock/Unconstrained/Rosenbrock.cpp 

CPP_DEPS += \
./Tests/Optimization/Rosenbrock/Unconstrained/Rosenbrock.d 

OBJS += \
./Tests/Optimization/Rosenbrock/Unconstrained/Rosenbrock.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Rosenbrock/Unconstrained/%.o: ../Tests/Optimization/Rosenbrock/Unconstrained/%.cpp Tests/Optimization/Rosenbrock/Unconstrained/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-Unconstrained

clean-Tests-2f-Optimization-2f-Rosenbrock-2f-Unconstrained:
	-$(RM) ./Tests/Optimization/Rosenbrock/Unconstrained/Rosenbrock.d ./Tests/Optimization/Rosenbrock/Unconstrained/Rosenbrock.o

.PHONY: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-Unconstrained

