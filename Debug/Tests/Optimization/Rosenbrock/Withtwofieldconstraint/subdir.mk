################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Rosenbrock/Withtwofieldconstraint/Rosenbrock.cpp 

CPP_DEPS += \
./Tests/Optimization/Rosenbrock/Withtwofieldconstraint/Rosenbrock.d 

OBJS += \
./Tests/Optimization/Rosenbrock/Withtwofieldconstraint/Rosenbrock.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Rosenbrock/Withtwofieldconstraint/%.o: ../Tests/Optimization/Rosenbrock/Withtwofieldconstraint/%.cpp Tests/Optimization/Rosenbrock/Withtwofieldconstraint/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-Withtwofieldconstraint

clean-Tests-2f-Optimization-2f-Rosenbrock-2f-Withtwofieldconstraint:
	-$(RM) ./Tests/Optimization/Rosenbrock/Withtwofieldconstraint/Rosenbrock.d ./Tests/Optimization/Rosenbrock/Withtwofieldconstraint/Rosenbrock.o

.PHONY: clean-Tests-2f-Optimization-2f-Rosenbrock-2f-Withtwofieldconstraint

