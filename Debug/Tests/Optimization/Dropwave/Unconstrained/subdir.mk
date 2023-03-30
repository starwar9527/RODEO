################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Dropwave/Unconstrained/Dropwave.cpp 

CPP_DEPS += \
./Tests/Optimization/Dropwave/Unconstrained/Dropwave.d 

OBJS += \
./Tests/Optimization/Dropwave/Unconstrained/Dropwave.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Dropwave/Unconstrained/%.o: ../Tests/Optimization/Dropwave/Unconstrained/%.cpp Tests/Optimization/Dropwave/Unconstrained/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Dropwave-2f-Unconstrained

clean-Tests-2f-Optimization-2f-Dropwave-2f-Unconstrained:
	-$(RM) ./Tests/Optimization/Dropwave/Unconstrained/Dropwave.d ./Tests/Optimization/Dropwave/Unconstrained/Dropwave.o

.PHONY: clean-Tests-2f-Optimization-2f-Dropwave-2f-Unconstrained

