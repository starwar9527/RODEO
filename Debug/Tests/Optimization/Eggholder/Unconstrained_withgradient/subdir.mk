################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Optimization/Eggholder/Unconstrained_withgradient/eggholder.cpp 

CPP_DEPS += \
./Tests/Optimization/Eggholder/Unconstrained_withgradient/eggholder.d 

OBJS += \
./Tests/Optimization/Eggholder/Unconstrained_withgradient/eggholder.o 


# Each subdirectory must supply rules for building sources it contributes
Tests/Optimization/Eggholder/Unconstrained_withgradient/%.o: ../Tests/Optimization/Eggholder/Unconstrained_withgradient/%.cpp Tests/Optimization/Eggholder/Unconstrained_withgradient/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-Tests-2f-Optimization-2f-Eggholder-2f-Unconstrained_withgradient

clean-Tests-2f-Optimization-2f-Eggholder-2f-Unconstrained_withgradient:
	-$(RM) ./Tests/Optimization/Eggholder/Unconstrained_withgradient/eggholder.d ./Tests/Optimization/Eggholder/Unconstrained_withgradient/eggholder.o

.PHONY: clean-Tests-2f-Optimization-2f-Eggholder-2f-Unconstrained_withgradient

