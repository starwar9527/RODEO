################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/subdir.mk
-include build/src/subdir.mk
-include UnitTests/src/subdir.mk
-include UnitTests/subdir.mk
-include Tests/Optimization/Rosenbrock/Withtwofieldconstraint/subdir.mk
-include Tests/Optimization/Rosenbrock/WithtwoconstrainedUsingGradients/subdir.mk
-include Tests/Optimization/Rosenbrock/WithoneconstrainedUsingGradients/subdir.mk
-include Tests/Optimization/Rosenbrock/Withmixedconstraint/subdir.mk
-include Tests/Optimization/Rosenbrock/Withfieldconstraint/subdir.mk
-include Tests/Optimization/Rosenbrock/WithTwoConstraint/subdir.mk
-include Tests/Optimization/Rosenbrock/WithOneConstraint/subdir.mk
-include Tests/Optimization/Rosenbrock/UnconstrainedUsingGradients/subdir.mk
-include Tests/Optimization/Rosenbrock/Unconstrained/subdir.mk
-include Tests/Optimization/Himmelblau/WithoneconstrainedUsingGradients/subdir.mk
-include Tests/Optimization/Himmelblau/WithOneConstraint/subdir.mk
-include Tests/Optimization/Himmelblau/UnconstrainedUsingGradients/subdir.mk
-include Tests/Optimization/Himmelblau/Unconstrained/subdir.mk
-include Tests/Optimization/Eggholder/Unconstrained_withgradient/subdir.mk
-include Tests/Optimization/Eggholder/Unconstrained/subdir.mk
-include Tests/Optimization/Dropwave/UnconstrainedUsingGradients/subdir.mk
-include Tests/Optimization/Dropwave/Unconstrained/subdir.mk
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

OPTIONAL_TOOL_DEPS := \
$(wildcard ../makefile.defs) \
$(wildcard ../makefile.init) \
$(wildcard ../makefile.targets) \


BUILD_ARTIFACT_NAME := Rodeo
BUILD_ARTIFACT_EXTENSION :=
BUILD_ARTIFACT_PREFIX :=
BUILD_ARTIFACT := $(BUILD_ARTIFACT_PREFIX)$(BUILD_ARTIFACT_NAME)$(if $(BUILD_ARTIFACT_EXTENSION),.$(BUILD_ARTIFACT_EXTENSION),)

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: main-build

# Main-build Target
main-build: Rodeo

# Tool invocations
Rodeo: $(OBJS) $(USER_OBJS) makefile $(OPTIONAL_TOOL_DEPS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++  -o "Rodeo" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) Rodeo
	-@echo ' '

.PHONY: all clean dependents main-build

-include ../makefile.targets
