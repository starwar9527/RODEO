################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../UnitTests/Constraint2.cpp \
../UnitTests/himmelblau.cpp \
../UnitTests/himmelblauAdjWithOneConstraint.cpp \
../UnitTests/himmelblauAdjoint.cpp \
../UnitTests/himmelblauDoETest1.cpp \
../UnitTests/himmelblauDoETest2.cpp \
../UnitTests/himmelblauDoETest3.cpp \
../UnitTests/himmelblauDoETest4.cpp \
../UnitTests/himmelblauDoETest5.cpp \
../UnitTests/himmelblauDoETest6.cpp \
../UnitTests/himmelblauDoETestConstraint1.cpp \
../UnitTests/himmelblauWithOneConstraint.cpp \
../UnitTests/linearTestFunction.cpp 

CPP_DEPS += \
./UnitTests/Constraint2.d \
./UnitTests/himmelblau.d \
./UnitTests/himmelblauAdjWithOneConstraint.d \
./UnitTests/himmelblauAdjoint.d \
./UnitTests/himmelblauDoETest1.d \
./UnitTests/himmelblauDoETest2.d \
./UnitTests/himmelblauDoETest3.d \
./UnitTests/himmelblauDoETest4.d \
./UnitTests/himmelblauDoETest5.d \
./UnitTests/himmelblauDoETest6.d \
./UnitTests/himmelblauDoETestConstraint1.d \
./UnitTests/himmelblauWithOneConstraint.d \
./UnitTests/linearTestFunction.d 

OBJS += \
./UnitTests/Constraint2.o \
./UnitTests/himmelblau.o \
./UnitTests/himmelblauAdjWithOneConstraint.o \
./UnitTests/himmelblauAdjoint.o \
./UnitTests/himmelblauDoETest1.o \
./UnitTests/himmelblauDoETest2.o \
./UnitTests/himmelblauDoETest3.o \
./UnitTests/himmelblauDoETest4.o \
./UnitTests/himmelblauDoETest5.o \
./UnitTests/himmelblauDoETest6.o \
./UnitTests/himmelblauDoETestConstraint1.o \
./UnitTests/himmelblauWithOneConstraint.o \
./UnitTests/linearTestFunction.o 


# Each subdirectory must supply rules for building sources it contributes
UnitTests/%.o: ../UnitTests/%.cpp UnitTests/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-UnitTests

clean-UnitTests:
	-$(RM) ./UnitTests/Constraint2.d ./UnitTests/Constraint2.o ./UnitTests/himmelblau.d ./UnitTests/himmelblau.o ./UnitTests/himmelblauAdjWithOneConstraint.d ./UnitTests/himmelblauAdjWithOneConstraint.o ./UnitTests/himmelblauAdjoint.d ./UnitTests/himmelblauAdjoint.o ./UnitTests/himmelblauDoETest1.d ./UnitTests/himmelblauDoETest1.o ./UnitTests/himmelblauDoETest2.d ./UnitTests/himmelblauDoETest2.o ./UnitTests/himmelblauDoETest3.d ./UnitTests/himmelblauDoETest3.o ./UnitTests/himmelblauDoETest4.d ./UnitTests/himmelblauDoETest4.o ./UnitTests/himmelblauDoETest5.d ./UnitTests/himmelblauDoETest5.o ./UnitTests/himmelblauDoETest6.d ./UnitTests/himmelblauDoETest6.o ./UnitTests/himmelblauDoETestConstraint1.d ./UnitTests/himmelblauDoETestConstraint1.o ./UnitTests/himmelblauWithOneConstraint.d ./UnitTests/himmelblauWithOneConstraint.o ./UnitTests/linearTestFunction.d ./UnitTests/linearTestFunction.o

.PHONY: clean-UnitTests

