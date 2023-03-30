################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Rodeo.cpp \
../src/aggregation_model.cpp \
../src/auxiliary_functions.cpp \
../src/bounds.cpp \
../src/configkey.cpp \
../src/constraint_functions.cpp \
../src/correlation_functions.cpp \
../src/design.cpp \
../src/drivers.cpp \
../src/gek.cpp \
../src/kriging_training.cpp \
../src/lhs.cpp \
../src/linear_regression.cpp \
../src/matrix_vector_operations.cpp \
../src/metric.cpp \
../src/multi_level_model.cpp \
../src/objective_function.cpp \
../src/optimization.cpp \
../src/output.cpp \
../src/polynomials.cpp \
../src/random_functions.cpp \
../src/sgek.cpp \
../src/standard_test_functions.cpp \
../src/surrogate_model.cpp \
../src/surrogate_model_data.cpp \
../src/surrogate_model_tester.cpp \
../src/test_functions.cpp 

CPP_DEPS += \
./src/Rodeo.d \
./src/aggregation_model.d \
./src/auxiliary_functions.d \
./src/bounds.d \
./src/configkey.d \
./src/constraint_functions.d \
./src/correlation_functions.d \
./src/design.d \
./src/drivers.d \
./src/gek.d \
./src/kriging_training.d \
./src/lhs.d \
./src/linear_regression.d \
./src/matrix_vector_operations.d \
./src/metric.d \
./src/multi_level_model.d \
./src/objective_function.d \
./src/optimization.d \
./src/output.d \
./src/polynomials.d \
./src/random_functions.d \
./src/sgek.d \
./src/standard_test_functions.d \
./src/surrogate_model.d \
./src/surrogate_model_data.d \
./src/surrogate_model_tester.d \
./src/test_functions.d 

OBJS += \
./src/Rodeo.o \
./src/aggregation_model.o \
./src/auxiliary_functions.o \
./src/bounds.o \
./src/configkey.o \
./src/constraint_functions.o \
./src/correlation_functions.o \
./src/design.o \
./src/drivers.o \
./src/gek.o \
./src/kriging_training.o \
./src/lhs.o \
./src/linear_regression.o \
./src/matrix_vector_operations.o \
./src/metric.o \
./src/multi_level_model.o \
./src/objective_function.o \
./src/optimization.o \
./src/output.o \
./src/polynomials.o \
./src/random_functions.o \
./src/sgek.o \
./src/standard_test_functions.o \
./src/surrogate_model.o \
./src/surrogate_model_data.o \
./src/surrogate_model_tester.o \
./src/test_functions.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src

clean-src:
	-$(RM) ./src/Rodeo.d ./src/Rodeo.o ./src/aggregation_model.d ./src/aggregation_model.o ./src/auxiliary_functions.d ./src/auxiliary_functions.o ./src/bounds.d ./src/bounds.o ./src/configkey.d ./src/configkey.o ./src/constraint_functions.d ./src/constraint_functions.o ./src/correlation_functions.d ./src/correlation_functions.o ./src/design.d ./src/design.o ./src/drivers.d ./src/drivers.o ./src/gek.d ./src/gek.o ./src/kriging_training.d ./src/kriging_training.o ./src/lhs.d ./src/lhs.o ./src/linear_regression.d ./src/linear_regression.o ./src/matrix_vector_operations.d ./src/matrix_vector_operations.o ./src/metric.d ./src/metric.o ./src/multi_level_model.d ./src/multi_level_model.o ./src/objective_function.d ./src/objective_function.o ./src/optimization.d ./src/optimization.o ./src/output.d ./src/output.o ./src/polynomials.d ./src/polynomials.o ./src/random_functions.d ./src/random_functions.o ./src/sgek.d ./src/sgek.o ./src/standard_test_functions.d ./src/standard_test_functions.o ./src/surrogate_model.d ./src/surrogate_model.o ./src/surrogate_model_data.d ./src/surrogate_model_data.o ./src/surrogate_model_tester.d ./src/surrogate_model_tester.o ./src/test_functions.d ./src/test_functions.o

.PHONY: clean-src

