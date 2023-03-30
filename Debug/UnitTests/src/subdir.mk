################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../UnitTests/src/aggregation_test.cpp \
../UnitTests/src/auxiliary_functions_test.cpp \
../UnitTests/src/bounds_test.cpp \
../UnitTests/src/configkey_test.cpp \
../UnitTests/src/constraint_functions_test.cpp \
../UnitTests/src/design_test.cpp \
../UnitTests/src/driver_test.cpp \
../UnitTests/src/gek_test.cpp \
../UnitTests/src/kriging_test.cpp \
../UnitTests/src/lhs_test.cpp \
../UnitTests/src/matrix_vector_operations_test.cpp \
../UnitTests/src/metric_test.cpp \
../UnitTests/src/multi_level_model_test.cpp \
../UnitTests/src/objective_function_test.cpp \
../UnitTests/src/optimization_test.cpp \
../UnitTests/src/output_test.cpp \
../UnitTests/src/polynomials_test.cpp \
../UnitTests/src/standard_test_functions_test.cpp \
../UnitTests/src/surrogate_model_data_test.cpp \
../UnitTests/src/surrogate_model_tester_test.cpp \
../UnitTests/src/test_functions_test.cpp 

CPP_DEPS += \
./UnitTests/src/aggregation_test.d \
./UnitTests/src/auxiliary_functions_test.d \
./UnitTests/src/bounds_test.d \
./UnitTests/src/configkey_test.d \
./UnitTests/src/constraint_functions_test.d \
./UnitTests/src/design_test.d \
./UnitTests/src/driver_test.d \
./UnitTests/src/gek_test.d \
./UnitTests/src/kriging_test.d \
./UnitTests/src/lhs_test.d \
./UnitTests/src/matrix_vector_operations_test.d \
./UnitTests/src/metric_test.d \
./UnitTests/src/multi_level_model_test.d \
./UnitTests/src/objective_function_test.d \
./UnitTests/src/optimization_test.d \
./UnitTests/src/output_test.d \
./UnitTests/src/polynomials_test.d \
./UnitTests/src/standard_test_functions_test.d \
./UnitTests/src/surrogate_model_data_test.d \
./UnitTests/src/surrogate_model_tester_test.d \
./UnitTests/src/test_functions_test.d 

OBJS += \
./UnitTests/src/aggregation_test.o \
./UnitTests/src/auxiliary_functions_test.o \
./UnitTests/src/bounds_test.o \
./UnitTests/src/configkey_test.o \
./UnitTests/src/constraint_functions_test.o \
./UnitTests/src/design_test.o \
./UnitTests/src/driver_test.o \
./UnitTests/src/gek_test.o \
./UnitTests/src/kriging_test.o \
./UnitTests/src/lhs_test.o \
./UnitTests/src/matrix_vector_operations_test.o \
./UnitTests/src/metric_test.o \
./UnitTests/src/multi_level_model_test.o \
./UnitTests/src/objective_function_test.o \
./UnitTests/src/optimization_test.o \
./UnitTests/src/output_test.o \
./UnitTests/src/polynomials_test.o \
./UnitTests/src/standard_test_functions_test.o \
./UnitTests/src/surrogate_model_data_test.o \
./UnitTests/src/surrogate_model_tester_test.o \
./UnitTests/src/test_functions_test.o 


# Each subdirectory must supply rules for building sources it contributes
UnitTests/src/%.o: ../UnitTests/src/%.cpp UnitTests/src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-UnitTests-2f-src

clean-UnitTests-2f-src:
	-$(RM) ./UnitTests/src/aggregation_test.d ./UnitTests/src/aggregation_test.o ./UnitTests/src/auxiliary_functions_test.d ./UnitTests/src/auxiliary_functions_test.o ./UnitTests/src/bounds_test.d ./UnitTests/src/bounds_test.o ./UnitTests/src/configkey_test.d ./UnitTests/src/configkey_test.o ./UnitTests/src/constraint_functions_test.d ./UnitTests/src/constraint_functions_test.o ./UnitTests/src/design_test.d ./UnitTests/src/design_test.o ./UnitTests/src/driver_test.d ./UnitTests/src/driver_test.o ./UnitTests/src/gek_test.d ./UnitTests/src/gek_test.o ./UnitTests/src/kriging_test.d ./UnitTests/src/kriging_test.o ./UnitTests/src/lhs_test.d ./UnitTests/src/lhs_test.o ./UnitTests/src/matrix_vector_operations_test.d ./UnitTests/src/matrix_vector_operations_test.o ./UnitTests/src/metric_test.d ./UnitTests/src/metric_test.o ./UnitTests/src/multi_level_model_test.d ./UnitTests/src/multi_level_model_test.o ./UnitTests/src/objective_function_test.d ./UnitTests/src/objective_function_test.o ./UnitTests/src/optimization_test.d ./UnitTests/src/optimization_test.o ./UnitTests/src/output_test.d ./UnitTests/src/output_test.o ./UnitTests/src/polynomials_test.d ./UnitTests/src/polynomials_test.o ./UnitTests/src/standard_test_functions_test.d ./UnitTests/src/standard_test_functions_test.o ./UnitTests/src/surrogate_model_data_test.d ./UnitTests/src/surrogate_model_data_test.o ./UnitTests/src/surrogate_model_tester_test.d ./UnitTests/src/surrogate_model_tester_test.o ./UnitTests/src/test_functions_test.d ./UnitTests/src/test_functions_test.o

.PHONY: clean-UnitTests-2f-src

