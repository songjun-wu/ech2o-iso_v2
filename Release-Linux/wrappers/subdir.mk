################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../wrappers/CalculateBudgets.cpp \
../wrappers/CreateWorld.cpp \
../wrappers/CrunchWorld.cpp \
../wrappers/Report2maps.cpp \
../wrappers/Report2nc.cpp \
../wrappers/ReportAll2nc.cpp \
../wrappers/Report2screen.cpp \
../wrappers/SolveTimeStep.cpp \
../wrappers/Splash.cpp \
../wrappers/globals.cpp 

OBJS += \
./wrappers/CalculateBudgets.o \
./wrappers/CreateWorld.o \
./wrappers/CrunchWorld.o \
./wrappers/Report2maps.o \
./wrappers/Report2nc.o \
./wrappers/ReportAll2nc.o \
./wrappers/Report2screen.o \
./wrappers/SolveTimeStep.o \
./wrappers/Splash.o \
./wrappers/globals.o 

CPP_DEPS += \
./wrappers/CalculateBudgets.d \
./wrappers/CreateWorld.d \
./wrappers/CrunchWorld.d \
./wrappers/Report2maps.d \
./wrappers/Report2nc.d \
./wrappers/ReportAll2nc.d \
./wrappers/Report2screen.d \
./wrappers/SolveTimeStep.d \
./wrappers/Splash.d \
./wrappers/globals.d 


# Each subdirectory must supply rules for building sources it contributes
#for the UFZ-eve cluster, additional locations of libraries are specified in CFLAGS
wrappers/%.o: ../wrappers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


