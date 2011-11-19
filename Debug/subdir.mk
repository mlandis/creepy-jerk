################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../CondLikes.cpp \
../Expression.cpp \
../FileMgr.cpp \
../MbBitfield.cpp \
../MbRandom.cpp \
../Mcmc.cpp \
../Model.cpp \
../Parm.cpp \
../Parm_alpha.cpp \
../Parm_lambda.cpp \
../Parm_sigma.cpp \
../Parm_tau.cpp \
../Parm_tree.cpp \
../Patron.cpp \
../Settings.cpp \
../Table.cpp \
../main.cpp 

OBJS += \
./CondLikes.o \
./Expression.o \
./FileMgr.o \
./MbBitfield.o \
./MbRandom.o \
./Mcmc.o \
./Model.o \
./Parm.o \
./Parm_alpha.o \
./Parm_lambda.o \
./Parm_sigma.o \
./Parm_tau.o \
./Parm_tree.o \
./Patron.o \
./Settings.o \
./Table.o \
./main.o 

CPP_DEPS += \
./CondLikes.d \
./Expression.d \
./FileMgr.d \
./MbBitfield.d \
./MbRandom.d \
./Mcmc.d \
./Model.d \
./Parm.d \
./Parm_alpha.d \
./Parm_lambda.d \
./Parm_sigma.d \
./Parm_tau.d \
./Parm_tree.d \
./Patron.d \
./Settings.d \
./Table.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


