#export SLEPC_DIR=/home/guido/libs/slepc-3.7.1
#export PETSC_DIR=/home/guido/libs/petsc-3.7.2
#export PETSC_ARCH=arch-linux2-c-opt

PWD:= $(shell pwd)

SRC_DIR=./src
OBJ_DIR=./obj
DEP_DIR=${PWD}/inc

CFLAGS=-g -O0 
	
DEPS = ${DEP_DIR}/bathe.h              \
       ${DEP_DIR}/global.h             \
       ${DEP_DIR}/types.h              \
       ${DEP_DIR}/list.h               \
       ${DEP_DIR}/gmsh.h               \
       ${DEP_DIR}/mesh.h               \
       ${DEP_DIR}/fem.h                \
       ${DEP_DIR}/utils.h              \
       ${DEP_DIR}/fun.h           

OBJ  = ${OBJ_DIR}/bth_main.o           \
       ${OBJ_DIR}/lst2msh.o            \
       ${OBJ_DIR}/bth_parser.o         \
       ${OBJ_DIR}/bth_init.o           \
       ${OBJ_DIR}/bth_force.o          \
       ${OBJ_DIR}/bth_calc_k.o         \
       ${OBJ_DIR}/bth_property.o       \
       ${OBJ_DIR}/bth_fini.o           \
       ${OBJ_DIR}/bth_elemental.o      \
       ${OBJ_DIR}/bth_strain.o         \
       ${OBJ_DIR}/bth_boundary.o       \
       ${OBJ_DIR}/bth_evolute.o        \
       ${OBJ_DIR}/bth_newrap.o         \
       ${OBJ_DIR}/bth_arclength_1.o    \
       ${OBJ_DIR}/bth_arclength_2.o    \
       ${OBJ_DIR}/bth_arclength_3.o    \
       ${OBJ_DIR}/bth_arclength_4.o    \
       ${OBJ_DIR}/bth_output.o         \
       ${OBJ_DIR}/utils.o              \
       ${OBJ_DIR}/list.o               \
       ${OBJ_DIR}/gmsh.o               \
       ${OBJ_DIR}/mesh.o               \
       ${OBJ_DIR}/fem.o                \
       ${OBJ_DIR}/fun.o           
       

.PHONY: clean_
	
all: ${OBJ} 
	gcc -o bathe $^ ${SLEPC_EPS_LIB} -lgsl -lgslcblas -lm
	
${OBJ_DIR}/%.o: ${SRC_DIR}/%.c $(DEPS) 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I${DEP_DIR} 

clean_:	    
	rm -f $(OBJ) bathe

#include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common



