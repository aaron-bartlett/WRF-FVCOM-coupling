#Makefile include paths.mk

#Prepsource related

LIB_RAMS_PATH=../aux_src/utils
RAMS_PATH    =../aux_src/brams
PREPSOURCE_SRC=./
LIB_WPS_PATH=../aux_src/wps
ifeq ($(CHEM), RADM_WRF_FIM)
AER_DIR=SIMPLE
endif

