###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# MODULE_DIR       - Location of the MOOSE modules directory
# FRAMEWORK_DIR    - Location of the MOOSE framework
#
###############################################################################
THIS_MODULE_DIR    ?= $(shell dirname `pwd`)
MODULE_DIR         ?= $(shell dirname $(THIS_MODULE_DIR))
MOOSE_DIR          ?= $(shell dirname $(MODULE_DIR))
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
ADDITIONAL_CPPFLAGS += -Wall -Wextra -DMOOSE_UNIT_TEST
###############################################################################

# Extra stuff for GTEST
ADDITIONAL_INCLUDES	:= -isystem $(FRAMEWORK_DIR)/contrib/gtest
ADDITIONAL_LIBS 	:= $(FRAMEWORK_DIR)/contrib/gtest/libgtest.la

# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
SOLID_PROPERTIES := yes
include           $(MOOSE_DIR)/modules/modules.mk
###############################################################################

APPLICATION_DIR  := $(shell pwd)
APPLICATION_NAME := solid_properties-unit
BUILD_EXEC       := yes
app_BASE_DIR     :=
DEP_APPS    ?= $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include $(FRAMEWORK_DIR)/app.mk

# Find all the MOOSE unit test source files and include their dependencies.
solid_properties_unit_srcfiles := $(shell find $(APPLICATION_DIR)/src -name "*.C")
solid_properties_unit_deps := $(patsubst %.C, %.$(obj-suffix).d, $(solid_properties_unit_srcfiles))
-include $(solid_properties_unit_deps)

###############################################################################
# Additional special case targets should be added here
