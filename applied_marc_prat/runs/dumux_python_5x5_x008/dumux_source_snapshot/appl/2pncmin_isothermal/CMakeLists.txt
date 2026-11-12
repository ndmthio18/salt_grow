# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES params.input)

# Isothermal variants of the 2pncmin application (energy balance disabled).
# Same physics as appl/2pncmin otherwise (brine/water, precipitation, theta reg.).
dune_add_test(NAME pnm_2pncmin
              SOURCES main.cc
              COMPILE_DEFINITIONS ISOTHERMAL=1 USETHETAREGULARIZATION=0 USEBRINE=1)

dune_add_test(NAME pnm_2pncmin_reg
              SOURCES main.cc
              COMPILE_DEFINITIONS ISOTHERMAL=1 USETHETAREGULARIZATION=1 USEBRINE=1)

dune_add_test(NAME pnm_2pncmin_water_reg
              SOURCES main.cc
              COMPILE_DEFINITIONS ISOTHERMAL=1 USETHETAREGULARIZATION=1 USEBRINE=0)

set(CMAKE_BUILD_TYPE Debug)
