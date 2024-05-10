cmake_minimum_required(VERSION 3.5)

project(lab_5 LANGUAGES C)

add_executable(lab_5 matrix.c)

include(GNUInstallDirs)
install(TARGETS lab_5
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
