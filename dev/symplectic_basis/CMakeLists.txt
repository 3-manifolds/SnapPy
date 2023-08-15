# To compile project 
# - run 'cmake -S . -B build/' in the dev/symplectic_basis dir 
# - run 'cmake --build build/'

# For clangd copy 'dev/symplectic_basis/build/compile_commands.json' to 'build/'

cmake_minimum_required(VERSION 3.12)

project(SnapPea)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Set the name of your program
set(YOUR_PROGRAM symplectic_basis)

# Set the location of the SnapPea kernel code
set(SNAPPEA_KERNEL ${PROJECT_SOURCE_DIR}/../../kernel)

# Compiler options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -fanalyzer")

# Add include directories
include_directories(${SNAPPEA_KERNEL}/addl_code ${SNAPPEA_KERNEL}/headers ${SNAPPEA_KERNEL}/real_type ${SNAPPEA_KERNEL}/unix_kit ${SNAPPEA_KERNEL}/)

# Add source files
file(GLOB KERNEL_SOURCES ${SNAPPEA_KERNEL}/*/*.c)
file(GLOB HEADER_FILES ${SNAPPEA_KERNEL}/*/*.h)
file(COPY ${CMAKE_SOURCE_DIR}/CuspedCensusData DESTINATION ${CMAKE_SOURCE_DIR}/build)      # Copy census data to be used by program

# Create executable target
add_executable(${YOUR_PROGRAM} ${YOUR_PROGRAM}_main.c ${KERNEL_SOURCES} ${HEADER_FILES} ${CUSPED_CENSUS})

# Link math library
target_link_libraries(${YOUR_PROGRAM} m)

# Define custom target to create BuildDate file
#add_custom_target(BuildDate ALL
#        COMMAND date > ${CMAKE_SOURCE_DIR}/BuildDate
#        DEPENDS ${KERNEL_SOURCES} ${HEADER_FILES}
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        )

# testing binary
add_executable(${YOUR_PROGRAM}_test ${YOUR_PROGRAM}_test.c ${KERNEL_SOURCES} ${HEADER_FILES} ${CUSPED_CENSUS})
target_link_libraries(${YOUR_PROGRAM}_test m)

# enable testing functionality
enable_testing()

# define tests
add_test(
        NAME symplectic_test
        COMMAND $<TARGET_FILE:${YOUR_PROGRAM}_test>
)

# Clean target
add_custom_target(clean-all
        COMMAND ${CMAKE_COMMAND} -E remove ${YOUR_PROGRAM} ${YOUR_PROGRAM}.o
        COMMAND ${CMAKE_COMMAND} -E remove_directory KernelObjects
        COMMAND ${CMAKE_COMMAND} -E remove *.pyc
        )

