cmake_minimum_required(VERSION 3.20)

project(ReaktoroUsereos CXX)

find_package(Reaktoro REQUIRED)
find_package(ThermoFun REQUIRED)
find_package(nlohmann_json REQUIRED)

add_library(USEREOS SHARED src/reaktoro-usereos.cpp src/reaktoro-usereos.hpp)
target_link_libraries(USEREOS PUBLIC Reaktoro::Reaktoro ThermoFun::ThermoFun)
target_include_directories(USEREOS PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
                                          "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

add_custom_command(TARGET USEREOS POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:USEREOS>
                                   $<TARGET_FILE_DIR:USEREOS>
  COMMAND_EXPAND_LISTS
)

add_executable(test-usereos src/main.cpp)
target_link_libraries(test-usereos USEREOS)

add_custom_command(TARGET test-usereos POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:test-usereos>
                                   $<TARGET_FILE_DIR:test-usereos>
  COMMAND_EXPAND_LISTS
)
