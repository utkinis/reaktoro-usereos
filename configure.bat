@echo off

set Reaktoro_ROOT=C:/Users/krjac/anaconda3/envs/reaktoro/Library/lib/cmake/Reaktoro
set ThermoFun_ROOT=C:/Users/krjac/anaconda3/envs/reaktoro/Library/lib/cmake/ThermoFun
set Boost_ROOT=C:/Users/krjac/anaconda3/envs/reaktoro/Library/lib/cmake/Boost-1.70.0
set nlohmann_json_ROOT=C:/Users/krjac/anaconda3/envs/reaktoro/Library/lib/cmake/nlohmann_json

cmake -G "Visual Studio 16 2019" -A x64 -S . -B build