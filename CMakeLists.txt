cmake_minimum_required(VERSION 3.8)
project(FFTConvolver)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        test/Test.cpp
        AudioFFT.cpp
        AudioFFT.h
        FFTConvolver.cpp
        FFTConvolver.h
        TwoStageFFTConvolver.cpp
        TwoStageFFTConvolver.h
        Utilities.cpp
        Utilities.h BinauralFFTConvolver.cpp BinauralFFTConvolver.h)

add_executable(FFTConvolver ${SOURCE_FILES})