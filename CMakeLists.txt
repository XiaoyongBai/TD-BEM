cmake_minimum_required(VERSION 2.8)

#project name
PROJECT(TD_BEM)

message(STATUS "This is BINARY dir " ${TD_BEM_BINARY_DIR})
message(STATUS "This is SOURCE dir " ${TD_BEM_SOURCE_DIR})


#head file path
INCLUDE_DIRECTORIES(
/Users/xiaoyongbai/Workspace/petsc/include
/Users/xiaoyongbai/Workspace/petsc/arch-darwin-c-debug/include
)


#library file path
link_directories(
/Users/xiaoyongbai/Workspace/petsc/arch-darwin-c-debug/lib
)



#source directory
add_subdirectory(src)




