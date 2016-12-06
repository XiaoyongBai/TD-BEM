FILE(GLOB src_list *.cpp *.h)

add_executable(TD_BEM ${src_list})


target_link_libraries(TD_BEM libpetsc.dylib libmpi.dylib)




