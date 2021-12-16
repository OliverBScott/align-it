find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(py3lib "python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}")
find_package(Boost COMPONENTS ${py3lib} REQUIRED)

# Without this, any build libraries automatically have names "lib{x}.so"
set(CMAKE_SHARED_MODULE_PREFIX "")

python3_add_library(cpyalignit MODULE cpyalignit.cpp)
target_link_libraries(cpyalignit PRIVATE alignit_lib ${RDKIT_DYN_LIBRARIES} Boost::${py3lib})
set_target_properties(cpyalignit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../pyalignit)
