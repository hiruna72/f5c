# Find all the files which need to be compiled
include_directories(${PROJECT_SOURCE_DIR}/3rdparty/include)
include_directories(${PROJECT_SOURCE_DIR}/3rdparty/include/htslib)

# until making f5c a C project use g++
set_source_files_properties(main.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(meth_main.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(f5c.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(events.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(nanopolish_read_db.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(nanopolish_index.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(model.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(align.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(meth.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(hmm.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(freq.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(eventalign.c PROPERTIES LANGUAGE CXX)

# creating config.h file
# add_custom_command does not create a new target. You have to define targets explicitly
# by add_executable, add_library or add_custom_target in order to make them visible to make
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/config.h
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/config.h
    COMMAND echo "#define HAVE_HDF5_H 1" > ${CMAKE_CURRENT_BINARY_DIR}/config.h

    # Display the given message before the commands are executed at build time
    COMMENT "Creating config.h"
    VERBATIM
)
# target zoo is always built
add_custom_target(creating_config ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/config.h
    # to make quotes printable,for example
    VERBATIM
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})


set (MAIN main.c f5cmisc.h error.h)
set (METH_MAIN meth_main.c f5c.h fast5lite.h f5cmisc.h logsum.h)
set (F5C f5c.c f5c.h fast5lite.h f5cmisc.h)
set (EVENTS events.c f5c.h fast5lite.h f5cmisc.h nanopolish_read_db.h ksort.h)
set (NANOPOLISH_READ nanopolish_read_db.c nanopolish_read_db.h)
set (NANOPOLISH_INDEX nanopolish_index.c nanopolish_read_db.h fast5lite.h)
set (MODEL model.c model.h f5c.h fast5lite.h f5cmisc.h)
set (ALIGN align.c f5c.h fast5lite.h)
set (METH meth.c f5c.h fast5lite.h f5cmisc.h)
set (HMM hmm.c f5c.h fast5lite.h f5cmisc.h matrix.h logsum.h)
set (FREQ freq.c khash.h)
set (EVENT_ALIGN eventalign.c)


# Build an executable
# add_executable(f5c_exe ${MAIN} ${METH_MAIN} ${F5C} ${EVENTS} ${NANOPOLISH_READ} ${NANOPOLISH_INDEX} ${MODEL} ${ALIGN} ${METH} ${HMM} ${FREQ} ${EVENT_ALIGN})

# Build a shared lib
add_library(f5cshared SHARED ${MAIN} ${METH_MAIN} ${F5C} ${EVENTS} ${NANOPOLISH_READ} ${NANOPOLISH_INDEX} ${MODEL} ${ALIGN} ${METH} ${HMM} ${FREQ} ${EVENT_ALIGN})


# htslib dependency
set (htslib_static "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/libhts.a")

# hdf5 dependency
set (hdf5_static "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/libhdf5.a")

# link with dependencies
# target_link_libraries(f5c_exe -static ${htslib_static} ${hdf5_static} -lpthread -lz -ldl)

# link with dependencies
IF (${DEPLOY_PLATFORM} STREQUAL "x86")
	target_link_libraries(f5cshared ${htslib_static} ${hdf5_static} -lz -lpthread -ldl)
ELSE()
	find_library(ANDROID_LOG_LIB log)
	target_link_libraries(f5cshared ${htslib_static} ${hdf5_static} ${ANDROID_LOG_LIB} -lz -ldl)

ENDIF()
