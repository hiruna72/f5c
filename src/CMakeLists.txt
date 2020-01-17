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
set_source_files_properties(nanopolish_fast5_io.c PROPERTIES LANGUAGE CXX)

# include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/3rdparty/include)
include_directories(${PROJECT_SOURCE_DIR}/3rdparty/include/htslib)

set (MAIN main.c)# f5cmisc.h error.h)
set (METH_MAIN meth_main.c)# f5c.h fast5lite.h f5cmisc.h logsum.h)
set (F5C f5c.c)# f5c.h fast5lite.h f5cmisc.h)
set (EVENTS events.c)# f5c.h fast5lite.h f5cmisc.h nanopolish_read_db.h ksort.h)
set (NANOPOLISH_READ nanopolish_read_db.c)# nanopolish_read_db.h)
set (NANOPOLISH_INDEX nanopolish_index.c)# nanopolish_read_db.h fast5lite.h)
set (MODEL model.c)# model.h f5c.h fast5lite.h f5cmisc.h)
set (ALIGN align.c)# f5c.h fast5lite.h)
set (METH meth.c)# f5c.h fast5lite.h f5cmisc.h)
set (HMM hmm.c)# f5c.h fast5lite.h f5cmisc.h matrix.h logsum.h)
set (FREQ freq.c)# khash.h)
set (EVENT_ALIGN eventalign.c)#)
set (NANOPOLISH_FAST5_IO nanopolish_fast5_io.c)

# Build a shared lib
add_library(f5cshared SHARED ${MAIN} ${METH_MAIN} ${F5C} ${EVENTS} ${NANOPOLISH_READ} ${NANOPOLISH_INDEX} ${MODEL} ${ALIGN} ${METH} ${HMM} ${FREQ} ${EVENT_ALIGN} ${NANOPOLISH_FAST5_IO})


# htslib dependency
set (htslib_shared "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/libhts.so")

# hdf5 dependency
set (hdf5_static "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/libhdf5.a")

# link with dependencies
IF (${DEPLOY_PLATFORM} STREQUAL "x86")
	target_link_libraries(f5cshared ${htslib_shared} ${hdf5_static} -lz  -ldl -lpthread)
ELSE()
	target_link_libraries(f5cshared ${htslib_shared} ${hdf5_static}  -lz -ldl)

ENDIF()