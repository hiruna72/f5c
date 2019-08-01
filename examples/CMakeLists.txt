# Get all the directories with needed .hpp files
include_directories ("${PROJECT_SOURCE_DIR}/src")

set_source_files_properties(f5c_example.c PROPERTIES LANGUAGE CXX)


# Create the executable
add_executable(f5c_example f5c_example.c)

# Link the required libraries
target_link_libraries(f5c_example f5cshared)