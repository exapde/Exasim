# Install script for directory: /Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/programs

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/build/Darwin-arm64/programs/pm_ptest")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_ptest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_ptest")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/GKlib/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/METIS/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_ptest")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_ptest")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/build/Darwin-arm64/programs/pm_mtest")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_mtest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_mtest")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/GKlib/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/METIS/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_mtest")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_mtest")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/build/Darwin-arm64/programs/pm_parmetis")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_parmetis" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_parmetis")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/GKlib/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/METIS/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_parmetis")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_parmetis")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/build/Darwin-arm64/programs/pm_pometis")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_pometis" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_pometis")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/GKlib/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/METIS/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_pometis")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_pometis")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/build/Darwin-arm64/programs/pm_dglpart")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_dglpart" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_dglpart")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/GKlib/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/METIS/lib"
      -delete_rpath "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_dglpart")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pm_dglpart")
    endif()
  endif()
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/cuongnguyen/Documents/GitHub/Exasim/text2code/ParMETIS/build/Darwin-arm64/programs/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
