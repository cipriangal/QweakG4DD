# Find Subversion 
find_package(Subversion REQUIRED)
get_filename_component(REAL_SOURCE_DIR "${SOURCE_DIR}" REALPATH)
Subversion_WC_INFO(${REAL_SOURCE_DIR} QWGEANT4)
# Write version info to local text file
file(WRITE QweakSimVersion.txt "// This file is generated automatically.  All changes will be lost.\n")
file(APPEND QweakSimVersion.txt "#define QWEAKSIM_URL		\"${QWGEANT4_WC_URL}\"\n")
file(APPEND QweakSimVersion.txt "#define QWEAKSIM_ROOT		\"${QWGEANT4_WC_ROOT}\"\n")
file(APPEND QweakSimVersion.txt "#define QWEAKSIM_REVISION	\"${QWGEANT4_WC_REVISION}\"\n")
#file(APPEND QweakSimVersion.txt "#define QWEAKSIM_INFO		\"${QWGEANT4_WC_INFO}\"\n")
file(APPEND QweakSimVersion.txt "#define QWEAKSIM_LAST_CHANGED_AUTHOR	\"${QWGEANT4_WC_LAST_CHANGED_AUTHOR}\"\n")
file(APPEND QweakSimVersion.txt "#define QWEAKSIM_LAST_CHANGED_DATE	\"${QWGEANT4_WC_LAST_CHANGED_DATE}\"\n")
file(APPEND QweakSimVersion.txt "#define QWEAKSIM_LAST_CHANGED_REV	\"${QWGEANT4_WC_LAST_CHANGED_REV}\"\n")
#file(APPEND QweakSimVersion.txt "#define QWEAKSIM_LAST_CHANGED_LOG	\"${QWGEANT4_LAST_CHANGED_LOG}\"\n")
# Copy only when file changes
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
	QweakSimVersion.txt QweakSimVersion.hh)
