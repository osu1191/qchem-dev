file(GLOB SRC_LIBPEQS_C RELATIVE
    "${QCHEM_SOURCE_DIR}/libpeqs" "${QCHEM_SOURCE_DIR}/libpeqs/*.C")
file(GLOB SRC_LIBPEQS_F RELATIVE
    "${QCHEM_SOURCE_DIR}/libpeqs" "${QCHEM_SOURCE_DIR}/libpeqs/*.F")

add_library(libpeqs STATIC ${SRC_LIBPEQS_C} ${SRC_LIBPEQS_F})
add_dependencies(libpeqs generated_headers)
target_compile_definitions(libpeqs PRIVATE _USE_MATH_DEFINES)
set_target_properties(libpeqs PROPERTIES OUTPUT_NAME peqs)
