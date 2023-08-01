#ifndef __OPENCL_UTILS__
#define __OPENCL_UTILS__

//#define CL_TARGET_OPENCL_VERSION 120
#define CL_TARGET_OPENCL_VERSION 110

#ifdef __APPLE__
	#include <OpenCL/cl.hpp>
#else
	#include <CL/cl.hpp>
#endif

#ifndef CL_DEVICE_DOUBLE_FP_CONFIG 
    #error "Your OpenCL installation does not support the double data type"
#endif

// Format: (variable name (no quotes), script_path (with quotes))
#define include_opencl_file(script_name, script_path) \
    asm(#script_name ": .incbin \"" script_path "\""); \
    asm(#script_name "End: .byte 0x00"); \
    extern char script_name[];

#endif