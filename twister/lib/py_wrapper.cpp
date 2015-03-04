// Build using: python setup.py build

#include <Python.h>
#include <string>
#include "manifold.h"
#include "twister.h"
#include "global.h"

static char twister_build_bundle_doc[] = "Usage: construct_bundle(name, surface, monodromy, optimise, peripheral_curves, warnings, debugging_level). \n \
Returns the triangluation of the mapping torus with monodromy <monodromy> and fiber <surface> along with the contents of the GLOBAL_message_stream. If an error occurs (None, GLOBAL_message_stream) is returned. ";

#if PY_MAJOR_VERSION < 3
	#define string_convert PyString_AsString
	#define string_deconvert PyString_FromString
	
#else
	#define string_convert PyBytes_AsString
	#define string_deconvert PyBytes_FromString
#endif

extern "C" PyObject* twister_version(PyObject *self, PyObject *args)
{
	return string_deconvert((char *) version.c_str());
}

extern "C" PyObject* twister_build_bundle(PyObject *self, PyObject *args)
{
	const char *char_name, *char_surface, *char_monodromy;
	bool bool_optimise, bool_peripheral_curves, bool_warnings;
	int int_debugging_level;
	
	if (!PyArg_ParseTuple(args, "sssbbbi", &char_name, &char_surface, &char_monodromy, &bool_optimise, &bool_peripheral_curves, &bool_warnings, &int_debugging_level))
		return NULL;
	
	std::string manifold_name = char_name, surface_file_contents = char_surface, monodromy = char_monodromy;
	std::string manifold_contents = "";
	set_globals_to_defaults();
	GLOBAL_warnings = bool_warnings;
	GLOBAL_optimise = bool_optimise;
	GLOBAL_calculate_peripheral_curves = bool_peripheral_curves;
	GLOBAL_debugging_level = int_debugging_level;
	
	try
	{
		manifold M(manifold_name, bundle);  // Build a (blank) manifold with the correct name & type.
		construct_manifold(M, surface_file_contents, monodromy, "");
		manifold_contents = M.to_string();
	}
	catch (...)
	{
		return Py_BuildValue("ss", NULL, (char *) GLOBAL_message_stream.c_str());
	}
	
	return Py_BuildValue("ss", (char *) manifold_contents.c_str(), (char *) GLOBAL_message_stream.c_str());
}

static char twister_build_splitting_doc[] = "Usage: construct_splitting(name, surface, gluing, handles, optimise, peripheral_curves, warnings, debugging_level). \n \
Returns the triangluation of the Heegaard splitting over the surface <surface> made by using gluing <gluing> \
and attaching handles according to <handles> along with the contents of the GLOBAL_message_stream. If an error occurs (None, GLOBAL_message_stream) is returned. ";

extern "C" PyObject* twister_build_splitting(PyObject *self, PyObject *args)
{
	const char *char_name, *char_surface, *char_gluing, *char_handles;
	bool bool_optimise, bool_peripheral_curves, bool_warnings;
	int int_debugging_level;
	
	if (!PyArg_ParseTuple(args, "ssssbbbi", &char_name, &char_surface, &char_gluing, &char_handles, &bool_optimise, &bool_peripheral_curves, &bool_warnings, &int_debugging_level))
		return NULL;
	
	std::string manifold_name = char_name, surface_file_contents = char_surface, gluing = char_gluing, handles = char_handles;
	std::string manifold_contents = "";
	set_globals_to_defaults();
	GLOBAL_warnings = bool_warnings;
	GLOBAL_optimise = bool_optimise;
	GLOBAL_calculate_peripheral_curves = bool_peripheral_curves;
	GLOBAL_debugging_level = int_debugging_level;
	
	try
	{
		manifold M(manifold_name, splitting);  // Build a (blank) manifold with the correct name & type.
		construct_manifold(M, surface_file_contents, gluing, handles);
		manifold_contents = M.to_string();
	}
	catch (...)
	{
		return Py_BuildValue("ss", NULL, (char *) GLOBAL_message_stream.c_str());
	}
	
	return Py_BuildValue("ss", (char *) manifold_contents.c_str(), (char *) GLOBAL_message_stream.c_str());
}


static PyMethodDef twister_methods[] = {
	{"build_bundle", twister_build_bundle, METH_VARARGS, twister_build_bundle_doc},
	{"build_splitting", twister_build_splitting, METH_VARARGS, twister_build_splitting_doc},
	{"twister_version", twister_version, METH_NOARGS, "Version number as string."},
	{NULL, NULL}
};

static char twister_doc[] = "This module is for interfacing with Twister.  It provides two functions: build_bundle() and build_splitting().";

#if PY_MAJOR_VERSION < 3
	PyMODINIT_FUNC inittwister_core(void) 
	{
		Py_InitModule3("twister_core", twister_methods, twister_doc);
	}
#else
	static struct PyModuleDef twister_core_module = {PyModuleDef_HEAD_INIT, "twister_core", twister_doc, -1, twister_methods};
	// -1: the module keeps state in global variables.
	
	extern "C" PyMODINIT_FUNC PyInit_twister_core(void)
	{
		return PyModule_Create(&twister_core_module);
	}
#endif
