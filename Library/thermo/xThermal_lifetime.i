// Lifetime of non-owning results (Python).
//
// A pointer or reference returned by a method or by a member getter points into memory owned by the object it was
// obtained from: a field of a result struct (sw.UpdateState_TPX(T, p, X).Rho, sw.Slice_constP(p).lines), the water
// model owned by a mixture model (cH2ONaCl.get_pWater()), a reference to an internal table. SWIG wraps such results
// in proxies that do not own their memory, so without further action a proxy can outlive the object and then reads
// freed memory (crash or wrong values).
//
// Rule: the proxy of a non-owning result keeps its owner alive, by holding a reference to it (the same mechanism SWIG
// uses for references to container elements). Results that own their object (returned by value, constructors,
// %newobject) are not affected. Include this file in every module after its %{ ... %} header block and before the
// wrapped declarations.

#ifdef SWIGPYTHON
%fragment("xThermal_keep_owner", "header") {
// objs: the wrapper's argument array; objs[0] is the owner ("self" for methods and member getters).
SWIGINTERN void xThermal_keep_owner(PyObject* result, PyObject* const* objs) {
  PyObject* owner = objs ? objs[0] : 0;
  if (!result || !owner || result == Py_None || result == owner) return;
  // fast exit for plain Python values (e.g. a double returned by const reference from a vector element): only SWIG
  // objects and proxy instances (which have an instance dictionary) can wrap a pointer
  if (!SwigPyObject_Check(result) && Py_TYPE(result)->tp_dictoffset == 0) return;
  SwigPyObject* sobj = SWIG_Python_GetSwigThis(result);
  if (sobj && !(sobj->own & SWIG_POINTER_OWN)) {
    if (PyObject_SetAttrString(result, "_xThermal_owner", owner) != 0) PyErr_Clear();
  }
}
// Wrappers without arguments (default constructors) have no local swig_obj; this null fallback is used there.
// Wrappers with arguments declare their own swig_obj, which hides it.
static PyObject* const* const swig_obj = 0;
}

%typemap(ret, fragment="xThermal_keep_owner") SWIGTYPE *, SWIGTYPE &, const SWIGTYPE & %{
  xThermal_keep_owner($result, swig_obj);
%}
#endif
