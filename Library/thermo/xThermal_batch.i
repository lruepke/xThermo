// Array export of the fields of a ThermodynamicPropertiesVector (stage 5a): one copy into a bytearray, which the
// Python layer (xthermo.tpx / xthermo.hpx) wraps with numpy.frombuffer. Python C-API only, no NumPy at build time.
// Included by the three modules after their declaration of ThermodynamicPropertiesVector.
#ifdef SWIGPYTHON
%extend ThermodynamicPropertiesVector {
    PyObject* _field_bytes(const std::string& name) {
        const std::vector<double>* d = nullptr;
#define XT_FIELD(f) if (name == #f) d = &self->f;
        XT_FIELD(T) XT_FIELD(p) XT_FIELD(X) XT_FIELD(H) XT_FIELD(S_l) XT_FIELD(S_v) XT_FIELD(S_h) XT_FIELD(X_l) XT_FIELD(X_v)
        XT_FIELD(Rho_l) XT_FIELD(Rho_v) XT_FIELD(Rho_h) XT_FIELD(H_l) XT_FIELD(H_v) XT_FIELD(H_h) XT_FIELD(Cp_l) XT_FIELD(Cp_v)
        XT_FIELD(Cp_h) XT_FIELD(Mu_l) XT_FIELD(Mu_v) XT_FIELD(Rho) XT_FIELD(Cp) XT_FIELD(Mu) XT_FIELD(dRhodP) XT_FIELD(dRhodP_l)
        XT_FIELD(dRhodP_v) XT_FIELD(dRhodT) XT_FIELD(dRhodT_l) XT_FIELD(dRhodT_v) XT_FIELD(IsothermalCompressibility)
        XT_FIELD(IsothermalCompressibility_l) XT_FIELD(IsothermalCompressibility_v) XT_FIELD(IsobaricExpansivity)
        XT_FIELD(IsobaricExpansivity_l) XT_FIELD(IsobaricExpansivity_v)
#undef XT_FIELD
        if (d) return PyByteArray_FromStringAndSize(reinterpret_cast<const char*>(d->data()), (Py_ssize_t)(d->size() * sizeof(double)));
        const std::vector<int>* k = name == "phase" ? &self->phase : name == "status" ? &self->status : nullptr;
        if (k) return PyByteArray_FromStringAndSize(reinterpret_cast<const char*>(k->data()), (Py_ssize_t)(k->size() * sizeof(int)));
        PyErr_SetString(PyExc_KeyError, name.c_str());
        return nullptr;
    }
};
#endif
