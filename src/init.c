
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* TMB-generated native routines */
extern SEXP MakeADFunObject(SEXP, SEXP, SEXP, SEXP);
extern SEXP FreeADFunObject(SEXP);
extern SEXP InfoADFunObject(SEXP);
extern SEXP tmbad_print(SEXP, SEXP);
extern SEXP EvalADFunObject(SEXP, SEXP, SEXP);
extern SEXP TransformADFunObject(SEXP, SEXP);
extern SEXP MakeDoubleFunObject(SEXP, SEXP, SEXP, SEXP);
extern SEXP EvalDoubleFunObject(SEXP, SEXP, SEXP);
extern SEXP getParameterOrder(SEXP, SEXP, SEXP, SEXP);
extern SEXP MakeADGradObject(SEXP, SEXP, SEXP, SEXP);
extern SEXP MakeADHessObject2(SEXP, SEXP, SEXP, SEXP);
extern SEXP getFramework(void);
extern SEXP getSetGlobalPtr(SEXP);
extern SEXP TMBconfig(SEXP, SEXP);

/* Register .Call routines */
static const R_CallMethodDef CallEntries[] = {
    {"MakeADFunObject", (DL_FUNC) &MakeADFunObject, 4},
    {"FreeADFunObject", (DL_FUNC) &FreeADFunObject, 1},
    {"InfoADFunObject", (DL_FUNC) &InfoADFunObject, 1},
    {"tmbad_print", (DL_FUNC) &tmbad_print, 2},
    {"EvalADFunObject", (DL_FUNC) &EvalADFunObject, 3},
    {"TransformADFunObject", (DL_FUNC) &TransformADFunObject, 2},
    {"MakeDoubleFunObject", (DL_FUNC) &MakeDoubleFunObject, 4},
    {"EvalDoubleFunObject", (DL_FUNC) &EvalDoubleFunObject, 3},
    {"getParameterOrder", (DL_FUNC) &getParameterOrder, 4},
    {"MakeADGradObject", (DL_FUNC) &MakeADGradObject, 4},
    {"MakeADHessObject2", (DL_FUNC) &MakeADHessObject2, 4},
    {"getFramework", (DL_FUNC) &getFramework, 0},
    {"getSetGlobalPtr", (DL_FUNC) &getSetGlobalPtr, 1},
    {"TMBconfig", (DL_FUNC) &TMBconfig, 2},
    {NULL, NULL, 0}
};

void R_init_eDNAModel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
