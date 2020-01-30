/**
 * Form.hpp introduce struct FormGrade and acronyms Form, Hodge, and Derivative.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _FORM_HPP_INCLUDED_
#define _FORM_HPP_INCLUDED_

#include "Column.hpp"

namespace gfd
{

enum FormGrade
{
	fg_prim0, fg_dual0,
	fg_prim1, fg_dual1,
	fg_prim2, fg_dual2,
	fg_prim3, fg_dual3,
	fg_prim4, fg_dual4,
	fg_num
};
inline uint FormGradeDimension(const FormGrade grade) { return uint(grade >> 1); }
inline bool FormGradeIsPrim(const FormGrade grade) { return bool(~grade & 1); }
inline bool FormGradeIsDual(const FormGrade grade) { return bool(grade & 1); }
inline FormGrade FormGradeDual(const FormGrade grade) { return FormGrade(grade ^ 1); }
inline FormGrade FormGradeDeriv(const FormGrade grade) { return FormGrade(FormGradeIsPrim(grade) ? grade + 2 : grade - 2); }
inline uint FormGradeVectorDimension(const FormGrade grade, const uint dim) { 
	switch (FormGradeDimension(grade)) {
	case 1: return dim;
	case 2: return dim * (dim - 1) / 2;
	case 3: return dim * (dim - 1) * (dim - 2) / 6;
	default: return 1; // cases 0 and 4
	}
}

template <typename T> using Form = Column<T>;
template <typename T> using Hodge = Diagonal<T>;
using Derivative = Sparse<sign>;


}

#endif //_DEC_HPP_INCLUDED_
