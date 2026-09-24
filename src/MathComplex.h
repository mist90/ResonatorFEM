#ifndef MATHCOMPLEX_H
#define MATHCOMPLEX_H

/*****************************************************************************
 *        Класс для работы с комплексными числами. Типы: float и double      *
 ****************************************************************************/

#include <math.h>

template <class T>
class MathComplex
{
public:
    MathComplex(T re,T im = 0);
    MathComplex();
    void           setExponentialArg(T amplitude, T phase);
    MathComplex<T> operator+(const MathComplex<T> &c1);
    MathComplex<T> operator+(const T &c1);
    MathComplex<T> operator-(const MathComplex<T> &c1);
    MathComplex<T> operator-(const T &c1);
    MathComplex<T> operator-();
    MathComplex<T> operator*(const MathComplex<T> &c1);
    MathComplex<T> operator*(const T &c1);
    MathComplex<T> operator/(const MathComplex<T> &c1);
    MathComplex<T> operator/(const T &c1);
	void           operator=(const MathComplex<T> &c1);
    void           operator=(const T &c1);
	void           operator+=(const MathComplex<T> &c1);
	void           operator-=(const MathComplex<T> &c1);
	void           operator*=(const MathComplex<T> &c1);
	void           operator/=(const MathComplex<T> &c1);
    bool           operator==(const MathComplex<T> &c1);
    bool           operator!=(const MathComplex<T> &c1);
    bool           operator==(const T &c1);
    bool           operator!=(const T &c1);
    T&             re();
    T&             im();
    T              phase();
    MathComplex<T> powc(const MathComplex<T> &p);
    T              absc();
	MathComplex<T> sinc();
	MathComplex<T> cosc();
	MathComplex<T> tanc();
	MathComplex<T> sinhc();
	MathComplex<T> coshc();
	MathComplex<T> tanhc();
	MathComplex<T> asinc();
	MathComplex<T> acosc();
	MathComplex<T> atanc();
	MathComplex<T> expc();
	MathComplex<T> logc();
	MathComplex<T> log10c();
	MathComplex<T> sqrtc();
private:
    T _re,_im;
};

//--------------------------------------------------

template <class T>
MathComplex<T> pow(MathComplex<T> c, MathComplex<T> p);

template <class T>
T	 		   fabs(MathComplex<T> c);

template <class T>
MathComplex<T> sin(MathComplex<T> c);

template <class T>
MathComplex<T> cos(MathComplex<T> c);

template <class T>
MathComplex<T> tan(MathComplex<T> c);

template <class T>
MathComplex<T> sinh(MathComplex<T> c);

template <class T>
MathComplex<T> cosh(MathComplex<T> c);

template <class T>
MathComplex<T> tanh(MathComplex<T> c);

template <class T>
MathComplex<T> asin(MathComplex<T> c);

template <class T>
MathComplex<T> acos(MathComplex<T> c);

template <class T>
MathComplex<T> atan(MathComplex<T> c);

template <class T>
MathComplex<T> exp(MathComplex<T> c);

template <class T>
MathComplex<T> log(MathComplex<T> c);

template <class T>
MathComplex<T> log10(MathComplex<T> c);

template <class T>
MathComplex<T> sqrt(MathComplex<T> c);

#endif
//EOF
