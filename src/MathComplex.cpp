#include "MathComplex.h"


template<class T>
MathComplex<T>::MathComplex(T re,T im)
{
    _re=re;
    _im=im;
}

//--------------------------------------------------------

template<class T>
MathComplex<T>::MathComplex()
{
    _re=0;
    _im=0;
}

//--------------------------------------------------------

template<class T>
void MathComplex<T>::setExponentialArg(T amplitude, T phase)
{
    _re=amplitude*cos(phase);
    _im=amplitude*sin(phase);
}

//--------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::operator+(const MathComplex<T> &c1)
{
    MathComplex<T> s;
    s._re=c1._re+_re;
    s._im=c1._im+_im;
    return s;
}

//--------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::operator+(const T &c1)
{
    MathComplex<T> s;
    s._re=c1+_re;
    s._im=_im;
    return s;
}

//--------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::operator-(const MathComplex<T> &c1)
{
    MathComplex<T> s;
    s._re=_re-c1._re;
    s._im=_im-c1._im;
    return s;
}

//--------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::operator-(const T &c1)
{
    MathComplex<T> s;
    s._re=_re-c1;
    s._im=_im;
    return s;
}

//---------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::operator-()
{
    MathComplex<T> s;
    s._re=-_re;
    s._im=-_im;
    return s;
}

//---------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::operator*(const MathComplex<T> &c1)
{
    MathComplex<T> s;
    s._re=_re*c1._re-_im*c1._im;
    s._im=_re*c1._im+_im*c1._re;
    return s;
}

//---------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::operator*(const T &c1)
{
    MathComplex<T> s;
    s._re=_re*c1;
    s._im=_im*c1;
    return s;
}

//----------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::operator/(const MathComplex<T> &c1)
{
    MathComplex<T> s;
    s._re = (_re * c1._re +_im * c1._im) / (c1._re * c1._re + c1._im * c1._im);
    s._im = (c1._re * _im -_re * c1._im) / (c1._re * c1._re + c1._im * c1._im);
    return s;
}

//----------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::operator/(const T &c1)
{
    MathComplex<T> s;
    s._re = _re/c1;
    s._im = _im/c1;
    return s;
}

//----------------------------------------------------------

template<class T>
void MathComplex<T>::operator=(const MathComplex<T> &c1)
{
    _re = c1._re;
    _im = c1._im;
}

//----------------------------------------------------------

template<class T>
void MathComplex<T>::operator=(const T &c1)
{
    _re = c1;
    _im = 0;
}

//----------------------------------------------------------

template<class T>
void MathComplex<T>::operator+=(const MathComplex<T> &c1)
{
    _re = _re + c1._re;
    _im = _im + c1._im;
}

//----------------------------------------------------------

template<class T>
void MathComplex<T>::operator-=(const MathComplex<T> &c1)
{
    _re = _re - c1._re;
    _im = _im - c1._im;
}

//----------------------------------------------------------

template<class T>
void MathComplex<T>::operator*=(const MathComplex<T> &c1)
{
	T tmpRe = _re;
    _re=_re*c1._re-_im*c1._im;
    _im= tmpRe*c1._im+_im*c1._re;
}

//----------------------------------------------------------

template<class T>
void MathComplex<T>::operator/=(const MathComplex<T> &c1)
{
	T tmpRe = _re;
    _re = (_re * c1._re +_im * c1._im) / (c1._re * c1._re + c1._im * c1._im);
    _im = (c1._re * _im - tmpRe * c1._im) / (c1._re * c1._re + c1._im * c1._im);
}

//----------------------------------------------------------

template<class T>
bool MathComplex<T>::operator==(const MathComplex<T> &c1)
{
    if(_re == c1._re && _im == c1._im)
        return true;
    else
        return false;
}

//----------------------------------------------------------

template<class T>
bool MathComplex<T>::operator!=(const MathComplex<T> &c1)
{
    if(_re != c1._re || _im != c1._im)
        return true;
    else
        return false;
}

//----------------------------------------------------------

template<class T>
bool MathComplex<T>::operator==(const T &c1)
{
    if(_re == c1 && _im == c1)
        return true;
    else
        return false;
}

//-----------------------------------------------------------

template<class T>
bool MathComplex<T>::operator!=(const T &c1)
{
    if(_re != c1 && _im != c1)
        return true;
    else
        return false;
}

//-----------------------------------------------------------

template<class T>
T& MathComplex<T>::re()
{
    return _re;
}

//-----------------------------------------------------------

template<class T>
T& MathComplex<T>::im()
{
    return _im;
}

//-----------------------------------------------------------

template<class T>
T MathComplex<T>::phase()
{
    double s;
    float sf;
    if(_re == 0 && _im == 0) return 0;
    if(_re == 0 && _im > 0) return M_PI/2;
    if(_re == 0 && _im < 0) return -M_PI/2;
    if(sizeof(T) == sizeof(double))
    {
        s = fabs(atan(_im/_re));
        if(_re > 0)
        {
            if(_im > 0)
                return s;
            else
                return -s;
        }
        else
        {
            if(_im > 0)
                return -s + M_PI;
            else
                return s - M_PI;
        }
    }
    else
    {
        sf = fabsf(atanf(_im/_re));
        if(_re > 0)
        {
            if(_im > 0)
                return sf;
            else
                return -sf;
        }
        else
        {
            if(_im > 0)
                return -sf + M_PI;
            else
                return sf - M_PI;
        }
    }
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> MathComplex<T>::powc(const MathComplex<T> &p)
{
    MathComplex<T> s = MathComplex<T>();
    T A;
    if(sizeof(T) == sizeof(double))
    {
        A = sqrt(_re *_re + _im *_im);
        s.re() = cos(phase() * p._re + p._im * (phase() + log(A))) * pow(A, p._re);
        s.im() = sin(phase() * p._re + p._im * (phase() + log(A))) * pow(A, p._re);
    }
    else
    {
        A = sqrtf(_re *_re +_im *_im);
        s.re() = cosf(phase() * p._re + p._im * (phase() + logf(A))) * powf(A, p._re);
        s.im() = sinf(phase() * p._re + p._im * (phase() + logf(A))) * powf(A, p._re);
    }
    return s;
}

//-----------------------------------------------------------

template <class T>
T MathComplex<T>::absc()
{
    if(sizeof(double) == sizeof(T))
    {
        return sqrt(_re *_re + _im *_im);
    }
    else
    {
        return sqrtf(_re *_re + _im *_im);
    }
}

//-----------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::sinc()
{
	MathComplex<T> s = MathComplex<T>();
	if(sizeof(double) == sizeof(T))
	{
		s.re() = sin(_re) * (exp(-_im) + exp(_im))/2.0;
		s.im() = cos(_re) * (exp(_im) - exp(-_im))/2.0;
	}
	else
	{
		s.re() = sinf(_re) * (expf(-_im) + expf(_im))/2.0f;
		s.im() = cosf(_re) * (expf(_im) - expf(-_im))/2.0f;
	}
	return s;
}

//------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::cosc()
{
	MathComplex<T> s = MathComplex<T>();
	if(sizeof(double) == sizeof(T))
	{
		s.re() = cos(_re) * (exp(-_im) + exp(_im))/2.0;
		s.im() = sin(_re) * (exp(-_im) - exp(_im))/2.0;
	}
	else
	{
		s.re() = cosf(_re) * (expf(-_im) + expf(_im))/2.0f;
		s.im() = sinf(_re) * (expf(-_im) - expf(_im))/2.0f;
	}
	return s;
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::tanc()
{
	return sinc()/cosc();
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::sinhc()
{
	return ((*this).expc() - (*this * (-1)).expc()) / 2;
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::coshc()
{
	return ((*this).expc() + (*this * (-1)).expc()) / 2;
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::tanhc()
{
	return sinhc() / coshc();
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::asinc()
{
	MathComplex<T> s = *this;
	s = - (s * s) + 1;
	s = s.sqrtc();
	s = MathComplex<T>(0, 1) * (*this) + s;
	return MathComplex<T>(0, -1) * s.logc();
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::acosc()
{
	MathComplex<T> s = *this;
	s = - (s * s) + 1;
	s = s.sqrtc();
	s = MathComplex<T>(0, 1) * (*this) + s;
	return MathComplex<T>(0, 1) * s.logc() + M_PI/2;
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::atanc()
{
	MathComplex<T> s1 = MathComplex<T>(0, 1) * (*this) * (-1) + 1;
	MathComplex<T> s2 = MathComplex<T>(0, 1) * (*this) + 1;
	return MathComplex<T>(0, 1) * (s1.logc() - s2.logc()) / 2;
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::expc()
{
	MathComplex<T> s = MathComplex<T>();
	if(sizeof(double) == sizeof(T))
	{
		s.re() = exp(_re) * cos(_im);
		s.im() = exp(_re) * sin(_im);
	}
	else
	{
		s.re() = expf(_re) * cosf(_im);
		s.im() = expf(_re) * sinf(_im);
	}
	return s;
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::logc()
{
	MathComplex<T> s = MathComplex<T>();
	T A;
    if(sizeof(T) == sizeof(double))
    {
        A = sqrt(_re *_re + _im *_im);
        s.re() = log(A);
        s.im() = phase();
    }
    else
    {
        A = sqrtf(_re *_re +_im *_im);
        s.re() = log(A);
        s.im() = phase();
    }
    return s;
}

//-------------------------------------------------------------

template<class T>
MathComplex<T> MathComplex<T>::log10c()
{
	if(sizeof(T) == sizeof(double))
		return logc() / log(10.0);
	else
		return logc() / logf(10.0f);
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> MathComplex<T>::sqrtc()
{
    MathComplex<T> s = MathComplex<T>();
    T A;
    if(sizeof(T) == sizeof(double))
    {
        A = sqrt(_re *_re + _im *_im);
        s.re() = cos(phase() * 0.5) * sqrt(A);
        s.im() = sin(phase() * 0.5) * sqrt(A);
    }
    else
    {
        A = sqrtf(_re *_re +_im *_im);
        s.re() = cosf(phase() * 0.5f) * sqrtf(A);
        s.im() = sinf(phase() * 0.5f) * sqrtf(A);
    }
    return s;
}

/************************************************************************
 *                Перегрузка функций math.h                             *
 * **********************************************************************/

//-----------------------------------------------------------

template <class T>
MathComplex<T> pow(MathComplex<T> c, MathComplex<T> p)
{
	return c.powc(p);
}

//-----------------------------------------------------------

template <class T>
T fabs(MathComplex<T> c)
{
	return c.absc();
}

//----------------------------------------------------------

template <class T>
MathComplex<T> sin(MathComplex<T> c)
{
	return c.sinc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> cos(MathComplex<T> c)
{
	return c.cosc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> tan(MathComplex<T> c)
{
	return c.tanc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> sinh(MathComplex<T> c)
{
	return c.sinhc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> cosh(MathComplex<T> c)
{
	return c.coshc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> tanh(MathComplex<T> c)
{
	return c.tanhc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> asin(MathComplex<T> c)
{
	return c.asinc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> acos(MathComplex<T> c)
{
	return c.acosc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> atan(MathComplex<T> c)
{
	return c.atanc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> exp(MathComplex<T> c)
{
	return c.expc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> log(MathComplex<T> c)
{
	return c.logc();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> log10(MathComplex<T> c)
{
	return c.log10c();
}

//-----------------------------------------------------------

template <class T>
MathComplex<T> sqrt(MathComplex<T> c)
{
	return c.sqrtc();
}

//-----------------------------------------------------------



template class MathComplex<float>;
template class MathComplex<double>;

template MathComplex<double> pow<double>(MathComplex<double> c, MathComplex<double> p);
template double fabs(MathComplex<double> c);
template MathComplex<double> sin<double>(MathComplex<double> c);
template MathComplex<double> cos<double>(MathComplex<double> c);
template MathComplex<double> tan<double>(MathComplex<double> c);
template MathComplex<double> sinh<double>(MathComplex<double> c);
template MathComplex<double> cosh<double>(MathComplex<double> c);
template MathComplex<double> tanh<double>(MathComplex<double> c);
template MathComplex<double> asin<double>(MathComplex<double> c);
template MathComplex<double> acos<double>(MathComplex<double> c);
template MathComplex<double> atan<double>(MathComplex<double> c);
template MathComplex<double> exp<double>(MathComplex<double> c);
template MathComplex<double> log<double>(MathComplex<double> c);
template MathComplex<double> log10<double>(MathComplex<double> c);
template MathComplex<double> sqrt<double>(MathComplex<double> c);

template MathComplex<float> pow<float>(MathComplex<float> c, MathComplex<float> p);
template float fabs<float>(MathComplex<float> c);
template MathComplex<float> sin<float>(MathComplex<float> c);
template MathComplex<float> cos<float>(MathComplex<float> c);
template MathComplex<float> tan<float>(MathComplex<float> c);
template MathComplex<float> sinh<float>(MathComplex<float> c);
template MathComplex<float> cosh<float>(MathComplex<float> c);
template MathComplex<float> tanh<float>(MathComplex<float> c);
template MathComplex<float> asin<float>(MathComplex<float> c);
template MathComplex<float> acos<float>(MathComplex<float> c);
template MathComplex<float> atan<float>(MathComplex<float> c);
template MathComplex<float> exp<float>(MathComplex<float> c);
template MathComplex<float> log<float>(MathComplex<float> c);
template MathComplex<float> log10<float>(MathComplex<float> c);
template MathComplex<float> sqrt<float>(MathComplex<float> c);

//EOF
