#include "MathMatrix.h"
#include "MathComplex.h"
#include <stdio.h>
#include <stdlib.h>

template<class T>
MathMatrix<T>::MathMatrix(uint32_t width, uint32_t height)
{
    _width = _height = 0;
    error = false;
    noRecordValue = 0;
    _matrix = 0;
    setSize(width, height);
}

//---------------------------------------------------------

template<class T>
MathMatrix<T>::MathMatrix()
{
    _width = _height = 0;
    error = false;
    noRecordValue = 0;
    _matrix = 0;
}

//---------------------------------------------------------

template<class T>
MathMatrix<T>::MathMatrix(const MathMatrix<T> &matrix)
{
    uint32_t i, j;

    _width = _height = 0;
    error = false;
    noRecordValue = 0;
    _matrix = 0;
    setSize(matrix._width, matrix._height);
    for(i=0; i<_width; i++)
        for(j=0; j<_height; j++)
            _matrix[i][j] = matrix._matrix[i][j];
}

//---------------------------------------------------------

template<class T>
MathMatrix<T>::~MathMatrix()
{
    clear();
}

//---------------------------------------------------------

template<class T>
MathMatrix<T> MathMatrix<T>::operator+(const MathMatrix<T> &matrix)
{
    MathMatrix<T> s(_width, _height);
    uint32_t i, j;
    if(_width != matrix._width || _height != matrix._height)
    {
        s.error = true;
        return s;
    }
    for(i=0; i<_height; i++)
        for(j=0; j<_width; j++)
            s._matrix[j][i] = _matrix[j][i]+matrix._matrix[j][i];
    s._height = _height;
    s._width = _width;
    return s;
}

//---------------------------------------------------------

template<class T>
MathMatrix<T> MathMatrix<T>::operator*(const MathMatrix<T> &matrix)
{
    MathMatrix<T> s(matrix._width, _height);
    uint32_t i, j, k;
    if(_width != matrix._height || _height != matrix._width)
    {
        s.error = true;
        return s;
    }
    s._height = _height;
    s._width = matrix._width;
    for(i=0; i<s._height; i++)
        for(j=0; j<s._width; j++)
        {
            s._matrix[j][i] = 0;
            for(k=0; k<_width; k++) s._matrix[j][i] = s._matrix[j][i] + _matrix[k][i] * matrix._matrix[j][k];
        }
    return s;
}


//---------------------------------------------------------

template<class T>
void MathMatrix<T>::operator=(const MathMatrix<T> &matrix)
{
    uint32_t i, j;
    if(_width != matrix._width || _height != matrix._height) setSize(((MathMatrix<T>)matrix).width(), ((MathMatrix<T>)matrix).height());
    for(i=0; i<_width; i++)
        for(j=0; j<_height; j++)
            _matrix[i][j] = matrix._matrix[i][j];
}

//---------------------------------------------------------

template<class T>
T& MathMatrix<T>::element(uint32_t width, uint32_t height)
{
    if(width<_width && height<_height)
    {
        error = false;
        return _matrix[width][height];
    }
    else
    {
        printf("Error in MathMatrix<T>::element\n");
        exit(-1);
    }
}

//----------------------------------------------------------

template<class T>
void MathMatrix<T>::setElement(uint32_t width, uint32_t height, const T &value)
{
    element(width, height) = value;
}

//----------------------------------------------------------

template<class T>
void MathMatrix<T>::setSize(uint32_t width, uint32_t height)
{
    uint32_t i, j;

    if(_width == width && _height == height)
    {
        for(i=0; i<_width; i++)
            for(j=0; j<_height; j++)
                _matrix[i][j]=0;
        return;
    }
    clear();
    _width = width;
    _height = height;
    if(_width==0 || _height==0)
    {
        _width = _height = 0;
        return;
    }
    _matrix = new T*[_width];
    if(_matrix == 0)
	{
		printf("Error in new\n");
		exit(ERROR_IN_MEMORY);
	}
    for(i=0; i<_width; i++)
	{
        _matrix[i] = new T[_height];
        if(_matrix[i] == 0)
		{
			printf("Error in new\n");
			exit(ERROR_IN_MEMORY);
		}
	}
    for(i=0; i<_width; i++)
        for(j=0; j<_height; j++)
            _matrix[i][j]=0;
}

//----------------------------------------------------------

template<class T>
uint32_t MathMatrix<T>::width()
{
    return _width;
}

//----------------------------------------------------------

template<class T>
uint32_t MathMatrix<T>::height()
{
    return _height;
}

//----------------------------------------------------------

template<class T>
bool MathMatrix<T>::solveNoCopyMatrix(std::vector<T>& roots)
{
    if(_width != _height + 1)
    {
        error = true;
        return false;
    }
    uint32_t i, j;
    T mulDiag, r;
    T tmpValue;

    triangleUpMatrix();
    roots.clear();
    for(i=0; i<_height; i++) roots.push_back((T)0);
    mulDiag = 1;
    for(i=0; i<_height; i++)
        mulDiag = mulDiag * _matrix[i][i];
    if(mulDiag == 0)
    {
        error = true;
        return false;
    }
    for(i=_height; i>0; i--)
    {
        r = _matrix[_height][i-1];
        for(j=_height-1; j>0; j--)
        {
            tmpValue = _matrix[j][i-1] * roots[j];
            r = r - tmpValue;
        }
        roots[i-1] = r / _matrix[i-1][i-1];
    }
    error = false;
    return true;
}//solveNoCopyMatrix

//-------------------------------------------------------------------------

template<class T>
bool MathMatrix<T>::solve(std::vector<T>& roots)
{
    MathMatrix<T> tmpMatrix = *this;

    return tmpMatrix.solveNoCopyMatrix(roots);
}//solve

//-------------------------------------------------------------------------

template<class T>
int32_t MathMatrix<T>::triangleUpMatrix()
{
    uint32_t i, j, k, l;
    uint32_t countStringDown;   //счетчик строк, начиная снизу
    int32_t sign = 1;
    T mulString;
    T tmpValue;

    if(!_width || !_height) return sign;

    for(j=0; j<_height - 1; j++)       //j-я строка, которая вычитается из других строк
    {
        countStringDown = _height - 1;
        if(_matrix[j][j] == 0)          /* замена строк местами, если главный элемент равен 0 */
        {
            while ((_matrix[j][countStringDown] == 0) && (countStringDown != j)) countStringDown--;
            if(countStringDown != j)
            {
                sign = -sign;
                for(k=0; k<_width; k++)
                {
                    tmpValue = _matrix[k][countStringDown];
                    _matrix[k][countStringDown] = _matrix[k][j];
                    _matrix[k][j] = tmpValue;
                }
            }
        }
        if(countStringDown != j)
            for(i=j+1; i<_height; i++) //i-строка из которой вычитается строка j*mm
            {
                mulString = _matrix[j][i] / _matrix[j][j];
                _matrix[j][i] = (T)0;
                for(l=j + 1; l<_width; l++)
                    _matrix[l][i] = _matrix[l][i] - _matrix[l][j] * mulString;
            }
    }
    error = false;
    return sign;
}

//-------------------------------------------------------------------------

template<class T>
T MathMatrix<T>::determinant()
{
    MathMatrix<T> tmpMatrix = *this;

    return tmpMatrix.determinantNoCopyMatrix();
}//determinant

//--------------------------------------------------------------------------

template<class T>
T MathMatrix<T>::determinantNoCopyMatrix()
{
    uint32_t i;
    T det, sign = (T)1;

    if(_width != _height)
    {
        error = true;
        return (T)0;
    }
    sign = triangleUpMatrix();
    det = 1;
    for(i=0; i<_height; i++)
        det = det * _matrix[i][i];
    det = det * sign;
    error = false;
    return det;
}//determinantNoCopyMatrix

//----------------------------------------------------------------------

template<class T>
MathMatrix<T> MathMatrix<T>::eraseIJ(uint32_t width, uint32_t height)
{
    uint32_t i, j, i2=0, j2=0;
    MathMatrix<T> ret;

    if(_width <= 1 || _height <= 1) return ret;
    ret.setSize(_width - 1, _height - 1);

    for(i=0; i<ret._height; i++)
    {
        if(i2 == height) i2++;
        j2 = 0;
        for(j=0; j<ret._width; j++)
        {
            if(j2 == width) j2++;
            ret._matrix[j][i] = _matrix[j2][i2];
            j2++;
        }
        i2++;
    }

    return ret;
}//eraseIJ

//----------------------------------------------------------------

template<class T>
T MathMatrix<T>::algebraicAddition(uint32_t width, uint32_t height)
{
    T sign = (T)-1;
    uint32_t i;

    for(i=0; i<width + height + 1; i++) sign = -sign;
    return eraseIJ(width, height).determinantNoCopyMatrix()*sign;
}

//----------------------------------------------------------------

template<class T>
MathMatrix<T> MathMatrix<T>::inverseMatrix()
{
    MathMatrix<T> ret;
    T det;
    uint32_t i, j;

    if(_width != _height)
    {
        error = true;
        return ret;
    }
    ret.setSize(_width, _height);
    det = determinant();
    for(i=0; i<_height; i++)
        for(j=0; j<_width; j++)
            ret.element(j, i) = algebraicAddition(i, j)/det;
    return ret;
}

//----------------------------------------------------------------

template<class T>
MathMatrix<T> MathMatrix<T>::transp()
{
    uint32_t i, j;
    MathMatrix<T> ret;

    ret.setSize(height(), width());
    for(i=0; i<ret.height(); i++)
        for(j=0; j<ret.width(); j++)
            ret.element(j, i) = element(i, j);
    return ret;
}

//----------------------------------------------------------------

template<class T>
void MathMatrix<T>::clear()
{
    uint32_t i;
    if(_matrix != 0)
    {
        for(i=0; i<_width; i++)
            delete[] _matrix[i];
        delete[] _matrix;
    }
    noRecordValue = 0;
    _width = _height = 0;
    _matrix = 0;
    error = false;
}

//----------------------------------------------------------------

template<class T>
bool MathMatrix<T>::getError()
{
    return error;
}


//template class MathMatrix<int>;
template class MathMatrix<float>;
template class MathMatrix<double>;
template class MathMatrix<MathComplex<float> >;
template class MathMatrix<MathComplex<double> >;
//EOF
