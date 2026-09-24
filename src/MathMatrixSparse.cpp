#include "MathMatrixSparse.h"
#include <stdlib.h>
#include <stdio.h>

template<class T>
MathMatrixSparse<T>::MathMatrixSparse(uint32_t width, uint32_t height)
{
    _width = 0;
    error = false;
    setSize(width, height);
}


template<class T>
MathMatrixSparse<T>::MathMatrixSparse()
{
    _width = 0;
    error = false;
}


template<class T>
MathMatrixSparse<T>::MathMatrixSparse(const MathMatrixSparse<T> &matrix)
{
    _width = matrix._width;
    error = matrix.error;
    _matrix = matrix._matrix;
}

template<class T>
MathMatrixSparse<T>::~MathMatrixSparse()
{
    _matrix.clear();
}

template<class T>
void MathMatrixSparse<T>::operator=(const MathMatrixSparse<T> &matrix)
{
    _width = matrix._width;
    error = matrix.error;
    _matrix = matrix._matrix;
}


template<class T>
T MathMatrixSparse<T>::element(uint32_t width, uint32_t height)
{
    typename std::map<uint32_t,T>::iterator itCurrent;

    if(width >= _width || height >= _matrix.size())
    {
        printf("Error in MathMatrixSparse<T>::element\n");
        exit(-1);
    }
    error = false;
    itCurrent = _matrix[height].find(width);
    if(itCurrent == _matrix[height].end()) return 0;
    else return itCurrent->second;
}

template<class T>
void MathMatrixSparse<T>::setElement(uint32_t width, uint32_t height, const T &value)
{
    typename std::map<uint32_t,T>::iterator itCurrent;

    if(width >= _width || height >= _matrix.size())
    {
        error = true;
        return;
    }
    error = false;
    itCurrent = _matrix[height].find(width);
    if(itCurrent == _matrix[height].end() && ((T)(value)) == 0) return;   /* не нужно добавлять элемент, если он нулевой */
    if(itCurrent == _matrix[height].end()) _matrix[height].insert(std::pair<uint32_t,T>(width, value));
    else
    {
        if(((T)(value)) == 0) _matrix[height].erase(itCurrent);
        else itCurrent->second = value;
    }
}

template<class T>
void MathMatrixSparse<T>::setElementExt(uint32_t width, uint32_t height, const T &value)
{
    uint32_t delta = 0, delta2 = 0;
    if(width >= _width)
        delta = width - _width + 1;
    if(height >= _matrix.size())
        delta2 = height - _matrix.size() + 1;
    if(delta2 > delta) delta = delta2;
    if(delta) extensionMatrix(delta, delta);
    setElement(width, height, value);
}


template<class T>
void MathMatrixSparse<T>::setSize(uint32_t width, uint32_t height)
{
    typename std::map<uint32_t,T>::iterator itCurrent;
    typename std::vector<std::map<uint32_t,T> >::iterator itHeight = _matrix.begin();
    uint32_t i;

    if(_width == width && _matrix.size() == height)
    {
        for(itHeight = _matrix.begin(); itHeight != _matrix.end(); itHeight++)
            for(itCurrent = itHeight->begin(); itCurrent != itHeight->end(); itCurrent++)
                itCurrent->second = 0;
        return;
    }
    clear();
    _width = width;
    error = false;
    for(i=0; i<height; i++) _matrix.push_back(std::map<uint32_t, T>());
}


template<class T>
bool MathMatrixSparse<T>::getCSC(std::vector<T> &values, std::vector<uint32_t> &rows, std::vector<uint32_t> &columns)
{
    uint32_t i, j;
    T curElement;
    bool beginColumn;

    if(!_width) return false;
    values.clear();
    rows.clear();
    columns.clear();
    for(i=0; i<_width; i++)
    {
        beginColumn = true;
        for(j=0; j<_matrix.size(); j++)
        {
            curElement = element(i, j);
            if(curElement != 0)
            {
                values.push_back(curElement);
                rows.push_back(j);
                if(beginColumn)
                {
                    columns.push_back(values.size() - 1);
                    beginColumn = false;
                }
            }

        }
        if(beginColumn)
        {
            columns.push_back(values.size() - 1);
            beginColumn = false;
        }
    }
    columns.push_back(values.size());
    return true;
}

template<class T>
uint32_t MathMatrixSparse<T>::width()
{
    return _width;
}


template<class T>
uint32_t MathMatrixSparse<T>::height()
{
    return _matrix.size();
}


template<class T>
int32_t MathMatrixSparse<T>::triangleUpMatrix()
{
    uint32_t i, j, k;
    uint32_t countStringDown;   //счетчик строк, начиная снизу
    int32_t sign = 1;
    T mulString;
    T tmpValue;

    if(!_width || !_matrix.size()) return sign;

    for(j=0; j<_matrix.size() - 1; j++)       //j-я строка, которая вычитается из других строк
    {
        countStringDown = _matrix.size() - 1;
        if(element(j, j) == 0)          /* замена строк местами, если главный элемент равен 0 */
        {
            while ((element(j, countStringDown) == 0) && (countStringDown != j)) countStringDown--;
            if(countStringDown != j)
            {
                sign = -sign;
                for(k=0; k<_width; k++)
                {
                    tmpValue = element(k, countStringDown);
                    setElement(k, countStringDown, element(k, j));
                    setElement(k, j, tmpValue);
                }
            }
        }
        if(countStringDown != j)
            for(i=j+1; i<_matrix.size(); i++) //i-строка из которой вычитается строка j*mulString
            {
                mulString = element(j, i) / element(j, j);
                if(mulString != 0) addString(j, i, -mulString);
                setElement(j, i,(T)0);
            }
    }
    error = false;
    return sign;
}


template<class T>
void MathMatrixSparse<T>::eraseIJ(uint32_t width, uint32_t height)
{
    typename std::map<uint32_t,T>::iterator itCurrent, itDelete;
    typename std::vector<std::map<uint32_t,T> >::iterator itHeight;
    T tmpValue;
    uint32_t tmpFirst;

    if(width >= _width || height >= _matrix.size() || _width <= 1 || _matrix.size() <= 1)
    {
        error = true;
        return;
    }
    error = false;
    /* удаление столбца */
    for(itHeight = _matrix.begin(); itHeight != _matrix.end(); itHeight++)
    {
        itCurrent = itHeight->find(width);
        if(itCurrent != itHeight->end()) itHeight->erase(itCurrent);
        itCurrent = itHeight->begin();
        while(itCurrent != itHeight->end())
        {
            if(itCurrent->first > width)
            {
                tmpValue = itCurrent->second;
                tmpFirst = itCurrent->first;
                itDelete = itCurrent;
                itCurrent++;
                itHeight->erase(itDelete);
                itHeight->insert(std::pair<uint32_t,T>(tmpFirst - 1, tmpValue));
            }
            else itCurrent++;
        }

    }
    /* удаление строки */
    itHeight = _matrix.begin();
    std::advance(itHeight, height);
    _matrix.erase(itHeight);
    _width--;
}

template<class T>
bool MathMatrixSparse<T>::solveNoCopyMatrix(std::vector<T> &roots)
{
    if(_width != _matrix.size() + 1)
    {
        error = true;
        return false;
    }
    uint32_t i, j;
    T mulDiag, r;
    T tmpValue;

    triangleUpMatrix();
    roots.clear();
    mulDiag = 1;
    for(i=0; i<_matrix.size(); i++)
        mulDiag = mulDiag * element(i, i);
    if(mulDiag == 0)
    {
        error = true;
        return false;
    }
    for(i=0; i<_matrix.size(); i++) roots.push_back(0);
    for(i=_matrix.size(); i>0; i--)
    {
        r = element(_matrix.size(), i-1);
        for(j=_matrix.size()-1; j>0; j--)
        {
            tmpValue = element(j, i-1) * roots[j];
            r = r - tmpValue;
        }
        roots[i-1] = r / element(i-1, i-1);
    }
    error = false;
    return true;
}


template<class T>
bool MathMatrixSparse<T>::solve(std::vector<T> &roots)
{
    MathMatrixSparse<T> tmpMatrix = *this;

    return tmpMatrix.solveNoCopyMatrix(roots);
}


template<class T>
T MathMatrixSparse<T>::determinant()
{
    MathMatrixSparse<T> tmpMatrix = *this;

    return tmpMatrix.determinantNoCopyMatrix();
}//determinant


template<class T>
T MathMatrixSparse<T>::determinantNoCopyMatrix()
{
    uint32_t i;
    T det, sign = 1;

    if(_width != _matrix.size())
    {
        error = true;
        return 0;
    }
    sign = triangleUpMatrix();
    det = 1;
    for(i=0; i<_matrix.size(); i++)
        det = det * element(i, i);
    det = det * sign;
    error = false;
    return det;
}//determinantNoCopyMatrix


template<class T>
T MathMatrixSparse<T>::algebraicAddition(uint32_t width, uint32_t height)
{
    T sign = -1;
    uint32_t i;
    MathMatrixSparse matrix = *this;

    for(i=0; i<width + height + 1; i++) sign = -sign;
    matrix.eraseIJ(width, height);
    return matrix.determinantNoCopyMatrix()*sign;
}


template<class T>
MathMatrixSparse<T> MathMatrixSparse<T>::inverseMatrix()
{
    MathMatrixSparse<T> ret;
    T det;
    uint32_t i, j;

    if(_width != _matrix.size())
    {
        error = true;
        return ret;
    }
    ret.setSize(_width, _matrix.size());
    det = determinant();
    for(i=0; i<_matrix.size(); i++)
        for(j=0; j<_width; j++)
            ret.setElement(j, i, algebraicAddition(i, j)/det);
    return ret;
}

template<class T>
MathMatrixSparse<T> MathMatrixSparse<T>::transp()
{
    uint32_t i, j;
    MathMatrixSparse<T> ret;

    ret.setSize(height(), width());
    for(i=0; i<ret.height(); i++)
        for(j=0; j<ret.width(); j++)
            ret.setElement(j, i, element(i, j));
    return ret;
}

template<class T>
void MathMatrixSparse<T>::extensionWidth(uint32_t deltaWidth)
{
    _width += deltaWidth;
}

template<class T>
void MathMatrixSparse<T>::extensionHeight(uint32_t deltaHeight)
{
    uint32_t i;

    for(i=0; i<deltaHeight; i++) _matrix.push_back(std::map<uint32_t, T>());
}

template<class T>
void MathMatrixSparse<T>::extensionMatrix(uint32_t deltaWidth, uint32_t deltaHeight)
{
    extensionWidth(deltaWidth);
    extensionHeight(deltaHeight);
}

template<class T>
void MathMatrixSparse<T>::clear()
{
    _matrix.clear();
    _width = 0;
    error = false;
}


template<class T>
bool MathMatrixSparse<T>::getError()
{
    return error;
}

template<class T>
void MathMatrixSparse<T>::addString(uint32_t srcNum, uint32_t destNum, const T &mul)
{
    typename std::map<uint32_t,T>::iterator itSrcCurrent, itDestCurrent;

    for(itSrcCurrent = _matrix[srcNum].begin(); itSrcCurrent != _matrix[srcNum].end(); itSrcCurrent++)
    {
        itDestCurrent = _matrix[destNum].find(itSrcCurrent->first);
        if(itDestCurrent == _matrix[destNum].end())
        {
            _matrix[destNum].insert(std::pair<uint32_t,T>(itSrcCurrent->first, itSrcCurrent->second*mul));
        }
        else
        {
            itDestCurrent->second += itSrcCurrent->second*mul;
            if(itDestCurrent->second == 0)
            {
                _matrix[destNum].erase(itDestCurrent);
            }
        }
    }
}


template class MathMatrixSparse<float>;
template class MathMatrixSparse<double>;
template class MathMatrixSparse<MathComplex<float> >;
template class MathMatrixSparse<MathComplex<double> >;


//EOF
