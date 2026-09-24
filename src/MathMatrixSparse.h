#ifndef MATHMATRIXSPARSE
#define MATHMATRIXSPARSE

/****************************************************************************************
 *     Класс для работы с разреженными матрицами. Типы: float, double, MathComplex      *
 ****************************************************************************************/

#include <math.h>
#include <stdint.h>
#include <vector>
#include <list>
#include <map>
#include "MathComplex.h"

#define ERROR_IN_MEMORY     (-1)

template<class T>
class MathMatrixSparse
{
public:
    MathMatrixSparse(uint32_t width, uint32_t height);
    MathMatrixSparse();
    MathMatrixSparse(const MathMatrixSparse<T> &matrix);
    ~MathMatrixSparse();
    void                    operator=(const MathMatrixSparse<T> &matrix);
    T                       element(uint32_t width, uint32_t height);          /* считать элемент матрицы */
    void                    setElement(uint32_t width, uint32_t height, const T& value);          /* записать элемент матрицы */
    void                    setElementExt(uint32_t width, uint32_t height, const T& value);       /* записать элемент матрицы, при этом матрица при необходимости увеличивается */
    void                    setSize(uint32_t width, uint32_t height);
    bool                    getCSC(std::vector<T>& values, std::vector<uint32_t>& rows, std::vector<uint32_t>& columns);
    uint32_t                width();
    uint32_t                height();
    int32_t                 triangleUpMatrix();                           /* приводит матрицу к треугольному виду */
    void                    eraseIJ(uint32_t width, uint32_t height);   /* вырезать столбец width, строку height */
    bool                    solveNoCopyMatrix(std::vector<T> &roots);   /* найти корни линейного уравнения, текущая матрица приводится к треугольной */
    bool                    solve(std::vector<T> &roots);               /* найти корни линейного уравнения, текущая матрица сохраняется */
    T                       determinant();                              /* вычислить определитель, текущая матрица сохраняется */
    T                       determinantNoCopyMatrix();                  /* вычислить определитель, текущая матрица приводится к треугольной */
    T                       algebraicAddition(uint32_t width, uint32_t height); /* алгебраическое дополнение */
    MathMatrixSparse<T>     inverseMatrix();
    MathMatrixSparse<T>     transp();                                   /* транспонирование */
    void                    extensionWidth(uint32_t deltaWidth);
    void                    extensionHeight(uint32_t deltaHeight);
    void                    extensionMatrix(uint32_t deltaWidth, uint32_t deltaHeight); /* увеличение размера матрицы */
    void                    clear();
    bool                    getError();                                 /* ошибки чтения/записи элемента */
private:
    void                    addString(uint32_t srcNum, uint32_t destNum, const T& mul);
    std::vector<std::map<uint32_t, T> > _matrix;
    uint32_t _width;
    bool error;
};

#endif // MATHMATRIXSPARSE

