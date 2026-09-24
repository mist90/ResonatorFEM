#ifndef MATHMATRIX_H
#define MATHMATRIX_H

/***************************************************************************
 *     Класс для работы с матрицами. Типы: float, double, MathComplex      *
 **************************************************************************/

#include <math.h>
#include <stdint.h>
#include <vector>
#include "MathComplex.h"

#define ERROR_IN_MEMORY     (-1)

template<class T>
class MathMatrix
{
public:
    MathMatrix(uint32_t width, uint32_t height);
    MathMatrix();
    MathMatrix(const MathMatrix<T> &matrix);
    ~MathMatrix();
    MathMatrix<T>           operator+(const MathMatrix<T> &matrix);
    MathMatrix<T>           operator*(const MathMatrix<T> &matrix);
    void                    operator=(const MathMatrix<T> &matrix);
    T&                      element(uint32_t width, uint32_t height);          /* считать/записать элемент матрицы */
    void                    setElement(uint32_t width, uint32_t height, const T& value);
    void                    setSize(uint32_t width, uint32_t height);
    uint32_t                width();
    uint32_t                height();
    int32_t                 triangleUpMatrix();                           /* приводит матрицу к треугольному виду */
    T                       determinant();                              /* вычислить определитель, текущая матрица сохраняется */
    T                       determinantNoCopyMatrix();                  /* вычислить определитель, текущая матрица приводится к треугольной */
    MathMatrix<T>           eraseIJ(uint32_t width, uint32_t height);   /* вырезать столбец width, строку height */
    T                       algebraicAddition(uint32_t width, uint32_t height); /* алгебраическое дополнение */
    MathMatrix<T>           inverseMatrix();
    bool                    solveNoCopyMatrix(std::vector<T> &roots);   /* найти корни линейного уравнения, текущая матрица приводится к треугольной */
    bool                    solve(std::vector<T> &roots);               /* найти корни линейного уравнения, текущая матрица сохраняется */
    MathMatrix<T>           transp();                                   /* транспонирование */
    void                    clear();
    bool                    getError();                                 /* ошибки чтения/записи элемента */
private:
    T **_matrix;
    T noRecordValue;
    uint32_t _width, _height;
    bool error;
};


#endif
//EOF
