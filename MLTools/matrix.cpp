#include "stdafx.h"
#include "matrix.h"

template<class T>
Matrix<T>::Matrix(int rows, int columns, T defaultValue)
    :m_rows(rows), m_columns(columns)
{
     m_values = new T[m_rows *  m_columns];

     for(int i= 0 ;i < m_rows ;i++)
     {
		 for (int j = 0; j < m_columns; j++)
		 {
             m_values[i * m_columns + j] = defaultValue;
		 }
     }
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& matrix)
    :m_rows(matrix.m_rows), m_columns(matrix.m_columns)
{
     m_values = new T[m_rows *  m_columns];

     for(int i=0 ;i< m_rows ;i++)
     {
         for(int j= 0 ; j < m_columns ; j++)
         {
             m_values[i * m_columns + j] = matrix.m_values[i * m_columns + j];
         }
     }
}

template<class T>
Matrix<T>::~Matrix()
{
    delete[] m_values;
}

template<class T>
int Matrix<T>::rows() const
{
    return m_rows;
}

template<class T>
int Matrix<T>::columns() const
{
    return m_columns;
}

template<class T>
bool Matrix<T>::isRowVector() const
{
	if (m_rows == 1)
		return true;

	return false;
}

template<class T>
bool Matrix<T>::isColumnVector() const
{
	if (m_columns == 1)
		return true;

	return false;
}

template<class T>
bool Matrix<T>::isSquareMatrix() const
{
	if (m_columns == m_rows)
		return true;

	return false;
}

template<class T>
Matrix<T> Matrix<T>::transpose() const
{
	Matrix<T> matrix(m_columns, m_rows);

	for (int i = 0; i< m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
            matrix.m_values[j * m_rows + i] = m_values[i*m_columns + j];
		}
	}

	return matrix;
}

template<class T>
Matrix<T> Matrix<T>::identity() const
{
	assert(isSquareMatrix());

	Matrix<T> matrix(m_rows, m_columns);

	for (int i = 0; i< m_rows; i++)
	{
        matrix.m_values[i *  m_columns + i] = (T)1;
	}

	return matrix;
}

template<class T>
Matrix<T> Matrix<T>::getRow(int i) const
{
	Matrix<T> matrix(1, m_columns);

	for (int j = 0; j < m_columns; j++)
	{
        matrix.m_values[j] = m_values[i*m_columns+j];
	}

	return matrix;
}

template<class T>
void Matrix<T>::setRow(const Matrix<T>& row, int i)
{
    assert(row.isRowVector() && row.m_columns == m_columns);

	for (int j = 0; j < m_columns; j++)
	{
        m_values[i*m_columns + j] = row.m_values[j];
	}

}

template<class T>
void Matrix<T>::appendAsRows(const Matrix<T>& rows)
{
	assert(m_columns = rows.m_columns);
}

template<class T>
Matrix<T> Matrix<T>::getColumn(int j) const
{
	Matrix<T> matrix(m_rows, 1);

	for (int i = 0; i < m_rows; i++)
	{
        matrix.m_values[i] = m_values[i*m_columns+j];
	}

	return matrix;
}

template<class T>
void Matrix<T>::setColumn(const Matrix<T>& column , int j)
{
    assert(column.isColumnVector() && column.m_rows == m_rows);

	for (int i = 0; i < m_rows; i++)
	{
        m_values[i*m_columns+j] = column.m_values[i];
	}
}

template<class T>
void Matrix<T>::appendAsColumns(const Matrix<T>& columns)
{

}

template<class T>
T Matrix<T>::operator ()(int i, int j) const
{
    return m_values[i*m_columns+j];
}

template<class T>
T& Matrix<T>::operator ()(int i , int j)
{
    return m_values[i*m_columns + j];
}

template<class T>
Matrix<T>& Matrix<T>::operator =(const Matrix<T>& matrix)
{
	m_rows = matrix.m_rows;
	m_columns = matrix.m_columns;

	if (this != &matrix)
	{
        m_values = new T[m_rows*m_columns];

		for (int i = 0; i< m_rows; i++)
		{
			for (int j = 0; j < m_columns; j++)
			{
                m_values[i*m_columns+j] = matrix.m_values[i*m_columns+j];
			}
		}
	}

	return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator +(const Matrix<T>& matrix)
{
    assert(m_rows == matrix.m_rows && m_columns == matrix.m_columns);

	Matrix<T> mat(m_rows, m_columns, 0);

	for (int i = 0; i< m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
            mat.m_values[i * m_columns + j] = matrix.m_values[i * m_columns + j] + m_values[i * m_columns + j];
		}
	}
	return mat;
}

template<class T>
Matrix<T> Matrix<T>::operator +(T value)
{
	Matrix<T> mat(m_rows, m_columns, 0);

	for (int i = 0; i< m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
            mat.m_values[i * m_columns + j] =  m_values[i * m_columns + j] + value;
		}
	}
	return mat;
}

template<class T>
Matrix<T> Matrix<T>::operator -(const Matrix<T>& matrix)
{

    assert(m_rows == matrix.m_rows && m_columns == matrix.m_columns);

	Matrix<T> mat(m_rows, m_columns, 0);

	for (int i = 0; i< m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
            mat.m_values[i * m_columns + j] = m_values[i * m_columns + j] - matrix.m_values[i * m_columns + j];
		}
	}
	return mat;
}

template<class T>
Matrix<T> Matrix<T>::operator -(T value)
{
	Matrix<T> mat(m_rows, m_columns, 0);

	for (int i = 0; i< m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
            mat.m_values[i * m_columns + j] = m_values[i * m_columns + j] - value;
		}
	}
	return mat;
}

template<class T>
Matrix<T> Matrix<T>::operator *(const Matrix<T>& matrix)
{
	Matrix<T> mat(m_rows, matrix.m_columns, 0);

	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < matrix.m_columns; j++)
		{
			T value = 0;

			for (int k = 0; k < m_columns; k++)
			{
                value += m_values[i * m_columns + k] * matrix.m_values[k * matrix.m_columns + j];
			}

            mat.m_values[i * matrix.m_columns + j] = value;
		}
	}
	
	return mat;
}

template<class T>
Matrix<T> Matrix<T>::operator *(T value)
{
	Matrix<T> mat(m_rows, m_columns, 0);

	for (int i = 0; i< m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
            mat.m_values[i * m_columns + j] = m_values[i * m_columns + j] * value;
		}
	}
	return mat;
}

template<class T>
Matrix<T> Matrix<T>::operator /(T value)
{
	Matrix<T> mat(m_rows, m_columns, 0);

	for (int i = 0; i< m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
            mat.m_values[i * m_columns + j] = m_values[i * m_columns + j] / value;
		}
	}

	return mat;
}

template<class T>
Matrix<T> Matrix<T>::inverse() const
{
	Matrix<T> matrix(m_rows, m_columns, 0);

	return matrix;
}

template MLTOOLS_EXPORT class Matrix<double>;
template MLTOOLS_EXPORT class Matrix<float>;
template MLTOOLS_EXPORT class Matrix<int>;
