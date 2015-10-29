#ifndef MATRIX_H
#define MATRIX_H

#include "mltools_global.h"
#include <iostream>
#include <assert.h>

using namespace std;

template<class T>
class MLTOOLS_EXPORT Matrix
{

public:
    Matrix(int rows = 1, int columns = 1, T defaultValue = 0);

    Matrix(const Matrix<T>& matrix);

    ~Matrix();
	
    int rows() const;

    int columns() const;

	bool isRowVector() const;

	bool isColumnVector() const;

	bool isSquareMatrix() const;

	Matrix<T> transpose() const;

	Matrix<T> identity() const;

	Matrix<T> getRow(int i) const;

	void setRow(const Matrix<T>& row, int i);

	void appendAsRows(const Matrix<T>& rows);

	Matrix<T> getColumn(int j) const;

	void setColumn(const Matrix<T>& column, int j);

	void appendAsColumns(const Matrix<T>& columns);

	T operator ()(int i, int j) const;

	T& operator ()(int i, int j);

	Matrix<T>& operator = (const Matrix<T>& matrix);

	Matrix<T> operator + (const Matrix<T>& matrix);

	Matrix<T> operator + (T value);

	Matrix<T> operator - (const Matrix<T>& matrix);

	Matrix<T> operator - (T value);

	Matrix<T> operator * (const Matrix<T>& matrix);

	Matrix<T> operator * (T value);

	Matrix<T> operator / (T value);

	Matrix<T> inverse() const;

	static Matrix<T> ones(int rows, int cols, T value);

	static Matrix<T> diagonal(int size, T value);

	static Matrix<T> identity(int size);

	friend ostream& operator << (ostream& out, const Matrix<T>& matrix)
	{
		for (int i = 0; i < matrix.m_rows; i++)
		{
			for (int j = 0; j < matrix.m_columns; j++)
			{
				if (j == 0)
				{
                    out << matrix.m_values[i * matrix.m_columns + j];
				}
				else
				{
                    out << "\t" << matrix.m_values[i * matrix.m_columns + j];
				}

				if (j == matrix.m_columns - 1)
					out << "\n";
			}
		}

		return out;
	}

	friend istream& operator >> (istream& in, Matrix<T>& matrix)
	{
		for (int i = 0; i < matrix.m_rows; i++)
		{
			for (int j = 0; j < matrix.m_columns; j++)
			{
				in >> matrix.m_values[i][j];
			}
		}

		return in;
	}

	
private:
    T* m_values;
    int m_rows, m_columns;
};

#endif // MATRIX_H
