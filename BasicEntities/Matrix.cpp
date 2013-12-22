#include "Matrix.hpp"

Matrix::Matrix()
{
	this->col = 0;
	this->rows = 0;
	this->elements = NULL;
}

Matrix::Matrix(int rows, int col)
{
	if (rows < 1 || col < 1)
		return;
	this->rows = rows;
	this->col = col;
	this->SafeAlloc();
}

Matrix::Matrix(const Matrix &other) : rows(other.rows), col(other.col)
{
	this->SafeAlloc();
	for (int i = 0; i < rows * col; i++)
		elements[i] = other.elements[i];
}

Matrix Matrix::operator + (const Matrix &other) const
{
	if (this->rows != other.rows || this->col != other.col)
		throw string("Matrix error: Matrix dimensions must agree");

	Matrix result(*this);

	for (int i = 0; i < rows * col; i++)
		result.elements[i] += other.elements[i];

	return result;
}

Matrix& Matrix::operator += (const Matrix &other)
{
	if (this->rows != other.rows || this->col != other.col)
		throw string("Matrix error: Matrix dimensions must agree");

	for (int i = 0; i < rows * col; i++)
		this->elements[i] += other.elements[i];

	return *this;
}

Matrix Matrix::operator - (const Matrix &other) const
{
	if (this->rows != other.rows || this->col != other.col)
		throw string("Matrix error: Matrix dimensions must agree");

	Matrix result(*this);

	for (int i = 0; i < rows * col; i++)
		result.elements[i] -= other.elements[i];

	return result;
}

Matrix& Matrix::operator -= (const Matrix &other)
{
	if (this->rows != other.rows || this->col != other.col)
		throw string("Matrix error: Matrix dimensions must agree");

	for (int i = 0; i < rows * col; i++)
		this->elements[i] -= other.elements[i];

	return *this;
}

Matrix Matrix::operator * (const float &value) const
{
	Matrix result(*this);

	for (int i = 0; i < rows * col; i++)
		result.elements[i] *= value;

	return result;
}

Matrix& Matrix::operator *= (const float &value)
{
	for (int i = 0; i < rows * col; i++)
		this->elements[i] *= value;
	return *this;
}

Matrix Matrix::operator * (const Matrix &other) const
{
	int m1 = this->rows;
	int n1 = this->col;
	int m2 = other.rows;
	int n2 = other.col;

	if (n1 != m2)
		throw string("Matrix error: Matrix dimensions must agree");

	Matrix result(m1, n2);
	result.Fill(0.0f);

	float *C = result.elements;
	float *A = this->elements;
	float *B = other.elements;

	for (int i = 0; i < m1; i++)
	{
		for (int k = 0; k < m2; k++)
		{
			for (int j = 0; j < n2; j++)
			{
				C[i * n2 + j] += A[i * n1 + k] * B[k * n2 + j];
			}
		}
	}

	return result;
}

Matrix& Matrix::operator = (const Matrix &other)
{
	if ( this != &other && ( this->col != other.GetColsCount() || this->rows != other.GetRowsCount() ) )
	{
		delete [] this->elements;
		this->col = other.GetColsCount();
		this->rows = other.GetRowsCount();
		this->elements = new float[col * rows];
	}
	memcpy(this->elements, other.elements, rows * col * sizeof(float));
	return *this;
}

float& Matrix::operator[] (unsigned i)
{
	return elements[i];
}

float Matrix::Norm() const
{
	if (this->col != 1)
		throw string("Matrix error: Norm defined only for vector");

	float res = 0.0f;

	for (int i = 0; i < this->rows; i++)
		res += this->elements[i] * this->elements[i];

	return sqrt(res);
}

float Matrix::Dot(const Matrix &other) const 
{
	if (this->col != 1 || other.GetColsCount() != 1)
		throw string("Matrix error: Dot defined only for vectors");
	if (this->rows != other.GetRowsCount())
		throw string("Matrix error: Dot defined only for same length vectors");

	float result = 0.0f;
	for (int i = 0; i< this->rows; i++)
		result += this->elements[i] * other.elements[i];

	return result;
}

int Matrix::GetColsCount() const
{
	return this->col;
}

int Matrix::GetRowsCount() const
{
	return this->rows;
}

void Matrix::Show()
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < col; j++)
			cout << this->elements[i * col + j] << " ";
		cout << endl;
	}
	return;
}

void Matrix::Fill(float value)
{
	for (int i = 0; i < rows * col; i++)
		this->elements[i] = value;
	return;
}

void Matrix::AddToDiagonal(float value)
{
	if (rows != col)
		return;
	
	for (int i = 0; i < rows; i++)
		this->elements[i * col + i] += value;

	return;
}

Matrix::~Matrix(void)
{
	delete [] elements;
}

void Matrix::SafeAlloc()
{
	try
	{
		this->elements = new float [rows * col];
	}
	catch(std::bad_alloc&)
	{
		throw string("Not enogh memory for storage matrix");
	}
}
