#pragma once
#include <iostream>
using namespace std;
#include <math.h>

#define dllExport _declspec(dllexport)

class _declspec(dllexport) Matrix
{
private:
	int col;  //����� �����
	int rows; //����� ��������
	void SafeAlloc();
public:
	//�������� �������
	float *elements;

	//���������� ��������
	//��������
	Matrix operator + (const Matrix &other) const;
	Matrix& operator += (const Matrix &other);
	//���������
	Matrix operator - (const Matrix &other) const;
	Matrix& operator -= (const Matrix &other);
	//��������� �� ������
	Matrix operator * (const float &value) const;
	Matrix& operator *= (const float &value);
	//��������� �� �������
	Matrix operator * (const Matrix &other) const;
	//������������
	Matrix& operator = (const Matrix &other);
	//����������
	float& operator [] (unsigned i);

	//����� �������
	float Norm() const;
	//��������� ������������
	float Dot(const Matrix& other) const;

	//������������
	Matrix();
	Matrix(int rows, int col);
	Matrix(const Matrix &other);

	//��������� � �������� �����
	int GetColsCount() const;
	int GetRowsCount() const;

	//������ ������� �� �����
	void Show();

	//��������� ������� ����������� �������
	void Fill(float value);
	void AddToDiagonal(float value);

	Matrix Transpose();

	//����������
	~Matrix(void);
};
