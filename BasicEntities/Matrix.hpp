#pragma once
#include <iostream>
using namespace std;
#include <math.h>

#define dllExport _declspec(dllexport)

class _declspec(dllexport) Matrix
{
private:
	int col;  //Число строк
	int rows; //Число столбцов
	void SafeAlloc();
public:
	//Элементы матрицы
	float *elements;

	//Перегрузка операций
	//Сложение
	Matrix operator + (const Matrix &other) const;
	Matrix& operator += (const Matrix &other);
	//Вычитание
	Matrix operator - (const Matrix &other) const;
	Matrix& operator -= (const Matrix &other);
	//Умножение на скаляр
	Matrix operator * (const float &value) const;
	Matrix& operator *= (const float &value);
	//Умножение на матрицу
	Matrix operator * (const Matrix &other) const;
	//Присваивания
	Matrix& operator = (const Matrix &other);
	//Индексация
	float& operator [] (unsigned i);

	//Норма вектора
	float Norm() const;
	//Скалярное произведение
	float Dot(const Matrix& other) const;

	//Конструкторы
	Matrix();
	Matrix(int rows, int col);
	Matrix(const Matrix &other);

	//Аксессоры к закрытым полям
	int GetColsCount() const;
	int GetRowsCount() const;

	//Печать матрицы на экран
	void Show();

	//Заполнить матрицу одинаковыми числами
	void Fill(float value);
	void AddToDiagonal(float value);

	Matrix Transpose();

	//Деструктор
	~Matrix(void);
};
