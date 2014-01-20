#pragma once
#include "GridData.hpp"
#include "Matrix.hpp"
#include "Task.hpp"
#include <vector>
#include "MultilayerTask.hpp"

class DirectSolverBase
{
protected:
	//Для более короткой записи формул
	float coeff;
	//Параметры сетки
	GridParameters gridParams;
	//Ядра решателя прямой задачи гравиметрии с матрицей оператора на выходе и без
	virtual void SolveDirectGP(float *in, float *out, float dSigm, float H) = 0;
	virtual void SolveDirectGP(float *in, float *out, float *opOut, float dSigm, float H) = 0;
	//Ядра решателя прямой задачи магнитометрии с матрицей оператора на выходе и без
	virtual void SolveDirectMP(float *in, float *out, float dJ, float H) = 0;
	virtual void SolveDirectMP(float *in, float *out, float *opOut, float dJ, float H) = 0;
public:
	//Конструктор, инициализация решателя сеткой
	DirectSolverBase();
	DirectSolverBase(GridParameters params);
	//Решение прямой задачи, описаной в переменной task
	GridData SolveDirectTask(Task task);
	//То же, но со входами из матриц
	void SolveDirectTask(const Matrix &Z, Matrix &ans, Matrix &opDer, Task task);
	void SolveDirectTask(const Matrix &Z, Matrix &ans, Task task);
	//Решение для нескольких слоёв
	void SolveDirectMultilayerTask(MultilayerTask mTask, Matrix &Z, Matrix &A, Matrix &Ank);
	void SolveDirectMultilayerTask(MultilayerTask mTask, Matrix &Z, Matrix &A);
	GridData SolveDirectMultilayerTask(MultilayerTask mTask);
};

