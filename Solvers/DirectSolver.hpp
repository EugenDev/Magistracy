#pragma once
#include "GridData.hpp"
#include "Matrix.hpp"
#include "Task.hpp"
#include <vector>
#include "MultilayerTask.hpp"
#include "DirectSolverBase.hpp"
#include <omp.h>

class DirectSolver : public DirectSolverBase
{
private:
	//Ядра решателя прямой задачи гравиметрии с матрицей оператора на выходе и без
	void SolveDirectGP(float *in, float *out, float dSigm, float H);
	void SolveDirectGP(float *in, float *out, float *opOut, float dSigm, float H);
	//Ядра решателя прямой задачи магнитометрии с матрицей оператора на выходе и без
	void SolveDirectMP(float *in, float *out, float dSigm, float H);
	void SolveDirectMP(float *in, float *out, float *opOut, float dSigm, float H);
public:
	//Конструктор, инициализация решателя сеткой
	DirectSolver(GridParameters params);
};
