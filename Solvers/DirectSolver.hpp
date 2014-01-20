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
	//���� �������� ������ ������ ����������� � �������� ��������� �� ������ � ���
	void SolveDirectGP(float *in, float *out, float dSigm, float H);
	void SolveDirectGP(float *in, float *out, float *opOut, float dSigm, float H);
	//���� �������� ������ ������ ������������� � �������� ��������� �� ������ � ���
	void SolveDirectMP(float *in, float *out, float dSigm, float H);
	void SolveDirectMP(float *in, float *out, float *opOut, float dSigm, float H);
public:
	//�����������, ������������� �������� ������
	DirectSolver(GridParameters params);
};
