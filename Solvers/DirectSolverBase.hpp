#pragma once
#include "GridData.hpp"
#include "Matrix.hpp"
#include "Task.hpp"
#include <vector>
#include "MultilayerTask.hpp"

class DirectSolverBase
{
protected:
	//��� ����� �������� ������ ������
	float coeff;
	//��������� �����
	GridParameters gridParams;
	//���� �������� ������ ������ ����������� � �������� ��������� �� ������ � ���
	virtual void SolveDirectGP(float *in, float *out, float dSigm, float H) = 0;
	virtual void SolveDirectGP(float *in, float *out, float *opOut, float dSigm, float H) = 0;
	//���� �������� ������ ������ ������������� � �������� ��������� �� ������ � ���
	virtual void SolveDirectMP(float *in, float *out, float dJ, float H) = 0;
	virtual void SolveDirectMP(float *in, float *out, float *opOut, float dJ, float H) = 0;
public:
	//�����������, ������������� �������� ������
	DirectSolverBase();
	DirectSolverBase(GridParameters params);
	//������� ������ ������, �������� � ���������� task
	GridData SolveDirectTask(Task task);
	//�� ��, �� �� ������� �� ������
	void SolveDirectTask(const Matrix &Z, Matrix &ans, Matrix &opDer, Task task);
	void SolveDirectTask(const Matrix &Z, Matrix &ans, Task task);
	//������� ��� ���������� ����
	void SolveDirectMultilayerTask(MultilayerTask mTask, Matrix &Z, Matrix &A, Matrix &Ank);
	void SolveDirectMultilayerTask(MultilayerTask mTask, Matrix &Z, Matrix &A);
	GridData SolveDirectMultilayerTask(MultilayerTask mTask);
};

