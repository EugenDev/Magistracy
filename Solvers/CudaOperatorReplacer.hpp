#include <cublas.h>
#include <GridData.hpp>
#include <Task.hpp>
#include <MultilayerTask.hpp>
#include <device_launch_parameters.h>
#include "CudaHelpers.hpp"
#define GRAVY_CONST 6.673

class CudaOperatorReplacer
{
private: 
	//��������� �����
	GridParameters gp;
	//������
	Task task;
	MultilayerTask *mTask;
	//������� ������������ ������ (�������� �� �����, ������ ������ ������������ ��������������)
	bool isMultilayer;

	//����� ����� ��������� �������
	int nRows;
	int nCols;

	//��������� �� ������� ������ �� GPU ��� ��������� ��������������� ������ 
	float *heights;

	//��������� �� ������� ������ �� GPU ��� ��������� ������ ���������
	float *sigmas;

	//����� ����
	int layersCount;

public :
	CudaOperatorReplacer(Task t);
	CudaOperatorReplacer(MultilayerTask t);

	~CudaOperatorReplacer();

	//������� ��������� ������ �� �� ������� - ���������� ...
	void CalcAX(float *Z, float *X, float *result);
	// ... � ������� ��� ����������������� �������
	void CalcATX(float *Z, float *X, float *result);

	//������� ��������� �������
	void GetMatrix(Matrix &Res, float* Z);
	//���������� ���������������� �������
	void GetTransposedMatrix(Matrix &Res, float* Z);
};