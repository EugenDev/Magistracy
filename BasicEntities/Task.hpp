#pragma once
#include "Matrix.hpp"
#include "GridData.hpp"

enum TaskType 
{
	TASK_TYPE_MAGNITOMETRY, //������ �������������
	TASK_TYPE_GRAVIMETRY	  //������ �����������
};

enum ResidualType
{
	RESIDUAL_TYPE_RIGHTHAND,				//���� �� �������
	RESIDUAL_TYPE_EXACTSOLUTION			//���� �� �����������
};

struct Task
{
	GridData grid;				//�����, �� ������ ������ ������
	TaskType taskType;			//��� ������
	float asimptHeight;			//��������������� ���������
	float initialZ;				//��������� �����������
	float geltaSigm;			//C����� ��������� ����� �������
	float geltaJ;				//������ ������������ ������������ ������� ���������������
	float precision;			//��������� ��������
	GridData exactSolution;	    //������ �������
	ResidualType residualType;	//��� 
};
