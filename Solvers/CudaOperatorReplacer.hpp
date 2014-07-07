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
	//Параметры сетки
	GridParameters gp;
	//Задачи
	Task task;
	MultilayerTask *mTask;
	//Признак многослойной задачи (наверное не нужен, просто всегда предполагаем многослойность)
	bool isMultilayer;

	//Число строк истолбцов матрицы
	int nRows;
	int nCols;

	//Указатель на область памяти на GPU где храняться асимптотические высоты 
	float *heights;

	//Указатель на область памяти на GPU где храняться скачки плотности
	float *sigmas;

	//Число слоёв
	int layersCount;

public :
	CudaOperatorReplacer(Task t);
	CudaOperatorReplacer(MultilayerTask t);

	~CudaOperatorReplacer();

	//Функции умножения вектор на за матрицу - заменитель ...
	void CalcAX(float *Z, float *X, float *result);
	// ... и вариант для транспонированной матрицы
	void CalcATX(float *Z, float *X, float *result);

	//Функция заполняет матрицу
	void GetMatrix(Matrix &Res, float* Z);
	//Заполнение транспонированой матрицы
	void GetTransposedMatrix(Matrix &Res, float* Z);
};