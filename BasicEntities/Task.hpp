#pragma once
#include "Matrix.hpp"
#include "GridData.hpp"

enum TaskType 
{
	TASK_TYPE_MAGNITOMETRY, //Задача магнитометрии
	TASK_TYPE_GRAVIMETRY	  //Задача гравиметрии
};

enum ResidualType
{
	RESIDUAL_TYPE_RIGHTHAND,				//Счёт по невязке
	RESIDUAL_TYPE_EXACTSOLUTION			//Счёт по погрешности
};

struct Task
{
	GridData grid;				//Сетка, на котрой решаем задачу
	TaskType taskType;			//Тип задачи
	float asimptHeight;			//Асимптотическая плоскость
	float initialZ;				//Начальное приближение
	float geltaSigm;			//Cкачок плотности между средами
	float geltaJ;				//Скачок вертикальной составляющей вектора намагниченности
	float precision;			//Требуемая точность
	GridData exactSolution;	    //Точное решение
	ResidualType residualType;	//Тип 
};
