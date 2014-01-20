#pragma once
#include "Task.hpp"
#include <vector>

class MultilayerTask
{
private:
	vector<Task> _tasks;
	GridData generalField;
public:
	//Single layer task addition
	void AddTask(Task task);
	//Returns general grid parameter for all tasks
	GridParameters GetGeneralGridParameters();
	GridData GetGeneralField();
	//Returns count of tasks
	int GetLayersCount();
	//Indexer
	Task MultilayerTask::operator[] (unsigned i);
	//Constructor and destructor
	MultilayerTask(GridData generalField);
	~MultilayerTask(void);
};

