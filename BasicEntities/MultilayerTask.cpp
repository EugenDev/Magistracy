#include "MultilayerTask.hpp"

void MultilayerTask::AddTask(Task new_task)
{
	vector<Task>::iterator item;
	for (item = _tasks.begin(); item != _tasks.end(); item++)
	{
		if (new_task.grid.GetGridParameters() != item->grid.GetGridParameters())
		{
			throw "New task grid is incompatible with existing task grids";
		}
	}

	_tasks.push_back(new_task);

	return;
}

GridParameters MultilayerTask::GetGeneralGridParameters()
{
	if (_tasks.size() != 0)
	{
		return _tasks[0].grid.GetGridParameters();
	}
	else
	{			
		return GridParameters();
	}
}

int MultilayerTask::GetLayersCount()
{
	return _tasks.size();
}

GridData MultilayerTask::GetGeneralField()
{
	return this->generalField;
}

Task MultilayerTask::operator[] (unsigned i)
{
	return _tasks[i];
}

MultilayerTask::MultilayerTask(GridData generalField)
{
	this->generalField = generalField;
}

MultilayerTask::~MultilayerTask(void)
{
}
