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

void MultilayerTask::InitZ(Matrix &Z)
{
	GridParameters gp = this->generalField.GetGridParameters();
	int L = this->GetLayersCount();
	int M = gp.NX * gp.NY;
	if (Z.GetColsCount() * Z.GetRowsCount() != L * M) 
	{
		throw new string("Inapropriate matrix");
	}

	for(int l = 0; l < L; l++)
	{
		std::fill(Z.elements + l * M, Z.elements + (l + 1) * M, _tasks[l].initialZ);
	}

	return;
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
