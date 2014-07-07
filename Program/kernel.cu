#include <iostream>
using namespace std;
#include <Task.hpp>
#include <GridData.hpp>
#include <CudaDirectSolver.hpp>
#include <ReverseSolver.hpp>
#include <SLESolver.hpp>
#include <cublas.h>
#include <Matrix.hpp>
#include <FinalSolver.hpp>

void TaskOn1000()
{
	GridData layer1(string("D:\\Workspace\\Science\\grids\\mag500\\model1\\layer1.dat"));
	GridData layer2(string("D:\\Workspace\\Science\\grids\\mag500\\model1\\layer2.dat"));
	/*GridData layer1_field(string("D:\\Workspace\\Science\\grids\\mag1000\\model1\\layer1_field.dat"));
	GridData layer2_field(string("D:\\Workspace\\Science\\grids\\mag1000\\model1\\layer2_field.dat"));
	GridData sum(string("D:\\Workspace\\Science\\grids\\mag1000\\model1\\sum.dat"));*/

	GridParameters gp = layer1.GetGridParameters();
	CudaDirectSolver dslvr(gp, false);

	Task t1;
	t1.asimptHeight = 5.0f;
	t1.geltaSigm = 0.25f;
	t1.initialZ = 5.0f;
	t1.taskType = TASK_TYPE_GRAVIMETRY;
	t1.residualType = RESIDUAL_TYPE_EXACTSOLUTION;
	t1.precision = 0.005f;
	t1.grid = layer1;
	t1.exactSolution = layer1;
		
	//dslvr.SolveDirectTask(t1).SaveToFile("D:\\Workspace\\Science\\grids\\mag\\layer1_field", DAT_FILE_FORMAT);
	//dslvr.SolveDirectTask(t1).SaveToFile("D:\\Workspace\\Science\\grids\\mag500\\model1\\layer1_field", DAT_FILE_FORMAT);
		
	Task t2;	
	t2.asimptHeight = 20.0f;
	t2.geltaSigm = 0.3f;
	t2.initialZ = 20.0f;
	t2.taskType = TASK_TYPE_GRAVIMETRY;
	t2.residualType = RESIDUAL_TYPE_EXACTSOLUTION;
	t2.precision = 0.005f;
	t2.grid = layer2;
	t2.exactSolution = layer2;
		
	//dslvr.SolveDirectTask(t2).SaveToFile("D:\\Workspace\\Science\\grids\\mag\\layer2_field", DAT_FILE_FORMAT);
	//dslvr.SolveDirectTask(t2).SaveToFile("D:\\Workspace\\Science\\grids\\mag500\\model1\\layer2_field", DAT_FILE_FORMAT);

	MultilayerTask mTask(layer1);
	mTask.AddTask(t1);
	mTask.AddTask(t2);

	dslvr.SolveDirectMultilayerTask(mTask).SaveToFile("D:\\Workspace\\Science\\grids\\mag500\\model1\\sum", DAT_FILE_FORMAT);
}

int main(int argc, char* argv[])
{
	try
	{
		setlocale(LC_ALL, "Russian");
		cout << "Программа начала работу!!" << endl;
		
		//****************************************************Модель 2*********************************************

		GridData layer1(string("D:\\Workspace\\Science\\grids\\mag1000\\model2\\layer1.dat"));
		GridData layer2(string("D:\\Workspace\\Science\\grids\\mag1000\\model2\\layer2.dat"));
		/*GridData layer1_field(string("D:\\Workspace\\Science\\grids\\mag250\\model2\\layer1_field.dat"));
		GridData layer2_field(string("D:\\Workspace\\Science\\grids\\mag250\\model2\\layer2_field.dat"));
		GridData sum(string("D:\\Workspace\\Science\\grids\\mag250\\model2\\sum.dat"));*/

		GridParameters gp = layer1.GetGridParameters();
		CudaDirectSolver dslvr(gp, false);

		Task t1;
		t1.asimptHeight = 5.0f;
		t1.geltaSigm = 0.25f;
		t1.initialZ = 5.0f;
		t1.taskType = TASK_TYPE_GRAVIMETRY;
		t1.residualType = RESIDUAL_TYPE_EXACTSOLUTION;
		t1.precision = 0.005f;
		t1.grid = layer1;
		t1.exactSolution = layer1;
		
		dslvr.SolveDirectTask(t1).SaveToFile("D:\\Workspace\\Science\\grids\\mag1000\\model2\\layer1_field", DAT_FILE_FORMAT);
		
		Task t2;	
		t2.asimptHeight = 20.0f;
		t2.geltaSigm = 0.3f;
		t2.initialZ = 20.0f;
		t2.taskType = TASK_TYPE_GRAVIMETRY;
		t2.residualType = RESIDUAL_TYPE_EXACTSOLUTION;
		t2.precision = 0.005f;
		t2.grid = layer2;
		t2.exactSolution = layer2;
		
		dslvr.SolveDirectTask(t2).SaveToFile("D:\\Workspace\\Science\\grids\\mag1000\\model2\\layer2_field", DAT_FILE_FORMAT);

		MultilayerTask mTask(layer1);
		mTask.AddTask(t1);
		mTask.AddTask(t2);

		dslvr.SolveDirectMultilayerTask(mTask).SaveToFile("D:\\Workspace\\Science\\grids\\mag1000\\model2\\sum", DAT_FILE_FORMAT);
		
		/*vector<GridData> result = LightMultilayerLinearisedMinimalError(mTask, 0.4f);
		string outName("D:\\Workspace\\Science\\grids\\mag100\\model2\\ans_layer");
		string ss;
		char num;
		int L = result.size();
		for (int i = 1; i <= L; i++)
		{
			ss = outName;
			num = '0' + i;
			ss.append(string(&num).substr(0, 1));
			result[i - 1].SaveToFile(ss, DAT_FILE_FORMAT);
		}*/

		return 0;
		
		//****************************************************Модель 1*********************************************
		//TaskOn1000();
		//return 0;

		/*GridData layer1(string("D:\\Workspace\\Science\\grids\\mag\\layer1.dat"));
		GridData layer2(string("D:\\Workspace\\Science\\grids\\mag\\layer2.dat"));
		GridData layer1_field(string("D:\\Workspace\\Science\\grids\\mag\\layer1_field.dat"));
		GridData layer2_field(string("D:\\Workspace\\Science\\grids\\mag\\layer2_field.dat"));
		GridData sum(string("D:\\Workspace\\Science\\grids\\mag\\sum.dat"));*/

		/*GridData layer1(string("D:\\Workspace\\Science\\grids\\mag200\\model1\\layer1.dat"));
		GridData layer2(string("D:\\Workspace\\Science\\grids\\mag200\\model1\\layer2.dat"));
		GridData layer1_field(string("D:\\Workspace\\Science\\grids\\mag200\\model1\\layer1_field.dat"));
		GridData layer2_field(string("D:\\Workspace\\Science\\grids\\mag200\\model1\\layer2_field.dat"));
		GridData sum(string("D:\\Workspace\\Science\\grids\\mag200\\model1\\sum.dat"));*/

		//GridData layer1(string("D:\\Workspace\\Science\\grids\\mag500\\model1\\layer1.dat"));
		//GridData layer2(string("D:\\Workspace\\Science\\grids\\mag1000\\model1\\layer2.dat"));
		//GridData layer1_field(string("D:\\Workspace\\Science\\grids\\mag500\\model1\\layer1_field.dat"));
		//GridData layer2_field(string("D:\\Workspace\\Science\\grids\\mag1000\\model1\\layer2_field.dat"));
		//GridData sum(string("D:\\Workspace\\Science\\grids\\mag1000\\model1\\sum.dat"));
				

		//GridParameters gp = layer1.GetGridParameters();
		//GridParameters gp = sum.GetGridParameters();
		//CudaDirectSolver dslvr(gp, false);

		/*Task t1;
		t1.asimptHeight = 5.0f;
		t1.geltaSigm = 0.25f;
		t1.initialZ = 5.0f;
		t1.taskType = TASK_TYPE_GRAVIMETRY;
		t1.residualType = RESIDUAL_TYPE_EXACTSOLUTION;
		t1.precision = 0.005f;
		t1.grid = layer1_field;
		t1.exactSolution = layer1;*/
		
		//dslvr.SolveDirectTask(t1).SaveToFile("D:\\Workspace\\Science\\grids\\mag\\layer1_field", DAT_FILE_FORMAT);
		//dslvr.SolveDirectTask(t1).SaveToFile("D:\\Workspace\\Science\\grids\\mag1000\\model1\\layer1_field", DAT_FILE_FORMAT);
		
		/*Task t2;	
		t2.asimptHeight = 20.0f;
		t2.geltaSigm = 0.3f;
		t2.initialZ = 20.0f;
		t2.taskType = TASK_TYPE_GRAVIMETRY;
		t2.residualType = RESIDUAL_TYPE_EXACTSOLUTION;
		t2.precision = 0.005f;
		t2.grid = layer2_field;
		t2.exactSolution = layer2;*/
		
		//dslvr.SolveDirectTask(t2).SaveToFile("D:\\Workspace\\Science\\grids\\mag\\layer2_field", DAT_FILE_FORMAT);
		//dslvr.SolveDirectTask(t2).SaveToFile("D:\\Workspace\\Science\\grids\\mag1000\\model1\\layer2_field", DAT_FILE_FORMAT);
		//return 0;

		/*MultilayerTask mTask(sum);
		mTask.AddTask(t1);
		mTask.AddTask(t2);*/

		//dslvr.SolveDirectMultilayerTask(mTask).SaveToFile("D:\\Workspace\\Science\\grids\\mag\\sum", DAT_FILE_FORMAT);
		//dslvr.SolveDirectMultilayerTask(mTask).SaveToFile("D:\\Workspace\\Science\\grids\\mag100\\model1\\sum", DAT_FILE_FORMAT);
		//return;

		//LightLinearisedMinimalError(t1, 0.5f).SaveToFile("D:\\Workspace\\Science\\grids\\mag\\res1", DAT_FILE_FORMAT);
		//LightLinearisedSpeedDescent(t1, 0.5f).SaveToFile("D:\\Workspace\\Science\\grids\\mag\\res1", DAT_FILE_FORMAT);
		//LinearisedSpeedDescent(&dslvr, t1, 0.5f).SaveToFile("D:\\Workspace\\Science\\grids\\mag\\res1", DAT_FILE_FORMAT);
		//return;
		
		//LightLevenbergMarkvardt(t1, 0.5f).SaveToFile("D:\\Workspace\\Science\\grids\\mag\\result1", DAT_FILE_FORMAT);
		//return;
		//vector<GridData> result = LightMultilayerLevenbergMarkvardt(mTask, 1.0f);
		/*vector<GridData> result = LightMultilayerLinearisedMinimalError(mTask, 1.0f);
		string outName("D:\\Workspace\\Science\\grids\\mag\\ans_layer");
		string ss;
		char num;
		int L = result.size();
		for (int i = 1; i <= L; i++)
		{
			ss = outName;
			num = '0' + i;
			ss.append(string(&num).substr(0, 1));
			result[i - 1].SaveToFile(ss, DAT_FILE_FORMAT);
		}
		return 0;*/
	}
	catch(string s)
	{
		cout << s << endl;
	}
	catch (char *str)
	{
		cout << str << endl;
	}

	system("PAUSE");

	return 0;
}
