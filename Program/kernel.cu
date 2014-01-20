#include <iostream>
using namespace std;
#include <Task.hpp>
#include <GridData.hpp>
#include <CudaDirectSolver.hpp>
#include <ReverseSolver.hpp>

int main(int argc, char* argv[])
{
	GridData layer1(string("D:\\Workspace\\Science\\grids\\gravy\\small2\\gr45x50_layer1.dat"));
	GridData layer2(string("D:\\Workspace\\Science\\grids\\gravy\\small2\\gr45x50_layer2.dat"));
	GridData diploma(string("D:\\Workspace\\Science\\grids\\magnito\\diploma\\17.dat"));

	//GridData layer1_field(string("D:\\Workspace\\Science\\grids\\gravy\\small2\\gr45x50_layer1_field.dat"));
	//GridData layer2_field(string("D:\\Workspace\\Science\\grids\\gravy\\small2\\gr45x50_layer2_field.dat"));

	Task t;
	/*t.asimptHeight = 5.0f;
	t.geltaSigm = 0.3f;
	t.grid = layer1;
	t.initialZ = 5.0f;
	t.taskType = TASK_TYPE_GRAVIMETRY;

	t.asimptHeight = 10.0f;
	t.geltaSigm = 0.3f;
	t.grid = layer2;
	t.initialZ = 10.0f;
	t.taskType = TASK_TYPE_GRAVIMETRY;*/

	t.asimptHeight = 6.0f;
	t.geltaJ = 2000.0f;
	t.grid = diploma;
	t.exactSolution = diploma;
	t.initialZ = 6.0f;
	t.taskType = TASK_TYPE_MAGNITOMETRY;
	t.residualType = RESIDUAL_TYPE_RIGHTHAND;
	t.precision = 0.01f;

	/*CudaDirectSolver dslvr1(layer2.GetGridParameters(), false);
	dslvr1.SolveDirectTask(t).SaveToFile("D:\\Workspace\\Science\\grids\\gravy\\small2\\ans1", DAT_FILE_FORMAT);*/
	
	G_LinearisedMinimalError(t, 1.0f).SaveToFile("D:\\Workspace\\Science\\grids\\magnito\\diploma\\layer", DAT_FILE_FORMAT);

	return 0;
}
