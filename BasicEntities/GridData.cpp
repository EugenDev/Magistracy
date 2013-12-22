#include "GridData.hpp"

void GridData::ReadDSRBStream(ifstream &stream)
{
	int buf_val32 = 0;
	double buf_val64 = 0;
	
	//Header section length
	stream.read((char*)&buf_val32, sizeof(int));
	//Skip header section length
	stream.seekg(buf_val32 + stream.tellg());

	//Grid section header
	stream.read((char*)&buf_val32, sizeof(int));

	if (buf_val32 != 0x44495247)
	{
		//Error: Grid definition missing
		throw "Grid definifion missing";
	}
	
	//Grid definition section length
	stream.read((char*)&buf_val32, sizeof(int));

	//Rows and columns count
	stream.read((char*)&this->gp.NX, sizeof(int));
	stream.read((char*)&this->gp.NY, sizeof(int));

	//Axis initial values
	stream.read((char*)&buf_val64, sizeof(double));
	this->gp.MinX = (float)buf_val64;
	stream.read((char*)&buf_val64, sizeof(double));
	this->gp.MinY = (float)buf_val64;

	//Grid steps
	stream.read((char*)&buf_val64, sizeof(double));
	this->gp.dX = (float)buf_val64;
	stream.read((char*)&buf_val64, sizeof(double));
	this->gp.dY = (float)buf_val64;

	//Useless fields
	stream.read((char*)&zMin, sizeof(double));
	stream.read((char*)&zMax, sizeof(double));
	stream.read((char*)&blank, sizeof(double));
	stream.read((char*)&rotation, sizeof(double));

	//Total array elements
	int total_elements = gp.NX * gp.NY;

	//Data section header
	stream.read((char*)&buf_val32, sizeof(int));

	if (buf_val32 != 0x41544144)
	{
		throw "Grid data missing";
	}

	//Data section size
	int data_size = 0;
	stream.read((char*)&data_size, sizeof(int));
	
	//Reading data
	double* tmpData = new double[data_size];
	stream.read((char*)tmpData, data_size * sizeof(double));
	
	//Converting data to float
	this->data = new float[data_size];
	for(int i = 0; i < data_size; i++)
	{
		data[i] = (float)tmpData[i];
	}

	return;
}

void GridData::ReadGRDFile(string fileName)
{
	ifstream stream(fileName, ios::in | ios::binary);

	if (!stream)
	{
		//TODO: Приписать имя файла
		throw "Error loading grid file";
	}

	int file_version = 0;
	stream.read((char*)&file_version, sizeof(int));

	switch (file_version)
	{
		case 0x42525344:
			this->ReadDSRBStream(stream);
			break;

		default:
			throw "Invalid grid file";
			break;
	}
		
	stream.close();

	return;
}

void GridData::ReadDATFile(string fileName)
{
	ifstream inputFileStream(fileName, ios::binary);
	if (!inputFileStream.is_open())	
		throw ("GridData: No such file or directory");
	
	int linesCount = this->GetLinesCount(&inputFileStream);

	float *tmp1 = new float[linesCount];
	float *tmp2 = new float[linesCount];
	this->data = new float[linesCount];
	
	for (int i = 0; i < linesCount; i++)
		inputFileStream >> tmp1[i] >> tmp2[i] >> data[i];

	this->gp.MinX = tmp1[0];
	this->gp.MinY = tmp2[0];

	float *tmpp = tmp2;
	int i = 0;
	while(*(tmpp++) == *tmp2)
		i++;
	this->gp.NY = linesCount / i;
	this->gp.NX = i;

	this->gp.dX = tmp1[1] - tmp1[0];
	this->gp.dY = tmp2[this->gp.NX] - tmp2[0];

	delete [] tmp1;
	delete [] tmp2;

	inputFileStream.close();
	return;
}

int GridData::GetLinesCount(ifstream* file)
{
	int linesCount = 0;
	float tmp;
	while(!file->eof())
	{
		*file >> tmp >> tmp >> tmp;
		linesCount++;
	}
	linesCount--;
	file->clear();
	file->seekg(0, ios::beg);

	return linesCount;
}

void GridData::SaveToDSRBFile(string filePath)
{
	int buf_val32;
	double buf_val64;

	ofstream output(filePath, ios::binary);

	//Writing file header
	buf_val32 = 0x42525344;
	output.write((char*)&buf_val32, sizeof(int));

	//Header section length
	buf_val32 = 4;
	output.write((char*)&buf_val32, sizeof(int));
	buf_val32 = 0;
	output.write((char*)&buf_val32, sizeof(int));

	//Grid definition section
	buf_val32 = 0x44495247;
	output.write((char*)&buf_val32, sizeof(int));
	buf_val32 = 2 * sizeof(int) + 8 * sizeof(double);
	output.write((char*)&buf_val32, sizeof(int));
	output.write((char*)&gp.NX, sizeof(int));
	output.write((char*)&gp.NY, sizeof(int));
	buf_val64 = gp.MinX;
	output.write((char*)&buf_val64, sizeof(double));
	buf_val64 = gp.MinY;
	output.write((char*)&buf_val64, sizeof(double));
	buf_val64 = gp.dX;
	output.write((char*)&buf_val64, sizeof(double));
	buf_val64 = gp.dY;
	output.write((char*)&buf_val64, sizeof(double));

	//Useless variables
	zMin = 10000000;
	zMax = -10000000;
	for(int i = 0; i < gp.NX * gp.NY; i++)
	{
		if (data[i] > zMax)
			zMax = data[i];
		if (data[i] < zMin)
			zMin = data[i];
	}
	output.write((char*)&zMin, sizeof(double));
	output.write((char*)&zMax, sizeof(double));
	output.write((char*)&blank, sizeof(double));
	output.write((char*)&rotation, sizeof(double));

	//Data section
	int M = gp.NX * gp.NY;
	buf_val32 = 0x41544144;
	output.write((char*)&buf_val32, sizeof(int));
	buf_val32 = M * sizeof(double);
	output.write((char*)&buf_val32, sizeof(int));

	double *tmpBuf = new double[M];
	for(int i = 0; i < M; i++)
	{
		tmpBuf[i] = data[i];
	}
	output.write((char*)tmpBuf, buf_val32);
	delete [] tmpBuf;

	output.close();

	return;
}

void GridData::SaveToDSBBFile(string filePath)
{
	throw "Saving to Surfer 6 format not implemented yet";
	return;
}

void GridData::SaveToDATFile(string filePath)
{
	ofstream outputFile(filePath);

	for(int j = 0; j < gp.NY; j++)
		for(int i = 0; i < gp.NX; i++)
			outputFile << gp.MinX + gp.dX * i << " " 
						<< gp.MinY + gp.dY * j << " " 
						<< data[j * gp.NX + i] << endl;

	outputFile.close();

	return;
}

//-------------------------------Public class members -----------------------------

void GridData::ShowInfo()
{
	cout << "Grid Information" << endl;
	cout << "Rows: " << gp.NX << endl;
	cout << "Columns: " << gp.NY << endl;
	cout << "Initial X: " << gp.MinX << endl;
	cout << "Initial Y: " << gp.MinY << endl;
	cout << "dX: " << gp.dX << endl;
	cout << "dY: " << gp.dX << endl;
	return;
}

GridParameters GridData::GetGridParameters()
{
	return this->gp;
}

GridData::GridData(void)
{
	gp.NX = gp.NY = 1;
	gp.dX = gp.dY = 0.0f;
	gp.MinX = gp.MinY = 0.0f;
	this->data = new float[1];
	data[0] = 0.0f;
}

GridData::GridData(string inputFileName)
{
	int dotPos = inputFileName.find_last_of('.');
	string fileFormat = inputFileName.substr(dotPos + 1);

	if (fileFormat == "dat")
	{
		//.dat file
		this->ReadDATFile(inputFileName);
	}
	else
	{
		//.grd file
		this->ReadGRDFile(inputFileName);
	}
}

GridData::GridData(const GridData &idol)
{
	this->gp = idol.gp;
	this->data = new float [gp.NX * gp.NY];
	memcpy(this->data, idol.data, gp.NX * gp.NY * sizeof(float));
}

GridData::GridData(GridParameters gridParameters)
{
	this->gp = gridParameters;
	this->data = new float [gp.NX * gp.NY];
}

GridData::GridData(GridParameters gridParameters, Matrix matrix)
{
	gp = gridParameters;
	int M = gp.NX * gp.NY;
	if (M != matrix.GetColsCount() * matrix.GetRowsCount())
	{
		throw string("Incompatible matrix and grid parameters");
	}
	else
	{
		this->data = new float[M];
		memcpy(this->data, matrix.elements, M * sizeof(float));
	}
}

Matrix GridData::AsMatrix()
{
	int NX = gp.NX;
	int NY = gp.NY;

	Matrix res(NY, NX);

	memcpy(res.elements, this->data, NX * NY * sizeof(float));

	return res;
}

GridData& GridData::operator = (GridData &other)
{
	GridParameters otherGp = other.GetGridParameters();
	if (otherGp.NX <= 0 || otherGp.NY <= 0)
	{
		throw "How this can be?! GridData with dimensions below zero!!";
	}

	delete [] this->data;	
	this->gp = other.GetGridParameters();
	this->data = new float [gp.NX * gp.NY];
	memcpy(this->data, other.data, gp.NX * gp.NY * sizeof(float));
	return *this;
}

void GridData::FillMatrix(Matrix &m)
{
	if (gp.NX * gp.NY != m.GetColsCount() * m.GetRowsCount())
		throw("GridData: Matrix is too small for this data\n");

	memcpy(m.elements, this->data, gp.NX * gp.NY * sizeof(float));
	return;
}

GridData::~GridData(void)
{
	this->gp.NX = this->gp.NY = 0;
	delete [] data;
}

void GridData::SaveToFile(string filePath, SurferFileFormat fileFormat)
{
	switch(fileFormat)
	{
		case DAT_FILE_FORMAT:
			filePath.append(".dat");
			this->SaveToDATFile(filePath);
			break;
		case SURFER7_FILE_FORMAT:
			filePath.append(".grd");
			this->SaveToDSRBFile(filePath);
			break;
		case SURFER6_FILE_FORMAT:
			filePath.append(".grd");
			this->SaveToDATFile(filePath);
			break;
		default:
			throw "Unrecognized file format";
	}

	return;
}
