/*
Term Project 01		====================================
Hoon Seok Kim		====================================
STD #: 200452816	====================================
*/

// Headers
#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <bitset>
using namespace std;

// Image size definitions
#define ROWS 256
#define COLS 256

// Buffer for input images
unsigned char in_img[ROWS][COLS];

// Buffer for input images type change => complex<double>
complex<double> C_img[ROWS][COLS];

// Temporary Buffer for shuffling array
complex<double> S_img[ROWS][COLS];		// horizontal shuffle
complex<double> S_img2[ROWS][COLS];		// vertical shuffle

// Buffer for saving original DFT
complex<double> DFT_ori[ROWS][COLS];	// raw DFT values

// Temporary Buffer for FFT images
complex<double> DFT[ROWS][COLS];		// F(0) centered DFT values
complex<double> DFT1[ROWS][COLS];		// For FFT 1D
complex<double> DFT2[ROWS][COLS];		// For FFT 2D

// Buffer for output images
unsigned char out_img[ROWS][COLS];

// Definition of pi
double pi = 2 * asin(1);

// number of elements in a row/column
int N = ROWS;
//===============================================================//
//========================End of variables=======================//
//===============================================================//



//===============================================================//
//========================Sub-functions==========================//
//===============================================================//

// convert raw image to complex image::: in_img => C_img==================
void convert()
{
	for (int i = 0; i < ROWS; i++)
	{
		for (int j = 0; j < COLS; j++)
		{
			// unsigned char to complex<double>: (real, imaginary)
			C_img[i][j] = complex<double>(in_img[i][j], 0.0);
		}
	}
}

// Column index shuffle for every Row ==> output: S_img[][]=============== 
void Rshuffle(complex<double> Comimg[ROWS][COLS])
{
	string bin;				// string for binary number
	int indextemp = 0;		// initilizae shuffling index to 0.
	// for loop for entire rows
	for (int i = 0; i < ROWS; i++)
	{
		// checking each Column index
		for (int j = 0; j < COLS; j++)
		{
			// convert integer j(column index) to binary string 
			// and save log2(N)=8 bit binary to bin
			// ****number of bits for index 0-255 is 8
			bin = bitset<8>(j).to_string();
			// convert binary to decimal reversly 
			// and save recalculated decimal to indextemp
			for (int z = 0; z < 8; z++)
			{
				// If any binary number is 1, convert it to decimal value
				// And add them to indextemp
				if (bin[z] == '1')
				{
					// Binary to decimal in reverse direction
					// 2^0 for MSB ....... 2^7 for LSB
					indextemp = indextemp + pow(2, z);
				}
			}
			// Element at Re-esimated, shuffled column index
			// Goes to element of current scanning index
			S_img[i][j] = Comimg[i][indextemp];
			indextemp = 0;			// Reset indextemp for next column index
		}
	}
}

// Row index shuffle for every column=====================================
void Cshuffle(complex<double> Comimg[ROWS][COLS])
{
	string bin;				// string for binary number
	int indextemp = 0;		// initilizae shuffling index to 0.
	// for loop for entire columns
	for (int i = 0; i < COLS; i++)
	{
		// checking each row index
		for (int j = (ROWS-1); j >= 0; j--)
		{
			// convert row index=j to binary string and save to bin
			// and save log2(N)=8 bit binary to bin
			// ****number of bits for index 0-255 is 8
			bin = bitset<8>(j).to_string();
			// convert binary to decimal reversly 
			// and save recalculated decimal to indextemp
			for (int z = 0; z < 8; z++)
			{
				if (bin[z] == '1')
				{
					// Binary to decimal in reverse direction
					// 2^0 for MSB ....... 2^7 for LSB
					indextemp = indextemp + pow(2, z);
				}
			}
			// Element at Re-esimated, shuffled row index
			// Goes to element of current scanning row index
			S_img2[j][i] = Comimg[indextemp][i];
			indextemp = 0;			// Reset indextemp for next column index
		}
	}
}

// function for calculating W_n^u========================================
complex<double> W(int M, int u)
{
	// real part: cos(-pi*u/M) = cos(pi*u/M)
	// Imaginary part: sin(-pi*u/M) = -sin(pi*u/M)
	double real = cos((pi * u) / M);
	double imaginary = -sin((pi * u) / M);
	// If real calculated real/imaginary numbers are smaller than 0.0001
	// Consider it as 0 for calculation simplicity
	if (abs(real) < 0.0001)
	{
		real = 0.0;
	}
	if (abs(imaginary) < 0.0001)
	{
		imaginary = 0.0;
	}
	// Return euler fomula as complex<double>
	return(complex<double>(real, imaginary));
}

// 1D FFT function for single row =======================================
// Basically the same with lecture note==================================
void FFT1D(int z, complex<double> Comimg[ROWS][COLS])
{
	double n = 0.0;			// initialze integer n
	n = log2(N);			// n = Number of merging stages
							// N is array size for a single row
	int M = 1;				// initial subgroup length
	int pair = N / 2;		// initial number of pairs of subgroup

	// FFT Calculation 
	for (int i = 0; i < n; i++)				// Conducting merging n times.
	{
		for (int k = 0; k < pair; k++)		// Merge pairs at the current level
		{
			int start1 = k * 2 * M;			// start index for first group
			int start2 = (k * 2 + 1) * M;	// start index for second group
			// Fast fourier transform
			for (int u = 0; u < M; u++)
			{
				// Temporarily save intermediate FFT values to DFT1 
				DFT1[z][u] = 0.5 * (Comimg[z][start1 + u] + Comimg[z][start2 + u] * W(M, u));
				DFT1[z][u + M] = 0.5 * (Comimg[z][start1 + u] - Comimg[z][start2 + u] * W(M, u));
			}
			// update the FFT rusult at input Aarray
			for (int u = 0; u < 2 * M; u++)
			{
				Comimg[z][start1 + u] = DFT1[z][u];
			}
		}
		// double the subgroup length for next stage
		M = 2 * M;
		// reduce number of groups by half for next stage
		pair = pair / 2;
	}
}

// FFT function for single column =======================================
// updated from FFT1D ===================================================
void FFT2D(int z, complex<double> Comimg[ROWS][COLS])
{
	double n = 0.0;			// initialze integer n
	n = log2(N);			// Number of merge stages
							// N is array size for a single row
	int M = 1;				// initial subgroup length
	int pair = N / 2;		// initial number of pairs of subgroup

	// FFT Calculation 
	for (int i = 0; i < n; i++)				// Conducting merging n times.
	{
		for (int k = 0; k < pair; k++)		// Merge pairs at the current level
		{
			int start1 = k * 2 * M;			// start index for first group
			int start2 = (k * 2 + 1) * M;	// start index for second group
			// Fast fourier transform for row index
			for (int u = 0; u < M; u++)
			{
				// Temporarily save intermediate FFT values to DFT2 
				DFT2[u][z] = 0.5 * (Comimg[start1 + u][z] + Comimg[start2 + u][z] * W(M, u));
				DFT2[u + M][z] = 0.5 * (Comimg[start1 + u][z] - Comimg[start2 + u][z] * W(M, u));
			}
			// update the FFT rusult at input Aarray
			for (int u = 0; u < 2 * M; u++)
			{
				Comimg[start1 + u][z] = DFT2[u][z];
			}
		}
		// double the subgroup length for next stage
		M = 2 * M;
		// reduce number of groups by half for next stage
		pair = pair / 2;
	}
}

// Save FFT Result to DFT_ori[][] =========================================
// Scale FFT result and save to out_img[][] ===============================
void DFTresult(complex<double> Comimg[ROWS][COLS])
{	
	double max = 0.0;			// Initialize max value of input
	double min = 2000.0;		// Initialize min value of input
	// Scanning row
	for (int i = 0; i < ROWS; i++)
	{	// scannning column
		for (int j = 0; j < COLS; j++)
		{	// Compare current value with max and if needed
			if (max < abs(Comimg[i][j]))
			{
				max = abs(Comimg[i][j]);
			}
			// Compare current value with max and if needed
			if (min > abs(Comimg[i][j]))
			{
				min = abs(Comimg[i][j]);
				// If scanning value is too small < 0.0001
				// Set it as 0 for simplicity
				if (abs(Comimg[i][j])<0.0001)
				{
					min = 0.0;
				}
			}
			// Center the input FFT results and save to DFT[][]
			DFT[(i+128)%256][(j+128)%256] = Comimg[i][j];
			// Save original FFT reuslts to DFT_ori[][]
			DFT_ori[i][j] = Comimg[i][j];
		}
	}
	// Print out Max, min, and scaling factor
	cout << max << " " << min << " " << (255.0 / (max - min)) << endl;
	// save scaled DFT[][] to out_img[][] as unsigned char
	for (int i = 0; i < ROWS; i++)
	{
		for (int j = 0; j < COLS; j++)
		{
			out_img[i][j] = (unsigned char)(int)((255.0 / (max - min)) * (abs(DFT[i][j]) - min));
		}
	}
}
//===============================================================//
//---------------- end of sub-functions -------------------------//
//===============================================================//


//=================================================================
//=================================================================
//===== Main function =============================================
//=================================================================
int main()
{
//-------------------------------------------------------------------------
//--------------- File managing--------------------------------------------
//-------------------------------------------------------------------------
	FILE* fin, *fout;
	int n = 0;

	//-----------------------------------------------------------
	//*****************************Select image part*************
	//-----------------------------------------------------------

	// Open the input image file
//	if (fopen_s(&fin, "square256.raw", "rb") != 0) {			// Enable this line only for 1st image
	if (fopen_s(&fin, "car.raw", "rb") != 0) {					// Enable this line only for 2nd image
		fprintf(stderr, "ERROR: Cann't open input image file 'square256.raw.'\n");
	}
	// Open the output image file
//	if (fopen_s(&fout, "square256DFT.raw", "wb") != 0) {		// Enable this line only for 1st image
	if (fopen_s(&fout, "carDFT.raw", "wb") != 0) {				// Enable this line only for 2nd image
		fprintf(stderr, "ERROR: Cann't open output image file 'square256DFT.raw.'\n");
	}
	//*************************end of Select part****************
	//-----------------------------------------------------------

	// Load the input image
	printf("... Load input image1\n");
	n = fread(in_img, sizeof(char), ROWS * COLS, fin);
	if (n < ROWS * COLS * sizeof(char)) {
		fprintf(stderr, "ERROR: Read first input image file error\n");
		return 1;
	}

//---------------------------------------------------------------------------
//--------------- File managing finished-------------------------------------
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Image processing starts---------------------------------------------------
//---------------------------------------------------------------------------

	// Converting RAW images to complex images
	convert();				// in_img[][] => C_img[][]
	// Shuffle image column index
	Rshuffle(C_img);		// C_img[][] => S_img[][]
	// FFT for every row 
	for (int z = 0; z < ROWS; z++)
	{
		FFT1D(z, S_img);	// S_img[][] => S_img[][]
	}
	// Shuffle image row index
	Cshuffle(S_img);		// S_img[][] => S_img2[][]
	// FFT for every column
	for (int z = 0; z < COLS; z++)
	{
		FFT2D(z, S_img2);	// S_img2[][] => S_img2[][]
	}
	// Covert FFT results into several formats
	DFTresult(S_img2);		// S_img2[][] => out_img[][] ::: scaled raw image
							// S_img2[][] => DFT_ori[][] ::: Original DFT result

//---------------------------------------------------------------------------
// Image processing ends ----------------------------------------------------
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------
//--------------- output managing starts-------------------------------------
//---------------------------------------------------------------------------
	
	// Store out_img[][] to the output image 
	printf("... Save the output image\n");
	n = fwrite(out_img, sizeof(char), (ROWS) * (COLS), fout);
	if (n < (ROWS) * (COLS) * sizeof(char)) {
		fprintf(stderr, "ERROR: Write output image file output error)\n");
		return 2;
	}

	//-----------------------------------------------------------
	//*****************************Select image part*************
	//-----------------------------------------------------------

	// Store DFT_ori[][] to the output image
//	ofstream myfile("DFTsquare.txt");		// activate only for the 1st image
//	ofstream myfile2("DFTsquareC.txt");		// activate only for the 1st image

	ofstream myfile("DFTcar.txt");			// activate only for the 2nd image
	ofstream myfile2("DFTcarC.txt");		// activate only for the 2nd image
	
	//*************************end of Select part****************
	//-----------------------------------------------------------

	// write DFT_ori[][] to myfile
	if (myfile.is_open())
	{
		for (int countr = 0; countr < ROWS; countr++) 
		{
			for (int countc = 0; countc < COLS; countc++)
			{
				myfile << DFT_ori[countr][countc] << " ";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";
	// writhe DFT to myfile2
	if (myfile2.is_open())
	{
		for (int countr = 0; countr < ROWS; countr++)
		{
			for (int countc = 0; countc < COLS; countc++)
			{
				myfile2 << DFT[countr][countc] << " ";
			}
		}
		myfile2.close();
	}
	else cout << "Unable to open file";

	// Close file open write functions
	fclose(fin);
	fclose(fout);
	return 0;
}