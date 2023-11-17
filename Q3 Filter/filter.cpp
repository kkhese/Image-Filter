/*
Term Project 03		====================================
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

// Image size
#define ROWS 256
#define COLS 256

// ****************************Set cut-off frequency****************************
// Radius for cutoff frequency of LPF, HPF
#define R 5

// Buffer for input images
unsigned char in_img[ROWS][COLS];

// Buffer for saving FFT values from input file
complex<double> DFT[ROWS][COLS];

// Buffer for FFT after filter
complex<double> LFT[ROWS][COLS];	// DFT after LPF
complex<double> HFT[ROWS][COLS];	// DFT after HPF

// Temporary Buffer for shuffling array
complex<double> ODFT[ROWS][COLS];	// Centered/Decentered DFT
complex<double> RDFT[ROWS][COLS];	// After row index shuffle
complex<double> CDFT[ROWS][COLS];	// After col index Shuffle

// Buffer for FFT images
complex<double> DFT1[ROWS][COLS];	// temp buffer for 1D IFFT
complex<double> DFT2[ROWS][COLS];	// temp buffer for 2D IFFT

// Buffer for output images
unsigned char out_img[ROWS][COLS];

// Definition of pi
double pi = 2 * asin(1);

// number of elements in a row/column
int N = ROWS;

//--------------------------------------------------------------------------------
//-------------------- Sub-function: Fiters --------------------------------------
//--------------------------------------------------------------------------------

//==========================================================
// Ideal Lowpass filter 
//==========================================================
// Input is DFT[][] => result will be saved in LFT[][] array
void LPF(int radius)
{	// Scanning entire pixels
	for (int row = 0; row < N; row++)
	{
		for (int col = 0; col < N; col++)
		{
			// For quadrant I
			if (row <= 127 && col >= 128)
			{
				// pixels are inside of LPF radius, keep values
				if ( (radius*radius) >= ((128 - row)*(128 - row) + (col - 127)*(col - 127)) )
				{
					LFT[row][col] = ODFT[row][col];
				}
				// pixels are outside of LPF radius, DFT values becomes 0
				else
				{
					LFT[row][col] = (0.0,0.0);
				}
			}
			// for quadrant II
			else if (row <= 127 && col <= 127)
			{
				// pixels are inside of LPF radius, keep values
				if ( (radius*radius) >= ((128 - row)*(128 - row) + (128 - col)*(128 - col)) )
				{
					LFT[row][col] = ODFT[row][col];
				}
				// pixels are outside of LPF radius, DFT values becomes 0
				else
				{
					LFT[row][col] = (0.0, 0.0);
				}
			}
			// for quadrant III
			else if (row >= 128 && col <= 127)
			{
				// pixels are inside of LPF radius, keep values
				if ((radius*radius) >= ((row - 127)*(row - 127) + (128 - col)*(128 - col)))
				{
					LFT[row][col] = ODFT[row][col];
				}
				// pixels are outside of LPF radius, DFT values becomes 0
				else
				{
					LFT[row][col] = (0.0, 0.0);
				}
			}
			// for quadrant IV
			else if (row >= 128 && col >= 128)
			{
				// pixels are inside of LPF radius, keep values
				if ((radius*radius) >= ((row - 127)*(row - 127) + (col-127)*(col-127)))
				{
					LFT[row][col] = ODFT[row][col];
				}
				// pixels are outside of LPF radius, DFT values becomes 0
				else
				{
					LFT[row][col] = (0.0, 0.0);
				}
			}
		}
	}
}

//==========================================================
// Gaussian Lowpass filter 
//==========================================================
// Input is DFT[][] => result will be saved in LFT[][] array
void GLPF(int radius)
{	// Scanning the entire pixels
	for (int row = 0; row < N; row++)
	{
		for (int col = 0; col < N; col++)
		{
			// For quadrant I
			if (row <= 127 && col >= 128)
			{
				LFT[row][col] = exp(-((128 - row)*(128 - row) + (col - 127)*(col - 127))/(2*radius*radius))*ODFT[row][col];
			}
			// for quadrant II
			else if (row <= 127 && col <= 127)
			{
				LFT[row][col] = exp(-((128 - row)*(128 - row) + (128 - col)*(128 - col)) / (2 * radius*radius))*ODFT[row][col];
			}
			// for quadrant III
			else if (row >= 128 && col <= 127)
			{
				LFT[row][col] = exp(-((row-127)*(row-127) + (128 - col)*(128 - col)) / (2 * radius*radius))*ODFT[row][col];
			}
			// for quadrant IV
			else if (row >= 128 && col >= 128)
			{
				LFT[row][col] = exp(-((row - 127)*(row - 127) + (col-127)*(col-127)) / (2 * radius*radius))*ODFT[row][col];
			}
		}
	}
}

//==========================================================
// Butterworth LPF
//==========================================================
// Input is DFT[][] => result will be saved in LFT[][] array
void BLPF(int radius)
{
	int n = 10;			// Filtering factor for butterworth filter
	// scanning the entire pixels
	for (int row = 0; row < N; row++)
	{
		for (int col = 0; col < N; col++)
		{
			// For quadrant I
			if (row <= 127 && col >= 128)
			{
				LFT[row][col] = 1 / ( 1 + pow( (((128 - row)*(128 - row) + (col - 127)*(col - 127)) / (radius*radius)),n ))*ODFT[row][col];
			}
			// for quadrant II
			else if (row <= 127 && col <= 127)
			{
				LFT[row][col] = 1 / (1 + pow((((128 - row)*(128 - row) + (128 - col)*(128 - col)) / (radius*radius)), n))*ODFT[row][col];
			}
			// for quadrant III
			else if (row >= 128 && col <= 127)
			{
				LFT[row][col] = 1 / (1 + pow((((row-127)*(row-127) + (128 - col)*(128 - col)) / (radius*radius)), n))*ODFT[row][col];
			}
			// for quadrant IV
			else if (row >= 128 && col >= 128)
			{
				LFT[row][col] = 1 / (1 + pow((((row - 127)*(row - 127) + (col-127)*(col-127)) / (radius*radius)), n))*ODFT[row][col];
			}
		}
	}
}

//==========================================================
// Ideal high pass filter(HPF)
//==========================================================
// Input is DFT[][] => result will be saved in HFT[][] array
void HPF(int radius)
{	// scanning the entire pixels
	for (int row = 0; row < N; row++)
	{
		for (int col = 0; col < N; col++)
		{
			// For quadrant I
			if (row <= 127 && col >= 128)
			{
				// If pixels are outside of HPF radius, keep values
				if ( (radius*radius) < ((128 - row)*(128 - row) + (col - 127)*(col - 127)) )
				{
					HFT[row][col] = ODFT[row][col];
				}
				// If pixels are inside of HPF radius, DFT values becomes 0
				else
				{
					HFT[row][col] = (0.0,0.0);
				}
			}
			// for quadrant II
			else if (row <= 127 && col <= 127)
			{
				// If pixels are outside of HPF radius, keep values
				if ( (radius*radius) < ((128 - row)*(128 - row) + (128 - col)*(128 - col)) )
				{
					HFT[row][col] = ODFT[row][col];
				}
				// If pixels are inside of HPF radius, DFT values becomes 0
				else
				{
					HFT[row][col] = (0.0, 0.0);
				}
			}
			// for quadrant III
			else if (row >= 128 && col <= 127)
			{
				// pixels are outside of HPF radius, keep values
				if ((radius*radius) < ((row - 127)*(row - 127) + (128 - col)*(128 - col)))
				{
					HFT[row][col] = ODFT[row][col];
				}
				// pixels are inside of HPF radius, DFT values becomes 0
				else
				{
					HFT[row][col] = (0.0, 0.0);
				}
			}
			// for quadrant IV
			else if (row >= 128 && col >= 128)
			{
				// pixels are outside of HPF radius, keep values
				if ((radius*radius) < ((row - 127)*(row - 127) + (col-127)*(col-127)))
				{
					HFT[row][col] = ODFT[row][col];
				}
				// pixels are inside of HPF radius, DFT values becomes 0
				else
				{
					HFT[row][col] = (0.0, 0.0);
				}
			}
		}
	}
}

//==========================================================
// Gaussian high pass filter(GHPF)
//==========================================================
// Input is DFT[][] => result will be saved in HFT[][] array
void GHPF(int radius)
{	// scanning the entire pixels
	for (int row = 0; row < N; row++)
	{
		for (int col = 0; col < N; col++)
		{
			// For quadrant I
			if (row <= 127 && col >= 128)
			{
				HFT[row][col] = (1-exp(-((128 - row)*(128 - row) + (col - 127)*(col - 127)) / (2 * radius*radius)))*ODFT[row][col];
			}
			// for quadrant II
			else if (row <= 127 && col <= 127)
			{
				HFT[row][col] = (1-exp(-((128 - row)*(128 - row) + (128 - col)*(128 - col)) / (2 * radius*radius)))*ODFT[row][col];
			}
			// for quadrant III
			else if (row >= 128 && col <= 127)
			{
				HFT[row][col] = (1-exp(-((row - 127)*(row - 127) + (128 - col)*(128 - col)) / (2 * radius*radius)))*ODFT[row][col];
			}
			// for quadrant IV
			else if (row >= 128 && col >= 128)
			{
				HFT[row][col] = (1-exp(-((row - 127)*(row - 127) + (col - 127)*(col - 127)) / (2 * radius*radius)))*ODFT[row][col];
			}
		}
	}
}

//==========================================================
// // Butterworth HPF
//==========================================================
// Input is DFT[][] => result will be saved in HFT[][] array
void BHPF(int radius)
{
	int n = 10;			// filter factor 
	// Scanning the entire pixels
	for (int row = 0; row < N; row++)
	{
		for (int col = 0; col < N; col++)
		{
			// For quadrant I
			if (row <= 127 && col >= 128)
			{
				HFT[row][col] = 1 / (1 + pow(( (radius*radius)/ ((128 - row)*(128 - row) + (col - 127)*(col - 127))), n))*ODFT[row][col];
			}
			// for quadrant II
			else if (row <= 127 && col <= 127)
			{
				HFT[row][col] = 1 / (1 + pow(( (radius*radius)/ ((128 - row)*(128 - row) + (128 - col)*(128 - col))), n))*ODFT[row][col];
			}
			// for quadrant III
			else if (row >= 128 && col <= 127)
			{
				HFT[row][col] = 1 / (1 + pow(( (radius*radius)/ ((row - 127)*(row - 127) + (128 - col)*(128 - col))), n))*ODFT[row][col];
			}
			// for quadrant IV
			else if (row >= 128 && col >= 128)
			{
				HFT[row][col] = 1 / (1 + pow(( (radius*radius)/ ((row - 127)*(row - 127) + (col - 127)*(col - 127))), n))*ODFT[row][col];
			}
		}
	}
}

//--------------------------------------------------------------------------------
//-------------------- Sub-function: End of Fiters -------------------------------
//--------------------------------------------------------------------------------


//------------------------------------------------------------------------
//--------------------- the rest of FFT sub functions --------------------
//------------------------------------------------------------------------

// Center DFT values::: output ==> ODFT[][]===============================
void Centering(complex<double> Comimg[ROWS][COLS])
{
	for (int i = 0; i < ROWS; i++)
	{
		for (int j = 0; j < COLS; j++)
		{	// Center the first DFT result and resorting the array
			ODFT[(i + (N/2)) % N][(j + (N/2)) % N] = Comimg[i][j];
		}
	}

}

// Index shuffle for each Row ::: Output => RDFT[][]=====================
// The same function in FFT codes.=======================================
void Rshuffle(complex<double> Comimg[ROWS][COLS])
{
	string bin;				// string for binary number
	int indextemp = 0;		// initilizae shuffling index to 0
	// for loop for entire rows
	for (int i = 0; i < ROWS; i++)
	{	// checking each Column index
		for (int j = 0; j < COLS; j++)
		{
			// convert integer j(column index) to binary string 
			// and save log2(N)=8 bit binary to bin
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
			// Element at Re - esimated, shuffled column index
			// Goes to element of current scanning index
			RDFT[i][j] = Comimg[i][indextemp];
			indextemp = 0;		// Reset indextemp for next column index
		}
	}
}

// Index shuffle for each Row ::: Output => CDFT[][]=====================
// The same function in FFT codes.=======================================
void Cshuffle(complex<double> Comimg[ROWS][COLS])
{
	string bin;					// string for binary number
	int indextemp = 0;			// initilizae shuffling index to 0
	// for loop for entire columns
	for (int i = 0; i < COLS; i++)
	{	// checking each row index
		for (int j = (ROWS - 1); j >= 0; j--)
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
			CDFT[j][i] = Comimg[indextemp][i];
			indextemp = 0;		// Reset indextemp for next column index
		}
	}
}

// function for calculating W_n^u for inverse FFT ========================
// return complex<double>=================================================
complex<double> invW(int M, int u)
{
	// real part: cos(pi*u/M)
	// Imaginary part: sin(pi*u/M)
	double real = cos((pi * u) / M);
	double imaginary = sin((pi * u) / M);
	// if real and imaginary values are smaller than 0.0001
	// Set those values to 0 
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

// 1D inverse FFT function for single row ===============================
// similar with FFT =====================================================
void IFFT1D(int z, complex<double> Comimg[ROWS][COLS])
{
	double n = 0.0;			// initialze integer n
	n = log2(N);			// n = Number of merging stages
							// N is array size for a single row
	int M = 1;				// initial subgroup length
	int pair = N / 2;		// initial number of pairs of subgroup

	// IFFT Calculation 
	for (int i = 0; i < n; i++)				// Conducting merging n times.
	{
		for (int k = 0; k < pair; k++)		// Merge pairs at the current level
		{
			int start1 = k * 2 * M;			// start index for first group
			int start2 = (k * 2 + 1) * M;	// start index for second group
			// Fast inverse fourier transfor form
			for (int u = 0; u < M; u++)
			{
				// Temporarily save intermediate IFFT values to DFT1
				DFT1[z][u] = (Comimg[z][start1 + u] + Comimg[z][start2 + u] * invW(M, u));
				DFT1[z][u + M] = (Comimg[z][start1 + u] - Comimg[z][start2 + u] * invW(M, u));
			}
			// update the IFFT rusult at input Aarray
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

// 2D inverse FFT function for single column ============================
// similar with 1D inverse FFT ==========================================
void IFFT2D(int z, complex<double> Comimg[ROWS][COLS])
{
	double n = 0.0;			// initialze integer n
	n = log2(N);			// Number of merge stages
							// N is array size for a single row
	int M = 1;				// initial subgroup length
	int pair = N / 2;		// initial number of pairs of subgroup

	// IFFT Calculation 
	for (int i = 0; i < n; i++)				// Conducting merging n times.
	{
		for (int k = 0; k < pair; k++)		// Merge pairs at the current level
		{
			int start1 = k * 2 * M;			// start index for first group
			int start2 = (k * 2 + 1) * M;	// start index for second group
			// fast inverse fourier transfor form
			for (int u = 0; u < M; u++)
			{
				// Temporarily save intermediate IFFT values to DFT2 
				DFT2[u][z] = (Comimg[start1 + u][z] + Comimg[start2 + u][z] * invW(M, u));
				DFT2[u + M][z] = (Comimg[start1 + u][z] - Comimg[start2 + u][z] * invW(M, u));
			}
			// update the IFFT rusult at input Aarray
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

// Generate output with IFFT result =====================================
// output ==> (unsigned char) out_img[][] ===============================
void Filteredout(complex<double> Comimg[ROWS][COLS])
{
	// initial max/min value
	double max = 0.0;
	double min = 2000.0;
	// update abs max/min value by scanning input arrays
	for (int i = 0; i < ROWS; i++)
	{
		for (int j = 0; j < COLS; j++)
		{
			if (max < abs(Comimg[i][j]))
			{
				max = abs(Comimg[i][j]);
			}
			if (min > abs(Comimg[i][j]))
			{
				min = abs(Comimg[i][j]);
				// if min value is smaller than 0.0001, set min to 0
				if (abs(Comimg[i][j])<0.0001)
				{
					min = 0.0;
				}
			}
		}
	}
	// check out max and min
	cout << max << " " << min << " " << (255.0 / (max - min)) << endl;

	// set overflow pixels' level to 255
	for (int i = 0; i < ROWS; i++)
	{
		for (int j = 0; j < COLS; j++)
		{
			// If intensity level is greater than 255 after IFFT
			// Set those to 255
			if (abs(Comimg[i][j]) > 255)
			{
				out_img[i][j] = (unsigned char)(int)(255);
			}
			// If no overflow, send IFFT results to unsigned char out_img[][]
			else
			{
				out_img[i][j] = (unsigned char)(int)(abs(Comimg[i][j]));
			}
		}
	}
}

//****************************************************************************
//****************************************************************************
//------------------------------------Main starts-----------------------------
//****************************************************************************
//****************************************************************************

int main()
{
//----------------------------------------------------------------------------
//--------------------------------Reading File & create output file---------//

	//*************************************************************
	//*********************** select part**************************
	//*************************************************************
	// Reading DFT data from .txt files => DFT[][]
//	string inFileName = "DFTsquare.txt";			// activate only for 1st image
	string inFileName = "DFTcar.txt";				// activete only for 2nd image
	//*************************************************************
	//*********************** End of select part*******************
	//*************************************************************
	ifstream inFile;
	inFile.open(inFileName.c_str());
	if (inFile.is_open())
	{
		for (int row = 0; row < N; row++)
		{
			for (int col = 0; col < N; col++)
			{	// input file to DFT[][]
				inFile >> DFT[row][col];
			}
		}
		inFile.close();
	}
	else
	{
		cerr << "Can't find input file" << inFileName << endl;
	}
	FILE* fout;
	int n = 0;

	//*************************************************************
	//*********************** select part**************************
	//*************************************************************
	// Open the output image file
//	if (fopen_s(&fout, "squareBHPF.raw", "wb") != 0) {	// Enable this line for 1st image
	if (fopen_s(&fout, "carBLPF.raw", "wb") != 0) {		// Enable this line for 2nd image
	//*************************************************************
	//*********************** End of select part*******************
	//*************************************************************
		fprintf(stderr, "ERROR: Cann't open output image file \n");
	}
//------------------------ end of file management---------------------------//
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
//--------------------------- Applying filters -----------------------------//

	// Centering DFT values: DFT => ODFT
	Centering(DFT);

	//*************************************************************
	//*********************** select part**************************
	//*************************************************************

	// Apply LPF to DFT array and shuffling column index -- Activate next 2 out of 4 lines for LPF
//	LPF(R);				// Activate only for Ideal LPF, R is radius for filter: ODFT => LFT
//	GLPF(R);			// Activate only for Ideal GLPF, R is radius for filter: ODFT => LFT
	BLPF(R);			// Activate only for Ideal BLPF, R is radius for filter: ODFT => LFT
	Centering(LFT);		// Decentering LFT to original shape: LFT => ODFT

	// Apply HPF to DFT array and shuffling column index -- Activate next 2 out of 4 lines for HPF
//	HPF(R);				// Activate only for Ideal HPF, R is radius for filter: ODFT => HFT
//	GHPF(R);			// Activate only for Ideal GHPF, R is radius for filter: ODFT => HFT
//	BHPF(R);			// Activate only for Ideal BHPF, R is radius for filter: ODFT => HFT
//	Centering(HFT);		// Decentering HFT to original shape HFT => ODFT
	//*************************************************************
	//*********************** End of select part*******************
	//*************************************************************


//------------------------ end of filter------------------------------------//
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
//------------------------------- Inverse FFT --------------------------------
	// Shuffle column index order
	Rshuffle(ODFT);		// ODFT => RDFT	

	// IFFT 1D - each row::: RDFT => RDFT
	for (int z = 0; z < ROWS; z++)
	{
		IFFT1D(z, RDFT);
	}
	
	// Shuffling row index order
	Cshuffle(RDFT);		// RDFT => CDFT

	// IFFT 2D - each column::: CDFT => CDFT
	for (int z = 0; z < COLS; z++)
	{
		IFFT2D(z, CDFT);
	}

	// Send results to output::: CDFT => out_img[][]
	Filteredout(CDFT);
//------------------------ end of inverse FFT -----------------------------//
//----------------------------------------------------------------------------

	// Save the output image.
	printf("... Save the output image\n");
	n = fwrite(out_img, sizeof(char), (ROWS) * (COLS), fout);
	if (n < (ROWS) * (COLS) * sizeof(char)) {
		fprintf(stderr, "ERROR: Write output image file output error)\n");
		return 2;
	}
	fclose(fout);
	return 0;
}