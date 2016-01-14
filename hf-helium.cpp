/*
	Author: Giuseppe Accaputo
	Date:	14.01.2016

	Build:	g++ -o hf-helium hf-helium.cpp -llapacke
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <lapacke.h>
#include <math.h>

using namespace std;

/*
	Returns a matrix index using the row-major order
*/
int get_index(int N, int i, int j){
	return (i*N + j);
}

void print_matrix(int N, double* A)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			cout << A[get_index(N,i,j)] << "\t";
		}
		cout << endl;
	}
}

void print_vector(int n, double* v)
{
	for (int i = 0; i < n; ++i)
	{
		cout << v[i] << "\t";
	}
}

/*
	Calculate an entry of the non-interacting matrix T
*/
double calc_Tij(double* alpha, int i, int j){
	return 3.0 * (alpha[i] * alpha[j] * pow(M_PI, 1.5)) / (pow(alpha[i] + alpha[j], 2.5)) - 4.0*M_PI/(alpha[i] + alpha[j]);
}

/*	
	Calculate an entry of the Hartree matrix V
*/
double calc_Vijkl(double* alpha, int i, int j, int k, int l){
	return 2.0 * pow(M_PI, 2.5) / ((alpha[i] + alpha[j]) * (alpha[k] + alpha[l]) * sqrt(alpha[i] + alpha[j] + alpha[k] + alpha[l]));
}

/*
	Calculate an entry of the overlap matrix S
*/
double calc_Sij(double* alpha, int i, int j){
	return pow(M_PI / (alpha[i] + alpha[j]), 1.5);
}

/*
	Update the Fock matrix F
*/
void update_fock_matrix(int N, double* F, double* d, double* alpha){
	double sum_ddVijkl = 0.0f;

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{		
			sum_ddVijkl = 0;
			
			for (int k = 0; k < N; ++k)
			{
				for (int l = 0; l < N; ++l)
				{
					sum_ddVijkl += d[k] * d[l] * calc_Vijkl(alpha,i,j,k,l);
				}
			}

			F[get_index(N,i,j)] = calc_Tij(alpha,i,j) + sum_ddVijkl;
		}
	}
}

void normalize(int N, double* ev, double* alpha){
	double sum = 0.0f;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			sum += ev[i] * calc_Sij(alpha,i,j) * ev[j];
		}
	}

	double sqrt_sum = sqrt(sum);

	for (int i = 0; i < N; ++i)
	{
		ev[i] /= sqrt_sum;
	}
}

int main(int argc, char** argv)
{

	cout.precision(10);

	/*
		Size of the basis set
	*/
	int N = 4;

	/*
		Fock operator matrix
	*/
	double* F = (double*)malloc(N * N * sizeof(double));

	/*
		Overlap matrix
	*/
	double* S = (double*)malloc(N * N * sizeof(double));
	double* S_backup = (double*)malloc(N * N * sizeof(double));
	
	double* alpha = (double*)malloc(N * sizeof(double));
	
	alpha[0] = 0.297104;
	alpha[1] = 1.236745;
	alpha[2] = 5.749982;
	alpha[3] = 38.216677;

	double* d = (double*)malloc(N * sizeof(double));
	for (int i = 0; i < N; ++i)
	{
		d[i] = 1;
	}

	normalize(N, d, alpha);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			S[get_index(N, i, j)] = calc_Sij(alpha, i, j);
			S_backup[get_index(N,i,j)] = S[get_index(N, i, j)];
		}
	}

	int		itype = 1;
	char	jobz = 'V';
	char	uplo = 'U';
	int		ev_size = 3 * N - 1;
	double* eigenvalues = (double*)malloc(ev_size * sizeof(double));

	double tol = 1e-10;
	double eps = 0;
	double eps_old = 2;

	while(fabs(eps - eps_old) > tol){

		/*
			The second matrix argument (in our case the overlap matrix S)
			gets overwritten by the LAPACKE_dsygv routine, so we have to restore it after each step
		*/
		for (int i = 0; i < N * N; ++i)
		{
			S[i] = S_backup[i];
		}

		update_fock_matrix(N, F, d, alpha);

		LAPACKE_dsygv(LAPACK_ROW_MAJOR, itype, jobz, uplo, N, F, N, S, N, eigenvalues);

		/*
			LAPACKE_dsygv stores the eigenvalues in ascending order in the provided eigenvalues array,
			so the eigenvector appertaining to the smallest eigenvalue is stored in the first column
			of the eigenvector matrix Z (F gets overwritten by Z by  the LAPACKE_dsygv call)
		*/
		for (int i = 0; i < N; ++i)
		{
			d[i] = F[get_index(N, i, 0)];
		}

		normalize(N, d, alpha);

		eps_old = eps;
		eps = eigenvalues[0];
	}

	double e_0_1 = 0.0f;
	double e_0_2 = 0.0f;
	double di_dj = 0.0f;

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{	
			di_dj = d[i] * d[j];
			e_0_1 +=  di_dj * calc_Tij(alpha, i, j);

			for (int k = 0; k < N; ++k)
			{
				for (int l = 0; l < N; ++l)
				{
					e_0_2 += di_dj * d[k] * d[l] * calc_Vijkl(alpha, i, j, k, l);
				}
			}
		}
	}

	e_0_1 *= 2.0;

	cout << "The ground state energy E_0 [Hartree] is " << (e_0_1 + e_0_2) << endl;

	return EXIT_SUCCESS;
}
