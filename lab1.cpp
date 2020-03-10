#include <iostream>
#include <time.h>
#include "mpi.h"

int* GenArr(int& q)
{
	int* arr;
	std::cout << "Enter the number of array elements: " << std::endl;
	std::cin >> q;

	arr = new int[q];

	srand(time(0));

	for (int i = 0; i < q; i++)
		arr[i] = rand();


	if (q < 51)
	{
		for (int i = 0; i < q; i++)
			std::cout << arr[i] << " ";
		std::cout << std::endl;
	}

	return arr;
}

int Min(int* arr, int q)
{
	int min;

	min = arr[0];

	for (int i = 1; i < q; i++)
	{
		if (arr[i] < min)
			min = arr[i];
	}

	return min;
}

int main(int argc, char* argv[])
{
	int n;

	int errCode;

	MPI_Status Status;

	if ((errCode = MPI_Init(&argc, &argv)) != 0)
	{
		return errCode;
	}

	int ProcNum, ProcRank;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0)
	{
		int* arr, q;

		arr = GenArr(q);
		double start, end;
		start = MPI_Wtime();
		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Send(&q, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		n = q / ProcNum;

		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Send(arr + n * (i - 1), n, MPI_INT, i, 0, MPI_COMM_WORLD);
		}


		int min, tmp;

		min = arr[0];

		for (int i = n * (ProcNum - 1) - 1; i < q; i++)
		{
			if (arr[i] < min)
				min = arr[i];
		}

		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Recv(&tmp, 1, MPI_INT, MPI_ANY_SOURCE,
				MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			if (tmp < min)
				min = tmp;
		}
		end = MPI_Wtime();

		std::cout << "Parallel version minimum:" << min << std::endl;
		std::cout << "Parallel version time:" << end - start << std::endl;
		start = MPI_Wtime();
		min = Min(arr, q);
		end = MPI_Wtime();
		std::cout << "Sequence version minimum:" << min << std::endl;
		std::cout << "Sequence version time:" << end - start << std::endl;
		delete[] arr;
	}
	else
	{
		int* buf, q;

		MPI_Recv(&q, 1, MPI_INT, 0,
			0, MPI_COMM_WORLD, &Status);

		n = q / ProcNum;

		buf = new int[n];

		MPI_Recv(buf, n, MPI_INT, 0,
			0, MPI_COMM_WORLD, &Status);

		int min = buf[0];

		for (int i = 1; i < n; i++)
		{
			if (buf[i] < min)
				min = buf[i];
		}

		MPI_Send(&min, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

	}

	MPI_Finalize();



	return 0;
}
