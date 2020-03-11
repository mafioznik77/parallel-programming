#include <mpi.h> 
#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;
#define Main_Process 0
#define Port_Merge 2
#define Port_size_vec 4

double* Create_vec(int _size) {// генерация вектора
	if (_size < 1)
		return NULL;
	double a = 0, b = 0;
	cout << endl << "Enter range [a,b] for array" << endl;
	cin >> a >> b;
	double* _matr_as_vec = new double[_size];
	for (int i = 0; i < _size; i++)
		_matr_as_vec[i] = (long long)(rand() * rand()) % (long long)(b - a + 1) + a;
	return  _matr_as_vec;
}
void Show_vec(double* _vec, int _size) {// вывод массива в поток и в файл
	if (_vec == NULL || _size < 1)
		return;
	int delta = 7;
	if (_size < 20) {
		cout << "Array [" << _size << "]" << endl;
		for (int i = 0; i < _size; i++) {
			cout.width(delta);
			cout << _vec[i];
		}
	}
	cout << endl << "Output matrix to file? (y/n) ";
	char answer = 'n';
	cin >> answer;
	if (answer == 'y') {
		ofstream strm("Array.txt");
		strm << "Array [" << _size << "]" << endl;
		for (int i = 0; i < _size; i++) {
			strm.width(delta);
			strm << _vec[i];
		}
		strm.close();
	}
}
int Equality_test(double* _mas_seque, double* _mas_paral, int _size) {// проверка идентичности результатов
	for (int i = 0; i < _size; i++)
		if (_mas_seque[i] != _mas_paral[i])
			return 0;
	return 1;
}


void Merg(double* mas, double* tmp_mas, int l1, int r1, int l2, int r2) {
	int i = l1, j = l2, k = l1;
	while ((i <= r1) && (j <= r2)) {
		if (mas[i] < mas[j])
			tmp_mas[k++] = mas[i++];
		else
			tmp_mas[k++] = mas[j++];
	}
	if (i > r1)
		for (j; j <= r2; j++)
			tmp_mas[k++] = mas[j];
	else
		for (i; i <= r1; i++)
			tmp_mas[k++] = mas[i];
	for (i = l1; i <= r2; i++)
		mas[i] = tmp_mas[i];
}
void Quick_Sort(double* m, int l, int r)
{
	int i = l, j = r, e, t;
	e = m[(r + l) / 2];
	while (i <= j) {
		while (m[i] < e)
			i++;
		while (m[j] > e)
			j--;
		if (i <= j) {
			if (i < j) {
				t = m[i];
				m[i] = m[j];
				m[j] = t;
			}
			i++;
			j--;
		}
	}
	if (j > l)
		Quick_Sort(m, l, j);
	if (r > i)
		Quick_Sort(m, i, r);
}


int* Reverse_notation(int _n, int _k) { // (n=2^k)
	int higher = -1;
	int* _rev = new int[_n];
	_rev[0] = 0;
	for (int i = 1; i < _n; i++) {
		if ((i & (i - 1)) == 0)
			higher++;
		_rev[i] = _rev[i ^ (1 << higher)];
		_rev[i] |= 1 << (_k - higher - 1);
	}
	return _rev;
}
void Dive(int* _border, int id, int l, int r, int num_dive, int max_dive) {
	if (num_dive == max_dive) {
		_border[id] = l;
		return;
	}
	num_dive++;
	int l1, r1;
	l1 = (l + r) / 2;
	r1 = l1 + 1;
	Dive(_border, id, l, l1, num_dive, max_dive);
	Dive(_border, id + (1 << (max_dive - num_dive)), r1, r, num_dive, max_dive);
}
int* Create_list_work(int _size, int _n, int _k) { // n=2^k
	int* _border = new int[_n];
	Dive(_border, 0, 0, _size - 1, 0, _k);
	return _border;
}
int* Calculation_length_mas(int* _sendcounts, int _k, int _proc_num) {
	int* _size = new int[_proc_num];
	for (int i = 0; i < _proc_num; i++)
		_size[i] = _sendcounts[i];
	int mask;
	for (_k--; _k >= 0; _k--) {
		mask = 1 << _k;
		for (int i = 0; i < mask; i++) {
			if ((i | mask) < _proc_num)
				_size[i] += _size[i | mask];
		}
	}
	return _size;
}


int main(int argc, char* argv[]) {
	int proc_num;										// число процессов
	int proc_rank;										// номер текущего процесса

	double* mas = NULL;									// массив
	int size = 0;										// длина массива

	double* mas_seque = NULL;							// отсортированный последовательной версии алгоритма
	double* mas_paral = NULL;							// отсортированный параллельной версии алгоритма

	double time_seque = 0, time_paral = 0;				// время работы последовательной и параллельной версии алгоритма
	double time_start_seque = 0, time_start_paral = 0;	// время начала работы -//-
	double time_end_seque = 0, time_end_paral = 0;		// время конца работы -//-

	int* displs = NULL;									// массив смещений в векторе для начальной передачи данных
	int* sendcounts = NULL;								// массив длин подвекторов для начальной передачи данных
	int* size_buf = NULL;								// массив длин подвекторов для работы сортировки

	int n = 1, k = 0;	// n=2^k						// параметры бинарного деления массива
	double* part_mas = NULL;							// массив мортируемый каждым процессом
	double* tmp_mas = NULL;								// побочный массив для слияния массивов
	int part_size_vec, part_size_buf;					// размер вектора и его заполненость



	int err_code = 0;
	MPI_Status stat;							// структура атрибутов сообщений для "общения" процессов
	err_code = MPI_Init(&argc, &argv);			// передача всем процессам аргументов командной строки

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);	// функция определения числа процессов
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);	// функция определения номера текущего процесса



	if (proc_rank == Main_Process) {			// код главного процесса

		cout << endl << "Number process: " << proc_num << endl;
		cout << "Enter the size array: ";
		cin >> size;

		// генерация матрицы, вектора, а также представление матрицы в виде вектора
		srand((unsigned)time(NULL));
		mas = Create_vec(size);
		if (mas == NULL) {
			cout << "Error! Incorrect input data for array";
			MPI_Finalize();
			return -1;
		}

		while (n < proc_num) {
			n = n << 1;
			k++;
		}
		int* rev = Reverse_notation(n, k);
		int* border = Create_list_work(size, n, k);

		displs = new int[proc_num];
		sendcounts = new int[proc_num];

		for (int i = 0; i < proc_num; i++)
			displs[i] = border[rev[i]];

		int* tmp = new int[n];
		for (int i = 0; i < n - 1; i++)
			tmp[i] = border[i + 1] - border[i];
		tmp[n - 1] = size - border[n - 1];
		for (int i = 0; i < n; i++) {
			if (rev[i] < proc_num)
				sendcounts[rev[i]] = tmp[i];
			else
				sendcounts[rev[i - 1]] += tmp[i];
		}
		size_buf = Calculation_length_mas(sendcounts, k, proc_num);
		delete[] tmp;
		delete[] rev;
		delete[] border;
	}
	// НАЧАЛО ПАРАЛЛЕЛЬНОГО АЛГОРИТМА

		// отправка каждому процессу массива данных и их кол-во

	MPI_Bcast(&k, 1, MPI_INT, Main_Process, MPI_COMM_WORLD);
	MPI_Scatter(sendcounts, 1, MPI_INT, &part_size_vec, 1, MPI_INT, Main_Process, MPI_COMM_WORLD);

	MPI_Scatter(size_buf, 1, MPI_INT, &part_size_buf, 1, MPI_INT, Main_Process, MPI_COMM_WORLD);
	part_mas = new double[part_size_buf];
	tmp_mas = new double[part_size_buf];
	MPI_Scatterv(mas, sendcounts, displs, MPI_DOUBLE, part_mas, part_size_vec, MPI_DOUBLE, Main_Process, MPI_COMM_WORLD);

	if (proc_rank == Main_Process)
		time_start_paral = MPI_Wtime();

	Quick_Sort(part_mas, 0, part_size_vec - 1);

	int max_rank = 1 << (k - 1);
	int source = proc_rank ^ max_rank;

	if (k > 0) {
		if (proc_rank < max_rank) {
			if (source < proc_num) {
				int add_size;
				MPI_Recv(&add_size, 1, MPI_INT, source, 8, MPI_COMM_WORLD, &stat);
				MPI_Recv(part_mas + part_size_vec, add_size, MPI_DOUBLE, source, Port_Merge, MPI_COMM_WORLD, &stat);
				Merg(part_mas, tmp_mas, 0, part_size_vec - 1, part_size_vec, part_size_vec + add_size - 1);
				part_size_vec += add_size;
			}
		}
		else {
			MPI_Send(&part_size_vec, 1, MPI_INT, source, 8, MPI_COMM_WORLD);
			MPI_Send(part_mas, part_size_vec, MPI_DOUBLE, source, Port_Merge, MPI_COMM_WORLD);
		}
	}
	for (int i = k - 2; (i >= 0) && (proc_rank < max_rank); i--) {
		source = proc_rank ^ (1 << i);
		max_rank = max_rank >> 1;
		if (proc_rank < max_rank) {
			int add_size;
			MPI_Recv(&add_size, 1, MPI_INT, source, Port_size_vec, MPI_COMM_WORLD, &stat);
			MPI_Recv(part_mas + part_size_vec, add_size, MPI_DOUBLE, source, Port_Merge, MPI_COMM_WORLD, &stat);
			Merg(part_mas, tmp_mas, 0, part_size_vec - 1, part_size_vec, part_size_vec + add_size - 1);
			part_size_vec += add_size;
		}
		else {
			MPI_Send(&part_size_vec, 1, MPI_INT, source, Port_size_vec, MPI_COMM_WORLD);
			MPI_Send(part_mas, part_size_vec, MPI_DOUBLE, source, Port_Merge, MPI_COMM_WORLD);
		}
	}

	// КОНЕЦ ПАРАЛЛЕЛЬНОГО АЛГОРИТМА

	if (proc_rank == Main_Process) {
		time_end_paral = MPI_Wtime();
		time_paral = 1000 * (time_end_paral - time_start_paral); // посчет времени работы алгоритма в миллисекудах
		cout << endl << "Spend time algorithm (Parallel version): " << time_paral << "ms";

		// НАЧАЛО ПОСЛЕДОВАТЕЛЬНОГО АЛГОРИТМА
		time_start_seque = MPI_Wtime();
		Quick_Sort(mas, 0, size - 1);
		time_end_seque = MPI_Wtime();

		time_seque = 1000 * (time_end_seque - time_start_seque); // посчет времени работы алгоритма в миллисекудах
		cout << endl << "Spend time algorithm (Sequence version): " << time_seque << "ms" << endl;
		// КОНЕЦ ПОСЛЕДОВАТЕЛЬНОГО АЛГОРИТМА

		mas_paral = part_mas;
		mas_seque = mas;

		// вывод результатов работы алгоритмов
		if (time_seque < time_paral)
			cout << "Sequence version faster parallel" << endl;
		else
			cout << "Parallel version faster sequence" << endl;

		if (Equality_test(mas_seque, mas_paral, size))
			cout << "Result sequence and parallel version identical" << endl;
		else
			cout << "Result sequence and parallel version not identical" << endl;

		// вывод в поток матрицы небольших размеров + вывод мартицы в файл
		Show_vec(mas, size);
		// очистка памяти
		delete[] mas;
		delete[] displs;
		delete[] sendcounts;
		delete[] size_buf;
	}
	delete[] part_mas;
	delete[] tmp_mas;

	MPI_Finalize();
	return 0;
}