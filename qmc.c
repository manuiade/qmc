#define _CRT_SECURE_NO_WARNINGS
#include "qmc.h"

// Print '_' based on input configuration (0 or 1)
void upperTerm(int x, int pos, int num) {
	if (pos) {
		int z;
		for (z = num - 1; z >= 0; z--) {
			if (pos & (1 << z)) {
				if (x & (1 << z))
					printf("_");
				else
					printf(" ");
			}
		}
	}
}

//Print the letter corresponding to one of the inputs
void lowerTerm(int pos, int num) {
	if (pos) {
		int z;
		for (z = 0; z < num; z++) {
			if (pos & (1 << z)) {
				printf("%c", 'z' - (num - 1) + z);
			}
		}
	}
}

//Print all possible inputs
void outputTerm(int x, int pos, int num) {
	upperTerm(x, pos, num);
	if (pos) printf("\n");
	lowerTerm(pos, num);
}


void swap(int i, int j, int k,int num) {
	int tmp[ING];
	int tmp2;
	int tmp3;
	if (!k) {
		for (int x = 0; x < num; ++x) {
			tmp[x] = tab[i][x]; tab[i][x] = tab[j][x]; tab[j][x] = tmp[x];
		}
		tmp2 = min_terms[i]; min_terms[i] = min_terms[j]; min_terms[j] = tmp2;
		tmp3 = not_ones_1[i]; not_ones_1[i] = not_ones_1[j]; not_ones_1[j] = tmp3;
	}
	else {
		for (int x = 0; x < num; ++x) {
			tmp[x] = tabind[i][x]; tabind[i][x] = tabind[j][x]; tabind[j][x] = tmp[x];
		}
		tmp2 = indiff[i]; indiff[i] = indiff[j]; indiff[j] = tmp2;
		tmp3 = not_ones_2[i]; not_ones_2[i] = not_ones_2[j]; not_ones_2[j] = tmp3;
	}
}

// Order min_terms based on their occurencies of 1
void sort(int nmin, int nind,int num) {
	for (int i = 0; i < nmin - 1; ++i) {
		for (int j = i + 1; j < nmin; ++j) {
			if (not_ones_1[i] > not_ones_1[j])
				swap(i, j, 0,num);
		}
	}
	for (int i = 0; i < nind - 1; ++i) {
		for (int j = i + 1; j < nind; ++j) {
			if (not_ones_2[i] > not_ones_2[j])
				swap(i, j, 1,num);
		}
	}
}

// Apply combination property based on hamming score
void fillRange(int a, int b, int step, int *dimint) {
	for (int x = 0; x < (1 << (step - 1)); ++x)
		intmp[*dimint][x] = range[a][x];
	for (int x = (1 << (step - 1)), y = 0; x < (1 << step); ++x, ++y)
		intmp[*dimint][x] = range[b][y];
	++(*dimint);
}

// Calculate hemming distance
int hemming(int*tab1, int*tab2, int x, int a, int b, int step, int num, int *dimint) {
	int h = 0;
	int tmp[ING];
	for (int i = 0; i < num && h < 2; ++i) {
		if (tab1[i] + tab2[i] == 0)
			tmp[i] = 0;
		else {
			if (tab1[i] + tab2[i] == 2)
				tmp[i] = 1;
			else {
				if (tab1[i] == '-' && tab2[i] == '-')
					tmp[i] = '-';
				else {
					tmp[i] = '-';
					++h;
				}
			}
		}
	}
	if (h == 1) {
		if (flag[a] != -1 && flag[b] != -1) {
			flag[a] = 0;
			flag[b] = 0;
		}
		int i = 0;
		int check = 0;
		for (int z = 0; z < x; ++z) {
			for (i = 0; i < num; ++i) {
				if (tmp[i] == reduction[z][i])
					++check;
			}
			if (check == num)
				return 0;
			check = 0;
		}
		for (i = 0; i < num; ++i) {
			reduction[x][i] = tmp[i];
		}
		fillRange(a, b, step, dimint);
		return 1;
	}
	return 0;
}

// Copy non combined minterms
int copy_minterms(int *x, int nmin,int num) {
	for (int a = 0; a < nmin; ++a) {
		if (flag[a] == 1 || flag[a] == -1) {
			for (int y = 0; y < num; ++y) {
				reduction[*x][y] = tab[a][y];
				intmp[*x][y] = range[a][y];
			}
			++(*x);
		}
	}
	return 0;
}

int fill_matrix(int x, int nmin, int step) {
	int c; int check = 0; int ip = 0;
	for (int a = 0; a < x; ++a) {
		check = 0;
		for (int b = 0; b < (1 << step); ++b) {
			for (c = 0; c < nmin; ++c) {
				if (intmp[a][b] == 0 && b >0)
					break;
				if (intmp[a][b] == mincolumn[c]) {
					matrbin[a][c] = 1;
					++c;
					check = 1;
					break;
				}
			}
		}
		if (check)
			++ip;
	}
	return ip;
}

// Branch and bound recursive algorithm to calculate essential primes
void bab(int s, int nmin, int x, int best, int*curr) {
	int var = 0;
	for (int i = 0; i < nmin; ++i) {
		if (check[i] > 0)
			++var;
	}
	if (var == nmin) {
		if (best < (*curr)) {
			(*curr) = best;
			for (int i = 0; i < x; ++i) {
				if (prime_ess[i] == 1)
					solution[i] = 1;
				else
					solution[i] = 0;
			}
		}
		return;
	}
	if (s == x)
		return;
	for (int i = s; i < x; ++i) {
		int ess = 0;
		for (int j = 0; j < nmin; ++j) {
			if (matrbin[i][j] == 1) {
				check[j] += matrbin[i][j];
				prime_ess[i] = 1;
				ess = 1;
			}
		}
		if (ess)
			++best;
		s = i;
		bab(s + 1, nmin, x, best, curr);
		prime_ess[i] = 0;
		for (int j = 0; j < nmin; ++j) {
			if (matrbin[i][j] == 1)
				check[j]--;
		}
		best--;
	}
}

// Print result in SP format
void printsp(int x,int num) {
	for (int a = 0; a < x; ++a) {
		if (solution[a] == 1) {
			for (int b = 0; b < num; ++b) {
				if (reduction[a][b] == 0)
					printf("_");
				if (reduction[a][b] == 1)
					printf(" ");
			}
			printf("   ");
		}
	}
	printf("\n");
	int tmp = 0;
	for (int a = 0; a < x; ++a) {
		if (solution[a] == 1) {
			if (tmp) {
				printf(" + ");
			}
			tmp = 1;
			for (int b = 0; b < num; ++b) {
				if (reduction[a][b] == 45)
					continue;
				printf("%c", 'z' - (num - 1) + b);
			}
		}
	}
}

// Print results in PS format
void printps(int x, int num) {
	printf(" ");
	for (int a = 0; a < x; ++a) {
		if (solution[a] == 1) {
			for (int b = 0; b < num; ++b) {
				if (reduction[a][b] == 0)
					printf("_");
				if (reduction[a][b] == 1)
					printf(" ");
				if (reduction[a][b] == 45)
					continue;
				printf(" ");
			}
			printf(" ");
		}
	}
	printf("\n");
	for (int a = 0; a < x; ++a) {
		int tmp = 0;
		if (solution[a] == 1) {
			printf("(");
			for (int b = 0; b < num; ++b) {
				if (reduction[a][b] == 45)
					continue;
				if (tmp)
					printf("+");
				tmp = 1;
				printf("%c", 'z' - (num - 1) + b);
			}
			printf(")");
		}
	}
}


void printrid(int x, int passo,int num) {
	printf("Reduction Number %d:\n", passo);
	for (int a = 0; a < x; ++a) {
		int b = 0;
		int tmp = 0;
		for (; b < num; ++b) {
			if (reduction[a][b] == 45)
				printf("-");
			else
				printf("%d", reduction[a][b]);
		}
		printf("\t\t\t { ");
		b = 0;
		while (b == 0 || intmp[a][b] != 0)
		{
			if (tmp)
				printf(",");
			tmp = 1;
			printf("%d", intmp[a][b]);
			++b;
		}
		printf(" }\n\n");
	}
}