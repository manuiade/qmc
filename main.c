//ALGORITMO DI QUINE-MCCLUSKEY
//Realizzato da Iaderosa Manuel
#include "qmc.h"


int main(void) {
	int num = 0;
	int pos = 0;
	int x = 0;
	int y = 0;
	int nmin = 0;
	int nind = 0;
	int dh = 0;
	int control = TRUE;
	int step = 0;

	while (1) {
		printf("Choose number of input variables (from 1 to 6 included) -> ");
		scanf("%d", &num);
		if (num > 0 && num <= ING)
			break;
		else
			printf("Error! Not valid value..\n");
	}

	//Set number of possible configurations based on input variables
	pos = (1 << num); 
	int sint = 0;
	char*v = NULL;

	printf("Choose type of minimization to perform :\n1 -> SP\n2 -> PS\n");
	scanf("%d", &sint);
	while (sint != 1 && sint != 2) {
		printf("Not valid choice!!\nChoose type of minimization to perform :\n1 -> SP\n2 -> PS\n");
		scanf("%d", &sint);
	}

	printf("\nInsert results of function to minimize ( 0 , 1 or -1 for indiff) \n\n");

	//Insert of outputs
	for (x = 0; x < pos; x++) {
		int value = 0;
		outputTerm(x, pos - 1, num); //Print all possible inputs
		printf("   (%d)   = ", x);
		scanf(" %d", &value);
		while (value != 0 && value != 1 && value != -1) {
			printf("Not valid result. Insert valid result (0 , 1 or -1 for indiff) ->  ");
			scanf("%d", &value);
		}

		//Minterms are tabulated in binary form, a decimal value and a flag set to 1
		if (sint == 1) {
			if (value == 1) {
				for (int i = 0; i < num; i++) {
					tab[nmin][num - i - 1] = (x >> i) & 1;
					if (tab[nmin][num - i - 1])
						++(not_ones_1[nmin]);
				}
				min_terms[nmin] = x;
				mincolumn[nmin] = x;
				flag[nmin] = 1;
				nmin++;
			}
			printf("\n");
		}
		else
			if (!value) {
				for (int i = 0; i < num; i++) {
					tab[nmin][num - i - 1] = (x >> i) & 1;
					if (tab[nmin][num - i - 1])
						++(not_ones_1[nmin]);
				}
				min_terms[nmin] = x;
				mincolumn[nmin] = x;
				flag[nmin] = 1;
				nmin++;
			}
		if (value == -1) {
			for (int i = 0; i < num; i++) {
				tabind[nind][num - i - 1] = (x >> i) & 1;
				if (tabind[nind][num - i - 1])
					++(not_ones_2[nind]);
			}
			indiff[nind] = x;
			nind++;
		}
	}

	int mintmp = nmin + nind;
	//Sort minterms based on their 1 occurencies
	sort(nmin, nind,num);

	// indiff are put below min_terms table with a flag set to -1
	for (int i = nmin; i < mintmp; ++i)
		flag[i] = -1;
	for (int k = 0; k < nmin; ++k)
		range[k][0] = min_terms[k];
	for (int k = nmin; k < mintmp; ++k) {
		range[k][0] = indiff[k - nmin];
		not_ones_1[k] = not_ones_2[k - nmin];
		for (int l = 0; l < num; ++l)
			tab[k][l] = tabind[k - nmin][l];
	}

	int totmin = nmin;
	clock_t t = clock();
	if (!totmin) {
		printf("No exit available with the following configuration!");
		return 0;
	}
	if (totmin == (1<<num)) {
		printf("All exit available with the following configuration!");
		return 0;
	}

	//Step 1: Prime implicants
	int dimint = 0;

	while (control) {
		step++;
		x = 0;
		y = 0;
		// Hemming only for partition with adjacency
		for (int a = 0; a < nmin; ++a) {
			for (int b = a + 1; b < nmin; ++b) {
				if (pos) {
					if (((not_ones_1[b]) - (not_ones_1[a]) > 1))
						break;
					if (((not_ones_1[b]) - (not_ones_1[a]) == 0))
						continue;
				}
				// Hemming distance calculation
				dh = hemming(tab[a], tab[b], x, a, b, step,num, &dimint); 
				if (dh)
					++x;
			}
			for (int b = nmin; b < mintmp; ++b) {
				if (pos) {
					if (((not_ones_1[b]) - (not_ones_1[a]) > 1))
						break;
					if (((not_ones_1[b]) - (not_ones_1[a]) == 0))
						continue;
				}
				dh = hemming(tab[a], tab[b], x, a, b, step,num, &dimint);
				if (dh)
					++x;
			}
		}
		y = x;
		//Step 1 only to the indiff
		for (int a = nmin; a < mintmp - 1; ++a) {
			for (int b = a + 1; b < mintmp; ++b) {
				if (pos) {
					if (((not_ones_1[b]) - (not_ones_1[a]) > 1))
						break;
					if (((not_ones_1[b]) - (not_ones_1[a]) == 0))
						continue;
				}
				dh = hemming(tab[a], tab[b], y, a, b, step,num, &dimint);
				if (dh)
					++y;
			}
		}
		// Check if all min_terms have been combinated with other min_terms
		for (int z = 0; z < nmin; ++z) {
			if (flag[z] == 1) {
				control = copy_minterms(&y, mintmp,num); // Copy non-combined min_terms
				break;
			}
		}
		if (!control)
			break;
		// Set flag for new implicants and update tables
		for (int z = 0; z < x; ++z)
			flag[z] = 1;
		for (int z = x; z < y; ++z)
			flag[z] = -1;
		mintmp = y;
		nmin = x;
		for (int z = 0; z < mintmp; ++z) {
			for (int y = 0; y < num; ++y)
				tab[z][y] = reduction[z][y];
			for (int y = 0; y < (1 << step); ++y)
				range[z][y] = intmp[z][y];
		}
		pos = 0;
		printrid(y, step,num);
	}
	printrid(y, step,num);

	//Step 2: Cover and calculate prime implicants
	int ip = fill_matrix(y, totmin, step); // Fill matrix of prime implicants
	int curr = ip + 1;
	printf("Calcolo sintesi minima in corso...\n");

	bab(0, totmin, ip, 0, &curr);// Branch and bound for calculate prime implicants

	// Print results and time
	int index = 0;
	if (sint == 1) {
		printf("Minimization SP :\n");
		printsp(ip,num);
	}
	else {
		printf("Minimization PS :\n");
		printps(ip,num);
	}
	printf("\n");
	clock_t t2 = clock();
	printf("\nTime : %lf ms\n", (double)(t2 - t));
}
