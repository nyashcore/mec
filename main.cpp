#include <iostream>
#include <cmath>

#define NUMBER_OF_SB 50
#define NUMBER_OF_SW 50
#define NUMBER_OF_S 50
#define DIMENSION 4
#define SIGMA 0.0001
// rastrigin
#define BORDER 1.
// sphere
//#define BORDER 10.0
// rosenbrock
//#define BORDER 2.0
//ackley
//#define BORDER 5.0
#define A 10
#define PI 3.14

using namespace std;

double*** Sb = NULL;
double*** Sw = NULL;
double** Cb = NULL;
double** Cw = NULL;
int* CbInd = NULL;
int* CwInd = NULL;
int iterations = 0;
int gnuplot = 1;
double border = BORDER;
double sigma = SIGMA;

// Rosenbrock
double fros(double* point) {
    double result = 0.0;

    for(int j = 0; j < DIMENSION - 1; j++) {
        result += 100.0 * pow(1 - point[j], 2) + 100 * pow(point[j + 1] - pow(point[j], 2), 2);
    }

    return result;
}

// Sphere
double fs(double* point) {
    double result = 0.0;
    for(int j = 0; j < DIMENSION; j++) {
        result = result + pow(point[j], 2);
    }
    return result;
}

// Ackley
double f(double* point)
{
    double result = 0;
    double a = 20.0;
    double b = 0.2;
    double c = 3.14 * 2;
    double sum1 = 0, sum2 = 0;

    for(size_t i = 0; i < DIMENSION; i++) {
        sum1 += pow(point[i], 2);
        sum2 += cos(c*point[i]);
    }
    result += -a*exp(-b*sqrt(sum1/(double)DIMENSION));
    result += -exp(sum2/(double)DIMENSION);
    result += a;
    result += exp(1);
    return result;
}

// Rastrigin
double fras(double* point) {
    double result = A * DIMENSION;

    for(int i = 0; i < DIMENSION; i++) {
        result += pow(point[i], 2);
        result -= A * cos(2 * PI * point[i]);
    }

    return result;
}

double myRand() {
    double random;

    random = (rand() / (RAND_MAX + 1.0) * 6 * sigma) - 3 * sigma;

    return random;
}

double myRand2() {
    double random;

    random = (rand() / (RAND_MAX + 1.0) * border * 2.0) - border;

    return random;
}

void print() {
    for(int i = 0; i < NUMBER_OF_SB ; i++) {
        cout << "Sb: " << i << endl;
        for(int j = 0; j < NUMBER_OF_S; j++) {
            for(int k = 0; k < DIMENSION + 1; k++) {
                cout << Sb[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << "Winner: " << CbInd[i] << " : ";
        for(int j = 0; j < DIMENSION + 1; j++) {
            cout << Cb[i][j] << " ";
        }
        cout << endl;
    }
    for(int i = 0; i < NUMBER_OF_SW ; i++) {
        cout << "Sw: " << i << endl;
        for(int j = 0; j < NUMBER_OF_S; j++) {
            for(int k = 0; k < DIMENSION + 1; k++) {
                cout << Sw[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << "Winner: " << CwInd[i] << " : ";
        for(int j = 0; j < DIMENSION + 1; j++) {
            cout << Cw[i][j] << " ";
        }
        cout << endl;
    }

    return;
}

double** initialize() {
    double** newS;

    newS = (double**)malloc(NUMBER_OF_S * (DIMENSION + 1) * sizeof(double));
    for(int i = 0; i < NUMBER_OF_S; i++) {
        newS[i] = (double*)malloc((DIMENSION + 1) * sizeof(double));
    }

    for(int i = 0; i < DIMENSION; i++) {
        newS[0][i] = myRand2();
    }
    for(int i = 1; i < NUMBER_OF_S; i++) {
        for(int j = 0; j < DIMENSION; j++) {
            newS[i][j] = newS[0][j] + myRand();
        }
    }

    return newS;
}

void rebuildSw() {
    for(int i = 0; i < NUMBER_OF_SW; i++) {
        for(int k = 0; k < NUMBER_OF_S; k++) {
            if(k != CwInd[i]) {
                for (int j = 0; j < DIMENSION; j++) {
                    Sw[i][k][j] = Cw[i][j] + myRand();
                }
            }
        }
    }

    return;
}

void rebuildSb() {
    for(int i = 0; i < NUMBER_OF_SB; i++) {
        for(int k = 0; k < NUMBER_OF_S; k++) {
            if(k != CbInd[i]) {
                for(int j = 0; j < DIMENSION; j++) {
                    Sb[i][k][j] = Cb[i][j] + myRand();
                }
            }
        }
    }

    return;
}

void scores() {
    for(int i = 0; i < NUMBER_OF_SB; i++) {
        for(int j = 0; j < NUMBER_OF_S; j++) {
            Sb[i][j][DIMENSION] = f(Sb[i][j]);
        }
    }
    for(int i = 0; i < NUMBER_OF_SW; i++) {
        for(int j = 0; j < NUMBER_OF_S; j++) {
            Sw[i][j][DIMENSION] = f(Sw[i][j]);
        }
    }

    return;
}

int localWinners() {
    if(Cb == NULL) {
        Cb = (double**)malloc(NUMBER_OF_SB * (DIMENSION + 1) * sizeof(double));
        CbInd = (int*)malloc(NUMBER_OF_SB * sizeof(int));
        for(int i = 0; i < NUMBER_OF_SB; i++) {
            Cb[i] = Sb[i][0];
        }
    }
    if(Cw == NULL) {
        Cw = (double**)malloc(NUMBER_OF_SW * (DIMENSION + 1) * sizeof(double));
        CwInd = (int*)malloc(NUMBER_OF_SW * sizeof(int));
        for(int i = 0; i < NUMBER_OF_SW; i++) {
            Cw[i] = Sw[i][0];
        }
    }

    int flag = 0;
    for(int i = 0; i < NUMBER_OF_SB; i++) {
        for(int j = 0; j < NUMBER_OF_S; j++) {
            if(Sb[i][j][DIMENSION] < Cb[i][DIMENSION]) {
                Cb[i] = Sb[i][j];
                CbInd[i] = j;
                flag = 1;
            }
        }
    }
    if(flag == 0) {
        return 0;
    }
    for(int i = 0; i < NUMBER_OF_SW; i++) {
        for(int j = 0; j < NUMBER_OF_S; j++) {
            if(Sw[i][j][DIMENSION] < Cw[i][DIMENSION]) {
                Cw[i] = Sw[i][j];
                CwInd[i] = j;
            }
        }
    }

    return 1;
}

void global() {
    double* temp;
    for(int i = 0; i < NUMBER_OF_SB; i++) {
        for(int j = 0; j < NUMBER_OF_SW; j++) {
            if(Cb[i][DIMENSION] > Cw[j][DIMENSION]) {
                temp = Cb[i];
                Cb[i] = Cw[j];
                Cw[j] = temp;
                Sw[j][CwInd[j]] = Cw[j];
                Sb[i][CbInd[i]] = Cb[i];
                break;
            }
        }
    }
    int flag;
    for(int i = 0; i < NUMBER_OF_SW; i++) {
        flag = 0;
        for(int j = 0; j < NUMBER_OF_SB; j++) {
            if(Cw[i][DIMENSION] < Cb[j][DIMENSION]) {
                temp = Cb[j];
                Cb[j] = Cw[i];
                Cw[i] = temp;
                Sw[i][CwInd[i]] = Cw[i];
                Sb[j][CbInd[j]] = Cb[j];
                flag = 1;
                break;
            }
        }
        if(flag == 0) {
            free(Sw[i]);
            Sw[i] = initialize();
            Sw[i][CwInd[i]][0] = 0.0;
            Cw[i] = Sw[i][CwInd[i]];
        }
    }

    return;
}

void print_best() {
    double* best = Cb[0];


    for(int i = 1; i < NUMBER_OF_SB; i++) {
        if(Cb[i][DIMENSION] < best[DIMENSION]) {
            best = Cb[i];
        }
    }
    if(gnuplot) {
        cout << iterations << " " << best[DIMENSION] << endl;
        return;
    } else {
        cout << endl << "iter: " << iterations << endl;
        for (int i = 0; i < DIMENSION + 1; i++) {
            cout << best[i] << " " << endl;
        }
    }

    return;
}

int main() {
    srand(time(NULL));
    Sb = (double***)malloc(NUMBER_OF_SB * NUMBER_OF_S * (DIMENSION + 1) * sizeof(double));
    Sw = (double***)malloc(NUMBER_OF_SW * NUMBER_OF_S * (DIMENSION + 1) * sizeof(double));

    for(int i = 0; i < NUMBER_OF_SB; i++) {
        Sb[i] = initialize();
    }
    for(int i = 0; i < NUMBER_OF_SW; i++) {
        Sw[i] = initialize();
    }
    scores();

//    for(int i = 0; i < 20; i++) {
//        cout << myRand2() << endl;
//    }
    while(1) {
        localWinners();
        print_best();
        rebuildSb();
        rebuildSw();
        scores();
        if(localWinners() == 0) {
            break;
        }
        global();
        rebuildSb();
        rebuildSw();
        scores();
        iterations++;
    }

    cout << "Finished, number of iterations: " << iterations << endl;
    print_best();

    for(int i = 0; i < NUMBER_OF_SB; i++) {
        for(int j = 0; j < NUMBER_OF_S; j++) {
            free(Sb[i][j]);
        }
        free(Sb[i]);
    }
    free(Sb);
    for(int i = 0; i < NUMBER_OF_SW; i++) {
        for(int j = 0; j < NUMBER_OF_S; j++) {
            free(Sw[i][j]);
        }
        free(Sw[i]);
    }
    free(Sw);

//    double* point;
//    point = (double*)malloc((DIMENSION+1) * sizeof(double));
//    for(int i = 0; i < DIMENSION; i++) {
//        point[i] = 1.0;
//        cout << point[i] << " " << endl;
//    }
//    cout << Rosenbrock(point) << endl;

    return 0;
}