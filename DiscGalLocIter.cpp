// project
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const double pi = 3.1415926;

double U_l_t(double t)
{
    return 5 * exp(-t);
    //return 0;
}

double U_0_t(double t)
{
    return 0;
}

double U_x_0(double x)
{
    return 5. * x;
    // return sin(pi * x);
      /*if (x <= 0.5 && x >= 0)
          return x;
      else return 0;*/
}

double phi1(double x, double xj)
{
    return x - xj;
}

double phi2(double x, double xj)
{
    return (x - xj) * (x - xj);
}

double simpson(double a, double b, double (*pf)(double))
{
    return (b - a) / 6 * (pf(a) + 4 * pf((a + b) / 2) + pf(b));
}

double cotes4(double a, double b, double (*pf1)(double))
{
    double h = (b - a) / 4;
    return (4 * h / 90) * (7 * pf1(a) + 32 * pf1(a + h) + 12 * pf1(a + 2 * h)
        + 32 * pf1(a + 3 * h) + 7 * pf1(b));
}

double cotes4_for_mult(double a, double b, double (*pf1)(double))
{
    double xj = (a + b) / 2;
    double h = (b - a) / 4;
    return (4 * h / 90) * (7 * pf1(a) * phi1(a, xj) + 32 * pf1(a + h) * phi1(a + h, xj) + 12 * pf1(a + 2 * h) * phi1(a + 2 * h, xj)
        + 32 * pf1(a + 3 * h) * phi1(a + 3 * h, xj) + 7 * pf1(b) * phi1(b, xj));
}

double simpson_for_mult(double a, double b, double (*pf)(double))
{
    double xj = (a + b) / 2;
    return (b - a) / 6 * (pf(a) * phi1(a, xj) + 4 * pf((a + b) / 2) * phi1((a + b) / 2, xj) + pf(b) * phi1(b, xj));
}

double simpson_for_2mult(double a, double b, double (*pf)(double))
{
    double xj = (a + b) / 2;
    return (b - a) / 6 * (pf(a) * phi2(a, xj) + 4 * pf((a + b) / 2) * phi2((a + b) / 2, xj) + pf(b) * phi2(b, xj));
}

double cotes4_for_2mult(double a, double b, double (*pf1)(double))
{
    double xj = (a + b) / 2;
    double h = (b - a) / 4;
    return (4 * h / 90) * (7 * pf1(a) * phi2(a, xj) + 32 * pf1(a + h) * phi2(a + h, xj) + 12 * pf1(a + 2 * h) * phi2(a + 2 * h, xj)
        + 32 * pf1(a + 3 * h) * phi2(a + 3 * h, xj) + 7 * pf1(b) * phi2(b, xj));
}

double* dif(const double x[], const double y[], int n)
{
    double* vec = new double[n];
    for (int i = 0; i < n; ++i)
        vec[i] = x[i] - y[i];
    return vec;
}

double* matrix_x_vec(const double A[], const double X[], int n)
{
    double* temp = new double[n];
    for (int i = 0; i < n; ++i)
    {
        temp[i] = 0;
        for (int j = 0; j < n; ++j)
            temp[i] += A[i * n + j] * X[j];
    }
    return temp;
}

double norm(const double X[], int n)
{
    double temp = 0;
    for (int i = 0; i < n; ++i)
        temp += X[i] * X[i];
    return sqrt(temp);
}

void normir(double X[], int n)
{
    double N = norm(X, n);
    for (int i = 0; i < n; ++i)
        X[i] = X[i] / N;
}

template<typename T>
void equate(T u1[], const T u2[], int n)
{
    for (int i = 0; i < n; ++i)
        u1[i] = u2[i];
}

double max_eigenvalue(const double A[], int n)
{
    int k = 0;
    double xTx = 0;
    double xAx = 0;
    double eps = 0.001;
    double* X0 = new double[n];
    double* X1 = new double[n];
    for (int i = 0; i < n; ++i)
        X1[i] = 1 / sqrt(n);

    do
    {
        equate<double>(X0, X1, n);
        equate<double>(X1, matrix_x_vec(A, X0, n), n);
        normir(X1, n);
    } while (norm(dif(X1, X0, n), n) > eps);
    for (int i = 0; i < n; ++i)
    {
        xTx += X1[i] * X1[i];
        for (int j = 0; j < n; ++j)
            xAx += X1[i] * A[i * n + j] * X1[j];
    }
    delete[] X0;
    delete[] X1;
    return  xAx / xTx;
}

void matrix_fill(double A[], double a, double h, int n)
{
    double** B = new double* [3 * n];

    for (int i = 0; i < 3 * n; ++i)
    {
        B[i] = new double[3 * n];
        for (int j = 0; j < 3 * n; ++j)
            B[i][j] = 0;
    }

    double c = 1 / (h * h);
    double d = -1 / (h * h);

    B[0][0] = 18 * c;
    B[0][1] = 9 * d;

    B[0][n] = 2 * d;
    B[0][n + 1] = 9 * c;

    B[0][2 * n] = 12 * c;
    B[0][2 * n + 1] = 9 * d;

    for (int i = 1; i < n - 1; ++i)
    {
        B[i][i - 1] = 9 * d;
        B[i][i] = 18 * c;
        B[i][i + 1] = 9 * d;

        B[i][n + i - 1] = 7 * d;
        B[i][n + i] = 2 * d;
        B[i][n + i + 1] = 9 * c;

        B[i][2 * n + i - 1] = 3 * d;
        B[i][2 * n + i] = 12 * c;
        B[i][2 * n + i + 1] = 9 * d;
    }

    B[n - 1][n - 1] = 18 * c;
    B[n - 1][n - 2] = 9 * d;

    B[n - 1][2 * n - 1] = 2 * d;
    B[n - 1][2 * n - 2] = 7 * d;

    B[n - 1][3 * n - 1] = 3 * d;
    B[n - 1][3 * n - 2] = 12 * c;


    B[n][0] = 6 * d;
    B[n][1] = 21 * d;


    B[n][n] = 54 * c;
    B[n][n + 1] = 21 * c;

    B[n][2 * n] = 24 * d;
    B[n][2 * n + 1] = 21 * d;


    for (int i = n + 1; i < 2 * n - 1; ++i)
    {
        B[i][i - n - 1] = 27 * c;
        B[i][i - n] = 6 * d;
        B[i][i - n + 1] = 21 * d;

        B[i][i - 1] = 21 * c;
        B[i][i] = 54 * c;
        B[i][i + 1] = 21 * c;

        B[i][i + n - 1] = 9 * c;
        B[i][i + n] = 24 * d;
        B[i][i + n + 1] = 21 * d;
    }

    B[2 * n - 1][n - 1] = 21 * d;
    B[2 * n - 1][n - 2] = 6 * d;

    B[2 * n - 1][2 * n - 1] = 54 * c;
    B[2 * n - 1][2 * n - 2] = 21 * c;

    B[2 * n - 1][3 * n - 1] = 24 * d;
    B[2 * n - 1][3 * n - 2] = 9 * c;


    B[2 * n][0] = 60 * c;
    B[2 * n][1] = 15 * d;

    B[2 * n][n] = 40 * d;
    B[2 * n][n + 1] = 15 * c;

    B[2 * n][2 * n] = 90 * c;
    B[2 * n][2 * n + 1] = 15 * d;

    for (int i = 2 * n + 1; i < 3 * n - 1; ++i)
    {
        B[i][i - 2 * n - 1] = 45 * d;
        B[i][i - 2 * n] = 60 * c;
        B[i][i - 2 * n + 1] = 15 * d;

        B[i][i - n - 1] = 35 * d;
        B[i][i - n] = 40 * d;
        B[i][i - n + 1] = 15 * c;

        B[i][i - 1] = 15 * d;
        B[i][i] = 90 * c;
        B[i][i + 1] = 15 * d;
    }

    B[3 * n - 1][n - 2] = 45 * d;
    B[3 * n - 1][n - 1] = 60 * c;

    B[3 * n - 1][2 * n - 2] = 35 * d;
    B[3 * n - 1][2 * n - 1] = 40 * d;

    B[3 * n - 1][3 * n - 2] = 15 * d;
    B[3 * n - 1][3 * n - 1] = 90 * c;


    for (int i = 0; i < 3 * n; ++i)
    {
        for (int j = 0; j < 3 * n; ++j)
            A[i * 3 * n + j] = B[i][j];
    }


    for (int i = 0; i < n; ++i)
        delete[] B[i];
    delete[] B;
}

int* reg(int n)
{
    int* tetta = new int[n];
    int k = n;
    int m = 0;
    int* temp;
    if (k == 1)
    {
        temp = new int[1];
        temp[0] = 1;
        return temp;
    }
    do
    {
        ++m;
        if (k % 2 == 0) k /= 2;
        else --k;
    } while (k != 1);
    ++m;
    temp = new int[m];
    k = n;
    for (int j = m - 1; j >= 0; --j)
    {
        temp[j] = k;
        if (k % 2 == 0) k /= 2;
        else --k;
    }
    int* tetta1_temp = new int[1];
    int* tetta2_temp = new int[1];
    tetta1_temp[0] = 1;
    for (int i = 1; i < m - 1; ++i)
    {
        delete[] tetta2_temp;
        tetta2_temp = new int[temp[i - 1]];
        equate<int>(tetta2_temp, tetta1_temp, temp[i - 1]);
        delete[] tetta1_temp;
        tetta1_temp = new int[temp[i]];
        tetta1_temp[0] = 1;
        for (int j = 1; j < temp[i]; ++j)
        {
            if (temp[i] / temp[i - 1] == 2)
            {
                if (temp[i + 1] / temp[i] == 2 || i == m - 1)
                {
                    if ((j - 1) % 2 == 1)
                        tetta1_temp[j] = tetta2_temp[j / 2];
                    else tetta1_temp[j] = 4 * temp[i - 1] - tetta1_temp[j - 1];
                }
                if (temp[i + 1] % temp[i] == 1)
                {

                    if ((j - 1) % 2 == 1)
                        tetta1_temp[j] = tetta2_temp[j / 2];
                    else tetta1_temp[j] = 4 * temp[i - 1] + 2 - tetta1_temp[j - 1];
                }

            }
            if (temp[i] % temp[i - 1] == 1)
            {
                for (int r = 0; r < temp[i] - 1; ++r)
                    tetta1_temp[r] = tetta2_temp[r];
                tetta1_temp[temp[i] - 1] = temp[i];
            }
        }
    }
    for (int l = 0; l < temp[m - 2]; ++l)
        tetta[0] = 1;
    if (temp[m - 1] % 2 == 1)
    {
        for (int j = 1; j < n - 1; ++j)
            tetta[j] = tetta1_temp[j];
        tetta[n - 1] = temp[m - 1];
    }
    else
        for (int j = 1; j < n; ++j)
        {
            if ((j - 1) % 2 == 1)
                tetta[j] = tetta1_temp[j / 2];
            else tetta[j] = 4 * (temp[m - 1] / 2) - tetta[j - 1];
        }
    delete[] tetta1_temp;
    delete[] tetta2_temp;
    delete[] temp;
    return tetta;
}

void linear_iteration(double time, double a, double L0, double L, double tau, double T, int n, double (*u_0_x)(double), double (*u_t_0)(double), double (*u_t_l)(double))
{
    int re;
    int k = 0;
    double xi = L0;
    double pi = 3.1415926;
    ofstream file;
    ofstream file1;
    file.open("F:/Бауманка/диплом/res.txt");
    file1.open("F:/Бауманка/диплом/koef.txt");
    //file1 << n << " " << time << endl;
    double h = L / n;
    double** U1 = new double* [n];
    double** U2 = new double* [n];
    double** U_temp = new double* [n];
    double* M = new double[9 * n * n];
    matrix_fill(M, a, h, n);
    double eps = 0.000000001;
    //for (int i = 0; i < 2 * n; ++i)
    //{
    //    for (int j = 0; j < 2 * n; ++j)
    //    {
    //        cout << M[i * 2 * n + j] << " ";
    //        //file << A[i * 2 * n + j] << " ";
    //    }
    //    cout << endl;
    //    //file << endl;
    //}

    for (int i = 0; i < n; ++i)
    {

        U_temp[i] = new double[3];
        U2[i] = new double[3];
        U1[i] = new double[3];
        U1[i][0] = cotes4(xi, xi + h, u_0_x) / h;
        U1[i][1] = cotes4_for_mult(xi, xi + h, u_0_x) * 6. / h / h;
        U1[i][2] = cotes4_for_2mult(xi, xi + h, u_0_x) * 30. / h / h / h - 5. / 2 * cotes4(xi, xi + h, u_0_x) / h;
        xi += h;
        //cout << U1[i][0] << " " << U1[i][1] << " " << U1[i][2] << endl;
        //file << L0 + h * i << " " << U1[i][0] - U1[i][1] + U1[i][2] << endl;
        //file << L0 + h * (i + 1) << " " << U1[i][0] + U1[i][1] + U1[i][2] << endl;
    }

    double lambda_max = 1.2 * max_eigenvalue(M, 3 * n);
    //double lambda_max = 1200 / (h * h);

    int p = ceil(pi / 4 * sqrt(tau * lambda_max + 1));
    cout << lambda_max << " " << p << endl;
    double* B = new double[p];
    int* tetta = new int[p];
    equate<int>(tetta, reg(p), p);

    for (int i = 0; i < p; ++i)
        B[i] = cos(tetta[i] * pi / (2 * p));

    double* A = new double[p];

    for (int i = 0; i < p; ++i)
        A[i] = lambda_max / (1 + B[0]) * (B[0] - B[i]);

    while (k != T)
    {
        ++k;
        for (int s = 0; s < n; ++s)
            equate<double>(U2[s], U1[s], 3);
        //U2[0][0] = (4 * U_0_t(tau * k) + 9 * U1[0][0] + 4 * U1[0][1] + 9 * tau / (h * h) * (U1[1][0] - 2 * U1[0][0] + U1[1][1] - 2 * U1[0][1] + U_0_t(tau * (k - 1)))) / 13;
        //U2[0][1] = U2[0][0] - U_0_t(tau * k);
        // file << U2[0] << " ";
        // U2[n - 1][0] = (4 * U_l_t(tau * k) + 9 * U1[n - 1][0] - 4 * U1[n - 1][1] + 9 * tau / (h * h) * (U1[n - 2][0] - 2 * U1[n - 1][0] + 2 * U1[n - 1][1] - U1[n - 2][1] + U_l_t(tau * (k - 1)))) / 13;
        //U2[n - 1][1] = U_l_t(tau * k) - U2[n - 1][0];

       /* U2[0][0] = U_0_t(tau * k);
        U2[n - 1][0] = U_l_t(tau * k);*/

        for (int i = 0; i < n; ++i)
            equate<double>(U_temp[i], U1[i], 3);

        for (int i = 0; i < p; ++i)
        {

            /* U2[0][0] = 1 / (1 + tau * A[i]) * (U1[0][0] + tau * A[i] * U_temp[0][0] + 16./79*(U_0_t(tau * k) - U_0_t(tau * (k - 1))) +  tau * a / (h * h) * (53* U_0_t(tau * (k - 1))-1066*U_temp[0][0]
                 + 1013*U_temp[1][0] - 1039 * U_temp[0][1] - 934*U_temp[1][1]+29*U_temp[0][2]+200*U_temp[1][2])/79);*/
            U2[0][0] = 1 / (1 + tau * A[i]) * (U1[0][0] + tau * A[i] * U_temp[0][0] + 1. / 9 * (U_0_t(tau * k) - U_0_t(tau * (k - 1))) - tau * a / (h * h) * (78 * U_temp[0][0]
                - 78 * U_temp[1][0] + 78 * U_temp[0][1] + 78 * U_temp[1][1] - 18 * U_temp[0][2] - 78 * U_temp[1][2]) / 9);


            /* U2[0][1] = 1 / (1 + tau * A[i]) * (U1[0][1] + tau * A[i] * U_temp[0][1] - 192. / 316 * (U_0_t(tau * k) - U_0_t(tau * (k - 1))) - 3. / 316 * tau * a / (h * h) * (104* U_0_t(tau * (k - 1))+2688 * U_temp[0][0]
                 - 2584 * U_temp[1][0] + 2717 * U_temp[0][1] + 2347 * U_temp[1][1] + 1143 * U_temp[0][2] - 385 * U_temp[1][2]));*/
            U2[0][1] = 1 / (1 + tau * A[i]) * (U1[0][1] + tau * A[i] * U_temp[0][1] - 1. / 3 * (U_0_t(tau * k) - U_0_t(tau * (k - 1))) - 2 * tau * a / (h * h) * (11 * U_temp[0][0]
                - 11 * U_temp[1][0] + 11 * U_temp[0][1] + 11 * U_temp[1][1] + 9 * U_temp[0][2] - 11 * U_temp[1][2]));

            /*U2[0][2] = 1 / (1 + tau * A[i]) * (U1[0][2] + tau * A[i] * U_temp[0][2] + 60. / 316 * (U_0_t(tau * k) - U_0_t(tau * (k - 1))) + 5. / 316 * tau * a / (h * h) * (20 * U_0_t(tau * (k - 1)) - 760 * U_temp[0][0]
                + 740 * U_temp[1][0] - 799 * U_temp[0][1] - 661 * U_temp[1][1] - 709 * U_temp[0][2] + 71 * U_temp[1][2]));*/
            U2[0][2] = 1 / (1 + tau * A[i]) * (U1[0][2] + tau * A[i] * U_temp[0][2] + 5. / 9 * (U_0_t(tau * k) - U_0_t(tau * (k - 1))) + 4. / 3 * tau * a / (h * h) * (-2 * U_temp[0][0]
                + 2 * U_temp[1][0] - 2 * U_temp[0][1] - 2 * U_temp[1][1] - 3 * U_temp[0][2] + 2 * U_temp[1][2]));


            // U2[0][0] = U2[0][1] + U_0_t(tau * k );
            /* U2[0][0] = 1 / (1 + tau * A[i]) * (U1[0][0] + tau * A[i] * U_temp[0][0] - tau * a / (h * h) * (U_temp[0][0]
                 +2* U_temp[0][1] - U_temp[1][0]));
             U2[0][1] = 1 / (1 + tau * A[i]) * (U1[0][1] + tau * A[i] * U_temp[0][1] +9/4* tau * a / (h * h) * (U_temp[1][1]
                 - U_temp[0][1]));*/

                 /* U2[n - 1][0] = 1 / (1 + tau * A[i]) * (U1[n - 1][0] + tau * A[i] * U_temp[n - 1][0] + 16. / 79 * (U_l_t(tau * k) - U_l_t(tau * (k - 1))) + 1./79*tau * a / (h * h) * (53*u_t_l(tau*(k-1))+1013*U_temp[n - 2][0]
                      - 1066*U_temp[n - 1][0] + 934*U_temp[n - 2][1] + 1039 * U_temp[n - 1][1]+200*U_temp[n-2][2]+29*U_temp[n-1][2]));*/
            U2[n - 1][0] = 1 / (1 + tau * A[i]) * (U1[n - 1][0] + tau * A[i] * U_temp[n - 1][0] + 1. / 5 * (U_l_t(tau * k) - U_l_t(tau * (k - 1))) - tau * a / (h * h) * (-18 * U_temp[n - 2][0]
                + 18 * U_temp[n - 1][0] - 14 * U_temp[n - 2][1] - 22 * U_temp[n - 1][1] - 6 * U_temp[n - 2][2] - 18 * U_temp[n - 1][2]) / 5);

            /* U2[n - 1][1] = 1 / (1 + tau * A[i]) * (U1[n - 1][1] + tau * A[i] * U_temp[n - 1][1] + 192. / 316 * (U_l_t(tau * k) - U_l_t(tau * (k - 1))) + 3./316*tau * a / (h * h) * (104*u_t_l(tau * (k - 1))-2584*U_temp[n - 2][0]
                 + 2688*U_temp[n - 1][0] - 2347*U_temp[n - 2][1] - 2717 * U_temp[n - 1][1]-385*U_temp[n-2][2]+1143*U_temp[n-1][2]));*/
            U2[n - 1][1] = 1 / (1 + tau * A[i]) * (U1[n - 1][1] + tau * A[i] * U_temp[n - 1][1] + 7. / 15 * (U_l_t(tau * k) - U_l_t(tau * (k - 1))) - tau * a / (h * h) * (594 * U_temp[n - 2][0]
                - 594 * U_temp[n - 1][0] + 462 * U_temp[n - 2][1] + 726 * U_temp[n - 1][1] + 198 * U_temp[n - 2][2] - 906 * U_temp[n - 1][2]) / 15);

            /* U2[n - 1][2] = 1 / (1 + tau * A[i]) * (U1[n - 1][2] + tau * A[i] * U_temp[n - 1][2] + 60. / 316 * (U_l_t(tau * k) - U_l_t(tau * (k - 1))) + 5./316*tau * a / (h * h) * (20 * u_t_l(tau * (k - 1)) + 740 * U_temp[n - 2][0]
                 - 760 * U_temp[n - 1][0] + 661 * U_temp[n - 2][1] + 799 * U_temp[n - 1][1] + 71 * U_temp[n - 2][2] - 709 * U_temp[n - 1][2]));*/
            U2[n - 1][2] = 1 / (1 + tau * A[i]) * (U1[n - 1][2] + tau * A[i] * U_temp[n - 1][2] + 1. / 3 * (U_l_t(tau * k) - U_l_t(tau * (k - 1))) - tau * a / (h * h) * (-108 * U_temp[n - 2][0]
                + 108 * U_temp[n - 1][0] - 84 * U_temp[n - 2][1] - 132 * U_temp[n - 1][1] - 36 * U_temp[n - 2][2] + 192 * U_temp[n - 1][2]) / 3);
            //U2[n - 1][0] = U_l_t(tau * k) - U2[n - 1][1];
           /* U2[n-1][0] = 1 / (1 + tau * A[i]) * (U1[n-1][0] + tau * A[i] * U_temp[n-1][0] - tau * a / (h * h) * (U_temp[n-1][0]
                - 2*U_temp[n-1][1] +  U_temp[n-2][0]));
            U2[n-1][1] = 1 / (1 + tau * A[i]) * (U1[n-1][1] + tau * A[i] * U_temp[n-1][1] + 9 / 4 * tau * a / (h * h) * (U_temp[n-2][1]
                - U_temp[n-1][1]));*/

            for (int s = 1; s < n - 1; ++s)
            {
                U2[s][0] = 1 / (1 + tau * A[i]) * (U1[s][0] + tau * A[i] * U_temp[s][0] + tau * a / (h * h) * (9 * U_temp[s + 1][0]
                    - 18 * U_temp[s][0] + 9 * U_temp[s - 1][0] - 9 * U_temp[s + 1][1] + 2 * U_temp[s][1] + 7 * U_temp[s - 1][1] + 9 * U_temp[s + 1][2] - 12 * U_temp[s][2] + 3 * U_temp[s - 1][2]));

                U2[s][1] = 1 / (1 + tau * A[i]) * (U1[s][1] + tau * A[i] * U_temp[s][1] + 3. * tau * a / (h * h) * (7 * U_temp[s + 1][0]
                    + 2 * U_temp[s][0] - 9 * U_temp[s - 1][0] - 7 * U_temp[s + 1][1] - 18 * U_temp[s][1] - 7 * U_temp[s - 1][1] + 7 * U_temp[s + 1][2] + 8 * U_temp[s][2] - 3 * U_temp[s - 1][2]));

                U2[s][2] = 1 / (1 + tau * A[i]) * (U1[s][2] + tau * A[i] * U_temp[s][2] + 5. * tau * a / (h * h) * (3 * U_temp[s + 1][0]
                    - 12 * U_temp[s][0] + 9 * U_temp[s - 1][0] - 3 * U_temp[s + 1][1] + 8 * U_temp[s][1] + 7 * U_temp[s - 1][1] + 3 * U_temp[s + 1][2] - 18 * U_temp[s][2] + 3 * U_temp[s - 1][2]));

            }
            /* if (i == 2)
                 cout << 1 / (1 + tau * A[i]) * (U1[1][0] + tau * A[i] * U_temp[1][0] + tau * a / (h * h) * (U_temp[2][0]
                     - 2 * U_temp[1][0] + U_temp[0][0])) << endl;
             cout << "tau" << " | " << "A[i]" << " | " << "U1[0]" << " | " << "Utemp[0]" << " | " << "Utemp[1]" << " | " << "Utemp[2]" << " | " << "h" << " | " << endl;
             cout << tau << " | " << A[i] << " | " << U1[1][0] << " | " << U_temp[0][0] << " | " << U_temp[1][0] << " | " << U_temp[2][0] << " | " << h << " | " << endl;
             cout << U2[0][0] << endl;

             cin >> re;*/
            for (int s = 0; s < n; ++s)
                equate<double>(U_temp[s], U2[s], 3);

        }
        for (int s = 0; s < n; ++s)
            equate<double>(U1[s], U2[s], 3);
        if (fabs(tau * k - time) < eps)
        {
            //cout << tau * k << endl;
            file1 << n << " " << tau * k << endl;
            /* cout << "lambda:" << " " << lambda_max << endl;
             cout << "tetta[i]" << " | " << "U2[k-1]" << " | " << "U2[k]" << " | " << "U2[k+1]" << " | " << endl;
             cout << tetta[0] << " | " << "000" << " | " << U2[0][0] << " | " << U2[1][0] << " | " << endl;
             for (int u = 1; u < n - 1; ++u)
             {
                 cout << tetta[u] << " | " << U2[u - 1][0] << " | " << U2[u][0] << " | " << U2[u + 1][0] << " | " << endl;
             }
             cout << tetta[n - 1] << " | " << U2[n - 2][0] << " | " << U2[n - 1][0] << " | " << "000" << " | " << endl;*/
            for (int i = 0; i < n; ++i)
            {
                file1.precision(12);
                //cout << U2[i][0] << "          " << U2[i][1] << endl;
                file << L0 + h * i << " " << U2[i][0] - U2[i][1] + U2[i][2] << endl;
                file << L0 + h * (i + 1) << " " << U2[i][0] + U2[i][1] + U2[i][2] << endl;
                file1 << U2[i][0] << " " << U2[i][1] << " " << U2[i][2] << endl;
                // cout << U2[i][0] << " " << U2[i][1] << " " << U2[i][2] << endl;
            }
            cin >> re;
        }
    }


    delete[] B;
    //delete[] a;
    delete[] U1;
    delete[] U2;
    delete[] A;
    delete[] tetta;

}

int main()
{
    double time = 0.1;
    ofstream fileE;
    fileE.open("F:/Бауманка/диплом/выход.txt");
    int n = 100;
    double a = 1;
    double L = 1;
    double h = L / n;
    double tau = 0.01;
    //double tau = 0.001;
    double T = 10000000;
    double* A = new double[9 * n * n];
    matrix_fill(A, a, h, n);
    /*for (int i = 0; i < 3*n; ++i)
    {
        for (int j = 0; j < 3 * n; ++j)
        {
            cout << A[i * 3 * n + j] << " ";
            fileE << A[i * 3 * n + j] << " ";
        }
        cout << endl;
        fileE << endl;
    }*/
    cout << max_eigenvalue(A, 3 * n) << endl;
    linear_iteration(time, a, 0, L, tau, T, n, U_x_0, U_0_t, U_l_t);
}

