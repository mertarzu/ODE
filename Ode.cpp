
#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

typedef struct _point
{
    double x;
    double y;
}point;

double fyprime(point* p)
{
    return -2 * p->x * p->y;
}

double fy(double x)
{

    return exp(-x * x);
}

point* exact(point init, int n, double dx)
{
    point* exactarray = new point[n];

    exactarray[0] = init;

    for (int i = 1; i < n; i++)
    {
        exactarray[i].x = exactarray[i - 1].x + dx;
        exactarray[i].y = fy(exactarray[i].x);
    }
    return exactarray;
}

point* euler(point init, int n, double dx)
{
    point* eulerarray = new point[n];

    eulerarray[0] = init;

    for (int i = 1; i < n; i++)
    {
        eulerarray[i].x = eulerarray[i - 1].x + dx;
        eulerarray[i].y = eulerarray[i - 1].y + dx * fyprime(&eulerarray[i - 1]);
    }
    return eulerarray;
}

point* cauchy(point init, int n, double dx) // midpoint
{
    point* cauchyarray = new point[n];

    cauchyarray[0] = init;

    for (int i = 1; i < n; i++)
    {
        point* point12 = new point{};
        point12->x = cauchyarray[i - 1].x + .5 * dx;
        point12->y = cauchyarray[i - 1].y + .5 * dx * fyprime(&cauchyarray[i - 1]);

        cauchyarray[i].x = cauchyarray[i - 1].x + dx;
        cauchyarray[i].y = cauchyarray[i - 1].y + dx * fyprime(point12);
    }
    return cauchyarray;
}

point* heun(point init, int n, double dx) // predictor-corrector
{
    point* heunyarray = new point[n];

    heunyarray[0] = init;

    for (int i = 1; i < n; i++)
    {
        point* point1p = new point{};
        point1p->x = heunyarray[i - 1].x + dx;
        point1p->y = heunyarray[i - 1].y + dx * fyprime(&heunyarray[i - 1]);

        heunyarray[i].x = heunyarray[i - 1].x + dx;
        heunyarray[i].y = heunyarray[i - 1].y + .5 * dx * (fyprime(point1p) + fyprime(&heunyarray[i - 1]));
    }
    return heunyarray;
}

point* rungekutta(point init, int n, double dx)
{
    point* rungekuttaarray = new point[n];

    rungekuttaarray[0] = init;

    for (int i = 1; i < n; i++)
    {
        rungekuttaarray[i].x = rungekuttaarray[i - 1].x + dx;

        double k1 = dx * fyprime(&rungekuttaarray[i - 1]);
        double k2 = dx * fyprime(new point{ rungekuttaarray[i - 1].x + .5 * dx, rungekuttaarray[i - 1].y + .5 * k1 });
        double k3 = dx * fyprime(new point{ rungekuttaarray[i - 1].x + .5 * dx, rungekuttaarray[i - 1].y + .5 * k2 });
        double k4 = dx * fyprime(new point{ rungekuttaarray[i - 1].x + dx, rungekuttaarray[i - 1].y + k3 });

        rungekuttaarray[i].y = rungekuttaarray[i - 1].y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    return rungekuttaarray;
}


int main()
{
    int n = 20;
    point init = { 0,1 };
    double dx = 0.1;

    point* eulerarray = euler(init, n, dx);
    point* cauchyarray = cauchy(init, n, dx);
    point* heunarray = heun(init, n, dx);
    point* rungekuttaarray = rungekutta(init, n, dx);
    point* exactarray = exact(init, n, dx);

    cout.setf(ios::fixed);
    cout.setf(ios::showpoint);
    cout << setw(4) << left << "i" << setw(6) << left << "x"
        << setw(12) << left << "y (euler)" << setw(12) << left << "y (cauchy)"
        << setw(12) << left << "y (heun)" << setw(16) << left << "y (runge-kutta)"
        << setw(12) << left << "y (exact)"
        << endl;
    for (int i = 0; i < n; i++)
    {
        cout << setw(4) << left << i << setw(6) << setprecision(2) << left << eulerarray[i].x << setw(12) << setprecision(8) << left << eulerarray[i].y
            << setw(12) << setprecision(8) << left << cauchyarray[i].y
            << setw(12) << setprecision(8) << left << heunarray[i].y
            << setw(16) << setprecision(8) << left << rungekuttaarray[i].y
            << setw(12) << setprecision(8) << left << exactarray[i].y
            << endl;
    }
}

