#include <iostream>
using namespace std;
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

double au (double x) {
    if (x <= -0.5) {
        return 0;
    } else if (x > -0.25) {
        return 1;
    } else {
        return 4*(x + 0.5);
    }
}

double F(double x) {
    return -0.5*x*x; //Для (не)линейного случая поменять эту строчку;
}

double V1(double x1, double x2, double g) {
    return x1 - g*(F(x1) - F(x2));
}

void U (vector<double> u0, vector<double>& u, double h, double t) {

    u[0] = 0;
    for(int i = 1; i < u.size() - 1; i++) {
        u[i] = 0.5*(u0[i] + V1(u0[i], u0[i-1], (t/h)) - (t/h)*(F(V1(u0[i+1], u0[i], (t/h) )) - F(V1(u0[i], u0[i-1], (t/h) ))));
    }
    u[u.size()-1] = 1;
}

int main () {
    FILE * fp;
    fp = fopen("ves1().txt", "w");

    double h = 0.01, t = 0.01;
    double start = -1, finish = 1;

    double dvl = 0, Dvl = 0, dvc = 0, Dvc = 0, max = 0, temp = 0, mod = 0, k = 0;

    vector<double> u0, u;
    u.resize((finish-start)/h + 1);

    for (int i = 0; i*h <= finish - start; i++) {
        if(i*h + start <= 0) {
            u0.push_back(0);
        } else if (i*h + start <= 0.25) {
            u0.push_back(4*(i*h + start));
        } else {
            u0.push_back(1);
        }
    }

    for (int i = 0; i <= 1/t; i++) {
        U( u0, u, h, t);
        u0 = u;
    }

    k = -1.0;
    for(int i = 0; i <= u.size() - 1; i ++) {
        
        temp = abs(au(k) - u[i]);
        if (temp > Dvc) {
            Dvc = temp;
        }

        temp = abs(u[i]);
        if (temp > max) {
            max = temp;
        }

        Dvl += abs(au(k) - u[i]);
        mod += abs(u[i]);
        k = k + h;
    } 
    dvc = Dvc/max;
    Dvl = h*Dvl;
    dvl = Dvl/h/mod;

    printf("Dvc: %.e\n", Dvc);
    printf("dvc: %.e\n", dvc);
    printf("Dvl: %.e\n", Dvl);
    printf("dvl: %.e\n", dvl);

    k = -1.0;
    for(int i = 0; i < u.size(); i++) {
        fprintf(fp, "%lf %lf %lf\n", i*h + start, u[i], au(k));
        k+=h;
    }

    fclose(fp);
    
    return 0;
}