#include <iostream>
using namespace std;
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

double au (double x) {
    if (x <= -3.0/8) {
        return 0;
    } else {
        return 1;
    }
}

void U (vector<double> u0, vector<double>& u, double h, double t) {
    u[0] = 0;
    u[u.size() - 1] = 1;

    for(int i = 1; i< u.size() - 1; i++) {
        u[i] = (u0[i]/t + 0.25*pow(u0[i+1],2)/h - 0.25*pow(u0[i],2)/h - u0[i-1]*u[i-1]/2/h)/(1.0/t - 1.0/2/h*u0[i-1]);
    }
}

int main() {

    FILE * fp;
    FILE * fp1;
    fp = fopen("ves3.txt", "w");
    fp1 = fopen("vspom.txt", "w");

    double h = 0.1, t = 0.001, k = 0;
    double start = -1, finish = 1;

    double dvl = 0, Dvl = 0, dvc = 0, Dvc = 0, max = 0, temp = 0, mod = 0;

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
    for(int i = 0; i < u.size(); i++) {
        fprintf(fp, "%lf %lf %lf\n", i*h + start, u[i], au(k));
        k+=h;
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

    
    for(k = -1.0; k <= 1; k+=0.0001) {
        fprintf(fp1, "%lf %lf\n" ,k , au(k));
    }

    fclose(fp);
    

    return 0;
}