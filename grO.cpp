#include <iostream>
using namespace std;
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

vector<double> RKfunc(vector<double> alpha, int mark, FILE * fp);
vector<double> dfunc( vector<double> alpha, FILE * fp, int pos);
double detfunc( vector<double> X1, vector<double> X2, vector<double> X3);
void func( FILE * fp);

int main()
{
    clock_t start, stop;
    start = clock ();

    FILE * fp;
    fp = fopen("22(10.0).txt", "w");
    func(fp);

    stop = clock();
    printf("Секунд: %lu\n", (stop - start) / CLOCKS_PER_SEC);
    return 0;
}

//метод стрельбы с модификацией Исаева-Сонина;
void func( FILE * fp)
{
    int k = 0, m = 0, iter = 0;
    double det, det1, det2, det3, eps = 1e-10;
    double h1, h2, h3, a1, a2, a3, norma;
    vector<double> alpha;
    vector<double> X, dX1, dX2, dX3;
    vector<double> ddX1, ddX2, ddX3;

    alpha.resize(3);
    X.resize(3);
    dX1.resize(3);
    dX2.resize(3);
    dX3.resize(3);
    ddX1.resize(3);
    ddX2.resize(3);
    ddX3.resize(3);

    alpha[0] = -0.37499998; //начальные приближения;
    alpha[1] = -0.49999998;
    alpha[2] = 4.01111;

    X = RKfunc(alpha, k, fp);

    while(pow(X[0]*X[0] + X[1]*X[1] + X[1]*X[1] ,0.5) > eps)
    {
        iter++;
        dX1 = dfunc(alpha, fp, 0); //первая строка матрицы производных;
        dX2 = dfunc(alpha, fp, 1); //вторая строка матрицы производных;
        dX3 = dfunc(alpha, fp, 2); //третья строка матрицы производных;

        //метод Крамера;
        det = detfunc(dX1, dX2, dX3); 

        ddX1[0] = -X[0];
        ddX1[1] = dX1[1];
        ddX1[2] = dX1[2];

        ddX2[0] = -X[1];
        ddX2[1] = dX2[1];
        ddX2[2] = dX2[2];

        ddX3[0] = -X[2];
        ddX3[1] = dX3[1];
        ddX3[2] = dX3[2];

        det1 = detfunc(ddX1, ddX2, ddX3);

        ddX1[0] = dX1[0];
        ddX1[1] = -X[0];
        ddX1[2] = dX1[2];

        ddX2[0] = dX2[0];
        ddX2[1] = -X[1];
        ddX2[2] = dX2[2];

        ddX3[0] = dX3[0];
        ddX3[1] = -X[2];
        ddX3[2] = dX3[2];

        det2 = detfunc(ddX1, ddX2, ddX3);

        ddX1[0] = dX1[0];
        ddX1[1] = dX1[1];
        ddX1[2] = -X[0];

        ddX2[0] = dX2[0];
        ddX2[1] = dX2[1];
        ddX2[2] = -X[1];

        ddX3[0] = dX3[0];
        ddX3[1] = dX3[1];
        ddX3[2] = -X[2];

        det3 = detfunc(ddX1, ddX2, ddX3);

        h1 = det1/det;
        h2 = det2/det;
        h3 = det3/det;

        a1 = alpha[0];
        a2 = alpha[1];
        a3 = alpha[2];

        ddX1 = RKfunc(alpha, k, fp);

        alpha[0] = a1 + h1; 
        alpha[1] = a2 + h2;
        alpha[2] = a3 + h3;

        X = RKfunc(alpha, k, fp);
        norma = pow(X[0]*X[0] + X[1]*X[1] + X[2]*X[2],0.5 );

        //модификация Исаева-Сонина;
        while (pow(X[0]*X[0] + X[1]*X[1] + X[2]*X[2],0.5 ) > pow(ddX1[0]*ddX1[0] + ddX1[1]*ddX1[1] + ddX1[2]*ddX1[2],0.5 ))
        {
            h1 = h1/2;
            h2 = h2/2;
            h3 = h3/2;

            alpha[0] = a1 + h1; 
            alpha[1] = a2 + h2;
            alpha[2] = a3 + h3;

            X = RKfunc(alpha, k, fp);
            m++;
            if(m > 32)
                break;
        }
        
    }
    k++;
    printf("alpha1: %.15lf\n", alpha[0]);
    printf("alpha2: %.15lf\n", alpha[1]);
    printf("alpha3: %.15lf\n", alpha[2]);
    X = RKfunc(alpha, k, fp);
    printf("Итераций: %d\n", iter);
}

//функция, считающая строку производных;
vector<double> dfunc( vector<double> alpha, FILE * fp, int pos)
{
    double delta = 1e-8;
    vector<double> dXi, dalpha;
    dXi.resize(3);
    dalpha.resize(3);

    dalpha[0] = alpha[0] + delta;
    dalpha[1] = alpha[1];
    dalpha[2] = alpha[2];

    dXi[0] = (RKfunc(dalpha, 0, fp)[pos] - RKfunc(alpha, 0, fp)[pos])/delta;

    dalpha[0] = alpha[0];
    dalpha[1] = alpha[1] + delta;
    dalpha[2] = alpha[2];

    dXi[1] = (RKfunc(dalpha, 0, fp)[pos] - RKfunc(alpha, 0, fp)[pos])/delta;

    dalpha[0] = alpha[0];
    dalpha[1] = alpha[1];
    dalpha[2] = alpha[2] + delta;

    dXi[2] = (RKfunc(dalpha, 0, fp)[pos] - RKfunc(alpha, 0, fp)[pos])/delta;

    return dXi;
}

//функция, считающая определитель 3х3;
double detfunc( vector<double> X1, vector<double> X2, vector<double> X3)
{
    return X1[0]*X2[1]*X3[2]  +  X1[2]*X2[0]*X3[1]  +  X1[1]*X2[2]*X3[0]  -  X1[2]*X2[1]*X3[0]  -  X1[0]*X2[2]*X3[1]  -  X1[1]*X2[0]*X3[2];
}

//метод Рунге-Кутты;
vector<double> RKfunc(vector<double> alpha, int mark, FILE * fp)
{
    double al = 0;
    double T = alpha[2];
    int st = 0;
    double ct = 0, tol = 1e-9;
    double step = 0.1, err = 0;
    double K1 = 0, K2 = 0, K3 = 0, K4 = 0, K0 = 0;

    vector<double> c, b, db, u, du, u0;
    vector<vector<double> > k;
    vector<vector<double> > a;

    a.resize(6);
    c.resize(7);
	b.resize(7);
    db.resize(7);
	k.resize(7);
    
    for(int i = 0; i < 7; i++)
    {
        k[i].resize(5);
    }

    u.resize(5);
    du.resize(5);
    u0.resize(5);

    u0[0] = 0;
    u0[1] = 0;
    u0[2] = alpha[0];
    u0[3] = alpha[1];
    u0[4] = 0;

    u = u0;
    du = u0;

    for(int i = 0; i < 6; i++)
	{
		a[i].resize(i+1);
	}

    //коэффициенты;
    a[0][0]=1.0/5; 

    a[1][0]=3.0/40; 
    a[1][1]=9.0/40; 

    a[2][0]=(44.0/45); 
    a[2][1]=-56.0/15;
    a[2][2]=32.0/9;

    a[3][0]=19372.0/6561;
    a[3][1]=-25360.0/2187;
    a[3][2]=64448.0/6561;
    a[3][3]=-212.0/729;

    a[4][0]=9017.0/3168;
    a[4][1]=-355.0/33;
    a[4][2]=46732.0/5247;
    a[4][3]=49.0/176;
    a[4][4]=-5103.0/18656;

    a[5][0]=35.0/384;
    a[5][1]=0;
    a[5][2]=500.0/1113;
    a[5][3]=125.0/192;
    a[5][4]=-2187.0/6784;
    a[5][5]=11.0/84;


    c[0] = 0; 
    c[1] = 1.0/5; 
    c[2] = 3.0/10; 
    c[3] = 4.0/5; 
    c[4] = 8.0/9;
    c[5] = 1; 
    c[6] = 1;


	b[0] = 35.0/384; 
    b[1] = 0; 
    b[2] = 500.0/1113; 
    b[3] = 125.0/192; 
    b[4] = -2187.0/6784; 
    b[5] = 11.0/84; 
    b[6] = 0;


    db[0] = 5179.0/57600; 
    db[1] = 0; 
    db[2] = 7571.0/16695; 
    db[3] = 393.0/640; 
    db[4] = -92097.0/339200; 
    db[5] = 187.0/2100; 
    db[6] = 1.0/40;

    if(mark == 1)
        fprintf(fp, "%.9lf %.9lf %.9lf %.9lf %.9lf %.9lf\n", ct, u[0], u[1], u[2], u[3], u[4]);

    while(T > ct)
    {
        if(T < ct + step)
            {
                step = T - ct;
            }
        
        for(int i = 1; i < 7; i++)
        {
            k[0][0] = u0[1];
            k[0][1] = u0[3]*(1+al*pow(sin(u0[0]),2));
            k[0][2] = -0.5*pow(u0[3],2)*al*sin(2*u0[0]);
            k[0][3] = -u0[2];
            k[0][4] = pow(u0[3],2)*(1+al*pow(sin(u0[0]),2));
            K0 = u0[0];
            K1 = u0[1];
            K2 = u0[2];
            K3 = u0[3];
            K4 = u0[4];
            for(int j = 1; j <= i; j++)
            {
                K0 += step*a[i-1][j-1]*k[j-1][0];
                K1 += step*a[i-1][j-1]*k[j-1][1];
                K2 += step*a[i-1][j-1]*k[j-1][2];
                K3 += step*a[i-1][j-1]*k[j-1][3];
                K4 += step*a[i-1][j-1]*k[j-1][4];
            }
            k[i][0] = K1;
            k[i][1] = K3*(1+al*pow(sin(K0),2));
            k[i][2] = -0.5*pow(K3,2)*al*sin(2*K0);
            k[i][3] = -K2;
            k[i][4] = pow(K3,2)*(1+al*pow(sin(K0),2));
        }
        
        u[0] = u0[0] + step*(b[0]*k[0][0] + b[1]*k[1][0] + b[2]*k[2][0] + b[3]*k[3][0] + b[4]*k[4][0] + b[5]*k[5][0] + b[6]*k[6][0]);
        u[1] = u0[1] + step*(b[0]*k[0][1] + b[1]*k[1][1] + b[2]*k[2][1] + b[3]*k[3][1] + b[4]*k[4][1] + b[5]*k[5][1] + b[6]*k[6][1]);
        u[2] = u0[2] + step*(b[0]*k[0][2] + b[1]*k[1][2] + b[2]*k[2][2] + b[3]*k[3][2] + b[4]*k[4][2] + b[5]*k[5][2] + b[6]*k[6][2]);
        u[3] = u0[3] + step*(b[0]*k[0][3] + b[1]*k[1][3] + b[2]*k[2][3] + b[3]*k[3][3] + b[4]*k[4][3] + b[5]*k[5][3] + b[6]*k[6][3]);
        u[4] = u0[4] + step*(b[0]*k[0][4] + b[1]*k[1][4] + b[2]*k[2][4] + b[3]*k[3][4] + b[4]*k[4][4] + b[5]*k[5][4] + b[6]*k[6][4]);

        du[0] = u0[0] + step*(db[0]*k[0][0] + db[1]*k[1][0] + db[2]*k[2][0] + db[3]*k[3][0] + db[4]*k[4][0] + db[5]*k[5][0] + db[6]*k[6][0]);
        du[1] = u0[1] + step*(db[0]*k[0][1] + db[1]*k[1][1] + db[2]*k[2][1] + db[3]*k[3][1] + db[4]*k[4][1] + db[5]*k[5][1] + db[6]*k[6][1]);
        du[2] = u0[2] + step*(db[0]*k[0][2] + db[1]*k[1][2] + db[2]*k[2][2] + db[3]*k[3][2] + db[4]*k[4][2] + db[5]*k[5][2] + db[6]*k[6][2]);
        du[3] = u0[3] + step*(db[0]*k[0][3] + db[1]*k[1][3] + db[2]*k[2][3] + db[3]*k[3][3] + db[4]*k[4][3] + db[5]*k[5][3] + db[6]*k[6][3]);
        du[4] = u0[4] + step*(db[0]*k[0][4] + db[1]*k[1][4] + db[2]*k[2][4] + db[3]*k[3][4] + db[4]*k[4][4] + db[5]*k[5][4] + db[6]*k[6][4]);

        err = pow( pow(u[0] - du[0], 2) + pow(u[1] - du[1], 2) + pow(u[2] - du[2], 2) + pow(u[3] - du[3], 2) + pow(u[4] - du[4], 2)   , 0.5);
        
        if (err <= tol) 
        {
            u0 = u;
            ct += step; 
            //проверка на "финальность" решения;
            if (mark == 1)
            {
                fprintf(fp, "%.9lf %.9lf %.9lf %.9lf %.9lf %.9lf\n", ct, u[0], u[1], u[2], u[3], u[4]);
            }
        }

        step = step*min( 1.3, max( 0.7,0.98*pow(tol/err,1.0/6)));
    }

    //возвращение вектор-функции невязок;
    vector<double> result;
    result.resize(3);
    result[0] = u[0];
    result[1] = u[1] - 1;
    result[2] = u[4] - 1;

    if(mark == 1)
    {
        printf("dx(T): %.e\n", fabs(u[0]));
        printf("dy(T): %.e\n", fabs(u[1] - 1.0));
        printf("dpx(T): %.e\n", fabs(u[2] + 0.375));
        printf("dpy(T): %.e\n", fabs(u[3] - 1));
        printf("dz(T): %.e\n", fabs(u[4] - 1));

        printf("x(T): %.15lf\n", fabs(u[0]));
        printf("y(T): %.15lf\n", u[1]);
        printf("px(T): %.15lf\n", u[2]);
        printf("py(T): %.15lf\n", u[3]);
        printf("z(T): %.15lf\n", u[4]);
    }

    return result;
}
