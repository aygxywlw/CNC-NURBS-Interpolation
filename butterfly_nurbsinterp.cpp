/* Note:Your choice is C IDE */

#include "stdio.h"
#include "math.h"
#include <iostream>
#include "engine.h"
#include <stdlib.h>
#include <string.h>
#pragma comment(lib,"libeng.lib")               
#pragma comment(lib,"libmx.lib")
#pragma comment(lib,"libmat.lib")

#define deltu 0.001
using namespace std;
double mx[1000], my[1000], mz[1000], mv[1000], mhigh[1000], ma[1000], mb[1000], mc[1000];



void deboor(double calresult[2][2],double coeff[][2], double knot[], double weight[], int degree, int l,  double d, int e)
{
	double pointx, pointy, dpointx, dpointy, unpoint, rpointx, rpointy;
	double t1, t2;//coeffa[9][2],weighta[9];
	int k, i, j, m, n;
	
	
	double **coeffa, *weighta, **dcoeffa, *dweighta, **dcoeff, *dweight;
	coeffa = new double*[degree + l];
	dcoeffa = new double*[degree + l];
	dcoeff = new double*[degree + l];
	int n3;
	for (n3 = 0; n3 < degree + l; n3++)
	{
		coeffa[n3] = new double[2];
		dcoeffa[n3] = new double[2];
		dcoeff[n3] = new double[2];
	}
		

	weighta = new double[degree + l];
	dweight = new double[degree + l];
	dweighta = new double[degree + l];
	
	for (m = e - degree; m <= e; m++)
	{
		weighta[m] = weight[m];
		dweight[m] = weight[m];
	
		for (i = 0; i <= 1; i++)
		{
			coeffa[m][i] = coeff[m][i] * weight[m];
			dcoeff[m][i] = coeff[m][i] * weight[m];
		}
	}
	for (k = 1; k <= degree; k++)
	{
		for (j = e; j >= e - degree + k; j--)
		{
			t1 = (knot[j + degree - k + 1] - d) / (knot[j + degree - k + 1] - knot[j]);
			t2 = 1.0 - t1;
			coeffa[j][0] = t1*coeffa[j - 1][0] + t2*coeffa[j][0];
			coeffa[j][1] = t1*coeffa[j - 1][1] + t2*coeffa[j][1];
			weighta[j] = t1*weighta[j - 1] + t2*weighta[j];
		}

	}

	for (n = e - degree + 1; n <= e; n++)
	{
		dweighta[n] = 3 * (dweight[n] - dweight[n - 1]) / (knot[n + 3] - knot[n]);
		for (i = 0; i <= 1; i++)
		{
			dcoeffa[n][i] = 3 * (dcoeff[n][i] - dcoeff[n - 1][i]) / (knot[n + 3] - knot[n]);
		}
	}


	for (k = 1; k < degree; k++)
	{
		for (j = e; j >= e - degree + k + 1; j--)
		{
			t1 = (knot[j + degree - k] - d) / (knot[j + degree - k] - knot[j]);
			t2 = 1.0 - t1;
			dcoeffa[j][0] = t1*dcoeffa[j - 1][0] + t2*dcoeffa[j][0];
			dcoeffa[j][1] = t1*dcoeffa[j - 1][1] + t2*dcoeffa[j][1];
			dweighta[j] = t1*dweighta[j - 1] + t2*dweighta[j];
		}
	}


	
	
	pointx = coeffa[j + 1][0];
	pointy = coeffa[j + 1][1];
	unpoint = weighta[j + 1];
	rpointx = pointx / unpoint;
	rpointy = pointy / unpoint;

	dpointx = (dcoeffa[j + 1][0] - dweighta[j + 1] * rpointx) / unpoint;
	dpointy = (dcoeffa[j + 1][1] - dweighta[j + 1] * rpointy) / unpoint;
	
	

	calresult[0][0] = rpointx;
	calresult[0][1] = rpointy;
	calresult[1][0] = dpointx;
	calresult[1][1] = dpointy;
}


void CalNurbs(double coeff[][2],double knot[],double weight[],int degree,int Section)
{
       
	double u = 0, u1 = 0.0, u2 = 0.0, mm = 0, nn=0;
	double result[2][2];
	double cx=0, cy=0,dcx=coeff[0][0],dcy=coeff[0][1];
	int mici, kk;
	static int countnew = 0;
	double midx, midy,midu,feedlen,midlenx,midleny,high,totallen=0;
	double speederr;
	double lastcx=coeff[0][0], lastcy=coeff[0][1];
	double lastdcx = coeff[0][0], lastdcy = coeff[0][1];
	mici = degree;
	kk = degree;
	Engine *ep;
	mxArray *TX = NULL, *TY = NULL, *TZ = NULL, *TXX = NULL, *TYY = NULL, *TZZ = NULL, *TV = NULL, *TT = NULL;
	


	if (!(ep = engOpen("\0")))
	{
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		//return EXIT_FAILURE;
	}
	TX = mxCreateDoubleMatrix(1, 1000, mxREAL);
	TY = mxCreateDoubleMatrix(1, 1000, mxREAL);
	/*TZ = mxCreateDoubleMatrix(1, 300, mxREAL);*/
	TXX = mxCreateDoubleMatrix(1, 1000, mxREAL);
	TYY = mxCreateDoubleMatrix(1, 1000, mxREAL);
	//TZZ = mxCreateDoubleMatrix(1, 300, mxREAL);
	//TV = mxCreateDoubleMatrix(1, 300, mxREAL);
	TT = mxCreateDoubleMatrix(1, 1000, mxREAL);
	
	while (u <= 1.0)
	{
/*		if (countnew == 0)
		{
			u = 0.0;
		}
		else
		{
			if (countnew == 1)
			{
				u = 0.01;
			}
			else if (countnew == 2)
			{
				u = 2 * u - u1;
			}
			else
				u = 2.5*u - 2 * u1 + 0.5*u2;
		}*/

		//if (u > 1)
			//u = 1;
		if (u <= knot[mici + 1])
		{
			if (knot[mici + 1] > knot[mici])
			{
				deboor(result,coeff, knot, weight, degree, Section, u, kk);
			}
		}
		else
		{
			if (mici<Section + degree)
			{
				mici++;
				kk++;
				if (knot[mici + 1]>knot[mici])
				{
					deboor(result,coeff, knot, weight, degree, Section, u, kk);
				}
			}
			else
				printf("out of range");
		}
		cx = result[0][0];
		cy = result[0][1];
		dcx = dcx+result[1][0]*deltu;
		dcy = dcy+result[1][1]*deltu; 
		printf("speederr=%f\n", sqrt((dcx - lastdcx)*(dcx - lastdcx) + (dcy - lastdcy)*(dcy - lastdcy)));
		feedlen = sqrt((cx- lastcx)*(cx - lastcx) + (cy - lastcy)*(cy - lastcy));//获得进给步长
		//printf("feedlen=%f\n", feedlen);
		midu = mm + (u - mm) / 2;
		/* printf("midu=%f",midu);*/
		if (midu <= knot[mici + 1])
		{
			if (knot[mici + 1] > knot[mici])
			{
				deboor(result,coeff, knot, weight, degree, Section, midu, kk);
			}
		}
		else
		{
			if (mici<Section + degree)
			{
				mici++;
				kk++;
				if (knot[mici + 1]>knot[mici])
				{
					deboor(result,coeff, knot, weight, degree, Section, midu, kk);
				}
			}
			else
				printf("超出范围");
		}
		midx = result[0][0];
		midy = result[0][1];
		
		

		midlenx = (cx + lastcx) / 2.0;
		midleny = (cy + lastcy) / 2.0;
		

		high = sqrt((midx - midlenx)*(midx - midlenx) + (midy - midleny)*(midy - midleny));

		/*u2 = nn;
		u1 = mm;
		nn = u1;*/
		mm = u;
		mx[countnew] = cx;
		my[countnew] = cy;

		ma[countnew] = dcx;
		mb[countnew] = dcy;
		mhigh[countnew] = high;
		printf("high=%f\n", high);
		countnew++;
		u += deltu;

		lastcx = cx;
		lastcy = cy;
		lastdcx = dcx;
		lastdcy = dcy;
		totallen += feedlen;
		
	}
	printf("totallen=%f\n", totallen);
	memcpy((void *)mxGetPr(TX), (void *)mx, 1000 * sizeof(double));
	memcpy((void *)mxGetPr(TY), (void *)my, 1000 * sizeof(double));
	//memcpy((void *)mxGetPr(TZ), (void *)mz, 300 * sizeof(double));
	memcpy((void *)mxGetPr(TXX), (void *)ma, 1000 * sizeof(double));
	memcpy((void *)mxGetPr(TYY), (void *)mb, 1000 * sizeof(double));
	/*memcpy((void *)mxGetPr(TZZ), (void *)mzz, 300 * sizeof(double));
	memcpy((void *)mxGetPr(TV), (void *)mv, 300 * sizeof(double));*/
	memcpy((void *)mxGetPr(TT), (void *)mhigh, 1000 * sizeof(double));
	engPutVariable(ep, "TX", TX);
	engPutVariable(ep, "TY", TY);
	//engPutVariable(ep, "TZ", TZ);
	//engPutVariable(ep, "TXX", TXX);
	//engPutVariable(ep, "TYY", TYY);
	//engPutVariable(ep, "TZZ", TZZ);
	//engPutVariable(ep, "TT", TT);
	//engPutVariable(ep, "TV", TV);
	//engEvalString(ep, "plot(TX,TT,'.b');");
	//engEvalString(ep, "hold on;");
	//engEvalString(ep, "figure;");

	//engEvalString(ep, "plot(TX,TY,'.r');");
	//engEvalString(ep, "plot3(TXX,TYY,TZZ,'.r');");

	//engEvalString(ep, "title('surface figure');");
	//engEvalString(ep, "xlabel('x axis');");
	//engEvalString(ep, "ylabel('y axis');");
	//engEvalString(ep, "zlabel('z axis');");
	//engEvalString(ep, "grid on;");
	//engEvalString(ep, "hold on;");

	fgetc(stdin);
	mxDestroyArray(TX);
	mxDestroyArray(TY);
	/*mxDestroyArray(TZ);*/
	mxDestroyArray(TXX);
	mxDestroyArray(TYY);
	//mxDestroyArray(TZZ);
	mxDestroyArray(TT);
	//mxDestroyArray(TV);*/

	engEvalString(ep, "close;");
	engClose(ep);


}



int main()
{  
    double coeff[9][2]={40,240,100,140,170,260,230,320,300,280,360,200,420,120,480,220,540,300}; /*控制点*/
    double knot[21]= {0,0,0,0,0.16,0.32,0.5,0.67,0.83,1,1,1,1}; /*节点*/
    double weight[17]= {1,1,1,1,1,1,1,1,1}; /*权值 */
	
	//double coeff[51][2] = { 54.493000,52.139000,55.507000,52.139000,56.082000,49.615000,56.780000,44.971000,69.575000,51.358000,77.786000,58.573000,90.526000,67.081000,105.973000,63.801000,100.400000,47.326000,94.567000,39.913000,92.369000,30.485000,83.440000,33.757000,91.892000,28.509000,89.444000,20.393000,83.218000,15.446000,87.621000,4.830000,80.945000,9.267000,79.834000,14.535000,76.074000,8.522000,70.183000,12.550000,64.171000,16.865000,59.993000,22.122000,55.680000,36.359000,56.925000,24.995000,59.765000,19.828000,54.493000,14.940000,49.220000,19.828000,52.060000,24.994000,53.305000,36.359000,48.992000,22.122000,44.814000,16.865000,38.802000,12.551000,32.911000,8.521000,29.152000,14.535000,28.040000,9.267000,21.364000,4.830000,25.768000,15.447000,19.539000,20.391000,17.097000,28.512000,25.537000,33.750000,16.602000,30.496000,14.199000,39.803000,8.668000,47.408000,3.000000,63.794000,18.465000,67.084000,31.197000,58.572000,39.411000,51.358000,52.204000,44.971000,52.904000,49.614000,53.478000,52.139000,54.492000,52.139000 }; /*控制点*/
	//double knot[55] = { 0.000000,0.000000,0.000000,0.000000,0.008286,0.014978,0.036118,0.085467,0.129349,0.150871,0.193075,0.227259,0.243467,0.256080,0.269242,0.288858,0.316987,0.331643,0.348163,0.355261,0.364853,0.383666,0.400499,0.426851,0.451038,0.465994,0.489084,0.499973,0.510862,0.533954,0.548910,0.573096,0.599447,0.616280,0.635094,0.644687,0.651784,0.668304,0.682958,0.711087,0.730703,0.743865,0.756479,0.772923,0.806926,0.849130,0.870652,0.914534,0.963883,0.985023,0.991714,1.000000,1.000000,1.000000,1.000000 }; /*节点*/
	//double weight[51] = { 1.000000,1.000000,1.000000,1.200000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,2.000000,1.000000,1.000000,5.000000,3.000000,1.000000,1.100000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.100000,1.000000,3.000000,5.000000,1.000000,1.000000,2.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.200000,1.000000,1.000000,1.000000 }; /*权值 */
	


    CalNurbs(coeff,knot,weight,3,48);
   
    getchar();
    
    return 0;
}
