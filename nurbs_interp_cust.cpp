/* Note:Your choice is C IDE */
/*画一条NURBS曲线*/
#include "stdio.h"
#include "math.h"
#include <graphics.h>
#include <time.h>
#include <dos.h>

void deboor(double calresult[2][2],double coeff[][2], double knot[], double weight[], int degree, int l, int dense, double d, int e)
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
	//circle(rpointx, rpointy, 1);
	
	
	
	
	calresult[0][0] = rpointx;
	calresult[0][1] = rpointy;
	calresult[1][0] = dpointx;
	calresult[1][1] = dpointy;
	
}


void CalNurbs(double coeff[][2],double knot[],double weight[],int degree,int Section,int dense)
{
       
	double u; 
	double positionx = coeff[0][0];
	double positiony = coeff[0][1];
	double result[2][2];
    for(int kk=degree; kk<Section+degree; kk++)
    {
        if(knot[kk+1]>knot[kk])
        {
            for(int ii=0; ii<dense; ii++)
            {
                u=knot[kk]+ii*(knot[kk+1]-knot[kk])/dense;
                deboor(result,coeff,knot,weight,degree,Section,dense,u,kk);
				//circle(result[0][0], result[0][1], 2);
				positionx = positionx + result[1][0] * 0.00167;
				positiony = positiony + result[1][1] * 0.00167;
				circle(positionx, positiony, 2);
            }
			
        }
    }
}



int main()
{
    
    clock_t start, end;  /*时间函数*/
  
    double coeff[9][2]={40,240,100,140,170,260,230,320,300,280,360,200,420,120,480,220,540,300}; /*控制点*/
    double knot[21]= {0,0,0,0,0.16,0.32,0.5,0.67,0.83,1,1,1,1}; /*节点*/
    double weight[17]= {1,1,1,1,1,1,1,1,1}; /*权值 */
	
    initgraph(800,600);
    circle(40,240,5);     /*画出控制点*/
    circle(100,140,5);
    circle(170,260,5);
    circle(230,320,5);
    circle(300,280,5);
    circle(360,200,5);
    circle(420,120,5);
    circle(480,220,5);
    circle(540,300,5);
	

    moveto(40,240);
    start=clock();
    CalNurbs(coeff,knot,weight,3,6,100);
    end=clock();
    printf("the time is%f ", (end-start)/CLK_TCK);
    getchar();
     closegraph();

    return 0;
}
