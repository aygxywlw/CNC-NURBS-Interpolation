
#include "stdio.h"
#include <math.h>
#include <iostream>
#include "engine.h"
#include <stdlib.h>
#include <string.h>
#pragma comment(lib,"libeng.lib")               
#pragma comment(lib,"libmx.lib")
#pragma comment(lib,"libmat.lib")
using namespace std;


/*double coeff[9][9][3] = {{{20,20,20},{50,20,40},{80,20,60},{110,20,80},{140,20,100},{170,20,80},{200,20,60},{230,20,40},{260,20,20}},
{ { 20,40,40 },{ 50,40,60 },{ 80,40,80 },{ 110,40,100 },{ 140,40,120 },{ 170,40,100 },{ 200,40,80 },{ 230,40,60 },{ 260,40,40 } },
{ { 20,60,60 },{ 50,60,80 },{ 80,60,100 },{ 110,60,120 },{ 140,60,140 },{ 170,60,120 },{ 200,60,100 },{ 230,60,80     },{ 260,60,60 } },
{ { 20,80,80 },{ 50,80,100 },{ 80,80,120 },{ 110,80,140 },{ 140,80,160 },{ 170,80,140 },{ 200,80,120 },{ 230,80,100 },{ 260,80,80 } },
{ { 20,100,100 },{ 50,100,120 },{ 80,100,140 },{ 110,100,160 },{ 140,100,180 },{ 170,100,160 },{ 200,100,140 },{ 230,100,120 },{ 260,100,100 } },
{ { 20,120,80 },{ 50,120,100 },{ 80,120,120 },{ 110,120,140 },{ 140,120,160 },{ 170,120,140 },{ 200,120,120 },{ 230,120,100 },{ 260,120,80 } },
{ { 20,140,60 },{ 50,140,80 },{ 80,140,100 },{ 110,140,120 },{ 140,140,140 },{ 170,140,140 },{ 200,140,100 },{ 230,140,80 },{ 260,140,60 } },
{ { 20,160,40 },{ 50,160,60 },{ 80,160,80 },{ 110,160,100 },{ 140,160,120 },{ 170,160,100 },{ 200,160,80 },{ 230,160,60 },{ 260,160,40 } },
{ { 20,180,20 },{ 50,180,40 },{ 80,180,60 },{ 110,180,80 },{ 140,180,100 },{ 170,180,80 },{ 200,180,60 },{ 230,180,40 },{ 260,180,20 } } }; //control points
/* double coeff[9][9][3]={{{30,20,20},{30,60,50},{30,100,90},{30,130,120},{30,170,160},{30,200,140},{30,260,100},{30,300,80}},
 {{60,20,50},{60,60,90},{60,100,140},{60,130,170},{60,180,200},{60,210,160},{60,250,120},{60,290,100}},
{{100,20,80},{100,60,100},{100,100,150},{100,130,180},{100,200,210},{100,240,180},{100,300,150},{100,330,110}},
{{140,20,100},{140,60,100},{140,100,160},{140,130,200},{140,170,220},{140,220,180},{140,260,150},{140,300,120}},
{{200,20,120},{205,60,160},{200,100,200},{195,130,260},{200,200,290},{190,240,260},{200,300,210},{205,350,160}},
{{240,20,100},{240,60,140},{245,100,190},{240,130,220},{235,200,270},{240,240,240},{245,300,200},{240,360,160}},
{{300,20,90},{295,60,130},{300,100,170},{305,130,200},{300,180,240},{300,220,210},{310,260,170},{300,300,140}},
{{350,20,70},{360,60,110},{350,100,150},{340,130,180},{350,200,220},{350,260,190},{360,300,160},{350,350,120}}
 }; //control points

double coeff[9][9][3] = {{{30,20,20},{30,60,50},{30,100,90},{30,130,120},{30,170,160},{30,200,140},{30,260,100},{30,300,80},{30,340,60}},
{{60,20,50},{60,60,90},{60,100,140},{60,130,170},{60,180,200},{60,210,160},{60,250,120},{60,290,100},{60,330,70}},
{{100,20,80},{100,60,100},{100,100,150},{100,130,180},{100,200,210},{100,240,180},{100,300,150},{100,330,110},{100,360,80}},
{{140,20,100},{140,60,100},{140,100,160},{140,130,200},{140,170,220},{140,220,180},{140,260,150},{140,300,120},{140,350,90}},
{{200,20,120},{205,60,160},{200,100,200},{195,130,260},{200,200,290},{190,240,260},{200,300,210},{205,350,160},{210,400,100}},
{{240,20,100},{240,60,140},{245,100,190},{240,130,220},{235,200,270},{240,240,240},{245,300,200},{240,360,160},{240,410,110}},
{{300,20,90},{295,60,130},{300,100,170},{305,130,200},{300,180,240},{300,220,210},{310,260,170},{300,300,140},{300,360,100}},
{{350,20,70},{360,60,110},{350,100,150},{340,130,180},{350,200,220},{350,260,190},{360,300,160},{350,350,120},{350,400,80}},
{{400,20,50},{400,60,90},{395,90,110},{400,130,140},{405,200,190},{400,240,150},{400,300,120},{405,350,100},{400,390,70}}}; ////control points
double knotu[13] = { 0,0,0,0,0.17,0.33,0.50,0.67,0.83,1,1,1,1 };  //knot vector
double knotw[13] = { 0,0,0,0,0.17,0.33,0.50,0.67,0.83,1,1,1,1 };
double weight[9][9] = { 1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1 }; //weigths
//double va[81] = { 40,25,30,25,20,30,20,25,25,60,65,60,55,60,65,60,65,55,100,95,100,105,100,105,100,95,100,140,135,140,145,140,140,140,145,140,200,205,200,195,200,190,200,205,210,240,240,245,240,235,240,245,240,240,300,295,300,305,300,300,310,300,300,350,360,350,340,350,350,360,350,350,400,400,395,400,405,400,400,405,400 };
//double vb[81] = { 20,60,100,130,170,200,260,300,340,30,70,110,150,180,210,250,290,330,20,60,100,140,200,240,300,330,360,25,55,90,130,170,220,260,300,350,30,60,110,160,200,240,300,350,400,20,50,100,140,200,240,300,360,410,25,60,110,140,180,220,260,300,360,30,70,120,160,200,260,300,350,400,20,50,90,150,200,240,300,350,390 };
//double vc[81] = { 20,50,90,120,160,140,100,80,60,50,90,140,170,200,160,120,100,70,80,100,150,180,210,180,150,110,80,100,100,160,200,220,180,150,120,90,120,160,200,260,290,260,210,160,100,100,140,190,220,270,240,200,160,110,90,130,170,200,240,210,170,140,100,70,110,150,180,220,190,160,120,80,50,90,110,140,190,150,120,100,70 };
//double mesha[9][9] = { { 40,25,30,25,20,30,20,25,25 },{ 60,65,60,55,60,65,60,65,55 },{ 100,95,100,105,100,105,100,95,100 },{ 140,135,140,145,140,140,140,145,140 },{ 200,205,200,195,200,190,200,205,210 },{ 240,240,245,240,235,240,245,240,240 },{ 300,295,300,305,300,300,310,300,300 },{ 350,360,350,340,350,350,360,350,350 },{ 400,400,395,400,405,400,400,405,400 } };
//double meshb[9][9] = { { 20,60,100,130,170,200,260,300,340 },{ 30,70,110,150,180,210,250,290,330 },{ 20,60,100,140,200,240,300,330,360 },{ 25,55,90,130,170,220,260,300,350 },{ 30,60,110,160,200,240,300,350,400 },{ 20,50,100,140,200,240,300,360,410 },{ 25,60,110,140,180,220,260,300,360 },{ 30,70,120,160,200,260,300,350,400 },{ 20,50,90,150,200,240,300,350,390 } };
//double meshc[9][9] = { { 20,50,90,120,160,140,100,80,60 },{ 50,90,140,170,200,160,120,100,70 },{ 80,100,150,180,210,180,150,110,80 },{ 100,100,160,200,220,180,150,120,90 },{ 120,160,200,260,290,260,210,160,100 },{ 100,140,190,220,270,240,200,160,110 },{ 90,130,170,200,240,210,170,140,100 },{ 70,110,150,180,220,190,160,120,80 },{ 50,90,110,140,190,150,120,100,70 } };

double x, y, z, ux, uy, uz;
double ucpnx, ucpny, ucpnz, uweight;   
double wcpnx, wcpny, wcpnz, wweight; 
double px1, py1, pz1; 
double d = 5.0;   
double radius = 50;  //radius of the cutter
double xx, yy, zz;
double dpointx, dpointy, dpointz; 
double dux, duy, duz, dwx, dwy, dwz;   
double ddpx, ddpy, ddpz;      
double ddwx, ddwy, ddwz;
double setl =3;       //the preset step length
double acc = 0;
double msetl = 3;
double setd; //tool path distance  
double seth = 0.02;//preset residual error
double high = 0.02;     //preset chord error
int   wmici = 3;   //NURBS order
int   mici = 3;    
int daoflag = 0;   
int countnew = 0;     
double u, u0;
double mm = 0;
double nn = 0;
double u1 = 0;
double u2 = 0;
int countw = 0;    
double mxx[300];
double myy[300];
double mzz[300];
double w;
double rl, velocity;
double deltw, deltw0;
						
				
long tick0 = 1, tick1 = 1, tick2 = 1;
double deriux0, deriuy0, deriuz0;
double setl0 = 0, setl1 = 0, setl2 = 0;
double jiasuleiji;

long tick00 = 1, tick01 = 1, tick02 = 1;
double setl00 = 0, setl01 = 0, setl02 = 0;
double jiansuleiji;
double lastcurlong = 0;


double w0;
double ll, vl;
double ucoeff[9][3]; 
double wcoeff[9][3]; 
double uwei[9];      
double wwei[9];    
int flag;          
int degree = 3, l = 6, dense = 50; 
double vecx, vecy, vecz; 
double anglea, anglec;  
void concurve();
void vercon();  
void direcw();  
void prepare();  
void XYZpoint(); 

Engine *ep;
mxArray *TX = NULL, *TY = NULL, *TZ = NULL, *TXX = NULL, *TYY = NULL, *TZZ = NULL, *TV = NULL, *TT = NULL;
double mx[30000], my[30000], mz[30000], mv[30000], mt[30000], ma[30000], mb[30000], mc[30000];
double meshx[9][9],meshy[9][9],meshz[9][9];
void deboor(double d, int e, double p[9][3], double wp[9], double knot[13]);
void udbor(double d, int e, int n);  
void udborf(double d, int e, int n);
void wdbor(double d, int e, int n);
void uccpnt(double w, int wmici, int k, int hh);  
void wccpnt(double u, int mici, int k);  
void main()
{

	double deriwx0, deriwy0, deriwz0;





	deriwx0 = 3 * weight[0][1] * (coeff[0][1][0] - coeff[0][0][0]) / weight[0][0];
	deriwy0 = 3 * weight[0][1] * (coeff[0][1][1] - coeff[0][0][1]) / weight[0][0];
	deriwz0 = 3 * weight[0][1] * (coeff[0][1][2] - coeff[0][0][2]) / weight[0][0];

	//w0 = 100 * knotw[4] / sqrt(deriwx0*deriwx0 + deriwy0*deriwy0 + deriwz0*deriwz0);

	direcw();

}

void direcw()
{
	int k;//kkk;
	int nofc;
	double wm = 0;
	double wn = 0;
	double w1 = 0;
	double w2 = 0;

	// double deriux0,deriuy0,deriuz0; 

	w0 = 0.01;

	//w = 0.1;
	if (!(ep = engOpen("\0")))
	{
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		//return EXIT_FAILURE;
	}
	TX = mxCreateDoubleMatrix(1, 30000, mxREAL);
	TY = mxCreateDoubleMatrix(1, 30000, mxREAL);
	TZ = mxCreateDoubleMatrix(1, 30000, mxREAL);
	TXX = mxCreateDoubleMatrix(1, 30000, mxREAL);
	TYY = mxCreateDoubleMatrix(1, 30000, mxREAL);
	TZZ = mxCreateDoubleMatrix(1, 30000, mxREAL);
	TV = mxCreateDoubleMatrix(1, 30000, mxREAL);
	TT = mxCreateDoubleMatrix(1, 30000, mxREAL);


	//memcpy((void *)mxGetPr(TXX), (void *)meshx, 81 * sizeof(double));
	//memcpy((void *)mxGetPr(TYY), (void *)meshy, 81 * sizeof(double));
	//memcpy((void *)mxGetPr(TZZ), (void *)meshz, 81 * sizeof(double));
	//engPutVariable(ep, "TXX", TXX);
	//engPutVariable(ep, "TYY", TYY);
	//engPutVariable(ep, "TZZ", TZZ);
	//engEvalString(ep, "figure;");
	//engEvalString(ep, "mesh(TXX,TYY,TZZ);");
	//engEvalString(ep, "axis([20,600,20,600,20,400]);");
	//engEvalString(ep, "hold on;");
	//engEvalString(ep, "figure;");
	//engEvalString(ep, "axis([20,600,20,600,20,400]);");*/
	/*  for(kkk=0;kkk<=80;kkk++)
	{
	ma[kkk]=va[kkk];
	mb[kkk]=vb[kkk];
	mc[kkk]=vc[kkk];
	// engEvalString(ep, "plot3(TXX,TYY,TZZ,'*r');");
	//engEvalString(ep, "hold on;");
	//cout<<ma[kkk]<<endl;
	//}
	memcpy((void *)mxGetPr(TXX), (void *)ma, 300*sizeof(double));
	memcpy((void *)mxGetPr(TYY), (void *)mb, 300*sizeof(double));
	memcpy((void *)mxGetPr(TZZ), (void *)mc, 300*sizeof(double));
	engPutVariable(ep, "TXX", TXX);
	engPutVariable(ep, "TYY", TYY);
	engPutVariable(ep, "TZZ", TZZ);
	engEvalString(ep, "axis(20,600,20,600,20,400);");
	//for(kkk=0;kkk<=80;kkk++){
	engEvalString(ep, "plot3(TXX,TYY,TZZ,'or');");

	engEvalString(ep, "hold on;");
	//engEvalString(ep, "mesh(va,vb,vc);");
	}*/
	//engEvalString(ep, "mesh(mesha,meshb,meshc);");



	//for(nofc=0; nofc<2; nofc++)
	while (w <=1.0)
	{
		printf("w=%f", w);
		lastcurlong = 0;
		rl = 0;
		velocity = 0;
		if (w>knotw[wmici + 1])
			wmici++;
		//cout << "wmici=" << wmici<< endl;

		for (k = 0; k <= 3; k++)
			uccpnt(w, wmici, k, countw % 2);

		mm = 0;
		nn = 0;
		u = 0;
		mici = 3;

		px1 = ucoeff[0][0];
		py1 = ucoeff[0][1];
		pz1 = ucoeff[0][2];



		deriux0 = 3 * uwei[1] * (ucoeff[1][0] - ucoeff[0][0]) / uwei[0];
		deriuy0 = 3 * uwei[1] * (ucoeff[1][1] - ucoeff[0][1]) / uwei[0];
		deriuz0 = 3 * uwei[1] * (ucoeff[1][2] - ucoeff[0][2]) / uwei[0];
		countnew = 0;
		
		// u0=0.01*knotu[4]*setl/sqrt(deriux0*deriux0+deriuy0*deriuy0+deriuz0*deriuz0);
		//printf("u0=%f ",u0);
		/*if (countw == 0)
		{
			w = w0;
		}
		else
		{
			if (countw == 1)
			{
				w = 2 * w - w1;
			}
			else
				w = 2.5*w - 2 * w1 + 0.5*w2;
		}
		w2 = wn;
		w1 = wm;
		wn = w1;
		wm = w;*/
		//deltw0 = 1;
		while (u<1.0)
		{
			XYZpoint();

		}
		printf("w=%.16lf \n", w);
		w = w +0.005;
		
		//printf("lastcurlong=%f \n", lastcurlong);



		memcpy((void *)mxGetPr(TX), (void *)mx, 30000 * sizeof(double));
		memcpy((void *)mxGetPr(TY), (void *)my, 30000 * sizeof(double));
		memcpy((void *)mxGetPr(TZ), (void *)mz, 30000 * sizeof(double));
		memcpy((void *)mxGetPr(TXX), (void *)mxx, 30000*sizeof(double));
		memcpy((void *)mxGetPr(TYY), (void *)myy, 30000*sizeof(double));
		memcpy((void *)mxGetPr(TZZ), (void *)mzz, 30000*sizeof(double));
		memcpy((void *)mxGetPr(TV), (void *)mv, 30000* sizeof(double));
		memcpy((void *)mxGetPr(TT), (void *)mt, 30000 * sizeof(double));
		engPutVariable(ep, "TX", TX);
		engPutVariable(ep, "TY", TY);
		engPutVariable(ep, "TZ", TZ);
		engPutVariable(ep, "TXX", TXX);
		engPutVariable(ep, "TYY", TYY);
		engPutVariable(ep, "TZZ", TZZ);
		engPutVariable(ep, "TT", TT);
		//engEvalString(ep, "plot(TT,TX,'b.');");
		//engPutVariable(ep, "TV", TV);
		engEvalString(ep, "plot3(TX,TY,TZ,'.b');");
		engEvalString(ep, "hold on;");
		engEvalString(ep, "plot3(TXX,TYY,TZZ,'.r');");
		engEvalString(ep, "hold on;");
		//engEvalString(ep, "figure;");
		//if(w>0.999){*/
		/*engEvalString(ep, "subplot(5,1,1)");
		engEvalString(ep, "plot(TT,TX,'b.');");
		engEvalString(ep, "subplot(5,1,2)");
		engEvalString(ep, "plot(TT,TY,'r.');");
		engEvalString(ep, "subplot(5,1,3)");
		engEvalString(ep, "plot(TT,TZ,'c.');");
		engEvalString(ep, "subplot(5,1,4)");
		engEvalString(ep, "plot(TT,TXX,'m.');");
		engEvalString(ep, "subplot(5,1,5)");
		engEvalString(ep, "plot(TT,TYY,'k.');");
		//engEvalString(ep, "hold on;");//}
		//engEvalString(ep, "subplot(2,1,2)");
	   // engEvalString(ep, "plot(TT,TX,'.r');");
		//engEvalString(ep, "hold on;");
		//engEvalString(ep, "plot(TT,TY,'.r');");
		//engEvalString(ep, "hold on;");
		//engEvalString(ep, "plot(TT,TZ,'+b');");
		//engEvalString(ep, "hold on;");
		//engEvalString(ep, "plot(TT,TXX,'.r');");
		//engEvalString(ep, "hold on;");
		//engEvalString(ep, "plot(TT,TYY,'squarer');");
		//engEvalString(ep, "hold on;");
		//engEvalString(ep, "plot(TT,TZZ,'^r');");
		//engEvalString(ep, "hold on;");*/
		//w = w + deltw0;
		countw++;

	}
	//printf("deltw0=%f \n ",deltw0);
	//engEvalString(ep, "hold off;");
	//engPutVariable(ep, "TXX", TXX);
	//     engPutVariable(ep, "TYY", TYY);
	//     engPutVariable(ep, "TZZ", TZZ);
	//engEvalString(ep, "plot3(TXX,TYY,TZZ,'*r');");
	//engEvalString(ep, "title('surface figure');");
	//engEvalString(ep, "xlabel('x axis');");
	//engEvalString(ep, "ylabel('y axis');");
	//engEvalString(ep, "zlabel('z axis');");
	//engEvalString(ep, "grid on;");
	//engEvalString(ep, "hold on;");

	fgetc(stdin);
	mxDestroyArray(TX);
	mxDestroyArray(TY);
	mxDestroyArray(TZ);
	mxDestroyArray(TXX);
	mxDestroyArray(TYY);
	mxDestroyArray(TZZ);
	mxDestroyArray(TT);
	mxDestroyArray(TV);

	engEvalString(ep, "close;");
	engClose(ep);

}
/*uccpnt is for getting control points reversely of the u-direction NURBS curve of the NURBS surface*/
void uccpnt(double w, int wmici, int k, int hh)
{


	if (hh == 1)
		udborf(w, wmici, k);
	else
		udbor(w, wmici, k);
	ucoeff[k][0] = ucpnx;
	ucoeff[k][1] = ucpny;
	ucoeff[k][2] = ucpnz;
	uwei[k] = uweight;

}

/*wccpnt is for getting control points reversely of the w-direction NURBS curve of the NURBS surface*/
void wccpnt(double u, int mici, int k)
{

	int i;
	for (i = k - 3; i <= k; i++)
	{
		wdbor(u, mici, i);
		wcoeff[i][0] = wcpnx;
		wcoeff[i][1] = wcpny;
		wcoeff[i][2] = wcpnz;
		wwei[i] = wweight;
	}
}
/*concurve()
{
int kk;
int ii;

for(kk=degree;kk<l+degree;kk++)
{/*printf("kk=%d,\n",kk);*/
/* 	if(knot[kk+1]>knot[kk])
{
for(ii=0;ii<dense;ii++)
{
u=knot[kk]+ii*(knot[kk+1]-knot[kk])/dense;
deboor(u,kk);
mx=ux;
my=uy;
circle(mx,my,2);

}
}
}
}

/*udbor is used to calculate the NURBS points for u-direction NURBS curve.*/
void udbor(double d, int e, int n)
{


	int k, i, j, m;
	double t1, t2;
	double coeffa[9][3];
	double weighta[9];
	/*printf("e=%d,\n",e); */
	for (m = e - degree; m <= e; m++)
	{
		weighta[m] = weight[n][m];
		for (i = 0; i <= 2; i++)
		{
			coeffa[m][i] = coeff[n][m][i] * weight[n][m];
		}
	}
	for (k = 1; k <= degree; k++)
	{
		for (j = e; j >= e - degree + k; j--)
		{
			t1 = (knotw[j + degree - k + 1] - d) / (knotw[j + degree - k + 1] - knotw[j]);
			t2 = 1.0 - t1;
			coeffa[j][0] = t1*coeffa[j - 1][0] + t2*coeffa[j][0];
			coeffa[j][1] = t1*coeffa[j - 1][1] + t2*coeffa[j][1];
			coeffa[j][2] = t1*coeffa[j - 1][2] + t2*coeffa[j][2];
			weighta[j] = t1*weighta[j - 1] + t2*weighta[j];
		}

	}
	ucpnx = coeffa[j + 1][0];
	ucpny = coeffa[j + 1][1];
	ucpnz = coeffa[j + 1][2];
	uweight = weighta[j + 1];


}

/*udborf is used to calculate the NURBS points reversely.*/
void udborf(double d, int e, int n)
{


	int k, i, j, m, r;
	double t1, t2;
	double coeffa[9][3];
	double weighta[9];
	/*printf("e=%d,\n",e); */
	r = 8 - n;
	for (m = e - degree; m <= e; m++)
	{
		weighta[m] = weight[r][m];
		for (i = 0; i <= 2; i++)
		{
			coeffa[m][i] = coeff[r][m][i] * weight[r][m];
		}
	}
	for (k = 1; k <= degree; k++)
	{
		for (j = e; j >= e - degree + k; j--)
		{
			t1 = (knotw[j + degree - k + 1] - d) / (knotw[j + degree - k + 1] - knotw[j]);
			t2 = 1.0 - t1;
			coeffa[j][0] = t1*coeffa[j - 1][0] + t2*coeffa[j][0];
			coeffa[j][1] = t1*coeffa[j - 1][1] + t2*coeffa[j][1];
			coeffa[j][2] = t1*coeffa[j - 1][2] + t2*coeffa[j][2];
			weighta[j] = t1*weighta[j - 1] + t2*weighta[j];
		}

	}
	ucpnx = coeffa[j + 1][0];
	ucpny = coeffa[j + 1][1];
	ucpnz = coeffa[j + 1][2];
	uweight = weighta[j + 1];


}
/*udbor is used to calculate the NURBS points for w-direction NURBS curve.*/
void wdbor(double d, int e, int n)
{


	int k, i, j, m;
	double t1, t2;
	double coeffa[9][3];
	double weighta[9];
	/*printf("e=%d,\n",e); */
	for (m = e - degree; m <= e; m++)
	{
		weighta[m] = weight[m][n];
		for (i = 0; i <= 2; i++)
		{
			coeffa[m][i] = coeff[m][n][i] * weight[m][n];
		}
	}
	for (k = 1; k <= degree; k++)
	{
		for (j = e; j >= e - degree + k; j--)
		{
			t1 = (knotu[j + degree - k + 1] - d) / (knotu[j + degree - k + 1] - knotu[j]);
			t2 = 1.0 - t1;
			coeffa[j][0] = t1*coeffa[j - 1][0] + t2*coeffa[j][0];
			coeffa[j][1] = t1*coeffa[j - 1][1] + t2*coeffa[j][1];
			coeffa[j][2] = t1*coeffa[j - 1][2] + t2*coeffa[j][2];
			weighta[j] = t1*weighta[j - 1] + t2*weighta[j];
		}

	}
	wcpnx = coeffa[j + 1][0];
	wcpny = coeffa[j + 1][1];
	wcpnz = coeffa[j + 1][2];
	wweight = weighta[j + 1];


}

/*deboor is used to calculate the NURBS point and its first and second derivatives.
parameter: 
d: the parameter value
e: NURBS order
p:control points
wp:weights
knot:knotvcetor*/

void deboor(double d, int e, double p[9][3], double wp[9], double knot[13])
{
	double pointx, pointy, pointz, unpoint;

	// double rpointx,rpointy,rpointz;	
	int    k, i, j, m, n, s;
	double t1, t2;
	double coeffb[9][3];
	double dcoeffb[9][3];
	double weightb[9];
	double dweightb[9];
	double dweight[9];
	double dcoeff[9][3];
	double ddcoeff[9][3];
	double ddweight[9];
	double ddcoefb[9][3];
	double ddweighb[9];
	for (m = e - degree; m <= e; m++)
	{
		weightb[m] = wp[m];
		dweight[m] = wp[m];
		ddweight[m] = wp[m];
		/*dweighta[m]=3*(weight[m]-weight[m-1])/(knot[m+3]-knot[m]); */
		for (i = 0; i <= 2; i++)
		{
			coeffb[m][i] = p[m][i];//*wp[m];
			dcoeff[m][i] = p[m][i];//*wp[m]; 
			ddcoeff[m][i] = p[m][i];//*wp[m];
									/*dcoeffa[m][i]=3*(coeff[m][i]*weight[m]-coeff[m-1][i]*weight[m])/(knot[m+3]-knot[m]);*/
		}
	}

	for (k = 1; k <= degree; k++)
	{
		for (j = e; j >= e - degree + k; j--)
		{
			t1 = (knot[j + degree - k + 1] - d) / (knot[j + degree - k + 1] - knot[j]);
			t2 = 1.0 - t1;
			coeffb[j][0] = t1*coeffb[j - 1][0] + t2*coeffb[j][0];
			coeffb[j][1] = t1*coeffb[j - 1][1] + t2*coeffb[j][1];
			coeffb[j][2] = t1*coeffb[j - 1][2] + t2*coeffb[j][2];
			weightb[j] = t1*weightb[j - 1] + t2*weightb[j];
		}

	}


	for (n = e - degree + 1; n <= e; n++)
	{
		dweightb[n] = 3 * (dweight[n] - dweight[n - 1]) / (knot[n + 3] - knot[n]);
		for (i = 0; i <= 2; i++)
		{
			dcoeffb[n][i] = 3 * (dcoeff[n][i] - dcoeff[n - 1][i]) / (knot[n + 3] - knot[n]);
		}
	}


	for (k = 1; k<degree; k++)
	{
		for (j = e; j >= e - degree + k + 1; j--)
		{
			t1 = (knot[j + degree - k] - d) / (knot[j + degree - k] - knot[j]);
			t2 = 1.0 - t1;
			dcoeffb[j][0] = t1*dcoeffb[j - 1][0] + t2*dcoeffb[j][0];
			dcoeffb[j][1] = t1*dcoeffb[j - 1][1] + t2*dcoeffb[j][1];
			dcoeffb[j][2] = t1*dcoeffb[j - 1][2] + t2*dcoeffb[j][2];
			dweightb[j] = t1*dweightb[j - 1] + t2*dweightb[j];
		}
	}

	for (s = e - degree + 2; s <= e; s++)
	{
		ddweighb[s] = 6 * (ddweight[s] - 2 * ddweight[s - 1] + ddweight[s - 2]) / ((knot[s + 3] - knot[s])*(knot[s + 2] - knot[s]));
		for (i = 0; i <= 2; i++)
		{
			//ddcoefb[s][i]=6*(ddcoeff[s][i]-2*ddcoeff[s-1][i]+ddcoeff[s-2][i])/((knot[s+3]-knot[s])*(knot[s+2]-knot[s]));
			ddcoefb[s][i] = 6 * ((ddcoeff[s][i] - ddcoeff[s - 1][i]) / ((knot[s + 3] - knot[s])*(knot[s + 2] - knot[s])) - (ddcoeff[s - 1][i] - ddcoeff[s - 2][i]) / ((knot[s + 2] - knot[s - 1])*(knot[s + 2] - knot[s])));
		}
	}

	
	for (k = 1; k<degree - 1; k++)
	{
		for (j = e; j >= e - degree + k + 2; j--)
		{
			t1 = (knot[j + degree - k - 1] - d) / (knot[j + degree - k - 1] - knot[j]);
			t2 = 1.0 - t1;
			ddcoefb[j][0] = t1*ddcoefb[j - 1][0] + t2*ddcoefb[j][0];
			ddcoefb[j][1] = t1*ddcoefb[j - 1][1] + t2*ddcoefb[j][1];
			ddcoefb[j][2] = t1*ddcoefb[j - 1][2] + t2*ddcoefb[j][2];
			ddweighb[j] = t1*ddweighb[j - 1] + t2*ddweighb[j];
		}
	}
	pointx = coeffb[j + 1][0];
	pointy = coeffb[j + 1][1];
	pointz = coeffb[j + 1][2];
	unpoint = weightb[j + 1];
	ux = pointx / unpoint;
	uy = pointy / unpoint;
	uz = pointz / unpoint;
	dpointx = (dcoeffb[j + 1][0] - dweightb[j + 1] * ux) / unpoint;
	dpointy = (dcoeffb[j + 1][1] - dweightb[j + 1] * uy) / unpoint;
	dpointz = (dcoeffb[j + 1][2] - dweightb[j + 1] * uz) / unpoint;
	ddpx = (ddcoefb[j + 1][0] - ddweighb[j + 1] * ux - 2 * dpointx*dweightb[j + 1]) / unpoint;
	ddpy = (ddcoefb[j + 1][1] - ddweighb[j + 1] * uy - 2 * dpointy*dweightb[j + 1]) / unpoint;
	ddpz = (ddcoefb[j + 1][2] - ddweighb[j + 1] * uz - 2 * dpointz*dweightb[j + 1]) / unpoint;
}



void XYZpoint()
{
	
	double qulv;
	double qulvr;
	double qulvx, qulvy, qulvz;
	double qux, quy, quz;
	double test;
	double cosa, sina;
	double midx, midy, midz;
	double midlx, midly, midlz;
	double midu;
	double h, ml,rh;
	double ddux, dduy, dduz;
	double actualh;
	double wx, wy, wz;
	double wwx, wwy, wwz;
	double asidestep; 


//for velocity and acceleration schematic
	//if(w==0)





/*if (lastcurlong < 140)
{
	if(setl<0.2*msetl)
	{
		setl=(1.0/650)*tick0*tick0;
		acc = (2.0 / 65) * tick0;
		tick0++;
		if(setl>=0.2*msetl)
		{
			setl=0.2*msetl;
			acc = 400.0 / 650;
			setl0=setl;
		}
	//	printf("tick0=%d\n", tick0);
	//	printf("setl=%f\n", setl);
	}
	else if(setl<0.8*msetl)
	{
		setl=setl0+0.08*tick1;
		acc = 400.0 / 650;
		tick1++;
		if(setl>=0.8*msetl)
		{
			setl=0.8*msetl;
			setl1=setl;
		}
		//printf("tick1=%d \n", tick1);
		//printf("setl=%f\n", setl);
	}
	else if(setl<msetl)
	{
		setl=setl1+ 0.061 *tick2-(1.0/650)*tick2*tick2;
		acc = 400.0 / 650 - (2.0 / 65) * tick2;
		tick2++;
		if(setl>msetl)
			setl=msetl;
		//printf("tick2=%d \n", tick2);
		//printf("setl=%f\n", setl);
		//printf("lastcurlong1=%f \n", lastcurlong);
	}
	else
	{ 
		setl = msetl;
		acc = 0;
	}
		
jiasuleiji+=rl;
}
//u0=knotu[4]*setl/sqrt(deriux0*deriux0+deriuy0*deriuy0+deriuz0*deriuz0);




if (lastcurlong > (408.402881 - 89.9))
{
	if(setl>0.8*msetl)
	{
		setl=msetl-(1.0/650)*tick00*tick00;
		acc = -(2.0 / 65) * tick00;
		tick00++;
		//printf("tick00=%d" ,tick00);
		if(setl<=0.8*msetl)
		{
			setl=0.8*msetl;
			acc = -400.0 / 650;
			setl00=setl;
		}
	}
	else if(setl>0.2*msetl)
	{
		setl=setl00-0.08*tick01;
		acc = -400.0 / 650;
		tick01++;
		//printf("tick01=%d" ,tick01);
		if(setl<=0.2*msetl)
		{
			setl=0.2*msetl;
			setl01=setl;
		}
	}
	else if(setl>0)
	{
		setl=setl01-0.05780*tick02+(1.0/650)*tick02*tick02;
		acc = -400.0 / 650 + (2.0 / 65) * tick02;;
		tick02++;
		//printf("lastcurlong=%f \n", lastcurlong);
//printf("tick02=%d ",tick02);
// if(setl<=0)
//  {
//	   u=1.0;
//	   rl=0;
// }
//if(setl>msetl)
// setl=msetl;
	}
	else
		setl=0;
}
// printf("lastcurlong=%f \n",lastcurlong);
 printf("setl=%f \n",setl);
//printf("rl=%f \n",ml);

 //u0=knotu[4]*setl/sqrt(deriux0*deriux0+deriuy0*deriuy0+deriuz0*deriuz0);
//printf("setl=%f ",setl);
//printf("uu=%f ",u0);
//printf("tick0=%d ",tick0);
//printf("tick1=%d ",tick1);
//printf("tick2=%d ",tick2);




//for deceleration
/*if (w>0.999)
{
if(lastcurlong>(408.897653-149.432667)){
if(setl>0.8*msetl)
{
setl=msetl-(1.0/400)*tick00*tick00;
tick00++;
printf("tick00=%d" ,tick00);
if(setl<=0.8*msetl){
setl=0.8*msetl;
setl00=setl;
}

}
else if(setl>0.2*msetl)
{

setl=setl00-0.15*tick01;
tick01++;
printf("tick01=%d" ,tick01);
if(setl<=0.2*msetl){
setl=0.2*msetl;
setl01=setl;
}

}
else if(setl>=0)
{

setl=setl01-0.1*tick02+(1.0/416.32)*tick02*tick02;
tick02++;
printf("tick02=%d ",tick02);
// if(setl<=0)
//  {
//	   u=1.0;
//	   rl=0;
// }
//if(setl>msetl)
// setl=msetl;
}
else
setl=0;
}
//printf("lastcurlong=%f \n",lastcurlong);
//printf("setl=%f \n",setl);
//printf("rl=%f \n",ml);
}*/
u0 = knotu[4] * setl / sqrt(deriux0*deriux0 + deriuy0*deriuy0 + deriuz0*deriuz0);
//printf("setl=%f ",setl);
//printf("uu=%f ",u0);
//printf("tick0=%d ",tick0);
//printf("tick1=%d ",tick1);
//printf("tick2=%d ",tick2);
//*/





	if (countnew == 0)
	{
		u = u0;
	}
	else
	{
		if (countnew == 1)
		{
			u = 2 * u - u1;
		}
		else
			u = 2.5*u - 2 * u1 + 0.5*u2;
	}
	printf("u=%f\n ",u);
	if (u>1)
		u = 1;
	if (u <= knotu[mici + 1])
	{
		if (knotu[mici + 1]>knotu[mici])
		{
			deboor(u, mici, ucoeff, uwei, knotu);
			wccpnt(u, mici, wmici);
		}
	}
	else
	{
		if (mici<l + degree)
		{
			mici++;
			uccpnt(w, wmici, mici, countw % 2);
			flag = 1;
			if (knotu[mici + 1]>knotu[mici])
			{
				deboor(u, mici, ucoeff, uwei, knotu);
				wccpnt(u, mici, wmici);
			}
		}
		else
			printf("out of range");
	}
	ddux = ddpx;
	dduy = ddpy;
	dduz = ddpz;
	x = ux;
	y = uy;
	z = uz;
	dux = dpointx;
	duy = dpointy;
	duz = dpointz;
	ml = sqrt((x - px1)*(x - px1) + (y - py1)*(y - py1) + (z - pz1)*(z - pz1));

	midu = mm + (u - mm) / 2;
	/* printf("midu=%f",midu);*/
	if (midu <= knotu[mici])
	{
		deboor(midu, mici - 1, ucoeff, uwei, knotu);
	}
	else
	{
		deboor(midu, mici, ucoeff, uwei, knotu);
	}
	midx = ux;
	midy = uy;
	midz = uz;



	midlx = (x + px1) / 2.0;
	midly = (y + py1) / 2.0;
	midlz = (z + pz1) / 2.0;

	h = sqrt((midx - midlx)*(midx - midlx) + (midy - midly)*(midy - midly) + (midz - midlz)*(midz - midlz));
	/*printf("h=%f\n",h);*/

	vl = ml*sqrt(high / h);    /*expected step*/
							   /*printf("vl=%f\n",vl);*/
	if (setl<vl)
		ll = setl;
	else
		ll = vl;

	if (fabs(ml - ll) / ll>0.01)
	{
		u = (ll / ml)*u + (1 - ll / ml)*mm;
		if (u>1)
			u = 1;
		if (flag == 1)
		{
			if (u <= knotu[mici])
			{
				deboor(u, mici - 1, ucoeff, uwei, knotu);
				wccpnt(u, mici - 1, wmici);
			}
			else
			{
				deboor(u, mici, ucoeff, uwei, knotu);
				wccpnt(u, mici, wmici);
			}
		}
		else
		{
			if (u <= knotu[mici + 1])
			{
				deboor(u, mici, ucoeff, uwei, knotu);
				wccpnt(u, mici, wmici);
			}
			else
			{
				deboor(u, mici + 1, ucoeff, uwei, knotu);
				wccpnt(u, mici + 1, wmici);
			}
		}

		x = ux;
		y = uy;
		z = uz;
		dux = dpointx;
		duy = dpointy;
		duz = dpointz;
		ddux = ddpx;
		dduy = ddpy;
		dduz = ddpz;
	}
	flag = 0;
	rl = sqrt((x - px1)*(x - px1) + (y - py1)*(y - py1) + (z - pz1)*(z - pz1));

	//actual chord error
	midu = mm + (u - mm) / 2;
	//printf("rl=%f",rl);
	if (midu <= knotu[mici])
	{
		deboor(midu, mici - 1, ucoeff, uwei, knotu);
	}
	else
	{
		deboor(midu, mici, ucoeff, uwei, knotu);
	}
	midx = ux;
	midy = uy;
	midz = uz;

	

	midlx = (x + px1) / 2.0;
	midly = (y + py1) / 2.0;
	midlz = (z + pz1) / 2.0;

	rh = sqrt((midx - midlx) * (midx - midlx) + (midy - midly) * (midy - midly) + (midz - midlz) * (midz - midlz));
	//printf("mu=%f ",u);
	//printf("ml=%f ",ml);
	//printf("rl=%f \n",rl);
	//  printf("h=%f\n",h);
	deboor(w, wmici, wcoeff, wwei, knotw);
	wx = ux;
	wy = uy;
	wz = uz;
	dwx = dpointx;
	dwy = dpointy; 
	dwz = dpointz;
	ddwx = ddpx;
	ddwy = ddpy;  
	ddwz = ddpz;
	qulv = sqrt((dwy*ddwz - dwz*ddwy)*(dwy*ddwz - dwz*ddwy) + (dwz*ddwx - dwx*ddwz)*(dwz*ddwx - dwx*ddwz) + (dwx*ddwy - dwy*ddwx)*(dwx*ddwy - dwy*ddwx)) / sqrt((dwx*dwx + dwy*dwy + dwz*dwz)*(dwx*dwx + dwy*dwy + dwz*dwz)*(dwx*dwx + dwy*dwy + dwz*dwz));
	//qulv=sqrt((duy*dduz-duz*dduy)*(duy*dduz-duz*dduy)+(duz*ddux-dux*dduz)*(duz*ddux-dux*dduz)+(dux*dduy-duy*ddux)*(dux*dduy-duy*ddux))/sqrt((dux*dux+duy*duy+duz*duz)*(dux*dux+duy*duy+duz*duz)*(dux*dux+duy*duy+duz*duz));
	qulvr = 1.0 / qulv;
	qulvx = dwy*ddwz - dwz*ddwy;
	qulvy = dwz*ddwx - dwx*ddwz;
	qulvz = dwx*ddwy - dwy*ddwx;
	qux = (dwy*qulvz - dwz*qulvy) / sqrt((dwy*qulvz - dwz*qulvy)*(dwy*qulvz - dwz*qulvy) + (dwz*qulvx - dwx*qulvz)*(dwz*qulvx - dwx*qulvz) + (dwx*qulvy - dwy*qulvx)*(dwx*qulvy - dwy*qulvx));
	quy = (dwz*qulvx - dwx*qulvz) / sqrt((dwy*qulvz - dwz*qulvy)*(dwy*qulvz - dwz*qulvy) + (dwz*qulvx - dwx*qulvz)*(dwz*qulvx - dwx*qulvz) + (dwx*qulvy - dwy*qulvx)*(dwx*qulvy - dwy*qulvx));
	quz = (dwx*qulvy - dwy*qulvx) / sqrt((dwy*qulvz - dwz*qulvy)*(dwy*qulvz - dwz*qulvy) + (dwz*qulvx - dwx*qulvz)*(dwz*qulvx - dwx*qulvz) + (dwx*qulvy - dwy*qulvx)*(dwx*qulvy - dwy*qulvx));
	
	cosa = ((qulvr + radius)*(qulvr + radius) + (qulvr + seth)*(qulvr + seth) - radius*radius) / (2 * (qulvr + radius)*(qulvr + seth));
	//printf("cosa=%f\n  ",qulvr);
	sina = sqrt(1 - cosa*cosa);
	setd = 2 * qulvr*sina;

	if (countw % 2 == 0)
	{
		vecx = (duy*dwz - dwy*duz) / sqrt((duy*dwz - dwy*duz)*(duy*dwz - dwy*duz) + (duz*dwx - dwz*dux)*(duz*dwx - dwz*dux) + (dux*dwy - dwx*duy)*(dux*dwy - dwx*duy));
		vecy = (duz*dwx - dwz*dux) / sqrt((duy*dwz - dwy*duz)*(duy*dwz - dwy*duz) + (duz*dwx - dwz*dux)*(duz*dwx - dwz*dux) + (dux*dwy - dwx*duy)*(dux*dwy - dwx*duy));
		vecz = (dux*dwy - dwx*duy) / sqrt((duy*dwz - dwy*duz)*(duy*dwz - dwy*duz) + (duz*dwx - dwz*dux)*(duz*dwx - dwz*dux) + (dux*dwy - dwx*duy)*(dux*dwy - dwx*duy));
	}
	else
	{
		vecx = -(duy*dwz - dwy*duz) / sqrt((duy*dwz - dwy*duz)*(duy*dwz - dwy*duz) + (duz*dwx - dwz*dux)*(duz*dwx - dwz*dux) + (dux*dwy - dwx*duy)*(dux*dwy - dwx*duy));
		vecy = -(duz*dwx - dwz*dux) / sqrt((duy*dwz - dwy*duz)*(duy*dwz - dwy*duz) + (duz*dwx - dwz*dux)*(duz*dwx - dwz*dux) + (dux*dwy - dwx*duy)*(dux*dwy - dwx*duy));
		vecz = -(dux*dwy - dwx*duy) / sqrt((duy*dwz - dwy*duz)*(duy*dwz - dwy*duz) + (duz*dwx - dwz*dux)*(duz*dwx - dwz*dux) + (dux*dwy - dwx*duy)*(dux*dwy - dwx*duy));
	}
	anglea = acos(vecz);
	anglec = atan(-vecx / vecy);
	test = qux*vecx + quy*vecy + quz*vecz;
	// printf("test=%f ",test);
	
	if (qulvr>500 * radius)
	{
		setd = 2 * sqrt(radius*radius - (radius - seth)*(radius - seth));
		//actualh = radius - sqrt(radius * radius - 0.25 * setd * setd);
	}
	else
	{
		if (test <= 0)
		{
			cosa = ((qulvr - radius)*(qulvr - radius) + (qulvr - seth)*(qulvr - seth) - radius*radius) / (2 * (qulvr - radius)*(qulvr - seth));
			
			sina = sqrt(1 - cosa*cosa);
			setd = 2 * qulvr*sina;
			//actualh =qulvr- (qulvr - radius) * cosa - sqrt(radius * radius - ((qulvr - radius) * (qulvr - radius) * sina * sina));
		}
		else
		{
			cosa = ((qulvr + radius)*(qulvr + radius) + (qulvr + seth)*(qulvr + seth) - radius*radius) / (2 * (qulvr + radius)*(qulvr + seth));
			
			sina = sqrt(1 - cosa*cosa);
			setd = 2 * qulvr*sina;
			//actualh = (qulvr + radius) * cosa - sqrt(radius * radius - ((qulvr + radius) * (qulvr + radius) * sina * sina))-qulvr;
		}
	}
	 deltw=setd/sqrt(dwx*dwx+dwy*dwy+dwz*dwz)-setd*setd*(dwx*ddwx+dwy*ddwy+dwz*ddwz)/(2*(dwx*dwx+dwy*dwy+dwz*dwz)*(dwx*dwx+dwy*dwy+dwz*dwz));   
	//deltw = setd / sqrt(dwx*dwx + dwy*dwy + dwz*dwz);
	 //printf("setd=%.16lf \n", setd);
	if (deltw<deltw0)
		deltw0 = deltw;
	
	//printf("detlw0=%.16lf \n  ", deltw0); 

	//w = w + deltw0;
	
	deboor(w + deltw0, wmici, wcoeff, wwei, knotw);//w+ 0.000445517
	wwx = ux;
	wwy = uy;
	wwz = uz;
	asidestep = sqrt((wwx - wx) * (wwx - wx) + (wwy - wy) * (wwy - wy) + (wwz - wz) * (wwz - wz));
	if (qulvr > 500 * radius)
	{
		actualh = radius - sqrt(radius * radius - 0.25 * asidestep * asidestep);
	}
	else
	{
		sina = (0.5 * asidestep )/ qulvr;
		//printf("qulvr=%f\n  ", qulvr);
		cosa = sqrt(1 - sina * sina);
		if (test <= 0)
		{
			actualh =qulvr- (qulvr - radius) * cosa - sqrt(radius * radius - ((qulvr - radius) * (qulvr - radius) * sina * sina));
		}
		else
		{
			actualh = (qulvr + radius) * cosa - sqrt(radius * radius - ((qulvr + radius) * (qulvr + radius) * sina * sina))-qulvr;
		}
	}
	
	
	//printf("qulvr=%f \n", qulvr);
	
	//printf("dux=%f ",dux/sqrt(dux*dux+duy*duy+duz*duz)); 
	//printf("duy=%f ",duy/sqrt(dux*dux+duy*duy+duz*duz));
	//printf("duz=%f ",duz/sqrt(dux*dux+duy*duy+duz*duz));
	//printf("dwx=%f ",ddwx); 
	//printf("dwy=%f ",ddwy);
	//printf("dwz=%f \n",ddwz);
	/*printf("x=%f ",x);
	printf("y=%f ",y);
	printf("z=%f ",z);
	printf("a=%f ",anglea);
	printf("c=%f\n",anglec);*/
	if (daoflag == 0)
	{
		xx = x + d*vecx;
		yy = y + d*vecy;
		zz = z + d*vecz;
	}
	else
	{
		xx = x - d*vecx;
		yy = y - d*vecy;
		zz = z - d*vecz;
	}
	mv[countnew] = rl;
	velocity = velocity +rl;
	mt[countnew] = velocity;
	//printf("mv[%d]=%f \n",countnew,mv[countnew]);
	//printf("mt[%d]=%f ",countnew,mt[countnew]);
	printf("x=%f   ", x);
	printf("y=%f   ", y);
	printf("z=%f \n", z);
	printf("xx=%f  ", xx);
	printf("yy=%f  ", yy);
	printf("zz=%f \n", zz);
	
	mx[countnew] =x;
	my[countnew] = y;
	mz[countnew] = z;
	mxx[countnew] = xx;
	myy[countnew] = yy;
	mzz[countnew] = zz;
	u2 = nn;
	u1 = mm;
	nn = u1;
	mm = u;


	px1 = x;
	py1 = y;
	pz1 = z;
	countnew++; 
	lastcurlong += rl;
	//printf("anglea=%f \n",anglea);
	// printf("anglec=%f \n",anglec);

}




