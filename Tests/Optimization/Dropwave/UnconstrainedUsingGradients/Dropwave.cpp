#include<stdio.h>
#include<math.h>

double DropwaveAdj(double *x, double *xb) {

	xb[0] = (x[0]*(cos(12*sqrt(x[0]*x[0] + x[1]*x[1])) + 1))/pow((x[0]*x[0]/2 + x[1]*x[1]/2 + 2),2) + (12*x[0]*sin(12*sqrt(x[0]*x[0] + x[1]*x[1])))/(sqrt(x[0]*x[0] + x[1]*x[1])*(x[0]*x[0]/2 + x[1]*x[1]/2 + 2));

	xb[1] = (x[1]*(cos(12*sqrt(x[0]*x[0] + x[1]*x[1])) + 1))/pow((x[0]*x[0]/2 + x[1]*x[1]/2 + 2),2) + (12*x[1]*sin(12*sqrt(x[0]*x[0] + x[1]*x[1])))/(sqrt(x[0]*x[0] + x[1]*x[1])*(x[0]*x[0]/2 + x[1]*x[1]/2 + 2));

	return -(1+cos(12*sqrt(x[0]*x[0]+x[1]*x[1])))/(0.5*(x[0]*x[0]+x[1]*x[1])+2);

}

int main(void){

double x[2];
double xb[2];

FILE *inp = fopen("dv.dat","r");
fscanf(inp,"%lf",&x[0]);
fscanf(inp,"%lf",&x[1]);
fclose(inp);

double result = DropwaveAdj(x, xb);
FILE *outp = fopen("objFunVal.dat","w");
fprintf(outp,"Dropwave_function = %15.10f\n",result);
fprintf(outp,"Dropwave_gradient = %15.10f, %15.10f\n",xb[0],xb[1]);
fclose(outp);

return 0;
}
