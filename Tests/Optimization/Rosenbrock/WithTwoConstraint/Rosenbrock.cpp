#include<stdio.h>
#include<math.h>
#include <unistd.h>
#include <iostream>
#include <armadillo>
using namespace arma;
using namespace std;

double Rosenbrock10(double *x) {

	double temp = 0.0;

	for(int i=0; i<9; i++){

		temp+= 100.0*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]) + (1-x[i])*(1-x[i]);
	}

	return temp;
}

double constraint(double *x){

	 double sum = 0;

	 for(int i=0; i<10; i++){
        sum = sum + x[i];
	 }

	 return sum;
}


int main(void){

  double x[10];

  FILE *inp = fopen("dv.dat","r");

  for (unsigned int k = 0; k < sizeof(x); k++ ){
		fscanf(inp,"%lf",&x[k]);
   }
  fclose(inp);

  double result = Rosenbrock10(x);
  double constraintValue = constraint(x);

  std::ofstream obj_value;
  obj_value.open("objFunVal.dat");
  obj_value << result << std::endl;
  obj_value.close();


  std::ofstream const_value;
  const_value.open("conFunVal.dat");
  const_value << constraintValue << std::endl;
  const_value.close();

   std::ofstream const_value1;
   const_value1.open("conFunVal1.dat");
   const_value1 << constraintValue+1 << std::endl;
   const_value1.close();


  return 0;
}

