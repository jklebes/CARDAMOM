#include "testC.h"
#include <stdio.h>

void main(){
	printf("running\n");
	printf("%f\n",PI);
	printf("%d\n",X);
  C_initialize_stresstest_circle();
  C_initialize_model(1);
}
