#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void rdf_(char*,float*);
void rdi_(char*,int*);

void rdf_(char *str, float *v) {
	sscanf(str,"%f",v);
}

void rdi_(char *str, int *v) {
	sscanf(str,"%d",v);
}
