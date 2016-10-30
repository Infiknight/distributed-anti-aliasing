#include "my_aa.h"
#include <stdio.h>

int main(int argc, char ** argv){
	if(argc == 8)
		return my_mpi_aa(argc, argv);
	else if(argc == 7)
		return serial_aa(argc, argv);
	else{
		printf("wrong num of args, 6 needed for serial and 7 for parallel comp\n");
		return 1;
	}
}
