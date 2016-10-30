#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <omp.h>

static const int lower_bound_convergence= 0.000001;

void s_update(
	int i,
	int no_channels,
	int ch,
	int height,
	int width,
	unsigned char * arr,
	unsigned char * arr_new)
{
	if(i == 0) //TOP LEFT
		arr_new[ch + no_channels * i]= (9 * arr[ch + no_channels * (i)] + 3 * arr[ch + no_channels * (i + 1)] + 3 * arr[ch + no_channels * (i + width)] + arr[ch + no_channels * (i + width + 1)]) / 16;
	else if(i == width - 1) //TOP RIGHT
		arr_new[ch + no_channels * i]= (9 * arr[ch + no_channels * (i)] + 3 * arr[ch + no_channels * (i - 1)] + 3 * arr[ch + no_channels * (i + width)] + arr[ch + no_channels * (i + width - 1)]) / 16;
	else if(i == (height - 1) * width) //BOTTOM LEFT
		arr_new[ch + no_channels * i]= (9 * arr[ch + no_channels * (i)] + 3 * arr[ch + no_channels * (i + 1)] + 3 * arr[ch + no_channels * (i - width)] + arr[ch + no_channels * (i - width + 1)]) / 16;
	else if(i == height * width - 1) //BOTTOM RIGHT
		arr_new[ch + no_channels * i]= (9 * arr[ch + no_channels * (i)] + 3 * arr[ch + no_channels * (i - 1)] + 3 * arr[ch + no_channels * (i - width)] + arr[ch + no_channels * (i - width - 1)]) / 16;
	else if(i < width) //TOP
		arr_new[ch + no_channels * i]= (6 * arr[ch + no_channels * (i)] + 3 * arr[ch + no_channels * (i - 1)] + 3 * arr[ch + no_channels * (i + 1)] + 2 * arr[ch + no_channels * (i + width)] + arr[ch + no_channels * (i + width - 1)] + arr[ch + no_channels * (i + width + 1)]) / 16;
	else if(i >= (height - 1) * width) //BOTTOM
		arr_new[ch + no_channels * i]= (6 * arr[ch + no_channels * (i)] + 3 * arr[ch + no_channels * (i - 1)] + 3 * arr[ch + no_channels * (i + 1)] + 2 * arr[ch + no_channels * (i - width)] + arr[ch + no_channels * (i - width - 1)] + arr[ch + no_channels * (i - width + 1)]) / 16;
	else if((i + 1) % width == 0) //RIGHT
		arr_new[ch + no_channels * i]= (6 * arr[ch + no_channels * (i)] + 3 * arr[ch + no_channels * (i - width)] + 3 * arr[ch + no_channels * (i + width)] + 2 * arr[ch + no_channels * (i - 1)] + arr[ch + no_channels * (i + width - 1)] + arr[ch + no_channels * (i - width - 1)]) / 16;
	else if(i % width == 0) //LEFT
		arr_new[ch + no_channels * i]= (6 * arr[ch + no_channels * (i)] + 3 * arr[ch + no_channels * (i - width)] + 3 * arr[ch + no_channels * (i + width)] + 2 * arr[ch + no_channels * (i + 1)] + arr[ch + no_channels * (i + width + 1)] + arr[ch + no_channels * (i - width + 1)]) / 16;
	else
		arr_new[ch + no_channels * i]= (4 * arr[ch + no_channels * (i)] + 2 * (arr[ch + no_channels * (i - 1)] + arr[ch + no_channels * (i + 1)] + arr[ch + no_channels * (i - width)] + arr[ch + no_channels * (i + width)]) + arr[ch + no_channels * (i - width - 1)] + arr[ch + no_channels * (i - width + 1)] + arr[ch + no_channels * (i + width - 1)] + arr[ch + no_channels * (i + width + 1)]) / 16;
}

int serial_aa(int argc, char ** argv){
	int maxthreads= omp_get_max_threads();	//number of threads on each CPU
	FILE * inFile;
	FILE * outFile;
	struct timeval start, stop;
	inFile= fopen(argv[1], "r");
	int height= atoi(argv[3]);
	int width= atoi(argv[4]);
	int no_channels= atoi(argv[6]);
	int no_iterations= atoi(argv[5]);
	int counter_1,
		bytes_read,
		i, k,
		ch;
	double sum_of_differences;
	counter_1= 0;
	unsigned char * arr= malloc( height * width * no_channels);
	unsigned char * arr_new= malloc( height * width * no_channels);
	unsigned char * temp_arr;
	while(counter_1 < (height * width) * no_channels){
		bytes_read= fread(arr + counter_1, 1, height * width * no_channels - counter_1, inFile);
		if(bytes_read == -1){
			fprintf(stderr, "%s\n", strerror(errno));
			return 1;
		}
		counter_1+= bytes_read;
	}
	fclose(inFile);
	gettimeofday(&start, NULL);
	#pragma omp parallel num_threads(maxthreads) \
		private(k)
	{

	for(k= 0; k < no_iterations; k++){
		#pragma omp for \
			private(ch)
		for(i= 0; i < height * width; i++){
			for(ch= 0; ch < no_channels; ch++){
				s_update(
					i,
					no_channels,
					ch,
					height,
					width,
					arr,
					arr_new);
			}
		}
		sum_of_differences= 0;
		#pragma omp for \
			reduction(+:sum_of_differences)
		for(i= 0; i < height * width * no_channels; i++){
			sum_of_differences+= ((double) abs((int) (arr[i] - arr_new[i]))) / 255;
		}
		#pragma omp single
		{
		sum_of_differences/= width * height;
		printf("Iteration no %d, difference between this and the previous array: %f%%\n", k, sum_of_differences);
		temp_arr= arr;
		arr= arr_new;
		arr_new= temp_arr;
		}//omp single ends
		if(sum_of_differences < lower_bound_convergence)
			break;
	}

	}//omp parallelization ends
	gettimeofday(&stop, NULL);
	printf("Total time excluding time it took to read and print: %lu microsecs\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);
	outFile= fopen(argv[2], "w");
	counter_1= 0;
	while(counter_1 < (height * width) * no_channels){
		bytes_read= fwrite(arr + counter_1, 1, height * width * no_channels - counter_1, outFile);
		if(bytes_read == -1){
			fprintf(stderr, "%s\n", strerror(errno));
			return 1;
		}
		counter_1+= bytes_read;
	}
	fclose(outFile);
	free(arr);
	free(arr_new);
	return 0;
}
