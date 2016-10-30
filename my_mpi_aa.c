#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <math.h>
#include "tools.h"
#include <omp.h>

static const int lower_bound_convergence= 0.000001;

//triple send (one for each intermediate array)
int tri_MPI_Isend( 
	int offset, 
	int count, 
	MPI_Datatype datatype, 
	int dest, 
	int tag, 
	MPI_Comm comm, 
	MPI_Request * request, 
	unsigned short * one, 
	unsigned short * two,
	unsigned short * three)
{
	MPI_Isend(
		one + offset,
		count,
		datatype,
		dest,
		tag,
		comm,
		request);
	MPI_Isend(
		two + offset,
		count,
		datatype,
		dest,
		tag,
		comm,
		request);
	MPI_Isend(
		three + offset,
		count,
		datatype,
		dest,
		tag,
		comm,
		request);
	return 0;
}

//triple receive (one for each intermediate array)
int tri_MPI_Irecv( 
	int offset, 
	int count, 
	MPI_Datatype datatype, 
	int source, 
	int tag, 
	MPI_Comm comm, 
	MPI_Request * request, 
	unsigned short * one, 
	unsigned short * two,
	unsigned short * three)
{
	MPI_Irecv(
		one + offset,
		count,
		datatype,
		source,
		tag,
		comm,
		request);
	MPI_Irecv(
		two + offset,
		count,
		datatype,
		source,
		tag,
		comm,
		request);
	MPI_Irecv(
		three + offset,
		count,
		datatype,
		source,
		tag,
		comm,
		request);
	return 0;
}

int my_mpi_aa(int argc, char ** argv){
	int self_rank,
		comm_size,
		cart_rank;
	int i, j, k,
		indx_1,
		indx_2,
		bytes_written,
		bytes_read,
		counter_1;
	int periods[2]= {1, 1},	//cartesian dimensions are periodic
		reorder= 1,			//cartesian ranks may be a reordering of teh inital ones
		coords[2],
		neighb_coords[8][2],
		neighb_ranks[8];	//ranks of process' neighbors
	unsigned char * arr= NULL;
	unsigned char * temp_ptr;
	unsigned long time_sum= 0;
	double sum_of_differences,
		difference_q;
	char done= 0;
	int maxthreads= omp_get_max_threads();	//number of threads on each CPU
	int threading_type;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &threading_type);
	//printf("threading type: %d\n", threading_type);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &self_rank);
	if(argc != 8){
		if(cart_rank == 0)
			printf("7 arguments needed: input_file output_file height width superblock_edge_size num_of_iterations num_of_colors\n");
		MPI_Finalize();
		return 1;
	}
	MPI_Comm cartcomm;
	MPI_Request reqs[16];
	MPI_Status stats[16];
	int superblock_edge= (int) sqrt(comm_size); //size of superblock's edge in blocks
	int dims[2]= {superblock_edge, superblock_edge};
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);
	MPI_Comm_rank(cartcomm, &cart_rank);
	MPI_Cart_coords(cartcomm, cart_rank, 2, coords);
	int iteration_no= 0;
	
	struct timeval start, stop;

	//calc ranks of neighboors
	for(i= 0; i < 8; i++){
		neighb_coords[i][0]= coords[0];
		neighb_coords[i][1]= coords[1];
	}
	neighb_coords[TOP_LEFT][0]+= -1;
	neighb_coords[TOP_LEFT][1]+= -1;
	neighb_coords[TOP][0]+= -1;
	neighb_coords[TOP_RIGHT][0]+= -1;
	neighb_coords[TOP_RIGHT][1]+= 1;
	neighb_coords[RIGHT][1]+= 1;
	neighb_coords[BOTTOM_RIGHT][0]+= 1;
	neighb_coords[BOTTOM_RIGHT][1]+= 1;
	neighb_coords[BOTTOM][0]+= 1;
	neighb_coords[BOTTOM_LEFT][0]+= 1;
	neighb_coords[BOTTOM_LEFT][1]+= -1;
	neighb_coords[LEFT][1]+= -1;
	for(i= 0; i < 8; i++){
		MPI_Cart_rank(cartcomm, neighb_coords[i], &neighb_ranks[i]);
	}
	//------------------------

	int block_edge= atoi(argv[5]) / superblock_edge; //size of block's edge in pixels
	int height= atoi(argv[3]);
	int width= atoi(argv[4]);
	int height_in_sb= height / (superblock_edge * block_edge);
	int width_in_sb= width / (superblock_edge * block_edge);
	unsigned char * local_arr= malloc( (height * width) / comm_size);
	unsigned char * local_arr_new= malloc( (height * width) / comm_size);	//the updated local_arr (we need both to test for convergence)
	unsigned short * one_arr= calloc( (block_edge + 2) * (block_edge + 2) * height_in_sb * width_in_sb, sizeof(short));	//intermediate arrays to hold multiples of 1, 2 and 4 of the local_arr's elements
	unsigned short * two_arr= calloc( (block_edge + 2) * (block_edge + 2) * height_in_sb * width_in_sb, sizeof(short));	//elements are size short so as to be able to hold the multiples that may not fit
	unsigned short * four_arr= calloc( (block_edge + 2) * (block_edge + 2) * height_in_sb * width_in_sb, sizeof(short));//in char variables

	FILE * outFile;	//output file
	FILE * inFile;	//input file
	MPI_Datatype block_temp, block_type,	//used for scattering and gathering between the reading/writing process and the whole set of processes
		vertical_vector,		//vector sued for transfering vertical sides of blocks (left n right, and for each, all of the sides at once)
		vertical_vector_2,		//like above but minutely shorter in 1 dimension to account for the extremities
		horizontal_vector,		//similarly for up and down
		horizontal_vector_2,	//like above but used in case of extremities
		diagonal_vector,					//for the single elements of up-left, down-right etc
		diagonal_vector_lacking_horizontal,	//shorter on 1 dimension
		diagonal_vector_lacking_vertical,	//shorter in the other dimension
		diagonal_vector_lacking_both;		//shorter in both dimensions
	MPI_Type_vector(block_edge, block_edge, width, MPI_BYTE, &block_temp);
	MPI_Type_create_resized(block_temp, 0, 1, &block_type);
	MPI_Type_commit(&block_type);

	int blocklengths_1[height_in_sb * width_in_sb * block_edge], 
		displacements_1[height_in_sb * width_in_sb * block_edge];
	for(i= 0; i < height_in_sb * width_in_sb; i++){
		for(j= 0; j < block_edge; j++){
			blocklengths_1[j + i * block_edge]= 1;
			displacements_1[j + i * block_edge]= j * (block_edge + 2) + i * (block_edge + 2) * (block_edge + 2);
		}
	}
	MPI_Type_indexed( height_in_sb * width_in_sb * block_edge, blocklengths_1, displacements_1, MPI_SHORT, &vertical_vector);
	MPI_Type_commit(&vertical_vector);
	for(i= 0; i < height_in_sb; i++){
		for(k= 0; k < width_in_sb - 1; k++){
			for(j= 0; j < block_edge; j++){
				displacements_1[j + k * block_edge + i * block_edge * (width_in_sb - 1)]
					= j * (block_edge + 2) + k * (block_edge + 2) * (block_edge + 2) + i * (block_edge + 2) * (block_edge + 2) * width_in_sb;
			}
		}
	}
	MPI_Type_indexed( height_in_sb * (width_in_sb - 1) * block_edge, blocklengths_1, displacements_1, MPI_SHORT, &vertical_vector_2);
	MPI_Type_commit(&vertical_vector_2);

	MPI_Type_vector( height_in_sb * width_in_sb, block_edge, (block_edge + 2) * (block_edge + 2), MPI_SHORT, &horizontal_vector);
	MPI_Type_commit(&horizontal_vector);
	MPI_Type_vector( (height_in_sb - 1) * width_in_sb, block_edge, (block_edge + 2) * (block_edge + 2), MPI_SHORT, &horizontal_vector_2);
	MPI_Type_commit(&horizontal_vector_2);

	MPI_Type_vector( height_in_sb * width_in_sb, 1, (block_edge + 2) * (block_edge +2), MPI_SHORT, &diagonal_vector);
	MPI_Type_commit(&diagonal_vector);
	MPI_Type_vector( (height_in_sb - 1) * width_in_sb, 1, (block_edge + 2) * (block_edge +2), MPI_SHORT, &diagonal_vector_lacking_horizontal);
	MPI_Type_commit(&diagonal_vector_lacking_horizontal);
	for(i= 0; i < height_in_sb; i++){
		for(j= 0; j < width_in_sb - 1; j++){
			displacements_1[j + i * (width_in_sb - 1)]= (j + i * width_in_sb) * (block_edge + 2) * (block_edge + 2);
		}
	}
	MPI_Type_indexed(height_in_sb * (width_in_sb - 1), blocklengths_1, displacements_1, MPI_SHORT, &diagonal_vector_lacking_vertical);
	MPI_Type_commit(&diagonal_vector_lacking_vertical);
	MPI_Type_indexed( (height_in_sb - 1) * (width_in_sb - 1), blocklengths_1, displacements_1, MPI_SHORT, &diagonal_vector_lacking_both);
	MPI_Type_commit(&diagonal_vector_lacking_both);


	int counts[comm_size];
	int displacements[comm_size];
	for(i= 0; i < comm_size; i++){
		counts[i]= 1;
		displacements[i]= (i % (int) sqrt(comm_size)) * block_edge + (i / (int) sqrt(comm_size)) * block_edge *  width;
	}
	//if RGB, seperate colors in 3 files
	int no_channels= atoi(argv[7]);
	FILE * color_file[3];
	char buffy_3[3];
	if(no_channels == 3 && cart_rank == 0){
		printf("constructing intermediate color files\n");
		inFile= fopen(argv[1], "r");
		color_file[0]= fopen(".red_temp", "w");
		color_file[1]= fopen(".green_temp", "w");
		color_file[2]= fopen(".blue_temp", "w");
		for(i= 0; i < width * height; i++){
			if( fread(buffy_3, 1, 3, inFile) != 3){
				fprintf(stderr, "%s %d\n", strerror(errno), __LINE__);
				MPI_Finalize();
				return 1;
			}
			fwrite(buffy_3, 1, 1, color_file[0]);
			fwrite(buffy_3 + 1, 1, 1, color_file[1]);
			fwrite(buffy_3 + 2, 1, 1, color_file[2]);
		}
		fclose(inFile);
		for(i= 0; i < no_channels; i++){
			fclose(color_file[i]);
		}
	}
	//----------------------------------
	int color_iter;
	#pragma omp parallel num_threads(maxthreads) \
		private(color_iter)
	{
	
	//outer iteration
	for(color_iter= 0; color_iter < no_channels; color_iter++){

	#pragma omp single
	{

	//read and scatter the contents of the input file
	if(cart_rank == 0){
		if(color_iter == 0)
			arr= malloc(superblock_edge * block_edge * width);
		if(no_channels == 1)
			inFile= fopen(argv[1], "r");
		else if(no_channels == 3){
			if(color_iter == 0)
				inFile= fopen(".red_temp", "r");
			else if(color_iter == 1)
				inFile= fopen(".green_temp", "r");
			else if(color_iter == 2)
				inFile= fopen(".blue_temp", "r");
		}
	}
	for(j= 0; j < height_in_sb; j++){
		switch(cart_rank){
		case 0:
			counter_1= 0;
			while(counter_1 < superblock_edge * block_edge * width){
				bytes_read= fread(arr + counter_1, 1, superblock_edge * block_edge * width - counter_1, inFile);
				if(bytes_read == -1){
					fprintf(stderr, "%s %d\n", strerror(errno), __LINE__);
					//MPI_Finalize();
					//return 1;
				}
				counter_1+= bytes_read;
			}
		default:
			for(k= 0; k < width_in_sb; k++){
				MPI_Scatterv( arr + (k * block_edge * superblock_edge), counts, displacements, block_type,
					local_arr + (k * block_edge * block_edge + j * (block_edge*block_edge) * width_in_sb),
					block_edge * block_edge, MPI_BYTE, 0, cartcomm);
			}
		}
	}
	if(cart_rank == 0)
		fclose(inFile);
	//-----------------------------------------------
	done= 0;	//initialize flag
	iteration_no= 1;

	//start timer 
	if(cart_rank == 0){
		gettimeofday(&start, NULL);
	}

	}//end of single thread execution

	//inner iteration
	while(done == 0){
	#pragma omp sections
	{
	//calculate, send and receive sides of neighbouring blocks
	//TOP_LEFT
	#pragma omp section
	{
	calc_side( TOP_LEFT, local_arr, one_arr, two_arr, four_arr, superblock_edge, block_edge, height_in_sb, width_in_sb);
	if( cart_rank == superblock_edge * superblock_edge - 1)
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 2) - 1, 1, diagonal_vector_lacking_both, neighb_ranks[BOTTOM_RIGHT], BOTTOM_RIGHT, cartcomm, &reqs[BOTTOM_RIGHT], one_arr, two_arr, four_arr);
	else if( cart_rank >= (superblock_edge + 2) * (superblock_edge + 1) )
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 2) - 1, 1, diagonal_vector_lacking_horizontal, neighb_ranks[BOTTOM_RIGHT], BOTTOM_RIGHT, cartcomm, &reqs[BOTTOM_RIGHT], one_arr, two_arr, four_arr);
	else if( (cart_rank + 1) % superblock_edge == 0 )
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 2) - 1, 1, diagonal_vector_lacking_vertical, neighb_ranks[BOTTOM_RIGHT], BOTTOM_RIGHT, cartcomm, &reqs[BOTTOM_RIGHT], one_arr, two_arr, four_arr);
	else
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 2) - 1, 1, diagonal_vector, neighb_ranks[BOTTOM_RIGHT], BOTTOM_RIGHT, cartcomm, &reqs[BOTTOM_RIGHT], one_arr, two_arr, four_arr);
	if( cart_rank == 0)
		tri_MPI_Isend( (block_edge + 2) * (block_edge + 2) * (width_in_sb + 1) + block_edge + 3, 1, diagonal_vector_lacking_both, neighb_ranks[TOP_LEFT], BOTTOM_RIGHT, cartcomm, &reqs[BOTTOM_RIGHT + 8], one_arr, two_arr, four_arr);
	else if( cart_rank < superblock_edge )
		tri_MPI_Isend( (block_edge + 2) * (block_edge + 2) * width_in_sb + block_edge + 3, 1, diagonal_vector_lacking_horizontal, neighb_ranks[TOP_LEFT], BOTTOM_RIGHT, cartcomm, &reqs[BOTTOM_RIGHT + 8], one_arr, two_arr, four_arr);
	else if( cart_rank % superblock_edge == 0)
		tri_MPI_Isend( (block_edge + 2) * (block_edge + 2) + block_edge + 3, 1, diagonal_vector_lacking_vertical, neighb_ranks[TOP_LEFT], BOTTOM_RIGHT, cartcomm, &reqs[BOTTOM_RIGHT + 8], one_arr, two_arr, four_arr);
	else
		tri_MPI_Isend( block_edge + 3, 1, diagonal_vector, neighb_ranks[TOP_LEFT], BOTTOM_RIGHT, cartcomm, &reqs[BOTTOM_RIGHT + 8], one_arr, two_arr, four_arr);
	}

	//TOP_RIGHT
	#pragma omp section
	{
	calc_side( TOP_RIGHT, local_arr, one_arr, two_arr, four_arr, superblock_edge, block_edge, height_in_sb, width_in_sb);
	if( cart_rank == superblock_edge * (superblock_edge - 1) )
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 2) + (block_edge + 2) * (block_edge + 1), 1, diagonal_vector_lacking_both, neighb_ranks[BOTTOM_LEFT], BOTTOM_LEFT, cartcomm, &reqs[BOTTOM_LEFT], one_arr, two_arr, four_arr);
	else if( cart_rank % superblock_edge == 0)
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 2) + (block_edge + 2) * (block_edge + 1), 1, diagonal_vector_lacking_vertical, neighb_ranks[BOTTOM_LEFT], BOTTOM_LEFT, cartcomm, &reqs[BOTTOM_LEFT], one_arr, two_arr, four_arr);
	else if( cart_rank >= superblock_edge *(superblock_edge - 1) )
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 1), 1, diagonal_vector_lacking_horizontal, neighb_ranks[BOTTOM_LEFT], BOTTOM_LEFT, cartcomm, &reqs[BOTTOM_LEFT], one_arr, two_arr, four_arr);
	else
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 1), 1, diagonal_vector, neighb_ranks[BOTTOM_LEFT], BOTTOM_LEFT, cartcomm, &reqs[BOTTOM_LEFT], one_arr, two_arr, four_arr);
	if( cart_rank == superblock_edge - 1)
		tri_MPI_Isend( 2 * (block_edge + 1) + (block_edge + 2) * (block_edge + 2) * width_in_sb, 1, diagonal_vector_lacking_both, neighb_ranks[TOP_RIGHT], BOTTOM_LEFT, cartcomm, &reqs[BOTTOM_LEFT + 8], one_arr, two_arr, four_arr);
	else if( cart_rank < superblock_edge )
		tri_MPI_Isend( 2 * (block_edge + 1) + (block_edge + 2) * (block_edge + 2) * width_in_sb, 1, diagonal_vector_lacking_horizontal, neighb_ranks[TOP_RIGHT], BOTTOM_LEFT, cartcomm, &reqs[BOTTOM_LEFT + 8], one_arr, two_arr, four_arr);
	else if( (cart_rank + 1) % superblock_edge == 0 )
		tri_MPI_Isend( 2 * (block_edge + 1), 1, diagonal_vector_lacking_vertical, neighb_ranks[TOP_RIGHT], BOTTOM_LEFT, cartcomm, &reqs[BOTTOM_LEFT + 8], one_arr, two_arr, four_arr);
	else
		tri_MPI_Isend( 2 * (block_edge + 1), 1, diagonal_vector, neighb_ranks[TOP_RIGHT], BOTTOM_LEFT, cartcomm, &reqs[BOTTOM_LEFT + 8], one_arr, two_arr, four_arr);
	}

	//BOTTOM_LEFT
	#pragma omp section
	{
	calc_side( BOTTOM_LEFT, local_arr, one_arr, two_arr, four_arr, superblock_edge, block_edge, height_in_sb, width_in_sb);
	if( cart_rank == superblock_edge - 1)
		tri_MPI_Irecv( block_edge + 1 + (block_edge + 2) * (block_edge + 2) * width_in_sb, 1, diagonal_vector_lacking_both, neighb_ranks[TOP_RIGHT], TOP_RIGHT, cartcomm, &reqs[TOP_RIGHT], one_arr, two_arr, four_arr);
	else if( cart_rank < superblock_edge )
		tri_MPI_Irecv( block_edge + 1 + (block_edge + 2) * (block_edge + 2) * width_in_sb, 1, diagonal_vector_lacking_horizontal, neighb_ranks[TOP_RIGHT], TOP_RIGHT, cartcomm, &reqs[TOP_RIGHT], one_arr, two_arr, four_arr);
	else if( (cart_rank + 1) % superblock_edge == 0 )
		tri_MPI_Irecv( block_edge + 1, 1, diagonal_vector_lacking_vertical, neighb_ranks[TOP_RIGHT], TOP_RIGHT, cartcomm, &reqs[TOP_RIGHT], one_arr, two_arr, four_arr);
	else
		tri_MPI_Irecv( block_edge + 1, 1, diagonal_vector, neighb_ranks[TOP_RIGHT], TOP_RIGHT, cartcomm, &reqs[TOP_RIGHT], one_arr, two_arr, four_arr);
	if( cart_rank == superblock_edge * (superblock_edge - 1) )
		tri_MPI_Isend( (block_edge + 2) * block_edge + 1 + (block_edge + 2) * (block_edge + 2), 1, diagonal_vector_lacking_both, neighb_ranks[BOTTOM_LEFT], TOP_RIGHT, cartcomm, &reqs[TOP_RIGHT + 8], one_arr, two_arr, four_arr);
	else if( cart_rank >= superblock_edge * (superblock_edge - 1) )
		tri_MPI_Isend( (block_edge + 2) * block_edge + 1, 1, diagonal_vector_lacking_horizontal, neighb_ranks[BOTTOM_LEFT], TOP_RIGHT, cartcomm, &reqs[TOP_RIGHT + 8], one_arr, two_arr, four_arr);
	else if( cart_rank % superblock_edge == 0 )
		tri_MPI_Isend( (block_edge + 2) * block_edge + 1 + (block_edge + 2) * (block_edge + 2), 1, diagonal_vector_lacking_vertical, neighb_ranks[BOTTOM_LEFT], TOP_RIGHT, cartcomm, &reqs[TOP_RIGHT + 8], one_arr, two_arr, four_arr);
	else
		tri_MPI_Isend( (block_edge + 2) * block_edge + 1, 1, diagonal_vector, neighb_ranks[BOTTOM_LEFT], TOP_RIGHT, cartcomm, &reqs[TOP_RIGHT + 8], one_arr, two_arr, four_arr);
	}

	//BOTTOM_RIGHT
	#pragma omp section
	{
	calc_side( BOTTOM_RIGHT, local_arr, one_arr, two_arr, four_arr, superblock_edge, block_edge, height_in_sb, width_in_sb);
	if( cart_rank == 0 )
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 2) * (width_in_sb + 1), 1, diagonal_vector_lacking_both, neighb_ranks[TOP_LEFT], TOP_LEFT, cartcomm, &reqs[TOP_LEFT], one_arr, two_arr, four_arr);
	else if( cart_rank < superblock_edge)
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 2) * width_in_sb, 1, diagonal_vector_lacking_horizontal, neighb_ranks[TOP_LEFT], TOP_LEFT, cartcomm, &reqs[TOP_LEFT], one_arr, two_arr, four_arr);
	else if( cart_rank % superblock_edge == 0 )
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 2), 1, diagonal_vector_lacking_vertical, neighb_ranks[TOP_LEFT], TOP_LEFT, cartcomm, &reqs[TOP_LEFT], one_arr, two_arr, four_arr);
	else
		tri_MPI_Irecv( 0, 1, diagonal_vector, neighb_ranks[TOP_LEFT], TOP_LEFT, cartcomm, &reqs[TOP_LEFT], one_arr, two_arr, four_arr);
	if( cart_rank == superblock_edge * superblock_edge - 1 )
		tri_MPI_Isend( (block_edge + 2) * (block_edge + 1) - 2, 1, diagonal_vector_lacking_both, neighb_ranks[BOTTOM_RIGHT], TOP_LEFT, cartcomm, &reqs[TOP_LEFT + 8], one_arr, two_arr, four_arr);
	else if( (cart_rank + 1) % superblock_edge == 0 )
		tri_MPI_Isend( (block_edge + 2) * (block_edge + 1) - 2, 1, diagonal_vector_lacking_vertical, neighb_ranks[BOTTOM_RIGHT], TOP_LEFT, cartcomm, &reqs[TOP_LEFT + 8], one_arr, two_arr, four_arr);
	else if( cart_rank >= superblock_edge * (superblock_edge - 1) )
		tri_MPI_Isend( (block_edge + 2) * (block_edge + 1) - 2, 1, diagonal_vector_lacking_horizontal, neighb_ranks[BOTTOM_RIGHT], TOP_LEFT, cartcomm, &reqs[TOP_LEFT + 8], one_arr, two_arr, four_arr);
	else
		tri_MPI_Isend( (block_edge + 2) * (block_edge + 1) - 2, 1, diagonal_vector, neighb_ranks[BOTTOM_RIGHT], TOP_LEFT, cartcomm, &reqs[TOP_LEFT + 8], one_arr, two_arr, four_arr);
	}

	}//end of "corners" sections

	#pragma omp sections
	{

	//TOP
	#pragma omp section
	{
	calc_side( TOP, local_arr, one_arr, two_arr, four_arr, superblock_edge, block_edge, height_in_sb, width_in_sb);
	if( cart_rank < neighb_ranks[BOTTOM] )
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 1) + 1, 1, horizontal_vector, neighb_ranks[BOTTOM], BOTTOM, cartcomm, &reqs[BOTTOM], one_arr, two_arr, four_arr);
	else
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 1) + 1, 1, horizontal_vector_2, neighb_ranks[BOTTOM], BOTTOM, cartcomm, &reqs[BOTTOM], one_arr, two_arr, four_arr);
	if( cart_rank > neighb_ranks[TOP] )
		tri_MPI_Isend( block_edge + 3, 1, horizontal_vector, neighb_ranks[TOP], BOTTOM, cartcomm, &reqs[BOTTOM + 8], one_arr, two_arr, four_arr);
	else
		tri_MPI_Isend( block_edge + 3 + width_in_sb * (block_edge + 2) * (block_edge + 2), 1, horizontal_vector_2, neighb_ranks[TOP], BOTTOM, cartcomm, &reqs[BOTTOM + 8], one_arr, two_arr, four_arr);
	}

	//BOTTOM
	#pragma omp section
	{
	calc_side( BOTTOM, local_arr, one_arr, two_arr, four_arr, superblock_edge, block_edge, height_in_sb, width_in_sb);
	if( cart_rank > neighb_ranks[TOP] )
		tri_MPI_Irecv( 1, 1, horizontal_vector, neighb_ranks[TOP], TOP, cartcomm, &reqs[TOP], one_arr, two_arr, four_arr);
	else
		tri_MPI_Irecv( 1 + (block_edge + 2) * (block_edge + 2) * width_in_sb, 1, horizontal_vector_2, neighb_ranks[TOP], TOP, cartcomm, &reqs[TOP], one_arr, two_arr, four_arr);
	if( cart_rank < neighb_ranks[BOTTOM] )
		tri_MPI_Isend( (block_edge + 2) * block_edge + 1, 1, horizontal_vector, neighb_ranks[BOTTOM], TOP, cartcomm, &reqs[TOP + 8], one_arr, two_arr, four_arr);
	else
		tri_MPI_Isend( (block_edge + 2) * block_edge + 1, 1, horizontal_vector_2, neighb_ranks[BOTTOM], TOP, cartcomm, &reqs[TOP + 8], one_arr, two_arr, four_arr);
	}

	//LEFT
	#pragma omp section
	{
	calc_side( LEFT, local_arr, one_arr, two_arr, four_arr, superblock_edge, block_edge, height_in_sb, width_in_sb);
	if( cart_rank < neighb_ranks[RIGHT] )
		tri_MPI_Irecv( (block_edge + 2) * 2 - 1, 1, vertical_vector, neighb_ranks[RIGHT], RIGHT, cartcomm, &reqs[RIGHT], one_arr, two_arr, four_arr);
	else
		tri_MPI_Irecv( (block_edge + 2) * 2 - 1, 1, vertical_vector_2, neighb_ranks[RIGHT], RIGHT, cartcomm, &reqs[RIGHT], one_arr, two_arr, four_arr);
	if( cart_rank > neighb_ranks[LEFT] )
		tri_MPI_Isend( block_edge + 3, 1, vertical_vector, neighb_ranks[LEFT], RIGHT, cartcomm, &reqs[RIGHT + 8], one_arr, two_arr, four_arr);
	else
		tri_MPI_Isend( block_edge + 3 + (block_edge + 2) * (block_edge + 2), 1, vertical_vector_2, neighb_ranks[LEFT], RIGHT, cartcomm, &reqs[RIGHT + 8], one_arr, two_arr, four_arr);	
	}

	//RIGHT
	#pragma omp section
	{
	calc_side( RIGHT, local_arr, one_arr, two_arr, four_arr, superblock_edge, block_edge, height_in_sb, width_in_sb);
	if( cart_rank > neighb_ranks[LEFT] )
		tri_MPI_Irecv( block_edge + 2, 1, vertical_vector, neighb_ranks[LEFT], LEFT, cartcomm, &reqs[LEFT], one_arr, two_arr, four_arr);
	else
		tri_MPI_Irecv( (block_edge + 2) * (block_edge + 3), 1, vertical_vector_2, neighb_ranks[LEFT], LEFT, cartcomm, &reqs[LEFT], one_arr, two_arr, four_arr);
	if( cart_rank < neighb_ranks[RIGHT] )
		tri_MPI_Isend( (block_edge + 2) * 2 - 2, 1, vertical_vector, neighb_ranks[RIGHT], LEFT, cartcomm, &reqs[LEFT + 8], one_arr, two_arr, four_arr);
	else
		tri_MPI_Isend( (block_edge + 2) * 2 - 2, 1, vertical_vector_2, neighb_ranks[RIGHT], LEFT, cartcomm, &reqs[LEFT + 8], one_arr, two_arr, four_arr);
	}

	}//end of "sides" sections

	//calculate blocks, except for sides since they have already been computed
	#pragma omp for \
		private(j, k, indx_1, indx_2)
	for(i= 0; i < height_in_sb * width_in_sb; i++){
		for(j= 1; j < block_edge - 1; j++){
			for(k= 1; k < block_edge - 1; k++){
				indx_1= block_edge + 3 + k + j * (block_edge + 2) + i * (block_edge + 2) * (block_edge + 2);
				indx_2= k + j * block_edge + i * block_edge * block_edge;
				one_arr[indx_1]= (unsigned short) local_arr[indx_2];
				two_arr[indx_1]= ((unsigned short) local_arr[indx_2])  * 2;
				four_arr[indx_1]= ((unsigned short) local_arr[indx_2]) * 4;
			}
		}
	}
	//------------------------------------------------------------------------

	#pragma omp single
	{
	MPI_Waitall(16, reqs, stats);
	}

	//compute outer frame of the intermediate arrays based on copies of neighboring elements of the edges of the initial array
	if( cart_rank < superblock_edge ){
		#pragma omp for \
			private(j, indx_1, indx_2)
		for(i= 0; i < width_in_sb; i++){
			for(j= 0; j < block_edge + 2; j++){
				indx_1= j + i * (block_edge + 2) * (block_edge + 2);
				indx_2= indx_1 + block_edge + 2;
				one_arr[indx_1]= one_arr[indx_2];
				two_arr[indx_1]= two_arr[indx_2];
			}
		}
	}
	if( cart_rank >= superblock_edge * (superblock_edge - 1) ){
		#pragma omp for \
			private(j, indx_1, indx_2)
		for(i= 0; i < width_in_sb; i++){
			for(j= 0; j < block_edge + 2; j++){
				indx_1= ((height_in_sb - 1) * width_in_sb * (block_edge + 2) * (block_edge + 2)) + (block_edge + 2) * (block_edge + 1) + j + i * (block_edge + 2) * (block_edge + 2);
				indx_2= indx_1 - (block_edge + 2);
				one_arr[indx_1]= one_arr[indx_2];
				two_arr[indx_1]= two_arr[indx_2];
			}
		}
	}
	if( cart_rank % superblock_edge == 0){
		#pragma omp for \
			private(j, indx_1, indx_2)
		for(i= 0; i < height_in_sb; i++){
			for(j= 0; j < block_edge + 2; j++){
				indx_1= j * (block_edge + 2) + i * width_in_sb * (block_edge + 2) * (block_edge + 2);
				indx_2= indx_1 + 1;
				one_arr[indx_1]= one_arr[indx_2];
				two_arr[indx_1]= two_arr[indx_2];
			}
		}
	}
	if( (cart_rank + 1) % superblock_edge == 0){
		#pragma omp for \
			private(j, indx_1, indx_2)
		for(i= 0; i < height_in_sb; i++){
			for(j= 0; j < block_edge + 2; j++){
				indx_1= (width_in_sb - 1) * (block_edge + 2) * (block_edge + 2) + block_edge + 1 + j * (block_edge + 2) + i * width_in_sb * (block_edge + 2) * (block_edge + 2);
				indx_2= indx_1 - 1;
				one_arr[indx_1]= one_arr[indx_2];
				two_arr[indx_1]= two_arr[indx_2];
			}
		}
	}
	//------------------------------------------------------------------------------------------------------------------------

	//calculate final array of current iteration
	#pragma omp for \
		private(j, k, indx_1, indx_2)
	for(i= 0; i < height_in_sb * width_in_sb; i++){
		for(j= 1; j < block_edge + 1; j++){
			for(k= 1; k < block_edge + 1; k++){
				indx_1= k - 1 + (j - 1) * block_edge + i * block_edge * block_edge;
				indx_2= k + j * (block_edge + 2) + i * (block_edge + 2) * (block_edge + 2);
				local_arr_new[indx_1]
					= (four_arr[indx_2] 
						+ two_arr[indx_2 - 1] + two_arr[indx_2 + 1] + two_arr[indx_2 - (block_edge + 2)] + two_arr[indx_2 + (block_edge + 2)]
						+ one_arr[indx_2 - 1 - (block_edge + 2)] + one_arr[indx_2 + 1 - (block_edge + 2)] 
						+ one_arr[indx_2 - 1 + (block_edge + 2)] + one_arr[indx_2 + 1 + (block_edge + 2)]) / 16;
			}
		}
	}
	//-------------------------------------------

	#pragma omp for \
		reduction(+:sum_of_differences)
	for(i= 0; i < block_edge * block_edge * width_in_sb * height_in_sb; i++){
		sum_of_differences+= ((double) abs((int) (local_arr[i] - local_arr_new[i]))) / 255;
	}
	#pragma omp single
	{//single thread execution starts
	sum_of_differences/= width * height;
	MPI_Reduce(&sum_of_differences, &difference_q, 1, MPI_DOUBLE, MPI_SUM, 0, cartcomm );
	if(cart_rank == 0){
		printf("Iteration no %d, difference between this and the previous array: %f%%\n", iteration_no, difference_q);
		if(difference_q < lower_bound_convergence || iteration_no >= atoi(argv[6])) 	//stop if difference between the last 2 iterations is negligible or if the max nm of iters
			done= 1;																	//has been reached
	}
	MPI_Bcast(&done, 1, MPI_CHAR, 0, cartcomm);
	if( done == 0){	//incase of remaining iterations, swap arrays
		temp_ptr= local_arr;
		local_arr= local_arr_new;
		local_arr_new= temp_ptr;
	}
	iteration_no++;
	}//single thread execution ends

	}//inner iteration

	#pragma omp single
	{//single thread execution starts

	//stop timer for current iteration
	if(cart_rank == 0){
		gettimeofday(&stop, NULL);
		time_sum+= (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
	}
	//----------


	//gather and print the array
	if(cart_rank == 0){
		if(no_channels == 1){
			outFile= fopen(argv[2], "w");
			printf("Printing output file.\n");
		}
		else if(no_channels == 3){
			if(color_iter == 0){
				outFile= fopen(".red_temp", "w");
				printf("Printing red_temp.\n\n");
			}
			else if(color_iter == 1){
				outFile= fopen(".green_temp", "w");
				printf("Printing green_temp.\n\n");
			}
			else if(color_iter == 2){
				outFile= fopen(".blue_temp", "w");
				printf("Printing blue_temp.\n\n");
			}
		}
	}
	for(j= 0; j < height_in_sb; j++){
		for(k= 0; k < width_in_sb; k++){
			MPI_Gatherv( local_arr_new + (k * block_edge * block_edge + j * (block_edge*block_edge) * width_in_sb),
				block_edge *block_edge, MPI_BYTE, arr + (k * block_edge *superblock_edge), counts, displacements,
				block_type, 0, cartcomm);
		}
		if(cart_rank == 0){
			counter_1= 0;
			while(counter_1 < superblock_edge * block_edge * width){
				bytes_written= fwrite(arr + counter_1, 1, superblock_edge * block_edge * width - counter_1, outFile);
				if(bytes_written == -1){
					fprintf(stderr, "%s\n %d", strerror(errno), __LINE__);
					//MPI_Finalize();
					//return 1;
				}
				counter_1+= bytes_written;
			}
		}
	}
	if(cart_rank == 0)
		fclose(outFile);
	//--------------------------

	}//single thread execution ends
	
	}//outer iteration

	}//omp parallelization ends

	//print total time
	if(cart_rank == 0){
		printf("Total time excluding time it took to read and print: %lu microsecs\n", time_sum);
	}
	//----------------

	//in case of RGB input file, combine the 3 output channel files to construct the final output file
	if(no_channels == 3 && cart_rank == 0){
		printf("Combining temps into final output.\n");
		outFile= fopen(argv[2], "w");
		color_file[0]= fopen(".red_temp", "r");
		color_file[1]= fopen(".green_temp", "r");
		color_file[2]= fopen(".blue_temp", "r");
		for(i= 0; i < width * height; i++){
			fread(arr, 1, 1, color_file[0]);
			fwrite(arr, 1, 1, outFile);
			fread(arr, 1, 1, color_file[1]);
			fwrite(arr, 1, 1, outFile);
			fread(arr, 1, 1, color_file[2]);
			fwrite(arr, 1, 1, outFile);
		}
		for(i= 0; i < no_channels; i++){
			fclose(color_file[i]);
		}
		fclose(outFile);
		unlink(".red_temp");
		unlink(".green_temp");
		unlink(".blue_temp");
	}
	MPI_Type_free(&block_type);
	MPI_Type_free(&horizontal_vector);
	MPI_Type_free(&horizontal_vector_2);
	MPI_Type_free(&vertical_vector);
	MPI_Type_free(&vertical_vector_2);
	MPI_Type_free(&diagonal_vector);
	MPI_Type_free(&diagonal_vector_lacking_vertical);
	MPI_Type_free(&diagonal_vector_lacking_horizontal);
	MPI_Type_free(&diagonal_vector_lacking_both);
	free(local_arr);
	free(local_arr_new);
	free(one_arr);
	free(two_arr);
	free(four_arr);
	if(cart_rank == 0)
		free(arr);
	//------------------------------------------------------------------------------------------------

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
