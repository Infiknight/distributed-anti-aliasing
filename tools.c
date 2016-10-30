#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include "tools.h"

int calc_horizontal_side(
	int side,
	unsigned char * local_arr,
	unsigned short * one_arr,
	unsigned short * two_arr,
	unsigned short * four_arr,
	int superblock_edge,
	int block_edge,
	int height_in_sb,
	int width_in_sb)
{
	int i, j, k;
	unsigned char current_byte;
	int starting_pos= block_edge + 3;
	if(side == BOTTOM){
		local_arr+= block_edge * (block_edge - 1);
		starting_pos= (block_edge + 2) * block_edge + 1;
	}
	else if(side != TOP)
		return 1;
	for(i= 0; i < height_in_sb * width_in_sb; i++){
		for(j= 1; j < block_edge - 1; j++){
			current_byte= *(local_arr + j + i * block_edge * block_edge);
			k= starting_pos + j + i * (block_edge + 2) * (block_edge + 2);
			one_arr[k]= (unsigned short) current_byte;
			two_arr[k]= ((unsigned short) current_byte) * 2;
			four_arr[k]= ((unsigned short) current_byte) * 4;
		}
	}
	return 0;
}

int calc_vertical_side(
	int side,
	unsigned char * local_arr,
	unsigned short * one_arr,
	unsigned short * two_arr,
	unsigned short * four_arr,
	int superblock_edge,
	int block_edge,
	int height_in_sb,
	int width_in_sb)
{
	int i, j, k;
	unsigned char current_byte;
	int starting_pos= block_edge + 3;
	if(side == RIGHT){
		local_arr+= block_edge - 1;
		starting_pos= 2*(block_edge + 2) - 2;
	}
	else if(side != LEFT)
		return 1;
	for(i= 0; i < height_in_sb * width_in_sb; i++){
		for(j= 1; j < block_edge - 1; j++){
			current_byte= *(local_arr + j * block_edge + i * block_edge * block_edge);
			k= starting_pos + j * (block_edge + 2) + i * (block_edge + 2) * (block_edge + 2);
			one_arr[k]= (unsigned short) current_byte;
			two_arr[k]= ((unsigned short) current_byte) * 2;
			four_arr[k]= ((unsigned short) current_byte) * 4;
		}
	}
	return 0;
}

int calc_corner(
	int side,
	unsigned char * local_arr,
	unsigned short * one_arr,
	unsigned short * two_arr,
	unsigned short * four_arr,
	int superblock_edge,
	int block_edge,
	int height_in_sb,
	int width_in_sb)
{
	int i, k;
	unsigned char current_byte;
	int starting_pos= block_edge + 3;
	if(side == TOP_RIGHT){
		local_arr+= block_edge - 1;
		starting_pos= 2*(block_edge + 2) - 2;
	}
	else if(side == BOTTOM_RIGHT){
		local_arr+= block_edge * block_edge - 1;
		starting_pos= (block_edge + 2) * (block_edge + 1) - 2;
	}
	else if(side == BOTTOM_LEFT){
		local_arr+= (block_edge - 1) * block_edge;
		starting_pos= (block_edge + 2) * block_edge + 1;
	}
	else if(side != TOP_LEFT)
		return 1;
	for(i= 0; i < height_in_sb * width_in_sb; i++){
		current_byte= *(local_arr + i * block_edge * block_edge);
		k= starting_pos + i * (block_edge + 2) * (block_edge + 2);
		one_arr[k]= (unsigned short) current_byte;
		two_arr[k]= ((unsigned short) current_byte) * 2;
		four_arr[k]= ((unsigned short) current_byte) * 4;
	}
	return 0;
}

int calc_side(
	int side,
	unsigned char * local_arr,
	unsigned short * one_arr,
	unsigned short * two_arr,
	unsigned short * four_arr,
	int superblock_edge,
	int block_edge,
	int height_in_sb,
	int width_in_sb)
{
	if(side == TOP || side == BOTTOM)
		return calc_horizontal_side(
			side,
			local_arr,
			one_arr,
			two_arr,
			four_arr,
			superblock_edge,
			block_edge,
			height_in_sb,
			width_in_sb);
	else if(side == RIGHT || side == LEFT)
		return calc_vertical_side(
			side,
			local_arr,
			one_arr,
			two_arr,
			four_arr,
			superblock_edge,
			block_edge,
			height_in_sb,
			width_in_sb);
	else
		return calc_corner(
			side,
			local_arr,
			one_arr,
			two_arr,
			four_arr,
			superblock_edge,
			block_edge,
			height_in_sb,
			width_in_sb);
}
