typedef enum Sides{
	TOP_LEFT= 0,
	TOP,
	TOP_RIGHT,
	RIGHT,
	BOTTOM_RIGHT,
	BOTTOM,
	BOTTOM_LEFT,
	LEFT
} Sides;



int calc_side(
	int side,
	unsigned char * local_arr,
	unsigned short * one_arr,
	unsigned short * two_arr,
	unsigned short * four_arr,
	int superblock_edge,
	int block_edge,
	int height_in_sb,
	int width_in_sb);
