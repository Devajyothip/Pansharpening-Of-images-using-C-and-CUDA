//by Devajyothi Potnuru
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int g_img_width;
int g_img_height;
short* g_img_data;

short* g_imgout1_data;
short* g_imgout2_data;
short* g_imgout3_data;

int img1_width = 0, img1_height = 0;
int img2_width = 0, img2_height = 0;
int img3_width = 0, img3_height = 0;

void read_TIFF(FILE *fp) {
	unsigned long fileLen = 0;
	unsigned long off_set = 0, next_offset = 0, pos_count = 0, tag_value[25] = { 0 }, strip_offset_val = 0, strip_offset = 0;
	int i, j, k, tag_id[25] = { 0 }, tag_type[25] = { 0 }, tag_count[25] = { 0 };
	int tagCount = 0, img_Width = 0, img_Len = 0, IFD_count = 0;
	short buffer[4] = { 0 };
	fseek(fp, 4, SEEK_SET);
	for (i = 3; i >= 0; i--) {
		buffer[i] = getc(fp);
	}
	off_set = ((int)buffer[0] << 24) | ((int)buffer[1] << 16) | ((int)buffer[2] << 8) | ((int)buffer[3]);

	fseek(fp, off_set, SEEK_SET);
	k = 1;
	while (k) {
		IFD_count++;                                            // No. of IFD's in this TIFF File
		buffer[1] = getc(fp);
		buffer[0] = getc(fp);
		tagCount = ((int)buffer[0] << 8) | ((int)buffer[1]);    // No. of entries in an IFD
		for (i = 0; i < tagCount; i++) {                           // Read all the entries of this IFD
			buffer[1] = getc(fp);
			buffer[0] = getc(fp);
			tag_id[i] = ((int)buffer[0] << 8) | ((int)buffer[1]);       // Tag ID

			for (j = 1; j >= 0; j--) {
				buffer[j] = getc(fp);
			}
			tag_type[i] = ((int)buffer[0] << 8) | ((int)buffer[1]);     // Type

			for (j = 3; j >= 0; j--) {
				buffer[j] = getc(fp);
			}
			tag_count[i] = ((int)buffer[0] << 24) | ((int)buffer[1] << 16) | ((int)buffer[2] << 8) | ((int)buffer[3]);
              // Gives no. of values for this Tag
			for (j = 3; j >= 0; j--) {
				buffer[j] = getc(fp);
			}
			tag_value[i] = ((int)buffer[0] << 24) | ((int)buffer[1] << 16) | ((int)buffer[2] << 8) | ((int)buffer[3]);
              // Gets the value if the above count is 1, else offset of the starting value
			if (tag_id[i] == 256)                                     // Tag ID 256 says about the image width
				img_Width = tag_value[i];
			if (tag_id[i] == 257)                                     // Tag ID 257 says about the image length
				img_Len = tag_value[i];
			if (tag_id[i] == 273)                                     // Tag ID 273 says about the offset which points to the offset of strips
				strip_offset_val = tag_value[i];
		}
		g_img_data = (short *)malloc(img_Len * img_Width * sizeof(short));    // Create a matrix of image size
	
		g_img_width = img_Width;
		g_img_height = img_Len;
		for (i = 0; i < img_Len; i++) {                                // Read the pixel values from image and store it in the matrix
			fseek(fp, (strip_offset_val + (i * 4)), SEEK_SET);              // Move to the offset of the current strip's offset

			for (j = 3; j >= 0; j--) {
				buffer[j] = getc(fp);
			}
			strip_offset = ((int)buffer[0] << 24) | ((int)buffer[1] << 16) | ((int)buffer[2] << 8) | ((int)buffer[3]);
			fseek(fp, strip_offset, SEEK_SET);    // Move to the offset of the current strip
			for (j = 0; j < img_Width; j++) {

				short tmp_c1 = getc(fp);
				short tmp_c2 = getc(fp);
				g_img_data[i*img_Width + j] = tmp_c2;                    // Read the current strip values
				
			}
		}
		pos_count = ((off_set + 2) + (tagCount * 12));
		fseek(fp, pos_count, SEEK_SET);
		for (i = 3; i >= 0; i--) {
			buffer[i] = getc(fp);
		}
		next_offset = ((int)buffer[0] << 24) | ((int)buffer[1] << 16) | ((int)buffer[2] << 8) | ((int)buffer[3]); // Next IFD offset
		if (next_offset != 0)
		{
			fseek(fp, next_offset, SEEK_SET);
		}
		else {
			k = 0;
		}
	}
}

void WriteHexString(FILE *fptr, char *s)
{
	unsigned int i, c;
	char hex[3];

	for (i = 0; i<strlen(s); i += 2) {
		hex[0] = s[i];
		hex[1] = s[i + 1];
		hex[2] = '\0';
		sscanf(hex, "%X", &c);
		putc(c, fptr);
	}
}

void write_TIFF(FILE *fptr){
	/* Write the header */
	WriteHexString(fptr, "4d4d002a");    /* Big endian & TIFF identifier */
	int nx = g_img_width;
	int ny = g_img_height;
	int offset = nx * ny * 3 + 8;
	putc((offset & 0xff000000) / 16777216, fptr);
	putc((offset & 0x00ff0000) / 65536, fptr);
	putc((offset & 0x0000ff00) / 256, fptr);
	putc((offset & 0x000000ff), fptr);

	/* Write the binary data */
	int i, j,p,q;
	short c_r, c_g, c_b, c_p;
	short c_i, c_h, c_s;
        short t1 =0,t2=0,t3=0;
	for(p=0;p<ny;p++)
	{
		for(q=0;q<nx;q++)
		{
		t1 = t1+g_imgout1_data[(int)(q / 2)*img1_width + (int)(p / 2)];
		t2 = t2+g_imgout2_data[(int)(q / 2)*img1_width + (int)(p / 2)];
		t3 = t3+g_imgout2_data[(int)(q / 2)*img1_width + (int)(p / 2)];
		}
	}
   
	for (j = 0; j<ny; j++) {
		for (i = 0; i<nx; i++) {
			// Getting RGB, P
			c_r = g_imgout1_data[(int)(j / 2)*img1_width + (int)(i / 2)];
			c_g = g_imgout2_data[(int)(j / 2)*img1_width + (int)(i / 2)];
			c_b = g_imgout3_data[(int)(j / 2)*img1_width + (int)(i / 2)];
			c_p = g_img_data[j*g_img_width + i];
                        
                        c_i =(short)((c_r/t1)*(c_p));
                        c_h =(short)((c_g/t1)*(c_p));
                        c_s =(short)((c_b/t1)*(c_p));
                       
			fputc(c_r, fptr);  //R
			fputc(c_g, fptr);  //G
			fputc(c_b, fptr);  //B
                                

		}
	}

	/* Write the footer */
	WriteHexString(fptr, "000e");  /* The number of directory entries (14) */

	/* Width tag, short int */
	WriteHexString(fptr, "0100000300000001");
	fputc((nx & 0xff00) / 256, fptr);    /* Image width */
	fputc((nx & 0x00ff), fptr);
	WriteHexString(fptr, "0000");

	/* Height tag, short int */
	WriteHexString(fptr, "0101000300000001");
	fputc((ny & 0xff00) / 256, fptr);    /* Image height */
	fputc((ny & 0x00ff), fptr);
	WriteHexString(fptr, "0000");

	/* Bits per sample tag, short int */
	WriteHexString(fptr, "0102000300000003");
	offset = nx * ny * 3 + 182;
	putc((offset & 0xff000000) / 16777216, fptr);
	putc((offset & 0x00ff0000) / 65536, fptr);
	putc((offset & 0x0000ff00) / 256, fptr);
	putc((offset & 0x000000ff), fptr);

	/* Compression flag, short int */
	WriteHexString(fptr, "010300030000000100010000");

	/* Photometric interpolation tag, short int */
	WriteHexString(fptr, "010600030000000100020000");

	/* Strip offset tag, long int */
	WriteHexString(fptr, "011100040000000100000008");

	/* Orientation flag, short int */
	WriteHexString(fptr, "011200030000000100010000");

	/* Sample per pixel tag, short int */
	WriteHexString(fptr, "011500030000000100030000");

	/* Rows per strip tag, short int */
	WriteHexString(fptr, "0116000300000001");
	fputc((ny & 0xff00) / 256, fptr);
	fputc((ny & 0x00ff), fptr);
	WriteHexString(fptr, "0000");

	/* Strip byte count flag, long int */
	WriteHexString(fptr, "0117000400000001");
	offset = nx * ny * 3;
	putc((offset & 0xff000000) / 16777216, fptr);
	putc((offset & 0x00ff0000) / 65536, fptr);
	putc((offset & 0x0000ff00) / 256, fptr);
	putc((offset & 0x000000ff), fptr);

	/* Minimum sample value flag, short int */
	WriteHexString(fptr, "0118000300000003");
	offset = nx * ny * 3 + 188;
	putc((offset & 0xff000000) / 16777216, fptr);
	putc((offset & 0x00ff0000) / 65536, fptr);
	putc((offset & 0x0000ff00) / 256, fptr);
	putc((offset & 0x000000ff), fptr);

	/* Maximum sample value tag, short int */
	WriteHexString(fptr, "0119000300000003");
	offset = nx * ny * 3 + 194;
	putc((offset & 0xff000000) / 16777216, fptr);
	putc((offset & 0x00ff0000) / 65536, fptr);
	putc((offset & 0x0000ff00) / 256, fptr);
	putc((offset & 0x000000ff), fptr);

	/* Planar configuration tag, short int */
	WriteHexString(fptr, "011c00030000000100010000");

	/* Sample format tag, short int */
	WriteHexString(fptr, "0153000300000003");
	offset = nx * ny * 3 + 200;
	putc((offset & 0xff000000) / 16777216, fptr);
	putc((offset & 0x00ff0000) / 65536, fptr);
	putc((offset & 0x0000ff00) / 256, fptr);
	putc((offset & 0x000000ff), fptr);

	/* End of the directory entry */
	WriteHexString(fptr, "00000000");

	/* Bits for each colour channel */
	WriteHexString(fptr, "000800080008");

	/* Minimum value for each component */
	WriteHexString(fptr, "000000000000");

	/* Maximum value per channel */
	WriteHexString(fptr, "00ff00ff00ff");

	/* Samples per pixel for each channel */
	WriteHexString(fptr, "000100010001");
}
int main(int argc, char** argv) {

	if (argc != 5){
		printf("ihs <R-filename1> <G-filename2> <B-filename3> <Pan-filename4>\n");
		return 1;
	}
	int rc = 0;
	FILE *fptr_r = fopen(argv[1], "rb");
	FILE *fptr_g = fopen(argv[2], "rb");
	FILE *fptr_b = fopen(argv[3], "rb");
	FILE *fptr_p = fopen(argv[4], "rb");

	// File Existing
	if (fptr_r == NULL || fptr_g == NULL || fptr_b == NULL || fptr_p == NULL) {
		printf("Failed to open file\n");
		return EXIT_FAILURE;
	}
	
	// Read R-File
	fseek(fptr_r, 0, SEEK_SET);
	rc = getc(fptr_r);
	if (rc == 73 || rc == 77) {
		printf("Reading R-Gray Image...\n");
		read_TIFF(fptr_r);

		img1_width = g_img_width;
		img1_height = g_img_height;
		g_imgout1_data = (short *)malloc(img1_width * img1_height * sizeof(short));    // Create a matrix of image size
		memccpy(g_imgout1_data, g_img_data, sizeof(short), img1_width * img1_height*sizeof(short));
	}
	else {
		printf("Invalid file.\n");
		return EXIT_FAILURE;
	}

	// Read G-File
	fseek(fptr_g, 0, SEEK_SET);
	rc = getc(fptr_g);
	if (rc == 73 || rc == 77) {
		printf("Reading G-Gray Image...\n");
		read_TIFF(fptr_g);

		img2_width = g_img_width;
		img2_height = g_img_height;
		g_imgout2_data = (short *)malloc(img2_width * img2_height * sizeof(short));    // Create a matrix of image size
		memccpy(g_imgout2_data, g_img_data, sizeof(short), img2_width * img2_height*sizeof(short));
	}
	else {
		printf("Invalid file.\n");
		return EXIT_FAILURE;
	}
	
	// Read B-File
	fseek(fptr_b, 0, SEEK_SET);
	rc = getc(fptr_b);
	if (rc == 73 || rc == 77) {
		printf("Reading B-Gray Image...\n");
		read_TIFF(fptr_b);

		img3_width = g_img_width;
		img3_height = g_img_height;
		g_imgout3_data = (short *)malloc(img3_width * img3_height * sizeof(short));    // Create a matrix of image size
		memccpy(g_imgout3_data, g_img_data, sizeof(short), img3_width * img3_height*sizeof(short));
	}
	else {
		printf("Invalid file.\n");
		return EXIT_FAILURE;
	}
	
	// Read P-File
	
	fseek(fptr_p, 0, SEEK_SET);
	rc = getc(fptr_p);
	if (rc == 73 || rc == 77) {
		printf("Reading P-Gray Image...\n");
		read_TIFF(fptr_p);
	}
	else {
		printf("Invalid file.\n");
		return EXIT_FAILURE;
	}
	
	if (img1_width != img2_width ||
		img2_width != img3_width ||
		img1_height != img2_height ||
		img2_height != img3_height
		){
		printf("Not Matched Image Size.\n");
		return EXIT_FAILURE;
	}

	printf("Reconstructing Color-Image...\n");
	
	FILE *fptr_out = fopen("ihs_browvey_output.tif", "wb");
	write_TIFF(fptr_out);
	printf("COMPLETED\n");

	fclose(fptr_r);
	fclose(fptr_g);
	fclose(fptr_b);
	fclose(fptr_p);
	fclose(fptr_out);
	return (EXIT_SUCCESS);
}
