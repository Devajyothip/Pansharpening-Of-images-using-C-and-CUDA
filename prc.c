#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Image structure
typedef struct {
	int width, height; // Image size
	int bytes_per_pixel; // 1 for grayscale image, 3 for rgb image
	unsigned long long total_bytes; // Total bytes in data, width * height * bytes_per_pixel
	unsigned char * data; // Image data - very large array of 8-bit values
	
	float mean, variance, stddev; // Image metrics, left unevaluated at the beginning
} image_t;

void alloc_image(image_t * image, int width, int height, int bytes_per_pixel) {
	// Allocate large chunk of memory for the image
	image->width = width;
	image->height = height;
	image->bytes_per_pixel = bytes_per_pixel;
	image->total_bytes = (unsigned long long) (width) * height * bytes_per_pixel;
	
	image->data = (unsigned char *) malloc(width * height * bytes_per_pixel);
	if (!image->data) {
		printf("Unable to allocate %llu bytes of memory!", image->total_bytes);
		exit(1); // Critical error for us
	}
	
	image->mean = 0.0f;
	image->variance = 0.0f;
	image->stddev = 0.0f;
	
	printf("%llu MiB allocated for the image\n",
	             (unsigned long long) (image->total_bytes + sizeof(image_t) + 1024 * 1024 - 1) / (1024 * 1024));
}

void dealloc_image(image_t * image) {
	// Free memory
	free(image->data);
	printf("%llu MiB deallocated\n",
	             (unsigned long long) (image->total_bytes + sizeof(image_t) + 1024 * 1024 - 1) / (1024 * 1024));
	
	image->data = NULL;
	image->width = 0;
	image->height = 0;
	image->bytes_per_pixel = 0;
	image->total_bytes = 0;
}

void clone_image(image_t * image, image_t * src) {
	// Create a full copy of an image src
	alloc_image(image, src->width, src->height, src->bytes_per_pixel);
	memcpy(image->data, src->data, src->total_bytes);
	
	image->mean = src->mean;
	image->variance = src->variance;
	image->stddev = src->stddev;
}

void read_image(image_t * image, const char * path) {
	FILE * fp = fopen(path, "rb");
	if (!fp) {
		printf("File %s not found or can't be opened! Exiting...", path); // Critical error for our program
		exit(1);
	}
	
	unsigned long off_set = 0, next_offset = 0, pos_count = 0, tag_value[25] = { 0 }, strip_offset_val = 0, strip_offset = 0;
	int i, j, k, tag_id[25] = { 0 }; //, tag_type[25] = { 0 }, tag_count[25] = { 0 };
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
		IFD_count++;														  // No. of IFD's in this TIFF File
		buffer[1] = getc(fp);
		buffer[0] = getc(fp);
		tagCount = ((int)buffer[0] << 8) | ((int)buffer[1]);	 // No. of entries in an IFD
		for (i = 0; i < tagCount; i++) {									// Read all the entries of this IFD
			buffer[1] = getc(fp);
			buffer[0] = getc(fp);
			tag_id[i] = ((int)buffer[0] << 8) | ((int)buffer[1]);		 // Tag ID

			for (j = 1; j >= 0; j--) {
				buffer[j] = getc(fp);
			}

			for (j = 3; j >= 0; j--) {
				buffer[j] = getc(fp);
			}
			
				  // Gives no. of values for this Tag
			for (j = 3; j >= 0; j--) {
				buffer[j] = getc(fp);
			}
			tag_value[i] = ((int)buffer[0] << 24) | ((int)buffer[1] << 16) | ((int)buffer[2] << 8) | ((int)buffer[3]);
				  // Gets the value if the above count is 1, else offset of the starting value
			if (tag_id[i] == 256)												 // Tag ID 256 says about the image width
				img_Width = tag_value[i];
			if (tag_id[i] == 257)												 // Tag ID 257 says about the image length
				img_Len = tag_value[i];
			if (tag_id[i] == 273)												 // Tag ID 273 says about the offset which points to the offset of strips
				strip_offset_val = tag_value[i];
		}	

		alloc_image(image, img_Width, img_Len, 1);
	
		for (i = 0; i < img_Len; i++) {										  // Read the pixel values from image and store it in the matrix
			fseek(fp, (strip_offset_val + (i * 4)), SEEK_SET);				  // Move to the offset of the current strip's offset

			for (j = 3; j >= 0; j--) {
				buffer[j] = getc(fp);
			}
			strip_offset = ((int)buffer[0] << 24) | ((int)buffer[1] << 16) | ((int)buffer[2] << 8) | ((int)buffer[3]);
			fseek(fp, strip_offset, SEEK_SET);	 // Move to the offset of the current strip
			for (j = 0; j < img_Width; j++) {

				getc(fp);
				short tmp_c2 = getc(fp);
				image->data[i * img_Width + j] = tmp_c2;
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
	
	fclose(fp);
	
	printf("Image %s loaded successfully\n", path);
}

void WriteHexString(FILE *fptr, char *s) {
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

void write_image(image_t * image, const char * path) {
	if ((image->bytes_per_pixel != 1) && (image->bytes_per_pixel != 3)) {
		printf("Only 1 and 3 bytes per pixel images are supported in write_image procedure");
		exit(1);
		return;
	}
	
	FILE * fptr = fopen(path, "wb");
	if (!fptr) {
		printf("File %s can't be opened for writing! Exiting...", path); // Critical error for our program
		exit(1);
	}
	
   /* Write the header */
   WriteHexString(fptr, "4d4d002a");    /* Big endian & TIFF identifier */
   int nx = image->width;
   int ny = image->height;
   int offset = nx * ny * 3 + 8;
   putc((offset & 0xff000000) / 16777216, fptr);
   putc((offset & 0x00ff0000) / 65536, fptr);
   putc((offset & 0x0000ff00) / 256, fptr);
   putc((offset & 0x000000ff), fptr);

   /* Write the binary data */
	unsigned long long i;
	
	if (image->bytes_per_pixel == 3) {
		// Just save the data "as is"
		for (i = 0; i < image->total_bytes; i++)
			fputc(image->data[i], fptr);
	} else {
		// Save each pixel three times as r, g, b component
		for (i = 0; i < image->total_bytes; i++) {
			fputc(image->data[i], fptr);
			fputc(image->data[i], fptr);
			fputc(image->data[i], fptr);
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
	
	fclose(fptr);
	printf("File %s written successfully.\n", path);
}

// Evaluate image mean and standart deviation
void eval_stats(image_t * image) {
	if (image->bytes_per_pixel != 1) {
		printf("Mean and standart deviation are only evaluated for grayscale images!\n");
		exit(1);
	}
	
	float mean = 0.0f, variance = 0.0f, stddev = 0.0f;
	unsigned long long i;
	
	// Evaluate mean
	for (i = 0; i < image->total_bytes; i++)
		mean += image->data[i];
	mean /= image->total_bytes;
	
	// Evaluate variance
	for (i = 0; i < image->total_bytes; i++) {
		float tmp = (float) (image->data[i]) - mean;
		variance += tmp * tmp;
	}
	variance /= image->total_bytes;
	
	stddev = sqrtf(variance);
	
	// Set up those values
	image->mean = mean;
	image->stddev = stddev;
}

unsigned char get_closest_point(image_t * src, unsigned long long idx, image_t * where) {
	int si = idx % src->width;
	int sj = idx / src->width;
	
	int wi = (si / (src->height - 1.0f)) * (where->height - 1.0f);
	int wj = (sj / (src->width - 1.0f))  * (where->width - 1.0f);
	
	return where->data[wj * where->width + wi];
}

// Resize image to new size
void resize(image_t * dst, image_t * src, int new_w, int new_h) {
	// Note that we can only resize grayscale images
	if (src->bytes_per_pixel != 1) {
		printf("Resizing is only implemented for grayscale images!\n");
		exit(1);
	}
	
	// First - allocate memory for the dst image
	alloc_image(dst, new_w, new_h, 1);
	
	// Aspect ratio should not be changed
	unsigned long long i;
	
	for (i = 0; i < dst->total_bytes; i++)
		dst->data[i] = get_closest_point(dst, i, src);
}

float covariance(image_t * b, image_t * gs) {
	// Estimate mean
	eval_stats(b);
	eval_stats(gs);
	
	float covariance = 0.0f;
	float variance = 0.0f;
	
	unsigned long long i;
	for (i = 0; i < b->total_bytes; i++) {
		covariance += (b->data[i] - b->mean) * (gs->data[i] - gs->mean);
		variance += (gs->data[i] - gs->mean) * (gs->data[i] - gs->mean);
	}
	
	// Both should be divided by N, but we're going to divide them anyway
	return covariance / variance;
}

// Core routine - Gram-Schmidt transformation
void GramSchmidtTransformation(image_t * gs, image_t * bands, float ** phi) {
	// Gram-Schmidt imlementation for 4 vectors
	
	// First gs element is the same as band 0, i.e. artificial low res pan image
	clone_image(gs + 0, bands + 0);
	alloc_image(gs + 1, bands[0].width, bands[0].height, bands[0].bytes_per_pixel);
	alloc_image(gs + 2, bands[0].width, bands[0].height, bands[0].bytes_per_pixel);
	alloc_image(gs + 3, bands[0].width, bands[0].height, bands[0].bytes_per_pixel);
	
	// For the rest three images, we need to follow modified Gram-Schmidt routine
	unsigned l, T;
	for (T = 1; T < 4; T++) {
		phi[T][0] = 0.0f;
		phi[T][1] = 0.0f;
		phi[T][2] = 0.0f;
		phi[T][3] = 0.0f;
		
		for (l = 0; l < T; l++) 
			phi[T][l] = covariance(bands + T, gs + l);
			
		unsigned long long i;
		for (i = 0; i < bands[T].total_bytes; i++) {
			float res = (bands[T].data[i] - bands[T].mean);
			for (l = 0; l < T; l++)
				res -= phi[T][l] * gs[l].data[i];
			
			gs[T].data[i] = roundf(res);
		}
	}
}

void GramSchmidtReverseTransformation(image_t * dst, image_t * gs, image_t * bands, float ** phi) {
	// Gram-Schmidt imlementation for 4 vectors
	
	// First gs element is the same as band 0, i.e. artificial low res pan image
	alloc_image(dst + 0, gs[0].width, gs[0].height, gs[0].bytes_per_pixel);
	alloc_image(dst + 1, gs[0].width, gs[0].height, gs[0].bytes_per_pixel);
	alloc_image(dst + 2, gs[0].width, gs[0].height, gs[0].bytes_per_pixel);
	alloc_image(dst + 3, gs[0].width, gs[0].height, gs[0].bytes_per_pixel);
	
	// For the rest three images, we need to follow modified Gram-Schmidt routine
	unsigned l, T;
	for (T = 0; T < 4; T++) {
		unsigned long long i;
		for (i = 0; i < dst[0].total_bytes; i++) {
			float res = (get_closest_point(dst + 0, i, gs + T) + bands[T].mean);
			for (l = 0; l < T; l++)
				res += phi[T][l] * get_closest_point(dst + 0, i, gs + l);
			
			dst[T].data[i] = roundf(res);
		}
	}
}

void coalesce(image_t * images, image_t * res) {
	alloc_image(res, images[1].width, images[1].height, 3);
	
	unsigned long long i;
	for (i = 0; i < images[1].total_bytes; i++) {
		res->data[3 * i + 0] = images[1].data[i];
		res->data[3 * i + 1] = images[2].data[i];
		res->data[3 * i + 2] = images[3].data[i];
	}
}

void stretch(image_t * dst, image_t * src) {
	// No need for dst and src to be of the same size
	eval_stats(dst);
	eval_stats(src);
	
	float gain = src->stddev / dst->stddev;
	float bias = src->mean - gain * dst->mean;
	
	unsigned long long i;
	for (i = 0; i < dst->total_bytes; i++)
		dst->data[i] = roundf(dst->data[i] * gain + bias);
}

int main(int argc, char * argv[]) {
	if (argc != 5){
		printf("%s <R-filename1> <G-filename2> <B-filename3> <Pan-filename4>\n", argv[0]);
		return 1;
	}
	
	// Original images
	image_t r, g, b, p;
	
	printf("--- Loading initial images...\n");
	
	read_image(&r, argv[1]); // r
	read_image(&g, argv[2]); // g
	read_image(&b, argv[3]); // b
	read_image(&p, argv[4]); // High res pan band
	
	if ((r.width != g.width) || (r.width != b.width) || (r.height != g.height) || (r.height != b.height)) {
		printf("red, green, blue images are not the same size!\n");
		return 2;
	}
	
	// Create low res pan band
	printf("--- Simulating low res pan band...\n");
	image_t sim_p;
	resize(&sim_p, &p, r.width, r.height);
	
	image_t gs[4];
	image_t bands[4]; // Original bands
	image_t out[4];
	bands[0] = sim_p;
	bands[1] = r;
	bands[2] = g;
	bands[3] = b;
	
	float ** phi;
	phi = malloc(4 * sizeof(float *));
	unsigned i;
	for (i = 0; i < 4; i++)
		phi[i] = calloc(4, sizeof(float));
	
	printf("--- Executing Gram-Schmidt transformation...\n");
	GramSchmidtTransformation(gs, bands, phi);
	
	printf("--- Stretching high res pan image...\n");
	stretch(&p, gs + 0);
	
	dealloc_image(gs + 0);
	gs[0] = p;
	
	printf("--- Executing inverse Gram-Schmidt transformation...\n");
	GramSchmidtReverseTransformation(out, gs, bands, phi);
	
	write_image(bands+0, "i_in_sim_p.tif");
	write_image(bands+1, "i_in_r.tif");
	write_image(bands+2, "i_in_g.tif");
	write_image(bands+3, "i_in_b.tif");
	
	write_image(gs+0, "i_gs0.tif");
	write_image(gs+1, "i_gs1.tif");
	write_image(gs+2, "i_gs2.tif");
	write_image(gs+3, "i_gs3.tif");
	
	// Deallocate all non-needed images here
	dealloc_image(bands + 0); // sim_p
	dealloc_image(bands + 1); // r
	dealloc_image(bands + 2); // g
	dealloc_image(bands + 3); // b
	dealloc_image(gs + 0); // Modified pan
	dealloc_image(gs + 1); // GS band 1
	dealloc_image(gs + 2); // GS band 2
	dealloc_image(gs + 3); // GS band 3
	
	for (i = 0; i < 4; i++)
		free(phi[i]);
	free(phi);
	
	printf("--- Coalescing image...\n");
	image_t res;
	coalesce(out, &res);
	
	write_image(out+0, "i_out_r.tif");
	write_image(out+1, "i_out_g.tif");
	write_image(out+2, "i_out_b.tif");
	write_image(out+3, "i_out_p.tif");
	
	dealloc_image(out + 0);
	dealloc_image(out + 1);
	dealloc_image(out + 2);
	dealloc_image(out + 3);
	
	printf("--- Writing image...\n");
	write_image(&res, "out.tif");
	
	dealloc_image(&res);
	printf("--- Everything is done!\n");
	
	return 0;
}
