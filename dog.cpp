#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
//#include <cmath>
#include<math.h>
#include <map>
#define PI 3.14
#define E 2.718281828

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below are two helper functions that overlay rectangles / circles 
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// Definition for a 2D pixel point.
// Add additional fields if necessary.
// 
struct Coordinate
{  
  int row, col;
};

// Definition for a descriptor.
// Add additional fields if necessary, and change DESCRIPTOR_LENGTH depending on the design of your descriptor.
// 
const int DESCRIPTOR_LENGTH = 128; 
struct Descriptor 
{
  double vector[DESCRIPTOR_LENGTH];
};

// Find corners in an image.
//
//%%%%%%%%%%%%%%%%%%%%%%
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  int rows=input.rows(), cols=input.cols(),i,j,k;
  SDoublePlane output(rows,cols);
  int kernel_size = row_filter.cols();
  int mid_row = floor(kernel_size/2);
  
  //create a temporary placeholder
  double ** temp_holder = new double *[rows];
  for(i=0;i<rows;i++)
	temp_holder[i] = new double[cols];


  //convolve the image with the 1-d vertical kernel
  for(i=0;i<rows;i++)
{
	for(j=0;j<cols;j++)
	{
		temp_holder[i][j] = 0;
		for(k=0;k<kernel_size;k++)
		{
			int mask_row = kernel_size - 1 - k;
			int image_row = i + k - mid_row;
			if(image_row >= 0 && image_row < rows && j >= 0 && j < cols)
			{
				temp_holder[i][j] += double(input[image_row][j] * col_filter[mask_row][0]);
			}
		}
		
  	}
  }

  //convolve the image with the 1-d horizontal kernel
  for(i=0;i<rows;i++){
	for(j=0;j<cols;j++){
		output[i][j] = 0;
		for(k=0;k<kernel_size;k++){
			int mask_col = kernel_size - 1 - k;
			int image_col = j + k - mid_row;
			if(image_col >= 0 && image_col < cols && i >= 0 && i < rows){
				output[i][j] += temp_holder[i][image_col] * row_filter[0][mask_col];				
			}
		}
	}
  }
  
  return output;
}
/////////////////////////////
SDoublePlane gaussian_filter(const SDoublePlane &input, double sigma,int GAUSS_KERNEL_SIZE)
{
  int kernel_size = GAUSS_KERNEL_SIZE, rows = input.rows(), cols = input.cols(),i,j,k,l;
  SDoublePlane output(rows,cols), row_filter(1,kernel_size), col_filter(kernel_size,1);


  double const1 = 1/(2 * PI * sigma * sigma);
  double const2 = -1/(2 * sigma * sigma);
  int mid_row = floor(kernel_size/2);
  
  
  
  int x = -mid_row;
  for(i=0;i<kernel_size;i++){
	col_filter[i][0] =  const1 * pow(E, (const2 * (pow(x,2) + pow(i+x,2))));
  }
  
  int y = 0;
  for(i=0;i<kernel_size;i++){
 	row_filter[0][i] =  const1 * pow(E, (const2 * (pow(i+x,2) + pow(y,2))));
  }

  double value = row_filter[0][0];
  for(i=0;i<kernel_size;i++){
	row_filter[0][i] = row_filter[0][i]/value;
  }
 
  output = convolve_separable(input,row_filter, col_filter);
  return output;
}


/////////////////////////////
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  int rows = input.rows(), cols = input.cols(),i,j,k;
  SDoublePlane sobel_row(1,3);
  SDoublePlane sobel_col(3,1);

  SDoublePlane output(rows, cols);

  //vertical or horizontal
  if (_gx==false){
/*	sobel_col[0][0] = -1;
	sobel_col[1][0] = -2;
	sobel_col[2][0] = -1;

	sobel_row[0][0] = -1;
	sobel_row[0][1] = 0;
	sobel_row[0][2] = 1;*/
sobel_col[0][0] = 1;
	sobel_col[1][0] = 1;
	sobel_col[2][0] = 1;

	sobel_row[0][0] = -1;
	sobel_row[0][1] = 0;
	sobel_row[0][2] = 1;

  }else{
/*	sobel_col[0][0] = -1;
	sobel_col[1][0] = 0;
	sobel_col[2][0] = 1;

	sobel_row[0][0] = -1;
	sobel_row[0][1] = -2;
	sobel_row[0][2] = -1;*/
sobel_col[0][0] = -1;
	sobel_col[1][0] = 0;
	sobel_col[2][0] = 1;

	sobel_row[0][0] = 1;
	sobel_row[0][1] = 1;
	sobel_row[0][2] = 1;
  
}

  //convolution function
  output = convolve_separable(input, sobel_row, sobel_col);
  return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &gradient_x, SDoublePlane &gradient_y, double thresh)
{
  int rows = gradient_x.rows(), cols = gradient_y.cols(), i, j, k;
  SDoublePlane output(rows, cols);

  for(i=0;i<rows;i++){
	for(int j=0;j<cols;j++){
		output[i][j]=sqrt(pow(gradient_x[i][j],2)+pow(gradient_y[i][j],2));
		if (output[i][j]>=thresh)
			output[i][j] = 0;
		else
			output[i][j] = 255;
	}
  }

  for(i=0;i<rows;i++){
	for(j=0;j<cols;j++){
		if(i<rows-1 && i>0 && j<cols-1 && j>0){
			if(output[i][j]>output[i+1][j+1] && output[i][j]>=output[i-1][j-1])
				output[i][j]=0;
			else
				output[i][j]=255;
		}
	}
  }

  return output;
}
vector<Coordinate> find_corners(const SDoublePlane &input,const SDoublePlane &input_x,const SDoublePlane &input_y)
{
//SImageIO::write_png_file("input.png", input , input , input);
Coordinate c;
c.row=3;
c.col=3;
SDoublePlane temp(input.rows(),input.cols());
SDoublePlane output(input.rows(),input.cols());
vector<Coordinate> points;
points.push_back(c);
//SDoublePlane ysquare(input.rows(),input.cols());
//double R=0;
//for(double sigma=0.25;sigma<=5;sigma+=0.25)
//{
temp = gaussian_filter(input,0.5,3);
SImageIO::write_png_file("temp.png", temp , temp , temp);
for(int i=0;i<input.rows();i++)
{
for(int j=0;j<input.cols();j++)
{
output[i][j]=(input[i][j]-temp[i][j]);
}
}


SImageIO::write_png_file("output.png", output , output , output);
return points;

//SImageIO::write_png_file("harris.png", gaussxsquare , gaussxsquare , gaussxsquare);
//SImageIO::write_png_file("harris1.png", gaussysquare , gaussysquare , gaussysquare);
//SImageIO::write_png_file("harris2.png", gaussxy , gaussxy , gaussxy);


cout<<"hi jayesh";
return points;
//return output;


}

// Compute an invarient descriptor for each feature point.
//
vector<Descriptor> invariant_descriptors(const SDoublePlane &input, const vector<Coordinate> &points)
{
  // (This placeholder code just returns a set of random descriptors.)
  vector<Descriptor> descriptors;
  for(int i = 0; i < points.size(); i++) {
    Descriptor descriptor;
    for(int j = 0; j < DESCRIPTOR_LENGTH; j++)
	  descriptor.vector[j] = double(rand()) / RAND_MAX;
      descriptors.push_back(descriptor);
  }

  return descriptors;
}

// Estimate a relative transation given two sets of interest points.
//
pair<int,int> translation_estimation(const SDoublePlane image1, const SDoublePlane image2, const vector<Descriptor> &descriptors1, const vector<Descriptor> &descriptors2)
{
  // Find correlated interest points from two point sets, 
  // then use RANSAC to estimate the homography matrix
  // (This placeholder returns a translation of 0)
  
  return make_pair(0,0);
}

// Match the query image with a retrieval set.
/*void stitch(SDoublePlane &image1, SDoublePlane &image2, string output_filename)
{
  // process the first image
  vector<Coordinate> image1_coordinates = find_corners(image1);
  vector<Descriptor> image1_descriptors = invariant_descriptors(image1, image1_coordinates);
  ofstream ofs(output_filename.c_str());
  
  // output a debugging image showing position of corners in the first image.
  SDoublePlane green_plane = image1;
  for(int i = 0; i<image1_coordinates.size(); i++) {
    int top = image1_coordinates[i].row, left = image1_coordinates[i].col;
    overlay_rectangle(green_plane, top, left, top+6, left+6, 255, 3);
  }
  SImageIO::write_png_file("debug.png", image1, green_plane, image1);
  
  // process the second image
  vector<Coordinate> image2_coordinates = find_corners(image2);
  vector<Descriptor> image2_descriptors = invariant_descriptors(image2, image2_coordinates);
	
  // estimate a relative translation
  pair<int,int> translation = translation_estimation(image1, image2, image1_descriptors, image2_descriptors);

  // Now stitch the two images together.
  // (This sample code just averages the two images together without applying a translation -- you'll
  //  need to fix this.)
  
  SDoublePlane output = image1;

  for(int i=0; i<image1.rows(); i++)
    for(int j=0; j<image1.cols(); j++)
      output[i][j] = (image1[i][j] + image2[i][j])/2;
  SImageIO::write_png_file(output_filename.c_str(), output, output, output);
}*/

int main(int argc, char *argv[])
{
float sigma_value,thresh_value;
int GAUSS_KERNEL_SIZE;
  if((argc < 4)) {
    cerr << "usage: " << argv[0] << " output_file image_file1 image_file2" << endl;
    return 1;
  }
double counter=0;
vector<Coordinate> output_points;
string image1_filename(argv[2]);
  string image2_filename(argv[3]);
  SDoublePlane image1= SImageIO::read_png_file(image1_filename.c_str());
  SDoublePlane image2= SImageIO::read_png_file(image2_filename.c_str());

  
  if(argc >= 5)
  	sigma_value = atof(argv[4]);
if(argc >= 6)
  	 GAUSS_KERNEL_SIZE = atoi(argv[5]);
if(argc >= 7)
  	thresh_value = atof(argv[6]);

  int image_count = argc - 3;
  string output_filename = argv[1];
   
  //out of the gray scale image
 // SImageIO::write_png_file("image1.png", image1, image1, image1);

  // subsample the input image
 cout<<"ok";
  //Call to Gaussian filter to smooth the image
  SDoublePlane gaussed = gaussian_filter(image1,sigma_value,GAUSS_KERNEL_SIZE);
  SImageIO::write_png_file("gaussed.png", gaussed , gaussed , gaussed);
  
  
  // compute gradient magnitude map of the input image
  //SDoublePlane gradient_x = sobel_gradient_filter(gaussed, true);
  //SImageIO::write_png_file("gradient_x.png", gradient_x, gradient_x, gradient_x);

  //SDoublePlane gradient_y = sobel_gradient_filter(gaussed, false);
  //SImageIO::write_png_file("gradient_y.png", gradient_y, gradient_y, gradient_y);
  
  // find edges in the input image
  //SDoublePlane edges = find_edges(gradient_x, gradient_y, thresh_value);
  //SImageIO::write_png_file("edges.png", edges, edges, edges);
  // read in query image and compute descriptors
cout<<"upto here";
find_corners(gaussed,gaussed,gaussed);

  
  // image matching
  //stitch(image1, image2, output_filename);
  
  return 0;
}
