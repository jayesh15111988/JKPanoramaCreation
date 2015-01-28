#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include<sys/time.h>
#define PI 3.14
#define E 2.718281828

using namespace std;

double harris_thresh = 0.05;
double harris_thresh2 = 100000;
int a=0;


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

//Convolve the filters with the image
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
  for(i=0;i<rows;i++){
	for(j=0;j<cols;j++){
		temp_holder[i][j] = 0;
		for(k=0;k<kernel_size;k++){
			int mask_row = kernel_size - 1 - k;
			int image_row = i + k - mid_row;
			if(image_row >= 0 && image_row < rows && j >= 0 && j < cols){
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

// Apply a Gaussian of the specified sigma to an image, and return the result
//
SDoublePlane gaussian_filter(const SDoublePlane &input, double sigma, int GAUSS_KERNEL_SIZE)
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

// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  int rows = input.rows(), cols = input.cols(),i,j,k;
  SDoublePlane sobel_row(1,3);
  SDoublePlane sobel_col(3,1);

  SDoublePlane output(rows, cols);

  //vertical or horizontal
  if (_gx==false){
	sobel_col[0][0] = 1;
	sobel_col[1][0] = 1;
	sobel_col[2][0] = 1;

	sobel_row[0][0] = -1;
	sobel_row[0][1] = 0;
	sobel_row[0][2] = 1;
  }else{
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
vector<Coordinate> find_corners(const SDoublePlane &input)
{
  /****************************************************************************
  This function will first smooth the image and then apply the Harris algorithm
  to detect the corners of the given image
  ****************************************************************************/
  int rows=input.rows(), cols=input.cols(), i,j,k;
  double R = 0;
  SDoublePlane xsquare(rows,cols);
  SDoublePlane ysquare(rows,cols);
  SDoublePlane xy(rows,cols);
  SDoublePlane gaussxsquare(rows,cols);
  SDoublePlane gaussysquare(rows,cols);
  SDoublePlane gaussxy(rows,cols);
  SDoublePlane output(rows,cols);
  vector<Coordinate> points;
  Coordinate c;

  //Call to Gaussian filter to smooth the image
  SDoublePlane gaussed = gaussian_filter(input,0.75,3);

  //compute x gradient magnitude map of the input image
  SDoublePlane input_x = sobel_gradient_filter(gaussed, true);

  //compute x gradient magnitude map of the input image
  SDoublePlane input_y = sobel_gradient_filter(gaussed, false);

  //compute Ix ^ 2, Iy ^ 2, Ixy
  for(i=0;i<rows;i++){
	for(j=0;j<cols;j++){
		xsquare[i][j] = pow(input_x[i][j],2);
		ysquare[i][j] = pow(input_y[i][j],2);
		xy[i][j] = (input_x[i][j]*input_y[i][j]);
	}
  }

  //convolve with the gaussian kernel
  gaussxsquare=gaussian_filter(xsquare,0.75,3);
  gaussysquare=gaussian_filter(ysquare,0.75,3);
  gaussxy=gaussian_filter(xy,0.75,3);

  //Corner detection with non maximal supression
  for(i=0;i<rows;i++){
	for(j=0;j<cols;j++){
		R = (gaussxsquare[i][j]*gaussysquare[i][j]) - (pow(gaussxy[i][j],2))- (harris_thresh * pow((gaussxsquare[i][j]+gaussysquare[i][j]),2));
		if(R > harris_thresh2){
			output[i][j]=sqrt(pow(input_x[i][j],2) + pow(input_y[i][j],2));
		}
	}
  }

  //Non maximal supression
  for(i=0;i<rows;i++){
	for(j=0;j<cols;j++){
		if(i<rows-1 && i>0 && j<cols-1 && j>0){
			if(output[i][j]<=output[i+1][j+1] && output[i][j]>output[i-1][j-1])
				output[i][j]=0;
			else{
				c.row=i;
				c.col=j;
				points.push_back(c);
			}
		}
	}
  }

  stringstream ss;
  ss<<"corners"<<++a<<".png";
  string str(ss.str());
  //SImageIO::write_png_file(str.c_str(), output , output , output);  
  return points;
}

//A function to help find the gradient
double* thresh(double r,double l, double t, double b,double *tempa){
  double threshold =0;
  double theta=0;
  tempa[0]=sqrt       (         pow((r-l),2)+pow((b-t),2)                 );
  tempa[1] = atan (      (r-l)/(b-t)          )* double(180/PI);
  return tempa;
}

// Compute an invarient descriptor for each feature point.
//
vector<Descriptor> invariant_descriptors(const SDoublePlane &input, const vector<Coordinate> &points)
{
  int arrsize=0;
  vector<Descriptor> descriptors;
  double *tempa = new double[2];
  double *tempb=new double[2];

  int count = 0;
  int counter=0;
  int rem=0;
  double theta = 0;
  double distance;
  double vectemp[8];
  int chy=0;

  for(int i=0;i<points.size();i++){    
	arrsize=0;
	Descriptor descriptor;        
	for(int e=points.at(i).row-8;e<=points.at(i).row+8;e+=4){
		if(e>points.at(i).row-1 && e<points.at(i).row+5){   
			e+=1;
	       }

		for(int r=points.at(i).col-8;r<=points.at(i).col+8;r+=4){
	       	if(r>points.at(i).col-1 && r<points.at(i).col+5){
				r+=1;
			}
			for(int k=0;k<8;k++){
				vectemp[k]=0;
			}

  for(int a=e;a<e+4;a++){
	for(int b=r;b<r+4;b++){
		if(a<0 || a>input.rows()-1 || b<0 || b>input.cols()-1){
	              tempb[0]=0;
       	       tempb[1]=0;
              }
		else if(a-1 <0 && b-1<0){
			tempb = thresh(input[a][b+1],0,0,input[a+1][b],tempa);
       	}
	       else if(a >input.rows()-2&& b>input.cols()-2){
		       tempb = thresh(0,input[a][b-1],input[a-1][b],0,tempa); 
       	}
	       else if(a-1 <0 && b>input.cols()-2){
		       tempb = thresh(0,input[a][b-1],0,input[a+1][b],tempa);
       	}
	       else if(a >input.rows()-2 && b-1<0){
		       tempb = thresh(input[a][b+1],0,input[a-1][b],0,tempa);
	       }
       	else if(a-1 <0){
		       tempb = thresh(input[a][b+1],input[a][b-1],0,input[a+1][b],tempa);
       	}
	       else if(a >input.rows()-2){
		       tempb = thresh(input[a][b+1],input[a][b-1],input[a-1][b],0,tempa); 
	       }
       	else if(b-1<0){
		       tempb = thresh(input[a][b+1],0,input[a-1][b],input[a+1][b],tempa);
	       }
	       else if(b>input.cols()-2){
		       tempb = thresh(0,input[a][b-1],input[a-1][b],input[a+1][b],tempa);
	       }
       	else if(b>=0 && b<input.cols()-1 && a>=0 && a<input.rows()-1){
		       tempb = thresh(input[a][b+1],input[a][b-1],input[a-1][b],input[a+1][b],tempa);
		}

	       if(tempb[1]<0){
	              tempb[1]+=360;             
              } 
		distance = sqrt       (         pow((points.at(i).col-b),2)+pow((points.at(i).row-a),2)                 );
		tempb[0]=tempb[0]/distance;

		if(tempb[0]!=0){
			rem = (int)(tempb[1]/45);
			vectemp[rem] += tempb[0];

		}
	}
  }

  			for(int l=0;l<8;l++){
				descriptor.vector[arrsize++]=vectemp[l];
  			}
  		}
  	}
       
        
descriptors.push_back(descriptor);
  }
  return descriptors;
}

// Estimate a relative transation given two sets of interest points.
//
pair<int,int> translation_estimation(const SDoublePlane image1, const SDoublePlane image2, const vector<Descriptor> &descriptors1, const vector<Descriptor> &descriptors2,const vector<Coordinate> &image1_coordinates, const vector<Coordinate> &image2_coordinates)
{
  int i,j,k,l;
  vector<Coordinate> possible_translations;
int translatecount[image1.rows()][image1.cols()];
int actual[image1.rows()][image1.cols()];


for(int q=0;q<image1.rows();q++)
{
for(int w=0;w<image1.cols();w++)
{
translatecount[q][w]=0;

}
}

  Coordinate translation;
  for(i=0;i<descriptors1.size();i++){
	for(j=0;j<descriptors2.size();j++){
		int d=0;
		for(k=0;k<128;k++){
			d+=sqrt(descriptors1.at(i).vector[k]-descriptors2.at(j).vector[k]);
			
		}
if(d<min)
{
min=d;
xtrans=image1_coordinates.at(i).col-image2_coordinates.at(j).col;
ytrans=image1_coordinates.at(i).row-image2_coordinates.at(j).row;
des1=i;
des2=j;
}

		
		}
translatecount[xtrans][ytrans]++;
actual[xtrans][ytrans]=i*10+j;
//transalation.row=i;
//translation.col=j;
//possible_translations.push_back(translation);

  }

for(int q=0;q<image1.rows();q++)
{
for(int w=0;w<image1.cols();w++)
{
if(translatecount[q][w]>max);
{
first=q;
second=w;
}


}
}

totalvalue=actual[first][second];
j=totalvalue%10;

i=totalvalue/10;





}

// Match the query image with a retrieval set.
void stitch(SDoublePlane &image1, SDoublePlane &image2, string output_filename)
{
  // process the first image
  vector<Coordinate> image1_coordinates = find_corners(image1);
  vector<Descriptor> image1_descriptors = invariant_descriptors(image1, image1_coordinates);
  ofstream ofs(output_filename.c_str());
  
 /* // output a debugging image showing position of corners in the first image.
  SDoublePlane green_plane = image1;
  for(int i = 0; i<image1_coordinates.size(); i++) {
    int top = image1_coordinates[i].row, left = image1_coordinates[i].col;
    overlay_rectangle(green_plane, top, left, top+6, left+6, 255, 3);
  }
  SImageIO::write_png_file("debug.png", image1, green_plane, image1);*/
  
  // process the second image
  vector<Coordinate> image2_coordinates = find_corners(image2);
  vector<Descriptor> image2_descriptors = invariant_descriptors(image2, image2_coordinates);

  // estimate a relative translation
  cout<<image1_descriptors.size()<<" "<<image1_coordinates.size()<<endl;
  exit(1);
  pair<int,int> translation = translation_estimation(image1, image2, image1_descriptors, image2_descriptors, image1_coordinates, image2_coordinates);

  // Now stitch the two images together.
  // (This sample code just averages the two images together without applying a translation -- you'll
  //  need to fix this.)
  
  SDoublePlane output = image1;

  for(int i=0; i<image1.rows(); i++)
    for(int j=0; j<image1.cols(); j++)
      output[i][j] = (image1[i][j] + image2[i][j])/2;
  SImageIO::write_png_file(output_filename.c_str(), output, output, output);
}

int main(int argc, char *argv[])
{
  if(!(argc == 4)) {
    cerr << "usage: " << argv[0] << " output_file image_file1 image_file2" << endl;
    return 1;
  }
  
  string output_filename = argv[1];

  string image_file1 = argv[2];
  string image_file2 = argv[3];
  vector<Coordinate> output_points1, output_points2;
  
  SDoublePlane image1= SImageIO::read_png_file(image_file1.c_str());
  SDoublePlane image2= SImageIO::read_png_file(image_file2.c_str());
  
  // image matching
  stitch(image1, image2, output_filename);
  
  return 0;
}
