#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include<math.h>
#include <map>
#define PI 3.14
#define E 2.718281828
using namespace std;
////
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

const int DESCRIPTOR_LENGTH = 128; 
struct Descriptor 
{
  double vector[DESCRIPTOR_LENGTH];
};
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
vector<Coordinate> find_corners(const SDoublePlane &input,const SDoublePlane &input_x,const SDoublePlane &input_y,double sigma,int size,string filename)
{
// Implement a Harris corner detector
  // (This placeholder code just returns a set of random points.)
  
vector<Coordinate> points;
     Coordinate c;
  
SDoublePlane xsquare(input.rows(),input.cols());
SDoublePlane ysquare(input.rows(),input.cols());
SDoublePlane xy(input.rows(),input.cols());
SDoublePlane gaussxsquare(input.rows(),input.cols());
SDoublePlane gaussysquare(input.rows(),input.cols());
SDoublePlane gaussxy(input.rows(),input.cols());
SDoublePlane output(input.rows(),input.cols());
double R=0;

for(int i=0;i<input.rows();i++)
{
for(int j=0;j<input.cols();j++)
{
xsquare[i][j]=pow(input_x[i][j],2);
ysquare[i][j]=pow(input_y[i][j],2);
xy[i][j]=(input_x[i][j]*input_y[i][j]);
}
}
gaussxsquare=gaussian_filter(xsquare,sigma,size);
gaussysquare=gaussian_filter(ysquare,sigma,size);
gaussxy=gaussian_filter(xy,sigma,size);
for(int i=0;i<input.rows();i++)
{
for(int j=0;j<input.cols();j++)
{
R=((gaussxsquare[i][j]*gaussysquare[i][j]) - (pow(gaussxy[i][j],2)))- (      (0.01)*     pow((gaussxsquare[i][j]+gaussysquare[i][j]),2)                );
if(R>3000)
{
c.row=i;
c.col=j;
points.push_back(c);
output[i][j]=sqrt(pow(input_x[i][j],2)+pow(input_y[i][j],2));

}
}
}
for(int i=0;i<output.rows();i++){
	for(int j=0;j<output.cols();j++){
		if(i<output.rows()-1 && i>0 && j<output.cols()-1 && j>0){
			if(output[i][j]>output[i+1][j+1] && output[i][j]<=output[i-1][j-1])
				{output[i][j]=output[i][j];}
			else
{		output[i][j]=0;}
		
}
	}
  }
SImageIO::write_png_file("filename.c_str", output , output , output);



SImageIO::write_png_file("harris.png", gaussxsquare , gaussxsquare , gaussxsquare);
SImageIO::write_png_file("harris1.png", gaussysquare , gaussysquare , gaussysquare);
SImageIO::write_png_file("harris2.png", gaussxy , gaussxy , gaussxy);


  
   
  return points;

}

// Compute an invarient descriptor for each feature point.

double* thresh(double r,double l, double t, double b,double *tempa)
{
         double threshold =0;
         double theta=0;
         tempa[0]=sqrt       (         pow((r-l),2)+pow((b-t),2)                 );
tempa[1] = atan (      (r-l)/(b-t)          )* double(180/PI);
//cout<<tempa[0]<<"  temps "<<tempa[1]<<endl;
//cout<<"now returning"<<endl;
return tempa;
         }
      
vector<Descriptor> invariant_descriptors(const SDoublePlane &input, const vector<Coordinate> &points)
{
//cout<<input.rows()<<" "<<input.cols()<<endl<<"  "<<points.size()<<endl;
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
Descriptor descriptor;   
for(int i=0;i<points.size();i++)
{
cout<<"point i "<<i<<endl;
arrsize=0;
    // if(points.at(i).row>8 && points.at(i).col >8 && points.at(i).row<input.rows()-8 && points.at(i).col<input.cols()-8)
if(i>1863)	
 {
goto data;
}
if(points.at(i).row>8 && points.at(i).col >8 && points.at(i).row<input.rows()-8 && points.at(i).col<input.cols()-8)
{
for(int e=points.at(i).row-8;e<=points.at(i).row+8;e+=4)
{
cout<<e-8<<endl;
      //if(e==points.at(i).row+5)
	  //{
	  //e=e-2;
	  //}
 if(e==points.at(i).row)
 {   
       e+=1;
       }

for(int r=points.at(i).col-8;r<=points.at(i).col+8;r+=4)
{
cout<<r-8<<endl;
 //if(r==points.at(i).col+3)
	//  {
	  //r=r-2;
	  //}
        if(r==points.at(i).col)
        {     r+=1;
              }
   

 for(int k=0;k<8;k++)
{
vectemp[k]=0;
}



for(int a=e;a<e+4;a++)
{
for(int b=r;b<r+4;b++)
{
if(a<0 || a>input.rows()-1 || b<0 || b>input.cols()-1)
       {
              tempb[0]=0;
              tempb[1]=0;
              }
else if(a-1 <0 && b-1<0)
{
cout<<a<<"a "<<b+1<<"b+1"<<a+1<<"a+1"<<b<<"b"<<endl;
 tempb = thresh(input[a][b+1],0,0,input[a+1][b],tempa);
          
       }
       else if(a >input.rows()-2&& b>input.cols()-2)
{
cout<<a<<"a "<<b-1<<"b-1"<<a-1<<"a-1"<<b<<"b"<<endl;

       tempb = thresh(0,input[a][b-1],input[a-1][b],0,tempa);
          
       }
       else if(a-1 <0 && b>input.cols()-2)
{
cout<<a<<"a "<<b-1<<"b-1"<<a+1<<"a+1"<<b<<"b"<<endl;

        tempb = thresh(0,input[a][b-1],0,input[a+1][b],tempa);
       }
       else if(a >input.rows()-2 && b-1<0)
{
cout<<a<<"a "<<b+1<<"b+1"<<a-1<<"a-1"<<b<<"b"<<endl;

        tempb = thresh(input[a][b+1],0,input[a-1][b],0,tempa);
          
       }
       else if(a-1 <0)
{
cout<<a<<"a "<<b-1<<"b-1"<<a+1<<"a+1"<<b<<"b"<<endl;

        tempb = thresh(input[a][b+1],input[a][b-1],0,input[a+1][b],tempa);
          
       }
       else if(a >input.rows()-2)
{
cout<<a<<"a "<<b+1<<"b+1"<<a-1<<"a-1"<<b<<"b"<<endl;

        tempb = thresh(input[a][b+1],input[a][b-1],input[a-1][b],0,tempa);
          
       }
       else if(b-1<0)
{
cout<<a-1<<"a-1 "<<b<<"b"<<a+1<<"a+1"<<b<<"b"<<endl;

        tempb = thresh(input[a][b+1],0,input[a-1][b],input[a+1][b],tempa);
          
       }
       
       else if(b>input.cols()-2)
{
cout<<a-1<<"a-1 "<<b<<"b"<<a+1<<"a+1"<<b<<"b"<<endl;

        tempb = thresh(0,input[a][b-1],input[a-1][b],input[a+1][b],tempa);
       }
       else if(b>=0 && b<input.cols()-1 && a>=0 && a<input.rows()-1)
       {
cout<<a<<"a "<<b+1<<"b+1 "<<a-1<<"a-1 "<<b<<"b "<<b-1<<"b-1 "<<a+1<<"a+1"<<endl;

           tempb = thresh(input[a][b+1],input[a][b-1],input[a-1][b],input[a+1][b],tempa);
}
       if(tempb[1]<0)
       {
                     tempb[1]+=360;
                    
                     }
 

      distance = sqrt       (         pow((points.at(i).col-b),2)+pow((points.at(i).row-a),2)                 );

tempb[0]=tempb[0]/distance;
if(tempb[0]!=0)
{
rem = (int)(tempb[1]/45);
cout<<rem<<"rem"<<endl;
vectemp[rem] += tempb[0];

}
}
}


for(int l=0;l<8;l++)
{
cout<<"upto this"<<endl;         
descriptor.vector[arrsize++]=vectemp[l];
  
 }
  
        }
        
        }
       
        
        
        }
descriptors.push_back(descriptor);

}
//cout<<"oops";
 data:
cout<<"retiurn";
 return descriptors;

}

// Estimate a relative transation given two sets of interest points.
////
struct imagedetails
{
int xtrans;
int ytrans;
int pointimage1;
int pointimage2;
int intensitydiff;
};
imagedetails translation_estimation(const SDoublePlane image1, const SDoublePlane image2, const vector<Descriptor> &descriptors1, const vector<Descriptor> &descriptors2,const vector<Coordinate> &image1_coordinates, const vector<Coordinate> &image2_coordinates)
{  

Coordinate array[image1.rows()][image1.cols()];
 int totalvalue=0;

int max=-10000,first,second;
int min=10000;
imagedetails ims;
  int i,j,k,l;
  vector<Coordinate> possible_translations;
int translatecount[image1.rows()][image1.cols()];
long int actual[image1.rows()][image1.cols()];


for(int q=0;q<image1.rows();q++)
{
for(int w=0;w<image1.cols();w++)
{
translatecount[q][w]=0;
}
}


int d;
int xtrans,ytrans,des1,des2,image1_x,image1_y,image2_x,image2_y;
  Coordinate translation;
  for(i=0;i<8000/*descriptors1.size()*/;i++){
	for(j=0;j<8000/*descriptors2.size()*/;j++){
		 d=0;



		for(k=0;k<128;k++){

			d+=abs(descriptors1.at(i).vector[k]-descriptors2.at(j).vector[k]);
	
		}

if(d<min )
{
min=d;
xtrans=abs(image1_coordinates.at(i).col-image2_coordinates.at(j).col);
ytrans=abs(image1_coordinates.at(i).row-image2_coordinates.at(j).row);
des1=i;
des2=j;
translatecount[xtrans][ytrans]++;
array[xtrans][ytrans].row=i;
array[xtrans][ytrans].col=j;
cout<<" "<<image1_coordinates.at(i).row<<" "<<image1_coordinates.at(i).col<<endl;
cout<<" "<<image2_coordinates.at(j).row<<" "<<image2_coordinates.at(j).col<<endl;
ims.xtrans=xtrans;
ims.ytrans=ytrans;
ims.pointimage1=des1;
ims.pointimage2=des2;
ims.intensitydiff=image1[image1_coordinates.at(i).row][image1_coordinates.at(i).col] - image2[image2_coordinates.at(j).row][image2_coordinates.at(j).col]; 
}

		
		}



 }

for(int q=0;q<image1.rows();q++)
{
for(int w=0;w<image1.cols();w++)
{
if(translatecount[q][w]>max)
{

max=translatecount[q][w];


first=q;
second=w;
}


}
}
cout<<"this is maximum"<<max<<" "<<first<<" "<<second<<endl;


int in=array[first][second].row;

int jn=array[first][second].col;
cout<<image1_coordinates.at(des1).row<<" "<<image1_coordinates.at(des1).col<<endl;
cout<<image2_coordinates.at(des2).row<<" "<<image2_coordinates.at(des2).col<<endl;


return ims;


}

// Match the query image with a retrieval set.
void stitch(SDoublePlane &image1, SDoublePlane &image2,SDoublePlane &gradient_x1,SDoublePlane &gradient_y1,SDoublePlane &gradient_x2,SDoublePlane &gradient_y2,string output_filename,string image1_filename,string image2_filename)
{
int xtrans,ytrans,ptimage1,ptimage2,xshift,yshift;
int image1xlow,image1ylow,image1yhigh,image1xhigh,image2xlow,image2ylow,image2yhigh,image2xhigh;
imagedetails ims;
int min=10000;
int max=-10000;
int left,right;
 int intdiff;
vector<Coordinate> image1_coordinates = find_corners(image1,gradient_x1,gradient_y1,1,3,image1_filename);
vector<Descriptor> image1_descriptors =invariant_descriptors(image1, image1_coordinates);
ofstream ofs(output_filename.c_str());
vector<Coordinate> image2_coordinates = find_corners(image2,gradient_x2,gradient_y2,1,3,image2_filename);
vector<Descriptor> image2_descriptors = invariant_descriptors(image2, image2_coordinates);
int variable=0;
int counter=0;
int maximum;
int count1=0;
double d;
cout<<"bala";
ims = translation_estimation(image1, image2, image1_descriptors, image2_descriptors, image1_coordinates, image2_coordinates);
  cout<<"bala";

xtrans=ims.xtrans;
ytrans=ims.ytrans;
ptimage1=ims.pointimage1;
ptimage2=ims.pointimage2;
intdiff=ims.intensitydiff;
cout<<xtrans<<" "<<ytrans<<" "<<ptimage1<<" "<<ptimage2<<endl;

xshift = abs(image1_coordinates.at(ptimage1).col-image2_coordinates.at(ptimage2).col);
yshift = abs(image1_coordinates.at(ptimage1).row-image2_coordinates.at(ptimage2).row);


if(image2_coordinates.at(ptimage2).row>image1_coordinates.at(ptimage1).row)
{

 image1ylow = 0; //to length-ytrans
 image1yhigh = image1.rows()-1-yshift;
 image2ylow=yshift;
 image2yhigh=image2.rows()-1;

}
else
{
 image2ylow = 0; //to length-ytrans
 image2yhigh = image2.rows()-1-yshift;
 image1ylow=yshift;
 image1yhigh=image1.rows()-1;


}


if(image1_coordinates.at(ptimage1).col>image2_coordinates.at(ptimage2).col)
{
//variable=image1xhigh;
 image1xlow = 0; //to length-ytrans
 image1xhigh = image1_coordinates.at(ptimage1).col-1;
int jayesh = image1_coordinates.at(ptimage1).col;
 image2xlow=image2_coordinates.at(ptimage2).col;

image2xhigh=image2.cols()-1;

}
else
{
//variable=image2xhigh;
 image1xlow =image1_coordinates.at(ptimage1).col; //to length-ytrans
 image1xhigh = image1.cols()-1;
 image2xlow=0;
 image2xhigh=image2_coordinates.at(ptimage2).col-1;

}



SDoublePlane stitched(image1yhigh-image1ylow+1,image2xhigh-image2xlow+1+image1xhigh-image1xlow+1);

int jayesh = image1_coordinates.at(ptimage1).col;
int jayesh1 = image2_coordinates.at(ptimage2).col;
if(image1_coordinates.at(ptimage1).row<image2_coordinates.at(ptimage2).row)
{
cout<<"ouch";
for(int i=image1ylow;i<=image1yhigh;i++)
{
for(int j=image1xlow;j<=image1xhigh;j++)
{

stitched[i][j]=image1[i][j];
}
}
}
else
{
for(int i=image2ylow;i<=image2yhigh;i++)
{
for(int j=image2xlow;j<=image2xhigh;j++)
{

stitched[i][j]=image2[i][j];
}
}
}
if(image1_coordinates.at(ptimage1).col>image2_coordinates.at(ptimage2).col)
{
for(int i=image2ylow;i<=image2yhigh;i++)
{
for(int j=0;j<=image2xhigh-image2xlow;j++)
{

stitched[i-yshift][jayesh+j]=image2[i][jayesh1+j];
}
}
}
else
{
for(int i=image1ylow;i<=image1yhigh;i++)
{
for(int j=0;j<=image1xhigh-image1xlow;j++)
{

stitched[i-yshift][jayesh1+j]=image1[i][jayesh+j];
}
}

}

SImageIO::write_png_file(output_filename.c_str(), stitched , stitched , stitched);



}

int main(int argc, char *argv[])
{
float sigma_value,thresh_value;
int GAUSS_KERNEL_SIZE;
  if((argc < 4)) {
    cerr << "usage: " << argv[0] << " output_file image_file1 image_file2" << endl;
    return 1;
  }
vector<Coordinate> output_points1;
vector<Coordinate> output_points2;

string image1_filename(argv[2]);
  string image2_filename(argv[3]);
  string output_filename(argv[1]);
SDoublePlane image1= SImageIO::read_png_file(image1_filename.c_str());
  SDoublePlane image2= SImageIO::read_png_file(image2_filename.c_str());

  
  if(argc >= 5)
  	sigma_value = atof(argv[4]);
if(argc >= 6)
  	 GAUSS_KERNEL_SIZE = atoi(argv[5]);
if(argc >= 7)
  	thresh_value = atof(argv[6]);

  int image_count = argc - 3;
//  string output_filename = argv[1];
   
   SImageIO::write_png_file("image1.png", image1, image1, image1);
  SImageIO::write_png_file("image2.png", image2, image2, image2);
  //Call to Gaussian filter to smooth the image
  SDoublePlane gaussed1 = gaussian_filter(image1,sigma_value,GAUSS_KERNEL_SIZE);
  SImageIO::write_png_file("gaussed1.png", gaussed1 , gaussed1 , gaussed1);
  SDoublePlane gaussed2 = gaussian_filter(image2,sigma_value,GAUSS_KERNEL_SIZE);
  SImageIO::write_png_file("gaussed2.png", gaussed2 , gaussed2 , gaussed2);
SDoublePlane google(image1.rows(),image1.cols());

  
  // compute gradient magnitude map of the input image
  SDoublePlane gradient_x1 = sobel_gradient_filter(gaussed1, true);
  SImageIO::write_png_file("gradient_x1.png", gradient_x1, gradient_x1, gradient_x1);

  SDoublePlane gradient_y1 = sobel_gradient_filter(gaussed1, false);
  SImageIO::write_png_file("gradient_y1.png", gradient_y1, gradient_y1, gradient_y1);
  
 SDoublePlane gradient_x2 = sobel_gradient_filter(gaussed2, true);
  SImageIO::write_png_file("gradient_x2.png", gradient_x2, gradient_x2, gradient_x2);

  SDoublePlane gradient_y2 = sobel_gradient_filter(gaussed2, false);
  SImageIO::write_png_file("gradient_y2.png", gradient_y2, gradient_y2, gradient_y2);

 stitch(image1, image2,gradient_x1,gradient_y1,gradient_x2,gradient_y2,output_filename,image1_filename,image2_filename);
  
  return 0;
}
