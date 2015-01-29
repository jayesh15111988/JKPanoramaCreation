#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <ctime>
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

const int DESCRIPTOR_LENGTH = 128;
string folderName = "";
string currentDateAndTimeValue = "";

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
//
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
  for(i=0;i<rows;i++)
  {
  for(j=0;j<cols;j++)
  {
    output[i][j] = 0;
    for(k=0;k<kernel_size;k++)
    {
    int mask_col = kernel_size - 1 - k;
    int image_col = j + k - mid_row;
    if(image_col >= 0 && image_col < cols && i >= 0 && i < rows)
      {
        output[i][j] += temp_holder[i][image_col] * row_filter[0][mask_col];        
      }
    }
  }
  }
  return output;
}


SDoublePlane gaussian_filter(const SDoublePlane &input, double sigma,int GAUSS_KERNEL_SIZE)
{
  int kernel_size = GAUSS_KERNEL_SIZE, rows = input.rows(), cols = input.cols(),i,j,k,l;
  SDoublePlane output(rows,cols), row_filter(1,kernel_size), col_filter(kernel_size,1);
  double const1 = 1/(2 * PI * sigma * sigma);
  double const2 = -1/(2 * sigma * sigma);
  int mid_row = floor(kernel_size/2);
  int x = -mid_row;
  for(i=0;i<kernel_size;i++)
  {
    col_filter[i][0] =  const1 * pow(E, (const2 * (pow(x,2) + pow(i+x,2))));
  }
  int y = 0;
  for(i=0;i<kernel_size;i++)
  {
    row_filter[0][i] =  const1 * pow(E, (const2 * (pow(i+x,2) + pow(y,2))));
  }
  double value = row_filter[0][0];
  for(i=0;i<kernel_size;i++)
  {
    row_filter[0][i] = row_filter[0][i]/value;
  }
  output = convolve_separable(input,row_filter, col_filter);
  return output;
}


SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  int rows = input.rows(), cols = input.cols(),i,j,k;
  SDoublePlane sobel_row(1,3);
  SDoublePlane sobel_col(3,1);
  SDoublePlane output(rows, cols);

  //vertical or horizontal
  if (_gx==false)
  {
    sobel_col[0][0] = 1;
    sobel_col[1][0] = 1;
    sobel_col[2][0] = 1;

    sobel_row[0][0] = -1;
    sobel_row[0][1] = 0;
    sobel_row[0][2] = 1;

  }
  else
  {
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

//tomasi corner detector
vector<Coordinate> find_corners_tomasi(const SDoublePlane &input,const SDoublePlane &input_x,const SDoublePlane &input_y, string imageSequenceNumber)
{
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
    
  gaussxsquare=gaussian_filter(xsquare,2,5);
  gaussysquare=gaussian_filter(ysquare,2,5);
  gaussxy=gaussian_filter(xy,2,5);
  double trace1;
  double det1;
  double lambda1;
  double lambda2;
  for(int i=0;i<input.rows();i++)
  {
    for(int j=0;j<input.cols();j++)
    {
      trace1=(gaussxsquare[i][j]+gaussysquare[i][j]);
      det1 = (gaussxsquare[i][j]*gaussysquare[i][j]) - (pow(gaussxy[i][j],2));
      lambda1=(0.5)*((trace1)+sqrt(pow(trace1,2)-(4*det1)));
      lambda2=(0.5)*((trace1)-sqrt(pow(trace1,2)-(4*det1)));
      R=min(lambda1,lambda2);
      if(R>20)
      {
        c.row=i;
        c.col=j;
        points.push_back(c);
        output[i][j]=255;//sqrt(pow(input_x[i][j],2)+pow(input_y[i][j],2));
      }
    }
  }
  for(int i=0;i<output.rows();i++)
  {
    for(int j=0;j<output.cols();j++)
    {
      if(i<output.rows()-1 && i>0 && j<output.cols()-1 && j>0)
      {
        if(output[i][j]<=output[i+1][j+1] && output[i][j]>output[i-1][j-1])
        //{
        //output[i][j]=0;
      //  }
      //  else
        {//       
    output[i][j]=0;
        }
      }
    }
  }
    string generalFileName = currentDateAndTimeValue + "tomasi";
    string tomasiConvertedFileName =  generalFileName + "_" + imageSequenceNumber +  ".png";
    string tomasiImageStorageFullPath = folderName + tomasiConvertedFileName;
    SImageIO::write_png_file(tomasiImageStorageFullPath.c_str(), output , output , output);
  return points;
}


vector<Coordinate> find_corners_harris(const SDoublePlane &input,const SDoublePlane &input_x,const SDoublePlane &input_y,double sigma,int size,string filename,string turn)
{
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
      R=((gaussxsquare[i][j]*gaussysquare[i][j]) - (pow(gaussxy[i][j],2)))- ((0.01)*pow((gaussxsquare[i][j]+gaussysquare[i][j]),2)                );
      if(R>1000)
      {
        c.row=i;
        c.col=j;
        points.push_back(c);
        output[i][j]=sqrt(pow(input_x[i][j],2)+pow(input_y[i][j],2));
      }
    }
  }
    
  for(int i=0;i<output.rows();i++)
  {
    for(int j=0;j<output.cols();j++)
    {
      if(i<output.rows()-1 && i>0 && j<output.cols()-1 && j>0)
      {
        if(output[i][j]>output[i+1][j+1] && output[i][j]<=output[i-1][j-1])
        {
          output[i][j]=output[i][j];
        }
        else
        {       
          output[i][j]=255;
        }
      }
    }
  }
    
string harrisCornerDetectorOutputFile = folderName + currentDateAndTimeValue + "harris_corner_detector" + turn + ".png";

SImageIO::write_png_file(harrisCornerDetectorOutputFile.c_str(), output , output , output);

  return points;
}

// Compute an invarient descriptor for each feature point.

double* thresh(double r,double l, double t, double b,double *tempa)
{
  double threshold =0;
    double theta=0;
    tempa[0]=sqrt(pow((r-l),2)+pow((b-t),2));
  tempa[1] = atan((r-l)/(b-t))* double(180/PI);
  return tempa;
}
      
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
    int totalNumberOfPoints = points.size();
  for(long i=0;i<totalNumberOfPoints;i++)
  {
    arrsize=0;
    Descriptor descriptor;        
    for(int e=points.at(i).row-8;e<=points.at(i).row+8;e+=4)
      {   
        if(e>points.at(i).row-1 && e<points.at(i).row+5)
      {   
        e+=1;
      }

    for(int r=points.at(i).col-8;r<=points.at(i).col+8;r+=4)
      {
        if(r>points.at(i).col-1 && r<points.at(i).col+5)
      {     
        r+=1;
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
          tempb = thresh(input[a][b+1],0,0,input[a+1][b],tempa);
        }
        else if(a >input.rows()-2&& b>input.cols()-2)
        {
          tempb = thresh(0,input[a][b-1],input[a-1][b],0,tempa);
        }
        else if(a-1 <0 && b>input.cols()-2)
        {
          tempb = thresh(0,input[a][b-1],0,input[a+1][b],tempa);
        }
        else if(a >input.rows()-2 && b-1<0)
        {
          tempb = thresh(input[a][b+1],0,input[a-1][b],0,tempa);
        }
        else if(a-1 <0)
        {
          tempb = thresh(input[a][b+1],input[a][b-1],0,input[a+1][b],tempa);
        }
        else if(a >input.rows()-2)
        {
          tempb = thresh(input[a][b+1],input[a][b-1],input[a-1][b],0,tempa);
        }
        else if(b-1<0)
        {
          tempb = thresh(input[a][b+1],0,input[a-1][b],input[a+1][b],tempa);
        }
        else if(b>input.cols()-2)
        {
          tempb = thresh(0,input[a][b-1],input[a-1][b],input[a+1][b],tempa);
        }
        else if(b>=0 && b<input.cols()-1 && a>=0 && a<input.rows()-1)
        {
          tempb = thresh(input[a][b+1],input[a][b-1],input[a-1][b],input[a+1][b],tempa);
        }
        if(tempb[1]<0)
        {
          tempb[1]+=360;
                }
 

        distance = sqrt(pow((points.at(i).col-b),2)+pow((points.at(i).row-a),2));

        tempb[0]=tempb[0]/distance;
        if(tempb[0]!=0)
        {
          rem = (int)(tempb[1]/45);
          vectemp[rem] += tempb[0];

        }
      }
    }

  for(int l=0;l<8;l++)
  {     
    descriptor.vector[arrsize++]=vectemp[l];
  }
      }
        }
      
      descriptors.push_back(descriptor);
  }
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
int limit1,limit2;
int xtrans,ytrans,des1,des2;
Coordinate translation;
if(descriptors1.size()>12000)
{
  limit1=10000;
}
else
{
  limit1=descriptors1.size();
}
if(descriptors2.size()>12000)
{
  limit2=10000;
}
else
{
  limit2=descriptors2.size();
}
for(i=0;i<limit1;i++)
{
  for(j=0;j<limit2;j++)
  {
    d=0;
    for(k=0;k<128;k++)
    {
      d+=abs(descriptors1.at(i).vector[k]-descriptors2.at(j).vector[k]);
    }
    if(d<min)
    {
    min=d;
    xtrans=abs(image1_coordinates.at(i).col-image2_coordinates.at(j).col);
    ytrans=abs(image1_coordinates.at(i).row-image2_coordinates.at(j).row);
    des1=i;
    des2=j;
    translatecount[xtrans][ytrans]++;
    array[xtrans][ytrans].row=i;
    array[xtrans][ytrans].col=j;
    ims.xtrans=xtrans;
    ims.ytrans=ytrans;
    ims.pointimage1=des1;
    ims.pointimage2=des2;
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
return ims;
}

void stitch(SDoublePlane &image1, SDoublePlane &image2,SDoublePlane &gradient_x1,SDoublePlane &gradient_y1,SDoublePlane &gradient_x2,SDoublePlane &gradient_y2,string output_filename,string image1_filename,string image2_filename)
{
  int xtrans,ytrans,ptimage1,ptimage2,xshift,yshift;
  int image1xlow,image1ylow,image1yhigh,image1xhigh,image2xlow,image2ylow,image2yhigh,image2xhigh;
  imagedetails ims;
  int min=10000;
  int max=-10000;
  int left,right;
 
    cout<<"Finding corners using Harris corner detection algorithm for first image\n\n";
  vector<Coordinate> image1_coordinates = find_corners_harris(image1,gradient_x1,gradient_y1,1,3,image1_filename,"1");
    
    vector<Descriptor> image1_descriptors =invariant_descriptors(image1, image1_coordinates);
    ofstream ofs(output_filename.c_str());
    cout<<"Finding corners using Harris corner detection algorithm for second image\n\n";

  vector<Coordinate> image2_coordinates = find_corners_harris(image2,gradient_x2,gradient_y2,1,3,image2_filename,"2");
  vector<Descriptor> image2_descriptors = invariant_descriptors(image2, image2_coordinates);

  int variable=0;
  int counter=0;
  int maximum;
  int count1=0;
  double d;
  ims = translation_estimation(image1, image2, image1_descriptors, image2_descriptors, image1_coordinates, image2_coordinates);

  xtrans=ims.xtrans;
  ytrans=ims.ytrans;
  ptimage1=ims.pointimage1;
  ptimage2=ims.pointimage2;
  //cout<<xtrans<<" "<<ytrans<<" "<<ptimage1<<" "<<ptimage2<<endl;

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
    //int jayesh = image1_coordinates.at(ptimage1).col;
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
cout<<"Stitch Begin\n\n";
  SDoublePlane stitched(image1yhigh-image1ylow+1,image2xhigh-image2xlow+1+image1xhigh-image1xlow+1);
    cout<<"Stitch End\n\n";

  int jayesh = image1_coordinates.at(ptimage1).col;
  int jayesh1 = image2_coordinates.at(ptimage2).col;
  if(image1_coordinates.at(ptimage1).row<image2_coordinates.at(ptimage2).row)
  {
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


//Get current date and time for us
const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%m-%d-%Y.%X", &tstruct);
    return buf;
}

bool doesStringContainCharacter(string inputString, string characterToFind) {
    return (inputString.find(characterToFind) != std::string::npos);
}

void split(vector<string>& tokens, const string &text, char sep) {
    int start = 0, end = 0;
    while ((end = text.find(sep, start)) != string::npos) {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
}

int main(int argc, char *argv[])
{
  float sigma_value,thresh_value;
  int GAUSS_KERNEL_SIZE;
  if((argc < 4)) 
  {
    cerr << "usage: " << argv[0] << "mode(stitch(1)/no stitch(0)) output_file image_file1 image_file2" << endl;
    return 1;
  }
    


  vector<Coordinate> output_points1;
  vector<Coordinate> output_points2;

  string image1_filename(argv[3]);
  string image2_filename(argv[4]);
  string output_filename(argv[2]);

  int toStitchImages = atoi(argv[1]);
    
  const char* firstInputImageName = image1_filename.c_str();
  const char* secondInputImageName = image2_filename.c_str();

    
  if(!doesStringContainCharacter(firstInputImageName,".png")) {
    cerr << "Please provide only png files as an input" << endl;
    return 1;
  }
    
    if(toStitchImages == 1) {
        if(!doesStringContainCharacter(secondInputImageName,".png")) {
            cerr << "Please provide only png files as an input" << endl;
            return 1;
        }
    }
    
    //Start the clock to measure execution time
    clock_t begin = clock();
    
    //We will check if image to be processed is taken from folder
    vector<string> tokensCollection;
    string fileName = "";
    
    if(doesStringContainCharacter(firstInputImageName, "/")) {
        split(tokensCollection, firstInputImageName, '/');
        int totalLength = tokensCollection.size();
        for (int i=0; i<totalLength - 1;i++){
            folderName.append(tokensCollection[i]);
            folderName.append("/");
        }
    }
    
    cout<<"Full folder path "<<folderName<<endl;

  SDoublePlane image1= SImageIO::read_png_file(firstInputImageName);
  
    int extraParametersIndex = 4 + toStitchImages;
    
  if(argc >= extraParametersIndex + 1) {
    sigma_value = atof(argv[extraParametersIndex]);
  }
  
  if(argc >= extraParametersIndex + 2) {
    GAUSS_KERNEL_SIZE = atoi(argv[extraParametersIndex + 1]);
  }

  if(argc >= extraParametersIndex + 3) {
    thresh_value = atof(argv[extraParametersIndex + 2]);
  }
    
    cout<<"Output image name : "<<output_filename<<"\nsigma value is : "<<sigma_value<<"\nGauss Kernel size is : "<<GAUSS_KERNEL_SIZE<<"\nThreshold value for detecting image corners is : "<<thresh_value<<endl<<endl;
    
    currentDateAndTimeValue = currentDateTime() + "_";
    
    cout<<"Converting first image to gray scale version\n\n";
    string firstBlackAndWhiteImage = folderName + currentDateAndTimeValue + "first_b_and_w_image.png";
    SImageIO::write_png_file(firstBlackAndWhiteImage.c_str(), image1, image1, image1);
   
    cout<<"Applying Gaussian filter to first input grayscale images\n\n";
    string firstGaussianImage = folderName + currentDateAndTimeValue + "first_gaussian_image.png";
    SDoublePlane gaussed1 = gaussian_filter(image1,1,3);
    SImageIO::write_png_file(firstGaussianImage.c_str(), gaussed1 , gaussed1 , gaussed1);
    
    cout<<"Applying Sobel filter in X and Y direction for first image\n\n";
    string firstGradientXImage = folderName + currentDateAndTimeValue +"first_gradient_X_image.png";
    SDoublePlane gradient_x1 = sobel_gradient_filter(gaussed1, true);
    SImageIO::write_png_file(firstGradientXImage.c_str(), gradient_x1, gradient_x1, gradient_x1);

    string firstGradientYImage = folderName + currentDateAndTimeValue +"first_gradient_Y_image.png";
    SDoublePlane gradient_y1 = sobel_gradient_filter(gaussed1, false);
    SImageIO::write_png_file(firstGradientYImage.c_str(), gradient_y1, gradient_y1, gradient_y1);

    cout<<"Applying Tomasi Corner detection algorithm for first image\n\n";
    find_corners_tomasi(image1,gradient_x1,gradient_y1, "1") ;

    if(toStitchImages == 1) {
    SDoublePlane image2= SImageIO::read_png_file(secondInputImageName);
    cout<<"Converting second image to gray scale version\n\n";
    string secondBlackAndWhiteImage = folderName + currentDateAndTimeValue +"second_b_and_w_image.png";
    SImageIO::write_png_file(secondBlackAndWhiteImage.c_str(), image2, image2, image2);
    
    cout<<"Applying Gaussian filter to second input grayscale images\n\n";
    string secondGaussianImage = folderName + currentDateAndTimeValue +"second_gaussian_image.png";
    SDoublePlane gaussed2 = gaussian_filter(image2,1,3);
    SImageIO::write_png_file(secondGaussianImage.c_str(), gaussed2 , gaussed2 , gaussed2);

  
    cout<<"Applying Sobel filter in X and Y direction for second image\n\n";
    string secondGradientXImage = folderName + currentDateAndTimeValue + "second_gradient_X_image.png";
    string secondGradientYImage = folderName + currentDateAndTimeValue + "second_gradient_Y_image.png";
    
    SDoublePlane gradient_x2 = sobel_gradient_filter(gaussed2, true);
    SImageIO::write_png_file(secondGradientXImage.c_str(), gradient_x2, gradient_x2, gradient_x2);

    SDoublePlane gradient_y2 = sobel_gradient_filter(gaussed2, false);
    SImageIO::write_png_file(secondGradientYImage.c_str(), gradient_y2, gradient_y2, gradient_y2);

    cout<<"Applying Tomasi Corner detection algorithm for second image\n\n";
    find_corners_tomasi(image2,gradient_x2,gradient_y2, "2") ;
    
    cout<<"Stitching two images together\n\n";
    stitch(image1, image2,gradient_x1,gradient_y1,gradient_x2,gradient_y2,output_filename,image1_filename,image2_filename);
    }
    else {
        cout<<"Applying harris corner detection to image\n\n";
        find_corners_harris(image1,gradient_x1,gradient_y1,1,3,image1_filename,"1");
        cout<<"Successfully applied harris corner detector to image\n\n";
    }
    
    //Stop the clock after measuring specific amount of time
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout<<"Program executed in "<<elapsed_secs<<" Seconds \n\n";
  
  return 0;
}
