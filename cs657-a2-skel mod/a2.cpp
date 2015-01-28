#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>

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
vector<Coordinate> find_corners(const SDoublePlane &input)
{
  // Implement a Harris corner detector
  // (This placeholder code just returns a set of random points.)
  
  vector<Coordinate> points;
  const int point_count = rand() % 20;
  for(int i = 0; i < point_count; i++) {
    Coordinate c;
    c.row = rand() % input.rows();
    c.col = rand() % input.cols();

    points.push_back(c);
  }
  return points;
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
void stitch(SDoublePlane &image1, SDoublePlane &image2, string output_filename)
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
}

int main(int argc, char *argv[])
{
  if(!(argc == 4)) {
    cerr << "usage: " << argv[0] << " output_file image_file1 image_file2" << endl;
    return 1;
  }

  int image_count = argc - 3;
  string output_filename = argv[1];
  
  // read in query image and compute descriptors
  string image1_filename(argv[2]);
  string image2_filename(argv[3]);
  SDoublePlane image1= SImageIO::read_png_file(image1_filename.c_str());
  SDoublePlane image2= SImageIO::read_png_file(image2_filename.c_str());
  
  // image matching
  stitch(image1, image2, output_filename);
  
  return 0;
}
