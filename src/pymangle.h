/*
Time-stamp: <pymangle.cc on Tuesday, 29 January, 2013 at 17:58:05 MST (philip)>
*/

#include <cstdio>
#include <cstdlib>
#include <cassert>

/*
  Three functions to take flat multidimensional arrays (arrays of one index) and
  their sets of dimensions from Python,
  returning a C pointer-to-pointer-to... object which can be referred to 
  as a multidimensional object, e.g. x[i][j][k]

  The only memory allocated is that holding the array pointers.
 */

// two-dimensional
template <class T>
T** pymangle(int nx, int ny, T* data_from_python);

template <class T>
void freepymangle(T** m);

// three-dimensional
template <class T>
T*** pymangle(int nx, int ny, int nz, T* data_from_python);

template <class T>
void freepymangle(T*** m);

// four-dimensional

template <class T>
T**** pymangle(int nx, int ny, int nz, int nq, T* data_from_python);

template <class T>
void freepymangle(T**** m);
