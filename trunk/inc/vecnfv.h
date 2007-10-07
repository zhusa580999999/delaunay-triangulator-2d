/*
===============================================================================

  FILE:  vecNfv.h
  
  CONTENTS:
  
    inlined functions for common vecNfv operations
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu  (with ideas stolen from Kenny Hoff)
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    28 June 2000 -- created for TET mesh property handling
  
===============================================================================
*/
#ifndef VECNFV_H
#define VECNFV_H

inline void VecZeroNfv(float v[], int n);

inline void VecCopyNfv(float v[], const float a[], int n);
inline void VecCopyNfv_swap_endian(float v[], const float a[], int n);

inline float VecMinNfv(const float v[], int n);
inline float VecMaxNfv(const float v[], int n);

inline void VecSelfScalarMultNfv(float v[], float s, int n);
inline void VecSelfScalarDivNfv(float v[], float s, int n);
inline void VecSelfNegateNfv(float v[], int n);

inline void VecSelfAddNfv(float v[], const float a[], int n);
inline void VecSelfSubtractNfv(float v[], const float a[], int n);

inline void VecAddNfv(float v[], const float a[], const float b[], int n);
inline void VecSubtractNfv(float v[], const float a[], const float b[], int n);
inline void VecAbsDiffNfv(float v[], const float a[], const float b[], int n);
inline void VecScalarDivNfv(float v[], const float a[], float s, int n);

inline void VecAddScalarMultNfv(float v[], const float a[], const float b[], float s, int n);
inline void VecSelfAddScalarMultNfv(float v[], float a[], float s, int n);

inline bool VecEqualNfv(const float a[], const float b[], int n);

inline void VecUpdateMinMaxNfv(float min[], float max[], const float v[], int n);
inline void VecUpdateMinNfv(float min[], const float v[], int n);
inline void VecUpdateMaxNfv(float max[], const float v[], int n);

inline void VecZeroNfv(float v[], int n)
{
  for (int i = 0; i < n; i++) v[i] = 0.0f;
}

inline void VecCopyNfv(float v[], const float a[], int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i];
}

inline void VecCopyNfv_swap_endian(float v[], const float a[], int n)
{
  for (int i = 0; i < n; i++)
  {
    int j = 4*i;
    ((char*)v)[j+0] = ((const char*)a)[j+3];
    ((char*)v)[j+1] = ((const char*)a)[j+2];
    ((char*)v)[j+2] = ((const char*)a)[j+1];
    ((char*)v)[j+3] = ((const char*)a)[j+0];
  }
}

inline float VecMinNfv(const float v[], int n)
{
  float min = v[0];
  for (int i = 1; i < n; i++) if (v[i] < min) min = v[i];
  return min;
}

inline float VecMaxNfv(const float v[], int n)
{
  float max = v[0];
  for (int i = 1; i < n; i++) if (v[i] > max) max = v[i];
  return max;
}

inline void VecSelfScalarMultNfv(float v[], float s, int n)
{
  for (int i = 0; i < n; i++) v[i] *= s;
}

inline void VecSelfScalarDivNfv(float v[], float s, int n)
{
  for (int i = 0; i < n; i++) v[i] /= s;
}

inline void VecSelfNegateNfv(float v[], int n)
{
  for (int i = 0; i < n; i++) v[i] = -v[i];
}

inline void VecSelfAddNfv(float v[], const float a[], int n)
{
  for (int i = 0; i < n; i++) v[i] += a[i];
}

inline void VecSelfSubtractNfv(float v[], const float a[], int n)
{
  for (int i = 0; i < n; i++) v[i] -= a[i];
}

inline void VecAddNfv(float v[], const float a[], const float b[], int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i] + b[i];
}

inline void VecSubtractNfv(float v[], const float a[], const float b[], int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i] - b[i];
}

inline void VecAbsDiffNfv(float v[], const float a[], const float b[], int n)
{
  for (int i = 0; i < n; i++)
  {
    v[i] = a[i] - b[i];
    if (v[i] < 0.0f) v[i] = -v[i];
  }
}

inline void VecScalarDivNfv(float v[], const float a[], float s, int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i]/s;
}

inline void VecAddScalarMultNfv(float v[], const float a[], const float b[], float s, int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i] + (s * b[i]);
}

inline void VecSelfAddScalarMultNfv(float v[], const float a[], float s, int n)
{
  for (int i = 0; i < n; i++) v[i] += (s * a[i]);
}

inline void VecUpdateMinMaxNfv(float min[], float max[], const float v[], int n)
{
  for (int i = 0; i < n; i++) if (v[i]<min[i]) min[i]=v[i]; else if (v[i]>max[i]) max[i]=v[i];
}

inline void VecUpdateMinNfv(float min[], const float v[], int n)
{
  for (int i = 0; i < n; i++) if (v[i]<min[i]) min[i]=v[i];
}

inline void VecUpdateMaxNfv(float max[], const float v[], int n)
{
  for (int i = 0; i < n; i++) if (v[i]>max[i]) max[i]=v[i];
}

inline bool VecEqualNfv(const float a[], const float b[], int n)
{
  for (int i = 0; i < n; i++) if (a[i]!=b[i]) return false;
  return true;
}

#endif
