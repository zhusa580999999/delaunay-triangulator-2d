/*
===============================================================================

  FILE:  vecNiv.h
  
  CONTENTS:
  
    inlined functions for common vecNiv operations
  
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
#ifndef VECNIV_H
#define VECNIV_H

inline void VecZeroNiv(int v[], int n);

inline void VecCopyNiv(int v[], const int a[], int n);
inline void VecCopyNiv_swap_endian(int v[], const int a[], int n);

inline int VecMinNiv(const int v[], int n);
inline int VecMaxNiv(const int v[], int n);

inline void VecSelfScalarMultNiv(int v[], int s, int n);
inline void VecSelfScalarDivNiv(int v[], int s, int n);
inline void VecSelfNegateNiv(int v[], int n);

inline void VecSelfAddNiv(int v[], const int a[], int n);
inline void VecSelfSubtractNiv(int v[], const int a[], int n);

inline void VecAddNiv(int v[], const int a[], const int b[], int n);
inline void VecSubtractNiv(int v[], const int a[], const int b[], int n);
inline void VecAbsDiffNiv(int v[], const int a[], const int b[], int n);
inline void VecScalarDivNiv(int v[], const int a[], int s, int n);

inline void VecAddScalarMultNiv(int v[], const int a[], const int b[], int s, int n);
inline void VecSelfAddScalarMultNiv(int v[], int a[], int s, int n);

inline bool VecEqualNiv(const int a[], const int b[], int n);

inline void VecUpdateMinMaxNiv(int min[], int max[], const int v[], int n);
inline void VecUpdateMinNiv(int min[], const int v[], int n);
inline void VecUpdateMaxNiv(int max[], const int v[], int n);

inline void VecZeroNiv(int v[], int n)
{
  for (int i = 0; i < n; i++) v[i] = 0;
}

inline void VecCopyNiv(int v[], const int a[], int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i];
}

inline void VecCopyNiv_swap_endian(int v[], const int a[], int n)
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

inline int VecMinNiv(const int v[], int n)
{
  int min = v[0];
  for (int i = 1; i < n; i++) if (v[i] < min) min = v[i];
  return min;
}

inline int VecMaxNiv(const int v[], int n)
{
  int max = v[0];
  for (int i = 1; i < n; i++) if (v[i] > max) max = v[i];
  return max;
}

inline void VecSelfScalarMultNiv(int v[], int s, int n)
{
  for (int i = 0; i < n; i++) v[i] *= s;
}

inline void VecSelfScalarDivNiv(int v[], int s, int n)
{
  for (int i = 0; i < n; i++) v[i] /= s;
}

inline void VecSelfNegateNiv(int v[], int n)
{
  for (int i = 0; i < n; i++) v[i] = -v[i];
}

inline void VecSelfAddNiv(int v[], const int a[], int n)
{
  for (int i = 0; i < n; i++) v[i] += a[i];
}

inline void VecSelfSubtractNiv(int v[], const int a[], int n)
{
  for (int i = 0; i < n; i++) v[i] -= a[i];
}

inline void VecAddNiv(int v[], const int a[], const int b[], int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i] + b[i];
}

inline void VecSubtractNiv(int v[], const int a[], const int b[], int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i] - b[i];
}

inline void VecAbsDiffNiv(int v[], const int a[], const int b[], int n)
{
  for (int i = 0; i < n; i++)
  {
    v[i] = a[i] - b[i];
    if (v[i] < 0.0f) v[i] = -v[i];
  }
}

inline void VecScalarDivNiv(int v[], const int a[], int s, int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i]/s;
}

inline void VecAddScalarMultNiv(int v[], const int a[], const int b[], int s, int n)
{
  for (int i = 0; i < n; i++) v[i] = a[i] + (s * b[i]);
}

inline void VecSelfAddScalarMultNiv(int v[], const int a[], int s, int n)
{
  for (int i = 0; i < n; i++) v[i] += (s * a[i]);
}

inline void VecUpdateMinMaxNiv(int min[], int max[], const int v[], int n)
{
  for (int i = 0; i < n; i++) if (v[i]<min[i]) min[i]=v[i]; else if (v[i]>max[i]) max[i]=v[i];
}

inline void VecUpdateMinNiv(int min[], const int v[], int n)
{
  for (int i = 0; i < n; i++) if (v[i]<min[i]) min[i]=v[i];
}

inline void VecUpdateMaxNiv(int max[], const int v[], int n)
{
  for (int i = 0; i < n; i++) if (v[i]>max[i]) max[i]=v[i];
}

inline bool VecEqualNiv(const int a[], const int b[], int n)
{
  for (int i = 0; i < n; i++) if (a[i]!=b[i]) return false;
  return true;
}

#endif
