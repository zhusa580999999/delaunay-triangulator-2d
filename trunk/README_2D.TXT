****************************************************************
      STREAMING COMPUTATION OF 2D DELAUNAY TRIANGULATIONS
****************************************************************

This document describes the software for streaming *2D* Delaunay
triangulation. There is another README_3D.txt that describes how
to use this software for Delaunay triangulation in 3D.

The concept of "spatial finalization" is implemented in form of
an adaptive quadtree whose cells are declared "finalized" if no
downstream points fall into them.

Only the x and y coordinate of the input points are considered so
that spatial finalization and Delaunay triangulation happen in the
x/y plane. The z coordinate is assumed the height during visualization.

!!!! IMPORTANT !!!!!

the option '-terrain' makes the finalizer create a 2D finalized point
cloud based on the x and y coordinate. without this option it creates
a 3D finalized point cloud using an octree. avoid that for terrains.

******************************************************************
          A QUICK GUIDE FOR GETTING YOUR HANDS DIRTY
******************************************************************

just copy and paste the following into the command line
 
spfinalize -i lidar_ncfloodmap.txt.gz -ispa -terrain -o lidar_ncfloodmap_finalized.spb

sp_viewer -i lidar_ncfloodmap_finalized.spb 

press <p>

sp_viewer -i lidar_ncfloodmap_finalized.spb -del2

press <p>

spdelaunay2d -i lidar_ncfloodmap_finalized.spb -o lidar_ncfloodmap.smb

sm_viewer -i lidar_ncfloodmap.smb -every 1000

press <p>

for a "larger" datasets try:

spfinalize -i noaa.spb -terrain -empty 3 -level 5 -o noaa_finalized.spb

sp_viewer -i noaa_finalized.spb

press <p>

sp_viewer -i noaa_finalized.spb -del2

press <p>

and without 'sprinkle'

spfinalize -i noaa.spb -terrain -empty 3 -level 5 -o noaa_finalized.spb -nosprinkle

sp_viewer -i noaa_finalized.spb -del2

press <p>

spdelaunay2d -i noaa_finalized.spb -o noaa_finalized.smb

sm_viewer -i noaa_finalized.smb -every 10000

press <p>

or PIPELINED

spfinalize -i noaa.spb -terrain -empty 3 -level 5 -ospb | spdelaunay2d -ispb -o noaa_finalized.smb

and a longer pipe

spfinalize -i noaa.spb -terrain -empty 3 -level 5 -ospb | spdelaunay2d -ispb -osmb | sm_viewer -ismb -every 10000

and now go and try your own data ...

              /*********************************************/
              /***            spfinalize.exe             ***/
              /*********************************************/

usage:

spfinalize -h

takes a point cloud as input and spatially finalizes it in two or
three passes.

we can also read *.obj files, *.ply files and *.sm[abc] files (optionally
gzipped) but use only their vertex data, which is treated as a 2D input
with z being the height. we can write *.spa (Streaming Point Ascii) and
*.spb (Streaming Point Binary) files. The binary representaion is much
more IO-efficient!!! 

              /*********************************************/
              /***             sp_viewer.exe             ***/
              /*********************************************/

usage:

sp_viewer -h

visualizes a streaming 2D point cloud and its spatial finalization quadtree.

option '-flatten' will set the z coordinate of all read point to zero.
option '-del2' will visualize the state of streaming delaunay triangulation

press 'p' to play/stop
press 's' to step

toggle render mode with
 'u'  -  unfinalized space
 'a'  -  active triangles
 'i'  -  infinite triangles

toggle interaction mode with <space>

change point size with -/=

change line width with _/+

move pointer above a triangle and press 'm'. shows you circumcircles that
intersect the grid (if computed already).

press 'q' or <ESC> to quit

the memory requirements are proportional to the maximal number of triangles
and vertices that are allocated concurrently (the actual number is output to
stderr in the end) and the maximal number of quadtree cells in the quadtree
that is used to maintain and query the finalized space.

if given a non-streaming point cloud nothing can get finalized and the memory
requirements are much higher (e.g. more than a standard delaunay triangulator).

              /*********************************************/
              /***           spdelaunay2d.exe            ***/
              /*********************************************/

usage:
spdelaunay2d -h

computes a delaunay triangulation of the input point cloud and outputs a streaming
mesh. for this to be efficient, the input cloud should be spatially finalized. if
the point cloud is not spatially finalized then we can not operate in a streaming
manner and incrementally build up the entire delaunay mesh in memory.

currently the streaming output mesh is not perfect, because not all vertices are
getting finalized (all those vertices that are incident to an infinite triangle).

              /*********************************************/
              /***             sm_viewer.exe             ***/
              /*********************************************/

usage:
sm_viewer -h

for visualizing the streaming output mesh (and the order in which it waqs created)
