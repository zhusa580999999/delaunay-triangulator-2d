****************************************************************
      STREAMING COMPUTATION OF 3D DELAUNAY TRIANGULATIONS
****************************************************************

This document describes the software for streaming *3D* Delaunay
triangulation. There is another README_2D.txt that describes how
to use this software for Delaunay triangulation in 2D.

The concept of "spatial finalization" is implemented in form of
an adaptive octree whose cells are declared "finalized" if no
downstream points fall into them.

Only the x, y, and z coordinates of the input points are considered.
Other scalar values are currently ignored.

******************************************************************
          A QUICK GUIDE FOR GETTING YOUR HANDS DIRTY
******************************************************************

first run ...

spfinalize -i spx_z.spb -o spx_z_finalized.spa -level 2

or better (produces more IO-efficient binary *.spb):

spfinalize -i spx_z.spb -o spx_z_finalized.spb -level 2

... then run ...

sp_viewer -i spx_z_finalized.spb -steps 20

... and press <p> and then run ...

sp_viewer -i spx_z_finalized.spb -del3

... and press <p> or for producing an actual output mesh run ...

spdelaunay3d -i spx_z_finalized.spb -o spx_z.svb

... and to look at the output run ...

sv_viewer -i spx_z.svb -every 2000

... and press <p> or try some PIPELINING again ...

spfinalize -i spx_z.spb -ospb -level 1 | spdelaunay3d -ispb -osvb | sv_viewer -isvb -every 1000

... and press <p> or for a bigger model try ...

sp_viewer -i torso_t.spb -del3 -steps 30

... and press <p> or for a rather difficult model try ...

sp_viewer -i  fighter_s_4m2.spb -del3

.. and ...

spdelaunay3d -i fighter_s_4m2.spb -o fighter_s_4m2.svb

to create the mesh that you can look at with ...

sv_viewer -i fighter_s_4m2.svb

... and press <p>.

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

the option '-terrain' makes the finalizer create a 2D finalized point
cloud based on the x and y coordinate that are used by our streaming
2D delaunay triangulator.


              /*********************************************/
              /***             sp_viewer.exe             ***/
              /*********************************************/

usage:

sp_viewer -h

visualizes a streaming 3D point cloud and its spatial finalization octtree.

option '-del3' will visualize the state of streaming delaunay tetrahedralization

press 'p' to play/stop
press 's' to single step

toggle render mode with
 'u'  -  unfinalized space
 'a'  -  active tetrahedra
 'i'  -  infinite tetrahedra (only their finite triangle is drawn)

toggle interaction mode with <space>

change point size with +/-

move pointer above a tetrahedron and press 'm'. shows you circumspheres that
intersect the grid (if computed already).

press 'q' or <ESC> to quit

the memory requirements are proportional to the maximal number of terahedra
and vertices that are allocated concurrently (the actual number is output to
stdout in the end) and the maximal number of octree cells in the octree that
is used to maintain and query the finalized space.

if given a non-streaming point cloud nothing can get finalized and the memory
requirements are much higher (e.g. like a standard delaunay tetrahedralizator).


              /*********************************************/
              /***           spdelaunay3d.exe            ***/
              /*********************************************/

usage:
spdelaunay3d -h

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
sv_viewer -h
sv_viewer mesh.sva -every 1000
sv_viewer mesh.svb

allows to visualize the meshes you created. the '-every <ntets>' option specifies that
the viewer shows a visualization step every <ntets> tetrahedra. the default is 10,000.
