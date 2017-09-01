from dolfin import *

import matplotlib.pyplot as plt
import matplotlib.tri as trii
import numpy as np

## TODO LIST: 
# (1) Add options to control view points - 3D options
# (2) 3D plot to deal with cell function
# (3) Add animation to quiver. 

#------------------------2D plotting--------------------#
# Inherit from Chris
# Add extra control option to the graph - Only on 2D work
# but also very slow
# import mpld3
# mpld3.enable_notebook()

# Generate triangles from mesh points
def mesh2triang(mesh):
	xy = mesh.coordinates()
	return trii.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

# 2D plot of function
def mplot_cellfunction(cellfn):
	C = cellfn.array()
	tri = mesh2triang(cellfn.mesh())
	return plt.tripcolor(tri, facecolors=C)


def mplot_function(f):
	mesh = f.function_space().mesh()
	if (mesh.geometry().dim() != 2):
		raise AttributeError('Mesh must be 2D')
	# DG0 cellwise function
	if f.vector().size() == mesh.num_cells():
		C = f.vector().array()
		return plt.tripcolor(mesh2triang(mesh), C)
	# Scalar function, interpolated to vertices
	elif f.value_rank() == 0:
		C = f.compute_vertex_values(mesh)
		return plt.tripcolor(mesh2triang(mesh), C, shading='gouraud')
	# Vector function, interpolated to vertices
	elif f.value_rank() == 1:
		w0 = f.compute_vertex_values(mesh)
		if (len(w0) != 2*mesh.num_vertices()):
			raise AttributeError('Vector field must be 2D')
		X = mesh.coordinates()[:, 0]
		Y = mesh.coordinates()[:, 1]
		U = w0[:mesh.num_vertices()]
		V = w0[mesh.num_vertices():]
		return plt.quiver(X,Y,U,V)

# Plot a generic dolfin object (if supported)
def matplot(obj):
	plt.gca().set_aspect('equal')
	if isinstance(obj, Function):
		return mplot_function(obj)
	elif isinstance(obj, CellFunctionSizet):
		return mplot_cellfunction(obj)
	elif isinstance(obj, CellFunctionDouble):
		return mplot_cellfunction(obj)
	elif isinstance(obj, CellFunctionInt):
		return mplot_cellfunction(obj)
	elif isinstance(obj, Mesh):
		if (obj.geometry().dim() != 2):
			raise AttributeError('Mesh must be 2D')
		return plt.triplot(mesh2triang(obj), color='#808080')
	raise AttributeError('Failed to plot %s'%type(obj))

## Plot in 3D
from mpl_toolkits.mplot3d import Axes3D

#---------------------------3D plotting------------------------#
# Usage: matplot3d(obj, elev, azim)
# elev stores the elevation angle in the z plane.
# azim stores the azimuth angle in the x,y plane.
## Plot a scalar function in 3D
def mplot_function3d(f, elev=None, azim=None):
	mesh = f.function_space().mesh()
	X = mesh.coordinates()[:, 0]
	Y = mesh.coordinates()[:, 1]
	# DG0 cellwise function
	if f.vector().size() == mesh.num_cells():
		C = f.vector().array()
		ax = plt.gca(projection='3d')
		ax.view_init(elev=elev, azim=azim)
		ax.plot_trisurf(X, Y, C, cmap=plt.cm.jet, linewidth=0.0)
		return ax
	# Scalar function, interpolated theo vertices
	elif f.value_rank() == 0:
		C = f.compute_vertex_values(mesh)
		ax = plt.gca(projection='3d')
		ax.view_init(elev=elev, azim=azim)
		ax.plot_trisurf(X, Y, C, cmap=plt.cm.jet, linewidth=0.0)
		return ax
	# Vector function, interpolated to vertices
	# Not really working...
	elif f.value_rank() == 1:
		w0 = f.compute_vertex_values(mesh)
		if (len(w0) != 2*mesh.num_vertices()):
			raise AttributeError('Vector field must be 2D')
		ax = plt.gca(projection='3d')
		X = mesh.coordinates()[:, 0]
		Y = mesh.coordinates()[:, 1]
		C = np.zeros(len(X)) # This is to ensure plotting in one plane
		U = w0[:mesh.num_vertices()]
		V = w0[mesh.num_vertices():]
		ax.quiver(X,Y,C,U,V,C, length=0.05, arrow_length_ratio=0.05, cmap=plt.cm.jet, linewidth=0.5)
		ax.view_init(elev=elev, azim=azim)
		return ax

def matplot3d(obj, elev=None, azim=None):
	plt.gca().set_aspect('equal')
	if isinstance(obj, Function):
		return mplot_function3d(obj, elev, azim)
	raise AttributeError('Not yet implemented type of %s'%type(obj))