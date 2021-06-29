"""
Copyright (c) 2020, Antony Orton.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of Antony Orton nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import sqlite3
import pandas as pd
import numpy as np
import os
import shapely.geometry as shp
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import re
import scipy.interpolate as sinterp
import scipy.ndimage as ndimg
import matplotlib.path as matpath
from itertools import islice
#import triangle
import matplotlib.tri as trimat
from pyproj import Proj
from pyproj import Transformer
import descartes as dc
import shapefile as pyshp
import scipy.spatial


def read_csv_input_files(file_directory = None):
	""" Read in all csv files in file_directory (or the current directory if = None) with 'hole' or 'geology' in their name
		and create return two dataframes [dfhole, dfgeology]
	
		'hole' csv: must contain columns ['borehole', 'x', 'y', 'top_rl', 'EOH_depth','dip','dip_dirn']
		
		'geology' csv: must contain columns ['borehole','fromDepth','material']
		
		returns two dataframes with common ID column 'borehole'
			dfhole: coordinate information
			dfgeology: geology information
	"""
	
	if file_directory != None:
		dir1 = os.getcwd()
		os.chdir(file_directory)
	
	
	#get names of all csv files in directory containing names 'hole' or 'geology'
	list1 = os.listdir()
	data_list=[l1 for l1 in list1 if l1[-3::]=='csv']
	data_hole = [l1 for l1 in data_list if ('hole' in l1 or 'Hole' in l1)]
	data_geology = [l1 for l1 in data_list if ('geology' in l1 or 'Geology' in l1)]


	#Read hole csv data (coords) to pandas dataframe
	dfhole = pd.DataFrame()
	for i in range(len(data_hole)):
		dfi = pd.read_csv(data_hole[i])
		if list(dfi)==['borehole', 'x', 'y', 'top_rl', 'EOH_depth','dip','dip_direction']:
			dfhole = pd.concat([dfhole,pd.read_csv(data_hole[i])])
			print('Hole data successfully loaded for file: '+data_hole[i])
		else:
			raise ValueError('Error in file '+data_hole[i]+' - columns are not labelled exactly as: [borehole, x, y, top_rl, EOH_depth, dip, dip_direction]')
			print('Data not loaded for file: '+data_hole[i])
	print('[WARNING]: In HOLES data, the following duplicated borehole entries have been dropped ...')
	print(dfhole[dfhole.duplicated(subset='borehole')])
	print('  ')
	dfhole.drop_duplicates(subset='borehole',inplace=True)
	dfhole.reset_index(inplace=True,drop=True)


	#Read geology csv data (coords) to pandas dataframe
	dfgeology = pd.DataFrame()
	for i in range(len(data_geology)):
		dfi = pd.read_csv(data_geology[i])
		if list(dfi)==['borehole','fromDepth','material']:
			dfgeology = pd.concat([dfgeology,pd.read_csv(data_geology[i])])
			print('Geology data successfully loaded for file: '+data_geology[i])
		else:
			raise ValueError('Error in file '+data_geology[i]+' - columns are not labelled exactly as: [borehole,fromDepth,material]')
			print('Data not loaded for file: '+data_geology[i])
	dfgeology.reset_index(inplace=True,drop=True)
	
	if file_directory != None:
		os.chdir(dir1)
	
	return [dfhole,dfgeology]
	
def shapefile_to_shapely(input_shapefile):
	
	"""reads a shapfile and converts to a list of shapely objects"""
	"""NOTE: re-written on 29-8-2020 to remove dependency on Fiona/Gdal"""
	
	
	sf = pyshp.Reader(input_shapefile)
	ftype = sf.shapeTypeName
	shapes = sf.shapes()
	
	list1=[]
	for i in range(len(sf)):
		coords =  shapes[i].points
		if ftype == 'POINT':
			list1.append(shp.Point(coords[0]))
		#NOTE: multipoint not done yet	
		#if ftype == 'MULTIPOINT':
		#	list1.append(shp.Point(coords[0]))
		elif ftype == 'POLYLINE' or ftype == 'POLYLINEZ':
			list1.append(shp.LineString(coords))
		elif ftype == 'POLYGON':
			list1.append(shp.Polygon(coords))
	
	print('shapfile successfully converted to list of shapely objects')
	print('shapefile records:')
	[print(sf.records()[i]) for i in range(len(sf))]
	
	print(list1)
	return list1

def shapely_to_shapefile(input,save_filename):

	"""Take an input list of shapely objects (LineString, Polygon, Point) and create and save as a shapefile
		NOTE: re-written on 29-8-2020 to remove dependency on Fiona/Gdal"""


	if '.shp' in save_filename:
		save_filename = save_filename[0:-4]

	if not type(input)==list:
		inputdata = [input]
	else:	
		inputdata = input[:]

	if not (type(inputdata[0])==shp.LineString or type(inputdata[0])==shp.MultiPoint or type(inputdata[0])==shp.Polygon or type(inputdata[0])==shp.Point):
		print('Error: Only input shapely Polygons, Linestrings, Points or MultiPoint are supported at present')
		return

	
	coords = []
	for i in range(len(inputdata)):
		coords.append(shp.mapping(inputdata[i])['coordinates'])
	
	w = pyshp.Writer(save_filename)
	w.field('name', 'C')
	
	#NOTE: may still be issues here with w.line, w.poly (data may need to be supplied as list)
	if type(inputdata[0])==shp.LineString:
		for i in range(len(inputdata)):
			print(i)
			w.record('Line'+str(i))
			w.line([coords[i]])
	elif type(inputdata[0])==shp.Polygon:
		for i in range(len(inputdata)):
			w.record('Polygon'+str(i))
			w.poly(coords[i])
	elif type(inputdata[0])==shp.Point:
		for i in range(len(inputdata)):
			w.record('Point'+str(i))
			w.point(coords[i][0],coords[i][1])			
	elif type(inputdata[0])==shp.MultiPoint:
		for i in range(len(inputdata)):
			w.record('Multipoint'+str(i))
			w.multipoint([coords[i]])		
	w.close()
	
	return
		
def get_nearby_extents(line,max_offset = 60):
	"""Get extents in relation to input section line
	   Input: line must be a shapely linestring 
	   
	   returns: extents = (xmin,ymin,xmax,ymax)
	   """
	if not type(line)==shp.LineString:
		print('Error: line must be a shapely linestring')
	 
	#centreL = np.array(line.centroid)
	#lengthL = line.length
	#xmin = centreL[0] - lengthL/2 - max_offset
	#xmax = centreL[0] + lengthL/2 + max_offset
	#ymin = centreL[1] - lengthL/2 - max_offset
	#ymax = centreL[1] + lengthL/2 + max_offset
	
	line1 = np.array(line)
	
	xmin = np.min(line1[:,0]) - max_offset
	xmax = np.max(line1[:,0]) + max_offset
	ymin = np.min(line1[:,1]) - max_offset
	ymax = np.max(line1[:,1]) + max_offset
	
	return (xmin,ymin,xmax,ymax)
	
def point_on_left_or_right(line,point):
	
	"""determine which side of 'line' a 'point' lies on. Both are shapely objects."""
	
	pt = point
	
	pt0 = line.interpolate(line.project(pt))
	pt1 = line.interpolate(line.project(pt)+1e-5)
	
	pt = np.array(pt)[0:2]
	pt0 = np.array(pt0)[0:2]  #ignore z coords
	pt1 = np.array(pt1)[0:2]  #ignore z coords
	
	cross1 = np.cross(pt - pt0, pt1-pt0)
	
	return np.sign(cross1)	
	
def get_points_near_ctrl_line(line,dfpoints,max_offset):
	
	"""
		line: a shapely linestring
		dfpoints: a dataframe with x,y columns to be checked
		max_offset: distance within control line to accept points
	
		returns: dataframe subset of dfpoints within max_offset of control line with 'chain' and 'offset' columns added
	"""
	
	if not (('x' in list(dfpoints)) and ('y' in list(dfpoints))):
		print('Error: Input dataframe does not contain both x and y columns')
		print('Error in function: get_points_near_ctrl_line()')
		return
	
	envelope1 = line.buffer(max_offset)  #get envelope
	list_points = [shp.Point(coords) for coords in np.array(dfpoints[['x','y']])] #create points
	list1 = [envelope1.contains(p1) for p1 in list_points] #check T/F if within max_offset of ctrl line
	
	new_holes = dfpoints[list1].copy()
	new_holes['offset'] = get_offset_from_ctrl_line(line,new_holes)
	new_holes['chain'] = get_chainage(line,new_holes)
	
	return new_holes
	
def get_offset_from_ctrl_line(line,dfpoints):

	"""
		line: a shapely linestring
		dfpoints: a dataframe with x,y columns
		
		returns: signed distance to line for all points in dfpoints
	"""
	
	if not (('x' in list(dfpoints)) and ('y' in list(dfpoints))):
		print('Error: Input dataframe does not contain both x and y columns')
		print('Error in function: get_offset_from_ctrl_line()')
		return

	list_points = [shp.Point(coords) for coords in np.array(dfpoints[['x','y']])] #create points
	list1 = [line.distance(p1)*point_on_left_or_right(line,p1) for p1 in list_points]
	
	return list1

def get_chainage(line,dfpoints,start_chainage=0,flag_points = True):

	"""
		line: a shapely linestring
		dfpoints: a dataframe with x,y columns to be checked or ARRAY with first two cols being xy coords
		start_chainage: starting chainage of line (eg may start at 200, 
						whereas in shapely it is always zero)
		
		
		returns: chainage along line for all points in dfpoints
	"""
	if type(dfpoints) == np.ndarray:
		dfpoints = pd.DataFrame(dfpoints[:,0:2],columns = ['x','y'])
		
	
	len1 = line.length

	start_point = line.coords[0]
	end_point = line.coords[-1]
	
	list_points = [shp.Point(coords) for coords in np.array(dfpoints[['x','y']])] #create points
	list1 = [line.project(p1) for p1 in list_points] #get chainage
	list2 = [l1 if (l1>0)&(l1<len1) else -1 for l1 in list1] #FLAG points outside the chainage range
	
	
	#include start and end point
	for i in range(len(list2)):
		if list2[i]==-1:
			if list_points[i].coords[:][0]==start_point[0:2]:
				list2[i]=0
			elif list_points[i].coords[:][0]==end_point[0:2]:
				list2[i]=len1
	
	
	return list2
	
def get_borehole_geology(borehole,dfgeology,dfholes):

	""" borehole: string - borehole name
		dfgeology: geology dataframe must have 'borehole' and 'fromDepth' columns
		dfholes: main hole coordinates dataframe to read off end of hole data, must have 'EOH_depth' column
		
		returns: dataframe with dfgeo['borehole']=borehole, sorted by depth, end_of_hole added
	"""
	
	bh1 = dfgeology[dfgeology['borehole'] == borehole].copy()
	bh1.sort_values('fromDepth',inplace=True)
	bh1.drop_duplicates(subset='fromDepth',inplace=True)

	s1 = {'borehole':borehole,'fromDepth':float(dfholes[dfholes['borehole']==borehole]['EOH_depth']),'material':'end_of_hole'}
	bh1 = bh1.append(s1,ignore_index=True)   # append EOH to end
	
	if np.abs(float(bh1.iloc[0]['fromDepth']))>0.001:   # add in Fill if zero depth unit not recorded
		s1 = {'borehole':borehole,'fromDepth':0.0,'material':'Fill'}  
		bh1 = pd.DataFrame(s1,index=[0]).append(bh1,ignore_index=True)
	
	bh1.reset_index(inplace=True,drop=True)
	
	return bh1
	
def get_geology_coordinates(bore_geology,dfholes,line,start_chainage = 0):
	""" add x,y,z and chainage columns to bore_geology dataframe
		
		bore_geology: intended to be for a single hole and contains 'borehole' and 'fromDepth' columns
		dfholes: is main holes coordinates dataframe with 'x','y','dip','dip_directon' columns
		line: is shapely linestring (section line for chainage)
		start_chainage: starting chainage on section line
		
		returns: new dataframe
	"""
	
	bore_new = bore_geology.copy()
	
	name1 = list(bore_geology['borehole'].drop_duplicates())
	if len(name1)>1:
		raise ValueError('Error: bore_geology input to get_geology_coordinates function contains more than one borehole')
		return
	else:
		name1 = name1[0]
   
	#get hole coordinates
	x0 = float(dfholes[dfholes['borehole']==name1]['x'])
	y0 = float(dfholes[dfholes['borehole']==name1]['y'])
	z0 = float(dfholes[dfholes['borehole']==name1]['top_rl'])
	dip = float(np.radians(dfholes[dfholes['borehole']==name1]['dip']))
	dip_dir = float(np.radians(dfholes[dfholes['borehole']==name1]['dip_direction']))
	
	X=[]
	Y=[]
	Z=[]
	
	#find coordinates for each top of stratum
	if np.abs(dip+np.radians(90))>1e-5:   #check for dip not -90 degrees
		for i in range(len(bore_geology)):
			depth1 = bore_geology.iloc[i]['fromDepth']
			X.append(x0 + depth1 * np.cos(dip) * np.sin(dip_dir))
			Y.append(y0 + depth1 * np.cos(dip) * np.cos(dip_dir))
			Z.append(z0 + depth1 * np.sin(dip))
	else:
		X = [x0]*len(bore_geology)
		Y = [y0]*len(bore_geology)
		for i in range(len(bore_geology)):
			depth1 = bore_geology.iloc[i]['fromDepth']
			Z.append(z0 - depth1)
   
	bore_new['x'] = X
	bore_new['y'] = Y
	bore_new['z'] = Z
	bore_new['chain'] = np.array(get_chainage(line,bore_new)) + start_chainage
	
	return bore_new
  
def get_colours(dfgeology,materials_list = '',colourmap = 'hsv'):
	""" provides a dictionary providing an rgba colour tuple for each material in dfgeology dataframe
	Good to use dfgeology containing all materials in the whole project to keep consistent colours 
	materials_list: List, optional, specify list of all materials in order, otherwise order will be chosen at random
	
	returns: dictionary
	"""
	
	if type(materials_list) == list:
		all_materials = materials_list
	else:
		all_materials = np.sort(list(dfgeology['material'].drop_duplicates()))  

	
	#make sure Fill is in list of materials, as used a material above the highest specified depth 
	all_materials = [item if item != 'fill' else 'Fill' for item in all_materials]  
	if 'Fill' not in all_materials:
		all_materials = np.hstack((['Fill'],all_materials))	  

	
	scalarMap = cmx.ScalarMappable(cmap=colourmap)
	all_colours = scalarMap.to_rgba(np.arange(len(all_materials)))
	
	all_colours = [tuple(colsx) for colsx in all_colours]
	
	cols_dict = dict(zip(all_materials, all_colours))
	
	return cols_dict
	
#PLOT BOREHOLES -  cycle through geology, each borehole and plot desired section
def plot_boreholes_on_section(dfholes, dfgeology, line, title = 'Geotechnical long section', materials_list = '', colourmap = 'RdYlGn', max_offset = 100, start_chainage = 0):

	""" plot boreholes on section 
		dfholes: dataframe with coordinate info
		dfgeology: dataframe with geology info
		line: Linstring control line
		materials_list: List, optional for colourmap, specify list of all geology materials in order, otherwise order will be chosen at random
		max_offset: maximum borehole offset from section line to include
		start_chainage: starting chainage on line
	"""
	
	print('plotting boreholes on section ...')
	#extract nearby holes (for plan plotting)
	holes1 = get_points_near_ctrl_line(line,dfholes,max_offset)

	#get holes for section: ie where chainage is >0 and less than len(line)
	sholes = holes1[holes1['chain']!=-1].copy()
	sholes.reset_index(inplace=True,drop=True)
	sholes['chain'] += start_chainage
	
	#######Section plot
	
	#get colours
	colour_dict = get_colours(dfgeology,colourmap = colourmap,materials_list = materials_list)
	#print(colour_dict)
	
	#get text height in map x coordinate units
	xtextheight = 6.6e-3*line.length
	
	#Plotting of borholes
	for i in range(len(sholes)):
		bh1 = get_borehole_geology(sholes.iloc[i]['borehole'],dfgeology,sholes)  #extract geology
		bh1 = get_geology_coordinates(bh1,sholes,line,start_chainage = start_chainage)	   #add coordinates
		for j in range(len(bh1)-1):
			#plot lines
			plt.plot(bh1.iloc[j:j+2]['chain'],bh1.iloc[j:j+2]['z'],color = colour_dict[str(bh1.iloc[j]['material'])])
			#plot dots
			plt.plot(bh1.iloc[j]['chain'],bh1.iloc[j]['z'],'o',color = colour_dict[str(bh1.iloc[j]['material'])])
			#plot text
			plt.text(0.5*(bh1.iloc[j]['chain']+bh1.iloc[j+1]['chain']),0.5*(bh1.iloc[j]['z']+bh1.iloc[j+1]['z']),str(bh1.iloc[j]['material']))
			
		plt.text(sholes.iloc[i]['chain'],sholes.iloc[i]['top_rl']+1,str(sholes.iloc[i]['borehole'])+' (o/s '+str(sholes.iloc[i]['offset'])[0:5]+'m)',rotation = 45,rotation_mode='anchor')
		plt.text(sholes.iloc[i]['chain']-3,bh1.iloc[-1]['z']-0.5,str(sholes.iloc[i]['EOH_depth'])[0:5]+'m')
	
	plt.xlabel('Chainage (m)')
	plt.ylabel('Elevation (m AHD)')
	plt.title(title)
	plt.grid('on')
	
#PLOT GROUND SURFACE
def plot_groundsurface_on_section (line,xyzgrid,point_space = 10,colour = 'r',start_chainage = 0,smoothing_window = 5,height_datum = 'AHD'):
	
	"""
	Plot elevation vs chainage along section line.
	
	line (LineString)
	point_spacing (float/Int): distance along section line to get surface elevation from xyzgrid
	xyzgrid: DataFrame with x,y,z, columns
	start_chainage int/float): chainage at start of line
	DEM_directory: string, directory location containing all DEM files
	set DEM_directory = '' to search current directory
	WARNING: if copying directory from windows must put as raw string eg DEM_directory = r'C:\..somelocation\..Desktop\..python_codes .... note the r in front of the string
	
	smoothing_window = (int) window size for pandas smoothing algorithm. Set as 0 for no smoothing
	
	Returns: dataframe with columns [chain,z]
	"""
	
	print('extracting ground surface from dem ...')
	
	#extract nearby points to save cost of griddata function in get_zvals_fom_grid
	extents = get_nearby_extents(line,max_offset = 200)
	xyzgridnew = crop_xyz_grid_by_extents(xyzgrid,extents)
	
	pts, chain = get_points_at_equal_chainages(line,point_space = point_space)
	chain = chain + start_chainage
	zvals = get_zvals_from_grid(pts,xyzgridnew)
	
	#get a moving average to smooth noise in the data
	if smoothing_window > 0:
		mov = pd.DataFrame(zvals).rolling(smoothing_window,center = True,min_periods = 1).mean()
	else:
		mov = zvals
	
	plt.plot(chain,zvals,color = 'k',linewidth = 0.2)
	plt.plot(chain,mov,color = colour)
	plt.xlabel('Chainage (m)')
	plt.ylabel('Elevation (m '+str(height_datum)+')')
	plt.grid('on')
	
	
	return pd.DataFrame(np.r_[[chain,zvals]].T,columns = ['chain','z'])

#click event
def onclick(event,xoffset,yoffset):
	"""onclick event to print coordinates from mouse click on plan"""
	print(event.xdata+xoffset,event.ydata+yoffset)

#PLOT SURFACE ELEVATION CONTOUR PLAN
def plot_xyz_grid_in_plan(xyzgrid,section_line = '',contour_levels = '',contour_space = 1,structured_grid = True, map_grid = 'MGA', height_datum = 'AHD'):

	"""Plots an x,y,z grid in plan view with pre set styles
		xyzgrid: DataFrame with 'x','y' and 'z' columns
		section_line: Shapely LineString to plot on plan
		contour levels: array like if required to exactly specify contours, otherwise use contour_space
		contour_space: desired contour spacing
		structured_grid (boolean): set to True if data forms a regular grid (equal spacing in x,y) for much faster plotting
				otherwise must be False and plotting uses tricontourf
		mga_grid: (string) for axis test
		datum: (string) for colorbar text
	"""   
	if 'component_array' not in list(xyzgrid):
		xyzgrid['component_array']=0
		
	#work out contours
	if len(contour_levels)==0:
		minlev=min(xyzgrid['z'][xyzgrid['z']!=min(xyzgrid['z'])])
		maxlev=max(xyzgrid['z'][xyzgrid['z']!=max(xyzgrid['z'])])
		minlev=minlev-np.mod(minlev,contour_space)-2*contour_space
		maxlev=maxlev-np.mod(maxlev,contour_space)+contour_space
		levels=np.arange(minlev,maxlev,contour_space)
	else:
		levels = contour_levels
	
	#translate to origin
	xoffset = min(xyzgrid['x'])
	yoffset = min(xyzgrid['y'])
	
	fig,ax = plt.subplots()
	
	if not structured_grid:
		#get triangulation
		triang = trimat.Triangulation(xyzgrid['x']-xoffset,xyzgrid['y']-yoffset)
		#plot using tricontour
		c=ax.tricontourf(xyzgrid['x']-xoffset,xyzgrid['y']-yoffset,xyzgrid['z'],triangles = triang.triangles,levels=levels,cmap='terrain')
		ax.tricontour(xyzgrid['x']-xoffset,xyzgrid['y']-yoffset,xyzgrid['z'],triangles = triang.triangles,levels=levels,linewidths=0.05,colors='k')
	else:
		for i in range(max(xyzgrid['component_array'])+1):
			#convert z value to array of size (numy, numx)
			(zvals,xvals,yvals) = zvals_to_raster_array(xyzgrid[xyzgrid['component_array']==i].copy(),column_id = 'z')
			#plot using contour (note, must reverse yvals)
			c = ax.contourf(xvals-xoffset,yvals[-1::-1]-yoffset,zvals,levels=levels,cmap = 'terrain',origin = 'upper')
			ax.contour(xvals-xoffset,yvals[-1::-1]-yoffset,zvals,levels=levels,linewidths=0.05,colors='k',origin = 'upper')
   
   
	if type(section_line) == shp.LineString:
		s1 = np.array(section_line)[:,0:2]
		ax.plot(s1[:,0] - xoffset, s1[:,1] - yoffset,'r')
		
	fig.canvas.mpl_connect('button_press_event',lambda event: onclick(event,xoffset,yoffset))
	print('click on map canvas to print coordinates')
	
	plt.axis('equal')
	plt.title('Surface Elevation (m '+str(height_datum)+')')
	plt.xlabel('	 '+str(map_grid)+' Eastings	NOTE: local map coordinates (0,0) = ( '+str(xoffset)+' , '+str(yoffset)+' )')
	plt.ylabel(str(map_grid)+' Northings')
	clb=plt.colorbar(c)
	clb.set_label('m '+str(height_datum))
	plt.axis('equal')
	plt.grid('on')

#PLOT BOREHOLES on plan
def plot_bh_on_plan(dfholes, line, xoffset = 0, yoffset = 0, max_offset = 100):
   
	print('plotting boreholes on plan ...')
	#extract nearby holes (for plan plotting)
	holes1 = get_points_near_ctrl_line(line,dfholes,max_offset)

	plt.plot(holes1['x']-xoffset,holes1['y']-yoffset,'og')
	for i in range(len(holes1)):
		plt.text(holes1.iloc[i]['x']-xoffset,holes1.iloc[i]['y']-yoffset,str(holes1.iloc[i]['borehole']))

#PLOT SLOPE GRADIENT
def plot_xy_slopegradient_in_plan (xyzgrid,section_line = '',contour_levels = '',structured_grid = True, map_grid = 'MGA'):

	""" Plots up slope gradient to pre set styles
		xyzgrid: dataframe with x,y and grade columns (grade being the slope gradient)
		section_line: Shapely LineString to plot on plan
		contour_levels: Set to '' for automatic or else input a list/array of contour values
		map_grid: String, information for axis labels only
	"""
	
	if 'component_array' not in list(xyzgrid):
		xyzgrid['component_array']=0
	
	if len(contour_levels)==0:
		levels = 0.05*np.arange(40)
	else:
		levels = contour_levels
		
	#translate to origin
	xoffset = min(xyzgrid['x'])
	yoffset = min(xyzgrid['y'])

	fig,ax = plt.subplots()
	
	if not structured_grid:
		#get triangulation
		triang = trimat.Triangulation(xyzgrid['x']-xoffset,xyzgrid['y']-yoffset)
		
		c=ax.tricontourf(xyzgrid['x'] - xoffset,xyzgrid['y'] - yoffset,xyzgrid['grade'],triangles = triang.triangles,levels=levels,cmap='terrain')
		ax.tricontour(xyzgrid['x'] - xoffset,xyzgrid['y'] - yoffset,xyzgrid['grade'],triangles = triang.triangles,levels=levels,linewidths=0.05,colors='k')
	else:
		for i in range(max(xyzgrid['component_array'])+1):
			#convert z value to array of size (numy, numx)
			(zvals,xvals,yvals) = zvals_to_raster_array(xyzgrid[xyzgrid['component_array']==i].copy(),column_id = 'grade')
			#plot using contour (note, must reverse yvals)
			c = ax.contourf(xvals-xoffset,yvals[-1::-1]-yoffset,zvals,levels=levels,cmap = 'terrain',origin = 'upper')
			ax.contour(xvals-xoffset,yvals[-1::-1]-yoffset,zvals,levels=levels,linewidths=0.05,colors='k',origin = 'upper')
   
	
	if type(section_line) == shp.LineString:
		s1 = np.array(section_line)[:,0:2]
		ax.plot(s1[:,0] - xoffset, s1[:,1] - yoffset,'r')
	
	fig.canvas.mpl_connect('button_press_event',lambda event: onclick(event,xoffset,yoffset))  
	print('click on map canvas to print coordinates')
	
	plt.axis('equal')
	plt.title('Slope gradient (V:H)')
	plt.xlabel('	 '+str(map_grid)+' Eastings	NOTE: local map coordinates (0,0) = ( '+str(xoffset)+' , '+str(yoffset)+' )')
	plt.ylabel(str(map_grid)+' Northings')
	clb=plt.colorbar(c)
	clb.set_label('Slope gradient (V:H)')
	plt.axis('equal')
	plt.grid('on')   
	
def find_DEMs_from_line(section_line, DEM_directory = '', num_splits = 4):

	"""goes through all .asc files in DEM_directory to find the first one containing the points generated by
	splitting the section_line into "num_splits" equally spaced points
	section_line: LineString or list of LineStrings
	num_splits: number of equally spaced points along section_line to check

	DEM_directory: string, directory location containing all DEM files
	set DEM_directory = '' to search current directory
	WARNING: if copying directory from windows must put as raw string eg DEM_directory = r'C:\..somelocation\..Desktop\..python_codes .... note the r in front of the string
	 
	Returns: list of DEM file names"""
	
	
	line1 = section_line
	
	if type(line1) == shp.LineString:
		line1 = [line1]
	
	if not type(line1)==list:
		raise ValueError('ERROR: Section_line must be a list of shapely linestrings')
	else:
		for i in range(len(line1)):		   
			splitslist = np.array([np.array(line1[i].interpolate(i1)) for i1 in np.linspace(0,line1[i].length,num=num_splits)])
			if i==0:
				filename = find_DEMs_from_points(splitslist, DEM_directory = DEM_directory) #Find valid DEM at equally spaced coordinates
			else:
				filename = np.hstack((filename,find_DEMs_from_points(splitslist, DEM_directory = DEM_directory)))
		filename = list(np.unique(filename))
	
	if len(filename)==1:
		filename = filename[0]
	
	return filename

def find_DEMs_from_points(points,DEM_directory = ''): 

	""" points: array like 2d coordinates
		DEM_directory: string, directory location containing all DEM files
		set DEM_directory = '' to search current directory
		WARNING: if copying directory from windows must put as raw string eg DEM_directory = r'C:\..somelocation\..Desktop\..python_codes .... note the r in front of the string
		"""

	
	
	if DEM_directory != '':
		dir1 = os.getcwd()
		os.chdir(DEM_directory)
	
	filelist = os.listdir()
	#extract only the DEM (.asc) files
	filelist = [filelist[i] for i in range(len(filelist)) if filelist[i][-3::]=='asc']

	filename = []
	found = np.zeros(points.shape[0])
	
	#check each filename against each of the points
	for name in filelist:
		extents = get_dem_extents(name)
		for i in range(len(points)):
			if found[i] == 0:  #check that point has not yet had its DEM located
				result = check_if_point_in_extents(points[i],extents)
				if result:
					filename.append(name)
					found[i]=1
   
	if DEM_directory != '':
		os.chdir(dir1)
	
	
	return np.unique(filename)
	 
def get_dem_extents(demname):
	
	"""returns a tuple (xmin,ymin,xmax,ymax) specifying the DEM extents
	   demname is filename of a .asc file in the current directory"""
	
	
	a=[]
	with open(demname) as file:
		for row in islice(file,6):
			a.append(tuple(row.strip('\n').split(' ')))
	a = dict(a)
			
	a['ncols'] = int(a['ncols'])
	a['nrows'] = int(a['nrows'])
	a['xllcorner'] = float(a['xllcorner'])
	a['yllcorner'] = float(a['yllcorner'])
	a['cellsize'] = float(a['cellsize'])
	
	xmin = a['xllcorner']
	xmax = a['xllcorner'] + (a['ncols']-1)*a['cellsize']
	ymin = a['yllcorner']
	ymax = a['yllcorner'] + (a['nrows']-1)*a['cellsize']
	
	extents = (xmin,ymin,xmax,ymax)
	return extents

def check_if_point_in_extents(point,extents):
	"""check if a 2D point lies within the bounding box extents = (xmin,ymin,xmax,ymax)
		returns: boolean"""
		
	if (point[0]>=extents[0]) & (point[0]<=extents[2]) & \
	   (point[1]>=extents[1]) & (point[1]<=extents[3]):
		return True
	else:
		return False

	return
	
def OLD_find_DEMs_from_points(points,DEM_directory = ''): 

	""" points: array like 2d coordinates
		DEM_directory: string, directory location containing all DEM files
		set DEM_directory = '' to search current directory
		WARNING: if copying directory from windows must put as raw string eg DEM_directory = r'C:\..somelocation\..Desktop\..python_codes .... note the r in front of the string
		"""
	print('finding DEM ..')
	
	
	if DEM_directory != '':
		dir1 = os.getcwd()
		os.chdir(DEM_directory)
	
	filelist = os.listdir()
	#extract only the DEM (.asc) files
	filelist = [filelist[i] for i in range(len(filelist)) if filelist[i][-3::]=='asc']

	filename = []
	
	for name in filelist:
		fopen = open(name)
		content = fopen.readlines()
		numx =int(re.findall('\d+',content[0])[0])
		numy = int(re.findall('\d+',content[1])[0])
		grid_space = int(re.findall('\d+',content[4])[0])
		
		xmin = int(re.findall('\d+',content[2])[0])
		ymin = int(re.findall('\d+',content[3])[0])
		xmax = xmin + numx*grid_space
		ymax = ymin + numy*grid_space
		
		for point in points:
			if (point[0]>=xmin)&(point[0]<=xmax):
				if (point[1]>=ymin)&(point[1]<=ymax):
					filename.append(name)
		
		fopen.close()
   
   
	if DEM_directory != '':
		os.chdir(dir1)
	
	
	return np.unique(filename)
 
def get_xyz_grid_from_DEM(filename, DEM_directory = '',grid_space = 1,fill_value=0.0):

	"""Outputs an x,y,z dataframe from an .asc format input file, 
	Input: filename (string) eg 'Elevations.asc'  
	grid_space (float): desired grid spacing  
	fill_value (float): value to fill when no_data is encountered


	
	DEM_directory: string, directory location containing all DEM files
	set DEM_directory = '' to search current directory
	WARNING: if copying directory from windows must put as raw string eg DEM_directory = r'C:\..somelocation\..Desktop\..python_codes .... note the r in front of the string """

	if DEM_directory != '':
		dir1 = os.getcwd()
		os.chdir(DEM_directory)

	file = open(filename)
	templist = [line.split() for line in file.readlines()]
	meta = templist[0:6]   #Meta data
	
	#Read meta data
	#nx = int(meta[0][1])
	#ny = int(meta[1][1])
	#xllcorner = round(float(meta[2][1]),3)
	#yllcorner = round(float(meta[3][1]),3)
	#cellsize = round(float(meta[4][1]),3)
	#nodata_value = float(meta[4][1])
	
	#Read meta data
	nx = int(meta[0][1])
	ny = int(meta[1][1])
	xllcorner = round(float(meta[2][1]),3)
	yllcorner = round(float(meta[3][1]),3)
	cellsize = round(float(meta[4][1]),6)
	nodata_value = float(meta[5][1])
	
	if grid_space == '':
		grid_space = cellsize
		
	#set coarse_factor
	coarse_factor = max(round(grid_space/cellsize),1)
	
	#Get z values
	templist = templist[6::]
	templist = [templist[i][0::coarse_factor] for i in np.arange(0,len(templist),coarse_factor)]
	datavals = np.array(templist).astype(float)
	file.close()	



	#create xy cordinates
	x = cellsize * np.arange(nx)[0::coarse_factor] + xllcorner
	y = cellsize * np.arange(ny)[-1::-1][0::coarse_factor] + yllcorner
	x, y = np.meshgrid(x,y)

	#create grid
	grid = pd.DataFrame()
	grid['x'] = x.flatten()
	grid['y'] = y.flatten()
	grid['z'] = datavals.flatten()
	
	grid['z'][grid['z']==nodata_value] = fill_value
	
	if DEM_directory != '':
		os.chdir(dir1)
	
	return grid
	
def get_xyz_grid_from_multipleDEM(filenames,DEM_directory = '',grid_space='',fill_value=0.0):
	
	""" get one combined xyz dataframe from multiple .asc files
		Outputs an x,y,z dataframe from an esri .asc format input file, 
		Input: filename (list of strings) eg ['Elevations1.asc', 'Elevations2.asc', ...]	
		grid_space (float): desired grid spacing  
		fill_value (float): value to fill when no_data is encountered
		it is good to have a small offset to avoid numerical rounding error and omission/doubling up of raster cell selection

		
		
		DEM_directory: string, directory location containing all DEM files
		set DEM_directory = '' to search current directory
		WARNING: if copying directory from windows must put as raw string eg DEM_directory = r'C:\..somelocation\..Desktop\..python_codes .... note the r in front of the string """
		   
	
	
	if type(filenames) == str or type(filenames) == np.str_:
		filenames = [filenames]
	
	coords = get_xyz_grid_from_DEM(filenames[0],DEM_directory = DEM_directory, grid_space=grid_space,fill_value=fill_value)
	coords['component_array']=0
	
	if len(filenames)>0:
		for i in range(1,len(filenames)):
			coords_new = get_xyz_grid_from_DEM(filenames[i],DEM_directory = DEM_directory, grid_space=grid_space,fill_value=fill_value)
			coords_new['component_array']=i
			coords = pd.concat([coords,coords_new])
	
	return coords

def get_points_at_equal_chainages(line,point_space = 10):
	""" line: shapely LineString
		outputs: 
		1. 2D ARRAY of xy coordinates equally spaced along line according to point_space
		2. 1D ARRAY chainage of each point
		points may not exactly match line due to corners not being located at point_space"""
 
	len1 = line.length
	dists = np.arange(0,len1,point_space)
	if dists[-1]!=len1:
		dists = np.append(dists,len1)

		
	coords = [line.interpolate(dists[i]) for i in range(len(dists))] #list of shapely points

	return (np.array(shp.LineString(coords))[:,0:2], np.array(dists))  #(xy coords, chainage)
	
def get_zvals_from_grid(line,xyzgrid):
	""" Input: line (LineString or ARRAY with first two columns the xy coords)
		get z values at every point in line from xyzgrid 
		xyzgrid: DataFrame with x,y,z columns or 3 column ARRAY
		
		returns: array of z values for each point on line  """
	

		
	#create dataframe if np.array input
	if type(xyzgrid) == np.ndarray:
		gridnew = pd.DataFrame(xyzgrid,columns = [['x','y','z']])
	else:
		gridnew = xyzgrid.copy()
	
	#get xy coords of line
	if type(line) == np.ndarray:
		xy = line[:,0:2]
	else:
		xy = np.array(line)[:,0:2]	 
		
	#reduce grid size for use by griddata function
	extents = (min(xy[:,0]),min(xy[:,1]),max(xy[:,0]),max(xy[:,1]))
	gridnew = crop_xyz_grid_by_extents(gridnew,extents) 
	

	#convert to np array to avoid dodgy indexing issues as results looked funny
	gridnew = gridnew[['x','y','z']].values

	zvals = sinterp.griddata(gridnew[:,0:2],gridnew[:,2],xy,method = 'nearest')   #get z coordinates - nearest seems to give most reliable results
	
	return zvals

def get_both_zvals_and_equal_chainages(line,xyzgrid,point_space = 10,smoothing_window = 3):
	"""
	input: line: (shapely linestring of array with first two columns xy coordinates)
		   xyzgrid: DataFrame with x,y,z columns or 3 column ARRAY
		   point_space (float): desired spacing of points along line
		   smoothing window (int): Uses pandas rolling window to take a moving average windown on the z values
	returns: Dataframe with columns [chain,x,y,z]
	"""
	
	chain = get_points_at_equal_chainages(line,point_space = point_space)
	zvals = get_zvals_from_grid(chain[0],xyzgrid)
	
	if smoothing_window > 0:
		zvals = pd.DataFrame(zvals).rolling(smoothing_window,center = True,min_periods = 1).mean()

	df = pd.DataFrame()
	df['chain'] = chain[1]
	df['x'] = chain[0][:,0]
	df['y'] = chain[0][:,1]
	df['z'] = zvals
	
	return df

def slope_gradient_from_xyz_grid(xyzgrid):
	
	#raise NotImplementedError('Do not use slopegradient: Needs further testing')
	
	"""input a regular spaced (ie constant grid spacing, equal spacing in x and y directions) 
		xyzgrid (DF), must have columns x,y and z. and (if applicable) a 'component_array' column to split up into such units
		Returns: A dataframe with x,y,z and "grade" columns with grade being the magnitude of the slope eg 0.1 for a 1 in 10 slope"""
	
	print('CAUTION: Hillshade function - z-units are assumed to match xy-units (eg metres)')
	print('CAUTION: Shows up poor quality horizontal ailiasing on the 1sec govt DEM. Use 5m or 1m DEM as a minimum')

	if 'component_array' not in list(xyzgrid):
		xyzgrid['component_array']=0
	
	grid_new = pd.DataFrame(columns = list(xyzgrid).append('grade'))
	
	
	
	for i in range(max(xyzgrid['component_array'])+1):
	
		grid_temp = xyzgrid[xyzgrid['component_array']==i].copy()
	
		#Get step size and number in x/y
		xvals = np.sort(np.unique(np.array(grid_temp['x'])))
		xstep = xvals[1]-xvals[0]	   
		yvals = np.sort(np.unique(np.array(grid_temp['y'])))
		ystep = yvals[1]-yvals[0]
		
		#Flag bad grid spacing and/or determine grid_space
		if np.abs(xstep -ystep) >0.01:
			raise ValueError("Error: Grid spacing not equal in x and y")
		else:
			grid_space = xstep
		
		numx = len(xvals)
		numy = len(yvals)
		
		print('x_space: '+str(xstep)+'  y_space: '+str(ystep)+'   num_x: '+str(numx)+'   num_y: '+str(numy))
		
		print('caluculating slope gradients...')
		
		
		####   Gradiant calculation
		zvals = np.array(grid_temp[['x','y','z']])[:,2].reshape((numy,numx))	#values in grid format
		grad1=np.sqrt(np.gradient(zvals,grid_space,grid_space)[0]**2+np.gradient(zvals,grid_space,grid_space)[1]**2)
		grad1=grad1.flatten()
		grid_temp['grade'] = grad1
		####
				
				

		grid_new = pd.concat([grid_new,grid_temp])
	
 
	
	print('NaN replaced with 0.0')
	grid_new = grid_new.fillna(value=0.0)
	
	return grid_new

def crop_xyz_grid_by_extents(xyzgrid,extents):
	
	""" crop grid in xy plane
	   xyzgrid: Dataframe with x,y,z columns or ARRAY with first 3 columns being x,y,z coordinates
	   extents: (xmin,ymin,xmax,ymax)
	   returns: cropped Dataframe
	"""
	
	if type(xyzgrid) == np.array:
		xyzgridnew = pd.DataFrame(xyzgrid,columns=['x','y','z'])
		xyzgridnew['component_array']=0
	else:
		xyzgridnew = xyzgrid.copy()
	
	flag = (xyzgridnew['x'] > extents[0]) & (xyzgridnew['x'] < extents[2]) & (xyzgridnew['y'] > extents[1]) & (xyzgridnew['y'] < extents[3])

	return xyzgridnew[flag]

def smooth_zvals_from_xyz_grid(xyzgrid, smooth_factor = 3, mode = 'nearest'):

	""" smooth z values in xyzgrid by averaging over a window
		xyzgrid = Well ordered (ie not random) rectangular grid (dataframe) with x,y,z columns, or, if not rectangular dataframe must have 'component_array' column which splits it into well ordered rectangular grids
		smooth_factor = number of adjacent cells to use in smoothing
		mode = scipy.ndimage.uniform_filter mode  (check documentation for that function)
		
		returns: Dataframe with z values smoothed
	"""
	
	if 'component_array' not in list(xyzgrid):
		xyzgrid['component_array']=0
	
	
	grid_new = pd.DataFrame(columns=list(xyzgrid))
	
	
	for i in range(max(xyzgrid['component_array'])+1):
	
		grid_temp = xyzgrid[xyzgrid['component_array']==i].copy()
		
		numx = len(np.unique(np.array(grid_temp['x'])))  
		numy = len(np.unique(np.array(grid_temp['y'])))		
		zvals = np.array(grid_temp['z'])

		zvals = ndimg.uniform_filter(zvals.reshape((numy,numx)),size = smooth_factor,mode = mode).flatten()
	
		grid_temp['z'] = zvals
		grid_new = pd.concat([grid_new,grid_temp])
		
	return grid_new

def grid_from_shapely_polygon(poly=[],extents=[],grid_space=100,plot_grid=True):
	"""Input a shapley polygon (or multi-polygon). The function returns a 2d numpy array 
		of grid coordinates which cover the convex hull of poly
		if extents are supplied extents=(xmin,ymin,xmax,ymax) a polygon will be created from those extents"""
	
	if len(extents)>0:
		print('creating shapely polygon from supplied extents')
		xmin,ymin,xmax,ymax=extents[0],extents[1],extents[2],extents[3]
		coords1=((xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax),(xmin,ymin))
		poly=shp.Polygon(coords1)
		conv=poly
		
	if len(extents)==0:
		conv=poly.convex_hull
		xmin,ymin,xmax,ymax=poly.bounds
	
	#####  OLD
	#xmax=np.ceil((xmax-xmin)/grid_space)*grid_space+xmin
	#ymax=np.ceil((ymax-ymin)/grid_space)*grid_space+ymin
	#xx=np.linspace(xmin,xmax,int((xmax-xmin)/grid_space)+1)
	#yy=np.linspace(ymin,ymax,int((ymax-ymin)/grid_space)+1)
	#####  End
	
	#rectangular grid extents
	xx=np.arange(xmin,xmax,grid_space)
	yy=np.arange(ymin,ymax,grid_space)
	x,y=np.meshgrid(xx,yy)
	grid=np.vstack((x.flatten(),y.flatten())).T 

	#cropping of grid to input polygon/extents
	print(str(len(grid))+ 'meshgrid points ....')  
	path1=matpath.Path(conv.boundary.coords)
	ids=path1.contains_points(grid,radius=0.001)
	grid=grid[ids]
	
	#add in points on boundary
	#bdy = conv.boundary
	#intvals = np.linspace(0,bdy.length,int(bdy.length/grid_space))
	#bpts = np.array([np.array(bdy.interpolate(item)) for item in intvals])
	#grid = np.vstack((grid,bpts))
	
	
	#extract points strictly interior to poly
	flag = [poly.contains(shp.Point(grid[i])) for i in range(len(grid))]
	grid = grid[flag]
	
	print(str(len(grid))+ 'grid points inside or on polygon.')

	if plot_grid==True:
		#PLOT GRID
		fig,ax=plt.subplots()
		patch1=dc.PolygonPatch(poly,facecolor='none',edgecolor='k',linewidth=0.7,zorder=9)
		#patch2=dc.PolygonPatch(conv,facecolor='none',edgecolor='m',linewidth=2,zorder=10)
		plt.plot(grid[:,0],grid[:,1],'.r')
		ax.add_patch(patch1)
		#ax.add_patch(patch2)
		plt.axis('equal')
		plt.show()
	
	#extents=conv.bounds
	#xmin=extents[0]
	#ymin=extents[1]
	#xmax=extents[2]
	#ymax=extents[3]
	
	return grid #(grid,(xmin,ymin,xmax,ymax))
 
def zvals_to_raster_array(xyzgrid,column_id = 'z'):

	"""Input: list of coordinates with z values (DataFrame with 'x','y',column_id columns)
			  column_id: name of the column containing the z values
	   Returns: [0] an array of the z values of size (numy x numx) which can be used in imshow plotting.
				[1] xvalues 1D array
				[2] yvalues 1D array
	Input coordinates must be a rectangular grid with equal grid spacings in x and y directions
	"""
	
	grid_temp = xyzgrid.copy()

	#Get step size and number in x/y
	xvals = np.sort(np.unique(np.array(grid_temp['x'])))
	xstep = xvals[1]-xvals[0]	   
	yvals = np.sort(np.unique(np.array(grid_temp['y'])))
	ystep = yvals[1]-yvals[0]
	
	#Flag bad grid spacing and/or determine grid_space
	if np.abs(xstep -ystep) >0.01:
		raise ValueError("Error: Grid spacing not equal in x and y")
	else:
		grid_space = xstep
	
	numx = len(xvals)
	numy = len(yvals)
	
	#print('x_space: '+str(xstep)+'  y_space: '+str(ystep)+'   num_x: '+str(numx)+'   num_y: '+str(numy))
	
	zvals = grid_temp[column_id].values.reshape((numy,numx))
	return (zvals,xvals,yvals)

def contours_to_shapefile(cs,filename):
	"""Extract polylines from matplotlib contour object and create shapefile
	input: cs must be of the form cs=plt.contour(x,y,z) (ie a matplotlib contour plot)
	writes to filename, filename must be form 'myshape.shp'. Make sure filename does not already exist
	"""

	
	#OLD METHOD USING FIONA - removed on 29-8-2020
	#set schema for shapefile
	#schema={'geometry':'LineString','properties':{'level':'float'}}
	#d={}
	#with fiona.open(filename,'w','ESRI Shapefile',schema) as layer:
	#	count=0
	#	for i in range(num_lev):
	#		 for j in range(np.shape(cs.collections[i].get_paths())[0]):
	#			 if cs.collections[i].get_paths()[j].vertices.shape[0]>1: #check more than one point in the contour
	#				 d['geometry']=shp.mapping(shp.LineString(cs.collections[i].get_paths()[j].vertices)) #points (LineString) on particular contour line
	#				 d['properties']={'level':cs.levels[i]}
	#				 layer.write(d)
	#				 count+=1
	#				 d={}	
	

	print('NOTE: Input contours must must be tricontour and NOT tricontourf input')

	num_lev=np.shape(cs.levels)[0]

	if '.shp' in filename:
		filename = filename[0:-4]

	w = pyshp.Writer(filename)
	w.field('level', 'N', decimal=10)
	for i in range(num_lev):
		for j in range(np.shape(cs.collections[i].get_paths())[0]):
			if (cs.collections[i].get_paths()[j].vertices).shape[0] > 1: #need more than one point for a line
				coords = [shp.mapping(shp.LineString(cs.collections[i].get_paths()[j].vertices))['coordinates']]		
				w.record(cs.levels[i])
				w.line(coords)
		
	return

def contours_to_smoothed_polygons(cs, tolerance, point_space, plot=False):
	"""takes an input matplotlib contour set 'cs'
		and returns:
			1. list of smoothed shapely polygons
			2. list of point arrays along polygon boundaries
			
		tolerance = shapely polygon smoothing tolerance
		point_space = desired point space along line
	"""
	
	temp = cs.collections[0].get_paths()
	allpolys = []
	allpoints = []
	for i in range(len(temp)):
		
		#create polygon
		poly1 = shp.Polygon(temp[i].vertices).simplify(tolerance)
		allpolys.append(poly1)
		
		#create points
		line1 = shp.LineString(poly1.boundary)
		interp_dist = np.linspace(0,line1.length,int(line1.length/point_space))
		pts = [np.array(line1.interpolate(interp_dist[j])) for j in range(len(interp_dist))]
		pts = np.array(pts)
		allpoints.append(pts)
		
	if plot:
		for i in range(len(allpolys)):
			poly1 = allpolys[i]
			plt.plot(np.array(poly1.boundary)[:,0],np.array(poly1.boundary)[:,1])
			plt.plot(allpoints[i][:,0],allpoints[i][:,1],'o')
			
		plt.axis('equal')
		plt.grid(True)
		plt.show()
	
	return allpolys, allpoints	
	
def remove_close_points(datapoints, min_dist = 2.5):

	"""removes close points from array of datapoints
		min_dist (float): minimum allowed distance between any two points
	"""
	print('removing close points ..')
	t2 = scipy.spatial.kdtree.KDTree(datapoints)
	excluded = []
	
	for i in range(len(datapoints)):
		if i not in excluded:
			nearpoints = t2.query_ball_point(datapoints[i],r = min_dist)
			if len(nearpoints)>1:
				for j in range(len(nearpoints)):
					if not nearpoints[j] == i:
						excluded.append(nearpoints[j])
	
	excluded = np.unique(excluded)
	
	included = []
	for i in range(len(datapoints)):
		if i not in excluded:
			included.append(i)
	
	print(len(datapoints),' total points')
	print(len(excluded),' excluded points')
	print(len(included),' retained points')

	return datapoints[included]
	
def DEPRECIATEDcoord_transform(x,y,inprojection = 'epsg:4326', outprojection = 'epsg:28355'):
	"""depreciated due to changes in pyproj"""
	"""x,y: arraylike input coordinates """
	
	#inProj = Proj(init=inprojection) #2020-11-06 depreciated 
	#outProj = Proj(init=outprojection) #2020-11-06 ddepreciated
	inProj = Proj(inprojection)
	outProj = Proj(outprojection)
	transformer = Transformer.from_proj(inProj,outProj)

	return np.array([transformer.transform(x[i],y[i]) for i in range(len(x))])

	
def coord_transform(x,y,inprojection = 4326, outprojection = 28355):

	"""x,y: arrays input coordinates """
	
	inProj = Proj(inprojection)
	outProj = Proj(outprojection)
	transformer = Transformer.from_crs(inprojection,outprojection,always_xy = True)

	return np.array([transformer.transform(x[i],y[i]) for i in range(len(x))])
	
	
	
def split_df_to_intervals(dfdata,x_col, y_col, other_y_col = '', min_interval = 250, max_range = 25, show_plot = False):

	"""Takes a dataframe and splits into intervals based on the maximum allowed range of y values over an interval
		INput: dfdata (dataframe)
			x_col (string): column name of x axis data
			y_col (string): column name of y axis data
			other_y_col (string): column name of an additional y axis data to include in output
		
		Returns:
			Dataframe with max, min and mean values over each of the intervals
	"""
	
	xvals = dfdata[x_col].values
	yvals = dfdata[y_col].values
	if len(other_y_col) >0:
		other_y_values = dfdata[other_y_col].values
	
	fromidx = 0
	a1 =[]
	
	while fromidx+1 < len(xvals):
		c0 = 0
		for idx in np.arange(fromidx+1,len(xvals)+1):
			yrange = np.max(yvals[fromidx:idx+1]) - np.min(yvals[fromidx:idx+1])
			if yrange > max_range:
				c0 = idx-1
				if len(other_y_col) >0:
					a1.append([xvals[fromidx],xvals[c0],np.max(yvals[fromidx:c0]),\
					np.min(yvals[fromidx:c0]),np.mean(yvals[fromidx:c0]),np.max(other_y_values[fromidx:c0]),\
					np.min(other_y_values[fromidx:c0]),np.mean(other_y_values[fromidx:c0])])
				else:
					a1.append([xvals[fromidx],xvals[c0],np.max(yvals[fromidx:c0]),\
					np.min(yvals[fromidx:c0]),np.mean(yvals[fromidx:c0])])
				fromidx = c0
				break
		
		if idx == len(xvals):
			if len(other_y_col) >0:
				a1.append([xvals[fromidx],xvals[-1],np.max(yvals[fromidx::]),\
				np.min(yvals[fromidx::]),np.mean(yvals[fromidx::]),\
				np.max(other_y_values[fromidx::]),np.min(other_y_values[fromidx::]),np.mean(other_y_values[fromidx::])])
			else:	
				a1.append([xvals[fromidx],xvals[-1],np.max(yvals[fromidx::]),\
				np.min(yvals[fromidx::]),np.mean(yvals[fromidx::])])
			break
	
	if len(other_y_col) >0:
		dataout = pd.DataFrame(a1,columns = ['from','to','max-'+y_col,'min-'+y_col,'mean-'+y_col, \
		'max-'+other_y_col,'min-'+other_y_col,'mean-'+other_y_col])
	else:
		dataout = pd.DataFrame(a1,columns = ['from','to','max-'+y_col,'min-'+y_col,'mean-'+y_col])

	if show_plot:
		cols = ['max-'+y_col,'min-'+y_col,'mean-'+y_col]
		for i in range(3):
			xdat = np.vstack((dataout['from'].values,dataout['to'].values)).T.flatten()
			ydat = np.vstack((dataout[cols[i]].values,dataout[cols[i]].values)).T.flatten()
			plt.plot(xdat,ydat)
		
		plt.plot(xvals,yvals)
		plt.show()
	
	return dataout

def smooth_data_with_rbf(xydata,vals,gridsize=100,function='multiquadric',epsilon=75,smooth=5):
	
	""" Smooth a 2D dataset using scipy.interpolate.Rbf
	
		xydata: 2D array of coordinates
		vals: 1D array of values
		gridsize: Float desired gridsize of output data
		
		---- scipy.interpolate.Rbf kwargs -----  (note: function = 'multiquadric',epsilon=75,smooth=5 works ok for groundwater around an open pit) 
		function : str or callable, optional
		The radial basis function, based on the radius, r, given by the norm
		(default is Euclidean distance); the default is 'multiquadric'::

			'multiquadric': sqrt((r/self.epsilon)**2 + 1)
			'inverse': 1.0/sqrt((r/self.epsilon)**2 + 1)
			'gaussian': exp(-(r/self.epsilon)**2)
			'linear': r
			'cubic': r**3
			'quintic': r**5
			'thin_plate': r**2 * log(r)

		If callable, then it must take 2 arguments (self, r). The epsilon
		parameter will be available as self.epsilon. Other keyword
		arguments passed in will be available as well.

		epsilon : float, optional
			Adjustable constant for gaussian or multiquadrics functions
			- defaults to approximate average distance between nodes (which is
			a good start).
		smooth : float, optional
			Values greater than zero increase the smoothness of the
			approximation. 0 is for interpolation (default), the function will
			always go through the nodal points in this case.	
			
		Returns: array [gridx,gridy,smoothed_vals]
		.
	"""
	
	#input data check (raise value errors)
	if not ((type(xydata) == np.ndarray) and (type(vals) == np.ndarray)):
		raise ValueError('Error: Input data must be numpy arrays')
	
	if not ((xydata.shape[-1] ==2) and (len(vals.shape) == 1)):
		raise ValueError('Error: xydata must be 2D, vals must be 1D')
	
	
	###### create grid
	coords = xydata[:,0:2]
	bdry = shp.MultiPoint(coords[:,0:2]).envelope #bounding rectangle
	bdrycoords = np.array(bdry.boundary.coords[:])
	xmin, xmax, ymin, ymax = min(bdrycoords[:,0]), max(bdrycoords[:,0]), min(bdrycoords[:,1]), max(bdrycoords[:,1]) #extents
	x1,y1 = np.meshgrid(np.linspace(xmin,xmax,int((xmax-xmin)/gridsize)),np.linspace(ymin,ymax,int((ymax-ymin)/gridsize)))
	coords = np.vstack((x1.flatten(),y1.flatten())).T
	
	## extract grid points inside convex hull of input xydata
	flag = np.zeros(len(coords),dtype=bool)
	poly1 = shp.MultiPoint(xydata[:,0:2]).convex_hull.buffer(gridsize)
	for j in range(len(coords)):
		flag[j] = shp.Point(coords[j]).intersects(poly1)
	coords = coords[flag]
	##
	
	###### END CREATE GRID
	
	# create rbf interpolation of input data
	print('create Rbf interpolator ..')
	rbfi = sinterp.Rbf(xydata[:,0],xydata[:,1], vals, function = function, epsilon = epsilon, smooth = smooth)
	
	#get smoothed values
	print('get interpolated values ..')
	smoothed_vals = rbfi(coords[:,0],coords[:,1])
	
	return np.hstack((coords,np.c_[smoothed_vals]))  #return array [x,y,smoothed_vals]
	
	
def dxf_to_csv(filename):
	""" extracts a dataframe [x,y,z] from a dxf mesh file
	
		Input:
			filename: string name of a dxf file (eg 'example.dxf')
		Returns:
			DataFrame [x,y,z]
	"""
	
	print('searching for AcDbFace in file ..')
	c=0 #look for AcDbFace
	startface=1e16
	checkset = range(2,25,2) #checkset to denote and append the lines positioned 2,4,6,.. after 'AcDbFace' is found
	coords=[]
	with open(filename, 'r') as f:
		for line in f:
			if np.mod(c,1000000)==0:
				print(c)
			if line == 'AcDbFace\n':  #TRIGGER text to start extract of coordinates
				startface=c
			if c-startface in checkset: #append the lines positioned 2,4,6,.. after 'TRIGGER' is found
				coords.append(line.strip('\n'))
			c+=1		

	print('len(coords) =',len(coords))
	if len(coords) == 0: #look for 3DFACE if AcDbFace not found
		print('searching for 3DFACE in file ..')
		c=0 #look for AcDbFace
		startface=1e16
		checkset = range(4,27,2) #checkset to denote and append the lines positioned 4,6,8,.. after '3DFACE' is found
		with open(filename, 'r') as f:
			for line in f:
				if np.mod(c,1000000)==0:
					print(c)
				if line == '3DFACE\n':  #TRIGGER text to start extract of coordinates
					startface=c
				if c-startface in checkset: #append the lines positioned 2,4,6,.. after 'TRIGGER' is found
					coords.append(line.strip('\n'))
				c+=1
	
	print('len(coords) =',len(coords))
	if len(coords) == 0: #look for AcDbPolyFaceMeshVertex if 3DFace not found
		print('searching for AcDbPolyFaceMeshVertex in file ..')
		c=0 #look for AcDbPolyFaceMeshVertex
		startface=1e16
		checkset = range(2,7,2) #checkset to denote and append the lines positioned 4,6,8,.. after 'AcDbPolyFaceMeshVertex' is found
		with open(filename, 'r') as f:
			for line in f:
				if np.mod(c,1000000)==0:
					print(c)
				if line == 'AcDbPolyFaceMeshVertex\n':  #TRIGGER text to start extract of coordinates
					startface=c
				if c-startface in checkset: #append the lines positioned 2,4,6,.. after 'TRIGGER' is found
					coords.append(line.strip('\n'))
				c+=1	
			
	
	xyzvals = np.array(coords).reshape(int(len(coords)/3),3).astype(float)

	xyzDF = pd.DataFrame(xyzvals,columns=['x','y','z'])
	xyzDF = xyzDF.drop_duplicates()

	return xyzDF	

	

#OLD -- def shapefile_to_shapely(input_shapefile):
#	"""input_shapefile: string shapefile filename
#		returns: list of shapely objects"""
#
#	file1 = fiona.open(input_shapefile)
#	list1 = []
#
#	for i in range(len(file1)):
#		ftype = file1[i]['geometry']['type']
#		coords = file1[i]['geometry']['coordinates']
#		if ftype == 'Point':
#			list1.append(shp.Point(coords))
#		if ftype == 'MultiPoint':
#			list1.append(shp.MultiPoint(coords))
#		elif ftype == 'LineString':
#			list1.append(shp.LineString(coords))
#		elif ftype == 'Polygon':
#			list1.append(shp.Polygon(coords[0]))
#		else:
#			'Error: shapefile contains elements which are not a "Point","MultiPoint","LineString" or "Polygon"'
#
#	return list1


#OLD -- def shapely_to_shapefile(input,save_filename):
#	
#	"""Take an input list of shapely objects (LineString, Polygon, Point) and create and save as a shapefile"""
#
#	if not type(input)==list:
#		inputdata = [input]
#	else:	
#		inputdata = input[:]
#
#	if not (type(inputdata[0])==shp.LineString or type(inputdata[0])==shp.MultiPoint or type(inputdata[0])==shp.Polygon):
#		print('Error: Only input shapely Linestrings or MultiPoint are supported at present')
#		return
#		
#	if type(inputdata[0])==shp.LineString:
#		schema={'geometry':'LineString','properties':{'level':'float'}}#,'properties':{'level':'float'}}
#	elif type(inputdata[0])==shp.MultiPoint:
#		schema={'geometry':'MultiPoint','properties':{'level':'float'}}#,'properties':{'level':'float'}}
#	elif type(inputdata[0])==shp.Polygon:
#		schema={'geometry':'Polygon','properties':{'level':'float'}}#,'properties':{'level':'float'}}
#	
#	d={}
#	with fiona.open(save_filename,'w','ESRI Shapefile',schema) as layer:
#		for i in range(len(inputdata)):
#			d['geometry']=shp.mapping(inputdata[i]) #points (LineString) on particular contour line
#			#print(d['geometry'])
#			d['properties']={'level':999.99}
#			#d['properties']={'level':cs.levels[i]}
#			layer.write(d)
#			d={}
#		
#	return
	
#def triangulate_xy_grid(xydata, min_angle = 20):
#	""""triangulates an unstructured grid of xy coords using J Schewchuck's algorithm.
#		input xydata is s dataframe containing 'x','y' columns or an array with first two columns the x and y coordinates
#		min_angle(float): minimum allowed angle in any triangle (ie generates a quality mesh)
#		returns the connectivity matrix (np.array) of triangles in the triangulation"""
	
#	if type(xydata) == pd.DataFrame:   
#		points = xydata[['x','y']].values
#	else:
#		points = xydata[:,0:2]
		
#	triang = triangle.triangulate({'vertices':points},opts='q'+str(min_angle))
#	
#	return triang['triangles']	


if __name__=="__main__":  
	
	#INPUT DATA

	line1 = geo.shapefile_to_shapely('RWB01cntrl.shp')[0]  
	start_chainage = 250
	
	demloc = r'C:\Users\Antony.Orton\Desktop\Python programs\DigElev\ARUP_DEM_section_tool\DEMs'
	DEMgrid_space = 2
	DEMsmoothing_factor = 3
	
	materials_list = ['Fill', 'Q1B', 'RSA', 'RSB', 'SHV', 'SHIV', 'SHIII','SHII']
	max_BH_offset = 50
	
	[dfholes,dfgeology] = geo.read_csv_input_files()  #Read borehole information
	
	
	#DATA PROCESS
	demnames = geo.find_DEMs_from_line(line1, DEM_directory = demloc)  #find DEM names
	grid1 = geo.get_xyz_grid_from_multipleDEM(demnames,DEM_directory=demloc,grid_space=DEMgrid_space)  #read surface elevation DEM to pandas dataframe (non-cropped, non-smoothed)
	
	#PLOTTING - section
	geo.plot_boreholes_on_section(dfholes, dfgeology, line1, materials_list = materials_list, colourmap = 'cool',\
	max_offset = max_BH_offset,start_chainage = start_chainage)
	zline = geo.plot_groundsurface_on_section(line1,grid1,smoothing_window = 7, point_space = DEMgrid_space,colour = 'r',start_chainage = start_chainage)
	
	plt.show()
	
	#PlOTTING - plan
	grid_for_plan = geo.crop_xyz_grid_by_extents(grid1,geo.get_nearby_extents(line1,max_offset = 400))
	grid_for_plan = geo.smooth_zvals_from_xyz_grid(grid_for_plan,smooth_factor = DEMsmoothing_factor)
	geo.plot_xyz_grid_in_plan(grid_for_plan,section_line = line1,contour_space = 1)
	
	plt.show()
	
	#PLOTTING - slope gradient
	grid_for_gradient = geo.crop_xyz_grid_by_extents(grid1,geo.get_nearby_extents(line1,max_offset = 400))
	grid_for_gradient = geo.smooth_zvals_from_xyz_grid(grid_for_gradient,smooth_factor = DEMsmoothing_factor)
	grid_for_gradient = geo.slope_gradient_from_xyz_grid(grid_for_gradient)   
	geo.plot_xy_slopegradient_in_plan(grid_for_gradient)
	
	plt.show()
 
	
	
	
	
	
	
	
	
	
	
   
