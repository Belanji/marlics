# trace generated using paraview version 5.6.0

#### import the simple module from the paraview
from paraview.simple import *
import math
import sys

if len(sys.argv)<5:
	print 'Please, use this script as \''+sys.argv[0]+' csvfilename [Nx] [Ny] [Nz]\''
	exit()

else:
	InputFile=str(sys.argv[1])
	Nx = int(sys.argv[2])
	Ny = int(sys.argv[3])
	Nz = int(sys.argv[4])

csvfile =  CSVReader(FileName=InputFile)

renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')

# set active view
SetActiveView(renderView1)

# Properties modified on tableToStructuredGrid1
tableToStructuredGrid1 = TableToStructuredGrid(Input=csvfile)
tableToStructuredGrid1.WholeExtent = [0, Nx, 0, Ny, 0, Nz]
tableToStructuredGrid1.XColumn = 'x'
tableToStructuredGrid1.YColumn = 'y'
tableToStructuredGrid1.ZColumn = 'z'

# reset view to fit data
renderView1.ResetCamera()

threshold1 = Threshold(Input=tableToStructuredGrid1)
threshold1.Scalars = ['POINTS', 'pt']
threshold1.ThresholdRange = [1.0,1.0]

bulkSurface= ExtractSurface(Input=threshold1)

bulkSmooth= Smooth(Input=bulkSurface)
bulkSmooth.NumberofIterations = 2000

Smin, Smax= threshold1.PointData.GetArray("S").GetRange()
Sdefect=(0.9*Smax)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Contour'
contour1 = Contour(Input=threshold1)
contour1.ComputeNormals = 0
contour1.ComputeScalars = 1
contour1.PointMergeMethod = 'Uniform Binning'
contour1.ContourBy = ['POINTS', 'S']
contour1.Isosurfaces = [Sdefect]

# show data in view
contour1Display = Show(contour1, renderView1)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = [None, '']
contour1Display.DiffuseColor = [1, .0, .0]
contour1Display.OSPRayScaleArray = 'S'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 5.0
contour1Display.SelectScaleArray = 'S'
contour1Display.GaussianRadius = 0.25
contour1Display.SetScaleArray = ['POINTS', 'S']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'S']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# show data in view
gridDisplay = Show(bulkSmooth, renderView1)

gridDisplay.Scale = [(Nx+2)/Nx, (Ny+2)/Ny, (Nz+2)/Nz]
gridDisplay.Position = [0.,0., 0.]
gridDisplay.DiffuseColor = [1.0, 1.0, 1.0]
gridDisplay.Opacity = 0.3
# trace defaults for the display properties.
gridDisplay.Representation = 'Surface'
gridDisplay.ColorArrayName = [None, '']
gridDisplay.OSPRayScaleArray = 'P'
gridDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
gridDisplay.SelectOrientationVectors = 'None'
gridDisplay.ScaleFactor = 5.0
gridDisplay.SelectScaleArray = 'None'
gridDisplay.GlyphType = 'Arrow'
gridDisplay.GlyphTableIndexArray = 'None'
gridDisplay.GaussianRadius = 0.25
gridDisplay.SetScaleArray = ['POINTS', 'P']
gridDisplay.ScaleTransferFunction = 'PiecewiseFunction'
gridDisplay.OpacityArray = ['POINTS', 'P']
gridDisplay.OpacityTransferFunction = 'PiecewiseFunction'
gridDisplay.DataAxesGrid = 'GridAxesRepresentation'
gridDisplay.PolarAxes = 'PolarAxesRepresentation'

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'S'
sPWF = GetOpacityTransferFunction('S')
sPWF.Points = [0.4000000059604645, 0.0, 0.5, 0.0, 0.45006102323532104, 1.0, 0.5, 0.0]
sPWF.ScalarRangeInitialized = 1

# Properties modified on contour1

# update the view to ensure updated data information
renderView1.Update()


# set active source
SetActiveSource(tableToStructuredGrid1)


Hide(tableToStructuredGrid1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

bounds = GetActiveSource().GetDataInformation().GetBounds()
# change solid color

# Properties modified on ptLUT
renderView1.Background = [0.7,0.7,0.7]
# current camera placement for renderView1
renderView1.CameraPosition = [3*(bounds[1]-bounds[0]), 3*(bounds[3]-bounds[2]),3*(bounds[5]-bounds[4])]
renderView1.CameraFocalPoint = [-1.5*(bounds[1]-bounds[0]), -1.5*(bounds[3]-bounds[2]), -1.5*(bounds[5]-bounds[4])]
renderView1.CameraViewUp = [0,0,1]
renderView1.ResetCamera()
renderView1.CameraParallelProjection = 1
renderView1.CameraParallelScale = bounds[5]/1.3
SaveScreenshot('iso_'+InputFile[:-4]+'.png', renderView1, ImageResolution=[1510, 1510])


renderView1.CameraPosition =   [0.5*bounds[1], -0.5*bounds[3],0.5*bounds[5]]
renderView1.CameraFocalPoint = [0.5*bounds[1], 3.5*bounds[3],0.5*bounds[5]]
renderView1.CameraViewUp = [0,0,1]
SaveScreenshot('iso_y_'+InputFile[:-4]+'.png', renderView1, ImageResolution=[1510, 1510])

renderView1.CameraPosition =   [0.5*bounds[1], 0.5*bounds[3],-0.5*bounds[5]]
renderView1.CameraFocalPoint = [0.5*bounds[1], 0.5*bounds[3], 3.5*bounds[5]]
renderView1.CameraViewUp = [1,0,0]
SaveScreenshot('iso_z_'+InputFile[:-4]+'.png', renderView1, ImageResolution=[1510, 1510])




