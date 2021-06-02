# trace generated using paraview version 5.6.0

#### import the simple module from the paraview
from paraview.simple import *
import math  
import sys


if len(sys.argv)==1:
	print 'Please, use this script as \''+sys.argv[0]+' csvfilename [px] [py] [pz]\''
	print 'csvfilename -> a input csv file containing the position, direction and pt of the director field.'
	print 'px, py and pz -> optional arguments to indicate the slice position of the 3 slices, the default value refers to the center point of the director field.'
	exit()
else:
	InputFile=str(sys.argv[1])
OutputFile= InputFile[0:-3]+'png'

paraview.simple._DisableFirstRenderCameraReset()
# create a new 'CSV Reader'
inputCSVFile = CSVReader(FileName=[InputFile])

# create a new 'Table To Structured Grid'
tableToStructuredGrid1 = TableToPoints(Input=inputCSVFile)
tableToStructuredGrid1.XColumn = 'x'
tableToStructuredGrid1.YColumn = 'y'
tableToStructuredGrid1.ZColumn = 'z'

## Evaluate the n vector
calculator1 = Calculator(Input=tableToStructuredGrid1)
calculator1.ResultArrayName = 'n'
calculator1.Function = 'nx*iHat+ny*jHat+nz*kHat'

## get data bounds
Show(calculator1)
Bounds = GetActiveSource().GetDataInformation().GetBounds()
Hide(calculator1)
if len(sys.argv)>4:
	px=sys.argv[2]
	py=sys.argv[3]
	pz=sys.argv[4]
else:
	px=Bounds[1]/2
	py=Bounds[3]/2
	pz=Bounds[5]/2

## clip the x slide ( the y and z slides will be clipped latter)
clip = Clip(Input=calculator1)
clip.ClipType = 'Box'
clip.ClipType.Scale = [1,Bounds[3],Bounds[5]]
clip.ClipType.Position = [px,0,0]

## Reescale the number and size of the cylinders 
oldNP = Bounds[3]*Bounds[5]
RescaledNP = 4000
if (oldNP > RescaledNP):
	NewSize=1.4*math.sqrt(oldNP/RescaledNP)
else:
	NewSize=0.9

## create cylinder representation
glyph1 = Glyph(Input=clip,
    GlyphType='Cylinder')
glyph1.OrientationArray = ['POINTS', 'n']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = NewSize
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'All Points'
glyph1.GlyphMode = 'Uniform Spatial Distribution'
glyph1.MaximumNumberOfSamplePoints = RescaledNP
glyph1.GlyphType.Radius = 0.25
glyph1.GlyphTransform.Rotate = [0.0, 0.0, 90.0]

## find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
SetActiveView(renderView1)
SetActiveSource(glyph1)

## Remove all 0 types
threshold1 = Threshold(Input=glyph1)
threshold1.Scalars = ['POINTS', 'pt']
threshold1.ThresholdRange = [1.0, 10.0]
threshold1.AllScalars = 1
threshold1.UseContinuousCellRange = 0


glyph1Display = Show(threshold1, renderView1)

## trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.DiffuseColor = [1.0, 0.6666666666666666, 1.0]
glyph1Display.OSPRayScaleArray = 'Normals'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 10.111634457111359
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 0.5055817228555679
glyph1Display.SetScaleArray = ['POINTS', 'Normals']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Normals']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.SelectionCellLabelFontFile = ''
glyph1Display.SelectionPointLabelFontFile = ''
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleFontFile = ''
glyph1Display.DataAxesGrid.YTitleFontFile = ''
glyph1Display.DataAxesGrid.ZTitleFontFile = ''
glyph1Display.DataAxesGrid.XLabelFontFile = ''
glyph1Display.DataAxesGrid.YLabelFontFile = ''
glyph1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleFontFile = ''
glyph1Display.PolarAxes.PolarAxisLabelFontFile = ''
glyph1Display.PolarAxes.LastRadialAxisTextFontFile = ''
glyph1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
glyph1Display.Interpolation = 'Flat'
glyph1.GlyphType.Resolution = 20

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'S'))
glyph1Display.RescaleTransferFunctionToDataRange(True, False)
sPWF = GetOpacityTransferFunction('S')
sPWF.Points = [0.4, 0.0, 0.5, 0.0, 0.4001, 1.0, 0.5, 0.0]
sPWF.ScalarRangeInitialized = 1
sPWF.RescaleTransferFunction(0.2, 0.6)
sLUT = GetColorTransferFunction('S')
sLUT.RGBPoints = [0.438, 0.705882, 0.0156863, 0.14902, 0.4381,0.8, 1.5, 1.2]
sLUT.ScalarRangeInitialized = 1.0
glyph1Display.SetScalarBarVisibility(renderView1, False)
renderView1.OrientationAxesVisibility = 0

#Printing the x slice
renderView1.CameraPosition =   [3*Bounds[1], Bounds[3]/2, Bounds[5]/2]
renderView1.CameraFocalPoint = [Bounds[1]/2, Bounds[3]/2, Bounds[5]/2]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelProjection = 1
renderView1.CameraParallelScale =  Bounds[5]/1.95
renderView1.Background=[0.6,0.6,0.6]
ratio=int(Bounds[3]/Bounds[5]+0.5)
# save screenshot
SaveScreenshot("x_"+OutputFile, renderView1, SaveAllViews=0,
    ImageResolution=[2052*ratio, 2052],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=0,
    SeparatorColor=[0.0, 0.0, 0.0],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0,
    ImageQuality=100)

clip.ClipType.Scale = [Bounds[1],1,Bounds[5]]
clip.ClipType.Position = [0,py,0]
oldNP = Bounds[3]*Bounds[5]
if (oldNP > RescaledNP):
	NewSize=1.4*math.sqrt(oldNP/RescaledNP)
else:
	NewSize=0.9
glyph1.ScaleFactor = NewSize

#Printing the y slice
renderView1.CameraPosition =   [Bounds[1]/2, 6*Bounds[3], Bounds[5]/2]
renderView1.CameraFocalPoint = [Bounds[1]/2, Bounds[3]/2, Bounds[5]/2]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale =  Bounds[5]/1.95
renderView1.CameraParallelProjection = 1
ratio=int(Bounds[1]/Bounds[5]+0.5)
# save screenshot
SaveScreenshot("y_"+OutputFile, renderView1, SaveAllViews=0,
    ImageResolution=[2052*ratio, 2052],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=0,
    SeparatorColor=[0.0, 0.0, 0.0],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0,
    ImageQuality=100)

clip.ClipType.Scale = [Bounds[1],Bounds[3],1]
clip.ClipType.Position = [0,0,pz]
oldNP = Bounds[3]*Bounds[5]
if (oldNP > RescaledNP):
	NewSize=1.4*math.sqrt(oldNP/RescaledNP)
else:
	NewSize=0.9
glyph1.ScaleFactor = NewSize

#Printing the z slice
renderView1.CameraPosition =   [Bounds[1]/2, Bounds[3]/2, 6*Bounds[5]]
renderView1.CameraFocalPoint = [Bounds[1]/2, Bounds[3]/2, Bounds[5]/2]
renderView1.CameraViewUp = [0.0, 1.0, 0.0]
renderView1.CameraParallelScale =  Bounds[3]/1.95
renderView1.CameraParallelProjection = 1
ratio=int(Bounds[1]/Bounds[3]+0.5)
# save screenshot
SaveScreenshot("z_"+OutputFile, renderView1, SaveAllViews=0,
    ImageResolution=[2052*ratio, 2052],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=0,
    SeparatorColor=[0.0, 0.0, 0.0],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0,
    ImageQuality=100)




