import os
import unittest
from __main__ import vtk, qt, ctk, slicer
import numpy as np



#
# CTP
#

class CTP:
  def __init__(self, parent):
    parent.title = "CT Perfusion" 
    parent.categories = ["Quantification"]
    parent.dependencies = []
    parent.contributors = ["Erol Cimen (POW)", "Ken Butcher (POW)"] # replace with "Firstname Lastname (Org)"
    parent.helpText = """
    This will perform mathematical processing on CT perfusion data.
    """
    parent.acknowledgementText = """This was originally created by Erol Cimen, POW""" # replace with organization, grant and thanks.
    self.parent = parent

    # TODO: 
    # Add this test to the SelfTest module's list for discovery when the module
    # is created.  Since this module may be discovered before SelfTests itself,
    # create the list if it doesn't already exist.

    try:
      slicer.selfTests
    except AttributeError:
      slicer.selfTests = {}
    slicer.selfTests['CTP'] = self.runTest

  def runTest(self):
    tester = CTPTest()
    tester.runTest()






#
# qCTPWidget 
#
# THIS CLASS SETS UP THE LAYOUT FOR THE MODULE, (USING QT)


class CTPWidget:


  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parent.show()

    self.widgets = []
    self.inputs = [] 
    self.widgetConnections = []

  def __del__(self): 
    self.widgetConnections = []
    self.widgets = []


  def setup(self):
    # Instantiate and connect widgets ...


    # make an instance of the logic for use by the slots
    self.logic = CTPLogic()



    #
    # RELOAD AND TEST AREA 
    #
    reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    reloadCollapsibleButton.text = "Debugging"
    self.layout.addWidget(reloadCollapsibleButton)
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)

    # RELOAD BUTTON 
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "CTP Reload"
    reloadFormLayout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)


    #
    #  IMPORT AREA 
    # 
    # 
    #  CT Perfusion - Import 4D 'Volume'
    #

    importCollapsibleButton = ctk.ctkCollapsibleButton()
    importCollapsibleButton.text = "Import CT Perfusion Data"
    self.layout.addWidget(importCollapsibleButton)
    importFormLayout = qt.QFormLayout(importCollapsibleButton)

    



    #testing 
    """
    self.importButton = qt.QPushButton("Get Data")
    self.importButton.toolTip = "Grab some data."
    self.importButton.name = "CTPData"
    importFormLayout.addWidget(self.importButton)
    self.importButton.connect('clicked()', self.logic.getInfo)
    """


    self.inputWidget = self.createInputWidget(0, nodeTypes = ["vtkMRMLMultiVolumeNode"]) #4D volume will be input '0' in self.inputs !
    name = "Input 4D CTP 'Volume': "
    inputSelectorLabel = qt.QLabel(name)
    self.widgets.append(inputSelectorLabel)
    importFormLayout.addRow(inputSelectorLabel, self.inputWidget)

    self.inputs.append(self.inputWidget.currentNode())
    print(self.inputs, "length: ", len(self.inputs))



    self.inputAIFwidget = self.createInputWidget(1, nodeTypes = ["vtkMRMLLabelMapVolumeNode"]) #AIF will be input '1' in self.inputs !
    AIF_label = "Input AIF (labelmapvolume): "
    inputSelectorLabel = qt.QLabel(AIF_label)
    self.widgets.append(inputSelectorLabel)
    importFormLayout.addRow(inputSelectorLabel, self.inputAIFwidget)

    self.inputs.append(self.inputAIFwidget.currentNode())




    #TESTING PLOTS  
    """
    self.plotButton = qt.QPushButton("Plot Data")
    self.plotButton.toolTip = "Plot some data."
    self.plotButton.name = "CTPData"
    importFormLayout.addWidget(self.plotButton)
    self.plotButton.connect('clicked()', self.logic.plottingDensity)
    """


    #
    # DO THE CALCULATIONS 
    #

    calculateCollapsibleButton = ctk.ctkCollapsibleButton()
    calculateCollapsibleButton.text = "Perform Calculations"
    self.layout.addWidget(calculateCollapsibleButton)
    # Layout within the scrolling collapsible button
    calculateFormLayout = qt.QFormLayout(calculateCollapsibleButton)


    #noncerebral tissue (sulci, arteries, veins) can be removed from CTP maps by image segmentation, using thresholds based on attenuation measurements 
    #Gray matter usually around 30 - 40 HU, white matter ~ 20 - 30 HU
    #remove all pixels < 0 HU or > 60-80 HU effectively eliminates bone, fat, and air from CT image

    #should form a segment, rather than "eliminate"

    """
    self.eliminateButton = qt.QPushButton("Eliminate Noncerebral Tissue")
    self.eliminateButton.toolTip = "Remove all noncerebral tissue (sulci, arteries, veins) voxels from the DICOM pixel array."
    self.eliminateButton.name = "CTP Eliminate"
    calculateFormLayout.addRow(self.eliminateButton)
    self.eliminateButton.connect('clicked()', self.logic.calculateRelativeValues)
    """


    self.mapsButton = qt.QPushButton("GenerateMaps")
    self.mapsButton.toolTip = "Generate CTP maps."
    self.mapsButton.name = "CTP Eliminate"
    calculateFormLayout.addRow(self.mapsButton)
    self.mapsButton.connect('clicked()', lambda input1=self.inputs[0], input2=self.inputs[1]: self.logic.generateMaps(input1, input2))




  def createInputWidget(self, n, noneEnabled=True, nodeTypes = ["vtkMRMLScalarVolumeNode", "vtkMRMLMultiVolumeNode", "vtkMRMLLabelMapVolumeNode"]):
    inputSelector = slicer.qMRMLNodeComboBox()
    self.widgets.append(inputSelector)
    inputSelector.nodeTypes = nodeTypes
    inputSelector.selectNodeUponCreation = False 
    inputSelector.addEnabled = False
    inputSelector.removeEnabled = False
    inputSelector.noneEnabled = noneEnabled
    inputSelector.showHidden = False
    inputSelector.showChildNodeTypes = False
    inputSelector.setMRMLScene( slicer.mrmlScene )
    inputSelector.setToolTip( "Pick the input to the algorithm." )

    # connect and verify parameters
    # inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", lambda node,i=n:self.onInputSelect(node,i))
    inputSelector.connect("nodeActivated(vtkMRMLNode*)", lambda node,i=n:self.onInputSelect(node,i))
    self.widgetConnections.append((inputSelector, "nodeActivated(vtkMRMLNode*)"))
    return inputSelector

  def onInputSelect(self, mrmlNode, n):
    # print("The length of inputs is", len(self.inputs))
    print("works", n, mrmlNode.GetID())
    # print(self.inputs[0].GetID(), self.inputs[1].GetID())
    self.inputs[n] = mrmlNode

    self.inputWidget.setMRMLScene(slicer.mrmlScene)
    self.inputAIFwidget.setMRMLScene(slicer.mrmlScene)

    self.mapsButton.connect('clicked()', lambda input1=self.inputs[0], input2=self.inputs[1]: self.logic.generateMaps(input1, input2))



  def check(self):
    print("was triggered")

  def cleanup(self):

    pass




  
  # Generic reload method for any scripted module.

  def onReload(self,moduleName="CTP"):

    import imp, sys, os, slicer

    slicer.app.pythonConsole().clear()

    widgetName = moduleName + "Widget"

    # reload the source code
    # - set source file path
    # - load the module to the global space
    filePath = eval('slicer.modules.%s.path' % moduleName.lower())
    p = os.path.dirname(filePath)
    if not sys.path.__contains__(p):
      sys.path.insert(0,p)
    fp = open(filePath, "r")
    globals()[moduleName] = imp.load_module(
        moduleName, fp, filePath, ('.py', 'r', imp.PY_SOURCE))
    fp.close()

    # rebuild the widget
    # - find and hide the existing widget
    # - create a new widget in the existing parent
    parent = slicer.util.findChildren(name='%s Reload' % moduleName)[0].parent().parent()
    for child in parent.children():
      try:
        child.hide()
      except AttributeError:
        pass
    # Remove spacer items
    item = parent.layout().itemAt(0)
    while item:
      parent.layout().removeItem(item)
      item = parent.layout().itemAt(0)

    # delete the old widget instance
    if hasattr(globals()['slicer'].modules, widgetName):
      getattr(globals()['slicer'].modules, widgetName).cleanup()

    # create new widget inside existing parent
    globals()[widgetName.lower()] = eval(
        'globals()["%s"].%s(parent)' % (moduleName, widgetName))
    globals()[widgetName.lower()].setup()
    setattr(globals()['slicer'].modules, widgetName, globals()[widgetName.lower()])








#
# CTPLogic
#

"""
This class should implement all the actual 
  computation done by your module.  The interface 
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget
"""



class CTPLogic:

  def __init__(self):
    pass

  def setLevels(self, map):

    displayNode = map.GetDisplayNode()
    displayNode.AutoWindowLevelOff()

    map_array = slicer.util.arrayFromVolume(map)

    #grab a slice and average it ! 

    num_axes = int(len(map_array)/2)
    width = int(len(map_array[0])/2)
    height = int(len(map_array[0][0])/2)

    start = int(height - 10)
    end = int(height + 10)

    random_slice = map_array[num_axes][width][start:end]
    print(random_slice)

    mean_random_slice = np.mean(random_slice)
    max_random_slice = np.amax(random_slice)


    displayNode.SetWindow(max_random_slice)
    displayNode.SetLevel(mean_random_slice)

    displayNode.AutoThresholdOff()
    displayNode.SetLowerThreshold(0.1)
    displayNode.ApplyThresholdOn()

    displayNode.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileColdToHotRainbow.txt')

    return map



  def generateMaps(self, inputVolume, inputAIF):

    try:
      print(inputVolume.GetID(), inputAIF.GetID())
    except:
      print("You haven't selected volumes!")
      print(inputAIF)

      #convert to segment then to labelmapvolume again? 

      inputAIF_array = slicer.util.arrayFromVolume(inputAIF)
      print(inputAIF_array.values())
      return

    """
    /home/billg/Work/Slicer-SuperBuild-Debug/CTP_Analysis-build/lib/Slicer-4.13/cli-modules/DSCMRIAnalysis --relaxivity 0.0039 --S0grad 15.0 --fTolerance 1e-4 --gTolerance 1e-4 --xTolerance 1e-5 --epsilon 1e-9 --maxIter 200 --aucTimeInterval 90 --aifMask /tmp/Slicer-billg/BGJBB_vtkMRMLLabelMapVolumeNodeC.nrrd --outputAUC /tmp/Slicer-billg/BGJBB_vtkMRMLScalarVolumeNodeB.nrrd --outputCBF /tmp/Slicer-billg/BGJBB_vtkMRMLScalarVolumeNodeC.nrrd --outputMTT /tmp/Slicer-billg/BGJBB_vtkMRMLScalarVolumeNodeD.nrrd --outputK2 /tmp/Slicer-billg/BGJBB_vtkMRMLScalarVolumeNodeE.nrrd /tmp/Slicer-billg/BGJBB_vtkMRMLMultiVolumeNodeB.nrrd
    """

    cliModule = slicer.modules.dscmrianalysis

    n=cliModule.cliModuleLogic().CreateNode()

    for groupIndex in range(n.GetNumberOfParameterGroups()):
      print(f'Group: {n.GetParameterGroupLabel(groupIndex)}')
      for parameterIndex in range(n.GetNumberOfParametersInGroup(groupIndex)):
        print('  {0} [{1}]: {2}'.format(n.GetParameterName(groupIndex, parameterIndex),
          n.GetParameterTag(groupIndex, parameterIndex),n.GetParameterLabel(groupIndex, parameterIndex)))


    # Set parameters
    parameters = {}
    parameters["InputFourDImageFileName"] = inputVolume
    parameters["AIFMaskFileName"] = inputAIF
    outputCBV = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "Output CBV")
    outputCBF = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "Output CBF")
    outputMTT = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "Output MTT")
    outputTTP = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "Output TTP")
    
    parameters["OutputAUCFileName"] = outputCBV
    parameters["OutputCBFFileName"] = outputCBF
    parameters["OutputMTTFileName"] = outputMTT
    parameters["OutputK2FileName"] = outputTTP


    outputMaps = []

    outputMaps.append(outputCBV)
    outputMaps.append(outputCBF)
    outputMaps.append(outputMTT)
    outputMaps.append(outputTTP)
    
    # Execute
    ctPerfusion = slicer.modules.dscmrianalysis
    cliNode = slicer.cli.runSync(ctPerfusion, None, parameters)

    # Process results
    if cliNode.GetStatus() & cliNode.ErrorsMask:
      # error
      errorText = cliNode.GetErrorText()
      slicer.mrmlScene.RemoveNode(cliNode)
      raise ValueError("CLI execution failed: " + errorText)

    # success
    slicer.mrmlScene.RemoveNode(cliNode)

    for map in outputMaps:
      self.setLevels(map)
  


    #plot the AIF as well!

    return outputTTP



  def getInfo(self):

    #import the CT DICOM files as a Multi-Volume Sequence 
    #either in 'advanced, examine --> select it' 
    #or go to application settings, DICOM, change MultiVolumeImporter default to 'Sequence'

    sequences = slicer.util.getNodesByClass('vtkMRMLSequenceNode')
    sequenceNode = sequences[0]

    print("Number of sequences: {}".format(len(sequences)))
    print("Names of those sequence nodes {}".format(sequences))


    # voxels = slicer.util.arrayFromVolume(node1) 
    #voxels = slicer.util.arrayFromVolume(sequenceNode.GetNthDataNode(13))
    #print(voxels) #len(voxels) = num slices, #len(voxels[0]) = num values 

    sequenceNode = sequenceNode.GetNthDataNode(0)

    def onMouseMoved(observer,eventid):  
      point_Ras=[0,0,0]
      crosshairNode.GetCursorPositionRAS(point_Ras)

      #point RAS coordinates 
      print("point RAS {}".format(point_Ras))
      

      # Get voxel coordinates from physical coordinates
      volumeRasToIjk = vtk.vtkMatrix4x4()
      sequenceNode.GetRASToIJKMatrix(volumeRasToIjk)
      point_Ijk = [0, 0, 0, 1]
      volumeRasToIjk.MultiplyPoint(np.append(point_Ras,1.0), point_Ijk)
      point_Ijk = [ int(round(c)) for c in point_Ijk[0:3] ]

      # point IJK coordinates!
      print("point IJK {}".format(point_Ijk))

      #value at those coordinates!
      voxels = slicer.util.arrayFromVolume(sequenceNode)
      voxelValue = voxels[point_Ijk[2], point_Ijk[1], point_Ijk[0]]  # note that numpy array index order is kji (not ijk)
      print(voxelValue)

      self.moving_plottingDensity(point_Ijk)

    crosshairNode=slicer.util.getNode('Crosshair') 
    crosshairNode.AddObserver(slicer.vtkMRMLCrosshairNode.CursorPositionModifiedEvent, onMouseMoved)


    return 



  def moving_plottingDensity(self, point_Ijk):
    
    #you can grab the cursor position and constantly update this table
    #that would be cool!

    # Create table with time and density (HU) columns 

    if (slicer.mrmlScene.GetNodesByClass("vtkMRMLTableNode").GetNumberOfItems() == 0):
      tableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
      table = tableNode.GetTable()

      arrX = vtk.vtkFloatArray()
      arrX.SetName("Time")
      table.AddColumn(arrX)

      arrY = vtk.vtkFloatArray()
      arrY.SetName("Density")
      table.AddColumn(arrY)
    else:
      tableNode = slicer.mrmlScene.GetNodeByID('vtkMRMLTableNode1')
      table = slicer.mrmlScene.GetNodeByID('vtkMRMLTableNode1').GetTable()



    # Fill in the table with values

    
    numPoints = slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNumberOfDataNodes() #TODO: replace with selector

    table.SetNumberOfRows(numPoints)


    first_instUids = slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNthDataNode(0).GetAttribute('DICOM.instanceUIDs').split()

    first_filename = slicer.dicomDatabase.fileForInstance(first_instUids[0])

    baseline = float(slicer.dicomDatabase.fileValue(first_filename,'0008,0032')) #(group, element) tag for "AcquisitionTime" in metadata


    for i in range(numPoints):

      instUids = slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNthDataNode(i).GetAttribute('DICOM.instanceUIDs').split()
      filename = slicer.dicomDatabase.fileForInstance(instUids[0])

      acquisition_time = float(slicer.dicomDatabase.fileValue(filename,'0008,0032')) - baseline

      table.SetValue(i, 0, acquisition_time)

      voxels = slicer.util.arrayFromVolume(slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNthDataNode(i))
      voxelValue = voxels[point_Ijk[2], point_Ijk[1], point_Ijk[0]]

      table.SetValue(i, 1, float(voxelValue))
      #for some reason doesn't like ints, has to be passed a float


    # Create a plot series nodes
    
   
    if (slicer.mrmlScene.GetNodesByClass("vtkMRMLPlotSeriesNode").GetNumberOfItems() == 0):
      plotSeriesNode1 = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", "ROI")
    else:

      while (slicer.mrmlScene.GetNodesByClass("vtkMRMLPlotSeriesNode").GetNumberOfItems() != 0):
        slicer.mrmlScene.RemoveNode(slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLPlotSeriesNode'))

      plotSeriesNode1 = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", "ROI")


    plotSeriesNode1.SetAndObserveTableNodeID(tableNode.GetID())
    plotSeriesNode1.SetXColumnName("Time")
    plotSeriesNode1.SetYColumnName("Density")
    plotSeriesNode1.SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatter)
    plotSeriesNode1.SetLineStyle(slicer.vtkMRMLPlotSeriesNode.LineStyleSolid)
    plotSeriesNode1.SetMarkerStyle(slicer.vtkMRMLPlotSeriesNode.MarkerStyleSquare)
    #plotSeriesNode1.SetUniqueColor()




    # Create plot chart node

    if (slicer.mrmlScene.GetNodesByClass("vtkMRMLPlotChartNode").GetNumberOfItems() == 0):

      plotChartNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode")

    else:

      while (slicer.mrmlScene.GetNodesByClass("vtkMRMLPlotChartNode").GetNumberOfItems() != 0):
        slicer.mrmlScene.RemoveNode(slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLPlotChartNode'))

      plotChartNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode")

    plotChartNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode")
    plotChartNode.AddAndObservePlotSeriesNodeID(plotSeriesNode1.GetID())
    plotChartNode.SetTitle('Time-Density Curve')
    plotChartNode.SetXAxisTitle('Time (s)')
    plotChartNode.SetYAxisTitle('Density (HU)')




    # Switch to a layout that contains a plot view to create a plot widget

    layoutManager = slicer.app.layoutManager()
    layoutWithPlot = slicer.modules.plots.logic().GetLayoutWithPlot(layoutManager.layout)
    layoutManager.setLayout(layoutWithPlot)

    # Select chart in plot view

    plotWidget = layoutManager.plotWidget(0)
    plotViewNode = plotWidget.mrmlPlotViewNode()
    plotViewNode.SetPlotChartNodeID(plotChartNode.GetID())


    return 



  def plottingDensity(self, Xcoord = 250, Ycoord = 250):
    
    #you can grab the cursor position and constantly update this table
    #that would be cool!

    # Create table with time and density (HU) columns 

    """
    tableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
    table = tableNode.GetTable()

    arrX = vtk.vtkFloatArray()
    arrX.SetName("Time")
    table.AddColumn(arrX)

    arrY = vtk.vtkFloatArray()
    arrY.SetName("Density")
    table.AddColumn(arrY)
    """

    #limit to a single table (only 1 new table node added to the mrml scene)
    if (slicer.mrmlScene.GetNodesByClass("vtkMRMLTableNode").GetNumberOfItems() == 0):
      tableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
      table = tableNode.GetTable()

      arrX = vtk.vtkFloatArray()
      arrX.SetName("Time")
      table.AddColumn(arrX)

      arrY = vtk.vtkFloatArray()
      arrY.SetName("Density")
      table.AddColumn(arrY)
    else:
      tableNode = slicer.mrmlScene.GetNodeByID('vtkMRMLTableNode1')
      table = slicer.mrmlScene.GetNodeByID('vtkMRMLTableNode1').GetTable()



    # Fill in the table with values

    axial_slice = 12
    
    numPoints = slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNumberOfDataNodes() #TODO: replace with selector

    table.SetNumberOfRows(numPoints)


    first_instUids = slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNthDataNode(0).GetAttribute('DICOM.instanceUIDs').split()

    first_filename = slicer.dicomDatabase.fileForInstance(first_instUids[0])

    baseline = float(slicer.dicomDatabase.fileValue(first_filename,'0008,0032')) #(group, element) tag for "AcquisitionTime" in metadata


    for i in range(numPoints):

      instUids = slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNthDataNode(i).GetAttribute('DICOM.instanceUIDs').split()
      filename = slicer.dicomDatabase.fileForInstance(instUids[0])

      acquisition_time = float(slicer.dicomDatabase.fileValue(filename,'0008,0032')) - baseline

      table.SetValue(i, 0, acquisition_time)
      table.SetValue(i, 1, float(slicer.util.arrayFromVolume(slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNthDataNode(i))[axial_slice][Xcoord][Ycoord]))
      #for some reason doesn't like ints, has to be passed a float


    # Create a plot series node


    if (slicer.mrmlScene.GetNodesByClass("vtkMRMLPlotSeriesNode").GetNumberOfItems() == 0):
      plotSeriesNode1 = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", "ROI")
    else:
      plotSeriesNode1 = slicer.mrmlScene.GetNodeByID('vtkMRMLPlotSeriesNode1')

    plotSeriesNode1.SetAndObserveTableNodeID(tableNode.GetID())
    plotSeriesNode1.SetXColumnName("Time")
    plotSeriesNode1.SetYColumnName("Density")
    plotSeriesNode1.SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatter)
    plotSeriesNode1.SetLineStyle(slicer.vtkMRMLPlotSeriesNode.LineStyleSolid)
    plotSeriesNode1.SetMarkerStyle(slicer.vtkMRMLPlotSeriesNode.MarkerStyleSquare)
    #plotSeriesNode1.SetUniqueColor()


    # Create plot chart node

    if (slicer.mrmlScene.GetNodesByClass("vtkMRMLPlotChartNode").GetNumberOfItems() == 0):
      plotChartNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode")
    else:
      plotChartNode = slicer.mrmlScene.GetNodeByID('vtkMRMLPlotChartNode1')

    plotChartNode.AddAndObservePlotSeriesNodeID(plotSeriesNode1.GetID())
    plotChartNode.SetTitle('Time-Density Curve')
    plotChartNode.SetXAxisTitle('Time (s)')
    plotChartNode.SetYAxisTitle('Density (HU)')




    # Switch to a layout that contains a plot view to create a plot widget

    layoutManager = slicer.app.layoutManager()
    layoutWithPlot = slicer.modules.plots.logic().GetLayoutWithPlot(layoutManager.layout)
    layoutManager.setLayout(layoutWithPlot)

    # Select chart in plot view

    plotWidget = layoutManager.plotWidget(0)
    plotViewNode = plotWidget.mrmlPlotViewNode()
    plotViewNode.SetPlotChartNodeID(plotChartNode.GetID())


    return 



  
  def calculateRelativeValues(self):

    #Take the time density curves FOR EACH VOXEL 

    
    time_to_peak = np.nan #time to peak 
    relative_CBF = np.nan #peak of curve
    relative_CBV = np.nan #Area under the curve 

    sequenceNode = slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0]

    numPoints = sequenceNode.GetNumberOfDataNodes() #TODO: replace with selector


    first_instUids = sequenceNode.GetNthDataNode(0).GetAttribute('DICOM.instanceUIDs').split()

    first_filename = slicer.dicomDatabase.fileForInstance(first_instUids[0])

    baseline = float(slicer.dicomDatabase.fileValue(first_filename,'0008,0032')) #(group, element) tag for "AcquisitionTime" in metadata

    num_slices = len(slicer.util.arrayFromVolume(sequenceNode.GetNthDataNode(0)))
    num_rows = len(slicer.util.arrayFromVolume(sequenceNode.GetNthDataNode(0))[0])
    num_columns = len(slicer.util.arrayFromVolume(sequenceNode.GetNthDataNode(0))[0][0])

    print(num_slices)
    print(num_rows)
    print(num_columns)

    for x in range(num_slices):
      for y in range(num_rows):
        for z in range(num_columns):

          densities = []
          acquisition_times = []

          for i in range(numPoints):

            instUids = slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNthDataNode(i).GetAttribute('DICOM.instanceUIDs').split()
            filename = slicer.dicomDatabase.fileForInstance(instUids[0])

            acquisition_time = float(slicer.dicomDatabase.fileValue(filename,'0008,0032')) - baseline

            acquisition_times.append(acquisition_time)
            densities.append(float(slicer.util.arrayFromVolume(slicer.util.getNodesByClass('vtkMRMLSequenceNode')[0].GetNthDataNode(i))[x][y][z]))
            #for some reason doesn't like ints, has to be passed a float
          
          #print(densities, acquisition_times)


    """
    #create new volume 

    imageData = vtk.vtkImageData()
    imageData.SetDimensions(24, 512, 512) #slicer.util.arrayFromVolume(slicer.mrmlScene.GetNodeByID('vtkMRMLScalarVolumeNode29')).shape
    imageData.AllocateScalars(vtk.VTK_SHORT, 1)

    volumeID = "timetopeak"
  
    node = slicer.vtkMRMLScalarVolumeNode()
    node.SetName(volumeID)
    slicer.mrmlScene.AddNode(node)
    node.CreateDefaultDisplayNodes()

    node.SetAndObserveImageData(imageData)
    #node.SetIJKToRASMatrix(ijkToRAS)

    pixels = np.array([[1, 2, 3, 5], [2, 2, 2, 5]])
    array = slicer.util.array(node.GetID())
    array[:] = pixels.reshape(array.shape)
    imageData.GetPointData().GetScalars().Modified()

    displayNode = node.GetDisplayNode()
    displayNode.ProcessMRMLEvents(displayNode, vtk.vtkCommand.ModifiedEvent, "")

    """

    print("WORKED")
    return 

  def volumeCount(self):
    return len(slicer.util.getNodes('vtkMRML*VolumeNode*'))

  def selectVolume(self,index):
    nodes = slicer.util.getNodes('vtkMRML*VolumeNode*')
    names = nodes.keys()
    names = sorted(names)
    selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    selectionNode.SetReferenceActiveVolumeID( nodes[names[index]].GetID() )
    slicer.app.applicationLogic().PropagateVolumeSelection(0)








# TODO: SET UP TESTS 

"""
  This is the test case for your scripted module.
"""


class CTPTest(unittest.TestCase):

  def delayDisplay(self,message,msec=1000):
    """This utility method displays a small dialog and waits.
    This does two things: 1) it lets the event loop catch up
    to the state of the test so that rendering and widget updates
    have all taken place before the test continues and 2) it
    shows the user/developer/tester the state of the test
    so that we'll know when it breaks.
    """
    print(message)
    self.info = qt.QDialog()
    self.infoLayout = qt.QVBoxLayout()
    self.info.setLayout(self.infoLayout)
    self.label = qt.QLabel(message,self.info)
    self.infoLayout.addWidget(self.label)
    qt.QTimer.singleShot(msec, self.info.close)
    self.info.exec_()

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_CTP1()

  def test_CTP1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests sould exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        print('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        print('Loading %s...\n' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading\n')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = CTPLogic()
    volumesLogic = slicer.modules.volumes.logic()

    blurLevelCount = 10
    for sigma in range(blurLevelCount):
      self.delayDisplay('Making blurred volume with sigma of %d\n' % sigma)
      outputVolume = volumesLogic.CloneVolume(slicer.mrmlScene, volumeNode, 'blur-%d' % sigma)
      parameters = {
          "inputVolume": slicer.util.getNode('FA'),
          "outputVolume": outputVolume,
          "sigma": sigma,
          }

      blur = slicer.modules.gaussianblurimagefilter
      slicer.cli.run(blur, None, parameters, wait_for_completion=True)

    slicer.modules.CTPWidget.onRefresh()
    self.delayDisplay('Selecting original volume')
    slicer.modules.CTPWidget.slider.value = 0
    self.delayDisplay('Selecting final volume')
    slicer.modules.CTPWidget.slider.value = blurLevelCount

    selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    selectedID = selectionNode.GetActiveVolumeID()
    lastVolumeID = outputVolume.GetID()
    if selectedID != lastVolumeID:
      raise Exception("Volume ID was not selected!\nExpected %s but got %s" % (lastVolumeID, selectedID))

    self.delayDisplay('Test passed!')