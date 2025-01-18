import qupath.lib.scripting.QP
import qupath.lib.roi.ROIs
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane
import qupath.lib.objects.PathDetectionObject
import qupath.ext.biop.cellpose.Cellpose2D
import qupath.lib.images.ImageData
import qupath.lib.measurements.MeasurementList
import qupath.lib.objects.PathCellObject
import qupath.lib.objects.PathObject
import qupath.lib.roi.RoiTools
import qupath.lib.roi.interfaces.ROI
import qupath.lib.gui.viewer.overlays.HierarchyOverlay
import qupath.lib.gui.images.servers.RenderedImageServer
import qupath.lib.regions.RegionRequest
import qupath.imagej.tools.IJTools


//basidc file-----------------------------
// Get project and image information
def project = getProject()
def projectName = project.getName().split('/')[0]

def imageName = getProjectEntry().getImageName()
def sample = imageName.split('_iWAT')[0]

// Define output directory (you can modify this path as needed)
def outputBaseDir = "F:/Scripps Research Dropbox/Yu Wang/A. LAB/B. Data/Yu Revision_Saba/prestain/qupath_export"  // Change this to your desired output location
def outputDir = buildFilePath(outputBaseDir, projectName)

// Create output directories
def dropletDir = buildFilePath(outputDir, sample, "droplets")
def tissueDir = buildFilePath(outputDir, sample, "tissues")
def screenshotDir = buildFilePath(outputDir, "screenshot")

// Ensure directories exist
mkdirs(dropletDir)
mkdirs(tissueDir)
mkdirs(screenshotDir)

//-----------------------------------------------------



// Clear existing objects
clearAllObjects()
runPlugin('qupath.imagej.detect.tissue.SimpleTissueDetection2', 
    '{"threshold": 50, "requestedPixelSizeMicrons": 4, "darkBackground": true, "dilateBoundaries": true, "singleAnnotation": true}')
selectAnnotations()
def imageData = getCurrentImageData()
def pathObjects = getSelectedObjects()

// Step 2: Run Cellpose detection

def cellpose = Cellpose2D.builder('cyto3_denoise')
    .pixelSize(0.6)
    .channels('Channel 1')
    .diameter(0)
    .measureShape()
    .measureIntensity()
    .classify("Droplets")
    .flowThreshold( 2.3 ) 
    .cellprobThreshold(1.7)
    .excludeEdges()
    .simplify( 0 ) 
    .build()

cellpose.detectObjects(imageData, pathObjects)

// Helper function for creating center ROI
def createCenterROI(centroidX, centroidY, centerDiameter) {
    try {
        return ROIs.createEllipseROI(
            centroidX - (centerDiameter/2),
            centroidY - (centerDiameter/2),
            centerDiameter,
            centerDiameter,
            ImagePlane.getDefaultPlane()
        )
    } catch (Exception e) {
        print("Error creating ROI: ${e.getMessage()}")
        return null
    }
}
// Step 3: First round of filtering
def firstRoundFilter() {
    def detections = getDetectionObjects()
    if (detections == null || detections.isEmpty()) {
        print("No detections found!")
        return
    }

    detections.each { droplet ->
        try {
            def measurements = droplet.getMeasurementList()
            def circularity = measurements.getMeasurementValue("Circularity")
            def area = measurements.getMeasurementValue("Area µm^2")
            
            // First round filtering criteria
            if (circularity < 0.3 || area > 10000 || area < 150) {
                droplet.setPathClass(null)
            }
        } catch (Exception e) {
            print("Error in first round filtering: ${e.getMessage()}")
        }
    }
}

// Step 4: Second round of filtering with center ROI measurements
def secondRoundFilter() {
    def detections = getDetectionObjects().findAll { it.getPathClass()?.getName() == "Droplets" }
    
    detections.each { droplet ->
        try {
            def measurements = droplet.getMeasurementList()
            def roi = droplet.getROI()
            
            // Get measurements with null checks
            def minDiameterUm = measurements.getMeasurementValue("Min diameter µm")
            if (minDiameterUm == null) {
                print("Warning: Min diameter measurement missing for droplet")
                return
            }
            
            // Get coordinates and create center ROI
            def centroidX = roi.getCentroidX()
            def centroidY = roi.getCentroidY()
            def centerDiameter = 0.7 * minDiameterUm
            
            def centerROI = createCenterROI(centroidX, centroidY, centerDiameter)
            if (centerROI == null) {
                print("Warning: Failed to create center ROI for droplet")
                return
            }
            
            // Measure center ROI intensity
            def centerAnnotation = PathObjects.createAnnotationObject(centerROI)
            addObject(centerAnnotation)
            selectObjects(centerAnnotation)
            runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', 
                '{"channel1": true, "doMean": true, "doStdDev": false, "doMinMax": false, "doMedian": false, "doHaralick": false}')
            
            def centerMeanIntensity = centerAnnotation.getMeasurementList().getMeasurementValue("ROI: 2.00 µm per pixel: Channel 1: Mean")
            def channelMean = measurements.getMeasurementValue("Channel 1: Mean")
            
            // Store center intensity measurement
            if (centerMeanIntensity != null) {
                measurements.putMeasurement("center_mean_intensity", centerMeanIntensity)
            }
            
            // Second round filtering criteria
            if (channelMean > 25000 && centerMeanIntensity != null && centerMeanIntensity > channelMean) {
                droplet.setPathClass(null)
            }
            
            removeObject(centerAnnotation, true)
            
        } catch (Exception e) {
            print("Error in second round filtering: ${e.getMessage()}")
        }
    }
}

// Run both rounds of filtering
firstRoundFilter()
secondRoundFilter()

// Step 5: Calculate area measurements for annotations
getAnnotationObjects().each { annotation ->
    def children = annotation.getChildObjects()
    def totalDropletArea = 0.0
    def dropletCount = 0
    
    children.each { child ->
        if (child.getPathClass()?.getName() == "Droplets") {
            totalDropletArea += child.getROI().getArea()
            dropletCount++
        }
    }
    
    def annotationArea = annotation.getROI().getArea()
    def areaPercentage = (totalDropletArea / annotationArea) * 100
    
    def measurements = annotation.getMeasurementList()
    measurements.putMeasurement("Total droplet area µm²", totalDropletArea)
    measurements.putMeasurement("Droplet count", dropletCount)
    measurements.putMeasurement("Droplet area percentage", areaPercentage)
}

def dropletFile = new File(dropletDir, imageName + "_droplet_measurements.csv")
saveDetectionMeasurements(dropletFile.getAbsolutePath())

def annotationFile = new File(tissueDir, imageName + "_annotation_measurements.csv")
saveAnnotationMeasurements(annotationFile.getAbsolutePath())

// Save screenshot
getCurrentHierarchy().getSelectionModel().clearSelection()
double downsample = 5
String screenshotPath = buildFilePath(screenshotDir, imageName + "_segment_rendered.png")
def viewer = getCurrentViewer()

def display = qupath.lib.display.ImageDisplay.create(imageData)
def server = new RenderedImageServer.Builder(imageData)
    .display(display)
    .downsamples(downsample)
    .layers(new HierarchyOverlay(viewer.getImageRegionStore(), viewer.getOverlayOptions(), imageData))
    .build()

// Write the rendered image
mkdirs(new File(screenshotPath).getParent())
writeImage(server, screenshotPath)
imageData.getServer().close()
print("Processing complete!")
print("Measurements exported to: ${outputDir}")
print("Screenshot saved as: ${screenshotPath}")
