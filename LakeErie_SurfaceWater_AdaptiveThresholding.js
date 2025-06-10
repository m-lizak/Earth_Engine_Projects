// GEOG 6060 Project
// Maciej Lizak
// Fall 2023 Semester


// The primary objective of this project is to leverage Google Earth Engine (GEE) for the automated generation of
// accurate surface water extent maps for a given month and year in the Lake Erie watershed using Sentinel-1
// Synthetic Aperture Radar (SAR) data. The project employs a sophisticated adaptive thresholding approach which
// should enhance the precision and reliability of surface water detection. When compared to Otsu thresholding,
// adaptive thresholding manages to provide this additional accuracy through the use of a buffered threshold sampling
// to generate a binary image for subsequent edge detection, addressing challenges associated with SAR speckle and other artifacts.
// The method further refines results by filtering defined edges based on length, which effectively eliminates small edges that could
// potentially distort histogram sampling.

// The goal of this project is to generate accurate binary classification maps depicting surface water extent in the Lake Erie watershed
// for any user specified month and year. This should offer valuable insights for the purposes water resource management and
// environmental monitoring. If successful, applying the methodology above to the entire Great Lakes & St. Lawrence basin should be
// feasible and is a realistic next step forward



//
//  DEFINE PARAMETERS
//


// Define user-specified season and year
var userSeason = 'Summer'; // Replace with user input (Fall, Winter, Spring, Summer)
var userYear = 2023;

// Define the resolution for this script (in metres, min = 10) (avoid running below 30m without exporting)
var setResolution = 30

// Specify Sentinel 1 Polarization Band of Interest ('VV', 'VH', 'HV', 'HH') - Default 'VV' for this project
var band = 'VV';



// Define parameters for the Adaptive Thresholding.

// Minimum Length of edges to be considered water edges (anything below this will be masked out as noise and not
// included in the buffer)
var edgeLength = 25;

// Buffer in meters to apply to edges - this buffered area is then sampled to obtain the Adaptive Threshold value.
var edgeBuffer = 60;

// The kernel neighbourhood size for the filter
var connectedPixels = 100;


// Define Parameters for Canny Edge Detection

// Pixels with gradient values above this threshold are confidently classified as edges
var cannyThreshold = 1; //Default = 1

// Sigma value for canny edge detection. (Sigma value for gaussian filter applied before edge detection. 0 means apply no filtering)
var cannySigma = 0.5; //Default = 1

// Pixels with gradient values below this lower threshold are immediately discarded and considered as non-edges.
var cannyLt = 0.05; // Default = 0.05


// Import Great Lakes Basin shapefile and filter to the Lake Erie Basin
var greatLakesBasin = ee.FeatureCollection('users/mlizak/GreatLakes_Basins');
var studyArea = greatLakesBasin.filter(ee.Filter.eq('merge', 'lk_huron')); // Can change study area to anywhere in Great Lakes here, choose from:
                                                                           // lk_erie,	lk_huron, lk_mich, lk_ont, lk_sup



// 
//  INITIALIZATION AND PREPROCESSING
//

// Define the start and end months based on the user-specified season
var seasonStartMonth;
var seasonEndMonth;

switch (userSeason) {
  case 'Fall':
    seasonStartMonth = 9; // September
    seasonEndMonth = 11; // November
    break;
  case 'Winter':
    seasonStartMonth = 12; // December
    seasonEndMonth = 2;  // February
    break;
  case 'Spring':
    seasonStartMonth = 3; // March
    seasonEndMonth = 5;  // May
    break;
  case 'Summer':
    seasonStartMonth = 6; // June
    seasonEndMonth = 8;  // August
    break;
  default:
    print('Invalid season input. Please choose Fall, Winter, Spring, or Summer.');
}


var startDate;
var endDate;

// Take into account the winter edge case (images are over 2 years)
if (userSeason === 'Winter' && seasonStartMonth === 12) {
  startDate = ee.Date.fromYMD(userYear - 1, seasonStartMonth, 1);
  endDate = ee.Date.fromYMD(userYear, seasonEndMonth, 1).advance(1, 'month');
} else {
  startDate = ee.Date.fromYMD(userYear, seasonStartMonth, 1);
  endDate = ee.Date.fromYMD(userYear, seasonEndMonth, 1).advance(1, 'month');
}




// Get the Sentinel-1 collection and filter by space/time for the specified season.
var s1Collection = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(studyArea)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.eq('instrumentMode', 'IW'));

// Define a common grid using a reference image.
var referenceImage = s1Collection.first();
var commonGrid = referenceImage.select(band).resample('bicubic');

// Resample the Sentinel-1 collection to 30m resolution for the band and warp to the common grid.
var s1CollectionResampled = s1Collection.map(function(image) {
  var resampledImage = image.select(band).resample('bicubic');
  var warpedImage = resampledImage
    .reproject({
      crs: commonGrid.projection(),
      scale: setResolution, // SET RESOLUTION HERE
    });
  return warpedImage;
});


// Separate ascending and descending orbit images into distinct collections.
var s1CollectionAsc = s1CollectionResampled.filter(
  ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
var s1CollectionDesc = s1CollectionResampled.filter(
  ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));

// Calculate means for for ascending and descending orbits.
var meanAsc = s1CollectionAsc.select(band).reduce(ee.Reducer.mean());
var meanDesc = s1CollectionDesc.select(band).reduce(ee.Reducer.mean());

// Add the means as bands to a new image.
var s1Mosaic = meanAsc.addBands(meanDesc);

// Calculate the overall mean across both bands.
var s1Mosaic = s1Mosaic.reduce(ee.Reducer.mean()).clip(studyArea);

// Update the band name
s1Mosaic = s1Mosaic.rename(band);



print('Projection:', commonGrid.projection());
// 
//  GLOBAL THRESHOLDING
//



// Define a reducer to calculate a histogram of values.
var histogramReducer = ee.Reducer.histogram(255, 0.1);

// Reduce all of the image values.
var globalHistogram = ee.Dictionary(
    s1Mosaic.select(band).reduceRegion({
        reducer: histogramReducer,
        geometry: studyArea.geometry(),  // Use studyArea's geometry here
        scale: setResolution,
        maxPixels: 1e10
    }).get(band)
);


// Create a function that takes in the above histogram as input, applies the algorithm, and returns a single value where Otsu’s method suggests
// breaking the histogram into two parts
function otsu(histogram) {
    // Make sure histogram is an ee.Dictionary object.
    histogram = ee.Dictionary(histogram);
    // Extract relevant values into arrays.
    var counts = ee.Array(histogram.get('histogram'));
    var means = ee.Array(histogram.get('bucketMeans'));
    // Calculate single statistics over arrays
    var size = means.length().get([0]);
    var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0])
        .get([0]);
    var mean = sum.divide(total);
    // Compute between sum of squares, where each mean partitions the data.
    var indices = ee.List.sequence(1, size);
    var bss = indices.map(function(i) {
        var aCounts = counts.slice(0, 0, i);
        var aCount = aCounts.reduce(ee.Reducer.sum(), [0])
            .get([0]);
        var aMeans = means.slice(0, 0, i);
        var aMean = aMeans.multiply(aCounts)
            .reduce(ee.Reducer.sum(), [0]).get([0])
            .divide(aCount);
        var bCount = total.subtract(aCount);
        var bMean = sum.subtract(aCount.multiply(aMean))
            .divide(bCount);
        return aCount.multiply(aMean.subtract(mean).pow(2))
            .add(
                bCount.multiply(bMean.subtract(mean).pow(2)));
    });
    // Return the mean value corresponding to the maximum BSS.
    return means.sort(bss).get([-1]);
}






// Now run otsu thresholding.
var globalThreshold = otsu(globalHistogram);
print('Global threshold value:', globalThreshold);


//
// MAKE A CHART
//

// Extract out the histogram buckets and counts per bucket.
var x = ee.List(globalHistogram.get('bucketMeans'));
var y = ee.List(globalHistogram.get('histogram'));

// Define a list of values to plot.
var dataCol = ee.Array.cat([x, y], 1).toList();

// Define the header information for data.
var columnHeader = ee.List([
    [
    {
        label: 'Backscatter',
        role: 'domain',
        type: 'number'
    },
    {
        label: 'Values',
        role: 'data',
        type: 'number'
    }, ]
]);

// Concat the header and data for plotting.
var dataTable = columnHeader.cat(dataCol);

// Convert the server-side table to a client-side table.
var chartData = dataTable.getInfo();

// Create list of empty strings that will be used for annotation.
var thresholdCol = ee.List.repeat('', x.length());
// Find the index where the bucketMean equals the threshold.
var threshIndex = x.indexOf(globalThreshold);
// Set the index to the annotation text.
thresholdCol = thresholdCol.set(threshIndex, 'Otsu Threshold');

// Redefine the column header information with annotation column.
var columnHeader = ee.List([
    [
    {
        label: 'Backscatter',
        role: 'domain',
        type: 'number'
    },
    {
        label: 'Values',
        role: 'data',
        type: 'number'
    },
    {
        label: 'Threshold',
        role: 'annotation',
        type: 'string'
    }]
]);

// Loop through the data rows and add the annotation column.
dataCol = ee.List.sequence(0, x.length().subtract(1)).map(function(
i) {
    i = ee.Number(i);
    var row = ee.List(dataCol.get(i));
    return row.add(ee.String(thresholdCol.get(i)));
});

// Concat the header and data for plotting.
dataTable = columnHeader.cat(dataCol);

// Create plot using the ui.Chart function with the dataTable.
// Use 'evaluate' to transfer the server-side table to the client.
// Define the chart and print it to the console.
dataTable.evaluate(function(dataTableClient) {
    // loop through the client-side table and set empty strings to null
    for (var i = 0; i < dataTableClient.length; i++) {
        if (dataTableClient[i][2] === '') {
            dataTableClient[i][2] = null;
        }
    }
    var chart = ui.Chart(dataTableClient)
        .setChartType('AreaChart')
        .setOptions({
            title: band +
                ' Global Threshold Histogram ' + userSeason + ' ' + userYear,
            hAxis: {
                title: 'Backscatter [dB]',
                viewWindow: {
                    min: -35,
                    max: 15
                }
            },
            vAxis: {
                title: 'Count'
            },
            annotations: {
                style: 'line'
            }
        });
    print(chart);
});



// Apply the threshold on the image to extract water.
var globalWater = s1Mosaic.select(band).lt(globalThreshold);







// //
// // ADAPTIVE THRESHOLDING
// //


// Use Global Threshold found earlier using Otsu as Initial Threshold
var initialThreshold = globalThreshold;

// Begin restraining histogram sampling
// Get preliminary water.
var binary = s1Mosaic.select(band).lt(initialThreshold)
    .rename('binary');

// Get projection information to convert buffer size to pixels.
var imageProj = commonGrid.projection();

// Get canny edges.
var canny = ee.Algorithms.CannyEdgeDetector({
    image: binary,
    threshold: cannyThreshold,
    sigma: cannySigma
});



// Process canny edges.

// Get the edges and length of edges.
var connected = canny.updateMask(canny).lt(cannyLt)
    .connectedPixelCount(connectedPixels, true);

// Mask short edges that can be noise.
var edges = connected.gte(edgeLength);

// Calculate the buffer in pixel size.
var edgeBufferPixel = ee.Number(edgeBuffer).divide(imageProj
    .nominalScale());

// Buffer the edges using a dilation
var bufferedEdges = edges.fastDistanceTransform().lt(edgeBufferPixel);

// Mask areas not within the buffer
var edgeImage = s1Mosaic.select(band).updateMask(bufferedEdges);






///
/// Make Another Chart
///

// We now have a more representative sample region defined by the buffer, with noise eliminated. 
// Now we can replot the histogram to see the difference the adaptive technique has made

// Reduce all of the image values.
var localHistogram = ee.Dictionary(
    edgeImage.reduceRegion({
        reducer: histogramReducer,
        geometry: studyArea.geometry(),
        scale: setResolution,
        maxPixels: 1e10
    }).get(band)
);

// Apply otsu thresholding.
var localThreshold = otsu(localHistogram);
print('Adaptive Threshold Value:', localThreshold);

// Extract out the histogram buckets and counts per bucket.
var x = ee.List(localHistogram.get('bucketMeans'));
var y = ee.List(localHistogram.get('histogram'));

// Define a list of values to plot.
var dataCol = ee.Array.cat([x, y], 1).toList();

// Concat the header and data for plotting.
var dataTable = columnHeader.cat(dataCol);

// Create list of empty strings that will be used for annotation.
var thresholdCol = ee.List.repeat('', x.length());
// Find the index that bucketMean equals the threshold.
var threshIndex = x.indexOf(localThreshold);
// Set the index to the annotation text.
thresholdCol = thresholdCol.set(threshIndex, 'Otsu Threshold');

// Redefine the column header information now with annotation col.
columnHeader = ee.List([
    [
    {
        label: 'Backscatter',
        role: 'domain',
        type: 'number'
    },
    {
        label: 'Values',
        role: 'data',
        type: 'number'
    },
    {
        label: 'Threshold',
        role: 'annotation',
        type: 'string'
    }]
]);

// Loop through the data rows and add the annotation col.
dataCol = ee.List.sequence(0, x.length().subtract(1)).map(function(
i) {
    i = ee.Number(i);
    var row = ee.List(dataCol.get(i));
    return row.add(ee.String(thresholdCol.get(i)));
});

// Concat the header and data for plotting.
dataTable = columnHeader.cat(dataCol);

// Create plot using the ui.Chart function with the dataTable.
// Use 'evaluate' to transfer the server-side table to the client.
// Define the chart and print it to the console.
dataTable.evaluate(function(dataTableClient) {
    // Loop through the client-side table and set empty strings to null.
    for (var i = 0; i < dataTableClient.length; i++) {
        if (dataTableClient[i][2] === '') {
            dataTableClient[i][2] = null;
        }
    }
    var chart = ui.Chart(dataTableClient)
        .setChartType('AreaChart')
        .setOptions({
            title: band +
                ' Adaptive Threshold Histogram ' + userSeason + ' ' + userYear,
            hAxis: {
                title: 'Backscatter [dB]',
                viewWindow: {
                    min: -35,
                    max: 15
                }
            },
            vAxis: {
                title: 'Count'
            },
            annotations: {
                style: 'line'
            }
        });
    print(chart);
});




// Apply the threshold on the image to extract water using the adaptive threshold.
var localWater = s1Mosaic.select(band).lt(localThreshold);






// Count the number of water pixels using global thresholding
var globalWaterCount = globalWater.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: studyArea.geometry(),
  scale: setResolution,
  maxPixels: 1e10
}).get(band);

print('Global Water Pixel Count:', globalWaterCount);

// Count the number of water pixels using adaptive thresholding
var localWaterCount = localWater.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: studyArea.geometry(),
  scale: setResolution,
  maxPixels: 1e10
}).get(band);

print('Adaptive Water Pixel Count:', localWaterCount);





  
  // /* Export images */
  
  // var exportMosaic = {
  //   image: s1Mosaic,
  //   scale: setResolution,
  //   region: studyArea,
  //   crs: 'EPSG:32617',
  //   folder: 'GEOG 6060 Exports/Base Images',
  //   fileNamePrefix: 's1Mosaic_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   description: 's1Mosaic_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   maxPixels: 1e13
  // };
  // // Start the export task
  // var exportTask = Export.image.toDrive(exportMosaic);
  
  
  
  
  // var exportGlobalThreshold = {
  //   image: globalWater,
  //   scale: setResolution,
  //   region: studyArea,
  //   crs: 'EPSG:32617',
  //   folder: 'GEOG 6060 Exports/Global Threshold Maps',
  //   fileNamePrefix: 'globalWater_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   description: 'globalWater_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   maxPixels: 1e13
  // };
  // // Start the export task
  // var exportTask = Export.image.toDrive(exportGlobalThreshold);
  
  
  
  
  
  // var exportLocalThreshold = {
  //   image: localWater,
  //   scale: setResolution,
  //   region: studyArea,
  //   crs: 'EPSG:32617',
  //   folder: 'GEOG 6060 Exports/Local Threshold Maps',
  //   fileNamePrefix: 'localWater_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   description: 'localWater_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   maxPixels: 1e13
  // };
  // // Start the export task
  // var exportTask = Export.image.toDrive(exportLocalThreshold);
  
  
  
  
  // var exportDetectedEdges = {
  //   image: edges,
  //   scale: setResolution,
  //   region: studyArea,
  //   crs: 'EPSG:32617',
  //   folder: 'GEOG 6060 Exports/Detected Water Edges',
  //   fileNamePrefix: 'edges_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   description: 'edges_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   maxPixels: 1e13
  // };
  // // Start the export task
  // var exportTask = Export.image.toDrive(exportDetectedEdges);
  
  
  
  
  // var exportBufferedEdges = {
  //   image: bufferedEdges,
  //   scale: setResolution,
  //   region: studyArea,
  //   crs: 'EPSG:32617',
  //   folder: 'GEOG 6060 Exports/Buffered Water Edges',
  //   fileNamePrefix: 'bufferedEdges_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   description: 'bufferedEdges_' + userSeason + '_' + userYear + '_' + setResolution + 'm',
  //   maxPixels: 1e13
  // };
  // // Start the export task
  // var exportTask = Export.image.toDrive(exportBufferedEdges);
  
  
  

//
// Visualization
//

// Center the map.
Map.centerObject(studyArea, 8);

// Add the Lake Erie Basin boundary to the map.
Map.addLayer(studyArea, {}, 'Lake Erie basin');

// Add the overall mean to the map.
Map.addLayer(s1Mosaic, {min: -25, max: 0}, 'Mean Image');

// Add the water image to the map and mask 0 (no-water) values.
Map.addLayer(globalWater.selfMask(),
    {
        palette: 'blue'
    },
    'Water (Global Threshold)');

// Add the detected edges and buffered edges to the map.
Map.addLayer(edges, {
    palette: 'red'
}, 'Detected Water Edges (red)');
var edgesVis = {
    palette: 'lime',
    opacity: 1
};
Map.addLayer(bufferedEdges.selfMask(), edgesVis,
    'Buffered Water Edges (green)');

// Add the water image to the map and mask 0 (no-water) values.
Map.addLayer(localWater.selfMask(),
    {
        palette: 'darkblue'
    },
    'Water (Adaptive Threshold)');