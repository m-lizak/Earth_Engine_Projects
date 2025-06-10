# Surface Water Mapping in the Great Lakes Basin using Sentinel-1 SAR and Adaptive Thresholding

**Author:** Maciej Lizak  
**Course:** GEOG 6060 – Fall 2023  
**Platform:** [Google Earth Engine](https://earthengine.google.com/)  
**Region of Interest:** Lake Erie Watershed

---

## Overview

This project automates the detection and mapping of surface water extent within the Lake Erie watershed using Sentinel-1 Synthetic Aperture Radar (SAR) data. Leveraging Google Earth Engine’s cloud computing capabilities, it implements an adaptive thresholding technique to generate reliable binary water masks for any user-specified month and year within a season.

The approach enhances traditional methods by using buffered adaptive threshold sampling combined with edge detection to mitigate SAR speckle noise and artifacts, improving accuracy in surface water classification.

---

## Features

- Automated surface water classification from Sentinel-1 SAR data  
- User-defined season and year input for flexible temporal analysis  
- Adaptive thresholding using buffered edge sampling for improved accuracy  
- Noise reduction through edge length filtering and morphological operations  
- Customizable spatial resolution and polarization band selection  
- Region selection within the Great Lakes Basin (default: Lake Erie)

---

## How It Works

1. **Data Selection**: Sentinel-1 SAR imagery filtered by season and year, over the chosen basin area.  
2. **Edge Detection**: Canny edge detection identifies strong edges in SAR data representing water boundaries.  
3. **Adaptive Thresholding**: Samples pixel intensities within buffered edge zones to determine optimal water/non-water thresholds.  
4. **Binary Water Mask Creation**: Generates a refined binary map of surface water extent for the specified period.  
5. **Noise Filtering**: Removes small edge components and speckle artifacts by filtering edges below a minimum length.

---

## Usage

1. Open the script in the Google Earth Engine Code Editor.  
2. Adjust parameters at the top of the script:  
   - `userSeason`: `'Spring'`, `'Summer'`, `'Fall'`, or `'Winter'`  
   - `userYear`: e.g., `2023`  
   - `setResolution`: spatial resolution in meters (default 30)  
   - `band`: Sentinel-1 polarization band (`'VV'`, `'VH'`, `'HH'`, `'HV'`)  
   - Other edge detection and thresholding parameters as needed  
3. Run the script to generate surface water extent maps for the specified time window.  
4. Export results if desired.

---

## Parameters

| Parameter          | Description                                  | Default Value |
|--------------------|----------------------------------------------|---------------|
| `userSeason`       | Season to analyze (Spring, Summer, Fall, Winter) | `'Summer'`    |
| `userYear`         | Year of interest                             | `2023`        |
| `setResolution`    | Spatial resolution in meters (minimum 10)   | `30`          |
| `band`             | Sentinel-1 polarization band                 | `'VV'`        |
| `edgeLength`       | Minimum edge length for water edges (pixels)| `25`          |
| `edgeBuffer`       | Buffer size (meters) around edges for sampling| `60`          |
| `connectedPixels`  | Kernel size for morphological filtering      | `100`         |
| `cannyThreshold`   | Upper threshold for Canny edge detection     | `1`           |
| `cannySigma`       | Gaussian smoothing sigma for edge detection  | `0.5`         |
| `cannyLt`          | Lower threshold for Canny edge detection     | `0.05`        |

---

## Extending the Study Area

The study area defaults to Lake Erie basin but can be changed to other Great Lakes basins by modifying the filter line:

```javascript
var greatLakesBasin = ee.FeatureCollection('users/mlizak/GreatLakes_Basins');
var studyArea = greatLakesBasin.filter(ee.Filter.eq('merge', 'lk_erie')); // Options: lk_erie, lk_huron, lk_mich, lk_ont, lk_sup
```

Note:
The greatLakesBasin asset refers to a shapefile of the Great Lakes basins uploaded to GEE.
You can replace this asset and the studyArea filter with any other area of interest by importing your own vector data.

## Contributors:
- Maciej 'Mac' Lizak

---

