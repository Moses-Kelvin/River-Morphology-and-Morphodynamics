# River Morphodynamics Analysis of NIger River Using RivMAP

## Overview

This repository contains MATLAB scripts for analyzing multi-year river channel masks using the **RivMAP** (River Morphodynamics Analysis Package) toolbox. The workflow extracts centerlines, channel widths, sinuosity, curvature, migration rates, erosion/accretion patterns, and cutoffs from binary river masks derived from satellite imagery.

The analysis was originally developed for the **Kainji–Jebba reach of the River Niger, Nigeria** (2001–2025), but the code is adaptable to any large river system with sufficient time-series mask data.

---

## Author

**Kelvin Akhere Moses**  
Geomatics Graduate, University of Benin, Nigeria

---

## Requirements

### Software

- MATLAB R2020b or later
- Image Processing Toolbox
- Signal Processing Toolbox (for `findpeaks`)
- Mapping Toolbox _(optional, for coordinate transformations)_

### RivMAP Toolbox

Download the RivMAP toolbox from the official source and add it to your MATLAB path:

```matlab
addpath('path/to/RivMAP');
```

---

## Input Data

- Binary river masks as `.tif` files _(one per year)_
- Recommended format: single-band GeoTIFF
  - `1 = channel`
  - `0 = non-channel`

### Example Naming Convention

```text
Segment1_activeChannel_mask_2001.tif
```

---

## Helper Functions

- `pix2latlon` – converts pixel coordinates to latitude/longitude _(included)_
- `bwareaopen`, `imfill`, `imclose`, `imopen` – native MATLAB functions

---

## File Structure

```text
├── RivMAP/                                   # RivMAP toolbox (must be added to path)
├── Segment1_activeChannel_mask_*.tif         # Binary masks for Segment 1
├── Segment2_activeChannel_mask_*.tif         # Binary masks for Segment 2
├── pix2latlon.m                              # Pixel to lat/lon conversion
├── NigerRiver_Morphodynamics.m                  # Main analysis script
└── README.md                                 # This file
```

---

## How to Run

### 1. Clone this Repository

```bash
git clone https://github.com/Moses-Kelvin/River-Morphology-and-Morphodynamics.git
```

### 2. Add RivMAP to Your MATLAB Path

```matlab
addpath(genpath('path/to/RivMAP'));
```

### 3. Set Your Pixel Size

```matlab
pixel_size = 30;   % Landsat resolution (adjust for your data)
```

### 4. Run the Script

```matlab
RivMAP_Analysis_Script
```

---

## Outputs

Outputs are saved automatically as:

- CSV tables _(width, sinuosity, migration, network metrics)_
- GeoTIFF maps _(total migration, erosion, accretion)_
- Figures _(width profiles, migration rates, cutoff analysis)_
- `Max_Migration_Points.csv` with coordinates of maximum migration

---

## Key Features

| Feature               | Description                                                  |
| --------------------- | ------------------------------------------------------------ |
| Centerline extraction | Uses `centerline_from_mask` from RivMAP                      |
| Width estimation      | From binary mask and centerline                              |
| Sinuosity             | Channel length ÷ valley length                               |
| Curvature             | Smoothed centerline curvature                                |
| Braiding index        | Number of significant channel threads                        |
| Migration rates       | Centerline-based and mask-based _(erosion/accretion)_        |
| Cutoff detection      | Identifies meander cutoffs between consecutive years         |
| Spatial migration     | Per-segment migration rates along the meander belt           |
| GeoTIFF export        | Saves results with CRS _(EPSG:32631 – WGS84 / UTM Zone 31N)_ |
| Tables                | Exports summary statistics as CSV files                      |

---

## Output Files

| File                                     | Description                                                          |
| ---------------------------------------- | -------------------------------------------------------------------- |
| `Total_Migration_Area.tif`               | Binary map of all migration areas                                    |
| `Total_Erosion.tif`                      | Binary map of erosion                                                |
| `Total_Accretion.tif`                    | Binary map of accretion                                              |
| `Migration_Intensity_m.tif`              | Migration intensity _(meters)_                                       |
| `Table_4_3_Channel_Width_Statistics.csv` | Yearly width statistics                                              |
| `Table_4_4_Sinuosity.csv`                | Sinuosity and channel/valley lengths                                 |
| `Table_4_5_Network_Metrics.csv`          | Node density, branching intensity                                    |
| `Table_4_6_Migration.csv`                | Migration rates per interval                                         |
| `Morphodynamics_Summary.csv`             | Yearly summary _(width, sinuosity, curvature, wavelength, braiding)_ |
| `Max_Migration_Points.csv`               | Coordinates _(pixel and lat/lon)_ of maximum migration per interval  |
| `width_profile_*.png`                    | Width profile figures per year                                       |
| `cutoff_analysis.png`                    | Cutoff map and bar chart                                             |

---

## Example Visualization

The script automatically generates:

- Meander belt with segmented centerline nodes
- Spatial variation of migration rate along the river
- Cumulative erosion vs accretion balance over time
- Time series of:
  - sinuosity
  - width
  - curvature
  - braiding index
  - wavelength
- Reach-averaged migration rates:
  - centerline
  - erosion
  - accretion
- Cutoff maps colored by interval and cutoff-area bar charts

---

## Adapting to a Different River

To apply this workflow to another river system:

1. Prepare binary masks _(GeoTIFF, `1 = channel`, `0 = non-channel`)_
2. Update pixel size to match imagery resolution
3. Modify EPSG code in `geotiffwrite` for the local projected CRS
4. Adjust `Wn_initial` scaling factor
   - `0.4 × median width` works for most rivers
5. Set flow direction:
   - `'WE'` for east–west flow
   - `'NS'` for north–south flow

---

## Notes on RivMAP

RivMAP was developed by Schwenk et al. (2017) and Schwenk & Foufoula-Georgiou (2016). The toolbox is designed for analyzing river planform dynamics from satellite-derived masks.

This repository does **not** redistribute RivMAP; users must download it separately.

---

## Known Limitations

- Requires binary masks as input _(no direct image classification included)_
- Works best for single-thread or moderately braided rivers
- Cutoff detection may miss very small or gradual cutoffs

---

## Related Publications

Schwenk, J., Lanzoni, S., & Foufoula-Georgiou, E. (2015). _The life of a meander bend: Connecting shape and dynamics via analysis of a numerical model._ Journal of Geophysical Research: Earth Surface, 120(4), 690–710.

Schwenk, J., & Foufoula-Georgiou, E. (2016). _Meander cutoffs nonlocally accelerate upstream and downstream migration and channel widening._ Geophysical Research Letters, 43(24), 12,437–12,445.

Schwenk, J., Khandelwal, A., Fratkin, M., Kumar, V., & Foufoula-Georgiou, E. (2017). _High spatiotemporal resolution of river planform dynamics from Landsat: The RivMAP toolbox and results from the Ucayali River._ Earth and Space Science, 4(2), 46–75.

---

## License

This code is provided for academic and research purposes. Please cite the RivMAP authors if you use their toolbox.

---

## Contact

For questions or collaboration inquiries:

**Kelvin Akhere Moses**  
📧 moseskelvin683@gmail.com
