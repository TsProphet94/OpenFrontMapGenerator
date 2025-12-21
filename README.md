# Map Generator Tool

A complete workflow for generating stylised, game-ready terrain maps from real-world data using QGIS and Python.

## Features

*   **Topographic Data Collection**: Automatically downloads high-resolution elevation (DEM) data for any selected region.
*   **Stylised Rendering**: Generates game-ready PNG maps using a custom palette (`OpenFront_Palette.qml`) to match specific game aesthetics.
*   **Nation & Region Data**: Fetches administrative boundaries to identify nations and regions within the selected area.
*   **Coordinate Mapping**: Calculates precise coordinates for provinces and nations.
*   **Flag Assignment**: Automatically assigns appropriate flags to nations based on the downloaded data.

## Project Structure

```
.
├── OpenFrontMapGenerator.py    # QGIS Python script for data extraction
├── OpenFront_Palette.qml       # QGIS styling palette
├── scripts/
│   └── map_generator.py        # Python script for processing assets
├── StylisedMaps/               # Input folder for QGIS exports
├── generated/                  # Output folder for game assets
└── requirements.txt            # Python dependencies
```

## Prerequisites

1.  **QGIS** (v3.40 or compatible).
2.  **Python 3.x** installed.
3.  Install required Python packages:
    ```bash
    pip install -r requirements.txt
    ```

## Configuration

### OpenTopography API Key
The QGIS script requires an API key to download elevation data.
1.  Create a free account at [OpenTopography](https://portal.opentopography.org/myopentopo).
2.  Go to **myOpenTopo** > **Authorizations and API Key** and request a key.
3.  Open `OpenFrontMapGenerator.py` in a text editor.
4.  Find the line `API_KEY = "YOUR_OPENTOPO_API_KEY"` (around line 40).
5.  Replace `"YOUR_OPENTOPO_API_KEY"` with your actual key.

## Workflow

### Step 1: Export from QGIS

This step captures real-world elevation and province data.

1.  Download and open QGIS.
2.  Open the Python Console (**Plugins > Python Console**).
3.  Load the script `OpenFrontMapGenerator.py`.
    *   *Note: The script automatically detects the project folder. Ensure `StylisedMaps` exists in the same directory.*
4.  **Draw an Extent** on the canvas for the area you want to capture.
5.  Run the script.
    *   It downloads DEM data, applies `OpenFront_Palette.qml`, overlays rivers, and calculates province centers.
6.  **Output:** Files are saved to `StylisedMaps/<MapName>/`.

### Step 2: Generate Game Assets

This step processes the raw map into binary formats optimized for the game engine.

1.  Open a terminal in the project root.
2.  Run the generator script:

```bash
python scripts/map_generator.py "MapName"
```

*Replace `"MapName"` with the name of the folder inside `StylisedMaps` (e.g., "United Kingdom").*

**Optional Arguments:**
*   `--input`: Custom input directory (default: `StylisedMaps`).
*   `--output`: Custom output directory (default: `generated/maps`).
*   `--test`: Skip removal of small islands/lakes (faster for debugging).

### Step 3: Output

The tool generates **two versions** of the map in `generated/maps/`.
*Note: Output folder names are automatically converted to lowercase with no spaces (e.g., "United Kingdom" -> "unitedkingdom").*

1.  **`<mapname>`**: Full resolution (1:1 scale with QGIS export).
2.  **`<mapname>small`**: Half resolution (1:2 scale).

**Output Files:**
*   `map.bin`: Main binary terrain data.
*   `map4x.bin`: 1/2 scale binary data.
*   `map16x.bin`: 1/4 scale binary data.
*   `thumbnail.png`: Preview image.
*   `manifest.json`: Metadata with scaled coordinates for nations/provinces.

## Troubleshooting

*   **"ModuleNotFoundError"**: Run `pip install -r requirements.txt`.
*   **QGIS Path Errors**: Ensure you are running the script from a saved file, or manually update the `PROJECT_ROOT` variable in the Python script if running directly from the console buffer.

