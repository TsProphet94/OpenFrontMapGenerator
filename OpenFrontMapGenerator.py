# -*- coding: utf-8 -*-
"""
Map Exporter – DEM + Major Rivers/Lakes PNG & Province Centers (pixels) JSON  [QGIS 3.40 safe-load]

What it does (per run):
  • You draw an extent.
  • Downloads Copernicus GLO-90 (OpenTopography) DEM, auto-tiling if needed; mosaics tiles.
  • Applies your QML palette to the DEM and extracts its water/ocean colour.
  • Overlays major rivers + lakes (Natural Earth 10m), styled with your water colour.
  • Selects up to N major provinces/states (Natural Earth admin-1) intersecting the canvas,
    computes an interior label point for each, and converts those to pixel coords for the PNG.
  • Writes two files into: <ProjectRoot>/StylisedMaps/<MapName>/
      - <MapName>.png
      - <MapName>.json   (pixel coordinates for province centers)
"""

import os, json, zipfile, tempfile, math, datetime
from urllib.parse import urlencode

from qgis.PyQt import QtCore
from qgis.PyQt.QtCore import QSize, QByteArray
from qgis.PyQt.QtGui import QColor, QImage
from qgis.PyQt.QtNetwork import QNetworkRequest, QNetworkAccessManager

from qgis.core import (
    QgsProcessing, QgsProcessingAlgorithm,
    QgsProcessingParameterExtent, QgsProcessingParameterCrs,
    QgsProcessingParameterNumber, QgsProcessingParameterBoolean,
    QgsProcessingParameterFolderDestination, QgsProcessingParameterString,
    QgsProcessingParameterEnum,
    QgsProcessingException, QgsCoordinateReferenceSystem, QgsCoordinateTransform,
    QgsProject, QgsRasterLayer, QgsVectorLayer, QgsFeature, QgsGeometry,
    QgsMapSettings, QgsMapRendererParallelJob,
    QgsUnitTypes, QgsSimpleLineSymbolLayer, QgsLineSymbol, QgsFillSymbol,
    QgsColorRampShader, QgsRasterShader, QgsSingleBandPseudoColorRenderer,
    QgsRasterBandStats, QgsSingleSymbolRenderer, QgsFeatureRequest,
    QgsRasterDataProvider
)
from qgis import processing  # gdal:merge

# ======================= USER CONSTANTS ==========================
# Get a free API Key from: https://portal.opentopography.org/myopentopo
API_KEY = "dd4e0a24f02cbbadff55e7d1c5171672"
# Embedded OpenFront Palette (Discrete)
DEM_COLOR_RAMP = [
    (0.0, "#00006a", "Water ≤ 0 m"),
    (30.0, "#bedc8c", "Plains ≤ 30 m"),
    (60.0, "#beda8e", "Plains ≤ 60 m"),
    (90.0, "#bed890", "Plains ≤ 90 m"),
    (120.0, "#bed692", "Plains ≤ 120 m"),
    (150.0, "#bed494", "Plains ≤ 150 m"),
    (180.0, "#bed296", "Plains ≤ 180 m"),
    (210.0, "#bed098", "Plains ≤ 210 m"),
    (240.0, "#bece9a", "Plains ≤ 240 m"),
    (270.0, "#becc9c", "Plains ≤ 270 m"),
    (300.0, "#beca9e", "Plains ≤ 300 m"),
    (420.0, "#dccba0", "Highlands ≤ 420 m"),
    (540.0, "#decda2", "Highlands ≤ 540 m"),
    (660.0, "#e0cfa4", "Highlands ≤ 660 m"),
    (780.0, "#e2d1a6", "Highlands ≤ 780 m"),
    (900.0, "#e4d3a8", "Highlands ≤ 900 m"),
    (1020.0, "#e6d5aa", "Highlands ≤ 1020 m"),
    (1140.0, "#e8d7ac", "Highlands ≤ 1140 m"),
    (1260.0, "#ead9ae", "Highlands ≤ 1260 m"),
    (1380.0, "#ecdbb0", "Highlands ≤ 1380 m"),
    (1500.0, "#eeddb2", "Highlands ≤ 1500 m"),
    (1818.0, "#f0f0b4", "Mountains ≤ 1818 m"),
    (2136.0, "#f0f0b6", "Mountains ≤ 2136 m"),
    (2455.0, "#f1f1b8", "Mountains ≤ 2455 m"),
    (2773.0, "#f2f2ba", "Mountains ≤ 2773 m"),
    (3091.0, "#f2f2bc", "Mountains ≤ 3091 m"),
    (3409.0, "#f2f2be", "Mountains ≤ 3409 m"),
    (3727.0, "#f3f3c0", "Mountains ≤ 3727 m"),
    (4045.0, "#f4f4c2", "Mountains ≤ 4045 m"),
    (4364.0, "#f4f4c4", "Mountains ≤ 4364 m"),
    (4682.0, "#f4f4c6", "Mountains ≤ 4682 m"),
    (5000.0, "#f5f5c8", "Mountains ≤ 5000 m"),
]

# Determine project root relative to this script
try:
    PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # Fallback if __file__ is not defined (e.g. running directly in console without saving)
    # PLEASE UPDATE THIS PATH TO YOUR PROJECT FOLDER IF RUNNING FROM CONSOLE
    PROJECT_ROOT = os.path.join(os.path.expanduser("~"), "Desktop", "OpenFrontMapGenerator")

EXPORT_DIR   = os.path.join(PROJECT_ROOT, "StylisedMaps")
FLAGS_DIR    = os.path.join(PROJECT_ROOT, "Examples", "flags")

# Water detail knobs:
MAJOR_WATER_SCALERANK_MAX = 3      # allow up to 3 so features appear in more places
MIN_LAKE_AREA_KM2          = 200.0 # skip small lakes
# Province selection knobs:
MIN_PROVINCE_AREA_KM2      = 1500.0 # skip tiny slivers/islands
# ================================================================


# ------------------------------- helpers --------------------------------------

def _download_bytes(url, headers=None, post_data=None, timeout_ms=120000):
    """Qt-friendly HTTP GET/POST; returns bytes or raises QgsProcessingException."""
    nam = QNetworkAccessManager()
    req = QNetworkRequest(QtCore.QUrl(url))
    if headers:
        for k, v in headers.items():
            req.setRawHeader(k.encode(), v.encode())
    if post_data is None:
        reply = nam.get(req)
    else:
        if isinstance(post_data, dict):
            post_data = urlencode(post_data).encode("utf-8")
        elif isinstance(post_data, str):
            post_data = post_data.encode("utf-8")
        reply = nam.post(req, QByteArray(post_data))
    loop = QtCore.QEventLoop()
    reply.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(timeout_ms, loop.quit)
    loop.exec_()
    if not reply.isFinished():
        reply.abort()
        raise QgsProcessingException(f"Network timeout for {url}")
    if reply.error():
        raise QgsProcessingException(f"Network error for {url}: {reply.errorString()}")
    return bytes(reply.readAll())

def _ensure_folder(path):
    os.makedirs(path, exist_ok=True)
    return path

# Global list to track temp files for cleanup
_TEMP_FILES = []

def _save_temp(b, suffix):
    fd, p = tempfile.mkstemp(suffix=suffix); os.close(fd)
    with open(p, "wb") as f: f.write(b)
    _TEMP_FILES.append(p)
    return p

def _cleanup_temp_files(feedback=None):
    """Deletes all tracked temporary files."""
    count = 0
    for p in _TEMP_FILES:
        if os.path.exists(p):
            try:
                os.remove(p)
                count += 1
            except Exception as e:
                if feedback: feedback.pushWarning(f"Failed to delete temp file {p}: {e}")
    if feedback and count > 0:
        feedback.pushInfo(f"Cleaned up {count} temporary files.")
    _TEMP_FILES.clear()

def _apply_dem_style(layer, ramp_list):
    """Applies the discrete color ramp to the DEM layer."""
    if not layer.isValid(): return
    
    ramp_items = []
    for val, color_hex, label in ramp_list:
        ramp_items.append(QgsColorRampShader.ColorRampItem(val, QColor(color_hex), label))
    
    fcn = QgsColorRampShader()
    fcn.setColorRampType(QgsColorRampShader.Discrete)
    fcn.setColorRampItemList(ramp_items)
    
    shader = QgsRasterShader()
    shader.setRasterShaderFunction(fcn)
    
    renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, shader)
    layer.setRenderer(renderer)
    layer.triggerRepaint()

def _deg_size_km(lat_deg):
    """~111 km per degree latitude; longitude shrinks by cos(lat)."""
    return 111.0, 111.0 * max(0.0, min(1.0, abs(math.cos(math.radians(lat_deg)))))

def _bbox_area_km2(south, west, north, east):
    lat_mid = 0.5 * (south + north)
    km_lat, km_lon = _deg_size_km(lat_mid)
    return max(0.0, north - south) * km_lat * max(0.0, east - west) * km_lon

def _download_ne_zip(ne_url, out_dir):
    """Download a Natural Earth zip and extract; return the first .shp path."""
    zbytes = _download_bytes(ne_url, timeout_ms=240000)
    zpath = _save_temp(zbytes, ".zip")
    with zipfile.ZipFile(zpath) as zf:
        zf.extractall(out_dir)
    for fname in os.listdir(out_dir):
        if fname.lower().endswith(".shp"):
            return os.path.join(out_dir, fname)
    return None

def _parse_int(val):
    if val is None: return None
    try:
        return int(val)
    except Exception:
        try:
            return int(float(val))
        except Exception:
            return None

def _best_attr(f, names):
    """First present stringy attribute from a list of field names."""
    for n in names:
        if n in f.fields().names():
            v = f[n]
            if v not in (None, ""):
                return str(v)
    return ""

def _geom_area_km2(geom, crs_authid):
    """Approx area by projecting to EPSG:3857 and measuring m² (→ km²)."""
    from qgis.core import QgsCoordinateTransform, QgsCoordinateReferenceSystem
    if not geom or geom.isEmpty(): return 0.0
    try:
        src = QgsCoordinateReferenceSystem(crs_authid)
        tgt = QgsCoordinateReferenceSystem("EPSG:3857")
        g = QgsGeometry(geom)
        g.transform(QgsCoordinateTransform(src, tgt, QgsProject.instance().transformContext()))
        return abs(g.area()) / 1_000_000.0
    except Exception:
        return 0.0


# -------------------------- Processing Algorithm ------------------------------

class MapExporterDemRiversProvPixels(QgsProcessingAlgorithm):
    P_EXTENT="EXTENT"; P_TARGET_CRS="TARGET_CRS"
    P_WIDTH_PX="WIDTH_PX"; P_DPI="DPI"; P_TRANSPARENT="TRANSPARENT_BG"
    P_RIVER_WIDTH="RIVER_WIDTH"; P_CACHE="CACHE_FOLDER"
    P_MAX_PROVINCES="MAX_PROVINCES"
    P_AUTO_RIVER_SCALE = "AUTO_RIVER_SCALE"
    P_MIN_LAKE_AREA_PX2 = "MIN_LAKE_AREA_PX2"
    P_MAP_NAME = "MAP_NAME"
    P_DYNAMIC_SCALE = "DYNAMIC_SCALE"
    P_DEM_SOURCE = "DEM_SOURCE"

    def name(self): return "map_exporter_dem_majorwaters_provincepixels_safe"
    def displayName(self): return "Map Exporter – DEM + Major Rivers/Lakes & Province Centers (pixels)"
    def group(self): return "Scripts"
    def groupId(self): return "scripts"

    def shortHelpString(self):
        return ("DEM (COP90) styled with your QML, overlay major rivers/lakes, "
                "export PNG + JSON (province centers as pixel coords). "
                "Outputs go to <ProjectRoot>/StylisedMaps/")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterExtent(self.P_EXTENT, "Draw extent (map CRS)"))
        self.addParameter(QgsProcessingParameterString(self.P_MAP_NAME, "Map Name (optional)", defaultValue="", optional=True))
        self.addParameter(QgsProcessingParameterCrs(self.P_TARGET_CRS, "Output CRS", defaultValue="EPSG:3857"))
        self.addParameter(QgsProcessingParameterNumber(
            self.P_WIDTH_PX, "PNG width (px)", QgsProcessingParameterNumber.Integer,
            defaultValue=4096, minValue=256, maxValue=32768
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.P_DPI, "PNG DPI", QgsProcessingParameterNumber.Integer,
            defaultValue=300, minValue=72, maxValue=1200
        ))
        self.addParameter(QgsProcessingParameterEnum(
            self.P_DEM_SOURCE, "DEM Source",
            options=["Auto (Adaptive)", "COP30 (High Res)", "COP90 (Medium Res)", "SRTM15+ (Low Res)"],
            defaultValue=0
        ))
        self.addParameter(QgsProcessingParameterBoolean(self.P_DYNAMIC_SCALE, "Dynamic elevation scaling (fit palette to local max height)", defaultValue=True))
        self.addParameter(QgsProcessingParameterBoolean(self.P_TRANSPARENT, "Transparent background", defaultValue=True))
        self.addParameter(QgsProcessingParameterNumber(
            self.P_RIVER_WIDTH, "River stroke width (px)", QgsProcessingParameterNumber.Double,
            defaultValue=2.0, minValue=0.1, maxValue=20.0
        ))
        self.addParameter(QgsProcessingParameterBoolean(self.P_AUTO_RIVER_SCALE, "Auto-scale river width by image width", defaultValue=True))
        self.addParameter(QgsProcessingParameterFolderDestination(self.P_CACHE, "Cache folder"))
        self.addParameter(QgsProcessingParameterNumber(
            self.P_MAX_PROVINCES, "Max provinces to output", QgsProcessingParameterNumber.Integer,
            defaultValue=20, minValue=1, maxValue=100
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.P_MIN_LAKE_AREA_PX2, "Min lake area (px²) for inclusion", QgsProcessingParameterNumber.Integer,
            defaultValue=16, minValue=0, maxValue=100000
        ))

    def processAlgorithm(self, params, context, feedback):
        # Build output filenames
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        map_name = self.parameterAsString(params, self.P_MAP_NAME, context).strip()
        
        if map_name:
            # Sanitize filename
            safe_name = "".join(c for c in map_name if c.isalnum() or c in (' ', '.', '_', '-')).strip()
            base_name = safe_name if safe_name else f"map_{ts}"
        else:
            base_name = f"map_{ts}"

        # Create subfolder for the map
        map_folder = os.path.join(EXPORT_DIR, base_name)
        _ensure_folder(map_folder)
        
        # Inputs
        extent = self.parameterAsExtent(params, self.P_EXTENT, context)
        out_crs = self.parameterAsCrs(params, self.P_TARGET_CRS, context)
        width_px = int(self.parameterAsInt(params, self.P_WIDTH_PX, context))
        dpi = int(self.parameterAsInt(params, self.P_DPI, context))
        dynamic_scale = self.parameterAsBool(params, self.P_DYNAMIC_SCALE, context)
        transparent = self.parameterAsBool(params, self.P_TRANSPARENT, context)
        river_width = float(self.parameterAsDouble(params, self.P_RIVER_WIDTH, context))
        auto_river_scale = bool(self.parameterAsBool(params, self.P_AUTO_RIVER_SCALE, context))
        cache_dir = self.parameterAsString(params, self.P_CACHE, context)
        max_provinces = int(self.parameterAsInt(params, self.P_MAX_PROVINCES, context))
        min_lake_area_px2 = int(self.parameterAsInt(params, self.P_MIN_LAKE_AREA_PX2, context))

        _ensure_folder(cache_dir)

        proj = QgsProject.instance()
        proj_crs = proj.crs() if proj.crs().isValid() else QgsCoordinateReferenceSystem("EPSG:3857")
        to4326 = QgsCoordinateTransform(proj_crs, QgsCoordinateReferenceSystem("EPSG:4326"), context.transformContext())
        to_out = QgsCoordinateTransform(proj_crs, out_crs, context.transformContext())

        # Extents & sizes
        extent_out = to_out.transformBoundingBox(extent)
        extent_ll = to4326.transformBoundingBox(extent)
        south, west, north, east = extent_ll.yMinimum(), extent_ll.xMinimum(), extent_ll.yMaximum(), extent_ll.xMaximum()
        extent_out_geom = QgsGeometry.fromRect(extent_out)

        aspect = extent_out.width() / max(1e-9, extent_out.height())
        height_px = max(1, int(width_px / max(1e-9, aspect)))
        xmin = extent_out.xMinimum(); xmax = extent_out.xMaximum()
        ymin = extent_out.yMinimum(); ymax = extent_out.yMaximum()
        width_units = xmax - xmin; height_units = ymax - ymin

        # Check resolution and warn if low for large areas
        if out_crs.mapUnits() == QgsUnitTypes.DistanceMeters:
             res_m = width_units / width_px
             if res_m > 1000:
                 feedback.pushWarning(f"Output resolution is low ({res_m:.0f} m/px). Consider increasing Width PX for large areas.")

        def world_to_pixel(x, y):
            px = (x - xmin) / width_units * width_px
            py = (ymax - y) / height_units * height_px  # origin = top-left
            return int(round(px)), int(round(py))

        # Compute effective river width (in pixels) with optional scaling for large images
        if auto_river_scale:
            # Scale sublinearly with image width; cap to reasonable range
            scale = max(1.0, (width_px / 2048.0) ** 0.6)
            river_width_eff = max(1.0, min(20.0, river_width * scale))
        else:
            river_width_eff = river_width
        feedback.pushInfo(f"River stroke width (effective): {river_width_eff:.2f} px")

        # ---------------- 1) DEM (OpenTopography) with auto-tiling ----------------
        dem_source_idx = self.parameterAsEnum(params, self.P_DEM_SOURCE, context)
        area_km2 = _bbox_area_km2(south, west, north, east)
        
        # Determine DEM type and tile size limit
        dem_type = "COP90"
        tile_max_km2 = 500_000.0
        
        if dem_source_idx == 1: # COP30
            dem_type = "COP30"
            tile_max_km2 = 50_000.0
        elif dem_source_idx == 2: # COP90
            dem_type = "COP90"
            tile_max_km2 = 500_000.0
        elif dem_source_idx == 3: # SRTM15+
            dem_type = "SRTM15Plus"
            tile_max_km2 = 10_000_000.0
        else: # Auto
            if area_km2 < 25000:
                dem_type = "COP30"
                tile_max_km2 = 50_000.0
            elif area_km2 < 2000000:
                dem_type = "COP90"
                tile_max_km2 = 500_000.0
            else:
                dem_type = "SRTM15Plus"
                tile_max_km2 = 10_000_000.0
        
        feedback.pushInfo(f"DEM Source: {dem_type} (Area: {area_km2:,.0f} km²)")

        # Calculate target resolution in degrees for downsampling optimization
        # COP90 is ~0.000833 deg/px. If we need much less, we downsample immediately.
        req_res_deg = min((east - west) / width_px, (north - south) / height_px)
        target_res_deg = req_res_deg / 2.0  # Nyquist safety factor
        should_downsample = target_res_deg > 0.002  # Only if > ~2x coarser than native (approx)

        def _download_dem_tile(s, w, n, e):
            qs = dict(demtype=dem_type, south=s, north=n, west=w, east=e, outputFormat="GTiff")
            if API_KEY: qs["API_Key"] = API_KEY
            url = "https://portal.opentopography.org/API/globaldem?" + urlencode(qs)
            headers = {"User-Agent": "QGIS Map Exporter Script/1.0"}
            return _save_temp(_download_bytes(url, headers=headers, timeout_ms=240000), ".tif")

        if area_km2 <= tile_max_km2:
            feedback.pushInfo("Downloading single DEM tile…")
            dem_path = _download_dem_tile(south, west, north, east)
            # Optional: Downsample single tile if huge? (Usually not needed if area is small)
        else:
            feedback.pushInfo("Area too large → tiling and mosaicking…")
            lat_span = north - south
            lon_span = east - west
            km_lat, km_lon = _deg_size_km(0.5*(south+north))
            km2_per_deg2 = km_lat * km_lon
            total_deg2 = max(1e-9, lat_span * lon_span)
            tiles_needed = max(1, math.ceil((total_deg2 * km2_per_deg2) / tile_max_km2))
            rows = cols = max(1, math.ceil(math.sqrt(tiles_needed)))
            feedback.pushInfo(f"Tiling into {rows}×{cols} = {rows*cols} tiles")
            if should_downsample:
                feedback.pushInfo(f"Optimization: Downsampling tiles to {target_res_deg:.5f} deg/px before merge.")

            dem_tiles = []
            for r in range(rows):
                s = south + (lat_span * r / rows)
                n = south + (lat_span * (r + 1) / rows)
                for c in range(cols):
                    w = west + (lon_span * c / cols)
                    e = west + (lon_span * (c + 1) / cols)
                    feedback.pushInfo(f"  → tile r{r+1}/c{c+1}: S{s:.4f}, W{w:.4f}, N{n:.4f}, E{e:.4f}")
                    try:
                        raw_tile = _download_dem_tile(s, w, n, e)
                        if should_downsample:
                            # Downsample immediately to save disk/memory for merge
                            ds_tile = _save_temp(b"", ".tif"); os.remove(ds_tile)
                            processing.run("gdal:warpreproject", {
                                'INPUT': raw_tile,
                                'SOURCE_CRS': QgsCoordinateReferenceSystem("EPSG:4326"),
                                'TARGET_CRS': QgsCoordinateReferenceSystem("EPSG:4326"),
                                'RESAMPLING': 1, # Bilinear
                                'NODATA': None,
                                'TARGET_RESOLUTION': target_res_deg,
                                'OPTIONS': "",
                                'DATA_TYPE': 0,
                                'EXTRA': "",
                                'OUTPUT': ds_tile
                            }, context=context, feedback=feedback)
                            
                            # Try to delete huge raw tile; if locked, skip (cleanup will catch it later)
                            try:
                                os.remove(raw_tile)
                            except OSError:
                                feedback.pushWarning(f"Could not delete intermediate tile {raw_tile} (locked).")
                                
                            dem_tiles.append(ds_tile)
                        else:
                            dem_tiles.append(raw_tile)
                    except Exception as ex:
                        raise QgsProcessingException(f"DEM tile r{r+1} c{c+1} failed: {ex}")
            feedback.pushInfo("Mosaicking DEM tiles…")
            merged = _save_temp(b"", ".tif"); os.remove(merged)
            processing.run(
                "gdal:merge",
                {"INPUT": dem_tiles, "PCT": False, "SEPARATE": False,
                 "NODATA_INPUT": None, "NODATA_OUTPUT": None, "OPTIONS": "",
                 "DATA_TYPE": 5, "EXTRA": "", "OUTPUT": merged},
                context=context, feedback=feedback
            )
            dem_path = merged

        # ---------------- 1b) Warp & Smooth DEM to Output Resolution ----------------
        # We warp the raw DEM to the target CRS and exact output pixel dimensions using Bilinear resampling.
        # This effectively smooths the coastlines and terrain to match the output map, preventing jagged edges.
        feedback.pushInfo("Warping and smoothing DEM to output resolution...")
        dem_warped = _save_temp(b"", ".tif"); os.remove(dem_warped)
        processing.run("gdal:warpreproject", {
            'INPUT': dem_path,
            'SOURCE_CRS': QgsCoordinateReferenceSystem("EPSG:4326"), # COP90 is usually 4326
            'TARGET_CRS': out_crs,
            'RESAMPLING': 1, # Bilinear for smoothing
            'NODATA': None,
            'TARGET_RESOLUTION': None,
            'OPTIONS': "",
            'DATA_TYPE': 0, # Use input type (Float)
            'EXTRA': f'-ts {width_px} {height_px}', # Force exact output size
            'OUTPUT': dem_warped
        }, context=context, feedback=feedback)
        
        dem = QgsRasterLayer(dem_warped, "DEM")
        if not dem.isValid():
            raise QgsProcessingException("Failed to load Warped DEM.")
        
        # Calculate dynamic ramp if requested
        current_ramp = DEM_COLOR_RAMP
        if dynamic_scale:
            feedback.pushInfo("Calculating DEM statistics for dynamic scaling...")
            # Calculate stats on the actual extent
            provider = dem.dataProvider()
            stats = provider.bandStatistics(1, QgsRasterBandStats.Max, extent_out, 0)
            max_z = stats.maximumValue
            feedback.pushInfo(f"Local Max Elevation: {max_z:.2f} m")
            
            if max_z > 0:
                # Scale the ramp so that the highest color (5000m) maps to max_z
                # But keep Water (<=0) fixed.
                scale_factor = max_z / 5000.0
                # Avoid extreme scaling for very flat areas (e.g. < 10m) to prevent noise amplification?
                # User asked for "nice palette" even for flat areas, so we scale fully.
                
                new_ramp = []
                for val, color, label in DEM_COLOR_RAMP:
                    if val <= 0:
                        new_ramp.append((val, color, label))
                    else:
                        new_val = val * scale_factor
                        new_ramp.append((new_val, color, f"{label} (scaled)"))
                current_ramp = new_ramp

        # Apply embedded style
        _apply_dem_style(dem, current_ramp)
        
        # Extract ocean color from the first item in our ramp (Water <= 0)
        ocean_hex = current_ramp[0][1] if current_ramp else "#3a78c2"
        ocean_qcolor = QColor(ocean_hex)
        
        feedback.pushInfo(f"Ocean color used: {ocean_hex}")
        QgsProject.instance().addMapLayer(dem)

        # ---------------- 2) Major rivers & lakes (Natural Earth) ----------------
        ne_dir = _ensure_folder(os.path.join(cache_dir, "natural_earth"))

        # Rivers: centerlines (as MultiLineString)
        rivers_mem = None
        ne_riv = os.path.join(ne_dir, "ne_10m_rivers_lake_centerlines.shp")
        if not os.path.exists(ne_riv):
            # Try multiple mirrors for Natural Earth
            ne_riv = _download_ne_zip("https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_rivers_lake_centerlines.zip", ne_dir)
            if not ne_riv or not os.path.exists(ne_riv):
                ne_riv = _download_ne_zip("https://naciscdn.org/naturalearth/10m/physical/ne_10m_rivers_lake_centerlines.zip", ne_dir)
        if not ne_riv or not os.path.exists(ne_riv):
            ne_riv = None
            feedback.pushWarning("Rivers shapefile not found or downloaded.")
        else:
            feedback.pushInfo(f"Rivers shapefile: {ne_riv}")
        if ne_riv:
            rivers_src = QgsVectorLayer(ne_riv, "NE Rivers", "ogr")
            feedback.pushInfo(f"Rivers layer valid: {rivers_src.isValid()}")
            if rivers_src.isValid():
                rivers_mem = QgsVectorLayer("MultiLineString?crs=" + out_crs.authid(), "Rivers (major)", "memory")
                dp = rivers_mem.dataProvider()
                dp.addAttributes(rivers_src.fields()); rivers_mem.updateFields()
                to_out_r = QgsCoordinateTransform(rivers_src.crs(), out_crs, context.transformContext())
                
                # Create filter rect in source CRS to optimize large datasets
                tr_r = QgsCoordinateTransform(out_crs, rivers_src.crs(), context.transformContext())
                filter_rect_r = tr_r.transformBoundingBox(extent_out)

                feats = []
                def _collect_river_feats(relax=False):
                    local = []
                    req = QgsFeatureRequest().setFilterRect(filter_rect_r)
                    for f in rivers_src.getFeatures(req):
                        sr = _parse_int(f["scalerank"] if "scalerank" in rivers_src.fields().names() else None)
                        # Include features when scalerank missing; if relax=True, ignore scalerank entirely
                        if not relax and (sr is not None and sr > MAJOR_WATER_SCALERANK_MAX):
                            continue
                        g = f.geometry()
                        if not g: continue
                        g = QgsGeometry(g); g.transform(to_out_r)
                        if not g.intersects(extent_out_geom): continue
                        g_clip = g.intersection(extent_out_geom)
                        if not g_clip or g_clip.isEmpty(): continue
                        nf = QgsFeature(rivers_mem.fields()); nf.setAttributes(f.attributes()); nf.setGeometry(g_clip); local.append(nf)
                    return local

                feats = _collect_river_feats(relax=False)
                feedback.pushInfo(f"Rivers features collected (initial): {len(feats)}")
                if not feats:
                    feedback.pushInfo("No rivers after initial filter → relaxing scalerank filter…")
                    feats = _collect_river_feats(relax=True)
                    feedback.pushInfo(f"Rivers features collected (relaxed): {len(feats)}")
                if feats:
                    dp.addFeatures(feats); rivers_mem.updateExtents()
                    line = QgsSimpleLineSymbolLayer(); line.setColor(ocean_qcolor); line.setWidth(river_width_eff); line.setWidthUnit(QgsUnitTypes.RenderPixels)
                    sym = QgsLineSymbol(); sym.changeSymbolLayer(0, line)
                    rivers_mem.setRenderer(QgsSingleSymbolRenderer(sym))
                    rivers_mem.triggerRepaint()
                    feedback.pushInfo(f"Rivers added: {len(feats)} features, color: {ocean_qcolor.name()}")
                    QgsProject.instance().addMapLayer(rivers_mem)
                else:
                    feedback.pushWarning("No river features to add.")

        # Lakes: polygons (as MultiPolygon)
        lakes_mem = None
        ne_lakes = os.path.join(ne_dir, "ne_10m_lakes.shp")
        if not os.path.exists(ne_lakes):
            ne_lakes = _download_ne_zip("https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_lakes.zip", ne_dir)
            if not ne_lakes or not os.path.exists(ne_lakes):
                ne_lakes = _download_ne_zip("https://naciscdn.org/naturalearth/10m/physical/ne_10m_lakes.zip", ne_dir)
        if ne_lakes and os.path.exists(ne_lakes):
            feedback.pushInfo(f"Lakes shapefile: {ne_lakes}")
            lakes_src = QgsVectorLayer(ne_lakes, "NE Lakes", "ogr")
            feedback.pushInfo(f"Lakes layer valid: {lakes_src.isValid()}")
            if lakes_src.isValid():
                lakes_mem = QgsVectorLayer("MultiPolygon?crs=" + out_crs.authid(), "Lakes (major)", "memory")
                dp = lakes_mem.dataProvider()
                dp.addAttributes(lakes_src.fields()); lakes_mem.updateFields()
                to_out_l = QgsCoordinateTransform(lakes_src.crs(), out_crs, context.transformContext())
                
                # Create filter rect in source CRS
                tr_l = QgsCoordinateTransform(out_crs, lakes_src.crs(), context.transformContext())
                filter_rect_l = tr_l.transformBoundingBox(extent_out)

                feats = []
                def _collect_lake_feats(relax=False):
                    local = []
                    req = QgsFeatureRequest().setFilterRect(filter_rect_l)
                    for f in lakes_src.getFeatures(req):
                        sr = _parse_int(f["scalerank"] if "scalerank" in lakes_src.fields().names() else None)
                        # Include when scalerank missing; if relax=True, ignore scalerank and area threshold
                        if not relax and (sr is not None and sr > MAJOR_WATER_SCALERANK_MAX):
                            continue
                        g = f.geometry()
                        if not g: continue
                        g = QgsGeometry(g); g.transform(to_out_l)
                        if not g.intersects(extent_out_geom): continue
                        g_clip = g.intersection(extent_out_geom)
                        if not g_clip or g_clip.isEmpty(): continue
                        if not relax:
                            # Dual threshold: km² or pixel area
                            area_km2 = _geom_area_km2(g_clip, out_crs.authid())
                            # Approximate pixel area using map extent scaling
                            try:
                                area_units = g_clip.area()
                                px_per_unit_x = width_px / max(1e-12, width_units)
                                px_per_unit_y = height_px / max(1e-12, height_units)
                                area_px2 = area_units * px_per_unit_x * px_per_unit_y
                            except Exception:
                                area_px2 = 0.0
                            if area_km2 < MIN_LAKE_AREA_KM2 and area_px2 < float(min_lake_area_px2):
                                continue
                        nf = QgsFeature(lakes_mem.fields()); nf.setAttributes(f.attributes()); nf.setGeometry(g_clip); local.append(nf)
                    return local

                feats = _collect_lake_feats(relax=False)
                feedback.pushInfo(f"Lakes features collected (initial): {len(feats)}")
                if not feats:
                    feedback.pushInfo("No lakes after initial filter → relaxing scalerank/area filters…")
                    feats = _collect_lake_feats(relax=True)
                    feedback.pushInfo(f"Lakes features collected (relaxed): {len(feats)}")
                if feats:
                    dp.addFeatures(feats); lakes_mem.updateExtents()
                    fill = QgsFillSymbol.createSimple({"color": ocean_qcolor.name(), "outline_style": "no"})
                    lakes_mem.setRenderer(QgsSingleSymbolRenderer(fill))
                    lakes_mem.triggerRepaint()
                    feedback.pushInfo(f"Lakes added: {len(feats)} features, color: {ocean_qcolor.name()}")
                    QgsProject.instance().addMapLayer(lakes_mem)
                else:
                    feedback.pushWarning("No lake features to add.")
        else:
            feedback.pushWarning("Lakes shapefile not found or downloaded.")

        # ---------------- 3) Province centers (select top N) ---------------------
        feedback.pushInfo("Selecting major provinces/states and computing pixel centers…")
        ne_admin1 = os.path.join(ne_dir, "ne_10m_admin_1_states_provinces.shp")
        if not os.path.exists(ne_admin1):
            ne_admin1 = _download_ne_zip("https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_1_states_provinces.zip", ne_dir)
            if not ne_admin1 or not os.path.exists(ne_admin1):
                ne_admin1 = _download_ne_zip("https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_1_states_provinces.zip", ne_dir)
        if not ne_admin1 or not os.path.exists(ne_admin1):
            raise QgsProcessingException("Failed to download Natural Earth admin-1 data.")
        admin1 = QgsVectorLayer(ne_admin1, "Admin-1", "ogr")
        if not admin1.isValid():
            raise QgsProcessingException("Failed to load Natural Earth admin-1 (states/provinces).")

        to_out_a1 = QgsCoordinateTransform(admin1.crs(), out_crs, context.transformContext())

        candidates = []
        for f in admin1.getFeatures():
            g = f.geometry()
            if not g: continue
            g = QgsGeometry(g); g.transform(to_out_a1)
            if not g.intersects(extent_out_geom): continue
            g_clip = g.intersection(extent_out_geom)
            if not g_clip or g_clip.isEmpty(): continue

            labelrank = _parse_int(f["labelrank"] if "labelrank" in admin1.fields().names() else None)
            if labelrank is None: labelrank = 9  # default less important
            area_km2 = _geom_area_km2(g_clip, out_crs.authid())
            
            # Relax area threshold for very important provinces (likely capitals)
            # labelrank <= 3 usually indicates high importance
            if area_km2 < MIN_PROVINCE_AREA_KM2 and labelrank > 3:
                continue

            center = g_clip.pointOnSurface()
            if center.isEmpty(): continue
            pt = center.asPoint()
            px, py = world_to_pixel(pt.x(), pt.y())

            candidates.append({
                "rank": (labelrank, -area_km2),
                "admin": _best_attr(f, ["admin", "ADM0_NAME", "ADM0NAME", "admin_name"]),
                "province": _best_attr(f, ["name", "name_en", "name_1", "name_local"]),
                "code": _best_attr(f, ["iso_3166_2", "adm1_code", "code_hasc"]),
                "pixel_x": px, "pixel_y": py
            })

        # Scan for available flags
        valid_flags = set()
        if os.path.exists(FLAGS_DIR):
            for fn in os.listdir(FLAGS_DIR):
                if fn.endswith(".svg"):
                    valid_flags.add(fn.rsplit('.', 1)[0].lower())
        else:
            feedback.pushWarning(f"Flags directory not found: {FLAGS_DIR}")

        # Sort, dedupe (admin+province), and select with country balancing
        # 1. Deduplicate
        candidates.sort(key=lambda d: d["rank"])
        unique_candidates = []
        seen = set()
        for c in candidates:
            key = (c["admin"], c["province"])
            if key in seen:
                continue
            seen.add(key)
            unique_candidates.append(c)

        # 2. Group by Country
        by_country = {}
        for c in unique_candidates:
            country = c["admin"]
            if country not in by_country:
                by_country[country] = []
            by_country[country].append(c)

        # 3. Selection Strategy:
        #    - Ensure at least 1 top province from each country (if possible)
        #    - Fill remaining slots from the global pool of remaining provinces
        
        final_selection = []
        remaining_pool = []
        
        # Sort countries by the rank of their best province to prioritize "important" countries first
        # if we have more countries than max_provinces (rare but possible)
        sorted_countries = sorted(by_country.keys(), key=lambda k: by_country[k][0]["rank"])
        
        # Pass 1: Take the best province from each country
        for country in sorted_countries:
            provinces = by_country[country]
            # Best one (already sorted by rank in candidates)
            final_selection.append(provinces[0])
            # Add rest to pool
            for p in provinces[1:]:
                remaining_pool.append(p)
        
        # Pass 2: Fill remaining slots from the pool
        remaining_pool.sort(key=lambda d: d["rank"])
        
        needed = max_provinces - len(final_selection)
        if needed > 0:
            final_selection.extend(remaining_pool[:needed])
        else:
            # If we have more countries than slots, we trim the selection.
            # Since we sorted countries by importance, we keep the most important ones.
            final_selection = final_selection[:max_provinces]
            
        # Sort final output by rank (importance)
        final_selection.sort(key=lambda d: d["rank"])

        points = []
        for c in final_selection:
            iso = c["code"]
            raw_flag = iso.split('-')[0].lower() if iso else "xx"
            flag_code = raw_flag if raw_flag in valid_flags else "xx"
            
            points.append({
                "name": c["province"],
                "flag": flag_code,
                "pixel_x": c["pixel_x"],
                "pixel_y": c["pixel_y"]
            })

        # ---------------- 4) Render PNG (DEM + water overlays) -------------------
        
        # --- Compose and render map (DEM at bottom, lakes, then rivers) ---
        ms = QgsMapSettings()
        # QgsMapSettings.setLayers() expects the list in Top-to-Bottom order (index 0 is top).
        layer_order = []
        if rivers_mem and rivers_mem.isValid():
            layer_order.append(rivers_mem)
        if lakes_mem and lakes_mem.isValid():
            layer_order.append(lakes_mem)
        layer_order.append(dem)
        
        # Ensure we actually include vectors in output
        if len(layer_order) == 1:
            feedback.pushWarning("No valid rivers/lakes layers to render; exporting DEM only.")
        else:
            feedback.pushInfo(f"Rendering layers (Top->Bottom): {[l.name() for l in layer_order]}")
        ms.setLayers(layer_order)
        ms.setDestinationCrs(out_crs)
        # Transparent background if requested; otherwise solid opaque backdrop
        ms.setBackgroundColor(QColor(0, 0, 0, 0 if transparent else 255))
        ms.setExtent(extent_out)
        ms.setOutputSize(QSize(width_px, height_px))
        ms.setOutputDpi(dpi)
        
        # Disable Antialiasing to ensure strict palette colors (no blending/blurring)
        ms.setFlag(QgsMapSettings.Antialiasing, False)
        ms.setFlag(QgsMapSettings.HighQualityImageTransforms, False)
        
        # Ensure DEM uses Nearest Neighbour for final render to avoid interpolation colors
        dem.dataProvider().setZoomedInResamplingMethod(QgsRasterDataProvider.ResamplingMethod.Nearest)
        dem.dataProvider().setZoomedOutResamplingMethod(QgsRasterDataProvider.ResamplingMethod.Nearest)

        job = QgsMapRendererParallelJob(ms)
        job.start(); job.waitForFinished()
        img = job.renderedImage()
        if img.isNull():
            raise QgsProcessingException("Renderer produced a null image.")
        
        out_png = os.path.join(map_folder, f"{base_name}.png")
        if not img.save(out_png, "png"):
            raise QgsProcessingException("Failed to save PNG.")
        feedback.pushInfo(f"Styled PNG written: {out_png}")

        # ---------------- 5) Write companion JSON --------------------------------
        out_prov_json = os.path.join(map_folder, f"{base_name}.json")
        
        out_obj = {
            "image": {"path": out_png, "width_px": width_px, "height_px": height_px, "dpi": dpi},
            "crs": out_crs.authid(),
            "extent": {"xmin": xmin, "ymin": ymin, "xmax": xmax, "ymax": ymax},
            "origin": "top-left",
            "points": points
        }
        with open(out_prov_json, "w", encoding="utf-8") as fjson:
            json.dump(out_obj, fjson, indent=2)
        feedback.pushInfo(f"Province centers JSON written: {out_prov_json} (count: {len(points)})")

        # Cleanup temporary files
        _cleanup_temp_files(feedback)

        return { "PNG": out_png, "PROVINCE_JSON": out_prov_json }

    # Some QGIS loaders call a method on the class:
    def createInstance(self):
        return MapExporterDemRiversProvPixels()

# Some QGIS loaders call a module-level function:
def createInstance():
    return MapExporterDemRiversProvPixels()
