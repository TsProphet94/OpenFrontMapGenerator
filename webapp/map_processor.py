"""
Map Processor - DEM Processing without QGIS
============================================
Uses GDAL/Rasterio to replicate the QGIS map generation workflow.
"""

import os
import json
import math
import tempfile
import requests
import zipfile
import numpy as np
from PIL import Image
from io import BytesIO
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.merge import merge
from rasterio.mask import mask
from pyproj import Transformer
from shapely.geometry import box, shape, Point
from shapely.ops import transform


# Game file constants
MIN_ISLAND_SIZE = 60
MIN_LAKE_SIZE = 200
TYPE_LAND = 0
TYPE_WATER = 1


# OpenFront Color Palette (Discrete)
DEM_COLOR_RAMP = [
    (0.0, (0, 0, 106), "Water"),           # #00006a
    (30.0, (190, 220, 140), "Plains"),      # #bedc8c
    (60.0, (190, 218, 142), "Plains"),      # #beda8e
    (90.0, (190, 216, 144), "Plains"),      # #bed890
    (120.0, (190, 214, 146), "Plains"),     # #bed692
    (150.0, (190, 212, 148), "Plains"),     # #bed494
    (180.0, (190, 210, 150), "Plains"),     # #bed296
    (210.0, (190, 208, 152), "Plains"),     # #bed098
    (240.0, (190, 206, 154), "Plains"),     # #bece9a
    (270.0, (190, 204, 156), "Plains"),     # #becc9c
    (300.0, (190, 202, 158), "Plains"),     # #beca9e
    (420.0, (220, 203, 160), "Highlands"),  # #dccba0
    (540.0, (222, 205, 162), "Highlands"),  # #decda2
    (660.0, (224, 207, 164), "Highlands"),  # #e0cfa4
    (780.0, (226, 209, 166), "Highlands"),  # #e2d1a6
    (900.0, (228, 211, 168), "Highlands"),  # #e4d3a8
    (1020.0, (230, 213, 170), "Highlands"), # #e6d5aa
    (1140.0, (232, 215, 172), "Highlands"), # #e8d7ac
    (1260.0, (234, 217, 174), "Highlands"), # #ead9ae
    (1380.0, (236, 219, 176), "Highlands"), # #ecdbb0
    (1500.0, (238, 221, 178), "Highlands"), # #eeddb2
    (1818.0, (240, 240, 180), "Mountains"), # #f0f0b4
    (2136.0, (240, 240, 182), "Mountains"), # #f0f0b6
    (2455.0, (241, 241, 184), "Mountains"), # #f1f1b8
    (2773.0, (242, 242, 186), "Mountains"), # #f2f2ba
    (3091.0, (242, 242, 188), "Mountains"), # #f2f2bc
    (3409.0, (242, 242, 190), "Mountains"), # #f2f2be
    (3727.0, (243, 243, 192), "Mountains"), # #f3f3c0
    (4045.0, (244, 244, 194), "Mountains"), # #f4f4c2
    (4364.0, (244, 244, 196), "Mountains"), # #f4f4c4
    (4682.0, (244, 244, 198), "Mountains"), # #f4f4c6
    (5000.0, (245, 245, 200), "Mountains"), # #f5f5c8
]

# Natural Earth data URLs
NE_RIVERS_URL = "https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_rivers_lake_centerlines.zip"
NE_LAKES_URL = "https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_lakes.zip"
NE_ADMIN1_URL = "https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_1_states_provinces.zip"
NE_ADMIN0_URL = "https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_countries.zip"


class MapProcessor:
    """Processes DEM data and generates styled terrain maps."""
    
    def __init__(self, api_key: str, output_dir: str):
        self.api_key = api_key
        self.output_dir = output_dir
        self.cache_dir = os.path.join(tempfile.gettempdir(), 'openfront_cache')
        os.makedirs(self.cache_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
    
    def generate(self, name: str, south: float, west: float, north: float, east: float,
                 width_px: int, height_px: int, dem_source: str = 'COP90') -> dict:
        """
        Generate a styled terrain map.
        
        Args:
            name: Map name
            south, west, north, east: Bounding box in WGS84
            dem_source: DEM source ('COP30', 'COP90', 'SRTM15+')
        
        Returns:
            dict with file paths and metadata
        """
        print(f"Generating map: {name}")
        print(f"Bounds: S={south}, W={west}, N={north}, E={east}")
        print(f"DEM Source: {dem_source}")
        print(f"Output size: {width_px} x {height_px} px (total: {width_px * height_px:,} px)")
        
        # Step 1: Download DEM
        print("Downloading DEM...")
        dem_path = self._download_dem(south, west, north, east, dem_source)
        
        # Step 2: Load and process DEM
        print("Processing DEM...")
        dem_array, dem_transform, dem_crs = self._load_dem(dem_path, width_px, height_px)
        
        # Step 3: Apply color palette
        print("Applying color palette...")
        styled_image = self._apply_palette(dem_array, dynamic_scale=True)
        
        # Step 4: Download and overlay water features
        print("Adding water features...")
        styled_image = self._add_water_features(
            styled_image, south, west, north, east, width_px, height_px
        )
        
        # Step 5: Get province/country points (deactivated for low-memory testing)
        print("[Skipping province/country points for testing]")
        points = []
        
        # Step 6: Save outputs
        print("Saving outputs...")
        base_name = name.lower().replace(' ', '')
        
        # Save PNG as image.png (game format)
        png_path = os.path.join(self.output_dir, "image.png")
        styled_image.save(png_path, 'PNG')
        
        # Save JSON with metadata
        json_path = os.path.join(self.output_dir, f"{base_name}.json")
        metadata = {
            "image": {
                "path": "image.png",
                "width_px": width_px,
                "height_px": height_px
            },
            "bounds": {
                "south": south,
                "west": west,
                "north": north,
                "east": east
            },
            "origin": "top-left",
            "points": points
        }
        with open(json_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # Step 7: Generate game files (bin files, manifest, thumbnail)
        print("Generating game files...")
        game_files = self._generate_game_files(styled_image, base_name, points)
        
        # Clean up temp DEM file
        if os.path.exists(dem_path):
            os.remove(dem_path)
        
        print(f"Map generated successfully!")
        
        all_files = ["image.png", f"{base_name}.json"] + game_files
        
        return {
            'files': all_files,
            'metadata': metadata
        }
    
    def _download_dem(self, south: float, west: float, north: float, east: float,
                      dem_source: str) -> str:
        """Download DEM from OpenTopography."""
        
        demtype_map = {
            'COP30': 'COP30',
            'COP90': 'COP90',
            'SRTM15+': 'SRTMGL1'
        }
        demtype = demtype_map.get(dem_source, 'COP90')
        
        url = "https://portal.opentopography.org/API/globaldem"
        params = {
            'demtype': demtype,
            'south': south,
            'north': north,
            'west': west,
            'east': east,
            'outputFormat': 'GTiff'
        }
        
        if self.api_key:
            params['API_Key'] = self.api_key
        
        print(f"Requesting DEM from OpenTopography...")
        response = requests.get(url, params=params, timeout=300)
        
        if response.status_code != 200:
            raise Exception(f"Failed to download DEM: {response.status_code} - {response.text[:200]}")
        
        # Save to temp file
        dem_path = os.path.join(self.cache_dir, f"dem_{south}_{west}_{north}_{east}.tif")
        with open(dem_path, 'wb') as f:
            f.write(response.content)
        
        return dem_path
    
    def _load_dem(self, dem_path: str, target_width: int, target_height: int) -> tuple:
        """Load DEM and resample to target size."""
        
        with rasterio.open(dem_path) as src:
            # Calculate new transform for target size
            transform_new, width_new, height_new = calculate_default_transform(
                src.crs, src.crs,
                src.width, src.height,
                *src.bounds,
                dst_width=target_width,
                dst_height=target_height
            )
            
            # Resample
            data = np.empty((target_height, target_width), dtype=np.float32)
            reproject(
                source=rasterio.band(src, 1),
                destination=data,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform_new,
                dst_crs=src.crs,
                resampling=Resampling.bilinear
            )
            
            return data, transform_new, src.crs
    
    def _apply_palette(self, dem_array: np.ndarray, dynamic_scale: bool = True) -> Image.Image:
        """Apply the OpenFront color palette to DEM data."""
        
        height, width = dem_array.shape
        
        # Create RGBA image
        rgba = np.zeros((height, width, 4), dtype=np.uint8)
        rgba[:, :, 3] = 255  # Fully opaque
        
        # Get ramp values
        ramp = DEM_COLOR_RAMP
        
        # Dynamic scaling if requested
        if dynamic_scale:
            max_elevation = np.nanmax(dem_array[dem_array > 0])
            if max_elevation > 0 and max_elevation < 5000:
                scale = max_elevation / 5000.0
                ramp = []
                for val, color, label in DEM_COLOR_RAMP:
                    if val <= 0:
                        ramp.append((val, color, label))
                    else:
                        ramp.append((val * scale, color, label))
        
        # Apply discrete color ramp
        for i in range(len(ramp)):
            threshold = ramp[i][0]
            color = ramp[i][1]
            
            if i == 0:
                # Water (everything <= 0)
                mask = dem_array <= threshold
            elif i == len(ramp) - 1:
                # Last color (everything above previous threshold)
                prev_threshold = ramp[i-1][0]
                mask = dem_array > prev_threshold
            else:
                prev_threshold = ramp[i-1][0]
                mask = (dem_array > prev_threshold) & (dem_array <= threshold)
            
            rgba[mask, 0] = color[0]  # R
            rgba[mask, 1] = color[1]  # G
            rgba[mask, 2] = color[2]  # B
        
        return Image.fromarray(rgba, 'RGBA')
    
    def _add_water_features(self, image: Image.Image, south: float, west: float,
                            north: float, east: float, width: int, height: int) -> Image.Image:
        """Add rivers and lakes to the image."""
        
        # Water color from palette
        water_color = (0, 0, 106, 255)  # Ocean blue
        
        # Try to get Natural Earth data
        try:
            # Download and cache rivers/lakes shapefiles
            rivers_path = self._get_ne_shapefile(NE_RIVERS_URL, 'rivers')
            lakes_path = self._get_ne_shapefile(NE_LAKES_URL, 'lakes')
            
            if rivers_path or lakes_path:
                from PIL import ImageDraw
                draw = ImageDraw.Draw(image)
                
                # Coordinate transform function
                def world_to_pixel(lon, lat):
                    px = int((lon - west) / (east - west) * width)
                    py = int((north - lat) / (north - south) * height)
                    return (px, py)
                
                # Draw lakes
                if lakes_path:
                    self._draw_polygons(draw, lakes_path, south, west, north, east,
                                       world_to_pixel, water_color[:3])
                
                # Draw rivers
                if rivers_path:
                    river_width = max(1, int(width / 1000))
                    self._draw_lines(draw, rivers_path, south, west, north, east,
                                    world_to_pixel, water_color[:3], river_width)
        
        except Exception as e:
            print(f"Warning: Could not add water features: {e}")
        
        return image
    
    def _get_ne_shapefile(self, url: str, name: str) -> str:
        """Download and cache a Natural Earth shapefile."""
        
        cache_dir = os.path.join(self.cache_dir, 'natural_earth', name)
        
        # Check if already cached
        if os.path.exists(cache_dir):
            for f in os.listdir(cache_dir):
                if f.endswith('.shp'):
                    return os.path.join(cache_dir, f)
        
        # Download
        try:
            os.makedirs(cache_dir, exist_ok=True)
            print(f"Downloading {name} from Natural Earth...")
            response = requests.get(url, timeout=120)
            
            if response.status_code == 200:
                # Extract ZIP
                zip_path = os.path.join(cache_dir, f'{name}.zip')
                with open(zip_path, 'wb') as f:
                    f.write(response.content)
                
                with zipfile.ZipFile(zip_path, 'r') as zf:
                    zf.extractall(cache_dir)
                
                os.remove(zip_path)
                
                # Find .shp file
                for f in os.listdir(cache_dir):
                    if f.endswith('.shp'):
                        return os.path.join(cache_dir, f)
        
        except Exception as e:
            print(f"Warning: Could not download {name}: {e}")
        
        return None
    
    def _draw_polygons(self, draw, shp_path: str, south: float, west: float,
                       north: float, east: float, world_to_pixel, color):
        """Draw polygon features from a shapefile."""
        try:
            import fiona
            
            bounds_box = box(west, south, east, north)
            
            with fiona.open(shp_path, 'r') as src:
                for feature in src:
                    try:
                        geom = shape(feature['geometry'])
                        if not geom.intersects(bounds_box):
                            continue
                        
                        clipped = geom.intersection(bounds_box)
                        if clipped.is_empty:
                            continue
                        
                        # Draw polygon(s)
                        self._draw_geometry(draw, clipped, world_to_pixel, color)
                    except Exception:
                        continue
        except ImportError:
            print("Warning: fiona not available for shapefile reading")
        except Exception as e:
            print(f"Warning: Error reading shapefile: {e}")
    
    def _draw_lines(self, draw, shp_path: str, south: float, west: float,
                    north: float, east: float, world_to_pixel, color, width: int):
        """Draw line features from a shapefile."""
        try:
            import fiona
            
            bounds_box = box(west, south, east, north)
            
            with fiona.open(shp_path, 'r') as src:
                for feature in src:
                    try:
                        # Filter by scalerank if available
                        props = feature.get('properties', {})
                        scalerank = props.get('scalerank', 0)
                        if scalerank is not None and scalerank > 3:
                            continue
                        
                        geom = shape(feature['geometry'])
                        if not geom.intersects(bounds_box):
                            continue
                        
                        clipped = geom.intersection(bounds_box)
                        if clipped.is_empty:
                            continue
                        
                        # Draw line(s)
                        self._draw_line_geometry(draw, clipped, world_to_pixel, color, width)
                    except Exception:
                        continue
        except ImportError:
            print("Warning: fiona not available for shapefile reading")
        except Exception as e:
            print(f"Warning: Error reading shapefile: {e}")
    
    def _draw_geometry(self, draw, geom, world_to_pixel, color):
        """Draw a shapely geometry as a filled polygon."""
        if geom.geom_type == 'Polygon':
            coords = [world_to_pixel(x, y) for x, y in geom.exterior.coords]
            if len(coords) >= 3:
                draw.polygon(coords, fill=color)
        elif geom.geom_type == 'MultiPolygon':
            for poly in geom.geoms:
                self._draw_geometry(draw, poly, world_to_pixel, color)
    
    def _draw_line_geometry(self, draw, geom, world_to_pixel, color, width):
        """Draw a shapely geometry as a line."""
        if geom.geom_type == 'LineString':
            coords = [world_to_pixel(x, y) for x, y in geom.coords]
            if len(coords) >= 2:
                draw.line(coords, fill=color, width=width)
        elif geom.geom_type == 'MultiLineString':
            for line in geom.geoms:
                self._draw_line_geometry(draw, line, world_to_pixel, color, width)
    
    def _get_province_points(self, south: float, west: float, north: float, east: float,
                             width: int, height: int, max_provinces: int = 20) -> list:
        """Get province/country center points."""
        
        points = []
        
        def world_to_pixel(lon, lat):
            px = int((lon - west) / (east - west) * width)
            py = int((north - lat) / (north - south) * height)
            return (px, py)
        
        try:
            # Try countries first (for world/continent maps)
            admin0_path = self._get_ne_shapefile(NE_ADMIN0_URL, 'admin0')
            if admin0_path:
                country_points = self._extract_points_from_shapefile(
                    admin0_path, south, west, north, east, world_to_pixel,
                    name_fields=['NAME', 'NAME_EN', 'ADMIN'],
                    flag_fields=['ISO_A2', 'ISO_A2_EH'],
                    is_country=True
                )
                
                if len(country_points) >= 3:
                    # Use countries
                    points = sorted(country_points, key=lambda x: x.get('area', 0), reverse=True)
                    points = points[:max_provinces]
        except Exception as e:
            print(f"Warning: Could not load countries: {e}")
        
        # Fall back to provinces/states
        if len(points) < 3:
            try:
                admin1_path = self._get_ne_shapefile(NE_ADMIN1_URL, 'admin1')
                if admin1_path:
                    province_points = self._extract_points_from_shapefile(
                        admin1_path, south, west, north, east, world_to_pixel,
                        name_fields=['name', 'name_en', 'name_1'],
                        flag_fields=['iso_3166_2', 'adm1_code'],
                        is_country=False
                    )
                    points = sorted(province_points, key=lambda x: x.get('area', 0), reverse=True)
                    points = points[:max_provinces]
            except Exception as e:
                print(f"Warning: Could not load provinces: {e}")
        
        # Clean up points for output
        cleaned_points = []
        for p in points:
            cleaned_points.append({
                'name': p['name'],
                'flag': p.get('flag', 'xx'),
                'pixel_x': p['pixel_x'],
                'pixel_y': p['pixel_y']
            })
        
        return cleaned_points
    
    def _extract_points_from_shapefile(self, shp_path: str, south: float, west: float,
                                        north: float, east: float, world_to_pixel,
                                        name_fields: list, flag_fields: list,
                                        is_country: bool) -> list:
        """Extract center points from a shapefile."""
        
        points = []
        
        try:
            import fiona
            
            bounds_box = box(west, south, east, north)
            
            with fiona.open(shp_path, 'r') as src:
                for feature in src:
                    try:
                        geom = shape(feature['geometry'])
                        if not geom.intersects(bounds_box):
                            continue
                        
                        clipped = geom.intersection(bounds_box)
                        if clipped.is_empty:
                            continue
                        
                        # Get center point
                        if hasattr(clipped, 'representative_point'):
                            center = clipped.representative_point()
                        else:
                            center = clipped.centroid
                        
                        if center.is_empty:
                            continue
                        
                        # Get name
                        props = feature.get('properties', {})
                        name = None
                        for field in name_fields:
                            val = props.get(field)
                            if val:
                                name = str(val)
                                break
                        
                        if not name:
                            continue
                        
                        # Get flag code
                        flag = 'xx'
                        for field in flag_fields:
                            val = props.get(field)
                            if val:
                                flag = str(val).lower().split('-')[0]
                                break
                        
                        # Calculate pixel coords
                        px, py = world_to_pixel(center.x, center.y)
                        
                        # Skip if outside image bounds
                        if px < 0 or px >= world_to_pixel(east, south)[0]:
                            continue
                        if py < 0 or py >= world_to_pixel(west, south)[1]:
                            continue
                        
                        # Estimate area
                        area = clipped.area
                        
                        points.append({
                            'name': name,
                            'flag': flag,
                            'pixel_x': px,
                            'pixel_y': py,
                            'area': area
                        })
                    
                    except Exception:
                        continue
        
        except ImportError:
            print("Warning: fiona not available")
        except Exception as e:
            print(f"Warning: Error extracting points: {e}")
        
        return points

    # =========================================================================
    # Game File Generation (map.bin, map4x.bin, map16x.bin, manifest.json, thumbnail.webp)
    # =========================================================================
    
    def _generate_game_files(self, styled_image: Image.Image, base_name: str, points: list) -> list:
        """
        Generate OpenFront game files from the styled image.
        
        Returns:
            List of generated file names
        """
        # Convert image to numpy array
        img = styled_image.convert('RGBA')
        width, height = img.size
        
        # Ensure even dimensions
        width -= width % 2
        height -= height % 2
        img = img.crop((0, 0, width, height))
        
        pixels = np.array(img)
        
        # Extract channels
        r = pixels[:, :, 0]
        g = pixels[:, :, 1]
        b = pixels[:, :, 2]
        a = pixels[:, :, 3]
        
        # Determine Type: Alpha < 20 or Blue == 106 -> Water
        terrain_type = np.where((a < 20) | (b == 106), TYPE_WATER, TYPE_LAND).astype(np.uint8)
        
        # Determine Magnitude: (Blue - 140) / 2, clamped 0-30
        mag_raw = np.clip(b.astype(float), 140, 200) - 140
        terrain_mag = mag_raw / 2.0
        
        terrain_shore = np.zeros((height, width), dtype=bool)
        terrain_ocean = np.zeros((height, width), dtype=bool)
        
        # Remove small islands
        terrain_type, terrain_mag = self._remove_small_areas(
            terrain_type, terrain_mag, TYPE_LAND, MIN_ISLAND_SIZE, TYPE_WATER
        )
        
        # Process water (identify oceans, remove small lakes, calc distances)
        terrain_type, terrain_mag, terrain_shore, terrain_ocean = self._process_water(
            terrain_type, terrain_mag, terrain_shore, terrain_ocean
        )
        
        # Create downscaled maps
        # L1 (1/2)
        l1_type, l1_mag, l1_shore, l1_ocean = self._downscale_terrain(
            terrain_type, terrain_mag, terrain_shore, terrain_ocean
        )
        
        # L2 (1/4)
        l2_type, l2_mag, l2_shore, l2_ocean = self._downscale_terrain(
            l1_type, l1_mag, l1_shore, l1_ocean
        )
        
        # L3 (1/8)
        l3_type, l3_mag, l3_shore, l3_ocean = self._downscale_terrain(
            l2_type, l2_mag, l2_shore, l2_ocean
        )
        
        # Pack data
        l0_data, l0_land = self._pack_terrain(terrain_type, terrain_mag, terrain_shore, terrain_ocean)
        l1_data, l1_land = self._pack_terrain(l1_type, l1_mag, l1_shore, l1_ocean)
        l2_data, l2_land = self._pack_terrain(l2_type, l2_mag, l2_shore, l2_ocean)
        l3_data, l3_land = self._pack_terrain(l3_type, l3_mag, l3_shore, l3_ocean)
        
        # Save game files (use small version: L1, L2, L3)
        generated_files = []
        
        # Write binary files
        with open(os.path.join(self.output_dir, "map.bin"), "wb") as f:
            f.write(l1_data)
        generated_files.append("map.bin")
        
        with open(os.path.join(self.output_dir, "map4x.bin"), "wb") as f:
            f.write(l2_data)
        generated_files.append("map4x.bin")
        
        with open(os.path.join(self.output_dir, "map16x.bin"), "wb") as f:
            f.write(l3_data)
        generated_files.append("map16x.bin")
        
        # Generate and save thumbnail
        thumbnail = self._create_thumbnail(l2_type, l2_mag, l2_shore, 0.5)
        thumbnail.save(os.path.join(self.output_dir, "thumbnail.webp"), "WEBP")
        generated_files.append("thumbnail.webp")
        
        # Create manifest
        l1_h, l1_w = l1_type.shape
        l2_h, l2_w = l2_type.shape
        l3_h, l3_w = l3_type.shape
        
        manifest = {
            "name": base_name,
            "map": {
                "width": width,
                "height": height,
                "num_land_tiles": l1_land
            },
            "map4x": {
                "width": l2_w,
                "height": l2_h,
                "num_land_tiles": l2_land
            },
            "map16x": {
                "width": l3_w,
                "height": l3_h,
                "num_land_tiles": l3_land
            },
            "nations": []
        }
        
        # Add nations from points
        for p in points:
            nation = {
                "name": p.get("name", "Unknown"),
                "flag": p.get("flag", "unknown"),
                "coordinates": [
                    int(p.get("pixel_x", 0) * 0.5),  # Scale for small version
                    int(p.get("pixel_y", 0) * 0.5)
                ]
            }
            manifest["nations"].append(nation)
        
        with open(os.path.join(self.output_dir, "manifest.json"), "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)
        generated_files.append("manifest.json")
        
        return generated_files
    
    def _remove_small_areas(self, t_type, t_mag, target_type, min_size, replace_with):
        """Remove small islands or lakes."""
        try:
            from scipy.ndimage import label
            mask = (t_type == target_type)
            labeled, n_components = label(mask)
            
            sizes = np.bincount(labeled.ravel())
            small_labels = np.where((sizes < min_size) & (sizes > 0))[0]
            
            remove_mask = np.isin(labeled, small_labels)
            t_type = t_type.copy()
            t_mag = t_mag.copy()
            t_type[remove_mask] = replace_with
            t_mag[remove_mask] = 0
            
            print(f"Removed {len(small_labels)} areas smaller than {min_size}")
            
        except ImportError:
            print("scipy not found, skipping small area removal")
        
        return t_type, t_mag
    
    def _process_water(self, t_type, t_mag, t_shore, t_ocean):
        """Process water bodies - identify oceans, remove small lakes, calc distances."""
        try:
            from scipy.ndimage import label, distance_transform_cdt, binary_dilation
            
            water_mask = (t_type == TYPE_WATER)
            labeled, n_components = label(water_mask)
            
            if n_components == 0:
                print("No water bodies found.")
                return t_type, t_mag, t_shore, t_ocean
            
            sizes = np.bincount(labeled.ravel())
            water_labels = np.arange(1, len(sizes))
            water_labels = water_labels[np.argsort(sizes[water_labels])[::-1]]
            
            # Largest is Ocean
            t_type = t_type.copy()
            t_mag = t_mag.copy()
            t_shore = t_shore.copy()
            t_ocean = t_ocean.copy()
            
            largest_label = water_labels[0]
            t_ocean[labeled == largest_label] = True
            print(f"Identified ocean with {sizes[largest_label]} tiles")
            
            # Remove small lakes
            small_labels = water_labels[sizes[water_labels] < MIN_LAKE_SIZE]
            remove_mask = np.isin(labeled, small_labels)
            t_type[remove_mask] = TYPE_LAND
            t_mag[remove_mask] = 0
            print(f"Removed {len(small_labels)} lakes smaller than {MIN_LAKE_SIZE}")
            
            water_mask = (t_type == TYPE_WATER)
            
            # Shoreline detection
            struct_4 = np.array([[0,1,0],[1,1,1],[0,1,0]], dtype=bool)
            land_mask = (t_type == TYPE_LAND)
            
            dilated_land = binary_dilation(land_mask, structure=struct_4)
            shore_water = dilated_land & water_mask
            
            dilated_water = binary_dilation(water_mask, structure=struct_4)
            shore_land = dilated_water & land_mask
            
            t_shore = shore_water | shore_land
            
            # Water magnitude (distance to land)
            dist = distance_transform_cdt(water_mask, metric='taxicab')
            water_mag = np.maximum(dist - 1, 0).astype(float)
            t_mag = np.where(water_mask, water_mag, t_mag)
            
        except ImportError:
            print("scipy not found, skipping water processing")
        
        return t_type, t_mag, t_shore, t_ocean
    
    def _downscale_terrain(self, t_type, t_mag, t_shore, t_ocean):
        """Downscale terrain by factor of 2."""
        h, w = t_type.shape
        if h % 2 != 0:
            h -= 1
        if w % 2 != 0:
            w -= 1
        
        t_type = t_type[:h, :w]
        t_mag = t_mag[:h, :w]
        t_shore = t_shore[:h, :w]
        t_ocean = t_ocean[:h, :w]
        
        # Get 2x2 blocks
        s00_type = t_type[0::2, 0::2]
        s01_type = t_type[1::2, 0::2]
        s10_type = t_type[0::2, 1::2]
        s11_type = t_type[1::2, 1::2]
        
        w00 = (s00_type == TYPE_WATER)
        w01 = (s01_type == TYPE_WATER)
        w10 = (s10_type == TYPE_WATER)
        
        # Start with s11
        mini_type = s11_type.copy()
        mini_mag = t_mag[1::2, 1::2].copy()
        mini_shore = t_shore[1::2, 1::2].copy()
        mini_ocean = t_ocean[1::2, 1::2].copy()
        
        # Priority: P00 > P01 > P10 > P11 (for Water)
        # P10
        mini_type[w10] = t_type[0::2, 1::2][w10]
        mini_mag[w10] = t_mag[0::2, 1::2][w10]
        mini_shore[w10] = t_shore[0::2, 1::2][w10]
        mini_ocean[w10] = t_ocean[0::2, 1::2][w10]
        
        # P01
        mini_type[w01] = t_type[1::2, 0::2][w01]
        mini_mag[w01] = t_mag[1::2, 0::2][w01]
        mini_shore[w01] = t_shore[1::2, 0::2][w01]
        mini_ocean[w01] = t_ocean[1::2, 0::2][w01]
        
        # P00
        mini_type[w00] = t_type[0::2, 0::2][w00]
        mini_mag[w00] = t_mag[0::2, 0::2][w00]
        mini_shore[w00] = t_shore[0::2, 0::2][w00]
        mini_ocean[w00] = t_ocean[0::2, 0::2][w00]
        
        return mini_type, mini_mag, mini_shore, mini_ocean
    
    def _pack_terrain(self, t_type, t_mag, t_shore, t_ocean):
        """Pack terrain into bytes."""
        # Bit 7: Land (1) / Water (0)
        # Bit 6: Shoreline
        # Bit 5: Ocean
        # Bits 0-4: Magnitude
        
        mag_byte = np.where(
            t_type == TYPE_LAND,
            np.minimum(np.ceil(t_mag), 31),
            np.minimum(np.ceil(t_mag / 2), 31)
        ).astype(np.uint8)
        
        packed = np.zeros_like(t_type, dtype=np.uint8)
        packed |= (t_type == TYPE_LAND).astype(np.uint8) << 7
        packed |= t_shore.astype(np.uint8) << 6
        packed |= t_ocean.astype(np.uint8) << 5
        packed |= mag_byte & 0x1F
        
        num_land = int(np.sum(t_type == TYPE_LAND))
        
        return packed.tobytes(), num_land
    
    def _create_thumbnail(self, t_type, t_mag, t_shore, quality):
        """Create thumbnail image from terrain data."""
        src_h, src_w = t_type.shape
        target_w = int(max(1, math.floor(src_w * quality)))
        target_h = int(max(1, math.floor(src_h * quality)))
        
        img = Image.new('RGBA', (target_w, target_h))
        pixels = img.load()
        
        for x in range(target_w):
            for y in range(target_h):
                src_x = int(min(math.floor(x / quality), src_w - 1))
                src_y = int(min(math.floor(y / quality), src_h - 1))
                
                tile_type = t_type[src_y, src_x]
                tile_mag = t_mag[src_y, src_x]
                tile_shore = t_shore[src_y, src_x]
                
                color = self._get_thumbnail_color(tile_type, tile_mag, tile_shore)
                pixels[x, y] = color
        
        return img
    
    def _get_thumbnail_color(self, t_type, t_mag, t_shore):
        """Get color for thumbnail pixel."""
        if t_type == TYPE_WATER:
            if t_shore:
                return (100, 143, 255, 0)
            
            water_adj = 11 - min(t_mag / 2, 10) - 10
            r = int(max(70 + water_adj, 0))
            g = int(max(132 + water_adj, 0))
            b = int(max(180 + water_adj, 0))
            return (r, g, b, 0)
        
        # Land
        if t_shore:
            return (204, 203, 158, 255)
        
        if t_mag < 10:
            adj = 220 - 2 * t_mag
            return (190, int(adj), 138, 255)
        elif t_mag < 20:
            adj = 2 * t_mag
            return (int(200 + adj), int(183 + adj), int(138 + adj), 255)
        else:
            adj = int(230 + t_mag / 2)
            return (adj, adj, adj, 255)


if __name__ == '__main__':
    # Test the processor
    processor = MapProcessor(
        api_key='',  # Add your key for testing
        output_dir='./test_output'
    )
    
    result = processor.generate(
        name='Cyprus Test',
        south=34.5,
        west=32.0,
        north=35.7,
        east=34.6,
        width_px=1024,
        dem_source='COP90'
    )
    
    print(result)
