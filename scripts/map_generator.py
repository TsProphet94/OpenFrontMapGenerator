import os
import json
import argparse
import numpy as np
from PIL import Image
import math
from collections import deque
import struct

# Constants
MIN_ISLAND_SIZE = 60
MIN_LAKE_SIZE = 200

# Terrain Types
TYPE_LAND = 0
TYPE_WATER = 1

class MapGenerator:
    def __init__(self, map_name, input_dir, output_dir, is_test=False):
        self.map_name = map_name
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.is_test = is_test
        
        self.image_path = os.path.join(input_dir, f"{map_name}.png")
        self.json_path = os.path.join(input_dir, f"{map_name}.json")
        
        # If input files don't exist with map_name, try generic names (like MapGenerator structure)
        if not os.path.exists(self.image_path):
            self.image_path = os.path.join(input_dir, "image.png")
        if not os.path.exists(self.json_path):
            self.json_path = os.path.join(input_dir, "info.json")

        self.width = 0
        self.height = 0
        self.terrain_type = None
        self.terrain_mag = None
        self.terrain_shore = None
        self.terrain_ocean = None

    def process(self):
        print(f"Processing map: {self.map_name}")
        
        # 1. Load Data
        img = Image.open(self.image_path).convert('RGBA')
        self.width, self.height = img.size
        
        # Ensure even dimensions
        self.width -= self.width % 2
        self.height -= self.height % 2
        img = img.crop((0, 0, self.width, self.height))
        
        pixels = np.array(img)
        
        # 2. Initialize Terrain
        # Extract channels
        r = pixels[:, :, 0]
        g = pixels[:, :, 1]
        b = pixels[:, :, 2]
        a = pixels[:, :, 3]
        
        # Determine Type: Alpha < 20 or Blue == 106 -> Water
        self.terrain_type = np.where((a < 20) | (b == 106), TYPE_WATER, TYPE_LAND).astype(np.uint8)
        
        # Determine Magnitude: (Blue - 140) / 2, clamped 0-30
        # Note: Go code uses float64 for magnitude, but packs it as byte.
        # mag = math.Min(200, math.Max(140, float64(blue))) - 140
        # terrain[x][y].Magnitude = mag / 2
        mag_raw = np.clip(b.astype(float), 140, 200) - 140
        self.terrain_mag = mag_raw / 2.0
        
        self.terrain_shore = np.zeros((self.height, self.width), dtype=bool)
        self.terrain_ocean = np.zeros((self.height, self.width), dtype=bool)

        # 3. Remove Small Islands
        if not self.is_test:
            self.remove_small_areas(TYPE_LAND, MIN_ISLAND_SIZE, replace_with=TYPE_WATER)

        # 4. Process Water (Identify Oceans, Remove Small Lakes, Calc Distances)
        self.process_water(remove_small=not self.is_test)

        # 5. Create Downscaled Maps
        # L1 (1/2)
        l1_type, l1_mag, l1_shore, l1_ocean = self.downscale_terrain(
            self.terrain_type, self.terrain_mag, self.terrain_shore, self.terrain_ocean)
            
        # L2 (1/4)
        l2_type, l2_mag, l2_shore, l2_ocean = self.downscale_terrain(
            l1_type, l1_mag, l1_shore, l1_ocean)

        # L3 (1/8) - New for Small version's map16x
        l3_type, l3_mag, l3_shore, l3_ocean = self.downscale_terrain(
            l2_type, l2_mag, l2_shore, l2_ocean)

        # 7. Pack Data
        # L0
        l0_data, l0_land = self.pack_terrain(self.terrain_type, self.terrain_mag, self.terrain_shore, self.terrain_ocean)
        # L1
        l1_data, l1_land = self.pack_terrain(l1_type, l1_mag, l1_shore, l1_ocean)
        # L2
        l2_data, l2_land = self.pack_terrain(l2_type, l2_mag, l2_shore, l2_ocean)
        # L3
        l3_data, l3_land = self.pack_terrain(l3_type, l3_mag, l3_shore, l3_ocean)

        # 8. Save Outputs
        # Normalize map name: lowercase and remove spaces
        clean_name = self.map_name.lower().replace(" ", "")
        
        # Big Version (L0, L1, L2)
        self.save_variant(clean_name + "big", 
                          l0_data, l1_data, l2_data,
                          l0_land, l1_land, l2_land,
                          (self.terrain_type, l1_type, l2_type),  # Terrain types for dimension calculation
                          (l1_type, l1_mag, l1_shore), # Thumbnail source (map4x of Big)
                          1.0)
                          
        # Small Version (L1, L2, L3)
        self.save_variant(clean_name,
                          l1_data, l2_data, l3_data,
                          l1_land, l2_land, l3_land,
                          (l1_type, l2_type, l3_type),  # Terrain types for dimension calculation
                          (l2_type, l2_mag, l2_shore), # Thumbnail source (map4x of Small)
                          0.5)

    def remove_small_areas(self, target_type, min_size, replace_with):
        # Using scipy.ndimage if available, else custom BFS
        try:
            from scipy.ndimage import label
            mask = (self.terrain_type == target_type)
            labeled, n_components = label(mask)
            
            sizes = np.bincount(labeled.ravel())
            # sizes[0] is background (where mask is False), ignore it
            
            # Find labels to remove
            small_labels = np.where((sizes < min_size) & (sizes > 0))[0] # sizes > 0 check is redundant but safe
            
            # Create a mask of pixels to flip
            remove_mask = np.isin(labeled, small_labels)
            
            # Flip
            self.terrain_type[remove_mask] = replace_with
            self.terrain_mag[remove_mask] = 0
            
            print(f"Removed {len(small_labels)} areas smaller than {min_size}")
            
        except ImportError:
            print("scipy not found, skipping small area removal (install scipy for this feature)")

    def process_water(self, remove_small):
        try:
            from scipy.ndimage import label, distance_transform_cdt
            
            # Identify Water Bodies
            water_mask = (self.terrain_type == TYPE_WATER)
            labeled, n_components = label(water_mask)
            
            if n_components == 0:
                print("No water bodies found.")
                return

            sizes = np.bincount(labeled.ravel())
            # sizes[0] is land
            
            # Get water body labels sorted by size (descending)
            # Skip 0
            water_labels = np.arange(1, len(sizes))
            water_labels = water_labels[np.argsort(sizes[water_labels])[::-1]]
            
            # Largest is Ocean
            largest_label = water_labels[0]
            self.terrain_ocean[labeled == largest_label] = True
            print(f"Identified ocean with {sizes[largest_label]} tiles")
            
            # Remove small lakes
            if remove_small:
                small_labels = water_labels[sizes[water_labels] < MIN_LAKE_SIZE]
                remove_mask = np.isin(labeled, small_labels)
                self.terrain_type[remove_mask] = TYPE_LAND
                self.terrain_mag[remove_mask] = 0
                print(f"Removed {len(small_labels)} lakes smaller than {MIN_LAKE_SIZE}")
                
                # Update water mask after removal
                water_mask = (self.terrain_type == TYPE_WATER)

            # Shoreline Detection
            # Land adjacent to Water OR Water adjacent to Land
            # Using binary dilation
            from scipy.ndimage import binary_dilation
            
            # Structure for 4-connectivity
            struct_4 = np.array([[0,1,0],[1,1,1],[0,1,0]], dtype=bool)
            
            land_mask = (self.terrain_type == TYPE_LAND)
            
            # Dilate Land -> overlap with Water is Shoreline Water
            dilated_land = binary_dilation(land_mask, structure=struct_4)
            shore_water = dilated_land & water_mask
            
            # Dilate Water -> overlap with Land is Shoreline Land
            dilated_water = binary_dilation(water_mask, structure=struct_4)
            shore_land = dilated_water & land_mask
            
            self.terrain_shore = shore_water | shore_land
            
            # Water Magnitude (Distance to Land)
            # distance_transform_cdt calculates distance to nearest 0
            # So we want distance in Water (1) to Land (0)
            # But we want shoreline water to be 0 magnitude?
            # Go code:
            # Shoreline waters added to queue with dist 0.
            # BFS adds 1 per step.
            # So shoreline water mag = 0.
            # Next to shoreline = 1.
            
            # distance_transform_cdt with metric='taxicab'
            # Input: Water=1, Land=0.
            # Result at shoreline water (adj to land) is 1.
            # We want 0. So subtract 1?
            
            dist = distance_transform_cdt(water_mask, metric='taxicab')
            # dist is 1 for shoreline, 2 for next...
            # Go code: shoreline water mag = 0.
            # So we subtract 1, and clip at 0.
            
            water_mag = np.maximum(dist - 1, 0).astype(float)
            
            # Update magnitude for water tiles
            self.terrain_mag = np.where(water_mask, water_mag, self.terrain_mag)
            
        except ImportError:
            print("scipy not found, skipping water processing (install scipy for this feature)")

    def downscale_terrain(self, t_type, t_mag, t_shore, t_ocean):
        # Ensure even dimensions for downscaling
        h, w = t_type.shape
        if h % 2 != 0:
            h -= 1
        if w % 2 != 0:
            w -= 1
            
        # Crop if needed
        if h != t_type.shape[0] or w != t_type.shape[1]:
            t_type = t_type[:h, :w]
            t_mag = t_mag[:h, :w]
            t_shore = t_shore[:h, :w]
            t_ocean = t_ocean[:h, :w]

        # Downscale by 2
        # Logic matches Go:
        # Iterate x, y.
        # If mini_map cell is not Water, overwrite with current.
        # This implies:
        # 1. If any pixel in 2x2 is Water, result is Water.
        # 2. If multiple are Water, the FIRST one visited wins.
        # 3. If all are Land, the LAST one visited wins.
        # Visit order: (0,0), (0,1), (1,0), (1,1) relative to block.
        
        # Slices for the 4 positions
        # Note: numpy is (y, x)
        # P00: x=0, y=0 -> (0, 0)
        # P01: x=0, y=1 -> (1, 0)  <-- Wait. Go loop: x then y.
        # x=2*mx, y=2*my     -> P00
        # x=2*mx, y=2*my+1   -> P01
        # x=2*mx+1, y=2*my   -> P10
        # x=2*mx+1, y=2*my+1 -> P11
        
        # In numpy (y, x):
        # P00: y=0, x=0
        # P01: y=1, x=0
        # P10: y=0, x=1
        # P11: y=1, x=1
        
        s00_type = t_type[0::2, 0::2]
        s01_type = t_type[1::2, 0::2]
        s10_type = t_type[0::2, 1::2]
        s11_type = t_type[1::2, 1::2]
        
        # Masks for Water
        w00 = (s00_type == TYPE_WATER)
        w01 = (s01_type == TYPE_WATER)
        w10 = (s10_type == TYPE_WATER)
        w11 = (s11_type == TYPE_WATER)
        
        # Initialize with S11 (Case: All Land -> Last one wins)
        # Also covers Case: P11 is Water and others are Land.
        mini_type = s11_type.copy()
        mini_mag = t_mag[1::2, 1::2].copy()
        mini_shore = t_shore[1::2, 1::2].copy()
        mini_ocean = t_ocean[1::2, 1::2].copy()
        
        # Apply Priority: P00 > P01 > P10 > P11 (for Water)
        # We apply in reverse order of priority so high priority overwrites.
        
        # 3. P10 (if Water)
        self._update_mini(mini_type, mini_mag, mini_shore, mini_ocean, 
                          t_type[0::2, 1::2], t_mag[0::2, 1::2], 
                          t_shore[0::2, 1::2], t_ocean[0::2, 1::2], w10)

        # 2. P01 (if Water)
        self._update_mini(mini_type, mini_mag, mini_shore, mini_ocean, 
                          t_type[1::2, 0::2], t_mag[1::2, 0::2], 
                          t_shore[1::2, 0::2], t_ocean[1::2, 0::2], w01)

        # 1. P00 (if Water)
        self._update_mini(mini_type, mini_mag, mini_shore, mini_ocean, 
                          t_type[0::2, 0::2], t_mag[0::2, 0::2], 
                          t_shore[0::2, 0::2], t_ocean[0::2, 0::2], w00)
        
        return mini_type, mini_mag, mini_shore, mini_ocean

    def _update_mini(self, m_type, m_mag, m_shore, m_ocean, s_type, s_mag, s_shore, s_ocean, mask):
        m_type[mask] = s_type[mask]
        m_mag[mask] = s_mag[mask]
        m_shore[mask] = s_shore[mask]
        m_ocean[mask] = s_ocean[mask]

    def create_thumbnail(self, terrain_type, terrain_mag, terrain_shore, quality):
        # Create thumbnail from mini-map (or terrain)
        # Go code: createMapThumbnail(miniTerrain, 0.5)
        # So it downscales the mini-map by another 0.5 (total 0.25 of original).
        
        src_h, src_w = terrain_type.shape
        target_w = int(max(1, math.floor(src_w * quality)))
        target_h = int(max(1, math.floor(src_h * quality)))
        
        img = Image.new('RGBA', (target_w, target_h))
        pixels = img.load()
        
        for x in range(target_w):
            for y in range(target_h):
                src_x = int(min(math.floor(x / quality), src_w - 1))
                src_y = int(min(math.floor(y / quality), src_h - 1))
                
                # Get attributes
                # Note: numpy is (row, col) -> (y, x)
                t_type = terrain_type[src_y, src_x]
                t_mag = terrain_mag[src_y, src_x]
                t_shore = terrain_shore[src_y, src_x]
                
                color = self.get_thumbnail_color(t_type, t_mag, t_shore)
                pixels[x, y] = color
                
        return img

    def get_thumbnail_color(self, t_type, t_mag, t_shore):
        if t_type == TYPE_WATER:
            if t_shore:
                return (100, 143, 255, 0) # Alpha 0? Go code says 0.
            
            water_adj = 11 - min(t_mag / 2, 10) - 10
            r = int(max(70 + water_adj, 0))
            g = int(max(132 + water_adj, 0))
            b = int(max(180 + water_adj, 0))
            return (r, g, b, 0)
        
        # Land
        if t_shore:
            return (204, 203, 158, 255)
        
        if t_mag < 10:
            # Plains
            adj = 220 - 2 * t_mag
            return (190, int(adj), 138, 255)
        elif t_mag < 20:
            # Highlands
            adj = 2 * t_mag
            return (int(200 + adj), int(183 + adj), int(138 + adj), 255)
        else:
            # Mountains
            adj = math.floor(230 + t_mag / 2)
            return (int(adj), int(adj), int(adj), 255)

    def pack_terrain(self, t_type, t_mag, t_shore, t_ocean):
        # Pack into bytes
        # Bit 7: Land (1) / Water (0)
        # Bit 6: Shoreline
        # Bit 5: Ocean
        # Bits 0-4: Magnitude
        
        # Prepare arrays
        # Note: numpy arrays are (y, x). We need to flatten in (x, y) order?
        # Go code:
        # for x := 0; x < width; x++ {
        #     for y := 0; y < height; y++ {
        #         packedData[y*width+x] = packedByte
        #     }
        # }
        # Wait, `packedData[y*width+x]` implies row-major order (y changes fastest? No).
        # Usually `y*width + x` is index for (x, y) in a row-major array where y is row, x is col.
        # But the loop is `for x ... for y`.
        # If it fills `y*width + x`, it's filling the array in standard row-major order (row 0, then row 1...).
        # So the order in the file is Row 0 (all x), Row 1 (all x)...
        # So we can just flatten the numpy array (which is row-major by default).
        
        # Calculate magnitude byte
        # Land: min(ceil(mag), 31)
        # Water: min(ceil(mag/2), 31)
        mag_byte = np.where(t_type == TYPE_LAND, 
                            np.minimum(np.ceil(t_mag), 31), 
                            np.minimum(np.ceil(t_mag / 2), 31)).astype(np.uint8)
        
        packed = np.zeros_like(t_type, dtype=np.uint8)
        
        # Set bits
        packed |= (t_type == TYPE_LAND).astype(np.uint8) << 7
        packed |= t_shore.astype(np.uint8) << 6
        packed |= t_ocean.astype(np.uint8) << 5
        packed |= mag_byte & 0x1F
        
        # Count land tiles
        num_land = np.sum(t_type == TYPE_LAND)
        
        return packed.tobytes(), int(num_land)

    def save_variant(self, variant_name, 
                     map_data, map4x_data, map16x_data, 
                     map_land, map4x_land, map16x_land,
                     dimension_terrains,  # tuple of (map_type, map4x_type, map16x_type) for dimension calculation
                     thumbnail_source_terrain, # tuple of (type, mag, shore) to generate thumbnail from
                     scale_factor):
        
        # Determine output directory based on original output_dir
        # self.output_dir is e.g. "generated/maps/United Kingdom"
        # We want "generated/maps/United KingdomBig"
        
        base_output_dir = os.path.dirname(self.output_dir)
        variant_output_dir = os.path.join(base_output_dir, variant_name)
        
        if not os.path.exists(variant_output_dir):
            os.makedirs(variant_output_dir)

        # Write binaries
        with open(os.path.join(variant_output_dir, "map.bin"), "wb") as f:
            f.write(map_data)
        with open(os.path.join(variant_output_dir, "map4x.bin"), "wb") as f:
            f.write(map4x_data)
        with open(os.path.join(variant_output_dir, "map16x.bin"), "wb") as f:
            f.write(map16x_data)

        # Generate and write thumbnail
        t_type, t_mag, t_shore = thumbnail_source_terrain
        thumbnail = self.create_thumbnail(t_type, t_mag, t_shore, 0.5)
        thumbnail.save(os.path.join(variant_output_dir, "thumbnail.webp"), "WEBP")

        # Manifest
        input_manifest = {}
        if os.path.exists(self.json_path):
            with open(self.json_path, "r") as f:
                input_manifest = json.load(f)
        
        # Calculate dimensions directly from the terrain arrays
        map_type, map4x_type, map16x_type = dimension_terrains
        map_h, map_w = map_type.shape
        map4x_h, map4x_w = map4x_type.shape
        map16x_h, map16x_w = map16x_type.shape
        
        manifest = {
            "name": variant_name,
            "map": {
                "width": map_w,
                "height": map_h,
                "num_land_tiles": map_land
            },
            "map4x": {
                "width": map4x_w,
                "height": map4x_h,
                "num_land_tiles": map4x_land
            },
            "map16x": {
                "width": map16x_w,
                "height": map16x_h,
                "num_land_tiles": map16x_land
            },
            "nations": []
        }

        # Nations
        if "points" in input_manifest:
            for p in input_manifest["points"]:
                nation = {
                    "name": p.get("name", "Unknown"),
                    "flag": p.get("flag", "unknown"),
                    "coordinates": [
                        int(p.get("pixel_x", 0) * scale_factor), 
                        int(p.get("pixel_y", 0) * scale_factor)
                    ]
                }
                manifest["nations"].append(nation)
        elif "nations" in input_manifest:
             for n in input_manifest["nations"]:
                nation = n.copy()
                nation["coordinates"] = [
                    int(n["coordinates"][0] * scale_factor),
                    int(n["coordinates"][1] * scale_factor)
                ]
                manifest["nations"].append(nation)

        with open(os.path.join(variant_output_dir, "manifest.json"), "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)
            
        print(f"Saved outputs to {variant_output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate map files from PNG and JSON")
    parser.add_argument("map_name", help="Name of the map (folder name)")
    parser.add_argument("--input", help="Input directory (containing map folder)", default=".")
    parser.add_argument("--output", help="Output directory", default="generated")
    parser.add_argument("--test", action="store_true", help="Test mode (skip small island removal)")
    
    args = parser.parse_args()
    
    # Construct paths
    # If input is "StylisedMaps", and map_name is "Korea", look in "StylisedMaps/Korea"
    input_path = os.path.join(args.input, args.map_name)
    output_path = os.path.join(args.output, args.map_name)
    
    generator = MapGenerator(args.map_name, input_path, output_path, args.test)
    generator.process()
