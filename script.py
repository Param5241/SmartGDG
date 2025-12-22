import ezdxf
from ezdxf.enums import TextEntityAlignment

# --- Constants ---
cell_width = 20000
cell_height = 3000
start_x = 140000
start_y = 285000
text_height = 550


# --- Data ---
bom_data = {
    "PIPE": ["Ø180", "Ø125", "Ø63", "Ø32", "Ø20"],
    "WARNING MAT": ["300x0.25MM", "300x1.0MM"],
    "HDPE PIPE 50ø": [],
    "HDPE PIPE 125ø": [],
    "COUPLER": ["Ø180", "Ø125", "Ø63", "Ø32", "Ø20"],
    "SADDLE": ["Ø125xØ63", "Ø125xØ32", "Ø63xØ32", "Ø32xØ20"],
    "REDUCER": ["Ø180xØ125", "Ø125xØ63", "Ø63xØ32", "Ø32xØ20"],
    "EQUAL TEE": ["Ø180", "Ø125", "Ø63", "Ø32", "Ø20"],
    "ELBOW": ["Ø125x90°", "Ø125x45°", "Ø63x90°"],
    "END CAP": ["Ø180","Ø125", "Ø63", "Ø32", "Ø20"],
    "TRANSITION FITTING WITH COUPLER": ["1/2\"xØ20", "1\"xØ32"],
    "BALL VALVE": ["Ø180", "Ø125", "Ø63", "Ø32"],
}
individual_merged_items = [
    "PLATE MARKER", "STONE MARKER TYPE–A", "STONE MARKER TYPE–B",
    "S.R.–200 SCMH", "I.V.–1/2\"", "I.V.–1\""
]

# --- Helper to center text ---
def add_centered_text(msp, text, x1, x2, y1, y2, height):
    cx, cy = (x1 + x2) / 2, (y1 + y2) / 2
    txt = msp.add_text(text, dxfattribs={"height": height})
    txt.set_placement((cx, cy), align=TextEntityAlignment.MIDDLE_CENTER)

# --- Load Template ---
template_path = "Template/FTP.dxf"
doc = ezdxf.readfile(template_path)
msp = doc.modelspace()

# --- BOM Header ---
add_centered_text(msp, "BILL OF MATERIALS", start_x, start_x + 3 * cell_width, start_y + 10, start_y + 7000, text_height+50)
add_centered_text(msp, "NEW DASHMESH COLONY, SECTOR–8, RAJPURA", start_x, start_x + 3 * cell_width, start_y, start_y + 3000, text_height+50)

# --- BOM Table Header ---
y = start_y - 5
msp.add_lwpolyline([(start_x, y), (start_x + cell_width, y), (start_x + cell_width, y - cell_height), (start_x, y - cell_height)], close=True)
add_centered_text(msp, "GROUP", start_x, start_x + cell_width, y, y - cell_height, text_height)

for i, header in enumerate(["SIZE", "LENGTH (Mtrs.)"]):
    x1 = start_x + cell_width + i * cell_width
    x2 = x1 + cell_width
    msp.add_lwpolyline([(x1, y), (x2, y), (x2, y - cell_height), (x1, y - cell_height)], close=True)
    add_centered_text(msp, header, x1, x2, y, y - cell_height, text_height)

y -= cell_height

# --- BOM Table Body ---
for group, items in bom_data.items():
    if group in ["HDPE PIPE 50ø", "HDPE PIPE 125ø"]:
        x1, x2 = start_x, start_x + 2 * cell_width
        x3 = start_x + 2 * cell_width
        msp.add_lwpolyline([(x1, y), (x2, y), (x2, y - cell_height), (x1, y - cell_height)], close=True)
        msp.add_lwpolyline([(x3, y), (x3 + cell_width, y), (x3 + cell_width, y - cell_height), (x3, y - cell_height)], close=True)
        add_centered_text(msp, group, x1, x2, y, y - cell_height, text_height)
        y -= cell_height
        continue

    num_rows = max(len(items), 1)
    y_top = y
    y_bot = y - num_rows * cell_height
    msp.add_lwpolyline([(start_x, y_top), (start_x + cell_width, y_top), (start_x + cell_width, y_bot), (start_x, y_bot)], close=True)
    add_centered_text(msp, group, start_x, start_x + cell_width, y_top, y_bot, text_height)

    for i in range(num_rows):
        # Size cell
        xs1, xs2 = start_x + cell_width, start_x + 2 * cell_width
        msp.add_lwpolyline([(xs1, y), (xs2, y), (xs2, y - cell_height), (xs1, y - cell_height)], close=True)
        # Length cell
        xl1, xl2 = xs2, xs2 + cell_width
        msp.add_lwpolyline([(xl1, y), (xl2, y), (xl2, y - cell_height), (xl1, y - cell_height)], close=True)
        if i < len(items):
            add_centered_text(msp, items[i], xs1, xs2, y, y - cell_height, text_height)
        y -= cell_height

# --- Final individual merged rows ---
for label in individual_merged_items:
    x1, x2 = start_x, start_x + 2 * cell_width
    x3, x4 = x2, x2 + cell_width
    y1, y2 = y, y - cell_height
    msp.add_lwpolyline([(x1, y1), (x2, y1), (x2, y2), (x1, y2)], close=True)
    msp.add_lwpolyline([(x3, y1), (x4, y1), (x4, y2), (x3, y2)], close=True)
    add_centered_text(msp, label, x1, x2, y1, y2, text_height)
    y -= cell_height

# ======================
# ✅ AS BUILT TABLE 
# ======================
asbuilt_start_y = y - 82000

asbuilt_cell_w = 22000
asbuilt_cell_h = 4000
asbuilt_col_widths = [8000, 12000, 8000, 12000, 8000, 12000, 8000, 12000, 8000, 20000]
col_x = [start_x - 43000]
for w in asbuilt_col_widths:
    col_x.append(col_x[-1] + w)

# Last 2 rows: 7 equal columns
total_width = sum(asbuilt_col_widths)
col_x_last = [start_x - 43000]
for _ in range(7):
    col_x_last.append(col_x_last[-1] + total_width / 7)

row_y = [asbuilt_start_y - i * asbuilt_cell_h for i in range(15)]

def draw_asbuilt_cell(c1, c2, r1, r2, cols, text=""):
    x1, x2 = cols[c1], cols[c2]
    y1, y2 = row_y[r1], row_y[r2]
    msp.add_lwpolyline([(x1, y1), (x2, y1), (x2, y2), (x1, y2)], close=True)
    add_centered_text(msp, text, x1, x2, y1, y2, text_height)

# Rows 1–3
for r in range(3):
    for c in range(10):
        draw_asbuilt_cell(c, c+1, r, r+1, col_x, f"R{r+1}C{c+1}")

# Row 4
for c in range(8):
    draw_asbuilt_cell(c, c+1, 3, 4, col_x, f"R4C{c+1}")
draw_asbuilt_cell(8, 9, 3, 5, col_x, "STATUS")
draw_asbuilt_cell(9, 10, 3, 4, col_x, "SUBJECT")

# Row 5
draw_asbuilt_cell(0, 2, 4, 5, col_x, "DRAWN")
draw_asbuilt_cell(2, 4, 4, 5, col_x, "CHECKED")
draw_asbuilt_cell(4, 6, 4, 5, col_x, "APPROVED")
draw_asbuilt_cell(6, 8, 4, 5, col_x, "SIGN")
draw_asbuilt_cell(9, 10, 4, 5, col_x, "NOTE")

# Rows 6–10
for r in range(5, 10):
    draw_asbuilt_cell(0, 10, r, r+1, col_x, f"Torrent Gas Pvt. Ltd. R{r+1}")

# Row 11
draw_asbuilt_cell(0, 3, 10, 11, col_x, "PROJECT")
draw_asbuilt_cell(3, 6, 10, 11, col_x, "DRAWING TITLE")
draw_asbuilt_cell(6, 10, 10, 11, col_x, "AS BUILT")

# Rows 12–13
for r in [11, 12]:
    for c in range(7):
        draw_asbuilt_cell(c, c+1, r, r+1, col_x_last, f"R{r+1}C{c+1}")

# === ADD VERTICAL "AS BUILT" COLUMN ===

# Define column position for the vertical text box (right of the last column)
x1 = col_x[-1]  # small gap after last column
x2 = x1 + 3000         # width of the vertical label

# Vertical span from first row to last row (rows 0 to 14)
y1 = row_y[0]
y2 = row_y[13]

# Draw vertical rectangle
msp.add_lwpolyline([(x1, y1), (x2, y1), (x2, y2), (x1, y2)], close=True)

# Add vertical text centered in the vertical cell
cx = (x1 + x2) / 2
cy = (y1 + y2) / 2

msp.add_text("AS BUILT", dxfattribs={
    "height": 1000,
    "rotation": 90
}).set_placement((cx, cy), align=TextEntityAlignment.MIDDLE_CENTER)


# ======================
# ✅ GENERAL NOTES + DIAMOND
# ======================
diamond_center_x = start_x - 0.5 * cell_width
diamond_center_y = row_y[-1] + 100000
diamond_half = 20000
text_height_main = 1000
text_height_notes = 800

top = (diamond_center_x, diamond_center_y + diamond_half)
right = (diamond_center_x + diamond_half, diamond_center_y)
bottom = (diamond_center_x, diamond_center_y - diamond_half)
left = (diamond_center_x - diamond_half, diamond_center_y)

msp.add_lwpolyline([top, right, bottom, left, top], close=True)
msp.add_line((left[0], left[1]), (right[0], right[1]))

def diamond_center_text(text, y_offset, height=text_height_main):
    msp.add_text(text, dxfattribs={"height": height}).set_placement(
        (diamond_center_x, diamond_center_y + y_offset),
        align=TextEntityAlignment.MIDDLE_CENTER
    )

# Diamond upper text
diamond_center_text("Ø32 MM", 6000)
diamond_center_text("&", 4000)
diamond_center_text("Ø20 MM PIPE", 2000)

# Diamond lower text
diamond_center_text("LAKAD MANDI RAJPURA", -2000)
diamond_center_text("SECTOR–8,", -4000)
diamond_center_text("RAJPURA,", -6000)
diamond_center_text("PUNJAB", -8000)

# General Notes below diamond
title_y = diamond_center_y - diamond_half - 10000
msp.add_text("GENERAL NOTES", dxfattribs={"height": 1200}).set_placement(
    (diamond_center_x, title_y),
    align=TextEntityAlignment.MIDDLE_CENTER
)

note_y = title_y - 3500
notes = [
    "1. ALL DIMENSIONS ARE IN METRES. UNLESS OTHERWISE SPECIFIED.",
    "2. PIPE SIZES ARE IN MILLIMETRES"
]

for note in notes:
    msp.add_text(note, dxfattribs={"height": text_height_notes}).set_placement(
        (diamond_center_x - 15000, note_y),
        align=TextEntityAlignment.LEFT
    )
    note_y -= 2000




# === Top Center Text ===
top_rect_x1 = 12000
top_rect_x2 = 30000
top_rect_y1 = 280000
top_rect_y2 = 305000

# Calculate center
top_center_x = (top_rect_x1 + top_rect_x2) / 2
top_center_y = (top_rect_y1 + top_rect_y2) / 2

# Add top center text
msp.add_text("Title", dxfattribs={
    "height": 4000
}).set_placement((top_center_x, top_center_y), align=TextEntityAlignment.MIDDLE_CENTER)

# === Add Compass Rose (Top-Right Corner) ===

# Center coordinates of the compass (adjust to fit your layout)
compass_center_x = start_x + 3 * cell_width   # top-right of BOM table
compass_center_y = start_y + 7000              # higher above the BOM title

compass_radius = 2000
line_length = 2500
text_height = 1000

# Draw outer circle
msp.add_circle(center=(compass_center_x, compass_center_y), radius=compass_radius)

# Draw 4 directional lines
msp.add_line(
    (compass_center_x, compass_center_y + compass_radius),
    (compass_center_x, compass_center_y + compass_radius + line_length)
)
msp.add_line(
    (compass_center_x, compass_center_y - compass_radius),
    (compass_center_x, compass_center_y - compass_radius - line_length)
)
msp.add_line(
    (compass_center_x - compass_radius, compass_center_y),
    (compass_center_x - compass_radius - line_length, compass_center_y)
)
msp.add_line(
    (compass_center_x + compass_radius, compass_center_y),
    (compass_center_x + compass_radius + line_length, compass_center_y)
)

# Add the "N" text
msp.add_text("N", dxfattribs={"height": text_height}).set_placement(
    (compass_center_x, compass_center_y),
    align=TextEntityAlignment.MIDDLE_CENTER
)

import time

# --- Save ---
output_path = "Template/FTP_FINAL.dxf"
doc.saveas(output_path)
print(f"✅ Final drawing saved: {output_path}")

# Wait 5 seconds before running the second script
time.sleep(5)









# ============ SCRIPT 2 ============









import os
import json
import base64
import shutil
import ezdxf
from ezdxf.math import Vec2
import requests
import matplotlib.pyplot as plt
from PIL import Image
import math
import numpy as np

# ================= CONFIG =================
TEXT_HEIGHT = 450
SYMBOL_SIZE = 4500
SYMBOL_FOLDER = "icons"
TEMP_IMAGE_PATH = "temp"
TEMPLATE_PATH = "Template/FTP_FINAL.dxf"

# Map background box in DXF units (final placement rectangle)
FIXED_MIN_X, FIXED_MAX_X = -170000, 20000
FIXED_MIN_Y, FIXED_MAX_Y = 19000, 280000

# Dimension arrow configuration
ARROW_SIZE = 1500  # Size of the arrow heads
DIMENSION_OFFSET = 3000  # Offset distance for dimension lines from the main line
DIMENSION_TEXT_HEIGHT = 500  # Text height for dimension labels

# Layer definitions (AutoCAD index colors)
TYPE_LAYER_COLOR_MAP = {
    "THIRTYTWO_MM": (1, "THIRTYTWO_MM"),
    "DOUBLE_ARROW": (3, "DOUBLE_ARROW"),
    "ONE_SIXTY_N_ONE_EIGHTY_MM": (151, "ONE_SIXTY_N_ONE_EIGHTY_MM"),
    "ONE_TWENTYFIVE_MM": (1, "ONE_TWENTYFIVE_MM"),
    "NINETY_MM": (5, "NINETY_MM"),
    "SIXTYTHREE_MM": (3, "SIXTYTHREE_MM"),
    "SERVICE_LINE": (30, "SERVICE_LINE"),
    "ALL_OTHER_VALUES": (134, "ALL_OTHER_VALUES"),
}

# ---------- Web Mercator (meters) ----------
def latlon_to_webmercator(lat, lon):
    R = 6378137.0
    x = R * math.radians(lon)
    y = R * math.log(math.tan(math.pi/4 + math.radians(lat)/2))
    return x, y

def webmercator_to_latlon(x, y):
    R = 6378137.0
    lon = math.degrees(x / R)
    lat = math.degrees(2 * math.atan(math.exp(y / R)) - math.pi/2)
    return lat, lon

def meters_to_degrees(lat, meters):
    """Convert meters to degrees at a given latitude"""
    lat_degree = meters / 111111
    lon_degree = meters / (111111 * math.cos(math.radians(lat)))
    return lat_degree, lon_degree

def get_osm_data_directly(min_lat, min_lon, max_lat, max_lon):
    """Get OSM data directly using Overpass API"""
    overpass_url = "http://overpass-api.de/api/interpreter"
    
    # Query for roads and buildings
    overpass_query = f"""
    [out:json][timeout:25];
    (
      way["highway"]({min_lat},{min_lon},{max_lat},{max_lon});
      way["building"]({min_lat},{min_lon},{max_lat},{max_lon});
    );
    out body;
    >;
    out skel qt;
    """
    
    try:
        print(f"🌐 Fetching OSM data from Overpass API for bbox: {min_lat:.6f}, {min_lon:.6f}, {max_lat:.6f}, {max_lon:.6f}")
        response = requests.post(overpass_url, data=overpass_query)
        response.raise_for_status()
        print("✅ OSM data fetched successfully")
        return response.json()
    except Exception as e:
        print(f"❌ Error fetching data from Overpass API: {e}")
        return None

def extract_map_features_from_osm_data(osm_data):
    """Extract road lines and building polygons from OSM data"""
    features = {
        'roads': [],
        'buildings': []
    }
    
    if not osm_data or 'elements' not in osm_data:
        return features
    
    # Group elements by type and collect node coordinates
    nodes = {}
    ways = []
    
    for element in osm_data['elements']:
        if element['type'] == 'node':
            nodes[element['id']] = (element['lat'], element['lon'])
        elif element['type'] == 'way':
            ways.append(element)
    
    # Process ways
    for way in ways:
        try:
            # Get coordinates for all nodes in this way
            coords = []
            for node_id in way.get('nodes', []):
                if node_id in nodes:
                    lat, lon = nodes[node_id]
                    coords.append((lon, lat))  # Store as (lon, lat) for consistency
            
            if len(coords) < 2:
                continue
            
            # Check if it's a road or building
            tags = way.get('tags', {})
            if 'highway' in tags:
                features['roads'].append({
                    'type': tags.get('highway', 'road'),
                    'coordinates': coords
                })
            elif 'building' in tags:
                features['buildings'].append({
                    'type': 'building',
                    'coordinates': coords
                })
                
        except Exception as e:
            print(f"⚠️ Error processing way {way.get('id', 'unknown')}: {e}")
            continue
    
    return features

def encode_pngs_to_base64(folder_path):
    encoded = {}
    for filename in os.listdir(folder_path):
        if filename.lower().endswith(".png"):
            symbol_name = filename[:-4].upper()
            with open(os.path.join(folder_path, filename), "rb") as img_file:
                encoded[symbol_name] = base64.b64encode(img_file.read()).decode("utf-8")
    return encoded

def write_base64_to_file(symbol_name, base64_data):
    os.makedirs(TEMP_IMAGE_PATH, exist_ok=True)
    output_path = os.path.join(TEMP_IMAGE_PATH, f"{symbol_name}.png")
    with open(output_path, "wb") as f:
        f.write(base64.b64decode(base64_data))
    return output_path

def calculate_optimal_scale_and_offset(line_data, text_data, symbol_data, map_features):
    """Calculate scale and offset to center the pipeline data within the fixed box"""
    # Collect all pipeline data coordinates (lines, text, symbols)
    all_pipeline_points = []
    
    # Add line endpoints
    for line in line_data:
        all_pipeline_points.append((line["starting_point"]["x"], line["starting_point"]["y"]))
        all_pipeline_points.append((line["ending_point"]["x"], line["ending_point"]["y"]))
    
    # Add text and symbol points
    for text in text_data:
        all_pipeline_points.append((text["coordinates"]["x"], text["coordinates"]["y"]))
    
    for symbol in symbol_data:
        all_pipeline_points.append((symbol["coordinates"]["x"], symbol["coordinates"]["y"]))
    
    if not all_pipeline_points:
        raise ValueError("No pipeline data coordinates found")
    
    # Calculate pipeline data bounds
    pipeline_lats = [p[0] for p in all_pipeline_points]
    pipeline_lons = [p[1] for p in all_pipeline_points]
    
    min_lat, max_lat = min(pipeline_lats), max(pipeline_lats)
    min_lon, max_lon = min(pipeline_lons), max(pipeline_lons)
    
    print(f"📏 Pipeline data bounds: {min_lat:.6f}, {min_lon:.6f} to {max_lat:.6f}, {max_lon:.6f}")
    
    # Calculate center of pipeline data
    center_lat = (min_lat + max_lat) / 2
    center_lon = (min_lon + max_lon) / 2
    
    # Add buffer around pipeline data (100 meters)
    lat_buffer, lon_buffer = meters_to_degrees(center_lat, 100)
    min_lat -= lat_buffer
    max_lat += lat_buffer
    min_lon -= lon_buffer
    max_lon += lon_buffer
    
    # Calculate bounds in Web Mercator for pipeline data
    tl_lat, tl_lon = max_lat, min_lon
    br_lat, br_lon = min_lat, max_lon
    
    tl_xm, tl_ym = latlon_to_webmercator(tl_lat, tl_lon)
    br_xm, br_ym = latlon_to_webmercator(br_lat, br_lon)
    
    pipeline_min_xm = min(tl_xm, br_xm)
    pipeline_min_ym = min(tl_ym, br_ym)
    pipeline_max_xm = max(tl_xm, br_xm)
    pipeline_max_ym = max(tl_ym, br_ym)
    
    # Calculate pipeline data dimensions
    pipeline_width_m = pipeline_max_xm - pipeline_min_xm
    pipeline_height_m = pipeline_max_ym - pipeline_min_ym
    
    # Calculate fixed box dimensions
    box_width = FIXED_MAX_X - FIXED_MIN_X
    box_height = FIXED_MAX_Y - FIXED_MIN_Y
    
    # Calculate scale factors based on pipeline data
    scale_x = box_width / pipeline_width_m if pipeline_width_m > 0 else 1.0
    scale_y = box_height / pipeline_height_m if pipeline_height_m > 0 else 1.0
    
    # Use the smaller scale to ensure pipeline data fits perfectly
    scale = min(scale_x, scale_y) * 2.5  # 90% scale to add margin
    
    # Calculate center of pipeline data in scaled coordinates
    pipeline_center_x = (pipeline_min_xm + pipeline_max_xm) / 2 * scale
    pipeline_center_y = (pipeline_min_ym + pipeline_max_ym) / 2 * scale
    
    # Calculate center of fixed box
    box_center_x = (FIXED_MIN_X + FIXED_MAX_X) / 2
    box_center_y = (FIXED_MIN_Y + FIXED_MAX_Y) / 2
    
    # Calculate offset to center the pipeline data
    offset_x = box_center_x - pipeline_center_x
    offset_y = box_center_y - pipeline_center_y
    
    print(f"🎯 Pipeline center: ({pipeline_center_x:.2f}, {pipeline_center_y:.2f})")
    print(f"🎯 Box center: ({box_center_x:.2f}, {box_center_y:.2f})")
    
    return scale, offset_x, offset_y

def is_point_in_bounds(x, y):
    """Check if a point is within the fixed bounding box"""
    return (FIXED_MIN_X <= x <= FIXED_MAX_X and 
            FIXED_MIN_Y <= y <= FIXED_MAX_Y)

def clip_line_to_bounds(p1, p2):
    """Clip a line segment to the bounding box using Cohen-Sutherland algorithm"""
    def compute_outcode(x, y):
        code = 0
        if x < FIXED_MIN_X: code |= 1  # left
        elif x > FIXED_MAX_X: code |= 2  # right
        if y < FIXED_MIN_Y: code |= 4  # bottom
        elif y > FIXED_MAX_Y: code |= 8  # top
        return code

    x1, y1 = p1
    x2, y2 = p2
    
    outcode1 = compute_outcode(x1, y1)
    outcode2 = compute_outcode(x2, y2)
    
    accept = False
    
    while True:
        if not (outcode1 | outcode2):
            # Both points inside the clip window
            accept = True
            break
        elif outcode1 & outcode2:
            # Both points share an outside zone -> completely outside
            break
        else:
            # Partially inside the clip window
            x, y = 0.0, 0.0
            outcode_out = outcode1 if outcode1 else outcode2
            
            if outcode_out & 8:  # top
                x = x1 + (x2 - x1) * (FIXED_MAX_Y - y1) / (y2 - y1)
                y = FIXED_MAX_Y
            elif outcode_out & 4:  # bottom
                x = x1 + (x2 - x1) * (FIXED_MIN_Y - y1) / (y2 - y1)
                y = FIXED_MIN_Y
            elif outcode_out & 2:  # right
                y = y1 + (y2 - y1) * (FIXED_MAX_X - x1) / (x2 - x1)
                x = FIXED_MAX_X
            elif outcode_out & 1:  # left
                y = y1 + (y2 - y1) * (FIXED_MIN_X - x1) / (x2 - x1)
                x = FIXED_MIN_X
            
            if outcode_out == outcode1:
                x1, y1 = x, y
                outcode1 = compute_outcode(x1, y1)
            else:
                x2, y2 = x, y
                outcode2 = compute_outcode(x2, y2)
    
    if accept:
        return [(x1, y1), (x2, y2)]
    else:
        return None

def clip_polygon_to_bounds(points):
    """Clip a polygon to the bounding box using Sutherland-Hodgman algorithm"""
    def inside(p, edge):
        x, y = p
        if edge == 'left': return x >= FIXED_MIN_X
        elif edge == 'right': return x <= FIXED_MAX_X
        elif edge == 'bottom': return y >= FIXED_MIN_Y
        elif edge == 'top': return y <= FIXED_MAX_Y
    
    def intersection(p1, p2, edge):
        x1, y1 = p1
        x2, y2 = p2
        
        if edge == 'left':
            x = FIXED_MIN_X
            if x2 != x1:
                y = y1 + (y2 - y1) * (FIXED_MIN_X - x1) / (x2 - x1)
            else:
                y = y1
        elif edge == 'right':
            x = FIXED_MAX_X
            if x2 != x1:
                y = y1 + (y2 - y1) * (FIXED_MAX_X - x1) / (x2 - x1)
            else:
                y = y1
        elif edge == 'bottom':
            y = FIXED_MIN_Y
            if y2 != y1:
                x = x1 + (x2 - x1) * (FIXED_MIN_Y - y1) / (y2 - y1)
            else:
                x = x1
        elif edge == 'top':
            y = FIXED_MAX_Y
            if y2 != y1:
                x = x1 + (x2 - x1) * (FIXED_MAX_Y - y1) / (y2 - y1)
            else:
                x = x1
                
        return (x, y)
    
    # Clip against each edge
    edges = ['left', 'right', 'bottom', 'top']
    output_list = points
    
    for edge in edges:
        if len(output_list) == 0:
            break
            
        input_list = output_list
        output_list = []
        
        s = input_list[-1]
        for p in input_list:
            if inside(p, edge):
                if not inside(s, edge):
                    output_list.append(intersection(s, p, edge))
                output_list.append(p)
            elif inside(s, edge):
                output_list.append(intersection(s, p, edge))
            s = p
    
    return output_list

def create_dimension_arrow(msp, start_point, end_point, length_text, layer_name="DIMENSIONS"):
    """Create a dimension line with arrows and text in the center"""
    # Calculate line vector and properties
    line_vector = Vec2(end_point) - Vec2(start_point)
    line_length = line_vector.magnitude
    line_angle = math.atan2(line_vector.y, line_vector.x)
    
    # Calculate perpendicular offset direction
    perpendicular_angle = line_angle + math.pi/2
    offset_vector = Vec2(math.cos(perpendicular_angle), math.sin(perpendicular_angle)) * DIMENSION_OFFSET
    
    # Calculate dimension line points (parallel to main line but offset)
    dim_start = Vec2(start_point) + offset_vector
    dim_end = Vec2(end_point) + offset_vector
    
    # Draw dimension line
    msp.add_line(dim_start, dim_end, dxfattribs={
        "layer": layer_name,
        "color": 2,  # Yellow color for dimensions
        "lineweight": 20
    })
    
    # Calculate arrow positions
    arrow1_start = dim_start
    arrow1_dir = Vec2(math.cos(line_angle), math.sin(line_angle))
    arrow2_start = dim_end
    arrow2_dir = -arrow1_dir
    
    # Create arrow heads
    def create_arrow_head(position, direction, size=ARROW_SIZE):
        angle = math.atan2(direction.y, direction.x)
        arrow_angle1 = angle + math.pi * 0.75  # 135 degrees
        arrow_angle2 = angle - math.pi * 0.75  # -135 degrees
        
        point1 = position + Vec2(math.cos(arrow_angle1), math.sin(arrow_angle1)) * size
        point2 = position + Vec2(math.cos(arrow_angle2), math.sin(arrow_angle2)) * size
        
        # Draw arrow head lines
        msp.add_line(position, point1, dxfattribs={
            "layer": layer_name,
            "color": 2,
            "lineweight": 20
        })
        msp.add_line(position, point2, dxfattribs={
            "layer": layer_name,
            "color": 2,
            "lineweight": 20
        })
    
    # Draw arrow heads
    create_arrow_head(arrow1_start, arrow1_dir)
    create_arrow_head(arrow2_start, arrow2_dir)
    
    # Draw extension lines from main line to dimension line
    msp.add_line(start_point, dim_start, dxfattribs={
        "layer": layer_name,
        "color": 2,
        "lineweight": 15,
        "linetype": "DASHED"
    })
    msp.add_line(end_point, dim_end, dxfattribs={
        "layer": layer_name,
        "color": 2,
        "lineweight": 15,
        "linetype": "DASHED"
    })
    
    # Calculate text position (center of dimension line)
    text_position = (dim_start + dim_end) / 2
    
    # Adjust text position to be slightly offset from the dimension line for better visibility
    text_offset = offset_vector.normalize() * (DIMENSION_OFFSET * 0.3)
    text_position += text_offset
    
    # Add dimension text
    msp.add_mtext(length_text, dxfattribs={
        "char_height": DIMENSION_TEXT_HEIGHT,
        "layer": layer_name,
        "color": 2,
        "rotation": math.degrees(line_angle)
    }).set_location(text_position)

def generate_dxf_from_json(json_file, output_file):
    if os.path.exists(TEMP_IMAGE_PATH):
        shutil.rmtree(TEMP_IMAGE_PATH)

    with open(json_file) as f:
        data = json.load(f)
    line_data = data.get("lines", [])
    text_data = data.get("texts", [])
    symbol_data = data.get("symbols", [])

    png_db = encode_pngs_to_base64(SYMBOL_FOLDER)

    # Read template DXF
    doc = ezdxf.readfile(TEMPLATE_PATH)
    msp = doc.modelspace()

    # Ensure layers exist
    for color, layer in TYPE_LAYER_COLOR_MAP.values():
        if layer not in doc.layers:
            doc.layers.add(layer, color=color)
    for l in ["TEXT", "SYMBOLS", "MAP_ROADS", "MAP_BUILDINGS", "BOUNDING_BOX", "LINE_LENGTHS", "DIMENSIONS"]:
        if l not in doc.layers:
            doc.layers.add(l, color=2)

    # -------- Build bbox from JSON for OSM data --------
    lat_vals, lon_vals = [], []
    for p in line_data:
        lat_vals.extend([p["starting_point"]["x"], p["ending_point"]["x"]])
        lon_vals.extend([p["starting_point"]["y"], p["ending_point"]["y"]])
    for p in text_data + symbol_data:
        lat_vals.append(p["coordinates"]["x"])
        lon_vals.append(p["coordinates"]["y"])

    if not lat_vals or not lon_vals:
        raise ValueError("No coordinates found in JSON.")

    min_lat, max_lat = min(lat_vals), max(lat_vals)
    min_lon, max_lon = min(lon_vals), max(lon_vals)

    # Store original bounds for reference
    json_min_lat, json_max_lat = min_lat, max_lat
    json_min_lon, json_max_lon = min_lon, max_lon

    # Calculate center point of the data
    center_lat = (min_lat + max_lat) / 2
    center_lon = (min_lon + max_lon) / 2

    # Add 200 meters buffer in all directions for OSM data
    lat_buffer, lon_buffer = meters_to_degrees(center_lat, 200)
    
    osm_min_lat = min_lat - lat_buffer
    osm_max_lat = max_lat + lat_buffer
    osm_min_lon = min_lon - lon_buffer
    osm_max_lon = max_lon + lon_buffer

    print(f"📏 Pipeline data bounds: {json_min_lat:.6f}, {json_min_lon:.6f} to {json_max_lat:.6f}, {json_max_lon:.6f}")
    print(f"🗺️  OSM query bounds (200m buffer): {osm_min_lat:.6f}, {osm_min_lon:.6f} to {osm_max_lat:.6f}, {osm_max_lon:.6f}")

    # -------- Get OSM data using Overpass API --------
    print("🌐 Fetching OSM data from Overpass API...")
    osm_data = get_osm_data_directly(osm_min_lat, osm_min_lon, osm_max_lat, osm_max_lon)
    
    if osm_data is None:
        print("⚠️ Could not fetch OSM data, proceeding without map background")
        map_features = {'roads': [], 'buildings': []}
    else:
        print("📊 Processing OSM data...")
        map_features = extract_map_features_from_osm_data(osm_data)
        print(f"✅ Extracted {len(map_features['roads'])} road segments and {len(map_features['buildings'])} buildings")

    # Calculate scale and offset based on PIPELINE DATA (not OSM data)
    scale, offset_x, offset_y = calculate_optimal_scale_and_offset(
        line_data, text_data, symbol_data, map_features
    )

    print(f"📏 Scale factor: {scale:.6f}")
    print(f"📍 Offset: X={offset_x:.2f}, Y={offset_y:.2f}")

    # ========== DRAW OSM MAP FEATURES (CLIPPED) ==========
    print("🗺️ Drawing OSM features (clipped to bounds)...")
    
    # Draw roads (clipped)
    road_count = 0
    for road in map_features['roads']:
        # Convert all points to DXF coordinates
        dxf_points = []
        for lon, lat in road['coordinates']:
            xm, ym = latlon_to_webmercator(lat, lon)
            dxf_x = xm * scale + offset_x
            dxf_y = ym * scale + offset_y
            dxf_points.append((dxf_x, dxf_y))
        
        # Clip the road line by line segments
        clipped_segments = []
        for i in range(len(dxf_points) - 1):
            p1 = dxf_points[i]
            p2 = dxf_points[i + 1]
            clipped_line = clip_line_to_bounds(p1, p2)
            if clipped_line:
                clipped_segments.append(clipped_line)
        
        # Draw the clipped segments
        for segment in clipped_segments:
            if len(segment) == 2:
                msp.add_line(segment[0], segment[1], dxfattribs={
                    "layer": "MAP_ROADS", 
                    "color": 7,
                    "lineweight": 25
                })
                road_count += 1

    # Draw buildings (clipped)
    building_count = 0
    for building in map_features['buildings']:
        dxf_points = []
        for lon, lat in building['coordinates']:
            xm, ym = latlon_to_webmercator(lat, lon)
            dxf_x = xm * scale + offset_x
            dxf_y = ym * scale + offset_y
            dxf_points.append((dxf_x, dxf_y))
        
        # Clip polygon to bounds
        if len(dxf_points) > 2:  # Need at least 3 points for a polygon
            clipped_polygon = clip_polygon_to_bounds(dxf_points)
            
            if len(clipped_polygon) > 2:
                # Close the polygon
                clipped_polygon.append(clipped_polygon[0])
                msp.add_lwpolyline(clipped_polygon, dxfattribs={
                    "layer": "MAP_BUILDINGS", 
                    "color": 8,
                    "lineweight": 35
                })
                building_count += 1

    print(f"✅ OSM features drawn: {road_count} road segments, {building_count} buildings")

    # ========== DRAW LINES WITH DIMENSION ARROWS ==========
    line_count = 0
    dimension_count = 0
    
    for line in line_data:
        sp_geo = latlon_to_webmercator(line["starting_point"]["x"], line["starting_point"]["y"])
        ep_geo = latlon_to_webmercator(line["ending_point"]["x"], line["ending_point"]["y"])
        sp = (sp_geo[0] * scale + offset_x, sp_geo[1] * scale + offset_y)
        ep = (ep_geo[0] * scale + offset_x, ep_geo[1] * scale + offset_y)
        
        # Calculate actual length in meters
        dx = ep_geo[0] - sp_geo[0]
        dy = ep_geo[1] - sp_geo[1]
        length_m = math.hypot(dx, dy)
        
        # Clip line to bounds
        clipped_line = clip_line_to_bounds(sp, ep)
        
        if clipped_line:
            color, layer_name = TYPE_LAYER_COLOR_MAP.get(
                line.get("type", "ALL_OTHER_VALUES"),
                TYPE_LAYER_COLOR_MAP["ALL_OTHER_VALUES"],
            )
            msp.add_line(clipped_line[0], clipped_line[1], dxfattribs={"layer": layer_name})
            line_count += 1

            # Add dimension arrows if line is fully visible
            if is_point_in_bounds(*sp) and is_point_in_bounds(*ep):
                try:
                    create_dimension_arrow(msp, clipped_line[0], clipped_line[1], f"{length_m:.1f} m")
                    dimension_count += 1
                except Exception as e:
                    print(f"⚠️ Could not create dimension for line: {e}")
                    # Fallback: add simple text
                    mid_x, mid_y = (sp[0] + ep[0]) / 2.0, (sp[1] + ep[1]) / 2.0
                    if is_point_in_bounds(mid_x, mid_y):
                        msp.add_mtext(f"{length_m:.1f} m", dxfattribs={
                            "char_height": DIMENSION_TEXT_HEIGHT,
                            "layer": "LINE_LENGTHS",
                            "color": 2
                        }).set_location(Vec2(mid_x, mid_y))

    print(f"✅ Pipeline lines drawn: {line_count} lines")
    print(f"✅ Dimension arrows drawn: {dimension_count} dimensions")

    # ========== DRAW TEXT (ONLY IF IN BOUNDS) ==========
    text_count = 0
    for text in text_data:
        tx, ty = latlon_to_webmercator(text["coordinates"]["x"], text["coordinates"]["y"])
        pos_x = tx * scale + offset_x
        pos_y = ty * scale + offset_y
        
        # Only draw text if it's within bounds
        if is_point_in_bounds(pos_x, pos_y):
            msp.add_mtext(text["text"], dxfattribs={
                "char_height": TEXT_HEIGHT,
                "layer": "TEXT",
                "color": 7
            }).set_location(Vec2(pos_x, pos_y))
            text_count += 1

    # ========== DRAW SYMBOLS (ONLY IF IN BOUNDS) ==========
    missing_symbols = set()
    symbol_count = 0
    for symbol in symbol_data:
        sx, sy = latlon_to_webmercator(symbol["coordinates"]["x"], symbol["coordinates"]["y"])
        pos_x = sx * scale + offset_x
        pos_y = sy * scale + offset_y
        
        # Only draw symbol if it's within bounds
        if is_point_in_bounds(pos_x, pos_y):
            symbol_name = symbol["symbol"].upper()
            if symbol_name not in png_db:
                missing_symbols.add(symbol_name)
                continue
            png_path = write_base64_to_file(symbol_name, png_db[symbol_name])
            image_def = doc.add_image_def(png_path, size_in_pixel=(96, 96))
            msp.add_image(
                image_def,
                insert=(pos_x, pos_y),
                size_in_units=(SYMBOL_SIZE, SYMBOL_SIZE),
                rotation=symbol.get("rotation", 0),
                dxfattribs={"layer": "SYMBOLS"},
            )
            symbol_count += 1

    print(f"✅ Text elements drawn: {text_count}")
    print(f"✅ Symbols drawn: {symbol_count}")

    # Draw the bounding box for reference
    # msp.add_lwpolyline([
    #     (FIXED_MIN_X, FIXED_MIN_Y),
    #     (FIXED_MAX_X, FIXED_MIN_Y),
    #     (FIXED_MAX_X, FIXED_MAX_Y),
    #     (FIXED_MIN_X, FIXED_MAX_Y),
    #     (FIXED_MIN_X, FIXED_MIN_Y)
    # ], dxfattribs={
    #     "layer": "BOUNDING_BOX",
    #     "color": 1,
    #     "lineweight": 50
    # })

    # Save DXF
    doc.saveas(output_file)
    print(f"\n✅ DXF saved at: {output_file}")
    print(f"📐 Fixed box: X={FIXED_MIN_X} to {FIXED_MAX_X}, Y={FIXED_MIN_Y} to {FIXED_MAX_Y}")
    print(f"🎯 Pipeline data perfectly centered in bounding box")
    print(f"🔒 All elements clipped to bounds - no overflow")
    print(f"📏 Dimension arrows added to {dimension_count} lines")
    if missing_symbols:
        print("⚠️ Missing PNGs for symbols:", ", ".join(sorted(missing_symbols)))

# ========== RUN EXAMPLE ==========
if __name__ == "__main__":
    generate_dxf_from_json(
        json_file="Data/SMRTCGD6856.json",
        # json_file="/Users/apple/Downloads/Azoca/Torent Gas/Primary/Final/Data/SMRTCGD9450.json",
        output_file="Output/SmartCGD.dxf"
    )