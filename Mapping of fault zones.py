import os
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
import geopandas as gpd
from shapely.geometry import LineString, Polygon, Point
from shapely.ops import unary_union, nearest_points

def read_multiline_txt(path):
    lines = []
    coords = []
    with open(path, 'r', encoding='utf-8') as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith('>'):
                if coords:
                    lines.append(LineString(coords))
                    coords = []
                continue
            parts = [p for p in s.replace(',', ' ').split() if p]
            if len(parts) < 2:
                continue
            lon, lat = map(float, parts[:2])
            coords.append((lon, lat))
    if coords:
        lines.append(LineString(coords))
    return gpd.GeoDataFrame(geometry=lines, crs="EPSG:4326")

txt_file = "duance.txt"  
fault_lines = read_multiline_txt(txt_file)

lon_min, lon_max = 100.0, 118.0
lat_min, lat_max = 31.0, 43.0
outer_poly = Polygon([
    (lon_min, lat_min),
    (lon_max, lat_min),
    (lon_max, lat_max),
    (lon_min, lat_max),
    (lon_min, lat_min),
])
outer_boundary = outer_poly.boundary

all_lines = list(fault_lines.geometry) + [outer_boundary]
union_lines = unary_union(all_lines)

fig, ax = plt.subplots(figsize=(8, 8))
fault_lines.plot(ax=ax, color="red", linewidth=0.8, label="Faults")
gpd.GeoSeries([outer_boundary], crs="EPSG:4326").plot(
    ax=ax, linestyle="--", color="black", linewidth=1.0, label="Boundary"
)

ax.set_xlim(lon_min, lon_max)
ax.set_ylim(lat_min, lat_max)
ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
ax.set_title("Hold down the left button to draw a block shape; release to generate\nPress 'u' to undo the last block, and press 'q' to exit")
ax.set_aspect("equal", adjustable="box")

block_id = 1

drawn_blocks = []



def onselect(verts):
    """
    verts: [(x1,y1), (x2,y2), ...] 
    """
    global block_id, drawn_blocks
    snapped = []
    for x, y in verts:
        p = Point(x, y)
        _, p_on_line = nearest_points(p, union_lines)
        snapped.append((p_on_line.x, p_on_line.y))
    if len(snapped) < 3:
        print("There are too few points to form a polygon. This drawing will be ignored.")
        return
    poly = Polygon(snapped)
    xs, ys = poly.exterior.xy
    line_artist, = ax.plot(xs, ys, color="orange", linewidth=1.5)
    cx, cy = poly.representative_point().x, poly.representative_point().y
    text_artist = ax.text(cx, cy, str(block_id), color="purple",
                          ha="center", va="center", fontsize=12)
    fig.canvas.draw_idle()
    out_name = f"block_{block_id}.txt"
    with open(out_name, "w", encoding="utf-8") as f:
        for x, y in poly.exterior.coords:
            f.write(f"{x:.6f} {y:.6f}\n")
    print(f"The block with ID {block_id} has been saved as {out_name}")
    drawn_blocks.append((line_artist, text_artist, out_name))

    block_id += 1

lasso = LassoSelector(ax, onselect=onselect)


def on_key(event):
    global block_id, drawn_blocks
    if event.key == "u":
        if not drawn_blocks:
            print("Currently, there are no blocks that can be undone.。")
            return

        line_artist, text_artist, filename = drawn_blocks.pop()
        line_artist.remove()
        text_artist.remove()
        if os.path.exists(filename):
            os.remove(filename)
            print(f"Has been revoked and deleted {filename}")
        else:
            print(f"The block has been deleted, but the file was not found. {filename}")
        block_id -= 1

        fig.canvas.draw_idle()

    elif event.key == "q":
        plt.close(fig)

fig.canvas.mpl_connect("key_press_event", on_key)

plt.legend(loc="lower left")
plt.show()
