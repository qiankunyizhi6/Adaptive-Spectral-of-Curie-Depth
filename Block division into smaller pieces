import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import unary_union, polygonize
from shapely.errors import GEOSException
from shapely.validation import make_valid   
import shapely
import matplotlib.pyplot as plt

def read_block_txt(path):
    xs, ys = [], []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            parts = [p for p in s.replace(",", " ").split() if p]
            if len(parts) < 2:
                continue
            lon, lat = map(float, parts[:2])
            xs.append(lon)
            ys.append(lat)
    if len(xs) < 3:
        raise ValueError(f"{path} point less than3，can't get a Polygon")
    return np.column_stack([xs, ys])

def clean_geom(g):
    if g.is_empty:
        return g
    try:
        g2 = make_valid(g)
    except Exception:
        g2 = g.buffer(0)
    return g2

def force_polygon(g):
    g = clean_geom(g)
    if g.is_empty:
        return None

    if g.geom_type == "Polygon":
        return g

    if g.geom_type == "MultiPolygon":
        polys = list(g.geoms)
        return max(polys, key=lambda p: p.area)

    try:
        polys = list(polygonize(g))
    except Exception:
        polys = []

    if polys:
        return max(polys, key=lambda p: p.area)
    buf = g.buffer(0)
    if buf.geom_type == "Polygon":
        return buf
    if buf.geom_type == "MultiPolygon":
        polys = list(buf.geoms)
        return max(polys, key=lambda p: p.area)
    return None

block_files = [f"block_{i}.txt" for i in range(1, 10)]

geoms = []
ids = []
for i, fname in enumerate(block_files, start=1):
    coords = read_block_txt(fname)
    poly = Polygon(coords)
    poly = clean_geom(poly)
    geoms.append(poly)
    ids.append(i)

blocks = gpd.GeoDataFrame({"Block_ID": ids, "geometry": geoms}, crs="EPSG:4326")
print("Original block：", len(blocks))

blocks = blocks.sort_values("Block_ID").reset_index(drop=True)

clean_geoms = []
current_union = None

for idx, row in blocks.iterrows():
    geom = clean_geom(row.geometry)

    if current_union is not None and not current_union.is_empty:
        cu = clean_geom(current_union)
        try:
            geom = geom.difference(cu)
        except GEOSException as e:
            print(f"Block {row['Block_ID']} difference ", e)
            geom = clean_geom(geom).difference(clean_geom(cu))

    clean_geoms.append(geom)

    if current_union is None or current_union.is_empty:
        current_union = geom
    else:
        try:
            current_union = unary_union([clean_geom(current_union), clean_geom(geom)])
        except GEOSException as e:
            print("union，buffer(0) ", e)
            current_union = unary_union([clean_geom(current_union).buffer(0),
                                         clean_geom(geom).buffer(0)])

blocks["geometry"] = clean_geoms
blocks = blocks[~blocks.geometry.is_empty].copy().reset_index(drop=True)
print("The number of effective blocks after removing overlaps：", len(blocks))

all_bounds = np.array([g.bounds for g in blocks.geometry])
minx, miny = all_bounds[:, 0].min(), all_bounds[:, 1].min()
maxx, maxy = all_bounds[:, 2].max(), all_bounds[:, 3].max()

outer_poly = Polygon([
    (minx, miny),
    (maxx, miny),
    (maxx, maxy),
    (minx, maxy),
    (minx, miny),
])
outer_poly = clean_geom(outer_poly)

union_clean = unary_union([clean_geom(g) for g in blocks.geometry])
gap_geom = outer_poly.difference(clean_geom(union_clean))

gap_polys = []
if gap_geom.is_empty:
    print("no white")
else:
    if gap_geom.geom_type == "Polygon":
        gap_polys = [gap_geom]
    elif gap_geom.geom_type == "MultiPolygon":
        gap_polys = list(gap_geom.geoms)
    else:
        if isinstance(gap_geom, shapely.geometry.collection.GeometryCollection):
            for gg in gap_geom.geoms:
                if gg.geom_type in ["Polygon", "MultiPolygon"]:
                    if gg.geom_type == "Polygon":
                        gap_polys.append(gg)
                    else:
                        gap_polys.extend(list(gg.geoms))

print("Number of blank area blocks：", len(gap_polys))

for gap in gap_polys:
    gap = clean_geom(gap)
    cen = gap.representative_point()
    distances = [cen.distance(clean_geom(g)) for g in blocks.geometry]
    nearest_idx = int(np.argmin(distances))
    blocks.at[nearest_idx, "geometry"] = clean_geom(
        blocks.geometry.iloc[nearest_idx]
    ).union(gap)

blocks = blocks.dissolve(by="Block_ID", as_index=False)
print("The number of blocks after filling the blanks：", len(blocks))

blocks["geometry"] = blocks["geometry"].apply(force_polygon)
blocks = blocks[blocks.geometry.notnull()].copy().reset_index(drop=True)
print("The number of blocks after being forcibly converted to Polygon：", len(blocks))

fig, ax = plt.subplots(figsize=(8, 8))

colors = [
    "#1f77b4", "#ff7f0e", "#2ca02c",
    "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22"
]

for _, row in blocks.iterrows():
    geom = row.geometry
    bid = int(row["Block_ID"])
    xs, ys = geom.exterior.xy
    ax.fill(xs, ys,
            facecolor=colors[(bid - 1) % len(colors)],
            edgecolor="black", alpha=0.5, linewidth=1.0)

    rp = geom.representative_point()
    ax.text(rp.x, rp.y, str(bid),
            ha="center", va="center",
            fontsize=12, color="black")

ax.set_xlim(minx, maxx)
ax.set_ylim(miny, maxy)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
ax.set_title("Nine Structural Blocks (fixed)")

plt.tight_layout()
plt.show()

for idx, row in blocks.iterrows():
    bid = int(row["Block_ID"])
    geom = row.geometry

    poly = force_polygon(geom)
    if poly is None:
        print(f"The block {bid} cannot be converted into a polygon; it will be skipped.")
        continue

    out_name = f"block_{bid}_fixed.txt"
    with open(out_name, "w", encoding="utf-8") as f:
        for x, y in poly.exterior.coords:
            f.write(f"{x:.6f} {y:.6f}\n")

    print(f"The revised block has been written. {bid} 为 {out_name}")

