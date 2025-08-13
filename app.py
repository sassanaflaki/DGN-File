import os
import io
import json
import tempfile
from collections import Counter, defaultdict

import streamlit as st

# Optional but recommended if installed in your env
try:
    from osgeo import ogr, osr
except Exception as e:
    st.error("GDAL/OGR Python bindings are required. Install via conda: 'conda install -c conda-forge gdal'.\nIf you're on pip-only, ensure wheels match your Python/OS.")
    st.stop()

st.set_page_config(page_title="DGN Viewer & Converter", layout="wide")
st.title("🗺️ DGN Viewer & Converter (Streamlit)")
st.caption("Upload a MicroStation .dgn file, inspect layers/attributes, filter by Level, preview on a map, and export to GeoPackage/GeoJSON.")

# ------------------------------
# Helpers
# ------------------------------

def save_upload_to_tmp(uploaded_file) -> str:
    suffix = os.path.splitext(uploaded_file.name)[1]
    fd, tmp_path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return tmp_path


def get_driver_and_layers(ds):
    return ds.GetDriver().GetName(), [ds.GetLayerByIndex(i) for i in range(ds.GetLayerCount())]


def layer_fields(lyr):
    defn = lyr.GetLayerDefn()
    return [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]


def summarize_field(lyr, field_name: str, max_unique=1000):
    idx = lyr.FindFieldIndex(field_name, True)
    if idx == -1:
        return None
    counts = Counter()
    unique_seen = 0
    for f in lyr:
        val = f.GetField(idx)
        counts[val] += 1
        unique_seen = len(counts)
        if unique_seen > max_unique:
            break
    lyr.ResetReading()
    return counts


def compute_bbox(lyr):
    # overall extent (minx, maxx, miny, maxy)
    try:
        ext = lyr.GetExtent()
        # OGR returns (minx, maxx, miny, maxy)
        return ext
    except Exception:
        return None


def make_transform(source_epsg: int | None, target_epsg: int | None):
    if not source_epsg or not target_epsg or source_epsg == target_epsg:
        return None
    src = osr.SpatialReference(); src.ImportFromEPSG(int(source_epsg))
    dst = osr.SpatialReference(); dst.ImportFromEPSG(int(target_epsg))
    return osr.CoordinateTransformation(src, dst)


def ogr_feature_to_geojson_feature(f, transform=None):
    g = f.GetGeometryRef()
    if g is None:
        return None
    g = g.Clone()
    if transform is not None:
        try:
            g.Transform(transform)
        except Exception:
            return None
    geom_json = json.loads(g.ExportToJson())
    # properties
    props = {}
    defn = f.GetDefnRef()
    for i in range(defn.GetFieldCount()):
        name = defn.GetFieldDefn(i).GetName()
        props[name] = f.GetField(i)
    return {"type": "Feature", "geometry": geom_json, "properties": props}


def layer_to_feature_collection(lyr, level_filter=None, transform=None, max_features=50000):
    features = []
    level_idx = lyr.FindFieldIndex("Level", True)
    n = 0
    for f in lyr:
        if level_filter is not None and level_idx != -1:
            lvl = f.GetField(level_idx)
            if lvl not in level_filter:
                continue
        feat = ogr_feature_to_geojson_feature(f, transform)
        if feat:
            features.append(feat)
            n += 1
            if n >= max_features:
                break
    lyr.ResetReading()
    return {"type": "FeatureCollection", "features": features}


def write_geopackage_from_ds(ds, out_path, source_epsg=None, target_epsg=None, level_filters=None):
    drv = ogr.GetDriverByName("GPKG")
    try:
        drv.DeleteDataSource(out_path)
    except Exception:
        pass
    out_ds = drv.CreateDataSource(out_path)

    transform = make_transform(source_epsg, target_epsg)
    out_srs = None
    if transform is not None and target_epsg:
        out_srs = osr.SpatialReference(); out_srs.ImportFromEPSG(int(target_epsg))
    elif source_epsg:
        out_srs = osr.SpatialReference(); out_srs.ImportFromEPSG(int(source_epsg))

    for i in range(ds.GetLayerCount()):
        in_lyr = ds.GetLayerByIndex(i)
        in_defn = in_lyr.GetLayerDefn()
        out_lyr_name = in_lyr.GetName() or f"layer_{i}"
        out_lyr = out_ds.CreateLayer(out_lyr_name, srs=out_srs, geom_type=in_lyr.GetGeomType())
        for j in range(in_defn.GetFieldCount()):
            out_lyr.CreateField(in_defn.GetFieldDefn(j))
        out_defn = out_lyr.GetLayerDefn()

        level_idx = in_lyr.FindFieldIndex("Level", True)
        for f in in_lyr:
            if level_filters and level_idx != -1:
                lvl = f.GetField(level_idx)
                if out_lyr_name in level_filters and len(level_filters[out_lyr_name]) > 0 and lvl not in level_filters[out_lyr_name]:
                    continue
            g = f.GetGeometryRef()
            if g is None:
                continue
            g = g.Clone()
            if transform is not None:
                try:
                    g.Transform(transform)
                except Exception:
                    continue
            of = ogr.Feature(out_defn)
            of.SetGeometry(g)
            for j in range(out_defn.GetFieldCount()):
                nm = out_defn.GetFieldDefn(j).GetName()
                of.SetField(nm, f.GetField(nm))
            out_lyr.CreateFeature(of)
        in_lyr.ResetReading()

    out_ds = None
    return out_path


# ------------------------------
# Sidebar Controls
# ------------------------------
with st.sidebar:
    st.header("Settings")
    source_epsg = st.text_input("Source EPSG (if DGN has no CRS)", value="", help="e.g., 2248 (NAD83 / Maryland), 2283 (NAD83 / Virginia North), 26918 (UTM18N), etc.")
    target_epsg = st.text_input("Target EPSG for map/export", value="4326", help="Default 4326 (WGS84 lat/lon) for web maps")
    max_preview = st.number_input("Max preview features per layer", min_value=1000, max_value=200000, value=50000, step=1000)
    st.markdown("---")
    st.caption("Tip: If you see odd geometry, check the correct source EPSG.")

uploaded = st.file_uploader("Upload a .dgn file", type=["dgn"]) 
if not uploaded:
    st.info("Drop a DGN above to begin.")
    st.stop()

# Save and open
path = save_upload_to_tmp(uploaded)
try:
    ds = ogr.Open(path)
except Exception as e:
    st.exception(e)
    st.stop()

if ds is None:
    st.error("Could not open DGN. Ensure the file is valid and GDAL supports your DGN version (V7 best supported).")
    st.stop()

# Dataset summary
with st.expander("Dataset Summary", expanded=True):
    dr_name, layers = get_driver_and_layers(ds)
    st.write({"Driver": dr_name, "Layer count": len(layers)})

# Build a per-layer Level index for UI filters
layer_level_values = {}
for i, lyr in enumerate(layers):
    lvl_counts = summarize_field(lyr, "Level")
    if lvl_counts:
        layer_level_values[lyr.GetName() or f"layer_{i}"] = sorted([k for k in lvl_counts.keys() if k is not None])

# Layer browser
tabs = st.tabs([ (l.GetName() or f"Layer {i}") for i, l in enumerate(layers) ])

# Create a transform for preview/export
src_epsg_val = int(source_epsg) if source_epsg.strip().isdigit() else None
tgt_epsg_val = int(target_epsg) if target_epsg.strip().isdigit() else None
transform = make_transform(src_epsg_val, tgt_epsg_val)

export_level_filters = {}

for i, (tab, lyr) in enumerate(zip(tabs, layers)):
    with tab:
        name = lyr.GetName() or f"layer_{i}"
        st.subheader(f"{name}")

        # Fields
        fields = layer_fields(lyr)
        cols = st.columns([1,1,1,1])
        with cols[0]:
            st.write("**Geom Type:**", ogr.GeometryTypeToName(lyr.GetGeomType()))
        with cols[1]:
            st.write("**Feature Count (raw):**", lyr.GetFeatureCount())
        with cols[2]:
            st.write("**Fields:**", ", ".join(fields[:12]) + (" …" if len(fields) > 12 else ""))
        with cols[3]:
            ext = compute_bbox(lyr)
            if ext:
                minx, maxx, miny, maxy = ext
                st.write("**Extent:**", f"[{minx:.3f}, {miny:.3f}] – [{maxx:.3f}, {maxy:.3f}]")
            else:
                st.write("**Extent:**", "(unknown)")

        # Level filter UI
        if name in layer_level_values:
            sel_levels = st.multiselect(
                "Filter by Level (optional)", options=layer_level_values[name], default=[], key=f"levels_{i}"
            )
            export_level_filters[name] = set(sel_levels)
        else:
            st.caption("No 'Level' field detected for this layer.")
            export_level_filters[name] = set()

        # Preview on map using GeoJSON + pydeck
        st.markdown("**Layer Preview** (first N features, transformed to target EPSG if provided)")
        fc = layer_to_feature_collection(lyr, level_filter=export_level_filters[name] if len(export_level_filters[name])>0 else None, transform=transform, max_features=max_preview)
        if len(fc["features"]) == 0:
            st.warning("No features matched the current filter (or no geometries).")
        else:
            # Compute a view state from the preview bbox (after transform via features)
            # Fallback to DC area if missing.
            xs = []
            ys = []
            for f in fc["features"]:
                g = f["geometry"]
                # Collect a few representative coords
                def add_coords(obj):
                    if obj is None:
                        return
                    if isinstance(obj[0], (int, float)) and len(obj) >= 2:
                        xs.append(obj[0]); ys.append(obj[1])
                    else:
                        for sub in obj:
                            add_coords(sub)
                coords = g.get("coordinates")
                add_coords(coords)
            if xs and ys:
                center = (sum(xs)/len(xs), sum(ys)/len(ys))
                zoom_guess = 12
            else:
                center = (-77.0365, 38.8977)  # DC fallback
                zoom_guess = 9

            try:
                import pydeck as pdk
                layer = pdk.Layer(
                    "GeoJsonLayer",
                    data=fc,
                    pickable=True,
                    auto_highlight=True,
                    extruded=False,
                )
                view_state = pdk.ViewState(longitude=center[0], latitude=center[1], zoom=zoom_guess)
                r = pdk.Deck(layers=[layer], initial_view_state=view_state, tooltip={"text": "{Level} {ColorIndex}"})
                st.pydeck_chart(r)
            except Exception as e:
                st.info("Install pydeck to see a map preview: 'pip install pydeck'. Showing raw GeoJSON below.")
                st.json(fc if len(fc["features"]) < 2000 else {"type": "FeatureCollection", "features": fc["features"][:2000]})

        # Quick attribute summaries
        with st.expander("Quick Summaries (by Level, ColorIndex, Weight)"):
            for fld in ("Level", "ColorIndex", "Weight"):
                counts = summarize_field(lyr, fld)
                if counts:
                    st.write(f"**{fld}**: {dict(counts.most_common(20))}{' …' if len(counts)>20 else ''}")
                else:
                    st.caption(f"No {fld} field.")

# ------------------------------
# Exports
# ------------------------------
st.markdown("---")
colA, colB = st.columns([1,1])

with colA:
    st.subheader("Export GeoPackage (.gpkg)")
    if st.button("Build GeoPackage", type="primary"):
        try:
            out_path = tempfile.mktemp(suffix=".gpkg")
            # Map layer names to filters only for layers with filters selected
            cleaned_filters = {k: v for k, v in export_level_filters.items() if len(v) > 0}
            write_geopackage_from_ds(ds, out_path, source_epsg=src_epsg_val, target_epsg=tgt_epsg_val, level_filters=cleaned_filters)
            with open(out_path, "rb") as f:
                st.download_button("Download .gpkg", data=f, file_name=(uploaded.name.rsplit('.',1)[0] + ".gpkg"), mime="application/octet-stream")
        except Exception as e:
            st.exception(e)

with colB:
    st.subheader("Export current tab preview as GeoJSON")
    tab_idx = st.number_input("Tab index to export (0-based)", min_value=0, max_value=len(layers)-1, value=0, step=1)
    lyr = layers[int(tab_idx)]
    name = lyr.GetName() or f"layer_{tab_idx}"
    level_filter = export_level_filters.get(name, set())
    fc = layer_to_feature_collection(lyr, level_filter=level_filter if len(level_filter)>0 else None, transform=transform, max_features=1_000_000)
    geojson_bytes = io.BytesIO(json.dumps(fc).encode("utf-8"))
    st.download_button("Download GeoJSON", data=geojson_bytes, file_name=f"{name}.geojson", mime="application/geo+json")

st.markdown(
    """
**Notes**
- DGN V7 is best supported by GDAL. For complex V8 files (true arcs, cells, rich symbology), consider converting to DXF/GeoPackage with MicroStation/ODA, then use this app.
- Most DGN files lack embedded CRS. Set **Source EPSG** in the sidebar if coordinates look off. Default export/preview CRS is WGS84 (EPSG:4326).
- The map preview limits features for performance. Increase the limit in *Settings* if needed.
"""
)
