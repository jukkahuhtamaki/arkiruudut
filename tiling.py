import sys
import pdb
from lxml import etree
from shapely.geometry import Polygon, MultiLineString, LineString
from shapely.geometry import shape
from shapely.geometry import box
from kml2geojson import convert
import matplotlib.pyplot as plt
import mercantile
import simplekml
from shapely.ops import unary_union

# -----------------------------------

def kml_to_shapes(kml_file):
    # Read KML as XML
    tree = etree.parse(kml_file)
    root = tree.getroot()

    # Find all placemark elements
    placemarks = root.findall(".//{http://www.opengis.net/kml/2.2}Placemark")

    polygon_storage = []

    for pm in placemarks:
        # find the polygons
        polygons = pm.findall(".//{http://www.opengis.net/kml/2.2}Polygon")
        for polygon in polygons:
            # find the outer boundary => 1 pcs

            outer_boundary = polygon.find(".//{http://www.opengis.net/kml/2.2}outerBoundaryIs")
            ob_coordinates = outer_boundary.find(".//{http://www.opengis.net/kml/2.2}coordinates")
            ob_coordinates = [x.strip() for x in ob_coordinates.text.split(' ') if x]
            ob_coordinates = [[float(y) for y in x.split(',')] for x in ob_coordinates]

            hole_coordinates = []

            inner_boundaries = polygon.findall(".//{http://www.opengis.net/kml/2.2}innerBoundaryIs")
            for hole in inner_boundaries:
                in_coordinates = hole.find(".//{http://www.opengis.net/kml/2.2}coordinates")
                in_coordinates = [x.strip() for x in in_coordinates.text.split(' ') if x]
                in_coordinates = [[float(y) for y in x.split(',')] for x in in_coordinates]
                hole_coordinates.append(in_coordinates)

            # Create the polygon
            polygon_storage.append(Polygon(ob_coordinates, holes = hole_coordinates))
    return polygon_storage

# --------------------------

def shapes_to_osm(shapes, multi_grid, output_file="output.osm"):
    osm = etree.Element("osm", version="0.6", generator="Python")
    
    node_id = 1
    way_id = 1
    nodes = {}  # Track nodes to avoid duplicates
    
    # Helper function to add nodes safely
    def add_node(lat, lon):
        nonlocal node_id
        node_key = (lat, lon)
        if node_key not in nodes:
            node = etree.SubElement(osm, "node", id=str(node_id), lat=str(lat), lon=str(lon), version="1")
            nodes[node_key] = node_id
            node_id += 1
        return nodes[node_key]

    # Process Polygon and LineString shapes
    for shape in shapes:
        if isinstance(shape, Polygon):
            way = etree.SubElement(osm, "way", id=str(way_id), version="1")
            for coord in shape.exterior.coords:
                node_ref = add_node(coord[1], coord[0])  # Store unique nodes
                etree.SubElement(way, "nd", ref=str(node_ref))
            way_id += 1

            for interior in shape.interiors:
                way = etree.SubElement(osm, "way", id=str(way_id), version="1")  # Separate way for holes
                for coord in interior.coords:
                    node_ref = add_node(coord[1], coord[0])
                    etree.SubElement(way, "nd", ref=str(node_ref))
                way_id += 1
        
        elif isinstance(shape, LineString):
            way = etree.SubElement(osm, "way", id=str(way_id), version="1")
            for coord in shape.coords:
                node_ref = add_node(coord[1], coord[0])
                etree.SubElement(way, "nd", ref=str(node_ref))
            way_id += 1
    
    # ✅ Process MultiLineString (Grid Lines)
    if isinstance(multi_grid, MultiLineString):
        for line in multi_grid.geoms:
            way = etree.SubElement(osm, "way", id=str(way_id), version="1")
            for coord in line.coords:
                node_ref = add_node(coord[1], coord[0])
                etree.SubElement(way, "nd", ref=str(node_ref))
            way_id += 1
    
    # Save to file
    tree = etree.ElementTree(osm)
    tree.write(output_file, encoding="utf-8", xml_declaration=True)
    print(f"OSM file saved as {output_file}")



###############

def OLD_shapes_to_osm(shapes, multi_grid, output_file="output.osm"):
    osm = etree.Element("osm", version="0.6", generator="Python")

    node_id = 1
    way_id = 1
    nodes = {}

    # Process Polygon and LineString shapes
    
    for shape in shapes:
        if isinstance(shape, Polygon):
            way = etree.SubElement(osm, "way", id=str(way_id), version="1")  # 🔹 Added version="1"
            for coord in shape.exterior.coords:
                lat, lon = coord[1], coord[0]
                node = etree.SubElement(osm, "node", id=str(node_id), lat=str(lat), lon=str(lon), version="1")  # 🔹 Added version="1"
                nodes[node_id] = node
                etree.SubElement(way, "nd", ref=str(node_id))
                node_id += 1
            way_id += 1

            for interior in shape.interiors:
                for coord in interior.coords:
                    lat, lon = coord[1], coord[0]
                    node = etree.SubElement(osm, "node", id=str(node_id), lat=str(lat), lon=str(lon), version="1")  # 🔹 Added version="1"
                    nodes[node_id] = node
                    etree.SubElement(way, "nd", ref=str(node_id))
                    node_id += 1
            way_id += 1

        elif isinstance(shape, LineString):
            way = etree.SubElement(osm, "way", id=str(way_id), version="1")  # 🔹 Added version="1"
            for coord in shape.coords:
                lat, lon = coord[1], coord[0]
                node = etree.SubElement(osm, "node", id=str(node_id), lat=str(lat), lon=str(lon), version="1")  # 🔹 Added version="1"
                nodes[node_id] = node
                etree.SubElement(way, "nd", ref=str(node_id))
                node_id += 1
            way_id += 1
    
    
    # ✅ Process MultiLineString (Grid Lines)
    if isinstance(multi_grid, MultiLineString):
        for line in multi_grid.geoms:
            way = etree.SubElement(osm, "way", id=str(way_id), version="1")  # 🔹 Added version="1"
            for coord in line.coords:
                lat, lon = coord[1], coord[0]
                node = etree.SubElement(osm, "node", id=str(node_id), lat=str(lat), lon=str(lon), version="1")  # 🔹 Added version="1"
                nodes[node_id] = node
                etree.SubElement(way, "nd", ref=str(node_id))
                node_id += 1
            way_id += 1
    

    # Save to file
    tree = etree.ElementTree(osm)
    tree.write(output_file, encoding="utf-8", xml_declaration=True)
    print(f"OSM file saved as {output_file}")

# -----------------------------------    

def create_tile_grid(minx, miny, maxx, maxy, zoom=14):
    grid_lines = []

    # Get tile coordinates for bounding box
    tile_min = mercantile.tile(minx, miny, zoom)
    tile_max = mercantile.tile(maxx, maxy, zoom)

    # Generate vertical lines (Longitude-based)
    if tile_min.x > tile_max.x +1:
        step = -1
    else:
        step = 1
    for x in range(tile_min.x, tile_max.x + 1, step):
        lon_left, lat_top = mercantile.bounds(x, tile_min.y, zoom)[:2]
        lon_right, lat_bottom = mercantile.bounds(x, tile_max.y, zoom)[:2]
        grid_lines.append(LineString([(lon_left, lat_bottom), (lon_left, lat_top)]))

    # Generate horizontal lines (Latitude-based)
    if tile_min.y > tile_max.y +1:
        step = -1
    else:
        step = 1
    for y in range(tile_min.y, tile_max.y + 1, step):
        lon_left, lat_top = mercantile.bounds(tile_min.x, y, zoom)[:2]
        lon_right, lat_bottom = mercantile.bounds(tile_max.x, y, zoom)[:2]
        grid_lines.append(LineString([(lon_left, lat_top), (lon_right, lat_top)]))

    return MultiLineString(grid_lines)


# --------------------------

def geometry_to_kml(shapes, multi_grid, output_file="output.kml"):
    kml = simplekml.Kml()
    for shape in shapes:
        if isinstance(shape, Polygon):  
            # Add polygon with outer boundary
            poly = kml.newpolygon(outerboundaryis=list(shape.exterior.coords))

            # Add ALL inner boundaries (holes)
            poly.innerboundaryis = [list(interior.coords) for interior in shape.interiors]
            
        elif isinstance(shape, MultiLineString):  
            for line in geom.geoms:
                kml.newlinestring(coords=list(line.coords))

        elif isinstance(shape, LineString):  
            kml.newlinestring(coords=list(geom.coords))

    # add the multigrid
    if isinstance(multi_grid, MultiLineString):
        for line in multi_grid.geoms:  # Iterate through each LineString in MultiLineString
            kml.newlinestring(coords=list(line.coords))

    kml.save(output_file)
    print(f"KML saved as {output_file}")

#------------------------------

#kml_file='polygon04.kml'
#bounding_box = [2.0, 39, 3.6, 40]

kml_file='squadrats-2025-05-10.kml'
bounding_box = [23.55,61.4, 23.92, 61.55]
#bounding_box = [0,0, 90, 90]

print("parsitaan kml alueisiin")
shapes = kml_to_shapes(kml_file)

print("leikataan alueet bounding boxilla")
bounding_shape = box(*bounding_box)

#new_shapes = []
#for shape in shapes:
#    if bounding_shape.contains(shape):
#        new_shapes.append(shape)
#    else:
#        new_shapes.append(bounding_shape.intersection(shape))

new_shapes = shapes


# Generate grid lines
grid_lines = []
step = .1  # Grid spacing
minx=bounding_box[0]
maxx=bounding_box[2]
miny=bounding_box[1]
maxy=bounding_box[3]



print("gridviivaston luonti")
multi_grid = create_tile_grid(minx, miny, maxx, maxy, zoom=17)

# vähennetään käydyt alueet gridiviivastosta
print(f"leikataan gridiviivat: {len(shapes)} leikkausta")
#for shape in shapes:
#    multi_grid = multi_grid.difference(shape) 

# 🔥 Optimized: Union all shapes first
print(f"Kasataan {len(shapes)} aluetta yhdeksi leikkausgeometriaksi")
merged_shapes = unary_union(shapes)  # Combine all shapes into one

# 🔥 Perform difference in one step
print("Leikataan gridiviivat yhdellä operaatiolla")
multi_grid = multi_grid.difference(merged_shapes)



print("kasataan KML file")
geometry_to_kml(shapes, multi_grid, output_file="output.kml")

print("kasataan OSM file")
shapes_to_osm(shapes, multi_grid)



"""
gpsbabel -i kml -f output.kml -o osm,tag=highway:primary -F missing_tiles.osm

mkgmap --mapname-format="p20" --polygon-size-limits=0:1 -c config.txt typ.txt missing_tiles.osm 

java -jar /usr/share/mkgmap/mkgmap.jar --keep-going --verbose --debug --log-level=DEBUG -c config.txt typ.txt missing_tiles.osm 


echo "# arkiruudut" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:Myrtillus/arkiruudut.git
git push -u origin main




"""