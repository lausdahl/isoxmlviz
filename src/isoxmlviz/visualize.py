import argparse
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
import sys
import pymap3d
import shapely.geometry
import shapely.geometry as SHP

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from shapely.geometry import LineString, JOIN_STYLE, MultiLineString, MultiPoint
import math
from isoxmlviz.LineStringUtil import extract_lines_within, get_coordinates
from isoxmlviz.webmap import WebMap


def PolygonPatch(polygon, **kwargs):
    """Constructs a matplotlib patch from a geometric object

    The `polygon` may be a Shapely or GeoJSON-like object with or without holes.
    The `kwargs` are those supported by the matplotlib.patches.Polygon class
    constructor. Returns an instance of matplotlib.patches.PathPatch.

    Example (using Shapely Point and a matplotlib axes):

      >>> b = Point(0, 0).buffer(1.0)
      >>> patch = PolygonPatch(b, fc='blue', ec='blue', alpha=0.5)
      >>> axis.add_patch(patch)

    """
    if type(polygon) is shapely.geometry.MultiPolygon:
        return PatchCollection([PolygonPatch(p, **kwargs) for p in polygon.geoms])
    # from descartes but no longer maintained so inspired by https://github.com/geopandas/geopandas/issues/1039
    from matplotlib.path import Path
    from matplotlib.patches import Polygon
    import numpy as np
    return Polygon(Path.make_compound_path(Path(np.asarray(polygon.exterior.coords)[:, :2]),
                                           *[Path(np.asarray(ring.coords)[:, :2]) for ring in
                                             polygon.interiors]).vertices, **kwargs)


ell_wgs84 = pymap3d.Ellipsoid.from_name('wgs84')

default_propagation_num = 100


def generate_web_safe_colors():
    color_levels = [0, 51, 102, 153, 204, 255]
    web_safe_colors = [
        "#{:02X}{:02X}{:02X}{:02X}".format(r, g, b, a)
        for a in [255]  # color_levels
        for r in color_levels
        for g in color_levels
        for b in color_levels
    ]
    return web_safe_colors


colour_map = {
    0: 'black',
    1: 'white',
    2: 'green',
    3: 'teal',
    4: 'maroon',
    5: 'purple',
    6: 'olive',
    7: 'silver',
    8: 'grey',
    9: 'blue',
    10: 'lime',
    11: 'cyan',
    12: 'red',
    13: 'magenta',
    14: 'yellow',
    15: 'navy',
}

for idx, c in enumerate(generate_web_safe_colors()):
    colour_map[idx + 16] = c


def pnt_to_pair(element: ET):
    # C=lat
    return float(element.attrib.get("C")), float(element.attrib.get("D"))


def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a


class WebGroups:

    def __init__(self) -> None:
        super().__init__()

        self.line_type_groups = {}
        self.polygon_type_groups = {}
        self.field_name_group = None


def main():
    options = argparse.ArgumentParser(prog="isoxmlviz")
    options.add_argument("-file", dest="file", type=str, required=True, help='Path to a isoxml task file XML or ZIP')
    options.add_argument("-p", "--pdf", dest="pdf", action="store_true", required=False, help='Write figure to pdf')
    options.add_argument("-hide", "--hide", dest="hide", action="store_true", required=False, help='Hide plot')
    options.add_argument("-svg", "--svg", dest="svg", action="store_true", required=False, help='Write figure to svg')
    options.add_argument("-html", "--html", dest="html", action="store_true", required=False,
                         help='Write figure to html')
    options.add_argument("-output-name", "-output-name", dest="output_base_name", type=str, default=None,
                         help='Base name to be used for output')

    options.add_argument("-compact", "--compact-subplot", dest="compact", action="store_true", required=False,
                         help='Compact plot using a subplot for each part field')
    options.add_argument("-vf", "--version-filter", dest="version_prefix", required=False, help='Filter on version')
    options.add_argument("-gpn", "--gpn-filter", dest="gpn_filter", required=False, help='Filter on GPN', type=str,
                         nargs='+', default=None)
    args = options.parse_args()

    save_pdf = False
    if args.pdf:
        save_pdf = True

    save_svg = False
    if args.svg:
        save_svg = True

    hide_plot = False
    if args.hide:
        hide_plot = True

    web_map = WebMap(0, 0)
    web_groups = WebGroups()

    if args.file:
        if args.file.endswith(".zip"):
            with zipfile.ZipFile(args.file, 'r') as zip:
                for fname in zip.namelist():
                    if fname.endswith("TASKDATA.XML") or fname.endswith("TASKDATA.xml"):
                        if fname.endswith("TASKDATA.xml") or fname != 'TASKDATA/TASKDATA.XML':
                            print("Invalid case in filename: '%s'" % fname, file=sys.stderr)
                        with zip.open(fname) as f:
                            print(fname)
                            tree = ET.parse(f)
                            show_task_file(args.version_prefix, fname.replace('/', '_'), tree, save_pdf,
                                           gpn_filter=args.gpn_filter, save_svg=save_svg, hide=hide_plot,
                                           use_subplot=args.compact, web_map=web_map,
                                           output_base_name=args.output_base_name, groups=web_groups)
        else:
            if args.file.endswith("TASKDATA.xml"):
                print("Invalid case in filename: '%s'" % args.file, file=sys.stderr)
            tree = ET.parse(args.file)
            show_task_file(args.version_prefix, Path(args.file).name, tree, save_pdf, gpn_filter=args.gpn_filter,
                           save_svg=save_svg, hide=hide_plot, use_subplot=args.compact, web_map=web_map,
                           output_base_name=args.output_base_name, groups=web_groups)

        if args.html:
            web_map.save("map.html")


def show_task_file(version_prefix, name, tree, save_pdf: bool = False, save_svg: bool = False, hide=False,
                   gpn_filter=None, use_subplot=False, web_map=None, output_base_name=None, groups: WebGroups = None):
    root = tree.getroot()

    if version_prefix is not None:
        if version_prefix not in (root.attrib.get("VersionMajor") + "." + root.attrib.get("VersionMinor")):
            return

    parent_map = {c: p for p in tree.iter() for c in p}
    ref_point_element = root.find(".//PNT").iter().__next__()
    ref = pnt_to_pair(ref_point_element)

    if web_map is not None:
        web_map.set_refernce(ref[0], ref[1])
    else:
        web_map = WebMap(ref[0], ref[1])

    if not groups:
        groups = WebGroups()

    part_fields = root.findall(".//PFD") + root.findall(".//TSK")

    if not use_subplot or len(part_fields) == 1:
        fig = plt.figure()
        axs = plt.gca()
        part_fields_ax = zip([axs for i in range(0, len(part_fields))], part_fields)
    else:

        n = len(part_fields)
        if n > 2:
            cols = 3
        else:
            cols = 2
        fig, axes = plt.subplots(nrows=round(len(part_fields) / cols), ncols=cols)
        part_fields_ax = zip(axes.flat, part_fields)

    if not groups.field_name_group:
        groups.field_name_group = web_map.create_group("Field names")

    for (ax, pfd) in part_fields_ax:
        if use_subplot:
            ax.title.set_text(pfd.attrib.get("C"))
        plot_all_pln(ax, parent_map, web_map, ref, pfd, groups.polygon_type_groups)
        plot_all_lsg(ax, parent_map, web_map, ref, pfd, gpn_filter=gpn_filter, line_type_groups=groups.line_type_groups)
        plot_all_pnt(ax, parent_map, web_map, ref, pfd)
        plot_center_name(name, pfd, ref, web_map, group=groups.field_name_group)

        ax.axis("equal")
        ax.axis("off")
    fig.tight_layout()
    # plt.legend(loc="upper left")

    if not output_base_name:
        output_base_name = name

    if save_pdf:
        plt.savefig(output_base_name + ".pdf")

    if save_svg:
        plt.savefig(output_base_name + ".svg")

    if not hide:
        plt.show()

    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def plot_center_name(name, pfd, ref, web_map, group=None):
    points_elements = pfd.findall('.//PNT')
    if len(points_elements) == 0:
        return
    point_data = [pnt_to_pair(pelement) for pelement in points_elements]
    points = [pymap3d.geodetic2enu(p[0], p[1], 0, ref[0], ref[1], 0, ell=ell_wgs84, deg=True) for p in
              point_data]
    field_name = str(pfd.attrib.get("C"))
    web_map.add_marker(name + '<br/> ' + field_name, MultiPoint(points).centroid, group=group)


def get_line_points(ref, line: ET.ElementTree):
    points_elements = line.findall('./PNT')
    if len(points_elements) == 0:
        return []
    point_data = [pnt_to_pair(pelement) for pelement in points_elements]
    points = [pymap3d.geodetic2enu(p[0], p[1], 0, ref[0], ref[1], 0, ell=ell_wgs84, deg=True) for p in
              point_data]
    return points


def get_polygon(ref, pln):
    if not pln:
        return None
    exterior_line = pln.find('./LSG[@A="1"]')

    points = get_line_points(ref, exterior_line)
    if len(points) == 0:
        # not implemented for polygons with no exterior
        return None
    interior_points = [l for l in [get_line_points(ref, l) for l in pln.findall('./LSG[@A="2"]')] if l is not None]
    exterior_points = [(p[0], p[1]) for p in points]

    interiors = [SHP.Polygon([[p[0], p[1]] for p in g]) for g in interior_points]

    return SHP.Polygon(exterior_points, [p.exterior.coords for p in interiors])


def plot_all_pln(ax, parent_map, web_map, ref, root, polygon_type_groups):
    for pln in root.findall(".//PLN"):
        designator = pln.attrib.get("C")
        print("Processing line '%s'" % designator)

        polygon = get_polygon(ref, pln)
        if polygon is None:
            continue

        polygon_type = int(pln.attrib.get("A"))
        group = None

        group_names = {1: 'Boundary', 2: 'Treatment zone', 3: 'Water', 4: 'PrimaryArea', 5: 'Road', 6: 'Obstacles',
                       9: 'Mainfield', 8: 'Other', 10: 'Headland', 11: 'BufferZone', 12: 'Windbreak'}

        if polygon_type in [1, 2, 3, 6, 8, 9, 10]:
            if polygon_type not in polygon_type_groups:
                polygon_type_groups[polygon_type] = web_map.create_group(group_names[polygon_type])
            group = polygon_type_groups[polygon_type]
        else:
            if polygon_type not in polygon_type_groups:
                polygon_type_groups[polygon_type] = web_map.create_group("Other polygons")
            group = polygon_type_groups[polygon_type]

        if polygon_type == 1:  # boundary
            patch = PolygonPatch(polygon.buffer(0), alpha=1, zorder=2, facecolor="black", linewidth=2, fill=False)
            web_map.add(polygon, tooltip=designator, style={'color': 'black', 'fillOpacity': '0'}, group=group)
        elif polygon_type == 2:  # treatmentzone
            patch = PolygonPatch(polygon.buffer(0), linewidth=2, facecolor="gray", alpha=0.1,
                                 hatch="...", fill=False)
            if 'B' in parent_map[pln].attrib:
                designator = parent_map[pln].attrib['B']

            color = get_color(parent_map[pln].attrib, 'C', 'gray')

            web_map.add(polygon, tooltip=designator, style={'color': color, 'fillOpacity': '0.1'}, group=group)
        elif polygon_type == 3:  # water
            patch = PolygonPatch(polygon, linewidth=2, facecolor="blue", alpha=0.1)
            web_map.add(polygon, tooltip=designator, style={'color': 'blue', 'fillOpacity': '0.1'}, group=group)
        elif polygon_type == 6:  # obstacle
            patch = PolygonPatch(polygon.buffer(0), alpha=0.1, zorder=2, facecolor="red")
            web_map.add(polygon, tooltip=designator, style={'color': 'red'}, group=group)
        elif polygon_type == 8:  # other
            patch = PolygonPatch(polygon.buffer(0), alpha=0.1, zorder=2, facecolor="gray")
            web_map.add(polygon, tooltip=designator, style={'color': 'gray', 'fillOpacity': '0.5'}, group=group)
        elif polygon_type == 9:  # mainland
            patch = PolygonPatch(polygon.buffer(0), alpha=0.2, zorder=2, facecolor="forestgreen", linewidth=0)
            web_map.add(polygon, tooltip=designator, style={'color': 'forestgreen', 'fillOpacity': '0.2'}, group=group)
        elif polygon_type == 10:  # headland
            patch = PolygonPatch(polygon.buffer(0), alpha=0.1, zorder=2, facecolor="springgreen")
            web_map.add(polygon, tooltip=designator,
                        style={'color': 'springgreen', 'opacity': '0.3', 'fillOpacity': '0.1', 'z-index': '2'},
                        group=group)
        else:
            patch = PolygonPatch(polygon.buffer(0), alpha=0.1, zorder=2, facecolor="violet")
            web_map.add(polygon, tooltip=designator, style={'color': 'violet', 'fillOpacity': '0.1'}, group=group)

        if type(patch) is PatchCollection:
            ax.add_collection(patch)
        else:
            ax.add_patch(patch)


def get_color(attrib, key, default_colour):
    if key in attrib:
        c = int(attrib.get(key))
        if c in colour_map:
            return colour_map[c]
    return default_colour


def get_pnts(node, ref, pnt_type_filter=None) -> [shapely.geometry.Point]:
    points_elements = node.findall("./PNT") if not pnt_type_filter else node.findall(
        "./PNT[@A='%s']" % pnt_type_filter)
    point_data = [pnt_to_pair(pelement) for pelement in points_elements]
    return [pymap3d.geodetic2enu(p[0], p[1], 0, ref[0], ref[1], 0, ell=ell_wgs84, deg=True) for p in
            point_data]

def plot_all_pnt(ax, parent_map, web_map, ref, root, line_type_groups=None, gpn_filter=None):
    for pnt in root.findall("./PNT"):
        print("Processing line '%s'" % pnt.attrib.get("B") if 'B' in pnt.attrib else 'No Designator')
        p_type = int(pnt.attrib.get("A"))
        if p_type not in [1,2,3,4,5,10,11]:
            continue

        point_data = [pnt_to_pair(pelement) for pelement in [pnt]]
        point=SHP.Point( [pymap3d.geodetic2enu(p[0], p[1], 0, ref[0], ref[1], 0, ell=ell_wgs84, deg=True) for p in
                point_data][0])

        web_map.add_point(pnt.attrib.get("B") if 'B' in pnt.attrib else 'Type %d' % p_type,
                            point, style={'color': get_color(pnt.attrib, 'F', 'black')}, group=None)
        ax.plot(point.x, point.y)



def plot_all_lsg(ax, parent_map, web_map, ref, root, line_type_groups, gpn_filter=None):
    for line in root.findall(".//LSG"):
        print("Processing line '%s'" % line.attrib.get("B"))

        points = get_pnts(line, ref)

        type = int(line.attrib.get("A"))
        designator = line.attrib.get("B")

        if type == 1 or type == 2:  # polygon exterior and interior
            pass
        elif type == 5:
            # this is guidance so check if replication

            parent = parent_map[line]

            if type not in line_type_groups:
                line_type_groups[type] = web_map.create_group("GuidanceLines")
            guidance_plot_group = line_type_groups[type]

            if parent.tag == "GPN":

                if gpn_filter is not None and parent.attrib.get("A") not in gpn_filter:
                    continue

                if not "C" in parent.attrib.keys():
                    print("Invalid GPN: " + str(parent))
                    continue

                gpn_type = int(parent.attrib.get("C"))

                if not gpn_type in [1, 2, 3, 4, 5]:
                    print("GPN %d not implemented" % gpn_type)
                    continue

                boundary_polygons = [get_polygon(ref, p) for p in parent_map[parent].findall('./PLN') if p is not None]
                boundary_polygons += [get_polygon(ref, p) for p in parent.findall('.//PLN') if p is not None]
                designator = parent.attrib.get('B')

                if gpn_type == 4:
                    center_pnts = get_pnts(line, ref, '8')
                    a_pnts = get_pnts(line, ref, '6')
                    b_pnts = get_pnts(line, ref, '7')

                    # if len(a_pnts)>0:
                    #     web_map.add_marker(designator + '-A-pivot', shapely.geometry.Point(a_pnts[0]), group=guidance_plot_group)
                    # if len(b_pnts)>0:
                    #     web_map.add_marker(designator + '-B-pivot', shapely.geometry.Point(b_pnts[0]), group=guidance_plot_group)

                    if not 'H' in parent.attrib.keys() or len(center_pnts) == 0:
                        print('Invalid pivot')
                        continue

                    pivot_cutout = None
                    if len(a_pnts) > 0 and len(b_pnts) > 0:
                        # we have a cut out to lets add it as a boundary limitation
                        from shapely.geometry import Polygon
                        pivot_cutout = Polygon(b_pnts + center_pnts + a_pnts + b_pnts)

                    radius = int(parent.attrib.get("H")) / 1000.0

                    center = shapely.geometry.Point(center_pnts[0])

                    # outer perimiter of pivot
                    poly = center.buffer(radius)

                    def calculate_angle(point1, point2):
                        x1, y1 = point1.x, point1.y
                        x2, y2 = point2.x, point2.y

                        angle = math.atan2(y2 - y1, x2 - x1)
                        angle_degrees = math.degrees(angle)

                        return angle_degrees

                    def is_angle_between(angle, start_angle, end_angle):
                        # Normalize angles to be within the range [0, 360)
                        angle = angle % 360
                        start_angle = start_angle % 360
                        end_angle = end_angle % 360

                        # Handle the case where end_angle is smaller than start_angle
                        if start_angle <= end_angle:
                            return start_angle <= angle <= end_angle
                        else:
                            # Angle range wraps around the circle (e.g., start_angle = 330, end_angle = 30)
                            return start_angle <= angle or angle <= end_angle

                    def filter_pivot_line_ab(line):
                        if len(a_pnts) > 0 and len(b_pnts) > 0:
                            a_angle = calculate_angle(center, shapely.geometry.Point(a_pnts[0]))
                            b_angle = calculate_angle(center, shapely.geometry.Point(b_pnts[0]))
                            # print('A %f B %f' % (a_angle, b_angle))
                            # print([calculate_angle(center, shapely.geometry.Point(p)) for p in line.coords     ])
                            # We should probably add the a and b points to the line as its will make sure it terminates at the a and b angles
                            return [p for p in line.coords if
                                    not is_angle_between(calculate_angle(center, shapely.geometry.Point(p)), a_angle,
                                                         b_angle)]
                        else:
                            return line

                    base_line_string = LineString(filter_pivot_line_ab(poly.exterior))

                    patch = LineCollection([base_line_string.coords], linewidths=1.5,
                                           edgecolors="black", zorder=7)
                    ax.add_collection(patch)

                    web_map.add(base_line_string, tooltip=designator + '-pivot-boundary', style={'color': 'black'},
                                group=guidance_plot_group)

                    if "C" in line.attrib.keys():
                        width = int(line.attrib.get("C")) / 1000 if "C" in line.attrib.keys() else 1

                        lines = []

                        for offset in range(1, int((radius / width) + 1)):
                            offset_line = LineString(
                                filter_pivot_line_ab(center.buffer(offset * width - (width / 2.0)).exterior))
                            if isinstance(offset_line, MultiLineString):
                                for line in offset_line:
                                    lines.append(line)
                            else:
                                lines.append(offset_line)

                        if len(lines) > 0:
                            guidance_lines = [extract_lines_within(line, boundary_polygons) for line in lines]
                            # if pivot_cutout:
                            #     guidance_lines = [extract_lines_within(line,[ pivot_cutout],invert=True) for line in lines]

                            patchc = LineCollection(
                                get_coordinates([item for sublist in guidance_lines for item in sublist]),
                                linewidths=1,
                                edgecolors="purple", zorder=5, alpha=0.5)

                            ax.add_collection(patchc)

                            g = web_map.create_group('GuidanceLines-' + str(designator) + '_replicated_GPNs')
                            for trimmed_line in [item for sublist in guidance_lines for item in sublist]:
                                web_map.addPoly(trimmed_line, tooltip=designator,
                                                style={'color': 'purple', 'z-index': '0', 'opacity': '0.3'},
                                                group=g)

                    """GuidancePatternOptions D
                    GuidancePatternRadius H
                    
                    pnt 8 for center
                    6 for A and 7 fot B"""

                    continue

                if gpn_type == 2:
                    # calculate second point
                    if "G" not in parent.attrib.keys():
                        print("A+ missing angle")
                        continue
                    angle = float(parent.attrib.get("G"))
                    length = 10
                    endy = points[0][1] + length * math.sin(math.radians(angle))
                    endx = points[0][0] + length * math.cos(math.radians(angle))
                    p = (endx, endy)
                    points.append(p)

                base_line_string = LineString([(p[0], p[1]) for p in points])

                if 'F' in parent.attrib.keys():
                    # ok extensions are enabled
                    extension_type = int(parent.attrib.get("F"))

                    # we dont know by how much
                    def getExtrapoledLine(p1, p2, ext_length, from_start=False):
                        print("Extending line  from start %s" % str(from_start))

                        direction_vector = (p2[0] - p1[0], p2[1] - p1[1])

                        if from_start:
                            ext_length = -ext_length

                        t = ext_length / (direction_vector[0] ** 2 + direction_vector[1] ** 2) ** 0.5

                        base = p1 if from_start else p2

                        point_at_distance = shapely.geometry.Point(
                            base[0] + t * direction_vector[0],
                            base[1] + t * direction_vector[1]
                        )

                        return LineString([p1, point_at_distance] if from_start else [p2, point_at_distance])

                    extension_length = 50  # base_line_string.length / 10.0
                    ext_a_length = extension_length
                    ext_b_length = extension_length

                    # maybe we have a boundary which sounds like a good idea to extend it to
                    if parent_map[parent_map[parent]].tag == "PFD":
                        boundary = get_polygon(ref, parent_map[parent_map[parent]].find('.//PLN[@A="1"]'))
                        if boundary:
                            base_line_string.intersection(boundary)
                        # if boundary:
                        #     first, last = base_line_string.boundary
                        #     ext_a_length = first.distance(boundary.exterior)
                        #     ext_b_length = last.distance(boundary.exterior)
                    from shapely.ops import linemerge
                    if extension_type == 1 or extension_type == 2:
                        # extend both ends
                        base_line_string = linemerge(
                            [getExtrapoledLine(*base_line_string.coords[:2], ext_a_length, from_start=True),
                             base_line_string])
                        # base_line_string = extend_line(base_line_string, ext_a_length, extend_start=True)
                    if extension_type == 1 or extension_type == 3:
                        # extend both ends
                        # base_line_string = extend_line(base_line_string, ext_b_length, extend_start=False)
                        base_line_string = linemerge([base_line_string,
                                                      getExtrapoledLine(*base_line_string.coords[-2:], ext_b_length)])

                number_of_swaths_left = 0
                number_of_swaths_right = 0

                propagation = parent.attrib.get("E")
                if propagation:
                    propagation_num = int(propagation)
                    if propagation_num == 1:
                        number_of_swaths_left = default_propagation_num
                        number_of_swaths_right = default_propagation_num
                    elif propagation_num == 2:
                        number_of_swaths_left = default_propagation_num
                    elif propagation_num == 3:
                        number_of_swaths_right = default_propagation_num

                if parent.attrib.get("N"):
                    number_of_swaths_left = int(parent.attrib.get("N"))
                if parent.attrib.get("O"):
                    number_of_swaths_right = int(parent.attrib.get("O"))

                lines = []

                if "C" in line.attrib.keys():
                    width = int(line.attrib.get("C")) / 1000 if "C" in line.attrib.keys() else 1

                    if number_of_swaths_left > 0:
                        for offset in range(1, number_of_swaths_left + 1):
                            offset_line = base_line_string.parallel_offset(width * offset, 'left',
                                                                           join_style=JOIN_STYLE.mitre)

                            if isinstance(offset_line, MultiLineString):
                                for line in offset_line:
                                    lines.append(line)
                            else:
                                lines.append(offset_line)
                    if number_of_swaths_right > 0:
                        for offset in range(1, number_of_swaths_right + 1):
                            offset_line = base_line_string.parallel_offset(width * offset * -1, 'left',
                                                                           join_style=JOIN_STYLE.mitre)
                            if isinstance(offset_line, MultiLineString):
                                for line in offset_line:
                                    lines.append(line)
                            else:
                                lines.append(offset_line)

                if len(lines) > 0:
                    guidance_lines = [extract_lines_within(line, boundary_polygons) for line in lines]

                    patchc = LineCollection(get_coordinates([item for sublist in guidance_lines for item in sublist]),
                                            linewidths=1,
                                            edgecolors="purple", zorder=5, alpha=0.5)

                    ax.add_collection(patchc)

                    g = web_map.create_group('GuidanceLines-' + str(designator) + '_replicated_GPNs')
                    for trimmed_line in [item for sublist in guidance_lines for item in sublist]:
                        web_map.addPoly(trimmed_line, tooltip=designator,
                                        style={'color': 'purple', 'z-index': '0', 'opacity': '0.3'},
                                        group=g)

                # plot the baseline
                # https://stackoverflow.com/questions/19877666/add-legends-to-linecollection-plot
                patch = LineCollection(get_coordinates(extract_lines_within(base_line_string, boundary_polygons)),
                                       linewidths=1.5,
                                       edgecolors="goldenrod", zorder=7)
                ax.add_collection(patch)

                for l in extract_lines_within(base_line_string, boundary_polygons):
                    web_map.add(l, tooltip=designator, style={'color': 'goldenrod'}, group=guidance_plot_group)

            else:
                if len(points) < 2:
                    print("Too few points in line skipping: " + str(line))
                    continue
                    # isoxml 3 guidance line

                base_line_string = LineString([(p[0], p[1]) for p in points])
                patch = LineCollection(get_coordinates([base_line_string]), linewidths=1.5,
                                       edgecolors="goldenrod", zorder=7)
                ax.add_collection(patch)

                web_map.add(base_line_string, tooltip=designator, style={'color': 'goldenrod'},
                            group=guidance_plot_group)
        elif type == 9:  # obstacle
            # its a polygon
            # polygon = SHP.Polygon(points)
            # patch = PatchCollection([polygon], linewidths=1, edgecolors="none", facecolor="red", alpha=0.3)
            # ax.add_patch(patch)
            base_line_string = LineString([(p[0], p[1]) for p in points])

            if type not in line_type_groups:
                line_type_groups[type] = web_map.create_group("Obstacles")
            group_obstacles = line_type_groups[type]

            if "C" in line.attrib.keys():
                # it is a polygon as a width is given
                width = int(line.attrib.get("C")) / 1000 if "C" in line.attrib.keys() else 1

                poly = shapely.concave_hull(base_line_string.parallel_offset(width / 2.0, 'left',
                                                        join_style=JOIN_STYLE.mitre).union(
                    base_line_string.parallel_offset(width / 2.0, 'right',
                                                     join_style=JOIN_STYLE.mitre)))

                web_map.add(poly, tooltip=designator, style={'color': get_color(line.attrib, 'E', 'red')},
                            group=group_obstacles)

                patch = PolygonPatch(poly, alpha=0.1, zorder=2, facecolor="red")
                ax.add_patch(patch)


            else:
                # it is a line
                patch = LineCollection(get_coordinates([base_line_string]), linewidths=1.5,
                                       edgecolors="red", zorder=7)

                web_map.add(base_line_string, tooltip=designator, style={'color': get_color(line.attrib, 'E', 'red')},
                            group=group_obstacles)
                ax.add_collection(patch)
        elif type == 3:  # tramlines
            color = get_color(line.attrib, 'E', 'brown')
            designator = line.attrib.get("B")
            base_line_string = LineString([(p[0], p[1]) for p in points])
            print(base_line_string.length)
            if designator:
                ax.plot([p[0] for p in points], [p[1] for p in points], label=designator, color=color)
            else:
                ax.plot([p[0] for p in points], [p[1] for p in points], color=color)

            if type not in line_type_groups:
                line_type_groups[type] = web_map.create_group("Tramlines")
            group_tramlines = line_type_groups[type]
            web_map.add(base_line_string, tooltip=designator, style={'color': color}, group=group_tramlines)
        else:
            color = get_color(line.attrib, 'E', 'gray')
            if designator:
                ax.plot([p[0] for p in points], [p[1] for p in points], label=designator, color=color)
            else:
                ax.plot([p[0] for p in points], [p[1] for p in points], color=color)
            web_map.add(LineString([(p[0], p[1]) for p in points]), tooltip=designator, style={'color': color})


if __name__ == '__main__':
    main()
