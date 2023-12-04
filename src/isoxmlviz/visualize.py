import argparse
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
import sys
import pymap3d
import shapely.geometry as SHP
from descartes.patch import PolygonPatch
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from shapely.geometry import LineString, JOIN_STYLE, MultiLineString, MultiPoint
import math
from isoxmlviz.LineStringUtil import extract_lines_within
from isoxmlviz.webmap import WebMap

ell_wgs84 = pymap3d.Ellipsoid.from_name('wgs84')

default_propagation_num = 100

def generate_web_safe_colors():
    color_levels = [0, 51, 102, 153, 204, 255]
    web_safe_colors = [
        "#{:02X}{:02X}{:02X}{:02X}".format( r, g, b,a)
        for a in [255]#color_levels
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

for idx,c in enumerate(generate_web_safe_colors()):
    colour_map[idx+16] = c


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

        self.line_type_groups= {}
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

    part_fields = root.findall(".//PFD")

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

        group_names = {1: 'Boundary', 2: 'Treatment zone', 3: 'Water',4:'PrimaryArea',5:'Road', 6: 'Obstacles',9:'Mainfield', 8: 'Other', 10: 'Headland',11:'BufferZone',12:'Windbreak'}

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
            patch = PolygonPatch([polygon], linewidth=2, facecolor="blue", alpha=0.1)
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
        ax.add_patch(patch)


def get_color(attrib, key, default_colour):
    if key in attrib:
        c = int(attrib.get(key))
        if c in colour_map:
            return colour_map[c]
    return default_colour


def plot_all_lsg(ax, parent_map, web_map, ref, root, line_type_groups, gpn_filter=None):
    for line in root.findall(".//LSG"):
        print("Processing line '%s'" % line.attrib.get("B"))
        points_elements = line.findall("./PNT")
        point_data = [pnt_to_pair(pelement) for pelement in points_elements]
        points = [pymap3d.geodetic2enu(p[0], p[1], 0, ref[0], ref[1], 0, ell=ell_wgs84, deg=True) for p in
                  point_data]

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

                if not gpn_type in [1, 2, 3, 5]:
                    print("GPN %d not implemented" % gpn_type)
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

                boundary_polygons = [get_polygon(ref, p) for p in parent_map[parent].findall('.//PLN') if p is not None]

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
                    trimmed_lines = [extract_lines_within(line, boundary_polygons) for line in lines]

                    patchc = LineCollection([item for sublist in trimmed_lines for item in sublist], linewidths=1,
                                            edgecolors="purple", zorder=5, alpha=0.5)

                    ax.add_collection(patchc)

                    g = web_map.create_group('GuidanceLines-'+str(designator) + '_replicated_GPNs')
                    for trimmed_line in [item for sublist in trimmed_lines for item in sublist]:
                        web_map.addPoly(trimmed_line, tooltip=designator,
                                        style={'color': 'purple', 'z-index': '0', 'opacity': '0.3'},
                                        group=g)

                    # if len(boundary_polygons) > 0:
                    #     for patch in [PolygonPatch(ggg, alpha=0.1, zorder=6, facecolor="pink", linewidth=2, fill=False,
                    #                                hatch="...")
                    #                   for ggg in
                    #                   boundary_polygons]:
                    #         ax.add_patch(patch)
                #
                # designator = line.attrib.get("B")
                # if designator:
                #     ax.plot([p[0] for p in points], [p[1] for p in points], color="goldenrod", label=designator)
                # else:
                #     ax.plot([p[0] for p in points], [p[1] for p in points], color="goldenrod")

                # if len(boundary_polygons) > 0:
                #     for patch in [PolygonPatch(ggg, alpha=0.1, zorder=6, facecolor="pink", linewidth=2, fill=False,
                #                                hatch="...")
                #                   for ggg in
                #                   boundary_polygons]:
                #         ax.add_patch(patch)

                # https://stackoverflow.com/questions/19877666/add-legends-to-linecollection-plot
                patch = LineCollection(extract_lines_within(base_line_string, boundary_polygons), linewidths=1.5,
                                       edgecolors="goldenrod", zorder=7)
                ax.add_collection(patch)

                for l in extract_lines_within(base_line_string, boundary_polygons):
                    web_map.add(l, tooltip=designator, style={'color': 'goldenrod'}, group=guidance_plot_group)
                # patch = LineCollection([base_line_string], linewidths=1.5,
                #                        edgecolors="black", zorder=6,linestyle="--")
                # ax.add_collection(patch)
            else:
                if len(points) < 2:
                    print("Too few points in line skipping: " + str(line))
                    continue
                    # isoxml 3 guidance line

                base_line_string = LineString([(p[0], p[1]) for p in points])
                patch = LineCollection([base_line_string], linewidths=1.5,
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

                poly =base_line_string.parallel_offset(width /2.0, 'left',
                                                 join_style=JOIN_STYLE.mitre).union(base_line_string.parallel_offset(width /2.0, 'right',
                                                 join_style=JOIN_STYLE.mitre)).convex_hull

                web_map.add(poly, tooltip=designator, style={'color': get_color(line.attrib, 'E', 'red')},
                            group=group_obstacles)

                patch = PolygonPatch(poly, alpha=0.1, zorder=2, facecolor="red")
                ax.add_patch(patch)


            else:
                #it is a line
                patch = LineCollection([base_line_string], linewidths=1.5,
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
            web_map.add([p[0] for p in points], [p[1] for p in points], tooltip=designator, style={'lineColor': color})


if __name__ == '__main__':
    main()
