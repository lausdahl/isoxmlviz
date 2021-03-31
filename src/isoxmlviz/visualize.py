import argparse
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pymap3d
import shapely.geometry as SHP
from descartes.patch import PolygonPatch
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from shapely.geometry import LineString, JOIN_STYLE

from isoxmlviz.LineStringUtil import extract_lines_within

ell_wgs84 = pymap3d.Ellipsoid('wgs84')

default_propagation_num = 100


def pnt_to_pair(element: ET):
    # C=lat
    return float(element.attrib.get("C")), float(element.attrib.get("D"))


def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a


def main():
    options = argparse.ArgumentParser(prog="isoxmlviz")
    options.add_argument("-file", dest="file", type=str, required=True, help='Path to a isoxml task file XML or ZIP')
    options.add_argument("-p", "--pdf", dest="pdf", action="store_true", required=False, help='Write figure to pdf')
    args = options.parse_args()

    save_pdf = False
    if args.pdf:
        save_pdf = True

    if args.file:
        if args.file.endswith(".zip"):
            with zipfile.ZipFile(args.file, 'r') as zip:
                for fname in zip.namelist():
                    if fname.endswith("TASKDATA.XML"):
                        with zip.open(fname) as f:
                            print(fname)
                            tree = ET.parse(f)
                            show_task_file(fname.replace('/', '_'), tree, save_pdf)
        else:
            tree = ET.parse(args.file)
            show_task_file(Path(args.file).name, tree, save_pdf)


def show_task_file(name, tree, save_pdf: bool):
    root = tree.getroot()
    parent_map = {c: p for p in tree.iter() for c in p}
    ref_point_element = root.find(".//PNT").iter().__next__()
    ref = pnt_to_pair(ref_point_element)
    ax = plt.gca()
    plot_all_pln(ax, parent_map, ref, root)
    plot_all_lsg(ax, parent_map, ref, root)

    plt.legend(loc="upper left")
    ax.axis("equal")
    ax.axis("off")
    if save_pdf:
        plt.savefig(name + ".pdf")
    plt.show()


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


def plot_all_pln(ax, parent_map, ref, root):
    for pln in root.findall(".//PLN"):
        print("Processing line '%s'" % pln.attrib.get("B"))

        polygon = get_polygon(ref, pln)
        if polygon is None:
            continue

        polygon_type = int(pln.attrib.get("A"))

        if polygon_type == 1:  # boundary
            patch = PolygonPatch(polygon.buffer(0), alpha=1, zorder=2, facecolor="black", linewidth=2, fill=False)
        elif polygon_type == 2:  # treatmentzone
            patch = PolygonPatch([polygon], linewidth=2, facecolor="gray", alpha=0.1,
                                 hatch="...", fill=False)
        elif polygon_type == 3:  # water
            patch = PolygonPatch([polygon], linewidth=2, facecolor="blue", alpha=0.1)
        elif polygon_type == 6:  # obstacle
            patch = PolygonPatch(polygon.buffer(0), alpha=0.1, zorder=2, facecolor="red")
        elif polygon_type == 8:  # other
            patch = PolygonPatch(polygon.buffer(0), alpha=0.1, zorder=2, facecolor="gray")
        elif polygon_type == 9:  # mainland
            patch = PolygonPatch(polygon.buffer(0), alpha=0.2, zorder=2, facecolor="forestgreen", linewidth=0)
        elif polygon_type == 10:  # headland
            patch = PolygonPatch(polygon.buffer(0), alpha=0.1, zorder=2, facecolor="springgreen")
        else:
            patch = PolygonPatch(polygon.buffer(0), alpha=0.1, zorder=2, facecolor="violet")
        ax.add_patch(patch)


def plot_all_lsg(ax, parent_map, ref, root):
    for line in root.findall(".//LSG"):
        print("Processing line '%s'" % line.attrib.get("B"))
        points_elements = line.findall("./PNT")
        point_data = [pnt_to_pair(pelement) for pelement in points_elements]
        points = [pymap3d.geodetic2enu(p[0], p[1], 0, ref[0], ref[1], 0, ell=ell_wgs84, deg=True) for p in
                  point_data]

        type = int(line.attrib.get("A"))

        if type == 1 or type == 2:  # polygon exterior and interior
            pass
        elif type == 5:
            # this is guidance so check if replication

            parent = parent_map[line]

            base_line_string = LineString([(p[0], p[1]) for p in points])

            if parent.tag == "GPN":
                width = int(line.attrib.get("C")) / 1000
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

                boundary_polygons = [get_polygon(ref, p) for p in parent.findall('./PLN') if p is not None]

                lines = []

                if number_of_swaths_left > 0:
                    for offset in range(1, number_of_swaths_left + 1):
                        offset_line = base_line_string.parallel_offset(width * offset, 'left',
                                                                       join_style=JOIN_STYLE.mitre)
                        lines.append(offset_line)
                if number_of_swaths_right > 0:
                    for offset in range(1, number_of_swaths_right + 1):
                        offset_line = base_line_string.parallel_offset(width * offset * -1, 'left',
                                                                       join_style=JOIN_STYLE.mitre)
                        lines.append(offset_line)

                if len(lines) > 0:
                    trimmed_lines = [extract_lines_within(line, boundary_polygons) for line in lines]

                    patchc = LineCollection([item for sublist in trimmed_lines for item in sublist], linewidths=1,
                                            edgecolors="purple", zorder=5, alpha=0.5)

                    ax.add_collection(patchc)

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

                # patch = LineCollection([base_line_string], linewidths=1.5,
                #                        edgecolors="black", zorder=6,linestyle="--")
                # ax.add_collection(patch)
        elif type == 9:  # obstacle
            # its a polygon
            polygon = SHP.Polygon(points)
            patch = PatchCollection([polygon], linewidths=1, edgecolors="none", facecolor="red", alpha=0.3)
            ax.add_patch(patch)
        else:
            designator = line.attrib.get("B")
            if designator:
                ax.plot([p[0] for p in points], [p[1] for p in points], label=designator)
            else:
                ax.plot([p[0] for p in points], [p[1] for p in points])


if __name__ == '__main__':
    main()
