import argparse
import xml.etree.ElementTree as ET
import zipfile

import matplotlib
import numpy as np
import pymap3d
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Polygon
from shapely.geometry import LineString, JOIN_STYLE
from pathlib import Path

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
            tree = ET.parse(Path(args.file).name, args.file, save_pdf)
            show_task_file(tree)


def show_task_file(name, tree, save_pdf: bool):
    root = tree.getroot()
    parent_map = {c: p for p in tree.iter() for c in p}
    ref_point_element = root.find(".//PNT").iter().__next__()
    ref = pnt_to_pair(ref_point_element)
    for line in root.findall(".//LSG"):
        print("Processing line '%s'" % line.attrib.get("B"))
        points_elements = line.findall("./PNT")
        point_data = [pnt_to_pair(pelement) for pelement in points_elements]
        points = [pymap3d.geodetic2enu(p[0], p[1], 0, ref[0], ref[1], 0, ell=ell_wgs84, deg=True) for p in
                  point_data]

        type = int(line.attrib.get("A"))
        ax = plt.gca()

        if type == 1 or type == 2:
            # its a polygon
            patches = []
            polygon = Polygon([(p[0], p[1]) for p in points], True)
            patches.append(polygon)
            p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.1)
            ax.add_collection(p)
        elif type == 5:
            # this is guidance so check if replication

            parent = parent_map[line]

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

                line_string = LineString([(p[0], p[1]) for p in points])
                lines = []

                if number_of_swaths_left > 0:
                    for offset in range(1, number_of_swaths_left + 1):
                        offset_line = line_string.parallel_offset(width * offset, 'left',
                                                                  join_style=JOIN_STYLE.mitre)
                        lines.append(offset_line)
                if number_of_swaths_right > 0:
                    for offset in range(1, number_of_swaths_right + 1):
                        offset_line = line_string.parallel_offset(width * offset * -1, 'left',
                                                                  join_style=JOIN_STYLE.mitre)
                        lines.append(offset_line)
                lineseg = [[totuple(np.array(ls.xy)[:, i]) for i in range(0, int(np.array(ls.xy).size / 2))] for ls
                           in lines]
                # lineseg=[totuple(np.array(line_string.xy)[:,0]),totuple(np.array(line_string.xy)[:,1])]#[[(l.xy[0][0],l.xy[0][1]),(l.xy[1][0],l.xy[1][1])] for l in lines]

                if len(lines) > 0:
                    line_segments = LineCollection(lineseg, color="pink")
                    ax.add_collection(line_segments)

            designator = line.attrib.get("B")
            if designator:
                ax.plot([p[0] for p in points], [p[1] for p in points], color="goldenrod", label=designator)
            else:
                ax.plot([p[0] for p in points], [p[1] for p in points], color="goldenrod")
        elif type == 9:
            # its a polygon
            patches_obstacles = []
            polygon = Polygon([(p[0], p[1]) for p in points], True)
            patches_obstacles.append(polygon)
            p = PatchCollection(patches_obstacles, alpha=0.3, color="red")
            ax.add_collection(p)
        else:
            designator = line.attrib.get("B")
            if designator:
                ax.plot([p[0] for p in points], [p[1] for p in points], label=designator)
            else:
                ax.plot([p[0] for p in points], [p[1] for p in points])
        plt.legend(loc="upper left")
        ax.axis("equal")
        ax.axis("off")
    if save_pdf:
        plt.savefig(name + ".pdf")
    plt.show()


if __name__ == '__main__':
    main()
