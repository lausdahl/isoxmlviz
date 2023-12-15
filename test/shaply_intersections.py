import unittest

import matplotlib.pyplot as plt
import shapely.geometry as SHP

from matplotlib.collections import LineCollection

from isoxmlviz.LineStringUtil import is_point_on_line, extract_line_within
from isoxmlviz.visualize import PolygonPatch


class ShapelyIntersectionTest(unittest.TestCase):

    def test_point_on_line(self):
        a = [0, 0]
        b = [1, 1]
        p = [0.5, 0.5]
        self.assertTrue(is_point_on_line(a, b, p))
        self.assertFalse(is_point_on_line(a, b, [0.4, 5.5]))

    def test_inter(self):
        polygon = SHP.Polygon([[0, 0], [0, 10], [5, 6], [10, 10], [10, 0]])
        line = SHP.LineString([(-5, -5), (5, 5)])
        fig, ax = plt.subplots()
        ax.add_collection(LineCollection([line], edgecolors="black"))
        ax.add_collection(
            LineCollection(extract_line_within(line, polygon), linewidths=10, edgecolors="black", alpha=0.3, zorder=3))

        line2s = SHP.LineString([(-5, 8), (15, 8)])
        ax.add_collection(LineCollection([line2s], edgecolors="black"))
        ax.add_collection(
            LineCollection(extract_line_within(line2s, polygon), linewidths=10, edgecolors="black", alpha=0.3,
                           zorder=3))

        line2s = SHP.LineString([(-5, 8), (8, 7)])
        ax.add_collection(LineCollection([line2s], edgecolors="black"))
        ax.add_collection(
            LineCollection(extract_line_within(line2s, polygon), linewidths=10, edgecolors="black", alpha=0.3,
                           zorder=3))

        line2s = SHP.LineString([(6, 5), (15, 7)])
        ax.add_collection(LineCollection([line2s], edgecolors="black"))
        ax.add_collection(
            LineCollection(extract_line_within(line2s, polygon), linewidths=10, edgecolors="black", alpha=0.3,
                           zorder=3))

        line2s = SHP.LineString([(2,2), (15, 2)])
        ax.add_collection(LineCollection([line2s], edgecolors="black"))
        ax.add_collection(
            LineCollection(extract_line_within(line2s, polygon), linewidths=10, edgecolors="black", alpha=0.3,
                           zorder=3))

        line2s = SHP.LineString([(2, -1), (15, 3)])
        ax.add_collection(LineCollection([line2s], edgecolors="black"))
        ax.add_collection(
            LineCollection(extract_line_within(line2s, polygon), linewidths=10, edgecolors="black", alpha=0.3,
                           zorder=3))

        patch = PolygonPatch(polygon.buffer(0), alpha=0.1)
        ax.add_patch(patch)
        # ax.add_collection(LineCollection([line, line2s], edgecolors="black"))
        #
        # inter = line.intersection(polygon.boundary)
        #
        # inter2 = line2s.intersection(polygon.boundary)
        #
        # lines = []
        #
        # c = []
        # for p in inter2:
        #     c.append(p)
        #
        # for a, b in pairwise(inter2):
        #     line = SHP.LineString([a, b])
        #
        #     mid_point = line.interpolate(0.5, normalized=True)
        #     if mid_point.within(polygon):
        #         lines.append(line)

        # cur = []
        # for p in inter2:
        #
        #     if len(cur) == 0:
        #         cur.append(p)
        #     else:
        #         l = SHP.LineString(cur + [p])
        #
        #         if not l.within(polygon.boundary):
        #             if len(cur)>1:
        #                 lines.append(SHP.LineString(cur ))
        #             cur=[]
        #         else:
        #             cur.append(p)
        # if len(cur)>1:
        #     lines.append(SHP.LineString(cur))

        # ax.add_collection(LineCollection(lines, linewidths=10, edgecolors="black", zorder=6))

        ax.set_xlim(-5, 15)
        ax.set_ylim(-5, 15)
        plt.show()
