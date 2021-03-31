from itertools import tee

# import shapely.geometry as SHP
from shapely.geometry import Point, Polygon, MultiPoint, LineString

'''utility function to iterate a list in pairs'''


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


'''
https://stackoverflow.com/questions/11907947/how-to-check-if-a-point-lies-on-a-line-between-2-other-points
'''


def is_point_on_line(a, b, point):
    dxl = b[0] - a[0]
    dyl = b[1] - a[1]
    dxp = point[0] - a[0]
    dyp = point[1] - a[1]
    c = ((dxp * dxl + dyp * dyl) / (dxl * dxl + dyl * dyl))
    return 0 < c < 1


'''
extract line segments from a line and a list of polygon restricting the line to within the polygon. May result in many line segments
'''


def extract_lines_within(base_line: LineString, polygons: list[Polygon]):
    lines = [base_line]
    if len(polygons) == 0:
        return lines

    for p in polygons:
        lines = [extract_line_within(line, p) for line in lines]
    return [item for sublist in lines for item in sublist]


'''
extract line segments from a line and a polygon restricting the line to within the polygon. May result in many line segments
'''


def extract_line_within(base_line: LineString, polygon: Polygon):
    inter2 = base_line.intersection(polygon.boundary)
    lines = []

    # from matplotlib import pyplot as plt
    # plt.scatter([p.coords.xy[0] for p in inter2], [p.coords.xy[1] for p in inter2], color="red")

    expanded_points = []

    interpolated_points = []

    if type(inter2) is Point:
        interpolated_points.append(inter2)
    elif type(inter2) is MultiPoint:
        for p in inter2:
            interpolated_points.append(p)

    for a, b in pairwise(base_line.coords):
        expanded_points.append(a)

        # check for points between and add them with key as distance from a
        points_between = [(Point(a).distance(Point(p)), p) for p in interpolated_points if
                          is_point_on_line(a, b, p.coords[0]) and LineString((a, b)).distance(Point(p)) < 0.01]

        # sort the points between with distance from a
        points_between.sort(key=lambda x: x[0])

        for p in points_between:
            # add point between in the order of distance between
            expanded_points.append(p[1])
            interpolated_points = [g for g in interpolated_points if g != p[1]]

        expanded_points.append(b)

    for a, b in pairwise(expanded_points):
        line = LineString([a, b])

        mid_point = line.interpolate(0.5, normalized=True)
        if mid_point.within(polygon):
            lines.append(line)
    return lines
