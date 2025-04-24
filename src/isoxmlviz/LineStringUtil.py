from itertools import tee

# import shapely.geometry as SHP
from shapely.geometry import Point, Polygon, MultiPoint, LineString
from shapely.ops import linemerge
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

def get_coordinates(lines: list[LineString]):
    return [l.coords for l in lines]

def extract_lines_within(base_line: LineString, polygons: list[Polygon]):
    lines = [base_line]
    if len(polygons) == 0:
        return  lines

    for p in polygons:
        lines = [extract_line_within(line, p) for line in lines]
        lines = [item for sublist in lines for item in sublist]
    return lines


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
        for p in inter2.geoms:
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


def insert_intersection_in_polygon(polygon: Polygon, line: LineString):
    # Get the intersection points
    intersection = polygon.intersection(line)

    # Collect all intersection points
    if intersection.is_empty:
        print("No intersection")
        new_polygon = polygon
    else:
        intersection_points = []
        intersection_points.append(Point(intersection.coords[-1]))
        # if intersection.geom_type == "Point":
        #     intersection_points = [intersection]
        # elif intersection.geom_type == "MultiPoint":
        #     intersection_points = list(intersection.geoms)
        # elif intersection.geom_type in ["LineString", "MultiLineString"]:
        #     # Convert LineString(s) to Points
        #     merged = linemerge(intersection)
        #     if merged.geom_type == "LineString":
        #         intersection_points = [Point(merged.coords[0]), Point(merged.coords[-1])]
        #     elif merged.geom_type == "MultiLineString":
        #         intersection_points = []
        #         for ls in merged.geoms:
        #             intersection_points.append(Point(ls.coords[0]))
        #             intersection_points.append(Point(ls.coords[-1]))
        #     else:
        #         raise ValueError("Unexpected merged type:", merged.geom_type)
        # else:
        #     raise ValueError("Unexpected intersection type:", intersection.geom_type)

        # Original polygon coordinates (excluding the closing point)
        coords = list(polygon.exterior.coords[:-1])

        # Insert intersection points into the right place
        new_coords = []
        for i in range(len(coords)):
            a = Point(coords[i])
            b = Point(coords[(i + 1) % len(coords)])
            segment = LineString([a, b])
            new_coords.append((a.x, a.y))

            for ip in intersection_points:
                if segment.distance(ip) < 1e-8 and segment.project(ip) > 0 and segment.project(ip) < segment.length:
                    # Only insert if truly on the segment and not a vertex
                    if (ip.x, ip.y) not in new_coords:
                        new_coords.append((ip.x, ip.y))

        # Re-close the polygon
        new_coords.append(new_coords[0])
        new_polygon = Polygon(new_coords)
    return new_polygon


def densify_linestring(linestring, max_segment_length):
    coords = list(linestring.coords)
    new_coords = []

    for i in range(len(coords) - 1):
        start = Point(coords[i])
        end = Point(coords[i + 1])
        segment = LineString([start, end])
        length = segment.length

        num_segments = max(int(length // max_segment_length), 1)

        for j in range(num_segments):
            fraction = j / num_segments
            point = segment.interpolate(fraction, normalized=True)
            new_coords.append((point.x, point.y))

    new_coords.append(coords[-1])  # ensure the final point is added
    return LineString(new_coords)