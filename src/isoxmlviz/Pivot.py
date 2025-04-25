from shapely.geometry import Point
import numpy as np
import math


def calculate_angle(point1, point2):
    x1, y1 = point1.x, point1.y
    x2, y2 = point2.x, point2.y

    angle = math.atan2(y2 - y1, x2 - x1)
    angle_degrees = math.degrees(angle)
    normalized_degrees = (angle_degrees + 360) % 360
    return normalized_degrees


def angle_from_center(center: Point, point: Point) -> float:
    """Return angle in degrees from center to point, in [0, 360)."""
    dx = point.x - center.x
    dy = point.y - center.y
    return np.degrees(np.arctan2(dy, dx)) % 360


def points_on_arc(center: Point, A: Point, B: Point, radius: float, delta_deg: float = 10.0) -> list:
    """Generate points on the arc from vector CA to CB, all equidistant from C."""
    r = radius  # center.distance(A)  # assume A and B are at same distance from C
    angle_A = angle_from_center(center, A)
    angle_B = angle_from_center(center, B)

    # Determine direction and sweep angle
    # clockwise
    angle_diff = (angle_A - angle_B) % 360
    clockwise = angle_diff > 180

    if clockwise:
        # Go clockwise from A to B
        angle_diff = (angle_A - angle_B) % 360
        num_steps = int(angle_diff // delta_deg)
        angles = [(angle_A - i * delta_deg) % 360 for i in range(num_steps + 1)]
    else:
        # Go counterclockwise from A to B
        num_steps = int(angle_diff // delta_deg)
        angles = [(angle_A + i * delta_deg) % 360 for i in range(num_steps + 1)]

    # Ensure the last point is exactly at B's angle
    if angles[-1] != angle_B:
        angles.append(angle_B)

    # Generate points at each angle
    points = [
        Point(center.x + r * np.cos(np.radians(angle)),
              center.y + r * np.sin(np.radians(angle)))
        for angle in angles
    ]
    return points


def pivot_arch_line_ab(center, a_pnts, b_pnts, line):
    if len(a_pnts) > 0 and len(b_pnts) > 0:

        arch_radius = center.distance(Point(line.coords[0]))

        a_point = Point(a_pnts[0])
        b_point = Point(b_pnts[0])

        delta_deg = 10
        line = points_on_arc(center, a_point, b_point, radius=arch_radius, delta_deg=delta_deg)
        while len(line) < 10:
            delta_deg = delta_deg / 2
            line = points_on_arc(center, a_point, b_point, radius=arch_radius, delta_deg=delta_deg)
        return line
    else:
        return line
