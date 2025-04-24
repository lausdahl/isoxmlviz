import folium
import pymap3d

ell_wgs84 = pymap3d.Ellipsoid.from_name('wgs84')


def flatten(matrix):
    return [item for row in matrix for item in row]


def transform(shape, func):
    ''' Apply a function to every coordinate in a geometry.
    '''
    from shapely.geometry import Point, LineString, Polygon, GeometryCollection
    construct = shape.__class__

    if shape.geom_type.startswith('Multi'):
        parts = [transform(geom, func) for geom in shape.geoms]
        return construct(parts)

    if shape.geom_type in ('LineString'):
        return construct(map(func, shape.coords))
    if shape.geom_type in ('Point'):
        return construct(flatten(map(func, shape.coords)))

    if shape.geom_type == 'Polygon':
        exterior = map(func, shape.exterior.coords)
        rings = [list(map(func, ring.coords)) for ring in shape.interiors]
        return construct(exterior, rings)

    if shape.geom_type == 'GeometryCollection':
        return construct()

    raise ValueError('Unknown geometry type, "%s"' % shape.geom_type)


class WebMap:
    def __init__(self, ref_lat, ref_lng) -> None:
        super().__init__()
        self.ref_lat = ref_lat
        self.ref_lng = ref_lng
        self.m = folium.Map(location=[40.730610, -73.935242],
                            zoom_start=12, control_scale=True, prefer_canvas=True, max_zoom=22)

        folium.TileLayer('openstreetmap').add_to(self.m)
        tile = folium.TileLayer(
            tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
            attr='Esri',
            name='Esri Satellite',
            overlay=False,
            control=True,
            max_zoom=22
        ).add_to(self.m)

        # Add custom basemaps to folium
        basemaps = {
            'Google Maps': folium.TileLayer(
                tiles='https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}',
                attr='Google',
                name='Google Maps',
                overlay=False,
                control=True,
                show=False
            ),
            'Google Satellite': folium.TileLayer(
                tiles='https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
                attr='Google',
                name='Google Satellite',
                overlay=True,
                control=True,
                show= False
            ),
            'Google Terrain': folium.TileLayer(
                tiles='https://mt1.google.com/vt/lyrs=p&x={x}&y={y}&z={z}',
                attr='Google',
                name='Google Terrain',
                overlay=True,
                control=True
            ),
            'Google Satellite Hybrid': folium.TileLayer(
                tiles='https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}',
                attr='Google',
                name='Google Satellite',
                overlay=False,
                control=True,
                show=False
            ),
            'Esri Satellite': folium.TileLayer(
                tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
                attr='Esri',
                name='Esri Satellite',
                overlay=True,
                control=True
            )
        }

        def add_ee_layer(self, ee_object, vis_params, name):

            try:
                # display ee.Image()
                if isinstance(ee_object, ee.image.Image):
                    map_id_dict = ee.Image(ee_object).getMapId(vis_params)
                    folium.raster_layers.TileLayer(
                        tiles=map_id_dict['tile_fetcher'].url_format,
                        attr='Google Earth Engine',
                        name=name,
                        overlay=True,
                        control=True
                    ).add_to(self)
                # display ee.ImageCollection()
                elif isinstance(ee_object, ee.imagecollection.ImageCollection):
                    ee_object_new = ee_object.mosaic()
                    map_id_dict = ee.Image(ee_object_new).getMapId(vis_params)
                    folium.raster_layers.TileLayer(
                        tiles=map_id_dict['tile_fetcher'].url_format,
                        attr='Google Earth Engine',
                        name=name,
                        overlay=True,
                        control=True
                    ).add_to(self)
                # display ee.Geometry()
                elif isinstance(ee_object, ee.geometry.Geometry):
                    folium.GeoJson(
                        data=ee_object.getInfo(),
                        name=name,
                        overlay=True,
                        control=True
                    ).add_to(self)
                # display ee.FeatureCollection()
                elif isinstance(ee_object, ee.featurecollection.FeatureCollection):
                    ee_object_new = ee.Image().paint(ee_object, 0, 2)
                    map_id_dict = ee.Image(ee_object_new).getMapId(vis_params)
                    folium.raster_layers.TileLayer(
                        tiles=map_id_dict['tile_fetcher'].url_format,
                        attr='Google Earth Engine',
                        name=name,
                        overlay=True,
                        control=True
                    ).add_to(self)

            except:
                print("Could not display {}".format(name))

        # Add EE drawing method to folium.
        # folium.Map.add_ee_layer = add_ee_layer

        vis_params = {
            'min': 0,
            'max': 4000,
            'palette': ['006633', 'E5FFCC', '662A00', 'D8D8D8', 'F5F5F5']}

        # Create a folium map object.
        my_map = folium.Map(location=[20, 0], zoom_start=3, height=500)

        # Add custom basemaps
        basemaps['Google Maps'].add_to(self.m)
        basemaps['Google Satellite Hybrid'].add_to(self.m)

        # Add the elevation model to the map object.
        # self.m.add_ee_layer(dem.updateMask(dem.gt(0)), vis_params, 'DEM')



    def set_refernce(self, ref_lat, ref_lng):
        self.ref_lat = ref_lat
        self.ref_lng = ref_lng

    def add(self, geom, tooltip=None, style={'fillColor': 'red', 'lineColor': '#228B22', 'opacity': '1'}, group=None):
        mkt = lambda x: folium.features.Popup(tooltip) if tooltip else None

        ff = lambda c: [c[1], c[0], c[2]]
        folium.GeoJson(transform(geom, lambda p: ff(
            pymap3d.enu2geodetic(p[0], p[1], 0, self.ref_lat, self.ref_lng, 0, ell=ell_wgs84, deg=True))),
                       style_function=lambda x: style, tooltip=mkt(tooltip), popup=mkt(tooltip)).add_to(
            self.m if not group else group)

    def addPoly(self, geom, tooltip=None, style={'color': 'red', 'lineColor': '#228B22', 'opacity': '1'}, group=None):
        mkt = lambda x: folium.features.Popup(tooltip) if tooltip else None
        ff = lambda c: [c[0], c[1]]
        g = transform(geom, lambda p: ff(
            pymap3d.enu2geodetic(p[0], p[1], 0, self.ref_lat, self.ref_lng, 0, ell=ell_wgs84, deg=True)))

        ff = lambda c: [c[1], c[0], c[2]]
        # folium.GeoJson(transform(geom, lambda p: ff(
        #     pymap3d.enu2geodetic(p[0], p[1], 0, self.ref_lat, self.ref_lng, 0, ell=ell_wgs84, deg=True))),
        #                style_function=lambda x: style,tooltip=mkt(tooltip),popup=mkt(tooltip)).add_to(self.m)
        folium.PolyLine(locations=g.coords, tooltip=mkt(tooltip), color=style['color'],
                        opacity=style['opacity'] if 'opacity' in style else 1).add_to(self.m if not group else group)

    def add_marker(self, name, geom, group=None):
        ff = lambda c: [c[0], c[1], c[2]]
        g = transform(geom, lambda p: ff(
            pymap3d.enu2geodetic(p[0], p[1], 0, self.ref_lat, self.ref_lng, 0, ell=ell_wgs84, deg=True)))

        folium.Marker(location=[c[:-1] for c in g.coords][0], popup=folium.features.Popup(name)).add_to(
            self.m if not group else group)

    def add_point(self, name, geom, group=None, style={'color': 'red', 'lineColor': '#228B22', 'opacity': '1'}):
        ff = lambda c: [c[0], c[1], c[2]]
        g = transform(geom, lambda p: ff(
            pymap3d.enu2geodetic(p[0], p[1], 0, self.ref_lat, self.ref_lng, 0, ell=ell_wgs84, deg=True)))

        folium.CircleMarker(location=[c[:-1] for c in g.coords][0],
                            radius=2,
                            weight=5, color=style['color'],
                            opacity=style['opacity'] if 'opacity' in style else 1,
                            popup=folium.features.Popup(name)).add_to(self.m if not group else group)

    def create_group(self, name):
        return folium.FeatureGroup(name=name).add_to(self.m)

    def save(self, path):
        folium.LayerControl().add_to(self.m)
        self.m.fit_bounds(self.m.get_bounds(), padding=(30, 30))
        self.m._repr_html_()
        self.m.save(path)
