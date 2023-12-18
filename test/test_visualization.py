import pytest
import glob
from pathlib import Path

from isoxmlviz.visualize import WebGroups, show_task_file
from isoxmlviz.webmap import WebMap
import xml.etree.ElementTree as ET


@pytest.mark.parametrize("filename",
                         [f for f in glob.glob(str(Path(__file__).parent / 'examples' / '*.xml'))],
                         ids=lambda params: f"{Path(params).name}")
def test_visualize(filename):
    output_dir = Path(__file__).parent / 'generated' / Path(filename).name
    output_dir.mkdir(exist_ok=True, parents=True)

    web_map = WebMap(0, 0)
    web_groups = WebGroups()

    tree = ET.parse(filename)
    show_task_file(None, str(Path(filename).name),
                   tree,
                   hide=False, save_pdf=True, save_svg=True, use_subplot=False,
                   web_map=web_map,
                   output_base_name=str(output_dir), groups=web_groups)
    web_map.save(str(output_dir / "map.html"))
