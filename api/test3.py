from ipywidgets.embed import embed_minimal_html

import nglview
from nglview.contrib.movie import MovieMaker
import mdtraj
trajectory_LCO = "/Users/sachin/Documents/Coding/Self-Assembly/api/temp/LCO_mlff.heat.lammpstrj"
topology_LCO = "/Users/sachin/Documents/Coding/Self-Assembly/api/temp/LCO_pristine.mol2"
trj_LCO = mdtraj.load_lammpstrj(trajectory_LCO, top=topology_LCO)
view = nglview.show_mdtraj(trj_LCO)
view.render_image()
view._display_image()
embed_minimal_html('export.html', views=[view], title='Widgets export')