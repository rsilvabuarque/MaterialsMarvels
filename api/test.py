import nglview
from nglview.contrib.movie import MovieMaker
import mdtraj
from time import sleep
trajectory_LCO = "/Users/sachin/Documents/Coding/Self-Assembly/api/temp/LCO_mlff.heat.lammpstrj"
topology_LCO = "/Users/sachin/Documents/Coding/Self-Assembly/api/temp/LCO_pristine.mol2"
trj_LCO = mdtraj.load_lammpstrj(trajectory_LCO, top=topology_LCO)
view = nglview.show_mdtraj(trj_LCO)
# view
view.render_image()
# view._display_image()
# sleep(5)
nglview.write_html('index.html', [view], frame_range=(0, view.max_frame, 1)) 

# movie = MovieMaker(view, output='test2.gif', in_memory=True)
# movie.make()