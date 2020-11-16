# %%
import MDAnalysis as mda
import nglview as nv
from funcs import unwrap_traj, add_elements_to_universe
# %%
u = mda.Universe('pack.lmps', 'traj.dcd', topology_format='DATA')
add_elements_to_universe(u, 'pack.lmps')
unwrap_traj(u)
nv.show_mdanalysis(u)
