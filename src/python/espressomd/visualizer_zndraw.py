"""
EspressoMD: Visualizer with zincware/ZnDraw

This is an extremely simple utilization of ZnDraw to visualize positional data
of Espresso. The speed is still extremely slow for large quantities of
particles.

-- Plans to be implemented
- automatic detection of position/size
- atomic bonds
- obstacles
- speed improvement

----------------
Known Problems

For some reason google chrome crashes the first time the visualizer is
initialized. To fix this simply run the simulation again while the browser
where it previously crashed is still open.

(Especially chrome)

----------------
-- How to use --                                                                                                                                  -Before the simulation loop, insert the setup line:                                                                                              vis = visualizer_zndraw.setup_vis()

-First you need to install zndraw on your machine. To do this run the following
-command in your console

pip install zndraw

-In the simulation where you want to visualize your solutions, you need to 
-import the module "visualizer_zndraw". 

import espressomd.visualizer_zndraw as visual

-Before the simulation loop, insert the setup line:

vis = visual.setup() 

-If you have a state you want to visualize use                                                                                                    

visual.update_single(vis, system, size)  
                                                                                                                                                                                                     
-The size parameter is optional, if left empty all particles appear in the same                                              
-color an size. If you want to customize the size for each particle, input an                                                
-numpy array where you assign each particle an integer size. The length of the size-array                                    
-must equal the number of particles or enter an integer size that all particles
-will recieve.                                                                                                                                                                                                                      
"""      

import numpy as np
import espressomd
from ase import Atoms
from zndraw import ZnDraw


def setup():
    """
    Sets up the environment for starting ZnDraw

    Returns
    -------
    visZn : ZnDraw()
        returns Visualizer object ZnDraw()

    """
    visZn = ZnDraw()
    visZn.socket.sleep(2)
    return visZn

def espresso_to_atoms(pos, size=None):
    """
    Converts positional and size data into an ase.atoms object

    Parameters
    ----------
    pos : np.ndarray
        position of particles
    size : np.array, optional
        size of particles. The default is None. If none then all particles
        have the same size

    Returns
    -------
    atoms_obj : ase.atoms 
        Data is saved as an ase.atoms object

    """
    
    if size == None:
        size = np.ones(np.shape(pos)[0])
    elif type(size) == int:
        size = size*np.ones(np.shape(pos)[0])
        
    atoms_obj = Atoms(numbers=size, positions=pos)
    return atoms_obj

def update_single(visualizer, system, size=None):
    """
    Appends the input data as the last frame in the visualizer

    Parameters
    ----------
    visualizer : ZnDraw
        Visualizer Object from ZnDraw
    system : espressomd.system
        system object where we draw positional information from
    size : np.array, optional
        size of particles. The default is None. If none then all particles
        have the same size
    """
    pos = system.part.all().pos
    vis_atoms = espresso_to_atoms(pos, size)
    visualizer.append(vis_atoms)
    
#def vis_update_all_after(visualizer, pos, size=None):
    



    
    