import openmc
import math
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

def build_model(radius, fuel):
   
   
    materials = openmc.Materials([fuel])
   
    # create sphere with radius parameter
    sphere_radius = openmc.Sphere(x0=0,y0=0,z0=0,r=radius, boundary_type='vacuum', name='BCM')
   
    # create core cell
    core_cell = openmc.Cell(name='Cylindrical Core')
    core_cell.fill = fuel
    core_cell.region = -sphere_radius
   
    # create universe geometry
    root_universe = openmc.Universe(name='root universe')
    root_universe.add_cells([core_cell])
   
    geometry = openmc.Geometry(root_universe)
    # define criticality settings
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue' # keff calculation
    settings.particles = 10000 # particles per batch
    settings.batches = 1050 # number of batches
    settings.inactive = 50 # inactive batches
   
    # settings.output = {'tallies': False}
   
   
    model = openmc.model.Model(geometry,materials,settings)
   
    return model

radius = 3
length = 3
