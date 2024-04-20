import openmc
import math
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import openmc.deplete
from IPython.display import Image
import xml.etree.ElementTree as et
from mpl_toolkits.mplot3d import Axes3D 
from pylab import *

# General Parameters

# Constants - via PNNL Mat'l Comp., Rev-2
ThO2OAtomRho = 0.044247 # Oxygen Atom Density in ThO2 [atom/b-cm]
ThO2ThAtomRho = 0.022124 # Thorium Atom Density in ThO2 [atom/b-cm]
ThMolW = 232. # Molecular Weight of Th [g/mol]
OMolW = 16. # Molecular Weight of O [g/mol]
U235MolW = 235. # Molecular Weight of U-235 [g/mol]
U238MolW = 238. # Molecular Weight of U-238 [g/mol]
ThO2MassRho = 9.7 # Mass Density of ThO2 [g/cm^3]
graphMassRho = 2.26 # Mass Density of Graphite [g/cm^3]
he_CoolMassRho = 0.000178 # Mass Density of He @ 293K [g/cm^3]
he_Cool_Z = 1.00652 # He Compression Factor @750C and 6 MPa
# MW*P/(Z*R*T) from https://www.cedengineering.com/userfiles/H02-011%20-%20Calculation%20of%20Gas%20Density%20and%20Viscosity.pdf
he_Cool_MW = 4.003 # g/mol
press_homog = 6000 #kPa
temp_homog = 750 + 273.15
R = 8.314472 #L*kPa/K*mol or m^3*Pa/K*mol
he_CoolMassRho = he_Cool_MW*press_homog/(he_Cool_Z*temp_homog)/1000000 # Mass Density of He @ 1073K and 6000 kPa [g/cm^3]
# Pressure density data from here https://backend.orbit.dtu.dk/ws/portalfiles/portal/52768509/ris_224.pdf
avo = 0.6022 # Avogadro's Number

# Volume Fraction Constants
vF_Pebbles_Core = 2./3 # Volume Fraction of Pebbles in Core (Max packing factor)
vF_TRISO_Pebbles = 0.15 # Volume Fraction of TRISO in Pebbles
vF_Fuel_TRISO = 0.15 # Volume Fraction of Fuel in TRISO
fuelVolFrac = vF_Pebbles_Core*vF_TRISO_Pebbles*vF_Fuel_TRISO # Volume Fraction of Fuel
graphVolFrac = vF_Pebbles_Core*(1 - (vF_TRISO_Pebbles*vF_Fuel_TRISO)) # Total Volume Fraction of Graphite
hel_CoolVolFrac = 1 - (fuelVolFrac + graphVolFrac) # Total Volume Fraction of Helium
totalVolFrac = fuelVolFrac + graphVolFrac + hel_CoolVolFrac # Total Volume

UO2MassRho = 10.45
enrichVal = 19.75/100
pctLEUVal = 60/100

# Potential 3D Plot
#x = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
#y = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
#z = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
#colo = [] 
  
# creating figures 
#fig = plt.figure(figsize=(10, 10)) 
#ax = fig.add_subplot(111, projection='3d') 
  
# setting color bar 
#color_map = cm.ScalarMappable(cmap=cm.Greens_r) 
#color_map.set_array(colo) 
  
# creating the heatmap 
#img = ax.scatter(x, y, z, marker='s', 
#                 s=200, color='green') 
#plt.colorbar(color_map) 
  
# adding title and labels 
#ax.set_title("Keffective Heat Map as a function of Thorium Mass %, TRISO packing factor, and Uranium Volume Fraction") 
#ax.set_xlabel('Thorium %') 
#ax.set_ylabel('Packing Factor') 
#ax.set_zlabel('Uranium Volume Fraction') 
  
# displaying plot 
#plt.show() 


## Calculate Weight Percentages ##
ThO2MolW = ThMolW + 2*OMolW # Molecular Weight of ThO2 [g/mol]
pctTh = 1 - pctLEUVal
# Total Mixture Density
UThMixMassRho = (pctLEUVal/UO2MassRho + pctTh/ThO2MassRho)**-1 # [g/cm]
# UO2 Atom Densities
UO2UMolW = ((enrichVal/U235MolW) + ((1 - enrichVal)/U238MolW))**-1 # Molecular Weight of U in UO2 [g/mol]
UO2MolW = UO2UMolW + 2*16 # Molecular Weight of UO2 [g/mol]
UO2OAtomRho = (UO2MassRho/UO2MolW)*2*avo # O Atom Density in UO2 [atom/b-cm]
UO2U235AtomRho = (enrichVal*UO2MassRho*(UO2UMolW/UO2MolW)/U235MolW)*avo # U-235 Atom Density in UO2 [atom/b-cm]
UO2U238AtomRho = ((1 - enrichVal)*UO2MassRho*(UO2UMolW/UO2MolW)/U238MolW)*avo # U-238 Atom Density in UO2 [atom/b-cm]
UO2U235AtomRhoFrac = UO2U235AtomRho/(UO2U235AtomRho+UO2U238AtomRho)
UO2U238AtomRhoFrac = UO2U238AtomRho/(UO2U235AtomRho+UO2U238AtomRho)
# UO2 ThO2 Mixture Atom Densities
UThMixOAtomRho = (pctTh*UThMixMassRho/ThO2MolW)*2*avo + (pctLEUVal*UThMixMassRho/UO2MolW)*2*avo # Oxygen Atom Density in Fuel Mixture [atom/b-cm]
UThMixThAtomRho = (pctTh*UThMixMassRho/ThO2MolW)*avo # Uranium Atom Density in Fuel Mixture [atom/b-cm]
UThMixU235AtomRho = (pctLEUVal*UThMixMassRho/UO2MolW)*avo*UO2U235AtomRhoFrac # Uranium Atom Density in Fuel Mixture [atom/b-cm]
UThMixU238AtomRho = (pctLEUVal*UThMixMassRho/UO2MolW)*avo*UO2U238AtomRhoFrac # Uranium Atom Density in Fuel Mixture [atom/b-cm]
# UO2 ThO2 Mixture Atom Fractions
# Normalized to 1
UThMixOAtomFrac = UThMixOAtomRho/(UThMixOAtomRho + UThMixThAtomRho + UThMixU235AtomRho + UThMixU238AtomRho) # Oxygen Atom Fraction in Fuel Mixture [atom/b-cm]
UThMixThAtomFrac = UThMixThAtomRho/(UThMixOAtomRho + UThMixThAtomRho + UThMixU235AtomRho + UThMixU238AtomRho) # Uranium Atom Fraction in Fuel Mixture [atom/b-cm]
UThMixU235AtomFrac = UThMixU235AtomRho/(UThMixOAtomRho + UThMixThAtomRho + UThMixU235AtomRho + UThMixU238AtomRho) # Uranium Atom Fraction in Fuel Mixture [atom/b-cm]
UThMixU238AtomFrac = UThMixU238AtomRho/(UThMixOAtomRho + UThMixThAtomRho + UThMixU235AtomRho + UThMixU238AtomRho) # Uranium Atom Density in Fuel Mixture [atom/b-cm]
# Normalized to 3
UThMixOAtomFrac3 = UThMixOAtomFrac*3 # Oxygen Atom Fraction in Fuel Mixture [atom/b-cm]
UThMixU235AtomFrac3 = UThMixU235AtomFrac*3 # Uranium Atom Fraction in Fuel Mixture [atom/b-cm]
UThMixU238AtomFrac3 = UThMixU238AtomFrac*3 # Uranium Atom Fraction in Fuel Mixture [atom/b-cm]
UThMixThAtomFrac3 = UThMixThAtomFrac*3 # Uranium Atom Density in Fuel Mixture [atom/b-cm]
print("Thorium atom fraction:", UThMixThAtomFrac)
print("Oxygen atom fraction:", UThMixOAtomFrac)
print("U-235 atom fraction:", UThMixU235AtomFrac)
print("U-238 atom fraction:", UThMixU238AtomFrac)

# Volume Parameters
# Cylinder
cyl_Diameter = 2*100 # based on shipping container width of 8*8.6 feet (2 meters will comfortably fit)
cyl_Radius = cyl_Diameter/2
cyl_Length = 10*100 # based on shipping container length of 40 feet (10 meters will comfortably fit)
cylinderVol = math.pi*cyl_Radius**2*cyl_Length #[m^3]
clad_cyl_Radius = 1.1*cyl_Radius
clad_cyl_Length = 1.05*cyl_Length
clad_Thickness = (clad_cyl_Length - cyl_Length)/2.
cladVol = math.pi*clad_cyl_Radius**2*clad_cyl_Length-cylinderVol #[m^3]
clad_rho = 6.56 #[g/cm^3]

print("Volume Fraction of Fuel: " + "{:.2f}".format(fuelVolFrac*100) + "%")
graphVolFrac = vF_Pebbles_Core*(1 - (vF_TRISO_Pebbles*vF_Fuel_TRISO)) # Total Volume Fraction of Graphite
print("Volume Fraction of Graphite: " + "{:.2f}".format(graphVolFrac*100) + "%")
hel_CoolVolFrac = 1 - (fuelVolFrac + graphVolFrac) # Total Volume Fraction of Helium
print("Volume Fraction of Helium: " + "{:.2f}".format(hel_CoolVolFrac*100) + "%")
totalVolFrac = fuelVolFrac + graphVolFrac + hel_CoolVolFrac # Total Volume
print("Total Volume Fraction (should equal 100%): " + "{:.2f}".format(totalVolFrac*100) + "%")

# Volume Output
# Cylinder
fuelCylinderVol = fuelVolFrac*cylinderVol
graphCylinderVol = graphVolFrac*cylinderVol
hel_CoolCylinderVol = hel_CoolVolFrac*cylinderVol
totalCylinderVol = fuelCylinderVol + graphCylinderVol + hel_CoolCylinderVol

# Mass Output
#Cylinder
fuelCylinderMass = UThMixMassRho*fuelCylinderVol
graphCylinderMass = graphMassRho*graphCylinderVol
hel_CoolCylinderMass = he_CoolMassRho*hel_CoolCylinderVol
totalCylinderMass = fuelCylinderMass + graphCylinderMass + hel_CoolCylinderMass


## Define Materials ##
fuel = openmc.Material(name='fuel');
#The below ratios were calculated assuming 95%LEU-5%Th Oxide Fuel.
fuel.add_nuclide('Th232', UThMixThAtomFrac, 'ao')
fuel.add_nuclide('U235', UThMixU235AtomFrac, 'ao') #5% U-235 enrichment
fuel.add_nuclide('U238', UThMixU238AtomFrac, 'ao')
fuel.add_element('O', UThMixOAtomFrac, 'ao')
fuel.set_density('g/cm3', UThMixMassRho) # Based on assumption of fuel density within TRISO
fuel.temperature = temp_homog


# Establish Graphite Moderator material
graph = openmc.Material(name='graph')
# add nuclides to graph
graph.add_element('C', 1.00)
graph.set_density('g/cm3', graphMassRho)

# Establish Helium Coolant material
hel_Cool = openmc.Material(name='hel_Cool')
# add nuclides to hel_Cool
hel_Cool.add_nuclide('He3', 0.000002)
hel_Cool.add_nuclide('He4', 0.999998)
hel_Cool.set_density('g/cm3', he_CoolMassRho)
hel_Cool.temperature = temp_homog


#materials = openmc.Materials([fuel, graph, hel_Cool])
mixMat = openmc.Material.mix_materials([fuel,graph,hel_Cool],
                                   fracs=[fuelVolFrac,graphVolFrac,hel_CoolVolFrac],
                                  percent_type='vo',
                                  name = 'mixMat_Fuel')
mixMat.volume = totalCylinderVol

# Define Cladding Material
# Zircaloy 4
# ref 1: https://www.azom.com/article.aspx?ArticleID=7644
# ref 2: https://www.sciencedirect.com/topics/engineering/zircaloy#:~:text=Zircaloy%2D4%20was%20obtained%20by,%E2%80%93%20O%20(0.003%20wt%25).
cladMat = openmc.Material(name = 'cladMat')
cladMat.add_element('Zr', 0.98197)
cladMat.add_element('Sn', 0.015)
cladMat.add_element('Fe', 0.002)
cladMat.add_element('Cr', 0.001)
cladMat.add_element('O', 0.00003)      
cladMat.set_density('g/cm3', clad_rho) # Based on assumption of fuel density within TRISO
cladMat.volume = cladVol


materials = openmc.Materials()
#materials += [fuel, graph, hel_Cool]
materials += [mixMat, cladMat]

materials.export_to_xml()

## Define Universe Geometry

universeCylinder = openmc.model.RightCircularCylinder([0, 0, -0.25*cyl_Length], 1.5*cyl_Length, 1.5*cyl_Radius, axis='z')

insideCylinder = -universeCylinder
# outsideCylinder = +universeCylinder #unused

cell = openmc.Cell()
cell.region = insideCylinder

universe = openmc.Universe()
universe.add_cell(cell)

## Define Bounding Geometry ##
matCylinder = openmc.model.RightCircularCylinder([0, 0, 0], cyl_Length, cyl_Radius, axis='z')
cladCylinder = openmc.model.RightCircularCylinder([0, 0, -clad_Thickness], clad_cyl_Length, 
                                                  clad_cyl_Radius, axis='z',
                                                  boundary_type='vacuum') #arbitrarily decided cladding width (may have to adjust)
material_region = -matCylinder
clad_region = -cladCylinder & +matCylinder

# fuel
material_Geom = openmc.Cell(name='material_Geom')
material_Geom.fill = mixMat
material_Geom.region = material_region

# cladding
clad_Geom = openmc.Cell(name='clad_Geom')
clad_Geom.fill = cladMat
clad_Geom.region = clad_region


root_universe = openmc.Universe(cells=[material_Geom, clad_Geom])

geometry = openmc.Geometry()
geometry.root_universe = root_universe

geometry.export_to_xml()


## Cross Sections ##
## Source ##
# create a point source
point = openmc.stats.Point((0,0,0))
source = openmc.IndependentSource(space=point)

settings = openmc.Settings()
settings.source = source

n_batches =110
n_inactive = 10
n_particlesPerBatch = 10000

settings.batches = n_batches
settings.inactive = n_inactive
settings.particles = n_particlesPerBatch

settings.export_to_xml()

## Tallies ##

# arrays to hold tally data
efilter_values = [0.0, 10, 14.0e6]
numEbins = len(efilter_values)-1
th232_cap_data_th = np.zeros((1,1))
u238_cap_data_th = np.zeros_like(th232_cap_data_th)
th232_cap_data_fast = np.zeros_like(th232_cap_data_th)
u238_cap_data_fast = np.zeros_like(th232_cap_data_th)

cell_filter = openmc.CellFilter(material_Geom)

efilter = openmc.EnergyFilter(values=efilter_values)

#flux_tally = openmc.Tally(name='flux')
#flux_tally.scores =['flux']
#flux_tally.filters = [cell_filter, efilter]

th232_capture_tally = openmc.Tally(name='Th232_capture')
th232_capture_tally.scores = ['(n,gamma)']
th232_capture_tally.filters = [cell_filter, efilter]

u238_capture_tally = openmc.Tally(name='U238_captuer')
u238_capture_tally.scores = ['(n,gamma)']
u238_capture_tally.filters = [cell_filter, efilter]
                       


# this tally basically scores most interactions with U235
# but does not distinguish them
# maybe not terribly useful
#tally = openmc.Tally()
#tally.filters = [cell_filter]
#tally.nuclides = ['U235']
#tally.scores = ['total','fission','absorption','(n,gamma)']



tallies = openmc.Tallies([th232_capture_tally,
                          u238_capture_tally])

tallies.export_to_xml()

openmc.run(output=False)
sp_filename = 'statepoint.' + str(settings.batches) + '.h5'
sp = openmc.StatePoint(sp_filename);
keffVal = sp.keff
print("{:.2f}".format(vF_TRISO_Pebbles*100), 
      "% Packing Factor,", "{:.2f}".format(pctLEUVal*100), 
      "% Uranium (","{:.2f}".format((1-pctLEUVal)*100), 
      "% Thorium ):", keffVal)
sep = '+/-'

keffValFloat = float(str(keffVal).split(sep, 1)[0])
print(keffValFloat)
reactMatFloat = (keffValFloat-1)/keffValFloat



# get capture tally data
th232_cap = sp.get_tally(name=th232_capture_tally.name)
th232_cap_df = th232_cap.get_pandas_dataframe()
th232_cap_vals = th232_cap_df['mean'].to_numpy()

u238_cap = sp.get_tally(name=u238_capture_tally.name)
u238_cap_df = u238_cap.get_pandas_dataframe()
u238_cap_vals = u238_cap_df['mean'].to_numpy()


#for k in range(numEbins):
#    th232_cap_data[i,j,k] = th232_cap_vals[k]
#    u238_cap_data[i,j,k] = u238_cap_vals[k]

th232_cap_data_th = th232_cap_vals[0];
th232_cap_data_fast = th232_cap_vals[1];

u238_cap_data_th = u238_cap_vals[0];
u238_cap_data_fast = u238_cap_vals[1];


sp.close()

model = openmc.model.Model(geometry,materials,settings)
operator = openmc.deplete.CoupledOperator(model,"chain_endfb71_pwr.xml");

# typical PWR power density
full_pd = 20.0; # W/gHM, estimate
power_density = [full_pd,full_pd,full_pd,full_pd,full_pd,
                full_pd,full_pd,full_pd,full_pd,full_pd,
                0,0,0,0,0]; # power density W/gHM 
# power 0 after 4.5 years with cooldown steps of a day, week, month to 2 years
days = 24*3600;
time_steps = [0.5*days,0.5*days,1*days,5*days,
              23*days,150*days,365*days,365*days,
              365*days,365*days,
              1*days,6*days,23*days,335*days,365*days];
cecm = openmc.deplete.CECMIntegrator(operator,time_steps,power_density=power_density);

repeat_depletion = False


if(repeat_depletion):
    cecm.integrate()
    
    
# get depletion results to manipulate
r = openmc.deplete.Results('depletion_results.h5')
burned_mats = r.export_to_materials(burnup_index=15)
burned_mats.export_to_xml('BurnedMaterials15.xml')

print(burned_mats)

mat_tree = et.parse('BurnedMaterials15.xml')
root = mat_tree.getroot()
i=0
for child in root:
    if child.attrib['name']=='mixMat_Fuel':
        uo2_elem = root[i]
    i+=1
    
# create Material object from element in burned Materials object
uo2_elem.set('id',241824)
print(uo2_elem.items())
type(uo2_elem)
burned_uo2 = openmc.Material.from_xml_element(uo2_elem)
burned_uo2_mass = burned_uo2.get_mass()

#burned_uo2 = openmc.Material(name='burned_uo2')
#Burned_uo2 = burned_uo2.from_xml_element(uo2_elem)
print(burned_uo2)
print(burned_uo2_mass)

listnuc = burned_uo2.get_nuclides() # list of nuclides present in burned fuel


# get string with all Pu isotopes present in burned fuel
# isotopes that will be present after chemical processing
import re
Puiso = []
for nuclide in listnuc:
    if re.search('Pu.+', nuclide):
        Puiso.append(nuclide)
        
pu_mass =0.
for nuclide in Puiso:
    pu_mass+=burned_uo2.get_mass(nuclide=nuclide)
print(pu_mass)

pu_mass_fraction = pu_mass/burned_uo2_mass
print(pu_mass_fraction)

# create metallic Pu from separated Pu product in Burned Fuel
SepPu = openmc.Material(name='PuProduct')
SepPu.set_density('g/cc',19.84) # density used for all metallic Plutonium in PNNL Compendium
print(Puiso)
i = len(Puiso)
n = 0
BurnPuAo = []
while (n < i):
    BurnPu = burned_uo2.get_nuclide_atom_densities(Puiso[n])
    BurnPuAo.append(BurnPu)
    SepPu.add_nuclide(Puiso[n],BurnPu[Puiso[n]])
    n+=1
print(BurnPuAo)
print(SepPu)


def build_model(radius, fuel):
    
    
    materials = openmc.Materials([fuel])
    
    # create sphere with radius parameter
    sphere_radius = openmc.Sphere(x0=0,y0=0,z0=0,r=radius, boundary_type='vacuum', name='BCM')
    
    # create core cell
    core_cell = openmc.Cell(name='Bare Critical Sphere')
    core_cell.fill = fuel
    core_cell.region = -sphere_radius
    
    # create universe geometry
    root_universe = openmc.Universe(name='root universe')
    root_universe.add_cells([core_cell])
    
    geometry = openmc.Geometry(root_universe)
    # define criticality settings
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue' # keff calculation
    settings.particles = 5000 # particles per batch (mo betta)
    settings.batches = 250 # number of batches
    settings.inactive = 50 # inactive batches
    
    settings.output = {'tallies': False}
    
    model = openmc.model.Model(geometry,materials,settings)
    
    return model


crit_r, guesses, keffs = openmc.search_for_keff(build_model, bracket=[1,50],model_args={'fuel':SepPu},
                                                tol=1e-4, print_iterations=True,
                                               run_args={'output':False})
# print results and collect data
print('Burned Plutonium Critical Mass')
print('The bare critical sphere radius is %7.4f cm.' % crit_r)
crit_v = 4/3*math.pi*crit_r**3 # volume of critical sphere (cc)

BCM = crit_v * 19.84 /1000 # mass of critical radius (kg)
print('The bare critical mass is %7.3f kg.' % BCM)

BCMs = np.array(BCM)
print(BCMs,
      '\n')

# get activity from burned fuel
print('Target material activity is %5.3g Bq/g ' % burned_uo2.get_activity())
burnact = burned_uo2.get_activity(units='Bq/g',by_nuclide=True)
print(burnact)


net_weight_HALEU = BCM/pu_mass_fraction
print(net_weight_HALEU,' kg') # in kg only fuel material (no clad)

total_spec_act = sum(burnact.values()) 
totalact = total_spec_act*net_weight_HALEU/(3.7e7) # total activity from nuclear fuel required for one BCM (Ci)
print(totalact,' Ci')