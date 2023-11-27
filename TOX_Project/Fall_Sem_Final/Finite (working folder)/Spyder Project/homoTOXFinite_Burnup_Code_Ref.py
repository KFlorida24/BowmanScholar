import openmc
import openmc.deplete
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
from math import pi
import xml.etree.ElementTree as et

# Create Materials
# fuel uranium dioxode (PNNL)
fuel = openmc.Material(name='uo2')
fuel.set_density('g/cc',10.96)
fuel.add_nuclide('U234', 0.000090, 'ao')
fuel.add_nuclide('U235', 0.010124, 'ao')
fuel.add_nuclide('U236', 0.000046, 'ao')
fuel.add_nuclide('U238', 0.323072, 'ao')
fuel.add_element('O', 0.666667, 'ao')
#fuel.remove_nuclide('O17')
fuel.depletable = True

clad = openmc.Material(name='Zirc4')
clad.set_density('g/cc', 6.56)
clad.add_element('O', 0.006790, 'ao')
clad.add_element('Cr', 0.001741, 'ao')
clad.add_element('Fe', 0.003242, 'ao')
clad.add_element('Zr', 0.977549, 'ao')
clad.add_element('Sn', 0.010677, 'ao')
# clad.add_nuclide('U235', 0.010124e-19, 'ao')
# clad.remove_nuclide('O17')
# clad.depletable = True

water = openmc.Material(name='h2o')
water.set_density('g/cc', 0.712)
water.add_element('H', 2)
water.add_element('O', 1)
#water.remove_nuclide('O17')
water.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([fuel, clad, water])
materials.export_to_xml()

# Create Geometry
h_cell = 300; # height of pincell

r_fuel = 0.42; # fuel radius
r_pin = 0.45; # clad radius

P_D = 1.4; # pitch to diameter ratio
pitch = P_D*(2*r_pin);

fuel_temp = 900; # representaive fuel temperature (K)
mod_temp = 600; # moderator temp (K)

# fuel cylinder
fuel_cyl = openmc.model.RightCircularCylinder([0,0,-h_cell/2],h_cell,r_fuel);

fuel.volume = np.pi*(r_fuel**2)*h_cell;

# pin cylinder
pin_cyl = openmc.model.RightCircularCylinder([0,0,-(h_cell+(r_pin-r_fuel))/2],h_cell+(r_pin-r_fuel)*2,r_pin);

clad.volume = np.pi*((r_pin**2) - (r_fuel**2))*h_cell;

# pin cell container

core_cell = openmc.model.RectangularParallelepiped(-pitch/2,pitch/2,
                                                   -pitch/2,pitch/2,
                                                   -(h_cell+100)/2,(h_cell+100)/2,
                                                   boundary_type = "reflective");
fuel_cell = openmc.Cell();
fuel_cell.region = -fuel_cyl
fuel_cell.fill = fuel;
fuel_cell.temperature = fuel_temp;

clad_cell = openmc.Cell();
clad_cell.region = +fuel_cyl & -pin_cyl;
clad_cell.fill = clad;

mod_cell = openmc.Cell();
mod_cell.region = +pin_cyl & - core_cell;
mod_cell.fill = water;
mod_cell.temperature = mod_temp;

root_univ = openmc.Universe();
root_univ.add_cells([fuel_cell,clad_cell,mod_cell]);

geometry = openmc.Geometry();
geometry.root_universe = root_univ;

geometry.export_to_xml();

# Criticality Settings
settings = openmc.Settings();
settings.run_mode = 'eigenvalue';
settings.particles = 20000;
settings.batches = 250;
settings.inactive = 50;

box = openmc.stats.Box(lower_left = (-r_fuel,-r_fuel,-h_cell/2),
                      upper_right = (r_fuel,r_fuel,h_cell/2),
                      only_fissionable=True);
src = openmc.Source(space=box);

settings.source = src;
settings.temperature['method']='interpolation';
settings.export_to_xml();

# create color dictionary
colors = {}
colors[water]='blue';
colors[fuel]='yellow';
colors[clad]='gray';

# plot pin cell universe to inspect geometry
root_univ.plot(width=(pitch,pitch),color_by='material',colors=colors);

# Depletion

model = openmc.model.Model(geometry,materials,settings)
operator = openmc.deplete.CoupledOperator(model,"chain_endfb71_pwr.xml");

# typical PWR power density
power_density = [30.5,30.5,30.5,30.5,30.5,
                30.5,30.5,30.5,30.5,30.5,
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
burned_mats = r.export_to_materials(burnup_index=15) #consider adding this once updated version maybe  ,path='burnedmats15.xml'
burned_mats.export_to_xml('BurnedMaterials15.xml')

print(burned_mats)

mat_tree = et.parse('BurnedMaterials15.xml')
root = mat_tree.getroot()
i=0
for child in root:
    if child.attrib['name']=='uo2':
        uo2_elem = root[i]
    i+=1
    
 # create Material object from element in burned Materials object
uo2_elem.set('id',23)
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

#Conduct BCM Search

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
crit_v = 4/3*pi*crit_r**3 # volume of critical sphere (cc)

BCM = crit_v * 19.84 /1000 # mass of critical radius (kg)
print('The bare critical mass is %7.3f kg.' % BCM)

BCMs = np.array(BCM)
print(BCMs,
      '\n')


net_weight_LEU = BCM/pu_mass_fraction
print(net_weight_LEU,' kg') # in kg only fuel material (no clad)

# get activity from burned fuel
print('Target material activity is %5.3g Bq/g ' % burned_uo2.get_activity())
burnact = burned_uo2.get_activity(by_nuclide=True,units='Bq/g')
print(burnact)


# plot activities in pie chart
# end of cool down period (2 years)

#newBurnact = burnact.copy()
newBurnact = {}
thresh = 3.7e8
for i, j in burnact.items():
        if j >= thresh:
                newBurnact[i] = j



labels = []
sizes = []
for x, y in newBurnact.items():
    labels.append(x)
    sizes.append(y)
    
print('newburnact',len(newBurnact))
print('burnact',len(burnact)) 
plt.pie(sizes, labels=labels)
plt.axis('equal')
plt.show()
print(newBurnact)

# fuel uranium dioxode (PNNL)
fuel = openmc.Material(name='uo2')
fuel.set_density('g/cc',10.96)
fuel.add_nuclide('U234', 0.000090, 'ao')
fuel.add_nuclide('U235', 0.010124, 'ao')
fuel.add_nuclide('U236', 0.000046, 'ao')
fuel.add_nuclide('U238', 0.323072, 'ao')
fuel.add_element('O', 0.666667, 'ao')

clad = openmc.Material(name='Zirc4')
clad.set_density('g/cc', 6.56)
clad.add_element('O', 0.006790, 'ao')
clad.add_element('Cr', 0.001741, 'ao')
clad.add_element('Fe', 0.003242, 'ao')
clad.add_element('Zr', 0.977549, 'ao')
clad.add_element('Sn', 0.010677, 'ao')

water = openmc.Material(name='h2o')
water.set_density('g/cc', 0.712)
water.add_element('H', 2)
water.add_element('O', 1)
water.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([fuel, clad, water])
materials.export_to_xml()

print(burnact)

total_spec_act = sum(burnact.values()) # sum of the specific activity of each isotope (Bq/g)
print(total_spec_act,' Bq/g')
print(net_weight_LEU,' kg')


totalact = total_spec_act*net_weight_LEU/(3.7e7) # total activity from nuclear fuel required for one BCM (Ci)
print(totalact,' Ci')


print('BurnedMaterials15.xml')

    
# Mass of fuel required for a BCM and number of assemblies needed

pin_fuel_mass = np.pi*(r_fuel**2)*h_cell*fuel.density; # mass in grams
print('fuel mass of one pin',pin_fuel_mass,'g')


pin_clad_mass = (np.pi*((r_pin**2) - (r_fuel**2))*h_cell)*clad.density;
print('clad mass of one pin',pin_clad_mass,'g')

clad_to_fuel = pin_clad_mass/pin_fuel_mass

pin_total_mass = pin_clad_mass + pin_fuel_mass
print('total mass of one pin',pin_total_mass/1000,'kg')
assem_fuel_mass = pin_fuel_mass*264 # based on AP1000 with 264 pins per assembly
print('fuel mass of one assembly',assem_fuel_mass/1000,'kg')
assem_total_mass = pin_total_mass*264 # based on AP1000 with 264 pins per assembly
print('total mass of one assembly',assem_total_mass/1000,'kg')

# total number of assemblies 
num_assem = net_weight_LEU/(assem_fuel_mass/1000); # net weight of fuel required for one BCM (kg) / fuel mass in one assembly (kg)
print('the number of assemblies required to gather the material for one BCM is',num_assem);
total_net_wt_of_assem = np.ceil(num_assem)*(assem_total_mass/1000); # mass of whole number assemblies needed for one BCM (kg)
print('mass of whole number assemblies needed for one BCM',total_net_wt_of_assem,'kg');
