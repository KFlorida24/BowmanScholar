import openmc
import math
import numpy as np
#import pandas as pd
#from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import openmc.deplete
from IPython.display import Image
import xml.etree.ElementTree as et

#nitial General Disclaimer
print("Code written by MIDN 1/C Karl Florida, USN under advisement of CAPT Stuart Blair, USN")
print("Using endfb VII.1 Cross Section Library")
print("This code is unfinished and should not be used outside of test-running\n")

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
avo = 0.6022 # Avogadro's Number

# Volume Fraction Constants
vF_Pebbles_Core = 2./3 # Volume Fraction of Pebbles in Core (Max packing factor)
vF_TRISO_Pebbles = 0.4 # Volume Fraction of TRISO in Pebbles
vF_Fuel_TRISO = 0.15 # Volume Fraction of Fuel in TRISO
fuelVolFrac = vF_Pebbles_Core*vF_TRISO_Pebbles*vF_Fuel_TRISO # Volume Fraction of Fuel
graphVolFrac = vF_Pebbles_Core*(1 - (vF_TRISO_Pebbles*vF_Fuel_TRISO)) # Total Volume Fraction of Graphite
hel_CoolVolFrac = 1 - (fuelVolFrac + graphVolFrac) # Total Volume Fraction of Helium
totalVolFrac = fuelVolFrac + graphVolFrac + hel_CoolVolFrac # Total Volume

def matMixFunInput():

    # Run Mode
    while True:
        try:
            print('Please select whether you want to do one run with user inputs (1), a continuous run (2), or Burn-Up (3).\n')
            userInputTF = int(input() or 2)
        except ValueError:
            print("Please type in a valid number (1 or 2)")
            continue
        if userInputTF not in (1,2,3):
            print("Please type in a valid entry (1 or 2)")
            continue
        else:
            break
    if userInputTF == 1:
        print("Single Run, User Input")
    elif userInputTF == 2:
        print("Continuous Run")
    elif userInputTF == 3:
        print("Burn-Up")
   
    #Geometry Definitions
    while True:
        try:
            print("Please input diameter of core (cylinder) in m. (Press Enter key for default for default value of 2 m)\n")
            cyl_Diameter = float(input()*100 or 2*100) # based on shipping container width of 8*8.6 feet (2 meters will comfortably fit)
            cyl_Radius = cyl_Diameter/2
        except ValueError:
            print("Please type in a valid number")
            continue
        else:
            break
        
    while True:
        try:
            print('Please input length of core (cylinder) in m. (Press Enter key for default for default value of 10 m)\n')
            cyl_Length = float(input()*100 or 10*100) # based on shipping container length of 40 feet (10 meters will comfortably fit)
        except ValueError:
            print("Please type in a valid number")
            continue
        else:
            print("Cylinder with diameter of", cyl_Diameter/100, "m and length", cyl_Length/100, "m generated.")
            break


    if userInputTF == 1:
        while True:
            try:
                UO2MassRho = float(input('Desired UO2 Mass Density [g/cc]? (Press Enter key for default for default value of 10.45 g/cc)\n') or 10.45)
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            else:
                break
        while True:
            try:
                UO2Enrichment = float(input('Desired level of UO2 Enrichment (%)? (Press Enter key for default for default value of 15%)\n') or 15)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            else:
                break
        while True:
            try:
                pctLEU = float(input('Desired TOX Fuel Uranium Fraction (%)? (Press Enter key for default for default value of 40%)\n') or 40)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            else:
                break
        print('UO2 Mass Density: ', UO2MassRho, 'g/cc')
        print('UO2 Enrichment: ', UO2Enrichment*100, '%')
        print('Fuel UO2 Percent: ', pctLEU*100, '%')
        print('Fuel Thorium Percent: ', (1 - pctLEU)*100, '%')
        while True:
            try:
                n_batches = int(input('How many batches per run? (Press Enter key for default for default value of 110)\n') or 110)
            except ValueError:
                print("Please type in a valid integer (no characters)")
                continue
            if n_batches <= 0:
                print("Must be greater than 0")
                continue
            else:
                break
        while True:
            try:
                n_inactive = int(input('How many discarded batches per run? (Press Enter key for default for default value of 10)\n') or 10)
            except ValueError:
                print("Please type in a valid integer (no characters)")
                continue
            if n_batches <= n_inactive:
                print("Must be less than the number of batches run")
                continue
            else:
                break
        while True:
            try:
                n_particlesPerBatch = int(input('How many particles ber batch? (Press Enter key for default for default value of 10000\n') or 10000)
            except ValueError:
                print("Please type in a valid integer (no characters)")
                continue
            if n_particlesPerBatch <= 0:
                print("Must be greater than 0")
                continue
            else:
                break
        print('# of batches: ', n_batches)
        print('# of discarded batches: ', n_inactive)
        print('# of particles per batch: ', n_particlesPerBatch)
        
        ## Calculate Weight Percentages ##
        ThO2MolW = ThMolW + 2*OMolW # Molecular Weight of ThO2 [g/mol]
        
        # Results of Constants and Input Parameters
        # Weight percent of Th in LEU mixture
        pctTh = 1 - pctLEU 
        # Total Mixture Density
        UThMixMassRho = (pctLEU/UO2MassRho + pctTh/ThO2MassRho)**-1 # [g/cm]
        # UO2 Atom Densities
        UO2UMolW = ((UO2Enrichment/U235MolW) + ((1 - UO2Enrichment)/U238MolW))**-1 # Molecular Weight of U in UO2 [g/mol]
        UO2MolW = UO2UMolW + 2*16 # Molecular Weight of UO2 [g/mol]
        UO2OAtomRho = (UO2MassRho/UO2MolW)*2*avo # O Atom Density in UO2 [atom/b-cm]
        UO2U235AtomRho = (UO2Enrichment*UO2MassRho*(UO2UMolW/UO2MolW)/U235MolW)*avo # U-235 Atom Density in UO2 [atom/b-cm] 
        UO2U238AtomRho = ((1 - UO2Enrichment)*UO2MassRho*(UO2UMolW/UO2MolW)/U238MolW)*avo # U-238 Atom Density in UO2 [atom/b-cm]
        UO2U235AtomRhoFrac = UO2U235AtomRho/(UO2U235AtomRho+UO2U238AtomRho)
        UO2U238AtomRhoFrac = UO2U238AtomRho/(UO2U235AtomRho+UO2U238AtomRho)
        # UO2 ThO2 Mixture Atom Densities
        UThMixOAtomRho = (pctTh*UThMixMassRho/ThO2MolW)*2*avo + (pctLEU*UThMixMassRho/UO2MolW)*2*avo # Oxygen Atom Density in Fuel Mixture [atom/b-cm]
        UThMixThAtomRho = (pctTh*UThMixMassRho/ThO2MolW)*avo # Uranium Atom Density in Fuel Mixture [atom/b-cm] 
        UThMixU235AtomRho = (pctLEU*UThMixMassRho/UO2MolW)*avo*UO2U235AtomRhoFrac # Uranium Atom Density in Fuel Mixture [atom/b-cm]
        UThMixU238AtomRho = (pctLEU*UThMixMassRho/UO2MolW)*avo*UO2U238AtomRhoFrac # Uranium Atom Density in Fuel Mixture [atom/b-cm]
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
        
        # Volume Parameters
        # Cylinder
        cylinderVol = math.pi*cyl_Radius**2*cyl_Length #[m^3]
        clad_cyl_Radius = 1.1*cyl_Radius
        clad_cyl_Length = 1.05*cyl_Length
        clad_Thickness = (clad_cyl_Length - cyl_Length)/2.
        cladVol = math.pi*clad_cyl_Radius**2*clad_cyl_Length-cylinderVol #[m^3]
        clad_rho = 6.56 #[g/cm^3]
        
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

        #materials = openmc.Materials([fuel, graph, hel_Cool])
        mixMat = openmc.Material.mix_materials([fuel,graph,hel_Cool],
                                           fracs=[fuelVolFrac,graphVolFrac,hel_CoolVolFrac],
                                          percent_type='vo',
                                          )
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
        # ref: https://openmc.discourse.group/t/cylindrical-geometry-more-than-95-of-external-source-sites-sampled-were-rejected/1432
        source = openmc.Source()
        rad_src = openmc.stats.Uniform(a=0, b=cyl_Radius)
        phi_src = openmc.stats.Uniform(a=0, b=2*math.pi)
        z_src = openmc.stats.Uniform(a=0, b=cyl_Length)
        origin_src = (0.0, 0.0, cyl_Length/2)
        source.space = openmc.stats.CylindricalIndependent(r=rad_src,phi=phi_src,z=z_src,origin = origin_src)
        source.angle = openmc.stats.Isotropic()
        #source.energy = openmc.stats.Discrete([10.0e6], [1.0])
        #source.time = openmc.stats.Uniform(0, 1e-6)
        
        settings = openmc.Settings()
        settings.source = source
        #settings.run_modes = 'fixed source'
        settings.run_mode = 'eigenvalue';
        settings.temperature['method']='interpolation';
        settings.batches = n_batches
        settings.inactive = n_inactive
        settings.particles = n_particlesPerBatch
    
        
        cell_filter = openmc.CellFilter(material_Geom)
        
        tally = openmc.Tally()
        tally.filters = [cell_filter]
    
        tally.nuclides = ['U235']
        tally.scores = ['total','fission','absorption','(n,gamma)']
       
        tallies = openmc.Tallies([tally])
        tallies.export_to_xml()

        openmc.run()
        
        
        # Plot XY Geometry
        plotXY = openmc.Plot()
        plotXY.filename = 'homogeneousMatMixPlotXY'
        plotXY.width = (2*clad_cyl_Radius*1.1,2*clad_cyl_Radius*1.1)
        plotXY.pixels = (600,600)
        plotXY.color_by = 'material'
        plotXY.colors = {mixMat: 'yellow',cladMat: 'blue'}
        plotXY.basis = 'xy'
        plotXY.origin = (0,0,clad_cyl_Length/2)
        
        plotsXY = openmc.Plots([plotXY])
        plotsXY.export_to_xml()

        openmc.plot_geometry()
        # openmc.plot command documentation: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.Plot.html
        # openmc.plot command color documentation: https://www.w3.org/TR/SVG11/types.html#ColorKeywords
        Image('homogeneousMatMixPlotXY.png')
        plotXY.to_ipython_image()
        
        
        # Plo XZ Geometry
        plotXZ = openmc.Plot()
        plotXZ.filename = 'homogeneousMatMixPlotXZ'
        plotXZ.width = (2*clad_cyl_Radius*1.1,clad_cyl_Length*1.1)
        plotXZ.pixels = (600,round(600*(clad_cyl_Length*1.1)/(2*clad_cyl_Radius*1.1)))
        plotXZ.color_by = 'material'
        plotXZ.colors = {mixMat: 'yellow',cladMat: 'blue'}
        plotXZ.basis = 'xz'
        plotXZ.origin = (0,0,clad_cyl_Length/2)
        
        plotsXZ = openmc.Plots([plotXZ])
        plotsXZ.export_to_xml()

        openmc.plot_geometry()
        # openmc.plot command documentation: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.Plot.html
        # openmc.plot command color documentation: https://www.w3.org/TR/SVG11/types.html#ColorKeywords
        Image('homogeneousMatMixPlotXZ.png')
        plotXZ.to_ipython_image()
        


       
    elif userInputTF == 2:
        while True:
            try:
                UO2MassRho = float(input('Desired UO2 Mass Density [g/cc]? (Press Enter key for default for default value of 10.45 g/cc)') or 10.45)
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            else:
                break
        while True:
            try:
                enrichRunNum = int(input('Desired number of enrichments tested? (Press Enter key for default value of 2)') or 2)
            except ValueError:
                print("Please type in a valid number (must be an integer, no characters)")
                continue
            if enrichRunNum < 2:
                print("# of runs must be more than 1")
                continue
            else:
                break
        while True:
            try:
                fuelMixRunNum = int(input('Desired number of fuel mixtures tested? (Press Enter key for default value of 2)') or 2)
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            if fuelMixRunNum < 2:
                print("# of runs must be more than 1")
                continue
            else:
                break
        while True:
            try:
                UO2EnrichmentFirst = float(input('Starting level of UO2 Enrichment (%)? (Press Enter key for default value of 15%)') or 15)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            if UO2EnrichmentFirst < 0:
                print("Value must be greater than 0.")
                continue
            else:
                break
        while True:
            try:
                UO2EnrichmentLast = float(input('Final level of UO2 Enrichment (%)? (Press Enter key for default value of 19.75%)') or 19.75)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            if UO2EnrichmentLast <= UO2EnrichmentFirst:
                print("Final value must be greater than starting value.")
                continue
            else:
                break
        while True:
            try:
                pctLEUFirst = float(input('Starting TOX Fuel Uranium Fraction (%)? (Press Enter key for default for default value of 40%)') or 40)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            if pctLEUFirst < 0:
                print("Value must be greater than 0.")
                continue
            else:
                break
        while True:
            try:
                pctLEULast = float(input('Final TOX Fuel Uranium Fraction (%)? (Press Enter key for default for default value of 95%)') or 95)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            if pctLEULast <= pctLEUFirst:
                print("Final value must be greater than starting value.")
                continue
            else:
                break
        while True:
            try:
                n_batches = int(input('How many batches per run? (Press Enter key for default for default value of 110)') or 110)
            except ValueError:
                print("Please type in a valid integer (no characters)")
                continue
            if n_batches <= 0:
                print("Must be greater than 0")
                continue
            else:
                break
        while True:
            try:
                n_inactive = int(input('How many discarded batches per run? (Press Enter key for default for default value of 10)') or 10)
            except ValueError:
                print("Please type in a valid integer (no characters)")
                continue
            if n_batches <= n_inactive:
                print("Must be less than the number of batches run")
                continue
            else:
                break
        while True:
            try:
                n_particlesPerBatch = int(input('How many particles ber batch? (Press Enter key for default for default value of 10000)') or 10000)
            except ValueError:
                print("Please type in a valid integer (no characters)")
                continue
            if n_particlesPerBatch <= 0:
                print("Must be greater than 0")
                continue
            else:
                break
        ## Calculate Weight Percentages ##
        ThO2MolW = ThMolW + 2*OMolW # Molecular Weight of ThO2 [g/mol]
   
        i = 0
        j = 0
        enrichVals = np.linspace(UO2EnrichmentFirst, UO2EnrichmentLast, num=enrichRunNum, retstep=False)
        pctLEUVals = np.linspace(pctLEUFirst, pctLEULast, num=fuelMixRunNum, retstep=False)
        print("Enrichment Values (%):", enrichVals*100)
        print("Thorium %:", (1-pctLEUVals)*100)

        keffMatInit=np.zeros((enrichRunNum+1,fuelMixRunNum+1))
        keffMat = np.array((keffMatInit), dtype=float)
        keffMat[1:, 0] = enrichVals
        keffMat[0,1:] = pctLEUVals
        print(keffMat)
        reactMatInit=np.zeros((enrichRunNum+1,fuelMixRunNum+1))
        reactMat = np.array((reactMatInit), dtype=float)
        reactMat[1:,0] = enrichVals
        reactMat[0,1:] = pctLEUVals
        
        print('# of batches: ', n_batches)
        print('# of discarded batches: ', n_inactive)
        print('# of particles per batch: ', n_particlesPerBatch)

   
        while i < enrichRunNum:
            while j < fuelMixRunNum:
                # Results of Constants and Input Parameters
                # Weight percent of Th in LEU mixture
                pctTh = 1 - pctLEUVals[j]
                # Total Mixture Density
                UThMixMassRho = (pctLEUVals[j]/UO2MassRho + pctTh/ThO2MassRho)**-1 # [g/cm]
                # UO2 Atom Densities
                UO2UMolW = ((enrichVals[i]/U235MolW) + ((1 - enrichVals[i])/U238MolW))**-1 # Molecular Weight of U in UO2 [g/mol]
                UO2MolW = UO2UMolW + 2*16 # Molecular Weight of UO2 [g/mol]
                UO2OAtomRho = (UO2MassRho/UO2MolW)*2*avo # O Atom Density in UO2 [atom/b-cm]
                UO2U235AtomRho = (enrichVals[i]*UO2MassRho*(UO2UMolW/UO2MolW)/U235MolW)*avo # U-235 Atom Density in UO2 [atom/b-cm]
                UO2U238AtomRho = ((1 - enrichVals[i])*UO2MassRho*(UO2UMolW/UO2MolW)/U238MolW)*avo # U-238 Atom Density in UO2 [atom/b-cm]
                UO2U235AtomRhoFrac = UO2U235AtomRho/(UO2U235AtomRho+UO2U238AtomRho)
                UO2U238AtomRhoFrac = UO2U238AtomRho/(UO2U235AtomRho+UO2U238AtomRho)
                # UO2 ThO2 Mixture Atom Densities
                UThMixOAtomRho = (pctTh*UThMixMassRho/ThO2MolW)*2*avo + (pctLEUVals[j]*UThMixMassRho/UO2MolW)*2*avo # Oxygen Atom Density in Fuel Mixture [atom/b-cm]
                UThMixThAtomRho = (pctTh*UThMixMassRho/ThO2MolW)*avo # Uranium Atom Density in Fuel Mixture [atom/b-cm]
                UThMixU235AtomRho = (pctLEUVals[j]*UThMixMassRho/UO2MolW)*avo*UO2U235AtomRhoFrac # Uranium Atom Density in Fuel Mixture [atom/b-cm]
                UThMixU238AtomRho = (pctLEUVals[j]*UThMixMassRho/UO2MolW)*avo*UO2U238AtomRhoFrac # Uranium Atom Density in Fuel Mixture [atom/b-cm]
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
      
           
                # Parameters
                # Cylinder
                cylinderVol = math.pi*cyl_Radius**2*cyl_Length #[m^3]
                clad_cyl_Radius = 1.1*cyl_Radius
                clad_cyl_Length = 1.05*cyl_Length
                clad_Thickness = (clad_cyl_Length - cyl_Length)/2.
                cladVol = math.pi*clad_cyl_Radius**2*clad_cyl_Length-cylinderVol #[m^3]
                clad_rho = 6.56 #[g/cm^3]
                
           
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

                #materials = openmc.Materials([fuel, graph, hel_Cool])
                mixMat = openmc.Material.mix_materials([fuel,graph,hel_Cool],
                                                   fracs=[fuelVolFrac,graphVolFrac,hel_CoolVolFrac],
                                                  percent_type='vo',
                                                  )
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
                # ref: https://openmc.discourse.group/t/cylindrical-geometry-more-than-95-of-external-source-sites-sampled-were-rejected/1432
                source = openmc.Source()
                rad_src = openmc.stats.Uniform(a=0, b=cyl_Radius)
                phi_src = openmc.stats.Uniform(a=0, b=2*math.pi)
                z_src = openmc.stats.Uniform(a=0, b=cyl_Length)
                origin_src = (0.0, 0.0, cyl_Length/2)
                source.space = openmc.stats.CylindricalIndependent(r=rad_src,phi=phi_src,z=z_src,origin = origin_src)
                source.angle = openmc.stats.Isotropic()
                #source.energy = openmc.stats.Discrete([10.0e6], [1.0])
                #source.time = openmc.stats.Uniform(0, 1e-6)
            
                settings = openmc.Settings()
                settings.source = source
                #settings.run_modes = 'fixed source'
                settings.run_mode = 'eigenvalue'
                settings.temperature['method']='interpolation'
                settings.batches = n_batches
                settings.inactive = n_inactive
                settings.particles = n_particlesPerBatch
                
                
                cell_filter = openmc.CellFilter(material_Geom)
                
                tally = openmc.Tally()
                tally.filters = [cell_filter]
            
                tally.nuclides = ['U235']
                tally.scores = ['total','fission','absorption','(n,gamma)']
               
                tallies = openmc.Tallies([tally])
                tallies.export_to_xml()

                openmc.run(output=False)
               
                sp = openmc.StatePoint('statepoint.100.h5');
                keffVal = sp.keff
                print("{:.2f}".format(enrichVals[i]*100), "% Enrichment,", "{:.2f}".format(pctLEUVals[j]*100), "% Uranium (","{:.2f}".format((1-pctLEUVals[j])*100), "% Thorium ):", keffVal)
                sep = '+/-'
           
                keffValFloat = float(str(keffVal).split(sep, 1)[0])
                reactMatFloat = (keffValFloat-1)/keffValFloat
       
                keffMat[i+1,j+1]= keffValFloat
                reactMat[i+1,j+1] = reactMatFloat
               
                sp.close()
           
                j += 1
            j = 0
            i += 1
        keffMatValsOnly = keffMat[1:,1:]
        print(keffMat)
        print(keffMatValsOnly)
        reactMatValsOnly = reactMat[1:,1:]
        print(reactMat)
        print(reactMatValsOnly)
       
        sns.heatmap(reactMatValsOnly, center=0, cmap = "PiYG", xticklabels = list(map(lambda enrichVals :str(enrichVals) + '%',100*enrichVals.round(2))), yticklabels = list(map(lambda pctLEUVals :str(pctLEUVals) + '%',100*pctLEUVals.round(2))))
        plt.xlabel("Enrichment Values (%)")
        plt.ylabel("Uranium Fraction (%)")
        plt.title("Reactivity as a Function of Enrichment and Uranium Concentration in TOX Fuel")
        plt.savefig('heatmap.png',bbox_inches='tight')
        plt.show()
        
        # Plot XY Geometry
        plotXY = openmc.Plot()
        plotXY.filename = 'homogeneousMatMixPlotXY'
        plotXY.width = (2*clad_cyl_Radius*1.1,2*clad_cyl_Radius*1.1)
        plotXY.pixels = (600,600)
        plotXY.color_by = 'material'
        plotXY.colors = {mixMat: 'yellow',cladMat: 'blue'}
        plotXY.basis = 'xy'
        plotXY.origin = (0,0,clad_cyl_Length/2)
        
        plotsXY = openmc.Plots([plotXY])
        plotsXY.export_to_xml()

        openmc.plot_geometry()
        # openmc.plot command documentation: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.Plot.html
        # openmc.plot command color documentation: https://www.w3.org/TR/SVG11/types.html#ColorKeywords
        Image('homogeneousMatMixPlotXY.png')
        plotXY.to_ipython_image()
        
        
        # Plot XZ Geometry
        plotXZ = openmc.Plot()
        plotXZ.filename = 'homogeneousMatMixPlotXZ'
        plotXZ.width = (2*clad_cyl_Radius*1.1,clad_cyl_Length*1.1)
        plotXZ.pixels = (600,round(600*(clad_cyl_Length*1.1)/(2*clad_cyl_Radius*1.1)))
        plotXZ.color_by = 'material'
        plotXZ.colors = {mixMat: 'yellow',cladMat: 'blue'}
        plotXZ.basis = 'xz'
        plotXZ.origin = (0,0,clad_cyl_Length/2)
        
        plotsXZ = openmc.Plots([plotXZ])
        plotsXZ.export_to_xml()

        openmc.plot_geometry()
        # openmc.plot command documentation: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.Plot.html
        # openmc.plot command color documentation: https://www.w3.org/TR/SVG11/types.html#ColorKeywords
        Image('homogeneousMatMixPlotXZ.png')
        plotXZ.to_ipython_image()
        

        
    elif userInputTF == 3:
        while True:
            try:
                UO2MassRho = float(input('Desired UO2 Mass Density [g/cc]? (Press Enter key for default for default value of 10.45 g/cc)\n') or 10.45)
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            else:
                break
        while True:
            try:
                UO2Enrichment = float(input('Desired level of UO2 Enrichment (%)? (Press Enter key for default for default value of 15%)\n') or 15)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            else:
                break
        while True:
            try:
                pctLEU = float(input('Desired TOX Fuel Uranium Fraction (%)? (Press Enter key for default for default value of 40%)\n') or 40)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            else:
                break
        print('UO2 Mass Density: ', UO2MassRho, 'g/cc')
        print('UO2 Enrichment: ', UO2Enrichment*100, '%')
        print('Fuel UO2 Percent: ', pctLEU*100, '%')
        print('Fuel Thorium Percent: ', (1 - pctLEU)*100, '%')
        while True:
            try:
                n_batches = int(input('How many batches per run? (Press Enter key for default for default value of 110)\n') or 110)
            except ValueError:
                print("Please type in a valid integer (no characters)")
                continue
            if n_batches <= 0:
                print("Must be greater than 0")
                continue
            else:
                break
        while True:
            try:
                n_inactive = int(input('How many discarded batches per run? (Press Enter key for default for default value of 10)\n') or 10)
            except ValueError:
                print("Please type in a valid integer (no characters)")
                continue
            if n_batches <= n_inactive:
                print("Must be less than the number of batches run")
                continue
            else:
                break
        while True:
            try:
                n_particlesPerBatch = int(input('How many particles ber batch? (Press Enter key for default for default value of 10000\n') or 10000)
            except ValueError:
                print("Please type in a valid integer (no characters)")
                continue
            if n_particlesPerBatch <= 0:
                print("Must be greater than 0")
                continue
            else:
                break
        print('# of batches: ', n_batches)
        print('# of discarded batches: ', n_inactive)
        print('# of particles per batch: ', n_particlesPerBatch)
        
        ## Calculate Weight Percentages ##
        ThO2MolW = ThMolW + 2*OMolW # Molecular Weight of ThO2 [g/mol]
        
        # Results of Constants and Input Parameters
        # Weight percent of Th in LEU mixture
        pctTh = 1 - pctLEU 
        # Total Mixture Density
        UThMixMassRho = (pctLEU/UO2MassRho + pctTh/ThO2MassRho)**-1 # [g/cm]
        # UO2 Atom Densities
        UO2UMolW = ((UO2Enrichment/U235MolW) + ((1 - UO2Enrichment)/U238MolW))**-1 # Molecular Weight of U in UO2 [g/mol]
        UO2MolW = UO2UMolW + 2*16 # Molecular Weight of UO2 [g/mol]
        UO2OAtomRho = (UO2MassRho/UO2MolW)*2*avo # O Atom Density in UO2 [atom/b-cm]
        UO2U235AtomRho = (UO2Enrichment*UO2MassRho*(UO2UMolW/UO2MolW)/U235MolW)*avo # U-235 Atom Density in UO2 [atom/b-cm] 
        UO2U238AtomRho = ((1 - UO2Enrichment)*UO2MassRho*(UO2UMolW/UO2MolW)/U238MolW)*avo # U-238 Atom Density in UO2 [atom/b-cm]
        UO2U235AtomRhoFrac = UO2U235AtomRho/(UO2U235AtomRho+UO2U238AtomRho)
        UO2U238AtomRhoFrac = UO2U238AtomRho/(UO2U235AtomRho+UO2U238AtomRho)
        # UO2 ThO2 Mixture Atom Densities
        UThMixOAtomRho = (pctTh*UThMixMassRho/ThO2MolW)*2*avo + (pctLEU*UThMixMassRho/UO2MolW)*2*avo # Oxygen Atom Density in Fuel Mixture [atom/b-cm]
        UThMixThAtomRho = (pctTh*UThMixMassRho/ThO2MolW)*avo # Uranium Atom Density in Fuel Mixture [atom/b-cm] 
        UThMixU235AtomRho = (pctLEU*UThMixMassRho/UO2MolW)*avo*UO2U235AtomRhoFrac # Uranium Atom Density in Fuel Mixture [atom/b-cm]
        UThMixU238AtomRho = (pctLEU*UThMixMassRho/UO2MolW)*avo*UO2U238AtomRhoFrac # Uranium Atom Density in Fuel Mixture [atom/b-cm]
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
        
        # Volume Parameters
        # Cylinder
        cylinderVol = math.pi*cyl_Radius**2*cyl_Length #[m^3]
        clad_cyl_Radius = 1.1*cyl_Radius
        clad_cyl_Length = 1.05*cyl_Length
        clad_Thickness = (clad_cyl_Length - cyl_Length)/2.
        cladVol = math.pi*clad_cyl_Radius**2*clad_cyl_Length-cylinderVol #[m^3]
        clad_rho = 6.56 #[g/cm^3]
        
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

        #materials = openmc.Materials([fuel, graph, hel_Cool])
        mixMat = openmc.Material.mix_materials([fuel,graph,hel_Cool],
                                           fracs=[fuelVolFrac,graphVolFrac,hel_CoolVolFrac],
                                          percent_type='vo',
                                          )
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
        # ref: https://openmc.discourse.group/t/cylindrical-geometry-more-than-95-of-external-source-sites-sampled-were-rejected/1432
        source = openmc.Source()
        rad_src = openmc.stats.Uniform(a=0, b=cyl_Radius)
        phi_src = openmc.stats.Uniform(a=0, b=2*math.pi)
        z_src = openmc.stats.Uniform(a=0, b=cyl_Length)
        origin_src = (0.0, 0.0, cyl_Length/2)
        source.space = openmc.stats.CylindricalIndependent(r=rad_src,phi=phi_src,z=z_src,origin = origin_src)
        source.angle = openmc.stats.Isotropic()
        #source.energy = openmc.stats.Discrete([10.0e6], [1.0])
        #source.time = openmc.stats.Uniform(0, 1e-6)

        settings = openmc.Settings()
        settings.source = source
        #settings.run_modes = 'fixed source'
        settings.run_mode = 'eigenvalue'
        settings.temperature['method']='interpolation'
        settings.batches = n_batches
        settings.inactive = n_inactive
        settings.particles = n_particlesPerBatch
        
        # settings.particles = 10000 # particles per batch
        # settings.batches = 1050 # number of batches
        # settings.inactive = 50 # inactive batches
    
        
        cell_filter = openmc.CellFilter(material_Geom)
        
        tally = openmc.Tally()
        tally.filters = [cell_filter]
    
        tally.nuclides = ['U235']
        tally.scores = ['total','fission','absorption','(n,gamma)']
       
        tallies = openmc.Tallies([tally])
        tallies.export_to_xml()

        openmc.run(output=False)
        sp = openmc.StatePoint('statepoint.100.h5');
        keffVal = sp.keff
        print("K-effective =", keffVal)
        print("Commencing Burn-Up Analysis")
        sp.close()
        
        
        # Plot XY Geometry
        plotXY = openmc.Plot()
        plotXY.filename = 'homogeneousMatMixPlotXY'
        plotXY.width = (2*clad_cyl_Radius*1.1,2*clad_cyl_Radius*1.1)
        plotXY.pixels = (600,600)
        plotXY.color_by = 'material'
        plotXY.colors = {mixMat: 'yellow',cladMat: 'blue'}
        plotXY.basis = 'xy'
        plotXY.origin = (0,0,clad_cyl_Length/2)
        
        plotsXY = openmc.Plots([plotXY])
        plotsXY.export_to_xml()

        openmc.plot_geometry()
        # openmc.plot command documentation: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.Plot.html
        # openmc.plot command color documentation: https://www.w3.org/TR/SVG11/types.html#ColorKeywords
        Image('homogeneousMatMixPlotXY.png')
        plotXY.to_ipython_image()
        
        
        # Plo XZ Geometry
        plotXZ = openmc.Plot()
        plotXZ.filename = 'homogeneousMatMixPlotXZ'
        plotXZ.width = (2*clad_cyl_Radius*1.1,clad_cyl_Length*1.1)
        plotXZ.pixels = (600,round(600*(clad_cyl_Length*1.1)/(2*clad_cyl_Radius*1.1)))
        plotXZ.color_by = 'material'
        plotXZ.colors = {mixMat: 'yellow',cladMat: 'blue'}
        plotXZ.basis = 'xz'
        plotXZ.origin = (0,0,clad_cyl_Length/2)
        
        plotsXZ = openmc.Plots([plotXZ])
        plotsXZ.export_to_xml()

        openmc.plot_geometry()
        # openmc.plot command documentation: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.Plot.html
        # openmc.plot command color documentation: https://www.w3.org/TR/SVG11/types.html#ColorKeywords
        Image('homogeneousMatMixPlotXZ.png')
        plotXZ.to_ipython_image()
        
        
        ## Burn Up Section
        fuel.depletable = True
        
        ## Cross Sections ##
        ## Source ##
        # create a point source
        # ref: https://openmc.discourse.group/t/cylindrical-geometry-more-than-95-of-external-source-sites-sampled-were-rejected/1432
        source = openmc.Source()
        rad_src = openmc.stats.Uniform(a=0, b=cyl_Radius)
        phi_src = openmc.stats.Uniform(a=0, b=2*math.pi)
        z_src = openmc.stats.Uniform(a=0, b=cyl_Length)
        origin_src = (0.0, 0.0, cyl_Length/2)
        source.space = openmc.stats.CylindricalIndependent(r=rad_src,phi=phi_src,z=z_src,origin = origin_src)
        source.angle = openmc.stats.Isotropic()
        #source.energy = openmc.stats.Discrete([10.0e6], [1.0])
        #source.time = openmc.stats.Uniform(0, 1e-6)
    
        settings = openmc.Settings()
        settings.source = source
        #settings.run_modes = 'fixed source'
        settings.run_mode = 'eigenvalue'
        settings.temperature['method']='interpolation'
        settings.batches = n_batches
        settings.inactive = n_inactive
        settings.particles = n_particlesPerBatch
        
        # Depletion
        #Depletion Code is created using ENS Cullinan's code from his GitHub Repo:
            #Link: https://github.com/Toshi-23/Proliferation_Research/blob/main/BCM_Calculations/BurnMats/BU_extraction.ipynb

        model = openmc.model.Model(geometry,materials,settings)
        operator = openmc.deplete.CoupledOperator(model,"chain_endfb71_pwr.xml");
        
        # typical PWR power density (assumption)
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
        
        repeat_depletion = True


        if(repeat_depletion == True):
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
        
        #Start main difference with ENS Cullinan's code
        burn_model_geometry = openmc.Geometry(root_universe)
        burn_model_materials = openmc.Materials([fuel])
        burn_model = openmc.model.Model(burn_model_geometry, burn_model_materials, settings)

        crit_r, guesses, keffs = openmc.search_for_keff(burn_model, bracket=[1,50],model_args={'fuel':SepPu},
                                                        tol=1e-4, print_iterations=True,
                                                       run_args={'output':False})
        #End main difference
        
        # print results and collect data
        print('Burned Plutonium Critical Mass')
        print('The bare critical sphere radius is %7.4f cm.' % crit_r)
        crit_v = 4/3*math.pi*crit_r**3 # volume of critical sphere (cc)

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
        
        print(burnact)

        total_spec_act = sum(burnact.values()) # sum of the specific activity of each isotope (Bq/g)
        print(total_spec_act,' Bq/g')
        print(net_weight_LEU,' kg')


        totalact = total_spec_act*net_weight_LEU/(3.7e7) # total activity from nuclear fuel required for one BCM (Ci)
        print(totalact,' Ci')


        print('BurnedMaterials15.xml')

    

        # Mass of fuel required for a BCM and number of assemblies needed

        print('Assumed total homogenized fuel mass ',totalCylinderMass,'g')


        clad_Mass = (np.pi*((clad_cyl_Radius**2)*cyl_Length - (cyl_Radius**2)*clad_cyl_Length))*cladMat.density;
        print('clad mass of one pin',clad_Mass,'g')

        clad_to_fuel_Ratio = clad_Mass/totalCylinderMass

        Structure_total_mass = clad_Mass + totalCylinderMass
        print('total mass of rhe assembly',Structure_total_mass/1000,'kg')
        print('the mass of Uranium fuel required for one BCM is',net_weight_LEU);

        
        
        
        
