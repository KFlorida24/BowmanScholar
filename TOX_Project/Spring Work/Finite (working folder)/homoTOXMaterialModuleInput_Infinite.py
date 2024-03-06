import openmc
import math
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


def matMixFunInput():

    while True:
        try:
            print('Please select whether you want to do one run with user inputs (1) or a continuous run (2).')
            userInputTF = int(input() or 2)
        except ValueError:
            print("Please type in a valid number (1 or 2)")
            continue
        if userInputTF not in (1,2):
            print("Please type in a valid entry (1 or 2)")
            continue
        else:
            break
    if userInputTF == 1: 
        print("Single Run, User Input")
    elif userInputTF == 2:
        print("Continuous Run")

    if userInputTF == 1:
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
                UO2Enrichment = float(input('Desired level of UO2 Enrichment (%)? (Press Enter key for default for default value of 5%)') or 5)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            else:
                break
        while True:
            try:
                pctLEU = float(input('Desired TOX Fuel Uranium Fraction (%)? (Press Enter key for default for default value of 40%)') or 40)/100
            except ValueError:
                print("Please type in a valid number (no characters)")
                continue
            else:
                break
        print('UO2 Mass Density: ', UO2MassRho, 'g/cc')
        print('UO2 Enrichment: ', UO2Enrichment*100, '%')
        print('Fuel UO2 Percent: ', pctLEU*100, '%')
        print('Fuel Thorium Percent: ', (1 - pctLEU)*100, '%')
        
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
                enrichRunNum = int(input('Desired number packing factors tested? (Press Enter key for default value of 2)') or 2)
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
                UO2EnrichmentFirst = float(input('Starting level of UO2 Enrichment (%)? (Press Enter key for default value of 5%)') or 5)/100
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
                pctLEUFirst = float(input('Starting TOX Fuel Uranium Fraction (%)? (Press Enter key for default for default value of 5%)') or 5)/100
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

        

        ## Calculate Weight Percentages ##
        # Constants - via PNNL Mat'l Comp., Rev-2
        ThO2OAtomRho = 0.044247 # Oxygen Atom Density in ThO2 [atom/b-cm]
        ThO2ThAtomRho = 0.022124 # Thorium Atom Density in ThO2 [atom/b-cm]
        ThMolW = 232. # Molecular Weight of Th [g/mol]
        OMolW = 16. # Molecular Weight of O [g/mol]
        ThO2MolW = ThMolW + 2*OMolW # Molecular Weight of ThO2 [g/mol]
        U235MolW = 235. # Molecular Weight of U-235 [g/mol]
        U238MolW = 238. # Molecular Weight of U-238 [g/mol]
        ThO2MassRho = 9.7 # Mass Density of ThO2 [g/cm^3]
        graphMassRho = 2.26 # Mass Density of Graphite [g/cm^3]
        he_Cool_Z = 1.00652 # He Compression Factor @750C and 6 MPa
        # MW*P/(Z*R*T) from https://www.cedengineering.com/userfiles/H02-011%20-%20Calculation%20of%20Gas%20Density%20and%20Viscosity.pdf
        he_Cool_MW = 4.003 # g/mol
        press_homog = 6000 #kPa
        temp_homog = 750 + 273.15
        R = 8.314472 #L*kPa/K*mol or m^3*Pa/K*mol
        he_CoolMassRho = he_Cool_MW*press_homog/(he_Cool_Z*temp_homog)/1000000 # Mass Density of He @ 293K [g/cm^3]
        # Pressure density data from here https://backend.orbit.dtu.dk/ws/portalfiles/portal/52768509/ris_224.pdf
        avo = 0.6022 # Avogadro's Number
    
        i = 0
        j = 0
        enrichVals = np.linspace(UO2EnrichmentFirst, UO2EnrichmentLast, num=enrichRunNum, retstep=False)
        packFactorVals = np.linspace(0.4, 0.7, num=enrichRunNum, retstep=False)
        pctLEUVals = np.linspace(pctLEUFirst, pctLEULast, num=fuelMixRunNum, retstep=False)
        print("Enrichment Values (%):", enrichVals*100)
        print("Packing Factor Values: ", packFactorVals)
        print("Thorium %:", (1-pctLEUVals)*100)

        keffMatInit=np.zeros((enrichRunNum+1,fuelMixRunNum+1))
        keffMat = np.array((keffMatInit), dtype=float)
        keffMat[1:, 0] = packFactorVals
        keffMat[0,1:] = pctLEUVals
        print(keffMat)
        reactMatInit=np.zeros((enrichRunNum+1,fuelMixRunNum+1))
        reactMat = np.array((reactMatInit), dtype=float)
        reactMat[1:,0] = packFactorVals
        reactMat[0,1:] = pctLEUVals

    
        while i < enrichRunNum:
            while j < fuelMixRunNum:
                # Results of Constants and Input Parameters
                # Weight percent of Th in LEU mixture
                pctTh = 1 - pctLEUVals[j]
                # Total Mixture Density
                UThMixMassRho = (pctLEUVals[j]/UO2MassRho + pctTh/ThO2MassRho)**-1 # [g/cm]
                # UO2 Atom Densities
                UO2UMolW = ((.1975/U235MolW) + ((1 - .1975)/U238MolW))**-1 # Molecular Weight of U in UO2 [g/mol]
                UO2MolW = UO2UMolW + 2*16 # Molecular Weight of UO2 [g/mol]
                UO2OAtomRho = (UO2MassRho/UO2MolW)*2*avo # O Atom Density in UO2 [atom/b-cm]
                UO2U235AtomRho = (.1975*UO2MassRho*(UO2UMolW/UO2MolW)/U235MolW)*avo # U-235 Atom Density in UO2 [atom/b-cm] 
                UO2U238AtomRho = ((1 - .1975)*UO2MassRho*(UO2UMolW/UO2MolW)/U238MolW)*avo # U-238 Atom Density in UO2 [atom/b-cm]
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
        
        
                ## Calculate Vol % ##
                # Constants
                vF_Pebbles_Core = 0.64 # Volume Fraction of Pebbles in Core (Max packing factor)
                vF_TRISO_Pebbles = packFactorVals[i]# Volume Fraction of TRISO in Pebbles
                # vF_TRISO_Pebbles = 0.5 # Volume Fraction of TRISO in Pebbles
                vF_Fuel_TRISO = 0.15 # Volume Fraction of Fuel in TRISO
                fuelVolFrac = vF_Pebbles_Core*vF_TRISO_Pebbles*vF_Fuel_TRISO # Volume Fraction of Fuel
                graphVolFrac = vF_Pebbles_Core*(1 - (vF_TRISO_Pebbles*vF_Fuel_TRISO)) # Total Volume Fraction of Graphite
                hel_CoolVolFrac = 1 - (fuelVolFrac + graphVolFrac) # Total Volume Fraction of Helium
                totalVolFrac = fuelVolFrac + graphVolFrac + hel_CoolVolFrac # Total Volume
            
                # Parameters
                # Cube
                matCubeL = 1 # [cm]
                cubeVol = matCubeL**3 # [cm^3]
                # Sphere
                matSphereRad = 1 # [cm]
                sphereVol = (4/3)*math.pi*matSphereRad**3
            
                # Volume Output
                #Cube
                fuelCubeVol = fuelVolFrac*cubeVol
                graphCubeVol = graphVolFrac*cubeVol
                hel_CoolCubeVol = hel_CoolVolFrac*cubeVol
                totalCubeVol = fuelCubeVol + graphCubeVol + hel_CoolCubeVol
                #Sphere
                fuelSphereVol = fuelVolFrac*sphereVol
                graphSphereVol = graphVolFrac*sphereVol
                hel_CoolSphereVol = hel_CoolVolFrac*sphereVol
                totalSphereVol = fuelSphereVol + graphSphereVol + hel_CoolSphereVol
                # Mass Output
                #Cube
                fuelCubeMass = UThMixMassRho*graphCubeVol
                graphCubeMass = graphMassRho*graphCubeVol
                hel_CoolCubeMass = he_CoolMassRho*graphCubeVol
                totalCubeMass = fuelCubeMass + graphCubeMass + hel_CoolCubeMass
                #Sphere
                fuelSphereMass = UThMixMassRho*graphSphereVol
                graphSphereMass = graphMassRho*graphSphereVol
                hel_CoolSphereMass = he_CoolMassRho*graphSphereVol
                totalSphereMass = fuelSphereMass + graphSphereMass + hel_CoolSphereMass
        
            
                ## Define Materials ##
                fuel = openmc.Material(name='fuel');
                #The below ratios were calculated assuming 95%LEU-5%Th Oxide Fuel.
                fuel.add_nuclide('Th232', UThMixThAtomFrac, 'ao')
                fuel.add_nuclide('U235', UThMixU235AtomFrac, 'ao') #5% U-235 enrichment
                fuel.add_nuclide('U238', UThMixU238AtomFrac, 'ao')
                fuel.add_element('O', UThMixOAtomFrac, 'ao')
                fuel.set_density('g/cm3', UThMixMassRho) # Based on assumption of fuel density within TRISO
                fuel.volume = cubeVol
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
                                                   fracs=[fuelVolFrac,graphVolFrac,hel_CoolCubeVol],
                                                  percent_type='vo',
                                                  )
            
                materials = openmc.Materials()
                #materials += [fuel, graph, hel_Cool]
                materials += [mixMat]
            
                materials.export_to_xml()
        
        
                ## Define Universe Geometry ##
                l_cube = 2.0;
                universeCube = openmc.model.RectangularParallelepiped(-l_cube, l_cube, -l_cube, l_cube, -l_cube, l_cube)
            
                insideCube = -universeCube
                outsideCube = +universeCube
            
                cell = openmc.Cell()
                cell.region = insideCube
            
                universe = openmc.Universe()
                universe.add_cell(cell)
            

            
            
                ## Define Bounding Geometry ##
                matBox = openmc.model.RectangularParallelepiped(-matCubeL/2, matCubeL/2, 
                                                            -matCubeL/2, matCubeL/2, 
                                                            -matCubeL/2, matCubeL/2,
                                                            boundary_type='reflective')
            
                material_region = -matBox
            
                material_Geom = openmc.Cell(name='material_Geom')
                #material_Geom.fill = materials #<<(SRB) This probably isn't correct.  Needs to be a material.
                material_Geom.fill = mixMat
                material_Geom.region = material_region
            
                root_universe = openmc.Universe(cells=[material_Geom])
            
                geometry = openmc.Geometry()
                geometry.root_universe = root_universe
                geometry.export_to_xml()
            
            
                ## Cross Sections ##
            
            
                ## Source ##
                # create a point source
                point = openmc.stats.Point((0,0,0))
                source = openmc.Source(space=point)
            
                settings = openmc.Settings()
                settings.source = source
                settings.batches = 100
                settings.inactive = 10
                settings.particles = 1000
            
                settings.export_to_xml()

                ## Tallies ##
            
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
                print("{:.2f}".format(packFactorVals[i]*100), "% Enrichment,", "{:.2f}".format(pctLEUVals[j]*100), "% Uranium (","{:.2f}".format((1-pctLEUVals[j])*100), "% Thorium ):", keffVal)
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
        reactMatTrans = np.transpose(reactMatValsOnly)
        print(reactMat)
        print(reactMatValsOnly)


        heatmap_Plot = sns.heatmap(reactMatTrans, center=0, cmap = "PiYG", xticklabels = list(map(lambda pctLEUVals :str(pctLEUVals) + '%',100*pctLEUVals.round(2))), yticklabels = list(map(lambda packFactorVals :str(packFactorVals) + '%',100*packFactorVals.round(2))))
        for ind, label in enumerate(heatmap_Plot.get_xticklabels()):
            if ind % 5 == 0:  # every 5th label is kept
                label.set_visible(True)
            else:
                label.set_visible(False)
        for ind, label in enumerate(heatmap_Plot.get_yticklabels()):
            if ind % 5 == 0:  # every 5th label is kept
                label.set_visible(True)
            else:
                label.set_visible(False)
        plt.ylabel("TRISO Packing Factor Values (%)")
        plt.xlabel("Uranium Fraction (%)")
        plt.title("Reactivity as a Function of TRISO Packing Factor and Uranium Concentration in TOX Fuel")
        plt.savefig('heatmap.png',bbox_inches='tight')
        plt.show()