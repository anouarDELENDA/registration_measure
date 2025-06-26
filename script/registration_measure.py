#!/usr/bin/env python
# coding: utf-8

################################################################################################################################################################
#Author: Gabriel Francas (gabriel.francas@postgrad.curtin.edu.au)                                                                                              #
#Last modified: 23/06/2025                                                                                                                                     #
#                              
#   
# Update Note 4: Since previous ring finding algorithm polypy didn't work well with PBC's and non-orthogonal cells, we have found a new ring finding algorithm #
# to employ. (https://github.com/MorrowChem/RingsStatisticsMatter.jl)
#                                                                                                                                                              #
# Update Note 3: Added functionality to output mean registration measure for each atom in .xyz output by taking mean registration measure of each bond it is   #
# part of. This does not give information on the type of stacking, which remains only part of the bonds output in the .dump file.                              #                                                                                             #
#                                                                                                                                                              #
# Update Note 2: We decided that given it is necessary to assess each atom relative to its neighbours to determine the type of stacking it is in (AA or AB     #
# etc.), that it would be more sensible to define the registration measure for each set of neighbouring sp2 atoms, ie. for each bond. Thus this program now    #
# outputs a .dump file which has the format of a LAMMPS dump/local file, in that it can store arbitrary information about each bond. This can be read in by    #
# OVITO using the 'load trajectory' modifier after the atoms have already been loaded in (ie. the original .lmp or .xyz file).                                 #
#                                                                                                                                                              #
#Update Note 1: I have realised that the code made a grave oversiplification/assumption about using distance from nearest atom to calculate distance from 'ideal #
#   stacking position' (above either an atom or ring centre). To do this soundly we really need the ring centres mapped which we can do separately using       #
#   the ring finding algorithm of polypy and then import those results here. The ring centres will be added to the atoms data set as H atoms in order to       #
#   differentiate them. Then, our algorithm can find the nearest atom for the projection, which will be either a carbon or hydrogen and then that is           #
#   automatically the offset and by which atom type it is, determines if the atom was closer to aligning with a ring centre or carbon atom (beta or alpha type)#
#                                                                                                                                                              #    
#Purpose: This program is designed to analyse the stacking registration for a largely sp2 carbon based atomistic model. The system should be largerly          #
#   composed of layered graphenic sheets, which may be curved. This program will calculate the local registration for each sp2 atom in terms of the layers     #
#   directly above and below it as well as its three neighbours. A degree of registration is calculated ranging from 0 (completely unregistered,               #
#   ie. turbostratic) to 1 being perfect registration. The type of stacking will also be determined as either being more closely in line with AAA (hexagonal), #
#   ABB (mixed), ABA (bernal) or ABC (rhombohedral). The type of stacking will change the final output number 'registration-measure' such that in the range    #
#   [0-1] = AAA, (1-2] = ABB, (2-3] = ABA, (3-4] = ABC. Thus the 'registration-measure' value can be used to colour code the atoms for example in OVITO        #
#   using an appropriate colour bar (eg. as attached). The modulus of the 'registration-measure' by 1 will return the original 0 --> 1 value describing        #
#   purely the degree of registration without regard to the type of stacking. In the case where the stacking cannot be measured because the atom is not        #
#   sp2, it is an edge or it is missing an adjacent layer above or below it will be assigned a 'registration-measure' of -1.                                   #                     
################################################################################################################################################################




# This is a new upgraded version 28/11/24 which will have each atom examine the sheet above *and* below in order to allow the quantification of ABC stacking.
# Considering just the middle layer relative to the layers above and below, I have characterised 4 'perfect' stacking possibilities for 3 layers: 
#BBB, ABB, ABA, ABC so I will colour each atom in a discrete range based on which of these categories it most closely falls into. 
#The mismatch for each atom is calculated based on how misaligned it is with the nearest atom in the layers above and below and should vary from 0 
#to half an atomic bond length (0.71A). Average 'mismatch' in atomic alignment between atoms in adjacent layers will be scaled to a value 
#between 0 (maximally misaligned) and 1 (perfectly aligned). 



#Import the required packages
import numpy as np
import copy
import argparse
import time
from multiprocessing import Pool


import ase as ase # We can make use of the package ASE to do alot of functionality with atoms.
from ase import Atom
from ase.visualize import view
from ase.io.lammpsdata import read_lammps_data
from ase.io.lammpsdata import write_lammps_data
from ase.io.extxyz import write_extxyz
from ase.neighborlist import NewPrimitiveNeighborList
from ase.neighborlist import NeighborList
from ase.geometry import get_duplicate_atoms
from ase.geometry.analysis import Analysis


def add_centres(atomsC):
    cutoff = 1.85 #define standard cutoff for carbon-carbon bond
    rs, rings = ring_statistics(atomsC, cutoff=cutoff)
    rings_array = np.array(rings, dtype=object)
    gCom = atomsC.get_center_of_mass(scaled = False)

    #we will only count rings of size 5 - 8 (ring_array elements 4 --> 7)
    for i in range(4):
        for j in range(len(rings_array[i+4])):
            com = atomsC.get_center_of_mass(scaled=False, indices=rings_array[i+4][j])
            centre = Atom('X', com) #tentatively add centre of mass
            atomsC.append(centre)
            cIdx = len(atomsC)-1
            for a in rings_array[i+4][j]:
                if atomsC.get_distance(a,cIdx,mic=True,vector=False)>cutoff:
                    # print('found an issue with com pbcs')
                    initPos = atomsC[a].position
                    shift = np.array(gCom) - np.array(initPos)
                    atomsCopy = copy.deepcopy(atoms)
                    atomsCopy.translate(shift)
                    atomsCopy.wrap()

                    com2 = atomsCopy.get_center_of_mass(scaled=False, indices=rings_array[i+4][j])
                    atomsC.pop()
                    com2 = com2 -1*shift
                    c = Atom('X', com2)
                    atomsC.append(c) #add new proper centre
                    break
    return atomsC

def find_neighbours(nl,index):
    indices, off = nl.get_neighbors(index)
    return indices

def find_normal(atoms, id0, neighbourList):
    #must find plane of the four neighbouring atoms using cross product
    id1 = neighbourList[0]
    id2 = neighbourList[1]
    id3 = neighbourList[2]

    #These are vectors going from the central atom to two of its three neighbours
    vec1 = atoms.get_distance(id0,id1,mic=True,vector=True)
    vec2 = atoms.get_distance(id0,id2,mic=True,vector=True)
    
    cross1 = np.cross(vec1, vec2)
    
    #repeat with other pairs of neighbours
    vec1 = atoms.get_distance(id0,id2,mic=True,vector=True)
    vec2 = atoms.get_distance(id0,id3,mic=True,vector=True)
    cross2 = np.cross(vec1, vec2)
    
    vec1 = atoms.get_distance(id0,id3,mic=True,vector=True)
    vec2 = atoms.get_distance(id0,id1,mic=True,vector=True)
    cross3 = np.cross(vec1, vec2)
    
    arr = np.zeros([3,3])
    arr[:,0]=cross1
    arr[:,1]=cross2
    arr[:,2]=cross3
    
    #mean cross product defines the local plane of the atom
    meanCross = np.mean(arr, axis=1) 
    
    return meanCross


def make_norm_pos(meanCross):
    #This ensures we are looking in the same direction - if the z is less than 0 we reverse the whole vector
    if meanCross[2]<0:
        meanCross = -1*meanCross

    meanCross = meanCross/np.sqrt(np.dot(meanCross,meanCross)) #reduce mean cross product to a unit vector
    return meanCross

def make_norm_neg(meanCross):
    #This ensures we are looking in the same direction - if the z is greater than 0 we reverse the whole vector
    if meanCross[2]>0:
        meanCross = -1*meanCross

    meanCross = meanCross/np.sqrt(np.dot(meanCross,meanCross)) #reduce mean cross product to a unit vector
    return meanCross


def calc_distance(atoms,pointA, normalA, pointB):
    
    c = atoms.get_distance(pointA,pointB,mic=True,vector=True)
    dist = np.dot(c,normalA)
    return dist


def find_plane(atoms, index, normal, meshC, nlM):
    dist = 3.4 #this is a parameter we can tweak depending roughly on the interlayer spacing, probably good to work on this at some point to be more robust
    pos_adjacent = atoms[index].position + dist*normal #we step away from the atom of interest in the direction of the normal a distance of ~1 interlayer dist
    
    #to make things easier with distances in the pbcs, we just add a proxy atom in the position that the normal lands
    a = Atom('P', pos_adjacent)
    meshC.append(a)

    #we then find the single nearest neighbour of this proxy atom (either an atom or ring centre)
    indices = range(len(meshC)-1)
    proxy_idx = len(meshC)-1
    dist = meshC.get_distances(proxy_idx, indices, mic=True)
    dist_min=min(dist)
    neighbour = int(np.where(dist == dist_min)[0])
    neighbours = find_neighbours(nlM,neighbour)
    n_neigh = len(neighbours)
    

    # if neighbour == index: #maybe this??? or (n_neigh != 3 and n_neigh !=6)
    #     #There is no sp2 atom or hexagonal ring centre nearby in the adjacent layer
    #     neighbour = -1 #assign a non-sensible value for neighbour index
    #     interlayer_dist = -1 #assign non-sensible value for distance to later filter out

    if dist_min > 1.85: #there is no actual neighbour vertex present in adjacent sheet
        neighbour = -1 #assign a non-sensible value for neighbour index
        interlayer_dist = -1 #assign non-sensible value for distance to later filter out
    else:
        interlayer_dist = calc_distance(meshC,index,normal,neighbour)

    
    meshC.pop() #remove the proxy atom
    return neighbour, interlayer_dist
    


def find_offset(atoms, atom1, pos2, interlayer_dist, normal, meshC ):
    pos_adjacent = atoms[atom1].position + interlayer_dist*normal
    
    #to make things easier with distances in the pbcs, we just add a proxy atom in the position that the normal lands
    a = Atom('P', pos_adjacent)
    meshC.append(a)
    
    #we then calculate the distance between this proxy atom (atom1 projected into the neighbouring sheet) 
    #and the nearest stacked neighbour (atom2) - we call this the offset
    proxy_idx = len(meshC)-1
    offset_dist = meshC.get_distance(pos2,proxy_idx,mic=True)  
    meshC.pop()
    return offset_dist


def distance_calc(i, atoms, mesh, nl, nlM):
    neighbours = find_neighbours(nl,i)
    meshC = copy.deepcopy(mesh)
    problem = False
    if len(neighbours)==3:
        
        normal = find_normal(atoms,i,neighbours)
        normal_up = make_norm_pos(normal)
        normal_down = make_norm_neg(normal)
        
        layer2U, interlayer_distU = find_plane(atoms,i,normal_up, meshC, nlM)
        layer2D, interlayer_distD = find_plane(atoms,i,normal_down, meshC, nlM)
        # print(layer2atomD, interlayer_distD, layer2atomU, interlayer_distU)
        if layer2U >= 0 and np.abs(interlayer_distU) >0.1 and layer2D >= 0 and np.abs(interlayer_distD) > 0.1:
            o_minU = find_offset(atoms,i, layer2U,interlayer_distU, normal_up, meshC)
            o_minD = find_offset(atoms,i, layer2D,interlayer_distD, normal_down, meshC) 


            layer2UType = meshC[layer2U].symbol
            # print(layer2UType, offset_distU)
            layer2DType = meshC[layer2D].symbol
            posU = 0 if layer2UType == 'C' else 1
            posD = 0 if layer2DType == 'C' else 1


            distance = np.mean(np.array([np.abs(interlayer_distU),np.abs(interlayer_distD)]))
        else:
            #there was a problem with the adjacent layer(s)
            problem = True
            # print(f"There was an issue with adjacent layer(s) such that atom {i} can not be defined as within bulk graphite material")
        
    else:
        #non sp2 atom
        # print(f"Atom {i} not 3-fold coordinated")
        problem = True
    if problem:
        #set all output as nonsense -1
        distance = -1
        o_minU = -1
        o_minD = -1
        posU = -1
        posD = -1

    return {
        "distance": distance,
        "offset_distanceU": o_minU,
        "offset_distanceD": o_minD,
        "posU": posU,
        "posD": posD
    }

def find_distances(atoms,mesh,nl,nlM,num_processes):
    indices = range(len(atoms))
    
    args = [(i,atoms,mesh,nl,nlM) for i in indices]

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(distance_calc,args)

    #gather results
    distances = [res["distance"] for res in results]
    offset_distancesU = [res["offset_distanceU"] for res in results]
    offset_distancesD = [res["offset_distanceD"] for res in results]
    posUs = [res["posU"] for res in results]
    posDs = [res["posD"] for res in results]
    return distances, offset_distancesU, offset_distancesD, posUs, posDs


# In[9]:
def mean_bondL(atoms, index, nl):
    bonds = []
    neighbours = find_neighbours(nl,index)
    for n in neighbours:
        bonds.append(atoms.get_distance(n,index,mic=True))
    if len(bonds) > 0:
        meanB = sum(bonds) / len(bonds)
    else:
        meanB = 1.42
    return meanB

def calc_interlayer_mean(atoms):
    
    boxL = atoms.get_cell()[0][0]
    count=0
    distances = []
    for i in range(len(atoms)):
        #Remove select region of atoms based on narrow work this was originally applied to

        # distX1 = atoms[i].position[0]-boxL/2
        # distY1 = atoms[i].position[1]-np.sqrt(3)/6*boxL

        # distX2 = atoms[i].position[0]-boxL
        # distY2 = atoms[i].position[1]-np.sqrt(3)/3*boxL
        # if np.sqrt(distX1**2+distY1**2)>2 and np.sqrt(distX2**2+distY2**2)>2:
        d = atoms.arrays['interlayer_dist'][i]
        if d>0:
            distances.append(d)
            count+=1
    mean = np.sum(np.array(distances))/count
    return mean, distances


# In[10]:


def calc_2nd_neighbour_mean(atoms,nl):
    boxL = atoms.get_cell()[0][0]
    count=0
    distances = []
    for i in range(len(atoms)):
        first_neighbours = find_neighbours(nl,i)
        if len(first_neighbours) == 3:
            for n in first_neighbours:
                second_neighbours = find_neighbours(nl,n)
                if len(second_neighbours) == 3:
                    second_neighbours = second_neighbours[second_neighbours != i]
                    dist1 = atoms.get_distance(i,second_neighbours[0],mic=True)
                    dist2 = atoms.get_distance(i,second_neighbours[1],mic=True)
                    distances.append(dist1)
                    distances.append(dist2)
                    count+=2
    mean = np.sum(np.array(distances))/count
    return mean, distances


def find_stacking_type(bondRegistration, posUs, posDs, bonds):
    for i in range(len(bonds)):
        a1 = bonds[i,0]
        a2 = bonds[i,1]
        if bondRegistration[i] < 0.001 and bondRegistration[i]>=0: #Hard code perfectly turbostratic to actually be slightly above 0 so when we sort into correct 
            #category it doesn't end up colouring as if it was perfectly stacked
            bondRegistration[i] = 0.001
            # print("found zero") 
        if bondRegistration[i] > 0:
            if (posUs[a1] + posUs[a2]) == 1:
                if (posDs[a1] + posDs[a2]) == 0:
                    #ABB stacking
                    bondRegistration[i] += 1
                else:
                    #now we are either in ABA or ABC, must check further
                    if posDs[a1] == posUs[a1]:
                        #ABA stacking
                        bondRegistration[i] += 2
                    else:
                        #ABC stacking
                        bondRegistration[i] += 3
            elif (posUs[a1] + posUs[a2]) == 0:
                if (posDs[a1] + posDs[a2]) == 1:
                    #ABB stacking
                    bondRegistration[i] += 1
                #else AAA stacking so do not adjust mismatch (varies from - to 1)
                
    return bondRegistration


def find_bond_registration(offset_distancesU, offset_distancesD, bonds, bondRegistration):
    for i in range(len(bonds)):
        a1 = bonds[i,0]
        a2 = bonds[i,1]
        if offset_distancesU[a1] != -1 and offset_distancesD[a1] != -1 and offset_distancesU[a2] != -1 and offset_distancesD[a2] != -1:    
            a1Tot = offset_distancesU[a1]+offset_distancesD[a1]
            a2Tot = offset_distancesU[a2]+offset_distancesD[a2]
            bondRegistration[i] = np.mean([a1Tot, a2Tot])

    max_registration = max(bondRegistration)
    # max_upper = max([max_registration, 1.6396]) #For a bond length of 1.42 A the max distance from an atom OR ring centre within a layer is 0.8198A, since we add up both up and down directions, the max total is twice this.
    max_upper = 1.6396
    print(f"Max upper offset before scaling = {max_upper}")

    #and rescale to make it go from 0 to 1 (least to most perfectly registered)
    for i in range(len(bonds)):
        #only compute a value if we have valid offset
        if bondRegistration[i] >= 0:
            bondRegistration[i] = -1*(bondRegistration[i]/max_upper-1)
    return bondRegistration
# In[13]:



def write_registration(bonds, bondRegistration, outFile):
    head = f"ITEM: NUMBER OF ENTRIES\n{len(bonds)}\nITEM: ENTRIES index atom1 atom2 registration_measure"

    indexes = np.arange(1,len(bonds)+1,1)
    data = np.column_stack((indexes,bonds[:,0]+1,bonds[:,1]+1, bondRegistration))
    np.savetxt(outFile, data,fmt=("%d", "%d", "%d", "%.4f"), header=head, comments="")



def find_atom_registrations(atoms, nl, bonds, bondRegistration):
    registrations = []
    for i in range(len(atoms)):
        first_neighbours = find_neighbours(nl,i)
        if len(first_neighbours) == 3:
            atomPos = np.where(bonds[:,0] == i)[0]
            atomPos = np.append(atomPos,(np.where(bonds[:,1] == i)[0]))
            bondRegistrationTot = 0
            if len(atomPos) != 3:
                print("Error in finding neighbours in bond list.")
            else:
                skip = False
                for j in range(3):
                    if bondRegistration[atomPos[j]] >=0:
                        bondRegistrationTot += bondRegistration[atomPos[j]]%1
                    else:
                        atomReg = -1
                        skip = True
                if not skip:
                    atomReg = bondRegistrationTot/3
        else:
            atomReg = -1
        registrations.append(atomReg)
    return registrations



start = time.time()
parser = argparse.ArgumentParser(description="This script analyses the stacking order within a graphitic system")

#Add file argument
parser.add_argument("--input", "-i", type=str, help="Path to input file to analyse")

parser.add_argument("--output", "-o", type=str, default="<input>_registration-measure", help="Default =<input>_registration-measure. Option to change the output file names. Both files will have the same name, differentiated by extensions of .dump and .xyz.")
parser.add_argument("--processors", "-p", type=int, default=1, help="Default=1. How many cpu's to use for parallel loop")
# parser.add_argument("--detail", "-d", type=bool, default=False, help="Default=False. Write extra .xyz file with more detail including neighbourhood bondlengths per atom and normal vector?")
parser.add_argument("--detail", "-d", action="store_true", help="Default=False. Write extra .xyz file with more detail including neighbourhood bondlengths per atom and normal vector?")

args = parser.parse_args()

inFile = args.input


mean_as=[]
mean_ds = []
js=[]

procs = args.processors
import os
os.environ["JULIA_NUM_THREADS"] = f"{procs}"
import julia
from julia_rings.rings import ring_statistics


print(f"Using {procs} processors...")
detail = args.detail


#############################################################################################################################################################################################
#Read in the data and set up the atomic and mesh systems
atoms=read_lammps_data(inFile, read_image_flags=False)
atoms.set_pbc(True)
atoms.center()
stM = time.time() #mesh start

mesh = copy.deepcopy(atoms)
mesh = add_centres(mesh) #mesh contains the set of ring centre vertices as 'atoms' of type X as well as all C atoms
# write_extxyz(inFile.replace('.lmp','_testWrap.xyz'), atoms, columns=['symbols','positions'])
# write_extxyz(inFile.replace('.lmp','_ringCentres.xyz'), mesh, columns=['symbols','positions'])

etM = time.time() #mesh end

print(f"Time to add ring centres: {(etM-stM):.2f}")


atom_count = len(atoms)
cutoffs = (1.85/2-0.3)*np.ones(len(atoms)) #-0.3 to account for the random 'skin' factor in the ase source code...it just works like this idk
nl = NewPrimitiveNeighborList(cutoffs, self_interaction=False, bothways=True)
nl.update(atoms.pbc, atoms.get_cell(), atoms.positions)

cutoffsM = (1.85/2-0.3)*np.ones(len(mesh)) #-0.3 to account for the random 'skin' factor in the ase source code...it just works like this idk
nlM = NewPrimitiveNeighborList(cutoffsM, self_interaction=False, bothways=True)
nlM.update(mesh.pbc, mesh.get_cell(), mesh.positions) #neighbour list for mesh vertices

################################################################################################################################################################################################################
#Actually run the analysis:

#Find offset distances and alpha/beta status for each atom in up and down directions:
distances, offset_distancesU, offset_distancesD, posUs, posDs = find_distances(atoms, mesh, nl, nlM, procs)


#Create list of bonds (atom pairs)
nl2 = NeighborList(cutoffs, self_interaction=False, bothways=True)
nl2.update(atoms)
atomAnalysis = Analysis(atoms,nl2)
atomAnalysis.adjacency_matrix
bonds = np.array(atomAnalysis.get_bonds('C', 'C')[0])

bondRegistration = -1* np.ones(len(bonds))

bondRegistration = find_bond_registration(offset_distancesU, offset_distancesD, bonds, bondRegistration)


#Finally must check what category to place into based on the posPairs in up and down directions
bondRegistration = find_stacking_type(bondRegistration, posUs, posDs,bonds)

inFile = args.input.split('.')[0] + '_registration-measure'

if args.output != "<input>_registration-measure":
    dumpOut = args.output + '.dump'
    xyzOut = args.output + '.xyz'
else:
    dumpOut = inFile + '.dump'
    xyzOut = inFile + '.xyz'

    
write_registration(bonds, bondRegistration,dumpOut)



#For extra functionality we can compute a per-atom registration which doesn't tell any info about the type of stacking, but just how good it is
#We will take the mean registration measure of each bond involving the atom after find the modulus by 1 which will reduce it to the registration
#measure from 0 to 1.

atom_registrations = find_atom_registrations(atoms, nl, bonds, bondRegistration)
atoms.arrays['registration_mean'] = np.array(atom_registrations)
atoms.arrays['interlayer_dist']=np.array(distances)

mean_interlayer, dList = calc_interlayer_mean(atoms)
mean_a, aList = calc_2nd_neighbour_mean(atoms, nl)
print(f"Interlayer mean: {mean_interlayer:.3f} A, a lattice value mean: {mean_a:.3f} A")


write_extxyz(xyzOut, atoms, columns=['symbols','positions','registration_mean','interlayer_dist'])

print(atom_count, " atoms")

#optional extra file write of more details from calculation
if detail:
    bondLs = []
    normalXs = []
    normalYs = []
    normalZs = []
    thetas = []


    for i in range(atom_count):
        neighbours = find_neighbours(nl,i)
        if len(neighbours)==3:
            bondLs = np.append(bondLs,mean_bondL(atoms,i,nl))
            normal = find_normal(atoms,i,neighbours)
            normal_up = make_norm_pos(normal)
            normX = normal_up[0]
            normY = normal_up[1]
            normZ = normal_up[2]


            normalXs = np.append(normalXs,normX)
            normalYs = np.append(normalYs,normY)
            normalZs = np.append(normalZs,normZ)

        else:
            # normals = normals.append([0,0,-1])
            bondLs = bondLs.append(-1)

    atoms.arrays['neighbourhood_bondLs']=np.array(bondLs)
    atoms.arrays['normal_x']=np.array(normalXs)
    atoms.arrays['normal_y']=np.array(normalYs)
    atoms.arrays['normal_z']=np.array(normalZs)


    outFile=xyzOut.replace('.xyz','_detailed.xyz')
    normals = np.transpose(np.array([normX,normY,normZ]))
    #TODO: add angle deviation from c-axis
    # for i in range(np.shape(normals)[0]):
    #     theta = np.rad2deg(np.arccos(np.dot(normals[i],[0,0,1])/np.linalg.norm(normals[i])))
    #     thetas.append(theta)
    # atoms.arrays['theta']=np.array(thetas)

    write_extxyz(outFile, atoms, columns=['symbols','positions','interlayer_dist', 'neighbourhood_bondLs', 'normal_x', 'normal_y', 'normal_z'])

end = time.time()
elapsed = end- start
print(f"Elapsed time = {elapsed:.2f} seconds")
