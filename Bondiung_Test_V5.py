########################################################################################################################
#   Bonding_Test_V5.py ----- A code to build the .data file for a LAMMPS Simulation ; WIth new bonding Scheme
                                # this one hopefully fixes some file formatting errors and reading by lammps/ovito

                                # Trying to add the cartesian to spherical coordinate arguments to guarantee optimal reults going forward 
#
#   General Description - New and Improved Methods implemented and now\
#                           build new bonding scheme
#
#This
########################################################################################################################

##      Import Extra Modules / Tools        ##

import numpy as np  # Numpy is used for Array Creation and Manipulation Mainly
from copy import deepcopy  # DeepCopy Is used to Copy Arrays
import random as random


########################################################################################################################
##      Self Defined Functions      ##
def find_atom_id(x_find, y_find, z_find):
    for iV2, posV2 in enumerate(molecule_1):
        if (molecule_1[posV2].get('X') == x_find) \
                and (molecule_1[posV2].get('Y') == y_find) \
                and (molecule_1[posV2].get('Z') == z_find):
            return molecule_1[posV2].get('Atom_ID')

def packing_roll(coeff):
    temp_roll_func = random.randint(1, 100)

    if temp_roll_func in range(0, coeff):
        return 'Place Atom'

    else:
        return 'Leave Empty'
    

def cart2sphere(x4, y4, z4):
    rho = np.sqrt((x4 ** 2) + (y4 ** 2) + (z4 ** 2))
    theta = np.arctan((y4 / x4))
    phi = np.arccos(z4 / (np.sqrt((x4 ** 2) + (y4 ** 2) + (z4 ** 2))))

    return [rho, theta, phi]

def sphere2cart(rho, theta, phi):
    x5 = rho * np.sin(phi) * np.cos(theta)
    y5 = rho * np.sin(phi) * np.sin(theta)
    z5 = rho * np.cos(phi)

    return [x5, y5, z5]

########################################################################################################################

###      Set Variables Used For Later       ###


#   Output File Options #

filename = input("Enter a name for the Output File (remember .data at the end!) :")  # String Name of File to be Created

#   Simulation Box Settings   #

sim_box_side_length = 100  # Each Side of the Simulation Box Will have sides of these lengths

SIM_SPACE = np.zeros((sim_box_side_length, sim_box_side_length, sim_box_side_length))
# Creates an Array Called SIM_SPACE (is a 3-D array (x,y,z of each simulation space position)

SIM_SPACE_PUREcopy = deepcopy(SIM_SPACE)  # Extra Copy of the Sim Space kept because its "pure"

x0 = int(np.floor(SIM_SPACE.shape[2] / 2))  # Finds the Center of the Sim Space ( x0, y0, z0 )
y0 = int(np.floor(SIM_SPACE.shape[1] / 2))
z0 = int(np.floor(SIM_SPACE.shape[0] / 2))

#   Molecule Geometry Creation Settings #
molecule_ids = 1  # LAMMPS Track molecules by their ID number ; this is the number of molecules created

molecule_1_thickness = 15  # The thickness is an extra range past the radius to include in the molecules

molecule_1_radius = 10  # Radius of each molecule

num_atoms_molecule_1 = 0  # Placeholder for total number of atoms created as part of the molecule

molecule_1_bond_length = 15  # Functions Similar to molecule radius but on the bond scale

molecule_1 = {}  # initializes the dictionary for molecule 1

molecule_1_COPY = {}  # Copy of Molecule dictionary

molecule_1_bonds = {}  # dictionary for bonds

Membrane_xyz = {}

Total_xyz = {}

num_molecule_1_bonds = 0

angle_1_bonds = {}

Membrane_Atom_Type =1

Membrane_Molecule_Type =1

num_angle_1_bonds = 0

atom_place_limit = 90

Membrane_Number_Atoms = 0

Total_Number_Atoms = 0
# Constants for spherical coords

rho_range_end = sim_box_side_length / 2

theta_range = np.linspace(0, 2 * np.pi)

phi_range = np.linspace(0, 2 * np.pi)

rho_range = np.linspace(0, rho_range_end)

#   LAMMPS Data Inserts #

atom_types = 1  # LAMMPS Data Field ; Number of different types of Atoms

bond_types = 1  # LAMMPS Data Field ; Number of different types of Bonds

#   Program Performance Tracking

positions_checked = 0  # Keeps Count of How many points are iterated through ; easy way to check it iterating correctly

########################################################################################################################

###     Creates 3-D Set of Points      ###

# Iterates through the spherical coordinates

num = 0

for Theta in theta_range:
    for Phi in phi_range:
        for Rho in rho_range:
            if Rho in np.linspace(15, 16):
                temp_x_y_z = sphere2cart(Rho, Theta, Phi)

                temp_x = int(temp_x_y_z[0])
                temp_y = int(temp_x_y_z[1])
                temp_z = int(temp_x_y_z[2])

                temp_x = temp_x + x0
                temp_y = temp_y + y0
                temp_z = temp_z + z0

                if (temp_x and temp_y and temp_z) in range(0, atom_place_limit):
                    SIM_SPACE[temp_x, temp_y, temp_z] = Membrane_Atom_Type

                    Total_Number_Atoms = Total_Number_Atoms + 1

                    Total_xyz[Total_Number_Atoms] = {'Atom_ID': Total_Number_Atoms,
                                                     'Molecule_ID': Membrane_Molecule_Type,
                                                     'Atom_Type': Membrane_Atom_Type,
                                                     'X': temp_x,
                                                     'Y': temp_y,
                                                     'Z': temp_z,
                                                     'Positions String': '{} {} {}'.format(temp_x,
                                                                                           temp_y,
                                                                                           temp_z)}
                    counter=0               
                    for x1 in range(temp_x - molecule_1_bond_length, temp_x + molecule_1_bond_length):
                        for y1 in range(temp_y - molecule_1_bond_length, temp_y + molecule_1_bond_length):
                            for z1 in range(temp_z - molecule_1_bond_length, temp_z + molecule_1_bond_length):
                                Atom_ID_2 = find_atom_id(x1, y1, z1)
                                if Atom_ID_2 is not None:
                                    if Atom_ID_2 != num_atoms_molecule_1:
                                        num_molecule_1_bonds = num_molecule_1_bonds + 1
                                        molecule_1_bonds[num_molecule_1_bonds] = {'Bond_ID': num_molecule_1_bonds,
                                                                                'Bond_Type': bond_types,
                                                                                'atom1': num_atoms_molecule_1,
                                                                                'atom2': Atom_ID_2}
                                        
                                else:
                                   # print('nothing here folks')
                                    counter=counter + 1
########################################################################################################################

###     Prints Sim Info To Make Sure Its Running Correct        ###

print("Number of Points in Simulation Space Iterated thorough: {} \n".format(positions_checked))

print("number of counts: ")
print(counter)

#print("Number of Bonds Created: {} \n".format(molecule_1_bonds))

#print("Sim Space Array: {} \n".format(SIM_SPACE))

##########################################################################################################################

###     Write LAMMPS Data File      ###

with open(filename, 'w+') as fdata:  # opens a text file named a for the 'filename' variable
    fdata.write('{}\n\n'.format(filename))  # First line is a comment line

    ##     Header of Data File     ##

    #   Atoms Header #

    fdata.write('{} atoms\n'.format(Total_Number_Atoms))  # Specify number of atoms
    fdata.write('{} atom types\n'.format(atom_types))

    # Bonds Header  #
    fdata.write('{} bonds\n'.format(num_molecule_1_bonds))  # Specify number of atoms
    fdata.write('{} bond types\n'.format(bond_types))

    #   specify box dimensions      #

    fdata.write('{} {} xlo xhi\n'.format(0.0, sim_box_side_length))  # Writes X Position
    fdata.write('{} {} ylo yhi\n'.format(0.0, sim_box_side_length))  # Writes Y Position
    fdata.write('{} {} zlo zhi\n'.format(0.0, sim_box_side_length))  # Writes Z Position

    fdata.write('\n')  # Skips the Next line for Formatting

    ##      Data File Body      ##

    #       Atoms section       #
    fdata.write('Atoms\n\n')
    for i, pos in enumerate(Total_xyz):
        fdata.write('{} {} {} {} \n'.format(Total_xyz[pos].get('Atom_ID'),
                                            Total_xyz[pos].get('Molecule_ID'),
                                            Total_xyz[pos].get('Atom_Type'),
                                            Total_xyz[pos].get('Positions String')))
    fdata.write('\n')
    fdata.write('\n')

    #       Bonds Section       #
    fdata.write('Bonds\n\n')
    for iV2, posV2 in enumerate(molecule_1_bonds):
        fdata.write('{} {} {} {} \n'.format(molecule_1_bonds[posV2].get('Bond_ID'),
                                            molecule_1_bonds[posV2].get('Bond_Type'),
                                            molecule_1_bonds[posV2].get('atom1'),
                                            molecule_1_bonds[posV2].get('atom2')))

########################################################################################################################

###         Finished!!!     ###

print(' Data File Created Successfully!!! ; File name => {} '.format(filename))

########################################################################################################################
