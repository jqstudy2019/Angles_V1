########################################################################################################################
#   REGION_TEST_V1.py ----- A code to build the .data file for a LAMMPS Simulation
#
#       8/28/2020
#
########################################################################################################################
##      Import Extra Modules / Tools        ##
import random as random

import numpy as np  # Numpy is used for Array Creation and Manipulation Mainly


########################################################################################################################
##      Self Defined Functions      ##
def find_atom_id(x_find, y_find, z_find):
    for iV3, posV3 in enumerate(Total_xyz):
        if (Total_xyz[posV3].get('X') == x_find) \
                and (Total_xyz[posV3].get('Y') == y_find) \
                and (Total_xyz[posV3].get('Z') == z_find):
            return Total_xyz[posV3].get('Atom_ID')


##--------------------------------------------------------------------------------------------------------------------##
def find_atom_type(atom_id):
    for iV5, posV5 in enumerate(Total_xyz):
        if Total_xyz[posV5].get('Atom_ID') == atom_id:
            temp_id = Total_xyz[posV5].get('Atom_Type')
            return temp_id


##--------------------------------------------------------------------------------------------------------------------##
def packing_roll(coeff):
    temp_roll_func = random.randint(1, 100)

    if temp_roll_func in range(0, coeff):
        return 'Place Atom'

    else:
        return 'Leave Empty'


##--------------------------------------------------------------------------------------------------------------------##
def cart2sphere(x4, y4, z4):
    rho = np.sqrt((x4 ** 2) + (y4 ** 2) + (z4 ** 2))
    theta = np.arctan((y4 / x4))
    phi = np.arccos(z4 / (np.sqrt((x4 ** 2) + (y4 ** 2) + (z4 ** 2))))

    return [rho, theta, phi]


##--------------------------------------------------------------------------------------------------------------------##
def sphere2cart(rho, theta, phi):
    x5 = rho * np.sin(phi) * np.cos(theta)
    y5 = rho * np.sin(phi) * np.sin(theta)
    z5 = rho * np.cos(phi)

    return [x5, y5, z5]


##--------------------------------------------------------------------------------------------------------------------##
########################################################################################################################
##      Set Constants and Variables Needed      ##
##--------------------------------------------------------------------------------------------------------------------##
#   Output File Options   #

filename = input("Enter a name for the Output File (remember .data at the end!) :")  # String Name of File to be Created

##--------------------------------------------------------------------------------------------------------------------##
#   Sim Box Settings    #
sim_box_side_length = 100  # Each Side of the Simulation Box Will have sides of these lengths

SIM_SPACE = np.zeros((sim_box_side_length, sim_box_side_length, sim_box_side_length))
# Creates an Array Called SIM_SPACE (is a 3-D array (x,y,z of each simulation space position)


x0 = int(np.floor(SIM_SPACE.shape[2] / 2))  # Finds the Center of the Sim Space ( x0, y0, z0 )
y0 = int(np.floor(SIM_SPACE.shape[1] / 2))
z0 = int(np.floor(SIM_SPACE.shape[0] / 2))

Sim_Box_Longest_Length = int(np.sqrt(((sim_box_side_length / 2) ** 2) + ((sim_box_side_length / 2) ** 2)))
# This should be the longest Rho Value Possible / The Distance from the center of the box to one of the corners

atom_place_limit = 90  # limits where atoms maybe placed based on box dimensions

##--------------------------------------------------------------------------------------------------------------------##
#   Membrane Settings   #
Membrane_Atom_Type = 1

Membrane_Molecule_Type = 1

Membrane_Thickness = 1

Membrane_Radius = 5

Membrane_Diameter = Membrane_Radius * 2

Membrane_Number_Atoms = 0

Membrane_xyz = {}

Membrane_Bonds = {}

Membrane_Number_Bonds = 0

Membrane_Bond_Type = 1

Membrane_Bond_1_Length = 1

##--------------------------------------------------------------------------------------------------------------------##
#   Water Settings  #

Water_Number_Atoms = 0

Water_xyz = {}

Water_Atom_Type = 2

Water_Molecule_Type = 2

Water_Packing_Coeff = 5

Water_Offset = 2

##--------------------------------------------------------------------------------------------------------------------##
#   Plasma Settings #

Plasma_Number_Atoms = 0

Plasma_xyz = {}

Plasma_Atom_Type = 3

Plasma_Molecule_Type = 3

Plasma_Packing_Coeff = 5

Plasma_Offset = 2

##--------------------------------------------------------------------------------------------------------------------##
#   Total / Aggregate Variables     #

Total_Number_Atoms = 0

Total_xyz = {}

Total_Atom_Types = 3

Total_Molecule_Type = 3

##--------------------------------------------------------------------------------------------------------------------##
#   Calculated Constants    #

rho_range_end = sim_box_side_length / 2

theta_range = np.linspace(0, 2 * np.pi)

phi_range = np.linspace(0, 2 * np.pi)

rho_range = np.linspace(0, rho_range_end)

# membrane_range = np.linspace(Membrane_Radius, Membrane_Radius + Membrane_Thickness)


total_bond_types = 1

total_angle_types = 1

total_dihedral_types = 1


bonds_per_atom_limit = 5000

########################################################################################################################
##      Geometry Calculations       ##

num = 0

for Theta in theta_range:
    for Phi in phi_range:
        for Rho in rho_range:

            # Membrane Section
            if Rho in np.linspace(0, 10):

                temp_roll = packing_roll(Plasma_Packing_Coeff)

                if temp_roll == 'Place Atom':
                    temp_x_y_z = sphere2cart(Rho, Theta, Phi)

                    temp_x = int(temp_x_y_z[0])
                    temp_y = int(temp_x_y_z[1])
                    temp_z = int(temp_x_y_z[2])

                    temp_x = temp_x + x0
                    temp_y = temp_y + y0
                    temp_z = temp_z + z0

                    if (temp_x and temp_y and temp_z) in range(0, atom_place_limit):
                        SIM_SPACE[temp_x, temp_y, temp_z] = Plasma_Atom_Type

                        Plasma_Number_Atoms = Plasma_Number_Atoms + 1
                        Total_Number_Atoms = Total_Number_Atoms + 1

                        Plasma_xyz[Plasma_Number_Atoms] = {'Atom_ID': Total_Number_Atoms,
                                                           'Molecule_ID': Plasma_Molecule_Type,
                                                           'Atom_Type': Plasma_Atom_Type,
                                                           'X': temp_x,
                                                           'Y': temp_y,
                                                           'Z': temp_z,
                                                           'Positions String': '{} {} {}'.format(temp_x,
                                                                                                 temp_y,
                                                                                                 temp_z)}

                        Total_xyz[Total_Number_Atoms] = {'Atom_ID': Total_Number_Atoms,
                                                         'Molecule_ID': Plasma_Molecule_Type,
                                                         'Atom_Type': Plasma_Atom_Type,
                                                         'X': temp_x,
                                                         'Y': temp_y,
                                                         'Z': temp_z,
                                                         'Positions String': '{} {} {}'.format(temp_x,
                                                                                               temp_y,
                                                                                               temp_z)}
            # Membrane Section
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

                    Membrane_Number_Atoms = Membrane_Number_Atoms + 1
                    Total_Number_Atoms = Total_Number_Atoms + 1

                    Membrane_xyz[Membrane_Number_Atoms] = {'Atom_ID': Total_Number_Atoms,
                                                           'Molecule_ID': Membrane_Molecule_Type,
                                                           'Atom_Type': Membrane_Atom_Type,
                                                           'X': temp_x,
                                                           'Y': temp_y,
                                                           'Z': temp_z,
                                                           'Positions String': '{} {} {}'.format(temp_x,
                                                                                                 temp_y,
                                                                                                 temp_z)}

                    Total_xyz[Total_Number_Atoms] = {'Atom_ID': Total_Number_Atoms,
                                                     'Molecule_ID': Membrane_Molecule_Type,
                                                     'Atom_Type': Membrane_Atom_Type,
                                                     'X': temp_x,
                                                     'Y': temp_y,
                                                     'Z': temp_z,
                                                     'Positions String': '{} {} {}'.format(temp_x,
                                                                                           temp_y,
                                                                                           temp_z)}

                    for Rho1 in range(Rho - Membrane_Bond_1_Length, Rho + Membrane_Bond_1_Length):
                        Bond_Neighbors = find_neighbors(Total_Number_Atoms)
                        Membrane_Number_Bonds = Membrane_Number_Bonds + 1
                        Membrane_Bonds[Membrane_Number_Bonds] = {'Bond_ID': Membrane_Number_Bonds,
                        'Bond_Type': Membrane_Bond_Type,
                        'atom1': Temp_Atom_ID_1,
                        'atom2': Temp_Atom_ID_2,
                        'bond_length': molecule_1_bond_type_1_length}

            # Water Section
            if Rho in np.linspace(20, 40):

                temp_roll = packing_roll(Water_Packing_Coeff)

                if temp_roll == 'Place Atom':
                    temp_x_y_z = sphere2cart(Rho, Theta, Phi)

                    temp_x = int(temp_x_y_z[0])
                    temp_y = int(temp_x_y_z[1])
                    temp_z = int(temp_x_y_z[2])

                    temp_x = temp_x + x0
                    temp_y = temp_y + y0
                    temp_z = temp_z + z0

                    if (temp_x and temp_y and temp_z) in range(0, atom_place_limit):
                        SIM_SPACE[temp_x, temp_y, temp_z] = Water_Atom_Type

                        Water_Number_Atoms = Water_Number_Atoms + 1
                        Total_Number_Atoms = Total_Number_Atoms + 1

                        Membrane_xyz[Membrane_Number_Atoms] = {'Atom_ID': Total_Number_Atoms,
                                                               'Molecule_ID': Water_Molecule_Type,
                                                               'Atom_Type': Water_Atom_Type,
                                                               'X': temp_x,
                                                               'Y': temp_y,
                                                               'Z': temp_z,
                                                               'Positions String': '{} {} {}'.format(temp_x,
                                                                                                     temp_y,
                                                                                                     temp_z)}

                        Total_xyz[Total_Number_Atoms] = {'Atom_ID': Total_Number_Atoms,
                                                         'Molecule_ID': Water_Molecule_Type,
                                                         'Atom_Type': Water_Atom_Type,
                                                         'X': temp_x,
                                                         'Y': temp_y,
                                                         'Z': temp_z,
                                                         'Positions String': '{} {} {}'.format(temp_x,
                                                                                               temp_y,
                                                                                               temp_z)}
                else:
                    pass
########################################################################################################################
###     Write LAMMPS Data File      ###
with open(filename, 'w+') as fdata:  # opens a text file named a for the 'filename' variable
    fdata.write('{}\n\n'.format(filename))  # First line is a comment line

    ##     Header of Data File     ##

    #   Atoms Header #

    fdata.write('{} atoms\n'.format(Total_Number_Atoms))  # Specify number of atoms
    fdata.write('{} atom types\n'.format(Total_Atom_Types))

    fdata.write('{} bonds per atom\n'.format(bonds_per_atom_limit))

    # Bonds Header  #
    # fdata.write('{} bonds\n'.format(num_molecule_1_bonds))  # Specify number of atoms

    fdata.write('{} bond types\n'.format(total_bond_types))
    fdata.write('{} angle types\n'.format(total_angle_types))
    fdata.write('{} dihedral types\n'.format(total_dihedral_types))

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

########################################################################################################################
###         Finished!!!     ###
print(' Data File Created Successfully!!! ; File name => {} '.format(filename))
