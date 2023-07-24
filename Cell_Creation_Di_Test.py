########################################################################################################################
#   Cell_Creation.py -- a Scrpti to make a LAMMPS Data file 
#
#       Author: John Quinn  
#
#       Created : 7/24/23
#
#       Last Mooified: 7/24/23
#
#   Refrences: 
#       LAMMPS Documentation 
#           https://docs.lammps.org/Manual.html
#
#       LAMMPS Python API Documentaion 
#           https://docs.lammps.org/Python_module.html
#
#       Ovito API 
#           https://www.ovito.org/docs/current/python/
#
########################################################################################################################
##      Import Extra Modules / Tools        ##
import random as random
import numpy as np  # Numpy is used for Array Creation and Manipulation Mainly

from ovito.io import import_file
from ovito.vis import Viewport

########################################################################################################################
##      Self Defined Functions      ##
def find_atom_id(x_find, y_find, z_find,verbose= False):
    for iV3, posV3 in enumerate(Total_xyz):
        if (Total_xyz[posV3].get('X') == x_find) \
                and (Total_xyz[posV3].get('Y') == y_find) \
                and (Total_xyz[posV3].get('Z') == z_find):
            return Total_xyz[posV3].get('Atom_ID')
        elif verbose == True:
            print('No Atom Found in that position!')
        else: 
            pass
        
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
def cart2sphere(x, y, z):
    rho = np.sqrt((x ** 2) + (y ** 2) + (z ** 2))
    theta = np.arctan((y / x))
    phi = np.arccos(z / (np.sqrt((x ** 2) + (y ** 2) + (z ** 2))))

    return [rho, theta, phi]

##--------------------------------------------------------------------------------------------------------------------##
def sphere2cart(rho, theta, phi):
    x = rho * np.sin(phi) * np.cos(theta)
    y = rho * np.sin(phi) * np.sin(theta)
    z = rho * np.cos(phi)

    return [x, y, z]

##--------------------------------------------------------------------------------------------------------------------##
########################################################################################################################
##      Set Constants and Variables Needed      ##
##--------------------------------------------------------------------------------------------------------------------##
#   Output File Options   #

filename = input("Enter a name for the Output File (.data will be added) :") + '.data' # String Name of File to be Created

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

Membrane_Bond_1_Length = 2

##--------------------------------------------------------------------------------------------------------------------##
#   Total / Aggregate Variables     #

Total_Number_Atoms = 0

Total_xyz = {}

Total_Atom_Types = 1

Total_Molecule_Type = 1

##--------------------------------------------------------------------------------------------------------------------##
#   Calculated Constants    #

rho_range_end = sim_box_side_length / 2

theta_range = np.linspace(0, np.pi) # Mofified per Nate and Chat GPT 

phi_range = np.linspace(0, 2 * np.pi)

rho_range = np.linspace(0, rho_range_end)

total_bond_types = 1

total_angle_types = 1

total_dihedral_types = 1

##--------------------------------------------------------------------------------------------------------------------##
#   Angle Creation Testing    #

Membrane_Angles = {}

Membrane_Num_Angles = 0

Membrane_Angle_Creation_Dist = 1

Membrane_Angle_Type = 1

##--------------------------------------------------------------------------------------------------------------------##
#   Dihederal Creation Testing    #

Membrane_Dihederals = {}

Membrane_Numer_Dihederals = 0

Membrane_Dihederal_Creation_Dist = 1

Membrane_Dihederal_Type = 1

########################################################################################################################
##      Geometry Calculations       ##
for Theta in theta_range:
    for Phi in phi_range:
        for Rho in rho_range:

            # Membrane Section
            if Rho in np.linspace(15, 16):
                temp_x_y_z = sphere2cart(Rho, Theta, Phi)

                temp_x = int(temp_x_y_z[0]) + x0
                temp_y = int(temp_x_y_z[1]) + y0
                temp_z = int(temp_x_y_z[2]) + z0

                if (temp_x and temp_y and temp_z) in range(0, atom_place_limit):
                    SIM_SPACE[temp_x, temp_y, temp_z] = Membrane_Atom_Type

                    Membrane_Number_Atoms = Membrane_Number_Atoms + 1
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
                    # Bonds Creation #
                    for x1 in range(temp_x - Membrane_Bond_1_Length,temp_x + Membrane_Bond_1_Length):
                        for y1 in range(temp_y - Membrane_Bond_1_Length,temp_y + Membrane_Bond_1_Length):
                            for z1 in range(temp_z - Membrane_Bond_1_Length,temp_z + Membrane_Bond_1_Length):
                                bond_atom_ID = find_atom_id(x1,y1,z1)
                                if bond_atom_ID is not None:
                                    if bond_atom_ID != Total_Number_Atoms:
                                        Membrane_Number_Bonds = Membrane_Number_Bonds + 1
                                        Membrane_Bonds[Membrane_Number_Bonds] = {'Bond_ID':Membrane_Number_Bonds,
                                                            'Bond Type': Membrane_Bond_Type,
                                                            'atom1':Total_Number_Atoms,
                                                            'atom2':bond_atom_ID}
                    
                                        #   Angles Creation #
                                        for x2 in range(x1 - Membrane_Angle_Creation_Dist,x1 + Membrane_Angle_Creation_Dist):
                                            for y2 in range(y1 - Membrane_Angle_Creation_Dist,y1 + Membrane_Angle_Creation_Dist):
                                                for z2 in range(z1 - Membrane_Angle_Creation_Dist,z1 + Membrane_Angle_Creation_Dist):
                                                    angle_atom_ID = find_atom_id(x2,y2,z2)
                                                    if angle_atom_ID is not None:
                                                        if angle_atom_ID != Total_Number_Atoms:
                                                            if angle_atom_ID != bond_atom_ID:
                                                                Membrane_Num_Angles = Membrane_Num_Angles + 1
                                                                Membrane_Angles[Membrane_Num_Angles] = {'Angle_ID':Membrane_Num_Angles,
                                                                                    'Angle Type': Membrane_Angle_Type,
                                                                                    'atom1':Total_Number_Atoms,
                                                                                    'atom2':bond_atom_ID,
                                                                                    'atom3':angle_atom_ID}
                                                        
                                                                #   Dihederals Creation #
                                                                for x3 in range(x2 - Membrane_Dihederal_Creation_Dist,x2 + Membrane_Dihederal_Creation_Dist):
                                                                    for y3 in range(y2 - Membrane_Dihederal_Creation_Dist,y2 + Membrane_Dihederal_Creation_Dist):
                                                                        for z3 in range(z2 - Membrane_Dihederal_Creation_Dist,z2 + Membrane_Dihederal_Creation_Dist):
                                                                            Dihederal_atom_ID = find_atom_id(x3,y3,z3)
                                                                            if Dihederal_atom_ID is not None:
                                                                                if Dihederal_atom_ID != Total_Number_Atoms:
                                                                                    if Dihederal_atom_ID != bond_atom_ID:
                                                                                        if Dihederal_atom_ID != angle_atom_ID:
                                                                                            Membrane_Numer_Dihederals = Membrane_Numer_Dihederals + 1
                                                                                            Membrane_Dihederals[Membrane_Numer_Dihederals] = {'Angle_ID':Membrane_Numer_Dihederals,
                                                                                                                'Angle Type': Membrane_Dihederal_Type,
                                                                                                                'atom1':Total_Number_Atoms,
                                                                                                                'atom2':bond_atom_ID,
                                                                                                                'atom3':angle_atom_ID,
                                                                                                                'atom4':Dihederal_atom_ID}

                else:
                    pass
########################################################################################################################
###     Write LAMMPS Data File      ###
with open(filename, 'w+') as fdata:  # opens a text file named a for the 'filename' variable
    fdata.write('{}\n\n'.format(filename))  # First line is a comment line

    ##     Header of Data File     ##

    #   Atoms Header #

    fdata.write('{} atoms\n'.format(Total_Number_Atoms))  # Specify number of atoms
    fdata.write('{} bonds\n'.format(Membrane_Number_Bonds))
    fdata.write('{} angles\n'.format(Membrane_Num_Angles))
    fdata.write('{} dihederals\n'.format(Membrane_Numer_Dihederals))
    
    fdata.write('\n')
    fdata.write('{} atom types\n'.format(Total_Atom_Types))
    fdata.write('{} bond types\n'.format(total_bond_types))
    fdata.write('{} angle types\n'.format(total_angle_types))
    fdata.write('{} dihedral types\n'.format(total_dihedral_types))

    #   specify box dimensions      #
    fdata.write('\n')
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

    #       Bonds section       #
    fdata.write('Bonds\n\n')
    for i, pos in enumerate(Membrane_Bonds):
        fdata.write('{} {} {} {} \n'.format(Membrane_Bonds[pos].get('Bond_ID'),
                                            Membrane_Bonds[pos].get('Bond Type'),
                                            Membrane_Bonds[pos].get('atom1'),
                                            Membrane_Bonds[pos].get('atom2')))    
    fdata.write('\n')

    #       Angles section       #
    fdata.write('Angles\n\n')
    for i, pos in enumerate(Membrane_Angles):
        fdata.write('{} {} {} {} {} \n'.format(Membrane_Angles[pos].get('Angle_ID'),
                                            Membrane_Angles[pos].get('Angle Type'),
                                            Membrane_Angles[pos].get('atom1'),
                                            Membrane_Angles[pos].get('atom2'),
                                            Membrane_Angles[pos].get('atom3')))    
    fdata.write('\n')

    #       Dihederals section       #
    fdata.write('Dihedreals\n\n')
    for i, pos in enumerate(Membrane_Dihederals):
        fdata.write('{} {} {} {} {} {} \n'.format(Membrane_Dihederals[pos].get('Angle_ID'),
                                            Membrane_Dihederals[pos].get('Angle Type'),
                                            Membrane_Dihederals[pos].get('atom1'),
                                            Membrane_Dihederals[pos].get('atom2'),
                                            Membrane_Dihederals[pos].get('atom3'),
                                            Membrane_Dihederals[pos].get('atom4')))    
    fdata.write('\n')
########################################################################################################################
###         Ovito Visulalization and Analyis        ###

# pipeline = import_file(filename, atom_style='molecular')
# pipeline.add_to_scene()

# vp = Viewport()
# vp.type = Viewport.Type.Perspective
# vp.camera_pos = (75,75,75)
# vp.camera_dir = (2,3,-3)
# vp.fov = (60.0)

# image_filename = filename + 'imgae.png'
# vp.render_image(filename=image_filename,size=(1920,1080))

########################################################################################################################
###         Finished!!!     ###
print(' Data File Created Successfully!!! ; File name => {} '.format(filename))
