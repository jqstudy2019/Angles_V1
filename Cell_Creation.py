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
from numpy import nan
import sys
import time as sleep 
import math

from ovito.io import import_file
from ovito.vis import Viewport

########################################################################################################################
##      Self Defined Functions      ##
##--------------------------------------------------------------------------------------------------------------------##
def angle_3_pts(pointA=list, pointB=list, pointC=list):
    # Convert the points to NumPy arrays
    pointA = np.array(pointA)
    pointB = np.array(pointB)
    pointC = np.array(pointC)
    
    # Calculate vectors AB and BC
    vector_AB = pointB - pointA
    vector_BC = pointC - pointB

    # Calculate dot product of AB and BC
    dot_product = np.dot(vector_AB, vector_BC)
    #print('Dot Prod:', dot_product)

    # Calculate magnitudes of AB and BC
    magnitude_AB = np.linalg.norm(vector_AB)
    magnitude_BC = np.linalg.norm(vector_BC)
    #print('Magnitues:',magnitude_AB,magnitude_BC)
    
# Calculate the cosine of the angle
    cos_theta = dot_product / 1
    a= ((magnitude_AB * magnitude_BC) + 1e-6)
    print(a)
    # Calculate the angle in radians
    theta_radians = np.arccos(cos_theta)
    
    # Convert the angle to degrees
    theta_degrees = np.degrees(theta_radians)

    return float(round(theta_degrees,2))

   


##--------------------------------------------------------------------------------------------------------------------##
def remove_NaNs(oldlist):
    newlist = list()
    for element in oldlist:
        if str(element) != 'nan':
            newlist.append(element)
    
    return newlist

##--------------------------------------------------------------------------------------------------------------------##
def are_sets_equal(set1, set2):
    tuple1 = tuple(sorted(set1))
    tuple2 = tuple(sorted(set2))
    return tuple1 == tuple2

##--------------------------------------------------------------------------------------------------------------------##
def check_list_for_duplicates(set_list, target_set):
    for s in set_list:
        if are_sets_equal(s, target_set):
            return True  # Found a duplicate set
    return False  # No duplicates found

##--------------------------------------------------------------------------------------------------------------------##
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
theta_range_len = len(theta_range)

phi_range = np.linspace(0, 2 * np.pi)
phi_range_len = len(phi_range)

rho_range = np.linspace(0, rho_range_end)
rho_range_len = len(rho_range)

total_bond_types = 1

total_angle_types = 1

total_dihedral_types = 1

max_iter = theta_range_len * phi_range_len

iter_counter = 0 
##--------------------------------------------------------------------------------------------------------------------##
#   Angle Creation Testing    #

Membrane_Angles = {}

Membrane_Num_Angles = 0

Membrane_Angle_Creation_Dist = 1

Membrane_Angle_Type = 1

angles_aggre = []

avg_angle_prev  =  141.27

angle_create_range = 5

max_angle_per_atom  = 5 

checklist = []

print('\n')
########################################################################################################################
##      Geometry Calculations       ##
for Theta in theta_range:
    for Phi in phi_range:
        sys.stdout.write('\r')
        iter_counter = iter_counter + 1
        percent_done = (round((iter_counter/max_iter),4))*100
        sys.stdout.write("\rGeometry Creation Progress: {:.2f}% done".format(percent_done))
        sys.stdout.flush()
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
                                                    angles_in_atom = 0 
                                                    angle_atom_ID = find_atom_id(x2,y2,z2)
                                                    angle_theta = angle_3_pts(temp_x_y_z,[x1,y1,z1],[x2,y2,z2])
                                                    angles_aggre.append(angle_theta)
                                                    if not math.isnan(angle_theta):
                                                        if avg_angle_prev - angle_create_range <= angle_theta <= avg_angle_prev + angle_create_range:
                                                            if angle_atom_ID is not None:
                                                                if bond_atom_ID is not None:
                                                                    if angle_atom_ID != Total_Number_Atoms:
                                                                        if angle_atom_ID != bond_atom_ID:
                                                                            if not check_list_for_duplicates(checklist,[Total_Number_Atoms,bond_atom_ID,angle_atom_ID]):
                                                                                checklist.append([Total_Number_Atoms,bond_atom_ID,angle_atom_ID])
                                                                                angles_in_atom = angles_in_atom +1
                                                                                Membrane_Num_Angles = Membrane_Num_Angles + 1
                                                                                Membrane_Angles[Membrane_Num_Angles] = {'Angle_ID':Membrane_Num_Angles,
                                                                                                                        'Angle Type': Membrane_Angle_Type,
                                                                                                                        'atom1':Total_Number_Atoms,
                                                                                                                        'atom2':bond_atom_ID,
                                                                                                                        'atom3':angle_atom_ID,
                                                                                                                        'Angle Theta': angle_theta}
                                                                                if angles_in_atom >= max_angle_per_atom:
                                                                                    break
                                                if angles_in_atom >= max_angle_per_atom:
                                                    break
                                            if angles_in_atom >= max_angle_per_atom:
                                                break
                                                                                                 
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

    fdata.write('{} atom types\n'.format(Total_Atom_Types))

    # Bonds Header  #
    # fdata.write('{} bonds\n'.format(num_molecule_1_bonds))  # Specify number of atoms

    fdata.write('{} bond types\n'.format(total_bond_types))
    fdata.write('{} angle types\n'.format(total_angle_types))
   # fdata.write('{} dihedral types\n'.format(total_dihedral_types))

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

# ########################################################################################################################
###         Finished!!!     ###
angles_aggre = remove_NaNs(angles_aggre)
avg_ang = sum(angles_aggre)/len(angles_aggre)

thetas = [entry['Angle Theta'] for entry in Membrane_Angles.values()]

print('\n\n Data File Created Successfully!!! ; File name => {} '.format(filename))

print('\n\n\n####  Angles Debugging  #####',
      '\n  Average Angle for all Particles (Includes non created angles): ', round(avg_ang,2),
      '\n  Average Angle of created Angles (only saved angles) ', round(sum(thetas)/len(thetas),2),
      #'\n  Average Angles Per Atom : ', sum(angles_per_atom)/len(angles_per_atom),
      '\n\n  ##  Angle Creation Params  ##',
      '\n # of Angles Considered: ',len(angles_aggre),
      '\n # of Angles Created: ', len(Membrane_Angles))
