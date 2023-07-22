# Python code to generate discretised morphing structure
# This code was written by Hirak Kansara 
# Corresponding author: Wei Tan (wei.tan@qmul.ac.uk), Liu Mingchao (mingchao.liu@ntu.edu.sg)

from __future__ import division
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
import math
import numpy as np
import os
from scipy import integrate

# Local directory
setPath = r'C:\Temp\Morphing_Structure\MorhComp\Morph_Gen'# 'C:\Temp\Morphing_Structure\MorhComp\Morph_Gen\shell_morph' # 'C:\Users\hirak\AppData\Local\Temp\Temp1_InverseTessellation_Hemiellipsoid.zip' # 'C:\Temp\Morphing_Structure\MorhComp\Morph_Gen'
os.chdir(setPath)

# HPC directory
cwd = os.getcwd()
os.chdir(cwd)

N = 8
eps = 0.1
Depth = 2.0 # mm
Width = 40 # mm
Length = 40    # mm
Square_Size = 0.5  # mm size of each cube based

patterns = ['rand-dist', 'rand-dist2', 'clumped-top', 'clumped-bot', 'half-half', 'half-half-mid', 'thickness-dist']

pattern = patterns[-1]  # random.randint(0, 4)

numSamp = Length/ Square_Size  # number of slices / number of samples
 # inner circle radius


Em = 3.14 # (MPa) Modulus of matrix i.e. Modulus of Tangoblack
Ef = 12.6# 2033.33
Ef_actual = Ef
Ei = 0.2*Ef

## thermal conductivity
k_soft = 0.3155
k_rigid = 2.01

## electrical conductivity
sigma_soft = 1e-7
sigma_rigid = 2.86e-23

def read_coeffs(a, b):
# Parameters - Bigger dimensions better the resolution

    # Fourier series coefficients for modulus fit
    str1 = 'vf_2mat' # 'vf_a_b_05'# 'E_fit_eps_'+str(eps)+'_a_'+str(a)+'_b_'+str(b)+'_N_'+str(N) # vf_newshape_negkappa vf_newshape_varkappa vf_newshape_nosecone
    str1 = str1.replace(".", "")
    with open(str1+'.txt') as f: # with open(str1+'.txt') as f:
        list_of_coeffs = [line.split(',') for line in f.readlines()]
        flat_list = [item for sublist in list_of_coeffs for item in sublist]

    extended = []
    for i in range(len(flat_list)):
        extended.extend(flat_list[i].rstrip('\n').split(','))

    e_coeffs = [float(i) for i in extended]

    # Fourier series coefficients for width fit
# w_newshape_negkappa w_newshape_varkappa w_newshape_nosecone
    str2 = 'w_fit_eps_01_a_1_b_1_N_8'# 'w_a_b_05'# 'w_fit_eps_'+str(eps)+'_a_'+str(a)+'_b_'+str(b)+'_N_'+str(N)# 'w_newshape_negkappa'#  ## 'w_fit_eps_'+str(eps)+'_a_'+str(a)+'_b_'+str(b)+'_N_'+str(N)
    str2 = str2.replace(".", "")

    with open(str2+'.txt') as f1: # with open(str2+'.txt') as f1:
        list_of_coeffs1 = [line.split(',') for line in f1.readlines()]
        flat_list1 = [item for sublist in list_of_coeffs1 for item in sublist]

    extended1 = []
    for i in range(len(flat_list1)):
        extended1.extend(flat_list1[i].rstrip('\n').split(','))

    w_coeffs = [float(i) for i in extended1]
    return e_coeffs, w_coeffs


x = list(np.linspace(0, 1, num=int(numSamp), endpoint=True))  # choosing intervals of x based on the number of slices

a_tot = 1
b_tot = 1
e_coeffs, w_coeffs = read_coeffs(a_tot, b_tot)

def f(x):

    # a0 = e_coeffs[0]
    # a1 = e_coeffs[1]
    # b1 = e_coeffs[2]
    # a2 = e_coeffs[3]
    # b2 = e_coeffs[4]
    # a3 = e_coeffs[5]
    # b3 = e_coeffs[6]
    # a4 = e_coeffs[7]
    # b4 = e_coeffs[8]
    # a5 = e_coeffs[9]
    # b5 = e_coeffs[10]
    # a6 = e_coeffs[11]
    # b6 = e_coeffs[12]
    # Omega = e_coeffs[13]

    a1 =       e_coeffs[0]
    b1 =      e_coeffs[1]
    c1 =      e_coeffs[2]
    a2 =    e_coeffs[3]
    b2 =      e_coeffs[4]
    c2 =     e_coeffs[5]
    a3 =     e_coeffs[6]
    b3 =      e_coeffs[7]
    c3 =     e_coeffs[8]
    a4 =   e_coeffs[9]
    b4 =     e_coeffs[10]
    c4 =     e_coeffs[11]
    a5 =   e_coeffs[12]
    b5 =      e_coeffs[13]
    c5 =     e_coeffs[14]
    a6 =      e_coeffs[15]
    b6 =    e_coeffs[16]
    c6 =     e_coeffs[17]
    a7 =    e_coeffs[18]
    b7 =     e_coeffs[19]
    c7 =    e_coeffs[20]
    a8 =      e_coeffs[21]
    b8 =       e_coeffs[22]
    c8 =     e_coeffs[23]
#     # nose cone
#     # a1 =0
#     # b1=0.969024267367378
#     # c1=0.0919312434846282
#     # a2=11763326939062.1
#     # b2=-2.50664267940204
#     # c2=0.456723812535920
#     # a3=-0.529596607698589
#     # b3=0.377058715219084
#     # c3=0.0282517876062011
#     # a4=0.436174984360824
#     # b4=0.371162522241522
#     # c4=0.0243594683862297
#     # a5=0
#     # b5=0.347092817332407
#     # c5=2.22044604925031e-14
#     # a6=0.0218735357848556
#     # b6=0.336489579949035
#     # c6=0.00880888846562960
#     # a7=0
#     # b7=1.71677957210599
#     # c7=0.000273688213538033
#     # a8=7.48471519753611
#     # b8=0.888000495638263
#     # c8=0.342824569160347
#
#     # a1 = 0.005412
#     # b1 = 0.4293
#     # c1 = 0.209
#     # a2 = 0.2549
#     # b2 = 0.5006
#     # c2 = 112.6
#     # a3 = -7.552e-06
#     # b3 = 0.2994
#     # c3 = 2.373e-05
#     # a4 = 3.451e-05
#     # b4 = 0.2962
#     # c4 = 0.01269
#     # a5 = -7.549e-06
#     # b5 = 0.2995
#     # c5 = 1.869e-05
#     # a6 = 0
#     # b6 = 0.299
#     # c6 = 2.22e-14
#     # a7 = 0.3389
#     # b7 = 0.2552
#     # c7 = 0.4197
#     # a8 = 0.2145
#     # b8 = 0.6476
#     # c8 = 0.3755
#
#     # a0 = 26066878.94
#     # a1 = -38086143.12
#     # b1 = -23755987.57
#     # a2 = 12507577.53
#     # b2 = 25538189.02
#     # a3 = 1319184.715
#     # b3 = -12859575.83
#     # a4 = -2453816.054
#     # b4 = 3162803.977
#     # a5 = 711390.0667
#     # b5 = -262429.0038
#     # a6 = -65069.14895
#     # b6 = -13467.81536
#     # w = 0.009610714
#
#     a0 = 62026674.3371860
#     a1= -91076668.1074805
#     b1=-55954930.4420836
#     a2=30674869.5279175
#     b2=60543571.7091132
#     a3=2535322.44707436
#     b3=-30867236.3553288
#     a4=-5702853.48805207
#     b4=7774777.00246328
#     a5=1703017.91593619
#     b5=-694109.851695349
#     a6=-160360.343129793
#     b6=-26522.7825813694
#     Omega=2.66472507858493
#
#     # if x == 0 or x == 1:
#     #     y = 1
#     # else:
#     y = a0 + a1 * cos(x*Omega) + b1 * sin(x*Omega) + a2 * cos(2 * x*Omega) + b2 * sin(2 * x*Omega) + a3 * cos(
#         3 * x*Omega) + b3 * sin(3 * x*Omega) + a4 * cos(4 * x*Omega) + b4 * sin(4 * x*Omega) + a5 * cos(
#         5 * x*Omega) + b5 * sin(5 * x*Omega) + a6 * cos(6 * x*Omega) + b6 * sin(6 * x*Omega)
#
    y = a1 * math.exp(-((x - b1) / c1) ** 2) + a2 * math.exp(-((x - b2) / c2) ** 2) + a3 * math.exp(-((x - b3) / c3) ** 2) + a4 * math.exp(-((x - b4) / c4) ** 2) + a5 * math.exp(-((x - b5) / c5) ** 2) + a6 * math.exp(-((x - b6) / c6) ** 2) + a7 * math.exp(-((x - b7) / c7) ** 2) + a8 * math.exp(-((x - b8) / c8) ** 2)
#         # y = a0 + a1*cos(X*w) + b1*sin(X*w) + a2*cos(2*X*w) + b2*sin(2*X*w) + a3*cos(3*X*w) + b3*sin(3*X*w) +  a4*cos(4*X*w) + b4*sin(4*X*w) + a5*cos(5*X*w) + b5*sin(5*X*w) + a6*cos(6*X*w) + b6*sin(6*X*w)
    return y


def width_dist(z): # inner circle radius of 0.5mm = Square_size

    a0 = w_coeffs[0]
    a1 = w_coeffs[1]
    b1 = w_coeffs[2]
    a2 = w_coeffs[3]
    b2 = w_coeffs[4]
    a3 = w_coeffs[5]
    b3 = w_coeffs[6]
    Omega = w_coeffs[7]

    # nose cone

    # a0 = 0.127296892557192
    # a1 = 0.349770299185691
    # b1=0.0708229683631302
    # a2=-0.0813172359048639
    # b2=-0.0657512907097582
    # a3=0.00726508591768930
    # b3=0.0180036913895488
    # Omega=2.36581484518248


    # cylinder
    # a0 = 0.999999999996847
    # a1 = 4.03677170631604e-12
    # b1 = 2.56746046942382e-12
    # a2 = -8.39039055547287e-13
    # b2 = -1.79383291375409e-12
    # a3 = -4.53210232501044e-14
    # b3 = 3.46694809253220e-13
    # Omega = 1.04719755119660

    # a0 =      0.2889  #(0.2889, 0.2889)
    # a1 =  -1.463e-09  #(-1.496e-09, -1.43e-09)
    # b1 =   3.019e-09  #(2.951e-09, 3.086e-09)
    # a2 =   2.055e-09  #(2.009e-09, 2.1e-09)
    # b2 =   2.568e-09  #(2.511e-09, 2.625e-09)
    # a3 =      0.1077 # (0.1077, 0.1077)
    # b3 =    -0.02749 # (-0.02749, -0.02749)
    # Omega =       4.189 # (4.189, 4.189)


    a = a0 + a1 * cos(z*Omega) + b1 * sin(z*Omega) + a2 * cos(2 * z*Omega) + b2 * sin(2 * z*Omega) + a3 * cos(
            3 * z*Omega) + b3 * sin(3 * z*Omega)
    # a = 1
    return a


a1 = list(np.linspace(0, 1, num=int(Length/Square_Size), endpoint=True))

def rounded(x):  # rounding up to the nearest Square_Size value
    a = round(x * 1 / Square_Size) / (1 / Square_Size)
    if a == 0:
        b = Square_Size
    else:
        b = round(x * 1 / Square_Size) / (1 / Square_Size)
    return b


l_half_width = []
width_dists = []
for i in range(len(a1)):
    w = width_dist(a1[i]) * Width
    a = rounded(w)
    width_dists.append(a)
    l_half_width.append(-1*a)

# vf = [0,1/4,2/4,3/4,4/4]
# random.shuffle(vf)

vf = e_coeffs
xi = []
k_composite = []
sigma_composite = []
E_composite = []
for i in range(len(vf)):
    k_composite.append(k_rigid*vf[i] + k_soft*(1-vf[i]))
    sigma_composite.append(sigma_rigid*vf[i] + sigma_soft*(1-vf[i]))
    E_composite.append(Ef * vf[i] + Em * (1 - vf[i]))
    xi.append(i*0.5/40)
k_comp = np.trapz(k_composite, xi)
sigma_comp = np.trapz(sigma_composite, xi)
E_comp = np.trapz(E_composite, xi)


# a = 0
# def vol_frac(x):
#     v = (x - Em/Ef)/(1-Em/Ef)
#     return v
# #
# vf = []
# for i in range(len(x)):
#     vf.append(vol_frac(f(x[i])))  # Volume fraction of fibre
#     # print(x[i],f(x[i]),vf[i])
#     # print(f(x[i]), Em, Ef, Em/Ef, (f(x[i]) - (Em/Ef)) / ((Ef/Ef) - (Em/Ef)))
#     a += (f(x[i]) - (Em/Ef)) / ((Ef/Ef) - (Em/Ef))
#     # print(f(x[i]))


if Square_Size > 1 or Square_Size < 0.25:
    raise ValueError('Check Square_Size')


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]
# Parameters

# Strain rate
l_measure = Length
strain_rate = 1
strain_applied = 0.001

total_time = strain_applied / strain_rate
Disp_BC = strain_applied * l_measure

origin_x = 0
origin_y = 0
origin_z = 0

A = Width / Square_Size
B = Length / Square_Size
C = Depth / Square_Size

m = mdb.models['Model-1']

# Used to sketch the lines of squares and find midpoint of the lines
# Width
W = []
w_midpoint = []
W = np.linspace(0, Width, A + 1)
w = W[1:-1]  # Select all but first and last elements in the list

# Length
L = []
l_midpoint = []
L = np.linspace(0, Length, B + 1)
l = L[1:-1]

# Depth
D = []
d_midpoint = []
D = np.linspace(0, Depth, C + 1)
d = D[1:-1]

# Sketch
m.ConstrainedSketch(name='__profile__', sheetSize=200.0)

start_point = eps*Length
# m.sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(Length, 0.0))
m.sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(0.0, width_dists[0]))
m.sketches['__profile__'].Line(point1=(Length, 0.0), point2=(Length, width_dists[-1]))
m.sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(0.0, l_half_width[0]))
m.sketches['__profile__'].Line(point1=(Length, 0.0), point2=(Length, l_half_width[-1]))

for i in range(len(a1)):
    m.sketches['__profile__'].Line(point1=(Square_Size*i, width_dists[i]), point2=(Square_Size*(i+1), width_dists[i]))
    m.sketches['__profile__'].Line(point1=(Square_Size * i, l_half_width[i]),
                                   point2=(Square_Size * (i + 1), l_half_width[i]))

for i in range(len(a1)-1):
    if width_dists[i+1] - width_dists[i] != 0:
        m.sketches['__profile__'].Line(point1=(Square_Size*(i+1), width_dists[i]),
                                       point2=(Square_Size*(i+1), width_dists[i+1]))
        m.sketches['__profile__'].Line(point1=(Square_Size*(i+1), l_half_width[i]),
                                       point2=(Square_Size*(i+1), l_half_width[i+1]))


m.Part(dimensionality=THREE_D, name='Part-1', type=DEFORMABLE_BODY)
m.parts['Part-1'].BaseSolidExtrude(depth=Depth, sketch=
m.sketches['__profile__'])

mp = m.parts['Part-1']
mpe = mp.edges
mpf = mp.faces


# Partition using Datum Plane
def datum_plane(offset, principalPlane):
    plane = mp.DatumPlaneByPrincipalPlane(offset=offset, principalPlane=principalPlane)
    all_cells = mp.cells.getByBoundingBox(origin_x, -Width, -Depth, Length, Width, Depth)
    mp.PartitionCellByDatumPlane(cells=all_cells, datumPlane=mp.datums[plane.id])


for i in range(len(l)):
    datum_plane(l[i], YZPLANE)

for i in range(len(d)):
    datum_plane(d[i], XYPLANE)

w1 = np.linspace(l_half_width[0], width_dists[0], ((width_dists[0]/Square_Size))*2+1, endpoint=True)

w2 = w1[1:-1]

for i in range(len(w2)):
    datum_plane(w2[i], XZPLANE)

# Material definition
m.Material(name='Tangoblack')
m.materials['Tangoblack'].Elastic(table=((Em, 0.3),))
m.HomogeneousSolidSection(material='Tangoblack', name='Section-1', thickness=None)

m.Material(name='Veroclear')
m.materials['Veroclear'].Elastic(table=((Ef_actual, 0.3),))
m.HomogeneousSolidSection(material='Veroclear', name='Section-2', thickness=None)

# m.Material(name='Intermediate')
# m.materials['Intermediate'].Elastic(table=((Ei, 0.3),))
# m.HomogeneousSolidSection(material='Intermediate', name='Section-3', thickness=None)


mdb.models['Model-1'].materials['Tangoblack'].Density(table=((1.145e-09,),))
mdb.models['Model-1'].materials['Veroclear'].Density(table=((1.17e-09,),))
# mdb.models['Model-1'].materials['Intermediate'].Density(table=((1.17e-09,),))

l_new = np.linspace(0.0, Length, (Length / Square_Size)+1)
l_mid = [sum(l_new[i:i + 2]) / 2 for i in range(len(l_new) - 2 + 1)]

d_new = np.linspace(0.0, Depth, (Depth/Square_Size)+1)
d_mid = [sum(d_new[i:i + 2]) / 2 for i in range(len(d_new) - 2 + 1)]

w_new = []
w_lists = []
for i in range(len(width_dists)):
    a = np.linspace(l_half_width[i], width_dists[i], ((width_dists[i]/Square_Size))*2+1)
    mid = [sum(a[i:i + 2]) / 2 for i in range(len(a) - 2 + 1)]
    w_new.append(a)
    w_lists.append(mid)

# # Creating a list of cells
cells = []
for i in range(len(l_mid)):
    sliced_cells = []
    for j in range(len(d_mid)):
        for k in range(len(w_lists[i])):
            sliced_cells.append(mp.cells.findAt(((l_mid[i], w_lists[i][k], d_mid[j]),), ))
    cells.append(sliced_cells)

# cells = []
# for i in range(len(l_mid)):
#     sliced_cells = []
#     for j in range(len(w_lists[i])):
#         for k in range(len(d_mid)):
#             sliced_cells.append(mp.cells.findAt(((l_mid[i], w_lists[i][j], d_mid[k]),), ))
#     cells.append(sliced_cells)

# randlist = []
# binarylist = []
# #
# for i in range(len(cells)):
#     l1 = len(cells[i])
#     n = random.randint(0, 2 ** l1)
#     randlist.append(n)
#     m = '{0:0{l1}b}'.format(n, l1=len(cells[i]))
#     a = [int(d) for d in str(m)]
#     binarylist.append(a)
#
# cell_list = []
# for i in range(len(cells)):
#     zipped = zip(binarylist[i], cells[i])
#     zipped.sort(key=lambda s: s[0])
#     new_bi, new = zip(*zipped)
#     cell_list.append(list(new))
cell_list = cells

# Material assignment
# 2 materials
soft = 0
rigid = 0
for q in range(len(cell_list)):

    l = len(cell_list[q])
    if pattern == 'rand-dist':

    ## pattern 1 ## random distribution

        # for sublist in cell_list:  # shuffle items in nested list
        #     random.Random(1).shuffle(sublist)
        # print('tangoblack %d cells and veroclear %d' % (l - int(round(l * vf[q])), int(round(l * vf[q]))))

        set_material_1 = cell_list[q][0:int(round(l * (1-vf[q])))]  # last item in array # tangoblack
        set_material_2 = cell_list[q][int(round(l * (1-vf[q]))):l]  # all but last item in array # veroclear

        for i in range(len(set_material_1)):  # tangoblack
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                                 FROM_SECTION)

        for i in range(len(set_material_2)):  # veroclear
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                 FROM_SECTION)
        # soft += len(set_material_1)
        # rigid += len(set_material_2)

    elif pattern == 'rand-dist2':

    ## pattern 1 ## random distribution

        for sublist in cell_list:  # shuffle items in nested list
            random.Random(2).shuffle(sublist)
        # print('tangoblack %d cells and veroclear %d' % (l - int(round(l * vf[q])), int(round(l * vf[q]))))

        set_material_1 = cell_list[q][0:int(round(l * (1-vf[q])))]  # last item in array # tangoblack
        set_material_2 = cell_list[q][int(round(l * (1-vf[q]))):l]  # all but last item in array # veroclear

        for i in range(len(set_material_1)):  # tangoblack
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                                 FROM_SECTION)

        for i in range(len(set_material_2)):  # veroclear
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                 FROM_SECTION)
        soft += len(set_material_1)
        rigid += len(set_material_2)
    elif pattern == 'clumped-top': ## causes distortion at large thickness ##

        ## pattern 2 ## tangoblack clumped at the top

        set_material_1 = cell_list[q][int(round(l * vf[q])):l]  # last item in array # tangoblack
        set_material_2 = cell_list[q][0:int(round(l * vf[q]))]  # all but last item in array # veroclear

        for i in range(len(set_material_1)):  # tangoblack
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                                 FROM_SECTION)

        for i in range(len(set_material_2)):  # veroclear
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                 FROM_SECTION)

    elif pattern == 'clumped-bot': ## causes distortion at large thickness ##

        ## pattern 3 ## tangoblack clumped at the bottom

        set_material_1 = cell_list[q][0:int(round(l * (1-vf[q])))]  # last item in array # tangoblack
        set_material_2 = cell_list[q][int(round(l * (1-vf[q]))):l]  # all but last item in array # veroclear

        for i in range(len(set_material_1)):  # tangoblack
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                                 FROM_SECTION)

        for i in range(len(set_material_2)):  # veroclear
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                 FROM_SECTION)

    elif pattern == 'half-half':

    ## pattern 4 ## tangoblack split into half; half at the top, half at the bottom

        if int(round(l * vf[q])) >= 1:
            set_material_1 = cell_list[q][0:int(round(l * (1-vf[q])/2))]
            set_material_2 = cell_list[q][int(round(l * (1-vf[q])/2)):int(round(l * (vf[q]+(1-vf[q])/2)))]
            set_material_1.extend(cell_list[q][int(round(l * (vf[q]+(1-vf[q])/2))):l])

            for i in range(len(set_material_1)):  # tangoblack
                mp.SectionAssignment(offset=0.0,
                                     offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                        cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                                     FROM_SECTION)

            for i in range(len(set_material_2)):  # veroclear
                mp.SectionAssignment(offset=0.0,
                                     offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                        cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                     FROM_SECTION)

        else:
            set_material_1 = cell_list[q][int(round(l * vf[q])):l]  # last item in array # tangoblack
            set_material_2 = cell_list[q][0:int(round(l * vf[q]))] # all but last item in array # veroclear

            for i in range(len(set_material_1)): # tangoblack
                mp.SectionAssignment(offset=0.0,
                                     offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                        cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                                     FROM_SECTION)

            for i in range(len(set_material_2)): # veroclear
                mp.SectionAssignment(offset=0.0,
                                     offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                        cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                         FROM_SECTION)
        soft += len(set_material_1)
        rigid += len(set_material_2)
    elif pattern == 'half-half-mid':

        ## pattern 5 ## veroclear split into half

        if int(round(l * vf[q])) >= 1:
            set_material_2 = cell_list[q][0:int(round(l * (vf[q])/2))]
            set_material_1 = cell_list[q][int(round(l * (vf[q])/2)):int(round(l * (vf[q]/2+(1-vf[q]))))]
            set_material_2.extend(cell_list[q][int(round(l * (vf[q]/2+(1-vf[q])))):l])
            print(len(cell_list[q]), len(set_material_1), len(set_material_2))
            for i in range(len(set_material_1)):  # tangoblack
                mp.SectionAssignment(offset=0.0,
                                     offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                        cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                                     FROM_SECTION)

            for i in range(len(set_material_2)):  # veroclear
                mp.SectionAssignment(offset=0.0,
                                     offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                        cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                     FROM_SECTION)

        else:
            set_material_1 = cell_list[q][int(round(l * vf[q])):l]  # last item in array # tangoblack
            set_material_2 = cell_list[q][0:int(round(l * vf[q]))] # all but last item in array # veroclear

            for i in range(len(set_material_1)): # tangoblack
                mp.SectionAssignment(offset=0.0,
                                     offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                        cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                                     FROM_SECTION)

            for i in range(len(set_material_2)): # veroclear
                mp.SectionAssignment(offset=0.0,
                                     offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                        cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                         FROM_SECTION)
        soft += len(set_material_1)
        rigid += len(set_material_2)
    elif pattern == 'thickness-dist':
        l = len(cell_list[q])
        set_material_1 = []
        set_material_2 = []
        N_s = round(l * (1 - vf[q]))  # number of voxels with soft phase
        layers = int(Depth / Square_Size)

        no_per_layer = int(len(cell_list[q]) / layers)
        N_s_counter = N_s

        stripped_x_y = [cell_list[q][i * no_per_layer:(i + 1) * no_per_layer] for i in range((len(cell_list[q]) + no_per_layer - 1) // no_per_layer )]
        for i in range(len(stripped_x_y)):
            if round(N_s_counter) >= len(stripped_x_y[i]):
                set_material_1.extend(stripped_x_y[i])
                N_s_counter -= len(stripped_x_y[i])
            elif round(N_s_counter) == 0:
                set_material_2.extend(stripped_x_y[i])

            elif round(N_s_counter) < len(stripped_x_y[i]):
                difference = len(stripped_x_y[i]) - N_s_counter
                a = stripped_x_y[i][int(difference//2):int(difference//2+N_s_counter)]
                set_material_1.extend(a)

                set_material_2.extend(stripped_x_y[i][0:int(difference//2)])
                set_material_2.extend(stripped_x_y[i][int(difference//2+N_s_counter):l])
                N_s_counter -= len(a)

        for i in range(len(set_material_1)):  # tangoblack
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                                 FROM_SECTION)

        for i in range(len(set_material_2)):  # veroclear
            mp.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                    cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                 FROM_SECTION)

        # if N_s <= no_per_layer:
        #     set_material_1 = []
        #     set_material_2 = []
        #     strt_idx = round(len(cell_list[q])// 2) - (N_s // 2)
        #     end_idx = round(len(cell_list[q])// 2) + (N_s // 2)
        #
        #     for idx in range(len(cell_list[q])):
        #         # checking for elements in range
        #
        #         if idx >= strt_idx and idx <= end_idx:
        #             set_material_1.append(cell_list[q][idx])
        #
        #     for element in cell_list[q]:
        #         if element not in set_material_1:
        #             set_material_2.append(element)
        #
        #     for i in range(len(set_material_1)):  # tangoblack
        #         mp.SectionAssignment(offset=0.0,
        #                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        #                 cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
        #                              FROM_SECTION)
        #
        #     for i in range(len(set_material_2)):  # veroclear
        #         mp.SectionAssignment(offset=0.0,
        #                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        #                 cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
        #                              FROM_SECTION)
        #     print(len(set_material_1), len(set_material_2), len(cell_list[q]))
        # else:
        #     no_per_layer = int(len(cell_list[q]) / layers)
        #     a = int(N_s / no_per_layer)
        #     set_material_1 = cell_list[q][0:a*no_per_layer]
        #
        #     left_over = N_s - a*no_per_layer
        #     left_over_cells = cell_list[q][a*no_per_layer:len(cell_list[q])]
        #
        #     print(cell_list[q])
        #
        #     set_material_2 = []
        #
        #     strt_idx = round(len(left_over_cells)// 2) - (N_s // 2)
        #     end_idx = round(len(left_over_cells)// 2) + (N_s // 2)
        #
        #     for idx in range(len(left_over_cells)):
        #         # checking for elements in range
        #
        #         if idx >= strt_idx and idx <= end_idx:
        #             set_material_1.append(left_over_cells[idx])
        #
        #     for element in cell_list[q]:
        #         if element not in set_material_1:
        #             set_material_2.append(element)
        #
        #     for i in range(len(set_material_1)):  # tangoblack
        #         mp.SectionAssignment(offset=0.0,
        #                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        #                 cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
        #                              FROM_SECTION)
        #
        #     for i in range(len(set_material_2)):  # veroclear
        #         mp.SectionAssignment(offset=0.0,
        #                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        #                 cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
        #                              FROM_SECTION)
        #
        #
        #
        # set_material_1 = cell_list[q][0:int(round(l * (1-vf[q])))]  # last item in array # tangoblack
        # set_material_2 = cell_list[q][int(round(l * (1-vf[q]))):l]  #
# print(soft, rigid)
# Creating a list of cells
# cells = []
# for i in range(len(d_mid)):
#     sliced_cells = []
#     for j in range(len(l_mid)):
#         for k in range(len(w_lists[j])):
#             sliced_cells.append(mp.cells.findAt(((l_mid[j], w_lists[j][k], d_mid[i]),), ))
#     cells.append(sliced_cells)


# # Creating a list of cells
# cells = []
# for i in range(len(l_mid)):
#     sliced_cells = []
#     for j in range(len(w_lists[i])):
#         for k in range(len(d_mid)):
#             sliced_cells.append(mp.cells.findAt(((l_mid[i], w_lists[i][j], d_mid[k]),), ))
#     cells.append(sliced_cells)
#
# cell_list = cells


# for q in range(len(cell_list)):

   # all but last item in array # veroclear
    # for sublist in cell_list:  # shuffle items in nested list
    #     random.Random(1).shuffle(sublist)

    # print('tangoblack %d cells and veroclear %d' % (l - int(round(l * vf[q])), int(round(l * vf[q]))))


    # for i in range(len(w_lists)):
    #     l = len(w_lists[i])
    #     sliced_cell_list = cell_list[q][i*l:(i+1)*l]
    #     set_material_1 = sliced_cell_list[0:int(round(l * (1-vf[i])))]  # last item in array # tangoblack
    #     set_material_2 = sliced_cell_list[int(round(l * (1-vf[i]))):l]  # all but last item in array # veroclear
    #     print(l)
    #
    #     for i in range(len(set_material_1)):  # tangoblack
    #         mp.SectionAssignment(offset=0.0,
    #                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    #                 cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
    #                              FROM_SECTION)
    #
    #     for i in range(len(set_material_2)):  # veroclear
    #         mp.SectionAssignment(offset=0.0,
    #                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    #                 cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
    #                              FROM_SECTION)

# 3 materials
# for q in range(len(cell_list)):
#
#     l = len(cell_list[q])
#     for sublist in cell_list:  # shuffle items in nested list
#         random.shuffle(sublist)
#     # if q <= np.where(vf == 0.0179974439188619)[0] or q > np.where(vf == 0.0304377478361769)[0]:
#
#     set_material_3 = cell_list[q][0:int(round(l * 0.2))]  # last item in array # intermediate
#     set_material_2 = cell_list[q][int(round(l * 0.2)):int(round(l * (vf[q]+0.2)))] # all but last item in array # veroclear
#     set_material_1 = cell_list[q][int(round(l *(vf[q]+0.2))):l]  # last item in array # tangoblack
#     print(len(set_material_3), len(set_material_2), len(set_material_1))
#     for i in range(len(set_material_1)): # tangoblack
#         mp.SectionAssignment(offset=0.0,
#                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
#                 cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
#                              FROM_SECTION)
#
#     for i in range(len(set_material_2)): # veroclear
#         mp.SectionAssignment(offset=0.0,
#                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
#                 cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
#                                  FROM_SECTION)
#
#     for i in range(len(set_material_3)): # intermediate
#         mp.SectionAssignment(offset=0.0,
#                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
#                 cells=set_material_3[i]), sectionName='Section-3', thicknessAssignment=
#                                  FROM_SECTION)

# Assembly
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1', part=p, dependent=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)

gap = eps*Length
mdb.models['Model-1'].rootAssembly.translate(instanceList=(
    'Part-1', ), vector=(-gap, 0.0, 0.0))

# Set definition

a.Set(faces=a.instances['Part-1'].faces.getByBoundingBox(origin_x-gap, -Width, -Depth,
                                                            origin_x-gap, Width, Depth), name='Z-face')


a.Set(faces=a.instances['Part-1'].faces.getByBoundingBox(Length-gap, -Width, -Depth,
                                                         Length-gap, Width, Depth), name='negZ-face')


a.Set(edges=a.instances['Part-1'].edges.getByBoundingBox(Length-gap, -Width, Depth,
                                                         Length-gap, Width, Depth), name='mid')

a.Set(edges=a.instances['Part-1'].edges.getByBoundingBox(origin_x-gap, origin_y, Depth,
                                                        Length-gap, origin_y, Depth), name='top')
a.Set(edges=a.instances['Part-1'].edges.getByBoundingBox(origin_x-gap, origin_y, origin_z,
                                                        Length-gap, origin_y, origin_z), name='bot')
a.Set(edges=a.instances['Part-1'].edges.getByBoundingBox(-(Length+gap), -Width, Depth,
                                                        Length-gap, Width, Depth), name='top-surf')

mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(origin_x-gap, origin_y, Depth/2))

mdb.models['Model-1'].rootAssembly.Set(name='m_Set-12', referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[10], )) #9 # 28
mdb.models['Model-1'].Coupling(controlPoint=
    mdb.models['Model-1'].rootAssembly.sets['m_Set-12'], couplingType=KINEMATIC
    , influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-1',
    surface=mdb.models['Model-1'].rootAssembly.sets['Z-face'], u1=ON, u2=ON,
    u3=ON, ur1=ON, ur2=ON, ur3=ON)
# Steps
mdb.models['Model-1'].StaticStep(initialInc=total_time, name='Step-1', previous=
'Initial')
mdb.models['Model-1'].steps['Step-1'].setValues(adaptiveDampingRatio=None,
                                                continueDampingFactors=False, initialInc=0.001,
                                                maxNumInc=10000,
                                                minInc=1e-15, nlgeom=ON,
                                                stabilizationMagnitude=1e-7,
                                                stabilizationMethod=DAMPING_FACTOR)
mdb.models['Model-1'].StaticStep(adaptiveDampingRatio=None,
                                 continueDampingFactors=False, initialInc=0.001, maxNumInc=10000,
                                 minInc=1e-15,
                                 name=
                                 'Step-2', previous='Step-1', stabilizationMagnitude=1e-7,
                                 stabilizationMethod=DAMPING_FACTOR)

mdb.models['Model-1'].StaticStep(adaptiveDampingRatio=None,
                                 continueDampingFactors=False, initialInc=0.001, maxNumInc=10000,
                                 minInc=1e-15,
                                 name=
                                 'Step-3', previous='Step-2', stabilizationMagnitude=1e-7,
                                 stabilizationMethod=DAMPING_FACTOR)

# indenter

del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0, 0.0), point1=(0.0, 10.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-3', type=
    DISCRETE_RIGID_SURFACE)
mdb.models['Model-1'].parts['Part-3'].BaseSolidExtrude(depth=10.0, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].parts['Part-3'].RemoveCells(cellList=
    mdb.models['Model-1'].parts['Part-3'].cells.getSequenceFromMask(mask=(
    '[#1 ]', ), ))

mdb.models['Model-1'].parts['Part-3'].ReferencePoint(point=
    mdb.models['Model-1'].parts['Part-3'].InterestingPoint(
    mdb.models['Model-1'].parts['Part-3'].edges[1], CENTER))

mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-3-1',
    part=mdb.models['Model-1'].parts['Part-3'])

mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-3-1', ),
    vector=(36.0, 0.0, 30.0))

mdb.models['Model-1'].rootAssembly.Set(name='Set-8', referencePoints=(
    mdb.models['Model-1'].rootAssembly.instances['Part-3-1'].referencePoints[3],
    ))

mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-3',
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'BC-1', region=mdb.models['Model-1'].rootAssembly.sets['Set-8'], u1=0.0,
    u2=0.0, u3=-15.0, ur1=0.0, ur2=0.0, ur3=0.0)
mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-3', name=
    'H-Output-2', rebar=EXCLUDE, region=
    mdb.models['Model-1'].rootAssembly.sets['Set-8'], sectionPoints=DEFAULT,
    variables=('U3', 'RF3'), frequency = 1)

mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-3', name=
'H-Output-1', rebar=EXCLUDE, region=
                                           mdb.models['Model-1'].rootAssembly.sets['m_Set-12'],
                                           sectionPoints=DEFAULT,
                                           variables=('U3', 'RF3'), frequency = 1)

mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'CDISP', 'CF', 'CSTRESS', 'LE', 'PE', 'PEEQ', 'PEMAG', 'RF', 'S', 'U',
    'STH'), frequency = 10)


tot_disp = ((1 - w_coeffs[8])+eps) * Length # e_coeffs[13] for fourier 6 and e_coeffs[24] for gauss 8 ((1 - e_coeffs[24])+0.1) * Length
# print(tot_disp, w_coeffs[8])
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                     createStepName='Step-1',
                                     distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                     name='BC-1',
                                     region=mdb.models['Model-1'].rootAssembly.sets['m_Set-12'],
                                     u1=0.1 * tot_disp, u2=0.0, u3=0.0, ur2=-0.571)  # 0.1*tot_disp


mdb.models['Model-1'].DisplacementBC(createStepName='Step-1',
                                     distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                     name='BC-4', region=mdb.models['Model-1'].rootAssembly.sets['mid'],
                                     u3=1 * Depth)

mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                     createStepName='Step-1',
                                     distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                     name='BC-5',
                                     region=mdb.models['Model-1'].rootAssembly.sets['negZ-face'],
                                     u1=0.0, u2=0.0)  #

mdb.models['Model-1'].boundaryConditions['BC-1'].setValuesInStep(stepName=
                                                                 'Step-2', u1=1 * tot_disp, u2=0.0, u3=0.0, ur2=-1.571)
mdb.models['Model-1'].boundaryConditions['BC-5'].setValuesInStep(stepName=
                                                                 'Step-2', u1=0.0, u2=0.0)  #
mdb.models['Model-1'].boundaryConditions['BC-4'].deactivate('Step-2')
mdb.models['Model-1'].rootAssembly.Set(name='Set-12', referencePoints=(
    mdb.models['Model-1'].rootAssembly.instances['Part-3-1'].referencePoints[3],
    ))
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'Step-3', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None
    , name='BC-14', region=
    mdb.models['Model-1'].rootAssembly.sets['Set-12'], u1=0.0, u2=0.0,
    u3=-15.0, ur1=0.0, ur2=0.0, ur3=0.0)

mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    allowSeparation=ON, constraintEnforcementMethod=DEFAULT,
    pressureOverclosure=HARD)
mdb.models['Model-1'].rootAssembly.Surface(name='m_Surf-1', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-3-1'].faces.getSequenceFromMask(
    ('[#4 ]', ), ))
mdb.models['Model-1'].SurfaceToSurfaceContactStd(adjustMethod=NONE,
    clearanceRegion=None, createStepName='Step-3', datumAxis=None,
    initialClearance=OMIT, interactionProperty='IntProp-1', master=
    mdb.models['Model-1'].rootAssembly.surfaces['m_Surf-1'], name='Int-1',
    slave=mdb.models['Model-1'].rootAssembly.sets['top-surf'], sliding=FINITE,
    thickness=ON)

a.seedPartInstance(deviationFactor=0.1, minSizeFactor=0.1,
                   regions=(a.instances['Part-1'],), size=Square_Size / 2)
a.generateMesh(regions=(a.instances['Part-1'],))
mdb.models['Model-1'].parts['Part-3'].seedPart(deviationFactor=0.1,
    minSizeFactor=0.1, size=2.0)
mdb.models['Model-1'].parts['Part-3'].generateMesh()

mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF,
        multiprocessingMode=DEFAULT, name='a_b_1_'+pattern, nodalOutputPrecision=FULL,
        numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
        ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

# mdb.jobs['a_b_1_'+pattern].writeInput()
# mdb.jobs['a_b_1_'+pattern].submit(consistencyChecking=ON)
# mdb.jobs['a_b_1_'+pattern].waitForCompletion()



