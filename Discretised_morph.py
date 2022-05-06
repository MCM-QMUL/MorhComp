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
# setPath = r'C:\Temp\Morphing_Structure\Hemiellipsoid'
# os.chdir(setPath)


# Parameters - Bigger dimensions better the resolution

Depth = 1  # mm
Width = 50  # mm
Length = 50  # mm
Square_Size = 0.5  # mm size of each cube based

numSamp = Length/ Square_Size  # number of slices / number of samples
eps = 0.1 # inner circle radius
a = 1
b = 4
a_b = a/b
N = 8
Em = 0.62  # (MPa) Modulus of matrix i.e. Modulus of Tangoblack
Ef = 2033.33
Ef_actual = 0.62*4

# # Fourier series coefficients for modulus fit
# str1 = 'E_fit_eps_'+str(eps)+'_a_'+str(a)+'_b_'+str(b)+'_N_'+str(N)
# str1 = str1.replace(".", "")
# with open(str1+'.txt') as f:
#     list_of_coeffs = [line.split(',') for line in f.readlines()]
#     flat_list = [item for sublist in list_of_coeffs for item in sublist]
#
# extended = []
# for i in range(len(flat_list)):
#     extended.extend(flat_list[i].rstrip('\n').split(','))
#
# e_coeffs = [float(i) for i in extended]
#
# # Fourier series coefficients for width fit
#
# str2 ='w_fit_eps_'+str(eps)+'_a_'+str(a)+'_b_'+str(b)+'_N_'+str(N)
# str2 = str2.replace(".", "")
#
# with open(str2+'.txt') as f1:
#     list_of_coeffs1 = [line.split(',') for line in f1.readlines()]
#     flat_list1 = [item for sublist in list_of_coeffs1 for item in sublist]
#
# extended1 = []
# for i in range(len(flat_list1)):
#     extended1.extend(flat_list1[i].rstrip('\n').split(','))
#
# w_coeffs = [float(i) for i in extended1]

x = list(np.linspace(0, 1, num=int(numSamp), endpoint=True))  # choosing intervals of x based on the number of slices


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

    # a0 = 1.604908610713824e+06# (-4.803e+08, 5.305e+08)
    # a1 = -2.449655390777114e+06#(-8.216e+08, 7.453e+08)
    # b1 = -1.282385958979823e+06#(-4.012e+08, 3.602e+08)
    # a2 = 9.994086624623607e+05#(-3.268e+08, 3.571e+08)
    # b2 = 1.441696003349646e+06#(-4.121e+08, 4.58e+08)
    # a3 = -9.855626409633849e+04#(-5.482e+07, 5.256e+07)
    # b3 = -7.933522785375235e+05#(-2.594e+08, 2.344e+08)
    # a4 = -8.734760579149603e+04#(-2.008e+07, 1.701e+07)
    # b4 = 2.325929969720600e+05#(-7.294e+07, 8.014e+07)
    # a5 = 3.528777062039608e+04#(-8.979e+06, 1.015e+07)
    # b5 = -3.147730362833932e+04#(-1.208e+07, 1.114e+07)
    # a6 = -4.044782930188850e+03#(-1.275e+06, 1.145e+06)
    # b6 = 1.010661375527444e+03#(-5.443e+05, 5.682e+05)
    # Omega = 1.047198564199498#(-0.7115, 2.806)


    # a0 = -34153138.9354597# (-4.803e+08, 5.305e+08)
    # a1 = 51041737.8830233#(-8.216e+08, 7.453e+08)
    # b1 = 30323425.8629846#(-4.012e+08, 3.602e+08)
    # a2 = -18593609.9082897#(-3.268e+08, 3.571e+08)
    # b2 = -33906611.1633922#(-4.121e+08, 4.58e+08)
    # a3 = -481970.361745478#(-5.482e+07, 5.256e+07)
    # b3 = 18376129.208518#(-2.594e+08, 2.344e+08)
    # a4 = 3143967.17881367#(-2.008e+07, 1.701e+07)
    # b4 = -5160121.0297764#(-7.294e+07, 8.014e+07)
    # a5 = -1074323.60537138#(-8.979e+06, 1.015e+07)
    # b5 = 600827.522657738#(-1.208e+07, 1.114e+07)
    # a6 = 117357.752246658#(-1.275e+06, 1.145e+06)
    # b6 = -369.150568791272#(-5.443e+05, 5.682e+05)
    # Omega = 0.0479559072908973#(-0.7115, 2.806)

    a1 =       0.502
    b1 =      0.3141
    c1 =      0.2848
    a2 =    -0.01209
    b2 =      0.2338
    c2 =     0.09678
    a3 =           0
    b3 =      -9.614
    c3 =     0.00206
    a4 =   5.169e-05
    b4 =      0.4352
    c4 =     0.02899
    a5 =   -0.008395
    b5 =      0.1414
    c5 =     0.08515
    a6 =      0.5235
    b6 =    -0.02723
    c6 =       1.764
    a7 =     0.02358
    b7 =     -0.1152
    c7 =      0.1433
    a8 =        2293
    b8 =       2.205
    c8 =      0.4034



    if x == 0 or x == 1 or x == 0.5:
        y = 1
    else:
        # y = a0 + a1 * cos(X*Omega) + b1 * sin(X*Omega) + a2 * cos(2 * X*Omega) + b2 * sin(2 * X*Omega) + a3 * cos(
        #     3 * X*Omega) + b3 * sin(3 * X*Omega) + a4 * cos(4 * X*Omega) + b4 * sin(4 * X*Omega) + a5 * cos(
        #     5 * X*Omega) + b5 * sin(5 * X*Omega) + a6 * cos(6 * X*Omega) + b6 * sin(6 * X*Omega)

        y = a1 * exp(-((x - b1) / c1) ** 2) + a2 * exp(-((x - b2) / c2) ** 2) + a3 * exp(-((x - b3) / c3) ** 2) + a4 * exp(-((x - b4) / c4) ** 2) + a5 * exp(-((x - b5) / c5) ** 2) + a6 * exp(-((x - b6) / c6) ** 2) + a7 * exp(-((x - b7) / c7) ** 2) + a8 * exp(-((x - b8) / c8) ** 2)

    return y



def width_dist(z): # inner circle radius of 0.5mm = Square_size

    # a0 = w_coeffs[0]
    # a1 = w_coeffs[1]
    # b1 = w_coeffs[2]
    # a2 = w_coeffs[3]
    # b2 = w_coeffs[4]
    # a3 =  w_coeffs[5]
    # b3 = w_coeffs[6]
    # Omega = w_coeffs[7]

    a0 = 0.102463492672831
    a1 = 0.154740740968973
    b1 = -0.0146699253012995
    a2 = -0.0447558216578581
    b2 = 0.00360311732963929
    a3 =  0.00960274049900937
    b3 = 0.00150556989143305
    Omega = 2.0943951023932


    a = a0 + a1 * cos(z*Omega) + b1 * sin(z*Omega) + a2 * cos(2 * z*Omega) + b2 * sin(2 * z*Omega) + a3 * cos(
            3 * z*Omega) + b3 * sin(3 * z*Omega)

    return a


a1 = list(np.linspace(0, 1, num=int(Length/Square_Size), endpoint=True))

def rounded(x):  # rounding up to the nearest Square_Size value
    return round(x * 1/Square_Size) / (1/Square_Size)


l_half_width = []
width_dists = []
for i in range(len(a1)):
    w = width_dist(a1[i]) * Width
    a = rounded(w)
    width_dists.append(a)
    l_half_width.append(-1*a)

Ec = []
for i in range(len(x)):
    Ec.append(f(x[i]) * Ef)
    # print(x[i], f(x[i]), f(x[i]) * Ef)
vf = []
a = 0
for i in range(len(Ec)):
    vf.append((Ec[i] - Em) / (Ef - Em))  # Volume fraction of fibre
    a += (Ec[i] - Em) / (Ef - Em)

# print(a/len(vf))
if Square_Size > 1 or Square_Size < 0.25:
    raise ValueError('Check Square_Size')

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

# m.sketches['__profile__'].Line(point1=(Length, 0.0), point2=(Length+1.5, 0.0))
# m.sketches['__profile__'].Line(point1=(Length, width_dists[-1]), point2=(Length+1.5, width_dists[-1]))
# m.sketches['__profile__'].Line(point1=(Length+1.5, 0.0), point2=(Length+1.5, width_dists[-1]))

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

m.Material(name='Rigid')
m.materials['Rigid'].Elastic(table=((Ef_actual*4, 0.3),))
m.HomogeneousSolidSection(material='Rigid', name='Section-3', thickness=None)

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

# for i in range(len(w_lists)):
#     print(w_lists[i], width_dists[i])

# Creating a list of cells
cells = []
for i in range(len(l_mid)):
    sliced_cells = []
    for j in range(len(w_lists[i])):
        for k in range(len(d_mid)):
            sliced_cells.append(mp.cells.findAt(((l_mid[i], w_lists[i][j], d_mid[k]),), ))
    cells.append(sliced_cells)
#

randlist = []
binarylist = []

for i in range(len(cells)):
    l1 = len(cells[i])
    n = random.randint(0, 2 ** l1)
    randlist.append(n)
    m = '{0:0{l1}b}'.format(n, l1=len(cells[i]))
    a = [int(d) for d in str(m)]
    binarylist.append(a)

cell_list = []
for i in range(len(cells)):
    zipped = zip(binarylist[i], cells[i])
    zipped.sort(key=lambda s: s[0])
    new_bi, new = zip(*zipped)
    cell_list.append(list(new))


# Material assignment

for q in range(len(cell_list)):

    l = len(cell_list[q])
    for sublist in cell_list:  # shuffle items in nested list
        random.shuffle(sublist)

    set_material_1 = cell_list[q][int(round(l * vf[q])):l]  # last item in array
    set_material_2 = cell_list[q][0:int(round(l * vf[q]))] # all but last item in array

    for i in range(len(set_material_1)):
        mp.SectionAssignment(offset=0.0,
                             offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                cells=set_material_1[i]), sectionName='Section-1', thicknessAssignment=
                             FROM_SECTION)

    for i in range(len(set_material_2)):
        mp.SectionAssignment(offset=0.0,
                             offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
                                 FROM_SECTION)
    # elif vf[q] == 1:
    #     set_material_2 = cell_list[q][0:l]
    #
    #     for i in range(len(set_material_2)):
    #         mp.SectionAssignment(offset=0.0,
    #                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    #                 cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
    #                              FROM_SECTION)
    #
    # else:
    #     set_material_3 = cell_list[q][int(round(l * (1-vf[q]))):l]
    #     set_material_2 = cell_list[q][0:int(round(l * 1-(1-vf[q])))]
    #
    #     for i in range(len(set_material_2)):
    #         mp.SectionAssignment(offset=0.0,
    #                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    #                 cells=set_material_2[i]), sectionName='Section-2', thicknessAssignment=
    #                              FROM_SECTION)
    #
    #     for i in range(len(set_material_3)):
    #         mp.SectionAssignment(offset=0.0,
    #                              offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    #                 cells=set_material_3[i]), sectionName='Section-3', thicknessAssignment=
    #                              FROM_SECTION)

# Assembly
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1', part=p, dependent=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)

# Set definition

a.Set(faces=a.instances['Part-1'].faces.getByBoundingBox(origin_x, -Width, -Depth,
                                                            origin_x, Width, Depth), name='Z-face')


a.Set(faces=a.instances['Part-1'].faces.getByBoundingBox(Length, -Width, -Depth,
                                                         Length, Width, Depth), name='negZ-face')

# a.Set(faces=a.instances['Part-1'].faces.getByBoundingBox(origin_x, origin_y, origin_z,
#                                                          Length, origin_y, Depth), name='sym')

# a.Set(cells=a.instances['Part-1'].cells.getByBoundingBox(Length, origin_y, origin_z,
#                                                          Length, Width, Depth), name='stiff')

# a.Set(edges=a.instances['Part-1'].edges.getByBoundingBox(Width/2, origin_y, Depth,
#                                                          Width/2, Length, Depth), name='mid')

# minval = min(vf)
# min_index = vf.index(minval) + 1

a.Set(edges=a.instances['Part-1'].edges.getByBoundingBox(Length, -Width, Depth,
                                                         Length, Width, Depth), name='mid')

a.Set(edges=a.instances['Part-1'].edges.getByBoundingBox(origin_x, origin_y, Depth,
                                                        Length, origin_y, Depth), name='top')
a.Set(edges=a.instances['Part-1'].edges.getByBoundingBox(origin_x, origin_y, origin_z,
                                                        Length, origin_y, origin_z), name='bot')
## Rigid part

# a.Set(faces=a.instances['Part-1'].faces.getByBoundingBox(Length+1.5, origin_y, origin_z,
#                                                          Length+1.5, Width, Depth), name='negZ-face')
#
# a.Set(faces=a.instances['Part-1'].faces.getByBoundingBox(origin_x, origin_y, origin_z,
#                                                          Length+1.5, origin_y, Depth), name='sym')
#
# a.Set(cells=a.instances['Part-1'].cells.getByBoundingBox(Length+1.5, origin_y, origin_z,
#                                                          Length+1.5, Width, Depth), name='stiff')


mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(origin_x, origin_y, Depth/2))

# mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(Length, width_dists[-1]/2, Depth/2))


mdb.models['Model-1'].rootAssembly.Set(name='m_Set-12', referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[9], )) #28
mdb.models['Model-1'].Coupling(controlPoint=
    mdb.models['Model-1'].rootAssembly.sets['m_Set-12'], couplingType=KINEMATIC
    , influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-1',
    surface=mdb.models['Model-1'].rootAssembly.sets['Z-face'], u1=ON, u2=ON,
    u3=ON, ur1=ON, ur2=ON, ur3=ON)




# Buckling method

def buckling_method(method):
    tot_disp = (1 - 0.535814) * Length # (1 - 0.73662) * Length

    if method == 'imperfection':

        mdb.models['Model-1'].BuckleStep(name='Step-1', numEigen=3, previous='Initial', vectors=6)

        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                             createStepName='Step-1',
                                             distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                             name='BC-1',
                                             region=mdb.models['Model-1'].rootAssembly.sets['outer-edge'],
                                             u1=1 * tot_disp, u2=0, u3=0, ur2=1.571)

        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                             createStepName='Step-1',
                                             distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                             name='BC-2',
                                             region=mdb.models['Model-1'].rootAssembly.sets['sym'],
                                             u2=0.0)

        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                             createStepName='Step-1',
                                             distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                             name='BC-3',
                                             region=mdb.models['Model-1'].rootAssembly.sets['negZ-face'],
                                             u1=-Square_Size*2, u2=0.0)

        a.seedPartInstance(deviationFactor=0.1, minSizeFactor=0.1,
                           regions=(a.instances['Part-1'],), size=Square_Size/2)
        a.generateMesh(regions=(a.instances['Part-1'],))

        mdb.models['Model-1'].keywordBlock.synchVersions(storeNodesAndElements=False)
        mdb.models['Model-1'].keywordBlock.insert(200, '\n*node file\nu,')

        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
                explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
                memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF,
                multiprocessingMode=DEFAULT, name='RVE_slice_flower', nodalOutputPrecision=FULL,
                numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
                ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

        mdb.jobs['RVE_slice_flower'].writeInput()
        # mdb.jobs['RVE_flower'].submit(consistencyChecking=OFF)
        # mdb.jobs['RVE_flower'].waitForCompletion()

        mdb.Model(name='Model-2', objectToCopy=mdb.models['Model-1'])
        mdb.models['Model-2'].StaticStep(adaptiveDampingRatio=None,
                                         continueDampingFactors=False, initialInc=0.0001, maintainAttributes=True,
                                         maxNumInc=5000, minInc=1e-15, name='Step-1', nlgeom=ON, previous='Initial',
                                         stabilizationMagnitude=1e-7, stabilizationMethod=DAMPING_FACTOR)

        mdb.models['Model-2'].keywordBlock.synchVersions(storeNodesAndElements=False)
        mdb.models['Model-2'].keywordBlock.insert(280, '\n*imperfection, file=RVE_slice_flower, step=1\n1, 5.0e-2\n2, 0.0e-2\n3, 0.0e-2')

        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
                explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
                memory=90, memoryUnits=PERCENTAGE, model='Model-2', modelPrint=OFF,
                multiprocessingMode=DEFAULT, name='RVE_slice_flower_Buckled', nodalOutputPrecision=FULL,
                numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
                ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

        mdb.jobs['RVE_slice_flower_Buckled'].writeInput()
        # mdb.jobs['RVE_slice_flower_Buckled'].submit(consistencyChecking=OFF)
        # mdb.jobs['RVE_slice_flower_Buckled'].waitForCompletion()

    if method == '2step_rotation':

        mdb.models['Model-1'].StaticStep(initialInc=total_time, name='Step-1', previous=
        'Initial')
        mdb.models['Model-1'].steps['Step-1'].setValues(adaptiveDampingRatio=None,
                                                        continueDampingFactors=False, initialInc=0.001, maxNumInc=10000,
                                                        minInc=1e-15, nlgeom=ON,
                                                        stabilizationMagnitude=1e-7, stabilizationMethod=DAMPING_FACTOR)
        mdb.models['Model-1'].StaticStep(adaptiveDampingRatio=None,
                                         continueDampingFactors=False, initialInc=0.001, maxNumInc=10000, minInc=1e-15,
                                         name=
                                         'Step-2', previous='Step-1', stabilizationMagnitude=1e-7,
                                         stabilizationMethod=DAMPING_FACTOR)


        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                             createStepName='Step-1',
                                             distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                             name='BC-1',
                                             region=mdb.models['Model-1'].rootAssembly.sets['m_Set-12'],
                                             u1=0.1 * tot_disp, u2=0.0, u3=0.0, ur2=0.571)  # 0.1*tot_disp

        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                             createStepName='Step-1',
                                             distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                             name='BC-2',
                                             region=mdb.models['Model-1'].rootAssembly.sets['sym'],
                                             u2=0.0)

        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                             createStepName='Step-1',
                                             distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                             name='BC-3',
                                             region=mdb.models['Model-1'].rootAssembly.sets['negZ-face'],
                                             u1=0.0, u2=0.0)

        mdb.models['Model-1'].boundaryConditions['BC-1'].setValuesInStep(stepName=
                                                                         'Step-2', u1=tot_disp, u2=0.0, u3=0.0, ur2=1.571)

        a.seedPartInstance(deviationFactor=0.1, minSizeFactor=0.1,
                           regions=(a.instances['Part-1'],), size=Square_Size/2)
        a.generateMesh(regions=(a.instances['Part-1'],))

        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
                explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
                memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF,
                multiprocessingMode=DEFAULT, name='RVE_flower_2step_rotation_50mm', nodalOutputPrecision=FULL,
                numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
                ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

        mdb.jobs['RVE_flower_2step_rotation_50mm'].writeInput()
        mdb.jobs['RVE_flower_2step_rotation_50mm'].submit(consistencyChecking=OFF)
        mdb.jobs['RVE_flower_2step_rotation_50mm'].waitForCompletion()


    if method == '2step_out_of_plane':

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

            mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                                 createStepName='Step-1',
                                                 distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                                 name='BC-1',
                                                 region=mdb.models['Model-1'].rootAssembly.sets['m_Set-12'],
                                                 u1=0.1 * tot_disp, u2=0.0, u3=0.0, ur2=-0.571)  # 0.1*tot_disp

            mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-2', name=
            'H-Output-2', rebar=EXCLUDE, region=
                                                       mdb.models['Model-1'].rootAssembly.sets['top'],
                                                       sectionPoints=DEFAULT,
                                                       variables=('U3','U1'))
            mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-2', name=
            'H-Output-3', rebar=EXCLUDE, region=
                                                       mdb.models['Model-1'].rootAssembly.sets['bot'],
                                                       sectionPoints=DEFAULT,
                                                       variables=('U3','U1'))

            # mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
            #                                      createStepName='Step-1',
            #                                      distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
            #                                      name='BC-3',
            #                                      region=mdb.models['Model-1'].rootAssembly.sets['sym'],
            #                                      u2=0.0)

            mdb.models['Model-1'].DisplacementBC(createStepName='Step-1',
                                                 distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                                 name='BC-4', region=mdb.models['Model-1'].rootAssembly.sets['mid'],
                                                 u3=1 * Depth)

            mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING,
                                                 createStepName='Step-1',
                                                 distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                                 name='BC-5',
                                                 region=mdb.models['Model-1'].rootAssembly.sets['negZ-face'],
                                                 u1=-0.1*0.1*Length, u2=0.0)  #

            mdb.models['Model-1'].boundaryConditions['BC-1'].setValuesInStep(stepName=
                                                                             'Step-2', u1=1 * tot_disp, u2=0.0, u3=0.0, ur2=-1.571)
            mdb.models['Model-1'].boundaryConditions['BC-5'].setValuesInStep(stepName=
                                                                             'Step-2', u1=-0.1*Length, u2=0.0)  #
            mdb.models['Model-1'].boundaryConditions['BC-4'].deactivate('Step-2')

            a.seedPartInstance(deviationFactor=0.1, minSizeFactor=0.1,
                               regions=(a.instances['Part-1'],), size=Square_Size / 2)
            a.generateMesh(regions=(a.instances['Part-1'],))

            mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
                    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
                    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF,
                    multiprocessingMode=DEFAULT, name='a_1_b_4_eps_01_50mm_2mat', nodalOutputPrecision=FULL,
                    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
                    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

            # mdb.jobs['a_1_b_4_eps_01_50mm_2mat'].writeInput()
            # mdb.jobs['a_1_b_4_eps_01_50mm_2mat'].submit(consistencyChecking=ON)
            # mdb.jobs['a_1_b_4_eps_01_50mm_2mat'].waitForCompletion()

buckling_method(method='2step_out_of_plane')

# # Model 2
#
# mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-2')
# mdb.models['Model-2'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
# mdb.models['Model-2'].sketches['__profile__'].rectangle(point1=(-15.0, 15.0),
#     point2=(15.0, -15.0))
# mdb.models['Model-2'].Part(dimensionality=THREE_D, name='Part-1', type=
#     DISCRETE_RIGID_SURFACE)
# mdb.models['Model-2'].parts['Part-1'].BaseShell(sketch=
#     mdb.models['Model-2'].sketches['__profile__'])
# del mdb.models['Model-2'].sketches['__profile__']
# mdb.models['Model-2'].Part(name='Part-2', objectToCopy=
#     mdb.models['Model-2'].parts['Part-1'])
# mdb.models['Model-2'].PartFromOdb(frame=52, instance='PART-1', name='GC', odb=
#     session.openOdb(r'C:/Temp/Morphing_Structure/Hemiellipsoid/2step-test.odb')
#     , shape=DEFORMED, step=1)
# mdb.models['Model-2'].Material(name='Material-1')
# mdb.models['Model-2'].materials['Material-1'].Elastic(table=((1000.0, 0.3), ))
# mdb.models['Model-2'].materials['Material-1'].elastic.setValues(table=((2.48,
#     0.3), ))
# mdb.models['Model-2'].HomogeneousSolidSection(material='Material-1', name=
#     'Section-1', thickness=None)
# mdb.models['Model-2'].parts['GC'].Set(elements=
#     mdb.models['Model-2'].parts['GC'].elements.getSequenceFromMask(mask=(
#     '[#ffffffff:268 ]', ), ), name='Set-1073')
# mdb.models['Model-2'].parts['GC'].SectionAssignment(offset=0.0, offsetField='',
#     offsetType=MIDDLE_SURFACE, region=
#     mdb.models['Model-2'].parts['GC'].sets['Set-1073'], sectionName='Section-1'
#     , thicknessAssignment=FROM_SECTION)
# mdb.models['Model-2'].rootAssembly.DatumCsysByDefault(CARTESIAN)
# mdb.models['Model-2'].rootAssembly.Instance(dependent=ON, name='GC-1', part=
#     mdb.models['Model-2'].parts['GC'])
# mdb.models['Model-2'].rootAssembly.Instance(dependent=ON, name='Part-1-1',
#     part=mdb.models['Model-2'].parts['Part-1'])
# mdb.models['Model-2'].rootAssembly.Instance(dependent=ON, name='Part-2-1',
#     part=mdb.models['Model-2'].parts['Part-2'])
# mdb.models['Model-2'].rootAssembly.translate(instanceList=('Part-2-1', ),
#     vector=(15.0, 0.0, 0.0))
# mdb.models['Model-2'].rootAssembly.translate(instanceList=('Part-2-1', ),
#     vector=(-15.0, 0.0, 0.0))
# mdb.models['Model-2'].rootAssembly.translate(instanceList=('GC-1', ), vector=(
#     -22.0, 0.0, 0.0))
# mdb.models['Model-2'].rootAssembly.translate(instanceList=('GC-1', ), vector=(
#     2.0, 0.0, 0.0))
# mdb.models['Model-2'].rootAssembly.RadialInstancePattern(axis=(0.0, 0.0, 1.0),
#     instanceList=('GC-1', ), number=8, point=(0.0, 0.0, 0.0), totalAngle=360.0)
# mdb.models['Model-2'].rootAssembly._previewMergeMeshes(instances=(
#     mdb.models['Model-2'].rootAssembly.instances['GC-1'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-2'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-3'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-4'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-5'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-6'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-7'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-8']),
#     keepOnlyOrphanElems=True, nodeMergingTolerance=1e-06)
# mdb.models['Model-2'].rootAssembly.InstanceFromBooleanMerge(domain=BOTH,
#     instances=(mdb.models['Model-2'].rootAssembly.instances['GC-1'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-2'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-3'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-4'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-5'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-6'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-7'],
#     mdb.models['Model-2'].rootAssembly.instances['GC-1-rad-8']),
#     keepIntersections=ON, mergeNodes=BOUNDARY_ONLY, name='Part-3',
#     nodeMergingTolerance=1e-06, originalInstances=SUPPRESS)
# mdb.models['Model-2'].rootAssembly.translate(instanceList=('Part-2-1', ),
#     vector=(0.0, 0.0, 17.0))
# mdb.models['Model-2'].StaticStep(adaptiveDampingRatio=None,
#     continueDampingFactors=False, maxNumInc=1000, name='Step-1', nlgeom=ON,
#     previous='Initial', stabilizationMagnitude=0.0002, stabilizationMethod=
#     DAMPING_FACTOR)






