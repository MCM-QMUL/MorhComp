# Matlab code to generate discretised morphing structure
# This code was developed by Hirak Kansara 
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

# Local directory
setPath = r'C:\Temp\Morphing_Structure\MorhComp\Morph_Gen\new_morphing_structure'
os.chdir(setPath)

# HPC directory
# cwd = os.getcwd()
# os.chdir(cwd)

a_tot = 1
b_tot = 1

N = 8
eps = 0.1
Depth = 1  # mm
Width = 100  # mm
Length = 100  # mm
Square_Size = 4 # mm size of each cube based

numSamp = Length/ Square_Size  # number of slices / number of samples
 # inner circle radius


Em = 3.14  # (MPa) Modulus of matrix i.e. Modulus of Tangoblack
Ef = 12.6 # 2033.33
Ef_actual = Ef


def read_coeffs(a, b):
# Parameters - Bigger dimensions better the resolution



    # Fourier series coefficients for modulus fit
    str1 = 'T_s' #'E_fit_eps_'+str(eps)+'_a_'+str(a)+'_b_'+str(b)+'_N_'+str(N)
    str1 = str1.replace(".", "")
    with open(str1+'.txt') as f: # with open(str1+'.txt') as f:
        list_of_coeffs = [line.split(',') for line in f.readlines()]
        flat_list = [item for sublist in list_of_coeffs for item in sublist]

    extended = []
    for i in range(len(flat_list)):
        extended.extend(flat_list[i].rstrip('\n').split(','))

    e_coeffs = [float(i) for i in extended]

    # Fourier series coefficients for width fit

    str2 ='w_fit_eps_'+str(eps)+'_a_'+str(a)+'_b_'+str(b)+'_N_'+str(N)
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

e_coeffs, w_coeffs = read_coeffs(a_tot, b_tot)

mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-1')

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



    if x == 0 or x == 1:
        y = 1
    else:
        # y = a0 + a1 * cos(x*Omega) + b1 * sin(x*Omega) + a2 * cos(2 * x*Omega) + b2 * sin(2 * x*Omega) + a3 * cos(
        #     3 * x*Omega) + b3 * sin(3 * x*Omega) + a4 * cos(4 * x*Omega) + b4 * sin(4 * x*Omega) + a5 * cos(
        #     5 * x*Omega) + b5 * sin(5 * x*Omega) + a6 * cos(6 * x*Omega) + b6 * sin(6 * x*Omega)

        y = a1 * math.exp(-((x - b1) / c1) ** 2) + a2 * math.exp(-((x - b2) / c2) ** 2) + a3 * math.exp(-((x - b3) / c3) ** 2) + a4 * math.exp(-((x - b4) / c4) ** 2) + a5 * math.exp(-((x - b5) / c5) ** 2) + a6 * math.exp(-((x - b6) / c6) ** 2) + a7 * math.exp(-((x - b7) / c7) ** 2) + a8 * math.exp(-((x - b8) / c8) ** 2)

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

    a = a0 + a1 * cos(z*Omega) + b1 * sin(z*Omega) + a2 * cos(2 * z*Omega) + b2 * sin(2 * z*Omega) + a3 * cos(
            3 * z*Omega) + b3 * sin(3 * z*Omega)

    return a


a1 = list(np.linspace(0, 1, num=int(Length/Square_Size), endpoint=True))

width_dists = []
l_half_width = []

for i in range(len(a1)):
    w = width_dist(a1[i]) * Width
    width_dists.append(w)

width_dists = sorted(width_dists)
for i in range(len(width_dists)):
    l_half_width.append(-1 * width_dists[i])

start_point = eps*Length

y = np.linspace(start_point,Length+start_point, int(numSamp))
y_list = y.tolist()
np.array(width_dists, dtype=np.float64)
np.array(l_half_width, dtype=np.float64)
np.array(y_list, dtype=np.float64)


u_half_coords = zip(width_dists, y_list)  # coordinates of upper half of petals
l_half_coords = zip(l_half_width, y_list)  # coordinates of lower half of petals

origin_x = 0
origin_y = 0
origin_z = 0

A = Width / Square_Size
B = Length / Square_Size
C = Depth / Square_Size

pi = math.pi

def transform(x,y,i):

    x_trans = x*cos(2*pi*(i/N)) - y*sin(2*pi*(i/N))
    y_trans = x * sin(2 * pi * (i / N)) + y * cos(2 * pi * (i / N))

    return x_trans, y_trans


x_u_half = []
y_u_half = []
x_l_half = []
y_l_half = []
for i in range(N):
    for j in range(len(y_list)):

        a, b = transform(y_list[j], width_dists[j], i)
        x_u_half.append(a)
        y_u_half.append(b)
        c, d = transform(y_list[j], l_half_width[j], i)
        x_l_half.append(c)
        y_l_half.append(d)

m = mdb.models['Model-1']


# Sketch
m.ConstrainedSketch(name='__profile__', sheetSize=200.0)

n = int(numSamp)
x_u = [x_u_half[i:i + n] for i in range(0, len(x_u_half), n)]
y_u = [y_u_half[i:i + n] for i in range(0, len(y_u_half), n)]
x_l = [x_l_half[i:i + n] for i in range(0, len(x_l_half), n)]
y_l = [y_l_half[i:i + n] for i in range(0, len(y_l_half), n)]


for i in range(-1, N-1):
    y_u[i][0] = y_l[i + 1][0]
    x_u[i][0] = x_l[i + 1][0]

for j in range(N):

    for i in range(len(y_list) - 1):

        m.sketches['__profile__'].Line(point1=(y_u[j][i], x_u[j][i]),
                                       point2=(y_u[j][i+1], x_u[j][i+1]))

        if i % Length/numSamp == 0:
            m.sketches['__profile__'].Line(point1=(y_u[j][(200*(1*i))-1], x_u[j][(200*(1*i))-1]),
                                           point2=(y_l[j][(200*(1*i))-1], x_l[j][(200*(1*i))-1]))

    for i in range(len(y_list) - 1):
        m.sketches['__profile__'].Line(point1=(y_l[j][i], x_l[j][i]),
                                       point2=(y_l[j][i+1], x_l[j][i+1]))

mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseShell(sketch=mdb.models['Model-1'].sketches['__profile__'])


mdb.models['Model-1'].parts['Part-1'].Set(faces=
    mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#1 ]',
    ), ), name='Set-1')

for i in range(-1,N-1):
    mdb.models['Model-1'].parts['Part-1'].PartitionFaceByShortestPath(faces=mdb.models['Model-1'].parts['Part-1'].faces.findAt(((0, 0, 0),), ),
                                                    point1=(y_u[i][0],x_u[i][0],0), point2=(y_u[i+1][0],x_u[i+1][0],0))


mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Elastic(table=((0.1 * Em, 0.3, 0.1),
                                                             (Em, 0.3, 1.0), (2.0 * Em, 0.3, 2.0),
                                                             (3.0 * Em, 0.3, 3.0)),
                                                      temperatureDependency=ON)
mdb.models['Model-1'].materials['Material-1'].Density(table=((1.145e-09,),))

mdb.models['Model-1'].Material(name='Material-2')
mdb.models['Model-1'].materials['Material-2'].Elastic(table=((10000, 0.3),))

mdb.models['Model-1'].materials['Material-2'].Density(table=((1.145e-09,),))

mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION,
                                              integrationRule=SIMPSON, material='Material-1',
                                              name='Section-1',
                                              nodalThicknessField='', numIntPts=5,
                                              poissonDefinition=DEFAULT,
                                              preIntegrate=OFF, temperature=GRADIENT, thickness=2.0,
                                              thicknessField='',
                                              thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)

mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION,
                                              integrationRule=SIMPSON, material='Material-2',
                                              name='Section-2',
                                              nodalThicknessField='', numIntPts=5,
                                              poissonDefinition=DEFAULT,
                                              preIntegrate=OFF, temperature=GRADIENT, thickness=2.0,
                                              thicknessField='',
                                              thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)


mdb.models['Model-1'].parts['Part-1'].Set(faces=
    mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#ff ]',
    ), ), name='petals')
mdb.models['Model-1'].parts['Part-1'].Set(faces=
    mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#100 ]',
    ), ), name='Set-3')

mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0,
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-1'].sets['petals'], sectionName=
    'Section-1', thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0,
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-1'].sets['Set-3'], sectionName=
    'Section-2', thicknessAssignment=FROM_SECTION)

mdb.models['Model-1'].StaticStep(adaptiveDampingRatio=None, continueDampingFactors=False, initialInc=0.01,
                                 maxNumInc=5000, minInc=
                                 1e-15, name='Step-1', nlgeom=ON, previous='Initial',
                                 stabilizationMagnitude=1E-7, stabilizationMethod=DAMPING_FACTOR, timePeriod=1.0)
mdb.models['Model-1'].StaticStep(adaptiveDampingRatio=None, continueDampingFactors=False, initialInc=0.01,
                                 maxNumInc=5000, minInc=
                                 1e-15, name='Step-2', previous='Step-1', stabilizationMagnitude=1E-7,
                                 stabilizationMethod=DAMPING_FACTOR, timePeriod=1.0)

mdb.models['Model-1'].ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP,
                                           application=QUASI_STATIC, initialConditions=OFF, initialInc=0.001, maxNumInc=1000,
                                           minInc=
                                           1e-10, name='Step-3', nohaf=OFF, previous='Step-2',
                                           timePeriod=1.0)

# mdb.models['Model-1'].ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP,
#                                            application=QUASI_STATIC, initialConditions=OFF, maxNumInc=1000,
#                                            minInc=
#                                            1e-10, name='Step-4', nohaf=OFF, previous='Step-3',
#                                            timePeriod=1.0)


mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1',
                                            part=mdb.models['Model-1'].parts['Part-1'])

# mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
#     sheetSize=200.0)
# mdb.models['Model-1'].sketches['__profile__'].rectangle(
#     point1=(-100.0, -100.0), point2=(100.0, 100.0))
# mdb.models['Model-1'].Part(dimensionality=THREE_D, name=
#     'Part-2', type=DISCRETE_RIGID_SURFACE)
# mdb.models['Model-1'].parts['Part-2'].BaseShell(sketch=
#     mdb.models['Model-1'].sketches['__profile__'])
# mdb.models['Model-1'].parts['Part-2'].ReferencePoint(point=
#     mdb.models['Model-1'].parts['Part-2'].InterestingPoint(
#     mdb.models['Model-1'].parts['Part-2'].edges[3], MIDDLE))
# mdb.models['Model-1'].rootAssembly.Instance(dependent=ON,
#     name='Part-2-1', part=
#     mdb.models['Model-1'].parts['Part-2'])
# mdb.models['Model-1'].rootAssembly.Instance(dependent=ON,
#     name='Part-2-2', part=
#     mdb.models['Model-1'].parts['Part-2'])


del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0, 0.0), point1=(0.0, 25.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-3', type=
    DISCRETE_RIGID_SURFACE)
mdb.models['Model-1'].parts['Part-3'].BaseSolidExtrude(depth=20.0, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].parts['Part-3'].RemoveCells(cellList=
    mdb.models['Model-1'].parts['Part-3'].cells.getSequenceFromMask(mask=(
    '[#1 ]', ), ))

mdb.models['Model-1'].parts['Part-3'].ReferencePoint(point=
    mdb.models['Model-1'].parts['Part-3'].InterestingPoint(
    mdb.models['Model-1'].parts['Part-3'].edges[1], CENTER))

height = e_coeffs[25]*Length+5

mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-3-1',
    part=mdb.models['Model-1'].parts['Part-3'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-3-1', ),
    vector=(0.0, 0.0, height))

mdb.models['Model-1'].rootAssembly.Set(name='top', referencePoints=(
    mdb.models['Model-1'].rootAssembly.instances['Part-3-1'].referencePoints[3],
    ))
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-4', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-3-1'].faces.getSequenceFromMask(
    ('[#5 ]', ), ))

mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-3', name=
    'H-Output-2', rebar=EXCLUDE, region=
    mdb.models['Model-1'].rootAssembly.sets['top'], sectionPoints=DEFAULT,
    variables=('RF3', 'U3'), timeInterval=0.01)

# mdb.models['Model-1'].rootAssembly.Set(name='bot', referencePoints=(
#     mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].referencePoints[2],
# ))


for i in range(N):
    x = start_point + Length
    y = width_dists[-1]/2
    y_neg = l_half_width[-1]/2

    x_t = x*cos(2*pi*(i/N)) - y*sin(2*pi*(i/N))
    y_t = x*sin(2*pi*(i/N)) + y*cos(2*pi*(i/N))

    x_neg_t = x*cos(2*pi*(i/N)) - y_neg*sin(2*pi*(i/N))
    y_neg_t = x*sin(2*pi*(i/N)) + y_neg*cos(2*pi*(i/N))

    x_orig = x
    y_orig = 0

    x_orig_t = x_orig*cos(2*pi*((i+1)/N)) - y_orig*sin(2*pi*((i+1)/N))
    y_orig_t = x_orig*sin(2*pi*((i+1)/N)) + y_orig*cos(2*pi*((i+1)/N))

    x_perp = x_orig*cos(2*pi*((1)/N)) - y_orig*sin(2*pi*((1)/N))
    y_perp = x_orig*sin(2*pi*((1)/N)) + y_orig*cos(2*pi*((1)/N))

    x_perp_t = x_perp*cos(2*pi*((i+1)/N)) - y_perp*sin(2*pi*((i+1)/N))
    y_perp_t = x_perp*sin(2*pi*((i+1)/N)) + y_perp*cos(2*pi*((i+1)/N))

    x_face = start_point + y_list[-10]
    y_face = width_dists[-1]/2
    y_neg_face = l_half_width[-1]/2

    x_face_t = x_face*cos(2*pi*(i/N)) - y_face*sin(2*pi*(i/N))
    y_face_t = x_face*sin(2*pi*(i/N)) + y_face*cos(2*pi*(i/N))

    x_face_neg_t = x_face*cos(2*pi*(i/N)) - y_neg_face*sin(2*pi*(i/N))
    y_face_neg_t = x_face*sin(2*pi*(i/N)) + y_neg_face*cos(2*pi*(i/N))


    mdb.models['Model-1'].rootAssembly.DatumCsysByThreePoints(coordSysType=CARTESIAN,
                                                              name='Datum csys-%d' % (i + 2),
                                                              origin=(0.0, 0.0, 0.0),
                                                              point1=(x_orig_t, y_orig_t, 0),
                                                              point2=(x_perp_t, y_perp_t, 0))

    mdb.models['Model-1'].rootAssembly.Set(edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt(((x_t,
                   y_t, 0),), ((x_neg_t, y_neg_t, 0),),), name='Edge-%d' % (i+1))

    mdb.models['Model-1'].rootAssembly.Set(faces=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(((x_face_t,
                   y_face_t, 0),), ((x_face_neg_t, y_face_neg_t, 0),),), name='P-%d' % (i+1))


mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    allowSeparation=ON, constraintEnforcementMethod=DEFAULT,
    pressureOverclosure=HARD)
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-1', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
    ('[#1ff ]', ), ))
mdb.models['Model-1'].SurfaceToSurfaceContactStd(adjustMethod=NONE,
                                                 clearanceRegion=None, createStepName='Step-3',
                                                 datumAxis=None,
                                                 initialClearance=OMIT, interactionProperty='IntProp-1',
                                                 master=
                                                 mdb.models['Model-1'].rootAssembly.surfaces['Surf-4'],
                                                 name='Int-1',
                                                 slave=mdb.models['Model-1'].rootAssembly.surfaces[
                                                     'Surf-1'], sliding=FINITE
                                                 , thickness=ON)

# mdb.models['Model-1'].rootAssembly.Set(edges=
#     mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.getSequenceFromMask(
#     ('[#100000 #10000000 #0 #10 #1000 #100000 #10000000', ' #0 #10 #1000 ]'), )
#     , name='edges')
# mdb.models['Model-1'].rootAssembly.Surface(name='m_Surf-3', side1Faces=
# mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].faces.getSequenceFromMask(
#     ('[#1 ]',), ))
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(((0, 0, 0),), ), name='mid')
# mdb.models['Model-1'].SurfaceToSurfaceContactStd(adjustMethod=NONE,
#                                                  clearanceRegion=None, createStepName='Step-3',
#                                                  datumAxis=None, initialClearance=OMIT, enforcement=NODE_TO_SURFACE,
#                                                  interactionProperty='IntProp-1',
#                                                  master=mdb.models['Model-1'].rootAssembly.surfaces[
#                                                      'm_Surf-3'], name=
#                                                  'Int-2',
#                                                  slave=mdb.models['Model-1'].rootAssembly.sets['edges'],
#                                                  sliding=
#                                                  FINITE, smooth=0.2, surfaceSmoothing=NONE, thickness=ON)

# mdb.models['Model-1'].ContactStd(createStepName='Initial', name='Int-3')
# mdb.models['Model-1'].interactions['Int-3'].includedPairs.setValuesInStep(
#     addPairs=((mdb.models['Model-1'].rootAssembly.surfaces['Surf-1'], SELF), ),
#     stepName='Initial', useAllstar=OFF)
# mdb.models['Model-1'].interactions['Int-3'].contactPropertyAssignments.appendInStep(
#     assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Initial')

# mdb.models['Model-1'].StdInitialization(name='CInit-1')
#
# mdb.models['Model-1'].interactions['Int-3'].initializationAssignments.appendInStep(
#     assignments=((mdb.models['Model-1'].rootAssembly.surfaces['Surf-1'], SELF,
#     'CInit-1'), ), stepName='Initial')

tot_disp = ((1 - e_coeffs[24]) + 0.1) * Length

mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM,
                                     fieldName='', fixed=OFF, localCsys=None, name=
                                     'BC-1', region=mdb.models['Model-1'].rootAssembly.sets['Edge-1'],
                                     u1=-0.1*tot_disp, u2=0.0, u3=0.0, ur1=0.0, ur2=0.571, ur3=0.0)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM,
                                     fieldName='', fixed=OFF, localCsys=mdb.models['Model-1'].rootAssembly.datums[8], name=
                                     'BC-2', region=mdb.models['Model-1'].rootAssembly.sets['Edge-2'],
                                     u1=-0.1*tot_disp, u2=0.0, u3=0.0, ur1=0.0, ur2=0.571, ur3=0.0)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM,
                                     fieldName='', fixed=OFF, localCsys=mdb.models['Model-1'].rootAssembly.datums[11], name=
                                     'BC-3', region=mdb.models['Model-1'].rootAssembly.sets['Edge-3'],
                                     u1=-0.1*tot_disp, u2=0.0, u3=0.0, ur1=0.0, ur2=0.571, ur3=0.0)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM,
                                     fieldName='', fixed=OFF, localCsys=mdb.models['Model-1'].rootAssembly.datums[14], name=
                                     'BC-4', region=mdb.models['Model-1'].rootAssembly.sets['Edge-4'],
                                     u1=-0.1*tot_disp, u2=0.0, u3=0.0, ur1=0.0, ur2=0.571, ur3=0.0)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM,
                                     fieldName='', fixed=OFF, localCsys=mdb.models['Model-1'].rootAssembly.datums[17], name=
                                     'BC-5', region=mdb.models['Model-1'].rootAssembly.sets['Edge-5'],
                                     u1=-0.1*tot_disp, u2=0.0, u3=0.0, ur1=0.0, ur2=0.571, ur3=0.0)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM,
                                     fieldName='', fixed=OFF, localCsys=mdb.models['Model-1'].rootAssembly.datums[20], name=
                                     'BC-6', region=mdb.models['Model-1'].rootAssembly.sets['Edge-6'],
                                     u1=-0.1*tot_disp, u2=0.0, u3=0.0, ur1=0.0, ur2=0.571, ur3=0.0)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM,
                                     fieldName='', fixed=OFF, localCsys=mdb.models['Model-1'].rootAssembly.datums[23], name=
                                     'BC-7', region=mdb.models['Model-1'].rootAssembly.sets['Edge-7'],
                                     u1=-0.1*tot_disp, u2=0.0, u3=0.0, ur1=0.0, ur2=0.571, ur3=0.0)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM,
                                     fieldName='', fixed=OFF, localCsys=mdb.models['Model-1'].rootAssembly.datums[26], name=
                                     'BC-8', region=mdb.models['Model-1'].rootAssembly.sets['Edge-8'],
                                     u1=-0.1*tot_disp, u2=0.0, u3=0.0, ur1=0.0, ur2=0.571, ur3=0.0)

mdb.models['Model-1'].boundaryConditions['BC-1'].setValuesInStep(stepName='Step-2', u1=-tot_disp, ur2=1.571)
mdb.models['Model-1'].boundaryConditions['BC-2'].setValuesInStep(stepName='Step-2', u1=-tot_disp, ur2=1.571)
mdb.models['Model-1'].boundaryConditions['BC-3'].setValuesInStep(stepName='Step-2', u1=-tot_disp, ur2=1.571)
mdb.models['Model-1'].boundaryConditions['BC-4'].setValuesInStep(stepName='Step-2', u1=-tot_disp, ur2=1.571)
mdb.models['Model-1'].boundaryConditions['BC-5'].setValuesInStep(stepName='Step-2', u1=-tot_disp, ur2=1.571)
mdb.models['Model-1'].boundaryConditions['BC-6'].setValuesInStep(stepName='Step-2', u1=-tot_disp, ur2=1.571)
mdb.models['Model-1'].boundaryConditions['BC-7'].setValuesInStep(stepName='Step-2', u1=-tot_disp, ur2=1.571)
mdb.models['Model-1'].boundaryConditions['BC-8'].setValuesInStep(stepName='Step-2', u1=-tot_disp, ur2=1.571)

# mdb.models['Model-1'].EncastreBC(createStepName='Step-3', localCsys=None, name=
# 'BC-9', region=mdb.models['Model-1'].rootAssembly.sets['bot'])
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-3',
                                     distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None,
                                     name=
                                     'BC-10', region=mdb.models['Model-1'].rootAssembly.sets['top'], u1=0.0,
                                     u2=
                                     0.0, u3=-height*0.3, ur1=0.0, ur2=0.0, ur3=0.0)

a1 = e_coeffs[0]
b1 = e_coeffs[1]
c1 = e_coeffs[2]
a2 = e_coeffs[3]
b2 = e_coeffs[4]
c2 = e_coeffs[5]
a3 = e_coeffs[6]
b3 = e_coeffs[7]
c3 = e_coeffs[8]
a4 = e_coeffs[9]
b4 = e_coeffs[10]
c4 = e_coeffs[11]
a5 = e_coeffs[12]
b5 = e_coeffs[13]
c5 = e_coeffs[14]
a6 = e_coeffs[15]
b6 = e_coeffs[16]
c6 = e_coeffs[17]
a7 = e_coeffs[18]
b7 = e_coeffs[19]
c7 = e_coeffs[20]
a8 = e_coeffs[21]
b8 = e_coeffs[22]
c8 = e_coeffs[23]

mdb.models['Model-1'].ExpressionField(description='',
                                      expression='a1 * math.exp(-((X - b1) / c1) ** 2) + a2 * math.exp(-((X - b2) / c2) ** 2) + a3 * math.exp(-((X - b3) / c3) ** 2) + a4 * math.exp(-((X - b4) / c4) ** 2) + a5 * math.exp(-((X - b5) / c5) ** 2) + a6 * math.exp(-((X - b6) / c6) ** 2) + a7 * math.exp(-((X - b7) / c7) ** 2) + a8 * math.exp(-((X - b8) / c8) ** 2)',
                                      localCsys=mdb.models['Model-1'].rootAssembly.datums[1],
                                      name='AnalyticalField-1')
mdb.models['Model-1'].Temperature(createStepName='Initial',
                                  crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
                                  distributionType=FIELD
                                  , field='AnalyticalField-1', magnitudes=(1.0,), name='Predefined Field-1',
                                  region=mdb.models['Model-1'].rootAssembly.sets['P-1'])
for i in range(N-1):
    mdb.models['Model-1'].ExpressionField(description='',
                                          expression='a1 * math.exp(-((X - b1) / c1) ** 2) + a2 * math.exp(-((X - b2) / c2) ** 2) + a3 * math.exp(-((X - b3) / c3) ** 2) + a4 * math.exp(-((X - b4) / c4) ** 2) + a5 * math.exp(-((X - b5) / c5) ** 2) + a6 * math.exp(-((X - b6) / c6) ** 2) + a7 * math.exp(-((X - b7) / c7) ** 2) + a8 * math.exp(-((X - b8) / c8) ** 2)',
                                          localCsys=mdb.models['Model-1'].rootAssembly.datums[8+3*i],
                                          name='AnalyticalField-%d' % (i+2))
    mdb.models['Model-1'].Temperature(createStepName='Initial',
                                      crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
                                      distributionType=FIELD
                                      , field='AnalyticalField-%d' % (i+2), magnitudes=(1.0,),
                                      name='Predefined Field-%d' % (i+2),
                                      region=mdb.models['Model-1'].rootAssembly.sets['P-%d' % (i+2)])
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(((0, 0, 0),), ), name='mid')
mdb.models['Model-1'].Temperature(createStepName='Initial',
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(1.0, ), name='Predefined Field-9', region=
    mdb.models['Model-1'].rootAssembly.sets['mid'])


mdb.models['Model-1'].parts['Part-1'].seedPart(deviationFactor=0.1,
                                               minSizeFactor=0.1, size=5.0)
mdb.models['Model-1'].parts['Part-1'].setMeshControls(regions=
    mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#ff ]',
    ), ), technique=SWEEP)
mdb.models['Model-1'].parts['Part-1'].generateMesh()


mdb.models['Model-1'].parts['Part-3'].seedPart(deviationFactor=0.1,
    minSizeFactor=0.1, size=4.2)
mdb.models['Model-1'].parts['Part-3'].generateMesh()


# mdb.models['Model-1'].parts['Part-1'].PartitionEdgeByParam(edges=
#     mdb.models['Model-1'].parts['Part-1'].edges.getSequenceFromMask((
#     '[#1 #40000 #0 #10 #400000 #0 #100', ' #4000000 #0 #1000 #40000000 ]'), ),
#     parameter=0.5)
# mdb.models['Model-1'].parts['Part-1'].PartitionFaceByShortestPath(faces=
#     mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#100 ]',
#     ), ), point1=mdb.models['Model-1'].parts['Part-1'].vertices[51], point2=
#     mdb.models['Model-1'].parts['Part-1'].vertices[251])
# mdb.models['Model-1'].parts['Part-1'].PartitionFaceByShortestPath(faces=
#     mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#201 ]',
#     ), ), point1=mdb.models['Model-1'].parts['Part-1'].vertices[5], point2=
#     mdb.models['Model-1'].parts['Part-1'].vertices[156])
# mdb.models['Model-1'].parts['Part-1'].PartitionFaceByShortestPath(faces=
#     mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#3 ]',
#     ), ), point1=mdb.models['Model-1'].parts['Part-1'].vertices[4], point2=
#     mdb.models['Model-1'].parts['Part-1'].vertices[9])
# mdb.models['Model-1'].parts['Part-1'].PartitionFaceByShortestPath(faces=
#     mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask((
#     '[#2010 ]', ), ), point1=mdb.models['Model-1'].parts['Part-1'].vertices[12]
#     , point2=mdb.models['Model-1'].parts['Part-1'].vertices[207])
# mdb.models['Model-1'].parts['Part-1'].setMeshControls(regions=
#     mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask((
#     '[#807f ]', ), ), technique=SWEEP)

mdb.models['Model-1'].parts['Part-1'].generateMesh()

mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF,
    multiprocessingMode=DEFAULT, name='a_'+str(a_tot)+'_b_'+str(b_tot)+'_std', nodalOutputPrecision=FULL,
    numCpus=8, numDomains=8, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

#
# mdb.jobs['a_'+str(a_tot)+'_b_'+str(b_tot)].writeInput()
# mdb.jobs['a_'+str(a_tot)+'_b_'+str(b_tot)].submit(consistencyChecking=ON)
# mdb.jobs['a_'+str(a_tot)+'_b_'+str(b_tot)].waitForCompletion()
