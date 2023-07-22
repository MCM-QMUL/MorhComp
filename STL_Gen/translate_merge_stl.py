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


def shape(name):
    m = mdb.models['Model-1']
    if name == 'rec':

        Width = 50  # mm
        aspect_ratio = 100 / 30
        Length = 15  # mm
        Depth = 2  # mm

        # Voxel size
        w = l = d = Cube_Size = 0.5  # mm


        A = Width / Cube_Size
        B = Length / Cube_Size
        C = Depth / Cube_Size

        # Rectangular graded beam coordinates
        x = np.linspace(0, (Width - w), A)  # range(0, Width, w)
        y = np.linspace(0, (Length - l), B)  # range(0, Length, l)
        z = np.linspace(0, (Depth - d), C)


        # Sketch
        m.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        m.sketches['__profile__'].rectangle(point1=(0.0, 0.0),
                                            point2=(w, l))
        m.Part(dimensionality=THREE_D, name='Part-1', type=
        DEFORMABLE_BODY)
        m.parts['Part-1'].BaseSolidExtrude(depth=d, sketch=
        m.sketches['__profile__'])

        mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)

        # translating instance
        all_slices = []
        for i in range(len(x)):
            sliced = []
            n = 0
            for j in range(len(y)):
                for k in range(len(z)):
                    n += 1
                    sliced.append(mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-%d-%d' % (i+1, n),
                        part=mdb.models['Model-1'].parts['Part-1']))

                    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-%d-%d' % (i+1, n), ),
                        vector=(x[i], y[j], z[k]))
            all_slices.append(sliced)

        # implementing linear gradient

        def f(x):

            y = x

            # if x <= 1 / 2:
            #
            #     y = 2 * x
            #
            # else:
            #
            #     y = -2 * (x - 1)
            #
            return y


        x = list(np.linspace(0, 1, num=int(Width/Cube_Size), endpoint=True))  # choosing intervals of x based on the number of slices

        vf = []  # volume fraction of each slice (If vf = 0, material = TB)
        for i in range(len(x)):
            vf.append(abs(f(x[i])))

        TB = []
        VC = []

        # assigning TB or VC
        for i in range(len(vf)):
            for sublist in all_slices:  # shuffle items in nested list
                random.shuffle(sublist)
            l = len(all_slices[i])
            if vf[i] == 0:
                TB += all_slices[i][0:l]
            elif vf[i] == 1:
                VC += all_slices[i][0:l]
            else:
                TB += all_slices[i][int(l * vf[i]):l]
                VC += all_slices[i][0:int(l * vf[i])]

        mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY,
            instances=TB[0:len(TB)],
            keepIntersections=ON, name='TB', originalInstances=DELETE)

        mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY,
            instances=VC[0:len(VC)],
            keepIntersections=ON, name='VC', originalInstances=DELETE)

        # seeding
        mdb.models['Model-1'].parts['TB'].seedPart(deviationFactor=0.1, minSizeFactor=
            0.1, size=0.5)
        mdb.models['Model-1'].parts['TB'].generateMesh()
        mdb.models['Model-1'].parts['VC'].seedPart(deviationFactor=0.1, minSizeFactor=
            0.1, size=0.5)
        mdb.models['Model-1'].parts['VC'].generateMesh()
        mdb.models['Model-1'].rootAssembly.regenerate()

        mdb.Job(activateLoadBalancing=False, atTime=None, contactPrint=OFF,
            description='', echoPrint=OFF, explicitPrecision=SINGLE,
            getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=
            PERCENTAGE, model='Model-1', modelPrint=OFF, multiprocessingMode=DEFAULT,
            name='rec_stl_%dmm' % Length, nodalOutputPrecision=SINGLE, numCpus=1, numDomains=1,
            numGPUs=0, parallelizationMethodExplicit=DOMAIN, queue=None, resultsFormat=
            ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0,
            waitMinutes=0)

        mdb.jobs['rec_stl_%dmm' % Length].writeInput()
    if name == 'flower':

        Depth = 1  # mm
        Width = 100  # mm
        Length = 100  # mm

        # Voxel size
        w = l = d = Square_Size = 0.5  # mm

        # Width distribution

        def width_dist(z):  # inner circle radius of 0.5mm = Square_size

            a0 = 0.04142
            a1 = 0.2637
            b1 = 0.00000003216
            omega = 1.571
            a = a0 + a1 * cos(z * omega) + b1 * sin(z * omega)
            return a

        x = list(np.linspace(0, 1, num=int(Length / Square_Size), endpoint=True))

        def rounded(x):  # rounding up to the nearest Square_Size value
            return round(x * 1 / Square_Size) / (1 / Square_Size)

        width_dists = []
        for i in range(len(x)):
            w = width_dist(x[i]) * Width
            a = rounded(w)
            width_dists.append(a)

        origin_x = 0
        origin_y = 0
        origin_z = 0

        B = Length / l
        C = Depth / d
        y = []
        for i in range(len(width_dists)):
            y.append(list(np.linspace(0.0, width_dists[i], (width_dists[i] / Square_Size)+1)))
        print(y)
        x = np.linspace(0, (Length - l), B)  # range(0, Length, l)
        z = np.linspace(0, (Depth - d), C)


        # Sketch
        m.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        m.sketches['__profile__'].rectangle(point1=(0.0, 0.0),
                                            point2=(0.5, 0.5))
        m.Part(dimensionality=THREE_D, name='Part-1', type=
        DEFORMABLE_BODY)
        m.parts['Part-1'].BaseSolidExtrude(depth=0.5, sketch=
        m.sketches['__profile__'])

        mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)

        # translating instance
        all_slices = []
        for i in range(len(x)):
            sliced = []
            n = 0
            for j in range(len(y[i])):
                for k in range(len(z)):
                    n += 1
                    sliced.append(mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-%d-%d' % (i+1, n),
                        part=mdb.models['Model-1'].parts['Part-1']))
                    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-%d-%d' % (i+1, n), ),
                        vector=(x[i], y[i][j], z[k]))

                    # sliced.append(mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-neg%d-neg%d' % (i+1, n),
                    #     part=mdb.models['Model-1'].parts['Part-1']))
                    #
                    # mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-neg%d-neg%d' % (i+1, n), ),
                    #     vector=(x[i], -(y[i][j]+Square_Size), z[k]))
            all_slices.append(sliced)

        # Modulus distribution

        def f(X):

            a0 = 1.604908610713824e+06
            a1 = -2.449655390777114e+06
            b1 = -1.282385958979823e+06
            a2 = 9.994086624623607e+05
            b2 = 1.441696003349646e+06
            a3 = -9.855626409633849e+04
            b3 = -7.933522785375235e+05
            a4 = -8.734760579149603e+04
            b4 = 2.325929969720600e+05
            a5 = 3.528777062039608e+04
            b5 = -3.147730362833932e+04
            a6 = -4.044782930188850e+03
            b6 = 1.010661375527444e+03
            Omega = 1.047198564199498

            if X == 0 or X == 1:
                y = 1
            else:
                y = a0 + a1 * cos(X * Omega) + b1 * sin(X * Omega) + a2 * cos(2 * X * Omega) + b2 * sin(
                    2 * X * Omega) + a3 * cos(
                    3 * X * Omega) + b3 * sin(3 * X * Omega) + a4 * cos(4 * X * Omega) + b4 * sin(
                    4 * X * Omega) + a5 * cos(
                    5 * X * Omega) + b5 * sin(5 * X * Omega) + a6 * cos(6 * X * Omega) + b6 * sin(6 * X * Omega)

            return y

        x = list(np.linspace(0, 1, num=int(Length / Square_Size), endpoint=True))

        vf = []  # vf of TB
        for i in range(len(x)):
            vf.append(f(x[i]))

        TB = []
        VC = []
        print(vf)
        # assigning TB or VC
        for i in range(len(vf)):
            for sublist in all_slices:  # shuffle items in nested list
                random.shuffle(sublist)
            l = len(all_slices[i])
            if vf[i] == 1:
                TB += all_slices[i][0:l]
            else:
                TB += all_slices[i][0:int(l * vf[i])]
                VC += all_slices[i][int(l * vf[i]):l]

        mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY,
            instances=TB[0:len(TB)],
            keepIntersections=ON, name='TB', originalInstances=DELETE)

        mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY,
            instances=VC[0:len(VC)],
            keepIntersections=ON, name='VC', originalInstances=DELETE)

        # seeding
        mdb.models['Model-1'].parts['TB'].seedPart(deviationFactor=0.1, minSizeFactor=
            0.1, size=0.5)
        mdb.models['Model-1'].parts['TB'].generateMesh()
        mdb.models['Model-1'].parts['VC'].seedPart(deviationFactor=0.1, minSizeFactor=
            0.1, size=0.5)
        mdb.models['Model-1'].parts['VC'].generateMesh()
        mdb.models['Model-1'].rootAssembly.regenerate()

        mdb.Job(activateLoadBalancing=False, atTime=None, contactPrint=OFF,
            description='', echoPrint=OFF, explicitPrecision=SINGLE,
            getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=
            PERCENTAGE, model='Model-1', modelPrint=OFF, multiprocessingMode=DEFAULT,
            name='flower_stl_%dmm' % Length, nodalOutputPrecision=SINGLE, numCpus=1, numDomains=1,
            numGPUs=0, parallelizationMethodExplicit=DOMAIN, queue=None, resultsFormat=
            ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0,
            waitMinutes=0)

        mdb.jobs['flower_stl_%dmm' % Length].writeInput()

shape(name='flower')
