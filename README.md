# MorhComp
Inverse Design Framework for Shape Morphing Composite Materials

First author: Hirak Kansara

Corresponding authors: wei.tan@qmul.ac.uk (Wei Tan), mingchao.liu@ntu.edu.sg (Mingchao Liu)


## Graphical abstract
Inverse design and additive manufacturing of shape-morphing structures based on functionally graded composites

![graphical_abstract_d2](https://github.com/MCM-QMUL/MorhComp/assets/105063133/1c6632c2-930f-4762-9ba8-5344d22db57d)

## Abstract
Shape-morphing structures possess the ability to change their shapes from one state to another, and therefore, offer great potential for a broad range of applications. A typical paradigm of morphing is transforming from an initial two-dimensional (2D) flat configuration into a three-dimensional (3D) target structure. One popular fabrication method for these structures involves programming cuts in specific locations of a thin sheet material (i.e.~kirigami), forming a desired 3D shape upon application of external mechanical load. By adopting the non-linear beam equation, an inverse design strategy has been proposed to determine the 2D cutting patterns required to achieve an axisymmetric 3D target shape. Specifically, tailoring the localised variation of bending stiffness is the key requirement. In this paper, a novel inverse design strategy is proposed by modifying the bending stiffness via introducing distributed modulus in functionally graded composites (FGCs). To fabricate the FGC-based shape-morphing structures, we use a multi-material 3D printer to print graded composites with voxel-like building blocks. The longitudinal modulus of each cross-sectional slice can be controlled through the rule of mixtures according to the micro-mechanics model, hence matching the required modulus distribution along the elastic strip. Following the proposed framework, a diverse range of structures is obtained with different Gaussian curvatures in both numerical simulations and experiments. A very good agreement is achieved between the measured shapes of morphed structures and the targets. 

In addition, the compressive rigidity and specific energy absorption during compression of FGC-based hemi-ellipsoidal morphing structures with various aspect ratios were also examined numerically and validated against experiments. By conducting systematical numerical simulations, we also demonstrate the multifunctionality of the modulus-graded shape-morphing composites. For example, they are capable of blending the distinct advantages of two different materials, i.e. one with high thermal (but low electrical) conductivity, and the other is the other way around, to achieve combined effective properties in a single structure made by FGCs. This new inverse design framework provides an opportunity to create shape-morphing structures by utilising modulus-graded composite materials, which can be employed in a variety of applications involving multi-physical environments. Furthermore, this framework underscores the versatility of the approach, enabling precise control over material properties at a local level.

## Usage of the scripts

### Morph_Gen
To generate a FEM model using the framework:
1. Use MATLAB files to generate the volume fraction and width distribution as .txt files.
2. Run 'Shell_morph.py' or 'Discretised_morph.py' to generate either solid or shell-based geometries in ABAQUS.

### STL_Gen
.stl files can be generated using the 'translate_merge_stl.py' script. To generate the axisymmetric geometry, one must fit the width and modulus distribution using a function, ideally Fourier series and vary the magnitude of coefficients accordingly.

## Reference

If you find this work useful or use this code, please cite the reference below:

Kansara, H., Liu, M. * , He, Y., Tan, W. * (2023). Inverse design and additive manufacturing of shape-morphing structures based on functionally graded composites. J Mech. Phys. Solids, 105382.

## Acknowledgement
Wei Tan acknowledges the financial support from the EPSRC New Investigator Award (grant No. EP/V049259/1). Mingchao Liu acknowledges the support from the Nanyang Technological University via the Presidential Postdoctoral Fellowship.
