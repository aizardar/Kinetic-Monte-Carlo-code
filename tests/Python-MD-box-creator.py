#!/usr/bin/python2.7
#print '!!!!! Under Modification  !!! Feb 2014'

# MD atomic box creator for TiAl alloys. TiAl bulk, Ti3Al bulk, bi-crystals of G/G or A2/G and also multi lamellar samples are possible.  Orientation analysor and Burgers vector length calculator for  any crystal structure.
# It is also possible to apply a translational shift at mid layer of bulk structures or at interfaces of bi-crystals or lamellar samples.
# This code needs text input files with suitable format including the parameters of crystal structure and their orientation. For each crystal orientation of lamellar structure there must be a separate input file.
# Written by Mansour Kanani; ICAMS; Ruhr-University Bochum; Feb. 2014
# Edited by Aram Asaadi; ICAMS; Ruhr-University Bochum; Sep. 2015

import ase
from ase import *
from numpy import *
import os
from ase.visualize import *
from ase.utils.geometry import *  
import math
import shutil

X=array([1,0,0]);Y=array([0,1,0]);Z=array([0,0,1]) 
set_printoptions(precision=6)
set_printoptions(suppress=True)
Tol=0.1      # for interplanar distance calculation

SLAB_DIC={'G':'./box-TiAl-EAM.in','G_001':'./box-TiAl-EAM-001.in', 'PT':'./box-PT-EAM.in', 'RB':'./box-RB-EAM.in', 'TT':'./box-TT-EAM.in', 'A2':'./box-Ti3Al-EAM.in', 'A2_Gbased':'./box-Ti3AlGbased-EAM.in'}

SF_DIC={'CSF_P':[1/6.,[1,-2,1]],'CSF_N':[1/6.,[-2,1,1]],'SISF_P':[1/6.,[1,1,-2]],'APB_P':[1/2.,[1,0,-1]],'APB_N':[1/2.,[-1,0,1]],'OD_P':[1/2.,[-1,1,0]],'SF_none':[0,[0,0,0]],'A2_SISF_P':[1/4.,[-1,-1,0]] }

# ^^ to read the parameter file including orientation and structure information
def Initial_Cell_Reader(ParamName,):
    with open(ParamName) as f_in:
        data = filter(None, (line.rstrip() for line in f_in)) # to store data line by line without empty lines
    StruName=data[0].split()[0]     # structure name
    p1=float(data[1].split()[0])    # lattice constant
    p2=matrix([[float(data[i].split()[0].split(',')[j]) for j in range(0,3)] for i in range(2,5)]) # unit vectors
    p3=int(data[5].split()[0])      # number of given atoms in cell definition
    p4=array([[float(data[i].split()[0].split(',')[j]) for j in range(0,3)] for i in range(6,6+p3)]) # positions of atoms in scaled coordinate
    p5=data[6+p3].split()[0]      # suitable atomic symbol according to ASE standard
    p6=[float(data[6+p3+1].split()[0].split(',')[j]) for j in range(0,3)] # First crystal direction
    p7=[float(data[6+p3+2].split()[0].split(',')[j]) for j in range(0,3)] # Second crystal direction
    p8=[float(data[6+p3+3].split()[0].split(',')[j]) for j in range(0,3)] # Third crystal direction
    u=p1*p2
    n1n2n3=[p6,p7,p8]
    return u,p3,p4,p5,n1n2n3

def Miller_to_Cart(A,Vi):
    UA=A.cell
    V=array(Vi*matrix(UA))[0]
    return V
# ^^ to calculate interplanar distances along desired direction in miller format but in the inputed structure
def Inter_Planar_dist(A,Vi,Tol):
    A=A.repeat(30)
    NLay, dLay= get_layers(A,Vi, tolerance=Tol)
    d=[]
    for i in range(0,(len(dLay)-2)): # if the direction is to tilted try to avoid final layers i.e. increase the 2 here
        d.append(dLay[i+1]-dLay[i])
    D=average(array(d),axis=0)
    Nr=len(d)                   # if requested can be extracted
    return D

#^^ finding a nearest atomic position to the general center of mass and pick that as specific center of mass CM
def Center_Atom(A,ASymb):
    A=A.copy()
    cm=A.get_center_of_mass()
    A.append(Atom('Cu',cm))
    N0=A.get_number_of_atoms()
    if ASymb != 'none':
        I0=[[atom.index,A.get_distance(N0-1,atom.index)] for atom in A if (N0-1 != atom.index) and atom.symbol==ASymb]
    else:
        I0=[[atom.index,A.get_distance(N0-1,atom.index)] for atom in A if (N0-1 != atom.index)]
    I0.sort(key=lambda k: (k[1], k[0]), reverse=False)
    CM=A[I0[0][0]].position
    CMi=A[I0[0][0]].index
    return CM,CMi

#^^ orienting the input structure into a requested directions in parameter files
def Cell_Alaigner(S,V_in):
    CM,CMi=Center_Atom(S,'Ti') # finding a nearest atomic position to the general center of mass and pick that as specific center of mass CM
    # find the Cartesian of the vectors in initial stage and also the interplanar distances along them
    VF0=[]
    for v in V_in:
        VF0.append(Miller_to_Cart(S,v))

    S.rotate(VF0[0],'x',center=CM,rotate_cell=True) # first orientation step: rotating first direction along x

    # calculate the cartesian of parameter direction in the oriented system
    VF1=[]
    for v in V_in:
        VF1.append(Miller_to_Cart(S,v))
    

    phi=math.acos(dot(VF1[1],Y)/(linalg.norm(VF1[1])*linalg.norm(Y))) # the angle between second direction with y after first rotation
    if cross(VF1[1],Y)[0]/X[0] < 0: # check that it needs positive or negative rotation
        phi=-phi

    S.rotate('x',phi,center=CM,rotate_cell=True) # second orientation step: rotating second axis around x into y
    return S

#^^ to cut a 3d box from an Atmos object; Note: y boundary is z-dependent
def SLAB_CUT(slab,L1,L2,L3,yfact):
    del slab[[atom.index for atom in slab if atom.x < L1[0] or atom.x > L1[1]]]
    del slab[[atom.index for atom in slab if (atom.y < L2[0]+(yfact*atom.z)) or (atom.y > L2[1]+(yfact*atom.z))]]
    del slab[[atom.index for atom in slab if atom.z < L3[0] or atom.z > L3[1]]]
    return(slab)

#^^ to create a periodic suppercell with desired number of layers in 3d; ** the number of layers which creates periodicity should be determined manually; try and error
def Supercell_Creator(A, v_in_c, vd, TSh):
    print 'min layers to maintaine periodicity for Gamma [4,12,3], Gamma_001 [2,2,2] and for alpha2 [4,8,2]'
    #NL_Sup=eval(raw_input('Insert number of layers for supercell in x,y and z e.g. [2,2,2] :\n'))
    if param == 'A2' or param == 'A2_Gbased':  NL_Sup=[4,8,2]
    elif param == 'G_001':  NL_Sup=[2,2,2]
    else: NL_Sup=[4,12,3]         
    pbc_Sup=[1,1,1]
    v_in_c=array(v_in_c)
    v_in_unit=[v/linalg.norm(v) for v in v_in_c]

    A=A.repeat(10)
    #view(A)
    CM,CMi=Center_Atom(A,'Ti')
    O=[0,0,0]-A[CMi].position
    A.translate([O[0],O[1],O[2]])
    if TSh !='none': A.translate(TSh)
    v_fi_c=[(NL_Sup[i]-1+pbc_Sup[i])*vd[i]*v_in_unit[i] for i in range(0,3)]

    print "The supercell will have ", NL_Sup,"layers in each directions \nwith this equation:\n",matrix(v_fi_c)
    #print v_in_c , v_in_unit 
    Y_factor=v_fi_c[2][1]/v_fi_c[2][2]
    XLimit=[-vd[0]/2.,(vd[0]/2.)+v_fi_c[0][0]-(pbc_Sup[0]*vd[0])]
    YLimit=[-vd[1]/2.,(vd[1]/2.)+v_fi_c[1][1]-(pbc_Sup[1]*vd[1])]
    ZLimit=[-vd[2]/2.,(vd[2]/2.)+v_fi_c[2][2]-(pbc_Sup[2]*vd[2])]
    #print "===>>",Y_factor, XLimit,YLimit,ZLimit,vd
    Super_Cell=SLAB_CUT(A,XLimit,YLimit,ZLimit,Y_factor)
    Super_Cell.set_cell(v_fi_c)
    
    if raw_input("See the supercell? (y/n) :  ")=='y':
        view(Super_Cell)

    return Super_Cell

#^^ once the supercell has created, it could be transformed in a large box with periodic or non-periodic directions
def Box_Creator(S, vd):
    R_mode='y'
    while R_mode=='y':
        NL_Box=eval(raw_input("Insert desired number of layers along x,y,z for a slab or box: \n"))
        #NL_Box=[8,24,6]
        #pbc_Box=eval(raw_input("Periodicity of the new box? e.g. [1,1,0] : \n"))
        pbc_Box=[1,1,1]
        U_Sup=S.cell
        Uu_Sup=[u/linalg.norm(u) for u in U_Sup]
        v_fi_c=[(NL_Box[i]-1+pbc_Box[i])*vd[i]*Uu_Sup[i] for i in range(0,3)]
        Rep_Box=[int(math.ceil(linalg.norm(v_fi_c[i])/linalg.norm(U_Sup[i]))) for i in range(3)]
        Box=S.repeat(Rep_Box)

        Y_factor=v_fi_c[2][1]/v_fi_c[2][2]
        XLimit=[-vd[0]/2.,(vd[0]/2.)+v_fi_c[0][0]-(pbc_Box[0]*vd[0])]
        YLimit=[-vd[1]/2.,(vd[1]/2.)+v_fi_c[1][1]-(pbc_Box[1]*vd[1])]
        ZLimit=[-vd[2]/2.,(vd[2]/2.)+v_fi_c[2][2]-(pbc_Box[2]*vd[2])]
        Box=SLAB_CUT(Box,XLimit,YLimit,ZLimit,Y_factor)
        Box.set_cell(v_fi_c)

        print 70*"."
        print "The new box is created with this number of atoms  :  ", Box.get_number_of_atoms()
        print "and with following box equation: \n", matrix(Box.cell)
        print "and these number of layers in each direction : ",  NL_Box
        print 70*"."
        if raw_input("See the box? (y/n) :  ")=='y':
            view(Box)
        R_mode=raw_input("Another box? (y/n) :  ")

    return Box

#^^ to introduce a SF on the created box
def SF_Introducer(atoms,S,SF_name): 
    print 'the cell of the current box is:\n',atoms.cell
    Z_interface=eval(raw_input('insert the z component of interface location to implement shifts over that :\n'))
    Factor=SF_DIC[SF_name][0]
    Shift=SF_DIC[SF_name][1]
    VShift=Factor*Miller_to_Cart(S,Shift)
    for a in atoms:
        if a.z > Z_interface:
            atoms[a.index].position=a.position+VShift
    #view(atoms)
    print SF_name,'stacking fault along the', SF_DIC[SF_name],'direction equal to cartesian of',VShift,'has been implemented in the requsted structure'
    return atoms


def Tuning_shift(temp_box,V_z,dd,AllaignedStru):
    print '-'*70, '\n=====>>  Stacking tuning step'
    YF=V_z[1]/V_z[2]
    XL=[-8.5*dd[0],8.5*dd[0]]
    YL=[-dd[1],24.5*dd[1]]
    ZL=[V_z[2]-1.5*dd[2],V_z[2]+0.5*dd[2]]
    temp_box=SLAB_CUT(temp_box,XL,YL,ZL,YF)
    view(temp_box)

    MoreShift='y'
    while MoreShift=='y':
        print 70*'.'
        factor,V=map(str,raw_input("Insert the pre-factor and your desired direction in Miller indices format \n e.g. 1./6. [1,-2,1]    Don't forget to put . for fractional factors: \n").split( ))
        V=eval(V)
        Factor=eval(factor)
        VShift=Factor*Miller_to_Cart(AllaignedStru,V)
        L=linalg.norm(VShift)
        D=Inter_Planar_dist(AllaignedStru,V,Tol)
        print "your input vector is:              ", factor,V
        print "in new coordinate system it is:     ", VShift
        print "with length of:                     ", L
        print "interplanar ditances along that is: ", D

        atoms=temp_box.copy()
        for a in atoms:
            if a.z > V_z[2]-0.5*dd[2]:
                atoms[a.index].position=a.position+VShift
        view(atoms)

        MoreShift=raw_input('more shift or accept this shift or deny shifting? (y/a/d) : ')

    return MoreShift, VShift

#^^ the main function to create a box as one bulk or on grain or one lamella including box size and SF 
def Slab_Creator(ParameterPath, RefVD, RefSuperCell, RefSlabCell, TuneShift):
    print '^'*70
    print '=====>> Working structure is : ' , param
    U0,Nr0,P0_scaled,Symb0,V_initial = Initial_Cell_Reader(ParameterPath)
    Stru=Atoms(Symb0,scaled_positions=P0_scaled,cell=U0) # Atomic object definition according to input data
    #view(Stru)

    VD=[]
    for v in V_initial:
        VD.append(Inter_Planar_dist(Stru,v,Tol))
    print Stru.cell
    #****  Orienting in new system 
    print '-'*70, '\n=====>>  Cell orienting step'
    STRU=Cell_Alaigner(Stru,V_initial)
    # check the orientation and length of the main vectors after orientation
    VF=[];VFl=[]
    for v in V_initial:
        v_cartesian=Miller_to_Cart(STRU,v)
        VF.append(v_cartesian)
        VFl.append(linalg.norm(v_cartesian))
    #view(STRU)

    #****   Supercell creation 
    print '-'*70, '\n=====>>  Supercell creation step'
    SUPERCell=Supercell_Creator(STRU, VF, VD, TuneShift)
    #view(SUPERCell)
    if RefSuperCell != 'base':
        if raw_input('adjusting to the reference slab? (y/n) : ')=='y':
            print 'Sepercell before addoption: \n', SUPERCell.cell
            SUPERCell.set_cell(RefSuperCell,scale_atoms=True)
            VD=RefVD
            print 'Sepercell After addoption: \n', SUPERCell.cell

    #****  Create larger boxes with repeatation 
    print '-'*70, '\n=====>>  Box creation step'
    if RefSlabCell != 'base' :
        x_nl=RefSlabCell[0][0]/VD[0]
        y_nl=RefSlabCell[1][1]/VD[1]
        print 'to have consistency with previous slab take:\n x_Nr_Layer= ',x_nl,'and y_Nr_layer= ',y_nl

    BOX_STRU=Box_Creator(SUPERCell, VD)

    #***** Introduce stacking fault shifts 
    #print '-'*70, '\n=====>>  SF introducing step'
    # SF_NAME=raw_input('insert desired SF direction:  ')
    # if SF_NAME in SF_DIC:
    #     BOX_STRU=SF_Introducer(BOX_STRU,STRU,SF_NAME)
    # else:
    #     print 'Inserted SF name was not defined therefore it hasnt been applied'

    return BOX_STRU,VD, STRU, SUPERCell   #STRU can change to Stru if directions in original structure is desired

#*********** Cutting the model in Y direction **************
 
def SLAB_CUT_Y(slab,L2,yfact):
    del slab[[atom.index for atom in slab if (atom.y < L2[0]+(yfact*atom.z)) or (atom.y > L2[1]+(yfact*atom.z))]]
    return(slab)

######################################################################################
#                                     Main Body                                      #
######################################################################################
print 'H'*30,'Starting','H'*30
print 'H'*70,'\n'

###  INPUT 
# input data must be provided in given path including required data as template; e.g. ./S3Al.in
param=raw_input('insert BASE grain structure (G,G_001,PT,RB,TT,A2,A2_Gbased) : ')

### Creating first slab as reference
REF_SLAB,REF_VD,REF_A_STRU,REF_A_SUPERCELL=Slab_Creator(SLAB_DIC[param],'base','base','base','none')

### Adding more slab on the reference
FINAL_BOX=REF_SLAB.copy()
print  '\n', 'ss'*15,'Additional slabs  ','ss'*15
print '='*72
MoreSlab=raw_input('\n more slabe to add ? (y/n) : ')
while MoreSlab == 'y':
    param=raw_input('insert ADDITIONAL slab structure (G,PT,RB,TT,A2,A2_Gbased) : ')
    
    Tuning='t';Tune_Shift='none'
    while Tuning=='t' or Tuning=='a':
        if Tuning=='a': print '@ @ @ @ @  recreation with tuning shift  @ @ @ @ @'
        Slab, DD, A_STRU, A_SUPERCELL=Slab_Creator(SLAB_DIC[param],REF_VD,REF_A_SUPERCELL.cell,REF_SLAB.cell,Tune_Shift)
        Trans=FINAL_BOX.cell[2]
        Slab.translate(Trans)

        TEMP=FINAL_BOX.copy()
        TEMP=TEMP+Slab
        TEMP.cell[2]=TEMP.cell[2]+Slab.cell[2]
        print 'now the box is with number of atoms = ' ,  TEMP.get_number_of_atoms()
        print 'with this cell equation : \n',  TEMP.cell
        if raw_input('see the joined structure ? (y/n) : ')=='y':
            view(TEMP)
        Tuning=raw_input('Tune the stacking or deny? (t/d) : ')
        if  Tuning == 't':
            Tuning, Tune_Shift=Tuning_shift(TEMP.copy(),Trans,DD,REF_A_STRU)
        else: 
	    FINAL_BOX=TEMP.copy()

    print  '\n', 'ss'*15,'Additional slabs  ','ss'*15
    print '='*72
    MoreSlab=raw_input('more slab to add ? (y/n) : ')
# SF_NAME=raw_input('insert desired SF direction:  ')
# if SF_NAME in SF_DIC:
#     FINAL_BOX=SF_Introducer(FINAL_BOX,REF_A_STRU,SF_NAME)

#################### Cutting the incline of the box to create flat surface ####
'''print  '\n', '<>'*10,'Cutting the model  ','<>'*10
print '='*68
#Want_cut=raw_input('Want to cut the model in y axis ? (y/n) : ')
#while Want_cut== 'y':
YCUT=[]
d_y= 0.8 #interplannar distance in y direction for Gamma phase
#d_y= 1.235 #interplannar distance in y direction for Alpha2 phase
YCUT_min_Sug=float((FINAL_BOX.cell[0][1] + (2* d_y)))
print YCUT_min_Sug
YCUT_max_Sug=float(((FINAL_BOX.cell[1][1]+FINAL_BOX.cell[2][1]) - (2* d_y)))
print YCUT_max_Sug
YCUT.append(float(raw_input('Enter Y_min for cutting from Bottom (For Gamma : 0.8) : '))) #(Should be around'+YCUT_min_Sug+')
YCUT.append(float(raw_input('Enter Y_max for cutting from Top (For Gamma:(b22-|b32|) in box dim. : '))) #(Should be around'+YCUT_max_Sug+') 
FINAL_BOX=SLAB_CUT_Y(FINAL_BOX,YCUT,0)

'''
########################## Show Final shape of model#################
Write_mode=raw_input('Do you want to see the created box? y/n  : ')
if Write_mode =='y':
    view(FINAL_BOX)

#**** Write the final box in a position file for IMD 
print  'ss'*10,' Writing on output ','ss'*10
print '='*72
Write_mode=raw_input('Do you want to write the created box in a IMD format position file? y/n  : ')
if Write_mode =='y':
    OutName=raw_input('Insert output position full file name : \n')
    OutFile=open(OutName,'w')
    print>>OutFile,'#F A 1 1 1 3 0 0'+'\n','#C number type mass x y z'
    print>>OutFile, ' '.join(map(str,['#X ','%.6f' % FINAL_BOX.cell[0][0],'%.6f' % FINAL_BOX.cell[0][1],'%.6f' % FINAL_BOX.cell[0][2]]))
    print>>OutFile, ' '.join(map(str,['#Y ','%.6f' % FINAL_BOX.cell[1][0],'%.6f' % FINAL_BOX.cell[1][1],'%.6f' % FINAL_BOX.cell[1][2]]))
    print>>OutFile, ' '.join(map(str,['#Z ','%.6f' % FINAL_BOX.cell[2][0],'%.6f' % FINAL_BOX.cell[2][1],'%.6f' % FINAL_BOX.cell[2][2]]))
    print>>OutFile,'#E'
    for i, atom in enumerate(FINAL_BOX):
        if atom.symbol == 'Ti': AtomType=0; AtomMass=0.004961
        else: AtomType=1; AtomMass=0.002796
        print>>OutFile, ' '.join(map(str,[atom.index,AtomType,AtomMass,'%.6f' % atom.x,'%.6f' % atom.y,'%.6f' % atom.z]))
    OutFile.close()
    print 'The box configuration file is written into : ' , OutName

########################## Fix the bottom of the model #########################

Write_mode=raw_input('Do you want to Fix the bottom of created box? y/n  : ')
if Write_mode =='y':
	OutName_1=raw_input('Insert output file name : \n')
	OutFile_1=open(OutName_1,'w')
	with open(OutName) as fin, open(OutName_1,'w') as fout:
	 for line in fin:
		if line.startswith("#"):
	   		fout.write(line)
		else :
			fields = line.split()
			#T_Y=float(fields[-2])
			T_Z=float(fields[-1])
					
			#The cutoff distance is 5.678 for TiAl 
			if T_Z < 25 :
				fields[1] = str(int(fields[1])+2)
			fout.write(" ".join(fields) + "\n")
	OutFile_1.close()

# #******************** Printing on screen ********************
print  'ss'*10,' Input information from user','ss'*10
print '='*72 

def Initial_Cell_Reader(ParamName,):
    with open(ParamName) as f_in:
        data = filter(None, (line.rstrip() for line in f_in)) # to store data line by line without empty lines
    StruName=data[0].split()[0]     # structure name
    p1=float(data[1].split()[0])    # lattice constant
    p2=matrix([[float(data[i].split()[0].split(',')[j]) for j in range(0,3)] for i in range(2,5)]) # unit vectors
    p3=int(data[5].split()[0])      # number of given atoms in cell definition
    p4=array([[float(data[i].split()[0].split(',')[j]) for j in range(0,3)] for i in range(6,6+p3)]) # positions of atoms in scaled coordinate
    p5=data[6+p3].split()[0]      # suitable atomic symbol according to ASE standard
    p6=[float(data[6+p3+1].split()[0].split(',')[j]) for j in range(0,3)] # First crystal direction
    p7=[float(data[6+p3+2].split()[0].split(',')[j]) for j in range(0,3)] # Second crystal direction
    p8=[float(data[6+p3+3].split()[0].split(',')[j]) for j in range(0,3)] # Third crystal direction
    n1n2n3=[p6,p7,p8]
    return p1,p2,p3,p4,p5,n1n2n3
#print 'System name: %s  with atomic positions of %d atoms.' % (StruName, Nr0)
#print 'lattice constant %f' % (p1)
'''print "lattice vectors",'\n', p2
print "atom types: " ,p5
    print "Miller indices of the main vectors:",'\n',n1,"---> x",'\n',n2,"---> y",'\n',n3,"---> z"
    print "Cartesian of the main vectors initially:",'\n',matrix(VF0)
    print "Cartesian of the main vectors after orienting:",'\n',matrix(VF),'\n'
    print "Main vector lengths:\n", VFl,'\n'
    print "interplanar distances along the main vectors:",'\n',VD,'\n'
    print "New scaled cell equation: ",'\n',matrix(VF_scaled),'\n'''
