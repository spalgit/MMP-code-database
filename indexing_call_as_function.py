#!/usr/bin/env python
# coding: utf-8

# In[2]:


import sys
import re
from rdkit import Chem
# from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from optparse import OptionParser
import pandas as pd


def heavy_atom_count(smi):

    m = Chem.MolFromSmiles(smi)
    return m.GetNumAtoms()

def add_to_index(smi,attachments,cmpd_heavy, use_ratio = None, max_size = 10):

    result = False
    core_size = heavy_atom_count(smi) - attachments

    if(use_ratio):
        core_ratio = float(core_size) / float(cmpd_heavy)
        if(core_ratio <= ratio ):
            result = True
    else:
        if(core_size <= max_size):
            result = True

    return result

def get_symmetry_class(smi):

    symmetry = []

    m = Chem.MolFromSmiles(smi)

    #determine the symmetry class
    #see: http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01894.html
    #A thanks to Greg (and Alan)
    Chem.AssignStereochemistry(m,cleanIt=True,force=True,flagPossibleStereoCenters=True)

    #get the symmetry class of the attachements points
    #Note: 1st star is the zero index,
    #2nd star is first index, etc
    for atom in m.GetAtoms():
        if(atom.GetMass() == 0):
            symmetry.append(atom.GetProp('_CIPRank'))

    return symmetry

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    
    return mol


def get_context_smarts_string(atnum,deg,aromatic,implicitH,bondtype,mapatom):
    
    atomsymb = {}
    atomsymb["6A"] = "C"
    atomsymb["6a"] = "c"
    atomsymb["7a"] = "n"
    atomsymb["7A"] = "N"
    atomsymb["8A"] = "O"
    atomsymb["8a"] = "o"
    atomsymb["16A"] = "S"
    atomsymb["15A"] = "P"
    atomsymb["16a"] = "s"
    atomsymb["17A"] = "Cl"
    atomsymb["9A"] = "F"
    atomsymb["35A"] = "Br"
    atomsymb["53A"] = "I"
    atomsymb["5A"] = "B"
    atomsymb["1A"] = "H"
   
    if (aromatic):
        symb = str(atnum) + "a"
    else:
        symb = str(atnum) + "A"
    
    
    if(bondtype == Chem.BondType.DOUBLE):
        smarts = "=["+atomsymb[symb]+"D"+str(deg)+";"
    elif(bondtype == Chem.BondType.TRIPLE):
        smarts = "#["+atomsymb[symb]+"D"+str(deg)+";"
    else:
        smarts = "["+atomsymb[symb]+"D"+str(deg)+";"

    smarts = smarts + "H" + str(implicitH) + ":" + str(mapatom) + "]"

    return smarts


def add_context_to_smirks(lhs,rhs,context):
    
    
    smirk = "%s>>%s"  % (lhs,rhs) 

    mol = Chem.MolFromSmiles(context)
    mol = mol_with_atom_index(mol)
    
    bonds = []
    
    atomindex = []
    atoms = {}
    neigh = {}
    lv1_keys = []
    
    for atom in mol.GetAtoms():
        atoms[atom.GetIdx()] = atom
        neigh[atom.GetIdx()] = atom.GetNeighbors()

        if(atom.GetMass()==0):
            atomindex.append(atom.GetIdx())
            lv1_keys.append(atom.GetIdx())
            
    lv2_keys = []
    
    
    #Level 1 Smarts
    
    sma = ""
    sma_lv1 = []
    sma_lv1_left = []
    sma_lv1_right = []
    
    mapatom = 1
    for keys in lv1_keys:
        for values in neigh[keys]:
            sma = get_context_smarts_string(values.GetAtomicNum(), values.GetDegree(),
                                            values.GetIsAromatic(), values.GetImplicitValence(),"SINGLE",mapatom)
            sma_lv1.append(sma)
            matchObjlhs = re.search(r'\[\*\:1\]\[H\]', lhs)
            matchObjrhs = re.search(r'\[\*\:1\]\[H\]', rhs)
            if(matchObjlhs):
                smalhs = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree()-1,
                                                  values.GetIsAromatic(),values.GetImplicitValence()+1,"SINGLE",mapatom)
                sma_lv1_left.append(smalhs)
                sma_lv1_right.append(sma)
            elif(matchObjrhs):
                smarhs = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree()-1,
                                                  values.GetIsAromatic(),values.GetImplicitValence()+1,"SINGLE",mapatom)
                sma_lv1_right.append(smarhs)
                sma_lv1_left.append(sma)
            else:
                sma_lv1_right.append(sma)
                sma_lv1_left.append(sma)
            
            lv2_keys.append(values.GetIdx())
            mapatom = mapatom + 1
            
    lv3_keys = []
    sma_lv2 = []
    
    sma_lv3 = []
    smalhs_lv3 = []
    smarhs_lv3 = []
    
    cut = 0

    
    for keys in lv2_keys:
        degree_lv_1 = len(neigh[keys])
        sma = sma_lv1[cut];
        smalhs = sma_lv1_left[cut];
        smarhs = sma_lv1_right[cut];
        
        level2frags = 0
        for values in neigh[keys]:
            if (values.GetAtomicNum() > 0):
                keysn = values.GetIdx()
                degree_lv_2 = len(neigh[keysn])
      
                sma1 = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree(),
                                                values.GetIsAromatic(),values.GetImplicitValence(),
                                                mol.GetBondBetweenAtoms(keys,values.GetIdx()).GetBondType(),mapatom)
                mapatom = mapatom + 1
                sma = sma + "(" + sma1
                smalhs = smalhs + "(" + sma1
                smarhs = smarhs + "(" + sma1
                level3frags = 0
                for valuesn in neigh[keysn]:
                    if(valuesn.GetIdx() not in lv2_keys):
                        sma1 = get_context_smarts_string(valuesn.GetAtomicNum(),valuesn.GetDegree(),
                                                        valuesn.GetIsAromatic(),valuesn.GetImplicitValence(),
                                                        mol.GetBondBetweenAtoms(keysn,valuesn.GetIdx()).GetBondType(),mapatom)
                        mapatom = mapatom + 1
                        if(degree_lv_2 == 2):
                            sma = sma + sma1 +")"
                            smalhs = smalhs + sma1 + ")"
                            smarhs = smarhs + sma1 + ")"
                        elif (degree_lv_2 ==3):
                            if(level3frags == 0):
                                sma = sma + "(" + sma1 +")"
                                smalhs= smalhs + "(" + sma1 +")"
                                smarhs= smarhs + "(" + sma1 + ")"
                            elif(level3frags == 1):
                                sma = sma + sma1 + ")"
                                smalhs = smalhs + sma1 + ")"
                                smarhs = smarhs + sma1 + ")"
                        elif (degree_lv_2 ==4):
                            if(level3frags == 0 or level3frags ==1):
                                sma = sma + "(" + sma1+")"
                                smalhs= smalhs + "(" + sma1 +")"
                                smarhs= smarhs + "(" + sma1 + ")"
                            elif(level3frags == 2):
                                sma = sma + sma1 + ")"
                                smalhs = smalhs + sma1 + ")"
                                smarhs = smarhs + sma1 + ")"
                        level3frags = level3frags + 1
                if(level3frags == 0):
                    sma = sma + ")"
                    smalhs = smalhs + ")"
                    smarhs = smarhs + ")"
                level2frags = level2frags + 1
        cut = cut +1
        sma_lv3.append(sma)
        smalhs_lv3.append(smalhs)
        smarhs_lv3.append(smarhs)
        
        
        matchlhsHyd = re.search(r'\[H\]', lhs)
        matchrhsHyd = re.search(r'\[H\]', rhs)
        
        
        ii = 0
        
        for x in smalhs_lv3:
            matchObjlhs = re.search(r'\[\*\:([123])\]',lhs)
#             matchObjrhs = re.search(r'\[\*\:([123])\]',rhs)
            if (matchlhsHyd==None):
                if matchObjlhs:
                    if(ii==0):
                        lhs = re.sub(r'\[\*:1\]',x,lhs)
                    elif(ii ==1):
                        lhs= re.sub(r'\[\*:2\]',x,lhs)
                    elif(ii ==2):
                        lhs = re.sub(r'\[\*:3]',x,lhs)
                else:
                    lhs = re.sub(r'\*',x,lhs)
            else:
                lhs = x
                
            ii = ii + 1
         
        ii = 0
        
        for x in smarhs_lv3:
            matchObjrhs = re.search(r'\[\*\:([123])\]',rhs)   
            if(matchrhsHyd == None):
                if matchObjrhs:
                    if(ii==0):
                        rhs= re.sub(r'\[\*:1\]',x,rhs)
                    elif(ii ==1):
                        rhs= re.sub(r'\[\*:2\]',x,rhs)
                    elif(ii ==2):
                        rhs = re.sub(r'\[\*:3]',x,rhs)
                else:
                    rhs = re.sub(r'\*',x,rhs)
            else:
                rhs = x          
            ii = ii + 1 
            
        smirk = "%s>>%s" %  (lhs,rhs)
        
#         print("smirk= ", smirk)
        
        
    return smirk
            
            
# def canmirk(lhs,rhs,context):

def cansmirk(lhs,rhs,context):

    #cansmirk algorithm
    #1) cansmi the LHS.
    #2) For the LHS the 1st star will have label 1, 2nd star will have label 2 and so on
    #3) Do a symmetry check of lhs and rhs and use that to decide if the labels on
    #   RHS or/and context need to change.
    #4) For the rhs, if you have a choice (ie. two attachement points are symmetrically
    #   equivalent), always put the label with lower numerical value on the earlier
    #   attachement point on the cansmi-ed smiles

    #print "in: %s,%s" % (lhs,rhs)

    isotope_track={}
    #if the star count of lhs/context/rhs is 1, single cut
    stars = lhs.count("*")

    if(stars > 1):
        #get the symmetry class of stars of lhs and rhs
        lhs_sym = get_symmetry_class(lhs)
        rhs_sym = get_symmetry_class(rhs)

    #deal with double cuts
    if(stars == 2):
        #simple cases
        #unsymmetric lhs and unsymmetric rhs
        if( (lhs_sym[0] != lhs_sym[1]) and (rhs_sym[0] != rhs_sym[1]) ):
            #get 1st and 2nd labels and store the new label for it in isotope_track
            #structure: isotope_track[old_label]=new_label (as strings)
            isotope_track = build_track_dictionary(lhs,stars)

            #switch labels using isotope track
            lhs = switch_labels_on_position(lhs)
            rhs = switch_labels(isotope_track,stars,rhs)
            context = switch_labels(isotope_track,stars,context)
                        
        #symmetric lhs and symmetric rhs
        elif( (lhs_sym[0] == lhs_sym[1]) and (rhs_sym[0] == rhs_sym[1]) ):            
            #the points are all equivalent so change labels on lhs and rhs based on position
            #labels on context don't need to change
            lhs = switch_labels_on_position(lhs)
            rhs = switch_labels_on_position(rhs)
            
        #more difficult cases..
        #symmetric lhs and unsymmetric rhs
        elif( (lhs_sym[0] == lhs_sym[1]) and (rhs_sym[0] != rhs_sym[1]) ):            
            #switch labels lhs based on position
            lhs = switch_labels_on_position(lhs)
            #change labels on rhs based on position but need to record
            #the changes as need to appy them to the context
            isotope_track = build_track_dictionary(rhs,stars)
            rhs = switch_labels_on_position(rhs)
            context = switch_labels(isotope_track,stars,context)
            
        #unsymmetric lhs and symmetric rhs
        elif( (lhs_sym[0] != lhs_sym[1]) and (rhs_sym[0] == rhs_sym[1]) ):            
            #change labels on lhs based on position but need to record
            #the changes as need to appy them to the context
            isotope_track = build_track_dictionary(lhs,stars)
            lhs = switch_labels_on_position(lhs)
            context = switch_labels(isotope_track,stars,context)
            #as rhs is symmetric, positions are equivalent so change labels on position
            rhs = switch_labels_on_position(rhs)
            
    #deal with triple cut
    #unwieldy code but most readable I can make it
    elif(stars == 3):
        #simple cases
        #completely symmetric lhs and completely symmetric rhs
        if( ( (lhs_sym[0] == lhs_sym[1]) and (lhs_sym[1] == lhs_sym[2]) and (lhs_sym[0] == lhs_sym[2]) ) and
        ( (rhs_sym[0] == rhs_sym[1]) and (rhs_sym[1] == rhs_sym[2]) and (rhs_sym[0] == rhs_sym[2]) ) ):
            #the points are all equivalent so change labels on lhs and rhs based on position
            #labels on context don't need to change
            lhs = switch_labels_on_position(lhs)
            rhs = switch_labels_on_position(rhs)            
            
        #completely symmetric lhs and completely unsymmetric rhs
        elif( ( (lhs_sym[0] == lhs_sym[1]) and (lhs_sym[1] == lhs_sym[2]) and (lhs_sym[0] == lhs_sym[2]) ) and
        ( (rhs_sym[0] != rhs_sym[1]) and (rhs_sym[1] != rhs_sym[2]) and (rhs_sym[0] != rhs_sym[2]) ) ):
            #alter lhs in usual way
            lhs = switch_labels_on_position(lhs)
            #change labels on rhs based on position but need to record
            #the changes as need to appy them to the context
            isotope_track = build_track_dictionary(rhs,stars)
            rhs = switch_labels_on_position(rhs)
            context = switch_labels(isotope_track,stars,context)
            
        #completely unsymmetric lhs and completely unsymmetric rhs
        elif( ( (lhs_sym[0] != lhs_sym[1]) and (lhs_sym[1] != lhs_sym[2]) and (lhs_sym[0] != lhs_sym[2]) ) and 
        ( (rhs_sym[0] != rhs_sym[1]) and (rhs_sym[1] != rhs_sym[2]) and (rhs_sym[0] != rhs_sym[2]) ) ):
            #build the isotope track
            isotope_track = build_track_dictionary(lhs,stars)
            #alter lhs in usual way
            lhs = switch_labels_on_position(lhs)
            #change rhs and context based on isotope_track
            rhs = switch_labels(isotope_track,stars,rhs)
            context = switch_labels(isotope_track,stars,context)
                        
        #completely unsymmetric lhs and completely symmetric rhs
        elif( ( (lhs_sym[0] != lhs_sym[1]) and (lhs_sym[1] != lhs_sym[2]) and (lhs_sym[0] != lhs_sym[2]) ) and
        ( (rhs_sym[0] == rhs_sym[1]) and (rhs_sym[1] == rhs_sym[2]) and (rhs_sym[0] == rhs_sym[2]) ) ):
            #build isotope trach on lhs
            isotope_track = build_track_dictionary(lhs,stars)
            #alter lhs in usual way
            lhs = switch_labels_on_position(lhs)
            #change labels on context
            context = switch_labels(isotope_track,stars,context)
            #all positions on rhs equivalent so add labels on position
            rhs = switch_labels_on_position(rhs)            

        #more difficult cases, partial symmetry
        #completely unsymmetric on lhs and partial symmetry on rhs
        elif( (lhs_sym[0] != lhs_sym[1]) and (lhs_sym[1] != lhs_sym[2]) and (lhs_sym[0] != lhs_sym[2]) ):            
            #build the isotope track
            isotope_track = build_track_dictionary(lhs,stars)
            #alter lhs in usual way
            lhs = switch_labels_on_position(lhs)
            #change rhs and context based on isotope_track
            rhs = switch_labels(isotope_track,stars,rhs)
            context = switch_labels(isotope_track,stars,context)

            #tweak positions on rhs based on symmetry
            #rhs 1,2 equivalent
            if(rhs_sym[0] == rhs_sym[1]):
                #tweak rhs position 1 and 2 as they are symmetric
                rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,1,2)

            #rhs 2,3 equivalent
            elif(rhs_sym[1] == rhs_sym[2]):
                #tweak rhs position 1 and 2 as they are symmetric
                rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,2,3)
                
            #rhs 1,3 equivalent - try for larger set in future
            elif(rhs_sym[0] == rhs_sym[2]):
                #tweak rhs position 1 and 2 as they are symmetric
                rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,1,3)
                
        #now we are left with things with partial symmetry on lhs and not completely symmetric or unsymmetric on rhs
        else:
            #lhs 1,2,3 equivalent and any sort of partial symmetry on rhs
            if( (lhs_sym[0] == lhs_sym[1]) and (lhs_sym[1] == lhs_sym[2]) and (lhs_sym[0] == lhs_sym[2]) ):
                
                #alter lhs in usual way
                lhs = switch_labels_on_position(lhs)
                #change labels on rhs based on position but need to record
                #the changes as need to appy them to the context
                isotope_track = build_track_dictionary(rhs,stars)
                rhs = switch_labels_on_position(rhs)
                context = switch_labels(isotope_track,stars,context)
                

            #now deal partial symmetry on lhs or rhs.
            #Cases where:
            #lhs 1,2 equivalent
            #lhs 2,3 equivalent
            #lhs 1,3 equivalent
            else:
                #build isotope track on lhs
                isotope_track = build_track_dictionary(lhs,stars)
                #alter lhs in usual way
                lhs = switch_labels_on_position(lhs)
                #change rhs and context based on isotope_track
                rhs = switch_labels(isotope_track,stars,rhs)
                context = switch_labels(isotope_track,stars,context)

                #tweak positions on rhs based on symmetry

                #lhs 1,2 equivalent
                if(lhs_sym[0] == lhs_sym[1]):
                    #tweak rhs position 1 and 2 as they are symmetric on lhs
                    rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,1,2)                    

                #lhs 2,3 equivalent
                elif(lhs_sym[1] == lhs_sym[2]):
                    #tweak rhs position 1 and 2 as they are symmetric on lhs
                    rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,2,3)
                    
                #lhs 1,3 equivalent - try for larger set in future
                elif(lhs_sym[0] == lhs_sym[2]):
                    #tweak rhs position 1 and 2 as they are symmetric on lhs
                    rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,1,3)

    smirk = "%s>>%s" % (lhs,rhs)

    return smirk,context

def switch_specific_labels_on_symmetry(smi,symmetry_class,a,b):

    #check if a and b positions are symmetrically equivalent
    #if equivalent, swap labels if the lower numerical label is not on the
    #1st symmetrically equivalent attachment points in the smi

    if(symmetry_class[a-1] == symmetry_class[b-1]):
        #what are the labels on a and b

        matchObj = re.search( r'\[\*\:([123])\].*\[\*\:([123])\].*\[\*\:([123])\]', smi )
        if matchObj:
            #if the higher label comes first, fix
            if(int(matchObj.group(a)) > int(matchObj.group(b))):
            #if(int(matchObj.group(1)) > int(matchObj.group(2))):
                smi = re.sub(r'\[\*\:'+matchObj.group(a)+'\]', '[*:XX' + matchObj.group(b) + 'XX]' , smi)
                smi = re.sub(r'\[\*\:'+matchObj.group(b)+'\]', '[*:XX' + matchObj.group(a) + 'XX]' , smi)
                smi = re.sub('XX', '' , smi)

    return smi

def switch_labels_on_position(smi):
    
    #move the labels in order of position
    smi = re.sub(r'\[\*\:[123]\]', '[*:XX1XX]' , smi, 1)
    smi = re.sub(r'\[\*\:[123]\]', '[*:XX2XX]' , smi, 1)
    smi = re.sub(r'\[\*\:[123]\]', '[*:XX3XX]' , smi, 1)
    smi = re.sub('XX', '' , smi)
    
    return smi

def switch_labels(track,stars,smi):

    #switch labels based on the input dictionary track
    if(stars > 1):
        #for k in track:
        #        print "old: %s, new: %s" % (k,track[k])

        if(track['1'] != '1'):
            smi = re.sub(r'\[\*\:1\]', '[*:XX' + track['1'] + 'XX]' , smi)

        if(track['2'] != '2'):
            smi = re.sub(r'\[\*\:2\]', '[*:XX' + track['2'] + 'XX]' , smi)

        if(stars == 3):
            if(track['3'] != '3'):
                smi = re.sub(r'\[\*\:3\]', '[*:XX' + track['3'] + 'XX]' , smi)

        #now remove the XX
        smi = re.sub('XX', '' , smi)

    return smi

def build_track_dictionary(smi,stars):

    isotope_track = {}

    #find 1st label, record it in isotope_track as key, with value being the
    #new label based on its position (1st star is 1, 2nd star 2 etc.)
    if(stars ==2):
        matchObj = re.search( r'\[\*\:([123])\].*\[\*\:([123])\]', smi )
        if matchObj:
            isotope_track[matchObj.group(1)] = '1'
            isotope_track[matchObj.group(2)] = '2'

    elif(stars ==3):
        matchObj = re.search( r'\[\*\:([123])\].*\[\*\:([123])\].*\[\*\:([123])\]', smi )
        if matchObj:
            isotope_track[matchObj.group(1)] = '1'
            isotope_track[matchObj.group(2)] = '2'
            isotope_track[matchObj.group(3)] = '3'

    return isotope_track

def index_hydrogen_change():
    #Algorithm details
    #have an index of common fragment(key) => fragments conected to it (values)
    #Need to add *-H to the values where appropriate - and its
    #appropriate when the key is what you would get if you chopped a H off a cmpd.
    #Therefore simply need to check if key with the * replaced with a H is
    #the same as any full smiles in the set
    #
    #Specific details:
    #1) Loop through keys of index
    #2) If key is the result of a single cut (so contains only 1 *) replace the * with H, and cansmi
    #3) If full smiles matches key in hash above, add *-H to that fragment index.

    for key in index:

        attachments = key.count('*')
        #print attachments

        if(attachments==1):

            smi = key

            #simple method
            smi = re.sub(r'\[\*\:1\]', '[H]' , smi)
            
            #now cansmi it
            temp = Chem.MolFromSmiles(smi)

            if(temp == None):
                sys.stderr.write('Error with key: %s, Added H: %s\n' %(key,smi) )
            else:
                c_smi = Chem.MolToSmiles( temp, isomericSmiles=True )

                if(c_smi in smi_to_id):
                    core = "[*:1][H]"
                    id = smi_to_id[c_smi]

                    value = "%s;t%s" % (id,core)
                    #add to index
                    index[key].append(value)


# In[1]:


def add_context_to_lhs(context):

    context1 = re.sub(r'\[\*\:1\]', '[U]' , context)
    context1 = re.sub(r'\[\*\:2\]', '[V]' , context1)
    context1 = re.sub(r'\[\*\:3\]', '[W]' , context1)

    mol = Chem.MolFromSmiles(context1)
    mol = mol_with_atom_index(mol)
    
    
    bonds = []
    
    atomindex = []
    atoms = {}
    neigh = {}
    lv1_keys = []
    
  # This bit needs a bit re-thinking.  
    
    for atom in mol.GetAtoms():
        atoms[atom.GetIdx()] = atom
        neigh[atom.GetIdx()] = atom.GetNeighbors()
        if(atom.GetSymbol()=="U"):
            atomindex.append(atom.GetIdx())
            lv1_keys.append(atom.GetIdx())
    if(re.search( r'V', context1 )):        
        for atom in mol.GetAtoms():      
            if(atom.GetSymbol()=="V"):
                atomindex.append(atom.GetIdx())
                lv1_keys.append(atom.GetIdx())
    if(re.search( r'W', context1 )):        
        for atom in mol.GetAtoms():
            if(atom.GetSymbol()=="W"):
                atomindex.append(atom.GetIdx())
                lv1_keys.append(atom.GetIdx())
    
    lv2_keys = []
    
    
    ##Level 1 Smarts
    
#     sma = ""
    sma_lv1 = []
    sma_lv1_H = []
    sma_lv1_left = []
    sma_lv1_right = []
    
    mapatom = 1
    for keys in lv1_keys:
        for values in neigh[keys]:
            sma = get_context_smarts_string(values.GetAtomicNum(), values.GetDegree(),
                                            values.GetIsAromatic(), values.GetImplicitValence(),"SINGLE",mapatom)
            sma_lv1.append(sma)
            sma = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree()-1,
                                                  values.GetIsAromatic(),values.GetImplicitValence()+1,"SINGLE",mapatom)
            sma_lv1_H.append(sma)
            
            lv2_keys.append(values.GetIdx())
            mapatom = mapatom + 1
  
    lv3_keys = []
    sma_lv2 = []
    
    sma_lv3 = []
#     sma_lv3 = []
    sma_lv3_H = []
    
    cut = 0

    
    for keys in lv2_keys:
        degree_lv_1 = len(neigh[keys])
        sma = sma_lv1[cut];
        smalhs = sma_lv1[cut];
        smarhs = sma_lv1_H[cut];
        
        level2frags = 0
        
        for values in neigh[keys]:
#             if (values.GetAtomicNum() > 0):
            if (values.GetSymbol() != "U" and values.GetSymbol() != "V" and values.GetSymbol() != "W"):
                keysn = values.GetIdx()
                degree_lv_2 = len(neigh[keysn])
                sma1 = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree(),
                                                values.GetIsAromatic(),values.GetImplicitValence(),
                                                mol.GetBondBetweenAtoms(keys,values.GetIdx()).GetBondType(),mapatom)
#                 print(values.GetSymbol(),values.GetImplicitValence())
                
                mapatom = mapatom + 1
                sma = sma + "(" + sma1
                smalhs = smalhs + "(" + sma1
                smarhs = smarhs + "(" + sma1
                level3frags = 0
                for valuesn in neigh[keysn]:
                    if(valuesn.GetIdx() not in lv2_keys):
                        sma1 = get_context_smarts_string(valuesn.GetAtomicNum(),valuesn.GetDegree(),
                                                        valuesn.GetIsAromatic(),valuesn.GetImplicitValence(),
                                                        mol.GetBondBetweenAtoms(keysn,valuesn.GetIdx()).GetBondType(),mapatom)
                        mapatom = mapatom + 1
                        if(degree_lv_2 == 2):
                            sma = sma + sma1 +")"
                            smalhs = smalhs + sma1 + ")"
                            smarhs = smarhs + sma1 + ")"
                        elif (degree_lv_2 ==3):
                            if(level3frags == 0):
                                sma = sma + "(" + sma1 +")"
                                smalhs= smalhs + "(" + sma1 +")"
                                smarhs= smarhs + "(" + sma1 + ")"
                            elif(level3frags == 1):
                                sma = sma + sma1 + ")"
                                smalhs = smalhs + sma1 + ")"
                                smarhs = smarhs + sma1 + ")"
                        elif (degree_lv_2 ==4):
                            if(level3frags == 0 or level3frags ==1):
                                sma = sma + "(" + sma1+")"
                                smalhs= smalhs + "(" + sma1 +")"
                                smarhs= smarhs + "(" + sma1 + ")"
                            elif(level3frags == 2):
                                sma = sma + sma1 + ")"
                                smalhs = smalhs + sma1 + ")"
                                smarhs = smarhs + sma1 + ")"
                        level3frags = level3frags + 1
                if(level3frags == 0):
                    sma = sma + ")"
                    smalhs = smalhs + ")"
                    smarhs = smarhs + ")"
                level2frags = level2frags + 1
        cut = cut +1
        sma_lv3.append(sma)
        sma_lv3_H.append(smarhs)
        
    
    yield (sma_lv3,sma_lv3_H)
    yield (sma_lv1,sma_lv1_H)
#     print("Level3 = ", sma_lv3,sma_lv3_H)
#     print("Level1 = ", sma_lv1,sma_lv1_H)

#     print("*******")
        
    
    
    


# In[3]:


def index_from_smiles_re(df):
    index={}
    smi_to_id={}
    id_to_smi={}
    for x in df.index:
        id = x
        cmpd_heavy = df.loc[id]['HeavyAtomCount']
        smi = df.loc[id]['smiles']
        smi_to_id[smi]=id
        for y in df.loc[x]['frags']:
            core = y[0]
            context = y[1]
    #         print(core,context)
     #deal with cmpds that have not been fragmented
            if(len(core) == 0) and (len(context) == 0):
                continue

    #deal with single cuts
            if(len(core) == 0):
                side_chains = context.split('.')

                #minus 1 for the attachement pt
                if(add_to_index(side_chains[1],1,cmpd_heavy)==True ):
                    context = side_chains[0]
                    core = side_chains[1]

                    value = "%s;t%s" % (id,core)

                    #add the array if no key exists
                    #add the context with id to index
                    index.setdefault(context, []).append(value)

                #minus 1 for the attachement pt
                if( add_to_index(side_chains[0],1,cmpd_heavy)==True ):
                    context = side_chains[1]
                    core = side_chains[0]

                    value = "%s;t%s" % (id,core)

                    #add the array if no key exists
                    #add the context with id to index
                    index.setdefault(context, []).append(value)

    #double or triple cut
            else:
                attachments = core.count('*')

                if( add_to_index(core,attachments,cmpd_heavy)==True ):
                    value = "%s;t%s" % (id,core)

                    #add the array if no key exists
                    #add the context with id to index
                    index.setdefault(context, []).append(value) 
                    
    for key in index:
        attachments = key.count('*')
        if(attachments==1):

            smi = key

            #simple method
            smi = re.sub(r'\[\*\:1\]', '[H]' , smi)

            #now cansmi it
            temp = Chem.MolFromSmiles(smi)

            if(temp == None):
                sys.stderr.write('Error with key: %s, Added H: %s\n' %(key,smi) )
            else:
                c_smi = Chem.MolToSmiles( temp, isomericSmiles=True )
                if(c_smi in smi_to_id):
                    core = "[*:1][H]"
                    id = smi_to_id[c_smi]
                    value = "%s;t%s" % (id,core)
                    #add to index
                    index[key].append(value)
    return index   


# In[4]:


def index_from_smiles(df):
    index={}
    smi_to_id={}
    id_to_smi={}
    for i in range (0,len(df)):
        id = df.loc[i]['Name']
        cmpd_heavy = df.loc[i]['HeavyAtomCount']
        smi = df.loc[i]['smiles']
        smi_to_id[smi]=id
        id_to_smi[id]=smi
        for y in df.loc[i]['frags']:
            core = y[0]
            context = y[1]
    #         print(core,context)
     #deal with cmpds that have not been fragmented
            if(len(core) == 0) and (len(context) == 0):
                continue

    #deal with single cuts
            if(len(core) == 0):
                side_chains = context.split('.')

                #minus 1 for the attachement pt
                if(add_to_index(side_chains[1],1,cmpd_heavy)==True ):
                    context = side_chains[0]
                    core = side_chains[1]

                    value = "%s;t%s" % (id,core)

                    #add the array if no key exists
                    #add the context with id to index
                    index.setdefault(context, []).append(value)

                #minus 1 for the attachement pt
                if( add_to_index(side_chains[0],1,cmpd_heavy)==True ):
                    context = side_chains[1]
                    core = side_chains[0]

                    value = "%s;t%s" % (id,core)

                    #add the array if no key exists
                    #add the context with id to index
                    index.setdefault(context, []).append(value)

    #double or triple cut
            else:
                attachments = core.count('*')

                if( add_to_index(core,attachments,cmpd_heavy)==True ):
                    value = "%s;t%s" % (id,core)

                    #add the array if no key exists
                    #add the context with id to index
                    index.setdefault(context, []).append(value) 
                    
    for key in index:
        attachments = key.count('*')
        if(attachments==1):

            smi = key

            #simple method
            smi = re.sub(r'\[\*\:1\]', '[H]' , smi)

            #now cansmi it
            temp = Chem.MolFromSmiles(smi)

            if(temp == None):
                sys.stderr.write('Error with key: %s, Added H: %s\n' %(key,smi) )
            else:
                c_smi = Chem.MolToSmiles( temp, isomericSmiles=True )
                if(c_smi in smi_to_id):
                    core = "[*:1][H]"
                    id = smi_to_id[c_smi]
                    value = "%s;t%s" % (id,core)
                    #add to index
                    index[key].append(value)
    return index


# In[5]:


def gen_queries(lhs,context):
#     smirk = "%s>>%s"  % (lhs,rhs) 

    mol = Chem.MolFromSmiles(context)
    mol = mol_with_atom_index(mol)
    
    bonds = []
    
    atomindex = []
    atoms = {}
    neigh = {}
    lv1_keys = []
    
    for atom in mol.GetAtoms():
        atoms[atom.GetIdx()] = atom
        neigh[atom.GetIdx()] = atom.GetNeighbors()

        if(atom.GetMass()==0):
            atomindex.append(atom.GetIdx())
            lv1_keys.append(atom.GetIdx())
            
    lv2_keys = []
    
    
    #Level 1 Smarts
    
    sma = ""
    sma_lv1 = []
    sma_lv1_left = []
    sma_lv1_right = []
    
    mapatom = 1
    for keys in lv1_keys:
        for values in neigh[keys]:
            sma = get_context_smarts_string(values.GetAtomicNum(), values.GetDegree(),
                                            values.GetIsAromatic(), values.GetImplicitValence(),"SINGLE",mapatom)
            sma_lv1.append(sma)
            matchObjlhs = re.search(r'\[\*\:1\]\[H\]', lhs)
#             matchObjrhs = re.search(r'\[\*\:1\]\[H\]', rhs)
            if(matchObjlhs):
                smalhs = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree()-1,
                                                  values.GetIsAromatic(),values.GetImplicitValence()+1,"SINGLE",mapatom)
                sma_lv1_left.append(smalhs)
                sma_lv1_right.append(sma)
#             elif(matchObjrhs):
#                 smarhs = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree()-1,
#                                                   values.GetIsAromatic(),values.GetImplicitValence()+1,"SINGLE",mapatom)
#                 sma_lv1_right.append(smarhs)
#                 sma_lv1_left.append(sma)
            else:
#                 sma_lv1_right.append(sma)
                sma_lv1_left.append(sma)
            
            lv2_keys.append(values.GetIdx())
            mapatom = mapatom + 1
            
    lv3_keys = []
    sma_lv2 = []
    
    sma_lv3 = []
    smalhs_lv3 = []
    smarhs_lv3 = []
    
    cut = 0

    
    for keys in lv2_keys:
        degree_lv_1 = len(neigh[keys])
        sma = sma_lv1[cut];
        smalhs = sma_lv1_left[cut];
#         smarhs = sma_lv1_right[cut];
        
        level2frags = 0
        for values in neigh[keys]:
            if (values.GetAtomicNum() > 0):
                keysn = values.GetIdx()
                degree_lv_2 = len(neigh[keysn])
                sma1 = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree(),
                                                values.GetIsAromatic(),values.GetImplicitValence(),
                                                mol.GetBondBetweenAtoms(keys,values.GetIdx()).GetBondType(),mapatom)
                mapatom = mapatom + 1
                sma = sma + "(" + sma1
                smalhs = smalhs + "(" + sma1
#                 smarhs = smarhs + "(" + sma1
                level3frags = 0
                for valuesn in neigh[keysn]:
                    if(valuesn.GetIdx() not in lv2_keys):
                        sma1 = get_context_smarts_string(valuesn.GetAtomicNum(),valuesn.GetDegree(),
                                                        valuesn.GetIsAromatic(),valuesn.GetImplicitValence(),
                                                        mol.GetBondBetweenAtoms(keysn,valuesn.GetIdx()).GetBondType(),mapatom)
                        mapatom = mapatom + 1
                        if(degree_lv_2 == 2):
                            sma = sma + sma1 +")"
                            smalhs = smalhs + sma1 + ")"
#                             smarhs = smarhs + sma1 + ")"
                        elif (degree_lv_2 ==3):
                            if(level3frags == 0):
                                sma = sma + "(" + sma1 +")"
                                smalhs= smalhs + "(" + sma1 +")"
#                                 smarhs= smarhs + "(" + sma1 + ")"
                            elif(level3frags == 1):
                                sma = sma + sma1 + ")"
                                smalhs = smalhs + sma1 + ")"
#                                 smarhs = smarhs + sma1 + ")"
                        elif (degree_lv_2 ==4):
                            if(level3frags == 0 or level3frags ==1):
                                sma = sma + "(" + sma1+")"
                                smalhs= smalhs + "(" + sma1 +")"
#                                 smarhs= smarhs + "(" + sma1 + ")"
                            elif(level3frags == 2):
                                sma = sma + sma1 + ")"
                                smalhs = smalhs + sma1 + ")"
#                                 smarhs = smarhs + sma1 + ")"
                        level3frags = level3frags + 1
                if(level3frags == 0):
                    sma = sma + ")"
                    smalhs = smalhs + ")"
#                     smarhs = smarhs + ")"
                level2frags = level2frags + 1
        cut = cut +1
        sma_lv3.append(sma)
        smalhs_lv3.append(smalhs)
#         smarhs_lv3.append(smarhs)
        
        
        matchlhsHyd = re.search(r'\[H\]', lhs)
#         matchrhsHyd = re.search(r'\[H\]', rhs)
        
        
        ii = 0
        
        for x in smalhs_lv3:
#             print(x)
            matchObjlhs = re.search(r'\[\*\:([123])\]',lhs)
#             matchObjrhs = re.search(r'\[\*\:([123])\]',rhs)
            if (matchlhsHyd==None):
                if matchObjlhs:
                    if(ii==0):
                        lhs = re.sub(r'\[\*:1\]',x,lhs)
                    elif(ii ==1):
                        lhs= re.sub(r'\[\*:2\]',x,lhs)
                    elif(ii ==2):
                        lhs = re.sub(r'\[\*:3]',x,lhs)
                else:
                    lhs = re.sub(r'\*',x,lhs)
            else:
                lhs = x
                
            ii = ii + 1
         
        ii = 0
        
#         for x in smarhs_lv3:
#             matchObjrhs = re.search(r'\[\*\:([123])\]',rhs)   
#             if(matchrhsHyd == None):
#                 if matchObjrhs:
#                     if(ii==0):
#                         rhs= re.sub(r'\[\*:1\]',x,rhs)
#                     elif(ii ==1):
#                         rhs= re.sub(r'\[\*:2\]',x,rhs)
#                     elif(ii ==2):
#                         rhs = re.sub(r'\[\*:3]',x,rhs)
#                 else:
#                     rhs = re.sub(r'\*',x,rhs)
#             else:
#                 rhs = x          
#             ii = ii + 1 
            
#         smirk = "%s>>%s" %  (lhs,rhs)
        
#         print("smirk= ", smirk)
        
        
    return lhs
    
    


# In[1]:


def gen_queries_for_H_transform(lhs,context):

    mol = Chem.MolFromSmiles(context)
    mol = mol_with_atom_index(mol)

    bonds = []
    
    atomindex = []
    atoms = []
    neigh = {}
    lv1_keys = []
    
    for atom in mol.GetAtoms():
        neigh[atom.GetIdx()] = atom.GetNeighbors()
        if(atom.GetNumImplicitHs()>=1):
            atomindex.append(atom.GetIdx())
            lv1_keys.append(atom.GetIdx())
            atoms.append(atom)
  
    lv2_keys = []

    #Level 1 Smarts
    
    sma = ""
    sma_lv1 = []
    sma_lv1_left = []
    sma_lv1_right = []
    
    mapatom = 1
    
    for values in atoms:
        sma = get_context_smarts_string(values.GetAtomicNum(), values.GetDegree(),
                                            values.GetIsAromatic(), values.GetImplicitValence(),"SINGLE",mapatom)
        sma_lv1.append(sma)
        lv2_keys.append(values.GetIdx())
    
#     print("lv2_keys =", lv2_keys)
#     print(sma_lv1)
    lv3_keys = []
    sma_lv2 = []
    
    sma_lv3 = []
    
    cut = 0
    
#     print(sma_lv1)
    
#     for keys in lv2_keys:
#         degree_lv_1 = len(neigh[keys])
#         sma = sma_lv1[cut];
#         level2frags = 0
#         mapatom = 2
#         print("0 = ", keys, len(neigh[keys]), sma)
#         for values in neigh[keys]:
#             keysn = values.GetIdx()
#             degree_lv_2 = len(neigh[keysn])
#             sma1 = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree(),
#                                             values.GetIsAromatic(),values.GetImplicitValence(),
#                                             mol.GetBondBetweenAtoms(keys,values.GetIdx()).GetBondType(),mapatom)
#             mapatom = mapatom + 1
#             sma = sma + "(" + sma1
#             level3frags = 0
#             print("1 = ", sma1)
#             for valuesn in neigh[keysn]:
#                 print("valuesn.GetIdx() = ", valuesn.GetIdx())
# #                 if(valuesn.GetIdx() not in lv2_keys):
#                 sma1 = get_context_smarts_string(valuesn.GetAtomicNum(),valuesn.GetDegree(),
#                                                 valuesn.GetIsAromatic(),valuesn.GetImplicitValence(),
#                                                 mol.GetBondBetweenAtoms(keysn,valuesn.GetIdx()).GetBondType(),mapatom)
#                 print("2 = ", sma1)
            
    
    
    for keys in lv2_keys:
        degree_lv_1 = len(neigh[keys])
        sma = sma_lv1[cut];
        level2frags = 0
        mapatom = 2
        for values in neigh[keys]:
            keysn = values.GetIdx()
            degree_lv_2 = len(neigh[keysn])
            sma1 = get_context_smarts_string(values.GetAtomicNum(),values.GetDegree(),
                                            values.GetIsAromatic(),values.GetImplicitValence(),
                                            mol.GetBondBetweenAtoms(keys,values.GetIdx()).GetBondType(),mapatom)
            mapatom = mapatom + 1
            sma = sma + "(" + sma1
            level3frags = 0
            for valuesn in neigh[keysn]:
                if(valuesn.GetIdx() != keys):
                    sma1 = get_context_smarts_string(valuesn.GetAtomicNum(),valuesn.GetDegree(),
                                                    valuesn.GetIsAromatic(),valuesn.GetImplicitValence(),
                                                    mol.GetBondBetweenAtoms(keysn,valuesn.GetIdx()).GetBondType(),mapatom)
                    mapatom = mapatom + 1
                    if(degree_lv_2 == 2):
                        sma = sma + sma1 +")"
                    elif (degree_lv_2 ==3):
                        if(level3frags == 0):
                            sma = sma + "(" + sma1 +")"
                        elif(level3frags == 1):
                            sma = sma + sma1 + ")"
                    elif (degree_lv_2 ==4):
                        if(level3frags == 0 or level3frags ==1):
                            sma = sma + "(" + sma1+")"
                        elif(level3frags == 2):
                            sma = sma + sma1 + ")"
                    level3frags = level3frags + 1
            if(level3frags == 0):
                sma = sma + ")"
            level2frags = level2frags + 1
        cut = cut +1
        sma_lv3.append(sma)
    
    return(sma_lv3,sma_lv1)
#     yield sma_lv3
#     yield sma_lv1
#     return sma_lv3


# In[ ]:


def atom_conn_to_H(comp):
    mol = Chem.MolFromSmiles(context)
    mol = mol_with_atom_index(mol)

    for atom in mol.GetAtoms():
        if(atom.GetNumImplicitHs()>=1):
            atomindex.append(atom.GetIdx())
    
    
    
    


# In[7]:


if __name__=='__main__':

    #note max heavy atom count does not
    #include the attachement points (*)
    max_size = 10
    ratio = 0.3
    use_ratio = False

    index={}
    smi_to_id={}
    id_to_smi={}

    id_to_heavy={}

    #set up the command line options
    #parser = OptionParser()
#     parser = OptionParser(description="Program to generate MMPs")
#     parser.add_option('-s', '--symmetric', default=False, action='store_true', dest='sym',
#                       help='Output symmetrically equivalent MMPs, i.e output both cmpd1,cmpd2, SMIRKS:A>>B and cmpd2,cmpd1, SMIRKS:B>>A')
#     parser.add_option('-m','--maxsize',action='store', dest='maxsize', type='int',
#                       help='Maximum size of change (in heavy atoms) allowed in matched molecular pairs identified. DEFAULT=10. \
#                       Note: This option overrides the ratio option if both are specified.')
#     parser.add_option('-r','--ratio',action='store', dest='ratio', type='float',
#                       help='Maximum ratio of change allowed in matched molecular pairs identified. The ratio is: size of change / \
#                       size of cmpd (in terms of heavy atoms). DEFAULT=0.3. Note: If this option is used with the maxsize option, the maxsize option will be used.')

#     #parse the command line options
#     (options, args) = parser.parse_args()

#     #print options
#     if(options.maxsize != None):
#         max_size = options.maxsize
#     elif(options.ratio != None):
#         ratio = options.ratio
#         if(ratio >= 1):
#             print "Ratio specified: %s. Ratio needs to be less than 1."
#             sys.exit(1)
#         use_ratio = True

    #read the STDIN
    
#     datafile = "/home/spal/MMP_source_code/Mol_Fragments.txt"
    datafile = "C:/Users/sande/DataFiles/Approved_Drugs_Fragments.txt.txt"

#     args = sys.argv
#     datafile = open(args[1], "r")
    
    for line in open(datafile):

        line = line.rstrip()
        smi,id,core,context = line.split(',')

        #fill in dictionaries
        smi_to_id[smi]=id
        id_to_smi[id]=smi
        
#         print("smi = ", smi)

        #if using the ratio option, check if heavy atom
        #of mol already calculated. If not, calculate and store
        cmpd_heavy = None
        if(use_ratio):
            if( (id in id_to_heavy) == False):
                id_to_heavy[id] = heavy_atom_count(smi)

            cmpd_heavy = id_to_heavy[id]

        #deal with cmpds that have not been fragmented
        if(len(core) == 0) and (len(context) == 0):
            continue

        #deal with single cuts
        if(len(core) == 0):
            side_chains = context.split('.')

            #minus 1 for the attachement pt
            if( add_to_index(side_chains[1],1,cmpd_heavy)==True ):
                context = side_chains[0]
                core = side_chains[1]

                value = "%s;t%s" % (id,core)

                #add the array if no key exists
                #add the context with id to index
                index.setdefault(context, []).append(value)

            #minus 1 for the attachement pt
            if( add_to_index(side_chains[0],1,cmpd_heavy)==True ):
                context = side_chains[1]
                core = side_chains[0]

                value = "%s;t%s" % (id,core)

                #add the array if no key exists
                #add the context with id to index
                index.setdefault(context, []).append(value)

        #double or triple cut
        else:

            attachments = core.count('*')

            if( add_to_index(core,attachments,cmpd_heavy)==True ):
                value = "%s;t%s" % (id,core)

                #add the array if no key exists
                #add the context with id to index
                index.setdefault(context, []).append(value)

    #index the H change
    index_hydrogen_change()

    #Now index is ready

    #loop through the index
    for key in index:

        total = len(index[key])

        #check if have more than one value
        if(total == 1):
            continue

        for xa in range(total):

            for xb in range(xa, total):

                if(xa != xb):
                    #now generate the pairs

                    id_a,core_a = index[key][xa].split(";t")
                    id_b,core_b = index[key][xb].split(";t")

                    #make sure pairs are not same molecule
                    if(id_a != id_b):

                        #make sure LHS and RHS of SMIRKS are not the same
                        if(core_a != core_b):

                            smirks,context = cansmirk(core_a,core_b,key)
                            
                            lhs,rhs = smirks.split(">>")
                            
                            
#                             print("lhs =", lhs,"rhs = ", rhs,"context = ",context)
                            
                            newsmirks = add_context_to_smirks(lhs,rhs,context)
                            
                            
#                             print(newsmirks, smirks,id_to_smi[id_a],id_to_smi[id_b],id_a,id_b)
#                             mol1 = Chem.MolFromSmarts(newsmirks)
#                             molSizes=(450,450)
#                             mc = Chem.Mol(mol1.ToBinary())
#                             drawer = rdMolDraw2D.MolDraw2DVG()
    
                    
                            print ("%s,%s,%s,%s,%s,%s,%s,%s" % (lhs,id_to_smi[id_a], id_to_smi[id_b], id_a, id_b,newsmirks,smirks,context))

                            newlhs,newrhs = newsmirks.split(">>")
                
                            smirks = "%s>>%s" %(rhs,lhs)
                            newsmirks = "%s>>%s" %(newrhs,newlhs)
                            
                            print ("%s,%s,%s,%s,%s,%s,%s,%s" % (rhs,id_to_smi[id_b], id_to_smi[id_a], id_b, id_a,newsmirks,smirks,context))
                            

