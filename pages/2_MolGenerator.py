import streamlit as st
#from multiapp import MultiApp
from stmol import showmol

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole

from PIL import Image

import py3Dmol

from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np

import leafmap.foliumap as leafmap


st.set_page_config(layout="wide")

# Customize the sidebar

markdown = """
This tool allows you to convert SMILES into 2D, and 3D structures and allows to analyze the basic properties of the molecule. Inspired and code adopted from [Nápoles Duarte's work](https://github.com/napoles-uach)
"""
markdown2 = """
This web app is maintained by [Suneel Kumar BVS](https://www.linkedin.com/in/suneelbommisetty). You can follow me on social media: [GitHub](https://github.com/suneelbvs) | [Twitter](https://twitter.com/suneel_bvs) | [WebSite](sbvs.me).
Streamlit template design inspired from [Qiusheng Wu](https://wetlands.io)
"""

st.sidebar.title("About")
st.sidebar.info(markdown)
st.sidebar.info(markdown2)

def app():
    st.title("MolGenerator")

def image(smi):
    m = Chem.MolFromSmiles(compound_smiles)
    image=Draw.MolToImage(m)
    return image 

def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def render_mol(xyz):
    xyzview = py3Dmol.view()#(width=400,height=400)
    xyzview.addModel(xyz,'mol')
    xyzview.setStyle({'stick':{}})
    xyzview.setBackgroundColor('grey')
    xyzview.zoomTo()
    showmol(xyzview,height=500,width=500)

#st.subheader('Molecular Generator (2D/3D), and its Property Profile')
compound_smiles=st.text_input('Enter SMILES & Press Enter:','OCSCN')
blk=makeblock(compound_smiles)
blk2=image(compound_smiles)

col1, col2 = st.columns((4,4))
with col1:
    col1.subheader("2D Structure of the compound")
    st.image(blk2)
with col2:
    col2.subheader("3D Structure of the compound")
    render_mol(blk)


######################
## Calculate molecular descriptors
def generate(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData= np.arange(1,1)
    i=0
    for mol in moldata:

        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_HBA = Descriptors.NumHAcceptors(mol)
        desc_HBD = Descriptors.NumHDonors(mol)

        row = np.array([desc_MolWt,
        desc_MolLogP,
        desc_HBA,
        desc_HBD,
        ])

        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1

    columnNames=["MolLogP","MolWt","#Acceptors","#Donors"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)

    return descriptors

######################

## Calculate molecular descriptors
st.subheader('Computed Molecular Properties')
X_desc = generate(compound_smiles)
X_desc[1:] # Skips the dummy first item
