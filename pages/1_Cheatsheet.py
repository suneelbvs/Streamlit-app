import streamlit as st
from pathlib import Path
#import base64

import leafmap.foliumap as leafmap

st.set_page_config(layout="wide")

# Customize the sidebar

markdown = """
We have cheatsheets for Python, ML, Tensorflow, Deep learning. Why can't we have cheatsheet for rdkit. Here is one - This is the first version of rdkit cheatsheet. I have curated basic functions of rdkit (such as., converters, descriptors, alignment, file handling, colab installations). Im working on the improvement of this cheatsheet and suggestions are welcome.
Pls DM your suggestions to [My Twitter handle](https://twitter.com/suneel_bvs).
"""
markdown2 = """
This web app is maintained by [Suneel Kumar BVS](https://www.linkedin.com/in/suneelbommisetty). You can follow me on social media: [GitHub](https://github.com/suneelbvs) | [Twitter](https://twitter.com/suneel_bvs) | [WebSite](sbvs.me).
Streamlit template design inspired from [Qiusheng Wu](https://wetlands.io)
"""

st.sidebar.title("About")
st.sidebar.info(markdown)
st.sidebar.info(markdown2)

def app():
    st.title='rdKit cheat sheet',
#    layout="wide"
#    initial_sidebar_state="expanded"

# Thanks to streamlitopedia for the following code snippet

def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded

##########################
# Main body of cheat sheet
##########################

#def cs_body():
    # Magic commands

st.title='rdKit cheat sheet'

col1, col2 = st.columns(2)

col1.subheader('Installation')
col1.code('''# Magic commands
>> conda install -c rdkit rdkit

# Check rdKit version
>> rdkit.__version__

#Import packages
>> import rdkit
>> from rdkit import Che
''')

# Display text
col1.subheader('Google Colab')
col1.code('''
#rdKit installation in Colab:
>> !wget -c https://repo.continuum.io/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh
>> !chmod +x Miniconda3-py37_4.8.3-Linux-x86_64.sh
>> !time bash ./Miniconda3-py37_4.8.3-Linux-x86_64.sh -b -f -p /usr/local
>> !time conda install -q -y -c conda-forge rdkit
import sys
sys.path.append('/usr/local/lib/python3.7/site-packages/')

# try to install condacolab first, to avoid any errors
>> !pip install -q condacolab
>> import condacolab
>> condacolab.install()

''')


# Display data

col1.subheader('Working with Formats')
col1.code('''
# SMILES to Mol
>> mol1 = Chem.MolFromSmiles('COc1ccc2c(c1)[nH]c(n2)[S@@](=O)Cc1ncc(c(c1C)OC)C')
>> mol1

# Converting Mol to Inchi
>> Chem.MolToInchiKey(mol)

# MolBlock (coordinative representation)
>> m2 = Chem.MolFromSmiles('COc1ccc2c(c1)[nH]c(n2)[S@@](=O)Cc1ncc(c(c1C)OC)C')
>> print(Chem.MolToMolBlock(m2))

# Writing SD
>> w1 = Chem.SDWriter('./data/data_set.sdf')
>> w1.write(m2)
>> w1.close()

# Writing SMILES
>> w2 = Chem.SmilesWriter('./data/data_set.smi')
>> w2.write(m2)
>> w2.close()

''')

# Display charts
col1.subheader('Fingerprint Calculation')
col1.code('''

#Calculate FingerPrints for two Mols:
Load Libraries:
>> from rdkit import Chem
>> from rdkit import DataStructs
>> from rdkit.Chem.Draw import IPythonConsole
>> from rdkit.Chem import AllChem
>> from rdkit.six import StringIO

#  define mols
>> mol1 = Chem.MolFromSmiles("Cc1ccccc1")
>> mol2 = Chem.MolFromSmiles("Clc1ccccc1")

Calculate Fingerprints:
>> fp1 = AllChem.GetMorganFingerprint(mol1, 2)
>> fp2 = AllChem.GetMorganFingerprint(mol2, 2)
>> DataStructs.TanimotoSimilarity(fp1, fp2)

# Calculate FingerPrints for datasets
>> esol_data = pd.read_csv(Path for dataset)
>> PandasTools.AddMoleculeColumnToFrame(esol_data, smilesCol='smiles')
>> esol_data.head(1)
# It generates pandas table with ROMol.
# Now, lets calculate Morgan Fingerprint (ECFPx):
# radius: no default value, usually set 2 for similarity search and 3 for machine learning.
# nBits: number of bits, default is 2048. 1024 is also widely used.

>> radius=3
>> nBits=1024

>> ECFP6 = [AllChem.GetMorganFingerprintAsBitVect(x,radius=radius, nBits=nBits) for x in esol_data['ROMol']]
>> ECFP6[0]
>> len(ECFP6[0])

Save as a .csv file for futher use (e.g., machine learning). I usually save (1) SMILES as index and (2) each bit as a column to the csv file.

>> ecfp6_name = [f'Bit_{i}' for i in range(nBits)]
>> ecfp6_bits = [list(l) for l in ECFP6]
>. df_morgan = pd.DataFrame(ecfp6_bits, index = esol_data.smiles, columns=ecfp6_name)
>. df_morgan.head(1)
## Source: Xinhao Li
''')

# Display media

col1.subheader('Compute Descriptors')
col1.code('''
#Import Libraries
>> from rdkit import Chem
>> from rdkit.Chem.Draw import IPythonConsole
>> from rdkit.Chem import Draw
# Load Mols
>> mol1 = Chem.MolFromSmiles('COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1')
>> print ('Heavy atoms:', Descriptors.HeavyAtomCount(mol2))
#output : Heavy atoms: 29

>> print ('H-bond acceptors:', Descriptors.NumHAcceptors(mol2))
# H-bond donors: 1

## Compute selected Descriptors
# Calculate some selected descriptors only
>> des_list = ['MolWt', 'NumHAcceptors', 'NumHDonors', 'MolLogP', 'NumRotatableBonds']
>> calculator = MoleculeDescriptors.MolecularDescriptorCalculator(des_list)
>> calculator.CalcDescriptors(mol)
>> des_list = [x[0] for x in Descriptors._descList]
>> calculator = MoleculeDescriptors.MolecularDescriptorCalculator(des_list)
>> d = calculator.CalcDescriptors(mol)
>> d #to print the calculated descriptors
## Source: John Mommers

''')

# Display interactive widgets
col1.subheader('Play with Substructure')
col1.code('''
>> from rdkit import Chem
>> from rdkit.Chem.Draw import IPythonConsole

m = Chem.MolFromSmiles('c1cc(C(=O)O)c(OC(=O)C)cc1')

substructure = Chem.MolFromSmarts('C(=O)O')
print(m.GetSubstructMatches(substructure))
m
#Source: Greg Landrum

#Generic (“Markush”) queries in substructure matching
>>> q = Chem.MolFromSmarts('OC* |$;;ARY$|')
>>> Chem.SetGenericQueriesFromProperties(q)
>>> Chem.MolFromSmiles('C1CCCCC1CO').HasSubstructMatch(q)
## True

''')

col2.subheader('Align molecules')
col2.code('''
>>> with Chem.SDMolSupplier('data/cdk2.sdf') as suppl:
...   ms = [x for x in suppl if x is not None]
>>> for m in ms: tmp=AllChem.Compute2DCoords(m)
>>> from rdkit.Chem import Draw
>>> Draw.MolToFile(ms[0],'images/cdk2_mol1.o.png')    
>>> Draw.MolToFile(ms[1],'images/cdk2_mol2.o.png')  

>>> p = Chem.MolFromSmiles('[nH]1cnc2cncnc21')
>>> subms = [x for x in ms if x.HasSubstructMatch(p)]
>>> len(subms)
14
>>> AllChem.Compute2DCoords(p)
0
>>> for m in subms:
...   _ = AllChem.GenerateDepictionMatching2DStructure(m,p)
>>> img=Draw.MolsToGridImage(subms,molsPerRow=4,subImgSize=(200,200),legends=[x.GetProp("_Name") for x in subms])    
>>> img.save('images/cdk2_molgrid.aligned.o.png')  

#Source: rdkit.org
''')

# Control flow
col2.subheader('Play with Pandas')
col2.code('''
# PandasTools enables using RDKit molecules as columns of a Pandas Dataframe
>> import pandas as pd
>> from rdkit.Chem import PandasTools

# reading CSV file (with SMILES) 
>> df = pd.read_csv('./data/EGFR_compounds.csv', index_col=0)
>> print(df.shape)
>> df.head()

# Convert SMILES into ROMol
> PandasTools.AddMoleculeColumnToFrame(df, smilesCol='SMILES')

# Display the ROMol using GridImage
>> PandasTools.FrameToGridImage(df.head(8), legendsCol='name', molsPerRow=4)

# Write SDF using Pandas 
>> PandasTools.WriteSDF(df, '../data/file.sdf', idName='ID', properties=df.columns)

# Delete empty rows --> rule of five
>> filtered_df = df[df['rule_of_five_conform']=='yes']

''')

# Lay out your app

col2.subheader('Handling Larger files')
col2.code('''
# Reading SMI file without Pandas Support
>> file_name = 'somedata.smi'
>> with open(file_name, "r") as ins:
>>    smiles = []
>>    for line in ins:
>>        smiles.append(line.split('\n')[0])
>> print('# of SMILES:', len(smiles))

# SD file reading
>>> suppl = Chem.SDMolSupplier('data/5ht3ligs.sdf')
>>> for mol in suppl:
...   print(mol.GetNumAtoms())

# Reading Compressed files:
>>> import gzip
>>> inf = gzip.open('data/actives_5ht3.sdf.gz')
>>> with Chem.ForwardSDMolSupplier(inf) as gzsuppl:
>>>    ms = [x for x in gzsuppl if x is not None]
>>> len(ms)
''')

col2.subheader(' Substructure Searching')
col2.code('''
#SMARTS based Substructure matching 
>>> m = Chem.MolFromSmiles('c1ccccc1O')
>>> patt = Chem.MolFromSmarts('ccO')
>>> m.HasSubstructMatch(patt)
True
>>> m.GetSubstructMatch(patt)
(0, 5, 6)

#'alkyl’ query
>>> p.GetAtomWithIdx(4).SetProp("queryType", "alkyl")
>>> checker = SidechainChecker(p)
>>> params.setExtraFinalCheck(checker)
>>> m.GetSubstructMatches(p,params)
((5, 17, 11, 6, 7), (6, 5, 17, 11, 12), (6, 11, 17, 5, 4))
#Source:rdkit.org

''')

# Display code

#col2.subheader('Display code')
#col2.code('''
#st.echo()
#>>> with st.echo():
#>>>     st.write('Code will be executed and printed')
#''')

# Display progress and status

col2.subheader('StereoAnnotations')
col2.code('''
>>> from rdkit import Chem
>>> from rdkit.Chem import Draw
>>> from rdkit.Chem.Draw import IPythonConsole
>>> IPythonConsole.drawOptions.addAtomIndices = False
>>> IPythonConsole.drawOptions.addStereoAnnotation = True

# Default Representation uses legacy FindMolChiralCenters() code
>>> m1 = Chem.MolFromSmiles('C1CC1[C@H](F)C1CCC1')
>>> m2 = Chem.MolFromSmiles('F[C@H]1CC[C@H](O)CC1')
>>> Draw.MolsToGridImage((m1,m2), subImgSize=(250,250))

# new stereochemistry code with more accurate CIP labels, 2020.09 release
>>> from rdkit.Chem import rdCIPLabeler
>>> rdCIPLabeler.AssignCIPLabels(m1)
>>> rdCIPLabeler.AssignCIPLabels(m2)
>>> Draw.MolsToGridImage((m1,m2), subImgSize=(250,250))

Source: rdkitbook
''')

# rdKit, IPython, and functions

col2.subheader('IPython Console Commands')
col2.code('''
# adjust molsize 
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = 250,250
#m = Chem.MolFromSmiles('c1ncncc1C(=O)[O-]')
#AllChem.ComputeGasteigerCharges(m)
#m

from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.drawOptions.addAtomIndices = False
IPythonConsole.drawOptions.addStereoAnnotation = True
#m2 = Chem.MolFromSmiles('F[C@H]1CC[C@H](O)CC1')
#Draw.MolsToGridImage((m1,m2), subImgSize=(250,250))

# Find chiral centers and double bond stereochemistry.

from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.drawOptions.addStereoAnnotation = False
IPythonConsole.molSize = 200,200

m = Chem.MolFromSmiles("C[C@H]1CCC[C@@H](C)[C@@H]1Cl")
m
''')