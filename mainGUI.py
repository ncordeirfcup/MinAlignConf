import tkinter as tk
from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk
import os
from tkinter.filedialog import askopenfilename
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SaltRemover
from align_ligands import align_set_of_ligands
from conformer import process
import sys

form = tk.Tk()
form.title("ChEMBL_autocurator")
form.geometry("650x380")

tab_parent = ttk.Notebook(form)

tab1 = tk.Frame(tab_parent, background="#00ffff")
tab_parent.add(tab1, text="Data preparation")

initialdir=os.getcwd()


def datatr():
    global filename1
    filename1 = askopenfilename(initialdir=initialdir,title = "Select downloaded ChEMBL csv file")
    firstEntryTabOne.delete(0, END)
    firstEntryTabOne.insert(0, filename1)
    global c_
    c_,d_=os.path.splitext(filename1)
    global file1
    dm=Criterion0.get()
    print(dm)
    if dm=="com":
       file1 = pd.read_csv(filename1)
    elif dm=="tab":
       file1 = pd.read_csv(filename1, delimiter='\t')
    elif dm=='sem':
       file1 = pd.read_csv(filename1, delimiter=';')
    
def smitosmi(smi):
    mol=Chem.MolFromSmiles(smi)
    remover = SaltRemover.SaltRemover();
    m = remover.StripMol(mol)
    canSmi= Chem.MolToSmiles(m)
    if len(canSmi.split('.'))>1:
       f=len(canSmi.split('.')[0])
       l=len(canSmi.split('.')[1])
       if f>l:
          canSmi=canSmi.split('.')[0]
       else:
          canSmi=canSmi.split('.')[1]
    elif len(canSmi.split('.'))==1:
         pass
    return canSmi

def smitosdf():
    df=file1
    ls3,ls4=[],[]
    df['CanSmi']=df.apply(lambda x:smitosmi(x.SMILES), axis=1)
    if len(df['CanSmi'].unique().tolist())<len(df['SMILES'].tolist()):
       print("Duplicate molecules were found")
       messagebox.showinfo('Warning', 'Duplicate structures are present')
    else:
       print("No duplicate molecule was found")
    ls3=list(set(df['CanSmi']))
    #print(len(ls3))
    ls4=list(df['CanSmi'])
    #print(len(ls4))
    x,y=[],[]
    for i in ls4:
        if i not in x:
           x.append(i)
        else:
           y.append(i)
    #print(len(x),len(y))
    df['DUPLICATE']=df.apply(lambda x:'Duplicate' if x['CanSmi'] in y else 'Not Duplicate', axis=1)
    df.to_csv(str(c_)+'_processed.csv', index=False)
    mols,ls=[],[]
    for smi in df.CanSmi:
        try:
          mols.append(Chem.MolFromSmiles(smi))
        except:
          print(smi)
          continue
    hmols = [Chem.AddHs(m) for m in mols]
    for mol  in hmols:
        AllChem.EmbedMolecule(mol,AllChem.ETKDG())
        try:
           #AllChem.MMFFOptimizeMolecule(mol)
           AllChem.UFFOptimizeMolecule(mol,int(sEntryTabFive_1.get()))
        except ValueError:
           ls.append(Chem.MolToSmiles(mol))
    smiles = list(df.CanSmi)
    #sid = list(df.SOURCE_ID)
    libs = df[df.columns[1]]
    act= df[df.columns[2]]
    smiles=df[df.columns[0]]
    #an=act.columns
    writer = Chem.SDWriter(str(c_)+'_minimized'+'.sdf')
    f=pd.DataFrame(ls, columns=['failed'])
    f.to_csv('failed_SMILES.csv')

    for n in range(len(libs)):
        hmols[n].SetProp("_Name","%s"%libs[n])
        #hmols[n].SetProp("_Name","%s"%sid[n])
        hmols[n].SetProp(df.columns[2],"%s"%act[n])
        hmols[n].SetProp("Can_SMILES","%s"%smiles[n])
        writer.write(hmols[n])
    writer.close()
    align=Criterion.get()
    if align:
       mols = [m for m in Chem.SDMolSupplier(str(c_)+'_minimized'+'.sdf') if m != None]
       aligned_molecules, crippen_score=align_set_of_ligands(mols, int(aEntryTabFive_1.get()), int(aEntryTabFive_2.get()))
       outn=str(c_)+'_aligned'+'.sdf'
       writer = Chem.SDWriter(outn)
       for n in range(len(aligned_molecules)):
           aligned_molecules[n].SetProp(df.columns[2],"%s"%act[n])
           writer.write(aligned_molecules[n])
       writer.close()
    else:
       pass

    conf=Criterion2.get()
    if conf:
       writer = Chem.SDWriter(str(c_)+'_rdkitConf.sdf')
       process(filename1,str(c_)+'_rdkitConf.sdf',int(seventhEntryTabFive_1.get()),int(seventhEntryTabFive_2.get()),float(seventhEntryTabFive_3.get()),int(seventhEntryTabFive_4.get()),writer)
    else:
       pass

def disable_clvar():
    seventhEntryTabFive_1['state']='disabled'
    seventhEntryTabFive_2['state']='disabled'
    seventhEntryTabFive_3['state']='disabled'
    seventhEntryTabFive_4['state']='disabled'

def enable_clvar():
    seventhEntryTabFive_1['state']='normal'
    seventhEntryTabFive_2['state']='normal'
    seventhEntryTabFive_3['state']='normal'
    seventhEntryTabFive_4['state']='normal'

def disable_al():
    aEntryTabFive_1['state']='disabled'
    aEntryTabFive_2['state']='disabled'
    

def enable_al():
    aEntryTabFive_1['state']='normal'
    aEntryTabFive_2['state']='normal'


Criterion_Label0 = ttk.Label(tab1, text="File delimiter",font=("Helvetica", 12),anchor=W, justify=LEFT)
Criterion0 = StringVar()
Criterion0.set('com')
Criterion_xc = ttk.Radiobutton(tab1, text='comma', variable=Criterion0, value='com')
Criterion_xt = ttk.Radiobutton(tab1, text='tab', variable=Criterion0, value='tab')
Criterion_xs = ttk.Radiobutton(tab1, text='Semicolon', variable=Criterion0, value='sem')
Criterion_Label0.place(x=150,y=25)
Criterion_xc.place(x=250,y=25)
Criterion_xt.place(x=330,y=25)
Criterion_xs.place(x=380,y=25)

firstLabelTabOne = tk.Label(tab1, text="Select file",font=("Helvetica", 12))
firstEntryTabOne = tk.Entry(tab1,text='',width=50)
firstLabelTabOne.place(x=110,y=70)
firstEntryTabOne.place(x=200,y=73, height=25)
b5=tk.Button(tab1,text='Browse', command=datatr,font=("Helvetica", 10))
b5.place(x=520,y=75)


sLabelTabFive_1 = tk.Label(tab1, text="Minimization steps",font=("Helvetica", 12),anchor=W, justify=LEFT)
sLabelTabFive_1.place(x=200,y=110)
sEntryTabFive_1 = tk.Entry(tab1)
sEntryTabFive_1.place(x=350,y=113, width=50)


Criterion_Label = ttk.Label(tab1, text="Alignment",font=("Helvetica", 12))
Criterion = BooleanVar()
Criterion.set(False)
Criterion_Gini = ttk.Radiobutton(tab1, text='Yes', variable=Criterion, value=True, command=enable_al)
Criterion_Entropy = ttk.Radiobutton(tab1, text='No', variable=Criterion, value=False, command=disable_al)
Criterion_Label.place(x=240,y=150)
Criterion_Gini.place(x=320,y=150)
Criterion_Entropy.place(x=370,y=150)

aLabelTabFive_1 = tk.Label(tab1, text="Number of conformers",font=("Helvetica", 12),anchor=W, justify=LEFT)
aLabelTabFive_1.place(x=20,y=180)
aEntryTabFive_1 = tk.Entry(tab1, state='disabled', textvariable=IntVar(tab1, value=100))
aEntryTabFive_1.place(x=200,y=183, width=50)

aLabelTabFive_2 = tk.Label(tab1, text="Seed",font=("Helvetica", 12),anchor=W, justify=LEFT)
aLabelTabFive_2.place(x=400,y=180)
aEntryTabFive_2 = tk.Entry(tab1, state='disabled', textvariable=IntVar(tab1, value=42))
aEntryTabFive_2.place(x=470,y=183, width=50)



Criterion_Label2 = ttk.Label(tab1, text="Conformer generation",font=("Helvetica", 12))
Criterion2 = BooleanVar()
Criterion2.set(False)
Criterion_Gini2 = ttk.Radiobutton(tab1, text='Yes', variable=Criterion2, value=True, command=enable_clvar)
Criterion_Entropy2 = ttk.Radiobutton(tab1, text='No', variable=Criterion2, value=False, command=disable_clvar)
Criterion_Label2.place(x=150,y=230)
Criterion_Gini2.place(x=320,y=230)
Criterion_Entropy2.place(x=370,y=230)

seventhLabelTabFive_1 = tk.Label(tab1, text="Number of conformers",font=("Helvetica", 12),anchor=W, justify=LEFT)
seventhLabelTabFive_1.place(x=20,y=260)
seventhEntryTabFive_1 = tk.Entry(tab1, state='disabled', textvariable=IntVar(tab1, value=100))
seventhEntryTabFive_1.place(x=200,y=263, width=50)

seventhLabelTabFive_2 = tk.Label(tab1, text="Energy",font=("Helvetica", 12),anchor=W, justify=LEFT)
seventhLabelTabFive_2.place(x=260,y=260)
seventhEntryTabFive_2 = tk.Entry(tab1, state='disabled', textvariable=IntVar(tab1, value=50))
seventhEntryTabFive_2.place(x=320,y=263, width=50)

seventhLabelTabFive_3 = tk.Label(tab1, text="RMS",font=("Helvetica", 12),anchor=W, justify=LEFT)
seventhLabelTabFive_3.place(x=390,y=260)
seventhEntryTabFive_3 = tk.Entry(tab1,state='disabled', textvariable=StringVar(tab1, value=0.5))
seventhEntryTabFive_3.place(x=440,y=263, width=50)

seventhLabelTabFive_4 = tk.Label(tab1, text="Seed",font=("Helvetica", 12),anchor=W, justify=LEFT)
seventhLabelTabFive_4.place(x=520,y=260)
seventhEntryTabFive_4 = tk.Entry(tab1, state='disabled', textvariable=IntVar(tab1, value=42))
seventhEntryTabFive_4.place(x=570,y=263, width=50)


b2=Button(tab1, text='Generate datasets', command=smitosdf,bg="green",font=("Helvetica", 10),anchor=W, justify=LEFT)
b2.place(x=310,y=300)

tab_parent.pack(expand=1, fill='both')

form.mainloop()

