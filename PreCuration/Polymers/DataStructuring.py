import pandas as pd
import numpy as np 

import pandas as pd
import os


PATH = os.getcwd() +'Polymers/'
# this code will take in our  unstructured data and make it structured for curation

# read in the unstructured data 

df=pd.read_csv('PolymerData_uncurated_unstructured.csv')


# functions for gahtering data and structuring 

# trying to make a function so that all I have to do is input the name



def settle_property(prop):
    a = prop[0]
    b=prop[1]
    if b is not '':
        return float(b)
    elif a is not '':
         if '-' in a:
            c,d= a.split('-')
            c = c.strip()
            c = float(c)
            d = d.strip()
            d = float (d)
            avg = np.average(c,d)
            return avg
         else:
              return float(a)
    else:
         return None



def gather_prop(name, source):
    holder = []
    for j in df[source]:
    
        try:

                j = str(j)
                
                data = j.replace(":",",")
                data = data.replace("{","[")
                data = data.replace("}","]")

                Data = eval(data)
                for k in Data:
                    if k[0] == name:
                        f = settle_property(k[2:])
                        holder.append(f)

                    else:
                        next
        except:
             holder.append(None)

    return holder


# run fucntions for experimental properties 
# gathering experimental properties
Exp_properties=['Molar Volume Vm',
 'Density ρ',
 'Solubility Parameter δ',
 'Molar Cohesive Energy Ecoh',
 'Glass Transition Temperature Tg',
 'Molar Heat Capacity Cp',
 'Entanglement Molecular Weight Me',
 'Index of Refraction n']



Molar_volume  =gather_prop(Exp_properties[0],'Thermo-Physical Properties: Experimental / Literature Data')
Density = gather_prop(Exp_properties[1],'Thermo-Physical Properties: Experimental / Literature Data')
Solubility_Parameter = gather_prop(Exp_properties[2],'Thermo-Physical Properties: Experimental / Literature Data')
Molar_cohesive_E = gather_prop(Exp_properties[3],'Thermo-Physical Properties: Experimental / Literature Data')
Glass_trans = gather_prop(Exp_properties[4],'Thermo-Physical Properties: Experimental / Literature Data')
Molar_heat_cap = gather_prop( Exp_properties[5],'Thermo-Physical Properties: Experimental / Literature Data')
Entaglemeent_MolWeight =gather_prop( Exp_properties[6],'Thermo-Physical Properties: Experimental / Literature Data')
N =gather_prop(Exp_properties[7],'Thermo-Physical Properties: Experimental / Literature Data')

# now for computed properties 


# gather computed properties
CompProperties = [
 'Molecular Weight of Repeat unit',
 'Van-der-Waals Volume VvW',
 'Molar Volume Vm',
 'Density ρ',
 'Solubility Parameter δ',
 'Molar Cohesive Energy Ecoh',
 'Glass Transition Temperature Tg',
 'Molar Heat Capacity Cp',
 'Entanglement Molecular Weight Me',
 'Index of Refraction n']




Molar_Weight_Repeat  =gather_prop(CompProperties[0],'Thermo-Physical Properties: Calculated Data')
VanDerWaal =gather_prop(CompProperties[1],'Thermo-Physical Properties: Calculated Data')
MolarVolume=gather_prop(CompProperties[2],'Thermo-Physical Properties: Calculated Data')
Density_Comp=gather_prop(CompProperties[3],'Thermo-Physical Properties: Calculated Data')
Solubility_Parameter_comp=gather_prop(CompProperties[4],'Thermo-Physical Properties: Calculated Data')
Molar_cohesive_E_comp=gather_prop(CompProperties[5],'Thermo-Physical Properties: Calculated Data')
Glass_trans_comp=gather_prop(CompProperties[6],'Thermo-Physical Properties: Calculated Data')
Molar_heat_cap_comp=gather_prop(CompProperties[7],'Thermo-Physical Properties: Calculated Data')
Entaglemeent_MolWeight_comp=gather_prop(CompProperties[8],'Thermo-Physical Properties: Calculated Data')
N_comp=gather_prop(CompProperties[9],'Thermo-Physical Properties: Calculated Data')


#for the identifiesrs and such 

def gather_IDs(name, source):
    holder = []
    for j in df[source]:
    
        try:

            j = str(j)
            
            data = j.replace(":",",")
            data = data.replace("{","[")
            data = data.replace("}","]")

            Data = eval(data)
           
 
            for k in Data:
                if k[0] == name:
                    
                    holder.append(k[1:])
        except:
            holder.append(None)



    return holder


#finally, get ID information

Pol_ID_list = [
'POLYMER CLASS',
 'COMMON NAMES',
 'STRUCTURE BASED NAME',
 'CAS #',
 'PSMILES']

Polymer_Class  =gather_IDs(Pol_ID_list[0],'NAMES AND IDENTIFIERS OF POLYMER')
Common_names=gather_IDs(Pol_ID_list[1],'NAMES AND IDENTIFIERS OF POLYMER')
Structural_name=gather_IDs(Pol_ID_list[2],'NAMES AND IDENTIFIERS OF POLYMER')
CAS=gather_IDs(Pol_ID_list[3],'NAMES AND IDENTIFIERS OF POLYMER')
PSMILES=gather_IDs(Pol_ID_list[4],'NAMES AND IDENTIFIERS OF POLYMER')


Mon_ID_list=[
 'CAS #',
 'COMMON NAMES',
 'SMILES',
 'Std. InChI',
 'Std. InChIKey',
 'STRUCTURE'
]

CAS_Monomer =gather_IDs(Mon_ID_list[0],'IDENTIFIERS OF MONOMER(S)')
Common_names_Monomer =gather_IDs(Mon_ID_list[1],'IDENTIFIERS OF MONOMER(S)')
SMILES_Monomer=gather_IDs(Mon_ID_list[2],'IDENTIFIERS OF MONOMER(S)')
Std_InChI_Monomer=gather_IDs(Mon_ID_list[3],'IDENTIFIERS OF MONOMER(S)')
Std_InChIKey_Monomer=gather_IDs(Mon_ID_list[4],'IDENTIFIERS OF MONOMER(S)')


# linking the data together
# try with dataframe

data = [CAS_Monomer, 
Common_names_Monomer, 
SMILES_Monomer,
Std_InChI_Monomer,
Std_InChIKey_Monomer,

Polymer_Class,
Common_names,
Structural_name,
CAS,
PSMILES,

Molar_Weight_Repeat, 
VanDerWaal, 
MolarVolume,
Density_Comp,
Solubility_Parameter_comp,
Molar_cohesive_E_comp,
Glass_trans_comp,
Molar_heat_cap_comp,
Entaglemeent_MolWeight_comp,
N_comp,

Molar_volume,
Density, 
Solubility_Parameter,
Molar_cohesive_E,
Glass_trans,
Molar_heat_cap,
Entaglemeent_MolWeight,
N
]

cols=[
    'CAS_Monomer',
    'Common_names_Monomer',
    'SMILES_Monomer',
    'Std_InChI_Monomer',
    'Std_InChIKey_Monomer',
    'Polymer_Class',
    
    'Common_names',
    'Structural_name',
    'CAS',
    'PSMILES',
    'Molar_Weight_Repeat',
    'VanDerWaal',
    'MolarVolume',
    'Density_Comp',
    'Solubility_Parameter_comp',
    'Molar_cohesive_E_comp',
    'Glass_trans_comp',
    'Molar_heat_cap_comp',
    'Entaglemeent_MolWeight_comp',
    'N_comp',
    'Molar_volume',
    'Density',
    'Solubility_Parameter',
    'Molar_cohesive_E',
    'Glass_trans',
    'Molar_heat_cap',
    'Entaglemeent_MolWeight',
    'N'

]


df1 = pd.DataFrame(data)

df1=df1.T

df1.columns=cols

# fix the psmiles

psmiles = df1['PSMILES']


res = []
for i in psmiles:
    if i is not None:
        s = i[0]
        res.append(s)
    else:
        res.append(None)


df1['PSMILES'] = df['PSMILES']

df_final=pd.concat([df,df1], keys=['PSMILES'])

# save to Precuration folder 
df_final.to_csv("PolymerData_uncurated.csv",index=False)

#save to DATA FOLDER
df_final.to_csv(os.path.join('..', '..')+"/Data/Polymers/PolymerData_uncurated.csv",index=False)
