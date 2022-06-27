#!/usr/bin/env python
# coding: utf-8

# In[1]:


#THE PROPAGATION AND STABILITY OF THE DECOMPOSITION FLAME OF HYDRAZINE


# In[2]:


import cantera as ct #run this with CANTERA 2.5.0
import numpy as np
import pandas as pd

print(f"Running Cantera Version: {ct.__version__}")


# In[3]:


# Import plotting modules and define plotting preference
get_ipython().run_line_magic('config', 'InlineBackend.figure_formats = ["svg"]')
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pylab as plt

plt.rcParams["axes.labelsize"] = 14
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 10
plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["figure.dpi"] = 120

# Get the best of both ggplot and seaborn
plt.style.use("ggplot")
plt.style.use("seaborn-deep")

plt.rcParams["figure.autolayout"] = True


# In[4]:


# Inlet Temperature in Kelvin and Inlet Pressure in Pascals
# In this case we are setting the inlet T and P to room temperature conditions
To = 335
Po = 10666


# Define the gas-mixutre and kinetics
# In this case, we are choosing a GRI3.0 gas
gas = ct.Solution("/home/hydrazine/runs/xh1005/cantera/chem_annotated.cti")

# Create a stoichiometric CH4/Air premixed mixture
gas.TPX = To, Po, {"H4N2(1)": 1.0}


# In[5]:


# Domain width in metres
width = 0.014

# Create the flame object
flame = ct.FreeFlame(gas, width=width)

# Define tolerances for the solver
flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)

# Define logging level
loglevel = 1


# In[6]:


flame.flame.set_transient_tolerances(default=(1e-15, 1e-20))


# In[7]:


flame.solve(loglevel=loglevel, auto=True)

Su0 = flame.u[0]
print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")


# In[8]:


for i, specie in enumerate(gas.species()):
    print(f"{i}. {specie}")


# In[9]:


# Extract concentration data
X_H4N2 = flame.X[2]
X_NH3 = flame.X[8]
X_N2 = flame.X[0]
X_H2 = flame.X[4]

plt.figure()

plt.plot(flame.grid * 100, X_H4N2, "-o", label="$H4N2_{4}$")
plt.plot(flame.grid * 100, X_NH3, "-s", label="$NH3_{2}$")
plt.plot(flame.grid * 100, X_N2 , "-<", label="$N2_{2}O$")
plt.plot(flame.grid * 100, X_H2, "-o", label="$H2_{4}O$")

plt.legend(loc=2)
plt.xlabel("Distance (cm)")
plt.ylabel("MoleFractions");


# In[10]:


# Create a dataframe to store sensitivity-analysis data
sensitivities = pd.DataFrame(data=[], index=gas.reaction_equations(range(gas.n_reactions)))


# In[11]:


# Set the value of the perturbation
dk = 1e-2


# In[14]:


# Create an empty column to store the sensitivities data
sensitivities["baseCase"] = ""
for m in range(gas.n_reactions):
    gas.set_multiplier(1.0) # reset all multipliers
    gas.set_multiplier(1+dk, m) # perturb reaction m    
    
    # Always force loglevel=0 for this
    # Make sure the grid is not refined, otherwise it won't strictly
    # be a small perturbation analysis
    flame.solve(loglevel=0, refine_grid=False)    
    
    # The new flame speed
    Su = flame.velocity[1]    
    
    sensitivities["baseCase"][m] = (Su-Su0)/(Su0*dk)

# This step is essential, otherwise the mechanism will have been altered
gas.set_multiplier(1.0)
sensitivities.head()
# Reaction mechanisms can contains thousands of elementary steps. Choose a threshold
# to see only the top few
threshold = 0.546
firstColumn = sensitivities.columns[0]


# In[26]:


# Reaction mechanisms can contains thousands of elementary steps. Choose a threshold
# to see only the top few
threshold = 0.03

# For plotting, collect only those steps that are above the threshold
# Otherwise, the y-axis gets crowded and illegible
sensitivities_subset = sensitivities[sensitivities["baseCase"].abs() > threshold]
reactions_above_threshold = (
    sensitivities_subset.abs().sort_values(by="baseCase", ascending=False).index
)
sensitivities_subset.loc[reactions_above_threshold].plot.barh(
    title="Sensitivities for Hydrazine", legend=None
)
plt.gca().invert_yaxis()

plt.rcParams.update({"axes.labelsize": 20})
plt.xlabel(r"Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$");

# Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
# plt.savefig('sensitivityPlot', dpi=300)








