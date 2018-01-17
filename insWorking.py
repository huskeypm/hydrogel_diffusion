
#builds mesh using the insMesh.py


import numpy as np

# In[66]:

import MC
import miscutil as mu
import insMesh


def GenRandomizedMesh(numInclusions, widthInclusion, widthBox,\
    randomize=True,
    inclusionResolution=0.5,
    outerResolution = 1.,
    multiple = True):
    if randomize: 
      locs, rads = randomizePoints(numInclusions, widthInclusion, widthBox)
      # take first solution
      locs = locs[0]
    else: 
      locs=np.zeros((2,2)); 
      locs[0,:] = [5,5]
      locs[1,:] = [80,80]
      rads = np.zeros(2); rads[:] = 1.
    if multiple==True:
    # define name based on num inclus. 
    	filePrefix = "mesh_%d_%d_ins"%(numInclusions,widthBox) # added second flag to allow for variations in volFrac based on numInclusions and ratio of widthInclusion to widthBox
    	xmlName = insMesh.GenMany(
      	  locs=locs,    
      	  widthInclusion = widthInclusion,
      	  inclRes = inclusionResolution,
      	  outerRes = outerResolution,
      	  widthBox = widthBox, # nm
      	  filePrefix = filePrefix
      	  )
    else:       
 	filePrefix = "mesh_%d_ins"%(widthBox) 
        xmlName = insMesh.GenMany(
      	  locs=locs,
      	  widthInclusion = widthInclusion,
      	  inclRes = inclusionResolution,
      	  outerRes = outerResolution,
      	  widthBox = widthBox, # nm
      	  filePrefix = filePrefix
     	  )
    return xmlName 



#pass in a number of inclusions per side, how many inclusions expected per singular dimension
def randomizePoints(numInclusions,widthInclusion,widthBox): #numInclusions, widthInclusion(nm), widthBox(nm)
    length = widthBox
    ptsPerSide = False
    ptRad = widthInclusion*sqrt(2)
    keep = 3  # randomizations to keep 
    radScale = 0.  # keep radii fixed
    trials = 1e4
    comDist = np.float(length/numInclusions)
    passingLocs, passingRads = MC.GenLocsRads(
      length=length,
      ptsPerSide = ptsPerSide,
      comDist = comDist,
      ptRad = ptRad,
      keep = keep,  # randomizations to keep 
      trials = trials,
      radScale =radScale
    )
    return passingLocs, passingRads


# In[67]:



# In[68]:

from dolfin import *
import logging 
logging.getLogger("FFC").setLevel(logging.ERROR) # suppress errors



# In[69]:

def getCoords(xmlName):
    name = xmlName.replace(".xml","")
    mesh = Mesh(name+".xml")
    coords = mesh.coordinates()
    outerBoxMin, outerBoxMax = np.min(coords,axis=0), np.max(coords,axis=0)
    return name, mesh, coords, outerBoxMin, outerBoxMax


# In[72]:

inclusionLineMarker = 101
inclusionMarker = 88
outerMarker = 99



# In[17]:

#test




# In[77]:

import diffusionSolver as iS


# In[78]:

def getOutput(
  xmlName,
  outerBoxMin, 
  outerBoxMax,
  lineMarkers,
  domainMarkers,
  condInclusion,
  condBulk,
  voltageL,
  voltageR):
    
  output = iS.myfunc(
  xmlName=xmlName,
  outerBoxMin=outerBoxMin, 
  outerBoxMax=outerBoxMax,
  lineMarkers = [inclusionLineMarker],
  domainMarkers = [inclusionMarker,outerMarker],
  condInclusion = condInclusion,
  condBulk = condBulk,
  voltageL = voltageL,
  voltageR = voltageR)
  
  return output


# In[79]:



# ### Test

# In[80]:



# #### Compute J based on conc. gradient on RHS 

# In[81]:

# Find normals for all 'surfaces' 


# #### Compute grad(V) based on V and LHS and RHS over distance

# In[82]:



# #### get sigma eff 
# $ \sigma  = \frac{\langle J_x \rangle}{\nabla V} $

# In[83]:




