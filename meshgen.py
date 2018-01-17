import numpy as np

class empty:pass
##
## PARAMS
##
inclusionLineMarker = 101
inclusionMarker = 88
outerMarker = 99


# In[142]:


# In[142]:




# In[143]:

def inclusion(loc,width,ptNum=5): 
    if ptNum<5:
        raise RuntimeErorr("ptNum must be >=5")
    
    entry = ""
    margin = 0.05
    # Inclusion
    entry+="// Inclusion \n"
    entry+="Point(%d) = {%f,%f,0,lc};\n"%(ptNum, loc[0]-(width-margin),loc[1]-(width-margin)); 
    entry+="Point(%d) = {%f,%f,0,lc};\n"%(ptNum+1, loc[0]+(width-margin),loc[1]-(width-margin)); 
    entry+="Point(%d) = {%f,%f,0,lc};\n"%(ptNum+2, loc[0]+(width-margin),loc[1]+(width-margin)); 
    entry+="Point(%d) = {%f,%f,0,lc};\n"%(ptNum+3, loc[0]-(width-margin),loc[1]+(width-margin)); 
    entry+="Line(%d)  = {%d,%d};\n"      %(ptNum,ptNum,ptNum+1); 
    entry+="Line(%d)  = {%d,%d};\n"      %(ptNum+1,ptNum+1,ptNum+2); 
    entry+="Line(%d)  = {%d,%d};\n"      %(ptNum+2,ptNum+2,ptNum+3); 
    entry+="Line(%d)  = {%d,%d};\n"      %(ptNum+3,ptNum+3,ptNum); 
    entry+="Line Loop(%d) = {%d,%d,%d,%d};\n"  %(ptNum,ptNum,ptNum+1,ptNum+2,ptNum+3)
    entry+="Plane Surface(%d) = {%d};\n"   % (ptNum,ptNum)
    # boundary
    ptNumOuter = ptNum+4
    entry+="// Boundary \n"
    entry+="Point(%d) = {%f,%f,0,lc};\n"%(ptNumOuter, loc[0]-(width),loc[1]-(width)); 
    entry+="Point(%d) = {%f,%f,0,lc};\n"%(ptNumOuter+1, loc[0]+(width),loc[1]-(width)); 
    entry+="Point(%d) = {%f,%f,0,lc};\n"%(ptNumOuter+2, loc[0]+(width),loc[1]+(width)); 
    entry+="Point(%d) = {%f,%f,0,lc};\n"%(ptNumOuter+3, loc[0]-(width),loc[1]+(width)); 
    entry+="Line(%d)  = {%d,%d};\n"      %(ptNumOuter,ptNumOuter,ptNumOuter+1);
    entry+="Line(%d)  = {%d,%d};\n"      %(ptNumOuter+1,ptNumOuter+1,ptNumOuter+2);
    entry+="Line(%d)  = {%d,%d};\n"      %(ptNumOuter+2,ptNumOuter+2,ptNumOuter+3);
    entry+="Line(%d)  = {%d,%d};\n"      %(ptNumOuter+3,ptNumOuter+3,ptNumOuter);
    entry+="Line Loop(%d) = {%d,%d,%d,%d};\n"  %(ptNumOuter,ptNumOuter,ptNumOuter+1,ptNumOuter+2,ptNumOuter+3)
    
    newPtNum = ptNum + 8
    inclusionPhysicalSurface = ptNum
    boundaryLineLoop =ptNumOuter
    return entry, inclusionPhysicalSurface,boundaryLineLoop,newPtNum    
    


# Function to generate a gmsh-based mesh given a set of inclusion coordinates
# and widths
def GenMany(
  locs=np.array([[0.5,0.5],[2.5,2.5],[4.5,4.5],[6.5,6.5]]),
  widthInclusion = 0.5,# units?
  widthBox = 10., # 
  filePrefix = "test",
  outerRes = 0.5,
  inclRes = 0.5
  ): 
  ## Define outer box
  header="""
lc = %f; 
lco =%f; 
// outer
Point(1) = {%f, %f, -0, lco};
Point(2) = { 0, %f, -0, lco};
Point(3) = { 0,  0, -0, lco};
Point(4) = {%f,  0, -0, lco};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(1) = {3, 4, 1, 2};
// 
"""%(inclRes,outerRes,widthBox,widthBox,widthBox,widthBox)

  nLocs = np.shape(locs)[0]

  # loop over coordinates 
  entries = []
  ptNum = 5 # Always start from 5, since outer box covers pts 1-4 
  for i in np.arange(nLocs):
    loc = locs[i,:]
    print loc
    
    entry = empty()
    entry.line,entry.inclusionPhysicalSurface,entry.boundaryLineLoop,newPtNum = inclusion(loc=loc,width=widthInclusion,ptNum=ptNum)
    entries.append(entry)
    ptNum = newPtNum


  # define outer line loops 
  outerLineLoops="Plane Surface(17) = {1,"
  for entry in entries:
      outerLineLoops+="%d, "%entry.boundaryLineLoop    
  outerLineLoops+="};\n";  outerLineLoops = outerLineLoops.replace(", }","}")
  
  # physical line (placeholder)
  #inclusionPhysicalLine="//inclusion line\nPhysical Line(%d) = {%d};\n" % (inclusionLineMarker,entries[0].inclusionPhysicalSurface)
  inclusionPhysicalLine="Physical Line(%d) = {5,6,7,8};\n"%inclusionLineMarker
  
  inclusionPhysicalSurface="//inclusion surface\nPhysical Surface(%d) = {" % inclusionMarker
  for entry in entries:
      #print entry.inclusionLineLoop
      inclusionPhysicalSurface+="%d, "%entry.inclusionPhysicalSurface
      #print entry.boundaryLineLoop
  inclusionPhysicalSurface+="};\n";  inclusionPhysicalSurface = inclusionPhysicalSurface.replace(", }","}")
  
  
  outerPhysicalSurface="""
//outer surface
Physical Surface(%d) = {17};
  
  """ % outerMarker
  
  
  #print outerLineLoops
  #print inclusionPhysicalSurface
  #print outerPhysicalSurface
  
  
  # In[148]:
  
  fileLines = []
  fileLines.append(header)
  
  for entry in entries:
      fileLines.append(entry.line)
      
  
  fileLines.append(outerLineLoops)    
  fileLines.append(inclusionPhysicalLine)    
  fileLines.append(inclusionPhysicalSurface)    
  fileLines.append(outerPhysicalSurface) 
  
  
  # In[149]:
  
  fileName = filePrefix+".geo"
  with open(fileName, "w") as text_file:
      for fileLine in fileLines:
        text_file.write(fileLine)
  
  
  # In[150]:
  
  import os
  mshName = fileName.replace(".geo",".msh")
  xmlName = fileName.replace(".geo",".xml")
  line = "gmsh -2 %s"%fileName
  os.system(line)
  line="dolfin-convert %s %s"%(mshName,xmlName)
  os.system(line)
  return xmlName 
  

