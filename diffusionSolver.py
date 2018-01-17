import numpy as np 
import numpy 
from dolfin import *

class empty:pass

# In[179]:
boxDepth = 50
class Front(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 0)

class Back(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 50)


# In[192]:

def myfunc(
  xmlName,
  condInclusion,
  condBulk,
  voltageL,
  voltageR,
  outerBoxMin=[],
  outerBoxMax=[],
  lineMarkers = [],# [inclusionLineMarker],
  domainMarkers = [], #[inclusionMarker,outerMarker]
  ):

  # hack 
  inclusionMarker,outerMarker =domainMarkers[0],domainMarkers[1]

  # Create mesh and define function space
  name = xmlName.replace(".xml","")  
  mesh = Mesh(name+".xml")
  subdomains = MeshFunction("size_t", mesh, name+"_physical_region.xml")
  boundaries = MeshFunction("size_t", mesh, name+"_facet_region.xml")

    
    
  V = FunctionSpace(mesh, "Lagrange", 1)

    
  # Define boundary condition
  frontBoundary = Front()
  frontBoundary.outerBoxMin = outerBoxMin   
  frontmarker = 2
  frontBoundary.mark(boundaries,frontmarker)
  backBoundary = Back()
  backBoundary.outerBoxMax = outerBoxMax  
  backmarker = 3
  backBoundary.mark(boundaries,backmarker)

  ## Mark domains according to their 'boundaries definitions' 
  dx = Measure("dx",subdomain_data=subdomains,domain=mesh)
  ds = Measure("ds",subdomain_data=boundaries,domain=mesh)
  
  totalArea = assemble(Constant(1.)*ds)  
  #print "Total Area ", totalArea
  for i, marker in enumerate(lineMarkers):
    area = assemble(Constant(1.)*ds(marker))
    #print " Area ", area, " for marker ", marker

  totalVol = assemble(Constant(1.)*dx)
  #print "Total Vol ", totalVol
  for i, marker in enumerate(domainMarkers):
    vol = assemble(Constant(1.)*dx(marker))
   # print " Vol ", vol, " for marker ", marker
    #print i
    #print marker

  # Define input data
  a0 = Constant(condInclusion)
  a1 = Constant(condBulk)
  g_L = Constant(voltageL) #mV
  g_R = Constant(voltageR) #mV

  dx = Measure("dx",subdomain_data=subdomains,domain=mesh)
  ds = Measure("ds",subdomain_data=boundaries,domain=mesh)
    
  bcs = []
  bcs.append(DirichletBC(V, g_L, boundaries,frontmarker))
  bcs.append(DirichletBC(V, g_R, boundaries,backmarker))
  
  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  
  f = Constant(0.)
  g = Constant(0.)
  form = a0 * inner(grad(u), grad(v))*dx(inclusionMarker, domain=mesh)
  form += a1 * inner(grad(u), grad(v))*dx(outerMarker, domain=mesh)
  form += g*v*ds #(frontmarker,domain=subdomains)
  
# Compute solution
  a = lhs(form)
  L = rhs(form)
  x = Function(V)
  solve(a == L, x, bcs)


  File("out.pvd") << x
  #import numpy as np 
  #print "L2 norm ", norm(x,"l2")

  results = empty()
  results.x = x
  results.dx = dx
  results.ds = ds
  results.mesh = mesh 
  
  n = FacetNormal(mesh)
  m1 = dot(grad(x), n)*ds(backmarker)  # evaluate  grad(u) on RHS 
  v1 = assemble(m1)
  sa = assemble(Constant(1.)*ds(backmarker,domain=mesh)) # get SA 
  jEff =  v1/sa
  results.sa = sa
  results.jEff = jEff


  xs = mesh.coordinates() # determine length of mesh along x direction
  #mins = np.min(xs,axis=0)
  #maxs = np.max(xs,axis=0)
  dist = 50 #maxs[0] - mins[0] # dist along x 
  voltageL = 100. # mV
  voltageR = 0. # mV
  gradV = (voltageL - voltageR)/dist
  
  sig = jEff / gradV
  sig *= -1 # since I think normals point into domain   
  results.sig = sig
  print results.sig
  return results

  
  


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

#
# MAIN routine executed when launching this script from command line 
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"


