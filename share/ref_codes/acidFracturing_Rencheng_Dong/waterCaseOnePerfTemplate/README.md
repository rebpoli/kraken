### This is a template case for running dissolFoam.
### Main features:
1. Pad fluid is water.
2. There is a perforation zone in the middle of the inlet boundary. The injection velocity is nonzero only in the perforation zone. 
3. The velocity profile in the perforation zone is uniform.
4. Since the inlet boundary is expanding along the fracture width direction, the velocity in the perforation zone is updated at each time step. Injection velocity = Injection volumetrc rate / Perforation zone surface area.
5. The perforation zone can be set by defining x coordinates of perforation ranges in Zero/U file.
### This case is ready to run with OpenFoam v1706 and dissolFoam installed. Procedures:
1. ./Mesh > log.Mesh & // create mesh
2. ./Run > log.Run & // run the case with single processor