20171219

1.Find that we can use the VTK to read into the NIFITY label image, and run successfully
2.Find that the coordinate of the mesh is the physical coordinate.


20171220
1.Use these coordinates to segment the scar and the nodeID change into the original ID (not the mesh ID), but the ID seems not be segmented.

2.The z axis seems be turned, so I changed it be the "size[2]-k"


20171220
1.Change the t-link weight with the wright version
2.Change the n-link weight to keep the currEdgeStrength in the range of "0~1".This is very important!
So when next days, I'll train the n-link in the range of "0~1", which is not difficult.

20180105
1.Remove the PV and Mitral valve
2.Using the Gold truth to generate n-link.

20180123
1. tlink changed into "hard constrain"
2. Add gaussian noise in all the prob (* 0<=prob<=1)
3. Add one cut and using Additional edges
