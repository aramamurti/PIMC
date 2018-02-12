import sys

if(len(sys.argv) == 1):
    folder_path = 'parameters'
elif(len(sys.argv) == 2):
    folder_path = sys.argv[1]

for i in range(0,36):
    f = open(folder_path+"/parameters_{}.cfg".format(i),'w')
    f.write("[Simulation Parameters]\n")
    f.write(" particles                              =           32\n")
    f.write(" periodic_boundary_conditions           =           true\n")
    f.write(" grand_canonical_ensemble               =           false\n")
    f.write(" particle_type                          =           boson_coulomb\n")
    f.write(" temperature                            =           {}\n".format(1.4+.1*i))
    f.write(" mu                                     =           0\n")
    f.write(" NN_grid_size                           =           1\n")
    f.write(" C0                                     =           1\n")
    f.write(" end_step                               =           50000\n")
    f.write(" equilibration                          =           5000\n")
    f.write(" time_slices                            =           32\n")
    f.write(" coupling                               =           0\n")
    f.write(" dimensions                             =           3\n")
    f.write(" kb                                     =           1.0\n")
    f.close();

