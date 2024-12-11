import os

#distance = 0.05
det = "H"
distance_values = list(float(i)/1000 for i in range(5,151,2))

for distance in distance_values:
    inputdir = f"{os.getcwd()}/Frames{distance}"
    label = f"{det}1_SSSM_m1_0.012_m2_0.012_dL_{distance}_tc_1000023915"
    frametype = label
    tseg = 4096
    start_gpstime = 1000000000
    Ndata = 5

    framecache = f"{inputdir}/framecache{distance}"

    with open(framecache, "w+") as f:
        for k in range(Ndata):
            filename = f"file://localhost{inputdir}/{label}_{start_gpstime + k * tseg:d}-{tseg:d}.gwf"
            f.write(f"{det} {frametype} {start_gpstime + k*tseg} {tseg:d} {filename}\n")

    #print("Done: make framecache.")
