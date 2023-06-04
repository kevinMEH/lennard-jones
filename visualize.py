import matplotlib.pyplot as plot

coordinates = """
0.317950669,-0.143832896,-0.636066046
0.695342839,-1.894777575,-0.456724788
-0.120408848,-0.513323092,0.379616912
0.805494464,-0.895830904,-0.044332331
0.823533386,-0.834989673,1.069841777
0.825040667,0.036071666,0.394663288
0.844941126,-0.975138854,-1.103523884
-0.165130917,-1.171786295,-0.540786627
0.166226685,-1.610719077,0.456690264
1.694208684,-1.352822481,-0.487363404
1.299489243,-1.658152059,0.505338560
1.768206560,-0.656905353,0.441433776
1.439632887,-0.175852916,-0.548569233
"""

coordinates = coordinates.strip().split("\n")
coordinates = [ line.split(",") for line in coordinates ]
coordinates = [ (float(line[0]), float(line[1]), float(line[2])) for line in coordinates ]

figure, axes = plot.subplots(subplot_kw={"projection": "3d"})

for x, y, z in coordinates:
    axes.scatter([x], [y], [z], c="tab:blue")
    for x2, y2, z2 in coordinates:
        distance2 = (x - x2)**2 + (y - y2)**2 + (z - z2)**2
        if distance2 < 1.5:
            axes.plot([x, x2], [y, y2], [z, z2], c="tab:blue")

axes.set_xlabel("x")
axes.set_ylabel("y")
axes.set_zlabel("z")

plot.show()
