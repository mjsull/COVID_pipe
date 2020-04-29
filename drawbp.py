import matplotlib.pyplot as plt
import sys
import random

zeros = 0
first = True
for zefile in sys.argv[2:]:
    with open(zefile) as f:
        if first:
            primers = f.readline().split()[1:]
            columns = len(primers)
            y = [[] for j in range(columns)]
            y_ave = [[] for j in range(columns)]
            first = False
            xscatters = []
            yscatters = []
            yscattersave = []
        else:
            f.readline()
        xscatter = []
        yscatter = []
        yscatterave = []
        for line in f:
            thesum, thenum = 0, 0
            for i in line.split()[1:]:
                thesum += float(i)
                thenum += 1
            average = thesum /thenum
            if average == 0.0:
                zeros += 1
                continue
            for num, i in enumerate(line.split()[1:]):
                y[num].append(float(i))
                y_ave[num].append(float(i)/average)
                yscatter.append(float(i))
                yscatterave.append(float(i)/average)
                xscatter.append(random.random()/2 -0.25+ num)
        xscatters.append(xscatter)
        yscatters.append(yscatter)
        yscattersave.append(yscatterave)
    x = [q for q in range(columns)]


fig, ax = plt.subplots()
markers = ['r.', 'bx', 'g+']
if sys.argv[1] == 'normal':
    bp = ax.boxplot(y_ave, positions=x, showfliers=False)
    ax.set_ylabel('median coverage (normalized)')
    for num in range(len(xscatters))    :
        xscatter = xscatters[num]
        yscatterave = yscattersave[num]
        marker = markers[num%len(markers)]
        ax.plot(xscatter, yscatterave, marker, alpha=0.2)
else:
    bp = ax.boxplot(y, positions=x)
    ax.set_ylabel('median coverage')
    ax.plot(xscatter, yscatter, 'r.', alpha=0.2)
ax.set_xlabel('primer')
plt.setp(bp['fliers'], marker='x', markersize=1)
ax.set_xticklabels(primers,
                    rotation=45, fontsize=8)
plt.show()