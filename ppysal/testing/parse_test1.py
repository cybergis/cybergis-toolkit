import numpy as np
import collections

d= {}
with open('test1.log', 'r') as f:
    for l in f.readlines():
        i = l.split(' ')
        if len(i) == 1:
            continue
        else:
            node = i[4]
            time = i[-2]
            d[node] = time
od = collections.OrderedDict(sorted(d.items()))
with open('test1.nodelist', 'r') as f:
    nodes = f.readlines()

live_nodes = set([n.rstrip('\n') for n in nodes])
tested_nodes = set([k for k in od.iterkeys()])

print
print "######WARNING######"
print "The following nodes were either intentionally omitted or never responded:"
for n in live_nodes.difference(tested_nodes):
    print n
print

times = [float(v)/8.0 for v in od.itervalues()]
meantime = np.mean(times)
stdtime = np.std(times)

for k, v in od.iteritems():
    #z-test
    v = float(v) / 8.0
    if  (v-meantime) / stdtime >= 1.96:
        print "Node {} required significantly more processing time than expected, using {} seconds.".format(k,v)
