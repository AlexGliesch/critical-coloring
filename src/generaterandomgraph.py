/*
* A new heuristic for finding verifiable k-vertex-critical subgraphs
* 
* Copyright (c) 2022 Alex Gliesch, Marcus Ritt
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPY lRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/
#!/usr/bin/python

# generates uniform graphs in DIMACS format 

import sys, random

if len(sys.argv) != 4:
  print ('Usage: ./makerandomgraph.py n d s,\nwhere n>0 is the number of ' + 
         'vertices, p \in [0,1] is the probability of each edge existing, ' + 
         'and s is a random seed. ' +
         'Parameters out of range will cause a crash.')
  exit()

n = int(sys.argv[1])
d = float(sys.argv[2])
s = int(sys.argv[3])
e = []
m = 0
random.seed(s)
for v in xrange(n):
  for w in xrange(v+1,n):
    if random.uniform(0.0,1.0) < d:
      e.append((v,w))
      m = m+1
print 'p edge {} {}'.format(n, m)
for i in e:
  print 'e {} {}'.format(i[0]+1, i[1]+1)
