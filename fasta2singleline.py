#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 12:41:23 2017

@author: rachel
"""

infile=open("h37ra.fasta",'r')
out=open("h37ra.txt",'w')
for line in infile:
    line = line.strip("\n")
    out.write('%s' %(line))
    
#infile=open("genome-test.txt", 'r')
#out=open("g-4t.txt",'w')
#line=line.strip("TTTT")