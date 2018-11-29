#!/usr/bin/env python

import sys

args = sys.argv

f = open(args[1])
of = open(args[2], "w")

keep = set(["gene", "transcript", "exon"])
for line in f:
    if line.startswith("#"):
        continue
    cols = line.split("\t")
    if cols[2] == "mRNA":
        cols[2] = "transcript"
    if cols[2] not in keep:
        continue

    attribs = cols[8].split(";")
    attribs = [x.replace("=", "\t").split("\t") for x in attribs]
    newAttribs = {}
    for t in attribs:
        if t[0] == "ID":
            _, ID = t[1].split(":")
            if cols[2] == "gene":
                newAttribs['gene_id'] = ID
            elif cols[2] == "transcript":
                newAttribs['transcript_id'] = ID
                x = ID.rindex(".")
                newAttribs['gene_id'] = ID[:x]
            elif cols[2] == "exon":
                x = ID.rindex(".")
                ID = ID[:x]
                newAttribs['transcript_id'] = ID
                x = ID.rindex(".")
                ID = ID[:x]
                newAttribs['gene_id'] = ID
    cols[8] = " ".join("{} \"{}\";".format(x, y) for x, y in newAttribs.items())
    of.write("\t".join(cols))
    of.write("\n")
f.close()
of.close()
