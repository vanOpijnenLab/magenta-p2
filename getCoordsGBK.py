#Margaret Antonio 15.10.18

#Use this script to get gene names and coordinates from a GBK file
#Output file can be inputted into Gviz as an annotation track for visualization

from Bio import SeqIO
import sys

gbkFile = sys.argv[1]
customFile = sys.argv[2]

print (gbkFile)
print (customFile)

out = open(customFile, "w")
out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("seqnames","start","end","width","strand","id"))


for seq_record in SeqIO.parse(gbkFile, "genbank") :
    print "Extracting gene coordinates from GBK file: %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS":
            assert len(seq_feature.qualifiers['translation'])==1
            width=seq_feature.location.end-seq_feature.location.start

            out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    "NC_003028",
                    seq_feature.location.start,
                    seq_feature.location.end,
                    width,
                    seq_feature.strand,
                    seq_feature.qualifiers['locus_tag'][0]))

out.close()
print "Done"
