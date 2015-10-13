#Margaret Antonio modified from

from Bio import SeqIO
import sys

gbk_filename = sys.argv[1]
faa_filename = sys.argv[2]

print (gbk_filename)
print (faa_filename)

output_handle = open(faa_filename, "w")
output_handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("seqnames","start","end","width","strand","id"))


for seq_record in SeqIO.parse(gbk_filename, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS":
            assert len(seq_feature.qualifiers['translation'])==1
            width=seq_feature.location.end-seq_feature.location.start

            output_handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    "NC_003028",
                    seq_feature.location.start,
                    seq_feature.location.end,
                    width,
                    seq_feature.strand,
                    seq_feature.qualifiers['locus_tag'][0]))

output_handle.close()
print "Done"
