from Bio import SeqIO

gbk_filename = "NC_003028.gbk"
faa_filename = "tigr4_convereted.faa"

output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(gbk_filename, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS":
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s from %s\n%s\n" % (
                                                       seq_feature.qualifiers['locus_tag'][0],
                                                       seq_record.name,
                                                       seq_feature.qualifiers['translation'][0]))
            print(int(seq_feature.location.start),int(seq_feature.location.end),seq_feature.strand)

output_handle.close()
print "Done"
