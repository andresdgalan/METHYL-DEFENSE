# In this script we will create a mapping txt to know the position in the genome of the start and end of each contig
exec > >(tee -i DENOVO8_genome_position_map.log) 2>&1

awk '
BEGIN {seq_len=0; contig_start=1}
/^>/ {
    if(NR>1){
        contig_end = seq_len
        print contig_id "\t" contig_start "\t" contig_end
        contig_start = seq_len + 1
    }
    contig_id = substr($0,2)
    next
}
{
    seq_len += length($0)
}
END {
    contig_end = seq_len
    print contig_id "\t" contig_start "\t" contig_end
}' pseudogenome_contigs.fa > Plantago_genome_position_map.txt
