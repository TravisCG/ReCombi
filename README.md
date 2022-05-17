# ReCombi
Idetifying recombination events and parental haplotypes

The codes found here can identify recombination events
and parental haplotypes using parental and sibling genotypes.
The whole information can be packed toghether into a VCF file.

If you do not have any parent-sibling VCF, you can generate one
with my simulation toolkit: 
[https://github.com/TravisCG/BUG/tree/main/simvcf](https://github.com/TravisCG/BUG/tree/main/simvcf)

## recombipos
The script reads the VCF file and identify recombination events in
siblings.

Usage:

python3 recodet.py input.vcf.gz motherid fatherid windowsize

input.vcf.gz: Gzipped input VCF file with parents and siblings only
motherid: ID of the mother in the VCF file
fatherid: ID of the father in the VCF file
windowsize: default value is 200. Use this window for the statistical calculation
