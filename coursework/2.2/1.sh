wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.MT.fa mt.fa
wget https://ftp.ensembl.org/pub/current_fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.primary_assembly.4.fa.gz
gzip -d Drosophila_melanogaster.BDGP6.54.dna.primary_assembly.4.fa.gz
mv Drosophila_melanogaster.BDGP6.54.dna.primary_assembly.4.fa Dr.fa
wget https://ftp.ensembl.org/pub/current_fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.chromosome.I.fa.gz
gzip -d Caenorhabditis_elegans.WBcel235.dna.chromosome.I.fa.gz
mv  Caenorhabditis_elegans.WBcel235.dna.chromosome.I.fa elegans.fa
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo.fa



