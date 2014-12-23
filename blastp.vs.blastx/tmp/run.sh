
echo "start"
date
blastall -p blastp -d sce.protein.faa -i query.protein.faa -F F -o _out.prot.txt
echo "blastp done"
date
blastall -p blastx -d sce.protein.faa -i query.cds.faa -F F -o _out.dna.txt
echo "blastx done"
date


