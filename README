# Projet Simple_ORF_finder

## Options de Simple_ORF_finder
usage: Simple_ORF_finder.py [-h] -i INFASTA -c COUNTTABLE -t CODONSTABLE -o OUTGFF -f OUTFASTA [-r MATRICERBS] [-s MINORFSIZE]
                            [-u ORFTHRLD] [-b RBSMOTIFTHRLD] [-e RBSREGION RBSREGION] [-d DISPLAYORF DISPLAYORF]

Programme de recherche d'ORF basé sur l'usage du biais des codon

options:
  -h, --help            show this help message and exit
  -i INFASTA, --infasta INFASTA
                        fichier fasta du génome à analyser
  -c COUNTTABLE, --counttable COUNTTABLE
                        table de comptage des codons issue de régions codantes connues
  -t CODONSTABLE, --codonstable CODONSTABLE
                        code génétique utilisé pour la traduction
  -o OUTGFF, --outgff OUTGFF
                        nom du fichier de sortie au format GFF
  -f OUTFASTA, --outfasta OUTFASTA
                        séquences des ORFs au format fasta
  -r MATRICERBS, --matricerbs MATRICERBS
                        matrice RBS
  -s MINORFSIZE, --minorfsize MINORFSIZE
                        taille minimum des ORFs (valeur par défaut = 300)
  -u ORFTHRLD, --orfthrld ORFTHRLD
                        Seuil du score des ORF à conserver (par défaut 0)
  -b RBSMOTIFTHRLD, --rbsmotifthrld RBSMOTIFTHRLD
                        Score seuil pour la recherche du motif RBS (valeur par défaut = 4.5)
  -e RBSREGION RBSREGION, --rbsregion RBSREGION RBSREGION
                        Région de recherche du RBS en amont du codon start (2 valeurs séparées par un espace: min max) par
                        défaut [5,15]
  -d DISPLAYORF DISPLAYORF, --displayorf DISPLAYORF DISPLAYORF
                        Affiche les ORF dans un interval donné (2 valeurs séparées par un espace: début fin)


## Ligne de commande de test:
```
./Simple_ORF_finder.py -i datas/GCF_000009045.1_ASM904v1_partial_genomic.fna -c datas/E_coli_K12_MG1655_CDS_table -f found_ORF_SimpleORF.fasta -t datas/code_genetique_bacteries -r datas/matrice_GGAGA_RBS_E_coli -o test2.gff -d 550 12000
```
