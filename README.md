# Simple_ORF_finder.py

## Description
Programme de recherche d'ORF basé sur l'usage du biais des codons.

## Utilisation
Ce programme permet de rechercher des ORFs (Open Reading Frames) dans un fichier fasta de génome en utilisant une table de comptage des codons issue de régions codantes connues et un code génétique pour la traduction.

## Arguments
- `-i`, `--infasta` : Fichier fasta du génome à analyser (obligatoire).
- `-c`, `--counttable` : Table de comptage des codons issue de régions codantes connues (obligatoire).
- `-t`, `--codonstable` : Code génétique utilisé pour la traduction (obligatoire).
- `-o`, `--outgff` : Nom du fichier de sortie au format GFF (obligatoire).
- `-f`, `--outfasta` : Séquences des ORFs au format fasta (obligatoire).
- `-r`, `--matricerbs` : Matrice RBS (optionnel).
- `-s`, `--minorfsize` : Taille minimum des ORFs (valeur par défaut = 300) (optionnel).
- `-u`, `--orfthrld` : Seuil du score des ORF à conserver (par défaut 0) (optionnel).
- `-b`, `--rbsmotifthrld` : Score seuil pour la recherche du motif RBS (valeur par défaut = 4.5) (optionnel).
- `-e`, `--rbsregion` : Région de recherche du RBS en amont du codon start (2 valeurs séparées par un espace: min max) par défaut [5,15] (optionnel).
- `-d`, `--displayorf` : Affiche les ORF dans un interval donné (2 valeurs séparées par un espace: début fin) (optionnel).

## Exemple de commande
```bash
python Simple_ORF_finder.py -i genome.fasta -c count_table.txt -t codon_table.txt -o output.gff -f orfs.fasta
```

## Ligne de commande avec les données de test:
```
./Simple_ORF_finder.py -i datas/GCF_000009045.1_ASM904v1_partial_genomic.fna -c datas/E_coli_K12_MG1655_CDS_table -f found_ORF_SimpleORF.fasta -t datas/code_genetique_bacteries -r datas/matrice_GGAGA_RBS_E_coli -o test2.gff -d 550 12000
```
