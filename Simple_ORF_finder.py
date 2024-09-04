#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Programme de recherche d'ORF et identification de gènes pour génome procaryote
Usage:
======
python Simple_ORF_finder.py -h
pour afficher l'aide

"""
__authors__ = ("Mathieu GENETE")
__contact__ = ("mathieu.genete@univ-lille.fr")
__copyright__ = "copyleft"
__date__ = "2022/01"
__version__= "1.0.0"

import SEQ
import sys
import os
import argparse
import numpy as np
from matplotlib import pyplot as plt


def main():
    """
    Boucle principale du programme
    """
    #sys.setrecursionlimit(1500)
    
    #Définition des paramètres du programme avec argparse
    description=""" Programme de recherche d'ORF basé sur l'usage du biais des codon """
    parser=argparse.ArgumentParser(prog="Simple_ORF_finder.py",description=description)
    parser.add_argument("-i","--infasta", help="fichier fasta du génome à analyser", required = True)
    parser.add_argument("-c","--counttable", help="table de comptage des codons issue de régions codantes connues", required = True)
    parser.add_argument("-t","--codonstable", help="code génétique utilisé pour la traduction", required = True)
    parser.add_argument("-o","--outgff", help="nom du fichier de sortie au format GFF", required = True)
    parser.add_argument("-f","--outfasta", help="séquences des ORFs au format fasta", required = True)
    parser.add_argument("-r","--matricerbs", help="matrice RBS")
    parser.add_argument("-s","--minorfsize", help="taille minimum des ORFs (valeur par défaut = 300)",type=int,default=300)
    parser.add_argument("-u","--orfthrld", help="Seuil du score des ORF à conserver (par défaut 0)",type=float,default=0.0)
    parser.add_argument("-b","--rbsmotifthrld", help="Score seuil pour la recherche du motif RBS (valeur par défaut = 4.5)",type=float,default=4.5)
    parser.add_argument("-e","--rbsregion", help="Région de recherche du RBS en amont du codon start (2 valeurs séparées par un espace: min max) par défaut [5,15]",nargs=2,type=int,default=[5,15])
    parser.add_argument("-d","--displayorf", help="Affiche les ORF dans un interval donné (2 valeurs séparées par un espace: début fin)",nargs=2,type=int)
    args = parser.parse_args(sys.argv[1:])
    
    infasta=args.infasta
    counttable=args.counttable
    outgff=args.outgff
    minorfsize=args.minorfsize
    outfasta=args.outfasta
    codonstable=args.codonstable
    rbsmotifthrld=args.rbsmotifthrld

    
    print("== Simple_ORF_finder v{ver} ==\nAutheur: {auth}\nmail: {mail}\n".format(ver=__version__,auth=__authors__,mail=__contact__))
    print("Commande: {} \n\n".format(" ".join(sys.argv[0:])))
    
    print("=== Parse le fichier fasta === ")
    #Parse le ficher fasta d'entrée et retourne une liste d'objets DnaSeq
    dna_seqs=parseFasta(infasta)
    print("{} séquence(s) récupérée(s)".format(len(dna_seqs)))
    for s in dna_seqs:
        print("\t{} : {} bp".format(s.id,s.length))
    
    print("\n=== Parse la table du code génétique ===")
    #Parse la table des codons, retourne 3 variables: dictionnaire tables des codons
    #liste des codons start, liste des codons stop
    codon_table,start_codons,stop_codons=SEQ.OrfSeq.parse_genetic_code(codonstable)
    Rscu_table=SEQ.OrfSeq.compute_RSCU(counttable,codon_table)
    for c in sorted(Rscu_table.keys()):
        print(c,Rscu_table[c])
    print("table utilisée: {}".format(os.path.basename(codonstable)))
    print("Codons start: {}\nCodons stop: {}".format(" , ".join(start_codons)," , ".join(stop_codons)))
    print("Liste des codons synonymes:")
    for aa,codons in sorted(codon_table.items()):
        print("{} : {}".format(aa," , ".join(codons)))
    
    print("\n=== Recherche des ORFs ===")
    #Recherche les ORF pour chaque objet DnaSeq
    #Les ajoute à une liste d'objets OrfSeq: found_ORF
    found_ORF=[]
    total_ORF=0
    for rec in dna_seqs:
        found_ORF+=rec.computeORF(minorfsize,start_codons,stop_codons)
        total_ORF+=rec.orfNbr
    
    print("Nombre total d'ORF: %i" % (total_ORF))
    print("\nNombre d'ORFs de taille >=%i : %i" %(minorfsize,len(found_ORF)))
    orf_len=[o.length for o in found_ORF]
    mean_ORF_len=np.mean(orf_len)
    min_ORF_len=min(orf_len)
    max_ORF_len=max(orf_len)
    print("Taille moyenne des ORFs : {} bp - [{},{}] ".format(round(mean_ORF_len,1),min_ORF_len,max_ORF_len))
    
    #Tri des ORF par position start croissante puis end..
    #Voir aussi fonction sort_ORF() tri à bulle
    #sorted_found_ORF=sorted(found_ORF, key= lambda x:x.start)
    sorted_found_ORF=sorted(found_ORF, key= lambda x:(x.start,x.end))
    
    if args.displayorf:
        debutI=args.displayorf[0]
        finI=args.displayorf[1]
        if debutI<finI:
            print("\n=== Affiche les ORF dans l'interval [{},{}] ===".format(args.displayorf[0],args.displayorf[1]))
            #Affiche les ORF dans un interval défini
            displayORF(sorted_found_ORF,debutI,finI)
    
    #Si une matrice pour les sites de fixation du ribosome est entrée
    if args.matricerbs:
        print("\n=== Filtrer les RBS avec la matrice: {} ===".format(os.path.basename(args.matricerbs)))
        print("paramètres utilisés:\n\tScore seuil pour la recherche du motif RBS: {}\n\tRégion de recherche du RBS en amont du codon start: {}".format(args.rbsmotifthrld,args.rbsregion))
        
        #Parse le fichier text de la matrice de pois
        matrice_RBS=SEQ.OrfSeq.parse_matrice_RBS(args.matricerbs)
        
        #Recherche les RBS dans la liste d'ORF précédemment triée
        #Si un motif est trouvé au bon endroit avec le bon seuil
        #ajout l'ORF à la liste sorted_found_ORF_RBS
        sorted_found_ORF_RBS=[]
        for orf in sorted_found_ORF:
            if orf.recherche_RBS(matrice_RBS,score_thrld=rbsmotifthrld,motif_range=args.rbsregion):
                sorted_found_ORF_RBS.append(orf)
                
        #On écrase la variable sorted_found_ORF avec les ORF avec un RBS
        #en amont
        sorted_found_ORF=sorted_found_ORF_RBS

        print("Nombre d'ORFs après filtrage: {}".format(len(sorted_found_ORF)))

    #Calcul le score pour chaque ORF. Les ORF avec un score <0 (valeur par défaut)
    #sont éliminés
    print("\n=== Compute ORF Score ===")
    putativ_ORF=computeORFscore(sorted_found_ORF,Rscu_table,score_thrld=args.orfthrld)
    print("Nombre d'ORFs avec un score >={}: {}".format(args.orfthrld,len(putativ_ORF)))
    write_GFF("putativ_ORF_score_up_zero.gff",putativ_ORF)

    #Elimine les ORF chevauchants et maximise le score des ORFs
    print("\n=== Maximise ORF ===")
    #maximiseORFscore()
    #Tri les ORFs par position de fin
    sorted_putativ_ORF=sorted(putativ_ORF, key= lambda x:x.end)

    #Récupère les scores des ORFs dans un vecteur
    ORF_scores=[v.score for v in sorted_putativ_ORF]
    
    #Génère le vacteur P des prédécesseurs
    p=predecessor(sorted_putativ_ORF)
    
    #Initialiser la matrice M pour récupérer la table des scores
    M=[-1]*(len(ORF_scores))
    
    #Calcule le score max pour l'ensemble des ORFs
    # et récupère le score max
    max_score=maximiseORFscore(len(ORF_scores)-1,ORF_scores,p,M)
    
    #Récupère les ORFs qui maximisent le score à  l'aide de la table M
    #et des prédécesseurs
    max_orf=[sorted_putativ_ORF[i] for i in return_ORFs_max(M,p)]
    
    print("Nombre d'ORFs non chevauchants après maximisation du score: {}".format(len(max_orf)))
    print("Score maximisé pour les {} ORFs: {}".format(len(max_orf),max_score))
    #Ecrit les séquences des ORF dans un fichier fasta
    writeORF_in_fasta(max_orf,outfasta)
    
    #Plot les ORFs sur un graphique pour situer les positions
    plot_ORFs(max_orf,"max_ORF_map.png")
    
    #Ecrit les coordonnées et annotation des ORF dans un fichier GFF
    write_GFF(outgff,max_orf)



#===================
#     Fonctions
#===================
    
def plot_ORFs(orf_list,outpng):
    """
    Trace les ORFs (Open Reading Frames) sur un graphique et enregistre l'image.

    Parameters
    ----------
    orf_list : list
        Liste d'objets ORF contenant les informations des ORFs.
    outpng : str
        Chemin du fichier de sortie pour enregistrer l'image.

    Returns
    -------
    None
    """
    f=5
    ymax=len(orf_list)*f
    fig, ax = plt.subplots()
    y=ymax
    n=1
    stcolors={"+":"green","-":"red"}
    for o in orf_list:
        ax.hlines(y,o.start,o.end,color=stcolors[o.strand])
        ax.text(o.start,y+1,"ORF_{}".format(n),fontsize=2)
        n+=1
        y+=f
    ax.get_yaxis().set_visible(False)
    fig.set_size_inches(10,10)
    fig.savefig(outpng, dpi=300)
    
def plot_RSCU(Rscu_table,DSeq,window=300):
    """
    Trace les valeurs RSCU (Relative Synonymous Codon Usage) sur des fenêtres glissantes de la séquence.

    Parameters
    ----------
    Rscu_table : dict
        Table des valeurs RSCU pour chaque codon.
    DSeq : SeqRecord
        Objet SeqRecord contenant la séquence d'ADN.
    window : int, optional
        Taille de la fenêtre glissante (par défaut 300).

    Returns
    -------
    None
    """
    seq=DSeq.seq
    RSCU_datas={-3:[],-2:[],-1:[],1:[],2:[],3:[]}
    xvals=[]
    for i in range(0,len(seq)-window-1):
        xvals.append(i)
        win_seq=seq[i:i+window]
        rev_win_seq=reverse_complement(win_seq)
        for c in range(0,3):
            win_codons=[win_seq[j:j+3] for j in range(c,len(win_seq)-c,3)]
            rev_win_codons=[rev_win_seq[j:j+3] for j in range(c,len(rev_win_seq)-c,3)]
            RSCU_datas[c+1].append(sum([Rscu_table[v]-1 for v in win_codons if len(v)==3]))
            RSCU_datas[-c-1].append(sum([Rscu_table[v]-1 for v in rev_win_codons if len(v)==3]))
    
    fig, axs = plt.subplots(6)
    fig.set_size_inches(100,20)
    ax_fig=0
    for x,y_val in RSCU_datas.items():
        axs[ax_fig].plot(xvals,y_val)
        ax_fig+=1
    fig.savefig('scatter.png', dpi=150)

def reverse_complement(sequence):
    """
    Retourne le complément inverse d'une séquence d'ADN.

    Parameters
    ----------
    sequence : str
        Séquence d'ADN.

    Returns
    -------
    str
        Complément inverse de la séquence.
    """
    comp={'A':'T','T':'A','C':'G','G':'C'}
    return "".join([comp[b.upper()] for b in sequence[::-1]])
    
def sort_ORF(orflist,orf_attr,reverse=False):
    """
    Réalise un tri à bulle sur la position start des ORF    
    Parameters
    ----------
    orflist : list
    liste dobjets OrfSeq
    
    orf_attr : string
    Attribut de l'objet OrfSeq sur lequel faire le tri
    
    reverse : boolean
    Tris croissant par défaut, si reverse=True, tri décroissant
    
    Returns
    -------
    list
    Retourne la liste des ORFs triés par la position start
    """
    try:
        l=len(orflist)
        while l>1:
            for i in range(0,l-1):
                if reverse:
                    test_val=getattr(orflist[i],orf_attr)<getattr(orflist[i+1],orf_attr)
                else:
                    test_val=getattr(orflist[i],orf_attr)>getattr(orflist[i+1],orf_attr)
                    
                if test_val:
                    tmp=orflist[i+1]
                    orflist[i+1]=orflist[i]
                    orflist[i]=tmp
            l-=1
    except AttributeError as err:
        print("Erreur fonction sort_ORF:",err)
        sys.exit()
    
    return orflist

def write_GFF(outgff,orflist):
    """
    Écrit les informations des ORFs dans un fichier GFF.

    Parameters
    ----------
    outgff : str
        Chemin du fichier de sortie GFF.
    orflist : list
        Liste d'objets OrfSeq.

    Returns
    -------
    None
    """
    sorted_orf=sorted(orflist,key=lambda x:x.start)
    with open(outgff,'w') as handle:
        handle.write("##gff-version 3\n")
        for orf in sorted_orf:
            attributes_dict={"ORF_id":orf.id,"Score_RBS":orf.score_rbs,"rbs_pos":orf.rbs_pos,"rbd_seq":orf.rbs_seq}
            attr=" ; ".join(['{} = {}'.format(k,v) for k,v in attributes_dict.items()])
            gff_line="{chr}\tSimple_ORF_finder\tgene\t{start}\t{end}\t{score}\t{brin}\t{phase}\t{attr}\n".format(chr=orf.refseq.id,start=orf.start+1,end=orf.end,score=orf.score,brin=orf.strand,phase=orf.frame,attr=attr)
            handle.write(gff_line)
            
def affiche_ORFs(orf_list):
    """
    Affiche les informations des ORFs.

    Parameters
    ----------
    orf_list : list
        Liste d'objets OrfSeq.

    Returns
    -------
    None
    """
    for orf in orf_list:
        print(orf.id,orf.strand,orf.start,orf.end,(orf.end-orf.start),orf.score,orf.score_rbs)
        
def computeORFscore(orf_list,Rscu_table,score_thrld=None):
    """
    Calcule le score des ORFs et filtre selon un seuil.

    Parameters
    ----------
    orf_list : list
        Liste d'objets OrfSeq.
    Rscu_table : dict
        Table des valeurs RSCU pour chaque codon.
    score_thrld : float, optional
        Seuil de score pour filtrer les ORFs (par défaut None).

    Returns
    -------
    list
        Liste des ORFs avec un score supérieur ou égal au seuil.
    """
    out_orf=[]
    for orf in orf_list:
        score=orf.scoreORF(Rscu_table)
        if score_thrld is None or score>=score_thrld:
            orf.score=score
            out_orf.append(orf)
    return out_orf

def recherher_dicot(sorted_ORF,DebutIntervalle):
    """
    Recherche dichotomique de l'intervalle de début dans une liste triée d'ORFs.

    Parameters
    ----------
    sorted_ORF : list
        Liste triée d'objets ORF.
    DebutIntervalle : int
        Position de début de l'intervalle à rechercher.

    Returns
    -------
    int
        Index de l'élément trouvé ou position d'insertion.
    """
    sorted_ORF_list=[o.start for o in sorted_ORF]
    a=0
    b=len(sorted_ORF_list)-1
    m=(a+b)//2
    while a<b:
        if sorted_ORF_list[m] == DebutIntervalle:
            return m
        elif sorted_ORF_list[m] > DebutIntervalle:
            b = m-1
        else:
            a = m+1
        m=(a+b)//2
    return a

def displayORF(sorted_ORF,DebutIntervalle,FinIntervalle):
    """
    Affiche les ORFs dans un intervalle donné.

    Parameters
    ----------
    sorted_ORF : list
        Liste triée d'objets ORF.
    DebutIntervalle : int
        Position de début de l'intervalle.
    FinIntervalle : int
        Position de fin de l'intervalle.

    Returns
    -------
    list
        Liste des ORFs affichés.
    """
    found_Debut=recherher_dicot(sorted_ORF,DebutIntervalle)
    out_orf=[]
    for p in range(found_Debut,len(sorted_ORF)):
        orf=sorted_ORF[p]
        if orf.start<=FinIntervalle:
            print("{chr}\tSimple_ORF_finder\tgene\t{start}\t{end}\t{score}\t{brin}\t{phase}".format(chr=orf.refseq.id,start=orf.start,end=orf.end,score=orf.score,brin=orf.strand,phase=orf.frame))
    return out_orf
        

def writeORF_in_fasta(ORF_list,outfasta):
    """
    Écrit les séquences des ORFs dans un fichier FASTA.

    Parameters
    ----------
    ORF_list : list
        Liste d'objets ORF.
    outfasta : str
        Chemin du fichier de sortie FASTA.

    Returns
    -------
    None
    """
    with open(outfasta,'w') as fhandle:
        for orf in ORF_list:
            fhandle.write(">{} {}_{}\n{}\n".format(orf.id,orf.length,orf.refseq.id,orf.seq))

def parseFasta(Fi):
    """
    Analyse un fichier FASTA et retourne une liste d'objets DnaSeq.

    Parameters
    ----------
    Fi : str
        Chemin du fichier FASTA à analyser.

    Returns
    -------
    list
        Liste d'objets DnaSeq.
    """
    out_DnaSeq=[]
    with open(Fi,'r') as infasta:
        seqid=""
        description=""
        seq=""
        for line in infasta:
            #vérifie si la ligne commence par >
            if line.startswith(">"):
                if seqid!="":
                    out_DnaSeq.append(SEQ.DnaSeq(seqid,description,seq))
                tmp=line[1:].strip().split()
                seqid=tmp[0]
                description=" ".join(tmp[1:])
                seq=""
            else:
                seq+=line.strip()
        out_DnaSeq.append(SEQ.DnaSeq(seqid,description,seq))
    return out_DnaSeq

def return_ORFs_max(M,p):
    """
    Retourne les indices des ORFs maximisant le score.

    Parameters
    ----------
    M : list
        Liste des scores maximisés.
    p : list
        Liste des prédécesseurs.

    Returns
    -------
    list
        Liste des indices des ORFs maximisant le score.
    """
    index=M.index(max(M))
    ORF_indexs=[]
    while index>0:
        if M[index]!=M[index-1]:
            ORF_indexs.append(index)
            index=p[index]
        else:
            index-=1
    return ORF_indexs
    
def predecessor(sorted_orf):
    """
    Calcule les prédécesseurs pour chaque ORF.

    Parameters
    ----------
    sorted_orf : list
        Liste triée d'objets ORF.

    Returns
    -------
    list
        Liste des prédécesseurs pour chaque ORF.
    """
    p=[0]*(len(sorted_orf))
    
    for i in range(0,len(sorted_orf)):
        for j in range(i,-1,-1):
            if sorted_orf[i].start>=sorted_orf[j].end:
                p[i]=j
                break
    return p


def maximiseORFscore(j,ORF_scores,p,M):
    """
    Calcule le score maximal des ORFs en utilisant la programmation dynamique.

    Parameters
    ----------
    j : int
        Index de l'ORF actuel.
    ORF_scores : list
        Liste des scores des ORFs.
    p : list
        Liste des prédécesseurs.
    M : list
        Liste des scores maximisés.

    Returns
    -------
    int
        Score maximal des ORFs.
    """
    if j==0:
        M[0]=0
        return 0
    if M[j] != -1:
        return M[j]
    T1=ORF_scores[j] + maximiseORFscore(p[j],ORF_scores,p,M)
    T2=maximiseORFscore(j-1,ORF_scores,p,M)
    M[j]=max(T1,T2)
    return M[j]
        
if __name__=="__main__":
    main()
