#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bibliothèque pour la manipulations de Séquence et d'ORF
Usage:
======
import SEQ
"""

__authors__ = ("Mathieu GENETE")
__contact__ = ("mathieu.genete@univ-lille.fr")
__copyright__ = "copyleft"
__date__ = "2022/01"
__version__= "0.5.0"

class DnaSeq:
    """
    Classe représentant une séquence d'ADN.

    Attributs :
    -----------
    id : str
        Identifiant de la séquence.
    description : str
        Description de la séquence.
    seq : str
        La séquence d'ADN.
    length : int
        Longueur de la séquence.
    orfNbr : int
        Nombre de cadres de lecture ouverts (ORF) trouvés.
    """
    
    def __init__(self,seqid,description,seq):
        """
        Initialise une nouvelle instance de la classe DnaSeq.

        Paramètres :
        ------------
        seqid : str
            Identifiant de la séquence.
        description : str
            Description de la séquence.
        seq : str
            La séquence d'ADN.
        """
        self.__id=seqid
        self.__description=description
        self.__seq=seq
        self.__length=len(seq)
        self.__orfNbr=0
        
    #===================
    #Getters Setters
    #===================
    
    @property
    def id(self):
        return self.__id
    
    @property
    def description(self):
        return self.__description
        
    @property
    def seq(self):
        return self.__seq
        
    @property
    def length(self):
        return self.__length
        
    @property
    def orfNbr(self):
        return self.__orfNbr
        
    #===================
    #Méthodes magiques
    #===================
    
    def __str__(self):
        """
        Retourne la séquence d'ADN sous forme de chaîne de caractères.

        Retour :
        --------
        str
            La séquence d'ADN.
        """
        return self.__seq
    
    #===================
    #Méthodes publiques
    #===================

    def computeORF(self,l,start_codons,stop_codons):
        """
        Calcule tous les cadres de lecture ouverts (ORF) dans la séquence et son complément inverse.

        Paramètres :
        ------------
        l : int
            Longueur minimale d'un ORF pour être considéré.
        start_codons : list
            Liste des codons de départ.
        stop_codons : list
            Liste des codons d'arrêt.

        Retour :
        --------
        list
            Liste des ORF trouvés.
        """
        out_ORF=[]
        out_ORF=self.__searchORF(self.__seq,start_codons,stop_codons,thrld=l)
        out_ORF+=self.__searchORF(self.reverse_complement(),start_codons,stop_codons,thrld=l,revcomp=True)
        
        return out_ORF
        
    def reverse_complement(self,sequence=None):
        """
        Calcule le complément inverse de la séquence d'ADN.

        Paramètres :
        ------------
        sequence : str, optionnel
            La séquence d'ADN à inverser et complémenter. Si None, utilise la séquence de l'instance.

        Retour :
        --------
        str
            Le complément inverse de la séquence.
        """
        if sequence is None:
            sequence = self.__seq
        comp={'A':'T','T':'A','C':'G','G':'C'}
        return "".join([comp[b.upper()] for b in sequence[::-1]])     
        
    #===================
    #Méthodes statiques
    #=================== 
    
    @staticmethod
    def parse_genetic_code(code_file):
        """
        Analyse un fichier de code génétique pour extraire les codons de départ et d'arrêt.

        Paramètres :
        ------------
        code_file : str
            Chemin vers le fichier de code génétique.

        Retour :
        --------
        tuple
            Un tuple contenant le code génétique, la liste des codons de départ et la liste des codons d'arrêt.
        """
        genetic_code={}
        start_codons=[]
        stop_codons=[]

        with open(code_file,'r') as handle:
            aa_line=handle.readline().strip().split('=')[1].strip()
            starts=handle.readline().strip().split('=')[1].strip()
            c1=handle.readline().strip().split('=')[1].strip()
            c2=handle.readline().strip().split('=')[1].strip()
            c3=handle.readline().strip().split('=')[1].strip()
        for i in range(len(aa_line)):
            aa=aa_line[i]
            cstart=starts[i]
            codon=c1[i]+c2[i]+c3[i]
            if aa not in genetic_code.keys():
                genetic_code[aa]={}
            if codon not in genetic_code[aa].keys():
                genetic_code[aa][codon]=0
            if cstart=='M':
                start_codons.append(codon)
            if cstart=='*':
                stop_codons.append(codon)
        return genetic_code,start_codons,stop_codons
        
    #===================
    #Méthodes privées
    #===================
    
    def __searchORF(self,sequence,start_codons,stop_codons,thrld=300,revcomp=False):
        """
        Recherche les cadres de lecture ouverts (ORF) dans la séquence.

        Paramètres :
        ------------
        sequence : str
            La séquence d'ADN à analyser.
        start_codons : list
            Liste des codons de départ.
        stop_codons : list
            Liste des codons d'arrêt.
        thrld : int, optionnel
            Longueur minimale d'un ORF pour être considéré (par défaut 300).
        revcomp : bool, optionnel
            Indique si la recherche doit être effectuée sur le complément inverse (par défaut False).

        Retour :
        --------
        list
            Liste des ORF trouvés.
        """
        hits = []
        
        for i in range(0,len(sequence)):
            cur_frame=i%3
            scodon = sequence[i:i+3]
            if scodon in start_codons:
                start_pos=i
                for j in range(i,len(sequence),3):
                    ecodon = sequence[j:j+3]
                    if ecodon in stop_codons:
                        stop_pos=j+3
                        self.__orfNbr+=1
                        if len(sequence[start_pos:stop_pos])>=thrld:
                            if revcomp:
                                Rend=len(sequence)-start_pos
                                Rstart=len(sequence)-stop_pos
                                rc_seq=self.reverse_complement(self.__seq[Rstart:Rend])
                                hits.append( OrfSeq("ORF_{}_Rev_{}_({},{})".format(len(hits),cur_frame,Rstart,Rend),"",rc_seq,self,Rstart,Rend,"-"))
                            else:
                                hits.append( OrfSeq("ORF_{}_Fwd_{}_({},{})".format(len(hits),cur_frame,start_pos,stop_pos),"",self.__seq[start_pos:stop_pos],self,start_pos,stop_pos,"+"))
                        break
        return hits
    
class OrfSeq(DnaSeq):
    """
    Classe représentant une séquence ORF (Open Reading Frame) dérivée d'une séquence d'ADN.

    Attributs :
    -----------
    refseq : DnaSeq
        Référence de la séquence d'ADN.
    start : int
        Position de départ de l'ORF.
    end : int
        Position de fin de l'ORF.
    strand : str
        Brin de la séquence ('+' ou '-').
    frame : int
        Cadre de lecture de l'ORF.
    score : float
        Score de l'ORF.
    score_rbs : float
        Score du site de liaison au ribosome (RBS).
    found_rbs : bool
        Indique si un site de liaison au ribosome a été trouvé.
    rbs_seq : str
        Séquence du site de liaison au ribosome.
    rbs_pos : int
        Position du site de liaison au ribosome.
    """
    
    def __init__(self,seqid,description,seq,refseq,start,end,strand):
        """
        Initialise une nouvelle instance de la classe OrfSeq.

        Paramètres :
        ------------
        seqid : str
            Identifiant de la séquence.
        description : str
            Description de la séquence.
        seq : str
            La séquence d'ADN.
        refseq : DnaSeq
            Référence de la séquence d'ADN.
        start : int
            Position de départ de l'ORF.
        end : int
            Position de fin de l'ORF.
        strand : str
            Brin de la séquence ('+' ou '-').
        """
        
        super().__init__(seqid,description,seq)
        self.__refseq=refseq
        self.__start=start
        self.__end=end
        self.__strand=strand
        if strand=="+":
            self.__frame=self.__start%3
        else:
            self.__frame=(refseq.length-self.__end)%3
        self.__score=0.0
        self.__score_rbs=0
        self.__found_rbs=False
        self.__rbs_seq=""
        self.__rbs_pos=0
        
    #===================
    #Getters Setters
    #===================
    
    @property
    def refseq(self):
        return self.__refseq
        
    @property
    def start(self):
        return self.__start
        
    @property
    def end(self):
        return self.__end
        
    @property
    def strand(self):
        return self.__strand
        
    @property
    def frame(self):
        return self.__frame
        
    @property
    def score(self):
        return self.__score
        
    @property
    def score_rbs(self):
        return self.__score_rbs
        
    @property
    def found_rbs(self):
        return self.__found_rbs
        
    @property
    def rbs_seq(self):
        return self.__rbs_seq
        
    @property
    def rbs_pos(self):
        return self.__rbs_pos
        
    @score.setter
    def score(self, value: float):
        self.__score = value
        
    #===================
    #Méthodes publiques
    #===================
    
    def recherche_RBS(self,matrice_poids,motif_range=(5,15),score_thrld=4.5):
        """
        Recherche un site de liaison au ribosome (RBS) dans la séquence.

        Paramètres :
        ------------
        matrice_poids : dict
            Matrice de poids pour le calcul du score RBS.
        motif_range : tuple, optionnel
            Intervalle de recherche du motif RBS (par défaut (5, 15)).
        score_thrld : float, optionnel
            Seuil de score pour considérer un RBS comme trouvé (par défaut 4.5).

        Retour :
        --------
        bool
            True si un RBS est trouvé, sinon False.
        """
        len_rbs_pattern=len(list(matrice_poids.values())[0])
        if self.__strand=='+':
            seq=self.__refseq.seq[self.__start-(motif_range[1]+len_rbs_pattern):self.__start-motif_range[0]]
        elif self.__strand=='-':
            seq=self.reverse_complement(self.__refseq.seq[self.__end+motif_range[0]:self.__end+(motif_range[1]+len_rbs_pattern)])
            
        tmpscore=0
        pos_pattern=0
        for i in range(0,len(seq)-len_rbs_pattern+1):
            sub_seq=seq[i:i+len_rbs_pattern]
            score=self.__calcul_score_RBS(matrice_poids,sub_seq)
                
            if score>=tmpscore:
                tmpscore=score
                pos_pattern=i+len_rbs_pattern-1
                self.__score_rbs=score
                self.__found_rbs=True
                self.__rbs_seq=sub_seq
                self.__rbs_pos=len(seq)-pos_pattern
                
        if self.__score_rbs>=score_thrld:
            return True

        return False
        
    def scoreORF(self,rscu_table):
        """
        Calcule le score de l'ORF basé sur la table RSCU.

        Paramètres :
        ------------
        rscu_table : dict
            Table RSCU (Relative Synonymous Codon Usage).

        Retour :
        --------
        float
            Score de l'ORF.
        """
        ri_values=[]
        for i in range(0,len(self.seq),3):
            codon=self.seq[i:i+3]
            if len(codon)==3:
                ri_values.append(rscu_table[codon]-1)
        return sum(ri_values)


    #===================
    #Méthodes privées
    #===================
    
    def __calcul_score_RBS(self,matrice_poids,seqRBS):
        """
        Calcule le score d'un site de liaison au ribosome (RBS).

        Paramètres :
        ------------
        matrice_poids : dict
            Matrice de poids pour le calcul du score RBS.
        seqRBS : str
            Séquence du site de liaison au ribosome.

        Retour :
        --------
        float
            Score du site de liaison au ribosome.
        """
        score=0
        for c,base in enumerate(seqRBS):
            score+=matrice_poids[base][c]
        return score
        
    #===================
    #Méthodes statiques
    #===================
    
    @staticmethod
    def parse_matrice_RBS(matricerbs_file):
        """
        Analyse un fichier de matrice de poids RBS.

        Paramètres :
        ------------
        matricerbs_file : str
            Chemin vers le fichier de matrice de poids RBS.

        Retour :
        --------
        dict
            Matrice de poids RBS.
        """
        rbs_mat={}
        with open(matricerbs_file,'r') as handle:
            bases=handle.readline().strip().split()
            rbs_mat={b:[] for b in bases}
            for line in handle:
                line=line.strip().split()
                for i,b in enumerate(bases):
                    rbs_mat[b].append(float(line[i]))
        return rbs_mat
    
    @staticmethod
    def compute_RSCU(count_table_file,codon_table):
        """
        Calcule la table RSCU (Relative Synonymous Codon Usage).

        Paramètres :
        ------------
        count_table_file : str
            Chemin vers le fichier de table de comptage des codons.
        codon_table : dict
            Table des codons.

        Retour :
        --------
        dict
            Table RSCU.
        """
        out_RSCU={}
        
        count_table=OrfSeq.parse_count_table(count_table_file)
        for aa,codons in codon_table.items():
            for c in codons:
                if c in count_table.keys():
                    codon_table[aa][c]=count_table[c]
        
        for aa,val in codon_table.items():
            sumCk=sum(val.values())
            dA=len(val.keys())
            for cod,Ci in val.items():
                out_RSCU[cod]=dA*(Ci/sumCk)
        return out_RSCU
        
    
    @staticmethod    
    def parse_count_table(count_table_file):
        """
        Analyse un fichier de table de comptage des codons.

        Paramètres :
        ------------
        count_table_file : str
            Chemin vers le fichier de table de comptage des codons.

        Retour :
        --------
        dict
            Table de comptage des codons.
        """
        out_count={}
        with open(count_table_file,'r') as handle:
            for line in handle:
                val=line.strip().split()
                out_count[val[0]]=int(val[1])
        return out_count
