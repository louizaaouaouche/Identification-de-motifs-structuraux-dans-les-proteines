# -*- coding: utf-8 -*
##Fonctions de Lecture de fichiers PDB

import math
import time
from dim2 import *
from dim1 import *

import re
import sys

class Atom:
	def __init__(self, x, y, z, numRes=-1, resType="X", sse="X"):
		self.x = x #
		self.y = y # positions de l'atome dans l'espace 
		self.z = z #
		self.numRes=numRes 	 #numéro PDB du résidu
		self.resType=resType #code 3 lettres de l'aa
		self.sse=sse 		 #Code 1 lettre de la Structure secondaire

def getCAAtom(ligne):
	"""str->str*Atom ou None
	ATOM   1272  CA  MET A 145      16.640  22.247  -9.383  1.00 56.65           C  """

	re_whitespaces = re.compile("[ ]+")
	resname = re_whitespaces.sub('',ligne[17:20])   #nom de la protéine ex: GLY, SER, HIS ...
	atomname = re_whitespaces.sub('',ligne[12:16])  #nom d'atome:Carbone, Oxygen ....
	chain = ligne[21]

	if chain == " ":
		chain ="-"
	resnumber = int(ligne[8:12])   # numéro de la protéine auquel appartient de l'atome 
	x = float(ligne[30:38])			#position x
	y = float(ligne[38:46])			#position y
	z = float(ligne[46:54])			#position z
	
	if atomname=="CA": #On ne retourne que les carbones alpha
		return chain, Atom(x,y,z,resnumber, resname)
	else:
		return chain, None

def getISeqNumHelix(ligne):
	"""str->[int]
	HELIX 1 HAGLYA 86 GLYA 94 1 9"""

	lesPosH=[]
	re_whitespaces = re.compile("[ ]+")
	
	initChainID = re_whitespaces.sub('',ligne[19])
	initSeqNum = int(ligne[21:25])
	endChainID = re_whitespaces.sub('',ligne[31])
	endSeqNum = int(ligne[33:37])
	if initChainID != endChainID:
 		print ("Problème deux chaines différentes ? ", initChainID, endChainID)
 		print (ligne)
	else:
		for i in range(initSeqNum, endSeqNum+1):
			lesPosH.append(i)
	return lesPosH


def getiSeqNumSheet(ligne):
	"""str->[int]"""

	lesPosS=[]
	re_whitespaces = re.compile("[ ]+")
	
	initChainID = re_whitespaces.sub('',ligne[21])
	initSeqNum = int(ligne[22:26])
	
	endResName = re_whitespaces.sub('',ligne[28:31])
	endChainID = re_whitespaces.sub('',ligne[32])
	endSeqNum = int(ligne[33:37])
	
	if initChainID != endChainID:
 		print ("Problème deux chaines différentes ? ", initChainID, endChainID)
 		print (ligne)
	else:
		for i in range(initSeqNum, endSeqNum+1):
			lesPosS.append(i)
	return lesPosS


def ajouterSSEAtomes(lesAtomes, lesH, lesS):
	"""[Atom]*[int]*[int]->[Atom]"""

	if(len(lesH)==0 and len(lesS)==0):
		return lesAtomes
	lesH.sort()
	lesS.sort()
	iAt=0
	iH=0
	iS=0
	#3)lesAtomes
	while iAt<len(lesAtomes) and (iH<len(lesH) or iS<len(lesS)):
		if iH<len(lesH) and lesAtomes[iAt].numRes==lesH[iH]:
			lesAtomes[iAt].sse="H"
			iH+=1
		elif iS<len(lesS) and lesAtomes[iAt].numRes==lesS[iS]:
			lesAtomes[iAt].sse="S"
			iS+=1
		iAt+=1
	return lesAtomes

def getNumModel(ligne):
	return int(ligne[10:14])

def readPDB(nomFi):
	"""str->[Atom]"""

	try:
		f = open(nomFi, "r")
	except IOError:
		print ("readPDB:: Fichier <%s> introuvable, arret du programme"%(nomFi))
		sys.exit(1) 

	re_resolution = re.compile("^REMARK {3}2 RESOLUTION.*ANGSTROMS")
	re_atom = re.compile("^ATOM")
	re_helix = re.compile("^HELIX")
	re_sheet = re.compile("^SHEET")
	re_model = re.compile("^MODEL")
	numModel=None
	lesResSheet=[]
	lesResHelix=[]
	lesChaines=[]
	lesAtomes=[]

	for ligne in f:
		if re_model.search(ligne):
			numModel=getNumModel(ligne)
		elif re_atom.search(ligne):
			c,a=getCAAtom(ligne)
			if (numModel==None):
			    numModel=0
			if (a!=None and numModel==0) :
				lesAtomes.append(a)
		elif re_helix.search(ligne):
			lesResHelix+=getISeqNumHelix(ligne)
		elif re_sheet.search(ligne):
			lesResSheet+=getiSeqNumSheet(ligne)

	lesAtomes=ajouterSSEAtomes(lesAtomes, lesResHelix, lesResSheet)

	return lesAtomes

## Calcul de matrices discretise.

#[Atom]*[Atom]->int
def calcul(at1,at2):
    return math.sqrt((at1.x-at2.x)**2+(at1.y-at2.y)**2+(at1.z - at2.z)**2)

#[Atom]->[Atom]
def matrice(lesCa):
    n=len(lesCa)
    mat=[[]]*n
    for num in range(0,n):
        mat[num]=[0]*n
        mat[num][num]=-1
    for a1 in range(0,n):
        for a2 in range(0,n):
            if a1 != a2:
                mat[a1][a2]=calcul(lesCa[a1],lesCa[a2])
            else :
                mat[a1][a2]=0

    return mat

def disct(mat,pas=5):
    n=len(mat)
    for a1 in range(0,n):
        for a2 in range(0,n):
            if mat[a1][a2]!=-1:
                mat[a1][a2]=(int)(mat[a1][a2]/pas)

    
#renvoie la borne inférieur de l'intervalle de la taille du pas au quel le nombre appartient
def disc_nb( nb , pas ) :
    ''' Number * Number -> Number 
    la taille de l'intervalle peut etre entière ou decimal du type (10^(-n)*diviseur_entier_de_10)'''
    res = 0
    if ( pas < 0 ):
        # on pourrait générer une exception
        print( "erreur le pas doit être positif")
    if (pas == 0):
        return nb
    if ( pas >= 1 ):
        return (nb//pas)*pas
    else :
        res = int(nb)
        nb = (nb % 1)
        cpt = 1
        pas = pas*10
        q = int(pas)
        while ( q == 0) :
            res += (nb//(10**cpt))/(10**cpt)
            nb = nb%(10**cpt)
            cpt += 1
            pas = pas *10
            q = pas // 1
                
        return res + (int(nb*(10**cpt)/pas)*pas)/(10**cpt)

def disc_mat(mat ,pas=1):
    n=len(mat)
    #on parcours la matrice de distance 
    for a1 in range(0,n):
        for a2 in range(0,n):
            mat[a1][a2]= disc_nb(mat[a1][a2],pas)

def calculmat(lesCa):
    tmp=matrice(lesCa)
    disc(tmp,1)
    return tmp

def affiMat(mat):
    for num in range(0,len(mat)):
        print (mat[num])

def readAA(filename):
        ''' str -> List[str] * intcette fonction nous sert à lire les fichiers .fasta '''
        f = open(filename) 
        lenADN=0
        lis=[]
        lines = f.readlines()
        for x in range(0,len(lines)):
                if (x%2==1):
                        stri=lines[x]
                        st=stri[:-1]
                        lenADN+=len(st)
                        lis.append(st)
        return lis,lenADN


def main() :
	if (len(sys.argv)<3):
		print ("USAGE: ", sys.argv[0], "<pdb file>","<dim>","<aff>")
		sys.exit(1)
	
	dim = (int)(sys.argv[2])
	filename = sys.argv[1]
	aff=(int)(sys.argv[3])

	if (dim==2):
		
		lesAtomes=readPDB(filename)
		print("la longueur de la liste d'atomes :" , len(lesAtomes))
                                
		mt=matrice(lesAtomes)
		disc_mat(mt)

		t0 = time.time()

		lp,lenMotif=comparerSequences_DIM2(mt)		

		print (time.time() - t0, "seconds")		
		
		if (aff>0):
			affichage_2d(lp,lesAtomes,lenMotif)
			print("lenMotif=",lenMotif)
			print ("Nombre des Atomes:",len(lesAtomes))
		

	elif (dim==1):

		al=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','*']
		lis,ladn=readAA(filename)

		t0 = time.time()
		
		comparerSequences(lis,al,sys.maxsize,1,aff)

		print ("\n",time.time() - t0, "seconds")
		print ("Nombre des Acides amines:",ladn)



if __name__ == '__main__':
	main()


