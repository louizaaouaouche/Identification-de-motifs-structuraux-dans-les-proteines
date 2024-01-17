# -*- coding: utf-8 -*
import sys
## Dimension 2

# fonction de création du 1er vecteur V à partir de la matrice de distance discrétisée ( dimension 2)

def creer_vv_de_base(Mdist , dimMdist):
    '''List[List[Number]]*int -> List[List[int]] 
      crée le premier vecteur V à partir de la matirice de distance discrétisée'''

    # vv : List[List[int]]
    vv=[]
    ens = set()

    for i in range (0,dimMdist):
        vv.append([])
        for j in range (0,dimMdist):
          vv[i].append(Mdist[i][j])
          if(not Mdist[i][j]==-1):
            ens.add(Mdist[i][j])
    d = dict.fromkeys(ens)
    i=1
    i=(int)(i)
    for key in d.keys():
      d[key]=i
      i=i+1
    print (d)

    for i in range (0,dimMdist):
      for j in range (0,dimMdist):
        if(not Mdist[i][j]==-1):
          vv[i][j]=((int)(d[(Mdist[i][j])]))

    return vv , (len)(ens)

# fonction de création du vecteur P en dimension 2

def creervp_2d (mat,nbSsListes):
    p=[]
    for i in range(0,nbSsListes):
        p.append([])
    haut=len(mat)
    #print ("haut",haut)
    larg=len(mat[0])
    #print ("larg",larg)
    #print "dim = ", dim
    for i in range(0,len(mat)):
        L=mat[i]
        for j in range(0, len(L)):
            tup=(i,j)
            if not mat[i][j]==-1 :
                ind= (int)(mat[i][j]-1)
                (p[ind]).append(tup)
    return p

# fonction de création du vectteur Q en dimension 2

def creervq_2d(mat,vp,k,flag):
    # flag 0: a droite
    # flag 1: en bas
    q=[]
    for i in range(0,len(vp)):
        q.append([])
    
    for i in range(0,len(vp)):
        
        tmp=[]
        for n in range(0,len(vp)):
            tmp.append([])
        
        ll = vp[i]     
        if len(ll)>1:
            for j in range(0,len(ll)):
                #print (ll[j])
                if flag==0:
                	if ((ll[j])[1]+k <=len(mat[ (ll[j])[0] ])-1) :      		
	                    ind=(int)(mat[ (ll[j])[0] ] [ (ll[j])[1]+k ]-1)
	                    if (ind >-1):
		                    tmp[ind].append(ll[j])
                elif flag==1:
                  if ((ll[j])[0]+k<=len(mat)-1):
                   ind=(int)(mat[ (ll[j])[0]+k ] [ (ll[j])[1] ]-1)
                   #print (ind+1)
                   if (ind >-1):
                     tmp[ind].append(ll[j])
        
        for m in range(0,len(q)):
            if len(tmp[m])>1:
                q[m].append(tmp[m])    
    return q

# la fonction de calcul du vecteur V en 2 dimensions 

def creervv_dim2(vq , v ,dimMdist):
    #List[List[int]] * List[List[int] -> List[List[int]]
    if(dimMdist == 0):
        printf("la matrice est vide")
        return []
    if(vq is None):
        return None 
        
    for i in range (0,dimMdist):
        for j in range (0,dimMdist):
            v[i][j] = -1
    
    cpt = 0
    for i in range ( 0, len(vq)):
        L = vq[i]
        for a in L:
            for j in a:
                (a,b)= j
                v[a][b]=cpt
            cpt=cpt+1
    return v,cpt
        

def videvq(vq):
    cpt=0
    for l in vq:
        for sl in l:
            for couple in sl:
                cpt=cpt+1

    if (cpt==0):
        #print ("vide")
        return cpt
    else:
        #print ("nbCouple:", cpt) 
        return cpt


# ALGORITHME DE COMPARAISON


#2020 07/03

import sys
def comparerSequences_DIM2(MDD,lenVoulue=sys.maxsize):
    ''' List[List[int]] * int -> List[List[List[tuple(int,int)]]] or None
    MDD : matrice de distance discrétisée '''
    
    if(MDD) :
        #dimMdist : int
        dimMdist = len(MDD)
        # flag = int 
        flag =0

        # vv : List[List[int]]
        # nbGroupes : int 
        #vp : List[List[int]]
        #vq : List[List[List[tuple(int,int)]]]


        vv,nbGroupes = creer_vv_de_base(MDD, dimMdist)
        lenMotif=0
        nb_couples_vq=1
        
        while ((lenMotif<lenVoulue) and vv != None) :
            if ( flag == 0 ) :
                flag = 1
            if (flag == 1) :
                flag =0

            vp = creervp_2d (vv,nbGroupes)

            if (vp != None):
                vq = creervq_2d (vv,vp,1,flag)
            else : 
                print('erreur lors de la création du vecteur P')
                return None 

            if (vq != None):
                vv ,nbGroupes= creervv_dim2(vq , vv ,dimMdist)
            else :
                print('erreur lors de la création du vecteur Q')
                return None
  

            nb_couples_vq = videvq(vq)

            if (nb_couples_vq != 0):
                res = vq
            else :
                break 

            lenMotif += 1
            print("la longueur du motif est à : ",lenMotif)
        # on revoie le dernier vq avant qu'il ne soit vide
        print("pour une longueur de motif : ",lenMotif, "on a ces résultats :")
        return res , lenMotif
    else:

    	return None



def affichage_2d(vq,listAtom,lenMotif):
    ''' List[List[List[tuple(int,int)]]] * List[Atom] -> None 
    Cette fonction va permettre un affichage des résultats de la fonction de comparaison en 2 dim 
    Et rcupérer la commande d'execution Pymol pour repérer visuellement les séquences similaires '''
      
    for l in vq :
        for sl in l : 
            cpt = 1
            for couple in sl :
                print("\n Groupe ", cpt," :  " )
                (a,b) = couple
                print("select sa,",end='')
                for i in range(0,lenMotif) :
                    print ("resi",(listAtom[a+i]).numRes,"or ",end='')
                    if (i==lenMotif-1) :
                        print("resi",(listAtom[b+i]).numRes)
                    else :
                        print("resi",(listAtom[b+i]).numRes,"or ",end='')
                cpt += 1 

        print ("--------------------------------- \n")
            
