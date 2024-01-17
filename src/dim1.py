# -*- coding: utf-8 -*
import sys

def init_alpha(alphabet):
  '''List[str]->Dictionary'''
  #D: Dictionary
  D= dict()

  #i:int
  i=0

  #lettre : str
  for lettre in alphabet:
    D[lettre]=i
    i=i+1
  
  return D

def creerTabBase (listeSeq,alphabet):
    '''List[List[str]]*Liste[str]-> List[int] or None'''
    
    if listeSeq and alphabet:

      #alpha_dic: dict()
      alpha_dic=init_alpha(alphabet)
      
      #tabRes:List[int] contiendra le tableau résultant
      tabRes=[]
  
      #on parcourt la liste de séquences
      for seq in listeSeq: 

        # a chaque lettre de la séquence, on affecte un indice conformément à l'indexaction de l'alphabet
        for c in seq:
          if alpha_dic.get(c.upper()) !=None:
            tabRes.append(alpha_dic.get(c.upper())) 
          else:
            print(c)
            print("ERREUR. L'alphabet n'est pas conforme pour les séquences à étudier")
            return None
              
      return tabRes
    else:
      return None

def creervp  (vecteurV,nbClasses):
    '''List[int] * int -> List[List[int]] or None'''
    
    # vérification de type
    if vecteurV and isinstance(nbClasses, int) and nbClasses>0:
      
      # p : List[List[int]]
      p=[[]for k in range(nbClasses)]
      
      # on parourt le vecteur v, on ajoute i au vecteur P à l'indice vecteurV[i]
      for i in range (0,len(vecteurV)) :
        if vecteurV[i]>=0 : 
          p[vecteurV[i]].append(i)
        
      return p
    else:
      return None

def creervq (vp , vv ,BordSeq, k ) :
    ''' List[List[int]] * List[int] *int  -> List[List[int]] or None'''

    #vérification 
    if vv and vp and isinstance(k, int) and k>0:
      #len_vv:int
      len_vv=len(vv)

      #q:list[list[int]]
      q=[[]  for k in range(len(vp))]
      
      # Pour chaque sous liste de vp on effectue des dépilements
      for i in range (0,len(vp)) :
          #L:list[int]
          L = vp[i]
          
          #séparateur à chaque nouvelle sous-liste
          for j in range (0,len(q)) :
              q[j].append(-1)

          #dépilement de la sous liste vp[i]
          for j in range (len(L)-1,-1,-1) :
              if (L[j]+k)<len_vv and L[j]+k>=0:
                if vv[L[j]+k]>=0:
                  if L[j] not in BordSeq:
                    q[vv[L[j]+k]].append(L.pop(j))
                    
                   
          #L.clear()
          del L[:]
      for j in range (0,len(q)) :
              q[j].append(-1) 
      return q  
    else: 
      return None

def creervv (vq , v):
    ''' List[List[int]] * List[int] ->int* List[int]*List[int] or None'''

    if vq and v:
      #vecteurV : List[int]
      vecteurV=[]*len(v)
      for i in range(len(v)):
          vecteurV.append(-2)#case vide

      #cpt:int , compteur de nombre groupes
      cpt=-1 
      #pos:List[nuplet] , contiendra la liste des listes des motifs de taille maximale
      pos=[[] for k in range(len(v))]
      #ipos:int , indice pour le tableau pos
      ipos=-1

      for i in range(0,len(vq)): # parcourt du Vecteur Q
          L=vq[i]
        
          for j in range (len(L)-2,0,-1) : 
            if L[j]!=-1:
              elemavant=L[j+1] #dépilement de droite à gauche ex: 4 5 6 (6 est avant 5 et 4 est après 5)
              elemapres=L[j-1]

              if (elemapres == -1 and elemavant==-1): # motif de taile 1 
                cpt=cpt+1
                vecteurV[L[j]]=cpt
                
              else:
                if elemavant==-1: # premier élément d'un motif de taille > 1, il appartiendra au groupe cpt+1
                  cpt=cpt+1
                  ipos=ipos+1

                vecteurV[L[j]]=cpt # autres éléments de motif de taille >1, il appartiennent au même groupe que le précedent cpt+1
                pos[ipos].append(L[j]) # on ajoute le motif de taille > 1 à la liste des positions de motifs
                

          #L.clear()
          del L[:]

      pos=list(filter(None,pos)) #enlève les listes vides
      return cpt+1,pos,vecteurV # nombre de groupes obtenus, listes des groupes de motifs à taille maximale,vecteur v

    else:
      return 0,None


def ListToStr(listeSeq,k):
  #Seq:str
  Seq=""
  #taille: int
  taille=len(listeSeq)
  #Ext:List[int]
  Ext=[]
  if listeSeq :
    
    for i in range (taille):
      Seq=Seq+listeSeq[i]

      for j in range (k):
        Ext.append(len(Seq)-1-j)


  return Ext,Seq

def comparerSequences(listeSeq,alphabet,lenVoulue=sys.maxsize , k=1,aff=1):
    ''' List[str]*Dictionary*int -> List[nuplet] or None'''
    
    if listeSeq and alphabet:
      #AllSeq:str , contient l'ensemble des sequences sous forme de chaine
      #BordSeq:List[int]
      BordSeq,AllSeq=ListToStr(listeSeq,k)

      #vv : List[int]
      vv=creerTabBase(listeSeq,alphabet)

      #nbGroupes:int
      nbGroupes=len(alphabet)

      #lenMotif:int
      lenMotif=1
      
      #listePos:List[int]
      listePos=[]
      
      while (vv!=None and lenMotif<lenVoulue) :
        vp = creervp (vv , nbGroupes) 
        
        if vp!=None:
          vq = creervq (vp , vv ,BordSeq, k)    
          
          
          if vq!=None:
            nbGroupes,listePos,vv = creervv (vq , vv)

            if listePos==[] : #Condition d'arrêt de la boucle: plus de motifs trouvés
              break

            lenMotif=lenMotif+1         

            # à chaque tour de boucle afficher les motifs similaires de longueur lenMotif
            if (aff):
              print("\n_____________________\nTAILLE DU MOTIF:",lenMotif,"\n---------------------\n")
              print("MOTIF->POSITION(S)\n---------------------")
            for i in range (0,len(listePos)):
              P=listePos[i]
              if (aff):
                print(AllSeq[P[0]:P[0]+lenMotif],"->",P)
            
            
                
      if lenMotif==1:
        print("Pas de motifs de taille supérieure à 1.\n")
       
    else:
      ("PARAMETRES NON COHERENTS DANS comparerSequences() \n")


  
