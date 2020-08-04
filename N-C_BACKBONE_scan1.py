# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 12:35:12 2018

@author: suman
"""
from math import *
import math
import re
import numpy as np
class atom:
    aid=0    
    atype='' 
    x=0.0  
    y=0.0    
    z=0.0    
    rid=0    
    rtype='' 
    model=0
    chainid=''

def getlen(atm1,atm2):
    dist=sqrt(pow(atm1.x-atm2.x,2)+pow(atm1.y-atm2.y,2)+pow(atm1.z-atm2.z,2)) 
    return dist

def getangle(atm1,atm2,atm3):
    dist1=sqrt(pow(atm1.x-atm2.x,2)+pow(atm1.y-atm2.y,2)+pow(atm1.z-atm2.z,2)) 
    dist2=sqrt(pow(atm3.x-atm2.x,2)+pow(atm3.y-atm2.y,2)+pow(atm3.z-atm2.z,2)) 
    dotp=(atm1.x-atm2.x)*(atm3.x-atm2.x)+(atm1.y-atm2.y)*(atm3.y-atm2.y)+(atm1.z-atm2.z)*(atm3.z-atm2.z) 
    angle=acos(dotp/(dist1*dist2))*180/pi 
    return angle

def getangledihedral(atm1,atm2,atm3,atm4):
    ab=np.zeros(3)
    bc=np.zeros(3)
    cd=np.zeros(3)
    p=[]
    q=[]
    ab[0]=atm2.x-atm1.x
    ab[1]=atm2.y-atm1.y
    ab[2]=atm2.z-atm1.z
    bc[0]=atm3.x-atm2.x
    bc[1]=atm3.y-atm2.y
    bc[2]=atm3.z-atm2.z
    cd[0]=atm4.x-atm3.x
    cd[1]=atm4.y-atm3.y
    cd[2]=atm4.z-atm3.z
    p.append(ab[1]*bc[2]-ab[2]*bc[1])
    p.append(ab[2]*bc[0]-ab[0]*bc[2])
    p.append(ab[0]*bc[1]-ab[1]*bc[0])
    q.append(bc[1]*cd[2]-bc[2]*cd[1])
    q.append(bc[2]*cd[0]-bc[0]*cd[2])
    q.append(bc[0]*cd[1]-bc[1]*cd[0])


    r1=0
    r2=0
    dp=0
    dpcd=0
    for i in range(0,3):
        r1 += math.pow(p[i],2)
        r2 += math.pow(q[i],2)
        dp += p[i]*q[i]
        dpcd += p[i]*cd[i]

    dih=(dpcd/abs(dpcd))*math.acos(dp/(math.sqrt(r1)*math.sqrt(r2)))*180/math.pi
    dih2=abs(dih)

    return dih2

n_c_i_l=1.8
n_c_i_u=3.0
a_l_1=0 
a_u_1=10
a_l_2=170
a_u_2=180

filetxt=open('filelist.txt') 
txt_lines=filetxt.read().split('\n') 
filetxt.close()
fileout=open('out_H-N...C=O_scan1.txt','w')
f1=open('error_H-N...C=O_scan1.txt','w')
intr=[]
lenlines=len(txt_lines)
for ppp in range(lenlines):
    filename=txt_lines[ppp]
    if filename=='':
        continue
    print('%.2f'%((ppp+1)*100.0/(lenlines-1))+'% ('+str(ppp+1)+'/'+str(lenlines-1)+')  Executing for:'+filename)
    file=open(filename,'r')
    lines=file.read().split('\n')
    file.close()
    
    N=[] 
    H=[] 
    C=[]
    O=[]
    modelno=0

 
    try:
        for ln in lines:
            if len(ln)>=6 and (ln[0:4]=='ATOM' or ln[0:6]=='HETATM'):
                atm=atom()
                atm.aid=int(ln[6:11]) 
                atm.atype=ln[12:16].strip() 
                atm.rtype=ln[17:20].strip() 
                atm.chainid=ln[21]
                atm.rid=int(ln[22:26]) 
                atm.x=float(ln[30:38]) 
                atm.y=float(ln[38:46]) 
                atm.z=float(ln[46:54]) 
                atm.model=modelno
                symb=ln[76:78].strip()
                if atm.atype=='H':
                    if atm.rtype=='ARG' or atm.rtype=='HIS' or atm.rtype=='LYS' or atm.rtype=='ASP' or atm.rtype=='GLU' or atm.rtype=='SER' or atm.rtype=='THR' or atm.rtype=='ASN' or atm.rtype=='GLN' or atm.rtype=='CYS' or atm.rtype=='SEC' or atm.rtype=='GLY' or atm.rtype=='PRO' or atm.rtype=='ALA' or atm.rtype=='MET' or atm.rtype=='ILE' or atm.rtype=='TYR' or atm.rtype=='VAL' or atm.rtype=='LEU' or atm.rtype=='PHE' or atm.rtype=='TRP' or atm.rtype=='MSE': 
                        H.append(atm)
                if atm.atype=='N':
                    if atm.rtype=='ARG' or atm.rtype=='HIS' or atm.rtype=='LYS' or atm.rtype=='ASP' or atm.rtype=='GLU' or atm.rtype=='SER' or atm.rtype=='THR' or atm.rtype=='ASN' or atm.rtype=='GLN' or atm.rtype=='CYS' or atm.rtype=='SEC' or atm.rtype=='GLY' or atm.rtype=='PRO' or atm.rtype=='ALA' or atm.rtype=='MET' or atm.rtype=='ILE' or atm.rtype=='TYR' or atm.rtype=='VAL' or atm.rtype=='LEU' or atm.rtype=='PHE' or atm.rtype=='TRP' or atm.rtype=='MSE': 
                        N.append(atm) 
                if atm.atype=='C':
                    if atm.rtype=='ARG' or atm.rtype=='HIS' or atm.rtype=='LYS' or atm.rtype=='ASP' or atm.rtype=='GLU' or atm.rtype=='SER' or atm.rtype=='THR' or atm.rtype=='ASN' or atm.rtype=='GLN' or atm.rtype=='CYS' or atm.rtype=='SEC' or atm.rtype=='GLY' or atm.rtype=='PRO' or atm.rtype=='ALA' or atm.rtype=='MET' or atm.rtype=='ILE' or atm.rtype=='TYR' or atm.rtype=='VAL' or atm.rtype=='LEU' or atm.rtype=='PHE' or atm.rtype=='TRP' or atm.rtype=='MSE': 
                        C.append(atm)
                if atm.atype=='O':
                    if atm.rtype=='ARG' or atm.rtype=='HIS' or atm.rtype=='LYS' or atm.rtype=='ASP' or atm.rtype=='GLU' or atm.rtype=='SER' or atm.rtype=='THR' or atm.rtype=='ASN' or atm.rtype=='GLN' or atm.rtype=='CYS' or atm.rtype=='SEC' or atm.rtype=='GLY' or atm.rtype=='PRO' or atm.rtype=='ALA' or atm.rtype=='MET' or atm.rtype=='ILE' or atm.rtype=='TYR' or atm.rtype=='VAL' or atm.rtype=='LEU' or atm.rtype=='PHE' or atm.rtype=='TRP' or atm.rtype=='MSE': 
                        O.append(atm)
            elif len(ln)>=5 and ln[0:5]=='MODEL':
                modelno=int(ln[12:])

    except:
        f1.write(filename+'\n')
    
    for j in range(len(H)):
        for k in range(len(N)):
            if H[j].rid==N[k].rid:
                for m in range(len(C)):
                    if getlen(N[k],C[m])>=n_c_i_l and getlen(N[k],C[m])<=n_c_i_u:
                        for n in range(len(O)):
                            if O[n].rid==C[m].rid:
                                a1=getangledihedral(H[j],N[k],C[m],O[n]) 
                                if (a1>=a_l_1 and a1<=a_u_1) or (a1>=a_l_2 and a1<=a_u_2):
                                    flag=True
                                    for i in range(0,len(intr)):
                                        if intr[i][0]==filename and intr[i][3]==H[j].rid and intr[i][7]==C[m].rid:
                                            flag=False
                                            if getlen(N[k],C[m])<intr[i][9]:
                                                intr[i][0]=filename                                            
                                                intr[i][1]=H[j].atype               
                                                intr[i][2]=H[j].rtype              
                                                intr[i][3]=H[j].rid                
                                                intr[i][4]=N[k].atype               
                                                intr[i][5]=C[m].atype                
                                                intr[i][6]=C[m].rtype                
                                                intr[i][7]=C[m].rid                  
                                                intr[i][8]=getlen(H[j],N[k])        
                                                intr[i][9]=getlen(C[m],N[k])         
                                                intr[i][10]=getlen(C[m],O[n])       
                                                intr[i][11]=a1                       
                                    if flag==True:
                                        intr.append([])
                                        intr[len(intr)-1].append(filename)                                        
                                        intr[len(intr)-1].append(H[j].atype)            
                                        intr[len(intr)-1].append(H[j].rtype)              
                                        intr[len(intr)-1].append(H[j].rid)                
                                        intr[len(intr)-1].append(N[k].atype)               
                                        intr[len(intr)-1].append(C[m].atype)                
                                        intr[len(intr)-1].append(C[m].rtype)              
                                        intr[len(intr)-1].append(C[m].rid)               
                                        intr[len(intr)-1].append(getlen(H[j],N[k]))        
                                        intr[len(intr)-1].append(getlen(C[m],N[k]))         
                                        intr[len(intr)-1].append(getlen(C[m],O[n]))      
                                        intr[len(intr)-1].append(a1)                      

    H=[]
    N=[]
    C=[]
    O=[]
    for line in intr:
        for xxd in line:
            fileout.write(str(xxd))
            fileout.write('\t')
        fileout.write('\n')
    intr=[]
    fileout.close()
    fileout=open('out_H-N...C=O_scan1.txt','a')
fileout.close()
f1.close()
