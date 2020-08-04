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
    model=[]
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
    

    return dih


def getdelta(atm1,atm2,atm3,atm4):
    ai=(atm2.y-atm1.y)*(atm3.z-atm2.z)-(atm3.y-atm2.y)*(atm2.z-atm1.z)
    bj=(atm2.z-atm1.z)*(atm3.x-atm2.x)-(atm2.x-atm1.x)*(atm3.z-atm2.z)
    ck=(atm2.x-atm1.x)*(atm3.y-atm2.y)-(atm2.y-atm1.y)*(atm3.x-atm2.x)
    pd=(ai*atm4.x+bj*atm4.y+ck*atm4.z-(ai*atm1.x+bj*atm1.y+ck*atm1.z))/sqrt(pow(ai,2)+pow(bj,2)+pow(ck,2))
    pd1=abs(pd)
    return pd1

c_o_c_o= 3.22 # O---C=O distance

a_l=99     # n-pi* interactions criteria
a_u=119     # n-pi* interactions criteria

filetxt=open('filelist.txt') 
txt_lines=filetxt.read().split('\n') 
filetxt.close()
fileout=open('out_n-pi-star.txt','w')
f1=open('error_out_n-pi-star.txt','w')
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
    O1=[]
    O2=[]
    C1=[] 
    CA=[]
    N=[]
    modelno=[]

 
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
                atm.z=float(ln[47:54]) 
                atm.model=modelno
                symb=ln[76:78].strip()
                #print(atm.atype)
                if atm.atype=='C' and atm.rid > 70 and atm.rid < 80:
                    C1.append(atm)
                if atm.atype=='O' and atm.rid > 70 and atm.rid < 80 :
                    O1.append(atm)
                    O2.append(atm)
                if atm.atype=='CA' and atm.rid > 70 and atm.rid < 80 :
                    CA.append(atm)
                if atm.atype=='N'and atm.rid > 70 and atm.rid < 80:
                    N.append(atm)
                        
               
            elif len(ln)>=5 and ln[0:5]=='MODEL':
                modelno=int(ln[12:])

    except:
        f1.write(filename+'\n')
    print(len(C1))
    for c1 in range(len(C1)):
        for o1 in range(len(O1)):
            if C1[c1].rid==O1[o1].rid and C1[c1].chainid==O1[o1].chainid and getlen(C1[c1],O1[o1])>=1.18  and getlen(C1[c1],O1[o1])<=1.38:
                for o2 in range(len(O2)):
                    if O2[o2].aid!=O1[o1].aid :
                        for ca in range(len(CA)):
                            for n in range(len(N)):
                                if C1[c1].rid==O1[o1].rid==CA[ca].rid==N[n].rid and C1[c1].chainid==O1[o1].chainid==CA[ca].chainid==N[n].chainid and getlen(O2[o2],C1[c1])<=c_o_c_o and getangle(O2[o2],C1[c1],O1[o1])>= a_l and getangle(O2[o2],C1[c1],O1[o1])<= a_u :
                                    intr.append([])
                                    intr[len(intr)-1].append(filename)                                                     
                                    intr[len(intr)-1].append(O2[o2].chainid )               
                                    intr[len(intr)-1].append(O2[o2].rtype)
                                    intr[len(intr)-1].append(O2[o2].rid )
                                    intr[len(intr)-1].append(O2[o2].atype)
                                    intr[len(intr)-1].append(O1[o1].atype)
                                    intr[len(intr)-1].append(O1[o1].aid)
                                    intr[len(intr)-1].append(O1[o1].rid)        
                                    intr[len(intr)-1].append(O1[o1].rtype)         
                                    intr[len(intr)-1].append(O1[o1].chainid)   
                                    intr[len(intr)-1].append(getlen(O2[o2],C1[c1]))
                                    intr[len(intr)-1].append(getangle(O2[o2],C1[c1],O1[o1]))
                                    intr[len(intr)-1].append(getdelta(CA[ca],N[n],O1[o1],C1[c1]))

    O1=[]
    O2=[]
    C1=[] 
    CA=[]
    N=[]
    for line in intr:
        for xxd in line:
            fileout.write(str(xxd))
            fileout.write('\t')
        fileout.write('\n')
    intr=[]
    fileout.close()
    fileout=open('out_n-pi-star.txt','a')
fileout.close()
f1.close()
