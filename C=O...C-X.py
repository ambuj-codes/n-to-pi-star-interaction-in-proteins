from math import *
import re
class atom:
    aid=0    #atomid initialisation 
    atype='' #atomtype initialisation
    x=0.0    #xcoord initialisation
    y=0.0    #ycoord initialisation
    z=0.0    #zcoord initialisation
    rid=0    #residueid initialisation
    rtype='' #residuetype initialisation

def getlen(atm1,atm2):
    dist=sqrt(pow(atm1.x-atm2.x,2)+pow(atm1.y-atm2.y,2)+pow(atm1.z-atm2.z,2)) #distance measurement between 1 and 2
    return dist #distance is returned

def getangle(atm1,atm2,atm3):
    dist1=sqrt(pow(atm1.x-atm2.x,2)+pow(atm1.y-atm2.y,2)+pow(atm1.z-atm2.z,2)) #distance measurement between 1 and 2
    dist2=sqrt(pow(atm3.x-atm2.x,2)+pow(atm3.y-atm2.y,2)+pow(atm3.z-atm2.z,2)) #distance measurement between 3 and 2
    dotp=(atm1.x-atm2.x)*(atm3.x-atm2.x)+(atm1.y-atm2.y)*(atm3.y-atm2.y)+(atm1.z-atm2.z)*(atm3.z-atm2.z) #BA.BC gives the angle 
    angle=acos(dotp/(dist1*dist2))*180/pi #angle measurement from dot product and changing to degree
    return angle #angle is returned

c_o_d_l=1.18 #carbon-oxygen double bond upper limit (residue i)
c_o_d_u=1.35 #carbon-oxygen double bond lower limit (residue i)
c_o_i_l=1.5 #oxygen-carbon bond lower limit (interaction between residue i and residue (j not equals i))
c_o_i_u=3.22 #oxygen-carbon bond upper limit (interaction between residue i and residue (j not equals i))
c_x_s_l=1.41 #carbon-x(x=c,o,n) double bond upper limit (residue j not equals i)
c_x_s_u=1.56 #carbon-x(x=c,o,n) double bond upper limit (residue j not equals i)
a_l=99 # c-o-c and o-c-x  minimum angle
a_u=119 # c-o-c and o-c-x maximum angle

filetxt=open('filelist.txt') #sample for test
txt_lines=filetxt.read().split('\n') #spaces are reduced
filetxt.close() #filelist edited and saved as filetxt
fileout=open('out_C_O.txt','w') #opened a text file for writing
f1=open('error_C_O.txt','w') #opened another text file for error writing

for filename in txt_lines:
    if filename=='':
        continue
    print('Executing for:',filename) #printing in output window the filenames the code is running for
    file=open(filename,'r') #file opened for reading
    lines=file.read().split('\n') #spaces reduced
    file.close() #file edited and closed
    
    C=[] #array containing all the carbon atoms
    O=[] #array containing all the oxygen atoms
    X=[] #array containing all the carbon, oxygen and the nitrogen atoms  

 
    try: #in our pdbs there are some anomalies, for normal file do as below or list error files
        for i in range(2,len(lines)): #first two lines of pdb are not needed
            linex=re.sub(' +',' ',lines[i]) 
            words=linex.split(' ')
            if words[0]=='ATOM' or words[0]=='HETATM':
                atm=atom()
                atm.aid=int(words[1]) #atom id
                atm.atype=words[2] #atom type
                atm.rtype=words[3] #residue type
                atm.rid=int(words[5]) #residue id
                atm.x=float(words[6]) #xcoord
                atm.y=float(words[7]) #ycoord
                atm.z=float(words[8]) #zcoord
                if words[11]=='C':
                    C.append(atm) #making the array of C
                    X.append(atm) #making the array of X
                elif words[11]=='O':
                    O.append(atm) #making the array of O
                    X.append(atm) #making the array of X
                elif words[11]=='N':
                    X.append(atm) #making the array of X
    except:
        f1.write(filename+'\n') #except does the left side thing and writes error files
        
    for i in range(len(C)): #first check of C
        for j in range(len(O)): #then check of oxygen
            if C[i].rid==O[j].rid and getlen(C[i],O[j])>=c_o_d_l and getlen(C[i],O[j])<=c_o_d_u: #interaction same residue and checks cond.
                for k in range(len(C)): #then another carbon
                    if C[k].rid!=O[j].rid and getlen(O[j],C[k])>=c_o_i_l and getlen(O[j],C[k])<=c_o_i_u: #interaction_diff_res and dif. cond.
                        for m in range(len(X)): #checks for O,C or N
                            if C[k].rid==X[m].rid and getlen(C[k],X[m])>=c_x_s_l and getlen(C[k],X[m])<=c_x_s_u: #same_res and cond. check
                                a1=getangle(C[i],O[j],C[k]) #c-o-c angle
                                a2=getangle(O[j],C[k],X[m]) #o-c-x angle
                                if a1>=a_l and a1<=a_u and a2>=a_l and a2<=a_u: #angle cond. check
                                    st='' #making a null string
                                    st+=filename+'\t' #storing a value to string st
                                    st+=O[j].atype+'\t' #adding atom type  
                                    st+=C[i].atype+'\t'  #adding atom type 
                                    st+=C[k].atype+'\t'  #adding atom type 
                                    st+=X[m].atype+'\t'  #adding atom type 
                                    st+=C[i].rtype+'\tA\t' #adding residue type 
                                    st+=str(C[i].rid)+'\t' # #adding residue id 
                                    st+=C[k].rtype+'\tA\t' # #adding residue type
                                    st+=str(C[k].rid)+'\t'  #adding residue id 
                                    st+=str(getlen(C[k],O[j]))+'\t' #getting length diff. residues
                                    st+=str(a1)+'\t' #getting c-o-c angle
                                    st+=str(a2)+'\n' #getting o-c-x angle
                                    fileout.write(st) #string st is completed and printed

fileout.close() #closing output file
f1.close() #closing error output file
