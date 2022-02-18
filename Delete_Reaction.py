from tkinter import *                      # Importing tkinter for designing the UI
from tkinter.messagebox import *           # Importing tkinter.messagebox for warning/cofirmation message
import os                                  # Importing os to run system command
import pandas                              # Importing pandas to for datastructure
import csv                                 # Importing csv to read/write csv file
from itertools import permutations         # Importing permutation for permuting list
import ast                                 # Importing ast to process string 
import webbrowser

def read_manual():
    webbrowser.open_new('Manual_ESMILES.pdf')

def check_stability(input_esmiles):          #Method to process products predicted by the algorithm
    
    poss_symbols=['-','=','~','/','#','(',')','[',']','*','.','+','{','}']
    
    for ed in input_esmiles:                  #Loop to process the educt string input given by the user
        

        if (ed.find(' + ') == -1):
            ed = [ed]
        else:
            ed = ed.split(' + ')

        ed_compounds = ed           #Storing the individual Educt compound in a separate variables
        s = ''
        stack=[]                    #Declaring an empty Stack
        top=-1

        total_ed_matrix=[]          #List for storing EBE matrices of individual educt compounds

        ed_elements=[]              #List for storing elements of educts

        
        rem_flag=0
        ed_compounds = ed                    #Storing the individual Educt compound in a separate variables
        s = ''
        stack=[]                            #Declaring an empty Stack
        top=-1

        total_ed_matrix=[]                  #List for storing EBE matrices of individual educt compounds

        ed_elements=[]


        for react in ed:
            s=react[1:len(react)-1]
            c=0
            for j in range(len(s)):
                
                if s[j] not in poss_symbols and s[j].isdigit()==0 and s[j].isalpha()==0: 
                    #print(s[j])      
                    return 0

                if (j+1)<len(s) and s[j].isupper() and s[j+1].islower():
                    c+=1
                elif s[j].isupper():
                    c+=1
                    
            

            for j in range(len(s)):

                if (j+1)<len(s) and s[j].isupper() and s[j+1].islower():
                    ed_elements.append(s[j]+s[j+1])
                elif s[j].isupper():
                    ed_elements.append(s[j])


            
            ed_matrix=[[0 for q in range(c)] for z in range(c)]
            top=-1
            stack=[-1 for i in range(2*c)]
            j=0
            k=0
            cnt=''
            chg_cnt_p=0
            chg_cnt_n=0
            flag=0
            add_flag=0
            add_done=0
            top=top+1
            stack[top]=j
            m=0
            cyc_hash=[]                              #List for cyclic compounds
            cyc_hash.append([])
            cyc_hash[0].append(0)
            cyc_bond=[]
            cyc_bond.append([])
            cyc_bond[0].append(0)


            while k<len(s):
                m=0

                if k+1<len(s) and s[k].isupper() and s[k+1].islower() and cnt=='':              #Identifying Chemical Symbol in educt string
                    cnt=s[k]+s[k+1]
                    k=k+1

                elif k+1<len(s) and s[k].isupper() and cnt=='':                                 #Identifying Chemical Symbol in educt string 
                    cnt=s[k]
                    k=k+1

                elif k+1<len(s) and s[k].isupper() and s[k+1].islower() and cnt!='':            #Identifying Chemical Symbol in educt string
                    j=j+1
                    top=top+1
                    stack[top]=j
                    
                    
                        

                    k=k+1

                elif k+1<len(s) and s[k].isupper() and cnt!='':                                 #Identifying Chemical Symbol in educt string
                    j=j+1
                    top=top+1
                    stack[top]=j
                    
                    

                    k=k+1
                    

                elif k+1<len(s) and s[k]=='.':                                                  #Identifying valency of an element
                     for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            m=m*10+int(s[l])
                        else:
                            break

                     t1=top
                     while stack[t1]==-100:
                         t1=t1-1
                     p=stack[t1]
                     ed_matrix[p][p]=m
                     k=k+1

                elif k+1<len(s) and (s[k]=='+' or s[k]=='-') and s[k+1].isdigit():              #Identifying charge of an element
                    cov='1'
                    for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            m=m*10+int(s[l])
                        else:
                            break

                        if flag==0:
                             if s[k]=='+':
                                 t1=top
                                 while stack[t1]==-100:
                                     t1=t1-1
                                 p=stack[t1]
                                 ed_matrix[p][p]=m*(10)+ed_matrix[p][p]

                             elif s[k]=='-':
                                 t1=top
                                 while stack[t1]==-100:
                                    t1=t1-1
                                 p=stack[t1]
                                 ed_matrix[p][p]=m*(-10)-ed_matrix[p][p]

                    k=k+1

                elif k+1<len(s) and s[k]=='[' and s[k+1].isalpha() and flag==0:                  #Identifying ionic compound
                    cnt=''
                    chg_cnt_p=0
                    chg_cnt_n=0
                    cov='1'

                    if flag==0:
                        flag=1

                    for d in range(k+1,len(s)):
                        if d+1<len(s) and s[d].isupper() and s[d+1].islower():
                            chg_cnt_p+=1
                        elif d+1<len(s) and s[d].isupper():
                            chg_cnt_p+=1
                        elif s[d]==']':
                            break



                    for q in range(chg_cnt_p):
                        ed_matrix[q][q]=10
                        for w in range(chg_cnt_p,c):
                            ed_matrix[q][w]=10



                    for t in range(d+1,len(s)):
                        if t+1<len(s) and s[t].isupper() and s[t+1].islower():
                            chg_cnt_n+=1
                        elif t+1<len(s) and s[t].isupper():
                            chg_cnt_n+=1
                        elif s[t]==']':
                            break


                    for q in range(chg_cnt_p,c):
                        ed_matrix[q][q]=-10
                        for w in range(chg_cnt_p):
                            ed_matrix[q][w]=-10



                    k=k+1


                elif k+1<len(s) and s[k]==']' and flag==1 and s[k+1]=='[':                      #Identifying Ionic compounds
                    m=0
                    flag=2
                    k=k+1
                    j=j+1
                    top=top+1
                    stack[top]=j
                    cnt=''

                elif k+1<len(s) and s[k]==']' and flag==2:                                      #Identifying Ionic compounds
                    m=0
                    flag=0
                    k=k+1
                    j=j+1
                    top=top+1
                    stack[top]=j
                    cnt=''

                elif k+1<len(s) and s[k]=='[' and s[k+1].isdigit():                             #Identifying Cyclic compounds

                    for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            m=m*10+int(s[l])
                        else:
                            break

                    if len(cyc_hash)<=m:
                         cyc_hash.append([])
                         cyc_hash[m].append(j)
                         if s[l]=='-':
                             cyc_bond.append([])
                             cyc_bond[m].append(1)
                         elif s[l]=='=':
                             cyc_bond.append([])
                             cyc_bond[m].append(2)
                         elif s[l]=='#':
                             cyc_bond.append([])
                             cyc_bond[m].append(3)
                         elif l+1<len(s) and s[l]=='~' and s[l+1]=='/':
                             cyc_bond.append([])
                             cyc_bond[m].append(200)

                         elif s[l]=='~':
                             cyc_bond.append([])
                             cyc_bond[m].append(100)


                    else:
                        cyc_hash[m].append(j)

                    k=k+1
                    flag=-1

                elif k+1<len(s) and s[k]==']' and flag==-1:                                 #Identifying cyclic compounds
                    flag=0
                    k=k+1





                elif s[k]=='-' and flag!=-1:                                                #Identifying single bond
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=1
                    ed_matrix[j][p]=1
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                   


                elif s[k]=='=' and flag!=-1:                                                #Identifying double bond
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=2
                    ed_matrix[j][p]=2
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                    
                              
                elif s[k]=='#' and flag!=-1:                                                #Identifying triple bond
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=3
                    ed_matrix[j][p]=3
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                    
                              
                elif s[k]=='~' and s[k+1]=='/' and flag!=-1:                                #Identifying coordinate bond
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=200
                    ed_matrix[j][p]=100
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                   
                elif s[k]=='~' and flag!=-1:                                            #Identifying coordinate bond                                 
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=100
                    ed_matrix[j][p]=200
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                    
                
                    

                elif s[k]=='(':                                                         #Identifying branched compound
                    
                    if s[k+1].isalpha():
                        j=j+1
                        add_flag=1
                        add_top=top
                        while stack[add_top]==-100:
                            add_top=add_top-1
                        
                        ed_matrix[0][stack[add_top]+1]=500
                        top=top+1
                        stack[top]=j
                    
                    else:                       
                        top=top+1
                        stack[top]=-100
                        
                        if add_flag>=1:
                            add_flag+=1
                            
                        
                        
                
                    k=k+1
                    cnt=''
                    

                elif k+1<len(s) and s[k]==')' and s[k+1].isalpha():
                    k=k+1

                elif s[k]==')':                                                         #Identifying branched compound


                    if add_flag==1:
                        add_flag=0
                        
                    elif add_flag!=1:                    
                    
                        while stack[top]!=-100:
                            top=top-1
                        top=top-1

                        if add_flag>=1:
                            add_flag-=1

                        
                        
                        
                    k=k+1
                    cnt=''
                        

                else:
                    k=k+1
                    
            if len(cyc_hash)>1:                                                         #Processing cyclic compounds
                for i in range(1,len(cyc_hash)):
                    p1=cyc_hash[i][0]
                    p2=cyc_hash[i][1]
                    bond=cyc_bond[i][0]

                    ed_matrix[p1][p2]=bond
                    ed_matrix[p2][p1]=bond

                    if bond==100:
                        ed_matrix[p1][p2]=100
                        ed_matrix[p2][p1]=200
                    elif bond==200:
                        ed_matrix[p1][p2]=200
                        ed_matrix[p2][p1]=100


            total_ed_matrix.append(ed_matrix)                                           #Storing the educt matrix in a list
        
        
        
       
        c = 0
        ed_element_slno=[]
        for i in ed_elements:
            c = c+1
            
            ed_element_slno.append(c)
            
        
        ed_len=len(ed_element_slno)
        ebe_matrix_educts = [[0 for i in range(ed_len+2)] for j in range(ed_len)]           #EBE matrix of educts
    
       
    
        c = 0
        for i in total_ed_matrix:                        #Loop to fill the EBE matrix of educts
            for j in range(len(i)):
                for k in range(len(i)):
                    x = ed_element_slno[c+j]-1
                    y = ed_element_slno[c+k]-1
                    ebe_matrix_educts[x][y] = i[j][k]
    
                    if ebe_matrix_educts[x][y] == 100:                      #For Coordinate bond
                        if ebe_matrix_educts[x][ed_len+1]!=0:
                            ebe_matrix_educts[x][ed_len+1]=-99
                        else:
                            ebe_matrix_educts[x][ed_len+1] = y+1
    
                        if ebe_matrix_educts[y][ed_len+1]!=0:
                            ebe_matrix_educts[y][ed_len+1]=-99
                        else:
                            ebe_matrix_educts[y][ed_len+1] = x+1
    
    
                    if ebe_matrix_educts[x][y] >= 10 and ebe_matrix_educts[x][y] <= 88:             #For cations
                        ebe_matrix_educts[x][ed_len] = ebe_matrix_educts[x][y]//10
                        if x != y:
                            ebe_matrix_educts[x][y] = 1000
                            ebe_matrix_educts[y][x] = 2000
                        else:
                            ebe_matrix_educts[x][x] = ebe_matrix_educts[x][x] % 10
    
                    if ebe_matrix_educts[x][y] >= -88 and ebe_matrix_educts[x][y] <= -10:           #For anions
                        ebe_matrix_educts[x][ed_len] = ebe_matrix_educts[x][y]//10
                        if x != y:
                            ebe_matrix_educts[x][y] = 2000
                            ebe_matrix_educts[y][x] = 1000
                        else:
                            ebe_matrix_educts[x][x] = (ebe_matrix_educts[x][x]*-1) % 10

                    
    
            c = c+len(i)

       
        
        sorted_educt_matrix=ebe_matrix_educts[:]                #Sorting the ebe matrix of educts based on chemical symbol 
      
    
        for i in range(ed_len-1):
            for j in range(ed_len-i-1):
                ch1=ed_elements[j][0]
                ch2=ed_elements[j+1][0]
    
                if(ch1>ch2):
                    temp=ed_elements[j]
                    ed_elements[j]=ed_elements[j+1]
                    ed_elements[j+1]=temp
    
                    sorted_educt_matrix[j], sorted_educt_matrix[j+1]=sorted_educt_matrix[j+1], sorted_educt_matrix[j]

    
                    temp_list=[sorted_educt_matrix[k][j] for k in range(ed_len)]
    
                    for k in range(ed_len):
                        sorted_educt_matrix[k][j]=sorted_educt_matrix[k][j+1]
    
                    for k in range(ed_len):
                        sorted_educt_matrix[k][j+1]=temp_list[k]
    
                        
        
        ebe_matrix_educts=sorted_educt_matrix[:]
    
        #print(ed)
        #print(ed_elements)
        #print(ebe_matrix_educts)
    
        max_bond=[]
        electronegativity=[]
        
        ele_ex=0
        
        for i in ed_elements:
            ele_ex=0
            for j in range(len(elements)):
                if i==elements[j]:
                    ele_ex=1
                    max_bond.append(elements_max_bond[j])
                    electronegativity.append(elements_eneg[j])
            if ele_ex==0:
                return 0
                
        
        
        
        #print(free_electrons)
        #print(stable_electrons)

        
       
        rem_flag=0
        for i in range(len(ebe_matrix_educts)):
           ionic=0
           shared=0
           
           for j in range(len(ebe_matrix_educts[i])-2):
               if ebe_matrix_educts[i][j]==1000:
                   ionic=1
                   break
               elif ebe_matrix_educts[i][j]==2000:
                   ionic=2
                   break
               
          
           for j in range(len(ebe_matrix_educts[i])-2):
               if i!=j:
                    a=ebe_matrix_educts[i][j]
                    if a==100:
                        a=1
                    
                    elif a==500 or a==1000 or a==2000 or a==200:
                        a=0
                
                    shared+=a
                
           if shared>max_bond[i]:
                #print(shared)
                rem_flag=1
                break
                       
        if rem_flag==1:
            return 0
        
    return 1
            
            
    



def delete_reaction():          #Method to delete reaction from Reaction database
    t1=''
    t2=''
    p=[]
    f="Reactions_"+str(cov)+str(parts)+str(bnd)+".csv"      #File containing the particular reaction
    fd=0
    
    if os.path.isfile(f):
        reader = open(f, "r")
        lines = reader.read().split("\n")       #Reading reactions from the file
        reader.close()
        for i in lines:
            if i:            
                s=list(ast.literal_eval(i))             
                
                if s[0]==educt.get() and s[3]==product.get():       #Matching educt and product with user input 
                    t1=s[4]
                    t2=s[5]
                    p=s[:]
                    fd=1
                    break
        
        if fd==1:
            del_temp=0
            writer = open(f,"w")                                
            for i in lines:                                         #Removing reaction when there is a match
                if i:            
                    s=list(ast.literal_eval(i))
                   
                    if p!=s:                         
                        if t1==s[4] and t2==s[5]:
                            del_temp=1
                    
                        writer.write(i+"\n")
                    
                    
            writer.close()
            
                    
            if del_temp==0:
                f="Templates_"+str(cov)+str(parts)+str(bnd)+".csv"        #Template file containing the template of the reaction
                if os.path.isfile(f):
                    reader = open(f, "r")
                    lines = reader.read().split("\n")
                    reader.close()
                    writer = open(f,"w")
    
                    for i in lines:                                      #Removing templates from templates
                        if i:            
                            s=list(ast.literal_eval(i))
                            if s[0]==t1 and s[4]==t2:
                                continue
                            else:
                                writer.write(i+"\n")
                    writer.close()
                
        if fd==0:
            showerror("Error", "No Chemical Reaction Found")
        elif fd==1:
            
            showinfo("Deletion Done", "Reaction has been deleted from the database")
            
                
    
    
    
def go_press():  #Function to process the Educts and Products entered in the textbox

    global ed, prod, ed_elements,prod_elements,total_ed_matrix,total_prod_matrix, ed_compounds,prod_compounds,cov,bnd,parts,filename
    cov='0'
    bnd='1'
    ed = educt.get()            #Educt input given by the user
    prod = product.get()        #Product input given by the user
    f1 = 0
    f2 = 0
    if not ed or not prod:
        showwarning("Warning", "Please Fill Both Educt and Product Part")       #Warning if Product or Educt is blank
    
    elif check_stability([ed,prod])!=1:
        showwarning("Warning", "Please give corrct educts or products")       #Warning if Educt or Product is unstable

    else:                               #Extracting Individual compound in Educts
        if (ed.find(' + ') == -1):
            parts='1'
            ed = [ed]
        else:
            ed = ed.split(' + ')
            if len(ed)==2:
                parts='2'
            if len(ed)==3:
                parts='3'
            if len(ed)>3:
                parts='4'

        ed_compounds = ed           #Storing the individual Educt compound in a separate variables
        s = ''
        stack=[]                    #Declaring an empty Stack
        top=-1

        total_ed_matrix=[]          #List for storing EBE matrices of individual educt compounds
        total_prod_matrix=[]        #List for storing EBE matrices of individual product compounds

        ed_elements=[]              #List for storing elements of educts
        prod_elements=[]            #List for storing elements of products



        for i in ed:                        # Loop to process the educt string input given by the user
            s=i[1:len(i)-1]
            c=0
            
            for j in range(len(s)):

                if (j+1)<len(s) and s[j].isupper() and s[j+1].islower():
                    c+=1
                elif s[j].isupper():
                    c+=1
                    
            

            for j in range(len(s)):

                if (j+1)<len(s) and s[j].isupper() and s[j+1].islower():
                    ed_elements.append(s[j]+s[j+1])
                elif s[j].isupper():
                    ed_elements.append(s[j])


            
            ed_matrix=[[0 for q in range(c)] for z in range(c)]
            top=-1
            stack=[-1 for i in range(2*c)]
            j=0
            k=0
            cnt=''
            chg_cnt_p=0
            chg_cnt_n=0
            flag=0
            top=top+1
            stack[top]=j
            m=0
            cyc_hash=[]                     #List for cyclic compounds
            cyc_hash.append([])
            cyc_hash[0].append(0)
            cyc_bond=[]
            cyc_bond.append([])
            cyc_bond[0].append(0)


            while k<len(s):
                m=0

                if k+1<len(s) and s[k].isupper() and s[k+1].islower() and cnt=='':              #Identifying Chemical Symbol in educt string
                    cnt=s[k]+s[k+1]
                    k=k+1

                elif k+1<len(s) and s[k].isupper() and cnt=='':                                 #Identifying Chemical Symbol in educt string 
                    cnt=s[k]
                    k=k+1

                elif k+1<len(s) and s[k].isupper() and s[k+1].islower() and cnt!='':            #Identifying Chemical Symbol in educt string
                    j=j+1
                    top=top+1
                    stack[top]=j

                    k=k+1

                elif k+1<len(s) and s[k].isupper() and cnt!='':                                 #Identifying Chemical Symbol in educt string
                    j=j+1
                    top=top+1
                    stack[top]=j

                    k=k+1

                elif k+1<len(s) and s[k]=='.':                                                  #Identifying valency of an element
                     for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            m=m*10+int(s[l])
                        else:
                            break

                     t1=top
                     while stack[t1]==-100:
                         t1=t1-1
                     p=stack[t1]
                     if ed_matrix[p][p]>0:                         
                         ed_matrix[p][p]+=m
                         
                     elif ed_matrix[p][p]<0:                         
                         ed_matrix[p][p]-=m
                     else:
                         ed_matrix[p][p]=m
                         
                         
                     k=k+1

                elif k+1<len(s) and (s[k]=='+' or s[k]=='-') and s[k+1].isdigit():              #Identifying charge of an element
                    cov='1'
                    for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            m=m*10+int(s[l])
                        else:
                            break

                        if flag==0:
                             if s[k]=='+':
                                 t1=top
                                 while stack[t1]==-100:
                                     t1=t1-1
                                 p=stack[t1]
                                 ed_matrix[p][p]=m*(10)+ed_matrix[p][p]

                             elif s[k]=='-':
                                 t1=top
                                 while stack[t1]==-100:
                                    t1=t1-1
                                 p=stack[t1]
                                 ed_matrix[p][p]=m*(-10)-ed_matrix[p][p]

                    k=k+1

                elif k+1<len(s) and s[k]=='[' and s[k+1].isalpha() and flag==0:                 #Identifying anionic compound
                    cnt=''
                    chg_cnt_p=0
                    chg_cnt_n=0
                    cov='1'

                    if flag==0:
                        flag=1

                    for d in range(k+1,len(s)):
                        if d+1<len(s) and s[d].isupper() and s[d+1].islower():
                            chg_cnt_p+=1
                        elif d+1<len(s) and s[d].isupper():
                            chg_cnt_p+=1
                        elif s[d]==']':
                            break



                    for q in range(chg_cnt_p):
                        ed_matrix[q][q]=10
                        for w in range(chg_cnt_p,c):
                            ed_matrix[q][w]=10



                    for t in range(d+1,len(s)):
                        if t+1<len(s) and s[t].isupper() and s[t+1].islower():
                            chg_cnt_n+=1
                        elif t+1<len(s) and s[t].isupper():
                            chg_cnt_n+=1
                        elif s[t]==']':
                            break


                    for q in range(chg_cnt_p,c):
                        ed_matrix[q][q]=-10
                        for w in range(chg_cnt_p):
                            ed_matrix[q][w]=-10



                    k=k+1


                elif k+1<len(s) and s[k]==']' and flag==1 and s[k+1]=='[':                      #Identifying cationic compounds
                    m=0
                    flag=2
                    k=k+1
                    j=j+1
                    top=top+1
                    stack[top]=j
                    cnt=''

                elif k+1<len(s) and s[k]==']' and flag==2:                                      #Identifying cationic compounds
                    m=0
                    flag=0
                    k=k+1
                    j=j+1
                    top=top+1
                    stack[top]=j
                    cnt=''

                elif k+1<len(s) and s[k]=='[' and s[k+1].isdigit():                             #Identifying Cyclic compounds

                    for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            m=m*10+int(s[l])
                        else:
                            break

                    if len(cyc_hash)<=m:
                         cyc_hash.append([])
                         cyc_hash[m].append(j)
                         if s[l]=='-':
                             cyc_bond.append([])
                             cyc_bond[m].append(1)
                         elif s[l]=='=':
                             cyc_bond.append([])
                             cyc_bond[m].append(2)
                         elif s[l]=='#':
                             cyc_bond.append([])
                             cyc_bond[m].append(3)
                         elif l+1<len(s) and s[l]=='~' and s[l+1]=='/':
                             cyc_bond.append([])
                             cyc_bond[m].append(200)

                         elif s[l]=='~':
                             cyc_bond.append([])
                             cyc_bond[m].append(100)


                    else:
                        cyc_hash[m].append(j)

                    k=k+1
                    flag=-1

                elif k+1<len(s) and s[k]==']' and flag==-1:                                     #Identifying cyclic compounds
                    flag=0
                    k=k+1





                elif s[k]=='-' and flag!=-1:                                                    #Identifying Single bond
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=1
                    ed_matrix[j][p]=1
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                    if bnd=='':
                        bnd='1'


                elif s[k]=='=' and flag!=-1:                                                    #Identifying Double bond
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=2
                    ed_matrix[j][p]=2
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                    if bnd=='' or bnd=='1':
                        bnd='2'
                              
                elif s[k]=='#' and flag!=-1:                                                    #Identifying Triple bond
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=3
                    ed_matrix[j][p]=3
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                    if bnd=='' or bnd=='1' or bnd=='2':
                        bnd='3'
                              
                elif s[k]=='~' and s[k+1]=='/' and flag!=-1:                                    #Identifying Coordinate bond
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=200
                    ed_matrix[j][p]=100
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                    if bnd=='' or bnd=='1' or bnd=='2' or bnd=='3':
                        bnd='4'

                elif s[k]=='~' and flag!=-1:                                                    #Identifying Coordinate bond
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[p][j]=100
                    ed_matrix[j][p]=200
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                    if bnd=='' or bnd=='1' or bnd=='2' or bnd=='3':
                        bnd='4'
                        
                elif s[k]=='*' and flag!=-1:                                                    #Identifying Addition compound
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    ed_matrix[1][j]=500
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''
                    if bnd=='' or bnd=='1' or bnd=='2' or bnd=='3':
                        bnd='4'


                elif s[k]=='(':                                                                 #Identifying Branched Compound
                    if s[k+1].isalpha():
                        add_flag=1
                        add_top=top
                        while stack[add_top]==-100:
                            add_top=add_top-1 
                        edd_matrix[0][stack[add_top]+1]=500
                        top=top+1
                        j=j+1
                        stack[top]=j
                    
                    else:                       
                        top=top+1
                        stack[top]=-100
                        
                        if add_flag>=1:
                            add_flag+=1
                    k=k+1
                    cnt=''

                elif k+1<len(s) and s[k]==')' and s[k+1].isalpha():
                    k=k+1

                elif s[k]==')':                                                                 #Identifying Branched compound
                    if add_flag==1:
                        add_flag=0
                        
                    elif add_flag!=1:                    
                    
                        while stack[top]!=-100:
                            top=top-1
                        top=top-1

                        if add_flag>=1:
                            add_flag-=1

                        
                        
                        
                    k=k+1
                    cnt=''


                else:
                    k=k+1
                    
            if len(cyc_hash)>1:                                                                 #Processing for cyclic compound
                for i in range(1,len(cyc_hash)):
                    p1=cyc_hash[i][0]
                    p2=cyc_hash[i][1]
                    bond=cyc_bond[i][0]

                    ed_matrix[p1][p2]=bond
                    ed_matrix[p2][p1]=bond

                    if bond==100:
                        ed_matrix[p1][p2]=100
                        ed_matrix[p2][p1]=200
                        
                    elif bond==200:
                        ed_matrix[p1][p2]=200
                        ed_matrix[p2][p1]=100



            total_ed_matrix.append(ed_matrix)                                                   #Storing the educt matrix in a list



        if (prod.find(' + ') == -1):                                                            #Extracting product compounds
            prod = [prod]
        else:
            prod = prod.split(' + ')



        prod_compounds = prod
        s = ''
        for i in prod:
            s=i[1:len(i)-1]
            c=0
            for j in range(len(s)):

                if (j+1)<len(s) and s[j].isupper() and s[j+1].islower():
                    c+=1
                elif s[j].isupper():
                    c+=1

        for i in prod:                                                                          #Processing esmiles notation of each product compound
            s=i[1:len(i)-1]
            c=0
            
            for j in range(len(s)):

                if (j+1)<len(s) and s[j].isupper() and s[j+1].islower():        
                    c+=1
                elif s[j].isupper():
                    c+=1
                    

            for j in range(len(s)):

                if (j+1)<len(s) and s[j].isupper() and s[j+1].islower():
                    prod_elements.append(s[j]+s[j+1])
                elif s[j].isupper():
                    prod_elements.append(s[j])



            prod_matrix=[[0 for q in range(c)] for z in range(c)]
            top=-1
            stack=[-1 for i in range(2*c)]
            j=0
            k=0
            cnt=''
            chg_cnt_p=0
            chg_cnt_n=0
            flag=0
            top=top+1
            stack[top]=j
            m=0
            cyc_hash=[]
            cyc_hash.append([])
            cyc_hash[0].append(0)
            cyc_bond=[]
            cyc_bond.append([])
            cyc_bond[0].append(0)
            add_flag=0


            while k<len(s):
                m=0

                if k+1<len(s) and s[k].isupper() and s[k+1].islower() and cnt=='':                  #Identifying Chemical Symbol in product string
                    cnt=s[k]+s[k+1]
                    k=k+1

                elif k+1<len(s) and s[k].isupper() and cnt=='':                                     #Identifying Chemical Symbol in product string
                    cnt=s[k]
                    k=k+1

                elif k+1<len(s) and s[k].isupper() and s[k+1].islower() and cnt!='':                #Identifying Chemical Symbol in product string
                    j=j+1
                    top=top+1
                    stack[top]=j

                    k=k+1

                elif k+1<len(s) and s[k].isupper() and cnt!='':                                     #Identifying Chemical Symbol in product string
                    j=j+1
                    top=top+1
                    stack[top]=j

                    k=k+1

                elif k+1<len(s) and s[k]=='.':                                                      #Identifying valency of an element in product string
                     for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            m=m*10+int(s[l])
                        else:
                            break
                     t1=top
                     while stack[t1]==-100:
                        t1=t1-1
                     p=stack[t1]
                     if prod_matrix[p][p]>0:
                         prod_matrix[p][p]+=m
                     elif prod_matrix[p][p]<0:
                         prod_matrix[p][p]-=m
                     else:
                         prod_matrix[p][p]=m
                         

                     k=k+1

                elif k+1<len(s) and (s[k]=='+' or s[k]=='-') and s[k+1].isdigit():                  #Identifying charge in product string

                    for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            m=m*10+int(s[l])
                        else:
                            break

                        if flag==0:
                             if s[k]=='+':
                                 t1=top
                                 while stack[t1]==-100:
                                     t1=t1-1
                                 p=stack[t1]
                                 prod_matrix[p][p]=m*(10)+prod_matrix[p][p]

                             elif s[k]=='-':
                                 t1=top
                                 while stack[t1]==-100:
                                    t1=t1-1
                                 p=stack[t1]
                                 prod_matrix[p][p]=m*(-10)-prod_matrix[p][p]

                    k=k+1

                elif k+1<len(s) and s[k]=='[' and s[k+1].isalpha() and flag==0:                     #Identifying Ionic bond in product string
                    cnt=''
                    chg_cnt_p=0
                    chg_cnt_n=0

                    if flag==0:
                        flag=1

                    for d in range(k+1,len(s)):
                        if d+1<len(s) and s[d].isupper() and s[d+1].islower():
                            chg_cnt_p+=1
                        elif d+1<len(s) and s[d].isupper():
                            chg_cnt_p+=1
                        elif s[d]==']':
                            break



                    for q in range(chg_cnt_p):
                        prod_matrix[q][q]=10
                        for w in range(chg_cnt_p,c):
                            prod_matrix[q][w]=10



                    for t in range(d+1,len(s)):
                        if t+1<len(s) and s[t].isupper() and s[t+1].islower():
                            chg_cnt_n+=1
                        elif t+1<len(s) and s[t].isupper():
                            chg_cnt_n+=1
                        elif s[t]==']':
                            break


                    for q in range(chg_cnt_p,c):
                        prod_matrix[q][q]=-10
                        for w in range(chg_cnt_p):
                            prod_matrix[q][w]=-10



                    k=k+1


                elif k+1<len(s) and s[k]==']' and flag==1 and s[k+1]=='[':                  #Identifying Ionic bond in product string
                    m=0
                    flag=2
                    k=k+1
                    j=j+1
                    top=top+1
                    stack[top]=j
                    cnt=''

                elif k+1<len(s) and s[k]==']' and flag==2:                                  #Identifying Ionic bond in product string
                    m=0
                    flag=0
                    k=k+1
                    cnt=''
                   

                elif k+1<len(s) and s[k]=='[' and s[k+1].isdigit():                         #Identifying cyclic bond in product string


                    for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            m=m*10+int(s[l])
                        else:
                            break

                    if len(cyc_hash)<=m:
                         cyc_hash.append([])
                         cyc_hash[m].append(j)
                         if s[l]=='-':
                             cyc_bond.append([])
                             cyc_bond[m].append(1)
                         elif s[l]=='=':
                             cyc_bond.append([])
                             cyc_bond[m].append(2)
                         elif s[l]=='#':
                             cyc_bond.append([])
                             cyc_bond[m].append(3)

                         elif l+1<len(s) and s[l]=='~' and s[l+1]=='/':
                             cyc_bond.append([])
                             cyc_bond[m].append(100)
                        

                         elif s[l]=='~':
                             cyc_bond.append([])
                             cyc_bond[m].append(100)


                    else:
                        cyc_hash[m].append(j)

                    k=k+1
                    flag=-1

                elif k+1<len(s) and s[k]==']' and flag==-1:                             #Identifying cyclic bond in product string
                    flag=0
                    k=k+1

                elif s[k]=='-' and flag!=-1:                                            #Identifying single bond in product string
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    prod_matrix[p][j]=1
                    prod_matrix[j][p]=1
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''


                elif s[k]=='=' and flag!=-1:                                            #Identifying double bond in product string
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    prod_matrix[p][j]=2
                    prod_matrix[j][p]=2
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''


                elif s[k]=='#' and flag!=-1:                                        #Identifying triple bond in product string
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    prod_matrix[p][j]=3
                    prod_matrix[j][p]=3
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''


                elif s[k]=='~' and s[k+1]=='/' and flag!=-1:                        #Identifying coordinate bond in product string
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    prod_matrix[p][j]=200
                    prod_matrix[j][p]=100
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''


                elif s[k]=='~' and flag!=-1:                                        #Identifying coordinate bond in product string
                    j=j+1
                    t1=top
                    while stack[t1]==-100:
                        t1=t1-1
                    p=stack[t1]
                    prod_matrix[p][j]=100
                    prod_matrix[j][p]=200
                    top=top+1
                    stack[top]=j
                    k=k+1
                    cnt=''


                elif s[k]=='(':                                                         #Identifying branched compound
                    
                    if s[k+1].isalpha():
                        add_flag=1
                        add_top=top
                        while stack[add_top]==-100:
                            add_top=add_top-1 
                        prod_matrix[0][stack[add_top]+1]=500
                        top=top+1
                        j=j+1
                        stack[top]=j
                    
                    else:                       
                        top=top+1
                        stack[top]=-100
                        
                        if add_flag>=1:
                            add_flag+=1
                            
                        
                        
                
                    k=k+1
                    cnt=''
                    
                elif k+1<len(s) and s[k]==')' and s[k+1].isalpha():
                    k=k+1


                elif s[k]==')':                                                         #Identifying branched compound


                    if add_flag==1:
                        add_flag=0
                        
                    elif add_flag!=1:                    
                    
                        while stack[top]!=-100:
                            top=top-1
                        top=top-1

                        if add_flag>=1:
                            add_flag-=1

                        
                        
                        
                    k=k+1
                    cnt=''

                else:
                    k=k+1
                    

            if len(cyc_hash)>1:                                                     #processing cyclic compounds for product
                for i in range(1,len(cyc_hash)):
                    p1=cyc_hash[i][0]
                    p2=cyc_hash[i][1]
                    bond=cyc_bond[i][0]

                    prod_matrix[p1][p2]=bond
                    prod_matrix[p2][p1]=bond

                    if bond==100:
                        prod_matrix[p1][p2]=100
                        prod_matrix[p2][p1]=200
                        
                    elif bond==200:
                        prod_matrix[p1][p2]=200
                        prod_matrix[p2][p1]=100


            total_prod_matrix.append(prod_matrix)
            

        if f1 == 0 and f2 == 0:
            delete_reaction()                                               #Method call to delete reaction from database
            
            
            
            
window = Tk()                       #First window for the user
ed = ''
prod = ''
educt = StringVar()                 #User input for educts
product = StringVar()               #User input for products
ed_types = []
prod_types = []
ed_compounds = []
prod_compounds = []
educt_matrix = []
total_ed_matrix = []
total_prod_matrix = []
ed_elements = []
prod_elements = []
ed_element_slno = []
prod_element_slno = []
wrng_ed = 0
wrng_prd = 0
v = StringVar()
educt_element = []
product_element = []
filename="Reactions_"           #
cov=''                          #
parts=''                        #
bnd=''                          # to identify the filaname to store reaction and templates





# First Window Screen that appears for the User
w, h = window.winfo_screenwidth(), window.winfo_screenheight()
window.geometry("%dx%d+0+0" % (w, h))
window.title("DELETION OF REACTION")

window.configure(bg="white")


Label(window, text="Welcome", font=("Arial Bold", 50), fg = "black", bg = "white").pack()

Label(window, text="Please Enter the Chemical Reactions in the text box in the ESMILES Notation which you want to delete:", font=("Arial Bold", 15),fg = "dark green", bg = "white").place(x=40, y=140)       #
                                                                                    
Label(window, text="Please click on the button to read the Manual and learn how to input chemical compounds and elements in ESMILES notation ", font=("Arial Bold", 10),fg = "dark green", bg = "white").place(x=60, y=200)                       #

b = Button(window, text='Read Manual', command=read_manual,fg = "blue", bg = "white").place(x=900, y=190)           # Calling "go_press" function after button click




Label(window, text="Educt", fg = "black", bg = "white").place(x=40, y=290)


E1 = Entry(window, textvariable=educt,width=70,font=("Arial Bold", 15))         # Textbox for Reading Educt Compounds
E1.place(x=100, y=290)


Label(window, text="Product", fg = "black", bg = "white").place(x=40, y=390)                        

E2 = Entry(window, textvariable=product,width=70,font=("Arial Bold", 15))       # Textbox for Reading Product Compounds
E2.place(x=100, y=390)


btn = Button(window, text='Delete', command=go_press,fg = "dark green", bg = "white").place(x=500, y=460)        # Calling "Delete" function after button click



ele = pandas.read_csv("elementlist.csv", header=None)                           #Reading Elemenets name and Valency from file

elements = ele[1].values.tolist()                                               #List containing all Elements symbol in Periodic table
elements_name=ele[2].values.tolist()                                            #List containing all Elements name in Periodic table
elements_valency = ele[3].values.tolist()                                       #List containing all Elements valency in Periodic table
elements_free_electons = ele[4].values.tolist()
elements_stable_electons = ele[5].values.tolist()
elements_max_bond = ele[6].values.tolist()
elements_eneg = ele[7].values.tolist()

window.mainloop()
