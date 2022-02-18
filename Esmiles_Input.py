from tkinter import *                      # Importing tkinter for designing the UI
from tkinter.messagebox import *           # Importing tkinter.messagebox for warning/cofirmation message
import os                                  # Importing os to run system command
import pandas                              # Importing pandas to for datastructure
import csv                                 # Importing csv to read/write csv file
from itertools import permutations         # Importing permutation for permuting list
import subprocess
import webbrowser

def read_manual():
    webbrowser.open_new('Manual_ESMILES.pdf')

def run(c):                                # Function to Run New Program
    c.destroy()
    window.destroy()


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
                    mp=0
                    for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            mp=mp*10+int(s[l])
                        else:
                            break

                        if flag==0:
                             if s[k]=='+':
                                 t1=top
                                 while stack[t1]==-100:
                                     t1=t1-1
                                 p=stack[t1]
                                 ed_matrix[p][p]=mp*(10)+ed_matrix[p][p]

                             elif s[k]=='-':
                                 t1=top
                                 while stack[t1]==-100:
                                    t1=t1-1
                                 p=stack[t1]
                                 ed_matrix[p][p]=mp*(-10)-ed_matrix[p][p]

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
                        ed_matrix[p1][p2]=100


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
    #print(ebe_matrix_educts)
    #print(ed_elements)
        
    return 1
            
            
                
    

    
def database_entry(cv_row,temp_row,cov,bnd,parts):      # Method to store reactions and templates in database/files
     
          
     with open("Reactions_"+str(cov)+str(parts)+str(bnd)+".csv", 'a') as csvFile:   #storing reaction in file
        csvFile.write(str(cv_row)+"\n")

     csvFile.close()
     
    

     with open("Templates_"+str(cov)+str(parts)+str(bnd)+".csv", 'a') as csvFile:   #storing templates in file
         csvFile.write(str(temp_row[0])+"\n")

     csvFile.close()
     
     showinfo("Insertion Done", "Reaction is inserted into the Database properly")  #confirmation after storing reeactions
     
     f="Reactions_"+str(cov)+str(parts)+str(bnd)+".csv"         #File  name where reactions has to be stored
     t="Templates_"+str(cov)+str(parts)+str(bnd)+".csv"         #File name where templates has to be stored


     if os.path.isfile(f):                              #removing duplicate reactions from file
        reader = open(f, "r")           
        
        lines = reader.read().split("\n")               #list containing reactions from the file
        reader.close()
        writer = open(f, "w")
        for line in set(lines):
            writer.write(line + "\n")                   #Removing duplicates and writing in the file
        writer.close()
        
     if os.path.isfile(t):                              #removing duplicate templates from file
        reader = open(t, "r")           
        
        lines = reader.read().split("\n")               #list containing templates from the file
        reader.close()
        writer = open(t, "w")
        for line in set(lines):
            writer.write(line + "\n")                   #Removing duplicates and writing in the file
        writer.close()

def reference(template,alpha):                          #Method to show legends for alphabets and symbols
    global ed_elements

    root = Tk()                                         #Window to display templates and legends
    
    root.config(bg="white")

    w, h = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry("%dx%d+0+0" % (w, h))
    root.title("Template Reference")
    

    Label(root, text="Template\n\n", font=("Arial Bold", 20),bg="white").pack()

    edt=StringVar(root,value=template)

    M3=Entry(root, textvariable=edt,state='readonly',width=70,font=("Arial Bold", 20))      
    M4=Scrollbar(root, orient='horizontal', command=M3.xview)
    M3.config(xscrollcommand=M4.set)


    M3.pack()
    M4.pack()
    
    
    a=20
    b=0
    for i in range(len(alpha)):                     #Loop to map alphabets to chemical symbols
        name=''
        if alpha[i]!='':
            for j in range(len(elements)):
                if elements[j]==ed_elements[i]:     #identifying the chemical symbol
                    name=elements_name[j]           #name of the chemical symbol
                    break

            Label(root, text=alpha[i] +" ---> "+ed_elements[i]+" ("+name+") " , font=("Arial Bold", 15), bg="white").place(x=90+b,y=200+a)
            a=a+50

            if (i+1)%9==0:
                a=20
                b=b+250

    #list of legends
    choices = ['Legends', '-- :    Single Bond','= :   Double Bond','# :  Triple Bond','~> :  Coordinate Bond','[ ] :   Charged Compounds','+2 :  Positive Charge with Value 2','-2 :  Negative Charge with Value 2', '.2 :  2 Valence Electrons' ]
    tkvar = StringVar(root)
    tkvar.set('Legends') 

    popupMenu = OptionMenu(root, tkvar, *choices)           #Pop up menu for legends
    popupMenu.place(x=40,y=40)

    nbtn = Button(root, text='EXIT', command=root.destroy,fg = "dark green", bg = "white").place(x=650, y=560) 



def ESMILES_to_EBE():  #Function to process the Educts and Products entered in the textbox

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
            
            for j in range(len(s)):         #Identifying chemical symbols

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
                    mp=0
                    for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            mp=mp*10+int(s[l])
                        else:
                            break

                        if flag==0:
                             if s[k]=='+':
                                 t1=top
                                 while stack[t1]==-100:
                                     t1=t1-1
                                 p=stack[t1]
                                 ed_matrix[p][p]=mp*(10)+ed_matrix[p][p]

                             elif s[k]=='-':
                                 t1=top
                                 while stack[t1]==-100:
                                    t1=t1-1
                                 p=stack[t1]
                                 ed_matrix[p][p]=mp*(-10)-ed_matrix[p][p]

                    k=k+1

                elif k+1<len(s) and s[k]=='[' and s[k+1].isalpha() and flag==0:                 #Identifying ionic compound
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
                        
                


                elif s[k]=='(':                                                         #Identifying branched compound
                    
                    if s[k+1].isalpha():
                        add_flag=1
                        add_top=top
                        while stack[add_top]==-100:
                            add_top=add_top-1 
                        ed_matrix[0][stack[add_top]+1]=500
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
            add_flag=0
            top=top+1
            stack[top]=j
            m=0
            cyc_hash=[]
            cyc_hash.append([])
            cyc_hash[0].append(0)
            cyc_bond=[]
            cyc_bond.append([])
            cyc_bond[0].append(0)



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
                    mp=0
                    
                    for l in range(k+1,len(s)):
                        if s[l].isdigit():
                            mp=mp*10+int(s[l])
                        else:
                            break

                        if flag==0:
                             if s[k]=='+':
                                 t1=top
                                 while stack[t1]==-100:
                                     t1=t1-1
                                 p=stack[t1]
                                 prod_matrix[p][p]=mp*(10)+prod_matrix[p][p]
                                 

                             elif s[k]=='-':
                                 t1=top
                                 while stack[t1]==-100:
                                    t1=t1-1
                                 p=stack[t1]
                                 prod_matrix[p][p]=mp*(-10)-prod_matrix[p][p]
                                 #print(p)

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
                             cyc_bond[m].append(200)          

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
            window.withdraw()
            build_template()                                        #Method call to get serial number of product elements





def setb(root):         #Method to store the serial numbers of elements in educts

    global prod_element_slno, prod_elements, ed_element_slno

    for i in range(len(prod_element_slno)):
        prod_element_slno[i] = int(prod_element_slno[i].get())

    print_template()  # Function Call To generate Reaction Matrix 

    root.destroy()


def build_template():  #Method to get the serial numbers of elements in products

    global total_ed_matrix, total_prd_matrix, ed_elements, prod_elements, ed_element_slno, prod_element_slno, ed_compounds, prod_compounds

    ed_string = ''
    prod_string= ''

    for i in range(len(ed_compounds)):
        if i==len(ed_compounds)-1:
            ed_string+=ed_compounds[i][1:len(ed_compounds[i])-1]
        else:
            ed_string+=ed_compounds[i][1:len(ed_compounds[i])-1] + " + "


    for i in range(len(prod_compounds)):
        if i==len(prod_compounds)-1:
            prod_string+=prod_compounds[i][1:len(prod_compounds[i])-1]
        else:
            prod_string+=prod_compounds[i][1:len(prod_compounds[i])-1] + " + "

    root = Tk()                                 #New window to enter serial nuber of product elemens by the user
    
    root.configure(bg="white")
    
    w, h = window.winfo_screenwidth(), window.winfo_screenheight()
    root.geometry("%dx%d+0+0" % (w, h))
    edt1=StringVar(root,value="Educts -> "+ed_string)
    edt2=StringVar(root,value="Products -> "+prod_string)

    M1=Entry(root, textvariable=edt1,state='readonly',width=40,font=("Arial Bold", 15))
    M2=Scrollbar(root, orient='horizontal', command=M1.xview)
    M1.config(xscrollcommand=M2.set)


    M1.place(x=45, y=40)
    M2.place(x=280, y=72)

    M3=Entry(root, textvariable=edt2,state='readonly',width=40,font=("Arial Bold", 15))         #Entry box to enter the product serial number
    M4=Scrollbar(root, orient='horizontal', command=M3.xview)
    M3.config(xscrollcommand=M4.set)


    M3.place(x=805, y=40)
    M4.place(x=1055, y=72)



    #Label(root, text="Educts -> "+ed_string, font=("Arial Bold", 15)).place(x=45, y=40)
    #Label(root, text="Products -> "+prod_string, font=("Arial Bold", 15)).place(x=805, y=40)

    m = 0
    n = 0
    c = 0
    x = 0
    f = 0
    for i in ed_elements:
        c = c+1
        if c % 17 == 0:
            m = m+150
            x = n
            f = 1
            n = 0
        Label(root, text=i+"  --> ", bg="white").place(x=45+m, y=100+n)
        Label(root, text=c, bg="white").place(x=90+m, y=100+n)
        ed_element_slno.append(c)
        n = n+30
        if f == 0:
            x = n

    n = 0
    m = 0
    c = 0
    for i in prod_elements:
        c = c+1
        if c % 17 == 0:
            m = m+180
            n = 0
        Label(root, text=i+"  --> ", bg="white").place(x=805+m, y=100+n)
        v = Entry(root, width=7)
        v.place(x=900+m, y=100+n)
        prod_element_slno.append(v)
        n = n+30

    B = Button(root, text='Set', command=lambda: setb(root), fg ="blue", bg="white")  #calling method to stre the serial number 
    B.place(x=805, y=100+x)

    Label(root, text="Enter Corresponding Sl. No. to the Product which Corresponds to same element in Educt Side", fg ="green", bg="white").place(x=45, y=110+x)

    

def print_template():  # To generate Ebe Matrix of Educts,Products, Reaction Matrix and Template
    
    global prod_element_slno, ed_element_slno, total_ed_matrix, total_prod_matrix, ed_compounds, prod_compounds,filename,cov,bnd,parts

    ed_len = len(ed_element_slno)           
    prod_len = len(prod_element_slno)

    reaction = ''
    
    for i in range(len(ed_compounds)):
        if i==len(ed_compounds)-1:
            reaction+=ed_compounds[i][1:len(ed_compounds[i])-1] + " ----> "
            
        else:
            reaction+=ed_compounds[i][1:len(ed_compounds[i])-1] + " + "
            


    for i in range(len(prod_compounds)):
        if i==len(prod_compounds)-1:
            reaction+=prod_compounds[i][1:len(prod_compounds[i])-1]
            
        else:
            reaction+=prod_compounds[i][1:len(prod_compounds[i])-1] + " + "
            


    #print(total_ed_matrix)
    #print(total_prod_matrix)


    
    ebe_matrix_educts = [[0 for i in range(ed_len+2)] for j in range(ed_len)]           # Ebe Matrix of Educts
    
    
    ebe_matrix_products = [[0 for i in range(prod_len+2)] for j in range(prod_len)]     # Ebe Matrix of Products

    c = 0
    for i in total_ed_matrix:                   #Loop to fill the EBE matrix of educts
        for j in range(len(i)):
            for k in range(len(i)):
                x = ed_element_slno[c+j]-1
                y = ed_element_slno[c+k]-1
                ebe_matrix_educts[x][y] = i[j][k]

                if ebe_matrix_educts[x][y] == 100:          #For Coordinate bond
                    #print(ebe_matrix_educts[x][ed_len+1])
                    if ebe_matrix_educts[x][ed_len+1]!=0:
                        ebe_matrix_educts[x][ed_len+1]=-99
                    else:
                        ebe_matrix_educts[x][ed_len+1] = y+1

                    if ebe_matrix_educts[y][ed_len+1]!=0:
                        ebe_matrix_educts[y][ed_len+1]=-99
                    else:
                        ebe_matrix_educts[y][ed_len+1] = x+1


                if ebe_matrix_educts[x][y] >= 10 and ebe_matrix_educts[x][y] <= 88:         #For Cations
                    ebe_matrix_educts[x][ed_len] = ebe_matrix_educts[x][y]//10
                    if x != y:
                        ebe_matrix_educts[x][y] = 1000
                        ebe_matrix_educts[y][x] = 2000
                    else:
                        ebe_matrix_educts[x][x] = ebe_matrix_educts[x][x] % 10

                if ebe_matrix_educts[x][y] >= -88 and ebe_matrix_educts[x][y] <= -10:       #For Anions
                    ebe_matrix_educts[x][ed_len] = ebe_matrix_educts[x][y]//10
                    if x != y:
                        ebe_matrix_educts[x][y] = 2000
                        ebe_matrix_educts[y][x] = 1000
                    else:
                        ebe_matrix_educts[x][x] = (ebe_matrix_educts[x][x]*-1) % 10

        c = c+len(i)

    c = 0
    for i in total_prod_matrix:                         #Loop to fill the EBE matrix of products
        for j in range(len(i)):
            for k in range(len(i)):
                x = prod_element_slno[c+j]-1
                y = prod_element_slno[c+k]-1
                ebe_matrix_products[x][y] = i[j][k]
                
                if ebe_matrix_products[x][y] == 100:                #For Coordinate bond
                    #print(ebe_matrix_products[x][ed_len+1])
                    if ebe_matrix_products[x][ed_len+1]!=0:
                        ebe_matrix_products[x][ed_len+1]=-99
                    else:
                        ebe_matrix_products[x][ed_len+1] = y+1

                    if ebe_matrix_products[y][ed_len+1]!=0:
                        ebe_matrix_products[y][ed_len+1]=-99
                    else:
                        ebe_matrix_products[y][ed_len+1] = x+1

                if ebe_matrix_products[x][y] >= 10 and ebe_matrix_products[x][y] <= 88:         #For Cations
                    ebe_matrix_products[x][ed_len] = ebe_matrix_products[x][y]//10
                    if x != y:
                        ebe_matrix_products[x][y] = 1000
                        ebe_matrix_products[y][x] = 2000
                    else:
                        ebe_matrix_products[x][x] = ebe_matrix_products[x][x] % 10

                if ebe_matrix_products[x][y] >= -88 and ebe_matrix_products[x][y] <= -10:           #For Anions
                    ebe_matrix_products[x][ed_len] = ebe_matrix_products[x][y]//10
                    if x != y:
                        ebe_matrix_products[x][y] = 2000
                        ebe_matrix_products[y][x] = 1000
                    else:
                        ebe_matrix_products[x][x] = (ebe_matrix_products[x][x]*-1) % 10

        c = c+len(i)



    sorted_educt_matrix=ebe_matrix_educts[:]        
    sorted_product_matrix=ebe_matrix_products[:]


    for i in range(ed_len-1):              #Sorting the ebe matrix of educts and prducts based on chemical symbol of elements
        for j in range(ed_len-i-1):
                ch1=ed_elements[j][0]
                ch2=ed_elements[j+1][0]

                if(ch1>ch2):
                    temp=ed_elements[j]
                    ed_elements[j]=ed_elements[j+1]
                    ed_elements[j+1]=temp

                    sorted_educt_matrix[j], sorted_educt_matrix[j+1]=sorted_educt_matrix[j+1], sorted_educt_matrix[j]
                    sorted_product_matrix[j], sorted_product_matrix[j+1]=sorted_product_matrix[j+1], sorted_product_matrix[j]

                    temp_list=[sorted_educt_matrix[k][j] for k in range(ed_len)]

                    for k in range(ed_len):
                        sorted_educt_matrix[k][j]=sorted_educt_matrix[k][j+1]

                    for k in range(ed_len):
                        sorted_educt_matrix[k][j+1]=temp_list[k]

                    temp_list=[sorted_product_matrix[k][j] for k in range(ed_len)]

                    for k in range(ed_len):
                        sorted_product_matrix[k][j]=sorted_product_matrix[k][j+1]

                    for k in range(ed_len):
                        sorted_product_matrix[k][j+1]=temp_list[k]



    prod_elements=ed_elements[:]
    csv.register_dialect('myDialect',
                         delimiter = '|',                                      
                         skipinitialspace=True)

    ebe_matrix_educts=sorted_educt_matrix[:]
    ebe_matrix_products=sorted_product_matrix[:]
    
    for x in range(ed_len):
        ebe_matrix_educts[x][ed_len+1]=0
        ebe_matrix_products[x][ed_len+1]=0
            
    
    for x in range(ed_len):                                     #Filling up the coordinate bond column correctly in educt matrix
        for y in range(ed_len):
            if ebe_matrix_educts[x][y] == 100:
                    #print(ebe_matrix_educts[x][ed_len+1])
                    if ebe_matrix_educts[x][ed_len+1]!=0:
                        ebe_matrix_educts[x][ed_len+1]=-99
                    else:
                        ebe_matrix_educts[x][ed_len+1] = y+1

                    if ebe_matrix_educts[y][ed_len+1]!=0:
                        ebe_matrix_educts[y][ed_len+1]=-99
                    else:
                        ebe_matrix_educts[y][ed_len+1] = x+1
            
    for x in range(ed_len):                                     #Filling up the coordinate bond column correctly in product matrix
        for y in range(ed_len):
            if ebe_matrix_products[x][y] == 100:
                    #print(ebe_matrix_products[x][ed_len+1])
                    if ebe_matrix_products[x][ed_len+1]!=0:
                        ebe_matrix_products[x][ed_len+1]=-99
                    else:
                        ebe_matrix_products[x][ed_len+1] = y+1

                    if ebe_matrix_products[y][ed_len+1]!=0:
                        ebe_matrix_products[y][ed_len+1]=-99
                    else:
                        ebe_matrix_products[y][ed_len+1] = x+1
        
    #print(ed_elements)
    #print(ebe_matrix_educts)
    #print(ebe_matrix_products)
    
    flx= [[0 for i in range(ed_len+2)] for j in range(ed_len)]
    for i in range((ed_len)):
        flx[i]=ebe_matrix_educts[i][:]
        
    cv_row=[educt.get(),ed_elements,flx,product.get()]          #List to stre in Reactions file

    
    '''with open("Reactions_"+str(cov)+str(parts)+str(bnd)+".csv", 'a') as csvFile:
        csvFile.write(str(cv_row)+"\n")

    csvFile.close()'''





    reaction_matrix = [[0 for i in range(ed_len)] for j in range(ed_len)]  # Reaction Matrix

    break_bond = [[] for i in range(ed_len)]  # List to store elements whose bond has broken
    make_bond = [[] for i in range(ed_len)]  # List to store elements whose new bond has been made

    for i in range(ed_len):                 #Loop to fill up reaction matrix
        for j in range(ed_len):
            a = ebe_matrix_educts[i][j]
            b = ebe_matrix_products[i][j]
            '''if a == 100 or a == 200 or a == 500:        #For coordinate bond and additiion compound
                a = 1
            if b == 100 or b == 200 or b == 500:
                b = 1

            if a == 1000 or a == 2000:                  #For Ionic compound
                if b == 1000 or b == 2000:
                    a = 1
                    b = 1
                elif b != 0:
                    a = 0
                    b = 1
                elif b == 0:
                    if ebe_matrix_products[i][ed_len]==0 and ebe_matrix_products[j][ed_len]==0:
                        a = 0
                    else:
                        a = 1
                    a = 1
                        
                    

            if b == 1000 or b == 2000:
                if a == 1000 or a == 2000:
                    a = 1
                    b = 1
                elif a != 0:
                    b = 0
                    a = 1
                elif a == 0:
                    if ebe_matrix_educts[i][ed_len]==0 and ebe_matrix_educts[j][ed_len]==0:
                        b = 0
                    else:
                        b = 1
                    b = 1
            '''
            reaction_matrix[i][j] = b-a

    

    #print(ebe_matrix_educts)
    #print("\n")
    #print(ebe_matrix_products)
    #print("\n")
    #print("\n")
    #print(reaction_matrix)
    #print("\n")

    connected_components_educts=[[] for i in range(ed_len)]

    for i in range(ed_len):
        for j in range(ed_len):
            if ebe_matrix_educts[i][j]!=0 and i!=j:
                
                for k in range(len(connected_components_educts[i])):
                    connected_components_educts[j].append(connected_components_educts[i][k])
                for k in range(len(connected_components_educts[j])):
                    connected_components_educts[i].append(connected_components_educts[j][k])
                connected_components_educts[j].append(i)
                connected_components_educts[i].append(j)

                connected_components_educts[i]=list(set(connected_components_educts[i]))
                connected_components_educts[j]=list(set(connected_components_educts[j]))

    connected_components_products=[[] for i in range(ed_len)]

    for i in range(ed_len):
        for j in range(ed_len):
            if ebe_matrix_products[i][j]!=0 and i!=j:
                
                for k in range(len(connected_components_products[i])):
                    connected_components_products[j].append(connected_components_products[i][k])
                for k in range(len(connected_components_products[j])):
                    connected_components_products[i].append(connected_components_products[j][k])
                connected_components_products[j].append(i)
                connected_components_products[i].append(j)

                connected_components_products[i]=list(set(connected_components_products[i]))
                connected_components_products[j]=list(set(connected_components_products[j]))
                
                


    elli_list=[]                        #list to store elements having zero rows and columns in reaction matrix
    for i in range(ed_len):
        for j in range(ed_len):
            flage=0
            if reaction_matrix[i][j]!=0:
                flage=1
                break

        if flage==0:
            elli_list.append(i)
                
                            
                        
    #print(elli_list)

    elli_list=list(set(elli_list))
    
    for i in elli_list:                             #removing the elements from reaction matrix

        for k in range(ed_len):
            reaction_matrix[i][k]=0
            reaction_matrix[k][i]=0


    #print(reaction_matrix)
    #print("\n")


    elli_list=[]                    #list to store elements having non zero row or columns in reaction matrix but does not take part in reaction
    for i in range(ed_len):
        for j in range(ed_len):
            if i!=j:
                flage=0
                for k in range(ed_len):
                    if reaction_matrix[i][k]!=0:
                        if reaction_matrix[j][k]==0:
                            flage=1
                            break            
                                    

                if flage==0:
                    if reaction_matrix[i][j]==0:
                        if j not in elli_list:
                            if i in connected_components_educts[j]:
                                if i in connected_components_products[j]:
                                    elli_list.append(i)
                                    break
                


    elli_list=list(set(elli_list))
    
    for i in elli_list:                 #removing the elements from reaction matrix
        for k in range(ed_len):
            reaction_matrix[i][k]=0
            reaction_matrix[k][i]=0
    
    '''
    elli_list=[]                        #list to store elements having zero rows and columns in reaction matrix
    for i in range(ed_len):
        for j in range(i+1,ed_len):
            flage=0
            if (ebe_matrix_educts[i][ed_len]==ebe_matrix_educts[j][ed_len] and ebe_matrix_educts[i][ed_len]!=0) or (ebe_matrix_products[i][ed_len]==ebe_matrix_products[j][ed_len] and ebe_matrix_products[i][ed_len]!=0):
                for k in range(ed_len):
                    if reaction_matrix[i][k]!=reaction_matrix[j][k]:
                        flage=1
                        break

                if flage==0:
                    if ebe_matrix_educts[i][j]!=0 or ebe_matrix_products[i][j]!=0:
                        if ebe_matrix_educts[i][ed_len]!=0 and ebe_matrix_products[i][ed_len]!=0 and (ebe_matrix_educts[j][ed_len]==0 or ebe_matrix_products[j][ed_len]==0):
                            elli_list.append(j)
                        elif ebe_matrix_educts[j][ed_len]!=0 and ebe_matrix_products[j][ed_len]!=0 and (ebe_matrix_educts[i][ed_len]==0 or ebe_matrix_products[i][ed_len]==0):
                            elli_list.append(i)
                        
                            
                        
    #print(elli_list)

    elli_list=list(set(elli_list))
    
    for i in elli_list:                             #removing the elements from reaction matrix

        for k in range(ed_len):
            reaction_matrix[i][k]=0
            reaction_matrix[k][i]=0


    #print(reaction_matrix)
    #print("\n")


    elli_list=[]                    #list to store elements having non zero row or columns in reaction matrix but does not take part in reaction
    for i in range(ed_len):
        for j in range(ed_len):
            flage=0
            if i!=j:
                for k in range(ed_len):
                    if reaction_matrix[j][k]!=0:
                        if reaction_matrix[j][k]!=reaction_matrix[i][k]:
                            flage=1
                            break            
                                    

                if flage==0 and ebe_matrix_educts[i][j]!=0:
                    if i not in elli_list:
                        elli_list.append(j)
            


    elli_list=list(set(elli_list))
    
    for i in elli_list:                 #removing the elements from reaction matrix
        for k in range(ed_len):
            reaction_matrix[i][k]=0
            reaction_matrix[k][i]=0
   
    
    '''
    '''
    flagZ=0            
    for i in range(ed_len):
        for j in range(ed_len):
            if reaction_matrix[i][j]!=0:
                flagZ=1
                break
        if flagZ==1:
            break
    

    if flagZ==0:
         for i in range(ed_len):
             if ebe_matrix_educts[i][ed_len] != 0:
                 break_bond[i].append(i)
                 
             if ebe_matrix_products[i][ed_len] != 0:
                 make_bond[i].append(i)
             
             for j in range(i, ed_len):
                 if ebe_matrix_educts[i][j] == 1000 or ebe_matrix_educts[i][j] == 2000:
                     break_bond[i].append(j)
                 if ebe_matrix_products[i][j] == 1000 or ebe_matrix_products[i][j] == 2000:
                     make_bond[i].append(j)'''
                 
                 
        

   
    #print(reaction_matrix)
    #print("\n")        
          
    for i in range(ed_len):             #Loop to fillup break_bond and make_bond list
        for j in range(ed_len):
            flag = 0

            if reaction_matrix[i][j] != 0:
                if ebe_matrix_educts[i][j] != 0:
                    break_bond[i].append(j)
                    if i!=j:
                        break_bond[j].append(i)

                if ebe_matrix_products[i][j] != 0:
                    make_bond[i].append(j)
                    if i!=j:
                        make_bond[j].append(i)

                if ebe_matrix_educts[i][ed_len] != 0:
                    break_bond[i].append(i)

                if ebe_matrix_educts[j][ed_len] != 0:
                    break_bond[j].append(j)

                if ebe_matrix_products[i][ed_len] != 0:
                    make_bond[i].append(i)

                if ebe_matrix_products[j][ed_len] != 0:
                    make_bond[j].append(j)

   
    
    char = 65
    f = 0

    alpha = ['' for i in range(ed_len)]         #List to store alphabets A to Z
    template = ''
    flag = 0

    for i in range(len(break_bond)):                        #Sorting break_bond list and removing duplicates
        res = list(set(break_bond[i]))
        break_bond[i] = res[:]
        break_bond[i].sort()
        if len( break_bond[i]) > 1:
            for k in range(len(break_bond[i])):
                if break_bond[i][k] == i:
                    r = break_bond[i][:k]+break_bond[i][k+1:]
                    r.append(i)
                    break_bond[i]=r[:]
                    break
                             
                             

    for i in range(len(make_bond)):                        #Sorting break_bond list and removing duplicates
        res = list(set(make_bond[i]))
        make_bond[i] = res[:]
        make_bond[i].sort()
        if len(make_bond[i]) > 1:
            for k in range(len(make_bond[i])):
                if make_bond[i][k] == i:
                    r = make_bond[i][:k]+make_bond[i][k+1:]
                    r.append(i)
                    make_bond[i]=r[:]
                    break

   
    #print(break_bond)
    #print("\n")
    #print(make_bond)
    #print("\n")
    
     
    
    
    '''for i in range(ed_len):
        cnt_ion=0
        if len(break_bond[i])!=0:
            for j in range(len(break_bond[i])):
                if break_bond[i][j]==1000 or break_bond[i][j]==2000:
                    cnt_ion+=1
            if len(break_bond[i])>cnt_ion+1:
                brack_flag[i]=1'''
                
    brack_flag=[0 for i in range(ed_len)]         #List to store elements having branches
                        
    visited=[[0 for i in range(ed_len)] for j in range(ed_len)]         #visited matrix for dfs algorithm
    stack_atom=[]
    st_top=-1       
    ed_temp=''                      #Variable to store educt template
    
    template=''                     #Variable to store current template
    
    for i in range(ed_len):         #Loop to generate template for educts
        if len(break_bond[i])!=0:    
            st_top=st_top+1
            stack_atom.append(i)
            k=stack_atom[st_top]
            e=0
            while True:
                while e<len(break_bond[k]):
                    if visited[k][break_bond[k][e]]==0:
                        ign=0
                        w=break_bond[k][e]
                        if template.find('[')!=-1 and template.find(']')!=-1 and (ebe_matrix_educts[k][w]==1000 or ebe_matrix_educts[k][w]==2000) and k!=w:
                            for z in break_bond[w]:
                                if visited[w][z]==0 and ebe_matrix_educts[w][z]!=0 and ebe_matrix_educts[w][z]!=1000 and ebe_matrix_educts[w][z]!=2000:
                                    ign=1
                                    break
                        if ign==1:
                            e=e+1
                            continue
                        
                        if k!=break_bond[k][e]:
                            st_top+=1
                            stack_atom.append(break_bond[k][e])
                            
                        j=stack_atom[st_top]
                        
                        if alpha[k] == '':
                            alpha[k] = chr(char)
                            char = char+1
                            
                        cnt_ion=0
                        for x in break_bond[k]:
                           if x!= k and visited[x][k]==0 and (ebe_matrix_educts[k][x]!=1000 or ebe_matrix_educts[k][x]!=2000): 
                               cnt_ion+=1
                        if cnt_ion>1:
                           brack_flag[k]=1  
                                  
                        if template=='' or template.find(alpha[k])==-1:     #For unvisited elements in educts
                            if alpha[j] == '':
                                alpha[j] = chr(char)
                                char = char+1
                            
                            if brack_flag[k]==1:            #For branched compounds
                                
                                 if k == j:                 #For valency and charge
                                    if ebe_matrix_educts[k][k] != 0 and ebe_matrix_educts[k][k] != 1000 and ebe_matrix_educts[k][k] != 2000:
                                        template += alpha[k] + "(."+str(ebe_matrix_educts[k][k])+")"
                                    elif ebe_matrix_educts[k][k] == 0:
                                        f1 = 0
                                        for t in range(ed_len):
                                            if ebe_matrix_educts[k][t] == 1000 or ebe_matrix_educts[k][t] == 2000:
                                                f1 = 1
                                                break
                                        if ebe_matrix_educts[k][ed_len] > 0 and f1 == 0:
                                            template = template + alpha[k] + "(+"+str(ebe_matrix_educts[k][ed_len])+")"
                                        elif ebe_matrix_educts[k][ed_len] < 0 and f1 == 0:
                                            template = template + alpha[k] + "("+str(ebe_matrix_educts[k][ed_len])+")"                       
                                                                                    
                                                   
                                 elif ebe_matrix_educts[k][j] == 1:                         #For Single bond
                                    template = template +alpha[k]+"(--"+alpha[j]+")"
                                 elif ebe_matrix_educts[k][j] == 2:                         #For Double bond
                                    template =template +alpha[k]+"(="+alpha[j]+")"
                                 elif ebe_matrix_educts[k][j] == 3:                         #For triple bond
                                    template = template+alpha[k]+"(#"+alpha[j]+")"
                                 elif ebe_matrix_educts[k][j] == 100:                       #For coordinate bond
                                    template =  template+alpha[k]+"(~>"+alpha[j]+")"
                                 elif ebe_matrix_educts[k][j] == 200:                       #For coordinate bond
                                    template =  template+alpha[k]+"(<~"+alpha[j]+")"
                                 elif ebe_matrix_educts[k][j] == 500:                       #For addition compound
                                    fladdb=0
                                    for addb in range(len(template)-1):
                                        if template[addb]=='.' and template[addb+1]=='(':
                                            template = template+alpha[k]+"("+alpha[j]+")"
                                            fladdb=1
                                            break
                                    if fladdb==0:
                                        template = template+alpha[k]+".("+alpha[j]+")"                             
                                            
                                        
                                 elif ebe_matrix_educts[k][j] == 1000:                      #For cations
                                    template=template+'['+alpha[k]+']'+'['+alpha[j]+']'
                                 elif ebe_matrix_educts[k][j] == 2000:                      #For anions
                                    template=template+'['+alpha[j]+']'+'['+alpha[k]+']'
                            
                            elif brack_flag[k]==0:      #For unbranched compunds
                                
                                 if k == j:         #For valency and charge
                                    if ebe_matrix_educts[k][k] != 0 and ebe_matrix_educts[k][k] != 1000 and ebe_matrix_educts[k][k] != 2000:
                                        template += alpha[k] + "(."+str(ebe_matrix_educts[k][k])+")"
                                    elif ebe_matrix_educts[k][k] == 0:
                                        f1 = 0
                                        for t in range(ed_len):
                                            if ebe_matrix_educts[k][t] == 1000 or ebe_matrix_educts[k][t] == 2000:
                                                f1 = 1
                                                break
                                        if ebe_matrix_educts[k][ed_len] > 0 and f1 == 0:
                                            template+= alpha[k] + "(+"+str(ebe_matrix_educts[k][ed_len])+")"
                                        elif ebe_matrix_educts[k][ed_len] < 0 and f1 == 0:
                                            template+= alpha[k] + "("+str(ebe_matrix_educts[k][ed_len])+")"
                                                                                                            
                                            
                                 elif ebe_matrix_educts[k][j] == 1:                         #For single bond
                                    template = template +alpha[k]+"--"+alpha[j]
                                 elif ebe_matrix_educts[k][j] == 2:                         #For double bond
                                    template =template +alpha[k]+"="+alpha[j]
                                 elif ebe_matrix_educts[k][j] == 3:                         #For triple bond
                                    template = template+alpha[k]+"#"+alpha[j]
                                 elif ebe_matrix_educts[k][j] == 100:                       #For coordinate bond
                                    template =  template+alpha[k]+"~>"+alpha[j]
                                 elif ebe_matrix_educts[k][j] == 200:                       #For coordinate bond
                                    template =  template+alpha[k]+"<~"+alpha[j]
                                 elif ebe_matrix_educts[k][j] == 500:                       #For addition compound
                                    fladdb=0
                                    for addb in range(len(template)-1):
                                        if template[addb]=='.' and template[addb+1]=='(':
                                            template = template+alpha[k]+"("+alpha[j]+")"
                                            fladdb=1
                                            break
                                    if fladdb==0:
                                        template = template+alpha[k]+".("+alpha[j]+")"     
                                 elif ebe_matrix_educts[k][j] == 1000:                      #For cations
                                    template=template+'['+alpha[k]+']'+'['+alpha[j]+']'
                                 elif ebe_matrix_educts[k][j] == 2000:                      #For anions
                                    template=template+'['+alpha[j]+']'+'['+alpha[k]+']'    
                        
                       
                        else:                                           #For visited elements in educts
                            l = len(template)
                            for t in range(l):
                                if template[t]==alpha[k]:
                                     s1 = template[:t]
                                     s2 = template[t+1:]
                                     s3 = template[t]
                                     
                                     if alpha[j] == '':
                                        alpha[j] = chr(char)
                                        char = char+1
                                         
                                     if brack_flag[k]==0:               #For unbranched compound
                                         
                                         if k == j:                     #For valency and charge
                                             
                                            if ebe_matrix_educts[k][k] != 0 and ebe_matrix_educts[k][k] != 1000 and ebe_matrix_educts[k][k] != 2000:
                                                s3 = alpha[k] + "(."+str(ebe_matrix_educts[k][k])+")"
                                            elif ebe_matrix_educts[k][k] == 0:
                                                f1 = 0
                                                for t in range(ed_len):
                                                    if ebe_matrix_educts[k][t] == 1000 or ebe_matrix_educts[k][t] == 2000:
                                                        f1 = 1
                                                        break
                                                if ebe_matrix_educts[k][ed_len] > 0 and f1 == 0:
                                                    s3 = alpha[k] + "(+"+str(ebe_matrix_educts[k][ed_len])+")"
                                                elif ebe_matrix_educts[k][ed_len] < 0 and f1 == 0:
                                                    s3 = alpha[k] + "("+str(ebe_matrix_educts[k][ed_len])+")"
                                                
                                                                                                                                     
                
                                         elif ebe_matrix_educts[k][j] == 1:         #For single bond
                                            s3 = s3+"--"+alpha[j]
                                         elif ebe_matrix_educts[k][j] == 2:         #For double bond
                                            s3 = s3+"="+alpha[j]
                                         elif ebe_matrix_educts[k][j] == 3:         #For tiple bond
                                            s3 = s3+"#"+alpha[j]
                                         elif ebe_matrix_educts[k][j] == 100:       #For coordinate bond
                                            s3 = s3+"~>"+alpha[j]
                                         elif ebe_matrix_educts[k][j] == 200:       #For coordinate bond
                                            s3 = s3+"<~"+alpha[j]
                                         elif ebe_matrix_educts[k][j] == 500:                       #For addition compound
                                            fladdb=0
                                            for addb in range(len(s3)-1):
                                                if s3[addb]=='.' and s3[addb+1]=='(':
                                                    s3 = s3+alpha[k]+"("+alpha[j]+")"
                                                    fladdb=1
                                                    break
                                            if fladdb==0:
                                                s3 = s3+alpha[k]+".("+alpha[j]+")"     
                                         elif ebe_matrix_educts[k][j] == 1000:      #For cations
                                            si=t
                                            ei=t
                                            temp=''
                                            while(True):
                                                if si==0:
                                                    break
                                                elif template[si]==' ':
                                                    si=si+1
                                                    break
                                                si-=1
                
                                            while(True):
                                                if ei==len(template)-1:
                                                    break
                                                elif template[ei]==' ':
                                                    ei-=1
                                                    break
                                                ei+=1
                
                                            for r in range(si,ei+1):
                                                temp=temp+template[r]
                
                                            if temp.find('[') != -1:
                                                if temp.find(alpha[j])==-1:
                                                    for p in range(len(temp)-1,-1,-1):
                                                        if temp[p].isalpha():
                                                            s3=template[:si]+temp[:p]+alpha[j]+temp[p:]+template[ei+1:]
                                                            break
                                                else:
                                                    s3=template
                
                                            else:
                                                s3 = template[:si]+"["+temp+"]"+"["+alpha[j]+"]"+template[ei+1:]
                
                                         elif ebe_matrix_educts[k][j] == 2000:                  #For anions
                                            si=t
                                            ei=t
                                            temp=''
                                            while(True):
                                                if si==0:
                                                    break
                                                elif template[si]==' ':
                                                    si=si+1
                                                    break
                                                si-=1
                
                                            while(True):
                                                if ei==len(template)-1:
                                                    break
                                                elif template[ei]==' ':
                                                    ei-=1
                                                    break
                                                ei+=1
                
                                            for r in range(si,ei+1):
                                                temp=temp+template[r]
                
                                            if temp.find('[') != -1:
                                                if temp.find(alpha[j])==-1:
                                                    for p in range(len(temp)):
                                                        if temp[p].isalpha():
                                                            s3=template[:si]+temp[:p]+alpha[j]+temp[p:]+template[ei+1:]
                                                            break
                                                else:
                                                    s3=template
                
                                            else:
                                                s3 = template[:si]+"["+alpha[j]+"]"+"["+temp+"]"+template[ei+1:]
                
                                         if ebe_matrix_educts[k][j] != 1000 and ebe_matrix_educts[k][j] != 2000:
                                            template = s1+s3+s2
                                         else:
                                            template = s3
                                    
                                        
                                        
                                     elif brack_flag[k]==1:                     #For branched compound                        
                                         if k == j:
                                            if ebe_matrix_educts[k][k] != 0 and ebe_matrix_educts[k][k] != 1000 and ebe_matrix_educts[k][k] != 2000:
                                                s3 = alpha[k] + "(."+str(ebe_matrix_educts[k][k])+")"
                                            elif ebe_matrix_educts[k][k] == 0:
                                                f1 = 0
                                                for t in range(ed_len):
                                                    if ebe_matrix_educts[k][t] == 1000 or ebe_matrix_educts[k][t] == 2000:
                                                        f1 = 1
                                                        break
                                                if ebe_matrix_educts[k][ed_len] > 0 and f1 == 0:
                                                    s3 = alpha[k] + "(+"+str(ebe_matrix_educts[k][ed_len])+")"
                                                elif ebe_matrix_educts[k][ed_len] < 0 and f1 == 0:
                                                    s3 = alpha[k] + "("+str(ebe_matrix_educts[k][ed_len])+")"
                                                
                                                                                                                                     
                
                                         elif ebe_matrix_educts[k][j] == 1:         #For single bond
                                            s3 = s3+"(--"+alpha[j]+")"
                                         elif ebe_matrix_educts[k][j] == 2:         #For double bond
                                            s3 = s3+"(="+alpha[j]+")"
                                         elif ebe_matrix_educts[k][j] == 3:         #For triple bond
                                            s3 = s3+"(#"+alpha[j]+")"
                                         elif ebe_matrix_educts[k][j] == 100:       #For coordinate bond
                                            s3 = s3+"(~>"+alpha[j]+")"
                                         elif ebe_matrix_educts[k][j] == 200:       #For coordinate bond
                                            s3 = s3+"(<~"+alpha[j]+")"
                                         elif ebe_matrix_educts[k][j] == 500:                       #For addition compound
                                            fladdb=0
                                            for addb in range(len(s3)-1):
                                                if s3[addb]=='.' and s3[addb+1]=='(':
                                                    s3 = s3+alpha[k]+"("+alpha[j]+")"
                                                    fladdb=1
                                                    break
                                            if fladdb==0:
                                                s3 = s3+alpha[k]+".("+alpha[j]+")"  
                                         elif ebe_matrix_educts[k][j] == 1000:      #For cations
                                            si=t
                                            ei=t
                                            temp=''
                                            while(True):
                                                if si==0:
                                                    break
                                                elif template[si]==' ':
                                                    si=si+1
                                                    break
                                                si-=1
                
                                            while(True):
                                                if ei==len(template)-1:
                                                    break
                                                elif template[ei]==' ':
                                                    ei-=1
                                                    break
                                                ei+=1
                
                                            for r in range(si,ei+1):
                                                temp=temp+template[r]
                
                                            if temp.find('[') != -1:
                                                if temp.find(alpha[j])==-1:
                                                    for p in range(len(temp)-1,-1,-1):
                                                        if temp[p].isalpha():
                                                            s3=template[:si]+temp[:p]+alpha[j]+temp[p:]+template[ei+1:]
                                                            break
                                                else:
                                                    s3=template
                
                                            else:
                                                s3 = template[:si]+"["+temp+"]"+"["+alpha[j]+"]"+template[ei+1:]
                
                                         elif ebe_matrix_educts[k][j] == 2000:                      #For anions
                                            si=t
                                            ei=t
                                            temp=''
                                            while(True):
                                                if si==0:
                                                    break
                                                elif template[si]==' ':
                                                    si=si+1
                                                    break
                                                si-=1
                
                                            while(True):
                                                if ei==len(template)-1:
                                                    break
                                                elif template[ei]==' ':
                                                    ei-=1
                                                    break
                                                ei+=1
                
                                            for r in range(si,ei+1):
                                                temp=temp+template[r]
                
                                            if temp.find('[') != -1:
                                                if temp.find(alpha[j])==-1:
                                                    for p in range(len(temp)):
                                                        if temp[p].isalpha():
                                                            s3=template[:si]+temp[:p]+alpha[j]+temp[p:]+template[ei+1:]
                                                            break
                                                else:
                                                    s3=template
                
                                            else:
                                                s3 = template[:si]+"["+alpha[j]+"]"+"["+temp+"]"+template[ei+1:]
                
                                         if ebe_matrix_educts[k][j] != 1000 and ebe_matrix_educts[k][j] != 2000:
                                            template = s1+s3+s2
                                         else:
                                            template = s3
                                    
                        
                        
                        visited[k][j]=1
                        visited[j][k]=1
                        k=j
                        e=0
                        continue
                            
                    e=e+1   
                    
                
                
                st_top-=1
                stack_atom.pop()
                if st_top==-1:   
                    if template!='':
                        if ed_temp=='':
                            ed_temp=template
                        else:
                            ed_temp=ed_temp+" + "+template
                        template=''
                    break
                
                
                             
                e=0
                k=stack_atom[st_top]        
                
                
                    
                    

    
    for i in range(len(alpha)):
        if alpha[i]!='' and ed_temp.find(alpha[i])==-1:
            alpha[i]=''




    brack_flag=[0 for i in range(ed_len)]   
    
    '''for i in range(ed_len):
        cnt_ion=0
        if len(make_bond[i])!=0:
            for j in range(len(make_bond[i])):
                if make_bond[i][j]==1000 or make_bond[i][j]==2000:
                    cnt_ion+=1
            if len(make_bond[i])>cnt_ion+1:
                brack_flag[i]=1'''
                        
    visited=[[0 for i in range(ed_len)] for j in range(ed_len)]         #Visited matrix for DFS algorithm
    stack_atom=[]
    st_top=-1       

    prod_temp=''    #Variable to store ptoduct template
    
    template=''    #variable to store current template
    
    for i in range(ed_len):             #Loop to generate the product template
        
        if len(make_bond[i])!=0:
        
            
            st_top=st_top+1
            stack_atom.append(i)
            k=stack_atom[st_top]
            e=0
            while True:                  
                    
                while e<len(make_bond[k]):
                    if visited[k][make_bond[k][e]]==0:
                        ign=0
                        w=make_bond[k][e]
                        if template.find('[')!=-1 and template.find(']')!=-1 and (ebe_matrix_products[k][w]==1000 or ebe_matrix_products[k][w]==2000) and k!=w:
                            for z in make_bond[w]:
                                if visited[w][z]==0 and ebe_matrix_products[w][z]!=0 and ebe_matrix_products[w][z]!=1000 and ebe_matrix_products[w][z]!=2000:
                                    ign=1
                                    break
                        if ign==1:
                            e=e+1
                            continue
                        
                        if k!=make_bond[k][e]:
                            st_top+=1
                            stack_atom.append(make_bond[k][e])
                        j=stack_atom[st_top]
                        if alpha[j] == '':
                            alpha[j] = chr(char)
                            char = char+1
                        if alpha[k] == '':
                            alpha[k] = chr(char)
                            char = char+1
                            
                        cnt_ion=0
                        for x in make_bond[k]:
                            if x!= k and visited[x][k]==0 and (ebe_matrix_products[k][x]!=1000 or ebe_matrix_products[k][x]!=2000): 
                                cnt_ion+=1
                        if cnt_ion>1:
                            brack_flag[k]=1    
                            
                            
                        if template=='' or template.find(alpha[k])==-1:          #For unvisited elemets in products 
                            
                            if brack_flag[k]==1:                                #For branched compound
                                
                                 if k == j:                                     #For valency and charge
                                     
                                    if ebe_matrix_products[k][k] != 0 and ebe_matrix_products[k][k] != 1000 and ebe_matrix_products[k][k] != 2000:
                                        template += alpha[k] + "(."+str(ebe_matrix_products[k][k])+")"
                                    elif ebe_matrix_products[k][k] == 0:
                                        f1 = 0
                                        for t in range(ed_len):
                                            if ebe_matrix_products[k][t] == 1000 or ebe_matrix_products[k][t] == 2000:
                                                f1 = 1
                                                break
                                        if ebe_matrix_products[k][ed_len] > 0 and f1 == 0:
                                           template += alpha[k] + "(+"+str(ebe_matrix_products[k][ed_len])+")"
                                        elif ebe_matrix_products[k][ed_len] < 0 and f1 == 0:
                                           template += alpha[k] + "("+str(ebe_matrix_products[k][ed_len])+")"                       
                                                                                    
                                                   
                                 elif ebe_matrix_products[k][j] == 1:                   #For Single bond 
                                    template = template +alpha[k]+"(--"+alpha[j]+")"
                                 elif ebe_matrix_products[k][j] == 2:                   #For Double bond 
                                    template =template +alpha[k]+"(="+alpha[j]+")"
                                 elif ebe_matrix_products[k][j] == 3:                   #For Triple bond 
                                    template = template+alpha[k]+"(#"+alpha[j]+")"
                                 elif ebe_matrix_products[k][j] == 100:                 #For coordinate bond 
                                    template =  template+alpha[k]+"(~>"+alpha[j]+")"
                                 elif ebe_matrix_products[k][j] == 200:                 #For coordinate bond         
                                    template =  template+alpha[k]+"(<~"+alpha[j]+")"
                                 elif ebe_matrix_products[k][j] == 500:                       #For addition compound
                                    fladdb=0
                                    for addb in range(len(template)-1):
                                        if template[addb]=='.' and template[addb+1]=='(':
                                            template =template+alpha[k]+"("+alpha[j]+")"
                                            fladdb=1
                                            break
                                    if fladdb==0:
                                        template = template+alpha[k]+".("+alpha[j]+")"  
                                 elif ebe_matrix_products[k][j] == 1000:                #For cations
                                    template=template+'['+alpha[k]+']'+'['+alpha[j]+']'
                                 elif ebe_matrix_products[k][j] == 2000:                #For anions
                                    template=template+'['+alpha[j]+']'+'['+alpha[k]+']'
                            
                            elif brack_flag[k]==0:                                      #For unbranched compound
                                
                                 if k == j:                                         #For valency and charge
                                     
                                    if ebe_matrix_products[k][k] != 0 and ebe_matrix_products[k][k] != 1000 and ebe_matrix_products[k][k] != 2000:
                                        template += alpha[k] + "(."+str(ebe_matrix_products[k][k])+")"
                                    elif ebe_matrix_products[k][k] == 0:
                                        f1 = 0
                                        for t in range(ed_len):
                                            if ebe_matrix_products[k][t] == 1000 or ebe_matrix_products[k][t] == 2000:
                                                f1 = 1
                                                break
                                        if ebe_matrix_products[k][ed_len] > 0 and f1 == 0:
                                            template += alpha[k] + "(+"+str(ebe_matrix_products[k][ed_len])+")"
                                        elif ebe_matrix_products[k][ed_len] < 0 and f1 == 0:
                                            template += alpha[k] + "("+str(ebe_matrix_products[k][ed_len])+")"
                                                                                                            
                                            
                                 elif ebe_matrix_products[k][j] == 1:                       #For single bond 
                                    template = template + alpha[k]+"--"+alpha[j]
                                 elif ebe_matrix_products[k][j] == 2:                       #For double bond 
                                    template =template +alpha[k]+"="+alpha[j]
                                 elif ebe_matrix_products[k][j] == 3:                       #For triple bond 
                                    template = template+alpha[k]+"#"+alpha[j]
                                 elif ebe_matrix_products[k][j] == 100:                     #For coordinate bond 
                                    template =  template+alpha[k]+"~>"+alpha[j]
                                 elif ebe_matrix_products[k][j] == 200:                     #For coordinate bond 
                                    template =  template+alpha[k]+"<~"+alpha[j]
                                 elif ebe_matrix_products[k][j] == 500:                       #For addition compound
                                    fladdb=0
                                    for addb in range(len(template)-1):
                                        if template[addb]=='.' and template[addb+1]=='(':
                                            template =template+alpha[k]+"("+alpha[j]+")"
                                            fladdb=1
                                            break
                                    if fladdb==0:
                                        template = template+alpha[k]+".("+alpha[j]+")"  
                                 elif ebe_matrix_products[k][j] == 1000:                    #For cations
                                    template=template+'['+alpha[k]+']'+'['+alpha[j]+']'
                                 elif ebe_matrix_products[k][j] == 2000:                    #For anions
                                    template=template+'['+alpha[j]+']'+'['+alpha[k]+']'    
                        
                       
                        else:                                       #For visited elements in products
                            l = len(template)
                            for t in range(l):
                                if template[t]==alpha[k]:
                                     s1 = template[:t]
                                     s2 = template[t+1:]
                                     s3 = template[t]
                                     
                                         
                                     if brack_flag[k]==0:           #For unbranched compound
                                         
                                         if k == j:                 #For valency and charge
                                             
                                            if ebe_matrix_products[k][k] != 0 and ebe_matrix_products[k][k] != 1000 and ebe_matrix_products[k][k] != 2000:
                                                s3 = alpha[k] + "(."+str(ebe_matrix_products[k][k])+")"
                                            elif ebe_matrix_products[k][k] == 0:
                                                f1 = 0
                                                for t in range(ed_len):
                                                    if ebe_matrix_products[k][t] == 1000 or ebe_matrix_products[k][t] == 2000:
                                                        f1 = 1
                                                        break
                                                if ebe_matrix_products[k][ed_len] > 0 and f1 == 0:
                                                    s3 = alpha[k] + "(+"+str(ebe_matrix_products[k][ed_len])+")"
                                                elif ebe_matrix_products[k][ed_len] < 0 and f1 == 0:
                                                    s3 = alpha[k] + "("+str(ebe_matrix_products[k][ed_len])+")"
                                                
                                                                                                                                     
                
                                         elif ebe_matrix_products[k][j] == 1:                   #For single bond        
                                            s3 = s3+"--"+alpha[j]
                                         elif ebe_matrix_products[k][j] == 2:                   #For double bond 
                                            s3 = s3+"="+alpha[j]
                                         elif ebe_matrix_products[k][j] == 3:                   #For triple bond 
                                            s3 = s3+"#"+alpha[j]
                                         elif ebe_matrix_products[k][j] == 100:                 #For coordinate bond 
                                            s3 = s3+"~>"+alpha[j]
                                         elif ebe_matrix_products[k][j] == 200:                 #For coordinate bond 
                                            s3 = s3+"<~"+alpha[j]
                                         elif ebe_matrix_products[k][j] == 500:                       #For addition compound
                                            fladdb=0
                                            for addb in range(len(template)-1):
                                                if template[addb]=='.' and template[addb+1]=='(':
                                                    template =template+alpha[k]+"("+alpha[j]+")"
                                                    fladdb=1
                                                    break
                                            if fladdb==0:
                                                template = template+alpha[k]+".("+alpha[j]+")"  
                                         elif ebe_matrix_products[k][j] == 1000:                #For cations
                                            si=t
                                            ei=t
                                            temp=''
                                            while(True):
                                                if si==0:
                                                    break
                                                elif template[si]==' ':
                                                    si=si+1
                                                    break
                                                si-=1
                
                                            while(True):
                                                if ei==len(template)-1:
                                                    break
                                                elif template[ei]==' ':
                                                    ei-=1
                                                    break
                                                ei+=1
                
                                            for r in range(si,ei+1):
                                                temp=temp+template[r]
                
                                            if temp.find('[') != -1:
                                                if temp.find(alpha[j])==-1:
                                                    for p in range(len(temp)-1,-1,-1):
                                                        if temp[p].isalpha():
                                                            s3=template[:si]+temp[:p]+alpha[j]+temp[p:]+template[ei+1:]
                                                            break
                                                else:
                                                    s3=template
                
                                            else:
                                                s3 = template[:si]+"["+temp+"]"+"["+alpha[j]+"]"+template[ei+1:]
                
                                         elif ebe_matrix_products[k][j] == 2000:                        #For anions
                                            si=t
                                            ei=t
                                            temp=''
                                            while(True):
                                                if si==0:
                                                    break
                                                elif template[si]==' ':
                                                    si=si+1
                                                    break
                                                si-=1
                
                                            while(True):
                                                if ei==len(template)-1:
                                                    break
                                                elif template[ei]==' ':
                                                    ei-=1
                                                    break
                                                ei+=1
                
                                            for r in range(si,ei+1):
                                                temp=temp+template[r]
                
                                            if temp.find('[') != -1:
                                                if temp.find(alpha[j])==-1:
                                                    for p in range(len(temp)):
                                                        if temp[p].isalpha():
                                                            s3=template[:si]+temp[:p]+alpha[j]+temp[p:]+template[ei+1:]
                                                            break
                                                else:
                                                    s3=template
                
                                            else:
                                                s3 = template[:si]+"["+alpha[j]+"]"+"["+temp+"]"+template[ei+1:]
                
                                         if ebe_matrix_products[k][j] != 1000 and ebe_matrix_products[k][j] != 2000:
                                            template = s1+s3+s2
                                         else:
                                            template = s3
                                         
                                        
                                        
                                     elif brack_flag[k]==1:                             #For branched compound
                                         
                                        if k == j:                                      #For valency and charge
                                            
                                            if ebe_matrix_products[k][k] != 0 and ebe_matrix_products[k][k] != 1000 and ebe_matrix_products[k][k] != 2000:
                                                s3 = alpha[k] + "(."+str(ebe_matrix_products[k][k])+")"
                                            elif ebe_matrix_products[k][k] == 0:
                                                f1 = 0
                                                for t in range(ed_len):
                                                    if ebe_matrix_products[k][t] == 1000 or ebe_matrix_products[k][t] == 2000:
                                                        f1 = 1
                                                        break
                                                if ebe_matrix_products[k][ed_len] > 0 and f1 == 0:
                                                    s3 = alpha[k] + "(+"+str(ebe_matrix_products[k][ed_len])+")"
                                                elif ebe_matrix_products[k][ed_len] < 0 and f1 == 0:
                                                    s3 = alpha[k] + "("+str(ebe_matrix_products[k][ed_len])+")"
                                                
                                                                                                                                     
                
                                        elif ebe_matrix_products[k][j] == 1:                    #For single bond 
                                            s3 = s3+"(--"+alpha[j]+")"
                                        elif ebe_matrix_products[k][j] == 2:                    #For double bond 
                                            s3 = s3+"(="+alpha[j]+")"
                                        elif ebe_matrix_products[k][j] == 3:                    #For triple bond 
                                            s3 = s3+"(#"+alpha[j]+")"
                                        elif ebe_matrix_products[k][j] == 100:                  #For coordinate bond 
                                            s3 = s3+"(~>"+alpha[j]+")"
                                        elif ebe_matrix_products[k][j] == 200:                  #For coordinate bond 
                                            s3 = s3+"(<~"+alpha[j]+")"
                                        elif ebe_matrix_products[k][j] == 500:                       #For addition compound
                                            fladdb=0
                                            for addb in range(len(template)-1):
                                                if template[addb]=='.' and template[addb+1]=='(':
                                                    template =template+alpha[k]+"("+alpha[j]+")"
                                                    fladdb=1
                                                    break
                                            if fladdb==0:
                                                template = template+alpha[k]+".("+alpha[j]+")"  
                                        elif ebe_matrix_products[k][j] == 1000:                 #For cations 
                                            si=t
                                            ei=t
                                            temp=''
                                            while(True):
                                                if si==0:
                                                    break
                                                elif template[si]==' ':
                                                    si=si+1
                                                    break
                                                si-=1
                
                                            while(True):
                                                if ei==len(template)-1:
                                                    break
                                                elif template[ei]==' ':
                                                    ei-=1
                                                    break
                                                ei+=1
                
                                            for r in range(si,ei+1):
                                                temp=temp+template[r]
                
                                            if temp.find('[') != -1:
                                                if temp.find(alpha[j])==-1:
                                                    for p in range(len(temp)-1,-1,-1):
                                                        if temp[p].isalpha():
                                                            s3=template[:si]+temp[:p]+alpha[j]+temp[p:]+template[ei+1:]
                                                            break
                                                else:
                                                    s3=template
                
                                            else:
                                                s3 = template[:si]+"["+temp+"]"+"["+alpha[j]+"]"+template[ei+1:]
                
                                        elif ebe_matrix_products[k][j] == 2000:                     #For anions
                                            si=t
                                            ei=t
                                            temp=''
                                            while(True):
                                                if si==0:
                                                    break
                                                elif template[si]==' ':
                                                    si=si+1
                                                    break
                                                si-=1
                
                                            while(True):
                                                if ei==len(template)-1:
                                                    break
                                                elif template[ei]==' ':
                                                    ei-=1
                                                    break
                                                ei+=1
                
                                            for r in range(si,ei+1):
                                                temp=temp+template[r]
                
                                            if temp.find('[') != -1:
                                                if temp.find(alpha[j])==-1:
                                                    for p in range(len(temp)):
                                                        if temp[p].isalpha():
                                                            s3=template[:si]+temp[:p]+alpha[j]+temp[p:]+template[ei+1:]
                                                            break
                                                else:
                                                    s3=template
                
                                            else:
                                                s3 = template[:si]+"["+alpha[j]+"]"+"["+temp+"]"+template[ei+1:]
                
                                        if ebe_matrix_products[k][j] != 1000 and ebe_matrix_products[k][j] != 2000:
                                            template = s1+s3+s2
                                        else:
                                            template = s3
                                            
                        
                        
                        visited[k][j]=1             #marking elements as visited 
                        visited[j][k]=1
                        k=j
                        e=0
                        continue
                            
                    e=e+1   
                
                st_top-=1
                stack_atom.pop()
                if st_top==-1:    
                    if template!='':
                        if prod_temp=='':
                            prod_temp=template
                        else:
                            prod_temp=prod_temp+" + "+template
                        template=''
                    break
                
                
                e=0
                k=stack_atom[st_top] 
                
                
    

    usles_list=[]               #List of elements not required in template of educts and products

    for i in range(1,len(ed_temp)-1):       #Eliminating useless symbols in educts
        if ed_temp[i].isalpha():
            if ed_temp[i+1].isalpha() and ed_temp[i-1].isalpha():
                usles_list.append(ed_temp[i])
            elif ed_temp[i+1].isalpha() and ed_temp[i-1]=='[':
                usles_list.append(ed_temp[i])
            elif ed_temp[i-1].isalpha() and ed_temp[i+1]==']':
                usles_list.append(ed_temp[i])


    s=prod_temp
    
    for i in range(1,len(s)-1):             #Eliminating useless symbols in products
        if s[i].isalpha():
            if s[i+1].isalpha() and s[i-1].isalpha() and usles_list.count(s[i])==1:
                ed_temp=ed_temp.replace(s[i],"")
                prod_temp=prod_temp.replace(s[i],"")
                for k in range(len(alpha)):
                    if alpha[k]==s[i]:
                        for p in range(ed_len):
                            reaction_matrix[p][k]=0
                            reaction_matrix[k][p]=0

            elif s[i+1].isalpha() and s[i-1]=='[' and usles_list.count(s[i])==1:
                ed_temp=ed_temp.replace(s[i],"")
                prod_temp=prod_temp.replace(s[i],"")
                for k in range(len(alpha)):
                    if alpha[k]==s[i]:
                        for p in range(ed_len):
                            reaction_matrix[p][k]=0
                            reaction_matrix[k][p]=0

                
            elif s[i-1].isalpha() and s[i+1]==']' and usles_list.count(s[i])==1:
                ed_temp=ed_temp.replace(s[i],"")
                prod_temp=prod_temp.replace(s[i],"")
                for k in range(len(alpha)):
                    if alpha[k]==s[i]:
                        for p in range(ed_len):
                            reaction_matrix[p][k]=0
                            reaction_matrix[k][p]=0




    els=[]                              #List of elements preent in educts but not in products
    for i in range(len(ed_temp)):
        s=ed_temp[i]
        if s.isalpha() and prod_temp.find(s)==-1:
            els.append(s)
            
            
    for i in els:
        ed_temp=ed_temp.replace(i,'')
        for j in range(len(alpha)):
            if alpha[j]==i:
                alpha[j]=''
        
    els=[]                              #List of elements preent in products but not in educts
    for i in range(len(prod_temp)):
        s=prod_temp[i]
        if s.isalpha() and ed_temp.find(s)==-1:
            els.append(s)
            
            
    for i in els:
        prod_temp=prod_temp.replace(i,'')
        for j in range(len(alpha)):
            if alpha[j]==i:
                alpha[j]=''
            
        
        
        
    #print(ed_temp)
    #print(prod_temp)

    remove_atom=[]                              
    for i in range(len(reaction_matrix)):               #Removing all useless elements from educt, product and reaction matrix
        react_flag=0
        for j in range(len(reaction_matrix[i])):
            if reaction_matrix[i][j]!=0:
                react_flag=1
                break
        if react_flag==0:
            remove_atom.append(i)
            
    remove_atom=list(set(remove_atom))

    for i in range(len(remove_atom)):
        remove_atom[i]=remove_atom[i]-i
            
    for i in remove_atom:
        del(reaction_matrix[i])
        del(ebe_matrix_educts[i])
        del(ebe_matrix_products[i])

    for i in range(len(reaction_matrix)):
        for j in remove_atom:
            del(reaction_matrix[i][j])
            del(ebe_matrix_educts[i][j])
            del(ebe_matrix_products[i][j])
    
    #print(remove_atom)
    #print("\n")
    #print(ebe_matrix_educts)
    #print("\n")
    #print(ebe_matrix_products)
    #print("\n")
    #print(reaction_matrix)
    #print("\n")
        
        
    template = ed_temp + " ----> " + prod_temp              # Final Template of Reaction

    template = template.replace('{', '(')
    template = template.replace('}', ')')
    
    ed_temp = ed_temp.replace('{', '(')
    ed_temp = ed_temp.replace('}', ')')
    
    prod_temp = prod_temp.replace('{', '(')
    prod_temp = prod_temp.replace('}', ')')
    
    
    
    
    #temp_row=[ed_temp,ebe_matrix_educts,reaction_matrix,ebe_matrix_products,prod_temp]
    
    for x in range(len(ebe_matrix_educts)):
        ebe_matrix_educts[x][len(ebe_matrix_educts)+1]=0
        ebe_matrix_products[x][len(ebe_matrix_products)+1]=0
            
    
    for x in range(len(ebe_matrix_educts)):                                 #For placing proper values in coordinate bond column in educt matrix
        for y in range(len(ebe_matrix_educts)):
            if ebe_matrix_educts[x][y] == 100:
                    if ebe_matrix_educts[x][len(ebe_matrix_educts)+1]!=0:
                        ebe_matrix_educts[x][len(ebe_matrix_educts)+1]=-99
                    else:
                        ebe_matrix_educts[x][len(ebe_matrix_educts)+1] = y+1

                    if ebe_matrix_educts[y][len(ebe_matrix_educts)+1]!=0:
                        ebe_matrix_educts[y][len(ebe_matrix_educts)+1]=-99
                    else:
                        ebe_matrix_educts[y][len(ebe_matrix_educts)+1] = x+1
            
    for x in range(len(ebe_matrix_products)):                               #For placing proper values in coordinate bond column in educt matrix
        for y in range(len(ebe_matrix_products)):
            if ebe_matrix_products[x][y] == 100:
                    if ebe_matrix_products[x][len(ebe_matrix_products)+1]!=0:
                        ebe_matrix_products[x][len(ebe_matrix_products)+1]=-99
                    else:
                        ebe_matrix_products[x][len(ebe_matrix_products)+1] = y+1

                    if ebe_matrix_products[y][len(ebe_matrix_products)+1]!=0:
                        ebe_matrix_products[y][len(ebe_matrix_products)+1]=-99
                    else:
                        ebe_matrix_products[y][len(ebe_matrix_products)+1] = x+1

    temp_alpha=[]
    for i in alpha:
        if i.isalpha():
            temp_alpha.append(i)

    
    temp_row=[]
    k=[x for x in range(len(ebe_matrix_educts))]
    
    '''
    b=[x for x in range(1)]
    
    permuted_educts=[[0 for x in range(len(ebe_matrix_educts)+2)] for y in range(len(ebe_matrix_educts))] 
    permuted_products=[[0 for x in range(len(ebe_matrix_educts)+2)] for y in range(len(ebe_matrix_educts))] 
    
    perm = permutations(b)                          #finding combinations matrix of all elements in template 
    
    for i in list(perm): 
        permuted_educts=[[0 for x in range(len(ebe_matrix_educts)+2)] for y in range(len(ebe_matrix_educts))] 
        permuted_products=[[0 for x in range(len(ebe_matrix_educts)+2)] for y in range(len(ebe_matrix_educts))] 
        permuted_reaction=[[0 for x in range(len(ebe_matrix_educts))] for y in range(len(ebe_matrix_educts))]
        
        k=list(i)  
        p=0
        q=0
        for x in k:
            q=0
            for y in k:
                permuted_educts[p][q]=ebe_matrix_educts[x][y]
                permuted_educts[p][len(ebe_matrix_educts)]=ebe_matrix_educts[x][len(ebe_matrix_educts)]
                permuted_educts[q][len(ebe_matrix_educts)]=ebe_matrix_educts[y][len(ebe_matrix_educts)]
                               
                permuted_products[p][q]=ebe_matrix_products[x][y]
                permuted_products[p][len(ebe_matrix_educts)]=ebe_matrix_products[x][len(ebe_matrix_educts)]
                permuted_products[q][len(ebe_matrix_educts)]=ebe_matrix_products[y][len(ebe_matrix_educts)]
                
                permuted_reaction[p][q]=reaction_matrix[x][y]
                
                q+=1
            p+=1
        
        for x in range(len(ebe_matrix_educts)):                                     #For placing proper values in coordinate bond column in educt matrix
            for y in range(len(ebe_matrix_educts)):
                if permuted_educts[x][y] == 100:
                        if permuted_educts[x][len(ebe_matrix_educts)+1]!=0:
                            permuted_educts[x][len(ebe_matrix_educts)+1]=-99
                        else:
                            permuted_educts[x][len(ebe_matrix_educts)+1] = y+1
    
                        if permuted_educts[y][len(ebe_matrix_educts)+1]!=0:
                            permuted_educts[y][len(ebe_matrix_educts)+1]=-99
                        else:
                            permuted_educts[y][len(ebe_matrix_educts)+1] = x+1
            
        for x in range(len(ebe_matrix_products)):                                   #For placing proper values in coordinate bond column in educt matrix
            for y in range(len(ebe_matrix_products)):
                if permuted_products[x][y] == 100:
                        if permuted_products[x][len(ebe_matrix_products)+1]!=0:
                            permuted_products[x][len(ebe_matrix_products)+1]=-99
                        else:
                            permuted_products[x][len(ebe_matrix_products)+1] = y+1
    
                        if permuted_products[y][len(ebe_matrix_products)+1]!=0:
                            permuted_products[y][len(ebe_matrix_products)+1]=-99
                        else:
                           permuted_products[y][len(ebe_matrix_products)+1] = x+1
        
        #print(permuted_educts)
        #print(permuted_products)

        '''
                           
    temp_row.append([ed_temp,ebe_matrix_educts,reaction_matrix,ebe_matrix_products,prod_temp,temp_alpha,k])   #List to be written in templates file
        
    
    
    
   
    
    cv_row.append(ed_temp)
    cv_row.append(prod_temp)

    rt = Tk()            

    rt.config(bg="white")                                                       #Window to view the output template and reaction

    w, h = rt.winfo_screenwidth(), rt.winfo_screenheight()
    rt.geometry("%dx%d+0+0" % (w, h))
    rt.title("Template")



    Label(rt, text="\n\nReaction\n\n", font=("Arial Bold", 20),bg="white").pack()          #label for the reaction

    edc=StringVar(rt,value=reaction)

    M1=Entry(rt, textvariable=edc,state='readonly',width=70,font=("Arial Bold", 20))
    M2=Scrollbar(rt, orient='horizontal', command=M1.xview)
    M1.config(xscrollcommand=M2.set)


    M1.pack()
    M2.pack()

    Label(rt, text="\n\n\n\n\n",bg="white").pack()


    Label(rt, text="Template\n\n", font=("Arial Bold", 20), bg="white").pack()          #label for the template

    edt=StringVar(rt,value=template)

    M3=Entry(rt, textvariable=edt,state='readonly',width=70,font=("Arial Bold", 20))
    M4=Scrollbar(rt, orient='horizontal', command=M3.xview)
    M3.config(xscrollcommand=M4.set)


    M3.pack()
    M4.pack()

    Label(rt, text="\n\n\n\n\n",bg="white").pack()

    but = Button(rt, text='EXIT', command=lambda: run(rt),bg="white")              #calling method to exut
    but.place(x=130, y=560)
    
    but3 = Button(rt, text='INSERT INTO DATABASE', command=lambda: database_entry(cv_row,temp_row,cov,bnd,parts),bg="white")       #calling method to insert the reaction and template in datavbase
    but3.place(x=590, y=560)


    but2 = Button(rt, text='REFERENCE', command=lambda: reference(template,alpha),bg="white")              #Method to call the reference window
    but2.place(x=1130, y=560)
    



window = Tk()                   #First window for the user
ed = ''
prod = ''
educt = StringVar()             #User input for educts
product = StringVar()           #User input for products
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
window.title("ESMILES INPUT")           #Title of the Window Screen

window.configure(bg="white")

Label(window, text="Welcome", font=("Arial Bold", 50), fg = "blue", bg = "white").pack()                                                                                   #
Label(window, text="Please Enter the Chemical Reactions in the text box in the ESMILES Notation:", font=("Arial Bold", 15),fg = "dark green", bg = "white").place(x=40, y=140)       #
                                                                                    
Label(window, text="Please click on the button to read the Manual and learn how to input chemical compounds and elements in ESMILES notation ", font=("Arial Bold", 10),fg = "dark green", bg = "white").place(x=60, y=200)                       #

b = Button(window, text='Read Manual', command=read_manual,fg = "blue", bg = "white").place(x=900, y=190)           # Calling "go_press" function after button click


Label(window, text="Educt",fg = "black", bg = "white").place(x=40, y=290)


E1 = Entry(window, textvariable=educt,width=70,font=("Arial Bold", 15))         # Textbox for Reading Educt Compounds
E1.place(x=100, y=290)


Label(window, text="Product",fg = "black", bg = "white").place(x=40, y=390)                        

E2 = Entry(window, textvariable=product,width=70,font=("Arial Bold", 15))       # Textbox for Reading Product Compounds
E2.place(x=100, y=390)


btn = Button(window, text='GO', command=ESMILES_to_EBE,fg = "dark green", bg = "white").place(x=500, y=460)           # Calling "go_press" function after button click

#print(os.name)

ele = pandas.read_csv("elementlist.csv", header=None)                           #Reading Elemenets name and Valency from file

elements = ele[1].values.tolist()                                               #List containing all Elements symbol in Periodic table
elements_name=ele[2].values.tolist()                                            #List containing all Elements name in Periodic table
elements_valency = ele[3].values.tolist()                                       #List containing all Elements valency in Periodic table
elements_free_electons = ele[4].values.tolist()
elements_stable_electons = ele[5].values.tolist()
elements_max_bond = ele[6].values.tolist()
elements_eneg = ele[7].values.tolist()


window.mainloop()                                                               #End of window




