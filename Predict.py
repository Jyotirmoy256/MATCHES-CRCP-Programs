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

def back(c):                                # Function to Run New Program
    c.destroy()
    os.system('python Predict.py')
    
def exit_os(c):                                # Function to Run New Program
    c.destroy()

def stable_products(est_products):          #Method to process products predicted by the algorithm
    
    global count_ed
    root = Tk()
    w, h = root.winfo_screenwidth(), root.winfo_screenheight()

    root.geometry("%dx%d+0+0" % (w, h))
    
    root.config(bg="white")
    
    root.title("Output Products")           #Title of the Window Screen
    
    Label(root, text="Predicted Products", font=("Arial Bold", 30), fg="dark green", bg="white").pack()
    
    scrollbar = Scrollbar(root)
    scrollbar.pack(side = RIGHT)

    mylist = Listbox(root, yscrollcommand = scrollbar.set, font=("Arial Bold", 20))
    
    
    
    output_products=est_products[:]
    vx=-1   
    
    for ed in est_products:                  #Loop to process the educt string input given by the user
        
        vx=vx+1

        if (ed.find(' + ') == -1):
            ed = [ed]
        else:
            ed = ed.split(' + ')

        ed_compounds = ed           #Storing the individual Educt compound in a separate variables
        s = ''
        stack=[]                    #Declaring an empty Stack
        top=-1

        total_ed_matrix=[]          #List for storing EBE matrices of individual educt compounds
        total_prod_matrix=[]        #List for storing EBE matrices of individual product compounds

        ed_elements=[]              #List for storing elements of educts
        prod_elements=[]            #List for storing elements of products

        
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
                            ebe_matrix_educts[x][x] = (
                                ebe_matrix_educts[x][x]*-1) % 10
    
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
        
        for i in ed_elements:
            for j in range(len(elements)):
                if i==elements[j]:
                    max_bond.append(elements_max_bond[j])
                    electronegativity.append(elements_eneg[j])
        
        
        
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
                elif a==200:
                    a=0
                elif a==500 or a==1000 or a==2000:
                    a=0
            
                shared+=a
                
           if shared>max_bond[i]:
                rem_flag=1
                break
             
           
                       
        if rem_flag==0:
            if len(ebe_matrix_educts)==count_ed:
                #print(output_products[vx])
                mylist.insert(END, output_products[vx])
                mylist.insert(END,"")
            
                
    mylist.pack(fill = BOTH, expand=0)
    
    scrollbar.config( command = mylist.yview )
    
    
    window.destroy()
    
    btn = Button(root, text='Back', command=lambda: back(root), fg="dark green", bg="white")
    btn.place(x=500,y=400)
    
    btn1 = Button(root, text='Exit', command=lambda: exit_os(root), fg="dark green", bg="white")
    btn1.place(x=700,y=400)
    

def print_products():                      # Method top print predicted products in the list
    corrected_esmiles=[]
    fail=0
    for i in list(set(product_esmiles)):        #Loop to eliminate products not having proper ESMILES notation
        fail=0
        for j in range(len(i)):
            if j>0 and i[j].isalpha() and (i[j-1]==')' or i[j-1]==']'):
                fail=1
                break
        if fail==0:
            corrected_esmiles.append(i)
  
    if corrected_esmiles==[]:
        showerror("Error", "No Products Predicted")
    else:        
        stable_products(corrected_esmiles)
    
            

def match(temp_mat,ebe_matrix_educts,prod_mat,matched_temp,prod_esm,temp_elements,permutes):           #Method to generate ebe matrix of products
       
    original_matched=[]
    
    for i in matched_temp:                          #Loop to match all elements of templates with ebe matrix of educts
        nf=0
        x=0
        for j in i:
            y=0
            for k in i:
                if temp_mat[x][y]!=ebe_matrix_educts[j][k]:
                    nf=1
                    break
                y+=1
            if nf==1:
                break
            x+=1
        if nf==0:
             original_matched.append(i)         #Storing the matched templates in a list
    


             
    for f in original_matched:                  #Generating ebe matrix of products
               
        
        ed_len=len(ed_elements)       

        predicted_product= [[0 for i in range(ed_len+2)] for j in range(ed_len)]        #EBE matrix of products
        for q in range(len(ed_elements)):
           predicted_product[q]=ebe_matrix_educts[q][:]
           
       
        x=0
        y=0
        for j in f:                                                               #Filling up the ebe matrix of products                  
            y=0
            for k in f:
                
                predicted_product[j][k]=prod_mat[x][y]
                    
                predicted_product[j][ed_len]=prod_mat[x][len(prod_mat)]
                predicted_product[k][ed_len]=prod_mat[y][len(prod_mat)]
                
                y+=1
            x+=1
            
            
        for x in range(len(ed_elements)):                                           #Filling up the coordinate bond column correctly in educt matrix properly
            for y in range(len(ed_elements)):
                if predicted_product[x][y] == 100:
                        if predicted_product[x][len(ed_elements)+1]!=0:
                           predicted_product[x][len(ed_elements)+1]=-99
                        else:
                            predicted_product[x][len(ed_elements)+1] = y+1
    
                        if predicted_product[y][len(ed_elements)+1]!=0:
                            predicted_product[y][len(ed_elements)+1]=-99
                        else:
                            predicted_product[y][len(ed_elements)+1] = x+1
                    
       
        
    
        alpha=[chr(65+i) for i in range(len(ed_elements))]          #List to store alphabets A-Z
        
        
        brack_flag=[0 for i in range(ed_len)]                       #List for branched compounds
                        
        visited=[[0 for i in range(len(ed_elements))] for j in range(len(ed_elements))]     #visited matrix for dfs algorithm
        stack_atom=[]
        st_top=-1       
        prod_temp=''
        template=''    
        
        bonded_elements=[[] for i in range(len(ed_elements))]            #List to store elements for generating product in esmiles notation
        
        
        
        
        for i in range(len(predicted_product)):                         #Loop to fill up the bonded_elements list       
            for j in range(len(predicted_product[i][:-2])):
                if predicted_product[i][j]!=0:
                    bonded_elements[i].append(j)
                    if i!=j:
                        bonded_elements[j].append(i)

                    
                        
        for i in bonded_elements:                                       #Sorting and removing duplicates in bonded_elements list
            res = list(set(i))
            i[:] = res[:]
            i.sort()
            if len(i) > 1:
                j = i[0]
                if j == bonded_elements.index(i):
                    r = i[1:]
                    r.append(j)
                    i[:] = r[:]
                        
       
        '''print(ed_elements)
        print(predicted_product)
        print("\n")
        print(bonded_elements)
        print("\n\n\n")'''
        
        for i in range(len(ed_elements)):        #Loop to generate ESMILES notation of product
            if len(bonded_elements[i])!=0:    
                st_top=st_top+1
                stack_atom.append(i)
                k=stack_atom[st_top]
                e=0
                while True:
                    while e<len(bonded_elements[k]):
                        if visited[k][bonded_elements[k][e]]==0:
                            ign=0
                            w=bonded_elements[k][e]
                            if template.find('[')!=-1 and template.find(']')!=-1 and (predicted_product[k][w]==1000 or predicted_product[k][w]==2000) and k!=w:
                                for z in bonded_elements[w]:
                                    if visited[w][z]==0 and predicted_product[w][z]!=0 and predicted_product[w][z]!=1000 and predicted_product[w][z]!=2000:
                                        ign=1
                                        break
                            if ign==1:
                                e=e+1
                                continue
                            
                            if k!=bonded_elements[k][e]:
                                st_top+=1
                                stack_atom.append(bonded_elements[k][e])
                                
                            j=stack_atom[st_top]
                            
                            
                                        
                                
                            cnt_ion=0
                            for x in bonded_elements[k]:
                               if x!= k and visited[x][k]==0 and (predicted_product[k][x]!=1000 or predicted_product[k][x]!=2000): 
                                   cnt_ion+=1
                            if cnt_ion>1:
                               brack_flag[k]=1  
                                      
                            if template=='' or template.find(alpha[k])==-1:        #For unvisited elements in products                  
                                
                                if brack_flag[k]==1:                               #For branched compounds
                                    
                                     if k == j:                                     #For valency or charge
                                         
                                        if predicted_product[k][k] != 0 and predicted_product[k][k] != 1000 and predicted_product[k][k] != 2000:
                                            template += alpha[k] + "."+str(predicted_product[k][k])
                                        elif predicted_product[k][k] == 0:
                                            f1 = 0
                                            for t in range(len(ed_elements)):
                                                if predicted_product[k][t] == 1000 or predicted_product[k][t] == 2000:
                                                    f1 = 1
                                                    break
                                            if predicted_product[k][len(ed_elements)] > 0 and f1 == 0:
                                                template = template + alpha[k] + "+"+str(predicted_product[k][len(ed_elements)])
                                            elif predicted_product[k][len(ed_elements)] < 0 and f1 == 0:
                                                template = template + alpha[k] +str(predicted_product[k][len(ed_elements)])                      
                                                                                        
                                                       
                                     elif predicted_product[k][j] == 1:                         #For Single bond
                                        template = template +alpha[k]+"(-"+alpha[j]+")"
                                     elif predicted_product[k][j] == 2:                         #For Double bond
                                        template =template +alpha[k]+"(="+alpha[j]+")"
                                     elif predicted_product[k][j] == 3:                         #For Triple bond
                                        template = template+alpha[k]+"(#"+alpha[j]+")"
                                     elif predicted_product[k][j] == 100:                       #For coordinate bond
                                        template =  template+alpha[k]+"(~"+alpha[j]+")"
                                     elif predicted_product[k][j] == 200:                       #For coordinate bond
                                        template =  template+alpha[k]+"(~/"+alpha[j]+")"
                                     elif predicted_product[k][j] == 500:                       #For addition compound
                                        template = template+alpha[k]+"*"+alpha[j]
                                     elif predicted_product[k][j] == 1000:                      #For cations
                                        template=template+'['+alpha[k]+']'+'['+alpha[j]+']'
                                     elif predicted_product[k][j] == 2000:                      #For anions
                                        template=template+'['+alpha[j]+']'+'['+alpha[k]+']'
                                
                                elif brack_flag[k]==0:                  #For unbranched compounds
                                    
                                     if k == j:                         #For valency or charge
                                         
                                        if predicted_product[k][k] != 0 and predicted_product[k][k] != 1000 and predicted_product[k][k] != 2000:
                                            template += alpha[k] + "."+str(predicted_product[k][k])
                                        elif predicted_product[k][k] == 0:
                                            f1 = 0
                                            for t in range(ed_len):
                                                if predicted_product[k][t] == 1000 or predicted_product[k][t] == 2000:
                                                    f1 = 1
                                                    break
                                            if predicted_product[k][ed_len] > 0 and f1 == 0:
                                                template+= alpha[k] + "+"+str(predicted_product[k][ed_len])
                                            elif predicted_product[k][ed_len] < 0 and f1 == 0:
                                                template+= alpha[k] +str(predicted_product[k][ed_len])
                                                                                                                
                                                
                                     elif predicted_product[k][j] == 1:                         #For Single bond
                                        template = template +alpha[k]+"-"+alpha[j]
                                     elif predicted_product[k][j] == 2:                         #For Double bond
                                        template =template +alpha[k]+"="+alpha[j]
                                     elif predicted_product[k][j] == 3:                         #For Triple bond
                                        template = template+alpha[k]+"#"+alpha[j]
                                     elif predicted_product[k][j] == 100:                       #For coordinate bond
                                        template =  template+alpha[k]+"~"+alpha[j]
                                     elif predicted_product[k][j] == 200:                       #For coordinate bond
                                        template =  template+alpha[k]+"~/"+alpha[j]
                                     elif predicted_product[k][j] == 500:                       #For addition compound
                                        template = template+alpha[k]+"*"+alpha[j]
                                     elif predicted_product[k][j] == 1000:                      #For cations
                                        template=template+'['+alpha[k]+']'+'['+alpha[j]+']'
                                     elif predicted_product[k][j] == 2000:                      #For anions
                                        template=template+'['+alpha[j]+']'+'['+alpha[k]+']'    
                            
                           
                            else:                           #For visited elements in products
                                
                                l = len(template)
                                for t in range(l):
                                    if template[t]==alpha[k]:
                                         s1 = template[:t]
                                         s2 = template[t+1:]
                                         s3 = template[t]
                                         
                                         
                                             
                                         if brack_flag[k]==0:               #For unbranched compounds
                                             
                                             if k == j:                     #For valency or charge
                                                 
                                                if predicted_product[k][k] != 0 and predicted_product[k][k] != 1000 and predicted_product[k][k] != 2000:
                                                    s3 = alpha[k] + "."+str(predicted_product[k][k])
                                                elif predicted_product[k][k] == 0:
                                                    f1 = 0
                                                    for t in range(ed_len):
                                                        if predicted_product[k][t] == 1000 or predicted_product[k][t] == 2000:
                                                            f1 = 1
                                                            break
                                                    if predicted_product[k][ed_len] > 0 and f1 == 0:
                                                        s3 = alpha[k] + "+"+str(predicted_product[k][ed_len])
                                                    elif predicted_product[k][ed_len] < 0 and f1 == 0:
                                                        s3 = alpha[k] + str(predicted_product[k][ed_len])
                                                    
                                                                                                                                         
                    
                                             elif predicted_product[k][j] == 1:             #For Single bond
                                                s3 = s3+"-"+alpha[j]
                                             elif predicted_product[k][j] == 2:             #For Double bond
                                                s3 = s3+"="+alpha[j]
                                             elif predicted_product[k][j] == 3:             #For Triple bond
                                                s3 = s3+"#"+alpha[j]
                                             elif predicted_product[k][j] == 100:           #For coordinate bond
                                                s3 = s3+"~"+alpha[j]
                                             elif predicted_product[k][j] == 200:           #For coordinate bond
                                                s3 = s3+"~/"+alpha[j]
                                             elif predicted_product[k][j] == 500:           #For addition compound
                                                s3 = s3+"*"+alpha[j]
                                             elif predicted_product[k][j] == 1000:          #For cations
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
                    
                                             elif predicted_product[k][j] == 2000:                              #For anions
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
                    
                                             if predicted_product[k][j] != 1000 and predicted_product[k][j] != 2000:
                                                template = s1+s3+s2
                                             else:
                                                template = s3
                                        
                                            
                                            
                                         elif brack_flag[k]==1:                                 #For branched compounds
                                             
                                             if k == j:                                         #For valency or charge                       
                                                     
                                                if predicted_product[k][k] != 0 and predicted_product[k][k] != 1000 and predicted_product[k][k] != 2000:
                                                    s3 = alpha[k] + "." + str(predicted_product[k][k])
                                                elif predicted_product[k][k] == 0:
                                                    f1 = 0
                                                    for t in range(ed_len):
                                                        if predicted_product[k][t] == 1000 or predicted_product[k][t] == 2000:
                                                            f1 = 1
                                                            break
                                                    if predicted_product[k][ed_len] > 0 and f1 == 0:
                                                        s3 = alpha[k] + "+"+str(predicted_product[k][ed_len])
                                                    elif predicted_product[k][ed_len] < 0 and f1 == 0:
                                                        s3 = alpha[k] +str(predicted_product[k][ed_len])
                                                    
                                                                                                                                         
                    
                                             elif predicted_product[k][j] == 1:             #For Single Bond
                                                s3 = s3+"(-"+alpha[j]+")"
                                             elif predicted_product[k][j] == 2:             #For Double Bond
                                                s3 = s3+"(="+alpha[j]+")"
                                             elif predicted_product[k][j] == 3:             #For Triple Bond
                                                s3 = s3+"(#"+alpha[j]+")"
                                             elif predicted_product[k][j] == 100:           #For Coordinate Bond
                                                s3 = s3+"(~"+alpha[j]+")"
                                             elif predicted_product[k][j] == 200:           #For Coordinate Bond
                                                s3 = s3+"(~/"+alpha[j]+")"
                                             elif predicted_product[k][j] == 500:           #For addition compound
                                                s3 = s3+"*"+alpha[j]
                                             elif predicted_product[k][j] == 1000:          #For cations
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
                    
                                             elif predicted_product[k][j] == 2000:                      #For anions
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
                    
                                             if predicted_product[k][j] != 1000 and predicted_product[k][j] != 2000:
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
                            if prod_temp=='':
                                prod_temp="{"+template+"}"
                            else:
                                prod_temp=prod_temp+" + "+ "{" + template + "}"
                            template=''
                        break
                    
                    
                                 
                    e=0
                    k=stack_atom[st_top]        
                    
                    
                        
                         
        
            
                
                
        
        
        
        for i in range(len(alpha)):
            if alpha[i]!='' and prod_temp.find(alpha[i])==-1:
                alpha[i]=''
                
        final_temp=''
        for i in prod_temp:
            if i in alpha:
                final_temp+=ed_elements[alpha.index(i)]
            else:
                final_temp+=i
            
            
       # print("\n"+"{"+final_temp+"}")   
           
        
        brackets_info=[0 for _ in range(1000)]
        bl=0
        for i in range(len(final_temp)):                                                              #Removing redundant brackets from ESMILES notation of product
            if i+2<len(final_temp) and final_temp[i]==']' and final_temp[i+1]=='[' and final_temp[i+2].isalpha():
                bl=100
                    
                
            if i>0 and final_temp[i]=='(' and brackets_info[i]==0 and final_temp[i-1]!=')':
                bl=bl+1
                brackets_info[i]=bl
            
                
            if i+1<len(final_temp) and final_temp[i]==')' and final_temp[i+1]!='(':
                brackets_info[i]=bl
                if bl!=0 and bl!=100:
                    bl=bl-1
                
                
            if i+1<len(final_temp) and final_temp[i]==')' and final_temp[i+1]=='(':
                
                for j in range(0,i+1):
                    if brackets_info[j]==bl:
                        brackets_info[j]=0
                        
                
                if bl!=0 and bl!=100:
                    bl=bl-1
          
        rem_brac=[]                                             
        for i in range(len(brackets_info)):                 
            if brackets_info[i]!=0 and brackets_info[i]!=100:
                rem_brac.append(i)
        for i in range(len(rem_brac)):
            final_temp=final_temp[:rem_brac[i]]+final_temp[rem_brac[i]+1:]
            for n in range(len(rem_brac)):
                rem_brac[n]=rem_brac[n]-1
            
                
        product_esmiles.append(final_temp)          #ESMILES notation of product stored in a list
        
        

def predict(ebe_matrix_educts,ed_elements,cov,parts):           #Method to perform LCS to generate EBE matrix of product

    global max_ind
    for file in range(1,5):
        t="Templates_"+str(cov)+str(parts)+str(file)+".csv"        #Files contating templates
        if os.path.isfile(t):
            reader = open(t,"r")
            lines = reader.read().split("\n")
            reader.close()
            for x in lines:
                if x:            
                    s=list(ast.literal_eval(x))
                    max_e=len(s[1])
                    e=0

                    if len(ed_elements)<len(s[1]):                        
                        continue

                    temp_string=s[0].split('+')

                    cp=0
                    for string in temp_string:
                        if cp==1:
                            break
                        max_temp_ind=0
                        for chars in string:
                            if chars.isalpha():
                                max_temp_ind+=1
                        if max_temp_ind>max_ind:
                            cp=1
                            break

                    if cp==1:
                        continue
                    
                    

                    perm = permutations(s[6])                          #finding combinations matrix of all elements in template 

                    for ii in list(perm): 
                        permuted_educts=[[0 for x in range(len(s[1])+2)] for y in range(len(s[1]))] 
                        permuted_products=[[0 for x in range(len(s[1])+2)] for y in range(len(s[1]))] 
                        
                        kk=list(ii)  
                        p=0
                        q=0
                        for x in kk:
                            q=0
                            for y in kk:
                                permuted_educts[p][q]=s[1][x][y]
                                permuted_educts[p][len(s[1])]=s[1][x][len(s[1])]
                                permuted_educts[q][len(s[1])]=s[1][y][len(s[1])]
                                
                                               
                                permuted_products[p][q]=s[3][x][y]
                                permuted_products[p][len(s[1])]=s[3][x][len(s[1])]
                                permuted_products[q][len(s[1])]=s[3][y][len(s[1])]
                                
                                
                                q+=1
                            p+=1
                        
                        for x in range(len(s[1])):                                     #For placing proper values in coordinate bond column in educt matrix
                            for y in range(len(s[1])):
                                if permuted_educts[x][y] == 100:
                                        if permuted_educts[x][len(s[1])+1]!=0:
                                            permuted_educts[x][len(s[1])+1]=-99
                                        else:
                                            permuted_educts[x][len(s[1])+1] = y+1

                                        if permuted_educts[y][len(s[1])+1]!=0:
                                            permuted_educts[y][len(s[1])+1]=-99
                                        else:
                                            permuted_educts[y][len(s[1])+1] = x+1
                            
                        for x in range(len(s[1])):                                   #For placing proper values in coordinate bond column in educt matrix
                            for y in range(len(s[1])):
                                if permuted_products[x][y] == 100:
                                        if permuted_products[x][len(s[1])+1]!=0:
                                            permuted_products[x][len(s[1])+1]=-99
                                        else:
                                            permuted_products[x][len(s[1])+1] = y+1

                                        if permuted_products[y][len(s[1])+1]!=0:
                                            permuted_products[y][len(s[1])+1]=-99
                                        else:
                                           permuted_products[y][len(s[1])+1] = x+1
                        
                                         



                        
                        short_string=''
                        v=permuted_educts[0]
                        for p in v[:-2]:                           #Converting integer values into string for LCS algorithm of template matrix
                            if p==0:
                                short_string=short_string+'a'
                            elif p==1:
                                short_string=short_string+'b'
                            elif p==2:
                                short_string=short_string+'c'
                            elif p==3:
                                short_string=short_string+'d'
                            elif p==100:
                                short_string=short_string+'e'
                            elif p==200:
                                short_string=short_string+'f'
                            elif p==500:
                                short_string=short_string+'g'
                            elif p==1000:
                                short_string=short_string+'h'
                            elif p==2000:
                                short_string=short_string+'i'
                            else:
                                short_string=short_string+'j'
                            
                        e=0        
                        for k in ebe_matrix_educts:
                            if len(ed_elements)-e<max_e:
                                break
                            
                            long_string=''
                            
                            for p in k[e:len(k)-2]:                         #Converting integer values into string for LCS algorithm of educt matrix
                                if p==0:
                                    long_string=long_string+'a'
                                elif p==1:
                                    long_string=long_string+'b'
                                elif p==2:
                                    long_string=long_string+'c'
                                elif p==3:
                                    long_string=long_string+'d'
                                elif p==100:
                                    long_string=long_string+'e'
                                elif p==200:
                                    long_string=long_string+'f'
                                elif p==500:
                                    long_string=long_string+'g'
                                elif p==1000:
                                    long_string=long_string+'h'
                                elif p==2000:
                                    long_string=long_string+'i'
                                else:
                                    long_string=long_string+'j'
                                
                            ls=len(long_string)
                            ss=len(short_string)
                            
                           
                            
                            #print(long_string,short_string)
                            c=[[0 for _ in range(ls+1)] for _ in range(ss+1)]           #Matrix for LCS algorithm
                            
                            for i in range(1,ss+1):                                     #LCS algorithm
                                for j in range(1,ls+1):
                                    if short_string[i-1]==long_string[j-1]:
                                        c[i][j]=c[i-1][j-1]+1
                                    elif c[i-1][j]>=c[i][j-1]:
                                        c[i][j]=c[i-1][j]             
                                    elif c[i-1][j]<c[i][j-1]:
                                        c[i][j]=c[i][j-1]     
                                        
                           # print(c[ss][ls])
                                        
                            if c[ss][ls]==ss and ss<=ls:                              #Extracting rows/columns form educt matrix after LCS algorithm

                                #print(k,long_string,v,short_string)
                                stack=[0 for _ in range(1000)]
                                top=-1
                                i=0
                                j=0
                                matched_temp=[]                                      #List containing extracted rows/columns of Educt matrix
                                while True:
                                    if j>=ls:		
                                        if top==-1:
                                            break
                                        else:
                                            i-=1
                                            j=stack[top]+1                   
                                            top-=1
                                                                       
                                            if j>=ls:
                                                j=stack[top]+1
                                                top=top-1
                                                i=i-1
                                                if i<0:
                                                    break
                                                continue
                                    
                                    #print(short_string,i)
                                    #print(long_string,j)
                                    if short_string[i]==long_string[j]:
                                        top+=1
                                        stack[top]=j
                                        i+=1
                                        j+=1
                                        if top==c[ss][ls]-1:
                                            j=stack[top]+1
                                            i-=1
                                            h=[x + e for x in stack[0:top+1]]
                                            matched_temp.append(h)                                  #Inserting extracted rows/columns of Educt matrix into the list
                                            top-=1
                                    else:
                                        j+=1
                            
                                
                                    
                                match(permuted_educts,ebe_matrix_educts,permuted_products,matched_temp,s[4],s[5],list(ii))      #Mathod call to match elements of extracted rows/columns with the template elements
                                
                            e=e+1 

def go_press():  # Method to process educts given by the user

    global ed, prod, ed_elements,prod_elements,total_ed_matrix,total_prod_matrix, ed_compounds,prod_compounds,cov,bnd,parts,filename,count_ed,max_ind
    cov='0'
    bnd='1'
    ed = educt.get()
    f1 = 0
    f2 = 0
    if not ed:
        showwarning("Warning", "Please Fill the Educt Part")

    else:                                           #Extracting educt compounds
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

        ed_elements=[] 



        for i in ed:                 # Loop to process the educt string input given by the user
            s=i[1:len(i)-1]
            c=0
            
            for j in range(len(s)):

                if (j+1)<len(s) and s[j].isupper() and s[j+1].islower():
                    c+=1
                elif s[j].isupper():
                    c+=1
            if c>max_ind:
                max_ind=c     
            

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
                    if bnd=='':
                        bnd='1'


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
                    if bnd=='' or bnd=='1':
                        bnd='2'
                              
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
                    if bnd=='' or bnd=='1' or bnd=='2':
                        bnd='3'
                              
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
                    if bnd=='' or bnd=='1' or bnd=='2' or bnd=='3':
                        bnd='4'

                elif s[k]=='~' and flag!=-1:                                                #Identifying coordinate bond                                 
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
                            ebe_matrix_educts[x][x] = (
                                ebe_matrix_educts[x][x]*-1) % 10
    
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
        
        count_ed=len(ebe_matrix_educts)
        
        #print(ed_elements)
        #print(ebe_matrix_educts)
        prod=''
        
        f="Reactions_"+str(cov)+str(parts)+str(bnd)+".csv"   # filename to serach reaction
        if os.path.isfile(f):
            reader = open(f, "r")
            lines = reader.read().split("\n")
            reader.close()
            for i in lines:
                if i:            
                    s=list(ast.literal_eval(i))
                    #print(s)
                    if s[1]==ed_elements and s[2]==ebe_matrix_educts:
                        
                        prod=s[3]
                        break
                    
        if prod=='':
                     
            predict(ebe_matrix_educts,ed_elements,cov,parts)            #Predicting products 
            
            print_products()                                            #Printing the products
            
            
        else:
            stable_products([prod])
            
            
 

window = Tk()       #First window for the user
ed = ''
prod = ''
educt = StringVar()         #User input for educts
ed_types = []
ed_compounds = []
educt_matrix = []
total_ed_matrix = []
ed_elements = []
ed_element_slno = []
wrng_ed = 0
v = StringVar()
educt_element = []
product_esmiles=[]              #List to store the esmiles notation of predicted products

filename="Reactions_"           #
cov=''                          #
parts=''                        #
bnd=''                          # to identify the filaname to sore reaction and templates

count_ed=0
max_ind=0

# First Window Screen that appears for the User
w, h = window.winfo_screenwidth(), window.winfo_screenheight()
window.geometry("%dx%d+0+0" % (w, h))
window.title("PREDICTION OF REACTION")

Label(window, text="Welcome", font=("Arial Bold", 50), fg = "blue", bg = "white").pack()
Label(window, text="Please Enter the Educt in the text box in ESMILES Notation:", font=("Arial Bold", 15),fg = "dark green", bg = "white").place(x=40, y=140)       #
                                                                                    
Label(window, text="Please click on the button to read the Manual and learn how to input chemical compounds and elements in ESMILES notation ", font=("Arial Bold", 10),fg = "dark green", bg = "white").place(x=60, y=200)                       #

b = Button(window, text='Read Manual', command=read_manual,fg = "blue", bg = "white").place(x=900, y=190)           

window.config(bg="white")


Label(window, text="Educt",font=("Arial Bold", 15), fg = "black", bg = "white").place(x=40, y=290)


E1 = Entry(window, textvariable=educt,width=70,font=("Arial Bold", 15))  # Textbox for Reading Educt Compounds
E1.place(x=140, y=290)





btn = Button(window, text='Predict', command=go_press, fg = "dark green", bg = "white").place(x=480, y=360)  # Calling Fuctions to predict on button press



# Reading Elemenets name and Valency from file
ele = pandas.read_csv("elementlist.csv", header=None)

elements = ele[1].values.tolist()
elements_name=ele[2].values.tolist()
elements_valency = ele[3].values.tolist()
elements_free_electons = ele[4].values.tolist()
elements_stable_electons = ele[5].values.tolist()
elements_max_bond = ele[6].values.tolist()
elements_eneg = ele[7].values.tolist()




window.mainloop()

    
