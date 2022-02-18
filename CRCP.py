from tkinter import *                      # Importing tkinter for designing the UI
from tkinter.messagebox import *           # Importing tkinter.messagebox for warning/cofirmation message
import os                                  # Importing os to run system command
import pandas                              # Importing pandas to for datastructure
import csv                                 # Importing csv to read/write csv file
from itertools import permutations         # Importing permutation for permuting list
import subprocess
import webbrowser

def read_manual():
    webbrowser.open_new('Manual_CRCP.pdf')

def enter_reaction():
    os.system("python Esmiles_Input.py")
    
def delete_reaction():
    os.system("python Delete_Reaction.py")
    
def predict_reaction():
    os.system("python Predict.py")
    
def run(c):                                # Function to Run New Program
    c.destroy()
    window.destroy()
    
    

    
window = Tk()                   #First window for the user  
w, h = window.winfo_screenwidth(), window.winfo_screenheight()

window.geometry("%dx%d+0+0" % (w, h))
window.title("CRCP SOFTWARE")           #Title of the Window Screen

window.configure(bg="white")

Label(window, text="Welcome To Chemical Reaction Classifier & Predictor Software", font=("Arial Bold", 30), fg = "blue", bg = "white").pack()       

Label(window, text="Please Enter Your Choice : ",fg = "black", bg = "white", font=("Arial", 20)).place(x=100, y=260)  

btn1 = Button(window, text='Enter New Reaction', command=enter_reaction,fg = "dark green", bg = "white").place(x=520, y=260)           # Calling function after button click

btn2 = Button(window, text='Delete Reaction', command=delete_reaction,fg = "dark green", bg = "white").place(x=720, y=260)           # Calling function after button click

btn3 = Button(window, text='Predict Reaction', command=predict_reaction,fg = "dark green", bg = "white").place(x=900, y=260)           # Calling function after button click

btn4 = Button(window, text='Exit', command=window.destroy,fg = "dark green", bg = "white").place(x=1090, y=260)           # Calling function after button click

Label(window, text="Please click on the button to read the Manual and learn how to use CRCP software ", font=("Arial Bold", 10), bg = "white").place(x=100, y=460)                       #

b = Button(window, text='Read Manual', command=read_manual,fg = "blue", bg = "white").place(x=950, y=455)  

window.mainloop()                                                               #End of window
