# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 20:50:04 2018

@author: faure
"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess
from math import *
#from colors import *

#parametre 
K = 1
x, y = 0.5,0.5 # point
Lcol = ['red', 'blue', 'green', 'black', 'yellow']

col = 0 # numero de couleur

#--definit l'application sur R2
def f(x,y):
    y2 = y +K*sin(x)
    x2 = x + y2
    y2 = y2 % (2*pi) # ramene dans [0,2 pi]
    x2 = x2 % (2*pi) # ramene dans [0,2 pi]
    return x2, y2
    
#---- calcule et dessin de la trajectoire
# entree: x,y,tmax (duree)
def dessin_trajectoire(x,y, tmax):
    Lx,Ly = [x],[y] #listes
    for t in range(tmax):
        x,y = f(x,y) # iteration
        Lx.append(x)
        Ly.append(y)
    plt.plot(Lx,Ly, linestyle='none', marker=',', color = Lcol[col])
   # plt.show()
    plt.pause(0.001) # montre le dessin sans bloquer
    return x,y
            
    

#---- si evenement souris (clik)
def on_click(event):
    print('click')
    global x,y, col
    x, y = event.xdata, event.ydata # coord souris en x,y
    #x, y = event.x, event.y # coord souris en pixels
    col = (col+1) % len(Lcol)
    if event.button == 1:
            tmax = 10000 # duree 
            x,y = dessin_trajectoire(x,y, tmax)
            
            
def on_key(event):
    print('key = ', event.key)
    global x,y
   # print('you pressed', event.key, event.xdata, event.ydata)
    if event.key == ' ': # barre espace
            tmax = 10000 # duree 
            x,y = dessin_trajectoire(x,y, tmax)



plt.axis([0,2*pi,0,2*pi]) # selectionne la vue
#plt.axis('equal') # pour avoir meme echelle en x,y
plt.xlabel('x1')
plt.ylabel('x2')
plt.title("Application standard, K=%2.f"%K +  ' (cliquer pour condition initiale, et barre espace)')
plt.connect('button_press_event', on_click) # associe la fonction aux evenements  souris
plt.connect('key_press_event', on_key)  # associe la fonction aux evenements clavier


plt.show()
    