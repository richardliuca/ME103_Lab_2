from numpy import *

A = array([0,1,2,3,4,5])

def change():
    global A
    A+=1


def test(a):
    change()
    print(a)

test(A)
