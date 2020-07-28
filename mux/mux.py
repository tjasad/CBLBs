import numpy as np
import pickle

from transfer import *
from data import *


with open('direct_params.pickle', 'rb') as handle:
    direct_params = pickle.load(handle)

def mux_macia(S, I, file = ""):
    ALD, aCa, EST, DEX = I #I_0, I_1, I_2, I_3
    PRO, DOX = S #S_0, S_1

    ALD = ALD*scales["ALD"]
    aCa = aCa*scales["aCa"]
    EST = EST*scales["EST"]
    DOX = DOX*scales["DOX"]
    PRO = PRO*scales["PRO"]
    DEX = DEX*scales["DEX"]

    """
    Psi_1 (I_0) chamber 1
    """
    
    # ID_PRO
    ch11 = T_f(PRO, *direct_params["ID_PRO"])

    # ID_DOX
    ch12 = T_f(DOX, *direct_params["ID_DOX"])

    # NOT_ALD
    ch13 = T_f(ALD, *direct_params["NOT_ALD"])



    # NOT_GFP
    ch1 = max([ch11, ch12, ch13])
    GFP1 =  T_f(ch1, *params[cells["NOT_GFP"]])

    """
    Psi_2 (I_1) chamber 2
    """
    
    # ID_PRO
    ch21 = T_f(PRO, *direct_params["ID_PRO"])

    # NOT_DOX
    ch22 = T_f(DOX, *direct_params["NOT_DOX"])

    # NOT_aCa
    ch23 = T_f(aCa, *direct_params["NOT_aCa"])

    # NOT_GFP
    ch2 = max([ch21, ch22, ch23])
    GFP2 =  T_f(ch2, *params[cells["NOT_GFP"]])

    """
    Psi_3 (I_2) chamber 3
    """
    
    # NOT_PRO
    ch31 = T_f(PRO, *direct_params["NOT_PRO"])

    # ID_DOX
    ch32 = T_f(DOX, *direct_params["ID_DOX"])

    # NOT_EST
    ch33 = T_f(EST, *direct_params["NOT_EST"])

    # NOT_GFP
    ch3 = max([ch31, ch32, ch33])
    GFP3 =  T_f(ch3, *params[cells["NOT_GFP"]])


    """
    Psi_4 (I_3) chamber 4
    """

    # NOT_PRO
    ch41 = T_f(PRO, *direct_params["NOT_PRO"])

    # NOT_DOX
    ch42 = T_f(DOX, *direct_params["NOT_DOX"])

    # NOT_DEX
    ch43 = T_f(DEX, *direct_params["NOT_DEX"])

    # NOT_GFP
    ch4 = max([ch41, ch42, ch43])
    GFP4 =  T_f(ch4, *params[cells["NOT_GFP"]])

    if not file:

        """
        output
        """
        print("******")
        print("S:", S)
        print("I:", I)
        print(f"ch1:{ch1:.2f}, ch2:{ch2:.2f}, ch3:{ch3:.2f}, ch4:{ch4:.2f}")
        print(f"GFP1:{GFP1:.2f}, GFP2:{GFP2:.2f}, GFP3:{GFP3:.2f}, GFP4:{GFP4:.2f}")
        print("******")
    else:
        f = open(file, "a", encoding="utf8")
        for s in S:
            print(s, end=",",file=f)
        

        for i in I:
            print(i, end=",", file=f)
        

        print(f"{ch1:.2f},{ch2:.2f},{ch3:.2f},{ch4:.2f},", file=f, end="")
        print(f"{GFP1:.2f},{GFP2:.2f},{GFP3:.2f},{GFP4:.2f}\n", file=f, end="")

        f.close()


def mux_general(S, I, file = ""):
    S = 1000 * np.array(S)
    I = 1000 * np.array(I)
    
    I_0, I_1, I_2, I_3 = I
    S_0, S_1 = S

    """
    Psi_1 (I_0) chamber 1
    S_0 or S_1 or not I_0 
    """
    
    # S_0
    ch11 = T_f(S_0, *params["IL_ID"])

    # S_1
    ch12 = T_f(S_1, *params["IL_ID"])

    # not I_0
    ch13 = T_f(I_0, *params["IL_NOT"])

    # NOT_GFP
    ch1 = max(ch11, ch12, ch13)
    GFP1 =  T_f(ch1, *params["OL_NOT"])

    """
    Psi_2 (I_1) chamber 2
    S_0 or not S_1 or not I_1
    """
    
    # S_0
    ch21 = T_f(S_0, *params["IL_ID"])

    # NOT S_1
    ch22 = T_f(S_1, *params["IL_NOT"])

    # NOT I_1
    ch23 = T_f(I_1, *params["IL_NOT"])

    # NOT_GFP
    ch2 = max([ch21, ch22, ch23])
    GFP2 =  T_f(ch2, *params["OL_NOT"])

    """
    Psi_3 (I_2) chamber 3
    not S_0 or S_1 or not I_2
    """
    
    # NOT S_0
    ch13 = T_f(S_0, *params["IL_NOT"])

    # S_1
    ch23 = T_f(S_1, *params["IL_ID"])

    # NOT I_2
    ch33 = T_f(I_2, *params["IL_NOT"])

    # NOT_GFP
    ch3 = max([ch13, ch23, ch33])
    GFP3 =  T_f(ch3, *params["OL_NOT"])


    """
    Psi_4 (I_3) chamber 4
    not S_0 or not S_1 or not I_3
    """

    # NOT S_0
    ch41 = T_f(S_0, *params["IL_NOT"])

    # NOT S_1
    ch42 = T_f(S_1, *params["IL_NOT"])

    # NOT I_3
    ch43 = T_f(I_3, *params["IL_NOT"])

    # NOT_GFP
    ch4 = max([ch41, ch42, ch43])
    GFP4 =  T_f(ch4, *params["OL_NOT"])

    if not file:

        """
        output
        """
        print("******")
        print("S:", S)
        print("I:", I)
        print(f"ch1:{ch1:.2f}, ch2:{ch2:.2f}, ch3:{ch3:.2f}, ch4:{ch4:.2f}")
        print(f"GFP1:{GFP1:.2f}, GFP2:{GFP2:.2f}, GFP3:{GFP3:.2f}, GFP4:{GFP4:.2f}")
        print("******")
    else:
        f = open(file, "a", encoding="utf8")
        for s in S:
            print(s, end=",",file=f)
        

        for i in I:
            print(i, end=",", file=f)
        

        print(f"{ch1:.2f},{ch2:.2f},{ch3:.2f},{ch4:.2f},", file=f, end="")
        print(f"{GFP1:.2f},{GFP2:.2f},{GFP3:.2f},{GFP4:.2f}\n", file=f, end="")

        f.close()


if __name__ == '__main__':

    get_bin = lambda x, n: tuple(map(int, format(x, 'b').zfill(n)))

    I_0, I_1, I_2, I_3 = 0,0,0,0
    S_0, S_1 = 0,0

    S = (S_0, S_1)
    I = (I_0, I_1, I_2, I_3)

    f = open('out_macia.csv', 'w')
    f.write("S_0,S_1,I_0,I_1,I_2,I_3,ch1,ch2,ch3,ch4,GFP1,GFP2,GFP3,GFP4\n")
    f.close()

    f = open('out_general.csv', 'w')
    f.write("S_0,S_1,I_0,I_1,I_2,I_3,ch1,ch2,ch3,ch4,GFP1,GFP2,GFP3,GFP4\n")
    f.close()

    for S in range(4):
        s = get_bin(S,2)
        
        for I in range(16):
            i = get_bin(I,4)
        
        
            mux_macia(s,i,'out_macia.csv')
            mux_general(s,i,'out_general.csv')


    
