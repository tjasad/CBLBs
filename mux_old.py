import numpy as np
from transfer import *
from data import *




def mux_macia(S, I, file = ""):
    ALD, aCa, EST, DOX = I #I_0, I_1, I_2, I_3
    PRO, DEX = S #S_0, S_1

    ALD = ALD*scales["ALD"]
    ACa = aCa*scales["aCa"]
    EST = EST*scales["EST"]
    DOX = DOX*scales["DOX"]
    PRO = PRO*scales["PRO"]
    DEX = DEX*scales["DEX"]

    """
    Psi_1 (I_0) chamber 1
    """
    ch1 = 0

    # ID_PRO
    ch1 += T_f(PRO, *params[cells["ID_PRO"]])

    # ID_DOX
    ch1 += T_f(DOX, *params[cells["ID_DOX"]])

    # NOT_ALD
    ch1 += T_f(ALD, *params[cells["NOT_ALD"]])

    # NOT_GFP
    GFP1 =  T_f(ch1*50, *params[cells["NOT_GFP"]])

    """
    Psi_2 (I_1) chamber 2
    """
    ch2 = 0

    # ID_PRO
    ch2 += T_f(PRO, *params[cells["ID_PRO"]])

    # NOT_DOX
    ch2 += T_f(DOX, *params[cells["NOT_DOX"]])

    # NOT_aCa
    ch2 += T_f(aCa, *params[cells["NOT_aCa"]])

    # NOT_GFP
    GFP2 =  T_f(ch2, *params[cells["NOT_GFP"]])

    """
    Psi_3 (I_2) chamber 3
    """
    ch3 = 0

    # NOT_PRO
    ch3 += T_f(PRO, *params[cells["NOT_PRO"]])

    # ID_DOX
    ch3 += T_f(DOX, *params[cells["ID_DOX"]])

    # NOT_EST
    ch3 += T_f(EST, *params[cells["NOT_EST"]])

    # NOT_GFP
    GFP3 =  T_f(ch3, *params[cells["NOT_GFP"]])


    """
    Psi_4 (I_3) chamber 4
    """

    ch4 = 0

    # NOT_PRO
    ch4 += T_f(PRO, *params[cells["NOT_PRO"]])

    # NOT_DOX
    ch4 += T_f(DOX, *params[cells["NOT_DOX"]])

    # NOT_DEX
    ch4 += T_f(DEX, *params[cells["NOT_DEX"]])

    # NOT_GFP
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
            print(s, end="\t",file=f)
        print(end="\t", file=f)

        for i in I:
            print(i, end="\t", file=f)
        print(end="\t", file=f)

        print(f"{ch1:.2f}\t{ch2:.2f}\t{ch3:.2f}\t{ch4:.2f}\t", file=f, end="")
        print(f"{GFP1:.2f}\t{GFP2:.2f}\t{GFP3:.2f}\t{GFP4:.2f}\n", file=f, end="")

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
    ch1 = 0

    # S_0
    ch1 += T_f(S_0, *params["IL_ID"])

    # S_1
    ch1 += T_f(S_1, *params["IL_ID"])

    # not I_0
    ch1 += T_f(I_0, *params["IL_NOT"])

    # NOT_GFP
    GFP1 =  T_f(ch1, *params["OL_NOT"])

    """
    Psi_2 (I_1) chamber 2
    S_0 or not S_1 or not I_1
    """
    ch2 = 0

    # S_0
    ch2 += T_f(S_0, *params["IL_ID"])

    # NOT S_1
    ch2 += T_f(S_1, *params["IL_NOT"])

    # NOT I_1
    ch2 += T_f(I_1, *params["IL_NOT"])

    # NOT_GFP
    GFP2 =  T_f(ch2, *params["OL_NOT"])

    """
    Psi_3 (I_2) chamber 3
    not S_0 or S_1 or not I_2
    """
    ch3 = 0

    # NOT S_0
    ch3 += T_f(S_0, *params["IL_NOT"])

    # S_1
    ch3 += T_f(S_1, *params["IL_ID"])

    # NOT I_2
    ch3 += T_f(I_2, *params["IL_NOT"])

    # NOT_GFP
    GFP3 =  T_f(ch3, *params["OL_NOT"])


    """
    Psi_4 (I_3) chamber 4
    not S_0 or not S_1 or not I_3
    """

    ch4 = 0

    # NOT S_0
    ch4 += T_f(S_0, *params["IL_NOT"])

    # NOT S_1
    ch4 += T_f(S_1, *params["IL_NOT"])

    # NOT I_3
    ch4 += T_f(I_3, *params["IL_NOT"])

    # NOT_GFP
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
            print(s, end="\t",file=f)
        print(end="\t", file=f)

        for i in I:
            print(i, end="\t", file=f)
        print(end="\t", file=f)

        print(f"{ch1:.2f}\t{ch2:.2f}\t{ch3:.2f}\t{ch4:.2f}\t", file=f, end="")
        print(f"{GFP1:.2f}\t{GFP2:.2f}\t{GFP3:.2f}\t{GFP4:.2f}\n", file=f, end="")

        f.close()


if __name__ == '__main__':

    get_bin = lambda x, n: tuple(map(int, format(x, 'b').zfill(n)))

    I_0, I_1, I_2, I_3 = 0,0,0,0
    S_0, S_1 = 0,0

    S = (S_0, S_1)
    I = (I_0, I_1, I_2, I_3)

    f = open('out_macia.csv', 'w')
    f.write("S_0\tS_1\tI_0\tI_1\tI_2\tI_3\tch1\tch2\tch3\tch4\tGFP1\tGFP2\tGFP3\tGFP4\n")
    f.close()

    f = open('out_general.csv', 'w')
    f.write("S_0\tS_1\tI_0\tI_1\tI_2\tI_3\tch1\tch2\tch3\tch4\tGFP1\tGFP2\tGFP3\tGFP4\n")
    f.close()

    for S in range(4):
        s = get_bin(S,2)
        
        for I in range(16):
            i = get_bin(I,4)
        
        
            mux_macia(s,i,'out_macia.csv')
            mux_general(s,i,'out_general.csv')


    
