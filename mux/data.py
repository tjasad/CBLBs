DOX_max = 4
PRO_max = 3
ALD_max = 5
aCA_max = 6
EST_max = 4
DEX_max = 4

#####
min_gamma = 0
max_gamma = 10000

min_alpha1 = 0
max_alpha1 = 100

min_alpha2 = -3
max_alpha2 = 1

min_omega1 = 1
max_omega1 = 10

min_omega2 = -15
max_omega2 = 15

min_n = 1
max_n = 10

#####


names = ('DOX', 'PRO', 'ALD', 'aCa', 'EST', 'DEX')
max_vals = (DOX_max, PRO_max, ALD_max, aCA_max, EST_max, DEX_max)

params = {"IL1": (65.3,0.26,5.66*10**(-11),3.7), #ID_DOX
          "IL2": (10,8.7,2.5*10**(-3),1.7), # ID_PRO
          "IL3": (97,0,8.2*10**(-5),3.5), # ID_ALD
          "IL4": (4.1,21.7,1.5*10**(-2),1.5), # ID_aCa
          "IL5": (93.5,0.0107,8.01*10**(-11),2.9), # ID_EST
          "IL6": (20.2,4.4,3.5*10**(-7),1.98), #ID_DEX
          "IL7": (86.5,0.208,1.52*10**(-6),1.4), # NOT_DOX
          "IL8": (9,9.4,4*10**(-10),2.9), # NOT_PRO
          "IL9": (89.7,0.0245,2.89*10**(-9),7.9), # NOT_ALD
          "IL10": (5.9,14.1,3.6*10**(-6),3.1), # NOT_aCa
          "IL11": (95,0,1.62*10**(-5),3.1), # NOT_EST
          "IL12": (2.1,36.6,1*10**(-3),1.45), # NOT_DEX
          "OL1_GFP": (98.01,0.0102,1.5*10**(-14),4.3), # NOT_GFP
          "OL1_GFP_au": (95,0.0105,4.8*10**(-14),4.3), # 
          "OL2": (95,0.105,3.39*10**(-11),2.65), # NOT_mCherry
          "OL3": (85,0.1335,1*10**(-5),1.6), # NOT_alphaCa
          "BL": (4,24.5,1.5*10**(-38),10), # produce GFP in the presence of the aCa
          "IL_ID": (4.1,21.7,1.5*10**(-2),1.5), # ID_aCa
          "IL_NOT": (95,0,1.62*10**(-5),3.1), # NOT_EST
          "OL_NOT": (95,0,1.62*10**(-5),3.1)} # NOT_EST
          


"""
cells = {"ID_DOX":"IL1",
         "ID_PRO":"IL2",
         "ID_ALD":"IL3",
         "ID_aCa":"IL4",
         "ID_EST":"IL5",
         "ID_DEX":"IL6",
         "NOT_DOX":"IL7",
         "NOT_PRO":"IL8",
         "NOT_ALD":"IL9",
         "NOT_aCa":"IL10",
         "NOT_EST":"IL11",
         "NOT_DEX":"IL12",
         "NOT_GFP":"OL1_GFP",
         "NOT_GFP_au":"OL1_GFP_au",
         "NOT_mCherry":"OL2",
         "NOT_alphaCa":"OL3",
         "BL":"BL"}
"""


cells = {"ID_DOX":"IL1",
         "NOT_DOX":"IL2",
         "ID_PRO":"IL3",
         "NOT_PRO":"IL4",
         "ID_ALD":"IL5",
         "NOT_ALD":"IL6",
         "ID_aCa":"IL7",
         "NOT_aCa":"IL8",
         "ID_EST":"IL9",
         "NOT_EST":"IL10",
         "ID_DEX":"IL11",
         "NOT_DEX":"IL12",
         "NOT_GFP":"OL1_GFP",
         "NOT_GFP_au":"OL1_GFP_au",
         "NOT_mCherry":"OL2",
         "NOT_alphaCa":"OL3",
         "BL":"BL"}

"""
cells = {"ID_DOX":"IL2",
         "NOT_DOX":"IL1",
         "ID_PRO":"IL4",
         "NOT_PRO":"IL3",
         "ID_ALD":"IL6",
         "NOT_ALD":"IL5",
         "ID_aCa":"IL8",
         "NOT_aCa":"IL7",
         "ID_EST":"IL10",
         "NOT_EST":"IL9",
         "ID_DEX":"IL12",
         "NOT_DEX":"IL11",
         "NOT_GFP":"OL1_GFP",
         "NOT_GFP_au":"OL1_GFP_au",
         "NOT_mCherry":"OL2",
         "NOT_alphaCa":"OL3",
         "BL":"BL"}
"""

scales = {"DOX": 10**4,
          "PRO": 10**3,
          "ALD": 10**5,
          "aCa": 10**6,
          "EST": 10**3,
          "DEX": 10**4}
