\* Layout_Problem_Model *\
Minimize
OBJ: SumCD + SumOmega + TLC
Subject To
_C1: 0 O_furnace + l_furnace = 1
_C10: d_flash + l_flash = 4.6
_C100: Wy_comp_furnace + Wy_furnace_comp = 1
_C101: CD_reactor_comp = 0
_C102: - L_reactor_comp + R_reactor_comp + X_comp - X_reactor = 0
_C103: A_reactor_comp - B_reactor_comp + Y_comp - Y_reactor = 0
_C104: R_reactor_comp - 1000 Wx_reactor_comp <= 0
_C105: L_reactor_comp + 1000 Wx_reactor_comp <= 1000
_C106: A_reactor_comp - 1000 Wy_reactor_comp <= 0
_C107: B_reactor_comp + 1000 Wy_reactor_comp <= 1000
_C108: - A_reactor_comp - B_reactor_comp + D_reactor_comp - L_reactor_comp
 - R_reactor_comp = 0
_C109: - D_comp_reactor + D_reactor_comp = 0
_C11: X_flash - 0.5 l_flash >= 0
_C110: 1000 E1_reactor_comp + 1000 E2_reactor_comp - X_comp + X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= 0
_C111: - 1000 E1_reactor_comp + 1000 E2_reactor_comp + X_comp - X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= -1000
_C112: 1000 E1_reactor_comp - 1000 E2_reactor_comp - Y_comp + Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -1000
_C113: - 1000 E1_reactor_comp - 1000 E2_reactor_comp + Y_comp - Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -2000
_C114: L_reactor_comp - R_comp_reactor = 0
_C115: - L_comp_reactor + R_reactor_comp = 0
_C116: A_reactor_comp - B_comp_reactor = 0
_C117: - A_comp_reactor + B_reactor_comp = 0
_C118: Wx_comp_reactor + Wx_reactor_comp = 1
_C119: Wy_comp_reactor + Wy_reactor_comp = 1
_C12: Y_flash - 0.5 d_flash >= 0
_C120: CD_flash_comp - 254.464800781 D_flash_comp = 0
_C121: - L_flash_comp + R_flash_comp + X_comp - X_flash = 0
_C122: A_flash_comp - B_flash_comp + Y_comp - Y_flash = 0
_C123: R_flash_comp - 1000 Wx_flash_comp <= 0
_C124: L_flash_comp + 1000 Wx_flash_comp <= 1000
_C125: A_flash_comp - 1000 Wy_flash_comp <= 0
_C126: B_flash_comp + 1000 Wy_flash_comp <= 1000
_C127: - A_flash_comp - B_flash_comp + D_flash_comp - L_flash_comp
 - R_flash_comp = 0
_C128: - D_comp_flash + D_flash_comp = 0
_C129: 1000 E1_flash_comp + 1000 E2_flash_comp - X_comp + X_flash - 0.5 l_comp
 - 0.5 l_flash >= 0
_C13: 0 O_comp + l_comp = 12
_C130: - 1000 E1_flash_comp + 1000 E2_flash_comp + X_comp - X_flash
 - 0.5 l_comp - 0.5 l_flash >= -1000
_C131: 1000 E1_flash_comp - 1000 E2_flash_comp - Y_comp + Y_flash - 0.5 d_comp
 - 0.5 d_flash >= -1000
_C132: - 1000 E1_flash_comp - 1000 E2_flash_comp + Y_comp - Y_flash
 - 0.5 d_comp - 0.5 d_flash >= -2000
_C133: L_flash_comp - R_comp_flash = 0
_C134: - L_comp_flash + R_flash_comp = 0
_C135: A_flash_comp - B_comp_flash = 0
_C136: - A_comp_flash + B_flash_comp = 0
_C137: Wx_comp_flash + Wx_flash_comp = 1
_C138: Wy_comp_flash + Wy_flash_comp = 1
_C139: CD_furnace_distil = 0
_C14: d_comp + l_comp = 24
_C140: - L_furnace_distil + R_furnace_distil + X_distil - X_furnace = 0
_C141: A_furnace_distil - B_furnace_distil + Y_distil - Y_furnace = 0
_C142: R_furnace_distil - 1000 Wx_furnace_distil <= 0
_C143: L_furnace_distil + 1000 Wx_furnace_distil <= 1000
_C144: A_furnace_distil - 1000 Wy_furnace_distil <= 0
_C145: B_furnace_distil + 1000 Wy_furnace_distil <= 1000
_C146: - A_furnace_distil - B_furnace_distil + D_furnace_distil
 - L_furnace_distil - R_furnace_distil = 0
_C147: - D_distil_furnace + D_furnace_distil = 0
_C148: 1000 E1_furnace_distil + 1000 E2_furnace_distil - X_distil + X_furnace
 - 0.5 l_distil - 0.5 l_furnace >= 0
_C149: - 1000 E1_furnace_distil + 1000 E2_furnace_distil + X_distil
 - X_furnace - 0.5 l_distil - 0.5 l_furnace >= -1000
_C15: X_comp - 0.5 l_comp >= 0
_C150: 1000 E1_furnace_distil - 1000 E2_furnace_distil - Y_distil + Y_furnace
 - 0.5 d_distil - 0.5 d_furnace >= -1000
_C151: - 1000 E1_furnace_distil - 1000 E2_furnace_distil + Y_distil
 - Y_furnace - 0.5 d_distil - 0.5 d_furnace >= -2000
_C152: L_furnace_distil - R_distil_furnace = 0
_C153: - L_distil_furnace + R_furnace_distil = 0
_C154: A_furnace_distil - B_distil_furnace = 0
_C155: - A_distil_furnace + B_furnace_distil = 0
_C156: Wx_distil_furnace + Wx_furnace_distil = 1
_C157: Wy_distil_furnace + Wy_furnace_distil = 1
_C158: CD_reactor_distil = 0
_C159: - L_reactor_distil + R_reactor_distil + X_distil - X_reactor = 0
_C16: Y_comp - 0.5 d_comp >= 0
_C160: A_reactor_distil - B_reactor_distil + Y_distil - Y_reactor = 0
_C161: R_reactor_distil - 1000 Wx_reactor_distil <= 0
_C162: L_reactor_distil + 1000 Wx_reactor_distil <= 1000
_C163: A_reactor_distil - 1000 Wy_reactor_distil <= 0
_C164: B_reactor_distil + 1000 Wy_reactor_distil <= 1000
_C165: - A_reactor_distil - B_reactor_distil + D_reactor_distil
 - L_reactor_distil - R_reactor_distil = 0
_C166: - D_distil_reactor + D_reactor_distil = 0
_C167: 1000 E1_reactor_distil + 1000 E2_reactor_distil - X_distil + X_reactor
 - 0.5 l_distil - 0.5 l_reactor >= 0
_C168: - 1000 E1_reactor_distil + 1000 E2_reactor_distil + X_distil
 - X_reactor - 0.5 l_distil - 0.5 l_reactor >= -1000
_C169: 1000 E1_reactor_distil - 1000 E2_reactor_distil - Y_distil + Y_reactor
 - 0.5 d_distil - 0.5 d_reactor >= -1000
_C17: 0 O_distil + l_distil = 5.5
_C170: - 1000 E1_reactor_distil - 1000 E2_reactor_distil + Y_distil
 - Y_reactor - 0.5 d_distil - 0.5 d_reactor >= -2000
_C171: L_reactor_distil - R_distil_reactor = 0
_C172: - L_distil_reactor + R_reactor_distil = 0
_C173: A_reactor_distil - B_distil_reactor = 0
_C174: - A_distil_reactor + B_reactor_distil = 0
_C175: Wx_distil_reactor + Wx_reactor_distil = 1
_C176: Wy_distil_reactor + Wy_reactor_distil = 1
_C177: CD_flash_distil - 39.2294860204 D_flash_distil = 0
_C178: - L_flash_distil + R_flash_distil + X_distil - X_flash = 0
_C179: A_flash_distil - B_flash_distil + Y_distil - Y_flash = 0
_C18: d_distil + l_distil = 11
_C180: R_flash_distil - 1000 Wx_flash_distil <= 0
_C181: L_flash_distil + 1000 Wx_flash_distil <= 1000
_C182: A_flash_distil - 1000 Wy_flash_distil <= 0
_C183: B_flash_distil + 1000 Wy_flash_distil <= 1000
_C184: - A_flash_distil - B_flash_distil + D_flash_distil - L_flash_distil
 - R_flash_distil = 0
_C185: - D_distil_flash + D_flash_distil = 0
_C186: 1000 E1_flash_distil + 1000 E2_flash_distil - X_distil + X_flash
 - 0.5 l_distil - 0.5 l_flash >= 0
_C187: - 1000 E1_flash_distil + 1000 E2_flash_distil + X_distil - X_flash
 - 0.5 l_distil - 0.5 l_flash >= -1000
_C188: 1000 E1_flash_distil - 1000 E2_flash_distil - Y_distil + Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -1000
_C189: - 1000 E1_flash_distil - 1000 E2_flash_distil + Y_distil - Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -2000
_C19: X_distil - 0.5 l_distil >= 0
_C190: L_flash_distil - R_distil_flash = 0
_C191: - L_distil_flash + R_flash_distil = 0
_C192: A_flash_distil - B_distil_flash = 0
_C193: - A_distil_flash + B_flash_distil = 0
_C194: Wx_distil_flash + Wx_flash_distil = 1
_C195: Wy_distil_flash + Wy_flash_distil = 1
_C196: CD_comp_distil = 0
_C197: - L_comp_distil + R_comp_distil - X_comp + X_distil = 0
_C198: A_comp_distil - B_comp_distil - Y_comp + Y_distil = 0
_C199: R_comp_distil - 1000 Wx_comp_distil <= 0
_C2: d_furnace + l_furnace = 2
_C20: Y_distil - 0.5 d_distil >= 0
_C200: L_comp_distil + 1000 Wx_comp_distil <= 1000
_C201: A_comp_distil - 1000 Wy_comp_distil <= 0
_C202: B_comp_distil + 1000 Wy_comp_distil <= 1000
_C203: - A_comp_distil - B_comp_distil + D_comp_distil - L_comp_distil
 - R_comp_distil = 0
_C204: D_comp_distil - D_distil_comp = 0
_C205: 1000 E1_comp_distil + 1000 E2_comp_distil + X_comp - X_distil
 - 0.5 l_comp - 0.5 l_distil >= 0
_C206: - 1000 E1_comp_distil + 1000 E2_comp_distil - X_comp + X_distil
 - 0.5 l_comp - 0.5 l_distil >= -1000
_C207: 1000 E1_comp_distil - 1000 E2_comp_distil + Y_comp - Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -1000
_C208: - 1000 E1_comp_distil - 1000 E2_comp_distil - Y_comp + Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -2000
_C209: L_comp_distil - R_distil_comp = 0
_C21: 0 O_store + l_store = 22.8
_C210: - L_distil_comp + R_comp_distil = 0
_C211: A_comp_distil - B_distil_comp = 0
_C212: - A_distil_comp + B_comp_distil = 0
_C213: Wx_comp_distil + Wx_distil_comp = 1
_C214: Wy_comp_distil + Wy_distil_comp = 1
_C215: CD_furnace_store - 37.3877487445 D_furnace_store = 0
_C216: - L_furnace_store + R_furnace_store - X_furnace + X_store = 0
_C217: A_furnace_store - B_furnace_store - Y_furnace + Y_store = 0
_C218: R_furnace_store - 1000 Wx_furnace_store <= 0
_C219: L_furnace_store + 1000 Wx_furnace_store <= 1000
_C22: d_store + l_store = 45.6
_C220: A_furnace_store - 1000 Wy_furnace_store <= 0
_C221: B_furnace_store + 1000 Wy_furnace_store <= 1000
_C222: - A_furnace_store - B_furnace_store + D_furnace_store - L_furnace_store
 - R_furnace_store = 0
_C223: D_furnace_store - D_store_furnace = 0
_C224: 1000 E1_furnace_store + 1000 E2_furnace_store + X_furnace - X_store
 - 0.5 l_furnace - 0.5 l_store >= 0
_C225: - 1000 E1_furnace_store + 1000 E2_furnace_store - X_furnace + X_store
 - 0.5 l_furnace - 0.5 l_store >= -1000
_C226: 1000 E1_furnace_store - 1000 E2_furnace_store + Y_furnace - Y_store
 - 0.5 d_furnace - 0.5 d_store >= -1000
_C227: - 1000 E1_furnace_store - 1000 E2_furnace_store - Y_furnace + Y_store
 - 0.5 d_furnace - 0.5 d_store >= -2000
_C228: L_furnace_store - R_store_furnace = 0
_C229: - L_store_furnace + R_furnace_store = 0
_C23: X_store - 0.5 l_store >= 0
_C230: A_furnace_store - B_store_furnace = 0
_C231: - A_store_furnace + B_furnace_store = 0
_C232: Wx_furnace_store + Wx_store_furnace = 1
_C233: Wy_furnace_store + Wy_store_furnace = 1
_C234: CD_reactor_store = 0
_C235: - L_reactor_store + R_reactor_store - X_reactor + X_store = 0
_C236: A_reactor_store - B_reactor_store - Y_reactor + Y_store = 0
_C237: R_reactor_store - 1000 Wx_reactor_store <= 0
_C238: L_reactor_store + 1000 Wx_reactor_store <= 1000
_C239: A_reactor_store - 1000 Wy_reactor_store <= 0
_C24: Y_store - 0.5 d_store >= 0
_C240: B_reactor_store + 1000 Wy_reactor_store <= 1000
_C241: - A_reactor_store - B_reactor_store + D_reactor_store - L_reactor_store
 - R_reactor_store = 0
_C242: D_reactor_store - D_store_reactor = 0
_C243: 1000 E1_reactor_store + 1000 E2_reactor_store + X_reactor - X_store
 - 0.5 l_reactor - 0.5 l_store >= 0
_C244: - 1000 E1_reactor_store + 1000 E2_reactor_store - X_reactor + X_store
 - 0.5 l_reactor - 0.5 l_store >= -1000
_C245: 1000 E1_reactor_store - 1000 E2_reactor_store + Y_reactor - Y_store
 - 0.5 d_reactor - 0.5 d_store >= -1000
_C246: - 1000 E1_reactor_store - 1000 E2_reactor_store - Y_reactor + Y_store
 - 0.5 d_reactor - 0.5 d_store >= -2000
_C247: L_reactor_store - R_store_reactor = 0
_C248: - L_store_reactor + R_reactor_store = 0
_C249: A_reactor_store - B_store_reactor = 0
_C25: CD_furnace_reactor - 390.264455871 D_furnace_reactor = 0
_C250: - A_store_reactor + B_reactor_store = 0
_C251: Wx_reactor_store + Wx_store_reactor = 1
_C252: Wy_reactor_store + Wy_store_reactor = 1
_C253: CD_flash_store = 0
_C254: - L_flash_store + R_flash_store - X_flash + X_store = 0
_C255: A_flash_store - B_flash_store - Y_flash + Y_store = 0
_C256: R_flash_store - 1000 Wx_flash_store <= 0
_C257: L_flash_store + 1000 Wx_flash_store <= 1000
_C258: A_flash_store - 1000 Wy_flash_store <= 0
_C259: B_flash_store + 1000 Wy_flash_store <= 1000
_C26: - L_furnace_reactor + R_furnace_reactor - X_furnace + X_reactor = 0
_C260: - A_flash_store - B_flash_store + D_flash_store - L_flash_store
 - R_flash_store = 0
_C261: D_flash_store - D_store_flash = 0
_C262: 1000 E1_flash_store + 1000 E2_flash_store + X_flash - X_store
 - 0.5 l_flash - 0.5 l_store >= 0
_C263: - 1000 E1_flash_store + 1000 E2_flash_store - X_flash + X_store
 - 0.5 l_flash - 0.5 l_store >= -1000
_C264: 1000 E1_flash_store - 1000 E2_flash_store + Y_flash - Y_store
 - 0.5 d_flash - 0.5 d_store >= -1000
_C265: - 1000 E1_flash_store - 1000 E2_flash_store - Y_flash + Y_store
 - 0.5 d_flash - 0.5 d_store >= -2000
_C266: L_flash_store - R_store_flash = 0
_C267: - L_store_flash + R_flash_store = 0
_C268: A_flash_store - B_store_flash = 0
_C269: - A_store_flash + B_flash_store = 0
_C27: A_furnace_reactor - B_furnace_reactor - Y_furnace + Y_reactor = 0
_C270: Wx_flash_store + Wx_store_flash = 1
_C271: Wy_flash_store + Wy_store_flash = 1
_C272: CD_comp_store = 0
_C273: - L_comp_store + R_comp_store - X_comp + X_store = 0
_C274: A_comp_store - B_comp_store - Y_comp + Y_store = 0
_C275: R_comp_store - 1000 Wx_comp_store <= 0
_C276: L_comp_store + 1000 Wx_comp_store <= 1000
_C277: A_comp_store - 1000 Wy_comp_store <= 0
_C278: B_comp_store + 1000 Wy_comp_store <= 1000
_C279: - A_comp_store - B_comp_store + D_comp_store - L_comp_store
 - R_comp_store = 0
_C28: R_furnace_reactor - 1000 Wx_furnace_reactor <= 0
_C280: D_comp_store - D_store_comp = 0
_C281: 1000 E1_comp_store + 1000 E2_comp_store + X_comp - X_store - 0.5 l_comp
 - 0.5 l_store >= 0
_C282: - 1000 E1_comp_store + 1000 E2_comp_store - X_comp + X_store
 - 0.5 l_comp - 0.5 l_store >= -1000
_C283: 1000 E1_comp_store - 1000 E2_comp_store + Y_comp - Y_store - 0.5 d_comp
 - 0.5 d_store >= -1000
_C284: - 1000 E1_comp_store - 1000 E2_comp_store - Y_comp + Y_store
 - 0.5 d_comp - 0.5 d_store >= -2000
_C285: L_comp_store - R_store_comp = 0
_C286: - L_store_comp + R_comp_store = 0
_C287: A_comp_store - B_store_comp = 0
_C288: - A_store_comp + B_comp_store = 0
_C289: Wx_comp_store + Wx_store_comp = 1
_C29: L_furnace_reactor + 1000 Wx_furnace_reactor <= 1000
_C290: Wy_comp_store + Wy_store_comp = 1
_C291: CD_distil_store - 198.966774782 D_distil_store = 0
_C292: - L_distil_store + R_distil_store - X_distil + X_store = 0
_C293: A_distil_store - B_distil_store - Y_distil + Y_store = 0
_C294: R_distil_store - 1000 Wx_distil_store <= 0
_C295: L_distil_store + 1000 Wx_distil_store <= 1000
_C296: A_distil_store - 1000 Wy_distil_store <= 0
_C297: B_distil_store + 1000 Wy_distil_store <= 1000
_C298: - A_distil_store - B_distil_store + D_distil_store - L_distil_store
 - R_distil_store = 0
_C299: D_distil_store - D_store_distil = 0
_C3: X_furnace - 0.5 l_furnace >= 0
_C30: A_furnace_reactor - 1000 Wy_furnace_reactor <= 0
_C300: 1000 E1_distil_store + 1000 E2_distil_store + X_distil - X_store
 - 0.5 l_distil - 0.5 l_store >= 0
_C301: - 1000 E1_distil_store + 1000 E2_distil_store - X_distil + X_store
 - 0.5 l_distil - 0.5 l_store >= -1000
_C302: 1000 E1_distil_store - 1000 E2_distil_store + Y_distil - Y_store
 - 0.5 d_distil - 0.5 d_store >= -1000
_C303: - 1000 E1_distil_store - 1000 E2_distil_store - Y_distil + Y_store
 - 0.5 d_distil - 0.5 d_store >= -2000
_C304: L_distil_store - R_store_distil = 0
_C305: - L_store_distil + R_distil_store = 0
_C306: A_distil_store - B_store_distil = 0
_C307: - A_store_distil + B_distil_store = 0
_C308: Wx_distil_store + Wx_store_distil = 1
_C309: Wy_distil_store + Wy_store_distil = 1
_C31: B_furnace_reactor + 1000 Wy_furnace_reactor <= 1000
_C310: - CD_comp_comp - CD_comp_distil - CD_comp_flash - CD_comp_furnace
 - CD_comp_reactor - CD_comp_store - CD_distil_comp - CD_distil_distil
 - CD_distil_flash - CD_distil_furnace - CD_distil_reactor - CD_distil_store
 - CD_flash_comp - CD_flash_distil - CD_flash_flash - CD_flash_furnace
 - CD_flash_reactor - CD_flash_store - CD_furnace_comp - CD_furnace_distil
 - CD_furnace_flash - CD_furnace_furnace - CD_furnace_reactor
 - CD_furnace_store - CD_reactor_comp - CD_reactor_distil - CD_reactor_flash
 - CD_reactor_furnace - CD_reactor_reactor - CD_reactor_store - CD_store_comp
 - CD_store_distil - CD_store_flash - CD_store_furnace - CD_store_reactor
 - CD_store_store + SumCD = 0
_C311: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + X_furnace + 0.5 l_furnace <= 0
_C312: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + Y_furnace + 0.5 d_furnace <= 0
_C313: Gn_1 + Gn_10 + Gn_11 + Gn_12 + Gn_13 + Gn_14 + Gn_15 + Gn_16 + Gn_17
 + Gn_18 + Gn_19 + Gn_2 + Gn_20 + Gn_21 + Gn_22 + Gn_23 + Gn_24 + Gn_25
 + Gn_26 + Gn_27 + Gn_28 + Gn_29 + Gn_3 + Gn_30 + Gn_31 + Gn_32 + Gn_33
 + Gn_34 + Gn_35 + Gn_36 + Gn_37 + Gn_38 + Gn_39 + Gn_4 + Gn_40 + Gn_41
 + Gn_42 + Gn_43 + Gn_44 + Gn_45 + Gn_46 + Gn_47 + Gn_48 + Gn_49 + Gn_5
 + Gn_50 + Gn_51 + Gn_52 + Gn_53 + Gn_54 + Gn_55 + Gn_56 + Gn_57 + Gn_58
 + Gn_59 + Gn_6 + Gn_7 + Gn_8 + Gn_9 = 1
_C314: - 3000 Gn_1 - 300000 Gn_10 - 363000 Gn_11 - 432000 Gn_12 - 507000 Gn_13
 - 588000 Gn_14 - 675000 Gn_15 - 768000 Gn_16 - 867000 Gn_17 - 972000 Gn_18
 - 1083000 Gn_19 - 12000 Gn_2 - 1200000 Gn_20 - 1323000 Gn_21 - 1452000 Gn_22
 - 1587000 Gn_23 - 1728000 Gn_24 - 1875000 Gn_25 - 2028000 Gn_26
 - 2187000 Gn_27 - 2352000 Gn_28 - 2523000 Gn_29 - 27000 Gn_3 - 2700000 Gn_30
 - 2883000 Gn_31 - 3072000 Gn_32 - 3267000 Gn_33 - 3468000 Gn_34
 - 3675000 Gn_35 - 3888000 Gn_36 - 4107000 Gn_37 - 4332000 Gn_38
 - 4563000 Gn_39 - 48000 Gn_4 - 4800000 Gn_40 - 5043000 Gn_41 - 5292000 Gn_42
 - 5547000 Gn_43 - 5808000 Gn_44 - 6075000 Gn_45 - 6348000 Gn_46
 - 6627000 Gn_47 - 6912000 Gn_48 - 7203000 Gn_49 - 75000 Gn_5 - 7500000 Gn_50
 - 7803000 Gn_51 - 8112000 Gn_52 - 8427000 Gn_53 - 8748000 Gn_54
 - 9075000 Gn_55 - 9408000 Gn_56 - 9747000 Gn_57 - 10092000 Gn_58
 - 10443000 Gn_59 - 108000 Gn_6 - 147000 Gn_7 - 192000 Gn_8 - 243000 Gn_9
 + TLC = 0
_C315: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + X_reactor + 0.5 l_reactor <= 0
_C316: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + Y_reactor + 0.5 d_reactor <= 0
_C317: Gn_1 + Gn_10 + Gn_11 + Gn_12 + Gn_13 + Gn_14 + Gn_15 + Gn_16 + Gn_17
 + Gn_18 + Gn_19 + Gn_2 + Gn_20 + Gn_21 + Gn_22 + Gn_23 + Gn_24 + Gn_25
 + Gn_26 + Gn_27 + Gn_28 + Gn_29 + Gn_3 + Gn_30 + Gn_31 + Gn_32 + Gn_33
 + Gn_34 + Gn_35 + Gn_36 + Gn_37 + Gn_38 + Gn_39 + Gn_4 + Gn_40 + Gn_41
 + Gn_42 + Gn_43 + Gn_44 + Gn_45 + Gn_46 + Gn_47 + Gn_48 + Gn_49 + Gn_5
 + Gn_50 + Gn_51 + Gn_52 + Gn_53 + Gn_54 + Gn_55 + Gn_56 + Gn_57 + Gn_58
 + Gn_59 + Gn_6 + Gn_7 + Gn_8 + Gn_9 = 1
_C318: - 3000 Gn_1 - 300000 Gn_10 - 363000 Gn_11 - 432000 Gn_12 - 507000 Gn_13
 - 588000 Gn_14 - 675000 Gn_15 - 768000 Gn_16 - 867000 Gn_17 - 972000 Gn_18
 - 1083000 Gn_19 - 12000 Gn_2 - 1200000 Gn_20 - 1323000 Gn_21 - 1452000 Gn_22
 - 1587000 Gn_23 - 1728000 Gn_24 - 1875000 Gn_25 - 2028000 Gn_26
 - 2187000 Gn_27 - 2352000 Gn_28 - 2523000 Gn_29 - 27000 Gn_3 - 2700000 Gn_30
 - 2883000 Gn_31 - 3072000 Gn_32 - 3267000 Gn_33 - 3468000 Gn_34
 - 3675000 Gn_35 - 3888000 Gn_36 - 4107000 Gn_37 - 4332000 Gn_38
 - 4563000 Gn_39 - 48000 Gn_4 - 4800000 Gn_40 - 5043000 Gn_41 - 5292000 Gn_42
 - 5547000 Gn_43 - 5808000 Gn_44 - 6075000 Gn_45 - 6348000 Gn_46
 - 6627000 Gn_47 - 6912000 Gn_48 - 7203000 Gn_49 - 75000 Gn_5 - 7500000 Gn_50
 - 7803000 Gn_51 - 8112000 Gn_52 - 8427000 Gn_53 - 8748000 Gn_54
 - 9075000 Gn_55 - 9408000 Gn_56 - 9747000 Gn_57 - 10092000 Gn_58
 - 10443000 Gn_59 - 108000 Gn_6 - 147000 Gn_7 - 192000 Gn_8 - 243000 Gn_9
 + TLC = 0
_C319: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + X_flash + 0.5 l_flash <= 0
_C32: - A_furnace_reactor - B_furnace_reactor + D_furnace_reactor
 - L_furnace_reactor - R_furnace_reactor = 0
_C320: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + Y_flash + 0.5 d_flash <= 0
_C321: Gn_1 + Gn_10 + Gn_11 + Gn_12 + Gn_13 + Gn_14 + Gn_15 + Gn_16 + Gn_17
 + Gn_18 + Gn_19 + Gn_2 + Gn_20 + Gn_21 + Gn_22 + Gn_23 + Gn_24 + Gn_25
 + Gn_26 + Gn_27 + Gn_28 + Gn_29 + Gn_3 + Gn_30 + Gn_31 + Gn_32 + Gn_33
 + Gn_34 + Gn_35 + Gn_36 + Gn_37 + Gn_38 + Gn_39 + Gn_4 + Gn_40 + Gn_41
 + Gn_42 + Gn_43 + Gn_44 + Gn_45 + Gn_46 + Gn_47 + Gn_48 + Gn_49 + Gn_5
 + Gn_50 + Gn_51 + Gn_52 + Gn_53 + Gn_54 + Gn_55 + Gn_56 + Gn_57 + Gn_58
 + Gn_59 + Gn_6 + Gn_7 + Gn_8 + Gn_9 = 1
_C322: - 3000 Gn_1 - 300000 Gn_10 - 363000 Gn_11 - 432000 Gn_12 - 507000 Gn_13
 - 588000 Gn_14 - 675000 Gn_15 - 768000 Gn_16 - 867000 Gn_17 - 972000 Gn_18
 - 1083000 Gn_19 - 12000 Gn_2 - 1200000 Gn_20 - 1323000 Gn_21 - 1452000 Gn_22
 - 1587000 Gn_23 - 1728000 Gn_24 - 1875000 Gn_25 - 2028000 Gn_26
 - 2187000 Gn_27 - 2352000 Gn_28 - 2523000 Gn_29 - 27000 Gn_3 - 2700000 Gn_30
 - 2883000 Gn_31 - 3072000 Gn_32 - 3267000 Gn_33 - 3468000 Gn_34
 - 3675000 Gn_35 - 3888000 Gn_36 - 4107000 Gn_37 - 4332000 Gn_38
 - 4563000 Gn_39 - 48000 Gn_4 - 4800000 Gn_40 - 5043000 Gn_41 - 5292000 Gn_42
 - 5547000 Gn_43 - 5808000 Gn_44 - 6075000 Gn_45 - 6348000 Gn_46
 - 6627000 Gn_47 - 6912000 Gn_48 - 7203000 Gn_49 - 75000 Gn_5 - 7500000 Gn_50
 - 7803000 Gn_51 - 8112000 Gn_52 - 8427000 Gn_53 - 8748000 Gn_54
 - 9075000 Gn_55 - 9408000 Gn_56 - 9747000 Gn_57 - 10092000 Gn_58
 - 10443000 Gn_59 - 108000 Gn_6 - 147000 Gn_7 - 192000 Gn_8 - 243000 Gn_9
 + TLC = 0
_C323: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + X_comp + 0.5 l_comp <= 0
_C324: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + Y_comp + 0.5 d_comp <= 0
_C325: Gn_1 + Gn_10 + Gn_11 + Gn_12 + Gn_13 + Gn_14 + Gn_15 + Gn_16 + Gn_17
 + Gn_18 + Gn_19 + Gn_2 + Gn_20 + Gn_21 + Gn_22 + Gn_23 + Gn_24 + Gn_25
 + Gn_26 + Gn_27 + Gn_28 + Gn_29 + Gn_3 + Gn_30 + Gn_31 + Gn_32 + Gn_33
 + Gn_34 + Gn_35 + Gn_36 + Gn_37 + Gn_38 + Gn_39 + Gn_4 + Gn_40 + Gn_41
 + Gn_42 + Gn_43 + Gn_44 + Gn_45 + Gn_46 + Gn_47 + Gn_48 + Gn_49 + Gn_5
 + Gn_50 + Gn_51 + Gn_52 + Gn_53 + Gn_54 + Gn_55 + Gn_56 + Gn_57 + Gn_58
 + Gn_59 + Gn_6 + Gn_7 + Gn_8 + Gn_9 = 1
_C326: - 3000 Gn_1 - 300000 Gn_10 - 363000 Gn_11 - 432000 Gn_12 - 507000 Gn_13
 - 588000 Gn_14 - 675000 Gn_15 - 768000 Gn_16 - 867000 Gn_17 - 972000 Gn_18
 - 1083000 Gn_19 - 12000 Gn_2 - 1200000 Gn_20 - 1323000 Gn_21 - 1452000 Gn_22
 - 1587000 Gn_23 - 1728000 Gn_24 - 1875000 Gn_25 - 2028000 Gn_26
 - 2187000 Gn_27 - 2352000 Gn_28 - 2523000 Gn_29 - 27000 Gn_3 - 2700000 Gn_30
 - 2883000 Gn_31 - 3072000 Gn_32 - 3267000 Gn_33 - 3468000 Gn_34
 - 3675000 Gn_35 - 3888000 Gn_36 - 4107000 Gn_37 - 4332000 Gn_38
 - 4563000 Gn_39 - 48000 Gn_4 - 4800000 Gn_40 - 5043000 Gn_41 - 5292000 Gn_42
 - 5547000 Gn_43 - 5808000 Gn_44 - 6075000 Gn_45 - 6348000 Gn_46
 - 6627000 Gn_47 - 6912000 Gn_48 - 7203000 Gn_49 - 75000 Gn_5 - 7500000 Gn_50
 - 7803000 Gn_51 - 8112000 Gn_52 - 8427000 Gn_53 - 8748000 Gn_54
 - 9075000 Gn_55 - 9408000 Gn_56 - 9747000 Gn_57 - 10092000 Gn_58
 - 10443000 Gn_59 - 108000 Gn_6 - 147000 Gn_7 - 192000 Gn_8 - 243000 Gn_9
 + TLC = 0
_C327: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + X_distil + 0.5 l_distil <= 0
_C328: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + Y_distil + 0.5 d_distil <= 0
_C329: Gn_1 + Gn_10 + Gn_11 + Gn_12 + Gn_13 + Gn_14 + Gn_15 + Gn_16 + Gn_17
 + Gn_18 + Gn_19 + Gn_2 + Gn_20 + Gn_21 + Gn_22 + Gn_23 + Gn_24 + Gn_25
 + Gn_26 + Gn_27 + Gn_28 + Gn_29 + Gn_3 + Gn_30 + Gn_31 + Gn_32 + Gn_33
 + Gn_34 + Gn_35 + Gn_36 + Gn_37 + Gn_38 + Gn_39 + Gn_4 + Gn_40 + Gn_41
 + Gn_42 + Gn_43 + Gn_44 + Gn_45 + Gn_46 + Gn_47 + Gn_48 + Gn_49 + Gn_5
 + Gn_50 + Gn_51 + Gn_52 + Gn_53 + Gn_54 + Gn_55 + Gn_56 + Gn_57 + Gn_58
 + Gn_59 + Gn_6 + Gn_7 + Gn_8 + Gn_9 = 1
_C33: D_furnace_reactor - D_reactor_furnace = 0
_C330: - 3000 Gn_1 - 300000 Gn_10 - 363000 Gn_11 - 432000 Gn_12 - 507000 Gn_13
 - 588000 Gn_14 - 675000 Gn_15 - 768000 Gn_16 - 867000 Gn_17 - 972000 Gn_18
 - 1083000 Gn_19 - 12000 Gn_2 - 1200000 Gn_20 - 1323000 Gn_21 - 1452000 Gn_22
 - 1587000 Gn_23 - 1728000 Gn_24 - 1875000 Gn_25 - 2028000 Gn_26
 - 2187000 Gn_27 - 2352000 Gn_28 - 2523000 Gn_29 - 27000 Gn_3 - 2700000 Gn_30
 - 2883000 Gn_31 - 3072000 Gn_32 - 3267000 Gn_33 - 3468000 Gn_34
 - 3675000 Gn_35 - 3888000 Gn_36 - 4107000 Gn_37 - 4332000 Gn_38
 - 4563000 Gn_39 - 48000 Gn_4 - 4800000 Gn_40 - 5043000 Gn_41 - 5292000 Gn_42
 - 5547000 Gn_43 - 5808000 Gn_44 - 6075000 Gn_45 - 6348000 Gn_46
 - 6627000 Gn_47 - 6912000 Gn_48 - 7203000 Gn_49 - 75000 Gn_5 - 7500000 Gn_50
 - 7803000 Gn_51 - 8112000 Gn_52 - 8427000 Gn_53 - 8748000 Gn_54
 - 9075000 Gn_55 - 9408000 Gn_56 - 9747000 Gn_57 - 10092000 Gn_58
 - 10443000 Gn_59 - 108000 Gn_6 - 147000 Gn_7 - 192000 Gn_8 - 243000 Gn_9
 + TLC = 0
_C331: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + X_store + 0.5 l_store <= 0
_C332: - Gn_1 - 10 Gn_10 - 11 Gn_11 - 12 Gn_12 - 13 Gn_13 - 14 Gn_14
 - 15 Gn_15 - 16 Gn_16 - 17 Gn_17 - 18 Gn_18 - 19 Gn_19 - 2 Gn_2 - 20 Gn_20
 - 21 Gn_21 - 22 Gn_22 - 23 Gn_23 - 24 Gn_24 - 25 Gn_25 - 26 Gn_26 - 27 Gn_27
 - 28 Gn_28 - 29 Gn_29 - 3 Gn_3 - 30 Gn_30 - 31 Gn_31 - 32 Gn_32 - 33 Gn_33
 - 34 Gn_34 - 35 Gn_35 - 36 Gn_36 - 37 Gn_37 - 38 Gn_38 - 39 Gn_39 - 4 Gn_4
 - 40 Gn_40 - 41 Gn_41 - 42 Gn_42 - 43 Gn_43 - 44 Gn_44 - 45 Gn_45 - 46 Gn_46
 - 47 Gn_47 - 48 Gn_48 - 49 Gn_49 - 5 Gn_5 - 50 Gn_50 - 51 Gn_51 - 52 Gn_52
 - 53 Gn_53 - 54 Gn_54 - 55 Gn_55 - 56 Gn_56 - 57 Gn_57 - 58 Gn_58 - 59 Gn_59
 - 6 Gn_6 - 7 Gn_7 - 8 Gn_8 - 9 Gn_9 + Y_store + 0.5 d_store <= 0
_C333: Gn_1 + Gn_10 + Gn_11 + Gn_12 + Gn_13 + Gn_14 + Gn_15 + Gn_16 + Gn_17
 + Gn_18 + Gn_19 + Gn_2 + Gn_20 + Gn_21 + Gn_22 + Gn_23 + Gn_24 + Gn_25
 + Gn_26 + Gn_27 + Gn_28 + Gn_29 + Gn_3 + Gn_30 + Gn_31 + Gn_32 + Gn_33
 + Gn_34 + Gn_35 + Gn_36 + Gn_37 + Gn_38 + Gn_39 + Gn_4 + Gn_40 + Gn_41
 + Gn_42 + Gn_43 + Gn_44 + Gn_45 + Gn_46 + Gn_47 + Gn_48 + Gn_49 + Gn_5
 + Gn_50 + Gn_51 + Gn_52 + Gn_53 + Gn_54 + Gn_55 + Gn_56 + Gn_57 + Gn_58
 + Gn_59 + Gn_6 + Gn_7 + Gn_8 + Gn_9 = 1
_C334: - 3000 Gn_1 - 300000 Gn_10 - 363000 Gn_11 - 432000 Gn_12 - 507000 Gn_13
 - 588000 Gn_14 - 675000 Gn_15 - 768000 Gn_16 - 867000 Gn_17 - 972000 Gn_18
 - 1083000 Gn_19 - 12000 Gn_2 - 1200000 Gn_20 - 1323000 Gn_21 - 1452000 Gn_22
 - 1587000 Gn_23 - 1728000 Gn_24 - 1875000 Gn_25 - 2028000 Gn_26
 - 2187000 Gn_27 - 2352000 Gn_28 - 2523000 Gn_29 - 27000 Gn_3 - 2700000 Gn_30
 - 2883000 Gn_31 - 3072000 Gn_32 - 3267000 Gn_33 - 3468000 Gn_34
 - 3675000 Gn_35 - 3888000 Gn_36 - 4107000 Gn_37 - 4332000 Gn_38
 - 4563000 Gn_39 - 48000 Gn_4 - 4800000 Gn_40 - 5043000 Gn_41 - 5292000 Gn_42
 - 5547000 Gn_43 - 5808000 Gn_44 - 6075000 Gn_45 - 6348000 Gn_46
 - 6627000 Gn_47 - 6912000 Gn_48 - 7203000 Gn_49 - 75000 Gn_5 - 7500000 Gn_50
 - 7803000 Gn_51 - 8112000 Gn_52 - 8427000 Gn_53 - 8748000 Gn_54
 - 9075000 Gn_55 - 9408000 Gn_56 - 9747000 Gn_57 - 10092000 Gn_58
 - 10443000 Gn_59 - 108000 Gn_6 - 147000 Gn_7 - 192000 Gn_8 - 243000 Gn_9
 + TLC = 0
_C335: Ve2_furnace_furnace = 0
_C336: - D_furnace_reactor + Din_furnace_reactor + Dout_furnace_reactor = 0
_C337: Din_furnace_reactor - 28.8 Psie_furnace_reactor <= 0
_C338: Dout_furnace_reactor + 28.8 Psie_furnace_reactor >= 28.8
_C339: Dout_furnace_reactor + 1000 Psie_furnace_reactor <= 1000
_C34: 1000 E1_furnace_reactor + 1000 E2_furnace_reactor + X_furnace
 - X_reactor - 0.5 l_furnace - 0.5 l_reactor >= 0
_C340: 41260.4166667 Din_furnace_reactor - 1188300 Psie_furnace_reactor
 + Ve2_furnace_reactor = 0
_C341: - D_furnace_flash + Din_furnace_flash + Dout_furnace_flash = 0
_C342: Din_furnace_flash - 28.8 Psie_furnace_flash <= 0
_C343: Dout_furnace_flash + 28.8 Psie_furnace_flash >= 28.8
_C344: Dout_furnace_flash + 1000 Psie_furnace_flash <= 1000
_C345: 25888.8888889 Din_furnace_flash - 745600 Psie_furnace_flash
 + Ve2_furnace_flash = 0
_C346: - D_furnace_comp + Din_furnace_comp + Dout_furnace_comp = 0
_C347: Din_furnace_comp - 28.8 Psie_furnace_comp <= 0
_C348: Dout_furnace_comp + 28.8 Psie_furnace_comp >= 28.8
_C349: Dout_furnace_comp + 1000 Psie_furnace_comp <= 1000
_C35: - 1000 E1_furnace_reactor + 1000 E2_furnace_reactor - X_furnace
 + X_reactor - 0.5 l_furnace - 0.5 l_reactor >= -1000
_C350: 95402.7777778 Din_furnace_comp - 2747600 Psie_furnace_comp
 + Ve2_furnace_comp = 0
_C351: - D_furnace_distil + Din_furnace_distil + Dout_furnace_distil = 0
_C352: Din_furnace_distil - 28.8 Psie_furnace_distil <= 0
_C353: Dout_furnace_distil + 28.8 Psie_furnace_distil >= 28.8
_C354: Dout_furnace_distil + 1000 Psie_furnace_distil <= 1000
_C355: 120750 Din_furnace_distil - 3477600 Psie_furnace_distil
 + Ve2_furnace_distil = 0
_C356: - D_furnace_store + Din_furnace_store + Dout_furnace_store = 0
_C357: Din_furnace_store - 28.8 Psie_furnace_store <= 0
_C358: Dout_furnace_store + 28.8 Psie_furnace_store >= 28.8
_C359: Dout_furnace_store + 1000 Psie_furnace_store <= 1000
_C36: 1000 E1_furnace_reactor - 1000 E2_furnace_reactor + Y_furnace
 - Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -1000
_C360: 162529.930556 Din_furnace_store - 4680862 Psie_furnace_store
 + Ve2_furnace_store = 0
_C361: - D_reactor_furnace + Din_reactor_furnace + Dout_reactor_furnace = 0
_C362: Din_reactor_furnace - 23.9 Psie_reactor_furnace <= 0
_C363: Dout_reactor_furnace + 23.9 Psie_reactor_furnace >= 23.9
_C364: Dout_reactor_furnace + 1000 Psie_reactor_furnace <= 1000
_C365: 31644.3514644 Din_reactor_furnace - 756300 Psie_reactor_furnace
 + Ve2_reactor_furnace = 0
_C366: Ve2_reactor_reactor = 0
_C367: - D_reactor_flash + Din_reactor_flash + Dout_reactor_flash = 0
_C368: Din_reactor_flash - 23.9 Psie_reactor_flash <= 0
_C369: Dout_reactor_flash + 23.9 Psie_reactor_flash >= 23.9
_C37: - 1000 E1_furnace_reactor - 1000 E2_furnace_reactor - Y_furnace
 + Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -2000
_C370: Dout_reactor_flash + 1000 Psie_reactor_flash <= 1000
_C371: 31196.6527197 Din_reactor_flash - 745600 Psie_reactor_flash
 + Ve2_reactor_flash = 0
_C372: - D_reactor_comp + Din_reactor_comp + Dout_reactor_comp = 0
_C373: Din_reactor_comp - 23.9 Psie_reactor_comp <= 0
_C374: Dout_reactor_comp + 23.9 Psie_reactor_comp >= 23.9
_C375: Dout_reactor_comp + 1000 Psie_reactor_comp <= 1000
_C376: 114962.343096 Din_reactor_comp - 2747600 Psie_reactor_comp
 + Ve2_reactor_comp = 0
_C377: - D_reactor_distil + Din_reactor_distil + Dout_reactor_distil = 0
_C378: Din_reactor_distil - 23.9 Psie_reactor_distil <= 0
_C379: Dout_reactor_distil + 23.9 Psie_reactor_distil >= 23.9
_C38: L_furnace_reactor - R_reactor_furnace = 0
_C380: Dout_reactor_distil + 1000 Psie_reactor_distil <= 1000
_C381: 145506.276151 Din_reactor_distil - 3477600 Psie_reactor_distil
 + Ve2_reactor_distil = 0
_C382: - D_reactor_store + Din_reactor_store + Dout_reactor_store = 0
_C383: Din_reactor_store - 23.9 Psie_reactor_store <= 0
_C384: Dout_reactor_store + 23.9 Psie_reactor_store >= 23.9
_C385: Dout_reactor_store + 1000 Psie_reactor_store <= 1000
_C386: 195851.966527 Din_reactor_store - 4680862 Psie_reactor_store
 + Ve2_reactor_store = 0
_C387: - D_distil_furnace + Din_distil_furnace + Dout_distil_furnace = 0
_C388: Din_distil_furnace - 35.4 Psie_distil_furnace <= 0
_C389: Dout_distil_furnace + 35.4 Psie_distil_furnace >= 35.4
_C39: - L_reactor_furnace + R_furnace_reactor = 0
_C390: Dout_distil_furnace + 1000 Psie_distil_furnace <= 1000
_C391: 21364.4067797 Din_distil_furnace - 756300 Psie_distil_furnace
 + Ve2_distil_furnace = 0
_C392: - D_distil_reactor + Din_distil_reactor + Dout_distil_reactor = 0
_C393: Din_distil_reactor - 35.4 Psie_distil_reactor <= 0
_C394: Dout_distil_reactor + 35.4 Psie_distil_reactor >= 35.4
_C395: Dout_distil_reactor + 1000 Psie_distil_reactor <= 1000
_C396: 33567.7966102 Din_distil_reactor - 1188300 Psie_distil_reactor
 + Ve2_distil_reactor = 0
_C397: - D_distil_flash + Din_distil_flash + Dout_distil_flash = 0
_C398: Din_distil_flash - 35.4 Psie_distil_flash <= 0
_C399: Dout_distil_flash + 35.4 Psie_distil_flash >= 35.4
_C4: Y_furnace - 0.5 d_furnace >= 0
_C40: A_furnace_reactor - B_reactor_furnace = 0
_C400: Dout_distil_flash + 1000 Psie_distil_flash <= 1000
_C401: 21062.1468927 Din_distil_flash - 745600 Psie_distil_flash
 + Ve2_distil_flash = 0
_C402: - D_distil_comp + Din_distil_comp + Dout_distil_comp = 0
_C403: Din_distil_comp - 35.4 Psie_distil_comp <= 0
_C404: Dout_distil_comp + 35.4 Psie_distil_comp >= 35.4
_C405: Dout_distil_comp + 1000 Psie_distil_comp <= 1000
_C406: 77615.819209 Din_distil_comp - 2747600 Psie_distil_comp
 + Ve2_distil_comp = 0
_C407: Ve2_distil_distil = 0
_C408: - D_distil_store + Din_distil_store + Dout_distil_store = 0
_C409: Din_distil_store - 35.4 Psie_distil_store <= 0
_C41: - A_reactor_furnace + B_furnace_reactor = 0
_C410: Dout_distil_store + 35.4 Psie_distil_store >= 35.4
_C411: Dout_distil_store + 1000 Psie_distil_store <= 1000
_C412: 132227.740113 Din_distil_store - 4680862 Psie_distil_store
 + Ve2_distil_store = 0
_C413: - D_store_furnace + Din_store_furnace + Dout_store_furnace = 0
_C414: Din_store_furnace - 45 Psie_store_furnace <= 0
_C415: Dout_store_furnace + 45 Psie_store_furnace >= 45
_C416: Dout_store_furnace + 1000 Psie_store_furnace <= 1000
_C417: 16806.6666667 Din_store_furnace - 756300 Psie_store_furnace
 + Ve2_store_furnace = 0
_C418: - D_store_reactor + Din_store_reactor + Dout_store_reactor = 0
_C419: Din_store_reactor - 45 Psie_store_reactor <= 0
_C42: Wx_furnace_reactor + Wx_reactor_furnace = 1
_C420: Dout_store_reactor + 45 Psie_store_reactor >= 45
_C421: Dout_store_reactor + 1000 Psie_store_reactor <= 1000
_C422: 26406.6666667 Din_store_reactor - 1188300 Psie_store_reactor
 + Ve2_store_reactor = 0
_C423: - D_store_flash + Din_store_flash + Dout_store_flash = 0
_C424: Din_store_flash - 45 Psie_store_flash <= 0
_C425: Dout_store_flash + 45 Psie_store_flash >= 45
_C426: Dout_store_flash + 1000 Psie_store_flash <= 1000
_C427: 16568.8888889 Din_store_flash - 745600 Psie_store_flash
 + Ve2_store_flash = 0
_C428: - D_store_comp + Din_store_comp + Dout_store_comp = 0
_C429: Din_store_comp - 45 Psie_store_comp <= 0
_C43: Wy_furnace_reactor + Wy_reactor_furnace = 1
_C430: Dout_store_comp + 45 Psie_store_comp >= 45
_C431: Dout_store_comp + 1000 Psie_store_comp <= 1000
_C432: 61057.7777778 Din_store_comp - 2747600 Psie_store_comp + Ve2_store_comp
 = 0
_C433: - D_store_distil + Din_store_distil + Dout_store_distil = 0
_C434: Din_store_distil - 45 Psie_store_distil <= 0
_C435: Dout_store_distil + 45 Psie_store_distil >= 45
_C436: Dout_store_distil + 1000 Psie_store_distil <= 1000
_C437: 77280 Din_store_distil - 3477600 Psie_store_distil + Ve2_store_distil
 = 0
_C438: Ve2_store_store = 0
_C439: - Ve2_furnace_comp - Ve2_furnace_distil - Ve2_furnace_flash
 - Ve2_furnace_furnace - Ve2_furnace_reactor - Ve2_furnace_store + Ve_furnace
 = 756300
_C44: CD_furnace_flash = 0
_C440: Omega0_furnace - 0.75 Ve_furnace = 0
_C441: - Ve2_reactor_comp - Ve2_reactor_distil - Ve2_reactor_flash
 - Ve2_reactor_furnace - Ve2_reactor_reactor - Ve2_reactor_store + Ve_reactor
 = 1188300
_C442: Omega0_reactor - 0.66 Ve_reactor = 0
_C443: - Ve2_distil_comp - Ve2_distil_distil - Ve2_distil_flash
 - Ve2_distil_furnace - Ve2_distil_reactor - Ve2_distil_store + Ve_distil
 = 3477600
_C444: Omega0_distil - 0.77 Ve_distil = 0
_C445: - Ve2_store_comp - Ve2_store_distil - Ve2_store_flash
 - Ve2_store_furnace - Ve2_store_reactor - Ve2_store_store + Ve_store
 = 4680862
_C446: Omega0_store - 0.82 Ve_store = 0
_C447: - Omega0_distil - Omega0_furnace - Omega0_reactor - Omega0_store
 + SumOmega = 0
_C45: - L_furnace_flash + R_furnace_flash + X_flash - X_furnace = 0
_C46: A_furnace_flash - B_furnace_flash + Y_flash - Y_furnace = 0
_C47: R_furnace_flash - 1000 Wx_furnace_flash <= 0
_C48: L_furnace_flash + 1000 Wx_furnace_flash <= 1000
_C49: A_furnace_flash - 1000 Wy_furnace_flash <= 0
_C5: 0 O_reactor + l_reactor = 2
_C50: B_furnace_flash + 1000 Wy_furnace_flash <= 1000
_C51: - A_furnace_flash - B_furnace_flash + D_furnace_flash - L_furnace_flash
 - R_furnace_flash = 0
_C52: - D_flash_furnace + D_furnace_flash = 0
_C53: 1000 E1_furnace_flash + 1000 E2_furnace_flash - X_flash + X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= 0
_C54: - 1000 E1_furnace_flash + 1000 E2_furnace_flash + X_flash - X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= -1000
_C55: 1000 E1_furnace_flash - 1000 E2_furnace_flash - Y_flash + Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -1000
_C56: - 1000 E1_furnace_flash - 1000 E2_furnace_flash + Y_flash - Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -2000
_C57: L_furnace_flash - R_flash_furnace = 0
_C58: - L_flash_furnace + R_furnace_flash = 0
_C59: A_furnace_flash - B_flash_furnace = 0
_C6: d_reactor + l_reactor = 4
_C60: - A_flash_furnace + B_furnace_flash = 0
_C61: Wx_flash_furnace + Wx_furnace_flash = 1
_C62: Wy_flash_furnace + Wy_furnace_flash = 1
_C63: CD_reactor_flash - 414.065742323 D_reactor_flash = 0
_C64: - L_reactor_flash + R_reactor_flash + X_flash - X_reactor = 0
_C65: A_reactor_flash - B_reactor_flash + Y_flash - Y_reactor = 0
_C66: R_reactor_flash - 1000 Wx_reactor_flash <= 0
_C67: L_reactor_flash + 1000 Wx_reactor_flash <= 1000
_C68: A_reactor_flash - 1000 Wy_reactor_flash <= 0
_C69: B_reactor_flash + 1000 Wy_reactor_flash <= 1000
_C7: X_reactor - 0.5 l_reactor >= 0
_C70: - A_reactor_flash - B_reactor_flash + D_reactor_flash - L_reactor_flash
 - R_reactor_flash = 0
_C71: - D_flash_reactor + D_reactor_flash = 0
_C72: 1000 E1_reactor_flash + 1000 E2_reactor_flash - X_flash + X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= 0
_C73: - 1000 E1_reactor_flash + 1000 E2_reactor_flash + X_flash - X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= -1000
_C74: 1000 E1_reactor_flash - 1000 E2_reactor_flash - Y_flash + Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -1000
_C75: - 1000 E1_reactor_flash - 1000 E2_reactor_flash + Y_flash - Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -2000
_C76: L_reactor_flash - R_flash_reactor = 0
_C77: - L_flash_reactor + R_reactor_flash = 0
_C78: A_reactor_flash - B_flash_reactor = 0
_C79: - A_flash_reactor + B_reactor_flash = 0
_C8: Y_reactor - 0.5 d_reactor >= 0
_C80: Wx_flash_reactor + Wx_reactor_flash = 1
_C81: Wy_flash_reactor + Wy_reactor_flash = 1
_C82: CD_furnace_comp - 474.372452168 D_furnace_comp = 0
_C83: - L_furnace_comp + R_furnace_comp + X_comp - X_furnace = 0
_C84: A_furnace_comp - B_furnace_comp + Y_comp - Y_furnace = 0
_C85: R_furnace_comp - 1000 Wx_furnace_comp <= 0
_C86: L_furnace_comp + 1000 Wx_furnace_comp <= 1000
_C87: A_furnace_comp - 1000 Wy_furnace_comp <= 0
_C88: B_furnace_comp + 1000 Wy_furnace_comp <= 1000
_C89: - A_furnace_comp - B_furnace_comp + D_furnace_comp - L_furnace_comp
 - R_furnace_comp = 0
_C9: 0 O_flash + l_flash = 2.3
_C90: - D_comp_furnace + D_furnace_comp = 0
_C91: 1000 E1_furnace_comp + 1000 E2_furnace_comp - X_comp + X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= 0
_C92: - 1000 E1_furnace_comp + 1000 E2_furnace_comp + X_comp - X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= -1000
_C93: 1000 E1_furnace_comp - 1000 E2_furnace_comp - Y_comp + Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -1000
_C94: - 1000 E1_furnace_comp - 1000 E2_furnace_comp + Y_comp - Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -2000
_C95: L_furnace_comp - R_comp_furnace = 0
_C96: - L_comp_furnace + R_furnace_comp = 0
_C97: A_furnace_comp - B_comp_furnace = 0
_C98: - A_comp_furnace + B_furnace_comp = 0
_C99: Wx_comp_furnace + Wx_furnace_comp = 1
Bounds
 X_comp free
 X_distil free
 X_flash free
 X_furnace free
 X_reactor free
 X_store = 11.4
 Y_comp free
 Y_distil free
 Y_flash free
 Y_furnace free
 Y_reactor free
 Y_store = 11.4
Binaries
E1_comp_distil
E1_comp_store
E1_distil_store
E1_flash_comp
E1_flash_distil
E1_flash_store
E1_furnace_comp
E1_furnace_distil
E1_furnace_flash
E1_furnace_reactor
E1_furnace_store
E1_reactor_comp
E1_reactor_distil
E1_reactor_flash
E1_reactor_store
E2_comp_distil
E2_comp_store
E2_distil_store
E2_flash_comp
E2_flash_distil
E2_flash_store
E2_furnace_comp
E2_furnace_distil
E2_furnace_flash
E2_furnace_reactor
E2_furnace_store
E2_reactor_comp
E2_reactor_distil
E2_reactor_flash
E2_reactor_store
Gn_1
Gn_10
Gn_11
Gn_12
Gn_13
Gn_14
Gn_15
Gn_16
Gn_17
Gn_18
Gn_19
Gn_2
Gn_20
Gn_21
Gn_22
Gn_23
Gn_24
Gn_25
Gn_26
Gn_27
Gn_28
Gn_29
Gn_3
Gn_30
Gn_31
Gn_32
Gn_33
Gn_34
Gn_35
Gn_36
Gn_37
Gn_38
Gn_39
Gn_4
Gn_40
Gn_41
Gn_42
Gn_43
Gn_44
Gn_45
Gn_46
Gn_47
Gn_48
Gn_49
Gn_5
Gn_50
Gn_51
Gn_52
Gn_53
Gn_54
Gn_55
Gn_56
Gn_57
Gn_58
Gn_59
Gn_6
Gn_7
Gn_8
Gn_9
O_comp
O_distil
O_flash
O_furnace
O_reactor
O_store
Psie_distil_comp
Psie_distil_flash
Psie_distil_furnace
Psie_distil_reactor
Psie_distil_store
Psie_furnace_comp
Psie_furnace_distil
Psie_furnace_flash
Psie_furnace_reactor
Psie_furnace_store
Psie_reactor_comp
Psie_reactor_distil
Psie_reactor_flash
Psie_reactor_furnace
Psie_reactor_store
Psie_store_comp
Psie_store_distil
Psie_store_flash
Psie_store_furnace
Psie_store_reactor
Wx_comp_distil
Wx_comp_flash
Wx_comp_furnace
Wx_comp_reactor
Wx_comp_store
Wx_distil_comp
Wx_distil_flash
Wx_distil_furnace
Wx_distil_reactor
Wx_distil_store
Wx_flash_comp
Wx_flash_distil
Wx_flash_furnace
Wx_flash_reactor
Wx_flash_store
Wx_furnace_comp
Wx_furnace_distil
Wx_furnace_flash
Wx_furnace_reactor
Wx_furnace_store
Wx_reactor_comp
Wx_reactor_distil
Wx_reactor_flash
Wx_reactor_furnace
Wx_reactor_store
Wx_store_comp
Wx_store_distil
Wx_store_flash
Wx_store_furnace
Wx_store_reactor
Wy_comp_distil
Wy_comp_flash
Wy_comp_furnace
Wy_comp_reactor
Wy_comp_store
Wy_distil_comp
Wy_distil_flash
Wy_distil_furnace
Wy_distil_reactor
Wy_distil_store
Wy_flash_comp
Wy_flash_distil
Wy_flash_furnace
Wy_flash_reactor
Wy_flash_store
Wy_furnace_comp
Wy_furnace_distil
Wy_furnace_flash
Wy_furnace_reactor
Wy_furnace_store
Wy_reactor_comp
Wy_reactor_distil
Wy_reactor_flash
Wy_reactor_furnace
Wy_reactor_store
Wy_store_comp
Wy_store_distil
Wy_store_flash
Wy_store_furnace
Wy_store_reactor
End
