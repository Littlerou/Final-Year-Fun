\* Layout_Problem_Model *\
Minimize
OBJ: SumCD + SumOmega
Subject To
_C1: 9 O_furnace + l_furnace = 16
_C10: d_flash + l_flash = 8
_C100: L_furnace_comp + 1000 Wx_furnace_comp <= 1000
_C101: A_furnace_comp - 1000 Wy_furnace_comp <= 0
_C102: B_furnace_comp + 1000 Wy_furnace_comp <= 1000
_C103: - A_furnace_comp - B_furnace_comp + D_furnace_comp - L_furnace_comp
 - R_furnace_comp = 0
_C104: - D_comp_furnace + D_furnace_comp = 0
_C105: 1000 E1_furnace_comp + 1000 E2_furnace_comp - X_comp + X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= 0
_C106: - 1000 E1_furnace_comp + 1000 E2_furnace_comp + X_comp - X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= -1000
_C107: 1000 E1_furnace_comp - 1000 E2_furnace_comp - Y_comp + Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -1000
_C108: - 1000 E1_furnace_comp - 1000 E2_furnace_comp + Y_comp - Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -2000
_C109: L_furnace_comp - R_comp_furnace = 0
_C11: X_flash - 0.5 l_flash >= 0
_C110: - L_comp_furnace + R_furnace_comp = 0
_C111: A_furnace_comp - B_comp_furnace = 0
_C112: - A_comp_furnace + B_furnace_comp = 0
_C113: Wx_comp_furnace + Wx_furnace_comp = 1
_C114: Wy_comp_furnace + Wy_furnace_comp = 1
_C115: CD_reactor_comp = 0
_C116: - L_reactor_comp + R_reactor_comp + X_comp - X_reactor = 0
_C117: A_reactor_comp - B_reactor_comp + Y_comp - Y_reactor = 0
_C118: R_reactor_comp - 1000 Wx_reactor_comp <= 0
_C119: L_reactor_comp + 1000 Wx_reactor_comp <= 1000
_C12: Y_flash - 0.5 d_flash >= 0
_C120: A_reactor_comp - 1000 Wy_reactor_comp <= 0
_C121: B_reactor_comp + 1000 Wy_reactor_comp <= 1000
_C122: - A_reactor_comp - B_reactor_comp + D_reactor_comp - L_reactor_comp
 - R_reactor_comp = 0
_C123: - D_comp_reactor + D_reactor_comp = 0
_C124: 1000 E1_reactor_comp + 1000 E2_reactor_comp - X_comp + X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= 0
_C125: - 1000 E1_reactor_comp + 1000 E2_reactor_comp + X_comp - X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= -1000
_C126: 1000 E1_reactor_comp - 1000 E2_reactor_comp - Y_comp + Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -1000
_C127: - 1000 E1_reactor_comp - 1000 E2_reactor_comp + Y_comp - Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -2000
_C128: L_reactor_comp - R_comp_reactor = 0
_C129: - L_comp_reactor + R_reactor_comp = 0
_C13: 0 O_comp + l_comp = 14
_C130: A_reactor_comp - B_comp_reactor = 0
_C131: - A_comp_reactor + B_reactor_comp = 0
_C132: Wx_comp_reactor + Wx_reactor_comp = 1
_C133: Wy_comp_reactor + Wy_reactor_comp = 1
_C134: CD_flash_comp - 842.70764367 D_flash_comp = 0
_C135: - L_flash_comp + R_flash_comp + X_comp - X_flash = 0
_C136: A_flash_comp - B_flash_comp + Y_comp - Y_flash = 0
_C137: R_flash_comp - 1000 Wx_flash_comp <= 0
_C138: L_flash_comp + 1000 Wx_flash_comp <= 1000
_C139: A_flash_comp - 1000 Wy_flash_comp <= 0
_C14: d_comp + l_comp = 28
_C140: B_flash_comp + 1000 Wy_flash_comp <= 1000
_C141: - A_flash_comp - B_flash_comp + D_flash_comp - L_flash_comp
 - R_flash_comp = 0
_C142: - D_comp_flash + D_flash_comp = 0
_C143: 1000 E1_flash_comp + 1000 E2_flash_comp - X_comp + X_flash - 0.5 l_comp
 - 0.5 l_flash >= 0
_C144: - 1000 E1_flash_comp + 1000 E2_flash_comp + X_comp - X_flash
 - 0.5 l_comp - 0.5 l_flash >= -1000
_C145: 1000 E1_flash_comp - 1000 E2_flash_comp - Y_comp + Y_flash - 0.5 d_comp
 - 0.5 d_flash >= -1000
_C146: - 1000 E1_flash_comp - 1000 E2_flash_comp + Y_comp - Y_flash
 - 0.5 d_comp - 0.5 d_flash >= -2000
_C147: L_flash_comp - R_comp_flash = 0
_C148: - L_comp_flash + R_flash_comp = 0
_C149: A_flash_comp - B_comp_flash = 0
_C15: X_comp - 0.5 l_comp >= 0
_C150: - A_comp_flash + B_flash_comp = 0
_C151: Wx_comp_flash + Wx_flash_comp = 1
_C152: Wy_comp_flash + Wy_flash_comp = 1
_C153: CD_furnace_distil = 0
_C154: - L_furnace_distil + R_furnace_distil + X_distil - X_furnace = 0
_C155: A_furnace_distil - B_furnace_distil + Y_distil - Y_furnace = 0
_C156: R_furnace_distil - 1000 Wx_furnace_distil <= 0
_C157: L_furnace_distil + 1000 Wx_furnace_distil <= 1000
_C158: A_furnace_distil - 1000 Wy_furnace_distil <= 0
_C159: B_furnace_distil + 1000 Wy_furnace_distil <= 1000
_C16: Y_comp - 0.5 d_comp >= 0
_C160: - A_furnace_distil - B_furnace_distil + D_furnace_distil
 - L_furnace_distil - R_furnace_distil = 0
_C161: - D_distil_furnace + D_furnace_distil = 0
_C162: 1000 E1_furnace_distil + 1000 E2_furnace_distil - X_distil + X_furnace
 - 0.5 l_distil - 0.5 l_furnace >= 0
_C163: - 1000 E1_furnace_distil + 1000 E2_furnace_distil + X_distil
 - X_furnace - 0.5 l_distil - 0.5 l_furnace >= -1000
_C164: 1000 E1_furnace_distil - 1000 E2_furnace_distil - Y_distil + Y_furnace
 - 0.5 d_distil - 0.5 d_furnace >= -1000
_C165: - 1000 E1_furnace_distil - 1000 E2_furnace_distil + Y_distil
 - Y_furnace - 0.5 d_distil - 0.5 d_furnace >= -2000
_C166: L_furnace_distil - R_distil_furnace = 0
_C167: - L_distil_furnace + R_furnace_distil = 0
_C168: A_furnace_distil - B_distil_furnace = 0
_C169: - A_distil_furnace + B_furnace_distil = 0
_C17: 0 O_distil + l_distil = 8
_C170: Wx_distil_furnace + Wx_furnace_distil = 1
_C171: Wy_distil_furnace + Wy_furnace_distil = 1
_C172: CD_reactor_distil = 0
_C173: - L_reactor_distil + R_reactor_distil + X_distil - X_reactor = 0
_C174: A_reactor_distil - B_reactor_distil + Y_distil - Y_reactor = 0
_C175: R_reactor_distil - 1000 Wx_reactor_distil <= 0
_C176: L_reactor_distil + 1000 Wx_reactor_distil <= 1000
_C177: A_reactor_distil - 1000 Wy_reactor_distil <= 0
_C178: B_reactor_distil + 1000 Wy_reactor_distil <= 1000
_C179: - A_reactor_distil - B_reactor_distil + D_reactor_distil
 - L_reactor_distil - R_reactor_distil = 0
_C18: d_distil + l_distil = 16
_C180: - D_distil_reactor + D_reactor_distil = 0
_C181: 1000 E1_reactor_distil + 1000 E2_reactor_distil - X_distil + X_reactor
 - 0.5 l_distil - 0.5 l_reactor >= 0
_C182: - 1000 E1_reactor_distil + 1000 E2_reactor_distil + X_distil
 - X_reactor - 0.5 l_distil - 0.5 l_reactor >= -1000
_C183: 1000 E1_reactor_distil - 1000 E2_reactor_distil - Y_distil + Y_reactor
 - 0.5 d_distil - 0.5 d_reactor >= -1000
_C184: - 1000 E1_reactor_distil - 1000 E2_reactor_distil + Y_distil
 - Y_reactor - 0.5 d_distil - 0.5 d_reactor >= -2000
_C185: L_reactor_distil - R_distil_reactor = 0
_C186: - L_distil_reactor + R_reactor_distil = 0
_C187: A_reactor_distil - B_distil_reactor = 0
_C188: - A_distil_reactor + B_reactor_distil = 0
_C189: Wx_distil_reactor + Wx_reactor_distil = 1
_C19: X_distil - 0.5 l_distil >= 0
_C190: Wy_distil_reactor + Wy_reactor_distil = 1
_C191: CD_flash_distil - 78.4818701125 D_flash_distil = 0
_C192: - L_flash_distil + R_flash_distil + X_distil - X_flash = 0
_C193: A_flash_distil - B_flash_distil + Y_distil - Y_flash = 0
_C194: R_flash_distil - 1000 Wx_flash_distil <= 0
_C195: L_flash_distil + 1000 Wx_flash_distil <= 1000
_C196: A_flash_distil - 1000 Wy_flash_distil <= 0
_C197: B_flash_distil + 1000 Wy_flash_distil <= 1000
_C198: - A_flash_distil - B_flash_distil + D_flash_distil - L_flash_distil
 - R_flash_distil = 0
_C199: - D_distil_flash + D_flash_distil = 0
_C2: d_furnace + l_furnace = 23
_C20: Y_distil - 0.5 d_distil >= 0
_C200: 1000 E1_flash_distil + 1000 E2_flash_distil - X_distil + X_flash
 - 0.5 l_distil - 0.5 l_flash >= 0
_C201: - 1000 E1_flash_distil + 1000 E2_flash_distil + X_distil - X_flash
 - 0.5 l_distil - 0.5 l_flash >= -1000
_C202: 1000 E1_flash_distil - 1000 E2_flash_distil - Y_distil + Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -1000
_C203: - 1000 E1_flash_distil - 1000 E2_flash_distil + Y_distil - Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -2000
_C204: L_flash_distil - R_distil_flash = 0
_C205: - L_distil_flash + R_flash_distil = 0
_C206: A_flash_distil - B_distil_flash = 0
_C207: - A_distil_flash + B_flash_distil = 0
_C208: Wx_distil_flash + Wx_flash_distil = 1
_C209: Wy_distil_flash + Wy_flash_distil = 1
_C21: - 10 O_store + l_store = 12
_C210: CD_comp_distil = 0
_C211: - L_comp_distil + R_comp_distil - X_comp + X_distil = 0
_C212: A_comp_distil - B_comp_distil - Y_comp + Y_distil = 0
_C213: R_comp_distil - 1000 Wx_comp_distil <= 0
_C214: L_comp_distil + 1000 Wx_comp_distil <= 1000
_C215: A_comp_distil - 1000 Wy_comp_distil <= 0
_C216: B_comp_distil + 1000 Wy_comp_distil <= 1000
_C217: - A_comp_distil - B_comp_distil + D_comp_distil - L_comp_distil
 - R_comp_distil = 0
_C218: D_comp_distil - D_distil_comp = 0
_C219: 1000 E1_comp_distil + 1000 E2_comp_distil + X_comp - X_distil
 - 0.5 l_comp - 0.5 l_distil >= 0
_C22: d_store + l_store = 34
_C220: - 1000 E1_comp_distil + 1000 E2_comp_distil - X_comp + X_distil
 - 0.5 l_comp - 0.5 l_distil >= -1000
_C221: 1000 E1_comp_distil - 1000 E2_comp_distil + Y_comp - Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -1000
_C222: - 1000 E1_comp_distil - 1000 E2_comp_distil - Y_comp + Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -2000
_C223: L_comp_distil - R_distil_comp = 0
_C224: - L_distil_comp + R_comp_distil = 0
_C225: A_comp_distil - B_distil_comp = 0
_C226: - A_distil_comp + B_comp_distil = 0
_C227: Wx_comp_distil + Wx_distil_comp = 1
_C228: Wy_comp_distil + Wy_distil_comp = 1
_C229: CD_furnace_store - 75.8045442191 D_furnace_store = 0
_C23: X_store - 0.5 l_store >= 0
_C230: - L_furnace_store + R_furnace_store - X_furnace + X_store = 0
_C231: A_furnace_store - B_furnace_store - Y_furnace + Y_store = 0
_C232: R_furnace_store - 1000 Wx_furnace_store <= 0
_C233: L_furnace_store + 1000 Wx_furnace_store <= 1000
_C234: A_furnace_store - 1000 Wy_furnace_store <= 0
_C235: B_furnace_store + 1000 Wy_furnace_store <= 1000
_C236: - A_furnace_store - B_furnace_store + D_furnace_store - L_furnace_store
 - R_furnace_store = 0
_C237: D_furnace_store - D_store_furnace = 0
_C238: 1000 E1_furnace_store + 1000 E2_furnace_store + X_furnace - X_store
 - 0.5 l_furnace - 0.5 l_store >= 0
_C239: - 1000 E1_furnace_store + 1000 E2_furnace_store - X_furnace + X_store
 - 0.5 l_furnace - 0.5 l_store >= -1000
_C24: Y_store - 0.5 d_store >= 0
_C240: 1000 E1_furnace_store - 1000 E2_furnace_store + Y_furnace - Y_store
 - 0.5 d_furnace - 0.5 d_store >= -1000
_C241: - 1000 E1_furnace_store - 1000 E2_furnace_store - Y_furnace + Y_store
 - 0.5 d_furnace - 0.5 d_store >= -2000
_C242: L_furnace_store - R_store_furnace = 0
_C243: - L_store_furnace + R_furnace_store = 0
_C244: A_furnace_store - B_store_furnace = 0
_C245: - A_store_furnace + B_furnace_store = 0
_C246: Wx_furnace_store + Wx_store_furnace = 1
_C247: Wy_furnace_store + Wy_store_furnace = 1
_C248: CD_reactor_store = 0
_C249: - L_reactor_store + R_reactor_store - X_reactor + X_store = 0
_C25: 0 O_ctrlroom + l_ctrlroom = 15
_C250: A_reactor_store - B_reactor_store - Y_reactor + Y_store = 0
_C251: R_reactor_store - 1000 Wx_reactor_store <= 0
_C252: L_reactor_store + 1000 Wx_reactor_store <= 1000
_C253: A_reactor_store - 1000 Wy_reactor_store <= 0
_C254: B_reactor_store + 1000 Wy_reactor_store <= 1000
_C255: - A_reactor_store - B_reactor_store + D_reactor_store - L_reactor_store
 - R_reactor_store = 0
_C256: D_reactor_store - D_store_reactor = 0
_C257: 1000 E1_reactor_store + 1000 E2_reactor_store + X_reactor - X_store
 - 0.5 l_reactor - 0.5 l_store >= 0
_C258: - 1000 E1_reactor_store + 1000 E2_reactor_store - X_reactor + X_store
 - 0.5 l_reactor - 0.5 l_store >= -1000
_C259: 1000 E1_reactor_store - 1000 E2_reactor_store + Y_reactor - Y_store
 - 0.5 d_reactor - 0.5 d_store >= -1000
_C26: d_ctrlroom + l_ctrlroom = 30
_C260: - 1000 E1_reactor_store - 1000 E2_reactor_store - Y_reactor + Y_store
 - 0.5 d_reactor - 0.5 d_store >= -2000
_C261: L_reactor_store - R_store_reactor = 0
_C262: - L_store_reactor + R_reactor_store = 0
_C263: A_reactor_store - B_store_reactor = 0
_C264: - A_store_reactor + B_reactor_store = 0
_C265: Wx_reactor_store + Wx_store_reactor = 1
_C266: Wy_reactor_store + Wy_store_reactor = 1
_C267: CD_flash_store = 0
_C268: - L_flash_store + R_flash_store - X_flash + X_store = 0
_C269: A_flash_store - B_flash_store - Y_flash + Y_store = 0
_C27: X_ctrlroom - 0.5 l_ctrlroom >= 0
_C270: R_flash_store - 1000 Wx_flash_store <= 0
_C271: L_flash_store + 1000 Wx_flash_store <= 1000
_C272: A_flash_store - 1000 Wy_flash_store <= 0
_C273: B_flash_store + 1000 Wy_flash_store <= 1000
_C274: - A_flash_store - B_flash_store + D_flash_store - L_flash_store
 - R_flash_store = 0
_C275: D_flash_store - D_store_flash = 0
_C276: 1000 E1_flash_store + 1000 E2_flash_store + X_flash - X_store
 - 0.5 l_flash - 0.5 l_store >= 0
_C277: - 1000 E1_flash_store + 1000 E2_flash_store - X_flash + X_store
 - 0.5 l_flash - 0.5 l_store >= -1000
_C278: 1000 E1_flash_store - 1000 E2_flash_store + Y_flash - Y_store
 - 0.5 d_flash - 0.5 d_store >= -1000
_C279: - 1000 E1_flash_store - 1000 E2_flash_store - Y_flash + Y_store
 - 0.5 d_flash - 0.5 d_store >= -2000
_C28: Y_ctrlroom - 0.5 d_ctrlroom >= 0
_C280: L_flash_store - R_store_flash = 0
_C281: - L_store_flash + R_flash_store = 0
_C282: A_flash_store - B_store_flash = 0
_C283: - A_store_flash + B_flash_store = 0
_C284: Wx_flash_store + Wx_store_flash = 1
_C285: Wy_flash_store + Wy_store_flash = 1
_C286: CD_comp_store = 0
_C287: - L_comp_store + R_comp_store - X_comp + X_store = 0
_C288: A_comp_store - B_comp_store - Y_comp + Y_store = 0
_C289: R_comp_store - 1000 Wx_comp_store <= 0
_C29: - 0.166666666667 X_comp - 0.166666666667 X_distil
 - 0.166666666667 X_flash - 0.166666666667 X_furnace
 - 0.166666666667 X_reactor - 0.166666666667 X_store + avgX = 0
_C290: L_comp_store + 1000 Wx_comp_store <= 1000
_C291: A_comp_store - 1000 Wy_comp_store <= 0
_C292: B_comp_store + 1000 Wy_comp_store <= 1000
_C293: - A_comp_store - B_comp_store + D_comp_store - L_comp_store
 - R_comp_store = 0
_C294: D_comp_store - D_store_comp = 0
_C295: 1000 E1_comp_store + 1000 E2_comp_store + X_comp - X_store - 0.5 l_comp
 - 0.5 l_store >= 0
_C296: - 1000 E1_comp_store + 1000 E2_comp_store - X_comp + X_store
 - 0.5 l_comp - 0.5 l_store >= -1000
_C297: 1000 E1_comp_store - 1000 E2_comp_store + Y_comp - Y_store - 0.5 d_comp
 - 0.5 d_store >= -1000
_C298: - 1000 E1_comp_store - 1000 E2_comp_store - Y_comp + Y_store
 - 0.5 d_comp - 0.5 d_store >= -2000
_C299: L_comp_store - R_store_comp = 0
_C3: X_furnace - 0.5 l_furnace >= 0
_C30: - 0.166666666667 Y_comp - 0.166666666667 Y_distil
 - 0.166666666667 Y_flash - 0.166666666667 Y_furnace
 - 0.166666666667 Y_reactor - 0.166666666667 Y_store + avgY = 0
_C300: - L_store_comp + R_comp_store = 0
_C301: A_comp_store - B_store_comp = 0
_C302: - A_store_comp + B_comp_store = 0
_C303: Wx_comp_store + Wx_store_comp = 1
_C304: Wy_comp_store + Wy_store_comp = 1
_C305: CD_distil_store - 451.755491306 D_distil_store = 0
_C306: - L_distil_store + R_distil_store - X_distil + X_store = 0
_C307: A_distil_store - B_distil_store - Y_distil + Y_store = 0
_C308: R_distil_store - 1000 Wx_distil_store <= 0
_C309: L_distil_store + 1000 Wx_distil_store <= 1000
_C31: - L2 + R2 - X_ctrlroom + avgX = 0
_C310: A_distil_store - 1000 Wy_distil_store <= 0
_C311: B_distil_store + 1000 Wy_distil_store <= 1000
_C312: - A_distil_store - B_distil_store + D_distil_store - L_distil_store
 - R_distil_store = 0
_C313: D_distil_store - D_store_distil = 0
_C314: 1000 E1_distil_store + 1000 E2_distil_store + X_distil - X_store
 - 0.5 l_distil - 0.5 l_store >= 0
_C315: - 1000 E1_distil_store + 1000 E2_distil_store - X_distil + X_store
 - 0.5 l_distil - 0.5 l_store >= -1000
_C316: 1000 E1_distil_store - 1000 E2_distil_store + Y_distil - Y_store
 - 0.5 d_distil - 0.5 d_store >= -1000
_C317: - 1000 E1_distil_store - 1000 E2_distil_store - Y_distil + Y_store
 - 0.5 d_distil - 0.5 d_store >= -2000
_C318: L_distil_store - R_store_distil = 0
_C319: - L_store_distil + R_distil_store = 0
_C32: A2 - B2 - Y_ctrlroom + avgY = 0
_C320: A_distil_store - B_store_distil = 0
_C321: - A_store_distil + B_distil_store = 0
_C322: Wx_distil_store + Wx_store_distil = 1
_C323: Wy_distil_store + Wy_store_distil = 1
_C324: CD_furnace_ctrlroom = 0
_C325: - L_furnace_ctrlroom + R_furnace_ctrlroom + X_ctrlroom - X_furnace = 0
_C326: A_furnace_ctrlroom - B_furnace_ctrlroom + Y_ctrlroom - Y_furnace = 0
_C327: R_furnace_ctrlroom - 1000 Wx_furnace_ctrlroom <= 0
_C328: L_furnace_ctrlroom + 1000 Wx_furnace_ctrlroom <= 1000
_C329: A_furnace_ctrlroom - 1000 Wy_furnace_ctrlroom <= 0
_C33: R2 - 1000 Wx2 <= 0
_C330: B_furnace_ctrlroom + 1000 Wy_furnace_ctrlroom <= 1000
_C331: - A_furnace_ctrlroom - B_furnace_ctrlroom + D_furnace_ctrlroom
 - L_furnace_ctrlroom - R_furnace_ctrlroom = 0
_C332: - D_ctrlroom_furnace + D_furnace_ctrlroom = 0
_C333: 1000 E1_furnace_ctrlroom + 1000 E2_furnace_ctrlroom - X_ctrlroom
 + X_furnace - 0.5 l_ctrlroom - 0.5 l_furnace >= 0
_C334: - 1000 E1_furnace_ctrlroom + 1000 E2_furnace_ctrlroom + X_ctrlroom
 - X_furnace - 0.5 l_ctrlroom - 0.5 l_furnace >= -1000
_C335: 1000 E1_furnace_ctrlroom - 1000 E2_furnace_ctrlroom - Y_ctrlroom
 + Y_furnace - 0.5 d_ctrlroom - 0.5 d_furnace >= -1000
_C336: - 1000 E1_furnace_ctrlroom - 1000 E2_furnace_ctrlroom + Y_ctrlroom
 - Y_furnace - 0.5 d_ctrlroom - 0.5 d_furnace >= -2000
_C337: L_furnace_ctrlroom - R_ctrlroom_furnace = 0
_C338: - L_ctrlroom_furnace + R_furnace_ctrlroom = 0
_C339: A_furnace_ctrlroom - B_ctrlroom_furnace = 0
_C34: L2 + 1000 Wx2 <= 1000
_C340: - A_ctrlroom_furnace + B_furnace_ctrlroom = 0
_C341: Wx_ctrlroom_furnace + Wx_furnace_ctrlroom = 1
_C342: Wy_ctrlroom_furnace + Wy_furnace_ctrlroom = 1
_C343: CD_reactor_ctrlroom = 0
_C344: - L_reactor_ctrlroom + R_reactor_ctrlroom + X_ctrlroom - X_reactor = 0
_C345: A_reactor_ctrlroom - B_reactor_ctrlroom + Y_ctrlroom - Y_reactor = 0
_C346: R_reactor_ctrlroom - 1000 Wx_reactor_ctrlroom <= 0
_C347: L_reactor_ctrlroom + 1000 Wx_reactor_ctrlroom <= 1000
_C348: A_reactor_ctrlroom - 1000 Wy_reactor_ctrlroom <= 0
_C349: B_reactor_ctrlroom + 1000 Wy_reactor_ctrlroom <= 1000
_C35: A2 - 1000 Wy2 <= 0
_C350: - A_reactor_ctrlroom - B_reactor_ctrlroom + D_reactor_ctrlroom
 - L_reactor_ctrlroom - R_reactor_ctrlroom = 0
_C351: - D_ctrlroom_reactor + D_reactor_ctrlroom = 0
_C352: 1000 E1_reactor_ctrlroom + 1000 E2_reactor_ctrlroom - X_ctrlroom
 + X_reactor - 0.5 l_ctrlroom - 0.5 l_reactor >= 0
_C353: - 1000 E1_reactor_ctrlroom + 1000 E2_reactor_ctrlroom + X_ctrlroom
 - X_reactor - 0.5 l_ctrlroom - 0.5 l_reactor >= -1000
_C354: 1000 E1_reactor_ctrlroom - 1000 E2_reactor_ctrlroom - Y_ctrlroom
 + Y_reactor - 0.5 d_ctrlroom - 0.5 d_reactor >= -1000
_C355: - 1000 E1_reactor_ctrlroom - 1000 E2_reactor_ctrlroom + Y_ctrlroom
 - Y_reactor - 0.5 d_ctrlroom - 0.5 d_reactor >= -2000
_C356: L_reactor_ctrlroom - R_ctrlroom_reactor = 0
_C357: - L_ctrlroom_reactor + R_reactor_ctrlroom = 0
_C358: A_reactor_ctrlroom - B_ctrlroom_reactor = 0
_C359: - A_ctrlroom_reactor + B_reactor_ctrlroom = 0
_C36: B2 + 1000 Wy2 <= 1000
_C360: Wx_ctrlroom_reactor + Wx_reactor_ctrlroom = 1
_C361: Wy_ctrlroom_reactor + Wy_reactor_ctrlroom = 1
_C362: CD_flash_ctrlroom = 0
_C363: - L_flash_ctrlroom + R_flash_ctrlroom + X_ctrlroom - X_flash = 0
_C364: A_flash_ctrlroom - B_flash_ctrlroom + Y_ctrlroom - Y_flash = 0
_C365: R_flash_ctrlroom - 1000 Wx_flash_ctrlroom <= 0
_C366: L_flash_ctrlroom + 1000 Wx_flash_ctrlroom <= 1000
_C367: A_flash_ctrlroom - 1000 Wy_flash_ctrlroom <= 0
_C368: B_flash_ctrlroom + 1000 Wy_flash_ctrlroom <= 1000
_C369: - A_flash_ctrlroom - B_flash_ctrlroom + D_flash_ctrlroom
 - L_flash_ctrlroom - R_flash_ctrlroom = 0
_C37: - A2 - B2 + D_road - L2 - R2 = 0
_C370: - D_ctrlroom_flash + D_flash_ctrlroom = 0
_C371: 1000 E1_flash_ctrlroom + 1000 E2_flash_ctrlroom - X_ctrlroom + X_flash
 - 0.5 l_ctrlroom - 0.5 l_flash >= 0
_C372: - 1000 E1_flash_ctrlroom + 1000 E2_flash_ctrlroom + X_ctrlroom
 - X_flash - 0.5 l_ctrlroom - 0.5 l_flash >= -1000
_C373: 1000 E1_flash_ctrlroom - 1000 E2_flash_ctrlroom - Y_ctrlroom + Y_flash
 - 0.5 d_ctrlroom - 0.5 d_flash >= -1000
_C374: - 1000 E1_flash_ctrlroom - 1000 E2_flash_ctrlroom + Y_ctrlroom
 - Y_flash - 0.5 d_ctrlroom - 0.5 d_flash >= -2000
_C375: L_flash_ctrlroom - R_ctrlroom_flash = 0
_C376: - L_ctrlroom_flash + R_flash_ctrlroom = 0
_C377: A_flash_ctrlroom - B_ctrlroom_flash = 0
_C378: - A_ctrlroom_flash + B_flash_ctrlroom = 0
_C379: Wx_ctrlroom_flash + Wx_flash_ctrlroom = 1
_C38: - 120 D_road + roadCost = 0
_C380: Wy_ctrlroom_flash + Wy_flash_ctrlroom = 1
_C381: CD_comp_ctrlroom = 0
_C382: - L_comp_ctrlroom + R_comp_ctrlroom - X_comp + X_ctrlroom = 0
_C383: A_comp_ctrlroom - B_comp_ctrlroom - Y_comp + Y_ctrlroom = 0
_C384: R_comp_ctrlroom - 1000 Wx_comp_ctrlroom <= 0
_C385: L_comp_ctrlroom + 1000 Wx_comp_ctrlroom <= 1000
_C386: A_comp_ctrlroom - 1000 Wy_comp_ctrlroom <= 0
_C387: B_comp_ctrlroom + 1000 Wy_comp_ctrlroom <= 1000
_C388: - A_comp_ctrlroom - B_comp_ctrlroom + D_comp_ctrlroom - L_comp_ctrlroom
 - R_comp_ctrlroom = 0
_C389: D_comp_ctrlroom - D_ctrlroom_comp = 0
_C39: CD_furnace_reactor - 1075.48971174 D_furnace_reactor = 0
_C390: 1000 E1_comp_ctrlroom + 1000 E2_comp_ctrlroom + X_comp - X_ctrlroom
 - 0.5 l_comp - 0.5 l_ctrlroom >= 0
_C391: - 1000 E1_comp_ctrlroom + 1000 E2_comp_ctrlroom - X_comp + X_ctrlroom
 - 0.5 l_comp - 0.5 l_ctrlroom >= -1000
_C392: 1000 E1_comp_ctrlroom - 1000 E2_comp_ctrlroom + Y_comp - Y_ctrlroom
 - 0.5 d_comp - 0.5 d_ctrlroom >= -1000
_C393: - 1000 E1_comp_ctrlroom - 1000 E2_comp_ctrlroom - Y_comp + Y_ctrlroom
 - 0.5 d_comp - 0.5 d_ctrlroom >= -2000
_C394: L_comp_ctrlroom - R_ctrlroom_comp = 0
_C395: - L_ctrlroom_comp + R_comp_ctrlroom = 0
_C396: A_comp_ctrlroom - B_ctrlroom_comp = 0
_C397: - A_ctrlroom_comp + B_comp_ctrlroom = 0
_C398: Wx_comp_ctrlroom + Wx_ctrlroom_comp = 1
_C399: Wy_comp_ctrlroom + Wy_ctrlroom_comp = 1
_C4: Y_furnace - 0.5 d_furnace >= 0
_C40: - L_furnace_reactor + R_furnace_reactor - X_furnace + X_reactor = 0
_C400: CD_distil_ctrlroom = 0
_C401: - L_distil_ctrlroom + R_distil_ctrlroom + X_ctrlroom - X_distil = 0
_C402: A_distil_ctrlroom - B_distil_ctrlroom + Y_ctrlroom - Y_distil = 0
_C403: R_distil_ctrlroom - 1000 Wx_distil_ctrlroom <= 0
_C404: L_distil_ctrlroom + 1000 Wx_distil_ctrlroom <= 1000
_C405: A_distil_ctrlroom - 1000 Wy_distil_ctrlroom <= 0
_C406: B_distil_ctrlroom + 1000 Wy_distil_ctrlroom <= 1000
_C407: - A_distil_ctrlroom - B_distil_ctrlroom + D_distil_ctrlroom
 - L_distil_ctrlroom - R_distil_ctrlroom = 0
_C408: - D_ctrlroom_distil + D_distil_ctrlroom = 0
_C409: 1000 E1_distil_ctrlroom + 1000 E2_distil_ctrlroom - X_ctrlroom
 + X_distil - 0.5 l_ctrlroom - 0.5 l_distil >= 0
_C41: A_furnace_reactor - B_furnace_reactor - Y_furnace + Y_reactor = 0
_C410: - 1000 E1_distil_ctrlroom + 1000 E2_distil_ctrlroom + X_ctrlroom
 - X_distil - 0.5 l_ctrlroom - 0.5 l_distil >= -1000
_C411: 1000 E1_distil_ctrlroom - 1000 E2_distil_ctrlroom - Y_ctrlroom
 + Y_distil - 0.5 d_ctrlroom - 0.5 d_distil >= -1000
_C412: - 1000 E1_distil_ctrlroom - 1000 E2_distil_ctrlroom + Y_ctrlroom
 - Y_distil - 0.5 d_ctrlroom - 0.5 d_distil >= -2000
_C413: L_distil_ctrlroom - R_ctrlroom_distil = 0
_C414: - L_ctrlroom_distil + R_distil_ctrlroom = 0
_C415: A_distil_ctrlroom - B_ctrlroom_distil = 0
_C416: - A_ctrlroom_distil + B_distil_ctrlroom = 0
_C417: Wx_ctrlroom_distil + Wx_distil_ctrlroom = 1
_C418: Wy_ctrlroom_distil + Wy_distil_ctrlroom = 1
_C419: CD_store_ctrlroom = 0
_C42: R_furnace_reactor - 1000 Wx_furnace_reactor <= 0
_C420: - L_store_ctrlroom + R_store_ctrlroom + X_ctrlroom - X_store = 0
_C421: A_store_ctrlroom - B_store_ctrlroom + Y_ctrlroom - Y_store = 0
_C422: R_store_ctrlroom - 1000 Wx_store_ctrlroom <= 0
_C423: L_store_ctrlroom + 1000 Wx_store_ctrlroom <= 1000
_C424: A_store_ctrlroom - 1000 Wy_store_ctrlroom <= 0
_C425: B_store_ctrlroom + 1000 Wy_store_ctrlroom <= 1000
_C426: - A_store_ctrlroom - B_store_ctrlroom + D_store_ctrlroom
 - L_store_ctrlroom - R_store_ctrlroom = 0
_C427: - D_ctrlroom_store + D_store_ctrlroom = 0
_C428: 1000 E1_store_ctrlroom + 1000 E2_store_ctrlroom - X_ctrlroom + X_store
 - 0.5 l_ctrlroom - 0.5 l_store >= 0
_C429: - 1000 E1_store_ctrlroom + 1000 E2_store_ctrlroom + X_ctrlroom
 - X_store - 0.5 l_ctrlroom - 0.5 l_store >= -1000
_C43: L_furnace_reactor + 1000 Wx_furnace_reactor <= 1000
_C430: 1000 E1_store_ctrlroom - 1000 E2_store_ctrlroom - Y_ctrlroom + Y_store
 - 0.5 d_ctrlroom - 0.5 d_store >= -1000
_C431: - 1000 E1_store_ctrlroom - 1000 E2_store_ctrlroom + Y_ctrlroom
 - Y_store - 0.5 d_ctrlroom - 0.5 d_store >= -2000
_C432: L_store_ctrlroom - R_ctrlroom_store = 0
_C433: - L_ctrlroom_store + R_store_ctrlroom = 0
_C434: A_store_ctrlroom - B_ctrlroom_store = 0
_C435: - A_ctrlroom_store + B_store_ctrlroom = 0
_C436: Wx_ctrlroom_store + Wx_store_ctrlroom = 1
_C437: Wy_ctrlroom_store + Wy_store_ctrlroom = 1
_C438: - CD_comp_comp - CD_comp_ctrlroom - CD_comp_distil - CD_comp_flash
 - CD_comp_furnace - CD_comp_reactor - CD_comp_store - CD_ctrlroom_comp
 - CD_ctrlroom_ctrlroom - CD_ctrlroom_distil - CD_ctrlroom_flash
 - CD_ctrlroom_furnace - CD_ctrlroom_reactor - CD_ctrlroom_store
 - CD_distil_comp - CD_distil_ctrlroom - CD_distil_distil - CD_distil_flash
 - CD_distil_furnace - CD_distil_reactor - CD_distil_store - CD_flash_comp
 - CD_flash_ctrlroom - CD_flash_distil - CD_flash_flash - CD_flash_furnace
 - CD_flash_reactor - CD_flash_store - CD_furnace_comp - CD_furnace_ctrlroom
 - CD_furnace_distil - CD_furnace_flash - CD_furnace_furnace
 - CD_furnace_reactor - CD_furnace_store - CD_reactor_comp
 - CD_reactor_ctrlroom - CD_reactor_distil - CD_reactor_flash
 - CD_reactor_furnace - CD_reactor_reactor - CD_reactor_store - CD_store_comp
 - CD_store_ctrlroom - CD_store_distil - CD_store_flash - CD_store_furnace
 - CD_store_reactor - CD_store_store + SumCD - roadCost = 0
_C439: X_furnace + 0.5 l_furnace <= 50
_C44: A_furnace_reactor - 1000 Wy_furnace_reactor <= 0
_C440: Y_furnace + 0.5 d_furnace <= 50
_C441: X_reactor + 0.5 l_reactor <= 50
_C442: Y_reactor + 0.5 d_reactor <= 50
_C443: X_flash + 0.5 l_flash <= 50
_C444: Y_flash + 0.5 d_flash <= 50
_C445: X_comp + 0.5 l_comp <= 50
_C446: Y_comp + 0.5 d_comp <= 50
_C447: X_distil + 0.5 l_distil <= 50
_C448: Y_distil + 0.5 d_distil <= 50
_C449: X_store + 0.5 l_store <= 50
_C45: B_furnace_reactor + 1000 Wy_furnace_reactor <= 1000
_C450: Y_store + 0.5 d_store <= 50
_C451: X_ctrlroom + 0.5 l_ctrlroom <= 50
_C452: Y_ctrlroom + 0.5 d_ctrlroom <= 50
_C453: Ve2_furnace_furnace = 0
_C454: - D_furnace_reactor + Din_furnace_reactor + Dout_furnace_reactor = 0
_C455: Din_furnace_reactor - 40.7293505963 Psie_furnace_reactor <= 0
_C456: Dout_furnace_reactor + 40.7293505963 Psie_furnace_reactor
 >= 40.7293505963
_C457: Dout_furnace_reactor + 1000 Psie_furnace_reactor <= 1000
_C458: 104190.711069 Din_furnace_reactor - 4243620 Psie_furnace_reactor
 + Ve2_furnace_reactor = 0
_C459: - D_furnace_flash + Din_furnace_flash + Dout_furnace_flash = 0
_C46: - A_furnace_reactor - B_furnace_reactor + D_furnace_reactor
 - L_furnace_reactor - R_furnace_reactor = 0
_C460: Din_furnace_flash - 40.7293505963 Psie_furnace_flash <= 0
_C461: Dout_furnace_flash + 40.7293505963 Psie_furnace_flash >= 40.7293505963
_C462: Dout_furnace_flash + 1000 Psie_furnace_flash <= 1000
_C463: 67858.6807679 Din_furnace_flash - 2763840 Psie_furnace_flash
 + Ve2_furnace_flash = 0
_C464: - D_furnace_comp + Din_furnace_comp + Dout_furnace_comp = 0
_C465: Din_furnace_comp - 40.7293505963 Psie_furnace_comp <= 0
_C466: Dout_furnace_comp + 40.7293505963 Psie_furnace_comp >= 40.7293505963
_C467: Dout_furnace_comp + 1000 Psie_furnace_comp <= 1000
_C468: 157788.914036 Din_furnace_comp - 6426640 Psie_furnace_comp
 + Ve2_furnace_comp = 0
_C469: - D_furnace_distil + Din_furnace_distil + Dout_furnace_distil = 0
_C47: D_furnace_reactor - D_reactor_furnace = 0
_C470: Din_furnace_distil - 40.7293505963 Psie_furnace_distil <= 0
_C471: Dout_furnace_distil + 40.7293505963 Psie_furnace_distil
 >= 40.7293505963
_C472: Dout_furnace_distil + 1000 Psie_furnace_distil <= 1000
_C473: 203996.378001 Din_furnace_distil - 8308640 Psie_furnace_distil
 + Ve2_furnace_distil = 0
_C474: - D_furnace_store + Din_furnace_store + Dout_furnace_store = 0
_C475: Din_furnace_store - 40.7293505963 Psie_furnace_store <= 0
_C476: Dout_furnace_store + 40.7293505963 Psie_furnace_store >= 40.7293505963
_C477: Dout_furnace_store + 1000 Psie_furnace_store <= 1000
_C478: 175592.487857 Din_furnace_store - 7151768 Psie_furnace_store
 + Ve2_furnace_store = 0
_C479: - D_furnace_ctrlroom + Din_furnace_ctrlroom + Dout_furnace_ctrlroom = 0
_C48: 1000 E1_furnace_reactor + 1000 E2_furnace_reactor + X_furnace
 - X_reactor - 0.5 l_furnace - 0.5 l_reactor >= 0
_C480: Din_furnace_ctrlroom - 40.7293505963 Psie_furnace_ctrlroom <= 0
_C481: Dout_furnace_ctrlroom + 40.7293505963 Psie_furnace_ctrlroom
 >= 40.7293505963
_C482: Dout_furnace_ctrlroom + 1000 Psie_furnace_ctrlroom <= 1000
_C483: 950837.89535 Din_furnace_ctrlroom - 38727010 Psie_furnace_ctrlroom
 + Ve2_furnace_ctrlroom = 0
_C484: - D_reactor_furnace + Din_reactor_furnace + Dout_reactor_furnace = 0
_C485: Din_reactor_furnace - 57.1907964624 Psie_reactor_furnace <= 0
_C486: Dout_reactor_furnace + 57.1907964624 Psie_reactor_furnace
 >= 57.1907964624
_C487: Dout_reactor_furnace + 1000 Psie_reactor_furnace <= 1000
_C488: 48588.5871834 Din_reactor_furnace - 2778820 Psie_reactor_furnace
 + Ve2_reactor_furnace = 0
_C489: Ve2_reactor_reactor = 0
_C49: - 1000 E1_furnace_reactor + 1000 E2_furnace_reactor - X_furnace
 + X_reactor - 0.5 l_furnace - 0.5 l_reactor >= -1000
_C490: - D_reactor_flash + Din_reactor_flash + Dout_reactor_flash = 0
_C491: Din_reactor_flash - 57.1907964624 Psie_reactor_flash <= 0
_C492: Dout_reactor_flash + 57.1907964624 Psie_reactor_flash >= 57.1907964624
_C493: Dout_reactor_flash + 1000 Psie_reactor_flash <= 1000
_C494: 48326.6569267 Din_reactor_flash - 2763840 Psie_reactor_flash
 + Ve2_reactor_flash = 0
_C495: - D_reactor_comp + Din_reactor_comp + Dout_reactor_comp = 0
_C496: Din_reactor_comp - 57.1907964624 Psie_reactor_comp <= 0
_C497: Dout_reactor_comp + 57.1907964624 Psie_reactor_comp >= 57.1907964624
_C498: Dout_reactor_comp + 1000 Psie_reactor_comp <= 1000
_C499: 112371.92691 Din_reactor_comp - 6426640 Psie_reactor_comp
 + Ve2_reactor_comp = 0
_C5: - 2.8 O_reactor + l_reactor = 5.8
_C50: 1000 E1_furnace_reactor - 1000 E2_furnace_reactor + Y_furnace
 - Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -1000
_C500: - D_reactor_distil + Din_reactor_distil + Dout_reactor_distil = 0
_C501: Din_reactor_distil - 57.1907964624 Psie_reactor_distil <= 0
_C502: Dout_reactor_distil + 57.1907964624 Psie_reactor_distil
 >= 57.1907964624
_C503: Dout_reactor_distil + 1000 Psie_reactor_distil <= 1000
_C504: 145279.319645 Din_reactor_distil - 8308640 Psie_reactor_distil
 + Ve2_reactor_distil = 0
_C505: - D_reactor_store + Din_reactor_store + Dout_reactor_store = 0
_C506: Din_reactor_store - 57.1907964624 Psie_reactor_store <= 0
_C507: Dout_reactor_store + 57.1907964624 Psie_reactor_store >= 57.1907964624
_C508: Dout_reactor_store + 1000 Psie_reactor_store <= 1000
_C509: 125051.029928 Din_reactor_store - 7151768 Psie_reactor_store
 + Ve2_reactor_store = 0
_C51: - 1000 E1_furnace_reactor - 1000 E2_furnace_reactor - Y_furnace
 + Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -2000
_C510: - D_reactor_ctrlroom + Din_reactor_ctrlroom + Dout_reactor_ctrlroom = 0
_C511: Din_reactor_ctrlroom - 57.1907964624 Psie_reactor_ctrlroom <= 0
_C512: Dout_reactor_ctrlroom + 57.1907964624 Psie_reactor_ctrlroom
 >= 57.1907964624
_C513: Dout_reactor_ctrlroom + 1000 Psie_reactor_ctrlroom <= 1000
_C514: 677154.584226 Din_reactor_ctrlroom - 38727010 Psie_reactor_ctrlroom
 + Ve2_reactor_ctrlroom = 0
_C515: - D_distil_furnace + Din_distil_furnace + Dout_distil_furnace = 0
_C516: Din_distil_furnace - 50.063160108 Psie_distil_furnace <= 0
_C517: Dout_distil_furnace + 50.063160108 Psie_distil_furnace >= 50.063160108
_C518: Dout_distil_furnace + 1000 Psie_distil_furnace <= 1000
_C519: 55506.2843417 Din_distil_furnace - 2778820 Psie_distil_furnace
 + Ve2_distil_furnace = 0
_C52: L_furnace_reactor - R_reactor_furnace = 0
_C520: - D_distil_reactor + Din_distil_reactor + Dout_distil_reactor = 0
_C521: Din_distil_reactor - 50.063160108 Psie_distil_reactor <= 0
_C522: Dout_distil_reactor + 50.063160108 Psie_distil_reactor >= 50.063160108
_C523: Dout_distil_reactor + 1000 Psie_distil_reactor <= 1000
_C524: 84765.3242593 Din_distil_reactor - 4243620 Psie_distil_reactor
 + Ve2_distil_reactor = 0
_C525: - D_distil_flash + Din_distil_flash + Dout_distil_flash = 0
_C526: Din_distil_flash - 50.063160108 Psie_distil_flash <= 0
_C527: Dout_distil_flash + 50.063160108 Psie_distil_flash >= 50.063160108
_C528: Dout_distil_flash + 1000 Psie_distil_flash <= 1000
_C529: 55207.0623196 Din_distil_flash - 2763840 Psie_distil_flash
 + Ve2_distil_flash = 0
_C53: - L_reactor_furnace + R_furnace_reactor = 0
_C530: - D_distil_comp + Din_distil_comp + Dout_distil_comp = 0
_C531: Din_distil_comp - 50.063160108 Psie_distil_comp <= 0
_C532: Dout_distil_comp + 50.063160108 Psie_distil_comp >= 50.063160108
_C533: Dout_distil_comp + 1000 Psie_distil_comp <= 1000
_C534: 128370.641928 Din_distil_comp - 6426640 Psie_distil_comp
 + Ve2_distil_comp = 0
_C535: Ve2_distil_distil = 0
_C536: - D_distil_store + Din_distil_store + Dout_distil_store = 0
_C537: Din_distil_store - 50.063160108 Psie_distil_store <= 0
_C538: Dout_distil_store + 50.063160108 Psie_distil_store >= 50.063160108
_C539: Dout_distil_store + 1000 Psie_distil_store <= 1000
_C54: A_furnace_reactor - B_reactor_furnace = 0
_C540: 142854.905375 Din_distil_store - 7151768 Psie_distil_store
 + Ve2_distil_store = 0
_C541: - D_distil_ctrlroom + Din_distil_ctrlroom + Dout_distil_ctrlroom = 0
_C542: Din_distil_ctrlroom - 50.063160108 Psie_distil_ctrlroom <= 0
_C543: Dout_distil_ctrlroom + 50.063160108 Psie_distil_ctrlroom
 >= 50.063160108
_C544: Dout_distil_ctrlroom + 1000 Psie_distil_ctrlroom <= 1000
_C545: 773563.033505 Din_distil_ctrlroom - 38727010 Psie_distil_ctrlroom
 + Ve2_distil_ctrlroom = 0
_C546: - D_store_furnace + Din_store_furnace + Dout_store_furnace = 0
_C547: Din_store_furnace - 63.6396103068 Psie_store_furnace <= 0
_C548: Dout_store_furnace + 63.6396103068 Psie_store_furnace >= 63.6396103068
_C549: Dout_store_furnace + 1000 Psie_store_furnace <= 1000
_C55: - A_reactor_furnace + B_furnace_reactor = 0
_C550: 43664.9436822 Din_store_furnace - 2778820 Psie_store_furnace
 + Ve2_store_furnace = 0
_C551: - D_store_reactor + Din_store_reactor + Dout_store_reactor = 0
_C552: Din_store_reactor - 63.6396103068 Psie_store_reactor <= 0
_C553: Dout_store_reactor + 63.6396103068 Psie_store_reactor >= 63.6396103068
_C554: Dout_store_reactor + 1000 Psie_store_reactor <= 1000
_C555: 66682.055084 Din_store_reactor - 4243620 Psie_store_reactor
 + Ve2_store_reactor = 0
_C556: - D_store_flash + Din_store_flash + Dout_store_flash = 0
_C557: Din_store_flash - 63.6396103068 Psie_store_flash <= 0
_C558: Dout_store_flash + 63.6396103068 Psie_store_flash >= 63.6396103068
_C559: Dout_store_flash + 1000 Psie_store_flash <= 1000
_C56: Wx_furnace_reactor + Wx_reactor_furnace = 1
_C560: 43429.5556914 Din_store_flash - 2763840 Psie_store_flash
 + Ve2_store_flash = 0
_C561: - D_store_comp + Din_store_comp + Dout_store_comp = 0
_C562: Din_store_comp - 63.6396103068 Psie_store_comp <= 0
_C563: Dout_store_comp + 63.6396103068 Psie_store_comp >= 63.6396103068
_C564: Dout_store_comp + 1000 Psie_store_comp <= 1000
_C565: 100984.904983 Din_store_comp - 6426640 Psie_store_comp + Ve2_store_comp
 = 0
_C566: - D_store_distil + Din_store_distil + Dout_store_distil = 0
_C567: Din_store_distil - 63.6396103068 Psie_store_distil <= 0
_C568: Dout_store_distil + 63.6396103068 Psie_store_distil >= 63.6396103068
_C569: Dout_store_distil + 1000 Psie_store_distil <= 1000
_C57: Wy_furnace_reactor + Wy_reactor_furnace = 1
_C570: 130557.681921 Din_store_distil - 8308640 Psie_store_distil
 + Ve2_store_distil = 0
_C571: Ve2_store_store = 0
_C572: - D_store_ctrlroom + Din_store_ctrlroom + Dout_store_ctrlroom = 0
_C573: Din_store_ctrlroom - 63.6396103068 Psie_store_ctrlroom <= 0
_C574: Dout_store_ctrlroom + 63.6396103068 Psie_store_ctrlroom
 >= 63.6396103068
_C575: Dout_store_ctrlroom + 1000 Psie_store_ctrlroom <= 1000
_C576: 608536.253024 Din_store_ctrlroom - 38727010 Psie_store_ctrlroom
 + Ve2_store_ctrlroom = 0
_C577: - 0.001 Ve2_furnace_comp - 0.001 Ve2_furnace_ctrlroom
 - 0.001 Ve2_furnace_distil - 0.001 Ve2_furnace_flash
 - 0.001 Ve2_furnace_furnace - 0.001 Ve2_furnace_reactor
 - 0.001 Ve2_furnace_store + Ve_furnace = 0
_C578: Omega0_furnace - 0.75 Ve_furnace = 0
_C579: - 0.001 Ve2_reactor_comp - 0.001 Ve2_reactor_ctrlroom
 - 0.001 Ve2_reactor_distil - 0.001 Ve2_reactor_flash
 - 0.001 Ve2_reactor_furnace - 0.001 Ve2_reactor_reactor
 - 0.001 Ve2_reactor_store + Ve_reactor = 0
_C58: CD_furnace_flash = 0
_C580: Omega0_reactor - 0.82 Ve_reactor = 0
_C581: - 0.001 Ve2_distil_comp - 0.001 Ve2_distil_ctrlroom
 - 0.001 Ve2_distil_distil - 0.001 Ve2_distil_flash - 0.001 Ve2_distil_furnace
 - 0.001 Ve2_distil_reactor - 0.001 Ve2_distil_store + Ve_distil = 0
_C582: Omega0_distil - 0.77 Ve_distil = 0
_C583: - 0.001 Ve2_store_comp - 0.001 Ve2_store_ctrlroom
 - 0.001 Ve2_store_distil - 0.001 Ve2_store_flash - 0.001 Ve2_store_furnace
 - 0.001 Ve2_store_reactor - 0.001 Ve2_store_store + Ve_store = 0
_C584: Omega0_store - 0.82 Ve_store = 0
_C585: - Omega0_distil - Omega0_furnace - Omega0_reactor - Omega0_store
 + SumOmega = 0
_C59: - L_furnace_flash + R_furnace_flash + X_flash - X_furnace = 0
_C6: d_reactor + l_reactor = 14.4
_C60: A_furnace_flash - B_furnace_flash + Y_flash - Y_furnace = 0
_C61: R_furnace_flash - 1000 Wx_furnace_flash <= 0
_C62: L_furnace_flash + 1000 Wx_furnace_flash <= 1000
_C63: A_furnace_flash - 1000 Wy_furnace_flash <= 0
_C64: B_furnace_flash + 1000 Wy_furnace_flash <= 1000
_C65: - A_furnace_flash - B_furnace_flash + D_furnace_flash - L_furnace_flash
 - R_furnace_flash = 0
_C66: - D_flash_furnace + D_furnace_flash = 0
_C67: 1000 E1_furnace_flash + 1000 E2_furnace_flash - X_flash + X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= 0
_C68: - 1000 E1_furnace_flash + 1000 E2_furnace_flash + X_flash - X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= -1000
_C69: 1000 E1_furnace_flash - 1000 E2_furnace_flash - Y_flash + Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -1000
_C7: X_reactor - 0.5 l_reactor >= 0
_C70: - 1000 E1_furnace_flash - 1000 E2_furnace_flash + Y_flash - Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -2000
_C71: L_furnace_flash - R_flash_furnace = 0
_C72: - L_flash_furnace + R_furnace_flash = 0
_C73: A_furnace_flash - B_flash_furnace = 0
_C74: - A_flash_furnace + B_furnace_flash = 0
_C75: Wx_flash_furnace + Wx_furnace_flash = 1
_C76: Wy_flash_furnace + Wy_furnace_flash = 1
_C77: CD_reactor_flash - 1231.60506417 D_reactor_flash = 0
_C78: - L_reactor_flash + R_reactor_flash + X_flash - X_reactor = 0
_C79: A_reactor_flash - B_reactor_flash + Y_flash - Y_reactor = 0
_C8: Y_reactor - 0.5 d_reactor >= 0
_C80: R_reactor_flash - 1000 Wx_reactor_flash <= 0
_C81: L_reactor_flash + 1000 Wx_reactor_flash <= 1000
_C82: A_reactor_flash - 1000 Wy_reactor_flash <= 0
_C83: B_reactor_flash + 1000 Wy_reactor_flash <= 1000
_C84: - A_reactor_flash - B_reactor_flash + D_reactor_flash - L_reactor_flash
 - R_reactor_flash = 0
_C85: - D_flash_reactor + D_reactor_flash = 0
_C86: 1000 E1_reactor_flash + 1000 E2_reactor_flash - X_flash + X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= 0
_C87: - 1000 E1_reactor_flash + 1000 E2_reactor_flash + X_flash - X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= -1000
_C88: 1000 E1_reactor_flash - 1000 E2_reactor_flash - Y_flash + Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -1000
_C89: - 1000 E1_reactor_flash - 1000 E2_reactor_flash + Y_flash - Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -2000
_C9: 0 O_flash + l_flash = 4
_C90: L_reactor_flash - R_flash_reactor = 0
_C91: - L_flash_reactor + R_reactor_flash = 0
_C92: A_reactor_flash - B_flash_reactor = 0
_C93: - A_flash_reactor + B_reactor_flash = 0
_C94: Wx_flash_reactor + Wx_reactor_flash = 1
_C95: Wy_flash_reactor + Wy_reactor_flash = 1
_C96: CD_furnace_comp - 1487.20883073 D_furnace_comp = 0
_C97: - L_furnace_comp + R_furnace_comp + X_comp - X_furnace = 0
_C98: A_furnace_comp - B_furnace_comp + Y_comp - Y_furnace = 0
_C99: R_furnace_comp - 1000 Wx_furnace_comp <= 0
Bounds
 X_comp free
 X_ctrlroom free
 X_distil free
 X_flash free
 X_furnace free
 X_reactor free
 X_store = 11
 Y_comp free
 Y_ctrlroom free
 Y_distil free
 Y_flash free
 Y_furnace free
 Y_reactor free
 Y_store = 6
Binaries
E1_comp_ctrlroom
E1_comp_distil
E1_comp_store
E1_distil_ctrlroom
E1_distil_store
E1_flash_comp
E1_flash_ctrlroom
E1_flash_distil
E1_flash_store
E1_furnace_comp
E1_furnace_ctrlroom
E1_furnace_distil
E1_furnace_flash
E1_furnace_reactor
E1_furnace_store
E1_reactor_comp
E1_reactor_ctrlroom
E1_reactor_distil
E1_reactor_flash
E1_reactor_store
E1_store_ctrlroom
E2_comp_ctrlroom
E2_comp_distil
E2_comp_store
E2_distil_ctrlroom
E2_distil_store
E2_flash_comp
E2_flash_ctrlroom
E2_flash_distil
E2_flash_store
E2_furnace_comp
E2_furnace_ctrlroom
E2_furnace_distil
E2_furnace_flash
E2_furnace_reactor
E2_furnace_store
E2_reactor_comp
E2_reactor_ctrlroom
E2_reactor_distil
E2_reactor_flash
E2_reactor_store
E2_store_ctrlroom
O_comp
O_ctrlroom
O_distil
O_flash
O_furnace
O_reactor
O_store
Psie_distil_comp
Psie_distil_ctrlroom
Psie_distil_flash
Psie_distil_furnace
Psie_distil_reactor
Psie_distil_store
Psie_furnace_comp
Psie_furnace_ctrlroom
Psie_furnace_distil
Psie_furnace_flash
Psie_furnace_reactor
Psie_furnace_store
Psie_reactor_comp
Psie_reactor_ctrlroom
Psie_reactor_distil
Psie_reactor_flash
Psie_reactor_furnace
Psie_reactor_store
Psie_store_comp
Psie_store_ctrlroom
Psie_store_distil
Psie_store_flash
Psie_store_furnace
Psie_store_reactor
Wx2
Wx_comp_ctrlroom
Wx_comp_distil
Wx_comp_flash
Wx_comp_furnace
Wx_comp_reactor
Wx_comp_store
Wx_ctrlroom_comp
Wx_ctrlroom_distil
Wx_ctrlroom_flash
Wx_ctrlroom_furnace
Wx_ctrlroom_reactor
Wx_ctrlroom_store
Wx_distil_comp
Wx_distil_ctrlroom
Wx_distil_flash
Wx_distil_furnace
Wx_distil_reactor
Wx_distil_store
Wx_flash_comp
Wx_flash_ctrlroom
Wx_flash_distil
Wx_flash_furnace
Wx_flash_reactor
Wx_flash_store
Wx_furnace_comp
Wx_furnace_ctrlroom
Wx_furnace_distil
Wx_furnace_flash
Wx_furnace_reactor
Wx_furnace_store
Wx_reactor_comp
Wx_reactor_ctrlroom
Wx_reactor_distil
Wx_reactor_flash
Wx_reactor_furnace
Wx_reactor_store
Wx_store_comp
Wx_store_ctrlroom
Wx_store_distil
Wx_store_flash
Wx_store_furnace
Wx_store_reactor
Wy2
Wy_comp_ctrlroom
Wy_comp_distil
Wy_comp_flash
Wy_comp_furnace
Wy_comp_reactor
Wy_comp_store
Wy_ctrlroom_comp
Wy_ctrlroom_distil
Wy_ctrlroom_flash
Wy_ctrlroom_furnace
Wy_ctrlroom_reactor
Wy_ctrlroom_store
Wy_distil_comp
Wy_distil_ctrlroom
Wy_distil_flash
Wy_distil_furnace
Wy_distil_reactor
Wy_distil_store
Wy_flash_comp
Wy_flash_ctrlroom
Wy_flash_distil
Wy_flash_furnace
Wy_flash_reactor
Wy_flash_store
Wy_furnace_comp
Wy_furnace_ctrlroom
Wy_furnace_distil
Wy_furnace_flash
Wy_furnace_reactor
Wy_furnace_store
Wy_reactor_comp
Wy_reactor_ctrlroom
Wy_reactor_distil
Wy_reactor_flash
Wy_reactor_furnace
Wy_reactor_store
Wy_store_comp
Wy_store_ctrlroom
Wy_store_distil
Wy_store_flash
Wy_store_furnace
Wy_store_reactor
End
