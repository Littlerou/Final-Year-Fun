\* Layout_Problem_Model *\
Minimize
OBJ: SumCD + SumOmega
Subject To
_C1: 9 O_furnace + l_furnace = 16
_C10: d_flash + l_flash = 8
_C100: - L_comp_furnace + R_furnace_comp = 0
_C101: A_furnace_comp - B_comp_furnace = 0
_C102: - A_comp_furnace + B_furnace_comp = 0
_C103: Wx_comp_furnace + Wx_furnace_comp = 1
_C104: Wy_comp_furnace + Wy_furnace_comp = 1
_C105: CD_reactor_comp = 0
_C106: - L_reactor_comp + R_reactor_comp + X_comp - X_reactor = 0
_C107: A_reactor_comp - B_reactor_comp + Y_comp - Y_reactor = 0
_C108: R_reactor_comp - 1000 Wx_reactor_comp <= 0
_C109: L_reactor_comp + 1000 Wx_reactor_comp <= 1000
_C11: X_flash - 0.5 l_flash >= 0
_C110: A_reactor_comp - 1000 Wy_reactor_comp <= 0
_C111: B_reactor_comp + 1000 Wy_reactor_comp <= 1000
_C112: - A_reactor_comp - B_reactor_comp + D_reactor_comp - L_reactor_comp
 - R_reactor_comp = 0
_C113: - D_comp_reactor + D_reactor_comp = 0
_C114: 1000 E1_reactor_comp + 1000 E2_reactor_comp - X_comp + X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= 0
_C115: - 1000 E1_reactor_comp + 1000 E2_reactor_comp + X_comp - X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= -1000
_C116: 1000 E1_reactor_comp - 1000 E2_reactor_comp - Y_comp + Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -1000
_C117: - 1000 E1_reactor_comp - 1000 E2_reactor_comp + Y_comp - Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -2000
_C118: L_reactor_comp - R_comp_reactor = 0
_C119: - L_comp_reactor + R_reactor_comp = 0
_C12: Y_flash - 0.5 d_flash >= 0
_C120: A_reactor_comp - B_comp_reactor = 0
_C121: - A_comp_reactor + B_reactor_comp = 0
_C122: Wx_comp_reactor + Wx_reactor_comp = 1
_C123: Wy_comp_reactor + Wy_reactor_comp = 1
_C124: CD_flash_comp - 842.70764367 D_flash_comp = 0
_C125: - L_flash_comp + R_flash_comp + X_comp - X_flash = 0
_C126: A_flash_comp - B_flash_comp + Y_comp - Y_flash = 0
_C127: R_flash_comp - 1000 Wx_flash_comp <= 0
_C128: L_flash_comp + 1000 Wx_flash_comp <= 1000
_C129: A_flash_comp - 1000 Wy_flash_comp <= 0
_C13: 0 O_comp + l_comp = 14
_C130: B_flash_comp + 1000 Wy_flash_comp <= 1000
_C131: - A_flash_comp - B_flash_comp + D_flash_comp - L_flash_comp
 - R_flash_comp = 0
_C132: - D_comp_flash + D_flash_comp = 0
_C133: 1000 E1_flash_comp + 1000 E2_flash_comp - X_comp + X_flash - 0.5 l_comp
 - 0.5 l_flash >= 0
_C134: - 1000 E1_flash_comp + 1000 E2_flash_comp + X_comp - X_flash
 - 0.5 l_comp - 0.5 l_flash >= -1000
_C135: 1000 E1_flash_comp - 1000 E2_flash_comp - Y_comp + Y_flash - 0.5 d_comp
 - 0.5 d_flash >= -1000
_C136: - 1000 E1_flash_comp - 1000 E2_flash_comp + Y_comp - Y_flash
 - 0.5 d_comp - 0.5 d_flash >= -2000
_C137: L_flash_comp - R_comp_flash = 0
_C138: - L_comp_flash + R_flash_comp = 0
_C139: A_flash_comp - B_comp_flash = 0
_C14: d_comp + l_comp = 28
_C140: - A_comp_flash + B_flash_comp = 0
_C141: Wx_comp_flash + Wx_flash_comp = 1
_C142: Wy_comp_flash + Wy_flash_comp = 1
_C143: CD_furnace_distil = 0
_C144: - L_furnace_distil + R_furnace_distil + X_distil - X_furnace = 0
_C145: A_furnace_distil - B_furnace_distil + Y_distil - Y_furnace = 0
_C146: R_furnace_distil - 1000 Wx_furnace_distil <= 0
_C147: L_furnace_distil + 1000 Wx_furnace_distil <= 1000
_C148: A_furnace_distil - 1000 Wy_furnace_distil <= 0
_C149: B_furnace_distil + 1000 Wy_furnace_distil <= 1000
_C15: X_comp - 0.5 l_comp >= 0
_C150: - A_furnace_distil - B_furnace_distil + D_furnace_distil
 - L_furnace_distil - R_furnace_distil = 0
_C151: - D_distil_furnace + D_furnace_distil = 0
_C152: 1000 E1_furnace_distil + 1000 E2_furnace_distil - X_distil + X_furnace
 - 0.5 l_distil - 0.5 l_furnace >= 0
_C153: - 1000 E1_furnace_distil + 1000 E2_furnace_distil + X_distil
 - X_furnace - 0.5 l_distil - 0.5 l_furnace >= -1000
_C154: 1000 E1_furnace_distil - 1000 E2_furnace_distil - Y_distil + Y_furnace
 - 0.5 d_distil - 0.5 d_furnace >= -1000
_C155: - 1000 E1_furnace_distil - 1000 E2_furnace_distil + Y_distil
 - Y_furnace - 0.5 d_distil - 0.5 d_furnace >= -2000
_C156: L_furnace_distil - R_distil_furnace = 0
_C157: - L_distil_furnace + R_furnace_distil = 0
_C158: A_furnace_distil - B_distil_furnace = 0
_C159: - A_distil_furnace + B_furnace_distil = 0
_C16: Y_comp - 0.5 d_comp >= 0
_C160: Wx_distil_furnace + Wx_furnace_distil = 1
_C161: Wy_distil_furnace + Wy_furnace_distil = 1
_C162: CD_reactor_distil = 0
_C163: - L_reactor_distil + R_reactor_distil + X_distil - X_reactor = 0
_C164: A_reactor_distil - B_reactor_distil + Y_distil - Y_reactor = 0
_C165: R_reactor_distil - 1000 Wx_reactor_distil <= 0
_C166: L_reactor_distil + 1000 Wx_reactor_distil <= 1000
_C167: A_reactor_distil - 1000 Wy_reactor_distil <= 0
_C168: B_reactor_distil + 1000 Wy_reactor_distil <= 1000
_C169: - A_reactor_distil - B_reactor_distil + D_reactor_distil
 - L_reactor_distil - R_reactor_distil = 0
_C17: 0 O_distil + l_distil = 8
_C170: - D_distil_reactor + D_reactor_distil = 0
_C171: 1000 E1_reactor_distil + 1000 E2_reactor_distil - X_distil + X_reactor
 - 0.5 l_distil - 0.5 l_reactor >= 0
_C172: - 1000 E1_reactor_distil + 1000 E2_reactor_distil + X_distil
 - X_reactor - 0.5 l_distil - 0.5 l_reactor >= -1000
_C173: 1000 E1_reactor_distil - 1000 E2_reactor_distil - Y_distil + Y_reactor
 - 0.5 d_distil - 0.5 d_reactor >= -1000
_C174: - 1000 E1_reactor_distil - 1000 E2_reactor_distil + Y_distil
 - Y_reactor - 0.5 d_distil - 0.5 d_reactor >= -2000
_C175: L_reactor_distil - R_distil_reactor = 0
_C176: - L_distil_reactor + R_reactor_distil = 0
_C177: A_reactor_distil - B_distil_reactor = 0
_C178: - A_distil_reactor + B_reactor_distil = 0
_C179: Wx_distil_reactor + Wx_reactor_distil = 1
_C18: d_distil + l_distil = 16
_C180: Wy_distil_reactor + Wy_reactor_distil = 1
_C181: CD_flash_distil - 78.4818701125 D_flash_distil = 0
_C182: - L_flash_distil + R_flash_distil + X_distil - X_flash = 0
_C183: A_flash_distil - B_flash_distil + Y_distil - Y_flash = 0
_C184: R_flash_distil - 1000 Wx_flash_distil <= 0
_C185: L_flash_distil + 1000 Wx_flash_distil <= 1000
_C186: A_flash_distil - 1000 Wy_flash_distil <= 0
_C187: B_flash_distil + 1000 Wy_flash_distil <= 1000
_C188: - A_flash_distil - B_flash_distil + D_flash_distil - L_flash_distil
 - R_flash_distil = 0
_C189: - D_distil_flash + D_flash_distil = 0
_C19: X_distil - 0.5 l_distil >= 0
_C190: 1000 E1_flash_distil + 1000 E2_flash_distil - X_distil + X_flash
 - 0.5 l_distil - 0.5 l_flash >= 0
_C191: - 1000 E1_flash_distil + 1000 E2_flash_distil + X_distil - X_flash
 - 0.5 l_distil - 0.5 l_flash >= -1000
_C192: 1000 E1_flash_distil - 1000 E2_flash_distil - Y_distil + Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -1000
_C193: - 1000 E1_flash_distil - 1000 E2_flash_distil + Y_distil - Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -2000
_C194: L_flash_distil - R_distil_flash = 0
_C195: - L_distil_flash + R_flash_distil = 0
_C196: A_flash_distil - B_distil_flash = 0
_C197: - A_distil_flash + B_flash_distil = 0
_C198: Wx_distil_flash + Wx_flash_distil = 1
_C199: Wy_distil_flash + Wy_flash_distil = 1
_C2: d_furnace + l_furnace = 23
_C20: Y_distil - 0.5 d_distil >= 0
_C200: CD_comp_distil = 0
_C201: - L_comp_distil + R_comp_distil - X_comp + X_distil = 0
_C202: A_comp_distil - B_comp_distil - Y_comp + Y_distil = 0
_C203: R_comp_distil - 1000 Wx_comp_distil <= 0
_C204: L_comp_distil + 1000 Wx_comp_distil <= 1000
_C205: A_comp_distil - 1000 Wy_comp_distil <= 0
_C206: B_comp_distil + 1000 Wy_comp_distil <= 1000
_C207: - A_comp_distil - B_comp_distil + D_comp_distil - L_comp_distil
 - R_comp_distil = 0
_C208: D_comp_distil - D_distil_comp = 0
_C209: 1000 E1_comp_distil + 1000 E2_comp_distil + X_comp - X_distil
 - 0.5 l_comp - 0.5 l_distil >= 0
_C21: - 10 O_store + l_store = 12
_C210: - 1000 E1_comp_distil + 1000 E2_comp_distil - X_comp + X_distil
 - 0.5 l_comp - 0.5 l_distil >= -1000
_C211: 1000 E1_comp_distil - 1000 E2_comp_distil + Y_comp - Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -1000
_C212: - 1000 E1_comp_distil - 1000 E2_comp_distil - Y_comp + Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -2000
_C213: L_comp_distil - R_distil_comp = 0
_C214: - L_distil_comp + R_comp_distil = 0
_C215: A_comp_distil - B_distil_comp = 0
_C216: - A_distil_comp + B_comp_distil = 0
_C217: Wx_comp_distil + Wx_distil_comp = 1
_C218: Wy_comp_distil + Wy_distil_comp = 1
_C219: CD_furnace_store - 75.8045442191 D_furnace_store = 0
_C22: d_store + l_store = 34
_C220: - L_furnace_store + R_furnace_store - X_furnace + X_store = 0
_C221: A_furnace_store - B_furnace_store - Y_furnace + Y_store = 0
_C222: R_furnace_store - 1000 Wx_furnace_store <= 0
_C223: L_furnace_store + 1000 Wx_furnace_store <= 1000
_C224: A_furnace_store - 1000 Wy_furnace_store <= 0
_C225: B_furnace_store + 1000 Wy_furnace_store <= 1000
_C226: - A_furnace_store - B_furnace_store + D_furnace_store - L_furnace_store
 - R_furnace_store = 0
_C227: D_furnace_store - D_store_furnace = 0
_C228: 1000 E1_furnace_store + 1000 E2_furnace_store + X_furnace - X_store
 - 0.5 l_furnace - 0.5 l_store >= 0
_C229: - 1000 E1_furnace_store + 1000 E2_furnace_store - X_furnace + X_store
 - 0.5 l_furnace - 0.5 l_store >= -1000
_C23: X_store - 0.5 l_store >= 0
_C230: 1000 E1_furnace_store - 1000 E2_furnace_store + Y_furnace - Y_store
 - 0.5 d_furnace - 0.5 d_store >= -1000
_C231: - 1000 E1_furnace_store - 1000 E2_furnace_store - Y_furnace + Y_store
 - 0.5 d_furnace - 0.5 d_store >= -2000
_C232: L_furnace_store - R_store_furnace = 0
_C233: - L_store_furnace + R_furnace_store = 0
_C234: A_furnace_store - B_store_furnace = 0
_C235: - A_store_furnace + B_furnace_store = 0
_C236: Wx_furnace_store + Wx_store_furnace = 1
_C237: Wy_furnace_store + Wy_store_furnace = 1
_C238: CD_reactor_store = 0
_C239: - L_reactor_store + R_reactor_store - X_reactor + X_store = 0
_C24: Y_store - 0.5 d_store >= 0
_C240: A_reactor_store - B_reactor_store - Y_reactor + Y_store = 0
_C241: R_reactor_store - 1000 Wx_reactor_store <= 0
_C242: L_reactor_store + 1000 Wx_reactor_store <= 1000
_C243: A_reactor_store - 1000 Wy_reactor_store <= 0
_C244: B_reactor_store + 1000 Wy_reactor_store <= 1000
_C245: - A_reactor_store - B_reactor_store + D_reactor_store - L_reactor_store
 - R_reactor_store = 0
_C246: D_reactor_store - D_store_reactor = 0
_C247: 1000 E1_reactor_store + 1000 E2_reactor_store + X_reactor - X_store
 - 0.5 l_reactor - 0.5 l_store >= 0
_C248: - 1000 E1_reactor_store + 1000 E2_reactor_store - X_reactor + X_store
 - 0.5 l_reactor - 0.5 l_store >= -1000
_C249: 1000 E1_reactor_store - 1000 E2_reactor_store + Y_reactor - Y_store
 - 0.5 d_reactor - 0.5 d_store >= -1000
_C25: 0 O_ctrlroom + l_ctrlroom = 15
_C250: - 1000 E1_reactor_store - 1000 E2_reactor_store - Y_reactor + Y_store
 - 0.5 d_reactor - 0.5 d_store >= -2000
_C251: L_reactor_store - R_store_reactor = 0
_C252: - L_store_reactor + R_reactor_store = 0
_C253: A_reactor_store - B_store_reactor = 0
_C254: - A_store_reactor + B_reactor_store = 0
_C255: Wx_reactor_store + Wx_store_reactor = 1
_C256: Wy_reactor_store + Wy_store_reactor = 1
_C257: CD_flash_store = 0
_C258: - L_flash_store + R_flash_store - X_flash + X_store = 0
_C259: A_flash_store - B_flash_store - Y_flash + Y_store = 0
_C26: d_ctrlroom + l_ctrlroom = 30
_C260: R_flash_store - 1000 Wx_flash_store <= 0
_C261: L_flash_store + 1000 Wx_flash_store <= 1000
_C262: A_flash_store - 1000 Wy_flash_store <= 0
_C263: B_flash_store + 1000 Wy_flash_store <= 1000
_C264: - A_flash_store - B_flash_store + D_flash_store - L_flash_store
 - R_flash_store = 0
_C265: D_flash_store - D_store_flash = 0
_C266: 1000 E1_flash_store + 1000 E2_flash_store + X_flash - X_store
 - 0.5 l_flash - 0.5 l_store >= 0
_C267: - 1000 E1_flash_store + 1000 E2_flash_store - X_flash + X_store
 - 0.5 l_flash - 0.5 l_store >= -1000
_C268: 1000 E1_flash_store - 1000 E2_flash_store + Y_flash - Y_store
 - 0.5 d_flash - 0.5 d_store >= -1000
_C269: - 1000 E1_flash_store - 1000 E2_flash_store - Y_flash + Y_store
 - 0.5 d_flash - 0.5 d_store >= -2000
_C27: X_ctrlroom - 0.5 l_ctrlroom >= 0
_C270: L_flash_store - R_store_flash = 0
_C271: - L_store_flash + R_flash_store = 0
_C272: A_flash_store - B_store_flash = 0
_C273: - A_store_flash + B_flash_store = 0
_C274: Wx_flash_store + Wx_store_flash = 1
_C275: Wy_flash_store + Wy_store_flash = 1
_C276: CD_comp_store = 0
_C277: - L_comp_store + R_comp_store - X_comp + X_store = 0
_C278: A_comp_store - B_comp_store - Y_comp + Y_store = 0
_C279: R_comp_store - 1000 Wx_comp_store <= 0
_C28: Y_ctrlroom - 0.5 d_ctrlroom >= 0
_C280: L_comp_store + 1000 Wx_comp_store <= 1000
_C281: A_comp_store - 1000 Wy_comp_store <= 0
_C282: B_comp_store + 1000 Wy_comp_store <= 1000
_C283: - A_comp_store - B_comp_store + D_comp_store - L_comp_store
 - R_comp_store = 0
_C284: D_comp_store - D_store_comp = 0
_C285: 1000 E1_comp_store + 1000 E2_comp_store + X_comp - X_store - 0.5 l_comp
 - 0.5 l_store >= 0
_C286: - 1000 E1_comp_store + 1000 E2_comp_store - X_comp + X_store
 - 0.5 l_comp - 0.5 l_store >= -1000
_C287: 1000 E1_comp_store - 1000 E2_comp_store + Y_comp - Y_store - 0.5 d_comp
 - 0.5 d_store >= -1000
_C288: - 1000 E1_comp_store - 1000 E2_comp_store - Y_comp + Y_store
 - 0.5 d_comp - 0.5 d_store >= -2000
_C289: L_comp_store - R_store_comp = 0
_C29: CD_furnace_reactor - 1075.48971174 D_furnace_reactor = 0
_C290: - L_store_comp + R_comp_store = 0
_C291: A_comp_store - B_store_comp = 0
_C292: - A_store_comp + B_comp_store = 0
_C293: Wx_comp_store + Wx_store_comp = 1
_C294: Wy_comp_store + Wy_store_comp = 1
_C295: CD_distil_store - 451.755491306 D_distil_store = 0
_C296: - L_distil_store + R_distil_store - X_distil + X_store = 0
_C297: A_distil_store - B_distil_store - Y_distil + Y_store = 0
_C298: R_distil_store - 1000 Wx_distil_store <= 0
_C299: L_distil_store + 1000 Wx_distil_store <= 1000
_C3: X_furnace - 0.5 l_furnace >= 0
_C30: - L_furnace_reactor + R_furnace_reactor - X_furnace + X_reactor = 0
_C300: A_distil_store - 1000 Wy_distil_store <= 0
_C301: B_distil_store + 1000 Wy_distil_store <= 1000
_C302: - A_distil_store - B_distil_store + D_distil_store - L_distil_store
 - R_distil_store = 0
_C303: D_distil_store - D_store_distil = 0
_C304: 1000 E1_distil_store + 1000 E2_distil_store + X_distil - X_store
 - 0.5 l_distil - 0.5 l_store >= 0
_C305: - 1000 E1_distil_store + 1000 E2_distil_store - X_distil + X_store
 - 0.5 l_distil - 0.5 l_store >= -1000
_C306: 1000 E1_distil_store - 1000 E2_distil_store + Y_distil - Y_store
 - 0.5 d_distil - 0.5 d_store >= -1000
_C307: - 1000 E1_distil_store - 1000 E2_distil_store - Y_distil + Y_store
 - 0.5 d_distil - 0.5 d_store >= -2000
_C308: L_distil_store - R_store_distil = 0
_C309: - L_store_distil + R_distil_store = 0
_C31: A_furnace_reactor - B_furnace_reactor - Y_furnace + Y_reactor = 0
_C310: A_distil_store - B_store_distil = 0
_C311: - A_store_distil + B_distil_store = 0
_C312: Wx_distil_store + Wx_store_distil = 1
_C313: Wy_distil_store + Wy_store_distil = 1
_C314: CD_furnace_ctrlroom = 0
_C315: - L_furnace_ctrlroom + R_furnace_ctrlroom + X_ctrlroom - X_furnace = 0
_C316: A_furnace_ctrlroom - B_furnace_ctrlroom + Y_ctrlroom - Y_furnace = 0
_C317: R_furnace_ctrlroom - 1000 Wx_furnace_ctrlroom <= 0
_C318: L_furnace_ctrlroom + 1000 Wx_furnace_ctrlroom <= 1000
_C319: A_furnace_ctrlroom - 1000 Wy_furnace_ctrlroom <= 0
_C32: R_furnace_reactor - 1000 Wx_furnace_reactor <= 0
_C320: B_furnace_ctrlroom + 1000 Wy_furnace_ctrlroom <= 1000
_C321: - A_furnace_ctrlroom - B_furnace_ctrlroom + D_furnace_ctrlroom
 - L_furnace_ctrlroom - R_furnace_ctrlroom = 0
_C322: - D_ctrlroom_furnace + D_furnace_ctrlroom = 0
_C323: 1000 E1_furnace_ctrlroom + 1000 E2_furnace_ctrlroom - X_ctrlroom
 + X_furnace - 0.5 l_ctrlroom - 0.5 l_furnace >= 0
_C324: - 1000 E1_furnace_ctrlroom + 1000 E2_furnace_ctrlroom + X_ctrlroom
 - X_furnace - 0.5 l_ctrlroom - 0.5 l_furnace >= -1000
_C325: 1000 E1_furnace_ctrlroom - 1000 E2_furnace_ctrlroom - Y_ctrlroom
 + Y_furnace - 0.5 d_ctrlroom - 0.5 d_furnace >= -1000
_C326: - 1000 E1_furnace_ctrlroom - 1000 E2_furnace_ctrlroom + Y_ctrlroom
 - Y_furnace - 0.5 d_ctrlroom - 0.5 d_furnace >= -2000
_C327: L_furnace_ctrlroom - R_ctrlroom_furnace = 0
_C328: - L_ctrlroom_furnace + R_furnace_ctrlroom = 0
_C329: A_furnace_ctrlroom - B_ctrlroom_furnace = 0
_C33: L_furnace_reactor + 1000 Wx_furnace_reactor <= 1000
_C330: - A_ctrlroom_furnace + B_furnace_ctrlroom = 0
_C331: Wx_ctrlroom_furnace + Wx_furnace_ctrlroom = 1
_C332: Wy_ctrlroom_furnace + Wy_furnace_ctrlroom = 1
_C333: CD_reactor_ctrlroom = 0
_C334: - L_reactor_ctrlroom + R_reactor_ctrlroom + X_ctrlroom - X_reactor = 0
_C335: A_reactor_ctrlroom - B_reactor_ctrlroom + Y_ctrlroom - Y_reactor = 0
_C336: R_reactor_ctrlroom - 1000 Wx_reactor_ctrlroom <= 0
_C337: L_reactor_ctrlroom + 1000 Wx_reactor_ctrlroom <= 1000
_C338: A_reactor_ctrlroom - 1000 Wy_reactor_ctrlroom <= 0
_C339: B_reactor_ctrlroom + 1000 Wy_reactor_ctrlroom <= 1000
_C34: A_furnace_reactor - 1000 Wy_furnace_reactor <= 0
_C340: - A_reactor_ctrlroom - B_reactor_ctrlroom + D_reactor_ctrlroom
 - L_reactor_ctrlroom - R_reactor_ctrlroom = 0
_C341: - D_ctrlroom_reactor + D_reactor_ctrlroom = 0
_C342: 1000 E1_reactor_ctrlroom + 1000 E2_reactor_ctrlroom - X_ctrlroom
 + X_reactor - 0.5 l_ctrlroom - 0.5 l_reactor >= 0
_C343: - 1000 E1_reactor_ctrlroom + 1000 E2_reactor_ctrlroom + X_ctrlroom
 - X_reactor - 0.5 l_ctrlroom - 0.5 l_reactor >= -1000
_C344: 1000 E1_reactor_ctrlroom - 1000 E2_reactor_ctrlroom - Y_ctrlroom
 + Y_reactor - 0.5 d_ctrlroom - 0.5 d_reactor >= -1000
_C345: - 1000 E1_reactor_ctrlroom - 1000 E2_reactor_ctrlroom + Y_ctrlroom
 - Y_reactor - 0.5 d_ctrlroom - 0.5 d_reactor >= -2000
_C346: L_reactor_ctrlroom - R_ctrlroom_reactor = 0
_C347: - L_ctrlroom_reactor + R_reactor_ctrlroom = 0
_C348: A_reactor_ctrlroom - B_ctrlroom_reactor = 0
_C349: - A_ctrlroom_reactor + B_reactor_ctrlroom = 0
_C35: B_furnace_reactor + 1000 Wy_furnace_reactor <= 1000
_C350: Wx_ctrlroom_reactor + Wx_reactor_ctrlroom = 1
_C351: Wy_ctrlroom_reactor + Wy_reactor_ctrlroom = 1
_C352: CD_flash_ctrlroom = 0
_C353: - L_flash_ctrlroom + R_flash_ctrlroom + X_ctrlroom - X_flash = 0
_C354: A_flash_ctrlroom - B_flash_ctrlroom + Y_ctrlroom - Y_flash = 0
_C355: R_flash_ctrlroom - 1000 Wx_flash_ctrlroom <= 0
_C356: L_flash_ctrlroom + 1000 Wx_flash_ctrlroom <= 1000
_C357: A_flash_ctrlroom - 1000 Wy_flash_ctrlroom <= 0
_C358: B_flash_ctrlroom + 1000 Wy_flash_ctrlroom <= 1000
_C359: - A_flash_ctrlroom - B_flash_ctrlroom + D_flash_ctrlroom
 - L_flash_ctrlroom - R_flash_ctrlroom = 0
_C36: - A_furnace_reactor - B_furnace_reactor + D_furnace_reactor
 - L_furnace_reactor - R_furnace_reactor = 0
_C360: - D_ctrlroom_flash + D_flash_ctrlroom = 0
_C361: 1000 E1_flash_ctrlroom + 1000 E2_flash_ctrlroom - X_ctrlroom + X_flash
 - 0.5 l_ctrlroom - 0.5 l_flash >= 0
_C362: - 1000 E1_flash_ctrlroom + 1000 E2_flash_ctrlroom + X_ctrlroom
 - X_flash - 0.5 l_ctrlroom - 0.5 l_flash >= -1000
_C363: 1000 E1_flash_ctrlroom - 1000 E2_flash_ctrlroom - Y_ctrlroom + Y_flash
 - 0.5 d_ctrlroom - 0.5 d_flash >= -1000
_C364: - 1000 E1_flash_ctrlroom - 1000 E2_flash_ctrlroom + Y_ctrlroom
 - Y_flash - 0.5 d_ctrlroom - 0.5 d_flash >= -2000
_C365: L_flash_ctrlroom - R_ctrlroom_flash = 0
_C366: - L_ctrlroom_flash + R_flash_ctrlroom = 0
_C367: A_flash_ctrlroom - B_ctrlroom_flash = 0
_C368: - A_ctrlroom_flash + B_flash_ctrlroom = 0
_C369: Wx_ctrlroom_flash + Wx_flash_ctrlroom = 1
_C37: D_furnace_reactor - D_reactor_furnace = 0
_C370: Wy_ctrlroom_flash + Wy_flash_ctrlroom = 1
_C371: CD_comp_ctrlroom = 0
_C372: - L_comp_ctrlroom + R_comp_ctrlroom - X_comp + X_ctrlroom = 0
_C373: A_comp_ctrlroom - B_comp_ctrlroom - Y_comp + Y_ctrlroom = 0
_C374: R_comp_ctrlroom - 1000 Wx_comp_ctrlroom <= 0
_C375: L_comp_ctrlroom + 1000 Wx_comp_ctrlroom <= 1000
_C376: A_comp_ctrlroom - 1000 Wy_comp_ctrlroom <= 0
_C377: B_comp_ctrlroom + 1000 Wy_comp_ctrlroom <= 1000
_C378: - A_comp_ctrlroom - B_comp_ctrlroom + D_comp_ctrlroom - L_comp_ctrlroom
 - R_comp_ctrlroom = 0
_C379: D_comp_ctrlroom - D_ctrlroom_comp = 0
_C38: 1000 E1_furnace_reactor + 1000 E2_furnace_reactor + X_furnace
 - X_reactor - 0.5 l_furnace - 0.5 l_reactor >= 0
_C380: 1000 E1_comp_ctrlroom + 1000 E2_comp_ctrlroom + X_comp - X_ctrlroom
 - 0.5 l_comp - 0.5 l_ctrlroom >= 0
_C381: - 1000 E1_comp_ctrlroom + 1000 E2_comp_ctrlroom - X_comp + X_ctrlroom
 - 0.5 l_comp - 0.5 l_ctrlroom >= -1000
_C382: 1000 E1_comp_ctrlroom - 1000 E2_comp_ctrlroom + Y_comp - Y_ctrlroom
 - 0.5 d_comp - 0.5 d_ctrlroom >= -1000
_C383: - 1000 E1_comp_ctrlroom - 1000 E2_comp_ctrlroom - Y_comp + Y_ctrlroom
 - 0.5 d_comp - 0.5 d_ctrlroom >= -2000
_C384: L_comp_ctrlroom - R_ctrlroom_comp = 0
_C385: - L_ctrlroom_comp + R_comp_ctrlroom = 0
_C386: A_comp_ctrlroom - B_ctrlroom_comp = 0
_C387: - A_ctrlroom_comp + B_comp_ctrlroom = 0
_C388: Wx_comp_ctrlroom + Wx_ctrlroom_comp = 1
_C389: Wy_comp_ctrlroom + Wy_ctrlroom_comp = 1
_C39: - 1000 E1_furnace_reactor + 1000 E2_furnace_reactor - X_furnace
 + X_reactor - 0.5 l_furnace - 0.5 l_reactor >= -1000
_C390: CD_distil_ctrlroom = 0
_C391: - L_distil_ctrlroom + R_distil_ctrlroom + X_ctrlroom - X_distil = 0
_C392: A_distil_ctrlroom - B_distil_ctrlroom + Y_ctrlroom - Y_distil = 0
_C393: R_distil_ctrlroom - 1000 Wx_distil_ctrlroom <= 0
_C394: L_distil_ctrlroom + 1000 Wx_distil_ctrlroom <= 1000
_C395: A_distil_ctrlroom - 1000 Wy_distil_ctrlroom <= 0
_C396: B_distil_ctrlroom + 1000 Wy_distil_ctrlroom <= 1000
_C397: - A_distil_ctrlroom - B_distil_ctrlroom + D_distil_ctrlroom
 - L_distil_ctrlroom - R_distil_ctrlroom = 0
_C398: - D_ctrlroom_distil + D_distil_ctrlroom = 0
_C399: 1000 E1_distil_ctrlroom + 1000 E2_distil_ctrlroom - X_ctrlroom
 + X_distil - 0.5 l_ctrlroom - 0.5 l_distil >= 0
_C4: Y_furnace - 0.5 d_furnace >= 0
_C40: 1000 E1_furnace_reactor - 1000 E2_furnace_reactor + Y_furnace
 - Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -1000
_C400: - 1000 E1_distil_ctrlroom + 1000 E2_distil_ctrlroom + X_ctrlroom
 - X_distil - 0.5 l_ctrlroom - 0.5 l_distil >= -1000
_C401: 1000 E1_distil_ctrlroom - 1000 E2_distil_ctrlroom - Y_ctrlroom
 + Y_distil - 0.5 d_ctrlroom - 0.5 d_distil >= -1000
_C402: - 1000 E1_distil_ctrlroom - 1000 E2_distil_ctrlroom + Y_ctrlroom
 - Y_distil - 0.5 d_ctrlroom - 0.5 d_distil >= -2000
_C403: L_distil_ctrlroom - R_ctrlroom_distil = 0
_C404: - L_ctrlroom_distil + R_distil_ctrlroom = 0
_C405: A_distil_ctrlroom - B_ctrlroom_distil = 0
_C406: - A_ctrlroom_distil + B_distil_ctrlroom = 0
_C407: Wx_ctrlroom_distil + Wx_distil_ctrlroom = 1
_C408: Wy_ctrlroom_distil + Wy_distil_ctrlroom = 1
_C409: CD_store_ctrlroom = 0
_C41: - 1000 E1_furnace_reactor - 1000 E2_furnace_reactor - Y_furnace
 + Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -2000
_C410: - L_store_ctrlroom + R_store_ctrlroom + X_ctrlroom - X_store = 0
_C411: A_store_ctrlroom - B_store_ctrlroom + Y_ctrlroom - Y_store = 0
_C412: R_store_ctrlroom - 1000 Wx_store_ctrlroom <= 0
_C413: L_store_ctrlroom + 1000 Wx_store_ctrlroom <= 1000
_C414: A_store_ctrlroom - 1000 Wy_store_ctrlroom <= 0
_C415: B_store_ctrlroom + 1000 Wy_store_ctrlroom <= 1000
_C416: - A_store_ctrlroom - B_store_ctrlroom + D_store_ctrlroom
 - L_store_ctrlroom - R_store_ctrlroom = 0
_C417: - D_ctrlroom_store + D_store_ctrlroom = 0
_C418: 1000 E1_store_ctrlroom + 1000 E2_store_ctrlroom - X_ctrlroom + X_store
 - 0.5 l_ctrlroom - 0.5 l_store >= 0
_C419: - 1000 E1_store_ctrlroom + 1000 E2_store_ctrlroom + X_ctrlroom
 - X_store - 0.5 l_ctrlroom - 0.5 l_store >= -1000
_C42: L_furnace_reactor - R_reactor_furnace = 0
_C420: 1000 E1_store_ctrlroom - 1000 E2_store_ctrlroom - Y_ctrlroom + Y_store
 - 0.5 d_ctrlroom - 0.5 d_store >= -1000
_C421: - 1000 E1_store_ctrlroom - 1000 E2_store_ctrlroom + Y_ctrlroom
 - Y_store - 0.5 d_ctrlroom - 0.5 d_store >= -2000
_C422: L_store_ctrlroom - R_ctrlroom_store = 0
_C423: - L_ctrlroom_store + R_store_ctrlroom = 0
_C424: A_store_ctrlroom - B_ctrlroom_store = 0
_C425: - A_ctrlroom_store + B_store_ctrlroom = 0
_C426: Wx_ctrlroom_store + Wx_store_ctrlroom = 1
_C427: Wy_ctrlroom_store + Wy_store_ctrlroom = 1
_C428: X_comp - X_ctrlroom + X_distil + X_flash + X_furnace + X_reactor
 + X_store + Y_comp - Y_ctrlroom + Y_distil + Y_flash + Y_furnace + Y_reactor
 + Y_store + roadCost = 0
_C429: - CD_comp_comp - CD_comp_ctrlroom - CD_comp_distil - CD_comp_flash
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
 - CD_store_reactor - CD_store_store + SumCD = 0
_C43: - L_reactor_furnace + R_furnace_reactor = 0
_C430: X_furnace + 0.5 l_furnace <= 150
_C431: Y_furnace + 0.5 d_furnace <= 150
_C432: X_reactor + 0.5 l_reactor <= 150
_C433: Y_reactor + 0.5 d_reactor <= 150
_C434: X_flash + 0.5 l_flash <= 150
_C435: Y_flash + 0.5 d_flash <= 150
_C436: X_comp + 0.5 l_comp <= 150
_C437: Y_comp + 0.5 d_comp <= 150
_C438: X_distil + 0.5 l_distil <= 150
_C439: Y_distil + 0.5 d_distil <= 150
_C44: A_furnace_reactor - B_reactor_furnace = 0
_C440: X_store + 0.5 l_store <= 150
_C441: Y_store + 0.5 d_store <= 150
_C442: X_ctrlroom + 0.5 l_ctrlroom <= 150
_C443: Y_ctrlroom + 0.5 d_ctrlroom <= 150
_C444: Ve2_furnace_furnace = 0
_C445: - D_furnace_reactor + Din_furnace_reactor + Dout_furnace_reactor = 0
_C446: Din_furnace_reactor - 28.8 Psie_furnace_reactor <= 0
_C447: Dout_furnace_reactor + 28.8 Psie_furnace_reactor >= 28.8
_C448: Dout_furnace_reactor + 1000 Psie_furnace_reactor <= 1000
_C449: 135010.416667 Din_furnace_reactor - 3888300 Psie_furnace_reactor
 + Ve2_furnace_reactor = 0
_C45: - A_reactor_furnace + B_furnace_reactor = 0
_C450: - D_furnace_flash + Din_furnace_flash + Dout_furnace_flash = 0
_C451: Din_furnace_flash - 28.8 Psie_furnace_flash <= 0
_C452: Dout_furnace_flash + 28.8 Psie_furnace_flash >= 28.8
_C453: Dout_furnace_flash + 1000 Psie_furnace_flash <= 1000
_C454: 88388.8888889 Din_furnace_flash - 2545600 Psie_furnace_flash
 + Ve2_furnace_flash = 0
_C455: - D_furnace_comp + Din_furnace_comp + Dout_furnace_comp = 0
_C456: Din_furnace_comp - 28.8 Psie_furnace_comp <= 0
_C457: Dout_furnace_comp + 28.8 Psie_furnace_comp >= 28.8
_C458: Dout_furnace_comp + 1000 Psie_furnace_comp <= 1000
_C459: 189152.777778 Din_furnace_comp - 5447600 Psie_furnace_comp
 + Ve2_furnace_comp = 0
_C46: Wx_furnace_reactor + Wx_reactor_furnace = 1
_C460: - D_furnace_distil + Din_furnace_distil + Dout_furnace_distil = 0
_C461: Din_furnace_distil - 28.8 Psie_furnace_distil <= 0
_C462: Dout_furnace_distil + 28.8 Psie_furnace_distil >= 28.8
_C463: Dout_furnace_distil + 1000 Psie_furnace_distil <= 1000
_C464: 245750 Din_furnace_distil - 7077600 Psie_furnace_distil
 + Ve2_furnace_distil = 0
_C465: - D_furnace_store + Din_furnace_store + Dout_furnace_store = 0
_C466: Din_furnace_store - 28.8 Psie_furnace_store <= 0
_C467: Dout_furnace_store + 28.8 Psie_furnace_store >= 28.8
_C468: Dout_furnace_store + 1000 Psie_furnace_store <= 1000
_C469: 197216.458333 Din_furnace_store - 5679834 Psie_furnace_store
 + Ve2_furnace_store = 0
_C47: Wy_furnace_reactor + Wy_reactor_furnace = 1
_C470: - D_furnace_ctrlroom + Din_furnace_ctrlroom + Dout_furnace_ctrlroom = 0
_C471: Din_furnace_ctrlroom - 28.8 Psie_furnace_ctrlroom <= 0
_C472: Dout_furnace_ctrlroom + 28.8 Psie_furnace_ctrlroom >= 28.8
_C473: Dout_furnace_ctrlroom + 1000 Psie_furnace_ctrlroom <= 1000
_C474: 1333755.20833 Din_furnace_ctrlroom - 38412150 Psie_furnace_ctrlroom
 + Ve2_furnace_ctrlroom = 0
_C475: - D_reactor_furnace + Din_reactor_furnace + Dout_reactor_furnace = 0
_C476: Din_reactor_furnace - 40.44 Psie_reactor_furnace <= 0
_C477: Dout_reactor_furnace + 40.44 Psie_reactor_furnace >= 40.44
_C478: Dout_reactor_furnace + 1000 Psie_reactor_furnace <= 1000
_C479: 63212.1661721 Din_reactor_furnace - 2556300 Psie_reactor_furnace
 + Ve2_reactor_furnace = 0
_C48: CD_furnace_flash = 0
_C480: Ve2_reactor_reactor = 0
_C481: - D_reactor_flash + Din_reactor_flash + Dout_reactor_flash = 0
_C482: Din_reactor_flash - 40.44 Psie_reactor_flash <= 0
_C483: Dout_reactor_flash + 40.44 Psie_reactor_flash >= 40.44
_C484: Dout_reactor_flash + 1000 Psie_reactor_flash <= 1000
_C485: 62947.5766568 Din_reactor_flash - 2545600 Psie_reactor_flash
 + Ve2_reactor_flash = 0
_C486: - D_reactor_comp + Din_reactor_comp + Dout_reactor_comp = 0
_C487: Din_reactor_comp - 40.44 Psie_reactor_comp <= 0
_C488: Dout_reactor_comp + 40.44 Psie_reactor_comp >= 40.44
_C489: Dout_reactor_comp + 1000 Psie_reactor_comp <= 1000
_C49: - L_furnace_flash + R_furnace_flash + X_flash - X_furnace = 0
_C490: 134708.209693 Din_reactor_comp - 5447600 Psie_reactor_comp
 + Ve2_reactor_comp = 0
_C491: - D_reactor_distil + Din_reactor_distil + Dout_reactor_distil = 0
_C492: Din_reactor_distil - 40.44 Psie_reactor_distil <= 0
_C493: Dout_reactor_distil + 40.44 Psie_reactor_distil >= 40.44
_C494: Dout_reactor_distil + 1000 Psie_reactor_distil <= 1000
_C495: 175014.836795 Din_reactor_distil - 7077600 Psie_reactor_distil
 + Ve2_reactor_distil = 0
_C496: - D_reactor_store + Din_reactor_store + Dout_reactor_store = 0
_C497: Din_reactor_store - 40.44 Psie_reactor_store <= 0
_C498: Dout_reactor_store + 40.44 Psie_reactor_store >= 40.44
_C499: Dout_reactor_store + 1000 Psie_reactor_store <= 1000
_C5: - 4 O_reactor + l_reactor = 5
_C50: A_furnace_flash - B_furnace_flash + Y_flash - Y_furnace = 0
_C500: 140450.890208 Din_reactor_store - 5679834 Psie_reactor_store
 + Ve2_reactor_store = 0
_C501: - D_reactor_ctrlroom + Din_reactor_ctrlroom + Dout_reactor_ctrlroom = 0
_C502: Din_reactor_ctrlroom - 40.44 Psie_reactor_ctrlroom <= 0
_C503: Dout_reactor_ctrlroom + 40.44 Psie_reactor_ctrlroom >= 40.44
_C504: Dout_reactor_ctrlroom + 1000 Psie_reactor_ctrlroom <= 1000
_C505: 949855.341246 Din_reactor_ctrlroom - 38412150 Psie_reactor_ctrlroom
 + Ve2_reactor_ctrlroom = 0
_C506: - D_distil_furnace + Din_distil_furnace + Dout_distil_furnace = 0
_C507: Din_distil_furnace - 35.4 Psie_distil_furnace <= 0
_C508: Dout_distil_furnace + 35.4 Psie_distil_furnace >= 35.4
_C509: Dout_distil_furnace + 1000 Psie_distil_furnace <= 1000
_C51: R_furnace_flash - 1000 Wx_furnace_flash <= 0
_C510: 72211.8644068 Din_distil_furnace - 2556300 Psie_distil_furnace
 + Ve2_distil_furnace = 0
_C511: - D_distil_reactor + Din_distil_reactor + Dout_distil_reactor = 0
_C512: Din_distil_reactor - 35.4 Psie_distil_reactor <= 0
_C513: Dout_distil_reactor + 35.4 Psie_distil_reactor >= 35.4
_C514: Dout_distil_reactor + 1000 Psie_distil_reactor <= 1000
_C515: 109838.983051 Din_distil_reactor - 3888300 Psie_distil_reactor
 + Ve2_distil_reactor = 0
_C516: - D_distil_flash + Din_distil_flash + Dout_distil_flash = 0
_C517: Din_distil_flash - 35.4 Psie_distil_flash <= 0
_C518: Dout_distil_flash + 35.4 Psie_distil_flash >= 35.4
_C519: Dout_distil_flash + 1000 Psie_distil_flash <= 1000
_C52: L_furnace_flash + 1000 Wx_furnace_flash <= 1000
_C520: 71909.6045198 Din_distil_flash - 2545600 Psie_distil_flash
 + Ve2_distil_flash = 0
_C521: - D_distil_comp + Din_distil_comp + Dout_distil_comp = 0
_C522: Din_distil_comp - 35.4 Psie_distil_comp <= 0
_C523: Dout_distil_comp + 35.4 Psie_distil_comp >= 35.4
_C524: Dout_distil_comp + 1000 Psie_distil_comp <= 1000
_C525: 153887.00565 Din_distil_comp - 5447600 Psie_distil_comp
 + Ve2_distil_comp = 0
_C526: Ve2_distil_distil = 0
_C527: - D_distil_store + Din_distil_store + Dout_distil_store = 0
_C528: Din_distil_store - 35.4 Psie_distil_store <= 0
_C529: Dout_distil_store + 35.4 Psie_distil_store >= 35.4
_C53: A_furnace_flash - 1000 Wy_furnace_flash <= 0
_C530: Dout_distil_store + 1000 Psie_distil_store <= 1000
_C531: 160447.288136 Din_distil_store - 5679834 Psie_distil_store
 + Ve2_distil_store = 0
_C532: - D_distil_ctrlroom + Din_distil_ctrlroom + Dout_distil_ctrlroom = 0
_C533: Din_distil_ctrlroom - 35.4 Psie_distil_ctrlroom <= 0
_C534: Dout_distil_ctrlroom + 35.4 Psie_distil_ctrlroom >= 35.4
_C535: Dout_distil_ctrlroom + 1000 Psie_distil_ctrlroom <= 1000
_C536: 1085088.98305 Din_distil_ctrlroom - 38412150 Psie_distil_ctrlroom
 + Ve2_distil_ctrlroom = 0
_C537: - D_store_furnace + Din_store_furnace + Dout_store_furnace = 0
_C538: Din_store_furnace - 45 Psie_store_furnace <= 0
_C539: Dout_store_furnace + 45 Psie_store_furnace >= 45
_C54: B_furnace_flash + 1000 Wy_furnace_flash <= 1000
_C540: Dout_store_furnace + 1000 Psie_store_furnace <= 1000
_C541: 56806.6666667 Din_store_furnace - 2556300 Psie_store_furnace
 + Ve2_store_furnace = 0
_C542: - D_store_reactor + Din_store_reactor + Dout_store_reactor = 0
_C543: Din_store_reactor - 45 Psie_store_reactor <= 0
_C544: Dout_store_reactor + 45 Psie_store_reactor >= 45
_C545: Dout_store_reactor + 1000 Psie_store_reactor <= 1000
_C546: 86406.6666667 Din_store_reactor - 3888300 Psie_store_reactor
 + Ve2_store_reactor = 0
_C547: - D_store_flash + Din_store_flash + Dout_store_flash = 0
_C548: Din_store_flash - 45 Psie_store_flash <= 0
_C549: Dout_store_flash + 45 Psie_store_flash >= 45
_C55: - A_furnace_flash - B_furnace_flash + D_furnace_flash - L_furnace_flash
 - R_furnace_flash = 0
_C550: Dout_store_flash + 1000 Psie_store_flash <= 1000
_C551: 56568.8888889 Din_store_flash - 2545600 Psie_store_flash
 + Ve2_store_flash = 0
_C552: - D_store_comp + Din_store_comp + Dout_store_comp = 0
_C553: Din_store_comp - 45 Psie_store_comp <= 0
_C554: Dout_store_comp + 45 Psie_store_comp >= 45
_C555: Dout_store_comp + 1000 Psie_store_comp <= 1000
_C556: 121057.777778 Din_store_comp - 5447600 Psie_store_comp + Ve2_store_comp
 = 0
_C557: - D_store_distil + Din_store_distil + Dout_store_distil = 0
_C558: Din_store_distil - 45 Psie_store_distil <= 0
_C559: Dout_store_distil + 45 Psie_store_distil >= 45
_C56: - D_flash_furnace + D_furnace_flash = 0
_C560: Dout_store_distil + 1000 Psie_store_distil <= 1000
_C561: 157280 Din_store_distil - 7077600 Psie_store_distil + Ve2_store_distil
 = 0
_C562: Ve2_store_store = 0
_C563: - D_store_ctrlroom + Din_store_ctrlroom + Dout_store_ctrlroom = 0
_C564: Din_store_ctrlroom - 45 Psie_store_ctrlroom <= 0
_C565: Dout_store_ctrlroom + 45 Psie_store_ctrlroom >= 45
_C566: Dout_store_ctrlroom + 1000 Psie_store_ctrlroom <= 1000
_C567: 853603.333333 Din_store_ctrlroom - 38412150 Psie_store_ctrlroom
 + Ve2_store_ctrlroom = 0
_C568: - 0.001 Ve2_furnace_comp - 0.001 Ve2_furnace_ctrlroom
 - 0.001 Ve2_furnace_distil - 0.001 Ve2_furnace_flash
 - 0.001 Ve2_furnace_furnace - 0.001 Ve2_furnace_reactor
 - 0.001 Ve2_furnace_store + Ve_furnace = 0
_C569: Omega0_furnace - 0.75 Ve_furnace = 0
_C57: 1000 E1_furnace_flash + 1000 E2_furnace_flash - X_flash + X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= 0
_C570: - 0.001 Ve2_reactor_comp - 0.001 Ve2_reactor_ctrlroom
 - 0.001 Ve2_reactor_distil - 0.001 Ve2_reactor_flash
 - 0.001 Ve2_reactor_furnace - 0.001 Ve2_reactor_reactor
 - 0.001 Ve2_reactor_store + Ve_reactor = 0
_C571: Omega0_reactor - 0.8 Ve_reactor = 0
_C572: - 0.001 Ve2_distil_comp - 0.001 Ve2_distil_ctrlroom
 - 0.001 Ve2_distil_distil - 0.001 Ve2_distil_flash - 0.001 Ve2_distil_furnace
 - 0.001 Ve2_distil_reactor - 0.001 Ve2_distil_store + Ve_distil = 0
_C573: Omega0_distil - 0.77 Ve_distil = 0
_C574: - 0.001 Ve2_store_comp - 0.001 Ve2_store_ctrlroom
 - 0.001 Ve2_store_distil - 0.001 Ve2_store_flash - 0.001 Ve2_store_furnace
 - 0.001 Ve2_store_reactor - 0.001 Ve2_store_store + Ve_store = 0
_C575: Omega0_store - 0.82 Ve_store = 0
_C576: - Omega0_distil - Omega0_furnace - Omega0_reactor - Omega0_store
 + SumOmega = 0
_C58: - 1000 E1_furnace_flash + 1000 E2_furnace_flash + X_flash - X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= -1000
_C59: 1000 E1_furnace_flash - 1000 E2_furnace_flash - Y_flash + Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -1000
_C6: d_reactor + l_reactor = 14
_C60: - 1000 E1_furnace_flash - 1000 E2_furnace_flash + Y_flash - Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -2000
_C61: L_furnace_flash - R_flash_furnace = 0
_C62: - L_flash_furnace + R_furnace_flash = 0
_C63: A_furnace_flash - B_flash_furnace = 0
_C64: - A_flash_furnace + B_furnace_flash = 0
_C65: Wx_flash_furnace + Wx_furnace_flash = 1
_C66: Wy_flash_furnace + Wy_furnace_flash = 1
_C67: CD_reactor_flash - 1231.60506417 D_reactor_flash = 0
_C68: - L_reactor_flash + R_reactor_flash + X_flash - X_reactor = 0
_C69: A_reactor_flash - B_reactor_flash + Y_flash - Y_reactor = 0
_C7: X_reactor - 0.5 l_reactor >= 0
_C70: R_reactor_flash - 1000 Wx_reactor_flash <= 0
_C71: L_reactor_flash + 1000 Wx_reactor_flash <= 1000
_C72: A_reactor_flash - 1000 Wy_reactor_flash <= 0
_C73: B_reactor_flash + 1000 Wy_reactor_flash <= 1000
_C74: - A_reactor_flash - B_reactor_flash + D_reactor_flash - L_reactor_flash
 - R_reactor_flash = 0
_C75: - D_flash_reactor + D_reactor_flash = 0
_C76: 1000 E1_reactor_flash + 1000 E2_reactor_flash - X_flash + X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= 0
_C77: - 1000 E1_reactor_flash + 1000 E2_reactor_flash + X_flash - X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= -1000
_C78: 1000 E1_reactor_flash - 1000 E2_reactor_flash - Y_flash + Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -1000
_C79: - 1000 E1_reactor_flash - 1000 E2_reactor_flash + Y_flash - Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -2000
_C8: Y_reactor - 0.5 d_reactor >= 0
_C80: L_reactor_flash - R_flash_reactor = 0
_C81: - L_flash_reactor + R_reactor_flash = 0
_C82: A_reactor_flash - B_flash_reactor = 0
_C83: - A_flash_reactor + B_reactor_flash = 0
_C84: Wx_flash_reactor + Wx_reactor_flash = 1
_C85: Wy_flash_reactor + Wy_reactor_flash = 1
_C86: CD_furnace_comp - 1487.20883073 D_furnace_comp = 0
_C87: - L_furnace_comp + R_furnace_comp + X_comp - X_furnace = 0
_C88: A_furnace_comp - B_furnace_comp + Y_comp - Y_furnace = 0
_C89: R_furnace_comp - 1000 Wx_furnace_comp <= 0
_C9: 0 O_flash + l_flash = 4
_C90: L_furnace_comp + 1000 Wx_furnace_comp <= 1000
_C91: A_furnace_comp - 1000 Wy_furnace_comp <= 0
_C92: B_furnace_comp + 1000 Wy_furnace_comp <= 1000
_C93: - A_furnace_comp - B_furnace_comp + D_furnace_comp - L_furnace_comp
 - R_furnace_comp = 0
_C94: - D_comp_furnace + D_furnace_comp = 0
_C95: 1000 E1_furnace_comp + 1000 E2_furnace_comp - X_comp + X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= 0
_C96: - 1000 E1_furnace_comp + 1000 E2_furnace_comp + X_comp - X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= -1000
_C97: 1000 E1_furnace_comp - 1000 E2_furnace_comp - Y_comp + Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -1000
_C98: - 1000 E1_furnace_comp - 1000 E2_furnace_comp + Y_comp - Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -2000
_C99: L_furnace_comp - R_comp_furnace = 0
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
