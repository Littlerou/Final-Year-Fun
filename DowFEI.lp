\* Layout_Problem_Model *\
Minimize
OBJ: SumCD + SumOmega
Subject To
_C1: 9 O_furnace + l_furnace = 16
_C10: d_flash + l_flash = 8
_C100: CD_furnace_comp - 1487.20883073 D_furnace_comp = 0
_C101: 1000 E1_furnace_comp + 1000 E2_furnace_comp - X_comp + X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= 0
_C102: - 1000 E1_furnace_comp + 1000 E2_furnace_comp + X_comp - X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= -1000
_C103: 1000 E1_furnace_comp - 1000 E2_furnace_comp - Y_comp + Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -1000
_C104: - 1000 E1_furnace_comp - 1000 E2_furnace_comp + Y_comp - Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -2000
_C105: - L_reactor_comp + R_reactor_comp + X_comp - X_reactor = 0
_C106: A_reactor_comp - B_reactor_comp + Y_comp - Y_reactor = 0
_C107: R_reactor_comp - 1000 Wx_reactor_comp <= 0
_C108: L_reactor_comp + 1000 Wx_reactor_comp <= 1000
_C109: A_reactor_comp - 1000 Wy_reactor_comp <= 0
_C11: X_flash - 0.5 l_flash >= 0
_C110: B_reactor_comp + 1000 Wy_reactor_comp <= 1000
_C111: - A_reactor_comp - B_reactor_comp + D_reactor_comp - L_reactor_comp
 - R_reactor_comp = 0
_C112: - D_comp_reactor + D_reactor_comp = 0
_C113: L_reactor_comp - R_comp_reactor = 0
_C114: - L_comp_reactor + R_reactor_comp = 0
_C115: A_reactor_comp - B_comp_reactor = 0
_C116: - A_comp_reactor + B_reactor_comp = 0
_C117: Wx_comp_reactor + Wx_reactor_comp = 1
_C118: Wy_comp_reactor + Wy_reactor_comp = 1
_C119: CD_reactor_comp = 0
_C12: Y_flash - 0.5 d_flash >= 0
_C120: 1000 E1_reactor_comp + 1000 E2_reactor_comp - X_comp + X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= 0
_C121: - 1000 E1_reactor_comp + 1000 E2_reactor_comp + X_comp - X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= -1000
_C122: 1000 E1_reactor_comp - 1000 E2_reactor_comp - Y_comp + Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -1000
_C123: - 1000 E1_reactor_comp - 1000 E2_reactor_comp + Y_comp - Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -2000
_C124: - L_flash_comp + R_flash_comp + X_comp - X_flash = 0
_C125: A_flash_comp - B_flash_comp + Y_comp - Y_flash = 0
_C126: R_flash_comp - 1000 Wx_flash_comp <= 0
_C127: L_flash_comp + 1000 Wx_flash_comp <= 1000
_C128: A_flash_comp - 1000 Wy_flash_comp <= 0
_C129: B_flash_comp + 1000 Wy_flash_comp <= 1000
_C13: 0 O_comp + l_comp = 14
_C130: - A_flash_comp - B_flash_comp + D_flash_comp - L_flash_comp
 - R_flash_comp = 0
_C131: - D_comp_flash + D_flash_comp = 0
_C132: L_flash_comp - R_comp_flash = 0
_C133: - L_comp_flash + R_flash_comp = 0
_C134: A_flash_comp - B_comp_flash = 0
_C135: - A_comp_flash + B_flash_comp = 0
_C136: Wx_comp_flash + Wx_flash_comp = 1
_C137: Wy_comp_flash + Wy_flash_comp = 1
_C138: CD_flash_comp - 842.70764367 D_flash_comp = 0
_C139: 1000 E1_flash_comp + 1000 E2_flash_comp - X_comp + X_flash - 0.5 l_comp
 - 0.5 l_flash >= 0
_C14: d_comp + l_comp = 28
_C140: - 1000 E1_flash_comp + 1000 E2_flash_comp + X_comp - X_flash
 - 0.5 l_comp - 0.5 l_flash >= -1000
_C141: 1000 E1_flash_comp - 1000 E2_flash_comp - Y_comp + Y_flash - 0.5 d_comp
 - 0.5 d_flash >= -1000
_C142: - 1000 E1_flash_comp - 1000 E2_flash_comp + Y_comp - Y_flash
 - 0.5 d_comp - 0.5 d_flash >= -2000
_C143: - L_furnace_distil + R_furnace_distil + X_distil - X_furnace = 0
_C144: A_furnace_distil - B_furnace_distil + Y_distil - Y_furnace = 0
_C145: R_furnace_distil - 1000 Wx_furnace_distil <= 0
_C146: L_furnace_distil + 1000 Wx_furnace_distil <= 1000
_C147: A_furnace_distil - 1000 Wy_furnace_distil <= 0
_C148: B_furnace_distil + 1000 Wy_furnace_distil <= 1000
_C149: - A_furnace_distil - B_furnace_distil + D_furnace_distil
 - L_furnace_distil - R_furnace_distil = 0
_C15: X_comp - 0.5 l_comp >= 0
_C150: - D_distil_furnace + D_furnace_distil = 0
_C151: L_furnace_distil - R_distil_furnace = 0
_C152: - L_distil_furnace + R_furnace_distil = 0
_C153: A_furnace_distil - B_distil_furnace = 0
_C154: - A_distil_furnace + B_furnace_distil = 0
_C155: Wx_distil_furnace + Wx_furnace_distil = 1
_C156: Wy_distil_furnace + Wy_furnace_distil = 1
_C157: CD_furnace_distil - 52.9583379413 D_furnace_distil = 0
_C158: 1000 E1_furnace_distil + 1000 E2_furnace_distil - X_distil + X_furnace
 - 0.5 l_distil - 0.5 l_furnace >= 0
_C159: - 1000 E1_furnace_distil + 1000 E2_furnace_distil + X_distil
 - X_furnace - 0.5 l_distil - 0.5 l_furnace >= -1000
_C16: Y_comp - 0.5 d_comp >= 0
_C160: 1000 E1_furnace_distil - 1000 E2_furnace_distil - Y_distil + Y_furnace
 - 0.5 d_distil - 0.5 d_furnace >= -1000
_C161: - 1000 E1_furnace_distil - 1000 E2_furnace_distil + Y_distil
 - Y_furnace - 0.5 d_distil - 0.5 d_furnace >= -2000
_C162: - L_reactor_distil + R_reactor_distil + X_distil - X_reactor = 0
_C163: A_reactor_distil - B_reactor_distil + Y_distil - Y_reactor = 0
_C164: R_reactor_distil - 1000 Wx_reactor_distil <= 0
_C165: L_reactor_distil + 1000 Wx_reactor_distil <= 1000
_C166: A_reactor_distil - 1000 Wy_reactor_distil <= 0
_C167: B_reactor_distil + 1000 Wy_reactor_distil <= 1000
_C168: - A_reactor_distil - B_reactor_distil + D_reactor_distil
 - L_reactor_distil - R_reactor_distil = 0
_C169: - D_distil_reactor + D_reactor_distil = 0
_C17: 0 O_distil + l_distil = 8
_C170: L_reactor_distil - R_distil_reactor = 0
_C171: - L_distil_reactor + R_reactor_distil = 0
_C172: A_reactor_distil - B_distil_reactor = 0
_C173: - A_distil_reactor + B_reactor_distil = 0
_C174: Wx_distil_reactor + Wx_reactor_distil = 1
_C175: Wy_distil_reactor + Wy_reactor_distil = 1
_C176: CD_reactor_distil = 0
_C177: 1000 E1_reactor_distil + 1000 E2_reactor_distil - X_distil + X_reactor
 - 0.5 l_distil - 0.5 l_reactor >= 0
_C178: - 1000 E1_reactor_distil + 1000 E2_reactor_distil + X_distil
 - X_reactor - 0.5 l_distil - 0.5 l_reactor >= -1000
_C179: 1000 E1_reactor_distil - 1000 E2_reactor_distil - Y_distil + Y_reactor
 - 0.5 d_distil - 0.5 d_reactor >= -1000
_C18: d_distil + l_distil = 16
_C180: - 1000 E1_reactor_distil - 1000 E2_reactor_distil + Y_distil
 - Y_reactor - 0.5 d_distil - 0.5 d_reactor >= -2000
_C181: - L_flash_distil + R_flash_distil + X_distil - X_flash = 0
_C182: A_flash_distil - B_flash_distil + Y_distil - Y_flash = 0
_C183: R_flash_distil - 1000 Wx_flash_distil <= 0
_C184: L_flash_distil + 1000 Wx_flash_distil <= 1000
_C185: A_flash_distil - 1000 Wy_flash_distil <= 0
_C186: B_flash_distil + 1000 Wy_flash_distil <= 1000
_C187: - A_flash_distil - B_flash_distil + D_flash_distil - L_flash_distil
 - R_flash_distil = 0
_C188: - D_distil_flash + D_flash_distil = 0
_C189: L_flash_distil - R_distil_flash = 0
_C19: X_distil - 0.5 l_distil >= 0
_C190: - L_distil_flash + R_flash_distil = 0
_C191: A_flash_distil - B_distil_flash = 0
_C192: - A_distil_flash + B_flash_distil = 0
_C193: Wx_distil_flash + Wx_flash_distil = 1
_C194: Wy_distil_flash + Wy_flash_distil = 1
_C195: CD_flash_distil - 78.4818701125 D_flash_distil = 0
_C196: 1000 E1_flash_distil + 1000 E2_flash_distil - X_distil + X_flash
 - 0.5 l_distil - 0.5 l_flash >= 0
_C197: - 1000 E1_flash_distil + 1000 E2_flash_distil + X_distil - X_flash
 - 0.5 l_distil - 0.5 l_flash >= -1000
_C198: 1000 E1_flash_distil - 1000 E2_flash_distil - Y_distil + Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -1000
_C199: - 1000 E1_flash_distil - 1000 E2_flash_distil + Y_distil - Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -2000
_C2: d_furnace + l_furnace = 23
_C20: Y_distil - 0.5 d_distil >= 0
_C200: - L_comp_distil + R_comp_distil - X_comp + X_distil = 0
_C201: A_comp_distil - B_comp_distil - Y_comp + Y_distil = 0
_C202: R_comp_distil - 1000 Wx_comp_distil <= 0
_C203: L_comp_distil + 1000 Wx_comp_distil <= 1000
_C204: A_comp_distil - 1000 Wy_comp_distil <= 0
_C205: B_comp_distil + 1000 Wy_comp_distil <= 1000
_C206: - A_comp_distil - B_comp_distil + D_comp_distil - L_comp_distil
 - R_comp_distil = 0
_C207: D_comp_distil - D_distil_comp = 0
_C208: L_comp_distil - R_distil_comp = 0
_C209: - L_distil_comp + R_comp_distil = 0
_C21: - 10 O_store + l_store = 12
_C210: A_comp_distil - B_distil_comp = 0
_C211: - A_distil_comp + B_comp_distil = 0
_C212: Wx_comp_distil + Wx_distil_comp = 1
_C213: Wy_comp_distil + Wy_distil_comp = 1
_C214: CD_comp_distil = 0
_C215: 1000 E1_comp_distil + 1000 E2_comp_distil + X_comp - X_distil
 - 0.5 l_comp - 0.5 l_distil >= 0
_C216: - 1000 E1_comp_distil + 1000 E2_comp_distil - X_comp + X_distil
 - 0.5 l_comp - 0.5 l_distil >= -1000
_C217: 1000 E1_comp_distil - 1000 E2_comp_distil + Y_comp - Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -1000
_C218: - 1000 E1_comp_distil - 1000 E2_comp_distil - Y_comp + Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -2000
_C219: - L_furnace_store + R_furnace_store - X_furnace + X_store = 0
_C22: d_store + l_store = 34
_C220: A_furnace_store - B_furnace_store - Y_furnace + Y_store = 0
_C221: R_furnace_store - 1000 Wx_furnace_store <= 0
_C222: L_furnace_store + 1000 Wx_furnace_store <= 1000
_C223: A_furnace_store - 1000 Wy_furnace_store <= 0
_C224: B_furnace_store + 1000 Wy_furnace_store <= 1000
_C225: - A_furnace_store - B_furnace_store + D_furnace_store - L_furnace_store
 - R_furnace_store = 0
_C226: D_furnace_store - D_store_furnace = 0
_C227: L_furnace_store - R_store_furnace = 0
_C228: - L_store_furnace + R_furnace_store = 0
_C229: A_furnace_store - B_store_furnace = 0
_C23: X_store - 0.5 l_store >= 0
_C230: - A_store_furnace + B_furnace_store = 0
_C231: Wx_furnace_store + Wx_store_furnace = 1
_C232: Wy_furnace_store + Wy_store_furnace = 1
_C233: CD_furnace_store - 75.8045442191 D_furnace_store = 0
_C234: 1000 E1_furnace_store + 1000 E2_furnace_store + X_furnace - X_store
 - 0.5 l_furnace - 0.5 l_store >= 0
_C235: - 1000 E1_furnace_store + 1000 E2_furnace_store - X_furnace + X_store
 - 0.5 l_furnace - 0.5 l_store >= -1000
_C236: 1000 E1_furnace_store - 1000 E2_furnace_store + Y_furnace - Y_store
 - 0.5 d_furnace - 0.5 d_store >= -1000
_C237: - 1000 E1_furnace_store - 1000 E2_furnace_store - Y_furnace + Y_store
 - 0.5 d_furnace - 0.5 d_store >= -2000
_C238: - L_reactor_store + R_reactor_store - X_reactor + X_store = 0
_C239: A_reactor_store - B_reactor_store - Y_reactor + Y_store = 0
_C24: Y_store - 0.5 d_store >= 0
_C240: R_reactor_store - 1000 Wx_reactor_store <= 0
_C241: L_reactor_store + 1000 Wx_reactor_store <= 1000
_C242: A_reactor_store - 1000 Wy_reactor_store <= 0
_C243: B_reactor_store + 1000 Wy_reactor_store <= 1000
_C244: - A_reactor_store - B_reactor_store + D_reactor_store - L_reactor_store
 - R_reactor_store = 0
_C245: D_reactor_store - D_store_reactor = 0
_C246: L_reactor_store - R_store_reactor = 0
_C247: - L_store_reactor + R_reactor_store = 0
_C248: A_reactor_store - B_store_reactor = 0
_C249: - A_store_reactor + B_reactor_store = 0
_C25: 0 O_ctrlroom + l_ctrlroom = 15
_C250: Wx_reactor_store + Wx_store_reactor = 1
_C251: Wy_reactor_store + Wy_store_reactor = 1
_C252: CD_reactor_store = 0
_C253: 1000 E1_reactor_store + 1000 E2_reactor_store + X_reactor - X_store
 - 0.5 l_reactor - 0.5 l_store >= 0
_C254: - 1000 E1_reactor_store + 1000 E2_reactor_store - X_reactor + X_store
 - 0.5 l_reactor - 0.5 l_store >= -1000
_C255: 1000 E1_reactor_store - 1000 E2_reactor_store + Y_reactor - Y_store
 - 0.5 d_reactor - 0.5 d_store >= -1000
_C256: - 1000 E1_reactor_store - 1000 E2_reactor_store - Y_reactor + Y_store
 - 0.5 d_reactor - 0.5 d_store >= -2000
_C257: - L_flash_store + R_flash_store - X_flash + X_store = 0
_C258: A_flash_store - B_flash_store - Y_flash + Y_store = 0
_C259: R_flash_store - 1000 Wx_flash_store <= 0
_C26: d_ctrlroom + l_ctrlroom = 30
_C260: L_flash_store + 1000 Wx_flash_store <= 1000
_C261: A_flash_store - 1000 Wy_flash_store <= 0
_C262: B_flash_store + 1000 Wy_flash_store <= 1000
_C263: - A_flash_store - B_flash_store + D_flash_store - L_flash_store
 - R_flash_store = 0
_C264: D_flash_store - D_store_flash = 0
_C265: L_flash_store - R_store_flash = 0
_C266: - L_store_flash + R_flash_store = 0
_C267: A_flash_store - B_store_flash = 0
_C268: - A_store_flash + B_flash_store = 0
_C269: Wx_flash_store + Wx_store_flash = 1
_C27: X_ctrlroom - 0.5 l_ctrlroom >= 0
_C270: Wy_flash_store + Wy_store_flash = 1
_C271: CD_flash_store = 0
_C272: 1000 E1_flash_store + 1000 E2_flash_store + X_flash - X_store
 - 0.5 l_flash - 0.5 l_store >= 0
_C273: - 1000 E1_flash_store + 1000 E2_flash_store - X_flash + X_store
 - 0.5 l_flash - 0.5 l_store >= -1000
_C274: 1000 E1_flash_store - 1000 E2_flash_store + Y_flash - Y_store
 - 0.5 d_flash - 0.5 d_store >= -1000
_C275: - 1000 E1_flash_store - 1000 E2_flash_store - Y_flash + Y_store
 - 0.5 d_flash - 0.5 d_store >= -2000
_C276: - L_comp_store + R_comp_store - X_comp + X_store = 0
_C277: A_comp_store - B_comp_store - Y_comp + Y_store = 0
_C278: R_comp_store - 1000 Wx_comp_store <= 0
_C279: L_comp_store + 1000 Wx_comp_store <= 1000
_C28: Y_ctrlroom - 0.5 d_ctrlroom >= 0
_C280: A_comp_store - 1000 Wy_comp_store <= 0
_C281: B_comp_store + 1000 Wy_comp_store <= 1000
_C282: - A_comp_store - B_comp_store + D_comp_store - L_comp_store
 - R_comp_store = 0
_C283: D_comp_store - D_store_comp = 0
_C284: L_comp_store - R_store_comp = 0
_C285: - L_store_comp + R_comp_store = 0
_C286: A_comp_store - B_store_comp = 0
_C287: - A_store_comp + B_comp_store = 0
_C288: Wx_comp_store + Wx_store_comp = 1
_C289: Wy_comp_store + Wy_store_comp = 1
_C29: - L_furnace_reactor + R_furnace_reactor - X_furnace + X_reactor = 0
_C290: CD_comp_store = 0
_C291: 1000 E1_comp_store + 1000 E2_comp_store + X_comp - X_store - 0.5 l_comp
 - 0.5 l_store >= 0
_C292: - 1000 E1_comp_store + 1000 E2_comp_store - X_comp + X_store
 - 0.5 l_comp - 0.5 l_store >= -1000
_C293: 1000 E1_comp_store - 1000 E2_comp_store + Y_comp - Y_store - 0.5 d_comp
 - 0.5 d_store >= -1000
_C294: - 1000 E1_comp_store - 1000 E2_comp_store - Y_comp + Y_store
 - 0.5 d_comp - 0.5 d_store >= -2000
_C295: - L_distil_store + R_distil_store - X_distil + X_store = 0
_C296: A_distil_store - B_distil_store - Y_distil + Y_store = 0
_C297: R_distil_store - 1000 Wx_distil_store <= 0
_C298: L_distil_store + 1000 Wx_distil_store <= 1000
_C299: A_distil_store - 1000 Wy_distil_store <= 0
_C3: X_furnace - 0.5 l_furnace >= 0
_C30: A_furnace_reactor - B_furnace_reactor - Y_furnace + Y_reactor = 0
_C300: B_distil_store + 1000 Wy_distil_store <= 1000
_C301: - A_distil_store - B_distil_store + D_distil_store - L_distil_store
 - R_distil_store = 0
_C302: D_distil_store - D_store_distil = 0
_C303: L_distil_store - R_store_distil = 0
_C304: - L_store_distil + R_distil_store = 0
_C305: A_distil_store - B_store_distil = 0
_C306: - A_store_distil + B_distil_store = 0
_C307: Wx_distil_store + Wx_store_distil = 1
_C308: Wy_distil_store + Wy_store_distil = 1
_C309: CD_distil_store - 70.4511911446 D_distil_store = 0
_C31: R_furnace_reactor - 1000 Wx_furnace_reactor <= 0
_C310: 1000 E1_distil_store + 1000 E2_distil_store + X_distil - X_store
 - 0.5 l_distil - 0.5 l_store >= 0
_C311: - 1000 E1_distil_store + 1000 E2_distil_store - X_distil + X_store
 - 0.5 l_distil - 0.5 l_store >= -1000
_C312: 1000 E1_distil_store - 1000 E2_distil_store + Y_distil - Y_store
 - 0.5 d_distil - 0.5 d_store >= -1000
_C313: - 1000 E1_distil_store - 1000 E2_distil_store - Y_distil + Y_store
 - 0.5 d_distil - 0.5 d_store >= -2000
_C314: - L_furnace_ctrlroom + R_furnace_ctrlroom + X_ctrlroom - X_furnace = 0
_C315: A_furnace_ctrlroom - B_furnace_ctrlroom + Y_ctrlroom - Y_furnace = 0
_C316: R_furnace_ctrlroom - 1000 Wx_furnace_ctrlroom <= 0
_C317: L_furnace_ctrlroom + 1000 Wx_furnace_ctrlroom <= 1000
_C318: A_furnace_ctrlroom - 1000 Wy_furnace_ctrlroom <= 0
_C319: B_furnace_ctrlroom + 1000 Wy_furnace_ctrlroom <= 1000
_C32: L_furnace_reactor + 1000 Wx_furnace_reactor <= 1000
_C320: - A_furnace_ctrlroom - B_furnace_ctrlroom + D_furnace_ctrlroom
 - L_furnace_ctrlroom - R_furnace_ctrlroom = 0
_C321: - D_ctrlroom_furnace + D_furnace_ctrlroom = 0
_C322: L_furnace_ctrlroom - R_ctrlroom_furnace = 0
_C323: - L_ctrlroom_furnace + R_furnace_ctrlroom = 0
_C324: A_furnace_ctrlroom - B_ctrlroom_furnace = 0
_C325: - A_ctrlroom_furnace + B_furnace_ctrlroom = 0
_C326: Wx_ctrlroom_furnace + Wx_furnace_ctrlroom = 1
_C327: Wy_ctrlroom_furnace + Wy_furnace_ctrlroom = 1
_C328: CD_furnace_ctrlroom = 0
_C329: 1000 E1_furnace_ctrlroom + 1000 E2_furnace_ctrlroom - X_ctrlroom
 + X_furnace - 0.5 l_ctrlroom - 0.5 l_furnace >= 0
_C33: A_furnace_reactor - 1000 Wy_furnace_reactor <= 0
_C330: - 1000 E1_furnace_ctrlroom + 1000 E2_furnace_ctrlroom + X_ctrlroom
 - X_furnace - 0.5 l_ctrlroom - 0.5 l_furnace >= -1000
_C331: 1000 E1_furnace_ctrlroom - 1000 E2_furnace_ctrlroom - Y_ctrlroom
 + Y_furnace - 0.5 d_ctrlroom - 0.5 d_furnace >= -1000
_C332: - 1000 E1_furnace_ctrlroom - 1000 E2_furnace_ctrlroom + Y_ctrlroom
 - Y_furnace - 0.5 d_ctrlroom - 0.5 d_furnace >= -2000
_C333: - L_reactor_ctrlroom + R_reactor_ctrlroom + X_ctrlroom - X_reactor = 0
_C334: A_reactor_ctrlroom - B_reactor_ctrlroom + Y_ctrlroom - Y_reactor = 0
_C335: R_reactor_ctrlroom - 1000 Wx_reactor_ctrlroom <= 0
_C336: L_reactor_ctrlroom + 1000 Wx_reactor_ctrlroom <= 1000
_C337: A_reactor_ctrlroom - 1000 Wy_reactor_ctrlroom <= 0
_C338: B_reactor_ctrlroom + 1000 Wy_reactor_ctrlroom <= 1000
_C339: - A_reactor_ctrlroom - B_reactor_ctrlroom + D_reactor_ctrlroom
 - L_reactor_ctrlroom - R_reactor_ctrlroom = 0
_C34: B_furnace_reactor + 1000 Wy_furnace_reactor <= 1000
_C340: - D_ctrlroom_reactor + D_reactor_ctrlroom = 0
_C341: L_reactor_ctrlroom - R_ctrlroom_reactor = 0
_C342: - L_ctrlroom_reactor + R_reactor_ctrlroom = 0
_C343: A_reactor_ctrlroom - B_ctrlroom_reactor = 0
_C344: - A_ctrlroom_reactor + B_reactor_ctrlroom = 0
_C345: Wx_ctrlroom_reactor + Wx_reactor_ctrlroom = 1
_C346: Wy_ctrlroom_reactor + Wy_reactor_ctrlroom = 1
_C347: CD_reactor_ctrlroom = 0
_C348: 1000 E1_reactor_ctrlroom + 1000 E2_reactor_ctrlroom - X_ctrlroom
 + X_reactor - 0.5 l_ctrlroom - 0.5 l_reactor >= 0
_C349: - 1000 E1_reactor_ctrlroom + 1000 E2_reactor_ctrlroom + X_ctrlroom
 - X_reactor - 0.5 l_ctrlroom - 0.5 l_reactor >= -1000
_C35: - A_furnace_reactor - B_furnace_reactor + D_furnace_reactor
 - L_furnace_reactor - R_furnace_reactor = 0
_C350: 1000 E1_reactor_ctrlroom - 1000 E2_reactor_ctrlroom - Y_ctrlroom
 + Y_reactor - 0.5 d_ctrlroom - 0.5 d_reactor >= -1000
_C351: - 1000 E1_reactor_ctrlroom - 1000 E2_reactor_ctrlroom + Y_ctrlroom
 - Y_reactor - 0.5 d_ctrlroom - 0.5 d_reactor >= -2000
_C352: - L_flash_ctrlroom + R_flash_ctrlroom + X_ctrlroom - X_flash = 0
_C353: A_flash_ctrlroom - B_flash_ctrlroom + Y_ctrlroom - Y_flash = 0
_C354: R_flash_ctrlroom - 1000 Wx_flash_ctrlroom <= 0
_C355: L_flash_ctrlroom + 1000 Wx_flash_ctrlroom <= 1000
_C356: A_flash_ctrlroom - 1000 Wy_flash_ctrlroom <= 0
_C357: B_flash_ctrlroom + 1000 Wy_flash_ctrlroom <= 1000
_C358: - A_flash_ctrlroom - B_flash_ctrlroom + D_flash_ctrlroom
 - L_flash_ctrlroom - R_flash_ctrlroom = 0
_C359: - D_ctrlroom_flash + D_flash_ctrlroom = 0
_C36: D_furnace_reactor - D_reactor_furnace = 0
_C360: L_flash_ctrlroom - R_ctrlroom_flash = 0
_C361: - L_ctrlroom_flash + R_flash_ctrlroom = 0
_C362: A_flash_ctrlroom - B_ctrlroom_flash = 0
_C363: - A_ctrlroom_flash + B_flash_ctrlroom = 0
_C364: Wx_ctrlroom_flash + Wx_flash_ctrlroom = 1
_C365: Wy_ctrlroom_flash + Wy_flash_ctrlroom = 1
_C366: CD_flash_ctrlroom = 0
_C367: 1000 E1_flash_ctrlroom + 1000 E2_flash_ctrlroom - X_ctrlroom + X_flash
 - 0.5 l_ctrlroom - 0.5 l_flash >= 0
_C368: - 1000 E1_flash_ctrlroom + 1000 E2_flash_ctrlroom + X_ctrlroom
 - X_flash - 0.5 l_ctrlroom - 0.5 l_flash >= -1000
_C369: 1000 E1_flash_ctrlroom - 1000 E2_flash_ctrlroom - Y_ctrlroom + Y_flash
 - 0.5 d_ctrlroom - 0.5 d_flash >= -1000
_C37: L_furnace_reactor - R_reactor_furnace = 0
_C370: - 1000 E1_flash_ctrlroom - 1000 E2_flash_ctrlroom + Y_ctrlroom
 - Y_flash - 0.5 d_ctrlroom - 0.5 d_flash >= -2000
_C371: - L_comp_ctrlroom + R_comp_ctrlroom - X_comp + X_ctrlroom = 0
_C372: A_comp_ctrlroom - B_comp_ctrlroom - Y_comp + Y_ctrlroom = 0
_C373: R_comp_ctrlroom - 1000 Wx_comp_ctrlroom <= 0
_C374: L_comp_ctrlroom + 1000 Wx_comp_ctrlroom <= 1000
_C375: A_comp_ctrlroom - 1000 Wy_comp_ctrlroom <= 0
_C376: B_comp_ctrlroom + 1000 Wy_comp_ctrlroom <= 1000
_C377: - A_comp_ctrlroom - B_comp_ctrlroom + D_comp_ctrlroom - L_comp_ctrlroom
 - R_comp_ctrlroom = 0
_C378: D_comp_ctrlroom - D_ctrlroom_comp = 0
_C379: L_comp_ctrlroom - R_ctrlroom_comp = 0
_C38: - L_reactor_furnace + R_furnace_reactor = 0
_C380: - L_ctrlroom_comp + R_comp_ctrlroom = 0
_C381: A_comp_ctrlroom - B_ctrlroom_comp = 0
_C382: - A_ctrlroom_comp + B_comp_ctrlroom = 0
_C383: Wx_comp_ctrlroom + Wx_ctrlroom_comp = 1
_C384: Wy_comp_ctrlroom + Wy_ctrlroom_comp = 1
_C385: CD_comp_ctrlroom = 0
_C386: 1000 E1_comp_ctrlroom + 1000 E2_comp_ctrlroom + X_comp - X_ctrlroom
 - 0.5 l_comp - 0.5 l_ctrlroom >= 0
_C387: - 1000 E1_comp_ctrlroom + 1000 E2_comp_ctrlroom - X_comp + X_ctrlroom
 - 0.5 l_comp - 0.5 l_ctrlroom >= -1000
_C388: 1000 E1_comp_ctrlroom - 1000 E2_comp_ctrlroom + Y_comp - Y_ctrlroom
 - 0.5 d_comp - 0.5 d_ctrlroom >= -1000
_C389: - 1000 E1_comp_ctrlroom - 1000 E2_comp_ctrlroom - Y_comp + Y_ctrlroom
 - 0.5 d_comp - 0.5 d_ctrlroom >= -2000
_C39: A_furnace_reactor - B_reactor_furnace = 0
_C390: - L_distil_ctrlroom + R_distil_ctrlroom + X_ctrlroom - X_distil = 0
_C391: A_distil_ctrlroom - B_distil_ctrlroom + Y_ctrlroom - Y_distil = 0
_C392: R_distil_ctrlroom - 1000 Wx_distil_ctrlroom <= 0
_C393: L_distil_ctrlroom + 1000 Wx_distil_ctrlroom <= 1000
_C394: A_distil_ctrlroom - 1000 Wy_distil_ctrlroom <= 0
_C395: B_distil_ctrlroom + 1000 Wy_distil_ctrlroom <= 1000
_C396: - A_distil_ctrlroom - B_distil_ctrlroom + D_distil_ctrlroom
 - L_distil_ctrlroom - R_distil_ctrlroom = 0
_C397: - D_ctrlroom_distil + D_distil_ctrlroom = 0
_C398: L_distil_ctrlroom - R_ctrlroom_distil = 0
_C399: - L_ctrlroom_distil + R_distil_ctrlroom = 0
_C4: Y_furnace - 0.5 d_furnace >= 0
_C40: - A_reactor_furnace + B_furnace_reactor = 0
_C400: A_distil_ctrlroom - B_ctrlroom_distil = 0
_C401: - A_ctrlroom_distil + B_distil_ctrlroom = 0
_C402: Wx_ctrlroom_distil + Wx_distil_ctrlroom = 1
_C403: Wy_ctrlroom_distil + Wy_distil_ctrlroom = 1
_C404: CD_distil_ctrlroom = 0
_C405: 1000 E1_distil_ctrlroom + 1000 E2_distil_ctrlroom - X_ctrlroom
 + X_distil - 0.5 l_ctrlroom - 0.5 l_distil >= 0
_C406: - 1000 E1_distil_ctrlroom + 1000 E2_distil_ctrlroom + X_ctrlroom
 - X_distil - 0.5 l_ctrlroom - 0.5 l_distil >= -1000
_C407: 1000 E1_distil_ctrlroom - 1000 E2_distil_ctrlroom - Y_ctrlroom
 + Y_distil - 0.5 d_ctrlroom - 0.5 d_distil >= -1000
_C408: - 1000 E1_distil_ctrlroom - 1000 E2_distil_ctrlroom + Y_ctrlroom
 - Y_distil - 0.5 d_ctrlroom - 0.5 d_distil >= -2000
_C409: - L_store_ctrlroom + R_store_ctrlroom + X_ctrlroom - X_store = 0
_C41: Wx_furnace_reactor + Wx_reactor_furnace = 1
_C410: A_store_ctrlroom - B_store_ctrlroom + Y_ctrlroom - Y_store = 0
_C411: R_store_ctrlroom - 1000 Wx_store_ctrlroom <= 0
_C412: L_store_ctrlroom + 1000 Wx_store_ctrlroom <= 1000
_C413: A_store_ctrlroom - 1000 Wy_store_ctrlroom <= 0
_C414: B_store_ctrlroom + 1000 Wy_store_ctrlroom <= 1000
_C415: - A_store_ctrlroom - B_store_ctrlroom + D_store_ctrlroom
 - L_store_ctrlroom - R_store_ctrlroom = 0
_C416: - D_ctrlroom_store + D_store_ctrlroom = 0
_C417: L_store_ctrlroom - R_ctrlroom_store = 0
_C418: - L_ctrlroom_store + R_store_ctrlroom = 0
_C419: A_store_ctrlroom - B_ctrlroom_store = 0
_C42: Wy_furnace_reactor + Wy_reactor_furnace = 1
_C420: - A_ctrlroom_store + B_store_ctrlroom = 0
_C421: Wx_ctrlroom_store + Wx_store_ctrlroom = 1
_C422: Wy_ctrlroom_store + Wy_store_ctrlroom = 1
_C423: CD_store_ctrlroom = 0
_C424: 1000 E1_store_ctrlroom + 1000 E2_store_ctrlroom - X_ctrlroom + X_store
 - 0.5 l_ctrlroom - 0.5 l_store >= 0
_C425: - 1000 E1_store_ctrlroom + 1000 E2_store_ctrlroom + X_ctrlroom
 - X_store - 0.5 l_ctrlroom - 0.5 l_store >= -1000
_C426: 1000 E1_store_ctrlroom - 1000 E2_store_ctrlroom - Y_ctrlroom + Y_store
 - 0.5 d_ctrlroom - 0.5 d_store >= -1000
_C427: - 1000 E1_store_ctrlroom - 1000 E2_store_ctrlroom + Y_ctrlroom
 - Y_store - 0.5 d_ctrlroom - 0.5 d_store >= -2000
_C428: - D_ctrlroom_furnace + D_road_ctrlroom_furnace = 0
_C429: - D_ctrlroom_reactor + D_road_ctrlroom_reactor = 0
_C43: CD_furnace_reactor - 1075.48971174 D_furnace_reactor = 0
_C430: - D_ctrlroom_flash + D_road_ctrlroom_flash = 0
_C431: - D_ctrlroom_comp + D_road_ctrlroom_comp = 0
_C432: - D_ctrlroom_distil + D_road_ctrlroom_distil = 0
_C433: - D_ctrlroom_store + D_road_ctrlroom_store = 0
_C434: - 15.428 D_road_ctrlroom_comp - 15.428 D_road_ctrlroom_ctrlroom
 - 15.428 D_road_ctrlroom_distil - 15.428 D_road_ctrlroom_flash
 - 15.428 D_road_ctrlroom_furnace - 15.428 D_road_ctrlroom_reactor
 - 15.428 D_road_ctrlroom_store + roadCost = 0
_C435: - CD_comp_comp - CD_comp_ctrlroom - CD_comp_distil - CD_comp_flash
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
_C436: - 0.8 X_furnace + Y_furnace + 0.5 d_furnace + 0.4 l_furnace <= 50
_C437: 0.8 X_furnace + Y_furnace + 0.5 d_furnace + 0.4 l_furnace <= 90
_C438: Y_furnace - 0.5 d_furnace >= 0
_C439: X_furnace + 0.5 l_furnace <= 50
_C44: 1000 E1_furnace_reactor + 1000 E2_furnace_reactor + X_furnace
 - X_reactor - 0.5 l_furnace - 0.5 l_reactor >= 0
_C440: X_furnace - 0.5 l_furnace >= 0
_C441: - 0.8 X_reactor + Y_reactor + 0.5 d_reactor + 0.4 l_reactor <= 50
_C442: 0.8 X_reactor + Y_reactor + 0.5 d_reactor + 0.4 l_reactor <= 90
_C443: Y_reactor - 0.5 d_reactor >= 0
_C444: X_reactor + 0.5 l_reactor <= 50
_C445: X_reactor - 0.5 l_reactor >= 0
_C446: - 0.8 X_flash + Y_flash + 0.5 d_flash + 0.4 l_flash <= 50
_C447: 0.8 X_flash + Y_flash + 0.5 d_flash + 0.4 l_flash <= 90
_C448: Y_flash - 0.5 d_flash >= 0
_C449: X_flash + 0.5 l_flash <= 50
_C45: - 1000 E1_furnace_reactor + 1000 E2_furnace_reactor - X_furnace
 + X_reactor - 0.5 l_furnace - 0.5 l_reactor >= -1000
_C450: X_flash - 0.5 l_flash >= 0
_C451: - 0.8 X_comp + Y_comp + 0.5 d_comp + 0.4 l_comp <= 50
_C452: 0.8 X_comp + Y_comp + 0.5 d_comp + 0.4 l_comp <= 90
_C453: Y_comp - 0.5 d_comp >= 0
_C454: X_comp + 0.5 l_comp <= 50
_C455: X_comp - 0.5 l_comp >= 0
_C456: - 0.8 X_distil + Y_distil + 0.5 d_distil + 0.4 l_distil <= 50
_C457: 0.8 X_distil + Y_distil + 0.5 d_distil + 0.4 l_distil <= 90
_C458: Y_distil - 0.5 d_distil >= 0
_C459: X_distil + 0.5 l_distil <= 50
_C46: 1000 E1_furnace_reactor - 1000 E2_furnace_reactor + Y_furnace
 - Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -1000
_C460: X_distil - 0.5 l_distil >= 0
_C461: - 0.8 X_store + Y_store + 0.5 d_store + 0.4 l_store <= 50
_C462: 0.8 X_store + Y_store + 0.5 d_store + 0.4 l_store <= 90
_C463: Y_store - 0.5 d_store >= 0
_C464: X_store + 0.5 l_store <= 50
_C465: X_store - 0.5 l_store >= 0
_C466: - 0.8 X_ctrlroom + Y_ctrlroom + 0.5 d_ctrlroom + 0.4 l_ctrlroom <= 50
_C467: 0.8 X_ctrlroom + Y_ctrlroom + 0.5 d_ctrlroom + 0.4 l_ctrlroom <= 90
_C468: Y_ctrlroom - 0.5 d_ctrlroom >= 0
_C469: X_ctrlroom + 0.5 l_ctrlroom <= 50
_C47: - 1000 E1_furnace_reactor - 1000 E2_furnace_reactor - Y_furnace
 + Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -2000
_C470: X_ctrlroom - 0.5 l_ctrlroom >= 0
_C471: Ve2_furnace_furnace = 0
_C472: - D_furnace_reactor + Din_furnace_reactor + Dout_furnace_reactor = 0
_C473: Din_furnace_reactor - 40.7293505963 Psie_furnace_reactor <= 0
_C474: Dout_furnace_reactor + 40.7293505963 Psie_furnace_reactor
 >= 40.7293505963
_C475: Dout_furnace_reactor + 1000 Psie_furnace_reactor <= 1000
_C476: 104190.711069 Din_furnace_reactor - 4243620 Psie_furnace_reactor
 + Ve2_furnace_reactor = 0
_C477: - D_furnace_flash + Din_furnace_flash + Dout_furnace_flash = 0
_C478: Din_furnace_flash - 40.7293505963 Psie_furnace_flash <= 0
_C479: Dout_furnace_flash + 40.7293505963 Psie_furnace_flash >= 40.7293505963
_C48: - L_furnace_flash + R_furnace_flash + X_flash - X_furnace = 0
_C480: Dout_furnace_flash + 1000 Psie_furnace_flash <= 1000
_C481: 67858.6807679 Din_furnace_flash - 2763840 Psie_furnace_flash
 + Ve2_furnace_flash = 0
_C482: - D_furnace_comp + Din_furnace_comp + Dout_furnace_comp = 0
_C483: Din_furnace_comp - 40.7293505963 Psie_furnace_comp <= 0
_C484: Dout_furnace_comp + 40.7293505963 Psie_furnace_comp >= 40.7293505963
_C485: Dout_furnace_comp + 1000 Psie_furnace_comp <= 1000
_C486: 157788.914036 Din_furnace_comp - 6426640 Psie_furnace_comp
 + Ve2_furnace_comp = 0
_C487: - D_furnace_distil + Din_furnace_distil + Dout_furnace_distil = 0
_C488: Din_furnace_distil - 40.7293505963 Psie_furnace_distil <= 0
_C489: Dout_furnace_distil + 40.7293505963 Psie_furnace_distil
 >= 40.7293505963
_C49: A_furnace_flash - B_furnace_flash + Y_flash - Y_furnace = 0
_C490: Dout_furnace_distil + 1000 Psie_furnace_distil <= 1000
_C491: 203996.378001 Din_furnace_distil - 8308640 Psie_furnace_distil
 + Ve2_furnace_distil = 0
_C492: - D_furnace_store + Din_furnace_store + Dout_furnace_store = 0
_C493: Din_furnace_store - 40.7293505963 Psie_furnace_store <= 0
_C494: Dout_furnace_store + 40.7293505963 Psie_furnace_store >= 40.7293505963
_C495: Dout_furnace_store + 1000 Psie_furnace_store <= 1000
_C496: 175592.487857 Din_furnace_store - 7151768 Psie_furnace_store
 + Ve2_furnace_store = 0
_C497: - D_furnace_ctrlroom + Din_furnace_ctrlroom + Dout_furnace_ctrlroom = 0
_C498: Din_furnace_ctrlroom - 40.7293505963 Psie_furnace_ctrlroom <= 0
_C499: Dout_furnace_ctrlroom + 40.7293505963 Psie_furnace_ctrlroom
 >= 40.7293505963
_C5: - 1.8 O_reactor + l_reactor = 5.8
_C50: R_furnace_flash - 1000 Wx_furnace_flash <= 0
_C500: Dout_furnace_ctrlroom + 1000 Psie_furnace_ctrlroom <= 1000
_C501: 950837.89535 Din_furnace_ctrlroom - 38727010 Psie_furnace_ctrlroom
 + Ve2_furnace_ctrlroom = 0
_C502: - D_reactor_furnace + Din_reactor_furnace + Dout_reactor_furnace = 0
_C503: Din_reactor_furnace - 57.1342279199 Psie_reactor_furnace <= 0
_C504: Dout_reactor_furnace + 57.1342279199 Psie_reactor_furnace
 >= 57.1342279199
_C505: Dout_reactor_furnace + 1000 Psie_reactor_furnace <= 1000
_C506: 48636.6946955 Din_reactor_furnace - 2778820 Psie_reactor_furnace
 + Ve2_reactor_furnace = 0
_C507: Ve2_reactor_reactor = 0
_C508: - D_reactor_flash + Din_reactor_flash + Dout_reactor_flash = 0
_C509: Din_reactor_flash - 57.1342279199 Psie_reactor_flash <= 0
_C51: L_furnace_flash + 1000 Wx_furnace_flash <= 1000
_C510: Dout_reactor_flash + 57.1342279199 Psie_reactor_flash >= 57.1342279199
_C511: Dout_reactor_flash + 1000 Psie_reactor_flash <= 1000
_C512: 48374.5051018 Din_reactor_flash - 2763840 Psie_reactor_flash
 + Ve2_reactor_flash = 0
_C513: - D_reactor_comp + Din_reactor_comp + Dout_reactor_comp = 0
_C514: Din_reactor_comp - 57.1342279199 Psie_reactor_comp <= 0
_C515: Dout_reactor_comp + 57.1342279199 Psie_reactor_comp >= 57.1342279199
_C516: Dout_reactor_comp + 1000 Psie_reactor_comp <= 1000
_C517: 112483.186244 Din_reactor_comp - 6426640 Psie_reactor_comp
 + Ve2_reactor_comp = 0
_C518: - D_reactor_distil + Din_reactor_distil + Dout_reactor_distil = 0
_C519: Din_reactor_distil - 57.1342279199 Psie_reactor_distil <= 0
_C52: A_furnace_flash - 1000 Wy_furnace_flash <= 0
_C520: Dout_reactor_distil + 57.1342279199 Psie_reactor_distil
 >= 57.1342279199
_C521: Dout_reactor_distil + 1000 Psie_reactor_distil <= 1000
_C522: 145423.160555 Din_reactor_distil - 8308640 Psie_reactor_distil
 + Ve2_reactor_distil = 0
_C523: - D_reactor_store + Din_reactor_store + Dout_reactor_store = 0
_C524: Din_reactor_store - 57.1342279199 Psie_reactor_store <= 0
_C525: Dout_reactor_store + 57.1342279199 Psie_reactor_store >= 57.1342279199
_C526: Dout_reactor_store + 1000 Psie_reactor_store <= 1000
_C527: 125174.842829 Din_reactor_store - 7151768 Psie_reactor_store
 + Ve2_reactor_store = 0
_C528: - D_reactor_ctrlroom + Din_reactor_ctrlroom + Dout_reactor_ctrlroom = 0
_C529: Din_reactor_ctrlroom - 57.1342279199 Psie_reactor_ctrlroom <= 0
_C53: B_furnace_flash + 1000 Wy_furnace_flash <= 1000
_C530: Dout_reactor_ctrlroom + 57.1342279199 Psie_reactor_ctrlroom
 >= 57.1342279199
_C531: Dout_reactor_ctrlroom + 1000 Psie_reactor_ctrlroom <= 1000
_C532: 677825.034309 Din_reactor_ctrlroom - 38727010 Psie_reactor_ctrlroom
 + Ve2_reactor_ctrlroom = 0
_C533: - D_distil_furnace + Din_distil_furnace + Dout_distil_furnace = 0
_C534: Din_distil_furnace - 50.063160108 Psie_distil_furnace <= 0
_C535: Dout_distil_furnace + 50.063160108 Psie_distil_furnace >= 50.063160108
_C536: Dout_distil_furnace + 1000 Psie_distil_furnace <= 1000
_C537: 55506.2843417 Din_distil_furnace - 2778820 Psie_distil_furnace
 + Ve2_distil_furnace = 0
_C538: - D_distil_reactor + Din_distil_reactor + Dout_distil_reactor = 0
_C539: Din_distil_reactor - 50.063160108 Psie_distil_reactor <= 0
_C54: - A_furnace_flash - B_furnace_flash + D_furnace_flash - L_furnace_flash
 - R_furnace_flash = 0
_C540: Dout_distil_reactor + 50.063160108 Psie_distil_reactor >= 50.063160108
_C541: Dout_distil_reactor + 1000 Psie_distil_reactor <= 1000
_C542: 84765.3242593 Din_distil_reactor - 4243620 Psie_distil_reactor
 + Ve2_distil_reactor = 0
_C543: - D_distil_flash + Din_distil_flash + Dout_distil_flash = 0
_C544: Din_distil_flash - 50.063160108 Psie_distil_flash <= 0
_C545: Dout_distil_flash + 50.063160108 Psie_distil_flash >= 50.063160108
_C546: Dout_distil_flash + 1000 Psie_distil_flash <= 1000
_C547: 55207.0623196 Din_distil_flash - 2763840 Psie_distil_flash
 + Ve2_distil_flash = 0
_C548: - D_distil_comp + Din_distil_comp + Dout_distil_comp = 0
_C549: Din_distil_comp - 50.063160108 Psie_distil_comp <= 0
_C55: - D_flash_furnace + D_furnace_flash = 0
_C550: Dout_distil_comp + 50.063160108 Psie_distil_comp >= 50.063160108
_C551: Dout_distil_comp + 1000 Psie_distil_comp <= 1000
_C552: 128370.641928 Din_distil_comp - 6426640 Psie_distil_comp
 + Ve2_distil_comp = 0
_C553: Ve2_distil_distil = 0
_C554: - D_distil_store + Din_distil_store + Dout_distil_store = 0
_C555: Din_distil_store - 50.063160108 Psie_distil_store <= 0
_C556: Dout_distil_store + 50.063160108 Psie_distil_store >= 50.063160108
_C557: Dout_distil_store + 1000 Psie_distil_store <= 1000
_C558: 142854.905375 Din_distil_store - 7151768 Psie_distil_store
 + Ve2_distil_store = 0
_C559: - D_distil_ctrlroom + Din_distil_ctrlroom + Dout_distil_ctrlroom = 0
_C56: L_furnace_flash - R_flash_furnace = 0
_C560: Din_distil_ctrlroom - 50.063160108 Psie_distil_ctrlroom <= 0
_C561: Dout_distil_ctrlroom + 50.063160108 Psie_distil_ctrlroom
 >= 50.063160108
_C562: Dout_distil_ctrlroom + 1000 Psie_distil_ctrlroom <= 1000
_C563: 773563.033505 Din_distil_ctrlroom - 38727010 Psie_distil_ctrlroom
 + Ve2_distil_ctrlroom = 0
_C564: - D_store_furnace + Din_store_furnace + Dout_store_furnace = 0
_C565: Din_store_furnace - 63.6396103068 Psie_store_furnace <= 0
_C566: Dout_store_furnace + 63.6396103068 Psie_store_furnace >= 63.6396103068
_C567: Dout_store_furnace + 1000 Psie_store_furnace <= 1000
_C568: 43664.9436822 Din_store_furnace - 2778820 Psie_store_furnace
 + Ve2_store_furnace = 0
_C569: - D_store_reactor + Din_store_reactor + Dout_store_reactor = 0
_C57: - L_flash_furnace + R_furnace_flash = 0
_C570: Din_store_reactor - 63.6396103068 Psie_store_reactor <= 0
_C571: Dout_store_reactor + 63.6396103068 Psie_store_reactor >= 63.6396103068
_C572: Dout_store_reactor + 1000 Psie_store_reactor <= 1000
_C573: 66682.055084 Din_store_reactor - 4243620 Psie_store_reactor
 + Ve2_store_reactor = 0
_C574: - D_store_flash + Din_store_flash + Dout_store_flash = 0
_C575: Din_store_flash - 63.6396103068 Psie_store_flash <= 0
_C576: Dout_store_flash + 63.6396103068 Psie_store_flash >= 63.6396103068
_C577: Dout_store_flash + 1000 Psie_store_flash <= 1000
_C578: 43429.5556914 Din_store_flash - 2763840 Psie_store_flash
 + Ve2_store_flash = 0
_C579: - D_store_comp + Din_store_comp + Dout_store_comp = 0
_C58: A_furnace_flash - B_flash_furnace = 0
_C580: Din_store_comp - 63.6396103068 Psie_store_comp <= 0
_C581: Dout_store_comp + 63.6396103068 Psie_store_comp >= 63.6396103068
_C582: Dout_store_comp + 1000 Psie_store_comp <= 1000
_C583: 100984.904983 Din_store_comp - 6426640 Psie_store_comp + Ve2_store_comp
 = 0
_C584: - D_store_distil + Din_store_distil + Dout_store_distil = 0
_C585: Din_store_distil - 63.6396103068 Psie_store_distil <= 0
_C586: Dout_store_distil + 63.6396103068 Psie_store_distil >= 63.6396103068
_C587: Dout_store_distil + 1000 Psie_store_distil <= 1000
_C588: 130557.681921 Din_store_distil - 8308640 Psie_store_distil
 + Ve2_store_distil = 0
_C589: Ve2_store_store = 0
_C59: - A_flash_furnace + B_furnace_flash = 0
_C590: - D_store_ctrlroom + Din_store_ctrlroom + Dout_store_ctrlroom = 0
_C591: Din_store_ctrlroom - 63.6396103068 Psie_store_ctrlroom <= 0
_C592: Dout_store_ctrlroom + 63.6396103068 Psie_store_ctrlroom
 >= 63.6396103068
_C593: Dout_store_ctrlroom + 1000 Psie_store_ctrlroom <= 1000
_C594: 608536.253024 Din_store_ctrlroom - 38727010 Psie_store_ctrlroom
 + Ve2_store_ctrlroom = 0
_C595: - 0.001 Ve2_furnace_comp - 0.001 Ve2_furnace_ctrlroom
 - 0.001 Ve2_furnace_distil - 0.001 Ve2_furnace_flash
 - 0.001 Ve2_furnace_furnace - 0.001 Ve2_furnace_reactor
 - 0.001 Ve2_furnace_store + Ve_furnace = 2778.82
_C596: Omega0_furnace - 0.75 Ve_furnace = 0
_C597: - 0.001 Ve2_reactor_comp - 0.001 Ve2_reactor_ctrlroom
 - 0.001 Ve2_reactor_distil - 0.001 Ve2_reactor_flash
 - 0.001 Ve2_reactor_furnace - 0.001 Ve2_reactor_reactor
 - 0.001 Ve2_reactor_store + Ve_reactor = 4243.62
_C598: Omega0_reactor - 0.82 Ve_reactor = 0
_C599: - 0.001 Ve2_distil_comp - 0.001 Ve2_distil_ctrlroom
 - 0.001 Ve2_distil_distil - 0.001 Ve2_distil_flash - 0.001 Ve2_distil_furnace
 - 0.001 Ve2_distil_reactor - 0.001 Ve2_distil_store + Ve_distil = 8308.64
_C6: d_reactor + l_reactor = 13.4
_C60: Wx_flash_furnace + Wx_furnace_flash = 1
_C600: Omega0_distil - 0.77 Ve_distil = 0
_C601: - 0.001 Ve2_store_comp - 0.001 Ve2_store_ctrlroom
 - 0.001 Ve2_store_distil - 0.001 Ve2_store_flash - 0.001 Ve2_store_furnace
 - 0.001 Ve2_store_reactor - 0.001 Ve2_store_store + Ve_store = 7151.768
_C602: Omega0_store - 0.82 Ve_store = 0
_C603: - Omega0_distil - Omega0_furnace - Omega0_reactor - Omega0_store
 + SumOmega = 0
_C61: Wy_flash_furnace + Wy_furnace_flash = 1
_C62: CD_furnace_flash = 0
_C63: 1000 E1_furnace_flash + 1000 E2_furnace_flash - X_flash + X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= 0
_C64: - 1000 E1_furnace_flash + 1000 E2_furnace_flash + X_flash - X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= -1000
_C65: 1000 E1_furnace_flash - 1000 E2_furnace_flash - Y_flash + Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -1000
_C66: - 1000 E1_furnace_flash - 1000 E2_furnace_flash + Y_flash - Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -2000
_C67: - L_reactor_flash + R_reactor_flash + X_flash - X_reactor = 0
_C68: A_reactor_flash - B_reactor_flash + Y_flash - Y_reactor = 0
_C69: R_reactor_flash - 1000 Wx_reactor_flash <= 0
_C7: X_reactor - 0.5 l_reactor >= 0
_C70: L_reactor_flash + 1000 Wx_reactor_flash <= 1000
_C71: A_reactor_flash - 1000 Wy_reactor_flash <= 0
_C72: B_reactor_flash + 1000 Wy_reactor_flash <= 1000
_C73: - A_reactor_flash - B_reactor_flash + D_reactor_flash - L_reactor_flash
 - R_reactor_flash = 0
_C74: - D_flash_reactor + D_reactor_flash = 0
_C75: L_reactor_flash - R_flash_reactor = 0
_C76: - L_flash_reactor + R_reactor_flash = 0
_C77: A_reactor_flash - B_flash_reactor = 0
_C78: - A_flash_reactor + B_reactor_flash = 0
_C79: Wx_flash_reactor + Wx_reactor_flash = 1
_C8: Y_reactor - 0.5 d_reactor >= 0
_C80: Wy_flash_reactor + Wy_reactor_flash = 1
_C81: CD_reactor_flash - 1231.60506417 D_reactor_flash = 0
_C82: 1000 E1_reactor_flash + 1000 E2_reactor_flash - X_flash + X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= 0
_C83: - 1000 E1_reactor_flash + 1000 E2_reactor_flash + X_flash - X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= -1000
_C84: 1000 E1_reactor_flash - 1000 E2_reactor_flash - Y_flash + Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -1000
_C85: - 1000 E1_reactor_flash - 1000 E2_reactor_flash + Y_flash - Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -2000
_C86: - L_furnace_comp + R_furnace_comp + X_comp - X_furnace = 0
_C87: A_furnace_comp - B_furnace_comp + Y_comp - Y_furnace = 0
_C88: R_furnace_comp - 1000 Wx_furnace_comp <= 0
_C89: L_furnace_comp + 1000 Wx_furnace_comp <= 1000
_C9: 0 O_flash + l_flash = 4
_C90: A_furnace_comp - 1000 Wy_furnace_comp <= 0
_C91: B_furnace_comp + 1000 Wy_furnace_comp <= 1000
_C92: - A_furnace_comp - B_furnace_comp + D_furnace_comp - L_furnace_comp
 - R_furnace_comp = 0
_C93: - D_comp_furnace + D_furnace_comp = 0
_C94: L_furnace_comp - R_comp_furnace = 0
_C95: - L_comp_furnace + R_furnace_comp = 0
_C96: A_furnace_comp - B_comp_furnace = 0
_C97: - A_comp_furnace + B_furnace_comp = 0
_C98: Wx_comp_furnace + Wx_furnace_comp = 1
_C99: Wy_comp_furnace + Wy_furnace_comp = 1
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
