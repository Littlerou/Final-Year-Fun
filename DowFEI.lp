\* Layout_Problem_Model *\
Minimize
OBJ: SumCD + SumOmega
Subject To
_C1: 9 O_furnace + l_furnace = 16
_C10: d_flash + l_flash = 8
_C100: - A_furnace_comp - B_furnace_comp + D_furnace_comp - L_furnace_comp
 - R_furnace_comp = 0
_C101: - D_comp_furnace + D_furnace_comp = 0
_C102: L_furnace_comp - R_comp_furnace = 0
_C103: - L_comp_furnace + R_furnace_comp = 0
_C104: A_furnace_comp - B_comp_furnace = 0
_C105: - A_comp_furnace + B_furnace_comp = 0
_C106: Wx_comp_furnace + Wx_furnace_comp = 1
_C107: Wy_comp_furnace + Wy_furnace_comp = 1
_C108: 1000 E1_furnace_comp + 1000 E2_furnace_comp - X_comp + X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= 0
_C109: - 1000 E1_furnace_comp + 1000 E2_furnace_comp + X_comp - X_furnace
 - 0.5 l_comp - 0.5 l_furnace >= -1000
_C11: X_flash - 0.5 l_flash >= 0
_C110: 1000 E1_furnace_comp - 1000 E2_furnace_comp - Y_comp + Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -1000
_C111: - 1000 E1_furnace_comp - 1000 E2_furnace_comp + Y_comp - Y_furnace
 - 0.5 d_comp - 0.5 d_furnace >= -2000
_C112: CD_reactor_comp = 0
_C113: - L_reactor_comp + R_reactor_comp + X_comp - X_reactor = 0
_C114: A_reactor_comp - B_reactor_comp + Y_comp - Y_reactor = 0
_C115: R_reactor_comp - 1000 Wx_reactor_comp <= 0
_C116: L_reactor_comp + 1000 Wx_reactor_comp <= 1000
_C117: A_reactor_comp - 1000 Wy_reactor_comp <= 0
_C118: B_reactor_comp + 1000 Wy_reactor_comp <= 1000
_C119: - A_reactor_comp - B_reactor_comp + D_reactor_comp - L_reactor_comp
 - R_reactor_comp = 0
_C12: Y_flash - 0.5 d_flash >= 0
_C120: - D_comp_reactor + D_reactor_comp = 0
_C121: L_reactor_comp - R_comp_reactor = 0
_C122: - L_comp_reactor + R_reactor_comp = 0
_C123: A_reactor_comp - B_comp_reactor = 0
_C124: - A_comp_reactor + B_reactor_comp = 0
_C125: Wx_comp_reactor + Wx_reactor_comp = 1
_C126: Wy_comp_reactor + Wy_reactor_comp = 1
_C127: 1000 E1_reactor_comp + 1000 E2_reactor_comp - X_comp + X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= 0
_C128: - 1000 E1_reactor_comp + 1000 E2_reactor_comp + X_comp - X_reactor
 - 0.5 l_comp - 0.5 l_reactor >= -1000
_C129: 1000 E1_reactor_comp - 1000 E2_reactor_comp - Y_comp + Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -1000
_C13: 0 O_comp + l_comp = 14
_C130: - 1000 E1_reactor_comp - 1000 E2_reactor_comp + Y_comp - Y_reactor
 - 0.5 d_comp - 0.5 d_reactor >= -2000
_C131: CD_flash_comp - 842.70764367 D_flash_comp = 0
_C132: - L_flash_comp + R_flash_comp + X_comp - X_flash = 0
_C133: A_flash_comp - B_flash_comp + Y_comp - Y_flash = 0
_C134: R_flash_comp - 1000 Wx_flash_comp <= 0
_C135: L_flash_comp + 1000 Wx_flash_comp <= 1000
_C136: A_flash_comp - 1000 Wy_flash_comp <= 0
_C137: B_flash_comp + 1000 Wy_flash_comp <= 1000
_C138: - A_flash_comp - B_flash_comp + D_flash_comp - L_flash_comp
 - R_flash_comp = 0
_C139: - D_comp_flash + D_flash_comp = 0
_C14: d_comp + l_comp = 28
_C140: L_flash_comp - R_comp_flash = 0
_C141: - L_comp_flash + R_flash_comp = 0
_C142: A_flash_comp - B_comp_flash = 0
_C143: - A_comp_flash + B_flash_comp = 0
_C144: Wx_comp_flash + Wx_flash_comp = 1
_C145: Wy_comp_flash + Wy_flash_comp = 1
_C146: 1000 E1_flash_comp + 1000 E2_flash_comp - X_comp + X_flash - 0.5 l_comp
 - 0.5 l_flash >= 0
_C147: - 1000 E1_flash_comp + 1000 E2_flash_comp + X_comp - X_flash
 - 0.5 l_comp - 0.5 l_flash >= -1000
_C148: 1000 E1_flash_comp - 1000 E2_flash_comp - Y_comp + Y_flash - 0.5 d_comp
 - 0.5 d_flash >= -1000
_C149: - 1000 E1_flash_comp - 1000 E2_flash_comp + Y_comp - Y_flash
 - 0.5 d_comp - 0.5 d_flash >= -2000
_C15: X_comp - 0.5 l_comp >= 0
_C150: CD_furnace_distil = 0
_C151: - L_furnace_distil + R_furnace_distil + X_distil - X_furnace = 0
_C152: A_furnace_distil - B_furnace_distil + Y_distil - Y_furnace = 0
_C153: R_furnace_distil - 1000 Wx_furnace_distil <= 0
_C154: L_furnace_distil + 1000 Wx_furnace_distil <= 1000
_C155: A_furnace_distil - 1000 Wy_furnace_distil <= 0
_C156: B_furnace_distil + 1000 Wy_furnace_distil <= 1000
_C157: - A_furnace_distil - B_furnace_distil + D_furnace_distil
 - L_furnace_distil - R_furnace_distil = 0
_C158: - D_distil_furnace + D_furnace_distil = 0
_C159: L_furnace_distil - R_distil_furnace = 0
_C16: Y_comp - 0.5 d_comp >= 0
_C160: - L_distil_furnace + R_furnace_distil = 0
_C161: A_furnace_distil - B_distil_furnace = 0
_C162: - A_distil_furnace + B_furnace_distil = 0
_C163: Wx_distil_furnace + Wx_furnace_distil = 1
_C164: Wy_distil_furnace + Wy_furnace_distil = 1
_C165: 1000 E1_furnace_distil + 1000 E2_furnace_distil - X_distil + X_furnace
 - 0.5 l_distil - 0.5 l_furnace >= 0
_C166: - 1000 E1_furnace_distil + 1000 E2_furnace_distil + X_distil
 - X_furnace - 0.5 l_distil - 0.5 l_furnace >= -1000
_C167: 1000 E1_furnace_distil - 1000 E2_furnace_distil - Y_distil + Y_furnace
 - 0.5 d_distil - 0.5 d_furnace >= -1000
_C168: - 1000 E1_furnace_distil - 1000 E2_furnace_distil + Y_distil
 - Y_furnace - 0.5 d_distil - 0.5 d_furnace >= -2000
_C169: CD_reactor_distil = 0
_C17: 0 O_distil + l_distil = 8
_C170: - L_reactor_distil + R_reactor_distil + X_distil - X_reactor = 0
_C171: A_reactor_distil - B_reactor_distil + Y_distil - Y_reactor = 0
_C172: R_reactor_distil - 1000 Wx_reactor_distil <= 0
_C173: L_reactor_distil + 1000 Wx_reactor_distil <= 1000
_C174: A_reactor_distil - 1000 Wy_reactor_distil <= 0
_C175: B_reactor_distil + 1000 Wy_reactor_distil <= 1000
_C176: - A_reactor_distil - B_reactor_distil + D_reactor_distil
 - L_reactor_distil - R_reactor_distil = 0
_C177: - D_distil_reactor + D_reactor_distil = 0
_C178: L_reactor_distil - R_distil_reactor = 0
_C179: - L_distil_reactor + R_reactor_distil = 0
_C18: d_distil + l_distil = 16
_C180: A_reactor_distil - B_distil_reactor = 0
_C181: - A_distil_reactor + B_reactor_distil = 0
_C182: Wx_distil_reactor + Wx_reactor_distil = 1
_C183: Wy_distil_reactor + Wy_reactor_distil = 1
_C184: 1000 E1_reactor_distil + 1000 E2_reactor_distil - X_distil + X_reactor
 - 0.5 l_distil - 0.5 l_reactor >= 0
_C185: - 1000 E1_reactor_distil + 1000 E2_reactor_distil + X_distil
 - X_reactor - 0.5 l_distil - 0.5 l_reactor >= -1000
_C186: 1000 E1_reactor_distil - 1000 E2_reactor_distil - Y_distil + Y_reactor
 - 0.5 d_distil - 0.5 d_reactor >= -1000
_C187: - 1000 E1_reactor_distil - 1000 E2_reactor_distil + Y_distil
 - Y_reactor - 0.5 d_distil - 0.5 d_reactor >= -2000
_C188: CD_flash_distil - 78.4818701125 D_flash_distil = 0
_C189: - L_flash_distil + R_flash_distil + X_distil - X_flash = 0
_C19: X_distil - 0.5 l_distil >= 0
_C190: A_flash_distil - B_flash_distil + Y_distil - Y_flash = 0
_C191: R_flash_distil - 1000 Wx_flash_distil <= 0
_C192: L_flash_distil + 1000 Wx_flash_distil <= 1000
_C193: A_flash_distil - 1000 Wy_flash_distil <= 0
_C194: B_flash_distil + 1000 Wy_flash_distil <= 1000
_C195: - A_flash_distil - B_flash_distil + D_flash_distil - L_flash_distil
 - R_flash_distil = 0
_C196: - D_distil_flash + D_flash_distil = 0
_C197: L_flash_distil - R_distil_flash = 0
_C198: - L_distil_flash + R_flash_distil = 0
_C199: A_flash_distil - B_distil_flash = 0
_C2: d_furnace + l_furnace = 23
_C20: Y_distil - 0.5 d_distil >= 0
_C200: - A_distil_flash + B_flash_distil = 0
_C201: Wx_distil_flash + Wx_flash_distil = 1
_C202: Wy_distil_flash + Wy_flash_distil = 1
_C203: 1000 E1_flash_distil + 1000 E2_flash_distil - X_distil + X_flash
 - 0.5 l_distil - 0.5 l_flash >= 0
_C204: - 1000 E1_flash_distil + 1000 E2_flash_distil + X_distil - X_flash
 - 0.5 l_distil - 0.5 l_flash >= -1000
_C205: 1000 E1_flash_distil - 1000 E2_flash_distil - Y_distil + Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -1000
_C206: - 1000 E1_flash_distil - 1000 E2_flash_distil + Y_distil - Y_flash
 - 0.5 d_distil - 0.5 d_flash >= -2000
_C207: CD_comp_distil = 0
_C208: - L_comp_distil + R_comp_distil - X_comp + X_distil = 0
_C209: A_comp_distil - B_comp_distil - Y_comp + Y_distil = 0
_C21: - 10 O_store + l_store = 12
_C210: R_comp_distil - 1000 Wx_comp_distil <= 0
_C211: L_comp_distil + 1000 Wx_comp_distil <= 1000
_C212: A_comp_distil - 1000 Wy_comp_distil <= 0
_C213: B_comp_distil + 1000 Wy_comp_distil <= 1000
_C214: - A_comp_distil - B_comp_distil + D_comp_distil - L_comp_distil
 - R_comp_distil = 0
_C215: D_comp_distil - D_distil_comp = 0
_C216: L_comp_distil - R_distil_comp = 0
_C217: - L_distil_comp + R_comp_distil = 0
_C218: A_comp_distil - B_distil_comp = 0
_C219: - A_distil_comp + B_comp_distil = 0
_C22: d_store + l_store = 34
_C220: Wx_comp_distil + Wx_distil_comp = 1
_C221: Wy_comp_distil + Wy_distil_comp = 1
_C222: 1000 E1_comp_distil + 1000 E2_comp_distil + X_comp - X_distil
 - 0.5 l_comp - 0.5 l_distil >= 0
_C223: - 1000 E1_comp_distil + 1000 E2_comp_distil - X_comp + X_distil
 - 0.5 l_comp - 0.5 l_distil >= -1000
_C224: 1000 E1_comp_distil - 1000 E2_comp_distil + Y_comp - Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -1000
_C225: - 1000 E1_comp_distil - 1000 E2_comp_distil - Y_comp + Y_distil
 - 0.5 d_comp - 0.5 d_distil >= -2000
_C226: CD_furnace_store - 75.8045442191 D_furnace_store = 0
_C227: - L_furnace_store + R_furnace_store - X_furnace + X_store = 0
_C228: A_furnace_store - B_furnace_store - Y_furnace + Y_store = 0
_C229: R_furnace_store - 1000 Wx_furnace_store <= 0
_C23: X_store - 0.5 l_store >= 0
_C230: L_furnace_store + 1000 Wx_furnace_store <= 1000
_C231: A_furnace_store - 1000 Wy_furnace_store <= 0
_C232: B_furnace_store + 1000 Wy_furnace_store <= 1000
_C233: - A_furnace_store - B_furnace_store + D_furnace_store - L_furnace_store
 - R_furnace_store = 0
_C234: D_furnace_store - D_store_furnace = 0
_C235: L_furnace_store - R_store_furnace = 0
_C236: - L_store_furnace + R_furnace_store = 0
_C237: A_furnace_store - B_store_furnace = 0
_C238: - A_store_furnace + B_furnace_store = 0
_C239: Wx_furnace_store + Wx_store_furnace = 1
_C24: Y_store - 0.5 d_store >= 0
_C240: Wy_furnace_store + Wy_store_furnace = 1
_C241: 1000 E1_furnace_store + 1000 E2_furnace_store + X_furnace - X_store
 - 0.5 l_furnace - 0.5 l_store >= 0
_C242: - 1000 E1_furnace_store + 1000 E2_furnace_store - X_furnace + X_store
 - 0.5 l_furnace - 0.5 l_store >= -1000
_C243: 1000 E1_furnace_store - 1000 E2_furnace_store + Y_furnace - Y_store
 - 0.5 d_furnace - 0.5 d_store >= -1000
_C244: - 1000 E1_furnace_store - 1000 E2_furnace_store - Y_furnace + Y_store
 - 0.5 d_furnace - 0.5 d_store >= -2000
_C245: CD_reactor_store = 0
_C246: - L_reactor_store + R_reactor_store - X_reactor + X_store = 0
_C247: A_reactor_store - B_reactor_store - Y_reactor + Y_store = 0
_C248: R_reactor_store - 1000 Wx_reactor_store <= 0
_C249: L_reactor_store + 1000 Wx_reactor_store <= 1000
_C25: 0 O_ctrlroom + l_ctrlroom = 15
_C250: A_reactor_store - 1000 Wy_reactor_store <= 0
_C251: B_reactor_store + 1000 Wy_reactor_store <= 1000
_C252: - A_reactor_store - B_reactor_store + D_reactor_store - L_reactor_store
 - R_reactor_store = 0
_C253: D_reactor_store - D_store_reactor = 0
_C254: L_reactor_store - R_store_reactor = 0
_C255: - L_store_reactor + R_reactor_store = 0
_C256: A_reactor_store - B_store_reactor = 0
_C257: - A_store_reactor + B_reactor_store = 0
_C258: Wx_reactor_store + Wx_store_reactor = 1
_C259: Wy_reactor_store + Wy_store_reactor = 1
_C26: d_ctrlroom + l_ctrlroom = 30
_C260: 1000 E1_reactor_store + 1000 E2_reactor_store + X_reactor - X_store
 - 0.5 l_reactor - 0.5 l_store >= 0
_C261: - 1000 E1_reactor_store + 1000 E2_reactor_store - X_reactor + X_store
 - 0.5 l_reactor - 0.5 l_store >= -1000
_C262: 1000 E1_reactor_store - 1000 E2_reactor_store + Y_reactor - Y_store
 - 0.5 d_reactor - 0.5 d_store >= -1000
_C263: - 1000 E1_reactor_store - 1000 E2_reactor_store - Y_reactor + Y_store
 - 0.5 d_reactor - 0.5 d_store >= -2000
_C264: CD_flash_store = 0
_C265: - L_flash_store + R_flash_store - X_flash + X_store = 0
_C266: A_flash_store - B_flash_store - Y_flash + Y_store = 0
_C267: R_flash_store - 1000 Wx_flash_store <= 0
_C268: L_flash_store + 1000 Wx_flash_store <= 1000
_C269: A_flash_store - 1000 Wy_flash_store <= 0
_C27: X_ctrlroom - 0.5 l_ctrlroom >= 0
_C270: B_flash_store + 1000 Wy_flash_store <= 1000
_C271: - A_flash_store - B_flash_store + D_flash_store - L_flash_store
 - R_flash_store = 0
_C272: D_flash_store - D_store_flash = 0
_C273: L_flash_store - R_store_flash = 0
_C274: - L_store_flash + R_flash_store = 0
_C275: A_flash_store - B_store_flash = 0
_C276: - A_store_flash + B_flash_store = 0
_C277: Wx_flash_store + Wx_store_flash = 1
_C278: Wy_flash_store + Wy_store_flash = 1
_C279: 1000 E1_flash_store + 1000 E2_flash_store + X_flash - X_store
 - 0.5 l_flash - 0.5 l_store >= 0
_C28: Y_ctrlroom - 0.5 d_ctrlroom >= 0
_C280: - 1000 E1_flash_store + 1000 E2_flash_store - X_flash + X_store
 - 0.5 l_flash - 0.5 l_store >= -1000
_C281: 1000 E1_flash_store - 1000 E2_flash_store + Y_flash - Y_store
 - 0.5 d_flash - 0.5 d_store >= -1000
_C282: - 1000 E1_flash_store - 1000 E2_flash_store - Y_flash + Y_store
 - 0.5 d_flash - 0.5 d_store >= -2000
_C283: CD_comp_store = 0
_C284: - L_comp_store + R_comp_store - X_comp + X_store = 0
_C285: A_comp_store - B_comp_store - Y_comp + Y_store = 0
_C286: R_comp_store - 1000 Wx_comp_store <= 0
_C287: L_comp_store + 1000 Wx_comp_store <= 1000
_C288: A_comp_store - 1000 Wy_comp_store <= 0
_C289: B_comp_store + 1000 Wy_comp_store <= 1000
_C29: - D_ctrlroom_furnace + D_road_ctrlroom_furnace = 0
_C290: - A_comp_store - B_comp_store + D_comp_store - L_comp_store
 - R_comp_store = 0
_C291: D_comp_store - D_store_comp = 0
_C292: L_comp_store - R_store_comp = 0
_C293: - L_store_comp + R_comp_store = 0
_C294: A_comp_store - B_store_comp = 0
_C295: - A_store_comp + B_comp_store = 0
_C296: Wx_comp_store + Wx_store_comp = 1
_C297: Wy_comp_store + Wy_store_comp = 1
_C298: 1000 E1_comp_store + 1000 E2_comp_store + X_comp - X_store - 0.5 l_comp
 - 0.5 l_store >= 0
_C299: - 1000 E1_comp_store + 1000 E2_comp_store - X_comp + X_store
 - 0.5 l_comp - 0.5 l_store >= -1000
_C3: X_furnace - 0.5 l_furnace >= 0
_C30: - D_ctrlroom_reactor + D_road_ctrlroom_reactor = 0
_C300: 1000 E1_comp_store - 1000 E2_comp_store + Y_comp - Y_store - 0.5 d_comp
 - 0.5 d_store >= -1000
_C301: - 1000 E1_comp_store - 1000 E2_comp_store - Y_comp + Y_store
 - 0.5 d_comp - 0.5 d_store >= -2000
_C302: CD_distil_store - 451.755491306 D_distil_store = 0
_C303: - L_distil_store + R_distil_store - X_distil + X_store = 0
_C304: A_distil_store - B_distil_store - Y_distil + Y_store = 0
_C305: R_distil_store - 1000 Wx_distil_store <= 0
_C306: L_distil_store + 1000 Wx_distil_store <= 1000
_C307: A_distil_store - 1000 Wy_distil_store <= 0
_C308: B_distil_store + 1000 Wy_distil_store <= 1000
_C309: - A_distil_store - B_distil_store + D_distil_store - L_distil_store
 - R_distil_store = 0
_C31: - D_ctrlroom_flash + D_road_ctrlroom_flash = 0
_C310: D_distil_store - D_store_distil = 0
_C311: L_distil_store - R_store_distil = 0
_C312: - L_store_distil + R_distil_store = 0
_C313: A_distil_store - B_store_distil = 0
_C314: - A_store_distil + B_distil_store = 0
_C315: Wx_distil_store + Wx_store_distil = 1
_C316: Wy_distil_store + Wy_store_distil = 1
_C317: 1000 E1_distil_store + 1000 E2_distil_store + X_distil - X_store
 - 0.5 l_distil - 0.5 l_store >= 0
_C318: - 1000 E1_distil_store + 1000 E2_distil_store - X_distil + X_store
 - 0.5 l_distil - 0.5 l_store >= -1000
_C319: 1000 E1_distil_store - 1000 E2_distil_store + Y_distil - Y_store
 - 0.5 d_distil - 0.5 d_store >= -1000
_C32: - D_ctrlroom_comp + D_road_ctrlroom_comp = 0
_C320: - 1000 E1_distil_store - 1000 E2_distil_store - Y_distil + Y_store
 - 0.5 d_distil - 0.5 d_store >= -2000
_C321: CD_furnace_ctrlroom = 0
_C322: - L_furnace_ctrlroom + R_furnace_ctrlroom + X_ctrlroom - X_furnace = 0
_C323: A_furnace_ctrlroom - B_furnace_ctrlroom + Y_ctrlroom - Y_furnace = 0
_C324: R_furnace_ctrlroom - 1000 Wx_furnace_ctrlroom <= 0
_C325: L_furnace_ctrlroom + 1000 Wx_furnace_ctrlroom <= 1000
_C326: A_furnace_ctrlroom - 1000 Wy_furnace_ctrlroom <= 0
_C327: B_furnace_ctrlroom + 1000 Wy_furnace_ctrlroom <= 1000
_C328: - A_furnace_ctrlroom - B_furnace_ctrlroom + D_furnace_ctrlroom
 - L_furnace_ctrlroom - R_furnace_ctrlroom = 0
_C329: - D_ctrlroom_furnace + D_furnace_ctrlroom = 0
_C33: - D_ctrlroom_distil + D_road_ctrlroom_distil = 0
_C330: L_furnace_ctrlroom - R_ctrlroom_furnace = 0
_C331: - L_ctrlroom_furnace + R_furnace_ctrlroom = 0
_C332: A_furnace_ctrlroom - B_ctrlroom_furnace = 0
_C333: - A_ctrlroom_furnace + B_furnace_ctrlroom = 0
_C334: Wx_ctrlroom_furnace + Wx_furnace_ctrlroom = 1
_C335: Wy_ctrlroom_furnace + Wy_furnace_ctrlroom = 1
_C336: 1000 E1_furnace_ctrlroom + 1000 E2_furnace_ctrlroom - X_ctrlroom
 + X_furnace - 0.5 l_ctrlroom - 0.5 l_furnace >= 0
_C337: - 1000 E1_furnace_ctrlroom + 1000 E2_furnace_ctrlroom + X_ctrlroom
 - X_furnace - 0.5 l_ctrlroom - 0.5 l_furnace >= -1000
_C338: 1000 E1_furnace_ctrlroom - 1000 E2_furnace_ctrlroom - Y_ctrlroom
 + Y_furnace - 0.5 d_ctrlroom - 0.5 d_furnace >= -1000
_C339: - 1000 E1_furnace_ctrlroom - 1000 E2_furnace_ctrlroom + Y_ctrlroom
 - Y_furnace - 0.5 d_ctrlroom - 0.5 d_furnace >= -2000
_C34: - D_ctrlroom_store + D_road_ctrlroom_store = 0
_C340: CD_reactor_ctrlroom = 0
_C341: - L_reactor_ctrlroom + R_reactor_ctrlroom + X_ctrlroom - X_reactor = 0
_C342: A_reactor_ctrlroom - B_reactor_ctrlroom + Y_ctrlroom - Y_reactor = 0
_C343: R_reactor_ctrlroom - 1000 Wx_reactor_ctrlroom <= 0
_C344: L_reactor_ctrlroom + 1000 Wx_reactor_ctrlroom <= 1000
_C345: A_reactor_ctrlroom - 1000 Wy_reactor_ctrlroom <= 0
_C346: B_reactor_ctrlroom + 1000 Wy_reactor_ctrlroom <= 1000
_C347: - A_reactor_ctrlroom - B_reactor_ctrlroom + D_reactor_ctrlroom
 - L_reactor_ctrlroom - R_reactor_ctrlroom = 0
_C348: - D_ctrlroom_reactor + D_reactor_ctrlroom = 0
_C349: L_reactor_ctrlroom - R_ctrlroom_reactor = 0
_C35: - 120 D_road_ctrlroom_comp - 120 D_road_ctrlroom_ctrlroom
 - 120 D_road_ctrlroom_distil - 120 D_road_ctrlroom_flash
 - 120 D_road_ctrlroom_furnace - 120 D_road_ctrlroom_reactor
 - 120 D_road_ctrlroom_store + roadCost = 0
_C350: - L_ctrlroom_reactor + R_reactor_ctrlroom = 0
_C351: A_reactor_ctrlroom - B_ctrlroom_reactor = 0
_C352: - A_ctrlroom_reactor + B_reactor_ctrlroom = 0
_C353: Wx_ctrlroom_reactor + Wx_reactor_ctrlroom = 1
_C354: Wy_ctrlroom_reactor + Wy_reactor_ctrlroom = 1
_C355: 1000 E1_reactor_ctrlroom + 1000 E2_reactor_ctrlroom - X_ctrlroom
 + X_reactor - 0.5 l_ctrlroom - 0.5 l_reactor >= 0
_C356: - 1000 E1_reactor_ctrlroom + 1000 E2_reactor_ctrlroom + X_ctrlroom
 - X_reactor - 0.5 l_ctrlroom - 0.5 l_reactor >= -1000
_C357: 1000 E1_reactor_ctrlroom - 1000 E2_reactor_ctrlroom - Y_ctrlroom
 + Y_reactor - 0.5 d_ctrlroom - 0.5 d_reactor >= -1000
_C358: - 1000 E1_reactor_ctrlroom - 1000 E2_reactor_ctrlroom + Y_ctrlroom
 - Y_reactor - 0.5 d_ctrlroom - 0.5 d_reactor >= -2000
_C359: CD_flash_ctrlroom = 0
_C36: CD_furnace_reactor - 1075.48971174 D_furnace_reactor = 0
_C360: - L_flash_ctrlroom + R_flash_ctrlroom + X_ctrlroom - X_flash = 0
_C361: A_flash_ctrlroom - B_flash_ctrlroom + Y_ctrlroom - Y_flash = 0
_C362: R_flash_ctrlroom - 1000 Wx_flash_ctrlroom <= 0
_C363: L_flash_ctrlroom + 1000 Wx_flash_ctrlroom <= 1000
_C364: A_flash_ctrlroom - 1000 Wy_flash_ctrlroom <= 0
_C365: B_flash_ctrlroom + 1000 Wy_flash_ctrlroom <= 1000
_C366: - A_flash_ctrlroom - B_flash_ctrlroom + D_flash_ctrlroom
 - L_flash_ctrlroom - R_flash_ctrlroom = 0
_C367: - D_ctrlroom_flash + D_flash_ctrlroom = 0
_C368: L_flash_ctrlroom - R_ctrlroom_flash = 0
_C369: - L_ctrlroom_flash + R_flash_ctrlroom = 0
_C37: - L_furnace_reactor + R_furnace_reactor - X_furnace + X_reactor = 0
_C370: A_flash_ctrlroom - B_ctrlroom_flash = 0
_C371: - A_ctrlroom_flash + B_flash_ctrlroom = 0
_C372: Wx_ctrlroom_flash + Wx_flash_ctrlroom = 1
_C373: Wy_ctrlroom_flash + Wy_flash_ctrlroom = 1
_C374: 1000 E1_flash_ctrlroom + 1000 E2_flash_ctrlroom - X_ctrlroom + X_flash
 - 0.5 l_ctrlroom - 0.5 l_flash >= 0
_C375: - 1000 E1_flash_ctrlroom + 1000 E2_flash_ctrlroom + X_ctrlroom
 - X_flash - 0.5 l_ctrlroom - 0.5 l_flash >= -1000
_C376: 1000 E1_flash_ctrlroom - 1000 E2_flash_ctrlroom - Y_ctrlroom + Y_flash
 - 0.5 d_ctrlroom - 0.5 d_flash >= -1000
_C377: - 1000 E1_flash_ctrlroom - 1000 E2_flash_ctrlroom + Y_ctrlroom
 - Y_flash - 0.5 d_ctrlroom - 0.5 d_flash >= -2000
_C378: CD_comp_ctrlroom = 0
_C379: - L_comp_ctrlroom + R_comp_ctrlroom - X_comp + X_ctrlroom = 0
_C38: A_furnace_reactor - B_furnace_reactor - Y_furnace + Y_reactor = 0
_C380: A_comp_ctrlroom - B_comp_ctrlroom - Y_comp + Y_ctrlroom = 0
_C381: R_comp_ctrlroom - 1000 Wx_comp_ctrlroom <= 0
_C382: L_comp_ctrlroom + 1000 Wx_comp_ctrlroom <= 1000
_C383: A_comp_ctrlroom - 1000 Wy_comp_ctrlroom <= 0
_C384: B_comp_ctrlroom + 1000 Wy_comp_ctrlroom <= 1000
_C385: - A_comp_ctrlroom - B_comp_ctrlroom + D_comp_ctrlroom - L_comp_ctrlroom
 - R_comp_ctrlroom = 0
_C386: D_comp_ctrlroom - D_ctrlroom_comp = 0
_C387: L_comp_ctrlroom - R_ctrlroom_comp = 0
_C388: - L_ctrlroom_comp + R_comp_ctrlroom = 0
_C389: A_comp_ctrlroom - B_ctrlroom_comp = 0
_C39: R_furnace_reactor - 1000 Wx_furnace_reactor <= 0
_C390: - A_ctrlroom_comp + B_comp_ctrlroom = 0
_C391: Wx_comp_ctrlroom + Wx_ctrlroom_comp = 1
_C392: Wy_comp_ctrlroom + Wy_ctrlroom_comp = 1
_C393: 1000 E1_comp_ctrlroom + 1000 E2_comp_ctrlroom + X_comp - X_ctrlroom
 - 0.5 l_comp - 0.5 l_ctrlroom >= 0
_C394: - 1000 E1_comp_ctrlroom + 1000 E2_comp_ctrlroom - X_comp + X_ctrlroom
 - 0.5 l_comp - 0.5 l_ctrlroom >= -1000
_C395: 1000 E1_comp_ctrlroom - 1000 E2_comp_ctrlroom + Y_comp - Y_ctrlroom
 - 0.5 d_comp - 0.5 d_ctrlroom >= -1000
_C396: - 1000 E1_comp_ctrlroom - 1000 E2_comp_ctrlroom - Y_comp + Y_ctrlroom
 - 0.5 d_comp - 0.5 d_ctrlroom >= -2000
_C397: CD_distil_ctrlroom = 0
_C398: - L_distil_ctrlroom + R_distil_ctrlroom + X_ctrlroom - X_distil = 0
_C399: A_distil_ctrlroom - B_distil_ctrlroom + Y_ctrlroom - Y_distil = 0
_C4: Y_furnace - 0.5 d_furnace >= 0
_C40: L_furnace_reactor + 1000 Wx_furnace_reactor <= 1000
_C400: R_distil_ctrlroom - 1000 Wx_distil_ctrlroom <= 0
_C401: L_distil_ctrlroom + 1000 Wx_distil_ctrlroom <= 1000
_C402: A_distil_ctrlroom - 1000 Wy_distil_ctrlroom <= 0
_C403: B_distil_ctrlroom + 1000 Wy_distil_ctrlroom <= 1000
_C404: - A_distil_ctrlroom - B_distil_ctrlroom + D_distil_ctrlroom
 - L_distil_ctrlroom - R_distil_ctrlroom = 0
_C405: - D_ctrlroom_distil + D_distil_ctrlroom = 0
_C406: L_distil_ctrlroom - R_ctrlroom_distil = 0
_C407: - L_ctrlroom_distil + R_distil_ctrlroom = 0
_C408: A_distil_ctrlroom - B_ctrlroom_distil = 0
_C409: - A_ctrlroom_distil + B_distil_ctrlroom = 0
_C41: A_furnace_reactor - 1000 Wy_furnace_reactor <= 0
_C410: Wx_ctrlroom_distil + Wx_distil_ctrlroom = 1
_C411: Wy_ctrlroom_distil + Wy_distil_ctrlroom = 1
_C412: 1000 E1_distil_ctrlroom + 1000 E2_distil_ctrlroom - X_ctrlroom
 + X_distil - 0.5 l_ctrlroom - 0.5 l_distil >= 0
_C413: - 1000 E1_distil_ctrlroom + 1000 E2_distil_ctrlroom + X_ctrlroom
 - X_distil - 0.5 l_ctrlroom - 0.5 l_distil >= -1000
_C414: 1000 E1_distil_ctrlroom - 1000 E2_distil_ctrlroom - Y_ctrlroom
 + Y_distil - 0.5 d_ctrlroom - 0.5 d_distil >= -1000
_C415: - 1000 E1_distil_ctrlroom - 1000 E2_distil_ctrlroom + Y_ctrlroom
 - Y_distil - 0.5 d_ctrlroom - 0.5 d_distil >= -2000
_C416: CD_store_ctrlroom = 0
_C417: - L_store_ctrlroom + R_store_ctrlroom + X_ctrlroom - X_store = 0
_C418: A_store_ctrlroom - B_store_ctrlroom + Y_ctrlroom - Y_store = 0
_C419: R_store_ctrlroom - 1000 Wx_store_ctrlroom <= 0
_C42: B_furnace_reactor + 1000 Wy_furnace_reactor <= 1000
_C420: L_store_ctrlroom + 1000 Wx_store_ctrlroom <= 1000
_C421: A_store_ctrlroom - 1000 Wy_store_ctrlroom <= 0
_C422: B_store_ctrlroom + 1000 Wy_store_ctrlroom <= 1000
_C423: - A_store_ctrlroom - B_store_ctrlroom + D_store_ctrlroom
 - L_store_ctrlroom - R_store_ctrlroom = 0
_C424: - D_ctrlroom_store + D_store_ctrlroom = 0
_C425: L_store_ctrlroom - R_ctrlroom_store = 0
_C426: - L_ctrlroom_store + R_store_ctrlroom = 0
_C427: A_store_ctrlroom - B_ctrlroom_store = 0
_C428: - A_ctrlroom_store + B_store_ctrlroom = 0
_C429: Wx_ctrlroom_store + Wx_store_ctrlroom = 1
_C43: - A_furnace_reactor - B_furnace_reactor + D_furnace_reactor
 - L_furnace_reactor - R_furnace_reactor = 0
_C430: Wy_ctrlroom_store + Wy_store_ctrlroom = 1
_C431: 1000 E1_store_ctrlroom + 1000 E2_store_ctrlroom - X_ctrlroom + X_store
 - 0.5 l_ctrlroom - 0.5 l_store >= 0
_C432: - 1000 E1_store_ctrlroom + 1000 E2_store_ctrlroom + X_ctrlroom
 - X_store - 0.5 l_ctrlroom - 0.5 l_store >= -1000
_C433: 1000 E1_store_ctrlroom - 1000 E2_store_ctrlroom - Y_ctrlroom + Y_store
 - 0.5 d_ctrlroom - 0.5 d_store >= -1000
_C434: - 1000 E1_store_ctrlroom - 1000 E2_store_ctrlroom + Y_ctrlroom
 - Y_store - 0.5 d_ctrlroom - 0.5 d_store >= -2000
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
_C436: X_furnace + 0.5 l_furnace <= 150
_C437: Y_furnace + 0.5 d_furnace <= 150
_C438: X_reactor + 0.5 l_reactor <= 150
_C439: Y_reactor + 0.5 d_reactor <= 150
_C44: D_furnace_reactor - D_reactor_furnace = 0
_C440: X_flash + 0.5 l_flash <= 150
_C441: Y_flash + 0.5 d_flash <= 150
_C442: X_comp + 0.5 l_comp <= 150
_C443: Y_comp + 0.5 d_comp <= 150
_C444: X_distil + 0.5 l_distil <= 150
_C445: Y_distil + 0.5 d_distil <= 150
_C446: X_store + 0.5 l_store <= 150
_C447: Y_store + 0.5 d_store <= 150
_C448: X_ctrlroom + 0.5 l_ctrlroom <= 150
_C449: Y_ctrlroom + 0.5 d_ctrlroom <= 150
_C45: L_furnace_reactor - R_reactor_furnace = 0
_C450: Ve2_furnace_furnace = 0
_C451: - D_furnace_reactor + Din_furnace_reactor + Dout_furnace_reactor = 0
_C452: Din_furnace_reactor - 40.7293505963 Psie_furnace_reactor <= 0
_C453: Dout_furnace_reactor + 40.7293505963 Psie_furnace_reactor
 >= 40.7293505963
_C454: Dout_furnace_reactor + 1000 Psie_furnace_reactor <= 1000
_C455: 104190.711069 Din_furnace_reactor - 4243620 Psie_furnace_reactor
 + Ve2_furnace_reactor = 0
_C456: - D_furnace_flash + Din_furnace_flash + Dout_furnace_flash = 0
_C457: Din_furnace_flash - 40.7293505963 Psie_furnace_flash <= 0
_C458: Dout_furnace_flash + 40.7293505963 Psie_furnace_flash >= 40.7293505963
_C459: Dout_furnace_flash + 1000 Psie_furnace_flash <= 1000
_C46: - L_reactor_furnace + R_furnace_reactor = 0
_C460: 67858.6807679 Din_furnace_flash - 2763840 Psie_furnace_flash
 + Ve2_furnace_flash = 0
_C461: - D_furnace_comp + Din_furnace_comp + Dout_furnace_comp = 0
_C462: Din_furnace_comp - 40.7293505963 Psie_furnace_comp <= 0
_C463: Dout_furnace_comp + 40.7293505963 Psie_furnace_comp >= 40.7293505963
_C464: Dout_furnace_comp + 1000 Psie_furnace_comp <= 1000
_C465: 157788.914036 Din_furnace_comp - 6426640 Psie_furnace_comp
 + Ve2_furnace_comp = 0
_C466: - D_furnace_distil + Din_furnace_distil + Dout_furnace_distil = 0
_C467: Din_furnace_distil - 40.7293505963 Psie_furnace_distil <= 0
_C468: Dout_furnace_distil + 40.7293505963 Psie_furnace_distil
 >= 40.7293505963
_C469: Dout_furnace_distil + 1000 Psie_furnace_distil <= 1000
_C47: A_furnace_reactor - B_reactor_furnace = 0
_C470: 203996.378001 Din_furnace_distil - 8308640 Psie_furnace_distil
 + Ve2_furnace_distil = 0
_C471: - D_furnace_store + Din_furnace_store + Dout_furnace_store = 0
_C472: Din_furnace_store - 40.7293505963 Psie_furnace_store <= 0
_C473: Dout_furnace_store + 40.7293505963 Psie_furnace_store >= 40.7293505963
_C474: Dout_furnace_store + 1000 Psie_furnace_store <= 1000
_C475: 175592.487857 Din_furnace_store - 7151768 Psie_furnace_store
 + Ve2_furnace_store = 0
_C476: - D_furnace_ctrlroom + Din_furnace_ctrlroom + Dout_furnace_ctrlroom = 0
_C477: Din_furnace_ctrlroom - 40.7293505963 Psie_furnace_ctrlroom <= 0
_C478: Dout_furnace_ctrlroom + 40.7293505963 Psie_furnace_ctrlroom
 >= 40.7293505963
_C479: Dout_furnace_ctrlroom + 1000 Psie_furnace_ctrlroom <= 1000
_C48: - A_reactor_furnace + B_furnace_reactor = 0
_C480: 950837.89535 Din_furnace_ctrlroom - 38727010 Psie_furnace_ctrlroom
 + Ve2_furnace_ctrlroom = 0
_C481: - D_reactor_furnace + Din_reactor_furnace + Dout_reactor_furnace = 0
_C482: Din_reactor_furnace - 57.1907964624 Psie_reactor_furnace <= 0
_C483: Dout_reactor_furnace + 57.1907964624 Psie_reactor_furnace
 >= 57.1907964624
_C484: Dout_reactor_furnace + 1000 Psie_reactor_furnace <= 1000
_C485: 48588.5871834 Din_reactor_furnace - 2778820 Psie_reactor_furnace
 + Ve2_reactor_furnace = 0
_C486: Ve2_reactor_reactor = 0
_C487: - D_reactor_flash + Din_reactor_flash + Dout_reactor_flash = 0
_C488: Din_reactor_flash - 57.1907964624 Psie_reactor_flash <= 0
_C489: Dout_reactor_flash + 57.1907964624 Psie_reactor_flash >= 57.1907964624
_C49: Wx_furnace_reactor + Wx_reactor_furnace = 1
_C490: Dout_reactor_flash + 1000 Psie_reactor_flash <= 1000
_C491: 48326.6569267 Din_reactor_flash - 2763840 Psie_reactor_flash
 + Ve2_reactor_flash = 0
_C492: - D_reactor_comp + Din_reactor_comp + Dout_reactor_comp = 0
_C493: Din_reactor_comp - 57.1907964624 Psie_reactor_comp <= 0
_C494: Dout_reactor_comp + 57.1907964624 Psie_reactor_comp >= 57.1907964624
_C495: Dout_reactor_comp + 1000 Psie_reactor_comp <= 1000
_C496: 112371.92691 Din_reactor_comp - 6426640 Psie_reactor_comp
 + Ve2_reactor_comp = 0
_C497: - D_reactor_distil + Din_reactor_distil + Dout_reactor_distil = 0
_C498: Din_reactor_distil - 57.1907964624 Psie_reactor_distil <= 0
_C499: Dout_reactor_distil + 57.1907964624 Psie_reactor_distil
 >= 57.1907964624
_C5: - 2.8 O_reactor + l_reactor = 5.8
_C50: Wy_furnace_reactor + Wy_reactor_furnace = 1
_C500: Dout_reactor_distil + 1000 Psie_reactor_distil <= 1000
_C501: 145279.319645 Din_reactor_distil - 8308640 Psie_reactor_distil
 + Ve2_reactor_distil = 0
_C502: - D_reactor_store + Din_reactor_store + Dout_reactor_store = 0
_C503: Din_reactor_store - 57.1907964624 Psie_reactor_store <= 0
_C504: Dout_reactor_store + 57.1907964624 Psie_reactor_store >= 57.1907964624
_C505: Dout_reactor_store + 1000 Psie_reactor_store <= 1000
_C506: 125051.029928 Din_reactor_store - 7151768 Psie_reactor_store
 + Ve2_reactor_store = 0
_C507: - D_reactor_ctrlroom + Din_reactor_ctrlroom + Dout_reactor_ctrlroom = 0
_C508: Din_reactor_ctrlroom - 57.1907964624 Psie_reactor_ctrlroom <= 0
_C509: Dout_reactor_ctrlroom + 57.1907964624 Psie_reactor_ctrlroom
 >= 57.1907964624
_C51: 1000 E1_furnace_reactor + 1000 E2_furnace_reactor + X_furnace
 - X_reactor - 0.5 l_furnace - 0.5 l_reactor >= 0
_C510: Dout_reactor_ctrlroom + 1000 Psie_reactor_ctrlroom <= 1000
_C511: 677154.584226 Din_reactor_ctrlroom - 38727010 Psie_reactor_ctrlroom
 + Ve2_reactor_ctrlroom = 0
_C512: - D_distil_furnace + Din_distil_furnace + Dout_distil_furnace = 0
_C513: Din_distil_furnace - 50.063160108 Psie_distil_furnace <= 0
_C514: Dout_distil_furnace + 50.063160108 Psie_distil_furnace >= 50.063160108
_C515: Dout_distil_furnace + 1000 Psie_distil_furnace <= 1000
_C516: 55506.2843417 Din_distil_furnace - 2778820 Psie_distil_furnace
 + Ve2_distil_furnace = 0
_C517: - D_distil_reactor + Din_distil_reactor + Dout_distil_reactor = 0
_C518: Din_distil_reactor - 50.063160108 Psie_distil_reactor <= 0
_C519: Dout_distil_reactor + 50.063160108 Psie_distil_reactor >= 50.063160108
_C52: - 1000 E1_furnace_reactor + 1000 E2_furnace_reactor - X_furnace
 + X_reactor - 0.5 l_furnace - 0.5 l_reactor >= -1000
_C520: Dout_distil_reactor + 1000 Psie_distil_reactor <= 1000
_C521: 84765.3242593 Din_distil_reactor - 4243620 Psie_distil_reactor
 + Ve2_distil_reactor = 0
_C522: - D_distil_flash + Din_distil_flash + Dout_distil_flash = 0
_C523: Din_distil_flash - 50.063160108 Psie_distil_flash <= 0
_C524: Dout_distil_flash + 50.063160108 Psie_distil_flash >= 50.063160108
_C525: Dout_distil_flash + 1000 Psie_distil_flash <= 1000
_C526: 55207.0623196 Din_distil_flash - 2763840 Psie_distil_flash
 + Ve2_distil_flash = 0
_C527: - D_distil_comp + Din_distil_comp + Dout_distil_comp = 0
_C528: Din_distil_comp - 50.063160108 Psie_distil_comp <= 0
_C529: Dout_distil_comp + 50.063160108 Psie_distil_comp >= 50.063160108
_C53: 1000 E1_furnace_reactor - 1000 E2_furnace_reactor + Y_furnace
 - Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -1000
_C530: Dout_distil_comp + 1000 Psie_distil_comp <= 1000
_C531: 128370.641928 Din_distil_comp - 6426640 Psie_distil_comp
 + Ve2_distil_comp = 0
_C532: Ve2_distil_distil = 0
_C533: - D_distil_store + Din_distil_store + Dout_distil_store = 0
_C534: Din_distil_store - 50.063160108 Psie_distil_store <= 0
_C535: Dout_distil_store + 50.063160108 Psie_distil_store >= 50.063160108
_C536: Dout_distil_store + 1000 Psie_distil_store <= 1000
_C537: 142854.905375 Din_distil_store - 7151768 Psie_distil_store
 + Ve2_distil_store = 0
_C538: - D_distil_ctrlroom + Din_distil_ctrlroom + Dout_distil_ctrlroom = 0
_C539: Din_distil_ctrlroom - 50.063160108 Psie_distil_ctrlroom <= 0
_C54: - 1000 E1_furnace_reactor - 1000 E2_furnace_reactor - Y_furnace
 + Y_reactor - 0.5 d_furnace - 0.5 d_reactor >= -2000
_C540: Dout_distil_ctrlroom + 50.063160108 Psie_distil_ctrlroom
 >= 50.063160108
_C541: Dout_distil_ctrlroom + 1000 Psie_distil_ctrlroom <= 1000
_C542: 773563.033505 Din_distil_ctrlroom - 38727010 Psie_distil_ctrlroom
 + Ve2_distil_ctrlroom = 0
_C543: - D_store_furnace + Din_store_furnace + Dout_store_furnace = 0
_C544: Din_store_furnace - 63.6396103068 Psie_store_furnace <= 0
_C545: Dout_store_furnace + 63.6396103068 Psie_store_furnace >= 63.6396103068
_C546: Dout_store_furnace + 1000 Psie_store_furnace <= 1000
_C547: 43664.9436822 Din_store_furnace - 2778820 Psie_store_furnace
 + Ve2_store_furnace = 0
_C548: - D_store_reactor + Din_store_reactor + Dout_store_reactor = 0
_C549: Din_store_reactor - 63.6396103068 Psie_store_reactor <= 0
_C55: CD_furnace_flash = 0
_C550: Dout_store_reactor + 63.6396103068 Psie_store_reactor >= 63.6396103068
_C551: Dout_store_reactor + 1000 Psie_store_reactor <= 1000
_C552: 66682.055084 Din_store_reactor - 4243620 Psie_store_reactor
 + Ve2_store_reactor = 0
_C553: - D_store_flash + Din_store_flash + Dout_store_flash = 0
_C554: Din_store_flash - 63.6396103068 Psie_store_flash <= 0
_C555: Dout_store_flash + 63.6396103068 Psie_store_flash >= 63.6396103068
_C556: Dout_store_flash + 1000 Psie_store_flash <= 1000
_C557: 43429.5556914 Din_store_flash - 2763840 Psie_store_flash
 + Ve2_store_flash = 0
_C558: - D_store_comp + Din_store_comp + Dout_store_comp = 0
_C559: Din_store_comp - 63.6396103068 Psie_store_comp <= 0
_C56: - L_furnace_flash + R_furnace_flash + X_flash - X_furnace = 0
_C560: Dout_store_comp + 63.6396103068 Psie_store_comp >= 63.6396103068
_C561: Dout_store_comp + 1000 Psie_store_comp <= 1000
_C562: 100984.904983 Din_store_comp - 6426640 Psie_store_comp + Ve2_store_comp
 = 0
_C563: - D_store_distil + Din_store_distil + Dout_store_distil = 0
_C564: Din_store_distil - 63.6396103068 Psie_store_distil <= 0
_C565: Dout_store_distil + 63.6396103068 Psie_store_distil >= 63.6396103068
_C566: Dout_store_distil + 1000 Psie_store_distil <= 1000
_C567: 130557.681921 Din_store_distil - 8308640 Psie_store_distil
 + Ve2_store_distil = 0
_C568: Ve2_store_store = 0
_C569: - D_store_ctrlroom + Din_store_ctrlroom + Dout_store_ctrlroom = 0
_C57: A_furnace_flash - B_furnace_flash + Y_flash - Y_furnace = 0
_C570: Din_store_ctrlroom - 63.6396103068 Psie_store_ctrlroom <= 0
_C571: Dout_store_ctrlroom + 63.6396103068 Psie_store_ctrlroom
 >= 63.6396103068
_C572: Dout_store_ctrlroom + 1000 Psie_store_ctrlroom <= 1000
_C573: 608536.253024 Din_store_ctrlroom - 38727010 Psie_store_ctrlroom
 + Ve2_store_ctrlroom = 0
_C574: - 0.001 Ve2_furnace_comp - 0.001 Ve2_furnace_ctrlroom
 - 0.001 Ve2_furnace_distil - 0.001 Ve2_furnace_flash
 - 0.001 Ve2_furnace_furnace - 0.001 Ve2_furnace_reactor
 - 0.001 Ve2_furnace_store + Ve_furnace = 0
_C575: Omega0_furnace - 0.75 Ve_furnace = 0
_C576: - 0.001 Ve2_reactor_comp - 0.001 Ve2_reactor_ctrlroom
 - 0.001 Ve2_reactor_distil - 0.001 Ve2_reactor_flash
 - 0.001 Ve2_reactor_furnace - 0.001 Ve2_reactor_reactor
 - 0.001 Ve2_reactor_store + Ve_reactor = 0
_C577: Omega0_reactor - 0.82 Ve_reactor = 0
_C578: - 0.001 Ve2_distil_comp - 0.001 Ve2_distil_ctrlroom
 - 0.001 Ve2_distil_distil - 0.001 Ve2_distil_flash - 0.001 Ve2_distil_furnace
 - 0.001 Ve2_distil_reactor - 0.001 Ve2_distil_store + Ve_distil = 0
_C579: Omega0_distil - 0.77 Ve_distil = 0
_C58: R_furnace_flash - 1000 Wx_furnace_flash <= 0
_C580: - 0.001 Ve2_store_comp - 0.001 Ve2_store_ctrlroom
 - 0.001 Ve2_store_distil - 0.001 Ve2_store_flash - 0.001 Ve2_store_furnace
 - 0.001 Ve2_store_reactor - 0.001 Ve2_store_store + Ve_store = 0
_C581: Omega0_store - 0.82 Ve_store = 0
_C582: - Omega0_distil - Omega0_furnace - Omega0_reactor - Omega0_store
 + SumOmega = 0
_C59: L_furnace_flash + 1000 Wx_furnace_flash <= 1000
_C6: d_reactor + l_reactor = 14.4
_C60: A_furnace_flash - 1000 Wy_furnace_flash <= 0
_C61: B_furnace_flash + 1000 Wy_furnace_flash <= 1000
_C62: - A_furnace_flash - B_furnace_flash + D_furnace_flash - L_furnace_flash
 - R_furnace_flash = 0
_C63: - D_flash_furnace + D_furnace_flash = 0
_C64: L_furnace_flash - R_flash_furnace = 0
_C65: - L_flash_furnace + R_furnace_flash = 0
_C66: A_furnace_flash - B_flash_furnace = 0
_C67: - A_flash_furnace + B_furnace_flash = 0
_C68: Wx_flash_furnace + Wx_furnace_flash = 1
_C69: Wy_flash_furnace + Wy_furnace_flash = 1
_C7: X_reactor - 0.5 l_reactor >= 0
_C70: 1000 E1_furnace_flash + 1000 E2_furnace_flash - X_flash + X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= 0
_C71: - 1000 E1_furnace_flash + 1000 E2_furnace_flash + X_flash - X_furnace
 - 0.5 l_flash - 0.5 l_furnace >= -1000
_C72: 1000 E1_furnace_flash - 1000 E2_furnace_flash - Y_flash + Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -1000
_C73: - 1000 E1_furnace_flash - 1000 E2_furnace_flash + Y_flash - Y_furnace
 - 0.5 d_flash - 0.5 d_furnace >= -2000
_C74: CD_reactor_flash - 1231.60506417 D_reactor_flash = 0
_C75: - L_reactor_flash + R_reactor_flash + X_flash - X_reactor = 0
_C76: A_reactor_flash - B_reactor_flash + Y_flash - Y_reactor = 0
_C77: R_reactor_flash - 1000 Wx_reactor_flash <= 0
_C78: L_reactor_flash + 1000 Wx_reactor_flash <= 1000
_C79: A_reactor_flash - 1000 Wy_reactor_flash <= 0
_C8: Y_reactor - 0.5 d_reactor >= 0
_C80: B_reactor_flash + 1000 Wy_reactor_flash <= 1000
_C81: - A_reactor_flash - B_reactor_flash + D_reactor_flash - L_reactor_flash
 - R_reactor_flash = 0
_C82: - D_flash_reactor + D_reactor_flash = 0
_C83: L_reactor_flash - R_flash_reactor = 0
_C84: - L_flash_reactor + R_reactor_flash = 0
_C85: A_reactor_flash - B_flash_reactor = 0
_C86: - A_flash_reactor + B_reactor_flash = 0
_C87: Wx_flash_reactor + Wx_reactor_flash = 1
_C88: Wy_flash_reactor + Wy_reactor_flash = 1
_C89: 1000 E1_reactor_flash + 1000 E2_reactor_flash - X_flash + X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= 0
_C9: 0 O_flash + l_flash = 4
_C90: - 1000 E1_reactor_flash + 1000 E2_reactor_flash + X_flash - X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= -1000
_C91: 1000 E1_reactor_flash - 1000 E2_reactor_flash - Y_flash + Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -1000
_C92: - 1000 E1_reactor_flash - 1000 E2_reactor_flash + Y_flash - Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -2000
_C93: CD_furnace_comp - 1487.20883073 D_furnace_comp = 0
_C94: - L_furnace_comp + R_furnace_comp + X_comp - X_furnace = 0
_C95: A_furnace_comp - B_furnace_comp + Y_comp - Y_furnace = 0
_C96: R_furnace_comp - 1000 Wx_furnace_comp <= 0
_C97: L_furnace_comp + 1000 Wx_furnace_comp <= 1000
_C98: A_furnace_comp - 1000 Wy_furnace_comp <= 0
_C99: B_furnace_comp + 1000 Wy_furnace_comp <= 1000
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
