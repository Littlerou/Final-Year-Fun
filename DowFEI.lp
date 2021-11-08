\* Layout_Problem_Model *\
Minimize
OBJ: SumCD + SumOmega + SumPZ
Subject To
_C1: - 5 O_reactor + l_reactor = 5
_C10: d_eoabs + l_eoabs = 15.36
_C100: - 1000 E1_reactor_hex2 + 1000 E2_reactor_hex2 + X_hex2 - X_reactor
 - 0.5 l_hex2 - 0.5 l_reactor >= -1000
_C101: 1000 E1_reactor_hex2 - 1000 E2_reactor_hex2 - Y_hex2 + Y_reactor
 - 0.5 d_hex2 - 0.5 d_reactor >= -1000
_C102: - 1000 E1_reactor_hex2 - 1000 E2_reactor_hex2 + Y_hex2 - Y_reactor
 - 0.5 d_hex2 - 0.5 d_reactor >= -2000
_C103: L_reactor_hex2 - R_hex2_reactor = 0
_C104: - L_hex2_reactor + R_reactor_hex2 = 0
_C105: A_reactor_hex2 - B_hex2_reactor = 0
_C106: - A_hex2_reactor + B_reactor_hex2 = 0
_C107: Wx_hex2_reactor + Wx_reactor_hex2 = 1
_C108: Wy_hex2_reactor + Wy_reactor_hex2 = 1
_C109: - 56 D_hex1_hex2 + deltaP_hex1_hex2 = 0
_C11: X_eoabs - 0.5 l_eoabs >= 0
_C110: CD_hex1_hex2 = 0
_C111: - L_hex1_hex2 + R_hex1_hex2 - X_hex1 + X_hex2 = 0
_C112: A_hex1_hex2 - B_hex1_hex2 - Y_hex1 + Y_hex2 = 0
_C113: R_hex1_hex2 - 1000 Wx_hex1_hex2 <= 0
_C114: L_hex1_hex2 + 1000 Wx_hex1_hex2 <= 1000
_C115: A_hex1_hex2 - 1000 Wy_hex1_hex2 <= 0
_C116: B_hex1_hex2 + 1000 Wy_hex1_hex2 <= 1000
_C117: - A_hex1_hex2 - B_hex1_hex2 + D_hex1_hex2 - L_hex1_hex2 - R_hex1_hex2
 = 0
_C118: D_hex1_hex2 - D_hex2_hex1 = 0
_C119: 1000 E1_hex1_hex2 + 1000 E2_hex1_hex2 + X_hex1 - X_hex2 - 0.5 l_hex1
 - 0.5 l_hex2 >= 0
_C12: Y_eoabs - 0.5 d_eoabs >= 0
_C120: - 1000 E1_hex1_hex2 + 1000 E2_hex1_hex2 - X_hex1 + X_hex2 - 0.5 l_hex1
 - 0.5 l_hex2 >= -1000
_C121: 1000 E1_hex1_hex2 - 1000 E2_hex1_hex2 + Y_hex1 - Y_hex2 - 0.5 d_hex1
 - 0.5 d_hex2 >= -1000
_C122: - 1000 E1_hex1_hex2 - 1000 E2_hex1_hex2 - Y_hex1 + Y_hex2 - 0.5 d_hex1
 - 0.5 d_hex2 >= -2000
_C123: L_hex1_hex2 - R_hex2_hex1 = 0
_C124: - L_hex2_hex1 + R_hex1_hex2 = 0
_C125: A_hex1_hex2 - B_hex2_hex1 = 0
_C126: - A_hex2_hex1 + B_hex1_hex2 = 0
_C127: Wx_hex1_hex2 + Wx_hex2_hex1 = 1
_C128: Wy_hex1_hex2 + Wy_hex2_hex1 = 1
_C129: - 56 D_eoabs_hex2 + deltaP_eoabs_hex2 = 0
_C13: 0 O_hex2 + l_hex2 = 8.48
_C130: CD_eoabs_hex2 - 27.3003837522 D_eoabs_hex2 = 0
_C131: - L_eoabs_hex2 + R_eoabs_hex2 - X_eoabs + X_hex2 = 0
_C132: A_eoabs_hex2 - B_eoabs_hex2 - Y_eoabs + Y_hex2 = 0
_C133: R_eoabs_hex2 - 1000 Wx_eoabs_hex2 <= 0
_C134: L_eoabs_hex2 + 1000 Wx_eoabs_hex2 <= 1000
_C135: A_eoabs_hex2 - 1000 Wy_eoabs_hex2 <= 0
_C136: B_eoabs_hex2 + 1000 Wy_eoabs_hex2 <= 1000
_C137: - A_eoabs_hex2 - B_eoabs_hex2 + D_eoabs_hex2 - L_eoabs_hex2
 - R_eoabs_hex2 = 0
_C138: D_eoabs_hex2 - D_hex2_eoabs = 0
_C139: 1000 E1_eoabs_hex2 + 1000 E2_eoabs_hex2 + X_eoabs - X_hex2
 - 0.5 l_eoabs - 0.5 l_hex2 >= 0
_C14: d_hex2 + l_hex2 = 16.96
_C140: - 1000 E1_eoabs_hex2 + 1000 E2_eoabs_hex2 - X_eoabs + X_hex2
 - 0.5 l_eoabs - 0.5 l_hex2 >= -1000
_C141: 1000 E1_eoabs_hex2 - 1000 E2_eoabs_hex2 + Y_eoabs - Y_hex2
 - 0.5 d_eoabs - 0.5 d_hex2 >= -1000
_C142: - 1000 E1_eoabs_hex2 - 1000 E2_eoabs_hex2 - Y_eoabs + Y_hex2
 - 0.5 d_eoabs - 0.5 d_hex2 >= -2000
_C143: L_eoabs_hex2 - R_hex2_eoabs = 0
_C144: - L_hex2_eoabs + R_eoabs_hex2 = 0
_C145: A_eoabs_hex2 - B_hex2_eoabs = 0
_C146: - A_hex2_eoabs + B_eoabs_hex2 = 0
_C147: Wx_eoabs_hex2 + Wx_hex2_eoabs = 1
_C148: Wy_eoabs_hex2 + Wy_hex2_eoabs = 1
_C149: - 56 D_reactor_co2abs + deltaP_reactor_co2abs = 0
_C15: X_hex2 - 0.5 l_hex2 >= 0
_C150: CD_reactor_co2abs - 27.3003837522 D_reactor_co2abs = 0
_C151: - L_reactor_co2abs + R_reactor_co2abs + X_co2abs - X_reactor = 0
_C152: A_reactor_co2abs - B_reactor_co2abs + Y_co2abs - Y_reactor = 0
_C153: R_reactor_co2abs - 1000 Wx_reactor_co2abs <= 0
_C154: L_reactor_co2abs + 1000 Wx_reactor_co2abs <= 1000
_C155: A_reactor_co2abs - 1000 Wy_reactor_co2abs <= 0
_C156: B_reactor_co2abs + 1000 Wy_reactor_co2abs <= 1000
_C157: - A_reactor_co2abs - B_reactor_co2abs + D_reactor_co2abs
 - L_reactor_co2abs - R_reactor_co2abs = 0
_C158: - D_co2abs_reactor + D_reactor_co2abs = 0
_C159: 1000 E1_reactor_co2abs + 1000 E2_reactor_co2abs - X_co2abs + X_reactor
 - 0.5 l_co2abs - 0.5 l_reactor >= 0
_C16: Y_hex2 - 0.5 d_hex2 >= 0
_C160: - 1000 E1_reactor_co2abs + 1000 E2_reactor_co2abs + X_co2abs
 - X_reactor - 0.5 l_co2abs - 0.5 l_reactor >= -1000
_C161: 1000 E1_reactor_co2abs - 1000 E2_reactor_co2abs - Y_co2abs + Y_reactor
 - 0.5 d_co2abs - 0.5 d_reactor >= -1000
_C162: - 1000 E1_reactor_co2abs - 1000 E2_reactor_co2abs + Y_co2abs
 - Y_reactor - 0.5 d_co2abs - 0.5 d_reactor >= -2000
_C163: L_reactor_co2abs - R_co2abs_reactor = 0
_C164: - L_co2abs_reactor + R_reactor_co2abs = 0
_C165: A_reactor_co2abs - B_co2abs_reactor = 0
_C166: - A_co2abs_reactor + B_reactor_co2abs = 0
_C167: Wx_co2abs_reactor + Wx_reactor_co2abs = 1
_C168: Wy_co2abs_reactor + Wy_reactor_co2abs = 1
_C169: - 56 D_hex1_co2abs + deltaP_hex1_co2abs = 0
_C17: 0 O_co2abs + l_co2abs = 5
_C170: CD_hex1_co2abs = 0
_C171: - L_hex1_co2abs + R_hex1_co2abs + X_co2abs - X_hex1 = 0
_C172: A_hex1_co2abs - B_hex1_co2abs + Y_co2abs - Y_hex1 = 0
_C173: R_hex1_co2abs - 1000 Wx_hex1_co2abs <= 0
_C174: L_hex1_co2abs + 1000 Wx_hex1_co2abs <= 1000
_C175: A_hex1_co2abs - 1000 Wy_hex1_co2abs <= 0
_C176: B_hex1_co2abs + 1000 Wy_hex1_co2abs <= 1000
_C177: - A_hex1_co2abs - B_hex1_co2abs + D_hex1_co2abs - L_hex1_co2abs
 - R_hex1_co2abs = 0
_C178: - D_co2abs_hex1 + D_hex1_co2abs = 0
_C179: 1000 E1_hex1_co2abs + 1000 E2_hex1_co2abs - X_co2abs + X_hex1
 - 0.5 l_co2abs - 0.5 l_hex1 >= 0
_C18: d_co2abs + l_co2abs = 10
_C180: - 1000 E1_hex1_co2abs + 1000 E2_hex1_co2abs + X_co2abs - X_hex1
 - 0.5 l_co2abs - 0.5 l_hex1 >= -1000
_C181: 1000 E1_hex1_co2abs - 1000 E2_hex1_co2abs - Y_co2abs + Y_hex1
 - 0.5 d_co2abs - 0.5 d_hex1 >= -1000
_C182: - 1000 E1_hex1_co2abs - 1000 E2_hex1_co2abs + Y_co2abs - Y_hex1
 - 0.5 d_co2abs - 0.5 d_hex1 >= -2000
_C183: L_hex1_co2abs - R_co2abs_hex1 = 0
_C184: - L_co2abs_hex1 + R_hex1_co2abs = 0
_C185: A_hex1_co2abs - B_co2abs_hex1 = 0
_C186: - A_co2abs_hex1 + B_hex1_co2abs = 0
_C187: Wx_co2abs_hex1 + Wx_hex1_co2abs = 1
_C188: Wy_co2abs_hex1 + Wy_hex1_co2abs = 1
_C189: - 56 D_eoabs_co2abs + deltaP_eoabs_co2abs = 0
_C19: X_co2abs - 0.5 l_co2abs >= 0
_C190: CD_eoabs_co2abs = 0
_C191: - L_eoabs_co2abs + R_eoabs_co2abs + X_co2abs - X_eoabs = 0
_C192: A_eoabs_co2abs - B_eoabs_co2abs + Y_co2abs - Y_eoabs = 0
_C193: R_eoabs_co2abs - 1000 Wx_eoabs_co2abs <= 0
_C194: L_eoabs_co2abs + 1000 Wx_eoabs_co2abs <= 1000
_C195: A_eoabs_co2abs - 1000 Wy_eoabs_co2abs <= 0
_C196: B_eoabs_co2abs + 1000 Wy_eoabs_co2abs <= 1000
_C197: - A_eoabs_co2abs - B_eoabs_co2abs + D_eoabs_co2abs - L_eoabs_co2abs
 - R_eoabs_co2abs = 0
_C198: - D_co2abs_eoabs + D_eoabs_co2abs = 0
_C199: 1000 E1_eoabs_co2abs + 1000 E2_eoabs_co2abs - X_co2abs + X_eoabs
 - 0.5 l_co2abs - 0.5 l_eoabs >= 0
_C2: d_reactor + l_reactor = 15
_C20: Y_co2abs - 0.5 d_co2abs >= 0
_C200: - 1000 E1_eoabs_co2abs + 1000 E2_eoabs_co2abs + X_co2abs - X_eoabs
 - 0.5 l_co2abs - 0.5 l_eoabs >= -1000
_C201: 1000 E1_eoabs_co2abs - 1000 E2_eoabs_co2abs - Y_co2abs + Y_eoabs
 - 0.5 d_co2abs - 0.5 d_eoabs >= -1000
_C202: - 1000 E1_eoabs_co2abs - 1000 E2_eoabs_co2abs + Y_co2abs - Y_eoabs
 - 0.5 d_co2abs - 0.5 d_eoabs >= -2000
_C203: L_eoabs_co2abs - R_co2abs_eoabs = 0
_C204: - L_co2abs_eoabs + R_eoabs_co2abs = 0
_C205: A_eoabs_co2abs - B_co2abs_eoabs = 0
_C206: - A_co2abs_eoabs + B_eoabs_co2abs = 0
_C207: Wx_co2abs_eoabs + Wx_eoabs_co2abs = 1
_C208: Wy_co2abs_eoabs + Wy_eoabs_co2abs = 1
_C209: - 56 D_hex2_co2abs + deltaP_hex2_co2abs = 0
_C21: 0 O_flash + l_flash = 2.68
_C210: CD_hex2_co2abs - 27.3003837522 D_hex2_co2abs = 0
_C211: - L_hex2_co2abs + R_hex2_co2abs + X_co2abs - X_hex2 = 0
_C212: A_hex2_co2abs - B_hex2_co2abs + Y_co2abs - Y_hex2 = 0
_C213: R_hex2_co2abs - 1000 Wx_hex2_co2abs <= 0
_C214: L_hex2_co2abs + 1000 Wx_hex2_co2abs <= 1000
_C215: A_hex2_co2abs - 1000 Wy_hex2_co2abs <= 0
_C216: B_hex2_co2abs + 1000 Wy_hex2_co2abs <= 1000
_C217: - A_hex2_co2abs - B_hex2_co2abs + D_hex2_co2abs - L_hex2_co2abs
 - R_hex2_co2abs = 0
_C218: - D_co2abs_hex2 + D_hex2_co2abs = 0
_C219: 1000 E1_hex2_co2abs + 1000 E2_hex2_co2abs - X_co2abs + X_hex2
 - 0.5 l_co2abs - 0.5 l_hex2 >= 0
_C22: d_flash + l_flash = 5.36
_C220: - 1000 E1_hex2_co2abs + 1000 E2_hex2_co2abs + X_co2abs - X_hex2
 - 0.5 l_co2abs - 0.5 l_hex2 >= -1000
_C221: 1000 E1_hex2_co2abs - 1000 E2_hex2_co2abs - Y_co2abs + Y_hex2
 - 0.5 d_co2abs - 0.5 d_hex2 >= -1000
_C222: - 1000 E1_hex2_co2abs - 1000 E2_hex2_co2abs + Y_co2abs - Y_hex2
 - 0.5 d_co2abs - 0.5 d_hex2 >= -2000
_C223: L_hex2_co2abs - R_co2abs_hex2 = 0
_C224: - L_co2abs_hex2 + R_hex2_co2abs = 0
_C225: A_hex2_co2abs - B_co2abs_hex2 = 0
_C226: - A_co2abs_hex2 + B_hex2_co2abs = 0
_C227: Wx_co2abs_hex2 + Wx_hex2_co2abs = 1
_C228: Wy_co2abs_hex2 + Wy_hex2_co2abs = 1
_C229: - 56 D_reactor_flash + deltaP_reactor_flash = 0
_C23: X_flash - 0.5 l_flash >= 0
_C230: CD_reactor_flash = 0
_C231: - L_reactor_flash + R_reactor_flash + X_flash - X_reactor = 0
_C232: A_reactor_flash - B_reactor_flash + Y_flash - Y_reactor = 0
_C233: R_reactor_flash - 1000 Wx_reactor_flash <= 0
_C234: L_reactor_flash + 1000 Wx_reactor_flash <= 1000
_C235: A_reactor_flash - 1000 Wy_reactor_flash <= 0
_C236: B_reactor_flash + 1000 Wy_reactor_flash <= 1000
_C237: - A_reactor_flash - B_reactor_flash + D_reactor_flash - L_reactor_flash
 - R_reactor_flash = 0
_C238: - D_flash_reactor + D_reactor_flash = 0
_C239: 1000 E1_reactor_flash + 1000 E2_reactor_flash - X_flash + X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= 0
_C24: Y_flash - 0.5 d_flash >= 0
_C240: - 1000 E1_reactor_flash + 1000 E2_reactor_flash + X_flash - X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= -1000
_C241: 1000 E1_reactor_flash - 1000 E2_reactor_flash - Y_flash + Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -1000
_C242: - 1000 E1_reactor_flash - 1000 E2_reactor_flash + Y_flash - Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -2000
_C243: L_reactor_flash - R_flash_reactor = 0
_C244: - L_flash_reactor + R_reactor_flash = 0
_C245: A_reactor_flash - B_flash_reactor = 0
_C246: - A_flash_reactor + B_reactor_flash = 0
_C247: Wx_flash_reactor + Wx_reactor_flash = 1
_C248: Wy_flash_reactor + Wy_reactor_flash = 1
_C249: - 56 D_hex1_flash + deltaP_hex1_flash = 0
_C25: 0 O_pump + l_pump = 2.4
_C250: CD_hex1_flash = 0
_C251: - L_hex1_flash + R_hex1_flash + X_flash - X_hex1 = 0
_C252: A_hex1_flash - B_hex1_flash + Y_flash - Y_hex1 = 0
_C253: R_hex1_flash - 1000 Wx_hex1_flash <= 0
_C254: L_hex1_flash + 1000 Wx_hex1_flash <= 1000
_C255: A_hex1_flash - 1000 Wy_hex1_flash <= 0
_C256: B_hex1_flash + 1000 Wy_hex1_flash <= 1000
_C257: - A_hex1_flash - B_hex1_flash + D_hex1_flash - L_hex1_flash
 - R_hex1_flash = 0
_C258: - D_flash_hex1 + D_hex1_flash = 0
_C259: 1000 E1_hex1_flash + 1000 E2_hex1_flash - X_flash + X_hex1
 - 0.5 l_flash - 0.5 l_hex1 >= 0
_C26: d_pump + l_pump = 4.8
_C260: - 1000 E1_hex1_flash + 1000 E2_hex1_flash + X_flash - X_hex1
 - 0.5 l_flash - 0.5 l_hex1 >= -1000
_C261: 1000 E1_hex1_flash - 1000 E2_hex1_flash - Y_flash + Y_hex1
 - 0.5 d_flash - 0.5 d_hex1 >= -1000
_C262: - 1000 E1_hex1_flash - 1000 E2_hex1_flash + Y_flash - Y_hex1
 - 0.5 d_flash - 0.5 d_hex1 >= -2000
_C263: L_hex1_flash - R_flash_hex1 = 0
_C264: - L_flash_hex1 + R_hex1_flash = 0
_C265: A_hex1_flash - B_flash_hex1 = 0
_C266: - A_flash_hex1 + B_hex1_flash = 0
_C267: Wx_flash_hex1 + Wx_hex1_flash = 1
_C268: Wy_flash_hex1 + Wy_hex1_flash = 1
_C269: - 56 D_eoabs_flash + deltaP_eoabs_flash = 0
_C27: X_pump - 0.5 l_pump >= 0
_C270: CD_eoabs_flash = 0
_C271: - L_eoabs_flash + R_eoabs_flash - X_eoabs + X_flash = 0
_C272: A_eoabs_flash - B_eoabs_flash - Y_eoabs + Y_flash = 0
_C273: R_eoabs_flash - 1000 Wx_eoabs_flash <= 0
_C274: L_eoabs_flash + 1000 Wx_eoabs_flash <= 1000
_C275: A_eoabs_flash - 1000 Wy_eoabs_flash <= 0
_C276: B_eoabs_flash + 1000 Wy_eoabs_flash <= 1000
_C277: - A_eoabs_flash - B_eoabs_flash + D_eoabs_flash - L_eoabs_flash
 - R_eoabs_flash = 0
_C278: D_eoabs_flash - D_flash_eoabs = 0
_C279: 1000 E1_eoabs_flash + 1000 E2_eoabs_flash + X_eoabs - X_flash
 - 0.5 l_eoabs - 0.5 l_flash >= 0
_C28: Y_pump - 0.5 d_pump >= 0
_C280: - 1000 E1_eoabs_flash + 1000 E2_eoabs_flash - X_eoabs + X_flash
 - 0.5 l_eoabs - 0.5 l_flash >= -1000
_C281: 1000 E1_eoabs_flash - 1000 E2_eoabs_flash + Y_eoabs - Y_flash
 - 0.5 d_eoabs - 0.5 d_flash >= -1000
_C282: - 1000 E1_eoabs_flash - 1000 E2_eoabs_flash - Y_eoabs + Y_flash
 - 0.5 d_eoabs - 0.5 d_flash >= -2000
_C283: L_eoabs_flash - R_flash_eoabs = 0
_C284: - L_flash_eoabs + R_eoabs_flash = 0
_C285: A_eoabs_flash - B_flash_eoabs = 0
_C286: - A_flash_eoabs + B_eoabs_flash = 0
_C287: Wx_eoabs_flash + Wx_flash_eoabs = 1
_C288: Wy_eoabs_flash + Wy_flash_eoabs = 1
_C289: - 56 D_hex2_flash + deltaP_hex2_flash = 0
_C29: - 56 D_reactor_hex1 + deltaP_reactor_hex1 = 0
_C290: CD_hex2_flash = 0
_C291: - L_hex2_flash + R_hex2_flash + X_flash - X_hex2 = 0
_C292: A_hex2_flash - B_hex2_flash + Y_flash - Y_hex2 = 0
_C293: R_hex2_flash - 1000 Wx_hex2_flash <= 0
_C294: L_hex2_flash + 1000 Wx_hex2_flash <= 1000
_C295: A_hex2_flash - 1000 Wy_hex2_flash <= 0
_C296: B_hex2_flash + 1000 Wy_hex2_flash <= 1000
_C297: - A_hex2_flash - B_hex2_flash + D_hex2_flash - L_hex2_flash
 - R_hex2_flash = 0
_C298: - D_flash_hex2 + D_hex2_flash = 0
_C299: 1000 E1_hex2_flash + 1000 E2_hex2_flash - X_flash + X_hex2
 - 0.5 l_flash - 0.5 l_hex2 >= 0
_C3: X_reactor - 0.5 l_reactor >= 0
_C30: CD_reactor_hex1 - 27.3003837522 D_reactor_hex1 = 0
_C300: - 1000 E1_hex2_flash + 1000 E2_hex2_flash + X_flash - X_hex2
 - 0.5 l_flash - 0.5 l_hex2 >= -1000
_C301: 1000 E1_hex2_flash - 1000 E2_hex2_flash - Y_flash + Y_hex2
 - 0.5 d_flash - 0.5 d_hex2 >= -1000
_C302: - 1000 E1_hex2_flash - 1000 E2_hex2_flash + Y_flash - Y_hex2
 - 0.5 d_flash - 0.5 d_hex2 >= -2000
_C303: L_hex2_flash - R_flash_hex2 = 0
_C304: - L_flash_hex2 + R_hex2_flash = 0
_C305: A_hex2_flash - B_flash_hex2 = 0
_C306: - A_flash_hex2 + B_hex2_flash = 0
_C307: Wx_flash_hex2 + Wx_hex2_flash = 1
_C308: Wy_flash_hex2 + Wy_hex2_flash = 1
_C309: - 56 D_co2abs_flash + deltaP_co2abs_flash = 0
_C31: - L_reactor_hex1 + R_reactor_hex1 + X_hex1 - X_reactor = 0
_C310: CD_co2abs_flash - 27.3003837522 D_co2abs_flash = 0
_C311: - L_co2abs_flash + R_co2abs_flash - X_co2abs + X_flash = 0
_C312: A_co2abs_flash - B_co2abs_flash - Y_co2abs + Y_flash = 0
_C313: R_co2abs_flash - 1000 Wx_co2abs_flash <= 0
_C314: L_co2abs_flash + 1000 Wx_co2abs_flash <= 1000
_C315: A_co2abs_flash - 1000 Wy_co2abs_flash <= 0
_C316: B_co2abs_flash + 1000 Wy_co2abs_flash <= 1000
_C317: - A_co2abs_flash - B_co2abs_flash + D_co2abs_flash - L_co2abs_flash
 - R_co2abs_flash = 0
_C318: D_co2abs_flash - D_flash_co2abs = 0
_C319: 1000 E1_co2abs_flash + 1000 E2_co2abs_flash + X_co2abs - X_flash
 - 0.5 l_co2abs - 0.5 l_flash >= 0
_C32: A_reactor_hex1 - B_reactor_hex1 + Y_hex1 - Y_reactor = 0
_C320: - 1000 E1_co2abs_flash + 1000 E2_co2abs_flash - X_co2abs + X_flash
 - 0.5 l_co2abs - 0.5 l_flash >= -1000
_C321: 1000 E1_co2abs_flash - 1000 E2_co2abs_flash + Y_co2abs - Y_flash
 - 0.5 d_co2abs - 0.5 d_flash >= -1000
_C322: - 1000 E1_co2abs_flash - 1000 E2_co2abs_flash - Y_co2abs + Y_flash
 - 0.5 d_co2abs - 0.5 d_flash >= -2000
_C323: L_co2abs_flash - R_flash_co2abs = 0
_C324: - L_flash_co2abs + R_co2abs_flash = 0
_C325: A_co2abs_flash - B_flash_co2abs = 0
_C326: - A_flash_co2abs + B_co2abs_flash = 0
_C327: Wx_co2abs_flash + Wx_flash_co2abs = 1
_C328: Wy_co2abs_flash + Wy_flash_co2abs = 1
_C329: - 56 D_reactor_pump + deltaP_reactor_pump = 0
_C33: R_reactor_hex1 - 1000 Wx_reactor_hex1 <= 0
_C330: CD_reactor_pump = 0
_C331: - L_reactor_pump + R_reactor_pump + X_pump - X_reactor = 0
_C332: A_reactor_pump - B_reactor_pump + Y_pump - Y_reactor = 0
_C333: R_reactor_pump - 1000 Wx_reactor_pump <= 0
_C334: L_reactor_pump + 1000 Wx_reactor_pump <= 1000
_C335: A_reactor_pump - 1000 Wy_reactor_pump <= 0
_C336: B_reactor_pump + 1000 Wy_reactor_pump <= 1000
_C337: - A_reactor_pump - B_reactor_pump + D_reactor_pump - L_reactor_pump
 - R_reactor_pump = 0
_C338: - D_pump_reactor + D_reactor_pump = 0
_C339: 1000 E1_reactor_pump + 1000 E2_reactor_pump - X_pump + X_reactor
 - 0.5 l_pump - 0.5 l_reactor >= 0
_C34: L_reactor_hex1 + 1000 Wx_reactor_hex1 <= 1000
_C340: - 1000 E1_reactor_pump + 1000 E2_reactor_pump + X_pump - X_reactor
 - 0.5 l_pump - 0.5 l_reactor >= -1000
_C341: 1000 E1_reactor_pump - 1000 E2_reactor_pump - Y_pump + Y_reactor
 - 0.5 d_pump - 0.5 d_reactor >= -1000
_C342: - 1000 E1_reactor_pump - 1000 E2_reactor_pump + Y_pump - Y_reactor
 - 0.5 d_pump - 0.5 d_reactor >= -2000
_C343: L_reactor_pump - R_pump_reactor = 0
_C344: - L_pump_reactor + R_reactor_pump = 0
_C345: A_reactor_pump - B_pump_reactor = 0
_C346: - A_pump_reactor + B_reactor_pump = 0
_C347: Wx_pump_reactor + Wx_reactor_pump = 1
_C348: Wy_pump_reactor + Wy_reactor_pump = 1
_C349: - 56 D_hex1_pump + deltaP_hex1_pump = 0
_C35: A_reactor_hex1 - 1000 Wy_reactor_hex1 <= 0
_C350: CD_hex1_pump = 0
_C351: - L_hex1_pump + R_hex1_pump - X_hex1 + X_pump = 0
_C352: A_hex1_pump - B_hex1_pump - Y_hex1 + Y_pump = 0
_C353: R_hex1_pump - 1000 Wx_hex1_pump <= 0
_C354: L_hex1_pump + 1000 Wx_hex1_pump <= 1000
_C355: A_hex1_pump - 1000 Wy_hex1_pump <= 0
_C356: B_hex1_pump + 1000 Wy_hex1_pump <= 1000
_C357: - A_hex1_pump - B_hex1_pump + D_hex1_pump - L_hex1_pump - R_hex1_pump
 = 0
_C358: D_hex1_pump - D_pump_hex1 = 0
_C359: 1000 E1_hex1_pump + 1000 E2_hex1_pump + X_hex1 - X_pump - 0.5 l_hex1
 - 0.5 l_pump >= 0
_C36: B_reactor_hex1 + 1000 Wy_reactor_hex1 <= 1000
_C360: - 1000 E1_hex1_pump + 1000 E2_hex1_pump - X_hex1 + X_pump - 0.5 l_hex1
 - 0.5 l_pump >= -1000
_C361: 1000 E1_hex1_pump - 1000 E2_hex1_pump + Y_hex1 - Y_pump - 0.5 d_hex1
 - 0.5 d_pump >= -1000
_C362: - 1000 E1_hex1_pump - 1000 E2_hex1_pump - Y_hex1 + Y_pump - 0.5 d_hex1
 - 0.5 d_pump >= -2000
_C363: L_hex1_pump - R_pump_hex1 = 0
_C364: - L_pump_hex1 + R_hex1_pump = 0
_C365: A_hex1_pump - B_pump_hex1 = 0
_C366: - A_pump_hex1 + B_hex1_pump = 0
_C367: Wx_hex1_pump + Wx_pump_hex1 = 1
_C368: Wy_hex1_pump + Wy_pump_hex1 = 1
_C369: - 56 D_eoabs_pump + deltaP_eoabs_pump = 0
_C37: - A_reactor_hex1 - B_reactor_hex1 + D_reactor_hex1 - L_reactor_hex1
 - R_reactor_hex1 = 0
_C370: CD_eoabs_pump = 0
_C371: - L_eoabs_pump + R_eoabs_pump - X_eoabs + X_pump = 0
_C372: A_eoabs_pump - B_eoabs_pump - Y_eoabs + Y_pump = 0
_C373: R_eoabs_pump - 1000 Wx_eoabs_pump <= 0
_C374: L_eoabs_pump + 1000 Wx_eoabs_pump <= 1000
_C375: A_eoabs_pump - 1000 Wy_eoabs_pump <= 0
_C376: B_eoabs_pump + 1000 Wy_eoabs_pump <= 1000
_C377: - A_eoabs_pump - B_eoabs_pump + D_eoabs_pump - L_eoabs_pump
 - R_eoabs_pump = 0
_C378: D_eoabs_pump - D_pump_eoabs = 0
_C379: 1000 E1_eoabs_pump + 1000 E2_eoabs_pump + X_eoabs - X_pump
 - 0.5 l_eoabs - 0.5 l_pump >= 0
_C38: - D_hex1_reactor + D_reactor_hex1 = 0
_C380: - 1000 E1_eoabs_pump + 1000 E2_eoabs_pump - X_eoabs + X_pump
 - 0.5 l_eoabs - 0.5 l_pump >= -1000
_C381: 1000 E1_eoabs_pump - 1000 E2_eoabs_pump + Y_eoabs - Y_pump
 - 0.5 d_eoabs - 0.5 d_pump >= -1000
_C382: - 1000 E1_eoabs_pump - 1000 E2_eoabs_pump - Y_eoabs + Y_pump
 - 0.5 d_eoabs - 0.5 d_pump >= -2000
_C383: L_eoabs_pump - R_pump_eoabs = 0
_C384: - L_pump_eoabs + R_eoabs_pump = 0
_C385: A_eoabs_pump - B_pump_eoabs = 0
_C386: - A_pump_eoabs + B_eoabs_pump = 0
_C387: Wx_eoabs_pump + Wx_pump_eoabs = 1
_C388: Wy_eoabs_pump + Wy_pump_eoabs = 1
_C389: - 56 D_hex2_pump + deltaP_hex2_pump = 0
_C39: 1000 E1_reactor_hex1 + 1000 E2_reactor_hex1 - X_hex1 + X_reactor
 - 0.5 l_hex1 - 0.5 l_reactor >= 0
_C390: CD_hex2_pump = 0
_C391: - L_hex2_pump + R_hex2_pump - X_hex2 + X_pump = 0
_C392: A_hex2_pump - B_hex2_pump - Y_hex2 + Y_pump = 0
_C393: R_hex2_pump - 1000 Wx_hex2_pump <= 0
_C394: L_hex2_pump + 1000 Wx_hex2_pump <= 1000
_C395: A_hex2_pump - 1000 Wy_hex2_pump <= 0
_C396: B_hex2_pump + 1000 Wy_hex2_pump <= 1000
_C397: - A_hex2_pump - B_hex2_pump + D_hex2_pump - L_hex2_pump - R_hex2_pump
 = 0
_C398: D_hex2_pump - D_pump_hex2 = 0
_C399: 1000 E1_hex2_pump + 1000 E2_hex2_pump + X_hex2 - X_pump - 0.5 l_hex2
 - 0.5 l_pump >= 0
_C4: Y_reactor - 0.5 d_reactor >= 0
_C40: - 1000 E1_reactor_hex1 + 1000 E2_reactor_hex1 + X_hex1 - X_reactor
 - 0.5 l_hex1 - 0.5 l_reactor >= -1000
_C400: - 1000 E1_hex2_pump + 1000 E2_hex2_pump - X_hex2 + X_pump - 0.5 l_hex2
 - 0.5 l_pump >= -1000
_C401: 1000 E1_hex2_pump - 1000 E2_hex2_pump + Y_hex2 - Y_pump - 0.5 d_hex2
 - 0.5 d_pump >= -1000
_C402: - 1000 E1_hex2_pump - 1000 E2_hex2_pump - Y_hex2 + Y_pump - 0.5 d_hex2
 - 0.5 d_pump >= -2000
_C403: L_hex2_pump - R_pump_hex2 = 0
_C404: - L_pump_hex2 + R_hex2_pump = 0
_C405: A_hex2_pump - B_pump_hex2 = 0
_C406: - A_pump_hex2 + B_hex2_pump = 0
_C407: Wx_hex2_pump + Wx_pump_hex2 = 1
_C408: Wy_hex2_pump + Wy_pump_hex2 = 1
_C409: - 56 D_co2abs_pump + deltaP_co2abs_pump = 0
_C41: 1000 E1_reactor_hex1 - 1000 E2_reactor_hex1 - Y_hex1 + Y_reactor
 - 0.5 d_hex1 - 0.5 d_reactor >= -1000
_C410: CD_co2abs_pump - 27.3003837522 D_co2abs_pump = 0
_C411: - L_co2abs_pump + R_co2abs_pump - X_co2abs + X_pump = 0
_C412: A_co2abs_pump - B_co2abs_pump - Y_co2abs + Y_pump = 0
_C413: R_co2abs_pump - 1000 Wx_co2abs_pump <= 0
_C414: L_co2abs_pump + 1000 Wx_co2abs_pump <= 1000
_C415: A_co2abs_pump - 1000 Wy_co2abs_pump <= 0
_C416: B_co2abs_pump + 1000 Wy_co2abs_pump <= 1000
_C417: - A_co2abs_pump - B_co2abs_pump + D_co2abs_pump - L_co2abs_pump
 - R_co2abs_pump = 0
_C418: D_co2abs_pump - D_pump_co2abs = 0
_C419: 1000 E1_co2abs_pump + 1000 E2_co2abs_pump + X_co2abs - X_pump
 - 0.5 l_co2abs - 0.5 l_pump >= 0
_C42: - 1000 E1_reactor_hex1 - 1000 E2_reactor_hex1 + Y_hex1 - Y_reactor
 - 0.5 d_hex1 - 0.5 d_reactor >= -2000
_C420: - 1000 E1_co2abs_pump + 1000 E2_co2abs_pump - X_co2abs + X_pump
 - 0.5 l_co2abs - 0.5 l_pump >= -1000
_C421: 1000 E1_co2abs_pump - 1000 E2_co2abs_pump + Y_co2abs - Y_pump
 - 0.5 d_co2abs - 0.5 d_pump >= -1000
_C422: - 1000 E1_co2abs_pump - 1000 E2_co2abs_pump - Y_co2abs + Y_pump
 - 0.5 d_co2abs - 0.5 d_pump >= -2000
_C423: L_co2abs_pump - R_pump_co2abs = 0
_C424: - L_pump_co2abs + R_co2abs_pump = 0
_C425: A_co2abs_pump - B_pump_co2abs = 0
_C426: - A_pump_co2abs + B_co2abs_pump = 0
_C427: Wx_co2abs_pump + Wx_pump_co2abs = 1
_C428: Wy_co2abs_pump + Wy_pump_co2abs = 1
_C429: - 56 D_flash_pump + deltaP_flash_pump = 0
_C43: L_reactor_hex1 - R_hex1_reactor = 0
_C430: CD_flash_pump - 27.3003837522 D_flash_pump = 0
_C431: - L_flash_pump + R_flash_pump - X_flash + X_pump = 0
_C432: A_flash_pump - B_flash_pump - Y_flash + Y_pump = 0
_C433: R_flash_pump - 1000 Wx_flash_pump <= 0
_C434: L_flash_pump + 1000 Wx_flash_pump <= 1000
_C435: A_flash_pump - 1000 Wy_flash_pump <= 0
_C436: B_flash_pump + 1000 Wy_flash_pump <= 1000
_C437: - A_flash_pump - B_flash_pump + D_flash_pump - L_flash_pump
 - R_flash_pump = 0
_C438: D_flash_pump - D_pump_flash = 0
_C439: 1000 E1_flash_pump + 1000 E2_flash_pump + X_flash - X_pump
 - 0.5 l_flash - 0.5 l_pump >= 0
_C44: - L_hex1_reactor + R_reactor_hex1 = 0
_C440: - 1000 E1_flash_pump + 1000 E2_flash_pump - X_flash + X_pump
 - 0.5 l_flash - 0.5 l_pump >= -1000
_C441: 1000 E1_flash_pump - 1000 E2_flash_pump + Y_flash - Y_pump
 - 0.5 d_flash - 0.5 d_pump >= -1000
_C442: - 1000 E1_flash_pump - 1000 E2_flash_pump - Y_flash + Y_pump
 - 0.5 d_flash - 0.5 d_pump >= -2000
_C443: L_flash_pump - R_pump_flash = 0
_C444: - L_pump_flash + R_flash_pump = 0
_C445: A_flash_pump - B_pump_flash = 0
_C446: - A_pump_flash + B_flash_pump = 0
_C447: Wx_flash_pump + Wx_pump_flash = 1
_C448: Wy_flash_pump + Wy_pump_flash = 1
_C449: - CD_co2abs_co2abs - CD_co2abs_eoabs - CD_co2abs_flash - CD_co2abs_hex1
 - CD_co2abs_hex2 - CD_co2abs_pump - CD_co2abs_reactor - CD_eoabs_co2abs
 - CD_eoabs_eoabs - CD_eoabs_flash - CD_eoabs_hex1 - CD_eoabs_hex2
 - CD_eoabs_pump - CD_eoabs_reactor - CD_flash_co2abs - CD_flash_eoabs
 - CD_flash_flash - CD_flash_hex1 - CD_flash_hex2 - CD_flash_pump
 - CD_flash_reactor - CD_hex1_co2abs - CD_hex1_eoabs - CD_hex1_flash
 - CD_hex1_hex1 - CD_hex1_hex2 - CD_hex1_pump - CD_hex1_reactor
 - CD_hex2_co2abs - CD_hex2_eoabs - CD_hex2_flash - CD_hex2_hex1
 - CD_hex2_hex2 - CD_hex2_pump - CD_hex2_reactor - CD_pump_co2abs
 - CD_pump_eoabs - CD_pump_flash - CD_pump_hex1 - CD_pump_hex2 - CD_pump_pump
 - CD_pump_reactor - CD_reactor_co2abs - CD_reactor_eoabs - CD_reactor_flash
 - CD_reactor_hex1 - CD_reactor_hex2 - CD_reactor_pump - CD_reactor_reactor
 + SumCD = 0
_C45: A_reactor_hex1 - B_hex1_reactor = 0
_C450: - 0.933333333333 X_reactor + Y_reactor + 0.5 d_reactor
 + 0.466666666667 l_reactor <= 21
_C451: 0.933333333333 X_reactor + Y_reactor + 0.5 d_reactor
 + 0.466666666667 l_reactor <= 49
_C452: Y_reactor - 0.5 d_reactor >= 0
_C453: X_reactor + 0.5 l_reactor <= 30
_C454: X_reactor - 0.5 l_reactor >= 0
_C455: - 0.933333333333 X_hex1 + Y_hex1 + 0.5 d_hex1 + 0.466666666667 l_hex1
 <= 21
_C456: 0.933333333333 X_hex1 + Y_hex1 + 0.5 d_hex1 + 0.466666666667 l_hex1
 <= 49
_C457: Y_hex1 - 0.5 d_hex1 >= 0
_C458: X_hex1 + 0.5 l_hex1 <= 30
_C459: X_hex1 - 0.5 l_hex1 >= 0
_C46: - A_hex1_reactor + B_reactor_hex1 = 0
_C460: - 0.933333333333 X_eoabs + Y_eoabs + 0.5 d_eoabs
 + 0.466666666667 l_eoabs <= 21
_C461: 0.933333333333 X_eoabs + Y_eoabs + 0.5 d_eoabs + 0.466666666667 l_eoabs
 <= 49
_C462: Y_eoabs - 0.5 d_eoabs >= 0
_C463: X_eoabs + 0.5 l_eoabs <= 30
_C464: X_eoabs - 0.5 l_eoabs >= 0
_C465: - 0.933333333333 X_hex2 + Y_hex2 + 0.5 d_hex2 + 0.466666666667 l_hex2
 <= 21
_C466: 0.933333333333 X_hex2 + Y_hex2 + 0.5 d_hex2 + 0.466666666667 l_hex2
 <= 49
_C467: Y_hex2 - 0.5 d_hex2 >= 0
_C468: X_hex2 + 0.5 l_hex2 <= 30
_C469: X_hex2 - 0.5 l_hex2 >= 0
_C47: Wx_hex1_reactor + Wx_reactor_hex1 = 1
_C470: - 0.933333333333 X_co2abs + Y_co2abs + 0.5 d_co2abs
 + 0.466666666667 l_co2abs <= 21
_C471: 0.933333333333 X_co2abs + Y_co2abs + 0.5 d_co2abs
 + 0.466666666667 l_co2abs <= 49
_C472: Y_co2abs - 0.5 d_co2abs >= 0
_C473: X_co2abs + 0.5 l_co2abs <= 30
_C474: X_co2abs - 0.5 l_co2abs >= 0
_C475: - 0.933333333333 X_flash + Y_flash + 0.5 d_flash
 + 0.466666666667 l_flash <= 21
_C476: 0.933333333333 X_flash + Y_flash + 0.5 d_flash + 0.466666666667 l_flash
 <= 49
_C477: Y_flash - 0.5 d_flash >= 0
_C478: X_flash + 0.5 l_flash <= 30
_C479: X_flash - 0.5 l_flash >= 0
_C48: Wy_hex1_reactor + Wy_reactor_hex1 = 1
_C480: - 0.933333333333 X_pump + Y_pump + 0.5 d_pump + 0.466666666667 l_pump
 <= 21
_C481: 0.933333333333 X_pump + Y_pump + 0.5 d_pump + 0.466666666667 l_pump
 <= 49
_C482: Y_pump - 0.5 d_pump >= 0
_C483: X_pump + 0.5 l_pump <= 30
_C484: X_pump - 0.5 l_pump >= 0
_C485: Ve2_reactor_reactor = 0
_C486: - D_reactor_hex1 + Din_reactor_hex1 + Dout_reactor_hex1 = 0
_C487: Din_reactor_hex1 - 40 Psie_reactor_hex1 <= 0
_C488: Dout_reactor_hex1 + 40 Psie_reactor_hex1 >= 40
_C489: Dout_reactor_hex1 + 1000 Psie_reactor_hex1 <= 1000
_C49: - 56 D_reactor_eoabs + deltaP_reactor_eoabs = 0
_C490: 12775 Din_reactor_hex1 - 511000 Psie_reactor_hex1 + Ve2_reactor_hex1
 = 0
_C491: - D_reactor_eoabs + Din_reactor_eoabs + Dout_reactor_eoabs = 0
_C492: Din_reactor_eoabs - 40 Psie_reactor_eoabs <= 0
_C493: Dout_reactor_eoabs + 40 Psie_reactor_eoabs >= 40
_C494: Dout_reactor_eoabs + 1000 Psie_reactor_eoabs <= 1000
_C495: 77675 Din_reactor_eoabs - 3107000 Psie_reactor_eoabs
 + Ve2_reactor_eoabs = 0
_C496: - D_reactor_hex2 + Din_reactor_hex2 + Dout_reactor_hex2 = 0
_C497: Din_reactor_hex2 - 40 Psie_reactor_hex2 <= 0
_C498: Dout_reactor_hex2 + 40 Psie_reactor_hex2 >= 40
_C499: Dout_reactor_hex2 + 1000 Psie_reactor_hex2 <= 1000
_C5: 0 O_hex1 + l_hex1 = 11.42
_C50: CD_reactor_eoabs = 0
_C500: 12600 Din_reactor_hex2 - 504000 Psie_reactor_hex2 + Ve2_reactor_hex2
 = 0
_C501: - D_reactor_co2abs + Din_reactor_co2abs + Dout_reactor_co2abs = 0
_C502: Din_reactor_co2abs - 40 Psie_reactor_co2abs <= 0
_C503: Dout_reactor_co2abs + 40 Psie_reactor_co2abs >= 40
_C504: Dout_reactor_co2abs + 1000 Psie_reactor_co2abs <= 1000
_C505: 77032.5 Din_reactor_co2abs - 3081300 Psie_reactor_co2abs
 + Ve2_reactor_co2abs = 0
_C506: - D_reactor_flash + Din_reactor_flash + Dout_reactor_flash = 0
_C507: Din_reactor_flash - 40 Psie_reactor_flash <= 0
_C508: Dout_reactor_flash + 40 Psie_reactor_flash >= 40
_C509: Dout_reactor_flash + 1000 Psie_reactor_flash <= 1000
_C51: - L_reactor_eoabs + R_reactor_eoabs + X_eoabs - X_reactor = 0
_C510: 12625 Din_reactor_flash - 505000 Psie_reactor_flash + Ve2_reactor_flash
 = 0
_C511: - D_reactor_pump + Din_reactor_pump + Dout_reactor_pump = 0
_C512: Din_reactor_pump - 40 Psie_reactor_pump <= 0
_C513: Dout_reactor_pump + 40 Psie_reactor_pump >= 40
_C514: Dout_reactor_pump + 1000 Psie_reactor_pump <= 1000
_C515: 12537.5 Din_reactor_pump - 501500 Psie_reactor_pump + Ve2_reactor_pump
 = 0
_C516: - D_eoabs_reactor + Din_eoabs_reactor + Dout_eoabs_reactor = 0
_C517: Din_eoabs_reactor - 21.8 Psie_eoabs_reactor <= 0
_C518: Dout_eoabs_reactor + 21.8 Psie_eoabs_reactor >= 21.8
_C519: Dout_eoabs_reactor + 1000 Psie_eoabs_reactor <= 1000
_C52: A_reactor_eoabs - B_reactor_eoabs + Y_eoabs - Y_reactor = 0
_C520: 152981.651376 Din_eoabs_reactor - 3335000 Psie_eoabs_reactor
 + Ve2_eoabs_reactor = 0
_C521: - D_eoabs_hex1 + Din_eoabs_hex1 + Dout_eoabs_hex1 = 0
_C522: Din_eoabs_hex1 - 21.8 Psie_eoabs_hex1 <= 0
_C523: Dout_eoabs_hex1 + 21.8 Psie_eoabs_hex1 >= 21.8
_C524: Dout_eoabs_hex1 + 1000 Psie_eoabs_hex1 <= 1000
_C525: 23440.3669725 Din_eoabs_hex1 - 511000 Psie_eoabs_hex1 + Ve2_eoabs_hex1
 = 0
_C526: Ve2_eoabs_eoabs = 0
_C527: - D_eoabs_hex2 + Din_eoabs_hex2 + Dout_eoabs_hex2 = 0
_C528: Din_eoabs_hex2 - 21.8 Psie_eoabs_hex2 <= 0
_C529: Dout_eoabs_hex2 + 21.8 Psie_eoabs_hex2 >= 21.8
_C53: R_reactor_eoabs - 1000 Wx_reactor_eoabs <= 0
_C530: Dout_eoabs_hex2 + 1000 Psie_eoabs_hex2 <= 1000
_C531: 23119.266055 Din_eoabs_hex2 - 504000 Psie_eoabs_hex2 + Ve2_eoabs_hex2
 = 0
_C532: - D_eoabs_co2abs + Din_eoabs_co2abs + Dout_eoabs_co2abs = 0
_C533: Din_eoabs_co2abs - 21.8 Psie_eoabs_co2abs <= 0
_C534: Dout_eoabs_co2abs + 21.8 Psie_eoabs_co2abs >= 21.8
_C535: Dout_eoabs_co2abs + 1000 Psie_eoabs_co2abs <= 1000
_C536: 141344.036697 Din_eoabs_co2abs - 3081300 Psie_eoabs_co2abs
 + Ve2_eoabs_co2abs = 0
_C537: - D_eoabs_flash + Din_eoabs_flash + Dout_eoabs_flash = 0
_C538: Din_eoabs_flash - 21.8 Psie_eoabs_flash <= 0
_C539: Dout_eoabs_flash + 21.8 Psie_eoabs_flash >= 21.8
_C54: L_reactor_eoabs + 1000 Wx_reactor_eoabs <= 1000
_C540: Dout_eoabs_flash + 1000 Psie_eoabs_flash <= 1000
_C541: 23165.1376147 Din_eoabs_flash - 505000 Psie_eoabs_flash
 + Ve2_eoabs_flash = 0
_C542: - D_eoabs_pump + Din_eoabs_pump + Dout_eoabs_pump = 0
_C543: Din_eoabs_pump - 21.8 Psie_eoabs_pump <= 0
_C544: Dout_eoabs_pump + 21.8 Psie_eoabs_pump >= 21.8
_C545: Dout_eoabs_pump + 1000 Psie_eoabs_pump <= 1000
_C546: 23004.587156 Din_eoabs_pump - 501500 Psie_eoabs_pump + Ve2_eoabs_pump
 = 0
_C547: - D_co2abs_reactor + Din_co2abs_reactor + Dout_co2abs_reactor = 0
_C548: Din_co2abs_reactor - 18.06 Psie_co2abs_reactor <= 0
_C549: Dout_co2abs_reactor + 18.06 Psie_co2abs_reactor >= 18.06
_C55: A_reactor_eoabs - 1000 Wy_reactor_eoabs <= 0
_C550: Dout_co2abs_reactor + 1000 Psie_co2abs_reactor <= 1000
_C551: 184662.236988 Din_co2abs_reactor - 3335000 Psie_co2abs_reactor
 + Ve2_co2abs_reactor = 0
_C552: - D_co2abs_hex1 + Din_co2abs_hex1 + Dout_co2abs_hex1 = 0
_C553: Din_co2abs_hex1 - 18.06 Psie_co2abs_hex1 <= 0
_C554: Dout_co2abs_hex1 + 18.06 Psie_co2abs_hex1 >= 18.06
_C555: Dout_co2abs_hex1 + 1000 Psie_co2abs_hex1 <= 1000
_C556: 28294.5736434 Din_co2abs_hex1 - 511000 Psie_co2abs_hex1
 + Ve2_co2abs_hex1 = 0
_C557: - D_co2abs_eoabs + Din_co2abs_eoabs + Dout_co2abs_eoabs = 0
_C558: Din_co2abs_eoabs - 18.06 Psie_co2abs_eoabs <= 0
_C559: Dout_co2abs_eoabs + 18.06 Psie_co2abs_eoabs >= 18.06
_C56: B_reactor_eoabs + 1000 Wy_reactor_eoabs <= 1000
_C560: Dout_co2abs_eoabs + 1000 Psie_co2abs_eoabs <= 1000
_C561: 172037.65227 Din_co2abs_eoabs - 3107000 Psie_co2abs_eoabs
 + Ve2_co2abs_eoabs = 0
_C562: - D_co2abs_hex2 + Din_co2abs_hex2 + Dout_co2abs_hex2 = 0
_C563: Din_co2abs_hex2 - 18.06 Psie_co2abs_hex2 <= 0
_C564: Dout_co2abs_hex2 + 18.06 Psie_co2abs_hex2 >= 18.06
_C565: Dout_co2abs_hex2 + 1000 Psie_co2abs_hex2 <= 1000
_C566: 27906.9767442 Din_co2abs_hex2 - 504000 Psie_co2abs_hex2
 + Ve2_co2abs_hex2 = 0
_C567: Ve2_co2abs_co2abs = 0
_C568: - D_co2abs_flash + Din_co2abs_flash + Dout_co2abs_flash = 0
_C569: Din_co2abs_flash - 18.06 Psie_co2abs_flash <= 0
_C57: - A_reactor_eoabs - B_reactor_eoabs + D_reactor_eoabs - L_reactor_eoabs
 - R_reactor_eoabs = 0
_C570: Dout_co2abs_flash + 18.06 Psie_co2abs_flash >= 18.06
_C571: Dout_co2abs_flash + 1000 Psie_co2abs_flash <= 1000
_C572: 27962.3477298 Din_co2abs_flash - 505000 Psie_co2abs_flash
 + Ve2_co2abs_flash = 0
_C573: - D_co2abs_pump + Din_co2abs_pump + Dout_co2abs_pump = 0
_C574: Din_co2abs_pump - 18.06 Psie_co2abs_pump <= 0
_C575: Dout_co2abs_pump + 18.06 Psie_co2abs_pump >= 18.06
_C576: Dout_co2abs_pump + 1000 Psie_co2abs_pump <= 1000
_C577: 27768.5492802 Din_co2abs_pump - 501500 Psie_co2abs_pump
 + Ve2_co2abs_pump = 0
_C578: - Ve2_reactor_co2abs - Ve2_reactor_eoabs - Ve2_reactor_flash
 - Ve2_reactor_hex1 - Ve2_reactor_hex2 - Ve2_reactor_pump
 - Ve2_reactor_reactor + Ve_reactor = 3335000
_C579: Omega0_reactor - 0.87 Ve_reactor = 0
_C58: - D_eoabs_reactor + D_reactor_eoabs = 0
_C580: - Ve2_eoabs_co2abs - Ve2_eoabs_eoabs - Ve2_eoabs_flash - Ve2_eoabs_hex1
 - Ve2_eoabs_hex2 - Ve2_eoabs_pump - Ve2_eoabs_reactor + Ve_eoabs = 3107000
_C581: Omega0_eoabs - 0.73 Ve_eoabs = 0
_C582: - Ve2_co2abs_co2abs - Ve2_co2abs_eoabs - Ve2_co2abs_flash
 - Ve2_co2abs_hex1 - Ve2_co2abs_hex2 - Ve2_co2abs_pump - Ve2_co2abs_reactor
 + Ve_co2abs = 3081300
_C583: Omega0_co2abs - 0.66 Ve_co2abs = 0
_C584: - 100000000 Z_reactor_1 + linOmega0_reactor_1 <= 0
_C585: - 100000000 Z_reactor_2 + linOmega0_reactor_2 <= 0
_C586: - 100000000 Z_reactor_3 + linOmega0_reactor_3 <= 0
_C587: - 100000000 Z_reactor_4 + linOmega0_reactor_4 <= 0
_C588: - 100000000 Z_reactor_5 + linOmega0_reactor_5 <= 0
_C589: - 100000000 Z_reactor_6 + linOmega0_reactor_6 <= 0
_C59: 1000 E1_reactor_eoabs + 1000 E2_reactor_eoabs - X_eoabs + X_reactor
 - 0.5 l_eoabs - 0.5 l_reactor >= 0
_C590: - 100000000 Z_eoabs_1 + linOmega0_eoabs_1 <= 0
_C591: - 100000000 Z_eoabs_2 + linOmega0_eoabs_2 <= 0
_C592: - 100000000 Z_eoabs_3 + linOmega0_eoabs_3 <= 0
_C593: - 100000000 Z_eoabs_4 + linOmega0_eoabs_4 <= 0
_C594: - 100000000 Z_eoabs_5 + linOmega0_eoabs_5 <= 0
_C595: - 100000000 Z_eoabs_6 + linOmega0_eoabs_6 <= 0
_C596: - 100000000 Z_co2abs_1 + linOmega0_co2abs_1 <= 0
_C597: - 100000000 Z_co2abs_2 + linOmega0_co2abs_2 <= 0
_C598: - 100000000 Z_co2abs_3 + linOmega0_co2abs_3 <= 0
_C599: - 100000000 Z_co2abs_4 + linOmega0_co2abs_4 <= 0
_C6: d_hex1 + l_hex1 = 22.84
_C60: - 1000 E1_reactor_eoabs + 1000 E2_reactor_eoabs + X_eoabs - X_reactor
 - 0.5 l_eoabs - 0.5 l_reactor >= -1000
_C600: - 100000000 Z_co2abs_5 + linOmega0_co2abs_5 <= 0
_C601: - 100000000 Z_co2abs_6 + linOmega0_co2abs_6 <= 0
_C602: Z_reactor_1 + Z_reactor_2 + Z_reactor_3 + Z_reactor_4 + Z_reactor_5
 + Z_reactor_6 = 1
_C603: Omega0_reactor - linOmega0_reactor_1 - linOmega0_reactor_2
 - linOmega0_reactor_3 - linOmega0_reactor_4 - linOmega0_reactor_5
 - linOmega0_reactor_6 = 0
_C604: Omega_reactor - linOmega0_reactor_1 - 0.9 linOmega0_reactor_2
 - 0.75 linOmega0_reactor_3 - 0.365 linOmega0_reactor_4
 - 0.292 linOmega0_reactor_5 - 0.117 linOmega0_reactor_6 = 0
_C605: - Omega_co2abs - Omega_eoabs - Omega_reactor + SumOmega = 0
_C606: SumPZ - 5000 Z_co2abs_2 - 20000 Z_co2abs_3 - 25000 Z_co2abs_4
 - 35000 Z_co2abs_5 - 55000 Z_co2abs_6 - 5000 Z_eoabs_2 - 20000 Z_eoabs_3
 - 25000 Z_eoabs_4 - 35000 Z_eoabs_5 - 55000 Z_eoabs_6 - 5000 Z_reactor_2
 - 15000 Z_reactor_3 - 40000 Z_reactor_4 - 60000 Z_reactor_5
 - 125000 Z_reactor_6 = 0
_C607: Z_eoabs_1 + Z_eoabs_2 + Z_eoabs_3 + Z_eoabs_4 + Z_eoabs_5 + Z_eoabs_6
 = 1
_C608: Omega0_eoabs - linOmega0_eoabs_1 - linOmega0_eoabs_2
 - linOmega0_eoabs_3 - linOmega0_eoabs_4 - linOmega0_eoabs_5
 - linOmega0_eoabs_6 = 0
_C609: Omega_eoabs - linOmega0_eoabs_1 - 0.9 linOmega0_eoabs_2
 - 0.76 linOmega0_eoabs_3 - 0.684 linOmega0_eoabs_4 - 0.612 linOmega0_eoabs_5
 - 0.465 linOmega0_eoabs_6 = 0
_C61: 1000 E1_reactor_eoabs - 1000 E2_reactor_eoabs - Y_eoabs + Y_reactor
 - 0.5 d_eoabs - 0.5 d_reactor >= -1000
_C610: - Omega_co2abs - Omega_eoabs - Omega_reactor + SumOmega = 0
_C611: SumPZ - 5000 Z_co2abs_2 - 20000 Z_co2abs_3 - 25000 Z_co2abs_4
 - 35000 Z_co2abs_5 - 55000 Z_co2abs_6 - 5000 Z_eoabs_2 - 20000 Z_eoabs_3
 - 25000 Z_eoabs_4 - 35000 Z_eoabs_5 - 55000 Z_eoabs_6 - 5000 Z_reactor_2
 - 15000 Z_reactor_3 - 40000 Z_reactor_4 - 60000 Z_reactor_5
 - 125000 Z_reactor_6 = 0
_C612: Z_co2abs_1 + Z_co2abs_2 + Z_co2abs_3 + Z_co2abs_4 + Z_co2abs_5
 + Z_co2abs_6 = 1
_C613: Omega0_co2abs - linOmega0_co2abs_1 - linOmega0_co2abs_2
 - linOmega0_co2abs_3 - linOmega0_co2abs_4 - linOmega0_co2abs_5
 - linOmega0_co2abs_6 = 0
_C614: Omega_co2abs - linOmega0_co2abs_1 - 0.9 linOmega0_co2abs_2
 - 0.76 linOmega0_co2abs_3 - 0.684 linOmega0_co2abs_4
 - 0.612 linOmega0_co2abs_5 - 0.465 linOmega0_co2abs_6 = 0
_C615: - Omega_co2abs - Omega_eoabs - Omega_reactor + SumOmega = 0
_C616: SumPZ - 5000 Z_co2abs_2 - 20000 Z_co2abs_3 - 25000 Z_co2abs_4
 - 35000 Z_co2abs_5 - 55000 Z_co2abs_6 - 5000 Z_eoabs_2 - 20000 Z_eoabs_3
 - 25000 Z_eoabs_4 - 35000 Z_eoabs_5 - 55000 Z_eoabs_6 - 5000 Z_reactor_2
 - 15000 Z_reactor_3 - 40000 Z_reactor_4 - 60000 Z_reactor_5
 - 125000 Z_reactor_6 = 0
_C62: - 1000 E1_reactor_eoabs - 1000 E2_reactor_eoabs + Y_eoabs - Y_reactor
 - 0.5 d_eoabs - 0.5 d_reactor >= -2000
_C63: L_reactor_eoabs - R_eoabs_reactor = 0
_C64: - L_eoabs_reactor + R_reactor_eoabs = 0
_C65: A_reactor_eoabs - B_eoabs_reactor = 0
_C66: - A_eoabs_reactor + B_reactor_eoabs = 0
_C67: Wx_eoabs_reactor + Wx_reactor_eoabs = 1
_C68: Wy_eoabs_reactor + Wy_reactor_eoabs = 1
_C69: - 56 D_hex1_eoabs + deltaP_hex1_eoabs = 0
_C7: X_hex1 - 0.5 l_hex1 >= 0
_C70: CD_hex1_eoabs - 27.3003837522 D_hex1_eoabs = 0
_C71: - L_hex1_eoabs + R_hex1_eoabs + X_eoabs - X_hex1 = 0
_C72: A_hex1_eoabs - B_hex1_eoabs + Y_eoabs - Y_hex1 = 0
_C73: R_hex1_eoabs - 1000 Wx_hex1_eoabs <= 0
_C74: L_hex1_eoabs + 1000 Wx_hex1_eoabs <= 1000
_C75: A_hex1_eoabs - 1000 Wy_hex1_eoabs <= 0
_C76: B_hex1_eoabs + 1000 Wy_hex1_eoabs <= 1000
_C77: - A_hex1_eoabs - B_hex1_eoabs + D_hex1_eoabs - L_hex1_eoabs
 - R_hex1_eoabs = 0
_C78: - D_eoabs_hex1 + D_hex1_eoabs = 0
_C79: 1000 E1_hex1_eoabs + 1000 E2_hex1_eoabs - X_eoabs + X_hex1 - 0.5 l_eoabs
 - 0.5 l_hex1 >= 0
_C8: Y_hex1 - 0.5 d_hex1 >= 0
_C80: - 1000 E1_hex1_eoabs + 1000 E2_hex1_eoabs + X_eoabs - X_hex1
 - 0.5 l_eoabs - 0.5 l_hex1 >= -1000
_C81: 1000 E1_hex1_eoabs - 1000 E2_hex1_eoabs - Y_eoabs + Y_hex1 - 0.5 d_eoabs
 - 0.5 d_hex1 >= -1000
_C82: - 1000 E1_hex1_eoabs - 1000 E2_hex1_eoabs + Y_eoabs - Y_hex1
 - 0.5 d_eoabs - 0.5 d_hex1 >= -2000
_C83: L_hex1_eoabs - R_eoabs_hex1 = 0
_C84: - L_eoabs_hex1 + R_hex1_eoabs = 0
_C85: A_hex1_eoabs - B_eoabs_hex1 = 0
_C86: - A_eoabs_hex1 + B_hex1_eoabs = 0
_C87: Wx_eoabs_hex1 + Wx_hex1_eoabs = 1
_C88: Wy_eoabs_hex1 + Wy_hex1_eoabs = 1
_C89: - 56 D_reactor_hex2 + deltaP_reactor_hex2 = 0
_C9: 0 O_eoabs + l_eoabs = 7.68
_C90: CD_reactor_hex2 = 0
_C91: - L_reactor_hex2 + R_reactor_hex2 + X_hex2 - X_reactor = 0
_C92: A_reactor_hex2 - B_reactor_hex2 + Y_hex2 - Y_reactor = 0
_C93: R_reactor_hex2 - 1000 Wx_reactor_hex2 <= 0
_C94: L_reactor_hex2 + 1000 Wx_reactor_hex2 <= 1000
_C95: A_reactor_hex2 - 1000 Wy_reactor_hex2 <= 0
_C96: B_reactor_hex2 + 1000 Wy_reactor_hex2 <= 1000
_C97: - A_reactor_hex2 - B_reactor_hex2 + D_reactor_hex2 - L_reactor_hex2
 - R_reactor_hex2 = 0
_C98: - D_hex2_reactor + D_reactor_hex2 = 0
_C99: 1000 E1_reactor_hex2 + 1000 E2_reactor_hex2 - X_hex2 + X_reactor
 - 0.5 l_hex2 - 0.5 l_reactor >= 0
Bounds
 X_co2abs free
 X_eoabs free
 X_flash free
 X_hex1 free
 X_hex2 free
 X_pump free
 X_reactor free
 Y_co2abs free
 Y_eoabs free
 Y_flash free
 Y_hex1 free
 Y_hex2 free
 Y_pump free
 Y_reactor free
Binaries
E1_co2abs_flash
E1_co2abs_pump
E1_eoabs_co2abs
E1_eoabs_flash
E1_eoabs_hex2
E1_eoabs_pump
E1_flash_pump
E1_hex1_co2abs
E1_hex1_eoabs
E1_hex1_flash
E1_hex1_hex2
E1_hex1_pump
E1_hex2_co2abs
E1_hex2_flash
E1_hex2_pump
E1_reactor_co2abs
E1_reactor_eoabs
E1_reactor_flash
E1_reactor_hex1
E1_reactor_hex2
E1_reactor_pump
E2_co2abs_flash
E2_co2abs_pump
E2_eoabs_co2abs
E2_eoabs_flash
E2_eoabs_hex2
E2_eoabs_pump
E2_flash_pump
E2_hex1_co2abs
E2_hex1_eoabs
E2_hex1_flash
E2_hex1_hex2
E2_hex1_pump
E2_hex2_co2abs
E2_hex2_flash
E2_hex2_pump
E2_reactor_co2abs
E2_reactor_eoabs
E2_reactor_flash
E2_reactor_hex1
E2_reactor_hex2
E2_reactor_pump
O_co2abs
O_eoabs
O_flash
O_hex1
O_hex2
O_pump
O_reactor
Psie_co2abs_eoabs
Psie_co2abs_flash
Psie_co2abs_hex1
Psie_co2abs_hex2
Psie_co2abs_pump
Psie_co2abs_reactor
Psie_eoabs_co2abs
Psie_eoabs_flash
Psie_eoabs_hex1
Psie_eoabs_hex2
Psie_eoabs_pump
Psie_eoabs_reactor
Psie_reactor_co2abs
Psie_reactor_eoabs
Psie_reactor_flash
Psie_reactor_hex1
Psie_reactor_hex2
Psie_reactor_pump
Wx_co2abs_eoabs
Wx_co2abs_flash
Wx_co2abs_hex1
Wx_co2abs_hex2
Wx_co2abs_pump
Wx_co2abs_reactor
Wx_eoabs_co2abs
Wx_eoabs_flash
Wx_eoabs_hex1
Wx_eoabs_hex2
Wx_eoabs_pump
Wx_eoabs_reactor
Wx_flash_co2abs
Wx_flash_eoabs
Wx_flash_hex1
Wx_flash_hex2
Wx_flash_pump
Wx_flash_reactor
Wx_hex1_co2abs
Wx_hex1_eoabs
Wx_hex1_flash
Wx_hex1_hex2
Wx_hex1_pump
Wx_hex1_reactor
Wx_hex2_co2abs
Wx_hex2_eoabs
Wx_hex2_flash
Wx_hex2_hex1
Wx_hex2_pump
Wx_hex2_reactor
Wx_pump_co2abs
Wx_pump_eoabs
Wx_pump_flash
Wx_pump_hex1
Wx_pump_hex2
Wx_pump_reactor
Wx_reactor_co2abs
Wx_reactor_eoabs
Wx_reactor_flash
Wx_reactor_hex1
Wx_reactor_hex2
Wx_reactor_pump
Wy_co2abs_eoabs
Wy_co2abs_flash
Wy_co2abs_hex1
Wy_co2abs_hex2
Wy_co2abs_pump
Wy_co2abs_reactor
Wy_eoabs_co2abs
Wy_eoabs_flash
Wy_eoabs_hex1
Wy_eoabs_hex2
Wy_eoabs_pump
Wy_eoabs_reactor
Wy_flash_co2abs
Wy_flash_eoabs
Wy_flash_hex1
Wy_flash_hex2
Wy_flash_pump
Wy_flash_reactor
Wy_hex1_co2abs
Wy_hex1_eoabs
Wy_hex1_flash
Wy_hex1_hex2
Wy_hex1_pump
Wy_hex1_reactor
Wy_hex2_co2abs
Wy_hex2_eoabs
Wy_hex2_flash
Wy_hex2_hex1
Wy_hex2_pump
Wy_hex2_reactor
Wy_pump_co2abs
Wy_pump_eoabs
Wy_pump_flash
Wy_pump_hex1
Wy_pump_hex2
Wy_pump_reactor
Wy_reactor_co2abs
Wy_reactor_eoabs
Wy_reactor_flash
Wy_reactor_hex1
Wy_reactor_hex2
Wy_reactor_pump
Z_co2abs_1
Z_co2abs_2
Z_co2abs_3
Z_co2abs_4
Z_co2abs_5
Z_co2abs_6
Z_eoabs_1
Z_eoabs_2
Z_eoabs_3
Z_eoabs_4
Z_eoabs_5
Z_eoabs_6
Z_reactor_1
Z_reactor_2
Z_reactor_3
Z_reactor_4
Z_reactor_5
Z_reactor_6
End
