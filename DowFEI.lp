\* Layout_Problem_Model *\
Minimize
OBJ: SumCD + SumOmega + SumPZ
Subject To
_C1: - 5 O_reactor + l_reactor = 5
_C10: d_eoabs + l_eoabs = 15.36
_C100: - L_hex2_reactor + R_reactor_hex2 = 0
_C101: A_reactor_hex2 - B_hex2_reactor = 0
_C102: - A_hex2_reactor + B_reactor_hex2 = 0
_C103: Wx_hex2_reactor + Wx_reactor_hex2 = 1
_C104: Wy_hex2_reactor + Wy_reactor_hex2 = 1
_C105: CD_hex1_hex2 = 0
_C106: - L_hex1_hex2 + R_hex1_hex2 - X_hex1 + X_hex2 = 0
_C107: A_hex1_hex2 - B_hex1_hex2 - Y_hex1 + Y_hex2 = 0
_C108: R_hex1_hex2 - 1000 Wx_hex1_hex2 <= 0
_C109: L_hex1_hex2 + 1000 Wx_hex1_hex2 <= 1000
_C11: X_eoabs - 0.5 l_eoabs >= 0
_C110: A_hex1_hex2 - 1000 Wy_hex1_hex2 <= 0
_C111: B_hex1_hex2 + 1000 Wy_hex1_hex2 <= 1000
_C112: - A_hex1_hex2 - B_hex1_hex2 + D_hex1_hex2 - L_hex1_hex2 - R_hex1_hex2
 = 0
_C113: D_hex1_hex2 - D_hex2_hex1 = 0
_C114: 1000 E1_hex1_hex2 + 1000 E2_hex1_hex2 + X_hex1 - X_hex2 - 0.5 l_hex1
 - 0.5 l_hex2 >= 0
_C115: - 1000 E1_hex1_hex2 + 1000 E2_hex1_hex2 - X_hex1 + X_hex2 - 0.5 l_hex1
 - 0.5 l_hex2 >= -1000
_C116: 1000 E1_hex1_hex2 - 1000 E2_hex1_hex2 + Y_hex1 - Y_hex2 - 0.5 d_hex1
 - 0.5 d_hex2 >= -1000
_C117: - 1000 E1_hex1_hex2 - 1000 E2_hex1_hex2 - Y_hex1 + Y_hex2 - 0.5 d_hex1
 - 0.5 d_hex2 >= -2000
_C118: L_hex1_hex2 - R_hex2_hex1 = 0
_C119: - L_hex2_hex1 + R_hex1_hex2 = 0
_C12: Y_eoabs - 0.5 d_eoabs >= 0
_C120: A_hex1_hex2 - B_hex2_hex1 = 0
_C121: - A_hex2_hex1 + B_hex1_hex2 = 0
_C122: Wx_hex1_hex2 + Wx_hex2_hex1 = 1
_C123: Wy_hex1_hex2 + Wy_hex2_hex1 = 1
_C124: CD_eoabs_hex2 - 2730.03837522 D_eoabs_hex2 = 0
_C125: - L_eoabs_hex2 + R_eoabs_hex2 - X_eoabs + X_hex2 = 0
_C126: A_eoabs_hex2 - B_eoabs_hex2 - Y_eoabs + Y_hex2 = 0
_C127: R_eoabs_hex2 - 1000 Wx_eoabs_hex2 <= 0
_C128: L_eoabs_hex2 + 1000 Wx_eoabs_hex2 <= 1000
_C129: A_eoabs_hex2 - 1000 Wy_eoabs_hex2 <= 0
_C13: 0 O_hex2 + l_hex2 = 8.48
_C130: B_eoabs_hex2 + 1000 Wy_eoabs_hex2 <= 1000
_C131: - A_eoabs_hex2 - B_eoabs_hex2 + D_eoabs_hex2 - L_eoabs_hex2
 - R_eoabs_hex2 = 0
_C132: D_eoabs_hex2 - D_hex2_eoabs = 0
_C133: 1000 E1_eoabs_hex2 + 1000 E2_eoabs_hex2 + X_eoabs - X_hex2
 - 0.5 l_eoabs - 0.5 l_hex2 >= 0
_C134: - 1000 E1_eoabs_hex2 + 1000 E2_eoabs_hex2 - X_eoabs + X_hex2
 - 0.5 l_eoabs - 0.5 l_hex2 >= -1000
_C135: 1000 E1_eoabs_hex2 - 1000 E2_eoabs_hex2 + Y_eoabs - Y_hex2
 - 0.5 d_eoabs - 0.5 d_hex2 >= -1000
_C136: - 1000 E1_eoabs_hex2 - 1000 E2_eoabs_hex2 - Y_eoabs + Y_hex2
 - 0.5 d_eoabs - 0.5 d_hex2 >= -2000
_C137: L_eoabs_hex2 - R_hex2_eoabs = 0
_C138: - L_hex2_eoabs + R_eoabs_hex2 = 0
_C139: A_eoabs_hex2 - B_hex2_eoabs = 0
_C14: d_hex2 + l_hex2 = 16.96
_C140: - A_hex2_eoabs + B_eoabs_hex2 = 0
_C141: Wx_eoabs_hex2 + Wx_hex2_eoabs = 1
_C142: Wy_eoabs_hex2 + Wy_hex2_eoabs = 1
_C143: CD_reactor_co2abs - 2730.03837522 D_reactor_co2abs = 0
_C144: - L_reactor_co2abs + R_reactor_co2abs + X_co2abs - X_reactor = 0
_C145: A_reactor_co2abs - B_reactor_co2abs + Y_co2abs - Y_reactor = 0
_C146: R_reactor_co2abs - 1000 Wx_reactor_co2abs <= 0
_C147: L_reactor_co2abs + 1000 Wx_reactor_co2abs <= 1000
_C148: A_reactor_co2abs - 1000 Wy_reactor_co2abs <= 0
_C149: B_reactor_co2abs + 1000 Wy_reactor_co2abs <= 1000
_C15: X_hex2 - 0.5 l_hex2 >= 0
_C150: - A_reactor_co2abs - B_reactor_co2abs + D_reactor_co2abs
 - L_reactor_co2abs - R_reactor_co2abs = 0
_C151: - D_co2abs_reactor + D_reactor_co2abs = 0
_C152: 1000 E1_reactor_co2abs + 1000 E2_reactor_co2abs - X_co2abs + X_reactor
 - 0.5 l_co2abs - 0.5 l_reactor >= 0
_C153: - 1000 E1_reactor_co2abs + 1000 E2_reactor_co2abs + X_co2abs
 - X_reactor - 0.5 l_co2abs - 0.5 l_reactor >= -1000
_C154: 1000 E1_reactor_co2abs - 1000 E2_reactor_co2abs - Y_co2abs + Y_reactor
 - 0.5 d_co2abs - 0.5 d_reactor >= -1000
_C155: - 1000 E1_reactor_co2abs - 1000 E2_reactor_co2abs + Y_co2abs
 - Y_reactor - 0.5 d_co2abs - 0.5 d_reactor >= -2000
_C156: L_reactor_co2abs - R_co2abs_reactor = 0
_C157: - L_co2abs_reactor + R_reactor_co2abs = 0
_C158: A_reactor_co2abs - B_co2abs_reactor = 0
_C159: - A_co2abs_reactor + B_reactor_co2abs = 0
_C16: Y_hex2 - 0.5 d_hex2 >= 0
_C160: Wx_co2abs_reactor + Wx_reactor_co2abs = 1
_C161: Wy_co2abs_reactor + Wy_reactor_co2abs = 1
_C162: CD_hex1_co2abs = 0
_C163: - L_hex1_co2abs + R_hex1_co2abs + X_co2abs - X_hex1 = 0
_C164: A_hex1_co2abs - B_hex1_co2abs + Y_co2abs - Y_hex1 = 0
_C165: R_hex1_co2abs - 1000 Wx_hex1_co2abs <= 0
_C166: L_hex1_co2abs + 1000 Wx_hex1_co2abs <= 1000
_C167: A_hex1_co2abs - 1000 Wy_hex1_co2abs <= 0
_C168: B_hex1_co2abs + 1000 Wy_hex1_co2abs <= 1000
_C169: - A_hex1_co2abs - B_hex1_co2abs + D_hex1_co2abs - L_hex1_co2abs
 - R_hex1_co2abs = 0
_C17: 0 O_co2abs + l_co2abs = 5
_C170: - D_co2abs_hex1 + D_hex1_co2abs = 0
_C171: 1000 E1_hex1_co2abs + 1000 E2_hex1_co2abs - X_co2abs + X_hex1
 - 0.5 l_co2abs - 0.5 l_hex1 >= 0
_C172: - 1000 E1_hex1_co2abs + 1000 E2_hex1_co2abs + X_co2abs - X_hex1
 - 0.5 l_co2abs - 0.5 l_hex1 >= -1000
_C173: 1000 E1_hex1_co2abs - 1000 E2_hex1_co2abs - Y_co2abs + Y_hex1
 - 0.5 d_co2abs - 0.5 d_hex1 >= -1000
_C174: - 1000 E1_hex1_co2abs - 1000 E2_hex1_co2abs + Y_co2abs - Y_hex1
 - 0.5 d_co2abs - 0.5 d_hex1 >= -2000
_C175: L_hex1_co2abs - R_co2abs_hex1 = 0
_C176: - L_co2abs_hex1 + R_hex1_co2abs = 0
_C177: A_hex1_co2abs - B_co2abs_hex1 = 0
_C178: - A_co2abs_hex1 + B_hex1_co2abs = 0
_C179: Wx_co2abs_hex1 + Wx_hex1_co2abs = 1
_C18: d_co2abs + l_co2abs = 10
_C180: Wy_co2abs_hex1 + Wy_hex1_co2abs = 1
_C181: CD_eoabs_co2abs = 0
_C182: - L_eoabs_co2abs + R_eoabs_co2abs + X_co2abs - X_eoabs = 0
_C183: A_eoabs_co2abs - B_eoabs_co2abs + Y_co2abs - Y_eoabs = 0
_C184: R_eoabs_co2abs - 1000 Wx_eoabs_co2abs <= 0
_C185: L_eoabs_co2abs + 1000 Wx_eoabs_co2abs <= 1000
_C186: A_eoabs_co2abs - 1000 Wy_eoabs_co2abs <= 0
_C187: B_eoabs_co2abs + 1000 Wy_eoabs_co2abs <= 1000
_C188: - A_eoabs_co2abs - B_eoabs_co2abs + D_eoabs_co2abs - L_eoabs_co2abs
 - R_eoabs_co2abs = 0
_C189: - D_co2abs_eoabs + D_eoabs_co2abs = 0
_C19: X_co2abs - 0.5 l_co2abs >= 0
_C190: 1000 E1_eoabs_co2abs + 1000 E2_eoabs_co2abs - X_co2abs + X_eoabs
 - 0.5 l_co2abs - 0.5 l_eoabs >= 0
_C191: - 1000 E1_eoabs_co2abs + 1000 E2_eoabs_co2abs + X_co2abs - X_eoabs
 - 0.5 l_co2abs - 0.5 l_eoabs >= -1000
_C192: 1000 E1_eoabs_co2abs - 1000 E2_eoabs_co2abs - Y_co2abs + Y_eoabs
 - 0.5 d_co2abs - 0.5 d_eoabs >= -1000
_C193: - 1000 E1_eoabs_co2abs - 1000 E2_eoabs_co2abs + Y_co2abs - Y_eoabs
 - 0.5 d_co2abs - 0.5 d_eoabs >= -2000
_C194: L_eoabs_co2abs - R_co2abs_eoabs = 0
_C195: - L_co2abs_eoabs + R_eoabs_co2abs = 0
_C196: A_eoabs_co2abs - B_co2abs_eoabs = 0
_C197: - A_co2abs_eoabs + B_eoabs_co2abs = 0
_C198: Wx_co2abs_eoabs + Wx_eoabs_co2abs = 1
_C199: Wy_co2abs_eoabs + Wy_eoabs_co2abs = 1
_C2: d_reactor + l_reactor = 15
_C20: Y_co2abs - 0.5 d_co2abs >= 0
_C200: CD_hex2_co2abs - 2730.03837522 D_hex2_co2abs = 0
_C201: - L_hex2_co2abs + R_hex2_co2abs + X_co2abs - X_hex2 = 0
_C202: A_hex2_co2abs - B_hex2_co2abs + Y_co2abs - Y_hex2 = 0
_C203: R_hex2_co2abs - 1000 Wx_hex2_co2abs <= 0
_C204: L_hex2_co2abs + 1000 Wx_hex2_co2abs <= 1000
_C205: A_hex2_co2abs - 1000 Wy_hex2_co2abs <= 0
_C206: B_hex2_co2abs + 1000 Wy_hex2_co2abs <= 1000
_C207: - A_hex2_co2abs - B_hex2_co2abs + D_hex2_co2abs - L_hex2_co2abs
 - R_hex2_co2abs = 0
_C208: - D_co2abs_hex2 + D_hex2_co2abs = 0
_C209: 1000 E1_hex2_co2abs + 1000 E2_hex2_co2abs - X_co2abs + X_hex2
 - 0.5 l_co2abs - 0.5 l_hex2 >= 0
_C21: 0 O_flash + l_flash = 2.68
_C210: - 1000 E1_hex2_co2abs + 1000 E2_hex2_co2abs + X_co2abs - X_hex2
 - 0.5 l_co2abs - 0.5 l_hex2 >= -1000
_C211: 1000 E1_hex2_co2abs - 1000 E2_hex2_co2abs - Y_co2abs + Y_hex2
 - 0.5 d_co2abs - 0.5 d_hex2 >= -1000
_C212: - 1000 E1_hex2_co2abs - 1000 E2_hex2_co2abs + Y_co2abs - Y_hex2
 - 0.5 d_co2abs - 0.5 d_hex2 >= -2000
_C213: L_hex2_co2abs - R_co2abs_hex2 = 0
_C214: - L_co2abs_hex2 + R_hex2_co2abs = 0
_C215: A_hex2_co2abs - B_co2abs_hex2 = 0
_C216: - A_co2abs_hex2 + B_hex2_co2abs = 0
_C217: Wx_co2abs_hex2 + Wx_hex2_co2abs = 1
_C218: Wy_co2abs_hex2 + Wy_hex2_co2abs = 1
_C219: CD_reactor_flash = 0
_C22: d_flash + l_flash = 5.36
_C220: - L_reactor_flash + R_reactor_flash + X_flash - X_reactor = 0
_C221: A_reactor_flash - B_reactor_flash + Y_flash - Y_reactor = 0
_C222: R_reactor_flash - 1000 Wx_reactor_flash <= 0
_C223: L_reactor_flash + 1000 Wx_reactor_flash <= 1000
_C224: A_reactor_flash - 1000 Wy_reactor_flash <= 0
_C225: B_reactor_flash + 1000 Wy_reactor_flash <= 1000
_C226: - A_reactor_flash - B_reactor_flash + D_reactor_flash - L_reactor_flash
 - R_reactor_flash = 0
_C227: - D_flash_reactor + D_reactor_flash = 0
_C228: 1000 E1_reactor_flash + 1000 E2_reactor_flash - X_flash + X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= 0
_C229: - 1000 E1_reactor_flash + 1000 E2_reactor_flash + X_flash - X_reactor
 - 0.5 l_flash - 0.5 l_reactor >= -1000
_C23: X_flash - 0.5 l_flash >= 0
_C230: 1000 E1_reactor_flash - 1000 E2_reactor_flash - Y_flash + Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -1000
_C231: - 1000 E1_reactor_flash - 1000 E2_reactor_flash + Y_flash - Y_reactor
 - 0.5 d_flash - 0.5 d_reactor >= -2000
_C232: L_reactor_flash - R_flash_reactor = 0
_C233: - L_flash_reactor + R_reactor_flash = 0
_C234: A_reactor_flash - B_flash_reactor = 0
_C235: - A_flash_reactor + B_reactor_flash = 0
_C236: Wx_flash_reactor + Wx_reactor_flash = 1
_C237: Wy_flash_reactor + Wy_reactor_flash = 1
_C238: CD_hex1_flash = 0
_C239: - L_hex1_flash + R_hex1_flash + X_flash - X_hex1 = 0
_C24: Y_flash - 0.5 d_flash >= 0
_C240: A_hex1_flash - B_hex1_flash + Y_flash - Y_hex1 = 0
_C241: R_hex1_flash - 1000 Wx_hex1_flash <= 0
_C242: L_hex1_flash + 1000 Wx_hex1_flash <= 1000
_C243: A_hex1_flash - 1000 Wy_hex1_flash <= 0
_C244: B_hex1_flash + 1000 Wy_hex1_flash <= 1000
_C245: - A_hex1_flash - B_hex1_flash + D_hex1_flash - L_hex1_flash
 - R_hex1_flash = 0
_C246: - D_flash_hex1 + D_hex1_flash = 0
_C247: 1000 E1_hex1_flash + 1000 E2_hex1_flash - X_flash + X_hex1
 - 0.5 l_flash - 0.5 l_hex1 >= 0
_C248: - 1000 E1_hex1_flash + 1000 E2_hex1_flash + X_flash - X_hex1
 - 0.5 l_flash - 0.5 l_hex1 >= -1000
_C249: 1000 E1_hex1_flash - 1000 E2_hex1_flash - Y_flash + Y_hex1
 - 0.5 d_flash - 0.5 d_hex1 >= -1000
_C25: 0 O_pump + l_pump = 2.4
_C250: - 1000 E1_hex1_flash - 1000 E2_hex1_flash + Y_flash - Y_hex1
 - 0.5 d_flash - 0.5 d_hex1 >= -2000
_C251: L_hex1_flash - R_flash_hex1 = 0
_C252: - L_flash_hex1 + R_hex1_flash = 0
_C253: A_hex1_flash - B_flash_hex1 = 0
_C254: - A_flash_hex1 + B_hex1_flash = 0
_C255: Wx_flash_hex1 + Wx_hex1_flash = 1
_C256: Wy_flash_hex1 + Wy_hex1_flash = 1
_C257: CD_eoabs_flash = 0
_C258: - L_eoabs_flash + R_eoabs_flash - X_eoabs + X_flash = 0
_C259: A_eoabs_flash - B_eoabs_flash - Y_eoabs + Y_flash = 0
_C26: d_pump + l_pump = 4.8
_C260: R_eoabs_flash - 1000 Wx_eoabs_flash <= 0
_C261: L_eoabs_flash + 1000 Wx_eoabs_flash <= 1000
_C262: A_eoabs_flash - 1000 Wy_eoabs_flash <= 0
_C263: B_eoabs_flash + 1000 Wy_eoabs_flash <= 1000
_C264: - A_eoabs_flash - B_eoabs_flash + D_eoabs_flash - L_eoabs_flash
 - R_eoabs_flash = 0
_C265: D_eoabs_flash - D_flash_eoabs = 0
_C266: 1000 E1_eoabs_flash + 1000 E2_eoabs_flash + X_eoabs - X_flash
 - 0.5 l_eoabs - 0.5 l_flash >= 0
_C267: - 1000 E1_eoabs_flash + 1000 E2_eoabs_flash - X_eoabs + X_flash
 - 0.5 l_eoabs - 0.5 l_flash >= -1000
_C268: 1000 E1_eoabs_flash - 1000 E2_eoabs_flash + Y_eoabs - Y_flash
 - 0.5 d_eoabs - 0.5 d_flash >= -1000
_C269: - 1000 E1_eoabs_flash - 1000 E2_eoabs_flash - Y_eoabs + Y_flash
 - 0.5 d_eoabs - 0.5 d_flash >= -2000
_C27: X_pump - 0.5 l_pump >= 0
_C270: L_eoabs_flash - R_flash_eoabs = 0
_C271: - L_flash_eoabs + R_eoabs_flash = 0
_C272: A_eoabs_flash - B_flash_eoabs = 0
_C273: - A_flash_eoabs + B_eoabs_flash = 0
_C274: Wx_eoabs_flash + Wx_flash_eoabs = 1
_C275: Wy_eoabs_flash + Wy_flash_eoabs = 1
_C276: CD_hex2_flash = 0
_C277: - L_hex2_flash + R_hex2_flash + X_flash - X_hex2 = 0
_C278: A_hex2_flash - B_hex2_flash + Y_flash - Y_hex2 = 0
_C279: R_hex2_flash - 1000 Wx_hex2_flash <= 0
_C28: Y_pump - 0.5 d_pump >= 0
_C280: L_hex2_flash + 1000 Wx_hex2_flash <= 1000
_C281: A_hex2_flash - 1000 Wy_hex2_flash <= 0
_C282: B_hex2_flash + 1000 Wy_hex2_flash <= 1000
_C283: - A_hex2_flash - B_hex2_flash + D_hex2_flash - L_hex2_flash
 - R_hex2_flash = 0
_C284: - D_flash_hex2 + D_hex2_flash = 0
_C285: 1000 E1_hex2_flash + 1000 E2_hex2_flash - X_flash + X_hex2
 - 0.5 l_flash - 0.5 l_hex2 >= 0
_C286: - 1000 E1_hex2_flash + 1000 E2_hex2_flash + X_flash - X_hex2
 - 0.5 l_flash - 0.5 l_hex2 >= -1000
_C287: 1000 E1_hex2_flash - 1000 E2_hex2_flash - Y_flash + Y_hex2
 - 0.5 d_flash - 0.5 d_hex2 >= -1000
_C288: - 1000 E1_hex2_flash - 1000 E2_hex2_flash + Y_flash - Y_hex2
 - 0.5 d_flash - 0.5 d_hex2 >= -2000
_C289: L_hex2_flash - R_flash_hex2 = 0
_C29: CD_reactor_hex1 - 2730.03837522 D_reactor_hex1 = 0
_C290: - L_flash_hex2 + R_hex2_flash = 0
_C291: A_hex2_flash - B_flash_hex2 = 0
_C292: - A_flash_hex2 + B_hex2_flash = 0
_C293: Wx_flash_hex2 + Wx_hex2_flash = 1
_C294: Wy_flash_hex2 + Wy_hex2_flash = 1
_C295: CD_co2abs_flash - 2730.03837522 D_co2abs_flash = 0
_C296: - L_co2abs_flash + R_co2abs_flash - X_co2abs + X_flash = 0
_C297: A_co2abs_flash - B_co2abs_flash - Y_co2abs + Y_flash = 0
_C298: R_co2abs_flash - 1000 Wx_co2abs_flash <= 0
_C299: L_co2abs_flash + 1000 Wx_co2abs_flash <= 1000
_C3: X_reactor - 0.5 l_reactor >= 0
_C30: - L_reactor_hex1 + R_reactor_hex1 + X_hex1 - X_reactor = 0
_C300: A_co2abs_flash - 1000 Wy_co2abs_flash <= 0
_C301: B_co2abs_flash + 1000 Wy_co2abs_flash <= 1000
_C302: - A_co2abs_flash - B_co2abs_flash + D_co2abs_flash - L_co2abs_flash
 - R_co2abs_flash = 0
_C303: D_co2abs_flash - D_flash_co2abs = 0
_C304: 1000 E1_co2abs_flash + 1000 E2_co2abs_flash + X_co2abs - X_flash
 - 0.5 l_co2abs - 0.5 l_flash >= 0
_C305: - 1000 E1_co2abs_flash + 1000 E2_co2abs_flash - X_co2abs + X_flash
 - 0.5 l_co2abs - 0.5 l_flash >= -1000
_C306: 1000 E1_co2abs_flash - 1000 E2_co2abs_flash + Y_co2abs - Y_flash
 - 0.5 d_co2abs - 0.5 d_flash >= -1000
_C307: - 1000 E1_co2abs_flash - 1000 E2_co2abs_flash - Y_co2abs + Y_flash
 - 0.5 d_co2abs - 0.5 d_flash >= -2000
_C308: L_co2abs_flash - R_flash_co2abs = 0
_C309: - L_flash_co2abs + R_co2abs_flash = 0
_C31: A_reactor_hex1 - B_reactor_hex1 + Y_hex1 - Y_reactor = 0
_C310: A_co2abs_flash - B_flash_co2abs = 0
_C311: - A_flash_co2abs + B_co2abs_flash = 0
_C312: Wx_co2abs_flash + Wx_flash_co2abs = 1
_C313: Wy_co2abs_flash + Wy_flash_co2abs = 1
_C314: CD_reactor_pump = 0
_C315: - L_reactor_pump + R_reactor_pump + X_pump - X_reactor = 0
_C316: A_reactor_pump - B_reactor_pump + Y_pump - Y_reactor = 0
_C317: R_reactor_pump - 1000 Wx_reactor_pump <= 0
_C318: L_reactor_pump + 1000 Wx_reactor_pump <= 1000
_C319: A_reactor_pump - 1000 Wy_reactor_pump <= 0
_C32: R_reactor_hex1 - 1000 Wx_reactor_hex1 <= 0
_C320: B_reactor_pump + 1000 Wy_reactor_pump <= 1000
_C321: - A_reactor_pump - B_reactor_pump + D_reactor_pump - L_reactor_pump
 - R_reactor_pump = 0
_C322: - D_pump_reactor + D_reactor_pump = 0
_C323: 1000 E1_reactor_pump + 1000 E2_reactor_pump - X_pump + X_reactor
 - 0.5 l_pump - 0.5 l_reactor >= 0
_C324: - 1000 E1_reactor_pump + 1000 E2_reactor_pump + X_pump - X_reactor
 - 0.5 l_pump - 0.5 l_reactor >= -1000
_C325: 1000 E1_reactor_pump - 1000 E2_reactor_pump - Y_pump + Y_reactor
 - 0.5 d_pump - 0.5 d_reactor >= -1000
_C326: - 1000 E1_reactor_pump - 1000 E2_reactor_pump + Y_pump - Y_reactor
 - 0.5 d_pump - 0.5 d_reactor >= -2000
_C327: L_reactor_pump - R_pump_reactor = 0
_C328: - L_pump_reactor + R_reactor_pump = 0
_C329: A_reactor_pump - B_pump_reactor = 0
_C33: L_reactor_hex1 + 1000 Wx_reactor_hex1 <= 1000
_C330: - A_pump_reactor + B_reactor_pump = 0
_C331: Wx_pump_reactor + Wx_reactor_pump = 1
_C332: Wy_pump_reactor + Wy_reactor_pump = 1
_C333: CD_hex1_pump = 0
_C334: - L_hex1_pump + R_hex1_pump - X_hex1 + X_pump = 0
_C335: A_hex1_pump - B_hex1_pump - Y_hex1 + Y_pump = 0
_C336: R_hex1_pump - 1000 Wx_hex1_pump <= 0
_C337: L_hex1_pump + 1000 Wx_hex1_pump <= 1000
_C338: A_hex1_pump - 1000 Wy_hex1_pump <= 0
_C339: B_hex1_pump + 1000 Wy_hex1_pump <= 1000
_C34: A_reactor_hex1 - 1000 Wy_reactor_hex1 <= 0
_C340: - A_hex1_pump - B_hex1_pump + D_hex1_pump - L_hex1_pump - R_hex1_pump
 = 0
_C341: D_hex1_pump - D_pump_hex1 = 0
_C342: 1000 E1_hex1_pump + 1000 E2_hex1_pump + X_hex1 - X_pump - 0.5 l_hex1
 - 0.5 l_pump >= 0
_C343: - 1000 E1_hex1_pump + 1000 E2_hex1_pump - X_hex1 + X_pump - 0.5 l_hex1
 - 0.5 l_pump >= -1000
_C344: 1000 E1_hex1_pump - 1000 E2_hex1_pump + Y_hex1 - Y_pump - 0.5 d_hex1
 - 0.5 d_pump >= -1000
_C345: - 1000 E1_hex1_pump - 1000 E2_hex1_pump - Y_hex1 + Y_pump - 0.5 d_hex1
 - 0.5 d_pump >= -2000
_C346: L_hex1_pump - R_pump_hex1 = 0
_C347: - L_pump_hex1 + R_hex1_pump = 0
_C348: A_hex1_pump - B_pump_hex1 = 0
_C349: - A_pump_hex1 + B_hex1_pump = 0
_C35: B_reactor_hex1 + 1000 Wy_reactor_hex1 <= 1000
_C350: Wx_hex1_pump + Wx_pump_hex1 = 1
_C351: Wy_hex1_pump + Wy_pump_hex1 = 1
_C352: CD_eoabs_pump = 0
_C353: - L_eoabs_pump + R_eoabs_pump - X_eoabs + X_pump = 0
_C354: A_eoabs_pump - B_eoabs_pump - Y_eoabs + Y_pump = 0
_C355: R_eoabs_pump - 1000 Wx_eoabs_pump <= 0
_C356: L_eoabs_pump + 1000 Wx_eoabs_pump <= 1000
_C357: A_eoabs_pump - 1000 Wy_eoabs_pump <= 0
_C358: B_eoabs_pump + 1000 Wy_eoabs_pump <= 1000
_C359: - A_eoabs_pump - B_eoabs_pump + D_eoabs_pump - L_eoabs_pump
 - R_eoabs_pump = 0
_C36: - A_reactor_hex1 - B_reactor_hex1 + D_reactor_hex1 - L_reactor_hex1
 - R_reactor_hex1 = 0
_C360: D_eoabs_pump - D_pump_eoabs = 0
_C361: 1000 E1_eoabs_pump + 1000 E2_eoabs_pump + X_eoabs - X_pump
 - 0.5 l_eoabs - 0.5 l_pump >= 0
_C362: - 1000 E1_eoabs_pump + 1000 E2_eoabs_pump - X_eoabs + X_pump
 - 0.5 l_eoabs - 0.5 l_pump >= -1000
_C363: 1000 E1_eoabs_pump - 1000 E2_eoabs_pump + Y_eoabs - Y_pump
 - 0.5 d_eoabs - 0.5 d_pump >= -1000
_C364: - 1000 E1_eoabs_pump - 1000 E2_eoabs_pump - Y_eoabs + Y_pump
 - 0.5 d_eoabs - 0.5 d_pump >= -2000
_C365: L_eoabs_pump - R_pump_eoabs = 0
_C366: - L_pump_eoabs + R_eoabs_pump = 0
_C367: A_eoabs_pump - B_pump_eoabs = 0
_C368: - A_pump_eoabs + B_eoabs_pump = 0
_C369: Wx_eoabs_pump + Wx_pump_eoabs = 1
_C37: - D_hex1_reactor + D_reactor_hex1 = 0
_C370: Wy_eoabs_pump + Wy_pump_eoabs = 1
_C371: CD_hex2_pump = 0
_C372: - L_hex2_pump + R_hex2_pump - X_hex2 + X_pump = 0
_C373: A_hex2_pump - B_hex2_pump - Y_hex2 + Y_pump = 0
_C374: R_hex2_pump - 1000 Wx_hex2_pump <= 0
_C375: L_hex2_pump + 1000 Wx_hex2_pump <= 1000
_C376: A_hex2_pump - 1000 Wy_hex2_pump <= 0
_C377: B_hex2_pump + 1000 Wy_hex2_pump <= 1000
_C378: - A_hex2_pump - B_hex2_pump + D_hex2_pump - L_hex2_pump - R_hex2_pump
 = 0
_C379: D_hex2_pump - D_pump_hex2 = 0
_C38: 1000 E1_reactor_hex1 + 1000 E2_reactor_hex1 - X_hex1 + X_reactor
 - 0.5 l_hex1 - 0.5 l_reactor >= 1
_C380: 1000 E1_hex2_pump + 1000 E2_hex2_pump + X_hex2 - X_pump - 0.5 l_hex2
 - 0.5 l_pump >= 0
_C381: - 1000 E1_hex2_pump + 1000 E2_hex2_pump - X_hex2 + X_pump - 0.5 l_hex2
 - 0.5 l_pump >= -1000
_C382: 1000 E1_hex2_pump - 1000 E2_hex2_pump + Y_hex2 - Y_pump - 0.5 d_hex2
 - 0.5 d_pump >= -1000
_C383: - 1000 E1_hex2_pump - 1000 E2_hex2_pump - Y_hex2 + Y_pump - 0.5 d_hex2
 - 0.5 d_pump >= -2000
_C384: L_hex2_pump - R_pump_hex2 = 0
_C385: - L_pump_hex2 + R_hex2_pump = 0
_C386: A_hex2_pump - B_pump_hex2 = 0
_C387: - A_pump_hex2 + B_hex2_pump = 0
_C388: Wx_hex2_pump + Wx_pump_hex2 = 1
_C389: Wy_hex2_pump + Wy_pump_hex2 = 1
_C39: - 1000 E1_reactor_hex1 + 1000 E2_reactor_hex1 + X_hex1 - X_reactor
 - 0.5 l_hex1 - 0.5 l_reactor >= -999
_C390: CD_co2abs_pump - 2730.03837522 D_co2abs_pump = 0
_C391: - L_co2abs_pump + R_co2abs_pump - X_co2abs + X_pump = 0
_C392: A_co2abs_pump - B_co2abs_pump - Y_co2abs + Y_pump = 0
_C393: R_co2abs_pump - 1000 Wx_co2abs_pump <= 0
_C394: L_co2abs_pump + 1000 Wx_co2abs_pump <= 1000
_C395: A_co2abs_pump - 1000 Wy_co2abs_pump <= 0
_C396: B_co2abs_pump + 1000 Wy_co2abs_pump <= 1000
_C397: - A_co2abs_pump - B_co2abs_pump + D_co2abs_pump - L_co2abs_pump
 - R_co2abs_pump = 0
_C398: D_co2abs_pump - D_pump_co2abs = 0
_C399: 1000 E1_co2abs_pump + 1000 E2_co2abs_pump + X_co2abs - X_pump
 - 0.5 l_co2abs - 0.5 l_pump >= 0
_C4: Y_reactor - 0.5 d_reactor >= 0
_C40: 1000 E1_reactor_hex1 - 1000 E2_reactor_hex1 - Y_hex1 + Y_reactor
 - 0.5 d_hex1 - 0.5 d_reactor >= -999
_C400: - 1000 E1_co2abs_pump + 1000 E2_co2abs_pump - X_co2abs + X_pump
 - 0.5 l_co2abs - 0.5 l_pump >= -1000
_C401: 1000 E1_co2abs_pump - 1000 E2_co2abs_pump + Y_co2abs - Y_pump
 - 0.5 d_co2abs - 0.5 d_pump >= -1000
_C402: - 1000 E1_co2abs_pump - 1000 E2_co2abs_pump - Y_co2abs + Y_pump
 - 0.5 d_co2abs - 0.5 d_pump >= -2000
_C403: L_co2abs_pump - R_pump_co2abs = 0
_C404: - L_pump_co2abs + R_co2abs_pump = 0
_C405: A_co2abs_pump - B_pump_co2abs = 0
_C406: - A_pump_co2abs + B_co2abs_pump = 0
_C407: Wx_co2abs_pump + Wx_pump_co2abs = 1
_C408: Wy_co2abs_pump + Wy_pump_co2abs = 1
_C409: CD_flash_pump - 2730.03837522 D_flash_pump = 0
_C41: - 1000 E1_reactor_hex1 - 1000 E2_reactor_hex1 + Y_hex1 - Y_reactor
 - 0.5 d_hex1 - 0.5 d_reactor >= -1999
_C410: - L_flash_pump + R_flash_pump - X_flash + X_pump = 0
_C411: A_flash_pump - B_flash_pump - Y_flash + Y_pump = 0
_C412: R_flash_pump - 1000 Wx_flash_pump <= 0
_C413: L_flash_pump + 1000 Wx_flash_pump <= 1000
_C414: A_flash_pump - 1000 Wy_flash_pump <= 0
_C415: B_flash_pump + 1000 Wy_flash_pump <= 1000
_C416: - A_flash_pump - B_flash_pump + D_flash_pump - L_flash_pump
 - R_flash_pump = 0
_C417: D_flash_pump - D_pump_flash = 0
_C418: 1000 E1_flash_pump + 1000 E2_flash_pump + X_flash - X_pump
 - 0.5 l_flash - 0.5 l_pump >= 20
_C419: - 1000 E1_flash_pump + 1000 E2_flash_pump - X_flash + X_pump
 - 0.5 l_flash - 0.5 l_pump >= -980
_C42: L_reactor_hex1 - R_hex1_reactor = 0
_C420: 1000 E1_flash_pump - 1000 E2_flash_pump + Y_flash - Y_pump
 - 0.5 d_flash - 0.5 d_pump >= -980
_C421: - 1000 E1_flash_pump - 1000 E2_flash_pump - Y_flash + Y_pump
 - 0.5 d_flash - 0.5 d_pump >= -1980
_C422: L_flash_pump - R_pump_flash = 0
_C423: - L_pump_flash + R_flash_pump = 0
_C424: A_flash_pump - B_pump_flash = 0
_C425: - A_pump_flash + B_flash_pump = 0
_C426: Wx_flash_pump + Wx_pump_flash = 1
_C427: Wy_flash_pump + Wy_pump_flash = 1
_C428: - CD_co2abs_co2abs - CD_co2abs_eoabs - CD_co2abs_flash - CD_co2abs_hex1
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
_C429: - 0.933333333333 X_reactor + Y_reactor + 0.5 d_reactor
 + 0.466666666667 l_reactor <= 21
_C43: - L_hex1_reactor + R_reactor_hex1 = 0
_C430: 0.933333333333 X_reactor + Y_reactor + 0.5 d_reactor
 + 0.466666666667 l_reactor <= 49
_C431: Y_reactor - 0.5 d_reactor >= 0
_C432: X_reactor + 0.5 l_reactor <= 30
_C433: - 0.933333333333 X_hex1 + Y_hex1 + 0.5 d_hex1 + 0.466666666667 l_hex1
 <= 21
_C434: 0.933333333333 X_hex1 + Y_hex1 + 0.5 d_hex1 + 0.466666666667 l_hex1
 <= 49
_C435: Y_hex1 - 0.5 d_hex1 >= 0
_C436: X_hex1 + 0.5 l_hex1 <= 30
_C437: - 0.933333333333 X_eoabs + Y_eoabs + 0.5 d_eoabs
 + 0.466666666667 l_eoabs <= 21
_C438: 0.933333333333 X_eoabs + Y_eoabs + 0.5 d_eoabs + 0.466666666667 l_eoabs
 <= 49
_C439: Y_eoabs - 0.5 d_eoabs >= 0
_C44: A_reactor_hex1 - B_hex1_reactor = 0
_C440: X_eoabs + 0.5 l_eoabs <= 30
_C441: - 0.933333333333 X_hex2 + Y_hex2 + 0.5 d_hex2 + 0.466666666667 l_hex2
 <= 21
_C442: 0.933333333333 X_hex2 + Y_hex2 + 0.5 d_hex2 + 0.466666666667 l_hex2
 <= 49
_C443: Y_hex2 - 0.5 d_hex2 >= 0
_C444: X_hex2 + 0.5 l_hex2 <= 30
_C445: - 0.933333333333 X_co2abs + Y_co2abs + 0.5 d_co2abs
 + 0.466666666667 l_co2abs <= 21
_C446: 0.933333333333 X_co2abs + Y_co2abs + 0.5 d_co2abs
 + 0.466666666667 l_co2abs <= 49
_C447: Y_co2abs - 0.5 d_co2abs >= 0
_C448: X_co2abs + 0.5 l_co2abs <= 30
_C449: - 0.933333333333 X_flash + Y_flash + 0.5 d_flash
 + 0.466666666667 l_flash <= 21
_C45: - A_hex1_reactor + B_reactor_hex1 = 0
_C450: 0.933333333333 X_flash + Y_flash + 0.5 d_flash + 0.466666666667 l_flash
 <= 49
_C451: Y_flash - 0.5 d_flash >= 0
_C452: X_flash + 0.5 l_flash <= 30
_C453: - 0.933333333333 X_pump + Y_pump + 0.5 d_pump + 0.466666666667 l_pump
 <= 21
_C454: 0.933333333333 X_pump + Y_pump + 0.5 d_pump + 0.466666666667 l_pump
 <= 49
_C455: Y_pump - 0.5 d_pump >= 0
_C456: X_pump + 0.5 l_pump <= 30
_C457: Ve2_reactor_reactor = 0
_C458: - D_reactor_hex1 + Din_reactor_hex1 + Dout_reactor_hex1 = 0
_C459: Din_reactor_hex1 - 40 Psie_reactor_hex1 <= 0
_C46: Wx_hex1_reactor + Wx_reactor_hex1 = 1
_C460: Dout_reactor_hex1 + 40 Psie_reactor_hex1 >= 40
_C461: Dout_reactor_hex1 + 1000 Psie_reactor_hex1 <= 1000
_C462: 12775 Din_reactor_hex1 - 511000 Psie_reactor_hex1 + Ve2_reactor_hex1
 = 0
_C463: - D_reactor_eoabs + Din_reactor_eoabs + Dout_reactor_eoabs = 0
_C464: Din_reactor_eoabs - 40 Psie_reactor_eoabs <= 0
_C465: Dout_reactor_eoabs + 40 Psie_reactor_eoabs >= 40
_C466: Dout_reactor_eoabs + 1000 Psie_reactor_eoabs <= 1000
_C467: 77675 Din_reactor_eoabs - 3107000 Psie_reactor_eoabs
 + Ve2_reactor_eoabs = 0
_C468: - D_reactor_hex2 + Din_reactor_hex2 + Dout_reactor_hex2 = 0
_C469: Din_reactor_hex2 - 40 Psie_reactor_hex2 <= 0
_C47: Wy_hex1_reactor + Wy_reactor_hex1 = 1
_C470: Dout_reactor_hex2 + 40 Psie_reactor_hex2 >= 40
_C471: Dout_reactor_hex2 + 1000 Psie_reactor_hex2 <= 1000
_C472: 12600 Din_reactor_hex2 - 504000 Psie_reactor_hex2 + Ve2_reactor_hex2
 = 0
_C473: - D_reactor_co2abs + Din_reactor_co2abs + Dout_reactor_co2abs = 0
_C474: Din_reactor_co2abs - 40 Psie_reactor_co2abs <= 0
_C475: Dout_reactor_co2abs + 40 Psie_reactor_co2abs >= 40
_C476: Dout_reactor_co2abs + 1000 Psie_reactor_co2abs <= 1000
_C477: 77032.5 Din_reactor_co2abs - 3081300 Psie_reactor_co2abs
 + Ve2_reactor_co2abs = 0
_C478: - D_reactor_flash + Din_reactor_flash + Dout_reactor_flash = 0
_C479: Din_reactor_flash - 40 Psie_reactor_flash <= 0
_C48: CD_reactor_eoabs = 0
_C480: Dout_reactor_flash + 40 Psie_reactor_flash >= 40
_C481: Dout_reactor_flash + 1000 Psie_reactor_flash <= 1000
_C482: 12625 Din_reactor_flash - 505000 Psie_reactor_flash + Ve2_reactor_flash
 = 0
_C483: - D_reactor_pump + Din_reactor_pump + Dout_reactor_pump = 0
_C484: Din_reactor_pump - 40 Psie_reactor_pump <= 0
_C485: Dout_reactor_pump + 40 Psie_reactor_pump >= 40
_C486: Dout_reactor_pump + 1000 Psie_reactor_pump <= 1000
_C487: 12537.5 Din_reactor_pump - 501500 Psie_reactor_pump + Ve2_reactor_pump
 = 0
_C488: - D_eoabs_reactor + Din_eoabs_reactor + Dout_eoabs_reactor = 0
_C489: Din_eoabs_reactor - 21.8 Psie_eoabs_reactor <= 0
_C49: - L_reactor_eoabs + R_reactor_eoabs + X_eoabs - X_reactor = 0
_C490: Dout_eoabs_reactor + 21.8 Psie_eoabs_reactor >= 21.8
_C491: Dout_eoabs_reactor + 1000 Psie_eoabs_reactor <= 1000
_C492: 152981.651376 Din_eoabs_reactor - 3335000 Psie_eoabs_reactor
 + Ve2_eoabs_reactor = 0
_C493: - D_eoabs_hex1 + Din_eoabs_hex1 + Dout_eoabs_hex1 = 0
_C494: Din_eoabs_hex1 - 21.8 Psie_eoabs_hex1 <= 0
_C495: Dout_eoabs_hex1 + 21.8 Psie_eoabs_hex1 >= 21.8
_C496: Dout_eoabs_hex1 + 1000 Psie_eoabs_hex1 <= 1000
_C497: 23440.3669725 Din_eoabs_hex1 - 511000 Psie_eoabs_hex1 + Ve2_eoabs_hex1
 = 0
_C498: Ve2_eoabs_eoabs = 0
_C499: - D_eoabs_hex2 + Din_eoabs_hex2 + Dout_eoabs_hex2 = 0
_C5: 0 O_hex1 + l_hex1 = 11.42
_C50: A_reactor_eoabs - B_reactor_eoabs + Y_eoabs - Y_reactor = 0
_C500: Din_eoabs_hex2 - 21.8 Psie_eoabs_hex2 <= 0
_C501: Dout_eoabs_hex2 + 21.8 Psie_eoabs_hex2 >= 21.8
_C502: Dout_eoabs_hex2 + 1000 Psie_eoabs_hex2 <= 1000
_C503: 23119.266055 Din_eoabs_hex2 - 504000 Psie_eoabs_hex2 + Ve2_eoabs_hex2
 = 0
_C504: - D_eoabs_co2abs + Din_eoabs_co2abs + Dout_eoabs_co2abs = 0
_C505: Din_eoabs_co2abs - 21.8 Psie_eoabs_co2abs <= 0
_C506: Dout_eoabs_co2abs + 21.8 Psie_eoabs_co2abs >= 21.8
_C507: Dout_eoabs_co2abs + 1000 Psie_eoabs_co2abs <= 1000
_C508: 141344.036697 Din_eoabs_co2abs - 3081300 Psie_eoabs_co2abs
 + Ve2_eoabs_co2abs = 0
_C509: - D_eoabs_flash + Din_eoabs_flash + Dout_eoabs_flash = 0
_C51: R_reactor_eoabs - 1000 Wx_reactor_eoabs <= 0
_C510: Din_eoabs_flash - 21.8 Psie_eoabs_flash <= 0
_C511: Dout_eoabs_flash + 21.8 Psie_eoabs_flash >= 21.8
_C512: Dout_eoabs_flash + 1000 Psie_eoabs_flash <= 1000
_C513: 23165.1376147 Din_eoabs_flash - 505000 Psie_eoabs_flash
 + Ve2_eoabs_flash = 0
_C514: - D_eoabs_pump + Din_eoabs_pump + Dout_eoabs_pump = 0
_C515: Din_eoabs_pump - 21.8 Psie_eoabs_pump <= 0
_C516: Dout_eoabs_pump + 21.8 Psie_eoabs_pump >= 21.8
_C517: Dout_eoabs_pump + 1000 Psie_eoabs_pump <= 1000
_C518: 23004.587156 Din_eoabs_pump - 501500 Psie_eoabs_pump + Ve2_eoabs_pump
 = 0
_C519: - D_co2abs_reactor + Din_co2abs_reactor + Dout_co2abs_reactor = 0
_C52: L_reactor_eoabs + 1000 Wx_reactor_eoabs <= 1000
_C520: Din_co2abs_reactor - 18.06 Psie_co2abs_reactor <= 0
_C521: Dout_co2abs_reactor + 18.06 Psie_co2abs_reactor >= 18.06
_C522: Dout_co2abs_reactor + 1000 Psie_co2abs_reactor <= 1000
_C523: 184662.236988 Din_co2abs_reactor - 3335000 Psie_co2abs_reactor
 + Ve2_co2abs_reactor = 0
_C524: - D_co2abs_hex1 + Din_co2abs_hex1 + Dout_co2abs_hex1 = 0
_C525: Din_co2abs_hex1 - 18.06 Psie_co2abs_hex1 <= 0
_C526: Dout_co2abs_hex1 + 18.06 Psie_co2abs_hex1 >= 18.06
_C527: Dout_co2abs_hex1 + 1000 Psie_co2abs_hex1 <= 1000
_C528: 28294.5736434 Din_co2abs_hex1 - 511000 Psie_co2abs_hex1
 + Ve2_co2abs_hex1 = 0
_C529: - D_co2abs_eoabs + Din_co2abs_eoabs + Dout_co2abs_eoabs = 0
_C53: A_reactor_eoabs - 1000 Wy_reactor_eoabs <= 0
_C530: Din_co2abs_eoabs - 18.06 Psie_co2abs_eoabs <= 0
_C531: Dout_co2abs_eoabs + 18.06 Psie_co2abs_eoabs >= 18.06
_C532: Dout_co2abs_eoabs + 1000 Psie_co2abs_eoabs <= 1000
_C533: 172037.65227 Din_co2abs_eoabs - 3107000 Psie_co2abs_eoabs
 + Ve2_co2abs_eoabs = 0
_C534: - D_co2abs_hex2 + Din_co2abs_hex2 + Dout_co2abs_hex2 = 0
_C535: Din_co2abs_hex2 - 18.06 Psie_co2abs_hex2 <= 0
_C536: Dout_co2abs_hex2 + 18.06 Psie_co2abs_hex2 >= 18.06
_C537: Dout_co2abs_hex2 + 1000 Psie_co2abs_hex2 <= 1000
_C538: 27906.9767442 Din_co2abs_hex2 - 504000 Psie_co2abs_hex2
 + Ve2_co2abs_hex2 = 0
_C539: Ve2_co2abs_co2abs = 0
_C54: B_reactor_eoabs + 1000 Wy_reactor_eoabs <= 1000
_C540: - D_co2abs_flash + Din_co2abs_flash + Dout_co2abs_flash = 0
_C541: Din_co2abs_flash - 18.06 Psie_co2abs_flash <= 0
_C542: Dout_co2abs_flash + 18.06 Psie_co2abs_flash >= 18.06
_C543: Dout_co2abs_flash + 1000 Psie_co2abs_flash <= 1000
_C544: 27962.3477298 Din_co2abs_flash - 505000 Psie_co2abs_flash
 + Ve2_co2abs_flash = 0
_C545: - D_co2abs_pump + Din_co2abs_pump + Dout_co2abs_pump = 0
_C546: Din_co2abs_pump - 18.06 Psie_co2abs_pump <= 0
_C547: Dout_co2abs_pump + 18.06 Psie_co2abs_pump >= 18.06
_C548: Dout_co2abs_pump + 1000 Psie_co2abs_pump <= 1000
_C549: 27768.5492802 Din_co2abs_pump - 501500 Psie_co2abs_pump
 + Ve2_co2abs_pump = 0
_C55: - A_reactor_eoabs - B_reactor_eoabs + D_reactor_eoabs - L_reactor_eoabs
 - R_reactor_eoabs = 0
_C550: - Ve2_reactor_co2abs - Ve2_reactor_eoabs - Ve2_reactor_flash
 - Ve2_reactor_hex1 - Ve2_reactor_hex2 - Ve2_reactor_pump
 - Ve2_reactor_reactor + Ve_reactor = 3335000
_C551: Omega0_reactor - 0.87 Ve_reactor = 0
_C552: - Ve2_eoabs_co2abs - Ve2_eoabs_eoabs - Ve2_eoabs_flash - Ve2_eoabs_hex1
 - Ve2_eoabs_hex2 - Ve2_eoabs_pump - Ve2_eoabs_reactor + Ve_eoabs = 3107000
_C553: Omega0_eoabs - 0.73 Ve_eoabs = 0
_C554: - Ve2_co2abs_co2abs - Ve2_co2abs_eoabs - Ve2_co2abs_flash
 - Ve2_co2abs_hex1 - Ve2_co2abs_hex2 - Ve2_co2abs_pump - Ve2_co2abs_reactor
 + Ve_co2abs = 3081300
_C555: Omega0_co2abs - 0.66 Ve_co2abs = 0
_C556: - 100000000 Z_reactor_1 + linOmega0_reactor_1 <= 0
_C557: - 100000000 Z_reactor_2 + linOmega0_reactor_2 <= 0
_C558: - 100000000 Z_reactor_3 + linOmega0_reactor_3 <= 0
_C559: - 100000000 Z_reactor_4 + linOmega0_reactor_4 <= 0
_C56: - D_eoabs_reactor + D_reactor_eoabs = 0
_C560: - 100000000 Z_reactor_5 + linOmega0_reactor_5 <= 0
_C561: - 100000000 Z_reactor_6 + linOmega0_reactor_6 <= 0
_C562: - 100000000 Z_eoabs_1 + linOmega0_eoabs_1 <= 0
_C563: - 100000000 Z_eoabs_2 + linOmega0_eoabs_2 <= 0
_C564: - 100000000 Z_eoabs_3 + linOmega0_eoabs_3 <= 0
_C565: - 100000000 Z_eoabs_4 + linOmega0_eoabs_4 <= 0
_C566: - 100000000 Z_eoabs_5 + linOmega0_eoabs_5 <= 0
_C567: - 100000000 Z_eoabs_6 + linOmega0_eoabs_6 <= 0
_C568: - 100000000 Z_co2abs_1 + linOmega0_co2abs_1 <= 0
_C569: - 100000000 Z_co2abs_2 + linOmega0_co2abs_2 <= 0
_C57: 1000 E1_reactor_eoabs + 1000 E2_reactor_eoabs - X_eoabs + X_reactor
 - 0.5 l_eoabs - 0.5 l_reactor >= 0
_C570: - 100000000 Z_co2abs_3 + linOmega0_co2abs_3 <= 0
_C571: - 100000000 Z_co2abs_4 + linOmega0_co2abs_4 <= 0
_C572: - 100000000 Z_co2abs_5 + linOmega0_co2abs_5 <= 0
_C573: - 100000000 Z_co2abs_6 + linOmega0_co2abs_6 <= 0
_C574: Z_reactor_1 + Z_reactor_2 + Z_reactor_3 + Z_reactor_4 + Z_reactor_5
 + Z_reactor_6 = 1
_C575: Omega0_reactor - linOmega0_reactor_1 - linOmega0_reactor_2
 - linOmega0_reactor_3 - linOmega0_reactor_4 - linOmega0_reactor_5
 - linOmega0_reactor_6 = 0
_C576: Omega_reactor - linOmega0_reactor_1 - 0.9 linOmega0_reactor_2
 - 0.75 linOmega0_reactor_3 - 0.365 linOmega0_reactor_4
 - 0.292 linOmega0_reactor_5 - 0.117 linOmega0_reactor_6 = 0
_C577: - Omega_co2abs - Omega_eoabs - Omega_reactor + SumOmega = 0
_C578: SumPZ - 5000 Z_co2abs_2 - 20000 Z_co2abs_3 - 25000 Z_co2abs_4
 - 35000 Z_co2abs_5 - 55000 Z_co2abs_6 - 5000 Z_eoabs_2 - 20000 Z_eoabs_3
 - 25000 Z_eoabs_4 - 35000 Z_eoabs_5 - 55000 Z_eoabs_6 - 5000 Z_reactor_2
 - 15000 Z_reactor_3 - 40000 Z_reactor_4 - 60000 Z_reactor_5
 - 125000 Z_reactor_6 = 0
_C579: Z_eoabs_1 + Z_eoabs_2 + Z_eoabs_3 + Z_eoabs_4 + Z_eoabs_5 + Z_eoabs_6
 = 1
_C58: - 1000 E1_reactor_eoabs + 1000 E2_reactor_eoabs + X_eoabs - X_reactor
 - 0.5 l_eoabs - 0.5 l_reactor >= -1000
_C580: Omega0_eoabs - linOmega0_eoabs_1 - linOmega0_eoabs_2
 - linOmega0_eoabs_3 - linOmega0_eoabs_4 - linOmega0_eoabs_5
 - linOmega0_eoabs_6 = 0
_C581: Omega_eoabs - linOmega0_eoabs_1 - 0.9 linOmega0_eoabs_2
 - 0.76 linOmega0_eoabs_3 - 0.684 linOmega0_eoabs_4 - 0.612 linOmega0_eoabs_5
 - 0.465 linOmega0_eoabs_6 = 0
_C582: - Omega_co2abs - Omega_eoabs - Omega_reactor + SumOmega = 0
_C583: SumPZ - 5000 Z_co2abs_2 - 20000 Z_co2abs_3 - 25000 Z_co2abs_4
 - 35000 Z_co2abs_5 - 55000 Z_co2abs_6 - 5000 Z_eoabs_2 - 20000 Z_eoabs_3
 - 25000 Z_eoabs_4 - 35000 Z_eoabs_5 - 55000 Z_eoabs_6 - 5000 Z_reactor_2
 - 15000 Z_reactor_3 - 40000 Z_reactor_4 - 60000 Z_reactor_5
 - 125000 Z_reactor_6 = 0
_C584: Z_co2abs_1 + Z_co2abs_2 + Z_co2abs_3 + Z_co2abs_4 + Z_co2abs_5
 + Z_co2abs_6 = 1
_C585: Omega0_co2abs - linOmega0_co2abs_1 - linOmega0_co2abs_2
 - linOmega0_co2abs_3 - linOmega0_co2abs_4 - linOmega0_co2abs_5
 - linOmega0_co2abs_6 = 0
_C586: Omega_co2abs - linOmega0_co2abs_1 - 0.9 linOmega0_co2abs_2
 - 0.76 linOmega0_co2abs_3 - 0.684 linOmega0_co2abs_4
 - 0.612 linOmega0_co2abs_5 - 0.465 linOmega0_co2abs_6 = 0
_C587: - Omega_co2abs - Omega_eoabs - Omega_reactor + SumOmega = 0
_C588: SumPZ - 5000 Z_co2abs_2 - 20000 Z_co2abs_3 - 25000 Z_co2abs_4
 - 35000 Z_co2abs_5 - 55000 Z_co2abs_6 - 5000 Z_eoabs_2 - 20000 Z_eoabs_3
 - 25000 Z_eoabs_4 - 35000 Z_eoabs_5 - 55000 Z_eoabs_6 - 5000 Z_reactor_2
 - 15000 Z_reactor_3 - 40000 Z_reactor_4 - 60000 Z_reactor_5
 - 125000 Z_reactor_6 = 0
_C59: 1000 E1_reactor_eoabs - 1000 E2_reactor_eoabs - Y_eoabs + Y_reactor
 - 0.5 d_eoabs - 0.5 d_reactor >= -1000
_C6: d_hex1 + l_hex1 = 22.84
_C60: - 1000 E1_reactor_eoabs - 1000 E2_reactor_eoabs + Y_eoabs - Y_reactor
 - 0.5 d_eoabs - 0.5 d_reactor >= -2000
_C61: L_reactor_eoabs - R_eoabs_reactor = 0
_C62: - L_eoabs_reactor + R_reactor_eoabs = 0
_C63: A_reactor_eoabs - B_eoabs_reactor = 0
_C64: - A_eoabs_reactor + B_reactor_eoabs = 0
_C65: Wx_eoabs_reactor + Wx_reactor_eoabs = 1
_C66: Wy_eoabs_reactor + Wy_reactor_eoabs = 1
_C67: CD_hex1_eoabs - 2730.03837522 D_hex1_eoabs = 0
_C68: - L_hex1_eoabs + R_hex1_eoabs + X_eoabs - X_hex1 = 0
_C69: A_hex1_eoabs - B_hex1_eoabs + Y_eoabs - Y_hex1 = 0
_C7: X_hex1 - 0.5 l_hex1 >= 0
_C70: R_hex1_eoabs - 1000 Wx_hex1_eoabs <= 0
_C71: L_hex1_eoabs + 1000 Wx_hex1_eoabs <= 1000
_C72: A_hex1_eoabs - 1000 Wy_hex1_eoabs <= 0
_C73: B_hex1_eoabs + 1000 Wy_hex1_eoabs <= 1000
_C74: - A_hex1_eoabs - B_hex1_eoabs + D_hex1_eoabs - L_hex1_eoabs
 - R_hex1_eoabs = 0
_C75: - D_eoabs_hex1 + D_hex1_eoabs = 0
_C76: 1000 E1_hex1_eoabs + 1000 E2_hex1_eoabs - X_eoabs + X_hex1 - 0.5 l_eoabs
 - 0.5 l_hex1 >= 1
_C77: - 1000 E1_hex1_eoabs + 1000 E2_hex1_eoabs + X_eoabs - X_hex1
 - 0.5 l_eoabs - 0.5 l_hex1 >= -999
_C78: 1000 E1_hex1_eoabs - 1000 E2_hex1_eoabs - Y_eoabs + Y_hex1 - 0.5 d_eoabs
 - 0.5 d_hex1 >= -999
_C79: - 1000 E1_hex1_eoabs - 1000 E2_hex1_eoabs + Y_eoabs - Y_hex1
 - 0.5 d_eoabs - 0.5 d_hex1 >= -1999
_C8: Y_hex1 - 0.5 d_hex1 >= 0
_C80: L_hex1_eoabs - R_eoabs_hex1 = 0
_C81: - L_eoabs_hex1 + R_hex1_eoabs = 0
_C82: A_hex1_eoabs - B_eoabs_hex1 = 0
_C83: - A_eoabs_hex1 + B_hex1_eoabs = 0
_C84: Wx_eoabs_hex1 + Wx_hex1_eoabs = 1
_C85: Wy_eoabs_hex1 + Wy_hex1_eoabs = 1
_C86: CD_reactor_hex2 = 0
_C87: - L_reactor_hex2 + R_reactor_hex2 + X_hex2 - X_reactor = 0
_C88: A_reactor_hex2 - B_reactor_hex2 + Y_hex2 - Y_reactor = 0
_C89: R_reactor_hex2 - 1000 Wx_reactor_hex2 <= 0
_C9: 0 O_eoabs + l_eoabs = 7.68
_C90: L_reactor_hex2 + 1000 Wx_reactor_hex2 <= 1000
_C91: A_reactor_hex2 - 1000 Wy_reactor_hex2 <= 0
_C92: B_reactor_hex2 + 1000 Wy_reactor_hex2 <= 1000
_C93: - A_reactor_hex2 - B_reactor_hex2 + D_reactor_hex2 - L_reactor_hex2
 - R_reactor_hex2 = 0
_C94: - D_hex2_reactor + D_reactor_hex2 = 0
_C95: 1000 E1_reactor_hex2 + 1000 E2_reactor_hex2 - X_hex2 + X_reactor
 - 0.5 l_hex2 - 0.5 l_reactor >= 0
_C96: - 1000 E1_reactor_hex2 + 1000 E2_reactor_hex2 + X_hex2 - X_reactor
 - 0.5 l_hex2 - 0.5 l_reactor >= -1000
_C97: 1000 E1_reactor_hex2 - 1000 E2_reactor_hex2 - Y_hex2 + Y_reactor
 - 0.5 d_hex2 - 0.5 d_reactor >= -1000
_C98: - 1000 E1_reactor_hex2 - 1000 E2_reactor_hex2 + Y_hex2 - Y_reactor
 - 0.5 d_hex2 - 0.5 d_reactor >= -2000
_C99: L_reactor_hex2 - R_hex2_reactor = 0
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
