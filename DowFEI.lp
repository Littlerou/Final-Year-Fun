\* Layout_Problem_Model *\
Minimize
OBJ: SumCD + SumOmega + SumPZ + TLC
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
_C124: CD_eoabs_hex2 - 98.4 D_eoabs_hex2 = 0
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
_C143: CD_reactor_co2abs - 98.4 D_reactor_co2abs = 0
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
_C17: 0 O_co2abs + l_co2abs = 7.68
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
_C18: d_co2abs + l_co2abs = 15.36
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
_C200: CD_hex2_co2abs - 98.4 D_hex2_co2abs = 0
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
_C29: CD_reactor_hex1 - 98.4 D_reactor_hex1 = 0
_C290: - L_flash_hex2 + R_hex2_flash = 0
_C291: A_hex2_flash - B_flash_hex2 = 0
_C292: - A_flash_hex2 + B_hex2_flash = 0
_C293: Wx_flash_hex2 + Wx_hex2_flash = 1
_C294: Wy_flash_hex2 + Wy_hex2_flash = 1
_C295: CD_co2abs_flash - 98.4 D_co2abs_flash = 0
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
 - 0.5 l_hex1 - 0.5 l_reactor >= 0
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
 - 0.5 l_hex1 - 0.5 l_reactor >= -1000
_C390: CD_co2abs_pump - 98.4 D_co2abs_pump = 0
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
 - 0.5 d_hex1 - 0.5 d_reactor >= -1000
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
_C409: CD_flash_pump - 98.4 D_flash_pump = 0
_C41: - 1000 E1_reactor_hex1 - 1000 E2_reactor_hex1 + Y_hex1 - Y_reactor
 - 0.5 d_hex1 - 0.5 d_reactor >= -2000
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
 - 0.5 l_flash - 0.5 l_pump >= 0
_C419: - 1000 E1_flash_pump + 1000 E2_flash_pump - X_flash + X_pump
 - 0.5 l_flash - 0.5 l_pump >= -1000
_C42: L_reactor_hex1 - R_hex1_reactor = 0
_C420: 1000 E1_flash_pump - 1000 E2_flash_pump + Y_flash - Y_pump
 - 0.5 d_flash - 0.5 d_pump >= -1000
_C421: - 1000 E1_flash_pump - 1000 E2_flash_pump - Y_flash + Y_pump
 - 0.5 d_flash - 0.5 d_pump >= -2000
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
_C429: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + X_reactor
 + 0.5 l_reactor <= 0
_C43: - L_hex1_reactor + R_reactor_hex1 = 0
_C430: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + Y_reactor
 + 0.5 d_reactor <= 0
_C431: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + X_hex1
 + 0.5 l_hex1 <= 0
_C432: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + Y_hex1
 + 0.5 d_hex1 <= 0
_C433: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + X_eoabs
 + 0.5 l_eoabs <= 0
_C434: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + Y_eoabs
 + 0.5 d_eoabs <= 0
_C435: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + X_hex2
 + 0.5 l_hex2 <= 0
_C436: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + Y_hex2
 + 0.5 d_hex2 <= 0
_C437: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + X_co2abs
 + 0.5 l_co2abs <= 0
_C438: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + Y_co2abs
 + 0.5 d_co2abs <= 0
_C439: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + X_flash
 + 0.5 l_flash <= 0
_C44: A_reactor_hex1 - B_hex1_reactor = 0
_C440: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + Y_flash
 + 0.5 d_flash <= 0
_C441: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + X_pump
 + 0.5 l_pump <= 0
_C442: - 2 Gn_1 - 20 Gn_10 - 200 Gn_100 - 202 Gn_101 - 204 Gn_102 - 206 Gn_103
 - 208 Gn_104 - 210 Gn_105 - 212 Gn_106 - 214 Gn_107 - 216 Gn_108 - 218 Gn_109
 - 22 Gn_11 - 220 Gn_110 - 222 Gn_111 - 224 Gn_112 - 226 Gn_113 - 228 Gn_114
 - 230 Gn_115 - 232 Gn_116 - 234 Gn_117 - 236 Gn_118 - 238 Gn_119 - 24 Gn_12
 - 240 Gn_120 - 242 Gn_121 - 244 Gn_122 - 246 Gn_123 - 248 Gn_124 - 250 Gn_125
 - 252 Gn_126 - 254 Gn_127 - 256 Gn_128 - 258 Gn_129 - 26 Gn_13 - 260 Gn_130
 - 262 Gn_131 - 264 Gn_132 - 266 Gn_133 - 268 Gn_134 - 270 Gn_135 - 272 Gn_136
 - 274 Gn_137 - 276 Gn_138 - 278 Gn_139 - 28 Gn_14 - 280 Gn_140 - 282 Gn_141
 - 284 Gn_142 - 286 Gn_143 - 288 Gn_144 - 290 Gn_145 - 292 Gn_146 - 294 Gn_147
 - 296 Gn_148 - 298 Gn_149 - 30 Gn_15 - 300 Gn_150 - 302 Gn_151 - 304 Gn_152
 - 306 Gn_153 - 308 Gn_154 - 310 Gn_155 - 312 Gn_156 - 314 Gn_157 - 316 Gn_158
 - 318 Gn_159 - 32 Gn_16 - 320 Gn_160 - 322 Gn_161 - 324 Gn_162 - 326 Gn_163
 - 328 Gn_164 - 330 Gn_165 - 332 Gn_166 - 334 Gn_167 - 336 Gn_168 - 338 Gn_169
 - 34 Gn_17 - 340 Gn_170 - 342 Gn_171 - 344 Gn_172 - 346 Gn_173 - 348 Gn_174
 - 350 Gn_175 - 352 Gn_176 - 354 Gn_177 - 356 Gn_178 - 358 Gn_179 - 36 Gn_18
 - 360 Gn_180 - 362 Gn_181 - 364 Gn_182 - 366 Gn_183 - 368 Gn_184 - 370 Gn_185
 - 372 Gn_186 - 374 Gn_187 - 376 Gn_188 - 378 Gn_189 - 38 Gn_19 - 380 Gn_190
 - 382 Gn_191 - 384 Gn_192 - 386 Gn_193 - 388 Gn_194 - 390 Gn_195 - 392 Gn_196
 - 394 Gn_197 - 396 Gn_198 - 398 Gn_199 - 4 Gn_2 - 40 Gn_20 - 400 Gn_200
 - 402 Gn_201 - 404 Gn_202 - 406 Gn_203 - 408 Gn_204 - 410 Gn_205 - 412 Gn_206
 - 414 Gn_207 - 416 Gn_208 - 418 Gn_209 - 42 Gn_21 - 420 Gn_210 - 422 Gn_211
 - 424 Gn_212 - 426 Gn_213 - 428 Gn_214 - 430 Gn_215 - 432 Gn_216 - 434 Gn_217
 - 436 Gn_218 - 438 Gn_219 - 44 Gn_22 - 440 Gn_220 - 442 Gn_221 - 444 Gn_222
 - 446 Gn_223 - 448 Gn_224 - 450 Gn_225 - 452 Gn_226 - 454 Gn_227 - 456 Gn_228
 - 458 Gn_229 - 46 Gn_23 - 460 Gn_230 - 462 Gn_231 - 464 Gn_232 - 466 Gn_233
 - 468 Gn_234 - 470 Gn_235 - 472 Gn_236 - 474 Gn_237 - 476 Gn_238 - 478 Gn_239
 - 48 Gn_24 - 480 Gn_240 - 482 Gn_241 - 484 Gn_242 - 486 Gn_243 - 488 Gn_244
 - 490 Gn_245 - 492 Gn_246 - 494 Gn_247 - 496 Gn_248 - 498 Gn_249 - 50 Gn_25
 - 52 Gn_26 - 54 Gn_27 - 56 Gn_28 - 58 Gn_29 - 6 Gn_3 - 60 Gn_30 - 62 Gn_31
 - 64 Gn_32 - 66 Gn_33 - 68 Gn_34 - 70 Gn_35 - 72 Gn_36 - 74 Gn_37 - 76 Gn_38
 - 78 Gn_39 - 8 Gn_4 - 80 Gn_40 - 82 Gn_41 - 84 Gn_42 - 86 Gn_43 - 88 Gn_44
 - 90 Gn_45 - 92 Gn_46 - 94 Gn_47 - 96 Gn_48 - 98 Gn_49 - 10 Gn_5 - 100 Gn_50
 - 102 Gn_51 - 104 Gn_52 - 106 Gn_53 - 108 Gn_54 - 110 Gn_55 - 112 Gn_56
 - 114 Gn_57 - 116 Gn_58 - 118 Gn_59 - 12 Gn_6 - 120 Gn_60 - 122 Gn_61
 - 124 Gn_62 - 126 Gn_63 - 128 Gn_64 - 130 Gn_65 - 132 Gn_66 - 134 Gn_67
 - 136 Gn_68 - 138 Gn_69 - 14 Gn_7 - 140 Gn_70 - 142 Gn_71 - 144 Gn_72
 - 146 Gn_73 - 148 Gn_74 - 150 Gn_75 - 152 Gn_76 - 154 Gn_77 - 156 Gn_78
 - 158 Gn_79 - 16 Gn_8 - 160 Gn_80 - 162 Gn_81 - 164 Gn_82 - 166 Gn_83
 - 168 Gn_84 - 170 Gn_85 - 172 Gn_86 - 174 Gn_87 - 176 Gn_88 - 178 Gn_89
 - 18 Gn_9 - 180 Gn_90 - 182 Gn_91 - 184 Gn_92 - 186 Gn_93 - 188 Gn_94
 - 190 Gn_95 - 192 Gn_96 - 194 Gn_97 - 196 Gn_98 - 198 Gn_99 + Y_pump
 + 0.5 d_pump <= 0
_C443: Gn_1 + Gn_10 + Gn_100 + Gn_101 + Gn_102 + Gn_103 + Gn_104 + Gn_105
 + Gn_106 + Gn_107 + Gn_108 + Gn_109 + Gn_11 + Gn_110 + Gn_111 + Gn_112
 + Gn_113 + Gn_114 + Gn_115 + Gn_116 + Gn_117 + Gn_118 + Gn_119 + Gn_12
 + Gn_120 + Gn_121 + Gn_122 + Gn_123 + Gn_124 + Gn_125 + Gn_126 + Gn_127
 + Gn_128 + Gn_129 + Gn_13 + Gn_130 + Gn_131 + Gn_132 + Gn_133 + Gn_134
 + Gn_135 + Gn_136 + Gn_137 + Gn_138 + Gn_139 + Gn_14 + Gn_140 + Gn_141
 + Gn_142 + Gn_143 + Gn_144 + Gn_145 + Gn_146 + Gn_147 + Gn_148 + Gn_149
 + Gn_15 + Gn_150 + Gn_151 + Gn_152 + Gn_153 + Gn_154 + Gn_155 + Gn_156
 + Gn_157 + Gn_158 + Gn_159 + Gn_16 + Gn_160 + Gn_161 + Gn_162 + Gn_163
 + Gn_164 + Gn_165 + Gn_166 + Gn_167 + Gn_168 + Gn_169 + Gn_17 + Gn_170
 + Gn_171 + Gn_172 + Gn_173 + Gn_174 + Gn_175 + Gn_176 + Gn_177 + Gn_178
 + Gn_179 + Gn_18 + Gn_180 + Gn_181 + Gn_182 + Gn_183 + Gn_184 + Gn_185
 + Gn_186 + Gn_187 + Gn_188 + Gn_189 + Gn_19 + Gn_190 + Gn_191 + Gn_192
 + Gn_193 + Gn_194 + Gn_195 + Gn_196 + Gn_197 + Gn_198 + Gn_199 + Gn_2 + Gn_20
 + Gn_200 + Gn_201 + Gn_202 + Gn_203 + Gn_204 + Gn_205 + Gn_206 + Gn_207
 + Gn_208 + Gn_209 + Gn_21 + Gn_210 + Gn_211 + Gn_212 + Gn_213 + Gn_214
 + Gn_215 + Gn_216 + Gn_217 + Gn_218 + Gn_219 + Gn_22 + Gn_220 + Gn_221
 + Gn_222 + Gn_223 + Gn_224 + Gn_225 + Gn_226 + Gn_227 + Gn_228 + Gn_229
 + Gn_23 + Gn_230 + Gn_231 + Gn_232 + Gn_233 + Gn_234 + Gn_235 + Gn_236
 + Gn_237 + Gn_238 + Gn_239 + Gn_24 + Gn_240 + Gn_241 + Gn_242 + Gn_243
 + Gn_244 + Gn_245 + Gn_246 + Gn_247 + Gn_248 + Gn_249 + Gn_25 + Gn_26 + Gn_27
 + Gn_28 + Gn_29 + Gn_3 + Gn_30 + Gn_31 + Gn_32 + Gn_33 + Gn_34 + Gn_35
 + Gn_36 + Gn_37 + Gn_38 + Gn_39 + Gn_4 + Gn_40 + Gn_41 + Gn_42 + Gn_43
 + Gn_44 + Gn_45 + Gn_46 + Gn_47 + Gn_48 + Gn_49 + Gn_5 + Gn_50 + Gn_51
 + Gn_52 + Gn_53 + Gn_54 + Gn_55 + Gn_56 + Gn_57 + Gn_58 + Gn_59 + Gn_6
 + Gn_60 + Gn_61 + Gn_62 + Gn_63 + Gn_64 + Gn_65 + Gn_66 + Gn_67 + Gn_68
 + Gn_69 + Gn_7 + Gn_70 + Gn_71 + Gn_72 + Gn_73 + Gn_74 + Gn_75 + Gn_76
 + Gn_77 + Gn_78 + Gn_79 + Gn_8 + Gn_80 + Gn_81 + Gn_82 + Gn_83 + Gn_84
 + Gn_85 + Gn_86 + Gn_87 + Gn_88 + Gn_89 + Gn_9 + Gn_90 + Gn_91 + Gn_92
 + Gn_93 + Gn_94 + Gn_95 + Gn_96 + Gn_97 + Gn_98 + Gn_99 = 1
_C444: - 12000 Gn_1 - 1200000 Gn_10 - 120000000 Gn_100 - 122412000 Gn_101
 - 124848000 Gn_102 - 127308000 Gn_103 - 129792000 Gn_104 - 132300000 Gn_105
 - 134832000 Gn_106 - 137388000 Gn_107 - 139968000 Gn_108 - 142572000 Gn_109
 - 1452000 Gn_11 - 145200000 Gn_110 - 147852000 Gn_111 - 150528000 Gn_112
 - 153228000 Gn_113 - 155952000 Gn_114 - 158700000 Gn_115 - 161472000 Gn_116
 - 164268000 Gn_117 - 167088000 Gn_118 - 169932000 Gn_119 - 1728000 Gn_12
 - 172800000 Gn_120 - 175692000 Gn_121 - 178608000 Gn_122 - 181548000 Gn_123
 - 184512000 Gn_124 - 187500000 Gn_125 - 190512000 Gn_126 - 193548000 Gn_127
 - 196608000 Gn_128 - 199692000 Gn_129 - 2028000 Gn_13 - 202800000 Gn_130
 - 205932000 Gn_131 - 209088000 Gn_132 - 212268000 Gn_133 - 215472000 Gn_134
 - 218700000 Gn_135 - 221952000 Gn_136 - 225228000 Gn_137 - 228528000 Gn_138
 - 231852000 Gn_139 - 2352000 Gn_14 - 235200000 Gn_140 - 238572000 Gn_141
 - 241968000 Gn_142 - 245388000 Gn_143 - 248832000 Gn_144 - 252300000 Gn_145
 - 255792000 Gn_146 - 259308000 Gn_147 - 262848000 Gn_148 - 266412000 Gn_149
 - 2700000 Gn_15 - 270000000 Gn_150 - 273612000 Gn_151 - 277248000 Gn_152
 - 280908000 Gn_153 - 284592000 Gn_154 - 288300000 Gn_155 - 292032000 Gn_156
 - 295788000 Gn_157 - 299568000 Gn_158 - 303372000 Gn_159 - 3072000 Gn_16
 - 307200000 Gn_160 - 311052000 Gn_161 - 314928000 Gn_162 - 318828000 Gn_163
 - 322752000 Gn_164 - 326700000 Gn_165 - 330672000 Gn_166 - 334668000 Gn_167
 - 338688000 Gn_168 - 342732000 Gn_169 - 3468000 Gn_17 - 346800000 Gn_170
 - 350892000 Gn_171 - 355008000 Gn_172 - 359148000 Gn_173 - 363312000 Gn_174
 - 367500000 Gn_175 - 371712000 Gn_176 - 375948000 Gn_177 - 380208000 Gn_178
 - 384492000 Gn_179 - 3888000 Gn_18 - 388800000 Gn_180 - 393132000 Gn_181
 - 397488000 Gn_182 - 401868000 Gn_183 - 406272000 Gn_184 - 410700000 Gn_185
 - 415152000 Gn_186 - 419628000 Gn_187 - 424128000 Gn_188 - 428652000 Gn_189
 - 4332000 Gn_19 - 433200000 Gn_190 - 437772000 Gn_191 - 442368000 Gn_192
 - 446988000 Gn_193 - 451632000 Gn_194 - 456300000 Gn_195 - 460992000 Gn_196
 - 465708000 Gn_197 - 470448000 Gn_198 - 475212000 Gn_199 - 48000 Gn_2
 - 4800000 Gn_20 - 480000000 Gn_200 - 484812000 Gn_201 - 489648000 Gn_202
 - 494508000 Gn_203 - 499392000 Gn_204 - 504300000 Gn_205 - 509232000 Gn_206
 - 514188000 Gn_207 - 519168000 Gn_208 - 524172000 Gn_209 - 5292000 Gn_21
 - 529200000 Gn_210 - 534252000 Gn_211 - 539328000 Gn_212 - 544428000 Gn_213
 - 549552000 Gn_214 - 554700000 Gn_215 - 559872000 Gn_216 - 565068000 Gn_217
 - 570288000 Gn_218 - 575532000 Gn_219 - 5808000 Gn_22 - 580800000 Gn_220
 - 586092000 Gn_221 - 591408000 Gn_222 - 596748000 Gn_223 - 602112000 Gn_224
 - 607500000 Gn_225 - 612912000 Gn_226 - 618348000 Gn_227 - 623808000 Gn_228
 - 629292000 Gn_229 - 6348000 Gn_23 - 634800000 Gn_230 - 640332000 Gn_231
 - 645888000 Gn_232 - 651468000 Gn_233 - 657072000 Gn_234 - 662700000 Gn_235
 - 668352000 Gn_236 - 674028000 Gn_237 - 679728000 Gn_238 - 685452000 Gn_239
 - 6912000 Gn_24 - 691200000 Gn_240 - 696972000 Gn_241 - 702768000 Gn_242
 - 708588000 Gn_243 - 714432000 Gn_244 - 720300000 Gn_245 - 726192000 Gn_246
 - 732108000 Gn_247 - 738048000 Gn_248 - 744012000 Gn_249 - 7500000 Gn_25
 - 8112000 Gn_26 - 8748000 Gn_27 - 9408000 Gn_28 - 10092000 Gn_29
 - 108000 Gn_3 - 10800000 Gn_30 - 11532000 Gn_31 - 12288000 Gn_32
 - 13068000 Gn_33 - 13872000 Gn_34 - 14700000 Gn_35 - 15552000 Gn_36
 - 16428000 Gn_37 - 17328000 Gn_38 - 18252000 Gn_39 - 192000 Gn_4
 - 19200000 Gn_40 - 20172000 Gn_41 - 21168000 Gn_42 - 22188000 Gn_43
 - 23232000 Gn_44 - 24300000 Gn_45 - 25392000 Gn_46 - 26508000 Gn_47
 - 27648000 Gn_48 - 28812000 Gn_49 - 300000 Gn_5 - 30000000 Gn_50
 - 31212000 Gn_51 - 32448000 Gn_52 - 33708000 Gn_53 - 34992000 Gn_54
 - 36300000 Gn_55 - 37632000 Gn_56 - 38988000 Gn_57 - 40368000 Gn_58
 - 41772000 Gn_59 - 432000 Gn_6 - 43200000 Gn_60 - 44652000 Gn_61
 - 46128000 Gn_62 - 47628000 Gn_63 - 49152000 Gn_64 - 50700000 Gn_65
 - 52272000 Gn_66 - 53868000 Gn_67 - 55488000 Gn_68 - 57132000 Gn_69
 - 588000 Gn_7 - 58800000 Gn_70 - 60492000 Gn_71 - 62208000 Gn_72
 - 63948000 Gn_73 - 65712000 Gn_74 - 67500000 Gn_75 - 69312000 Gn_76
 - 71148000 Gn_77 - 73008000 Gn_78 - 74892000 Gn_79 - 768000 Gn_8
 - 76800000 Gn_80 - 78732000 Gn_81 - 80688000 Gn_82 - 82668000 Gn_83
 - 84672000 Gn_84 - 86700000 Gn_85 - 88752000 Gn_86 - 90828000 Gn_87
 - 92928000 Gn_88 - 95052000 Gn_89 - 972000 Gn_9 - 97200000 Gn_90
 - 99372000 Gn_91 - 101568000 Gn_92 - 103788000 Gn_93 - 106032000 Gn_94
 - 108300000 Gn_95 - 110592000 Gn_96 - 112908000 Gn_97 - 115248000 Gn_98
 - 117612000 Gn_99 + TLC = 0
_C445: Ve2_reactor_reactor = 0
_C446: - D_reactor_hex1 + Din_reactor_hex1 + Dout_reactor_hex1 = 0
_C447: Din_reactor_hex1 - 40 Psie_reactor_hex1 <= 0
_C448: Dout_reactor_hex1 + 40 Psie_reactor_hex1 >= 40
_C449: Dout_reactor_hex1 + 1000 Psie_reactor_hex1 <= 1000
_C45: - A_hex1_reactor + B_reactor_hex1 = 0
_C450: 12775 Din_reactor_hex1 - 511000 Psie_reactor_hex1 + Ve2_reactor_hex1
 = 0
_C451: - D_reactor_eoabs + Din_reactor_eoabs + Dout_reactor_eoabs = 0
_C452: Din_reactor_eoabs - 40 Psie_reactor_eoabs <= 0
_C453: Dout_reactor_eoabs + 40 Psie_reactor_eoabs >= 40
_C454: Dout_reactor_eoabs + 1000 Psie_reactor_eoabs <= 1000
_C455: 77675 Din_reactor_eoabs - 3107000 Psie_reactor_eoabs
 + Ve2_reactor_eoabs = 0
_C456: - D_reactor_hex2 + Din_reactor_hex2 + Dout_reactor_hex2 = 0
_C457: Din_reactor_hex2 - 40 Psie_reactor_hex2 <= 0
_C458: Dout_reactor_hex2 + 40 Psie_reactor_hex2 >= 40
_C459: Dout_reactor_hex2 + 1000 Psie_reactor_hex2 <= 1000
_C46: Wx_hex1_reactor + Wx_reactor_hex1 = 1
_C460: 12600 Din_reactor_hex2 - 504000 Psie_reactor_hex2 + Ve2_reactor_hex2
 = 0
_C461: - D_reactor_co2abs + Din_reactor_co2abs + Dout_reactor_co2abs = 0
_C462: Din_reactor_co2abs - 40 Psie_reactor_co2abs <= 0
_C463: Dout_reactor_co2abs + 40 Psie_reactor_co2abs >= 40
_C464: Dout_reactor_co2abs + 1000 Psie_reactor_co2abs <= 1000
_C465: 77032.5 Din_reactor_co2abs - 3081300 Psie_reactor_co2abs
 + Ve2_reactor_co2abs = 0
_C466: - D_reactor_flash + Din_reactor_flash + Dout_reactor_flash = 0
_C467: Din_reactor_flash - 40 Psie_reactor_flash <= 0
_C468: Dout_reactor_flash + 40 Psie_reactor_flash >= 40
_C469: Dout_reactor_flash + 1000 Psie_reactor_flash <= 1000
_C47: Wy_hex1_reactor + Wy_reactor_hex1 = 1
_C470: 12625 Din_reactor_flash - 505000 Psie_reactor_flash + Ve2_reactor_flash
 = 0
_C471: - D_reactor_pump + Din_reactor_pump + Dout_reactor_pump = 0
_C472: Din_reactor_pump - 40 Psie_reactor_pump <= 0
_C473: Dout_reactor_pump + 40 Psie_reactor_pump >= 40
_C474: Dout_reactor_pump + 1000 Psie_reactor_pump <= 1000
_C475: 12537.5 Din_reactor_pump - 501500 Psie_reactor_pump + Ve2_reactor_pump
 = 0
_C476: - D_eoabs_reactor + Din_eoabs_reactor + Dout_eoabs_reactor = 0
_C477: Din_eoabs_reactor - 21.8 Psie_eoabs_reactor <= 0
_C478: Dout_eoabs_reactor + 21.8 Psie_eoabs_reactor >= 21.8
_C479: Dout_eoabs_reactor + 1000 Psie_eoabs_reactor <= 1000
_C48: CD_reactor_eoabs = 0
_C480: 152981.651376 Din_eoabs_reactor - 3335000 Psie_eoabs_reactor
 + Ve2_eoabs_reactor = 0
_C481: - D_eoabs_hex1 + Din_eoabs_hex1 + Dout_eoabs_hex1 = 0
_C482: Din_eoabs_hex1 - 21.8 Psie_eoabs_hex1 <= 0
_C483: Dout_eoabs_hex1 + 21.8 Psie_eoabs_hex1 >= 21.8
_C484: Dout_eoabs_hex1 + 1000 Psie_eoabs_hex1 <= 1000
_C485: 23440.3669725 Din_eoabs_hex1 - 511000 Psie_eoabs_hex1 + Ve2_eoabs_hex1
 = 0
_C486: Ve2_eoabs_eoabs = 0
_C487: - D_eoabs_hex2 + Din_eoabs_hex2 + Dout_eoabs_hex2 = 0
_C488: Din_eoabs_hex2 - 21.8 Psie_eoabs_hex2 <= 0
_C489: Dout_eoabs_hex2 + 21.8 Psie_eoabs_hex2 >= 21.8
_C49: - L_reactor_eoabs + R_reactor_eoabs + X_eoabs - X_reactor = 0
_C490: Dout_eoabs_hex2 + 1000 Psie_eoabs_hex2 <= 1000
_C491: 23119.266055 Din_eoabs_hex2 - 504000 Psie_eoabs_hex2 + Ve2_eoabs_hex2
 = 0
_C492: - D_eoabs_co2abs + Din_eoabs_co2abs + Dout_eoabs_co2abs = 0
_C493: Din_eoabs_co2abs - 21.8 Psie_eoabs_co2abs <= 0
_C494: Dout_eoabs_co2abs + 21.8 Psie_eoabs_co2abs >= 21.8
_C495: Dout_eoabs_co2abs + 1000 Psie_eoabs_co2abs <= 1000
_C496: 141344.036697 Din_eoabs_co2abs - 3081300 Psie_eoabs_co2abs
 + Ve2_eoabs_co2abs = 0
_C497: - D_eoabs_flash + Din_eoabs_flash + Dout_eoabs_flash = 0
_C498: Din_eoabs_flash - 21.8 Psie_eoabs_flash <= 0
_C499: Dout_eoabs_flash + 21.8 Psie_eoabs_flash >= 21.8
_C5: 0 O_hex1 + l_hex1 = 11.42
_C50: A_reactor_eoabs - B_reactor_eoabs + Y_eoabs - Y_reactor = 0
_C500: Dout_eoabs_flash + 1000 Psie_eoabs_flash <= 1000
_C501: 23165.1376147 Din_eoabs_flash - 505000 Psie_eoabs_flash
 + Ve2_eoabs_flash = 0
_C502: - D_eoabs_pump + Din_eoabs_pump + Dout_eoabs_pump = 0
_C503: Din_eoabs_pump - 21.8 Psie_eoabs_pump <= 0
_C504: Dout_eoabs_pump + 21.8 Psie_eoabs_pump >= 21.8
_C505: Dout_eoabs_pump + 1000 Psie_eoabs_pump <= 1000
_C506: 23004.587156 Din_eoabs_pump - 501500 Psie_eoabs_pump + Ve2_eoabs_pump
 = 0
_C507: - D_co2abs_reactor + Din_co2abs_reactor + Dout_co2abs_reactor = 0
_C508: Din_co2abs_reactor - 18.06 Psie_co2abs_reactor <= 0
_C509: Dout_co2abs_reactor + 18.06 Psie_co2abs_reactor >= 18.06
_C51: R_reactor_eoabs - 1000 Wx_reactor_eoabs <= 0
_C510: Dout_co2abs_reactor + 1000 Psie_co2abs_reactor <= 1000
_C511: 184662.236988 Din_co2abs_reactor - 3335000 Psie_co2abs_reactor
 + Ve2_co2abs_reactor = 0
_C512: - D_co2abs_hex1 + Din_co2abs_hex1 + Dout_co2abs_hex1 = 0
_C513: Din_co2abs_hex1 - 18.06 Psie_co2abs_hex1 <= 0
_C514: Dout_co2abs_hex1 + 18.06 Psie_co2abs_hex1 >= 18.06
_C515: Dout_co2abs_hex1 + 1000 Psie_co2abs_hex1 <= 1000
_C516: 28294.5736434 Din_co2abs_hex1 - 511000 Psie_co2abs_hex1
 + Ve2_co2abs_hex1 = 0
_C517: - D_co2abs_eoabs + Din_co2abs_eoabs + Dout_co2abs_eoabs = 0
_C518: Din_co2abs_eoabs - 18.06 Psie_co2abs_eoabs <= 0
_C519: Dout_co2abs_eoabs + 18.06 Psie_co2abs_eoabs >= 18.06
_C52: L_reactor_eoabs + 1000 Wx_reactor_eoabs <= 1000
_C520: Dout_co2abs_eoabs + 1000 Psie_co2abs_eoabs <= 1000
_C521: 172037.65227 Din_co2abs_eoabs - 3107000 Psie_co2abs_eoabs
 + Ve2_co2abs_eoabs = 0
_C522: - D_co2abs_hex2 + Din_co2abs_hex2 + Dout_co2abs_hex2 = 0
_C523: Din_co2abs_hex2 - 18.06 Psie_co2abs_hex2 <= 0
_C524: Dout_co2abs_hex2 + 18.06 Psie_co2abs_hex2 >= 18.06
_C525: Dout_co2abs_hex2 + 1000 Psie_co2abs_hex2 <= 1000
_C526: 27906.9767442 Din_co2abs_hex2 - 504000 Psie_co2abs_hex2
 + Ve2_co2abs_hex2 = 0
_C527: Ve2_co2abs_co2abs = 0
_C528: - D_co2abs_flash + Din_co2abs_flash + Dout_co2abs_flash = 0
_C529: Din_co2abs_flash - 18.06 Psie_co2abs_flash <= 0
_C53: A_reactor_eoabs - 1000 Wy_reactor_eoabs <= 0
_C530: Dout_co2abs_flash + 18.06 Psie_co2abs_flash >= 18.06
_C531: Dout_co2abs_flash + 1000 Psie_co2abs_flash <= 1000
_C532: 27962.3477298 Din_co2abs_flash - 505000 Psie_co2abs_flash
 + Ve2_co2abs_flash = 0
_C533: - D_co2abs_pump + Din_co2abs_pump + Dout_co2abs_pump = 0
_C534: Din_co2abs_pump - 18.06 Psie_co2abs_pump <= 0
_C535: Dout_co2abs_pump + 18.06 Psie_co2abs_pump >= 18.06
_C536: Dout_co2abs_pump + 1000 Psie_co2abs_pump <= 1000
_C537: 27768.5492802 Din_co2abs_pump - 501500 Psie_co2abs_pump
 + Ve2_co2abs_pump = 0
_C538: - Ve2_reactor_co2abs - Ve2_reactor_eoabs - Ve2_reactor_flash
 - Ve2_reactor_hex1 - Ve2_reactor_hex2 - Ve2_reactor_pump
 - Ve2_reactor_reactor + Ve_reactor = 3335000
_C539: Omega0_reactor - 0.87 Ve_reactor = 0
_C54: B_reactor_eoabs + 1000 Wy_reactor_eoabs <= 1000
_C540: - Ve2_eoabs_co2abs - Ve2_eoabs_eoabs - Ve2_eoabs_flash - Ve2_eoabs_hex1
 - Ve2_eoabs_hex2 - Ve2_eoabs_pump - Ve2_eoabs_reactor + Ve_eoabs = 3107000
_C541: Omega0_eoabs - 0.73 Ve_eoabs = 0
_C542: - Ve2_co2abs_co2abs - Ve2_co2abs_eoabs - Ve2_co2abs_flash
 - Ve2_co2abs_hex1 - Ve2_co2abs_hex2 - Ve2_co2abs_pump - Ve2_co2abs_reactor
 + Ve_co2abs = 3081300
_C543: Omega0_co2abs - 0.66 Ve_co2abs = 0
_C544: - 100000000 Z_reactor_1 + linOmega0_reactor_1 <= 0
_C545: - 100000000 Z_reactor_2 + linOmega0_reactor_2 <= 0
_C546: - 100000000 Z_reactor_3 + linOmega0_reactor_3 <= 0
_C547: - 100000000 Z_reactor_4 + linOmega0_reactor_4 <= 0
_C548: - 100000000 Z_reactor_5 + linOmega0_reactor_5 <= 0
_C549: - 100000000 Z_reactor_6 + linOmega0_reactor_6 <= 0
_C55: - A_reactor_eoabs - B_reactor_eoabs + D_reactor_eoabs - L_reactor_eoabs
 - R_reactor_eoabs = 0
_C550: - 100000000 Z_eoabs_1 + linOmega0_eoabs_1 <= 0
_C551: - 100000000 Z_eoabs_2 + linOmega0_eoabs_2 <= 0
_C552: - 100000000 Z_eoabs_3 + linOmega0_eoabs_3 <= 0
_C553: - 100000000 Z_eoabs_4 + linOmega0_eoabs_4 <= 0
_C554: - 100000000 Z_eoabs_5 + linOmega0_eoabs_5 <= 0
_C555: - 100000000 Z_eoabs_6 + linOmega0_eoabs_6 <= 0
_C556: - 100000000 Z_co2abs_1 + linOmega0_co2abs_1 <= 0
_C557: - 100000000 Z_co2abs_2 + linOmega0_co2abs_2 <= 0
_C558: - 100000000 Z_co2abs_3 + linOmega0_co2abs_3 <= 0
_C559: - 100000000 Z_co2abs_4 + linOmega0_co2abs_4 <= 0
_C56: - D_eoabs_reactor + D_reactor_eoabs = 0
_C560: - 100000000 Z_co2abs_5 + linOmega0_co2abs_5 <= 0
_C561: - 100000000 Z_co2abs_6 + linOmega0_co2abs_6 <= 0
_C562: Z_reactor_1 + Z_reactor_2 + Z_reactor_3 + Z_reactor_4 + Z_reactor_5
 + Z_reactor_6 = 1
_C563: Omega0_reactor - linOmega0_reactor_1 - linOmega0_reactor_2
 - linOmega0_reactor_3 - linOmega0_reactor_4 - linOmega0_reactor_5
 - linOmega0_reactor_6 = 0
_C564: Omega_reactor - linOmega0_reactor_1 - 0.9 linOmega0_reactor_2
 - 0.75 linOmega0_reactor_3 - 0.365 linOmega0_reactor_4
 - 0.292 linOmega0_reactor_5 - 0.117 linOmega0_reactor_6 = 0
_C565: - Omega_co2abs - Omega_eoabs - Omega_reactor + SumOmega = 0
_C566: SumPZ - 5000 Z_co2abs_2 - 20000 Z_co2abs_3 - 25000 Z_co2abs_4
 - 35000 Z_co2abs_5 - 55000 Z_co2abs_6 - 5000 Z_eoabs_2 - 20000 Z_eoabs_3
 - 25000 Z_eoabs_4 - 35000 Z_eoabs_5 - 55000 Z_eoabs_6 - 5000 Z_reactor_2
 - 15000 Z_reactor_3 - 40000 Z_reactor_4 - 60000 Z_reactor_5
 - 125000 Z_reactor_6 = 0
_C567: Z_eoabs_1 + Z_eoabs_2 + Z_eoabs_3 + Z_eoabs_4 + Z_eoabs_5 + Z_eoabs_6
 = 1
_C568: Omega0_eoabs - linOmega0_eoabs_1 - linOmega0_eoabs_2
 - linOmega0_eoabs_3 - linOmega0_eoabs_4 - linOmega0_eoabs_5
 - linOmega0_eoabs_6 = 0
_C569: Omega_eoabs - linOmega0_eoabs_1 - 0.9 linOmega0_eoabs_2
 - 0.76 linOmega0_eoabs_3 - 0.684 linOmega0_eoabs_4 - 0.612 linOmega0_eoabs_5
 - 0.465 linOmega0_eoabs_6 = 0
_C57: 1000 E1_reactor_eoabs + 1000 E2_reactor_eoabs - X_eoabs + X_reactor
 - 0.5 l_eoabs - 0.5 l_reactor >= 0
_C570: - Omega_co2abs - Omega_eoabs - Omega_reactor + SumOmega = 0
_C571: SumPZ - 5000 Z_co2abs_2 - 20000 Z_co2abs_3 - 25000 Z_co2abs_4
 - 35000 Z_co2abs_5 - 55000 Z_co2abs_6 - 5000 Z_eoabs_2 - 20000 Z_eoabs_3
 - 25000 Z_eoabs_4 - 35000 Z_eoabs_5 - 55000 Z_eoabs_6 - 5000 Z_reactor_2
 - 15000 Z_reactor_3 - 40000 Z_reactor_4 - 60000 Z_reactor_5
 - 125000 Z_reactor_6 = 0
_C572: Z_co2abs_1 + Z_co2abs_2 + Z_co2abs_3 + Z_co2abs_4 + Z_co2abs_5
 + Z_co2abs_6 = 1
_C573: Omega0_co2abs - linOmega0_co2abs_1 - linOmega0_co2abs_2
 - linOmega0_co2abs_3 - linOmega0_co2abs_4 - linOmega0_co2abs_5
 - linOmega0_co2abs_6 = 0
_C574: Omega_co2abs - linOmega0_co2abs_1 - 0.9 linOmega0_co2abs_2
 - 0.76 linOmega0_co2abs_3 - 0.684 linOmega0_co2abs_4
 - 0.612 linOmega0_co2abs_5 - 0.465 linOmega0_co2abs_6 = 0
_C575: - Omega_co2abs - Omega_eoabs - Omega_reactor + SumOmega = 0
_C576: SumPZ - 5000 Z_co2abs_2 - 20000 Z_co2abs_3 - 25000 Z_co2abs_4
 - 35000 Z_co2abs_5 - 55000 Z_co2abs_6 - 5000 Z_eoabs_2 - 20000 Z_eoabs_3
 - 25000 Z_eoabs_4 - 35000 Z_eoabs_5 - 55000 Z_eoabs_6 - 5000 Z_reactor_2
 - 15000 Z_reactor_3 - 40000 Z_reactor_4 - 60000 Z_reactor_5
 - 125000 Z_reactor_6 = 0
_C58: - 1000 E1_reactor_eoabs + 1000 E2_reactor_eoabs + X_eoabs - X_reactor
 - 0.5 l_eoabs - 0.5 l_reactor >= -1000
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
_C67: CD_hex1_eoabs - 98.4 D_hex1_eoabs = 0
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
 - 0.5 l_hex1 >= 0
_C77: - 1000 E1_hex1_eoabs + 1000 E2_hex1_eoabs + X_eoabs - X_hex1
 - 0.5 l_eoabs - 0.5 l_hex1 >= -1000
_C78: 1000 E1_hex1_eoabs - 1000 E2_hex1_eoabs - Y_eoabs + Y_hex1 - 0.5 d_eoabs
 - 0.5 d_hex1 >= -1000
_C79: - 1000 E1_hex1_eoabs - 1000 E2_hex1_eoabs + Y_eoabs - Y_hex1
 - 0.5 d_eoabs - 0.5 d_hex1 >= -2000
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
Gn_1
Gn_10
Gn_100
Gn_101
Gn_102
Gn_103
Gn_104
Gn_105
Gn_106
Gn_107
Gn_108
Gn_109
Gn_11
Gn_110
Gn_111
Gn_112
Gn_113
Gn_114
Gn_115
Gn_116
Gn_117
Gn_118
Gn_119
Gn_12
Gn_120
Gn_121
Gn_122
Gn_123
Gn_124
Gn_125
Gn_126
Gn_127
Gn_128
Gn_129
Gn_13
Gn_130
Gn_131
Gn_132
Gn_133
Gn_134
Gn_135
Gn_136
Gn_137
Gn_138
Gn_139
Gn_14
Gn_140
Gn_141
Gn_142
Gn_143
Gn_144
Gn_145
Gn_146
Gn_147
Gn_148
Gn_149
Gn_15
Gn_150
Gn_151
Gn_152
Gn_153
Gn_154
Gn_155
Gn_156
Gn_157
Gn_158
Gn_159
Gn_16
Gn_160
Gn_161
Gn_162
Gn_163
Gn_164
Gn_165
Gn_166
Gn_167
Gn_168
Gn_169
Gn_17
Gn_170
Gn_171
Gn_172
Gn_173
Gn_174
Gn_175
Gn_176
Gn_177
Gn_178
Gn_179
Gn_18
Gn_180
Gn_181
Gn_182
Gn_183
Gn_184
Gn_185
Gn_186
Gn_187
Gn_188
Gn_189
Gn_19
Gn_190
Gn_191
Gn_192
Gn_193
Gn_194
Gn_195
Gn_196
Gn_197
Gn_198
Gn_199
Gn_2
Gn_20
Gn_200
Gn_201
Gn_202
Gn_203
Gn_204
Gn_205
Gn_206
Gn_207
Gn_208
Gn_209
Gn_21
Gn_210
Gn_211
Gn_212
Gn_213
Gn_214
Gn_215
Gn_216
Gn_217
Gn_218
Gn_219
Gn_22
Gn_220
Gn_221
Gn_222
Gn_223
Gn_224
Gn_225
Gn_226
Gn_227
Gn_228
Gn_229
Gn_23
Gn_230
Gn_231
Gn_232
Gn_233
Gn_234
Gn_235
Gn_236
Gn_237
Gn_238
Gn_239
Gn_24
Gn_240
Gn_241
Gn_242
Gn_243
Gn_244
Gn_245
Gn_246
Gn_247
Gn_248
Gn_249
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
Gn_60
Gn_61
Gn_62
Gn_63
Gn_64
Gn_65
Gn_66
Gn_67
Gn_68
Gn_69
Gn_7
Gn_70
Gn_71
Gn_72
Gn_73
Gn_74
Gn_75
Gn_76
Gn_77
Gn_78
Gn_79
Gn_8
Gn_80
Gn_81
Gn_82
Gn_83
Gn_84
Gn_85
Gn_86
Gn_87
Gn_88
Gn_89
Gn_9
Gn_90
Gn_91
Gn_92
Gn_93
Gn_94
Gn_95
Gn_96
Gn_97
Gn_98
Gn_99
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
