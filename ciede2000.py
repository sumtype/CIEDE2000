import numpy as np

# Converts RGB pixel array to XYZ format.
# Implementation derived from http://www.easyrgb.com/en/math.php
def rgb2xyz(rgb):
    def format(c):
        c = c / 255.
        if c > 0.04045: c = ((c + 0.055) / 1.055) ** 2.4
        else: c = c / 12.92
        return c * 100
    rgb = list(map(format, rgb))
    xyz = [None, None, None]
    xyz[0] = rgb[0] * 0.4124 + rgb[1] * 0.3576 + rgb[2] * 0.1805
    xyz[1] = rgb[0] * 0.2126 + rgb[1] * 0.7152 + rgb[2] * 0.0722
    xyz[2] = rgb[0] * 0.0193 + rgb[1] * 0.1192 + rgb[2] * 0.9505
    return xyz

# Converts XYZ pixel array to LAB format.
# Implementation derived from http://www.easyrgb.com/en/math.php
def xyz2lab(xyz):
    def format(c):
        if c > 0.008856: c = c ** (1. / 3.)
        else: c = (7.787 * c) + (16. / 116.)
        return c
    xyz[0] = xyz[0] / 95.047
    xyz[1] = xyz[1] / 100.00
    xyz[2] = xyz[2] / 108.883
    xyz = list(map(format, xyz))
    lab = [None, None, None]
    lab[0] = (116. * xyz[1]) - 16.
    lab[1] = 500. * (xyz[0] - xyz[1])
    lab[2] = 200. * (xyz[1] - xyz[2])
    return lab

# Converts RGB pixel array into LAB format.
def rgb2lab(rgb):
    return xyz2lab(rgb2xyz(rgb))

# Returns CIEDE2000 comparison results of two LAB formatted colors.
# Implementation derived from the excel spreadsheet provided here: http://www2.ece.rochester.edu/~gsharma/ciede2000/
def ciede2000(lab1, lab2):
    L1 = lab1[0]
    A1 = lab1[1]
    B1 = lab1[2]
    L2 = lab2[0]
    A2 = lab2[1]
    B2 = lab2[2]
    C1 = np.sqrt((A1 ** 2.) + (B1 ** 2.))
    C2 = np.sqrt((A2 ** 2.) + (B2 ** 2.))
    aC1C2 = np.average([C1, C2])
    G = 0.5 * (1. - np.sqrt((aC1C2 ** 7.) / ((aC1C2 ** 7.) + (25. ** 7.))))
    a1P = (1. + G) * A1
    a2P = (1. + G) * A2
    c1P = np.sqrt((a1P ** 2.) + (B1 ** 2.))
    c2P = np.sqrt((a2P ** 2.) + (B2 ** 2.))
    if a1P == 0 and B1 == 0: h1P = 0
    else:
        if B1 >= 0: h1P = np.degrees(np.arctan2(B1, a1P))
        else: h1P = np.degrees(np.arctan2(B1, a1P)) + 360.
    if a2P == 0 and B2 == 0: h2P = 0
    else:
        if B2 >= 0: h2P = np.degrees(np.arctan2(B2, a2P))
        else: h2P = np.degrees(np.arctan2(B2, a2P)) + 360.
    dLP = L2 - L1
    dCP = c2P - c1P
    if h2P - h1P > 180: dhC = 1
    elif h2P - h1P < -180: dhC = 2
    else: dhC = 0
    if dhC == 0: dhP = h2P - h1P
    elif dhC == 1: dhP = h2P - h1P - 360.
    else: dhP = h2P + 360 - h1P
    dHP = 2. * np.sqrt(c1P * c2P) * np.sin(np.radians(dhP / 2.))
    aL = np.average([L1, L2])
    aCP = np.average([c1P, c2P])
    if c1P * c2P == 0: haC = 3
    elif np.absolute(h2P - h1P) <= 180: haC = 0
    elif h2P + h1P < 360: haC = 1
    else: haC = 2
    haP = np.average([h1P, h2P])
    if haC == 3: aHP = h1P + h2P
    elif haC == 0: aHP = haP
    elif haC == 1: aHP = haP + 180
    else: aHP = haP - 180
    lPa50 = (aL - 50) ** 2.
    sL = 1. + (0.015 * lPa50 / np.sqrt(20. + lPa50))
    sC = 1. + 0.045 * aCP
    T = 1. - 0.17 * np.cos(np.radians(aHP - 30.)) + 0.24 * np.cos(np.radians(2. * aHP)) + 0.32 * np.cos(np.radians(3. * aHP + 6.)) - 0.2 * np.cos(np.radians(4. * aHP - 63.))
    sH = 1. + 0.015 * aCP * T
    dTheta = 30. * np.exp(-1. * ((aHP - 275.) / 25.) ** 2.)
    rC = 2. * np.sqrt((aCP ** 7.) / ((aCP ** 7.) + (25. ** 7.)))
    rC = 2. * np.sqrt(aCP ** 7. / (aCP ** 7. + 25. ** 7.))
    rT = -np.sin(np.radians(2. * dTheta)) * rC
    fL = dLP / sL / 1.
    fC = dCP / sC / 1.
    fH = dHP / sH / 1.
    dE2000 = np.sqrt(fL ** 2. + fC ** 2. + fH ** 2. + rT * fC * fH)
    return dE2000
