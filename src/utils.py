# ####### The cube on the cubie level is described by the permutation and orientations of corners and edges ############

from enum import IntEnum
from random import randrange

# ######################################## Miscellaneous functions #####################################################

def rotate_right(arr, left, right):
    """"Rotate array arr right between left and right. right is included."""
    temp = arr[right]
    for i in range(right, left, -1):
        arr[i] = arr[i-1]
    arr[left] = temp


def rotate_left(arr, left, right):
    """"Rotate array arr left between left and right. right is included."""
    temp = arr[left]
    for i in range(left, right):
        arr[i] = arr[i+1]
    arr[right] = temp


def c_nk(n, k):
    """Binomial coefficient [n choose k]."""
    if n < k:
        return 0
    if k > n // 2:
        k = n - k
    s, i, j = 1, n, 1
    while i != n - k:
        s *= i
        s //= j
        i -= 1
        j += 1
    return s

# ################## The basic six cube moves described by permutations and changes in orientation #####################

# Setup
class Color(IntEnum):
    """ The possible colors of the cube facelets. Color U refers to the color of the U(p)-face etc.
    Also used to name the faces itself."""
    U = 0
    R = 1
    F = 2
    D = 3
    L = 4
    B = 5

class Co(IntEnum):
    """The names of the corner positions of the cube. Corner URF e.g. has an U(p), a R(ight) and a F(ront) facelet."""
    URF = 0
    UFL = 1
    ULB = 2
    UBR = 3
    DFR = 4
    DLF = 5
    DBL = 6
    DRB = 7


class Ed(IntEnum):
    """The names of the edge positions of the cube. Edge UR e.g. has an U(p) and R(ight) facelet."""
    UR = 0
    UF = 1
    UL = 2
    UB = 3
    DR = 4
    DF = 5
    DL = 6
    DB = 7
    FR = 8
    FL = 9
    BL = 10
    BR = 11

# Up-move
cpU = [Co.UBR, Co.URF, Co.UFL, Co.ULB, Co.DFR, Co.DLF, Co.DBL, Co.DRB]
coU = [0, 0, 0, 0, 0, 0, 0, 0]
epU = [Ed.UB, Ed.UR, Ed.UF, Ed.UL, Ed.DR, Ed.DF, Ed.DL, Ed.DB, Ed.FR, Ed.FL, Ed.BL, Ed.BR]
eoU = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Right-move
cpR = [Co.DFR, Co.UFL, Co.ULB, Co.URF, Co.DRB, Co.DLF, Co.DBL, Co.UBR]  # permutation of the corners
coR = [2, 0, 0, 1, 1, 0, 0, 2]  # changes of the orientations of the corners
epR = [Ed.FR, Ed.UF, Ed.UL, Ed.UB, Ed.BR, Ed.DF, Ed.DL, Ed.DB, Ed.DR, Ed.FL, Ed.BL, Ed.UR]  # permutation of the edges
eoR = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # changes of the permutations of the edges

# Front-move
cpF = [Co.UFL, Co.DLF, Co.ULB, Co.UBR, Co.URF, Co.DFR, Co.DBL, Co.DRB]
coF = [1, 2, 0, 0, 2, 1, 0, 0]
epF = [Ed.UR, Ed.FL, Ed.UL, Ed.UB, Ed.DR, Ed.FR, Ed.DL, Ed.DB, Ed.UF, Ed.DF, Ed.BL, Ed.BR]
eoF = [0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0]

# Down-move
cpD = [Co.URF, Co.UFL, Co.ULB, Co.UBR, Co.DLF, Co.DBL, Co.DRB, Co.DFR]
coD = [0, 0, 0, 0, 0, 0, 0, 0]
epD = [Ed.UR, Ed.UF, Ed.UL, Ed.UB, Ed.DF, Ed.DL, Ed.DB, Ed.DR, Ed.FR, Ed.FL, Ed.BL, Ed.BR]
eoD = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Left-move
cpL = [Co.URF, Co.ULB, Co.DBL, Co.UBR, Co.DFR, Co.UFL, Co.DLF, Co.DRB]
coL = [0, 1, 2, 0, 0, 2, 1, 0]
epL = [Ed.UR, Ed.UF, Ed.BL, Ed.UB, Ed.DR, Ed.DF, Ed.FL, Ed.DB, Ed.FR, Ed.UL, Ed.DL, Ed.BR]
eoL = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Back-move
cpB = [Co.URF, Co.UFL, Co.UBR, Co.DRB, Co.DFR, Co.DLF, Co.ULB, Co.DBL]
coB = [0, 0, 1, 2, 0, 0, 2, 1]
epB = [Ed.UR, Ed.UF, Ed.UL, Ed.BR, Ed.DR, Ed.DF, Ed.DL, Ed.BL, Ed.FR, Ed.FL, Ed.UB, Ed.DB]
eoB = [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1]
########################################################################################################################

CUBE_OK = True


class CubieCube:
    """Represent a cube on the cubie level with 8 corner cubies, 12 edge cubies and the cubie orientations.
    
    Is also used to represent:
    1. the 18 cube moves
    2. the 48 symmetries of the cube.
    """
    def __init__(self, cp=None, co=None, ep=None, eo=None):
        """
        Initializes corners and edges.
        :param cp: corner permutation
        :param co: corner orientation
        :param ep: edge permutation
        :param eo: edge orientation
        """
        if cp is None:
            self.cp = [Co(i) for i in range(8)]  # You may not put this as the default two lines above!
        else:
            self.cp = cp[:]
        if co is None:
            self.co = [0]*8
        else:
            self.co = co[:]
        if ep is None:
            self.ep = [Ed(i) for i in range(12)]
        else:
            self.ep = ep[:]
        if eo is None:
            self.eo = [0] * 12
        else:
            self.eo = eo[:]

    def __str__(self):
        """Print string for a cubie cube."""
        s = ''
        for i in Co:
            s = s + '(' + str(self.cp[i]) + ',' + str(self.co[i]) + ')'
        s += '\n'
        for i in Ed:
            s = s + '(' + str(self.ep[i]) + ',' + str(self.eo[i]) + ')'
        return s

    def __eq__(self, other):
        """Define equality of two cubie cubes."""
        if self.cp == other.cp and self.co == other.co and self.ep == other.ep and self.eo == other.eo:
            return True
        else:
            return False

    def corner_multiply(self, b):
        """Multiply this cubie cube with another cubie cube b, restricted to the corners. Does not change b."""
        c_perm = [0]*8
        c_ori = [0]*8
        ori = 0
        for c in Co:
            c_perm[c] = self.cp[b.cp[c]]
            ori_a = self.co[b.cp[c]]
            ori_b = b.co[c]
            if ori_a < 3 and ori_b < 3:  # two regular cubes
                ori = ori_a + ori_b
                if ori >= 3:
                    ori -= 3
            elif ori_a < 3 <= ori_b:  # cube b is in a mirrored state
                ori = ori_a + ori_b
                if ori >= 6:
                    ori -= 3  # the composition also is in a mirrored state
            elif ori_a >= 3 > ori_b:  # cube a is in a mirrored state
                ori = ori_a - ori_b
                if ori < 3:
                    ori += 3  # the composition is a mirrored cube
            elif ori_a >= 3 and ori_b >= 3:  # if both cubes are in mirrored states
                ori = ori_a - ori_b
                if ori < 0:
                    ori += 3  # the composition is a regular cube
            c_ori[c] = ori
        for c in Co:
            self.cp[c] = c_perm[c]
            self.co[c] = c_ori[c]

    def edge_multiply(self, b):
        """ Multiply this cubie cube with another cubiecube b, restricted to the edges. Does not change b."""
        e_perm = [0]*12
        e_ori = [0]*12
        for e in Ed:
            e_perm[e] = self.ep[b.ep[e]]
            e_ori[e] = (b.eo[e] + self.eo[b.ep[e]]) % 2
        for e in Ed:
            self.ep[e] = e_perm[e]
            self.eo[e] = e_ori[e]

    def multiply(self, b):
        self.corner_multiply(b)
        self.edge_multiply(b)

    def inv_cubie_cube(self, d):
        """Store the inverse of this cubie cube in d."""
        for e in Ed:
            d.ep[self.ep[e]] = e
        for e in Ed:
            d.eo[e] = self.eo[d.ep[e]]

        for c in Co:
            d.cp[self.cp[c]] = c
        for c in Co:
            ori = self.co[d.cp[c]]
            if ori >= 3:
                d.co[c] = ori
            else:
                d.co[c] = -ori
                if d.co[c] < 0:
                    d.co[c] += 3

    def corner_parity(self):
        """Give the parity of the corner permutation."""
        s = 0
        for i in range(Co.DRB, Co.URF, -1):
            for j in range(i - 1, Co.URF - 1, -1):
                if self.cp[j] > self.cp[i]:
                    s += 1
        return s % 2

    def edge_parity(self):
        """Give the parity of the edge permutation. A solvable cube has the same corner and edge parity."""
        s = 0
        for i in range(Ed.BR, Ed.UR, -1):
            for j in range(i - 1, Ed.UR - 1, -1):
                if self.ep[j] > self.ep[i]:
                    s += 1
        return s % 2

# ############################################ other usefull functions #################################################
    def randomize(self):
        """Generate a random cube. The probability is the same for all possible states."""
        def set_edges(idx):
            """The permutation of the 12 edges. 0 <= idx < 12!."""
            self.ep = [i for i in Ed]
            for j in Ed:
                k = idx % (j + 1)
                idx //= j + 1
                while k > 0:
                    rotate_right(self.ep, 0, j)
                    k -= 1
        set_edges(randrange(479001600))  # 12!
        p = self.edge_parity()
        while True:
            self.set_corners(randrange(40320))  # 8!
            if p == self.corner_parity():  # parities of edge and corner permutations must be the same
                break
        self.set_flip(randrange(2048))  # 2^11
        self.set_twist(randrange(2187))  # 3^7

    def verify(self):
        """Check if cubiecube is valid."""
        edge_count = [0]*12
        for i in Ed:
            edge_count[self.ep[i]] += 1
        for i in Ed:
            if edge_count[i] != 1:
                return 'Error: Some edges are undefined.'

        s = 0
        for i in Ed:
            s += self.eo[i]
        if s % 2 != 0:
            return 'Error: Total edge flip is wrong.'

        corner_count = [0] * 8
        for i in Co:
            corner_count[self.cp[i]] += 1
        for i in Co:
            if corner_count[i] != 1:
                return 'Error: Some corners are undefined.'

        s = 0
        for i in Co:
            s += self.co[i]
        if s % 3 != 0:
            return 'Error: Total corner twist is wrong.'

        if self.edge_parity() != self.corner_parity():
            return 'Error: Wrong edge and corner parity'

        return CUBE_OK
########################################################################################################################

# ################################## these cubes represent the basic cube moves ########################################


basicMoveCube = [CubieCube()] * 6
basicMoveCube[Color.U] = CubieCube(cpU, coU, epU, eoU)
basicMoveCube[Color.R] = CubieCube(cpR, coR, epR, eoR)
basicMoveCube[Color.F] = CubieCube(cpF, coF, epF, eoF)
basicMoveCube[Color.D] = CubieCube(cpD, coD, epD, eoD)
basicMoveCube[Color.L] = CubieCube(cpL, coL, epL, eoL)
basicMoveCube[Color.B] = CubieCube(cpB, coB, epB, eoB)
########################################################################################################################

# ################################# these cubes represent the all 18 cube moves ########################################

moveCube = [CubieCube()] * 18
for c1 in Color:
    cc = CubieCube()
    for k1 in range(3):
        cc.multiply(basicMoveCube[c1])
        moveCube[3 * c1 + k1] = CubieCube(cc.cp, cc.co, cc.ep, cc.eo)
########################################################################################################################