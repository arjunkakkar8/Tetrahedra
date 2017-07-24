from math import cos, sin, pi, sqrt

margin = 2 * 10**-9

pairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

cyclic_indices = [(1, 2, 3, 4), (2, 3, 4, 1), (3, 4, 1, 2), (4, 1, 2, 3)]

def heron(a, b, c):
    '''
    >>> heron(3, 4, 5)
    6.0
    '''
    s = (a + b + c)/2
    
    return sqrt(s*(s-a)*(s-b)*(s-c))

class Tetrahedron:
    def __init__(self, dihedrals, face_areas, name):
        self._dihedrals = {pairs[i]: dihedrals[i] for i in range(len(dihedrals))}
        self.face_areas = face_areas
        self.name = name

        # Compute edge lengths
        self._compute_edge_lengths()

        # Check edge lengths from a tetrahedron
        self._check_edge_lengths_determine_tetrahedron()
        
        # Check edge lengths determine face areas correctly
        self._check_edge_lengths_face_areas()

        # Check edge lengths determine diherals correctly
        self._check_edge_lengths_dihedrals()

        self._compute_volume()

    def dihedral(self, i, j):
        return self._dihedrals[tuple(sorted([i, j]))]

    @property
    def dihedrals(self):
        return [self._dihedrals[pairs[i]] for i in range(len(pairs))]
        
    def face_area(self, i):
        return self.face_areas[i-1]

    def edge_length(self, i, j):
        return self._edge_lengths[tuple(sorted([i, j]))]

    @property
    def edge_lengths(self):
        return [self._edge_lengths[pairs[i]] for i in range(len(pairs))]
    
    def _compute_edge_lengths(self):
        self._edge_lengths = {}

        theta = self.dihedral
        Delta = self.face_area
        
        for i, j, k, l in cyclic_indices:
            ratio_kl = Delta(j) * sin(theta(k, l))
            ratio_jl = Delta(k) * sin(theta(j, l))
            ratio_jk = Delta(l) * sin(theta(j, k))

            area_from_ratio = heron(ratio_kl, ratio_jl, ratio_jk)

            scaling = sqrt(Delta(i) / area_from_ratio)

            self._edge_lengths[tuple(sorted([k, l]))] = ratio_kl * scaling
            self._edge_lengths[tuple(sorted([j, l]))] = ratio_jl * scaling
            self._edge_lengths[tuple(sorted([j, k]))] = ratio_jk * scaling

    def _check_edge_lengths_determine_tetrahedron(self):
        e = self.edge_length
        
        # Triangle inequality
        for i, j, k, l in cyclic_indices:
            d1, d2, d3 = e(j, k), e(k, l), e(j, l)

            assert d1 < d2 + d3
            assert d2 < d1 + d3
            assert d3 < d1 + d2

        # Determinant condition is verified in another code
        
    def _check_edge_lengths_face_areas(self):
        Delta = self.face_area
        e = self.edge_length
        
        for i, j, k, l in cyclic_indices:
            area = heron(e(j, k), e(k, l), e(j, l))
            
            assert abs(Delta(i) - area) < margin

    def _check_edge_lengths_dihedrals(self): 
        Delta = self.face_area
        theta = self.dihedral
        e = self.edge_length
        
        for i, j in pairs:
            k, l = {1, 2, 3, 4}.difference({i, j})

            D_ij = -e(i, j) ** 4
            D_ij += (e(i, k) ** 2 + e(i, l) ** 2 + e(j, k) ** 2 + e(j, l) ** 2 \
                    - 2 * e(k, l) ** 2) * e(i, j) ** 2
            D_ij += (e(i, k) ** 2 - e(j, k) ** 2) * (e(j, l) ** 2 - e(i, l) ** 2)

            D_ijk = -16 * Delta(l) ** 2
            D_ijl = -16 * Delta(k) ** 2

            a = cos(theta(i, j))
            b = D_ij / sqrt(D_ijk * D_ijl)

            assert abs(a - b) < margin

    def _compute_volume(self):
        Delta = self.face_area
        theta = self.dihedral
        e = self.edge_length
        
        self.volume = 2/3 * Delta(1) * Delta(4) * sin(theta(2, 3)) / e(2, 3)
        
               
dihedrals_data = [
    [3, 4, 5, 10, 6, 6, 'Non-tile A'],
    [3, 5, 5, 10, 10, 4, 'Non-tile B'],
    [3, 5, 10, 10, 6, 4, 'Non-tile C'],
    [3, 6, 6, 8, 8, 4, 'Sommerville No. 3'],
    [3, 6, 10, 10, 10, 3, 'Non-tile D'],
    [4, 4, 4, 5, 6, 10, 'Non-tile E'],
    [4, 4, 4, 6, 6, 8, 'Sommerville No. 2'],
    [4, 4, 8, 8, 4, 6, 'First goldberg family alpha = 2pi/8'],
    [4, 5, 6, 5, 6, 5, 'Non-tile F'],
    [4, 5, 6, 10, 5, 4, 'First goldberg family alpha = 2pi/5'],
    [4, 6, 6, 6, 6, 4, 'Sommerville No. 1'],
]

face_areas_data = [
    [ -0.7236067977, -0.4472135955, -0.2763932023, -0.4472135955],
    [ -0.7236067977, -0.2763932023, -0.4472135955, -0.4472135955],
    [ -0.5437701715, -0.5437701715, -0.5437701715, -0.3360684481],
    [  -0.632455532, -0.4472135955, -0.4472135955, -0.4472135955],
    [  -0.601500955, -0.3717480345,  -0.601500955, -0.3717480345],
    [  0.7071067812,  0.5720614028,  0.3535533906,  0.2185080122],
    [ -0.7071067812,          -0.5, -0.3535533906, -0.3535533906],
    [ -0.5773502692, -0.5773502692, -0.4082482905, -0.4082482905],
    [ -0.5437701715, -0.5437701715, -0.5437701715, -0.3360684481],
    [   0.601500955,  0.3717480345,  0.3717480345,   0.601500955],
    [          -0.5,          -0.5,          -0.5,          -0.5],
]

tetrahedra = []

for i in range(len(dihedrals_data)):
    d = dihedrals_data[i]

    dihedrals = [2*pi/x for x in d[:-1]]
    
    f = face_areas_data[i]

    if f[0] < 0:
        f = [-x for x in f]

    # Check positive face areas
    for x in f:
        assert x > 0

    t = Tetrahedron(dihedrals, f, d[-1])

    tetrahedra.append(t)

for t in tetrahedra:
    print(t.name)
    print(sum(t.face_areas)/t.volume**(2/3))
    print()
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
