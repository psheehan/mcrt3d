from .mcrt3d import Gas
import numpy

def load_gas_properties_lamda(filename):
    f = open(filename)

    for i in range(3):
        f.readline()

    mass = float(f.readline())

    f.readline()
    nlev = int(f.readline())
    f.readline()

    levels = numpy.empty(nlev, dtype=int)
    J = numpy.empty(nlev, dtype="<U6")
    E = numpy.empty(nlev, dtype=float)
    g = numpy.empty(nlev, dtype=float)
    for i in range(nlev):
        levels[i], E[i], g[i], J[i] = tuple(f.readline().split())

    f.readline()
    ntrans = int(f.readline())
    f.readline()

    transition = numpy.empty(ntrans, dtype=int)
    J_u = numpy.empty(ntrans, dtype=int)
    J_l = numpy.empty(ntrans, dtype=int)
    A_ul = numpy.empty(ntrans, dtype=float)
    nu = numpy.empty(ntrans, dtype=float)
    E_u = numpy.empty(ntrans, dtype=float)
    for i in range(ntrans):
        transition[i], J_u[i], J_l[i], A_ul[i], nu[i], \
                E_u[i] = tuple(f.readline().split())

    nu *= 1.0e9

    """
    f.readline()
    npartners = int(f.readline())

    partners = []
    temp = []
    J_u_coll = []
    J_l_coll = []
    gamma = []
    for i in range(npartners):
        f.readline()
        partners.append(f.readline())
        f.readline()
        ncolltrans = int(f.readline())
        f.readline()
        ncolltemps = int(f.readline())
        f.readline()
        temp.append(numpy.array(f.readline().split(), dtype=float))
        f.readline()

        J_u_coll.append(numpy.empty(ncolltrans, dtype=int))
        J_l_coll.append(numpy.empty(ncolltrans, dtype=int))
        gamma.append(numpy.empty((ncolltrans,ncolltemps), \
                dtype=float))

        for j in range(ncolltrans):
            temp, J_u_coll[i][j], J_l_coll[i][j], temp2 = \
                    tuple(f.readline().split(None,3))
            gamma[i][j,:] = numpy.array(temp2.split())
    """

    f.close()

    g = Gas(mass, levels, E, g, J, transition, J_u, J_l, A_ul, nu, E_u)

    return g
