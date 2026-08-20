The dynamical law will be:

        dx/dt = -x
        dy/dt =  x - y

or, under the framework's residual/action convention,

        S(q) = [ x,
                 y - x ]

with:

        qdot = -S(q).

Initial state for later levels:

        q0 = [2, 0]^T

and eventual nominal step:

        h = 1/2.

Forward Euler will therefore later have the simple oracle:

        q0 = [2,   0]
        q1 = [1,   1]
        q2 = [1/2, 1]
        q3 = [1/4, 3/4]
        q4 = [1/8, 1/2].

The structural specimen is only:

    state coordinates:

        Q = {x,y}

    time instants:

        T = {t0,t1,t2,t3,t4}

    time steps:

        E = {e1,e2,e3,e4}

with:

        e1 : t0 -> t1
        e2 : t1 -> t2
        e3 : t2 -> t3
        e4 : t3 -> t4

        E = {t0 -> t1, t1 -> t2, t2 -> t3, t3 -> t4}

Q_T = Q(T)



Create a Level-0 fixture, preferably:

    common/time_sets_fixture.f90

earned at Level 0.

It constructs only:

    Q = index_set('state coordinates', 2)
    T = index_set('time instants', 5)
    E = index_set('time steps', 4)