'''
This tidal ellipse package is taken from the PyFVCOM module written Created by
Pierre Cazenave available at https://github.com/pwcazenave/PyFVCOM.
The author retains the copyright of this program.
'''

"""
Python translation of Zhigang Xu's tidal_ellipse MATLAB tools, available at
http://woodshole.er.usgs.gov/operations/sea-mat/tidal_ellipse-html/ap2ep.html

Converted to Python by Pierre Cazenave (Plymouth Marine Laboratory), October
2012. Email: pica@pml.ac.uk

Authorship Copyright:

   The author retains the copyright of this program, while you are welcome
to use and distribute it as long as you credit the author properly and
respect the program name itself. Particularly, you are expected to retain the
original author's name in this original version or any of its modified version
that you might make. You are also expected not to essentially change the name
of the programs except for adding possible extension for your own version you
might create, e.g. ap2ep_xx is acceptable.  Any suggestions are welcome and
enjoy my program(s)!


Author Info:
_______________________________________________________________________________

  Zhigang Xu, Ph.D.
  (pronounced as Tsi Gahng Hsu)
  Research Scientist
  Coastal Circulation
  Bedford Institute of Oceanography
  1 Challenge Dr.
  P.O. Box 1006                    Phone  (902) 426-2307 (o)
  Dartmouth, Nova Scotia           Fax    (902) 426-7827
  CANADA B2Y 4A2                   email xuz@dfo-mpo.gc.ca
_______________________________________________________________________________

Release Date: Nov. 2000, Revised on May. 2002 to adopt Foreman's northern semi
major axis convention.

"""

import numpy as np
import matplotlib.pyplot as plt


def ap2ep(Au, PHIu, Av, PHIv, plot_demo=False):
    """
    Convert tidal amplitude and phase lag (ap-) parameters into tidal ellipse
    (ep-) parameters. Please refer to ep2ap for its inverse function.

    Usage:

    SEMA, ECC, INC, PHA, w = ap2ep(Au, PHIu, Av, PHIv, plot_demo=False)

    Where:

        Au, PHIu, Av, PHIv are the amplitudes and phase lags (in degrees) of
        u- and v- tidal current components. They can be vectors or
        matrices or multidimensional arrays.

        plot_demo is an optional argument, when it is supplied as an array
        of indices, say [i j k l], the program will plot an ellipse
        corresponding to Au[i, j, k, l], PHIu[i, j, k, l], Av[i, j, k, l], and
        PHIv[i, j, k, l]. Defaults to False (i.e. no plot).

        Any number of dimensions are allowed as long as your computer
        resource can handle.

        SEMA: Semi-major axes, or the maximum speed.

        ECC:  Eccentricity, the ratio of semi-minor axis over the semi-major
              axis; its negative value indicates that the ellipse is traversed
              in clockwise direction.

        INC:  Inclination, the angles (in degrees) between the semi-major axes
              and u-axis.

        PHA:  Phase angles, the time (in angles and in degrees) when the tidal
              currents reach their maximum speeds,  (i.e.  PHA=omega*tmax).

              These four ep-parameters will have the same dimensionality (i.e.,
              vectors, or matrices) as the input ap-parameters.

        w:    A matrix whose rows allow for plotting ellipses and whose columns
              are for different ellipses corresponding columnwise to SEMA. For
              example, plot(np.real(w[0, :]), np.imag(w[0, :])) will let you
              see the first ellipse. You may need to use squeeze function when
              w is a more than two dimensional array. See example.py.

    Document:   tidal_ellipse.ps

    Revisions: May  2002, by Zhigang Xu,  --- adopting Foreman's northern semi
    major axis convention.

    For a given ellipse, its semi-major axis is undetermined by 180. If we
    borrow Foreman's terminology to call a semi major axis whose direction lies
    in a range of [0, 180) as the northern semi-major axis and otherwise as a
    southern semi major axis, one has freedom to pick up either northern or
    southern one as the semi major axis without affecting anything else.
    Foreman (1977) resolves the ambiguity by always taking the northern one as
    the semi-major axis. This revision is made to adopt Foreman's convention.
    Note the definition of the phase, PHA, is still defined as the angle
    between the initial current vector, but when converted into the maximum
    current time, it may not give the time when the maximum current first
    happens; it may give the second time that the current reaches the maximum
    (obviously, the 1st and 2nd maximum current times are half tidal period
    apart) depending on where the initial current vector happen to be and its
    rotating sense.

    Version 2, May 2002

    Converted to Python by Pierre Cazenave, October 2012.

    Authorship Copyright:

       The author retains the copyright of this program, while  you are welcome
    to use and distribute it as long as you credit the author properly and
    respect the program name itself. Particularly, you are expected to retain
    the original author's name in this original version or any of its modified
    version that you might make. You are also expected not to essentially
    change the name of the programs except for adding possible extension for
    your own version you might create, e.g. ap2ep_xx is acceptable.  Any
    suggestions are welcome and enjoy my program(s)!


    Author Info:
    _______________________________________________________________________
      Zhigang Xu, Ph.D.
      (pronounced as Tsi Gahng Hsu)
      Research Scientist
      Coastal Circulation
      Bedford Institute of Oceanography
      1 Challenge Dr.
      P.O. Box 1006                    Phone  (902) 426-2307 (o)
      Dartmouth, Nova Scotia           Fax    (902) 426-7827
      CANADA B2Y 4A2                   email xuz@dfo-mpo.gc.ca
    _______________________________________________________________________

    Release Date: Nov. 2000, Revised on May. 2002 to adopt Foreman's northern
    semi major axis convention.

    """

    # Assume the input phase lags are in degrees and convert them in radians.
    PHIu = PHIu / 180 * np.pi
    PHIv = PHIv / 180 * np.pi

    # Make complex amplitudes for u and v
    i = 1j
    u = Au * np.exp(-i * PHIu)
    v = Av * np.exp(-i * PHIv)

    # Calculate complex radius of anticlockwise and clockwise circles:
    wp = (u + i * v) / 2           # for anticlockwise circles
    wm = np.conj(u - i * v) / 2    # for clockwise circles
    # and their amplitudes and angles
    Wp = np.abs(wp)
    Wm = np.abs(wm)
    THETAp = np.angle(wp)
    THETAm = np.angle(wm)

    # calculate ep-parameters (ellipse parameters)
    SEMA = Wp + Wm                 # Semi Major Axis, or maximum speed
    SEMI = Wp - Wm                 # Semi Minor Axis, or minimum speed
    ECC = SEMI / SEMA              # Eccentricity

    PHA = (THETAm - THETAp) / 2    # Phase angle, the time (in angle) when
                                   # the velocity reaches the maximum
    INC = (THETAm + THETAp) / 2    # Inclination, the angle between the
                                   # semi major axis and x-axis (or u-axis).

    # convert to degrees for output
    PHA = PHA / np.pi*180
    INC = INC / np.pi*180
    THETAp = THETAp / np.pi*180
    THETAm = THETAm / np.pi*180

    # map the resultant angles to the range of [0, 360].
    PHA = np.mod(PHA + 360, 360)
    INC = np.mod(INC + 360, 360)

    # Mar. 2, 2002 Revision by Zhigang Xu    (REVISION_1)
    # Change the southern major axes to northern major axes to conform the tidal
    # analysis convention  (cf. Foreman, 1977, p. 13, Manual For Tidal Currents
    # Analysis Prediction, available in www.ios.bc.ca/ios/osap/people/foreman.htm)
    k = np.fix(INC / 180)
    INC = INC - k * 180
    PHA = PHA + k * 180
    PHA = np.mod(PHA, 360)

    if plot_demo:
        plot_ell(SEMA, ECC, INC, PHA, plot_demo)

    ndot = np.prod(np.shape(SEMA))
    dot = 2 * np.pi / ndot
    ot = np.arange(0, 2 * np.pi, dot)
    w = wp.flatten() * np.exp(i * ot) + wm.flatten() * np.exp(-i * ot)
    w = np.reshape(w, np.shape(wp))

    return SEMA, ECC, INC, PHA, w


def ep2ap(SEMA, ECC, INC, PHA, plot_demo=False):
    """
    Convert tidal ellipse parameters into amplitude and phase lag parameters.
    Its inverse is app2ep.m. Please refer to app2ep for the meaning of the
    inputs and outputs.

    Zhigang Xu
    Oct. 20, 2000

    Converted to Python by Pierre Cazenave, October 2012.

    Authorship Copyright:

        The author of this program retains the copyright of this program, while
    you are welcome to use and distribute this program as long as you credit
    the author properly and respect the program name itself. Particularly, you
    are expected to retain the original author's name in this original version
    of the program or any of its modified version that you might make.  You are
    also expected not to essentially change the name of the programs except for
    adding possible extension for your own version you might create, e.g.
    app2ep_xx is acceptable.  Any suggestions are welcome and enjoy my
    program(s)!


    Author Info:
    _______________________________________________________________________
      Zhigang Xu, Ph.D.
      (pronounced as Tsi Gahng Hsu)
      Research Scientist
      Coastal Circulation
      Bedford Institute of Oceanography
      1 Challenge Dr.
      P.O. Box 1006                    Phone  (902) 426-2307 (o)
      Dartmouth, Nova Scotia           Fax    (902) 426-7827
      CANADA B2Y 4A2                   email xuz@dfo-mpo.gc.ca
    _______________________________________________________________________

    Release Date: Nov. 2000

    """

    i = 1j

    Wp = (1 + ECC) / 2 * SEMA
    Wm = (1 - ECC) / 2 * SEMA
    THETAp = INC - PHA
    THETAm = INC + PHA

    # Convert degrees into radians
    THETAp = THETAp / 180 * np.pi
    THETAm = THETAm / 180 * np.pi

    # Calculate wp and wm.
    wp = Wp * np.exp(i * THETAp)
    wm = Wm * np.exp(i * THETAm)

    ndot = np.prod(np.shape(SEMA))
    dot = 2 * np.pi / ndot
    ot = np.arange(0, 2 * np.pi, dot)
    w = wp.flatten() * np.exp(i * ot) + wm.flatten() * np.exp(-i * ot)
    w = np.reshape(w, np.shape(wp))

    # Calculate cAu, cAv --- complex amplitude of u and v
    cAu = wp + np.conj(wm)
    cAv = -i * (wp-np.conj(wm))
    Au = np.abs(cAu)
    Av = np.abs(cAv)
    PHIu = -np.angle(cAu) * 180 / np.pi
    PHIv = -np.angle(cAv) * 180 / np.pi

    # flip angles in the range of [-180 0) to the range of [180 360).
    id = PHIu < 0
    PHIu[id] = PHIu[id] + 360
    id = PHIv < 0
    PHIv[id] = PHIv[id] + 360

    if plot_demo:
        plot_ell(SEMA, ECC, INC, PHA, plot_demo)

    return Au, PHIu, Av, PHIv, w


def cBEpm(g, f, sigma, nu, kappa, z, h):
    """
    Evaluate the theoretical vertical profiles (or Bottom Ekman spiral
    profiles) of tidal currents in the two rotary directions driven by
    half-unit of sea surface gradients in the two directions respectively. Eddy
    viscosity is assumed as vertically invariant. See tidal_ellipse.ps for more
    details.

    Inputs:

        g:      acceleration gravity
        f:      the Coriolis parameter
        nu:     the eddy viscosity
        kappa:  the bottom frictional coefficient
        z:      the vertical coordinates, can be a vector but must be
                within [0 -h];
        h:      the water depth, must be positive.

        Note: except for z, all other inputs must be scalars.

    Outputs:

        BEp and BEm, the same dimensions of z,  the outputs for the vertical
            velocity profiles driven respectively by a unit of sea surface
            slope in the positive rotation direction and negative rotation
            direction for when the eddy viscosity is vertically invariant. See
            the associated document for more details.

    Authorship Copyright:

       The author of this program retains the copyright of this program, while
    you are welcome to use and distribute this program as long as you credit
    the author properly and respect the program name itself. Particularly,
    you are expected to retain the original author's name in this original
    version of the program or any of its modified version that you might make.
    You are also expected not to essentially change the name of the programs
    except for adding possible extension for your own version you might create,
    e.g. ap2ep_xx is acceptable.  Any suggestions are welcome and enjoying my
    program(s)!


    Author Info:
    _______________________________________________________________________
      Zhigang Xu, Ph.D.
      (pronounced as Tsi Gahng Hsu)
      Research Scientist
      Coastal Circulation
      Bedford Institute of Oceanography
      1 Challenge Dr.
      P.O. Box 1006                    Phone  (902) 426-2307 (o)
      Dartmouth, Nova Scotia           Fax    (902) 426-7827
      CANADA B2Y 4A2                   email zhigangx@emerald.bio.dfo.ca
                                             zhigang_xu_98@yahoo.com
    _______________________________________________________________________

    Release Date: Nov. 2000

    """

    if (len(g) > 1) | (len(f) > 1) | (len(sigma) > 1) | \
            (len(nu) > 1) | (len(kappa) > 1) | (len(h) > 1):
        print('inputs of g, f, sigma, nu, kappa, and h should be all scalars!')
        raise

    if (any(z / h > 0)) | (any(z / h < -1)):
        print('z must be negative and must be within [0 -h]')

    delta_e = np.sqrt(2 * nu / f)  # Ekman depth
    alpha = (1 + 1j) / delta_e * np.sqrt(1 + sigma / f)
    beta = (1 + 1j) / delta_e * np.sqrt(1 - sigma / f)

    BEp = get_BE(g, alpha, h, z, nu, kappa)
    BEm = get_BE(g, beta, h, z, nu, kappa)

    return BEp, BEm


def get_BE(g, alpha, h, z, nu, kappa):
    """ Child function of cBEpm """

    z = z.flatten()
    z_h = z / h
    ah = alpha * h
    az = alpha * z
    ah2 = ah * 2
    anu_k = alpha * nu / kappa
    nu_kh = nu / (kappa * h)

    # Series solution
    if abs(ah) < 1:
        T = 10
        C = -g * h * h / (nu * (1 + anu_k * np.tanh(ah))) * 2
        A1 = (1 - z_h * z_h) / 2 + nu_kh
        B1 = np.exp(-ah) / (1 + np.exp(-ah2))
        B = B1
        series_sum = A1 * B1

        for t in np.arange(2, T):
            t2 = 2*t
            A = (1 - z_h**t2) / t2 + nu_kh
            B = B * ah * ah / (t2 - 1) / (t2 - 2)
            series_sum = series_sum + A * B

        BE = C*series_sum

    # Finite solution
    else:
        c = -g * h * h / nu
        denom = (np.exp(az - ah) + np.exp(-(az + ah))) / (1 + np.exp(-2 * ah))

        numer = 1 + anu_k * np.tanh(ah)
        BE = c * ((1 - denom / numer) / (ah * ah))

    return BE


def sub2ind(shape, pos):
    """
    Substitute of MATLAB's sub2ind function for NumPy.

    t = numpy.random.random([2, 4, 5, 2])
    n = sub2ind(numpy.shape(t), [1, 2, 4, 1])
    >>> n
    69

    From http://stackoverflow.com/questions/4114461

    """
    res = 0
    acc = 1
    for pi, si in zip(reversed(pos), reversed(shape)):
        res += pi * acc
        acc *= si

    return res


def plot_ell(SEMA, ECC, INC, PHA, IND=[1]):
    """
    An auxiliary function used in ap2ep and ep2ap for plotting tidal ellipse.
    The inputs, MA, ECC, INC and PHA are the output of ap2ep and IND is a
    vector for indices for plotting a particular ellipse, e.g., if IND=[2 3 1]
    the ellipse corresponding to the indices of [2,3,1] will be plotted.

    By default, the first ellipse is always plotted.

    Converted to Python by Pierre Cazenave, October 2012.

    """

    len_IND = len(IND)
    if IND:
        cmd = 'sub2ind(size_SEMA, '
        if len_IND == 1:
            titletxt = 'Ellipse '
        else:
            titletxt = 'Ellipse ('

        for k in range(len_IND):
            if k == 0:
                cmd = cmd + '[' + str(IND[k])
            else:
                cmd = cmd + ',' + str(IND[k])

            if k < len_IND-1:
                titletxt = titletxt + str(IND[k]) + ','
            elif len_IND == 1:
                titletxt = titletxt + str(IND[k])
            else:
                titletxt = titletxt + str(IND[k]) + ')'

        cmd = 'n = ' + cmd + '])'
        # This is pretty nasty, but it works.
        exec(cmd)

        plt.gcf()
        plt.clf()
        do_the_plot(SEMA.flatten()[n], ECC.flatten()[n], INC.flatten()[n], PHA.flatten()[n])
        titletxt = titletxt + ',  (red) green (anti-) clockwise component'
        plt.title(titletxt)
    elif len_IND:
        print('IND input contains zero element(s)!\nNo ellipse will be plotted.')


def do_the_plot(SEMA, ECC, INC, PHA):
    """
    Ellipse plot subfunction.

    Converted to Python by Pierre Cazenave, October 2012.

    """

    i = 1j

    SEMI = SEMA * ECC
    Wp = (1 + ECC) / 2 * SEMA
    Wm = (1 - ECC) / 2 * SEMA
    THETAp = INC - PHA
    THETAm = INC + PHA

    # Convert degrees into radians
    THETAp = THETAp / 180 * np.pi
    THETAm = THETAm / 180 * np.pi
    INC = INC / 180 * np.pi
    PHA = PHA / 180 * np.pi

    # Calculate wp and wm.
    wp = Wp * np.exp(i * THETAp)
    wm = Wm * np.exp(i * THETAm)

    dot = np.pi / 36
    ot = np.arange(0, 2 * np.pi, dot)
    a = wp * np.exp(i * ot)
    b = wm * np.exp(-i * ot)
    w = a + b

    wmax = SEMA * np.exp(i * INC)
    wmin = SEMI * np.exp(i * (INC + np.pi / 2))

    plt.plot(np.real(w), np.imag(w))
    plt.axis('equal')
    plt.hold('on')
    plt.plot([0, np.real(wmax)], [0, np.imag(wmax)], 'm')
    plt.plot([0, np.real(wmin)], [0, np.imag(wmin)], 'm')
    plt.xlabel('u')
    plt.ylabel('v')
    plt.plot(np.real(a), np.imag(a), 'r')
    plt.plot(np.real(b), np.imag(b), 'g')
    plt.plot([0, np.real(a[0])], [0, np.imag(a[0])], 'ro')
    plt.plot([0, np.real(b[0])], [0, np.imag(b[0])], 'go')
    plt.plot([0, np.real(w[0])], [0, np.imag(w[0])], 'bo')
    plt.plot(np.real(a[0]), np.imag(a[0]), 'ro')
    plt.plot(np.real(b[0]), np.imag(b[0]), 'go')
    plt.plot(np.real(w[0]), np.imag(w[0]), 'bo')
    plt.plot(np.real([a[0], a[0]+b[0]]), np.imag([a[0], a[0]+b[0]]), linestyle='--', color='g')
    plt.plot(np.real([b[0], a[0]+b[0]]), np.imag([b[0], a[0]+b[0]]), linestyle='--', color='r')

    for n in range(len(ot)):
        plt.hold('on')
        plt.plot(np.real(a[n]), np.imag(a[n]), 'ro')
        plt.plot(np.real(b[n]), np.imag(b[n]), 'go')
        plt.plot(np.real(w[n]), np.imag(w[n]), 'bo')

    plt.hold('off')
    plt.show()


def prep_plot(SEMA, ECC, INC, PHA):
    """
    Take the output of ap2ep (SEMA, ECC, INC, and PHA) and prepare it for
    plotting.

    This is extracted from do_the_plot above, but allows quicker access when
    all that is required is a plot of an ellipse, for which only w is really
    required.

    Returns w, wmin and wmax (w is used for plotting the ellipse, see
    plot_ell).

    """

    i = 1j

    SEMI = SEMA * ECC
    Wp = (1 + ECC) / 2 * SEMA
    Wm = (1 - ECC) / 2 * SEMA
    THETAp = INC - PHA
    THETAm = INC + PHA

    # Convert degrees into radians
    THETAp = THETAp / 180 * np.pi
    THETAm = THETAm / 180 * np.pi
    INC = INC / 180 * np.pi
    PHA = PHA / 180 * np.pi

    # Calculate wp and wm.
    wp = Wp * np.exp(i * THETAp)
    wm = Wm * np.exp(i * THETAm)

    dot = np.pi / 36
    ot = np.arange(0, 2 * np.pi, dot)
    a = wp * np.exp(i * ot)
    b = wm * np.exp(-i * ot)
    w = a + b

    # Repeat the first position in w so we close the ellipse.
    w = np.hstack((w, w[0]))

    wmax = SEMA * np.exp(i * INC)
    wmin = SEMI * np.exp(i * (INC + np.pi / 2))

    return w, wmin, wmax


if __name__ == '__main__':

    """
    Replicate the tidal ellipse example file from Zhigang Xu's tidal_ellipse
    MATLAB toolbox.

    Pierre Cazenave (Plymouth Marine Laboratory), October 2012.

    """

    # Demonstrate how to use ap2ep and ep2ap
    Au = np.random.random([4, 3, 2])           # so 4x3x2 multi-dimensional matrices
    Av = np.random.random([4, 3, 2])           # are used for the demonstration.
    Phi_v = np.random.random([4, 3, 2]) * 360  # phase lags inputs are expected to
    Phi_u = np.random.random([4, 3, 2]) * 360  # be in degrees.

    plt.figure(1)
    plt.clf()
    SEMA, ECC, INC, PHA, w = ap2ep(Au, Phi_u, Av, Phi_v, [2, 3, 1])
    plt.figure(2)
    plt.clf()
    rAu, rPhi_u, rAv, rPhi_v, rw = ep2ap(SEMA, ECC, INC, PHA, [2, 3, 1])

    # Check if ep2ap has recovered Au, Phi_u, Av, Phi_v
    print(np.max(np.abs(rAu - Au).flatten()))        # = 9.9920e-16, = 2.22044604925e-16
    print(np.max(np.abs(rAv - Av).flatten()))        # = 6.6613e-16, = 7.77156117238e-16
    print(np.max(np.abs(rPhi_u - Phi_u).flatten()))  # = 4.4764e-13, = 1.70530256582e-13
    print(np.max(np.abs(rPhi_v - Phi_v).flatten()))  # = 1.1369e-13, = 2.27373675443e-13
    print(np.max(np.max(np.abs(w - rw).flatten())))  # = 1.3710e-15, = 1.1322097734e-15
    # For the random realization I (Zhigang Xu) had, the differences are listed
    # on the right hand of the above column. I (Pierre Cazenave) got the second
    # column with the Python version. What are yours?

    # Zhigang Xu
    # Nov. 12, 2000
    # Pierre Cazenave
    # October, 2012
