import numpy as np
import scipy as sp

def wavenumber_spectrum(ExcitationFreq, freq, r, rho = 1.205, K = 1.41767e5):
    """
    Analytical Solution of Wavenumber Spectrum of 3D Point Source Free Radiation
    
    Input:

    ExcitationFreq ... 3D Point Source Excitation Frequency
    freq ... frequency corresponding to axial wavenumber kz, freq = kz*c0 / 2pi
    r ... Distance to point source on the 2D Plane r = sqrt(x^2 + y^2)

    Output:

    solution in wavenumber domain p(r(x,y),kz(freq))
    """
    from scipy.special import hankel2, kv
    import numpy as np

    c0 = np.sqrt(K / rho)
    k = ExcitationFreq*2*np.pi / c0
    kz = freq*2*np.pi / c0
    if np.isclose(k, np.abs(kz),rtol=1e-10):
        return np.NaN
    elif k > np.abs(kz):
        kxy = np.sqrt(k**2 - kz**2)
        return (-1j / 4) * hankel2(0, kxy*r + 0j)
    elif k < np.abs(kz):
        kxy = np.sqrt(kz**2 - k**2)
        return (1/(2*np.pi))*kv(0, kxy*r + 0j)

def SphericalWaveAnalytic(rhsValue, freq, r, result_type='Pressure', rho = 1.24, c0 = 343, sphere = 'full'):
    '''
    Helper function to generate analytical solution of a monopole sound source.

    Returns the requested acoustic results over a defined distance.
    
    Parameters
    ----------
    rhsValue: complex
        Equivalent to openCFS rhsValues source type, the volume flow rate Q is calculated as Q = rhsValue / (1j*2*pi*freq*rho).
    freq: float
        Excitation frequency of the sound source
    r: array_like
        An array of discrete distance to sound source. r can either be an array of float or an array of cartesian coordinates (x,y,z).
    result_type: str, optional:
        Acoustic result to be returned. Options are 'Pressure', 'Velocity', 'Impedance', 'Intensity', 'SoundPower'.
    rho: float, optional:
        Air density in kg/m^3.
    c0: float, optional:
        Speed of sound in m/s.
    sphere: str, optional:
        Directivity of the sound source. 'full' means no reflections, 'fourth' means two reflecting planes, 'eighth' means three reflecting planes.

    Returns
    -------
    r: ndarray
        Array containing the absolute distances to sound source.
    result:
        Array containing the requested result type.
    '''
    import numpy as np

    factor = 2.0 if sphere=='half' else 4.0 if sphere=='fourth' else 8.0 if sphere=='eighth' else 1.0 if sphere=='full' else 0.0
    assert factor != 0.0, "sphere must must be one of the following type: full, half, fourth, eighth"
    Q = rhsValue * factor / (1j*2*np.pi*freq*rho)
    k = 2*np.pi*freq / c0
    if r.ndim != 1:
        assert r.ndim == 2 and r.shape[-1] == 3, "Coordinate array must have shape (number of coordinate points, 3), current shape is %s" % str(r.shape)
        r = np.linalg.norm(r, axis=-1)
    match result_type:
        case 'Pressure':
            result = 1j*rho*freq*Q*np.exp(-1j*k*r) / 2 / r
        case 'Velocity':
            result =  (Q/4*np.pi)*((1j*k*r + 1)/r**2)*np.exp(-1j*k*r)
        case 'Impedance':
            result =  rho*c0*(1j*k*r/(1j*k*r+1))
        case 'Intensity':
            result =  (rho*c0*k**2*np.abs(Q)**2) / (32*np.pi**2*r**2)
        case 'SoundPower':
            result =  rho*c0*k**2*np.abs(Q)**2 / (8*np.pi)
        case _:
            raise Exception("result_type must be one of 'Pressure, Velocity, Impedance, Intensity, SoundPower'")
    return r, result
    
def FourierSpectrumSphericalWave(freq, kz, xy = [1.0,0.0],rhsValue = 1.0, upperBound = 100, epsilon = 0.1, NumdZ = 1000, dz = 0.01, rho = 1.24, c0=343.0, sphere = 'fourth', printHalfSpectrum = True, rule='simpson'):
    k = 2*np.pi*freq / c0
    factor = 2.0 if sphere=='half' else 4.0 if sphere=='fourth' else 8.0 if sphere=='eighth' else 1.0 if sphere=='full' else 0.0
    assert factor != 0.0, "sphere must must be one of the following type: full, half, fourth, eighth"
    Q = rhsValue * factor / (1j*2*np.pi*freq*rho)
    const = 1j*rho*freq*Q/2
    result = []
    if np.linalg.norm(xy) < 1e-16:
        z_negative = np.linspace(-upperBound,-epsilon,NumdZ)
        z_positive = np.linspace(epsilon,upperBound,NumdZ)
        for Kz in kz:
            # left integral
            integrand_left = (1/(-z_negative))*np.exp(1j*z_negative*(Kz+k))
            integral_left = np.trapz(integrand_left, z_negative)
            # right integral
            integrand_right = (1/(z_positive))*np.exp(1j*z_positive*(Kz-k))
            integral_right = np.trapz(integrand_right, z_positive)
            result.append(integral_left+integral_right)
    else:
        #z = np.linspace(-upperBound,upperBound,NumdZ*2)
        #z = np.linspace(0,upperBound,NumdZ)
        z = np.arange(dz,upperBound,dz)
        #t = np.arange(dz,1,dz)
        #z = t/(1-t)
        r = np.sqrt(xy[0]**2 + xy[1]**2 + z**2)
        for Kz in kz:
            integrand = 1/r * np.exp(-1j*k*r) * np.exp(1j*Kz*z)
            #integrand = 1/r * np.exp(-1j*k*r) * np.exp(1j*Kz*z) * (1/(1-t)**2)
            #result.append(2*np.trapz(integrand,z))
            #result.append(2*sp.integrate.simpson(y=integrand,x=t))
            result.append(2*sp.integrate.simpson(y=integrand,x=z))
    if printHalfSpectrum:
        return kz[np.argwhere(kz >= 0.0).ravel()], const * np.array(result)[np.argwhere(kz >= 0.0).ravel()]
    else:
        return kz, const * np.array(result)

def ExtractResultOnDefinedAxis(hdf_name, region, step, multistep, result = 'acouPressure', dim = '3d', axis = 'x', ShiftCoord = [0,0], ExtractSpectrum = False):
    """
        Helper function for extracting openCFS result along a defined axis (similar to plot over line function in paraview)
        
        The cartesian coordinate system is defined below:
                  ^x
                  |  .(x,y,z)
                  |_____>y
                  /
                 /
                z
        
        Parameters
        ----------
        hdf_name: h5py.File or str
            openCFS hdf5 data file
        region: str
            region name for a subset of coordinates
        step: integer, list, or string
            defining the step as single integer, list of integers or 'last' for laststep or 'all' for all steps
        multistep: int
            defining in which analysis step we should extract the result
        result: str, optional
            openCFS postprocessing result type
        dim: str, optional
            dimension of the problem, either '3d' or '2.5d'
        axis: str, optional
            define along which axis we should extract the results, possible options are 'x','y' or 'z'
        ShiftCoord: arraylike, optional
            Offset to the coordinate axis
        
        Returns
        -------
        AxisValue: arraylike
            Array containing discrete values on the defined axis e.g. for [x1,x2,x3,...]
        Cooridates: arrylike, optional
            Array of coordinates in cartesian system r = (x,y,z), will only be returned if (x,y) != (0,0)
        Result: arraylike
            Array of results
    """
    import hdf5_tools as h5t
    import numpy as np

    coords = h5t.get_coordinates(hdf5_file=hdf_name, region=region)
    result = h5t.get_result(hdf5_file=hdf_name,result=result,region=region,step=step,multistep=multistep)
    match (dim,axis):
        case ('3d','x') | ('2.5d','x') | ('3d','y') | ('2.5d','y'):
            Axis = 1 if axis=='x' else 0
            Axis1 = 1 if axis=='z' else 2
            Axis2 = 0 if axis=='x' else 1 if axis=='y' else 2
            idx1 = np.where(np.abs(coords[:,Axis]) < 1e-16)[0]
            idx2 = np.where(np.abs(coords[:,Axis1]) < 1e-16)[0]
            idx = np.intersect1d(idx1,idx2)
            factor = -1 if axis=='z' else 1
            idx = idx[np.argsort(factor*coords[idx][:,Axis2])]
            return factor*coords[idx][:,Axis2], result[idx]
        case ('3d','z'):
            closestCoord = coords[np.linalg.norm(coords[:,0:2] - ShiftCoord,axis=-1).argmin()]
            idx1 = np.where(np.abs(coords[:,0] - closestCoord[0]) <= 1e-3)[0]
            idx2 = np.where(np.abs(coords[:,1] - closestCoord[1]) <= 1e-3)[0]
            idx = np.intersect1d(idx1,idx2)
            idx = idx[np.argsort(-1*coords[idx][:,2])]
            if np.linalg.norm(ShiftCoord) > 1e-10:
                return coords[idx][:,2], coords[idx], result[idx]
            else:
                return coords[idx][:,2], result[idx]
        case ('2.5d','z'):
            try:
                idx = np.linalg.norm(coords[:,0:2] - ShiftCoord,axis=-1).argmin()
            except:
                print("Cannot find the closest node coordinate to (%f,%f)" % (ShiftCoord[0],ShiftCoord[1]))
            z_coord = h5t.get_step_values(hdf_name)[0] if ExtractSpectrum else h5t.get_step_values(hdf_name)[1]
            result = result[:,idx].ravel()
            if np.linalg.norm(ShiftCoord) > 1e-10:
                coord = np.column_stack((np.full(z_coord.shape, coords[idx,0]), np.full(z_coord.shape, coords[idx,1]), z_coord))
                return z_coord, coord ,result
            else:
                return z_coord, result
                
                
def calc_angle_range(alpha, angle_unit_in_deg = True):
    # calculate the allowed angle range of a complex valued surface impedance for a given normal absorption coefficient alpha
    angle_max = np.arctan( (2*np.sqrt(1-alpha)) / alpha )
    if angle_unit_in_deg:
        return np.rad2deg(angle_max)
    else:
        return angle_max

def calc_surface_impedance(alpha, angle = 0.0, angle_unit_in_deg = True, rho = 1.205, K = 1.41767e5, normalized=True):
    # calculate speed of sound based on air density and compression modulus
    c = np.sqrt(K / rho)
    Z0 = rho * c

    # check if angle is in degree
    if angle_unit_in_deg:
        print_angle_deg = angle
        angle = np.deg2rad(angle)
    else:
        print_angle_deg = np.rad2deg(angle)
    # check if angle exceeds valid range
    angle_range = np.arctan( (2*np.sqrt(1-alpha)) / alpha )
    # if yes, print error message
    assert angle_range - np.abs(angle) >= 1e-8, "Given angle (%.3f°/%.3f Rad) exceeds valid range, which is %.3f°/%.3f Rad" % (print_angle_deg, angle, np.rad2deg(angle_range), angle_range)

    # Coefficients a,b,c corresponding to "Große Lösungsformel"
    a = alpha * (1 + np.tan(angle)**2)
    b = (2 * alpha - 4)
    c = alpha

    # solve a * R**2 + b * R + c = 0 with Große Lösungsformel
    R1 = (-b + np.sqrt(b**2 - 4*a*c)) / (2 * a)
    R2 = (-b - np.sqrt(b**2 - 4*a*c)) / (2 * a)
    X1 = np.tan(angle) * R1
    X2 = np.tan(angle) * R2
    Z1 = R1 + 1j * X1
    Z2 = R2 + 1j * X2
    
    # normalized impedance is acoustic impedance normalzed by characteristic impedance of fluid Z0 (density * speed of sound)
    if normalized:
        return Z1, Z2
    else:
        return Z1 * Z0, Z2 * Z0
