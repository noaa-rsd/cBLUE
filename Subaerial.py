"""
cBLUE (comprehensive Bathymetric Lidar Uncertainty Estimator)
Copyright (C) 2019 
Oregon State University (OSU)
Center for Coastal and Ocean Mapping/Joint Hydrographic Center, University of New Hampshire (CCOM/JHC, UNH)
NOAA Remote Sensing Division (NOAA RSD)

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

Contact:
Christopher Parrish, PhD
School of Construction and Civil Engineering
101 Kearney Hall
Oregon State University
Corvallis, OR  97331
(541) 737-5688
christopher.parrish@oregonstate.edu

"""

import logging
from Merge import Merge
from sympy import lambdify, symbols, Matrix, cos, sin
import numpy as np
import numexpr as ne


class SensorModel:
    """This class is used to define and access the sensor model of a particular 
    lidar sensor, including the laser geolocation equation and any
    supporting information or parameters.  Currently, only a single 
    sensor, the Riegl VQ-880-G, is supported, but development plans
    include extending support to the Chiroptera II (or III or IV).

    TODO:  move the a, b uncertainty values here

    """


    eval_type = 'numexpr'

    def __init__(self, sensor):
        self.sensor = sensor
        self.R, self.fR = self.set_rotation_matrix_airplane() 
        self.M = self.set_rotation_matrix_scanning_sensor()
        self.obs_eq, self.obs_eq_pre_poly = self.define_obseration_equation()

    def set_rotation_matrix_airplane(self):
        """define rotation matrix for airplane

        This method first generates the airplane rotation matrix, R, using
        symbolic calculations.  The symbolic components of the matrix R are
        then "functionized", or "lambdified", for faster processing (because
        symbolic calculations are relatively slow).  The components of R are
        functionized separately from the general observation equation, which
        includes R, M (the sensor rotation matrix), and polynomial-surface-
        correction terms, because R is later used to estimate parameters
        describing the assumed scan pattern, which is an approximation of the
        manufacturer's proprietary scan pattern.

        .. image:: ../images/eq_AirplaneRotationMatrix.png

        :return: Matrix
        :return: List[lambdify functions]
        """

        r, p, h = symbols('r p h')
        R1 = Matrix([[1, 0, 0],
                     [0, cos(r), -sin(r)],
                     [0, sin(r), cos(r)]])

        R2 = Matrix([[cos(p), 0, sin(p)],
                     [0, 1, 0],
                     [-sin(p), 0, cos(p)]])

        R3 = Matrix([[cos(h), -sin(h), 0],
                     [sin(h), cos(h), 0],
                     [0, 0, 1]])

        R = R3 * R2 * R1

        # "functionize" the necessary R components for a and b estimation
        r00 = lambdify((h, p), R[0], self.eval_type)
        r01 = lambdify((r, p, h), R[1], self.eval_type)
        r10 = lambdify((h, p), R[3], self.eval_type)
        r11 = lambdify((r, p, h), R[4], self.eval_type)
        r20 = lambdify(p, R[6], self.eval_type)
        r21 = lambdify((r, p), R[7], self.eval_type)

        fR = [r00, r01, None,
              r10, r11, None,
              r20, r21, None]
        
        return R, fR

    @staticmethod
    def set_rotation_matrix_scanning_sensor():
        """define the lidar sensor rotation matrix

        This method generates the rotation matrix associated with the
        scanning sensor.  The variables a, b, and w describe the assumed
        scan pattern, which is an approximation of the manufacturer's
        proprietary scan pattern.

        a: the rotation in the YZ plane
        b: the rotation in the XZ plane

        .. image:: ../images/eq_SensorRotationMatrix.png

        :return Matrix M: the scanning sensor rotation matrix
        """
        a, b = symbols('a b')
        M1 = Matrix([[1, 0, 0],
                     [0, cos(a), -sin(a)],
                     [0, sin(a), cos(a)]])

        M2 = Matrix([[cos(b), 0, sin(b)],
                     [0, 1, 0],
                     [-sin(b), 0, cos(b)]])

        M = M2 * M1
        
        return M

    def define_obseration_equation(self):
        """define the lidar geolocation observation equation

        The inital observation equation is defined as follows:

        .. image:: ../images/eq_OriginalObsEq.png

        However, to account for the differences between the assumed sensor 
        model and the proprietary sensor model, the initial observation equation 
        is modified to include terms derived from polynomial surface fitting of 
        differences in the X, Y, and Z components of the LAS positions and the
        positions calculated from the intial cBLUE observation equation.

        .. image:: ../images/eq_ModifiedObsEq.png

        :return: (sympy object, sympy object, sympy object, function)
        """

        # create variables for symbolic computations
        a, b, r, p, h, x, y, z, rho, p00, p10, p01, p20, p11, p02, p21, p12, p03 \
            = symbols('a b r p h x y z rho p00 p10 p01 p20 p11 p02 p21 p12 p03')

        # define observation equations
        # [00, 01, 02      matrix       [0 1 2
        #  10, 11, 12   ---indices-->    3 4 5
        #  20, 21, 22]                   6 7 8]
        F1 = x - rho * (self.R[0] * self.M[2] + self.R[1] * self.M[5] + self.R[2] * self.M[8])
        F2 = y - rho * (self.R[3] * self.M[2] + self.R[4] * self.M[5] + self.R[5] * self.M[8])
        F3 = z - rho * (self.R[6] * self.M[2] + self.R[7] * self.M[5] + self.R[8] * self.M[8])

        # converting symbolic to function (for faster computations)
        fF1 = lambdify((a, b, h, p, r, rho, x), F1, self.eval_type)
        fF2 = lambdify((a, b, h, p, r, rho, y), F2, self.eval_type)
        fF3 = lambdify((a, b, p, r, rho, z), F3, self.eval_type)
        fF_orig = [fF1, fF2, fF3]

        # least squares adjustment mimics the Matlab "fit (poly23)" function
        polysurfcorr = p00 + \
                       p10 * a + \
                       p01 * b + \
                       p20 * a ** 2 + \
                       p11 * a * b + \
                       p02 * b ** 2 + \
                       p21 * a ** 2 * b + \
                       p12 * a * b ** 2 + \
                       p03 * b ** 3

        # Redefine obs eqs by adding polynomial surface fitting eq
        F1 += polysurfcorr
        F2 += polysurfcorr
        F3 += polysurfcorr
        
        return (F1, F2, F3,),  fF_orig

    def estimate_rho_a_b(self, data):
        """calculates estimates for rho, alpha, and beta

        This method calculates the estimated values for rho, alpha, and beta, which
        are the lidar range, angle in the YZ plane, and angle in the XZ plane,
        respectively (see the following image).

        .. image:: ../images/rho_alpha_beta.png

        Alpha and beta are used to model the scan pattern, as a substitute for the
        actual, unknown, proprietary scan pattern model implemented by the
        manufacturer.  Polynomial-surface error modeling is used to account for the
        positional differences resulting from the difference between the cBLUE scan
        model and the manufacturer scan model.

        :return: (list[], list[], list[], list[])
        """

        #0       t_sbet
        #1       t_las 
        #2       x_las 
        #3       y_las 
        #4       z_las 
        #5       x_sbet
        #6       y_sbet
        #7       z_sbet
        #8       r     
        #9       p     
        #10      h   

        # TODO: only pass needed columns

        x_las = data[2]  # x_las
        y_las = data[3]  # y_las
        z_las = data[4]  # z_las
        x_sbet = data[5]  # x_sbet
        y_sbet = data[6]  # y_sbet
        z_sbet = data[7]  # z_sbet

        rho_x = ne.evaluate("x_las - x_sbet")
        rho_y = ne.evaluate("y_las - y_sbet")
        rho_z = ne.evaluate("z_las - z_sbet")

        fR0 = self.fR[0](data[10], data[9])
        fR3 = self.fR[3](data[10], data[9])
        fR6 = self.fR[6](data[9])
        fR1 = self.fR[1](data[8], data[9], data[10])
        fR4 = self.fR[4](data[8], data[9], data[10])
        fR7 = self.fR[7](data[8], data[9])

        rho_est = ne.evaluate("sqrt(rho_x**2 + rho_y**2 + rho_z**2)")
        b_est = ne.evaluate("arcsin(((fR0 * rho_x) + (fR3 * rho_y) + (fR6 * rho_z)) / (-rho_est))")
        a_est = ne.evaluate("arcsin(((fR1 * rho_x) + (fR4 * rho_y) + (fR7 * rho_z)) / (rho_est * cos(b_est)))")

        return rho_est, a_est, b_est

    def calc_poly_surf_coeffs(self, a_est, b_est, dx, dy, dz, itv=10):
        """estimates error model using polynomial surface fitting

        This method calculates the coefficients of the polynomial-surface
        error model intended to account for the positional errors resulting
        from differences between the sensor model implemented in cBLUE and
        the unknown, proprietary manufacturer sensor model.

        The original Matlab research code used a 'fit' function with a
        'poly23' option, which is emulated here by using np.linalg.lstsq
        with terms for a, b, a^2, ab, b^2, a^2b, ab^2, and b^3.

        Only every itv-th point is used to calculate the polynomial surface
        coefficients, for small speed gains in the calculations of the 
        coefficients.

        :param a_est:
        :param b_est:
        :param dx:
        :param dy:
        :param dz:
        :param itv:
        :return: list[tuple, tuple, tuple] TODO: verify
        """

        B0 = b_est[::itv]
        A0 = a_est[::itv]

        A = np.vstack((
            ne.evaluate('A0 * 0 + 1'),
            ne.evaluate('A0'),
            ne.evaluate('B0'),
            ne.evaluate('A0 ** 2'),
            ne.evaluate('A0 * B0'),
            ne.evaluate('B0 ** 2'),
            ne.evaluate('A0 ** 2 * B0'),
            ne.evaluate('A0 * B0 ** 2'),
            ne.evaluate('B0 ** 3'))).T

        (coeffs_x, __, __, __) = np.linalg.lstsq(A, dx[::itv], rcond=None)
        (coeffs_y, __, __, __) = np.linalg.lstsq(A, dy[::itv], rcond=None)
        (coeffs_z, __, __, __) = np.linalg.lstsq(A, dz[::itv], rcond=None)

        return coeffs_x, coeffs_y, coeffs_z

    @staticmethod
    def calcRMSE(data):
        """calc root mean square error for input data  NOT CURRENTLY USED"""
        num_coords, num_points = data.shape
        AMDE = np.mean(np.abs(data), axis=1)  # average mean distance error
        RMSE = sqrt(sum(sum(np.square(data))) / num_points)  # root mean squares error
        logging.info('Mean Difference:\n'
                     'X: {:.3f}\n'
                     'Y: {:.3f}\n'
                     'Z: {:.3f}'.format(AMDE[0], AMDE[1], AMDE[2]))
        logging.info('RMSE: {:.3f}\n'.format(RMSE))

    def calc_diff(self, subaer_pos_pre, x_las, y_las, z_las):
        """calculate the difference between the las position and the initial cBLUE position

        This method calculates the difference between the x, y, and z components of the
        positions in the las file and the respective cBLUE-calculated position components.
        Ideally, a cBLUE-calculated position would identically match the corresponding
        position in the las file, but due to differences between the proprietary manufacturer
        sensor model and the sensor model used by cBLUE, the positions are not identical.
        The differences calculated by this method are used in the polynomial-surface error
        modeling process to correct for the errors caused by the sensor model discrepancies.

        :param subaer_pos_pre:
        :return:
        """

        # calc diff between true and est las xyz (aer_pos = Laser Estimates)
        aer_x_pre = subaer_pos_pre[0]
        aer_y_pre = subaer_pos_pre[1]
        aer_z_pre = subaer_pos_pre[2]

        dx = ne.evaluate("x_las - aer_x_pre")
        dy = ne.evaluate("y_las - aer_y_pre")
        dz = ne.evaluate("z_las - aer_z_pre")

        return dx, dy, dz

    def calc_aer_pos_pre(self, rho_est, a_est, b_est, data):
        """calculates the inital cBLUE aubaerial position

        This method calculates the inital cBLUE subaerial position using the
        'lambdified' geolocation equation (without the polynomial-surface
        error terms).

        The data parameter contains the following ndarrays:

        =====   =========   =======================
        Index   ndarray     description
        =====   =========   =======================
        0       t_sbet      sbet timestamps
        1       t_las       las timestamps
        2       x_las       las x coordinates
        3       y_las       las y coordinates
        4       z_las       las z coordinates
        5       x_sbet      sbet x coordinates
        6       y_sbet      sbet y coordinates
        7       z_sbet      sbet z coordinates
        8       r           sbet roll
        9       p           sbet pitch
        10      h           sbet heading
        =====   =========   =======================

        :param rho_est:
        :param a_est:
        :param b_est:
        :return:
        """


        aer_x_pre_poly = self.obs_eq_pre_poly[0](a_est, b_est, data[10], data[9], data[8], rho_est, data[5])
        aer_y_pre_poly = self.obs_eq_pre_poly[1](a_est, b_est, data[10], data[9], data[8], rho_est, data[6])
        aer_z_pre_poly = self.obs_eq_pre_poly[2](a_est, b_est, data[9], data[8], rho_est, data[7])
        
        return aer_x_pre_poly, aer_y_pre_poly, aer_z_pre_poly

    def calc_aer_pos(self, coeffs, a_est, b_est, aer_pos_pre):
        """calculates the final cBLUE subearial position

        This method calculates the final cBLUE subaerial position by adding
        a polynomial-surface modelled error term

        :param coeffs:
        :param a_est:
        :param b_est:
        :param aer_pos_pre:
        :return:
        """

        A = np.vstack((
            ne.evaluate('a_est * 0 + 1'),
            ne.evaluate('a_est'),
            ne.evaluate('b_est'),
            ne.evaluate('a_est ** 2'),
            ne.evaluate('a_est * b_est'),
            ne.evaluate('b_est ** 2'),
            ne.evaluate('a_est ** 2 * b_est'),
            ne.evaluate('a_est * b_est ** 2'),
            ne.evaluate('b_est ** 3'))).T

        aer_x_pre = aer_pos_pre[0]
        aer_y_pre = aer_pos_pre[1]
        aer_z_pre = aer_pos_pre[2]

        err_x = np.sum(A * coeffs[0], axis=1)
        err_y = np.sum(A * coeffs[1], axis=1)
        err_z = np.sum(A * coeffs[2], axis=1)

        aer_pos_x = ne.evaluate('aer_x_pre + err_x')
        aer_pos_y = ne.evaluate('aer_y_pre + err_y')
        aer_pos_z = ne.evaluate('aer_z_pre + err_z')

        return aer_pos_x, aer_pos_y, aer_pos_z

    def calc_aer_pos_err(self, aer_pos, data):
        """calculates the difference between the las and cBLUE positions

        This method calculates the differences between the x, y, and z
        components of the final cBLUE positions and the corresponding
        las file positions.

        The data parameter contains the following ndarrays:

        =====   =========   =======================
        Index   ndarray     description
        =====   =========   =======================
        0       t_sbet      sbet timestamps
        1       t_las       las timestamps
        2       x_las       las x coordinates
        3       y_las       las y coordinates
        4       z_las       las z coordinates
        5       x_sbet      sbet x coordinates
        6       y_sbet      sbet y coordinates
        7       z_sbet      sbet z coordinates
        8       r           sbet roll
        9       p           sbet pitch
        10      h           sbet heading
        =====   =========   =======================

        :param aer_pos:
        :return:
        """

        aer_x = aer_pos[0]
        aer_y = aer_pos[1]
        aer_z = aer_pos[2]

        x_las = data[2]
        y_las = data[3]
        z_las = data[4]

        aer_x_err = ne.evaluate('aer_x - x_las')
        aer_y_err = ne.evaluate('aer_y - y_las')
        aer_z_err = ne.evaluate('aer_z - z_las')

        return [aer_x_err, aer_y_err, aer_z_err]

class Jacobian:
    """This class is used to calculate and evaluate the Jacobian of a
    sensor model's laser geolocation equation.  The class Jacobian attempts
    to decouple a Jacobian and the data used to evaluate it.  For example, 
    the inputs to the lambdified Jacobian components are not hard-coded in 
    the function call, but are determined from accessing the 
    .__code__.co_varnames attribute of the Jacobian component and then 
    looking up the corresponding values in a dict.  Although this somewhat
    decouples the Jacobian from the data used to evaluate it, the dict 
    containing the corresponding values is manually created, separate from
    the sensor model.  Development plans for future versions include
    decoupling the Jacobian and the data to evaluate it even more, by 
    creating the dict based on the sensor model.

    Two key modules that are used throughout are sympy and numexpr:

    The module sympy is used to symbolically define the laser geolocation
    equation and the corresponding Jacobian and to numerically evalulate
    the Jacobian.  
    
    The module numexpr is used to accelerate calculations using large 
    numpy arrays (https://github.com/pydata/numexpr).  One characteristic 
    of numexpr is that numexpr expressions do not allow indexing of 
    variables, so what might normally be coded as, for example,
    *var = data[1] * 3* would require something like *data1 = data[1]* 
    before executing the numexpr expression *"var = data1 * 3"*.
   
    """

    def __init__(self, S):
        self.S = S
        self.OEx = S.obs_eq[0]
        self.OEy = S.obs_eq[1]
        self.OEz = S.obs_eq[2]
        self.Jx, self.Jy, self.Jz  = self.form_jacobian()
        self.lJx, self.lJy, self.lJz = self.lambdify_jacobian()

    def form_jacobian(self):
        """generate the jacobian of the specified geolocation equation

        This method generates the Jacobian (i.e., the matrix of partial
        derivatives with respect to component variables) of the specified
        geoloation equation using the sympy symbolic math package.  Using
        sympy to symbolically calculate the Jacobian simplifies the coding
        of what would otherwise be very long equations.

        .. image:: ../images/eq_Jacobian.png

        :return (Matrix, Matrix, Matrix): sympy matrices for x, y, and z Jacobian components
        """

        a, b, r, p, h, x, y, z, rho = symbols('a b r p h x y z rho')

        v = Matrix([a, b, r, p, h, x, y, z, rho])  # vector of unknowns

        Jx = Matrix([self.OEx]).jacobian(v)
        Jy = Matrix([self.OEy]).jacobian(v)
        Jz = Matrix([self.OEz]).jacobian(v)

        return Jx, Jy, Jz

    def lambdify_jacobian(self, eval_type='numexpr'):
        """turn the symbolic Jacobian into a function for faster computation

        This method "lambdifies" (or "functionizes") the Jacobian components, for
        faster calculations. Part of this lambdify process includes simplifying
        the numerous trigonometric calculations of the Jacobian by defining the
        Jacobian functions to be functions of the sines and cosines of the various
        parameters, instead of the parameters directly.

        Reference:
        https://docs.sympy.org/latest/modules/utilities/lambdify.html

        :param str eval_type: the eval type for sympy lambdification
        :return (function, function, function): lambdified x, y, and z Jacobian components
        """

        # create variables for symbolic computations
        a, b, r, p, h, rho, \
        p00, p10, p01, p20, p11, p02, p21, p12, p03, \
        sin_a, sin_b, sin_r, sin_p, sin_h, \
        cos_a, cos_b, cos_r, cos_p, cos_h \
            = symbols('a b r p h rho '
                      'p00 p10 p01 p20 p11 p02 p21 p12 p03 '
                      'sin_a sin_b sin_r sin_p sin_h '
                      'cos_a cos_b cos_r cos_p cos_h')

        trig_substitutions = [
            (sin(a), sin_a),
            (sin(b), sin_b),
            (sin(r), sin_r),
            (sin(p), sin_p),
            (sin(h), sin_h),
            (cos(a), cos_a),
            (cos(b), cos_b),
            (cos(r), cos_r),
            (cos(p), cos_p),
            (cos(h), cos_h)]

        # functionize the trig terms of the Jacobian components
        Jxsub = [j.subs(trig_substitutions) for j in self.Jx]
        Jysub = [j.subs(trig_substitutions) for j in self.Jy]
        Jzsub = [j.subs(trig_substitutions) for j in self.Jz]

        # functionize the Jacobian x, y, and z components
        # (9 terms in each Jacobian component correspond to a, b, r, p, h, x, y, z, and rho 
        lJx = [lambdify(Jxsub[i].free_symbols, Jxsub[i], eval_type) for i in range(9)]  
        lJy = [lambdify(Jysub[i].free_symbols, Jysub[i], eval_type) for i in range(9)]
        lJz = [lambdify(Jzsub[i].free_symbols, Jzsub[i], eval_type) for i in range(9)]

        return lJx, lJy, lJz

    def calc_trig_terms(self, a_est, b_est, r, p, h):
        """helper method to evaluate the trigonometric terms in the Jacobian

        This method aims to simplify evaluation of the Jacobian by pre-evaluating
        the trigonometic terms of the Jacobian.  The reasoning is that this speeds
        up the computations because the trigonometric terms are only evaluated
        once, instead of every time they show up in the Jacobian.

        :param a_est: a calculated from the data
        :param b_est: b calculated from the data
        :param r: roll data
        :param p: pitch data
        :param h: heave data
        :return tupe(ndarray): the evaluated trigonometric terms
        """

        sin_a = ne.evaluate("sin(a_est)")
        sin_b = ne.evaluate("sin(b_est)")

        cos_a = ne.evaluate("cos(a_est)")
        cos_b = ne.evaluate("cos(b_est)")

        sin_r = ne.evaluate("sin(r)")
        sin_p = ne.evaluate("sin(p)")
        sin_h = ne.evaluate("sin(h)")

        cos_r = ne.evaluate("cos(r)")
        cos_p = ne.evaluate("cos(p)")
        cos_h = ne.evaluate("cos(h)")

        return (sin_a, sin_b, sin_r, sin_p, sin_h,
                cos_a, cos_b, cos_r, cos_p, cos_h, )

    def get_calc_vals_for_J_eval(self, data):
        """calculatse and assembles the values needed to evaluate the Jacobian

        This methods calculates and assembles the values needed to evaluate the Jacobian.  

        1. estimate rho, a, and b from data
        2. use rho, a, and b estimates to calculate initial X, Y, and Z
        3. calculate difference betwween inital X, Y, and Z and LAS X, Y, and Z (dX, dY, and dZ)
        4. calculate polynomial surfae coefficients to account for dX, dY, and dZ
        5. precalculate sine and cosine of attitude data to simplify evaluation of Jacobian

        The returned dictionary contains the following data:

        =========   ========================
        data        description
        =========   ========================
        a           calculated a values
        b           calculated b values
        rho         calculated rho values
        p_coeffs    {'x':{coeff:value,...},'y':{coeff:value,...},'z':{coeff:value,...}}
        sin_a       calculated sin(a) values
        sin_b       calculated sin(b) values
        sin_r       calculated sin(r) values
        sin_p       calculated sin(p) values
        sin_h       calculated sin(h) values
        cos_a       calculated cos(a) values
        cos_b       calculated cos(b) values
        cos_r       calculated cos(r) values
        cos_p       calculated cos(p) values
        cos_h       calculated cos(h) values
        =========   ========================

        :param data
        :return dict: calcualted values used to evaluate Jacobian
        
        """

        # estimate rho, a, and b from data
        rho_est, a_est, b_est = self.S.estimate_rho_a_b(data)
        
        # use rho, a, and b estimates to calculate initial estimate of X, Y, Z
        aer_pos_pre = self.S.calc_aer_pos_pre(rho_est, a_est, b_est, data)
        
        # calculate differece between initial X, Y, and Z estimates and las X, Y, and Z
        dx, dy, dz = self.S.calc_diff(aer_pos_pre, data[2], data[3], data[4])
        
        # calculate polynomial surface coefficients to account for differences dx, dy, and dz
        p_coeffs_x, p_coeffs_y, p_coeffs_z = self.S.calc_poly_surf_coeffs(a_est, b_est, dx, dy, dz)
        
        # precalculate sin and cos of attitude data to simplify evaluation of Jacobian
        trig_subs = self.calc_trig_terms(a_est, b_est, data[8], data[9], data[10]) 

        # 
        p_coeffs_vars = ['p00', 'p10', 'p01', 'p20', 'p11', 'p02', 'p21', 'p12', 'p03']

        J_params = {
            'a': a_est,
            'b': b_est,
            'rho': rho_est,
            'p_coeffs': {
                'x': {k: p_coeffs_x[i] for i, k in enumerate(p_coeffs_vars)},  # e.g., {'p00': p00x, ...}
                'y': {k: p_coeffs_y[i] for i, k in enumerate(p_coeffs_vars)},
                'z': {k: p_coeffs_z[i] for i, k in enumerate(p_coeffs_vars)},
                },
            'sin_a': trig_subs[0],
            'sin_b': trig_subs[1],
            'sin_r': trig_subs[2],
            'sin_p': trig_subs[3],
            'sin_h': trig_subs[4],
            'cos_a': trig_subs[5],
            'cos_b': trig_subs[6],
            'cos_r': trig_subs[7],
            'cos_p': trig_subs[8],
            'cos_h': trig_subs[9],
            }        

        return J_params

    def get_J_term_values(self, J_comp, J_term, values_for_J_eval): 
        """gets the calculated values needed to evaluate the specified Jacobian component

        This method retrieves from the passed 'values_for_J_eval parameter the
        calculated values needed to evaluate the specified Jacobian component (i.e., 
        the x, y, or z component).

        :param J_comp:
        :param J_term:
        :param values_for_J_eval:
        :return vals:
        
        """

        vals = []
        term_vars = J_term.__code__.co_varnames
        for var in term_vars:
            if var[0] == 'p':
                vals.append(values_for_J_eval['p_coeffs'][J_comp][var])
            else:
                vals.append(values_for_J_eval[var])
        
        return vals

    def eval_jacobian(self, data):
        """evaluate the Jacobian of the modified laser geolocation equation

        This method evaluates the Jacobian by passing the relevant parameters
        to the lambdified functions representing the x, y, and z components
        of the Jacobian.

        To simplify the Jacobian evaluation, only the non-zero terms are kept.  
        Accordingly, the rows of variance/covariance matrix corresponding to the
        Jacobian zero terms are deleted.  Additionally, the Jacobian evaluation
        is simplied further by not calling get_J_term_values() for Jacobian
        terms equal to 1; rather, the corresponding row in the evaluated Jacobian
        array is set to all 1s.

        :param data:
        :return (ndarray, ndarray, ndarray): x, y, and z evaluated Jacobian components
        """

        J_param_values = self.get_calc_vals_for_J_eval(data)

        Jx = np.vstack(
            (self.lJx[0](*self.get_J_term_values('x', self.lJx[0], J_param_values)),
             self.lJx[1](*self.get_J_term_values('x', self.lJx[1], J_param_values)),
             self.lJx[2](*self.get_J_term_values('x', self.lJx[2], J_param_values)),
             self.lJx[3](*self.get_J_term_values('x', self.lJx[3], J_param_values)),
             self.lJx[4](*self.get_J_term_values('x', self.lJx[4], J_param_values)),
             np.ones(data[0].size),
             self.lJx[8](*self.get_J_term_values('x', self.lJx[8], J_param_values))))

        Jy = np.vstack(
            (self.lJy[0](*self.get_J_term_values('y', self.lJy[0], J_param_values)),
             self.lJy[1](*self.get_J_term_values('y', self.lJy[1], J_param_values)),
             self.lJy[2](*self.get_J_term_values('y', self.lJy[2], J_param_values)),
             self.lJy[3](*self.get_J_term_values('y', self.lJy[3], J_param_values)),
             self.lJy[4](*self.get_J_term_values('y', self.lJy[4], J_param_values)),
             np.ones(data[0].size),
             self.lJy[8](*self.get_J_term_values('y', self.lJy[8], J_param_values))))

        Jz = np.vstack(
            (self.lJz[0](*self.get_J_term_values('z', self.lJz[0], J_param_values)),
             self.lJz[1](*self.get_J_term_values('z', self.lJz[1], J_param_values)),
             self.lJz[2](*self.get_J_term_values('z', self.lJz[2], J_param_values)),
             self.lJz[3](*self.get_J_term_values('z', self.lJz[3], J_param_values)),
             np.ones(data[0].size),
             self.lJz[8](*self.get_J_term_values('z', self.lJz[8], J_param_values))))

        return (Jx, Jy, Jz, )


class Subaerial:
    """
    This class provides the functionality to calculate the subaerial
    portion of the total propagated uncertainty (TPU), given the Jacobian 
    of a laser geolocation equation, merged lidar/trajectory 
    data, and the standard deviations of the provided data.  
    
    The following table lists the contents of merged lidar/trajectory
    data array:

    =====   =========   =======================
    Index   ndarray     description
    =====   =========   =======================
    0       t_sbet      sbet timestamps
    1       t_las       las timestamps
    2       x_las       las x coordinates
    3       y_las       las y coordinates
    4       z_las       las z coordinates
    5       x_sbet      sbet x coordinates
    6       y_sbet      sbet y coordinates
    7       z_sbet      sbet z coordinates
    8       r           sbet roll
    9       p           sbet pitch
    10      h           sbet heading
    =====   =========   =======================
        
    The following table lists the contents of the array of standard
    deviations corresponding to the variables of the merged data
    array:

    =====   =========   =======================
    Index   ndarray     description
    =====   =========   =======================
    0       std_ang1    ang1 uncertainty
    1       std_ang2    ang2 uncertainty
    2       std_r       sbet roll uncertainty
    3       std_p       sbet pitch uncertainty
    4       std_h       sbet heading uncertainty
    5       stdx_sbet   sbet x uncertainty
    6       stdy_sbet   sbet y uncertainty
    7       stdz_sbet   sbet z uncertainty
    8       std_rho     ?
    =====   =========   =======================

    :param Jacobian J: Jacobian object
    :param ndarray: merged Lidar/Trajectory data
    :param ndarray: standard deviations of component variables
    """

    def __init__(self, J, merged_data, stddev):

        self.J = J  # Jacobian object
        self.merged_data = merged_data  # merged-data ndarray
        self.stddev = stddev  # nparray of standard deviations

    def propogate_uncertainty(self, J_eval):
        """propogates the subaerial uncertatinty

        This method propogates the uncertainty of the component uncertainties
        using the following equation:

        .. image:: ../images/eq_PropogateError.png

        Because only the non-zero terms of the Jacobian are kept (to
        simplify Jacobian evaluation), the rows in the variance/covariance
        matrix corresponding to the zero terms of the Jacobian are
        deleted.  The table below summarizes which terms of the Jacobian
        are zero (and one).

        =====   ========    ==  ==  ==
        index   variable    Jx  Jy  Jz
        =====   ========    ==  ==  ==
        0       a           .   .   .
        1       b           .   .   .
        2       r           .   .   .
        3       p           .   .   .
        4       h           .   .   0
        5       x           1   0   0       
        6       y           0   1   0 
        7       z           0   0   1
        8       rho         .   .   .   
        =====   ========    ==  ==  ==

        :param tuple(ndarray) J_eval:  evaluated Jacobian values for X, Y, and Z components
        :return (ndarray, ndarray, list[str]): subaerial THU, subaerial TVU, THU and TVU column headers
        """

        stddev = self.stddev
        V = ne.evaluate("stddev * stddev")  # variance = stddev**2

        # delete the rows corresponding to the Jacobian terms that equal 0
        Vx = np.delete(V, [6, 7], 0)
        Vy = np.delete(V, [5, 7], 0)
        Vz = np.delete(V, [4, 5, 6], 0)

        Jx = J_eval[0]
        Jy = J_eval[1]
        Jz = J_eval[2]

        sum_Jx = ne.evaluate("sum(Jx * Jx * Vx, axis=0)")
        sum_Jy = ne.evaluate("sum(Jy * Jy * Vy, axis=0)")
        sum_Jz = ne.evaluate("sum(Jz * Jz * Vz, axis=0)")

        sx = ne.evaluate("sqrt(sum_Jx)")
        sy = ne.evaluate("sqrt(sum_Jy)")
        aer_tvu = ne.evaluate("sqrt(sum_Jz)")
        aer_thu = ne.evaluate('sqrt(sx**2 + sy**2)')

        return aer_thu, aer_tvu, ['subaerial_thu', 'subaerial_tvu']

    def calc_subaerial_tpu(self):
        """calculates the subaerial uncertainty

        This method calculates the subaerial uncertainty through two major
        steps:

        1. EVALUATE JACOBIAN

            The eval_jacobian() method of the Jacobian object evaluates the
            Jacobian with the merged lidar/trajectory data passed to it.  Although
            the Jacobian is first calculated symbolically, it's evaluated as a
            collection of lambda functions, which are "used to calculate numerical
            values very fast." (https://docs.sympy.org/latest/modules/utilities/lambdify.html)

        2. PROPAGATE UNCERTAINTY

            Once the Jacobian is evaluated, uncertainty is progagated by multiplying the
            square of the Jacobian with the squares of the standard deviations defined in the
            stddev parameter.  The covariances are assumed to be zero.  TODO:  explain how this is 
            implemented differently than shown by the propagation equation because the covariances are 
            assumed to be 0.

        :return: (ndarray, ndarray, list[str])
        """

        # EVALUATE JACOBIAN
        J_eval = self.J.eval_jacobian(self.merged_data)

        # PROPAGATE UNCERTAINTY
        aer_thu, aer_tvu, aer_cols = self.propogate_uncertainty(J_eval)

        return aer_thu, aer_tvu, aer_cols


if __name__ == '__main__':
    pass
