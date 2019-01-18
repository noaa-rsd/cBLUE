import logging
logging.basicConfig(format='%(asctime)s:%(message)s', level=logging.INFO)
from Merge import Merge
from sympy import *
import numpy as np
import numexpr as ne

"""
This class provides the functionality to calculate the subaerial
portion of the total propagated uncertainty (TPU).
"""


class Subaerial:

    def __init__(self, D):
        """
        :param List[ndarray] D: the merged data

        The following table lists the contents of D:

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
        11      std_ang1    ang1 uncertainty
        12      std_ang2    ang2 uncertainty
        13      std_r       sbet roll uncertainty
        14      std_p       sbet pitch uncertainty
        15      std_h       sbet heading uncertainty
        16      stdx_sbet   sbet x uncertainty
        17      stdy_sbet   sbet y uncertainty
        18      stdz_sbet   sbet z uncertainty
        19      std_rho     ?
        =====   =========   =======================
        """
        self.is_empty = not D
        self.num_pts = None
        self.x_las = D[2]
        self.y_las = D[3]
        self.z_las = D[4]
        self.x_sbet = D[5]
        self.y_sbet = D[6]
        self.z_sbet = D[7]
        self.r0 = D[8]
        self.p0 = D[9]
        self.h0 = D[10]
        self.x0 = self.x_sbet
        self.y0 = self.y_sbet
        self.z0 = self.z_sbet
        self.C = np.asarray(D[11:])
        self.R, self.fR = self.set_rotation_matrix_airplane()
        self.M = self.set_rotation_matrix_scanning_sensor()
        self.F1, self.F2, self.F3, self.fF_orig = self.define_obseration_equation()
        self.fJ1, self.fJ2, self.fJ3 = self.calc_jacobian(self.F1, self.F2, self.F3)

    @staticmethod
    def set_rotation_matrix_airplane():
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

        :return:
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
        R = R3*R2*R1
        logging.info(R)

        # "functionize" the necessary R components for a, b, and w estimation
        # (http://docs.sympy.org/latest/modules/utilities/lambdify.html)
        r00 = lambdify((h, p), R[0], 'numexpr')
        r01 = lambdify((r, p, h), R[1], 'numexpr')
        r10 = lambdify((h, p), R[3], 'numexpr')
        r11 = lambdify((r, p, h), R[4], 'numexpr')
        r20 = lambdify(p, R[6], 'numexpr')
        r21 = lambdify((r, p), R[7], 'numexpr')
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
        w: ?

        :return Matrix M: the scanning sensor rotation matrix
        """
        a, b, w = symbols('a b w')
        M1 = Matrix([[1, 0, 0],
                     [0, cos(a), -sin(a)],
                     [0, sin(a), cos(a)]])
        M2 = Matrix([[cos(b), 0, sin(b)],
                     [0, 1, 0],
                     [-sin(b), 0, cos(b)]])
        M3 = Matrix([[cos(w), -sin(w), 0],
                     [sin(w), cos(w), 0],
                     [0, 0, 1]])
        M = M3*M2*M1
        logging.info('sensor rotation matrix: {}'.format(M))
        return M

    def define_obseration_equation(self):
        """define the lidar geolocation observation equation



        :return:
        """

        # create variables for symbolic computations
        a, b, w, r, p, h, x, y, z, rho, p00, p10, p01, p20, p11, p02, p21, p12, p03 \
            = symbols('a b w r p h x y z rho p00 p10 p01 p20 p11 p02 p21 p12 p03')

        # define observation equations
        # [00, 01, 02      matrix       [0 1 2
        #  10, 11, 12   ---indices-->    3 4 5
        #  20, 21, 22]                   6 7 8]
        F1 = x - rho * (self.R[0] * self.M[2] + self.R[1] * self.M[5] + self.R[2] * self.M[8])
        F2 = y - rho * (self.R[3] * self.M[2] + self.R[4] * self.M[5] + self.R[5] * self.M[8])
        F3 = z - rho * (self.R[6] * self.M[2] + self.R[7] * self.M[5] + self.R[8] * self.M[8])

        # converting symbolic to function (faster)
        fF1 = lambdify((a, b, h, p, r, rho, w, x), F1, 'numexpr')
        fF2 = lambdify((a, b, h, p, r, rho, w, y), F2, 'numexpr')
        fF3 = lambdify((a, b, p, r, rho, w, z), F3, 'numexpr')
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
        return F1, F2, F3, fF_orig

    def calc_jacobian(self, F1, F2, F3):
        """define observation equations and calculate jacobian
        (i.e., matrix of partial derivatives with respect to
        component variables) using sympy symbolic math package
        """
        a, b, r, p, h, x, y, z, rho = symbols('a b r p h x y z rho')

        v = Matrix([a, b, r, p, h, x, y, z, rho])  # vector of unknowns
        J1 = Matrix([F1]).jacobian(v)
        J2 = Matrix([F2]).jacobian(v)
        J3 = Matrix([F3]).jacobian(v)
        return self.lambdify_jacobian(J1, J2, J3)

    def lambdify_jacobian(self, J1, J2, J3):
        """lambdify the Jacobian components (faster than symbolic calculations) and
        simplify the numerous trigonometric calculations of the Jacobian by
        defining the Jacobian functions to be functions of the sines and cosines
        of the various parameters, instead of the parameters directly
        """
        # create variables for symbolic computations
        a, b, w, r, p, h, rho, \
        p00, p10, p01, p20, p11, p02, p21, p12, p03, \
        sin_a, sin_b, sin_w, sin_r, sin_p, sin_h, \
        cos_a, cos_b, cos_w, cos_r, cos_p, cos_h \
            = symbols('a b w r p h rho '
                      'p00 p10 p01 p20 p11 p02 p21 p12 p03 '
                      'sin_a sin_b sin_w sin_r sin_p sin_h '
                      'cos_a cos_b cos_w cos_r cos_p cos_h')

        J1sub = []
        for i, j in enumerate(J1):
            J1sub.append(j.subs([
                (sin(a), sin_a),
                (sin(b), sin_b),
                (sin(w), sin_w),
                (sin(r), sin_r),
                (sin(p), sin_p),
                (sin(h), sin_h),
                (cos(a), cos_a),
                (cos(b), cos_b),
                (cos(w), cos_w),
                (cos(r), cos_r),
                (cos(p), cos_p),
                (cos(h), cos_h)]))

        J2sub = []
        for i, j in enumerate(J2):
            J2sub.append(j.subs([
                (sin(a), sin_a),
                (sin(b), sin_b),
                (sin(w), sin_w),
                (sin(r), sin_r),
                (sin(p), sin_p),
                (sin(h), sin_h),
                (cos(a), cos_a),
                (cos(b), cos_b),
                (cos(w), cos_w),
                (cos(r), cos_r),
                (cos(p), cos_p),
                (cos(h), cos_h)]))

        J3sub = []
        for i, j in enumerate(J3):
            J3sub.append(j.subs([
                (sin(a), sin_a),
                (sin(b), sin_b),
                (sin(w), sin_w),
                (sin(r), sin_r),
                (sin(p), sin_p),
                (sin(h), sin_h),
                (cos(a), cos_a),
                (cos(b), cos_b),
                (cos(w), cos_w),
                (cos(r), cos_r),
                (cos(p), cos_p),
                (cos(h), cos_h)]))

        eval_type = 'numexpr'


        '''functionize the x component'''
        J1sub0_vars = (a, b, rho, p10, p21, p12, p11, p20,
                       sin_a, sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_p, cos_h)

        J1sub1_vars = (a, b, rho, p02, p03, p21, p12, p01, p11,
                       sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_p, cos_h)

        J1sub2_vars = (rho, sin_a, sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_h)

        J1sub3_vars = (rho, sin_a, sin_b, sin_w, sin_r, sin_p,
                       cos_a, cos_b, cos_w, cos_r, cos_p, cos_h)

        J1sub4_vars = (rho, sin_a, sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_p, cos_h)

        J1sub5_vars = ()
        J1sub6_vars = ()
        J1sub7_vars = ()

        J1sub8_vars = (sin_a, sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_p, cos_h)

        fJ1sub0 = lambdify(J1sub0_vars, J1sub[0], eval_type)
        fJ1sub1 = lambdify(J1sub1_vars, J1sub[1], eval_type)
        fJ1sub2 = lambdify(J1sub2_vars, J1sub[2], eval_type)
        fJ1sub3 = lambdify(J1sub3_vars, J1sub[3], eval_type)
        fJ1sub4 = lambdify(J1sub4_vars, J1sub[4], eval_type)
        fJ1sub5 = lambdify(J1sub5_vars, J1sub[5], eval_type)
        fJ1sub6 = lambdify(J1sub6_vars, J1sub[6], eval_type)
        fJ1sub7 = lambdify(J1sub7_vars, J1sub[7], eval_type)
        fJ1sub8 = lambdify(J1sub8_vars, J1sub[8], eval_type)

        '''functionize the y component'''
        J2sub0_vars = (a, b, rho, p10, p21, p12, p11, p20,
                       sin_a, sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_p, cos_h)

        J2sub1_vars = (a, b, rho, p02, p03, p21, p12, p01, p11,
                       sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_p, cos_h)

        J2sub2_vars = (rho, sin_a, sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_h)

        J2sub3_vars = (rho, sin_a, sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_p)

        J2sub4_vars = (rho, sin_a, sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_p, cos_h)

        J2sub5_vars = ()
        J2sub6_vars = ()
        J2sub7_vars = ()

        J2sub8_vars = (sin_a, sin_b, sin_w, sin_r, sin_p, sin_h,
                       cos_a, cos_b, cos_w, cos_r, cos_p, cos_h)

        fJ2sub0 = lambdify(J2sub0_vars, J2sub[0], eval_type)
        fJ2sub1 = lambdify(J2sub1_vars, J2sub[1], eval_type)
        fJ2sub2 = lambdify(J2sub2_vars, J2sub[2], eval_type)
        fJ2sub3 = lambdify(J2sub3_vars, J2sub[3], eval_type)
        fJ2sub4 = lambdify(J2sub4_vars, J2sub[4], eval_type)
        fJ2sub5 = lambdify(J2sub5_vars, J2sub[5], eval_type)
        fJ2sub6 = lambdify(J2sub6_vars, J2sub[6], eval_type)
        fJ2sub7 = lambdify(J2sub7_vars, J2sub[7], eval_type)
        fJ2sub8 = lambdify(J2sub8_vars, J2sub[8], eval_type)

        '''functionize the z component'''
        J3sub0_vars = (a, b, rho, p10, p21, p12, p11, p20,
                       sin_a, sin_b, sin_w, sin_r, sin_p,
                       cos_a, cos_b, cos_w, cos_r, cos_p)

        J3sub1_vars = (a, b, rho, p02, p03, p21, p12, p01, p11,
                       sin_b, sin_w, sin_r, sin_p,
                       cos_a, cos_b, cos_w, cos_r, cos_p)

        J3sub2_vars = (rho, sin_a, sin_b, sin_w, sin_r,
                       cos_a, cos_b, cos_w, cos_r, cos_p)

        J3sub3_vars = (rho, sin_a, sin_b, sin_w, sin_r, sin_p,
                       cos_a, cos_b, cos_w, cos_r, cos_p)

        J3sub4_vars = ()
        J3sub5_vars = ()
        J3sub6_vars = ()
        J3sub7_vars = ()

        J3sub8_vars = (sin_a, sin_b, sin_w, sin_r, sin_p,
                       cos_a, cos_b, cos_w, cos_r, cos_p)

        fJ3sub0 = lambdify(J3sub0_vars, J3sub[0], eval_type)
        fJ3sub1 = lambdify(J3sub1_vars, J3sub[1], eval_type)
        fJ3sub2 = lambdify(J3sub2_vars, J3sub[2], eval_type)
        fJ3sub3 = lambdify(J3sub3_vars, J3sub[3], eval_type)
        fJ3sub4 = lambdify(J3sub4_vars, J3sub[4], eval_type)
        fJ3sub5 = lambdify(J3sub5_vars, J3sub[5], eval_type)
        fJ3sub6 = lambdify(J3sub6_vars, J3sub[6], eval_type)
        fJ3sub7 = lambdify(J3sub7_vars, J3sub[7], eval_type)
        fJ3sub8 = lambdify(J3sub8_vars, J3sub[8], eval_type)

        '''group the functioned jacobian components'''
        fJ1 = [fJ1sub0, fJ1sub1, fJ1sub2, fJ1sub3, fJ1sub4, fJ1sub5, fJ1sub6, fJ1sub7, fJ1sub8]
        fJ2 = [fJ2sub0, fJ2sub1, fJ2sub2, fJ2sub3, fJ2sub4, fJ2sub5, fJ2sub6, fJ2sub7, fJ2sub8]
        fJ3 = [fJ3sub0, fJ3sub1, fJ3sub2, fJ3sub3, fJ3sub4, fJ3sub5, fJ3sub6, fJ3sub7, fJ3sub8]
        
        return fJ1, fJ2, fJ3

    @staticmethod
    def calcRMSE(data):
        """calc root mean square error for input data"""
        num_coords, num_points = data.shape
        AMDE = np.mean(np.abs(data), axis=1)  # average mean distance error
        RMSE = sqrt(sum(sum(np.square(data))) / num_points)  # root mean squares error
        logging.info('Mean Difference:\n'
                     'X: {:.3f}\n'
                     'Y: {:.3f}\n'
                     'Z: {:.3f}'.format(AMDE[0], AMDE[1], AMDE[2]))
        logging.info('RMSE: {:.3f}\n'.format(RMSE))

    def calc_subaerial(self):
        itv = 10  # surface polynomial fitting interval

        if self.is_empty:  # i.e., if D is empty (las and sbet not merged)
            logging.warning('SBET and LAS not merged because max delta '
                            'time exceeded acceptable threshold of {} sec(s).'.format(Merge.dt_threshold))
        else:
            num_points = self.x_las.shape

            # assign variables because numexpr variables don't handle indexing
            x_las = self.x_las
            y_las = self.y_las
            z_las = self.z_las
            x_sbet = self.x_sbet
            y_sbet = self.y_sbet
            z_sbet = self.z_sbet
            x_ = ne.evaluate("x_las - x_sbet")
            y_ = ne.evaluate("y_las - y_sbet")
            z_ = ne.evaluate("z_las - z_sbet")

            fR0 = self.fR[0](self.h0, self.p0)
            fR3 = self.fR[3](self.h0, self.p0)
            fR6 = self.fR[6](self.p0)
            fR1 = self.fR[1](self.r0, self.p0, self.h0)
            fR4 = self.fR[4](self.r0, self.p0, self.h0)
            fR7 = self.fR[7](self.r0, self.p0)

            rho_est = ne.evaluate("sqrt(x_**2 + y_**2 + z_**2)")
            b_est = ne.evaluate("arcsin(((fR0 * x_) + (fR3 * y_) + (fR6 * z_)) / (-rho_est))")
            a_est = ne.evaluate("arcsin(((fR1 * x_) + (fR4 * y_) + (fR7 * z_)) / (rho_est * cos(b_est)))")
            w_est = np.zeros(num_points)

            # calc diff between true and est las xyz (LE = Laser Estimates)
            LE_pre_x = self.fF_orig[0](a_est, b_est, self.h0, self.p0, self.r0, rho_est, w_est, self.x_sbet)
            LE_pre_y = self.fF_orig[1](a_est, b_est, self.h0, self.p0, self.r0, rho_est, w_est, self.y_sbet)
            LE_pre_z = self.fF_orig[2](a_est, b_est, self.p0, self.r0, rho_est, w_est, self.z0)

            Er1x = ne.evaluate("x_las - LE_pre_x")
            Er1y = ne.evaluate("y_las - LE_pre_y")
            Er1z = ne.evaluate("z_las - LE_pre_z")

            # estimate error model using polynomial surface fitting
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

            '''
            Originally, Jaehoon used the Matlab "fit" function with a "poly23" option,
            but I didn't find a matching function to use in Python, so I manually coded
            the same functionally using np.linalg.lstsq with terms for
            a, b, a^2, ab, b^2, a^2b, ab^2, and b^3
            '''
            (coeffs_x, __, __, __) = np.linalg.lstsq(A, Er1x[::itv], rcond=None)  # rcond=None is FutureWarning
            (coeffs_y, __, __, __) = np.linalg.lstsq(A, Er1y[::itv], rcond=None)  # rcond=None is FutureWarning
            (coeffs_z, __, __, __) = np.linalg.lstsq(A, Er1z[::itv], rcond=None)  # rcond=None is FutureWarning

            p00x, p10x, p01x, p20x, p11x, p02x, p21x, p12x, p03x = coeffs_x
            p00y, p10y, p01y, p20y, p11y, p02y, p21y, p12y, p03y = coeffs_y
            p00z, p10z, p01z, p20z, p11z, p02z, p21z, p12z, p03z = coeffs_z

            #  calc ert2
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

            err_x = np.sum(A * coeffs_x, axis=1)
            err_y = np.sum(A * coeffs_y, axis=1)
            err_z = np.sum(A * coeffs_z, axis=1)

            LE_post_x = ne.evaluate('LE_pre_x + err_x')
            LE_post_y = ne.evaluate('LE_pre_y + err_y')
            LE_post_z = ne.evaluate('LE_pre_z + err_z')

            poly_surf_err_x = ne.evaluate('LE_post_x - x_las')
            poly_surf_err_y = ne.evaluate('LE_post_y - y_las')
            poly_surf_err_z = ne.evaluate('LE_post_z - z_las')

            # LE_post_arr = np.asarray([self.x_las, self.y_las, self.z_las])
            # LE_post1_arr = np.concatenate((self.LE_post1_arr, LE_post_arr), axis=1)

            '''
            simplify the numerous trigonometric calculations of the Jacobian by
            defining the Jacobian functions to be functions of the calculated sines and cosines
            of the various parameters, instead of the sin and cos functions directly;
            for example, instead of calculating sin(a_est) every time it shows up
            in an equation (which is a lot), calculate sin(a_est) once and use the value
            in the equation;  (it complicates the code, but it saves a noticeable
            chunk of time computationally)
            '''
            r0 = self.r0
            p0 = self.p0
            h0 = self.h0
            x0 = self.x0
            y0 = self.y0
            z0 = self.z0
            sin_a0 = ne.evaluate("sin(a_est)")
            sin_b0 = ne.evaluate("sin(b_est)")
            sin_w0 = np.zeros(num_points)  # ne.evaluate("sin(w_est)") ... sin(0) = 0
            sin_r0 = ne.evaluate("sin(r0)")
            sin_p0 = ne.evaluate("sin(p0)")
            sin_h0 = ne.evaluate("sin(h0)")
            cos_a0 = ne.evaluate("cos(a_est)")
            cos_b0 = ne.evaluate("cos(b_est)")
            cos_w0 = np.ones(num_points)  # ne.evaluate("cos(w_est)") ... cos(0) = 1
            cos_r0 = ne.evaluate("cos(r0)")
            cos_p0 = ne.evaluate("cos(p0)")
            cos_h0 = ne.evaluate("cos(h0)")

            pJ1 = np.vstack(
                (self.fJ1[0](a_est, b_est, rho_est, p10x, p21x, p12x, p11x, p20x,
                             sin_a0, sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0, cos_h0),
                 self.fJ1[1](a_est, b_est, rho_est, p02x, p03x, p21x, p12x, p01x, p11x,
                             sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0, cos_h0),
                 self.fJ1[2](rho_est, sin_a0, sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_h0),
                 self.fJ1[3](rho_est, sin_a0, sin_b0, sin_w0, sin_r0, sin_p0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0, cos_h0),
                 self.fJ1[4](rho_est, sin_a0, sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0, cos_h0),
                 np.ones(num_points),
                 self.fJ1[8](sin_a0, sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0, cos_h0)))

            pJ2 = np.vstack(
                (self.fJ2[0](a_est, b_est, rho_est, p10y, p21y, p12y, p11y, p20y,
                             sin_a0, sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0, cos_h0),
                 self.fJ2[1](a_est, b_est, rho_est, p02y, p03y, p21y, p12y, p01y, p11y,
                             sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0, cos_h0),
                 self.fJ2[2](rho_est, sin_a0, sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_h0),
                 self.fJ2[3](rho_est, sin_a0, sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0),
                 self.fJ2[4](rho_est, sin_a0, sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0, cos_h0),
                 np.ones(num_points),
                 self.fJ2[8](sin_a0, sin_b0, sin_w0, sin_r0, sin_p0, sin_h0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0, cos_h0)))

            pJ3 = np.vstack(
                (self.fJ3[0](a_est, b_est, rho_est, p10z, p21z, p12z, p11z, p20z,
                             sin_a0, sin_b0, sin_w0, sin_r0, sin_p0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0),
                 self.fJ3[1](a_est, b_est, rho_est, p02z, p03z, p21z, p12z, p01z, p11z,
                             sin_b0, sin_w0, sin_r0, sin_p0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0),
                 self.fJ3[2](rho_est, sin_a0, sin_b0, sin_w0, sin_r0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0),
                 self.fJ3[3](rho_est, sin_a0, sin_b0, sin_w0, sin_r0, sin_p0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0),
                 np.ones(num_points),
                 self.fJ3[8](sin_a0, sin_b0, sin_w0, sin_r0, sin_p0,
                             cos_a0, cos_b0, cos_w0, cos_r0, cos_p0)))

            # prop error
            C = self.C
            C = ne.evaluate("C * C")  # variance = stddev**2

            # delete the rows corresponding to the Jacobian functions that equal 0
            Cx = np.delete(C, [6, 7], 0)
            Cy = np.delete(C, [5, 7], 0)
            Cz = np.delete(C, [4, 5, 6], 0)

            sum_pJ1 = ne.evaluate("sum(pJ1 * pJ1 * Cx, axis=0)")
            sum_pJ2 = ne.evaluate("sum(pJ2 * pJ2 * Cy, axis=0)")
            sum_pJ3 = ne.evaluate("sum(pJ3 * pJ3 * Cz, axis=0)")

            sx = ne.evaluate("sqrt(sum_pJ1)")
            sy = ne.evaluate("sqrt(sum_pJ2)")
            subaerial_tvu = ne.evaluate("sqrt(sum_pJ3)")

            subaerial_thu = ne.evaluate('sqrt(sx**2 + sy**2)')

            column_headers = ['cblue_x', 'cblue_y', 'cblue_z',
                              'subaerial_thu', 'subaerial_tvu']

            return np.vstack((LE_post_x, LE_post_y, LE_post_z,
                              subaerial_thu, subaerial_tvu)).T, column_headers


if __name__ == '__main__':
    pass
