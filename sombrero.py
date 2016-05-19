#! /usr/bin/env python

"""
File: sombrero.py
Copyright (c) 2016 Austin Ayers
License: MIT


Imported from CW 11


Course: PHYS227
Assignment: Plots sombrero and motion of ball in said sombrero.
Date: April 26, 2016
Email: ayers111@mail.chapman.edu
Name: Austin Ayers
Description: Numerically solves Newton's equations for a ball being shaken in a sombrero by sinusoidal driving force.
The pinnacle of groupwork, efficiency,
"""
from math import sin, cos
from unittest import TestCase
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class Sombrero():
    """
    Sombrero equation plotter and solver
    Solves the coupled equations:
    x'(t) = y(t)
    y'(t) = -delta * y(t) + x(t) - x^3(t) + Fcos(omega * t)
    """
    def __init__(self, F, x_0, y_0, w = 1.0, delta = 0.25, m = 1):
        """
        Takes initial conditions F, x_0, and y_0
        """

        self.delta = delta
        self.F = F
        self.w = w
        self.m = m

        self.rk4_output = [[0],[x_0],[y_0]]

    def equation_1(self, x, y, t):
        return y
    def equation_2(self, x, y, t):
        return ( -1 * self.delta * y + x - x ** 3 + self.F * cos(self.w * t) ) / self.m
    
    def calc_n_rk4_steps(self, n, h):
        """
        Completely reworked from austin's implementation
        
        """
        f1 = self.equation_1
        f2 = self.equation_2
        vt = np.asarray([self.rk4_output[0][-1]+(h*i) for i in xrange(n+1)])
        vx = np.zeros(n+1) + self.rk4_output[1][-1]
        vy = np.zeros(n+1) + self.rk4_output[2][-1]

        for i in xrange(1, n+1):
            k1_x = h*f1(vx[i-1], vy[i-1], vt[i-1])
            k1_y = h*f2(vx[i-1], vy[i-1], vt[i-1])

            k2_x = h*f1(vx[i-1] + 0.5*k1_x, vy[i-1] + 0.5*k1_y, (vt[i-1] + vt[i]) / 2.0)
            k2_y = h*f2(vx[i-1] + 0.5*k1_x, vy[i-1] + 0.5*k1_y, (vt[i-1] + vt[i]) / 2.0)

            k3_x = h*f1(vx[i-1] + 0.5*k2_x, vy[i-1] + 0.5*k2_y,(vt[i-1] + vt[i]) / 2.0)
            k3_y = h*f2(vx[i-1] + 0.5*k2_x, vy[i-1] + 0.5*k2_y, (vt[i-1] + vt[i]) / 2.0)

            k4_x = h*f1(vx[i-1] + k3_x, vy[i-1] + k3_y, vt[i])
            k4_y = h*f2(vx[i-1] + k3_x, vy[i-1] + k3_y, vt[i])


            vx[i] = vx[i-1] + (k1_x + k2_x + k2_x + k3_x + k3_x + k4_x)/6
            vy[i] = vy[i-1] + (k1_y + k2_y + k2_y + k3_y + k3_y + k4_y)/6
        self.rk4_output[0] += vt[1:]
        self.rk4_output[1] += vx[1:]
        self.rk4_output[2] += vy[1:]
        return vt, vx, vy

    def rk4(self, f1, f2, x0, y0, n, h):
        """
        Runge Kutta 4 function

        Arguments: f1 (first coupled function being differentiated)
                   f2 (second coupled function being differentiated)
                   x0 (initial x point)
                   y0 (initial y point)
                   n (number of points to be calculated in rk4)
                   h (delta of x)

        Outputs: vx_1 (array of x points, subdivided n times)
                 vy_1 (array of y points, answer array which is output of rk4 - tied to vx_1)
                 vx_2 (array of x points, subdivided n times)
                 vy_2 (array of y points, answer array which is output of rk4 - tied to vx_1)
        """
        vx = [0]*(n + 1)
        vy = [0]*(n + 1)
        vx[0] = x = x0
        vy[0] = y = y0
        t = 0
        for i in range(1, n + 1):
            k1_x = h*f1(x, y, t)
            k1_y = h*f2(x, y, t)

            k2_x = h*f1(x + 0.5*k1_x, y + 0.5*k1_y, t + h / 2.0)
            k2_y = h*f2(x + 0.5*k1_x, y + 0.5*k1_y, t + h / 2.0)

            k3_x = h*f1(x + 0.5*k2_x, y + 0.5*k2_y, t + h / 2.0)
            k3_y = h*f2(x + 0.5*k2_x, y + 0.5*k2_y, t + h / 2.0)

            k4_x = h*f1(x + k3_x, y + k3_y, t + h)
            k4_y = h*f2(x + k3_x, y + k3_y, t + h)


            vx[i] = x = x + (k1_x + k2_x + k2_x + k3_x + k3_x + k4_x)/6
            vy[i] = y = y + (k1_y + k2_y + k2_y + k3_y + k3_y + k4_y)/6

            t += h

        return vx, vy

    def generate_parametric_graph(self):
        """
        Generates images ('.png's) of the plots of the Runge Kutta's
        solutions as well as the exact solutions.
        """

        approx_values = self.rk4_output
        u_list_approx = approx_values[1]
        v_list_approx = approx_values[2]

        fig, ax = plt.subplots(nrows = 1, ncols = 1)
        ax.plot(u_list_approx, v_list_approx, 'r')
        ax.set_xlabel("x(t)")
        ax.set_ylabel("y(t)")
        ax.set_title("Sombrero approximation (parametric) for n = " + str(self.n) + " and F = " + str(self.F))
        ax.grid(True)
        fig.savefig(self.__class__.__name__ + '_parametric'+ '_%3d' % (self.F) + '_%3d' % (self.x_0)  +  '_%3d' % (self.y_0) +  '_%3d.png' % (self.n))
        plt.close(fig)

    def generate_x_graph(self):
        """
        Generates images ('.png's) of the plots of the Runge Kutta's
        solutions as well as the exact solutions.
        """

        approx_values = self.rk4_output
        x_list_approx = approx_values[1]
        t_list = approx_values[0]

        fig, ax = plt.subplots(nrows = 1, ncols = 1)
        ax.plot(t_list, x_list_approx, 'r')
        ax.set_xlabel("t")
        ax.set_ylabel("x(t)")
        ax.set_title("Sombrero approximation for x(t) for n = " + str(self.n) + " and F = " + str(self.F))
        ax.grid(True)
        fig.savefig(self.__class__.__name__ + '_x' + '_%3d' % (self.F) + '_%3d' % (self.x_0)  +  '_%3d' % (self.y_0) +  '_%3d.png' % (self.n))
        plt.close(fig)

    def generate_poincare_graph(self, step, theta = 0 , xLim = (-1.5, 1.5), yLim = (-1.1, 1.1)):
        """
        Generates images ('.png's) of the plots of the Runge Kutta's
        solutions as well as the exact solutions.
        """
        approx_values = self.rk4_output
        u_list_approx = approx_values[1][int(theta*1000) :: step]
        v_list_approx = approx_values[2][int(theta*1000) :: step]

        #This method is more elegant, but rounding errors keep getting in the way
        #value_mask = np.where(approx_values[0] == step + theta)
        #u_list_approx = approx_values[1][value_mask]
        #v_list_approx = approx_values[2][value_mask]

        fig, ax = plt.subplots(nrows = 1, ncols = 1)
        ax.plot(u_list_approx[0::3], v_list_approx[0::3], 'r.')
        ax.plot(u_list_approx[1::3], v_list_approx[1::3], 'b.')
        ax.plot(u_list_approx[2::3], v_list_approx[2::3], 'g.')
        ax.set_xlabel("x")
        ax.set_ylabel("x'")
        ax.set_ylim(yLim)
        ax.set_xlim(xLim)
        ax.set_title("Sombrero approximation for theta = ." + str(theta))
        ax.grid(True)
        fig.savefig('frame%02d.png' % (100*theta))
        plt.close(fig)

class test_sombrero(TestCase):
    def test_linear(self):
        thisone = Sombrero(0, 0, 0)
        def line1(x, y, t):
            return 2
        def line2(x, y, t):
            return 5
        y1, y2 = thisone.rk4(line1, line2, 0, 0, 1000)
        apt = (math.fabs(y1[5] - 0.001*5*2) < 1e-5)  and ((y2[5] - 0.001*5*5) < 1e-5)
        msg = 'Linear test failed.'
        assert apt, msg 

