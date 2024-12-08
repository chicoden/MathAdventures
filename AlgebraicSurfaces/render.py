from Polynomial import analysis, solve, prettify
from Polynomial.arithmetic import *
from PIL import Image
import numpy as np
import sys

# Choices: "Genus 2 Surface", "Chmutov Surface", "Sarti Dodecic",
# "Heart", "Lampret", "Klein Bottle", "Blob Mix", "Bumpy Sphere",
# "y=Re(1/(1+(x+zi)^2))", "Trefoil"
surface_choice = sys.argv[1]

def normalize(v):
    return v / np.linalg.norm(v)

# Wrapper class for brevity
class Poly(list):
    def __init__(self, *coeffs):
        super().__init__(coeffs)

    # Operator overloading
    def __add__(self, other):
        return Poly(*poly_add(self, other))

    def __iadd__(self, other):
        return self + other

    def __sub__(self, other):
        return Poly(*poly_sub(self, other))

    def __isub__(self, other):
        return self - other

    def __mul__(self, other):
        return Poly(*poly_mul(self, other))

    def __imul__(self, other):
        return self * other

    def __pow__(self, other):
        return Poly(*poly_pow(self, other))

    def __ipow__(self, other):
        return self**other

    # Pretty printing for debugging
    def __str__(self):
        return prettify.stringify_poly(self)

    def __repr__(self):
        return str(self)

# Automatic differentiation (for calculating gradients)
class DualNumber3(object):
    def __init__(self, val, grad):
        self.val = val
        self.grad = np.array(grad)
        super().__init__()

    # Operator overloading
    def __add__(self, other):
        val = self.val + other.val
        grad_x = self.grad[0] + other.grad[0]
        grad_y = self.grad[1] + other.grad[1]
        grad_z = self.grad[2] + other.grad[2]
        return DualNumber3(val, [grad_x, grad_y, grad_z])

    def __iadd__(self, other):
        return self + other

    def __sub__(self, other):
        val = self.val - other.val
        grad_x = self.grad[0] - other.grad[0]
        grad_y = self.grad[1] - other.grad[1]
        grad_z = self.grad[2] - other.grad[2]
        return DualNumber3(val, [grad_x, grad_y, grad_z])

    def __isub__(self, other):
        return self - other

    def __mul__(self, other):
        val = self.val * other.val
        grad_x = self.val * other.grad[0] + self.grad[0] * other.val
        grad_y = self.val * other.grad[1] + self.grad[1] * other.val
        grad_z = self.val * other.grad[2] + self.grad[2] * other.val
        return DualNumber3(val, [grad_x, grad_y, grad_z])

    def __imul__(self, other):
        return self * other

    def __pow__(self, other):
        res = DualNumber3(1, [0, 0, 0])
        for i in range(0, other):
            res *= self

        return res

    def __ipow__(self, other):
        return self**other

    # Pretty printing for debugging
    def __str__(self):
        return "{0} + {1}dx + {2}dy + {3}dz".format(self.val, *self.grad)

    def __repr__(self):
        return str(self)

def surface(**kwargs):
    # Setup
    if len(kwargs) == 2: # Use case: compute intersection
        x = Poly(kwargs["ray_origin"][0], kwargs["ray_dir"][0])
        y = Poly(kwargs["ray_origin"][1], kwargs["ray_dir"][1])
        z = Poly(kwargs["ray_origin"][2], kwargs["ray_dir"][2])
        const = lambda c: Poly(c)

    else: # Use case: compute gradient
        x = DualNumber3(kwargs["pos"][0], [1, 0, 0])
        y = DualNumber3(kwargs["pos"][1], [0, 1, 0])
        z = DualNumber3(kwargs["pos"][2], [0, 0, 1])
        const = lambda c: DualNumber3(c, [0, 0, 0])

    # Surface definition
    if surface_choice == "Genus 2 Surface":
        q1 = const(2) * y
        q1 *= y**2 - const(3) * x**2
        q1 *= const(1) - z**2
        q2 = (x**2 + y**2)**2
        q3 = const(9) * z**2 - const(1)
        q3 *= const(1) - z**2
        return q1 + q2 - q3

    if surface_choice == "Chmutov Surface":
        Tx = const(32) * x**6 - const(48) * x**4 + const(18) * x**2 - const(1)
        Ty = const(32) * y**6 - const(48) * y**4 + const(18) * y**2 - const(1)
        Tz = const(32) * z**6 - const(48) * z**4 + const(18) * z**2 - const(1)
        return Tx + Ty + Tz

    if surface_choice == "Sarti Dodecic":
        # Scale 2x
        x *= const(0.5)
        y *= const(0.5)
        z *= const(0.5)

        # Parameter for the surface
        w = 1

        Q12 = (x**2 + y**2 + z**2 + const(w**2))**6

        l1 = x**4 + y**4 + z**4 + const(w**4)
        l2 = x**2 * y**2 + z**2 * const(w**2)
        l3 = x**2 * z**2 + y**2 * const(w**2)
        l4 = x**2 * const(w**2) + y**2 * z**2
        l5 = x * y * z * const(w)

        s10 = l1 * (l2 * l3 + l2 * l4 + l3 * l4)
        s11 = l1**2 * (l2 + l3 + l4)
        s12 = l1 * (l2**2 + l3**2 + l4**2)
        s51 = l5**2 * (l2 + l3 + l4)
        s234 = l2**3 + l3**3 + l4**3
        s23_p = l2**2 * l3 + l2 * l3**2
        s23_m = l2**2 * l3 - l2 * l3**2
        s34_p = l3**2 * l4 + l3 * l4**2
        s34_m = l3**2 * l4 - l3 * l4**2
        s42_p = l4**2 * l2 + l4 * l2**2
        s42_m = l4**2 * l2 - l4 * l2**2

        S12 = const(33 * np.sqrt(5)) * (s23_m + s34_m + s42_m)
        S12 += const(19) * (s23_p + s34_p + s42_p)
        S12 += const(10) * s234 - const(14) * s10 + const(2) * s11
        S12 -= const(6) * s12 + const(352) * s51
        S12 += const(336) * l5**2 * l1 + const(48) * l2 * l3 * l4

        return const(243) * S12 - const(22) * Q12

    if surface_choice == "Heart":
        q = x**2 + const(9 / 4) * z**2 + y**2 - const(1)
        return q**3 - (x**2 + const(9 / 80) * z**2) * y**3

    if surface_choice == "Lampret":
        q = (const(2.92) * (x - const(1)) * x**2 * (x + const(1)) + const(1.7) * y**2)**2 * (y**2 - const(0.88))**2
        q += (const(2.92) * (y - const(1)) * y**2 * (y + const(1)) + const(1.7) * z**2)**2 * (z**2 - const(0.88))**2
        q += (const(2.92) * (z - const(1)) * z**2 * (z + const(1)) + const(1.7) * x**2)**2 * (x**2 - const(0.88))**2
        return q - const(0.04)

    if surface_choice == "Klein Bottle":
        # Scale (4/9)x
        x *= const(2.25)
        y *= const(2.25)
        z *= const(2.25)

        shift = const(2) * y
        d = x**2 + y**2 + z**2 - const(1)
        d1, d2 = d + shift, d - shift
        return d1 * (d2**2 - const(8) * z**2) + const(16) * x * z * d2

    if surface_choice == "Blob Mix":
        # Scale (2/3)x
        x *= const(1.5)
        y *= const(1.5)
        z *= const(1.5)

        blobs = const(1)
        k = 2 / np.sqrt(3)
        for px in (-k, k):
            for py in (-k, k):
                for pz in (-k, k):
                    blobs *= (x - const(px))**2 + (y - const(py))**2 + (z - const(pz))**2 - const(1)

        return blobs + (x**2 + y**2 + z**2 - const(1) - blobs) * const(0.99996)

    if surface_choice == "Bumpy Sphere":
        # Parameters
        R = 1
        a = 0.1

        r = x**2 + z**2
        r5 = r**5
        p = const(10) * x**8 - const(120) * x**6 * z**2 + const(252) * x**4 * z**4 - const(120) * x**2 * z**6 + const(10) * z**8
        return (r + y**2) * r5**2 - (const(R) * r5 + const(a) * p * x * z)**2

    if surface_choice == "y=Re(1/(1+(x+zi)^2))":
        a = const(1) + x**2 - z**2
        b = const(2) * x * z
        return (a**2 + b**2) * y - a

    if surface_choice == "Trefoil":
        x *= const(1.125)
        y *= const(1.125)
        z *= const(1.125)
        r2 = x**2 + y**2 + z**2
        w = r2 - const(1)
        s = x**3 - const(3) * x * y**2 + z**2 - w**2
        t = y**3 - const(3) * x**2 * y - const(2) * z * w
        return s**2 + t**2 - (const(0.75) + const(0.12) * r2)**8

def get_intersection(ray_origin, ray_dir):
    # Checking for intersections with the sphere bounding the
    # surface to speed up rendering and to help the polynomial solver
    a = np.dot(ray_dir, ray_dir)
    b = np.dot(ray_origin, ray_dir)
    c = np.dot(ray_origin, ray_origin) - 4
    dis = b * b - a * c
    if dis > 0: # The ray is intersecting the bounding sphere
        ray_poly = surface(ray_origin=ray_origin, ray_dir=ray_dir)

        #lower_bound = analysis.lower_real_root_bound(ray_poly)
        #upper_bound = analysis.upper_real_root_bound(ray_poly)

        # Compute the intersections with the bounding sphere
        dis = np.sqrt(dis)
        lower_bound = (-b - dis) / a
        upper_bound = (-b + dis) / a

        # Compute the intersections with the surface
        roots = solve.solve_poly(ray_poly, iters=10, min_bound=lower_bound, max_bound=upper_bound)
        roots = [solve.newtons_method(ray_poly, root, 2, 1e-6) for root in roots]
        if len(roots) == 0: # No intersections
            return -1
    
        else:
            # Find the closest intersection that is in front of the camera
            t = -1
            for root in roots:
                if (True if t < 0 else root > 0 and root < t):
                    t = root

        return t

    else: # Didn't even hit the bounding sphere
        return -1

def get_normal(pos):
    grad = surface(pos=pos).grad
    return normalize(grad)

# Image to render to
img = Image.new("RGB", (800, 600), color=(0, 0, 0))
hw, hh = img.size[0] / 2, img.size[1] / 2

# Define the camera
cam_pos = np.array([5.7735, 3.7735, 4.7735])

# Camera coordinate frame
cam_fd = -normalize(cam_pos)
cam_rt = normalize(np.array([-cam_fd[2], 0, cam_fd[0]]))
cam_up = np.cross(cam_rt, cam_fd)

for y in range(0, img.size[1]):
    for x in range(0, img.size[0]):
        # Viewport coordinates
        uv = ((x - hw) / img.size[1], (hh - y) / img.size[1])

        # Compute the ray direction in the camera's coordinate frame
        ray_dir = normalize(cam_rt * uv[0] + cam_up * uv[1] + cam_fd)

        t = get_intersection(cam_pos, ray_dir)
        if t > 0:
            ray_hit = cam_pos + ray_dir * t
            surf_normal = get_normal(ray_hit)

            red = int((0.5 + 0.5 * surf_normal[0]) * 255)
            green = int((0.5 + 0.5 * surf_normal[1]) * 255)
            blue = int((0.5 + 0.5 * surf_normal[2]) * 255)

            img.putpixel((x, y), (red, green, blue))

    print("Row {0} of {1} rendered.".format(y + 1, img.size[1]))
    if y == hh: # Just in case something goes wrong, save progress
        img.save("render.png")

img.save("render.png")
