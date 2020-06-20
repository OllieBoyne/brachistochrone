from scipy.optimize import fsolve
from scipy.signal import resample
import numpy as np
from matplotlib import pyplot as plt
from viseng.figure import FuncFig

g = 9.81
G = 6.67408e-11  # Gravitational constant
M_earth = 5.972e24  # Mass of the earth
R_earth = 6371e3  # Radius of the earth

def linear_resample(y, x_orig, x_new):
	"""Given a data set (x_orig, y), resamples to produce new set (x_new, y_new)"""
	return np.interp(x_new, x_orig, y)

def compute_time_uniform(y, x, sum = True):
	"""***Numerical approach***
	Given a function y(x), compute the time for an object to travel from (0,0) to (L,H),
	along the surface y. Assumes uniform gravity = g

	Uses suvat equation v = u + at, summed across each pair of adjacent points (segments)
	v and u derived from y by conservation of energy
	a derived from angle of segment.

	If sum, return only total time. Otherwise, return array of t corresponding to x values
	"""

	dy = y[1:] - y[:-1]
	dx = x[1:] - x[:-1]
	L = np.sqrt(dy ** 2 + dx ** 2)  # segment length

	sinthetas = dy / L  # sin theta of each segment, measured from left horizontal
	a = g * sinthetas  # accelerations parallel to each segment

	u = np.sqrt(2 * g * y[:-1])  # initial speed parallel to each segment
	v = np.sqrt(2 * g * y[1:])  # final speed parallel to each segment

	t = (v - u) / a

	# if any segments slope upwards with insufficient initial speed, time is infinite
	if np.isnan(t).any():
		return np.inf

	if sum:	return t.sum()

	else: return np.cumsum(t)

def get_brach_times(L, H, x_range):
	"""Given the width L, and height H, of path A->B, returns the analytical solution for a brach time.
	Given parametrics x = a(p - sin p), y = a(1 - cos p),
	the travel time T can be found as sqrt(a/g) * phi, where phi = p @ (L, H)

	returns the t values corresponding to the input x values, using the relation p = sqrt(g/a) * t"""

	phi_condition = lambda phi: ((phi - np.sin(phi)) / (1-np.cos(phi))) - L/H
	phi = fsolve(phi_condition, x0 = np.pi)
	a = H / (1-np.cos(phi))

	p_findfunc = lambda p: ( a*(p-np.sin(p)) - x_range )
	p_range = fsolve(p_findfunc, np.linspace(0, phi, x_range.size)) # find p values corresponding to x values

	t = np.sqrt(a/g) * p_range
	return t

def get_brach_x_y(L, H, mu=0.0, n_samples = 100):
	"""Given the width L, and height H, of the path A -> B, returns a parametric (proportional to theta, or time)
	sampling of the brach curve. Also allows for non-zero friction, mu"""

	x_par = lambda A, theta: A * (theta - np.sin(theta) + mu * (1-np.cos(theta)))
	y_par = lambda A, theta: A * (1-np.cos(theta) + mu * (theta+np.sin(theta)))

	# Solve for phi (maximum theta) first. set k=1, as phi is independent of k
	phi_eqn = lambda t: (y_par(1, t) / x_par(1, t)) - (H / L)
	phi = fsolve(phi_eqn, x0=np.pi)[0]

	A_eqn = lambda A: y_par(A,phi) - H
	A = fsolve(A_eqn, x0=np.array([0.5]))[0]

	theta = np.linspace(0, phi, n_samples)

	return x_par(A, theta), y_par(A, theta)

def get_linear_x_y(L, H, n_samples = 100):
	"""Given the width, L, and height H, of the path A -> B, produces a sampling of a linear path,
	that corresponds to linear time"""
	t_max = np.sqrt(2 * L ** 2 * ( (H/L) ** 2 + 1) / (g * H))
	t = np.linspace(0, t_max, n_samples)
	x = t ** 2 * (g * H) / (2 * L * ((H/L)**2 + 1))
	y = x * H/L

	return x,y

def get_drop_x_y(L, H, n_samples = 100):
	"""Given the width, L, and height H, of the path A -> B, produces a sampling of a drop path,
	that corresponds to linear time.
	Drop path is pure drop until hits minimum y, then pure x movement"""
	t_drop = np.sqrt(2 * H / g) # time to drop
	v_drop = t_drop * g # velocity at bottom of drop
	t_roll = L / v_drop # time to roll (assuming perfect transfer of velocity)

	t = np.linspace(0, t_drop + t_roll, n_samples)
	x, y = np.zeros_like(t), np.zeros_like(t)

	x[t > t_drop] = v_drop * (t[t > t_drop] - t_drop)  # rolling x position
	y[t < t_drop] = 0.5 * g * t[t < t_drop] ** 2
	y[t > t_drop] = H

	return x, y

def get_quad_x_y(L, H, n_samples = 100):
	"""Given the width, L, and height H, of the path A -> B, produces a sampling of a quad path,
	that corresponds to linear time.
	quad path is y = H x ^ 2 / L ^ 2"""

	# solve numerically
	dt = 0
	pass


def get_nonuniform_brach_x_y(L, H, W, Z, R = None, M = None):
	"""Given a path that starts at (W, Z) from the centre of a planet, and ends at a width L and height H
	from there, find the minimum path function.
	R, M are Earth by default"""

	if R is None: R = R_earth
	if M is None: M = M_earth

	X = np.linspace(0, H, 10)

	def f(x, y, y_p):
		"""Return equation in integral of total time"""
		r = np.sqrt( (W - x) ** 2 + (Z - y) ** 2)
		r0 = np.sqrt( W ** 2 + Z ** 2 )

		return np.sqrt( (y_p**2 + 1) * r * r0 / (r - r0) )

	e_lagrange = lambda y: 0 #  euler lagrange equation to minimise


def produce_uniform_animation():
	"""Produce animation of brach, straight line, and straight drop, between (0, 0) and (L, H).
	Save as .gif"""

	anim_time = 5  # animation time in seconds
	play_speed = 5e-2#0.25 # speed factor relative to realtime
	fig, ax = plt.subplots(FigureClass=FuncFig, animate=True, time = anim_time)

	T = np.linspace(0, anim_time, fig.n_frames)  # all sampled times

	L = 2
	H = 1
	rad = max(L, H) / 50


	x_brach, y_brach = get_brach_x_y(L, H, n_samples = T.size)
	t_brach = get_brach_times(L, H, x_brach).max() / play_speed

	x_lin, y_lin = get_linear_x_y(L, H, n_samples = T.size)
	t_lin = np.sqrt( 2 * L ** 2 * ( (H/L)**2 + 1) / (g * H)) / play_speed # analytical solution for time

	# downwards drop
	x_drop, y_drop = get_drop_x_y(L, H, n_samples = T.size)
	t_drop = ( (2*H/g)**.5 + L / ((2*H/g)**.5*g) ) / play_speed

	# quadratic
	x_quad = np.linspace(0, L, 100)
	A = -1 # coefficient for determining graph shape
	y_quad = A * x_quad ** 2 + ((H-A*L**2)/L) * x_quad
	t_quad = compute_time_uniform(y_quad, x_quad, sum=False)
	t_quad_linear = np.linspace(0, t_quad[-1], T.size*10) # resample for linear time
	x_quad = linear_resample(x_quad[1:], t_quad, t_quad_linear)
	y_quad = linear_resample(y_quad[1:], t_quad, t_quad_linear)
	t_quad_total = t_quad[-1] / play_speed # total time for whole sequence

	def plot_fall(x, y, t_finish, color="blue"):
		### if missing origin, add to data
		x, y = np.concatenate([[0], x]), np.concatenate([[0], y])

		if color=="purple": print(x[0], y[0])

		ax.plot(x, y, color=color)
		ax.add_patch("Circle", xy=(x[0], y[0]), radius = rad, fill=False, ec=color, lw=4,
					 animate_kwargs=dict( xy = np.dstack([x, y])[0]),
					 timings = [0, t_finish])

	plot_fall(x_brach, y_brach, t_brach, "blue")
	plot_fall(x_lin, y_lin, t_lin, "red")
	plot_fall(x_drop, y_drop, t_drop, "orange")
	plot_fall(x_quad, y_quad, t_quad_total, "purple")

	ax.axis("equal")
	ax.axis("off")
	ax.invert_yaxis()
	fig.render()

	fig.extend()
	# fig.save(title="brach_comparison", fmt="gif")
	plt.show()

if __name__ == "__main__":

	# produce_uniform_animation()

	x, y = get_brach_x_y(3, 2.7, 1.0, 100)
	plt.plot(x, y)
	plt.show()

